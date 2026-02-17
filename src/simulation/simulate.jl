function _kw(kwargs, name, default)
    return haskey(kwargs, name) ? kwargs[name] : default
end

function simulate(model::AbstractBarcodeModel, exp::ExperimentParams; kwargs...)
    de = _kw(kwargs, :drug_effect, model.params.drug_effect)
    return simulate_expansion_and_treatment_hybrid(
        exp.n0, model.params.b, model.params.d,
        model.params.rho, model.params.mu, model.params.sig, model.params.del, model.params.al,
        model.params.Dc, model.params.k, model.params.psi,
        exp.t_exp, exp.tmax, exp.t_Pass,
        exp.Nseed, exp.Nmax, exp.Cc,
        exp.treat_ons, exp.treat_offs,
        exp.t_keep, exp.Nswitch;
        save_at = _kw(kwargs, :save_at, exp.save_at),
        n_Pass = _kw(kwargs, :n_Pass, exp.n_Pass),
        n_rep = _kw(kwargs, :n_rep, exp.n_rep),
        drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment),
        drug_effect = de,
        full_sol = _kw(kwargs, :full_sol, exp.full_sol),
        run_IC = _kw(kwargs, :run_IC, exp.run_IC),
        IC_n0 = _kw(kwargs, :IC_n0, exp.IC_n0),
        IC_tmax = _kw(kwargs, :IC_tmax, exp.IC_tmax),
        IC_treat_on = _kw(kwargs, :IC_treat_on, exp.IC_treat_on),
        run_colony = _kw(kwargs, :run_colony, exp.run_colony),
        nCol = _kw(kwargs, :nCol, exp.nCol),
        tCol = _kw(kwargs, :tCol, exp.tCol),
        ColNmax = _kw(kwargs, :ColNmax, exp.ColNmax)
    )
end

simulate_expansion_and_treatment_hybrid(model::AbstractBarcodeModel, exp::ExperimentParams; kwargs...) =
    simulate(model, exp; kwargs...)

"""
simulate_expansion_and_treatment_hybrid(
    n0, b, d, rho, mu, sig, del, al,
    Dc, k, psi, t_exp, tmax, t_Pass,
    Nseed, Nmax, Cc, treat_ons, treat_offs,
    t_keep, Nswitch;
    save_at=0.5, n_Pass=1, n_rep=4,
    drug_effect="d", full_sol=false
) -> Dict

Simulates clonal expansion and pulsed drug treatment using a hybrid model combining ODEs and stochastic jump processes.

Each replicate begins with untreated expansion, followed by splitting into replicates exposed to drug pulses. Phenotype switching and logistic growth are included. Population size is updated deterministically or stochastically depending on compartment size.
"""
function simulate_expansion_and_treatment_hybrid(n0::Int64, b::Float64, d::Float64,
    rho::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t_exp::Union{Float64, Vector{Float64}},
    tmax::Float64,
    t_Pass::Union{Float64, Vector{Float64}},
    Nseed::Union{Int64, Vector{Int64}},
    Nmax::Int64, Cc::Int64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64},
    t_keep::Vector{Float64}, Nswitch::Int64;
    save_at::Float64 = 0.5,
    n_Pass::Int64 = 1, n_rep::Int64 = 4,
    drug_treatment::Bool = true,
    drug_effect::Union{String, Symbol} = "d",
    full_sol::Bool = false,
    run_IC::Bool = false,
    IC_n0::Int64 = 1000, IC_tmax::Float64 = 4.0, IC_treat_on::Float64 = 1.0,
    run_colony::Bool = false,
    nCol::Int64 = 1000, tCol::Float64 = 12.0, ColNmax::Int64 = 50)

    @assert !(run_IC && run_colony) "Cannot run IC and colony assays at the same time."

    n_pop_obsv = n_rep * ((run_IC * 2) + length(t_keep) + n_Pass) + (run_colony * 1)
    de = normalize_drug_effect(drug_effect)

    if (al + sig) > 1.0 || psi < 0.0
        fin_t_outs = repeat([-1.0], n_pop_obsv)
        fin_u_outs = repeat([-1.0], n_pop_obsv)
    else
        exp_cells = nothing

        if (t_exp isa AbstractVector) && (Nseed isa AbstractVector)
            @assert length(t_exp) == length(Nseed) "Length of t_exp and Nseed vectors must match."

            nR = Float64(round(rho * n0))
            nS = n0 - nR
            nE = 0.0

            for i in 1:(length(t_exp) - 1)
                N_total = Int64(round(nS + nR + nE))
                exp_cells = simulate_grow_kill(N_total, nS, nR, nE, b, d, mu, sig, del, al,
                                          Dc, k, psi, 0.0, t_exp[i], -1.0, Nmax, Cc,
                                          treat_ons, treat_offs, Nswitch,
                                          save_at = save_at, treat = false, n_Pass = 1,
                                          drug_effect = de)

                nS = Int64(round(pop_fun(last(exp_cells.u))[1]))
                nR = Int64(round(pop_fun(last(exp_cells.u))[2]))
                nE = Int64(round(pop_fun(last(exp_cells.u))[3]))

                total_cells = nS + nR + nE
                if total_cells < Nseed[i]
                    return Dict("t" => repeat([-1.0], n_pop_obsv),
                                "u" => repeat([-1.0], n_pop_obsv))
                end

                counts = multivariate_hypergeometric_draw([nS, nR, nE], Nseed[i])

                nS = Float64(counts[1])
                nR = Float64(counts[2])
                nE = Float64(counts[3])
            end

            N_total = Int64(round(nS + nR + nE))
            exp_cells = simulate_grow_kill(N_total, nS, nR, nE, b, d, mu, sig, del, al,
                                      Dc, k, psi, 0.0, t_exp[end], -1.0, Nmax, Cc,
                                      treat_ons, treat_offs, Nswitch,
                                      save_at = save_at, treat = false, n_Pass = 1,
                                      drug_effect = de)

            exp_nS = Int64(round(pop_fun(last(exp_cells.u))[1]))
            exp_nR = Int64(round(pop_fun(last(exp_cells.u))[2]))
            exp_nE = Int64(round(pop_fun(last(exp_cells.u))[3]))

        elseif (t_exp isa Real) && (Nseed isa Integer)
            nR = Float64(round(rho * n0))
            nS = n0 - nR
            nE = 0.0

            exp_cells = simulate_grow_kill(n0, nS, nR, 0.0,
                                      b, d, mu, sig, del, al, Dc, k, psi,
                                      0.0, t_exp, -1.0, Nmax, Cc,
                                      treat_ons, treat_offs, Nswitch,
                                      save_at = save_at, treat = false, n_Pass = 1,
                                      drug_effect = de)

            exp_nS = Int64(round(pop_fun(last(exp_cells.u))[1]))
            exp_nR = Int64(round(pop_fun(last(exp_cells.u))[2]))
            exp_nE = Int64(round(pop_fun(last(exp_cells.u))[3]))
        else
            error("t_exp and Nseed must either both be scalars or both be vectors of the same length.")
        end

        nseed_last = Nseed isa AbstractVector ? last(Nseed) : Nseed

        if (exp_nS + exp_nR + exp_nE) < (nseed_last * n_rep)
            fin_t_outs = repeat([-1.0], n_pop_obsv)
            fin_u_outs = repeat([-1.0], n_pop_obsv)
            return Dict("t" => fin_t_outs, "u" => fin_u_outs)
        else
            rep_phenos = Array{Float64}(undef, 3, n_rep)

            nS = exp_nS
            nR = exp_nR
            nE = exp_nE

            for i in 1:n_rep
                sampled_counts = multivariate_hypergeometric_draw([nS, nR, nE], nseed_last)

                rep_phenos[1, i] = Float64(sampled_counts[1])
                rep_phenos[2, i] = Float64(sampled_counts[2])
                rep_phenos[3, i] = Float64(sampled_counts[3])

                nS -= sampled_counts[1]
                nR -= sampled_counts[2]
                nE -= sampled_counts[3]
            end

            t0 = 0.0

            fin_t_outs = Float64[]
            fin_u_outs = Float64[]

            if full_sol == true
                sol_dfs = DataFrame[]
            end

            if run_colony
                tx1_cols = 0
                tx0_cols = 0

                col_nR = max(Float64(round(rho * n0)), 0.0)
                col_nS = max((n0 - col_nR), 0.0)
                col_nE = 0.0

                for j in 1:nCol
                    total_cells = max(col_nS + col_nR + col_nE, 1e-10)
                    cell_probs = [col_nS, col_nR, col_nE] ./ total_cells
                    if any(cell_probs .< 0) || sum(cell_probs) == 0
                        error("Invalid probability vector: $(cell_probs)")
                    end
                    sampled_phenos = zeros(Float64, 3)
                    sampled_phenos[rand(Categorical(cell_probs))] = 1.0

                    col_sol_tx1 = simulate_grow_kill(1,
                                                sampled_phenos[1],
                                                sampled_phenos[2],
                                                sampled_phenos[3],
                                                b, d, mu, sig, del, al, Dc, k, psi, 0.0,
                                                tCol, -1.0, ColNmax + 10, Cc,
                                                [1.0], [100.0],
                                                Nswitch, save_at = save_at, treat = true,
                                                n_Pass = 1,
                                                drug_effect = de)

                    col_sol_tx0 = simulate_grow_kill(1,
                                                sampled_phenos[1],
                                                sampled_phenos[2],
                                                sampled_phenos[3],
                                                b, d, mu, sig, del, al, Dc, k, psi, 0.0,
                                                tCol, -1.0, ColNmax + 10, Cc,
                                                [1.0], [100.0],
                                                Nswitch, save_at = save_at, treat = false,
                                                n_Pass = 1,
                                                drug_effect = de)

                    if round(sum(last(pop_fun.(col_sol_tx1.u)))) >= ColNmax
                        tx1_cols += 1
                    end
                    if round(sum(last(pop_fun.(col_sol_tx0.u)))) >= ColNmax
                        tx0_cols += 1
                    end
                end

                col_prop = tx1_cols / tx0_cols

                append!(fin_t_outs, tCol)
                append!(fin_u_outs, col_prop)
            end

            for i in 1:n_rep
                if run_IC == true
                    nR = Float64(round(rho * IC_n0))
                    nS = IC_n0 - nR
                    nE = 0.0

                    IC_sol_tx1 = simulate_grow_kill(IC_n0,
                                               nS,
                                               nR,
                                               nE,
                                               b, d, mu, sig, del, al, Dc, k, psi, t0,
                                               IC_tmax, -1.0, Nmax, Cc,
                                               [IC_treat_on], [100.0],
                                               Nswitch, save_at = save_at,
                                               treat = true,
                                               n_Pass = 1,
                                               drug_effect = de)

                    IC_sol_tx0 = simulate_grow_kill(IC_n0,
                                               nS,
                                               nR,
                                               nE,
                                               b, d, mu, sig, del, al, Dc, k, psi, t0,
                                               IC_tmax, -1.0, Nmax, Cc,
                                               [IC_treat_on], [100.0],
                                               Nswitch, save_at = save_at,
                                               treat = false,
                                               n_Pass = 1,
                                               drug_effect = de)
                end

                sol = simulate_grow_kill(nseed_last,
                                    rep_phenos[1, i],
                                    rep_phenos[2, i],
                                    rep_phenos[3, i],
                                    b, d, mu, sig, del, al, Dc, k, psi, t0,
                                    tmax, t_Pass, Nmax, Cc,
                                    treat_ons, treat_offs,
                                    Nswitch, save_at = save_at,
                                    treat = drug_treatment,
                                    n_Pass = n_Pass,
                                    drug_effect = de)

                if full_sol == true
                    sol_df = DataFrame(
                        t = sol.t,
                        nS = map(x -> x[NS_INDEX], sol.u),
                        nR = map(x -> x[NR_INDEX], sol.u),
                        nE = map(x -> x[NE_INDEX], sol.u),
                        N = map(x -> x[NS_INDEX], sol.u) .+
                            map(x -> x[NR_INDEX], sol.u) .+
                            map(x -> x[NE_INDEX], sol.u),
                        rep = i
                    )
                    push!(sol_dfs, sol_df)
                end

                t_outs = Vector{Float64}(undef, (run_IC * 2) + length(t_keep) + n_Pass)
                u_outs = Vector{Float64}(undef, (run_IC * 2) + length(t_keep) + n_Pass)

                if run_IC
                    t_outs[1] = last(IC_sol_tx1.t)
                    u_outs[1] = Float64(round(sum(last(pop_fun.(IC_sol_tx1.u)))))
                    t_outs[2] = last(IC_sol_tx0.t)
                    u_outs[2] = Float64(round(sum(last(pop_fun.(IC_sol_tx0.u)))))
                end

                for j in 1:length(t_keep)
                    if !(t_keep[j] in sol.t)
                        t_closest_pos = findmin(abs.(sol.t .- t_keep[j]))[2]
                        t_closest = sol.t[t_closest_pos]
                        t_outs[j + (run_IC * 2)] = t_closest
                        u_outs[j + (run_IC * 2)] = Float64(round(sum.(pop_fun.(sol.u))[t_closest_pos]))
                    else
                        t_realised_pos = findlast(sol.t .== t_keep[j])
                        t_realised = sol.t[t_realised_pos]
                        t_outs[j + (run_IC * 2)] = t_realised
                        u_outs[j + (run_IC * 2)] = Float64(round(sum.(pop_fun.(sol.u))[t_realised_pos]))
                    end
                end

                for j in 1:n_Pass
                    if sum(pass_fun.(sol.u) .== j) > 0
                        t_realised_pos = findlast(pass_fun.(sol.u) .== j)
                        t_outs[length(t_keep) + j + (run_IC * 2)] = sol.t[t_realised_pos]
                        u_outs[length(t_keep) + j + (run_IC * 2)] = Float64(round(sum.(pop_fun.(sol.u))[t_realised_pos]))
                    else
                        t_outs[length(t_keep) + j + (run_IC * 2)] = t_outs[length(t_keep) + j + ((run_IC * 2) - 1)]
                        u_outs[length(t_keep) + j + (run_IC * 2)] = u_outs[length(t_keep) + j + ((run_IC * 2) - 1)]
                    end
                end

                append!(fin_t_outs, t_outs)
                append!(fin_u_outs, u_outs)
            end

            fin_t_outs = Float64.(round.(fin_t_outs, digits = 0))
        end
    end

    if full_sol == true
        sol_df = vcat(sol_dfs...)
        return Dict("t" => fin_t_outs, "u" => fin_u_outs, "sol_df" => sol_df)
    else
        return Dict("t" => fin_t_outs, "u" => fin_u_outs)
    end
end

function model_meas_noise(u_outs::Array{Float64}, phi::Float64,
                          run_colony::Bool = false)

    if sum(u_outs .< 0) > 0
        noise_u_outs = u_outs
    else
        if run_colony
            noise_u_outs = rand.(Normal.(u_outs[2:length(u_outs)],
                                         u_outs[2:length(u_outs)] * phi))
        else
            noise_u_outs = rand.(Normal.(u_outs, u_outs * phi))
        end

        noise_u_outs[noise_u_outs .< 0] .= 0
        noise_u_outs = Int64.(round.(noise_u_outs, digits = -3))

        if run_colony
            insert!(noise_u_outs, 1, u_outs[1])
        end
    end

    return noise_u_outs
end
