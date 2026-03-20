"""Resistance + damage population agent-based (ABM) model."""
struct ResDmg_ABM <: ABMModel
    params::ResDmgParams
    abm::ABMParams
    function ResDmg_ABM(params::ResDmgParams, abm::ABMParams)
        validate_model_params(params)
        new(params, abm)
    end
end

"""Construct a `ResDmg_ABM` from parameters and ABM settings."""
ResDmg_ABM(params::ResDmgParams; abm::ABMParams = ABMParams()) = ResDmg_ABM(params, abm)

"""Construct a `ResDmg_ABM` from keyword arguments forwarded to `ResDmgParams`."""
ResDmg_ABM(; abm::ABMParams = ABMParams(), kwargs...) = ResDmg_ABM(ResDmgParams(; kwargs...), abm)

mutable struct ResDmgCell
    barcode::Float64
    D::Bool
    R::Bool
    E::Bool
    alive::Bool
end

mutable struct ResDmgPhenoCounts
    Scount::Int64
    Dcount::Int64
    Rcount::Int64
    Ecount::Int64
end

mutable struct ResDmgGrowOut
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    Svec::Vector{Int64}
    Dvec::Vector{Int64}
    Rvec::Vector{Int64}
    Evec::Vector{Int64}
    Pvec::Vector{Int64}
    fin_t::Float64
end

mutable struct ResDmgABMState
    cells::Vector{ResDmgCell}
end

make_dead_resdmg_cell() = ResDmgCell(0.0, false, false, false, false)

function seed_resdmg_cells(N::Int64, rho::Float64, Nbuff::Int64;
    skew_lib::Bool = false, bc_unif::Float64 = 0.0, Nbc::Int64 = 0)

    cells = Vector{ResDmgCell}(undef, max(N, Nbuff))

    if skew_lib
        bcs = collect(1:Nbc)
        bc_probs = generate_probabilities(bcs, bc_unif)
        samp_bcs = bcs[rand(Categorical(bc_probs), N)]
    else
        samp_bcs = collect(1:N)
    end

    for i in 1:N
        cells[i] = ResDmgCell(samp_bcs[i], false, false, false, true)
    end

    if rho > 0.0
        nR = Int64(round(rho * N))
        R_cells = sample(1:N, nR, replace = false)
        for i in 1:nR
            cells[R_cells[i]].R = true
        end
    end

    if Nbuff > N
        for i in (N + 1):Nbuff
            cells[i] = make_dead_resdmg_cell()
        end
    end

    return cells
end

function resdmg_birth_mutate_event!(cell_arr::Vector{ResDmgCell},
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].D = cell_arr[cell_pos].D
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].E = cell_arr[cell_pos].E
    cell_arr[birth_pos].alive = true

    mut_p = rand()

    if cell_arr[birth_pos].E
        phen_counts.Ecount += 1
    elseif cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        elseif sig <= mut_p < (sig + al)
            cell_arr[birth_pos].R = false
            cell_arr[birth_pos].E = true
            phen_counts.Ecount += 1
        else
            phen_counts.Rcount += 1
        end
    elseif cell_arr[birth_pos].D
        throw(ArgumentError("Cannot perform birth event attempt on a damaged cell."))
    else
        if mu > mut_p
            cell_arr[birth_pos].R = true
            phen_counts.Rcount += 1
        else
            phen_counts.Scount += 1
        end
    end
end

function resdmg_damage_event!(cell_arr::Vector{ResDmgCell},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    if !cell_arr[cell_pos].alive
        return
    end

    if cell_arr[cell_pos].D
        throw(ArgumentError("Cell is already damaged"))
    end

    if cell_arr[cell_pos].R
        cell_arr[cell_pos].R = false
        cell_arr[cell_pos].D = true
        phen_counts.Rcount -= 1
        phen_counts.Dcount += 1
    elseif cell_arr[cell_pos].E
        nothing
    else
        cell_arr[cell_pos].D = true
        phen_counts.Scount -= 1
        phen_counts.Dcount += 1
    end
end

function resdmg_repair_event!(cell_arr::Vector{ResDmgCell},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    if !cell_arr[cell_pos].alive
        return
    end

    if !cell_arr[cell_pos].R && !cell_arr[cell_pos].E && !cell_arr[cell_pos].D
        throw(ArgumentError("Cell is already sensitive"))
    end

    if cell_arr[cell_pos].R
        nothing
    elseif cell_arr[cell_pos].E
        nothing
    else
        cell_arr[cell_pos].D = false
        phen_counts.Dcount -= 1
        phen_counts.Scount += 1
    end
end

function resdmg_death_event!(cell_arr::Vector{ResDmgCell},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    cell_arr[cell_pos].alive = false

    if cell_arr[cell_pos].E
        phen_counts.Ecount -= 1
    elseif cell_arr[cell_pos].R
        phen_counts.Rcount -= 1
    elseif cell_arr[cell_pos].D
        phen_counts.Dcount -= 1
    else
        phen_counts.Scount -= 1
    end
end

function update_track_vec_resdmg!(kmc_out,
    Nvec::Vector{Int64},
    nS_vec::Vector{Int64}, nD_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64})

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Svec)
    append!(nD_vec, kmc_out.Dvec)
    append!(nR_vec, kmc_out.Rvec)
    append!(nE_vec, kmc_out.Evec)
    append!(tvec, kmc_out.tvec)
    append!(Pvec, kmc_out.Pvec)
end

function _core_grow_kill_abm!(
    cells::Vector{ResDmgCell},
    params::ResDmgParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    @assert 0.0 <= sim.t_frac <= 1.0 "t_frac must be between 0 and 1."
    @assert params.drug_effect in RESDMG_DRUG_EFFECTS "drug effect can only take :b, :d, or :c."
    @assert sim.R_real in ["b", "d", "l"] "R_real can only take 'b', 'd' or 'l' as a value."

    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    bmax = params.b
    omemax = params.ome
    zetmax = params.zet
    lam = bmax - params.d

    # Cost-aware no-treatment death envelope (maximum possible death rate without drug).
    dmax_no_treat = if params.del > 0.0 && sim.R_real == "d"
        params.d + (lam * params.del)
    else
        params.d
    end

    # Treatment-on death envelope adds maximum drug effect to no-treatment envelope.
    dmax_treat = dmax_no_treat + params.Dc

    drug_concs = nothing
    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
        bdmax = bmax + dmax_treat
    else
        bdmax = bmax + dmax_no_treat
    end
    bdzetmax = bdmax + zetmax
    bdzetomemax = bdmax + omemax + zetmax

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_resdmg_cell)
    live_vec = live_positions(cells)
    dead_vec = dead_positions(cells)

    Nt = n_alive(cells)
    Dcount = Int64(sum(cell -> cell.D && cell.alive, cells))
    Rcount = Int64(sum(cell -> cell.R && cell.alive, cells))
    Ecount = Int64(sum(cell -> cell.E && cell.alive, cells))
    Scount = Int64(sum(cell -> !cell.D && !cell.R && !cell.E && cell.alive, cells))

    phen_counts = ResDmgPhenoCounts(Scount, Dcount, Rcount, Ecount)

    Nvec = Int64[Nt]
    Svec = Int64[phen_counts.Scount]
    Dvec = Int64[phen_counts.Dcount]
    Rvec = Int64[phen_counts.Rcount]
    Evec = Int64[phen_counts.Ecount]
    Pvec = Int64[sim.Passage]

    while t <= sim.tmax
        Nt > 0 || error("ABM core loop encountered non-positive live population (Nt=$(Nt)).")
        any(!c.alive for c in cells) || error("ABM core loop has no available dead slots for births; increase Nbuff.")

        live_pos = rand(1:length(cells))
        dead_pos = rand(1:length(cells))

        while live_vec[live_pos] == 0
            live_pos = rand(1:length(cells))
        end

        while dead_vec[dead_pos] == 0
            dead_pos = rand(1:length(cells))
        end

        curr_rconc = 0.0
        if treat
            curr_rconc = curr_rc(t, drug_concs)
            dmax_window = curr_rconc > 0.0 ? dmax_treat : dmax_no_treat
            ran_max = curr_rconc > 0.0 ? bdzetomemax : bdzetmax
            ran = rand(Uniform(0, ran_max))
            denom = ran_max * Nt
            denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (rate*Nt=$(denom)).")
            dt = -log(rand()) / denom
            t += dt
        else
            dmax_window = dmax_no_treat
            ran = rand(Uniform(0, bdzetmax))
            denom = bdzetmax * Nt
            denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (rate*Nt=$(denom)).")
            dt = -log(rand()) / denom
            t += dt
        end

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Dvec, phen_counts.Dcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            t_rec += t_rec_change
        end

        if t >= sim.tmax
            t = sim.tmax
            break
        end

        if Nt >= sim.Nmax
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Dvec, phen_counts.Dcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            break
        end

        if params.del > 0.0
            if cells[live_pos].R
                cell_b, cell_d = if sim.R_real == "b"
                    (params.b - (lam * params.del), params.d)
                elseif sim.R_real == "d"
                    (params.b, params.d + (lam * params.del))
                else
                    (params.b * (1 - params.del), params.d * (1 - params.del))
                end
            else
                cell_b = params.b
                cell_d = params.d
            end
            cell_b = max(cell_b, 0.0)
            cell_d = max(cell_d, 0.0)
        else
            cell_b = params.b
            cell_d = params.d
        end

        if cells[live_pos].D
            cell_b = 0.0
            cell_d = params.d + (treat ? params.Dc : 0.0)
        elseif treat
            curr_dconc = curr_dc(t, drug_concs)
            if de == :d
                cell_d += if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :b
                cell_b -= if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :c
                cell_b -= if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
                if cell_b < 0.0
                    negative_surplus = abs(cell_b)
                    cell_b = 0.0
                    cell_d += negative_surplus
                end
            end
        end

        cell_b *= (1 - (Nt / sim.Cc))
        cell_d *= (1 - (Nt / sim.Cc))

        if treat
            if cells[live_pos].R
                cell_ome = params.ome * (1 - params.psi) * curr_rconc
            elseif cells[live_pos].D || cells[live_pos].E
                cell_ome = 0.0
            else
                cell_ome = params.ome * curr_rconc
            end
        else
            cell_ome = 0.0
        end

        cell_ome <= omemax || error("Calculated cell_ome value $(cell_ome) exceeds maximum possible value $(omemax). Check parameter values and drug effect settings.")

        if cells[live_pos].D
            cell_zet = params.zet
        else
            cell_zet = 0.0
        end

        cell_zet <= zetmax || error("Calculated cell_zet value $(cell_zet) exceeds maximum possible value $(zetmax). Check parameter values.")

        if ran < cell_b
            if treat
                al_scal = params.al * curr_rconc
                resdmg_birth_mutate_event!(cells, live_pos, dead_pos,
                                           params.mu, params.sig, al_scal,
                                           phen_counts)
            else
                resdmg_birth_mutate_event!(cells, live_pos, dead_pos,
                                           params.mu, params.sig, 0.0,
                                           phen_counts)
            end

            live_vec[dead_pos] = 1
            dead_vec[dead_pos] = 0
            Nt += 1
        end

        if bmax <= ran < bmax + cell_d
            resdmg_death_event!(cells, live_pos, phen_counts)
            live_vec[live_pos] = 0
            dead_vec[live_pos] = 1
            Nt -= 1
        end

        if (bmax + dmax_window) <= ran < (bmax + dmax_window + cell_zet)
            resdmg_repair_event!(cells, live_pos, phen_counts)
        end

        if treat && ((bmax + dmax_window + zetmax) <= ran < (bmax + dmax_window + zetmax + cell_ome))
            resdmg_damage_event!(cells, live_pos, phen_counts)
        end

        if Nt <= 0
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Dvec, phen_counts.Dcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            break
        end
    end

    fin_t = round(t, digits = 4)
    return ResDmgGrowOut(Nvec, tvec, Svec, Dvec, Rvec, Evec, Pvec, fin_t)
end

function run_model_core_abm(model::ResDmg_ABM, state::ResDmgABMState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm!(state.cells, model.params, sim; treat = treat)
end

function grow_kill_lin_kmc!(cells::Vector{ResDmgCell},
    b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64, al::Float64,
    ome::Float64, zet::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64}, dt_save_at::Float64;
    R_real::String = "b", t_frac::Float64 = 0.050, treat::Bool = false,
    Passage::Int64 = 1,
    drug_effect::Union{String, Symbol} = :d)

    model = ResDmg_ABM(ResDmgParams(
        b = b, d = d, mu = mu, sig = sig, del = del, al = al,
        ome = ome, zet = zet,
        Dc = Dc, k = k, psi = psi, drug_effect = drug_effect
    ))

    sim = ABMSimParams(
        t0 = t0, tmax = tmax, Nmax = Nmax, Cc = Cc,
        treat_ons = treat_ons, treat_offs = treat_offs, dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, Passage = Passage, drug_effect = drug_effect
    )

    state = ResDmgABMState(cells)
    return run_model_core_abm(model, state, sim; treat = treat)
end
