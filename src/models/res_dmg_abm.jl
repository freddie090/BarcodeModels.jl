"""Resistance + damage population agent-based (ABM) model."""
struct ResDmg_ABM <: ABMModel
    params::ResDmgParams
    abm::ABMParams
    function ResDmg_ABM(params::ResDmgParams, abm::ABMParams)
        validate_model_params_strict(params)
        new(params, abm)
    end
end

"""Construct a `ResDmg_ABM` from parameters and ABM settings."""
ResDmg_ABM(params::ResDmgParams; abm::ABMParams = ABMParams()) = ResDmg_ABM(params, abm)

"""Construct a `ResDmg_ABM` from keyword arguments forwarded to `ResDmgParams`."""
ResDmg_ABM(; abm::ABMParams = ABMParams(), kwargs...) = ResDmg_ABM(ResDmgParams(; kwargs...), abm)

function _core_grow_kill_abm!(
    cells::Vector{ResDmgCell},
    params::ResDmgParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    @assert 0.0 <= sim.t_frac <= 1.0 "t_frac must be between 0 and 1."
    @assert params.drug_effect in RESDMG_DRUG_EFFECTS "drug effect can only take :b, :d, or :c."
    @assert sim.R_real in ["b", "d", "l"] "R_real can only take 'b', 'd' or 'l' as a value."

    independent = params.transition_coupling === :independent_rate
    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    lam = params.b - params.d
    drug_concs = nothing
    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
    end

    # Fixed full-run bounds cover every phenotype, in either coupling mode.
    bR_bound = sim.R_real == "b" ? params.b - lam * params.del :
        (sim.R_real == "l" ? params.b * (1 - params.del) : params.b)
    bmax = max(params.b, bR_bound, 0.0)
    dR_bound = sim.R_real == "d" ? params.d + lam * params.del : params.d
    drug_death_bound = treat && de in (:d, :c) ? params.Dc * max(1.0, 1 - params.psi) : 0.0
    dmax_fixed = max(max(params.d, dR_bound, 0.0) + drug_death_bound, params.d + params.Dc)
    zetmax = max(params.zet_S, params.zet_R)
    omemax = treat ? max(params.ome, params.ome * (1 - params.psi)) : 0.0
    transitionmax = independent ? max(params.mu, params.sig) : 0.0
    eventmax = bmax + dmax_fixed + zetmax + omemax + transitionmax

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_resdmg_cell)
    live_vec = live_positions(cells)
    dead_vec = dead_positions(cells)

    Nt = n_alive(cells)
    DScount = Int64(sum(cell -> cell.DS && cell.alive, cells))
    DRcount = Int64(sum(cell -> cell.DR && cell.alive, cells))
    Rcount = Int64(sum(cell -> cell.R && cell.alive, cells))
    Scount = Int64(sum(cell -> !cell.DS && !cell.DR && !cell.R && cell.alive, cells))

    phen_counts = ResDmgPhenoCounts(Scount, DScount, DRcount, Rcount)

    Nvec = Int64[Nt]
    Svec = Int64[phen_counts.Scount]
    DSvec = Int64[phen_counts.DScount]
    DRvec = Int64[phen_counts.DRcount]
    Rvec = Int64[phen_counts.Rcount]
    Pvec = Int64[sim.Passage]

    if eventmax == 0.0
        t = sim.tmax
        push!(tvec, t)
        push!(Nvec, Nt)
        push!(Svec, phen_counts.Scount)
        push!(DSvec, phen_counts.DScount)
        push!(DRvec, phen_counts.DRcount)
        push!(Rvec, phen_counts.Rcount)
        push!(Pvec, sim.Passage)
    end

    while t <= sim.tmax && eventmax > 0.0
        Nt > 0 || error("ABM core loop encountered non-positive live population (Nt=$(Nt)).")

        live_pos = rand(1:length(cells))

        while live_vec[live_pos] == 0
            live_pos = rand(1:length(cells))
        end

        dmax_window = dmax_fixed
        ran = rand(Uniform(0, eventmax))
        denom = eventmax * Nt
        denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (rate*Nt=$(denom)).")
        dt = -log(rand()) / denom
        t += dt

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, min(t, sim.tmax))
            push!(Svec, phen_counts.Scount)
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
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
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Pvec, sim.Passage)
            break
        end

        curr_rconc = treat ? curr_rc(t, drug_concs) : 0.0

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

        if cells[live_pos].DS || cells[live_pos].DR
            cell_b = 0.0
            cell_d = params.d + params.Dc
        elseif treat
            curr_dconc = curr_dc(t, drug_concs)
            if de == :d
                cell_d += if cells[live_pos].R
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :b
                cell_b -= if cells[live_pos].R
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :c
                cell_b -= if cells[live_pos].R
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

        density_factor = max(1 - Nt / sim.Cc, 0.0)
        cell_b *= density_factor
        cell_d *= density_factor

        if treat
            if cells[live_pos].R
                cell_ome = params.ome * (1 - params.psi) * curr_rconc
            elseif cells[live_pos].DS || cells[live_pos].DR
                cell_ome = 0.0
            else
                cell_ome = params.ome * curr_rconc
            end
        else
            cell_ome = 0.0
        end

        cell_ome <= omemax || error("Calculated cell_ome value $(cell_ome) exceeds maximum possible value $(omemax). Check parameter values and drug effect settings.")

        if cells[live_pos].DS
            cell_zet = params.zet_S
        elseif cells[live_pos].DR
            cell_zet = params.zet_R
        else
            cell_zet = 0.0
        end

        cell_zet <= zetmax || error("Calculated cell_zet value $(cell_zet) exceeds maximum possible value $(zetmax). Check parameter values.")

        if ran < cell_b
            any(!c.alive for c in cells) || error("ABM core loop has no available dead slots for an accepted birth; increase Nbuff.")
            dead_pos = rand(1:length(cells))
            while dead_vec[dead_pos] == 0
                dead_pos = rand(1:length(cells))
            end
            resdmg_birth_mutate_event!(cells, live_pos, dead_pos,
                                       independent ? 0.0 : params.mu, independent ? 0.0 : params.sig,
                                       phen_counts)

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

        if independent
            # Fixed, disjoint intervals permit at most one event per proposal.
            transition_start = bmax + dmax_window + zetmax + omemax
            if !cells[live_pos].DS && !cells[live_pos].DR
                if cells[live_pos].R
                    if transition_start <= ran < transition_start + params.sig
                        resdmg_R_to_S_transition_event!(cells, live_pos, phen_counts)
                    end
                elseif transition_start <= ran < transition_start + params.mu
                    resdmg_S_to_R_transition_event!(cells, live_pos, phen_counts)
                end
            end
        end

        if Nt <= 0
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Pvec, sim.Passage)
            break
        end
    end

    fin_t = round(t, digits = 4)
    return ResDmgGrowOut(Nvec, tvec, Svec, DSvec, DRvec, Rvec, Pvec, fin_t)
end

function run_model_core_abm(model::ResDmg_ABM, state::ResDmgABMState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm!(state.cells, model.params, sim; treat = treat)
end

function grow_kill_lin_kmc!(cells::Vector{ResDmgCell},
    b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64,
    ome::Float64, zet_S::Float64, zet_R::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64}, dt_save_at::Float64;
    R_real::String = "b", t_frac::Float64 = 0.050, treat::Bool = false,
    Passage::Int64 = 1,
    drug_effect::Union{String, Symbol} = :d,
    transition_coupling::Union{String, Symbol} = :birth_coupled)

    model = ResDmg_ABM(ResDmgParams(
        b = b, d = d, mu = mu, sig = sig, del = del,
        ome = ome, zet_S = zet_S, zet_R = zet_R,
        Dc = Dc, k = k, psi = psi, drug_effect = drug_effect,
        transition_coupling = transition_coupling
    ))

    sim = ABMSimParams(
        t0 = t0, tmax = tmax, Nmax = Nmax, Cc = Cc,
        treat_ons = treat_ons, treat_offs = treat_offs, dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, Passage = Passage, drug_effect = drug_effect
    )

    state = ResDmgABMState(cells)
    return run_model_core_abm(model, state, sim; treat = treat)
end
