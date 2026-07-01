"""Resistance population in vivo ABM model with engraftment phenotype parameters."""
struct ResPopInVivo_ABM <: ABMModel
    params::ResPopInVivoParams
    abm::ABMParams
    function ResPopInVivo_ABM(params::ResPopInVivoParams, abm::ABMParams)
        validate_model_params_strict(params)
        new(params, abm)
    end
end

ResPopInVivo_ABM(params::ResPopInVivoParams; abm::ABMParams = ABMParams()) = ResPopInVivo_ABM(params, abm)
ResPopInVivo_ABM(; abm::ABMParams = ABMParams(), kwargs...) = ResPopInVivo_ABM(ResPopInVivoParams(; kwargs...), abm)

mutable struct InVivoCancerCell
    barcode::Float64
    R::Bool
    E::Bool
    EG::Bool
    alive::Bool
end

mutable struct ResPopInVivoABMState
    cells::Vector{InVivoCancerCell}
end

make_dead_cell_invivo() = InVivoCancerCell(0.0, false, false, false, false)

mutable struct ResPopInVivoGrowOut
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    Svec::Vector{Int64}
    Rvec::Vector{Int64}
    Evec::Vector{Int64}
    Pvec::Vector{Int64}
    fin_t::Float64
    S_EG0_vec::Vector{Int64}
    S_EG1_vec::Vector{Int64}
    R_EG0_vec::Vector{Int64}
    R_EG1_vec::Vector{Int64}
    E_EG0_vec::Vector{Int64}
    E_EG1_vec::Vector{Int64}
end

function _invivo_abm_pheno_eg_counts(cells::Vector{InVivoCancerCell})
    nS_EG0 = 0
    nS_EG1 = 0
    nR_EG0 = 0
    nR_EG1 = 0
    nE_EG0 = 0
    nE_EG1 = 0

    for cell in cells
        cell.alive || continue
        if cell.E
            if cell.EG
                nE_EG1 += 1
            else
                nE_EG0 += 1
            end
        elseif cell.R
            if cell.EG
                nR_EG1 += 1
            else
                nR_EG0 += 1
            end
        else
            if cell.EG
                nS_EG1 += 1
            else
                nS_EG0 += 1
            end
        end
    end

    return nS_EG0, nS_EG1, nR_EG0, nR_EG1, nE_EG0, nE_EG1
end

function _push_invivo_abm_pheno_eg_counts!(S_EG0_vec, S_EG1_vec, R_EG0_vec, R_EG1_vec, E_EG0_vec, E_EG1_vec, cells)
    nS_EG0, nS_EG1, nR_EG0, nR_EG1, nE_EG0, nE_EG1 = _invivo_abm_pheno_eg_counts(cells)
    push!(S_EG0_vec, nS_EG0)
    push!(S_EG1_vec, nS_EG1)
    push!(R_EG0_vec, nR_EG0)
    push!(R_EG1_vec, nR_EG1)
    push!(E_EG0_vec, nE_EG0)
    push!(E_EG1_vec, nE_EG1)
    nothing
end

function seed_invivo_cells(N::Int64, rho::Float64, fEG1::Float64, Nbuff::Int64;
    skew_lib::Bool = false, use_lib_probs::Bool = false, bc_unif::Float64 = 0.0,
    Nbc::Int64 = 0, bc_probs = Float64[])

    cells = Vector{InVivoCancerCell}(undef, max(N, Nbuff))

    samp_bcs = sample_initial_barcodes(N;
                                       skew_lib = skew_lib,
                                       use_lib_probs = use_lib_probs,
                                       bc_unif = bc_unif,
                                       Nbc = Nbc,
                                       bc_probs = bc_probs)

    nEG1 = Int64(round(N * fEG1))
    nEG1 = clamp(nEG1, 0, N)
    eg1_indices = nEG1 == 0 ? Int64[] : sample(1:N, nEG1, replace = false)
    eg1_mask = falses(N)
    for idx in eg1_indices
        eg1_mask[idx] = true
    end

    for i in 1:N
        cells[i] = InVivoCancerCell(samp_bcs[i], false, false, eg1_mask[i], true)
    end

    if rho > 0.0
        nR = Int64(round(rho * N))
        if nR > 0
            R_cells = sample(1:N, nR, replace = false)
            for i in 1:nR
                cells[R_cells[i]].R = true
            end
        end
    end

    if Nbuff > N
        for i in (N + 1):Nbuff
            cells[i] = make_dead_cell_invivo()
        end
    end

    return cells
end

function birth_mutate_event_invivo!(cell_arr,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts)

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].E = cell_arr[cell_pos].E
    cell_arr[birth_pos].EG = cell_arr[cell_pos].EG
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
    else
        if mu > mut_p
            cell_arr[birth_pos].R = true
            phen_counts.Rcount += 1
        else
            phen_counts.Scount += 1
        end
    end
end

function _core_grow_kill_abm_invivo!(
    cells::Vector{InVivoCancerCell},
    params::ResPopInVivoParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    bmax = params.b
    dmax = params.d
    lam = bmax - dmax

    c_bdmax = bmax + dmax
    if params.del > 0.0 && sim.R_real == "d"
        dmax = dmax + (lam * params.del)
        c_bdmax = bmax + dmax
    end

    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
        t_bdmax = if de == :d || de == :c
            bmax + params.d + maximum(drug_concs["dconc"])
        else
            bmax + params.d
        end
        bdmax = max(c_bdmax, t_bdmax)
    else
        drug_concs = nothing
        bdmax = c_bdmax
    end

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_cell_invivo)
    live_vec = live_positions(cells)
    dead_vec = dead_positions(cells)

    Nt = n_alive(cells)
    Rcount = Int64(sum(cell -> cell.R && cell.alive, cells))
    Ecount = Int64(sum(cell -> cell.E && cell.alive, cells))
    Scount = Int64(sum(cell -> !cell.R && !cell.E && cell.alive, cells))
    phen_counts = PhenoCounts(Scount, Rcount, Ecount)

    Nvec = Int64[Nt]
    Svec = Int64[phen_counts.Scount]
    Rvec = Int64[phen_counts.Rcount]
    Evec = Int64[phen_counts.Ecount]
    Pvec = Int64[sim.Passage]
    nS_EG0, nS_EG1, nR_EG0, nR_EG1, nE_EG0, nE_EG1 = _invivo_abm_pheno_eg_counts(cells)
    S_EG0_vec = Int64[nS_EG0]
    S_EG1_vec = Int64[nS_EG1]
    R_EG0_vec = Int64[nR_EG0]
    R_EG1_vec = Int64[nR_EG1]
    E_EG0_vec = Int64[nE_EG0]
    E_EG1_vec = Int64[nE_EG1]

    while t <= sim.tmax
        Nt > 0 || break

        live_pos = rand(1:length(cells))
        dead_pos = rand(1:length(cells))

        while live_vec[live_pos] == 0
            live_pos = rand(1:length(cells))
        end

        while dead_vec[dead_pos] == 0
            dead_pos = rand(1:length(cells))
        end

        ran = rand(Uniform(0, bdmax))
        dt = -log(rand()) / (bdmax * Nt)
        t += dt

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            _push_invivo_abm_pheno_eg_counts!(S_EG0_vec, S_EG1_vec,
                                              R_EG0_vec, R_EG1_vec,
                                              E_EG0_vec, E_EG1_vec, cells)
            t_rec += t_rec_change
        end

        if t >= sim.tmax || Nt >= sim.Nmax
            push!(Nvec, Nt)
            push!(tvec, min(t, sim.tmax))
            push!(Svec, phen_counts.Scount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            _push_invivo_abm_pheno_eg_counts!(S_EG0_vec, S_EG1_vec,
                                              R_EG0_vec, R_EG1_vec,
                                              E_EG0_vec, E_EG1_vec, cells)
            break
        end

        if params.del > 0.0 && cells[live_pos].R
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

        if treat
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
            else
                cell_b -= if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
                if cell_b < 0.0
                    cell_d += abs(cell_b)
                    cell_b = 0.0
                end
            end
        end

        cell_b *= (1 - (Nt / sim.Cc))
        cell_d *= (1 - (Nt / sim.Cc))

        if ran < cell_b
            if treat
                curr_rconc = curr_rc(t, drug_concs)
                al_scal = params.al * curr_rconc
                birth_mutate_event_invivo!(cells, live_pos, dead_pos,
                                           params.mu, params.sig, al_scal,
                                           phen_counts)
            else
                birth_mutate_event_invivo!(cells, live_pos, dead_pos,
                                           params.mu, params.sig, 0.0,
                                           phen_counts)
            end
            live_vec[dead_pos] = 1
            dead_vec[dead_pos] = 0
            Nt += 1
        end

        if bmax <= ran < bmax + cell_d
            death_event!(cells, live_pos, phen_counts)
            live_vec[live_pos] = 0
            dead_vec[live_pos] = 1
            Nt -= 1
        end
    end

    fin_t = round(min(t, sim.tmax), digits = 4)
    return ResPopInVivoGrowOut(Nvec, tvec, Svec, Rvec, Evec, Pvec, fin_t,
                               S_EG0_vec, S_EG1_vec,
                               R_EG0_vec, R_EG1_vec,
                               E_EG0_vec, E_EG1_vec)
end

function run_model_core_abm(model::ResPopInVivo_ABM, state::ResPopInVivoABMState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm_invivo!(state.cells, model.params, sim; treat = treat)
end
