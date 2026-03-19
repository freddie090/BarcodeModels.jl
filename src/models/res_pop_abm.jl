"""Resistance population agent-based (ABM) model (validates parameters on construction)."""
struct ResPop_ABM <: ABMModel
    params::ResPopParams
    abm::ABMParams
    function ResPop_ABM(params::ResPopParams, abm::ABMParams)
        validate_model_params(params)
        new(params, abm)
    end
end

"""Construct a `ResPop_ABM` from parameters and ABM settings."""
ResPop_ABM(params::ResPopParams; abm::ABMParams = ABMParams()) = ResPop_ABM(params, abm)

"""Construct a `ResPop_ABM` from keyword arguments forwarded to `ResPopParams`."""
ResPop_ABM(; abm::ABMParams = ABMParams(), kwargs...) = ResPop_ABM(ResPopParams(; kwargs...), abm)

mutable struct CancerCell
    barcode::Float64
    R::Bool
    E::Bool
    alive::Bool
end

mutable struct PhenoCounts
    Scount::Int64
    Rcount::Int64
    Ecount::Int64
end

mutable struct GrowOut
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    Svec::Vector{Int64}
    Rvec::Vector{Int64}
    Evec::Vector{Int64}
    Pvec::Vector{Int64}
    fin_t::Float64
end

mutable struct ResPopABMState
    cells::Vector{CancerCell}
end

struct ABMSimParams
    t0::Float64
    tmax::Float64
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    dt_save_at::Float64
    R_real::String
    t_frac::Float64
    Passage::Int64
    drug_effect::Symbol
end

function ABMSimParams(; t0 = 0.0, tmax, Nmax, Cc, treat_ons, treat_offs,
    dt_save_at = 0.1, R_real = "b", t_frac = 0.05, Passage = 1, drug_effect = :d)

    return ABMSimParams(
        Float64(t0),
        Float64(tmax),
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Float64(dt_save_at),
        String(R_real),
        Float64(t_frac),
        Int64(Passage),
        normalize_respop_drug_effect(drug_effect)
    )
end

make_dead_cell() = CancerCell(0.0, false, false, false)

function seed_cells(N::Int64, rho::Float64, Nbuff::Int64;
    skew_lib::Bool = false, bc_unif::Float64 = 0.0, Nbc::Int64 = 0)

    cells = Vector{CancerCell}(undef, max(N, Nbuff))

    if skew_lib
        bcs = collect(1:Nbc)
        bc_probs = generate_probabilities(bcs, bc_unif)
        samp_bcs = bcs[rand(Categorical(bc_probs), N)]
    else
        samp_bcs = collect(1:N)
    end

    for i in 1:N
        cells[i] = CancerCell(samp_bcs[i], false, false, true)
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
            cells[i] = make_dead_cell()
        end
    end

    return cells
end

function _core_grow_kill_abm!(
    cells::Vector{CancerCell},
    params::ResPopParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    @assert 0.0 <= sim.t_frac <= 1.0 "t_frac must be between 0 and 1."
    @assert params.drug_effect in RESPOP_DRUG_EFFECTS "drug effect can only take :b, :d, or :c."
    @assert sim.R_real in ["b", "d", "l"] "R_real can only take 'b', 'd' or 'l' as a value."

    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    bmax = params.b
    dmax = params.d
    lam = bmax - dmax

    if params.del > 0.0
        if sim.R_real == "d"
            dmax = dmax + (lam * params.del)
            c_bdmax = bmax + dmax
        else
            c_bdmax = bmax + dmax
        end
    else
        c_bdmax = bmax + dmax
    end

    drug_concs = nothing
    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
        if de == :d || de == :c
            dmax = params.d + maximum(drug_concs["dconc"])
            t_bdmax = bmax + dmax
        else
            t_bdmax = bmax + dmax
        end
    else
        t_bdmax = bmax + dmax
    end

    bdmax = maximum([c_bdmax, t_bdmax])

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_cell)
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

        ran = rand(Uniform(0, bdmax))
        denom = bdmax * Nt
        denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (bdmax*Nt=$(denom)).")
        dt = -log(rand()) / denom
        t += dt

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
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
            cell_b = max(cell_b, 0)
            cell_d = max(cell_d, 0)
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

        if ran < cell_b
            if treat
                curr_rconc = curr_rc(t, drug_concs)
                al_scal = params.al * curr_rconc
                birth_mutate_event!(cells, live_pos, dead_pos,
                                    params.mu, params.sig, al_scal,
                                    phen_counts)
            else
                birth_mutate_event!(cells, live_pos, dead_pos,
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

        if Nt <= 0
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            break
        end
    end

    fin_t = round(t, digits = 4)
    return GrowOut(Nvec, tvec, Svec, Rvec, Evec, Pvec, fin_t)
end

function run_model_core_abm(model::ResPop_ABM, state::ResPopABMState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm!(state.cells, model.params, sim; treat = treat)
end

function grow_kill_lin_kmc!(cells::Vector{CancerCell},
    b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64}, dt_save_at::Float64;
    R_real::String = "b", t_frac::Float64 = 0.050, treat::Bool = false,
    Passage::Int64 = 1,
    drug_effect::Union{String, Symbol} = :d)

    model = ResPop_ABM(ResPopParams(
        b = b, d = d, mu = mu, sig = sig, del = del, al = al,
        Dc = Dc, k = k, psi = psi, drug_effect = drug_effect
    ))

    sim = ABMSimParams(
        t0 = t0, tmax = tmax, Nmax = Nmax, Cc = Cc,
        treat_ons = treat_ons, treat_offs = treat_offs, dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, Passage = Passage, drug_effect = drug_effect
    )

    state = ResPopABMState(cells)
    return run_model_core_abm(model, state, sim; treat = treat)
end


