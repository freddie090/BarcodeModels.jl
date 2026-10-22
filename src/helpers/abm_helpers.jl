function drug_treat_concs(tmax::Float64, k::Float64, Dc::Float64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64}, dt_save_at::Float64)

    # Discrete callbacks do not apply the initial-time onset; initialise its slope.
    p = [Dc, (0.0 in treat_ons && !(0.0 in treat_offs)) ? k : 0.0]
    u0 = [0.0, 0.0]
    tspan = (0.0, tmax)

    function ode_fxn!(du, u, p, t)
        Dc_local, kp = p
        du[1] = kp
        du[2] = kp * Dc_local
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    condition1(u, t, integrator) = any(ton -> isapprox(t, ton; atol = 1e-10), treat_ons)
    affect1!(integrator) = (integrator.p[2] = k)

    condition2(u, t, integrator) = any(toff -> isapprox(t, toff; atol = 1e-10), treat_offs)
    affect2!(integrator) = (integrator.p[2] = -k)

    condition3a(u, t, integrator) = integrator.p[2] > 0.0 ? (u[1] - 1.0) : 1.0
    function affect3a!(integrator)
        integrator.u[1] = 1.0
        integrator.u[2] = integrator.p[1]
        integrator.p[2] = 0.0
    end

    condition3b(u, t, integrator) = integrator.p[2] < 0.0 ? u[1] : 1.0
    function affect3b!(integrator)
        integrator.u[1] = 0.0
        integrator.u[2] = 0.0
        integrator.p[2] = 0.0
    end

    cb1 = DiscreteCallback(condition1, affect1!, save_positions = (false, true))
    cb2 = DiscreteCallback(condition2, affect2!, save_positions = (false, true))
    cb3a = ContinuousCallback(condition3a, affect3a!, save_positions = (false, true))
    cb3b = ContinuousCallback(condition3b, affect3b!, save_positions = (false, true))

    cbs = CallbackSet(cb1, cb2, cb3a, cb3b)

    save_times = collect(0.0:dt_save_at:tmax)
    filter_times_in_span(times) = sort(unique(filter(t -> (t >= 0.0 && t <= tmax), times)))
    event_times = vcat(filter_times_in_span(treat_ons), filter_times_in_span(treat_offs))
    tstops = isempty(event_times) ? save_times : sort(unique(vcat(save_times, event_times)))

    sol = solve(prob, OrdinaryDiffEqRosenbrock.Rodas5(), callback = cbs, tstops = tstops)

    return Dict(
        "t" => sol.t,
        "dconc" => map(x -> x[1] * Dc, sol.u),
        "rconc" => map(x -> x[1], sol.u)
    )
end

curr_dc(curr_t::Float64, drug_concs::Dict) = drug_concs["dconc"][find_closest(curr_t, drug_concs["t"])]
curr_rc(curr_t::Float64, drug_concs::Dict) = drug_concs["rconc"][find_closest(curr_t, drug_concs["t"])]

function generate_probabilities(X::Vector, prob_unif::Float64)
    n = length(X)
    if prob_unif == 0.0
        return vcat([1.0], zeros(n - 1))
    end
    prob_unif_vector = fill(prob_unif, n)
    return rand(Dirichlet(prob_unif_vector))
end

"""Sample initial ABM barcode identities using the selected library settings."""
function _sample_initial_barcodes(N::Int64;
    skew_lib::Bool = false, use_lib_probs::Bool = false, bc_unif::Float64 = 0.0,
    Nbc::Int64 = 0, bc_probs = Float64[])

    if use_lib_probs
        bc_probs_vec = Vector{Float64}(bc_probs)
        samp_bcs = rand(Categorical(bc_probs_vec), N)
    elseif skew_lib
        bcs = collect(1:Nbc)
        bc_probs = generate_probabilities(bcs, bc_unif)
        samp_bcs = bcs[rand(Categorical(bc_probs), N)]
    else
        samp_bcs = collect(1:N)
    end

    return samp_bcs
end

function find_first_dead_index(cell_arr)::Int
    for i in 1:length(cell_arr)
        if !cell_arr[i].alive
            return i
        end
    end
    return 0
end

function order_cells!(cell_arr, make_dead_cell::Function)
    alive_index = 1
    dead_index = length(cell_arr)

    while alive_index < dead_index
        if cell_arr[alive_index].alive
            alive_index += 1
        else
            cell_arr[alive_index], cell_arr[dead_index] = cell_arr[dead_index], cell_arr[alive_index]
            dead_index -= 1
        end
    end

    for i in alive_index:dead_index
        if !cell_arr[i].alive
            cell_arr[i] = make_dead_cell()
        end
    end
end

function n_alive(cell_arr)::Int
    count = 0
    for cell in cell_arr
        count += cell.alive ? 1 : 0
    end
    return count
end

live_positions(cells) = [cell.alive ? 1 : 0 for cell in cells]
dead_positions(cells) = [cell.alive ? 0 : 1 for cell in cells]
alive_cells(cells) = filter(cell -> cell.alive, cells)

"""One main-replicate instantaneous treatment and its observed outcome."""
mutable struct InstantKillState
    settings::ABMParams
    applied::Bool
    mutation_id::Union{Missing, Int64}
    carrier_fraction::Union{Missing, Float64}
    N_before::Union{Missing, Int64}
    N_killed::Union{Missing, Int64}
    N_after::Union{Missing, Int64}
end

_instant_kill_state(abm::ABMParams) = abm.insta_kill ?
    InstantKillState(abm, false, missing, missing, missing, missing, missing) : nothing

_instant_kill_due(::Nothing, start_t, end_t) = false
function _instant_kill_due(event::InstantKillState, start_t, end_t)
    event.applied && return false
    time = event.settings.inst_kill_time
    time < start_t && error("Instantaneous killing time $(time) was missed before core restart at $(start_t).")
    return time <= end_t
end

function _validate_instant_kill_run(abm::ABMParams, drug_treatment::Bool, horizons)
    abm.insta_kill || return nothing
    drug_treatment && throw(ArgumentError("insta_kill=true requires drug_treatment=false; both treatment types cannot be enabled."))
    all(t -> abm.inst_kill_time <= t, horizons) ||
        throw(ArgumentError("inst_kill_time must lie within every main replicate's simulation horizon."))
    return nothing
end

"""Select the closest eligible live-cell carrier fraction, with random ties."""
function _instant_resistance_mutation(cells, abm::ABMParams)
    counts = Dict{Int64, Int64}()
    N = 0
    for cell in cells
        cell.alive || continue
        N += 1
        cell.muts === nothing && error("Mutation-based killing requires neutral mutations on every live cell.")
        for id in cell.muts
            counts[id] = get(counts, id, Int64(0)) + 1
        end
    end
    ids = sort!(collect(keys(counts)))
    fractions = [counts[id] / N for id in ids]
    distances = abs.(fractions .- abm.neut_mut_res_prop)
    # Allow floating-point roundoff at inclusive boundaries and equal-distance ties.
    roundoff = 4eps(Float64)
    eligible = findall(d -> d <= abm.neut_mut_res_prop_tol ||
        isapprox(d, abm.neut_mut_res_prop_tol; atol = roundoff, rtol = 0), distances)
    if isempty(eligible)
        available = isempty(fractions) ? "no observed mutations" : "$(minimum(fractions)) to $(maximum(fractions))"
        throw(ArgumentError("No mutation matches target $(abm.neut_mut_res_prop) +/- $(abm.neut_mut_res_prop_tol); available carrier fractions: $(available)."))
    end
    closest = minimum(distances[eligible])
    ties = filter(i -> isapprox(distances[i], closest; atol = roundoff, rtol = 0), eligible)
    chosen = length(ties) == 1 ? only(ties) : rand(ties)
    return ids[chosen], fractions[chosen]
end

"""Simulation-local neutral mutation IDs and bulk observation tables."""
mutable struct NeutralMutationTracker
    next_id::Int64
    poisson::Poisson{Float64}
    counts::Vector{DataFrame}
    sample_counts::Vector{DataFrame}
end

function _neutral_tracker(abm::ABMParams)
    abm.neutral_mut || return nothing
    return NeutralMutationTracker(abm.clonal_mut + 1, Poisson(abm.mut_rate), DataFrame[], DataFrame[])
end

"""Assign shared clonal IDs once, before any expansion or replicate splitting."""
function _seed_neutral_mutations!(cells, abm::ABMParams)
    abm.neutral_mut || return nothing
    for cell in cells
        cell.muts = cell.alive ? collect(Int64, 1:abm.clonal_mut) : nothing
    end
    return nothing
end

# A standalone core call can initialise fresh cells or resume an existing population.
# Experiment workflows instead share one tracker across all cores and replicates.
function _neutral_tracker(abm::ABMParams, cells)
    tracker = _neutral_tracker(abm)
    tracker === nothing && return nothing
    for cell in cells
        cell.alive || continue
        if cell.muts === nothing
            cell.muts = collect(Int64, 1:abm.clonal_mut)
        elseif !isempty(cell.muts)
            tracker.next_id = max(tracker.next_id, maximum(cell.muts) + 1)
        end
    end
    return tracker
end

_inherit_neutral_mutations!(parent, daughter, ::Nothing) = nothing

function _inherit_neutral_mutations!(parent, daughter, tracker::NeutralMutationTracker)
    inherited = parent.muts
    inherited === nothing && error("Neutral mutation tracking requires initialised live cells.")
    daughter.muts = copy(inherited)
    n_new = rand(tracker.poisson)
    n_new <= typemax(Int64) - tracker.next_id || error("Neutral mutation IDs exhausted Int64 capacity.")
    append!(daughter.muts, tracker.next_id:(tracker.next_id + n_new - 1))
    tracker.next_id += n_new
    return nothing
end

function extend_with_dead_cells!(cells, Nbuff::Int64, make_dead_cell::Function)
    initial_length = length(cells)
    if Nbuff > initial_length
        append!(cells, [make_dead_cell() for _ in 1:(Nbuff - initial_length)])
    end
end
