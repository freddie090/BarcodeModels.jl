const GAM_INDEX = 1
const NS_INDEX = 2
const NR_INDEX = 3
const NE_INDEX = 4
const PASS_INDEX = 5

# Population and logistic helpers
total_population(u) = u[NS_INDEX] + u[NR_INDEX] + u[NE_INDEX]
logistic_factor(N::Real, Cc::Real) = 1 - (N / Cc)
logistic_factor(u::AbstractVector, Cc::Real) = logistic_factor(total_population(u), Cc)

"""Select a behavior based on `drug_effect` (:d, :b, or :c)."""
function select_drug_effect(de, f_d, f_b, f_c)
    return de == :d ? f_d : de == :b ? f_b : de == :c ? f_c :
           error("Invalid drug_effect value: must be :b, :d, or :c")
end

# Drug-effect modifiers for ODE rates (no allocations inside ode_fxn!).
function apply_drug_effect_death(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_mod = b * lf
    d_mod = (d + ((d / d_ref) * Dc * gam * (1 - psi))) * lf
    return b_mod, d_mod
end

function apply_drug_effect_birth(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_mod = (b - ((b / b_ref) * Dc * gam * (1 - psi))) * lf
    d_mod = d * lf
    return b_mod, d_mod
end

function apply_drug_effect_combined(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_eff = b - ((b / b_ref) * Dc * gam * (1 - psi))
    if b_eff >= 0
        b_mod = b_eff * lf
        d_mod = d * lf
    else
        b_mod = 0.0
        d_mod = (d + abs(b_eff)) * lf
    end
    return b_mod, d_mod
end

"""Return the population slice (sensitive, resistant, escape)."""
pop_fun(x) = x[NS_INDEX:NE_INDEX]

"""Return the passage counter from a state vector."""
pass_fun(x) = x[PASS_INDEX]

"""Draw multivariate hypergeometric counts (uses RNG)."""
function multivariate_hypergeometric_draw(population_sizes, num_samples)
    remaining_samples = num_samples
    total_population = sum(population_sizes)
    draws = zeros(Int, length(population_sizes))

    for i in 1:(length(population_sizes) - 1)
        if remaining_samples > 0
            draws[i] = rand(Hypergeometric(population_sizes[i],
                                            total_population - population_sizes[i],
                                            remaining_samples))
            remaining_samples -= draws[i]
            total_population -= population_sizes[i]
        end
    end
    draws[end] = remaining_samples
    return draws
end

"""Resistance population (res_pop) model (validates parameters on construction)."""
struct ResPop <: AbstractBarcodeModel
    params::ModelParams
    function ResPop(params::ModelParams)
        validate_model_params(params)
        new(params)
    end
end

"""Construct a `ResPop` from keyword arguments forwarded to `ModelParams`."""
ResPop(; kwargs...) = ResPop(ModelParams(; kwargs...))

"""Build component-rate parameters and return them with baseline rates."""
function build_component_params(params::ModelParams)
    bS = params.b
    dS = params.d

    bR = params.b * (1 - params.del)
    dR = params.d * (1 - params.del)

    bE = params.b
    dE = params.d

    p = ComponentArray(
        bS_o = 0.0, dS_o = 0.0, bS_j = bS, dS_j = dS,
        bR_o = 0.0, dR_o = 0.0, bR_j = bR, dR_j = dR,
        bE_o = 0.0, dE_o = 0.0, bE_j = bE, dE_j = dE,
        Dc = params.Dc, kp = 0.0, psi = params.psi,
        mu_o = 0.0, sig_o = 0.0, al_o = 0.0,
        mu_j = params.mu, sig_j = params.sig, al_j = params.al
    )

    return p, bS, dS, bR, dR, bE, dE
end

"""
Run the hybrid growth/kill model with a model, state, and simulation settings.

Arguments: `model::ResPop`, `state::ModelState`, `sim::SimParams`.
Returns: `ODESolution`.
"""
function simulate_grow_kill(model::ResPop, state::ModelState, sim::SimParams)
    params = model.params
    de = params.drug_effect

    if de == :b
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")
    end

    p, bS, dS, bR, dR, bE, dE = build_component_params(params)

    u0 = to_componentarray(state)
    tspan = (sim.t0, sim.tmax)

    rate_modifier = select_drug_effect(de,
        apply_drug_effect_death,
        apply_drug_effect_birth,
        apply_drug_effect_combined
    )

    t_pass_vec = sim.t_Pass isa AbstractVector ? Vector{Float64}(sim.t_Pass) : fill(Float64(sim.t_Pass), sim.n_Pass)
    if length(t_pass_vec) < sim.n_Pass
        error("t_Pass length must be >= n_Pass")
    end

    # Helper builders are grouped here to keep repeated rate/callback logic centralized.

    # Build birth/death rate closures for a phenotype and drug mode.
    function build_birth_death_rates(de, idx, b_ref, d_ref, b_sym::Symbol, d_sym::Symbol; psi_fn = p -> 0.0)
        function logistic_scale(u)
            N = total_population(u)
            return logistic_factor(N, sim.Cc)
        end

        function birth_drug_adjustment(u, p)
            return ((getproperty(p, b_sym) / b_ref) * p.Dc * u[GAM_INDEX] * (1 - psi_fn(p)))
        end

        function death_drug_adjustment(u, p)
            return ((getproperty(p, d_sym) / d_ref) * p.Dc * u[GAM_INDEX] * (1 - psi_fn(p)))
        end

        function birth_d(u, p, t)
            return (u[idx] * getproperty(p, b_sym)) * logistic_scale(u)
        end

        function death_d(u, p, t)
            drug_interaction = death_drug_adjustment(u, p)
            return (u[idx] * (getproperty(p, d_sym) + drug_interaction)) * logistic_scale(u)
        end

        function birth_b(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            return (u[idx] * (getproperty(p, b_sym) - drug_interaction)) * logistic_scale(u)
        end

        function death_b(u, p, t)
            return (u[idx] * getproperty(p, d_sym)) * logistic_scale(u)
        end

        function birth_c(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            b_drug = getproperty(p, b_sym) - drug_interaction
            if b_drug >= 0.0
                return (u[idx] * b_drug) * logistic_scale(u)
            else
                return 0.0
            end
        end

        function death_c(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            b_drug = getproperty(p, b_sym) - drug_interaction
            if b_drug >= 0.0
                return (u[idx] * getproperty(p, d_sym)) * logistic_scale(u)
            else
                negative_surplus = abs(b_drug)
                return (u[idx] * (getproperty(p, d_sym) + negative_surplus)) * logistic_scale(u)
            end
        end

        birth = select_drug_effect(de, birth_d, birth_b, birth_c)
        death = select_drug_effect(de, death_d, death_b, death_c)
        return birth, death
    end

    # Build a switching-rate closure for a phenotype and drug mode.
    function build_switch_rate(de, idx, b_total_fn, rate_fn; gamma_factor=false)
        function rate_d(u, p, t)
            N = total_population(u)
            base = u[idx] * rate_fn(p) * b_total_fn(p)
            if gamma_factor
                base *= u[GAM_INDEX]
            end
            return base * logistic_factor(N, sim.Cc)
        end

        function rate_b(u, p, t)
            N = total_population(u)
            b_total = b_total_fn(p)
            drug_interaction = (b_total * p.Dc * u[GAM_INDEX])
            b_drug = b_total - drug_interaction
            base = u[idx] * b_drug * rate_fn(p)
            if gamma_factor
                base *= u[GAM_INDEX]
            end
            return base * logistic_factor(N, sim.Cc)
        end

        function rate_c(u, p, t)
            N = total_population(u)
            b_total = b_total_fn(p)
            drug_interaction = (b_total * p.Dc * u[GAM_INDEX])
            b_drug = b_total - drug_interaction
            if b_drug >= 0.0
                base = u[idx] * b_drug * rate_fn(p)
                if gamma_factor
                    base *= u[GAM_INDEX]
                end
                return base * logistic_factor(N, sim.Cc)
            else
                return 0.0
            end
        end

        return select_drug_effect(de, rate_d, rate_b, rate_c)
    end

    # ODE system for deterministic dynamics (mutates `du` in-place).
    function ode_fxn!(du, u, p, t)
        @unpack kp, psi, mu_o, sig_o, al_o = p
        @unpack gam, nS, nR, nE = u

        N = nS + nR + nE

        bS_o_mod, dS_o_mod = rate_modifier(bS, dS, p.bS_o, p.dS_o, p.Dc, gam, N, sim.Cc, 0.0)
        bR_o_mod, dR_o_mod = rate_modifier(bR, dR, p.bR_o, p.dR_o, p.Dc, gam, N, sim.Cc, psi)
        bE_o_mod, dE_o_mod = rate_modifier(bE, dE, p.bE_o, p.dE_o, p.Dc, gam, N, sim.Cc, psi)

        du.gam = kp
        du.nS = (bS_o_mod - dS_o_mod) * nS -
            mu_o * nS * bS_o_mod +
            sig_o * nR * bR_o_mod
        du.nR = (bR_o_mod - dR_o_mod) * nR +
            mu_o * nS * bS_o_mod -
            sig_o * nR * bR_o_mod -
            al_o * gam * nR * bR_o_mod
        du.nE = (bE_o_mod - dE_o_mod) * nE +
                al_o * gam * nR * bR_o_mod
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    ################################
    # Jump Problem Birth-Death Jumps
    ################################

    # Sensitive
    # Event handlers (mutate integrator state)
    function S_birth!(integrator)
        integrator.u[NS_INDEX] += 1
        nothing
    end

    function S_death!(integrator)
        integrator.u[NS_INDEX] -= 1
        nothing
    end

    Sb_rate, Sd_rate = build_birth_death_rates(de, NS_INDEX, bS, dS, :bS_j, :dS_j)

    Sb_jump = VariableRateJump(Sb_rate, S_birth!)
    Sd_jump = VariableRateJump(Sd_rate, S_death!)

    # Resistant
    function R_birth!(integrator)
        integrator.u[NR_INDEX] += 1
        nothing
    end

    function R_death!(integrator)
        integrator.u[NR_INDEX] -= 1
        nothing
    end

    Rb_rate, Rd_rate = build_birth_death_rates(de, NR_INDEX, bR, dR, :bR_j, :dR_j; psi_fn = p -> p.psi)

    Rb_jump = VariableRateJump(Rb_rate, R_birth!)
    Rd_jump = VariableRateJump(Rd_rate, R_death!)

    # Escape
    function E_birth!(integrator)
        integrator.u[NE_INDEX] += 1
        nothing
    end

    function E_death!(integrator)
        integrator.u[NE_INDEX] -= 1
        nothing
    end

    Eb_rate, Ed_rate = build_birth_death_rates(de, NE_INDEX, bE, dE, :bE_j, :dE_j; psi_fn = p -> p.psi)

    Eb_jump = VariableRateJump(Eb_rate, E_birth!)
    Ed_jump = VariableRateJump(Ed_rate, E_death!)

    ########################################
    # Jump Problem Phenotype Switching Jumps
    ########################################

    # Sensitive -> Resistant
    function SR_switch!(integrator)
        integrator.u[NS_INDEX] = integrator.u[NS_INDEX] - 1
        integrator.u[NR_INDEX] = integrator.u[NR_INDEX] + 1
        nothing
    end

    SR_switch_rate = build_switch_rate(de, NS_INDEX, p -> (p.bS_o + p.bS_j), p -> p.mu_j)
    SR_switch_jump = VariableRateJump(SR_switch_rate, SR_switch!)

    # Resistant -> Sensitive
    function RS_switch!(integrator)
        integrator.u[NR_INDEX] = integrator.u[NR_INDEX] - 1
        integrator.u[NS_INDEX] = integrator.u[NS_INDEX] + 1
        nothing
    end

    RS_switch_rate = build_switch_rate(de, NR_INDEX, p -> (p.bR_o + p.bR_j), p -> p.sig_j)
    RS_switch_jump = VariableRateJump(RS_switch_rate, RS_switch!)

    # Resistant -> Escape
    function RE_switch!(integrator)
        integrator.u[NR_INDEX] = integrator.u[NR_INDEX] - 1
        integrator.u[NE_INDEX] = integrator.u[NE_INDEX] + 1
        nothing
    end

    RE_switch_rate = build_switch_rate(de, NR_INDEX, p -> (p.bR_o + p.bR_j), p -> p.al_j; gamma_factor=true)
    RE_switch_jump = VariableRateJump(RE_switch_rate, RE_switch!)

    ###########
    # Callbacks
    ###########

    mu = params.mu
    sig = params.sig
    al = params.al
    n0 = sim.n0

    # ODE/jump toggle callbacks for phenotype birth/death rates.
    function build_bd_toggle_callbacks(idx, b_o_sym, b_j_sym, d_o_sym, d_j_sym, Nswitch)
        cond_to_ode(u, t, integrator) = integrator.u[idx] >= Nswitch
        function switch_to_ode!(integrator)
            if getproperty(integrator.p, b_o_sym) == 0.0
                setproperty!(integrator.p, b_o_sym, deepcopy(getproperty(integrator.p, b_j_sym)))
                setproperty!(integrator.p, b_j_sym, 0.0)
            end
            if getproperty(integrator.p, d_o_sym) == 0.0
                setproperty!(integrator.p, d_o_sym, deepcopy(getproperty(integrator.p, d_j_sym)))
                setproperty!(integrator.p, d_j_sym, 0.0)
            end
            nothing
        end
        cb1 = DiscreteCallback(cond_to_ode, switch_to_ode!,
                               save_positions = (false, true))

        cond_to_jump(u, t, integrator) = integrator.u[idx] < Nswitch
        function switch_to_jump!(integrator)
            integrator.u[idx] = round(integrator.u[idx])
            if getproperty(integrator.p, b_j_sym) == 0.0
                setproperty!(integrator.p, b_j_sym, deepcopy(getproperty(integrator.p, b_o_sym)))
                setproperty!(integrator.p, b_o_sym, 0.0)
            end
            if getproperty(integrator.p, d_j_sym) == 0.0
                setproperty!(integrator.p, d_j_sym, deepcopy(getproperty(integrator.p, d_o_sym)))
                setproperty!(integrator.p, d_o_sym, 0.0)
            end
            nothing
        end
        cb2 = DiscreteCallback(cond_to_jump, switch_to_jump!,
                               save_positions = (false, true))

        return cb1, cb2
    end

    # ODE/jump toggle callbacks for switch rates.
    function build_rate_toggle_callbacks(idx, rate_value, rate_o_sym, rate_j_sym, epsi)
        cond_to_ode(u, t, integrator) = (integrator.u[idx] * rate_value) >= epsi
        function switch_to_ode!(integrator)
            if getproperty(integrator.p, rate_o_sym) == 0.0
                setproperty!(integrator.p, rate_o_sym, deepcopy(getproperty(integrator.p, rate_j_sym)))
                setproperty!(integrator.p, rate_j_sym, 0.0)
            end
            nothing
        end
        cb1 = DiscreteCallback(cond_to_ode, switch_to_ode!,
                               save_positions = (false, true))

        cond_to_jump(u, t, integrator) = (integrator.u[idx] * rate_value) < epsi
        function switch_to_jump!(integrator)
            if getproperty(integrator.p, rate_j_sym) == 0.0
                setproperty!(integrator.p, rate_j_sym, deepcopy(getproperty(integrator.p, rate_o_sym)))
                setproperty!(integrator.p, rate_o_sym, 0.0)
            end
            nothing
        end
        cb2 = DiscreteCallback(cond_to_jump, switch_to_jump!,
                               save_positions = (false, true))

        return cb1, cb2
    end

    # Resample populations for a passage bottleneck (mutates state; uses RNG).
    function resample_passage_bottleneck!(integrator, n0)
        curr_nR = round(Int, integrator.u[NR_INDEX])
        curr_nS = round(Int, integrator.u[NS_INDEX])
        curr_nE = round(Int, integrator.u[NE_INDEX])

        counts = multivariate_hypergeometric_draw([curr_nR, curr_nS, curr_nE], n0)

        integrator.u[NR_INDEX] = Float64(counts[1])
        integrator.u[NS_INDEX] = Float64(counts[2])
        integrator.u[NE_INDEX] = Float64(counts[3])
        nothing
    end

    function attempt_passage!(integrator, n0; require_enough::Bool)
        if require_enough
            curr_nR = round(Int, integrator.u[NR_INDEX])
            curr_nS = round(Int, integrator.u[NS_INDEX])
            curr_nE = round(Int, integrator.u[NE_INDEX])
            if (curr_nR + curr_nS + curr_nE) < n0
                return false
            end
        end

        resample_passage_bottleneck!(integrator, n0)
        integrator.u[PASS_INDEX] += 1
        return true
    end

    # Non-negativity clamp callback for a population component.
    function build_clamp_callback(idx)
        condition(u, t, integrator) = u[idx] < 0
        function affect!(integrator)
            integrator.u[idx] = 0.0
        end
        return DiscreteCallback(condition, affect!,
                                save_positions = (false, true))
    end

    # Sensitive
    S_cb_switch1, S_cb_switch2 = build_bd_toggle_callbacks(
        NS_INDEX, :bS_o, :bS_j, :dS_o, :dS_j, sim.Nswitch
    )

    # Resistant
    R_cb_switch1, R_cb_switch2 = build_bd_toggle_callbacks(
        NR_INDEX, :bR_o, :bR_j, :dR_o, :dR_j, sim.Nswitch
    )

    # Escape
    E_cb_switch1, E_cb_switch2 = build_bd_toggle_callbacks(
        NE_INDEX, :bE_o, :bE_j, :dE_o, :dE_j, sim.Nswitch
    )

    # Sensitive -> Resistant
    SR_cb_switch1, SR_cb_switch2 = build_rate_toggle_callbacks(
        NS_INDEX, mu, :mu_o, :mu_j, sim.epsi
    )

    # Resistant -> Sensitive
    RS_cb_switch1, RS_cb_switch2 = build_rate_toggle_callbacks(
        NR_INDEX, sig, :sig_o, :sig_j, sim.epsi
    )

    # Resistant -> Escape
    RE_cb_switch1, RE_cb_switch2 = build_rate_toggle_callbacks(
        NR_INDEX, al, :al_o, :al_j, sim.epsi
    )

    # Treatment schedule callbacks (use preset times for robustness).
    function set_treatment_on!(integrator)
        integrator.p.kp = params.k
    end

    function set_treatment_off!(integrator)
        integrator.p.kp = -params.k
    end

    filter_times_in_span(times) = sort(unique(filter(t -> (t >= sim.t0 && t <= sim.tmax), times)))

    treat_on_times = filter_times_in_span(sim.treat_ons)
    treat_off_times = filter_times_in_span(sim.treat_offs)

    treat_on_cb = isempty(treat_on_times) ? nothing :
                  PresetTimeCallback(treat_on_times, set_treatment_on!; save_positions = (false, true))
    treat_off_cb = isempty(treat_off_times) ? nothing :
                   PresetTimeCallback(treat_off_times, set_treatment_off!; save_positions = (false, true))

    # Drug level clamping
    gam_above_max(u, t, integrator) = (u[GAM_INDEX] - 1.0)
    function clamp_gam_max!(integrator)
        integrator.u[GAM_INDEX] = 1.0
        integrator.p.kp = 0.0
    end

    gam_below_min(u, t, integrator) = (u[GAM_INDEX])
    function clamp_gam_min!(integrator)
        integrator.u[GAM_INDEX] = 0.0
        integrator.p.kp = 0.0
    end

    gam_max_cb = ContinuousCallback(gam_above_max, clamp_gam_max!,
                                    save_positions = (false, true))
    gam_min_cb = ContinuousCallback(gam_below_min, clamp_gam_min!,
                                    save_positions = (false, true))

    # Nmax, Extinction & Passage
    nmax_reached(u, t, integrator) = (total_population(integrator.u) >= sim.Nmax)

    function handle_nmax!(integrator)
        if integrator.u[PASS_INDEX] >= sim.n_Pass
            terminate!(integrator)
        else
            attempt_passage!(integrator, n0; require_enough = false)
        end
    end

    nmax_cb = DiscreteCallback(nmax_reached, handle_nmax!,
                               save_positions = (false, true))

    extinction(u, t, integrator) = (total_population(integrator.u) < 1.0)
    extinction_cb = DiscreteCallback(extinction, terminate!,
                                     save_positions = (false, false))

    function apply_scheduled_passage!(integrator)
        curr_pass = round(Int, integrator.u[PASS_INDEX])
        if curr_pass < sim.n_Pass
            expected_t = t_pass_vec[curr_pass]
            # Only attempt the passage scheduled for the current pass.
            if integrator.t == expected_t
                attempt_passage!(integrator, n0; require_enough = true)
            end
        end
    end

    passage_times = filter_times_in_span(t_pass_vec)
    passage_cb = isempty(passage_times) ? nothing :
                 PresetTimeCallback(passage_times, apply_scheduled_passage!; save_positions = (false, true))

    # Population Limits
    cb7a = build_clamp_callback(NS_INDEX)
    cb7b = build_clamp_callback(NR_INDEX)
    cb7c = build_clamp_callback(NE_INDEX)

    # Collect Callbacks
    cb_list = Any[
        S_cb_switch1, S_cb_switch2,
        R_cb_switch1, R_cb_switch2,
        E_cb_switch1, E_cb_switch2,
        SR_cb_switch1, SR_cb_switch2,
        RS_cb_switch1, RS_cb_switch2,
        RE_cb_switch1, RE_cb_switch2,
        nmax_cb, extinction_cb,
        cb7a, cb7b, cb7c
    ]

    if passage_cb !== nothing
        push!(cb_list, passage_cb)
    end

    if sim.treat == true
        if treat_on_cb !== nothing
            push!(cb_list, treat_on_cb)
        end
        if treat_off_cb !== nothing
            push!(cb_list, treat_off_cb)
        end
        push!(cb_list, gam_max_cb, gam_min_cb)
    end

    cbs = CallbackSet(cb_list...)

    sjm_prob = JumpProblem(prob, Direct(),
                           Sb_jump, Sd_jump,
                           Rb_jump, Rd_jump,
                           Eb_jump, Ed_jump,
                           SR_switch_jump,
                           RS_switch_jump,
                           RE_switch_jump)

    base_tstops = collect(sim.t0:sim.save_at:sim.tmax)
    event_times = passage_times
    if sim.treat == true
        event_times = vcat(event_times, treat_on_times, treat_off_times)
    end
    tstops = isempty(event_times) ? base_tstops : sort(unique(vcat(base_tstops, event_times)))

    sol = @suppress solve(sjm_prob, Tsit5(), callback = cbs,
                          tstops = tstops)

    return sol
end

"""
Run the hybrid growth/kill model using primitive arguments.

Arguments: see signature.
Returns: `ODESolution`.
"""
function simulate_grow_kill(n0::Int64, nS::Float64, nR::Float64, nE::Float64,
    b::Float64, d::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64,
    t_Pass::Union{Float64, Vector{Float64}},
    Nmax::Int64, Cc::Int64,
    treat_ons::Vector{Float64}, treat_offs::Vector{Float64},
    Nswitch::Int64; save_at::Float64 = 0.5, treat::Bool = false,
    n_Pass::Int64 = 1, epsi::Float64 = 100.0,
    drug_effect::Union{String, Symbol} = "d")

    model = ResPop(ModelParams(
        b = b, d = d, mu = mu, sig = sig, del = del, al = al,
        Dc = Dc, k = k, psi = psi, drug_effect = drug_effect
    ))

    state = ModelState(nS, nR, nE; gam = 0.0, pass_num = 1)

    sim = SimParams(
        n0 = n0, t0 = t0, tmax = tmax, t_Pass = t_Pass,
        Nmax = Nmax, Cc = Cc, Nswitch = Nswitch,
        treat_ons = treat_ons, treat_offs = treat_offs,
        save_at = save_at, treat = treat, n_Pass = n_Pass, epsi = epsi
    )

    return simulate_grow_kill(model, state, sim)
end
