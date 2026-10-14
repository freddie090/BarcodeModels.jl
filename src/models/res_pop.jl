"""Resistance population (res_pop) model (validates parameters on construction)."""
struct ResPop <: HybridModel
    params::ResPopParams
    function ResPop(params::ResPopParams)
        validate_model_params(params)
        new(params)
    end
end

"""Construct a `ResPop` from keyword arguments forwarded to `ResPopParams`."""
ResPop(; kwargs...) = ResPop(ResPopParams(; kwargs...))

"""Build full biological rates and event representation flags, with baseline rates."""
function build_component_params(params::ResPopParams)
    bS = params.b
    dS = params.d

    bR = params.b * (1 - params.del)
    dR = params.d * (1 - params.del)

    bE = params.b
    dE = params.d

    phenotype_modes = params.transition_coupling === :birth_coupled ? (
        S_asymmetric_R_birth_ode = 0.0,
        R_asymmetric_S_birth_ode = 0.0,
        R_asymmetric_E_birth_ode = 0.0,
    ) : (
        S_to_R_transition_ode = 0.0,
        R_to_S_transition_ode = 0.0,
        R_to_E_transition_ode = 0.0,
    )

    p = ComponentArray(;
        bS = bS, dS = dS, bR = bR, dR = dR, bE = bE, dE = dE,
        Dc = params.Dc, kp = 0.0, psi = params.psi,
        mu = params.mu, sig = params.sig, al = params.al,
        S_symmetric_birth_ode = 0.0,
        S_death_ode = 0.0,
        R_symmetric_birth_ode = 0.0,
        R_death_ode = 0.0,
        E_symmetric_birth_ode = 0.0,
        E_death_ode = 0.0,
        phenotype_modes...
    )

    return p, bS, dS, bR, dR, bE, dE
end

"""
Run the hybrid growth/kill model with a model, state, and simulation settings.

Arguments: `model::ResPop`, `state::ResPopState`, `sim::SimParams`.
Returns: `ODESolution`.
"""
function run_model_core_hybrid(model::ResPop, state::ResPopState, sim::SimParams;
    treat::Bool = sim.treat, rng = Random.default_rng(), _setup_only::Bool = false, solve_kwargs...)
    sim_eff = treat == sim.treat ? sim : SimParams(
        n0 = sim.n0,
        t0 = sim.t0,
        tmax = sim.tmax,
        t_Pass = sim.t_Pass,
        Nmax = sim.Nmax,
        Cc = sim.Cc,
        Nswitch = sim.Nswitch,
        treat_ons = sim.treat_ons,
        treat_offs = sim.treat_offs,
        save_at = sim.save_at,
        treat = treat,
        Aswitch = sim.Aswitch
    )

    params = model.params
    independent = params.transition_coupling === :independent_rate
    de = params.drug_effect

    if de == :b
        bR = params.b * (1 - params.del)
        psi_scale = 1 - params.psi
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")
        if psi_scale > 0.0
            params.Dc <= (bR / psi_scale) || error("When drug_effect == :b, Dc*(1-psi) must be <= b*(1-del) to keep resistant birth non-negative. Use drug_effect == :c if stronger drug effects are intended.")
        end
    end

    p, bS, dS, bR, dR, bE, dE = build_component_params(params)

    u0 = to_componentarray(state)
    tspan = (sim_eff.t0, sim_eff.tmax)

    normalize_t_pass(times) = begin
        if times isa AbstractVector
            Vector{Float64}(times)
        else
            t_val = Float64(times)
            t_val < 0.0 ? Float64[] : [t_val]
        end
    end

    t_pass_vec = normalize_t_pass(sim_eff.t_Pass)
    n_pass_eff = length(t_pass_vec) + 1

    # Helper builders are grouped here to keep repeated rate/callback logic centralized.

    # Build birth/death rate closures for a phenotype and drug mode.
    function build_birth_death_rates(de, idx, b_ref, d_ref, ::Val{B}, ::Val{D}; psi_fn = p -> 0.0) where {B,D}
        safe_ratio(num, den) = abs(den) > eps(Float64) ? (num / den) : 0.0

        function logistic_scale(u)
            N = respop_total_population(u)
            return max(logistic_factor(N, sim_eff.Cc), 0.0)
        end

        function birth_drug_adjustment(u, p)
            return (safe_ratio(getproperty(p, B), b_ref) * p.Dc * clamp(u[RESPOP_GAM_INDEX], 0.0, 1.0) * (1 - psi_fn(p)))
        end

        function death_drug_adjustment(u, p)
            return (safe_ratio(getproperty(p, D), d_ref) * p.Dc * clamp(u[RESPOP_GAM_INDEX], 0.0, 1.0) * (1 - psi_fn(p)))
        end

        function birth_d(u, p, t)
            return (max(u[idx], 0.0) * getproperty(p, B)) * logistic_scale(u)
        end

        function death_d(u, p, t)
            drug_interaction = death_drug_adjustment(u, p)
            return (max(u[idx], 0.0) * (getproperty(p, D) + drug_interaction)) * logistic_scale(u)
        end

        function birth_b(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            return (max(u[idx], 0.0) * (getproperty(p, B) - drug_interaction)) * logistic_scale(u)
        end

        function death_b(u, p, t)
            return (max(u[idx], 0.0) * getproperty(p, D)) * logistic_scale(u)
        end

        function birth_c(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            b_drug = getproperty(p, B) - drug_interaction
            if b_drug >= 0.0
                return (max(u[idx], 0.0) * b_drug) * logistic_scale(u)
            else
                return 0.0
            end
        end

        function death_c(u, p, t)
            drug_interaction = birth_drug_adjustment(u, p)
            b_drug = getproperty(p, B) - drug_interaction
            if b_drug >= 0.0
                return (max(u[idx], 0.0) * getproperty(p, D)) * logistic_scale(u)
            else
                negative_surplus = abs(b_drug)
                return (max(u[idx], 0.0) * (getproperty(p, D) + negative_surplus)) * logistic_scale(u)
            end
        end

        birth = select_drug_effect(de, birth_d, birth_b, birth_c)
        death = select_drug_effect(de, death_d, death_b, death_c)
        return birth, death
    end

    S_division_activity, S_death_activity = build_birth_death_rates(
        de, RESPOP_NS_INDEX, bS, dS, Val(:bS), Val(:dS))
    R_division_activity, R_death_activity = build_birth_death_rates(
        de, RESPOP_NR_INDEX, bR, dR, Val(:bR), Val(:dR); psi_fn = p -> p.psi)
    E_symmetric_birth_activity, E_death_activity = build_birth_death_rates(
        de, RESPOP_NE_INDEX, bE, dE, Val(:bE), Val(:dE); psi_fn = p -> p.psi)

    S_symmetric_birth_activity(u, p, t) = S_division_activity(u, p, t) * (independent ? 1.0 : 1 - p.mu)
    S_asymmetric_R_birth_activity(u, p, t) = S_division_activity(u, p, t) * p.mu
    R_symmetric_birth_activity(u, p, t) = R_division_activity(u, p, t) *
        (independent ? 1.0 : max(1 - p.sig - p.al * clamp(u[RESPOP_GAM_INDEX], 0.0, 1.0), 0.0))
    R_asymmetric_S_birth_activity(u, p, t) = R_division_activity(u, p, t) * p.sig
    R_asymmetric_E_birth_activity(u, p, t) = R_division_activity(u, p, t) *
        p.al * clamp(u[RESPOP_GAM_INDEX], 0.0, 1.0)

    S_to_R_transition_activity(u, p, t) = max(u[RESPOP_NS_INDEX], 0.0) * p.mu
    R_to_S_transition_activity(u, p, t) = max(u[RESPOP_NR_INDEX], 0.0) * p.sig
    R_to_E_transition_activity(u, p, t) = max(u[RESPOP_NR_INDEX], 0.0) * p.al * clamp(u[RESPOP_GAM_INDEX], 0.0, 1.0)

    # ODE equations remain explicit by phenotype; each term is a complete event.
    function ode_fxn!(du, u, p, t)
        du.gam = p.kp
        if independent
            du.nS = (p.S_symmetric_birth_ode * S_symmetric_birth_activity(u, p, t)
                - p.S_death_ode * S_death_activity(u, p, t)
                - p.S_to_R_transition_ode * S_to_R_transition_activity(u, p, t)
                + p.R_to_S_transition_ode * R_to_S_transition_activity(u, p, t))
            du.nR = (p.R_symmetric_birth_ode * R_symmetric_birth_activity(u, p, t)
                - p.R_death_ode * R_death_activity(u, p, t)
                + p.S_to_R_transition_ode * S_to_R_transition_activity(u, p, t)
                - p.R_to_S_transition_ode * R_to_S_transition_activity(u, p, t)
                - p.R_to_E_transition_ode * R_to_E_transition_activity(u, p, t))
            du.nE = (p.E_symmetric_birth_ode * E_symmetric_birth_activity(u, p, t)
                - p.E_death_ode * E_death_activity(u, p, t)
                + p.R_to_E_transition_ode * R_to_E_transition_activity(u, p, t))
        else
            du.nS = (p.S_symmetric_birth_ode * S_symmetric_birth_activity(u, p, t)
                - p.S_death_ode * S_death_activity(u, p, t)
                + p.R_asymmetric_S_birth_ode * R_asymmetric_S_birth_activity(u, p, t))
            du.nR = (p.S_asymmetric_R_birth_ode * S_asymmetric_R_birth_activity(u, p, t)
                + p.R_symmetric_birth_ode * R_symmetric_birth_activity(u, p, t)
                - p.R_death_ode * R_death_activity(u, p, t))
            du.nE = (p.R_asymmetric_E_birth_ode * R_asymmetric_E_birth_activity(u, p, t)
                + p.E_symmetric_birth_ode * E_symmetric_birth_activity(u, p, t)
                - p.E_death_ode * E_death_activity(u, p, t))
        end
        du.Pass_num = 0.0
        return nothing
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    ################################
    # Jump Problem: Division Outcomes, Death, Damage and Repair
    ################################

    function S_symmetric_birth!(integrator)
        integrator.u[RESPOP_NS_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NS_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_symmetric_birth_rate(u, p, t)
        p.S_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_NS_INDEX] > 0.0 || return 0.0
        return S_symmetric_birth_activity(u, p, t)
    end
    S_symmetric_birth_jump = VariableRateJump(S_symmetric_birth_rate, S_symmetric_birth!)

    function S_death!(integrator)
        integrator.u[RESPOP_NS_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NS_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_death_rate(u, p, t)
        p.S_death_ode == 0.0 || return 0.0
        u[RESPOP_NS_INDEX] >= 1.0 || return 0.0
        return S_death_activity(u, p, t)
    end
    S_death_jump = VariableRateJump(S_death_rate, S_death!)

    function R_symmetric_birth!(integrator)
        integrator.u[RESPOP_NR_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NR_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_symmetric_birth_rate(u, p, t)
        p.R_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] > 0.0 || return 0.0
        return R_symmetric_birth_activity(u, p, t)
    end
    R_symmetric_birth_jump = VariableRateJump(R_symmetric_birth_rate, R_symmetric_birth!)

    function R_death!(integrator)
        integrator.u[RESPOP_NR_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NR_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_death_rate(u, p, t)
        p.R_death_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] >= 1.0 || return 0.0
        return R_death_activity(u, p, t)
    end
    R_death_jump = VariableRateJump(R_death_rate, R_death!)

    function E_symmetric_birth!(integrator)
        integrator.u[RESPOP_NE_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NE_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_symmetric_birth_rate(u, p, t)
        p.E_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_NE_INDEX] > 0.0 || return 0.0
        return E_symmetric_birth_activity(u, p, t)
    end
    E_symmetric_birth_jump = VariableRateJump(E_symmetric_birth_rate, E_symmetric_birth!)

    function E_death!(integrator)
        integrator.u[RESPOP_NE_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NE_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_death_rate(u, p, t)
        p.E_death_ode == 0.0 || return 0.0
        u[RESPOP_NE_INDEX] >= 1.0 || return 0.0
        return E_death_activity(u, p, t)
    end
    E_death_jump = VariableRateJump(E_death_rate, E_death!)

    # Asymmetric daughter births.

    function S_asymmetric_R_birth!(integrator)
        integrator.u[RESPOP_NS_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NR_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_asymmetric_R_birth_rate(u, p, t)
        p.S_asymmetric_R_birth_ode == 0.0 || return 0.0
        u[RESPOP_NS_INDEX] > 0.0 || return 0.0
        return S_asymmetric_R_birth_activity(u, p, t)
    end
    function S_to_R_transition!(integrator)
        integrator.u[RESPOP_NS_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NS_INDEX] -= 1
        integrator.u[RESPOP_NR_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_to_R_transition_rate(u, p, t)
        p.S_to_R_transition_ode == 0.0 || return 0.0
        u[RESPOP_NS_INDEX] >= 1.0 || return 0.0
        return S_to_R_transition_activity(u, p, t)
    end
    S_to_R_transition_jump = independent ? VariableRateJump(S_to_R_transition_rate, S_to_R_transition!) : nothing
    S_asymmetric_R_birth_jump = independent ? nothing : VariableRateJump(S_asymmetric_R_birth_rate, S_asymmetric_R_birth!)

    function R_asymmetric_S_birth!(integrator)
        integrator.u[RESPOP_NR_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NS_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_asymmetric_S_birth_rate(u, p, t)
        p.R_asymmetric_S_birth_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] > 0.0 || return 0.0
        return R_asymmetric_S_birth_activity(u, p, t)
    end
    function R_to_S_transition!(integrator)
        integrator.u[RESPOP_NR_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NR_INDEX] -= 1
        integrator.u[RESPOP_NS_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_to_S_transition_rate(u, p, t)
        p.R_to_S_transition_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] >= 1.0 || return 0.0
        return R_to_S_transition_activity(u, p, t)
    end
    R_to_S_transition_jump = independent ? VariableRateJump(R_to_S_transition_rate, R_to_S_transition!) : nothing
    R_asymmetric_S_birth_jump = independent ? nothing : VariableRateJump(R_asymmetric_S_birth_rate, R_asymmetric_S_birth!)

    function R_asymmetric_E_birth!(integrator)
        integrator.u[RESPOP_NR_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_NE_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_asymmetric_E_birth_rate(u, p, t)
        p.R_asymmetric_E_birth_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] > 0.0 || return 0.0
        return R_asymmetric_E_birth_activity(u, p, t)
    end
    function R_to_E_transition!(integrator)
        integrator.u[RESPOP_NR_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_NR_INDEX] -= 1
        integrator.u[RESPOP_NE_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_to_E_transition_rate(u, p, t)
        p.R_to_E_transition_ode == 0.0 || return 0.0
        u[RESPOP_NR_INDEX] >= 1.0 || return 0.0
        return R_to_E_transition_activity(u, p, t)
    end
    R_to_E_transition_jump = independent ? VariableRateJump(R_to_E_transition_rate, R_to_E_transition!) : nothing
    R_asymmetric_E_birth_jump = independent ? nothing : VariableRateJump(R_asymmetric_E_birth_rate, R_asymmetric_E_birth!)

    ###########
    # Callbacks
    ###########

    n0 = sim_eff.n0

    # Decide each complete event from its source and full biological activity.
    event_flags = independent ? (
            :S_symmetric_birth_ode, :S_to_R_transition_ode, :S_death_ode,
            :R_symmetric_birth_ode, :R_to_S_transition_ode, :R_to_E_transition_ode,
            :R_death_ode, :E_symmetric_birth_ode, :E_death_ode
        ) : (
            :S_symmetric_birth_ode, :S_asymmetric_R_birth_ode, :S_death_ode,
            :R_symmetric_birth_ode, :R_asymmetric_S_birth_ode, :R_asymmetric_E_birth_ode,
            :R_death_ode, :E_symmetric_birth_ode, :E_death_ode
        )

    function classify_event_modes!(u, p, t)
        previous = Tuple(getproperty(p, name) for name in event_flags)
        p.S_symmetric_birth_ode = u[RESPOP_NS_INDEX] >= sim_eff.Nswitch &&
            S_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.S_to_R_transition_ode = u[RESPOP_NS_INDEX] >= sim_eff.Nswitch &&
                S_to_R_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.S_asymmetric_R_birth_ode = u[RESPOP_NS_INDEX] >= sim_eff.Nswitch &&
                S_asymmetric_R_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.S_death_ode = u[RESPOP_NS_INDEX] >= sim_eff.Nswitch &&
            S_death_activity(u, p, t) >= sim_eff.Aswitch
        p.R_symmetric_birth_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
            R_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.R_to_S_transition_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
                R_to_S_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_asymmetric_S_birth_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
                R_asymmetric_S_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        if independent
            p.R_to_E_transition_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
                R_to_E_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_asymmetric_E_birth_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
                R_asymmetric_E_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.R_death_ode = u[RESPOP_NR_INDEX] >= sim_eff.Nswitch &&
            R_death_activity(u, p, t) >= sim_eff.Aswitch
        p.E_symmetric_birth_ode = u[RESPOP_NE_INDEX] >= sim_eff.Nswitch &&
            E_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        p.E_death_ode = u[RESPOP_NE_INDEX] >= sim_eff.Nswitch &&
            E_death_activity(u, p, t) >= sim_eff.Aswitch
        return previous != Tuple(getproperty(p, name) for name in event_flags)
    end

    function refresh_event_modes!(u, p, t)
        changed = classify_event_modes!(u, p, t)
        # Evaluate all eligibility before changing any count; net drift is not sufficient.
        NS_has_ode = (p.S_symmetric_birth_ode == 1.0 && S_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.S_death_ode == 1.0 && S_death_activity(u, p, t) > 0.0)
        NR_has_ode = (p.R_symmetric_birth_ode == 1.0 && R_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.R_death_ode == 1.0 && R_death_activity(u, p, t) > 0.0)
        NE_has_ode = (p.E_symmetric_birth_ode == 1.0 && E_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.E_death_ode == 1.0 && E_death_activity(u, p, t) > 0.0)
        if independent
            S_to_R_transition_active = p.S_to_R_transition_ode == 1.0 && S_to_R_transition_activity(u, p, t) > 0.0
            R_to_S_transition_active = p.R_to_S_transition_ode == 1.0 && R_to_S_transition_activity(u, p, t) > 0.0
            R_to_E_transition_active = p.R_to_E_transition_ode == 1.0 && R_to_E_transition_activity(u, p, t) > 0.0
            NS_has_ode |= S_to_R_transition_active || R_to_S_transition_active
            NR_has_ode |= S_to_R_transition_active || R_to_S_transition_active || R_to_E_transition_active
            NE_has_ode |= R_to_E_transition_active
        else
            NS_has_ode |= p.R_asymmetric_S_birth_ode == 1.0 && R_asymmetric_S_birth_activity(u, p, t) > 0.0
            NR_has_ode |= p.S_asymmetric_R_birth_ode == 1.0 && S_asymmetric_R_birth_activity(u, p, t) > 0.0
            NE_has_ode |= p.R_asymmetric_E_birth_ode == 1.0 && R_asymmetric_E_birth_activity(u, p, t) > 0.0
        end
        projected = false
        if !NS_has_ode && 0.0 <= u[RESPOP_NS_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_NS_INDEX])
            lower = floor(u[RESPOP_NS_INDEX])
            u[RESPOP_NS_INDEX] = lower + (rand(rng) < u[RESPOP_NS_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NR_has_ode && 0.0 <= u[RESPOP_NR_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_NR_INDEX])
            lower = floor(u[RESPOP_NR_INDEX])
            u[RESPOP_NR_INDEX] = lower + (rand(rng) < u[RESPOP_NR_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NE_has_ode && 0.0 <= u[RESPOP_NE_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_NE_INDEX])
            lower = floor(u[RESPOP_NE_INDEX])
            u[RESPOP_NE_INDEX] = lower + (rand(rng) < u[RESPOP_NE_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        projected && classify_event_modes!(u, p, t)
        return changed || projected
    end

    function refresh_integrator!(integrator; modified = false)
        changed = refresh_event_modes!(integrator.u, integrator.p, integrator.t)
        if changed || modified
            u_modified!(integrator, true)
            reset_aggregated_jumps!(integrator; update_jump_params = false)
        end
        return nothing
    end

    # Resample populations for a passage bottleneck (mutates state; uses RNG).
    function resample_passage_bottleneck!(integrator, n0)
        curr_nR = round(Int, integrator.u[RESPOP_NR_INDEX])
        curr_nS = round(Int, integrator.u[RESPOP_NS_INDEX])
        curr_nE = round(Int, integrator.u[RESPOP_NE_INDEX])

        counts = multivariate_hypergeometric_draw([curr_nR, curr_nS, curr_nE], n0)

        integrator.u[RESPOP_NR_INDEX] = Float64(counts[1])
        integrator.u[RESPOP_NS_INDEX] = Float64(counts[2])
        integrator.u[RESPOP_NE_INDEX] = Float64(counts[3])
        nothing
    end

    function attempt_passage!(integrator, n0; require_enough::Bool)
        if require_enough
            curr_nR = round(Int, integrator.u[RESPOP_NR_INDEX])
            curr_nS = round(Int, integrator.u[RESPOP_NS_INDEX])
            curr_nE = round(Int, integrator.u[RESPOP_NE_INDEX])
            if (curr_nR + curr_nS + curr_nE) < n0
                return false
            end
        end

        resample_passage_bottleneck!(integrator, n0)
        integrator.u[RESPOP_PASS_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return true
    end

    # Non-negativity clamp callback for a population component.
    function build_clamp_callback(idx)
        condition(u, t, integrator) = u[idx] < 0
        function affect!(integrator)
            integrator.u[idx] = 0.0
            refresh_integrator!(integrator; modified = true)
        end
        return DiscreteCallback(condition, affect!,
                                save_positions = (false, true))
    end

    event_mode_cb = DiscreteCallback((u, t, integrator) -> true,
        integrator -> refresh_integrator!(integrator); save_positions = (false, false))

    # Treatment schedule callbacks (use preset times for robustness).
    function set_treatment_on!(integrator)
        integrator.p.kp = params.k
        refresh_integrator!(integrator; modified = true)
    end

    function set_treatment_off!(integrator)
        integrator.p.kp = -params.k
        refresh_integrator!(integrator; modified = true)
    end

    filter_times_in_span(times) = sort(unique(filter(t -> (t >= sim_eff.t0 && t <= sim_eff.tmax), times)))

    treat_on_times = filter_times_in_span(sim_eff.treat_ons)
    treat_off_times = filter_times_in_span(sim_eff.treat_offs)

    treat_on_cb = isempty(treat_on_times) ? nothing :
                  PresetTimeCallback(treat_on_times, set_treatment_on!; save_positions = (false, true))
    treat_off_cb = isempty(treat_off_times) ? nothing :
                   PresetTimeCallback(treat_off_times, set_treatment_off!; save_positions = (false, true))

    # Drug level clamping
    gam_above_max(u, t, integrator) = (u[RESPOP_GAM_INDEX] - 1.0)
    function clamp_gam_max!(integrator)
        integrator.u[RESPOP_GAM_INDEX] = 1.0
        integrator.p.kp = 0.0
        refresh_integrator!(integrator; modified = true)
    end

    gam_below_min(u, t, integrator) = (u[RESPOP_GAM_INDEX])
    function clamp_gam_min!(integrator)
        integrator.u[RESPOP_GAM_INDEX] = 0.0
        integrator.p.kp = 0.0
        refresh_integrator!(integrator; modified = true)
    end

    gam_max_cb = ContinuousCallback(gam_above_max, clamp_gam_max!,
                                    save_positions = (false, true))
    gam_min_cb = ContinuousCallback(gam_below_min, clamp_gam_min!,
                                    save_positions = (false, true))

    # Nmax, Extinction & Passage
    nmax_reached(u, t, integrator) = (respop_total_population(integrator.u) >= sim_eff.Nmax)

    function handle_nmax!(integrator)
        if integrator.u[RESPOP_PASS_INDEX] >= n_pass_eff
            terminate!(integrator)
        else
            attempt_passage!(integrator, n0; require_enough = false)
        end
    end

    nmax_cb = DiscreteCallback(nmax_reached, handle_nmax!,
                               save_positions = (false, true))

    extinction(u, t, integrator) = (respop_total_population(integrator.u) < 1.0)
    extinction_cb = DiscreteCallback(extinction, terminate!,
                                     save_positions = (false, false))

    function apply_scheduled_passage!(integrator)
        curr_pass = round(Int, integrator.u[RESPOP_PASS_INDEX])
        if curr_pass < n_pass_eff
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
    cb7a = build_clamp_callback(RESPOP_NS_INDEX)
    cb7b = build_clamp_callback(RESPOP_NR_INDEX)
    cb7c = build_clamp_callback(RESPOP_NE_INDEX)

    # Collect Callbacks
    cb_list = Any[
        event_mode_cb, nmax_cb, extinction_cb,
        cb7a, cb7b, cb7c
    ]

    if passage_cb !== nothing
        push!(cb_list, passage_cb)
    end

    if sim_eff.treat == true
        if treat_on_cb !== nothing
            push!(cb_list, treat_on_cb)
        end
        if treat_off_cb !== nothing
            push!(cb_list, treat_off_cb)
        end
        push!(cb_list, gam_max_cb, gam_min_cb)
    end

    refresh_event_modes!(u0, p, sim_eff.t0)
    if _setup_only
        # Internal test access to these model-local scientific definitions.
        return (p = p, u0 = u0, ode! = ode_fxn!,
            activities = (
                S_symmetric_birth_activity, (independent ? S_to_R_transition_activity : S_asymmetric_R_birth_activity), S_death_activity,
                R_symmetric_birth_activity, (independent ? R_to_S_transition_activity : R_asymmetric_S_birth_activity), (independent ? R_to_E_transition_activity : R_asymmetric_E_birth_activity),
                R_death_activity, E_symmetric_birth_activity, E_death_activity
            ),
            jump_rates = (
                S_symmetric_birth_rate, (independent ? S_to_R_transition_rate : S_asymmetric_R_birth_rate), S_death_rate,
                R_symmetric_birth_rate, (independent ? R_to_S_transition_rate : R_asymmetric_S_birth_rate), (independent ? R_to_E_transition_rate : R_asymmetric_E_birth_rate),
                R_death_rate, E_symmetric_birth_rate, E_death_rate
            ),
            flags = event_flags,
            sources = (
                RESPOP_NS_INDEX, RESPOP_NS_INDEX, RESPOP_NS_INDEX,
                RESPOP_NR_INDEX, RESPOP_NR_INDEX, RESPOP_NR_INDEX,
                RESPOP_NR_INDEX, RESPOP_NE_INDEX, RESPOP_NE_INDEX
            ),
            classify! = classify_event_modes!, refresh! = refresh_event_modes!)
    end

    cbs = CallbackSet(cb_list...)

    sjm_prob = JumpProblem(prob, Direct(),
        S_symmetric_birth_jump,
        (independent ? S_to_R_transition_jump : S_asymmetric_R_birth_jump),
        S_death_jump,
        R_symmetric_birth_jump,
        (independent ? R_to_S_transition_jump : R_asymmetric_S_birth_jump),
        (independent ? R_to_E_transition_jump : R_asymmetric_E_birth_jump),
        R_death_jump,
        E_symmetric_birth_jump,
        E_death_jump; rng = rng)

    event_times = passage_times
    if sim_eff.treat == true
        event_times = vcat(event_times, treat_on_times, treat_off_times)
    end
    tstops = sort(unique(event_times))

    sol = if get(ENV, "BARCODEMODELS_SUPPRESS_SOLVER", "true") == "false"
        solve(sjm_prob, Tsit5(); callback = cbs,
              tstops = tstops, saveat = sim_eff.save_at, solve_kwargs...)
    else
        @suppress solve(sjm_prob, Tsit5(); callback = cbs,
                        tstops = tstops, saveat = sim_eff.save_at, solve_kwargs...)
    end

    return sol
end

simulate_grow_kill(model::ResPop, state::ResPopState, sim::SimParams) =
    run_model_core_hybrid(model, state, sim; treat = sim.treat)

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
    Aswitch::Float64 = 1000.0,
    drug_effect::Union{String, Symbol} = "d",
    transition_coupling::Union{String, Symbol} = :birth_coupled)

    model = ResPop(ResPopParams(
        b = b, d = d, mu = mu, sig = sig, del = del, al = al,
        Dc = Dc, k = k, psi = psi, drug_effect = drug_effect,
        transition_coupling = transition_coupling
    ))

    state = ResPopState(nS, nR, nE; gam = 0.0, pass_num = 1)

    sim = SimParams(
        n0 = n0, t0 = t0, tmax = tmax, t_Pass = t_Pass,
        Nmax = Nmax, Cc = Cc, Nswitch = Nswitch,
        treat_ons = treat_ons, treat_offs = treat_offs,
        save_at = save_at, treat = treat, Aswitch = Aswitch
    )

    return run_model_core_hybrid(model, state, sim; treat = treat)
end
