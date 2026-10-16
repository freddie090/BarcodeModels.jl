"""Resistance population model with a fixed heritable engraftment phenotype."""
struct ResPopInVivo <: HybridModel
    params::ResPopInVivoParams
    function ResPopInVivo(params::ResPopInVivoParams)
        validate_model_params(params)
        new(params)
    end
end

"""Construct a `ResPopInVivo` from keyword arguments forwarded to `ResPopInVivoParams`."""
ResPopInVivo(; kwargs...) = ResPopInVivo(ResPopInVivoParams(; kwargs...))

function build_component_params(params::ResPopInVivoParams; R_real = "l")
    R_real = normalize_R_real(R_real)
    bS = params.b
    dS = params.d

    bR, dR = if R_real == "b"
        (params.b - (params.b - params.d) * params.del, params.d)
    elseif R_real == "d"
        (params.b, params.d + (params.b - params.d) * params.del)
    else
        (params.b * (1 - params.del), params.d * (1 - params.del))
    end

    bE = params.b
    dE = params.d

    phenotype_modes_EG0 = params.transition_coupling === :birth_coupled ? (
        S_EG0_asymmetric_R_birth_ode = 0.0,
        R_EG0_asymmetric_S_birth_ode = 0.0,
        R_EG0_asymmetric_E_birth_ode = 0.0,
    ) : (
        S_EG0_to_R_transition_ode = 0.0,
        R_EG0_to_S_transition_ode = 0.0,
        R_EG0_to_E_transition_ode = 0.0,
    )

    phenotype_modes_EG1 = params.transition_coupling === :birth_coupled ? (
        S_EG1_asymmetric_R_birth_ode = 0.0,
        R_EG1_asymmetric_S_birth_ode = 0.0,
        R_EG1_asymmetric_E_birth_ode = 0.0,
    ) : (
        S_EG1_to_R_transition_ode = 0.0,
        R_EG1_to_S_transition_ode = 0.0,
        R_EG1_to_E_transition_ode = 0.0,
    )

    p = ComponentArray(;
        bS = bS, dS = dS, bR = bR, dR = dR, bE = bE, dE = dE,
        Dc = params.Dc, kp = 0.0, psi = params.psi,
        mu = params.mu, sig = params.sig, al = params.al,
        S_EG0_symmetric_birth_ode = 0.0,
        S_EG0_death_ode = 0.0,
        R_EG0_symmetric_birth_ode = 0.0,
        R_EG0_death_ode = 0.0,
        E_EG0_symmetric_birth_ode = 0.0,
        E_EG0_death_ode = 0.0,
        S_EG1_symmetric_birth_ode = 0.0,
        S_EG1_death_ode = 0.0,
        R_EG1_symmetric_birth_ode = 0.0,
        R_EG1_death_ode = 0.0,
        E_EG1_symmetric_birth_ode = 0.0,
        E_EG1_death_ode = 0.0,
        phenotype_modes_EG0..., phenotype_modes_EG1...
    )

    return p, bS, dS, bR, dR, bE, dE
end

function run_model_core_hybrid(model::ResPopInVivo, state::ResPopInVivoState, sim::SimParams;
    treat::Bool = sim.treat, rng = Random.default_rng(), _setup_only::Bool = false,
    on_engraft = (passage, stats) -> nothing, solve_kwargs...)
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
        Aswitch = sim.Aswitch,
        R_real = sim.R_real
    )

    params = model.params
    independent = params.transition_coupling === :independent_rate
    de = params.drug_effect

    validate_resistant_birth_drug_effect(params, sim_eff.R_real)
    p, bS, dS, bR, dR, bE, dE = build_component_params(params; R_real = sim_eff.R_real)

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
    function build_birth_death_rates(de, idx, ::Val{B}, ::Val{D}; psi_fn = p -> 0.0) where {B, D}
        function logistic_scale(u)
            N = respop_invivo_total_population(u)
            return max(logistic_factor(N, sim_eff.Cc), 0.0)
        end

        function birth_drug_adjustment(u, p)
            return p.Dc * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0) * (1 - psi_fn(p))
        end

        function death_drug_adjustment(u, p)
            return p.Dc * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0) * (1 - psi_fn(p))
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

    S_EG0_division_activity, S_EG0_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NS_EG0_INDEX, Val(:bS), Val(:dS))
    R_EG0_division_activity, R_EG0_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NR_EG0_INDEX, Val(:bR), Val(:dR); psi_fn = p -> p.psi)
    E_EG0_symmetric_birth_activity, E_EG0_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NE_EG0_INDEX, Val(:bE), Val(:dE); psi_fn = p -> p.psi)

    S_EG0_symmetric_birth_activity(u, p, t) = S_EG0_division_activity(u, p, t) * (independent ? 1.0 : 1 - p.mu)
    S_EG0_asymmetric_R_birth_activity(u, p, t) = S_EG0_division_activity(u, p, t) * p.mu
    R_EG0_symmetric_birth_activity(u, p, t) = R_EG0_division_activity(u, p, t) *
        (independent ? 1.0 : max(1 - p.sig - p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0), 0.0))
    R_EG0_asymmetric_S_birth_activity(u, p, t) = R_EG0_division_activity(u, p, t) * p.sig
    R_EG0_asymmetric_E_birth_activity(u, p, t) = R_EG0_division_activity(u, p, t) *
        p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0)

    S_EG0_to_R_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NS_EG0_INDEX], 0.0) * p.mu
    R_EG0_to_S_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NR_EG0_INDEX], 0.0) * p.sig
    R_EG0_to_E_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NR_EG0_INDEX], 0.0) * p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0)


    S_EG1_division_activity, S_EG1_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NS_EG1_INDEX, Val(:bS), Val(:dS))
    R_EG1_division_activity, R_EG1_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NR_EG1_INDEX, Val(:bR), Val(:dR); psi_fn = p -> p.psi)
    E_EG1_symmetric_birth_activity, E_EG1_death_activity = build_birth_death_rates(
        de, RESPOP_INVIVO_NE_EG1_INDEX, Val(:bE), Val(:dE); psi_fn = p -> p.psi)

    S_EG1_symmetric_birth_activity(u, p, t) = S_EG1_division_activity(u, p, t) * (independent ? 1.0 : 1 - p.mu)
    S_EG1_asymmetric_R_birth_activity(u, p, t) = S_EG1_division_activity(u, p, t) * p.mu
    R_EG1_symmetric_birth_activity(u, p, t) = R_EG1_division_activity(u, p, t) *
        (independent ? 1.0 : max(1 - p.sig - p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0), 0.0))
    R_EG1_asymmetric_S_birth_activity(u, p, t) = R_EG1_division_activity(u, p, t) * p.sig
    R_EG1_asymmetric_E_birth_activity(u, p, t) = R_EG1_division_activity(u, p, t) *
        p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0)

    S_EG1_to_R_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NS_EG1_INDEX], 0.0) * p.mu
    R_EG1_to_S_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NR_EG1_INDEX], 0.0) * p.sig
    R_EG1_to_E_transition_activity(u, p, t) = max(u[RESPOP_INVIVO_NR_EG1_INDEX], 0.0) * p.al * clamp(u[RESPOP_INVIVO_GAM_INDEX], 0.0, 1.0)

    # ODE equations remain explicit by phenotype; each term is a complete event.
    function ode_fxn!(du, u, p, t)
        du.gam = p.kp
        if independent
            du.nS_EG0 = (p.S_EG0_symmetric_birth_ode * S_EG0_symmetric_birth_activity(u, p, t)
                - p.S_EG0_death_ode * S_EG0_death_activity(u, p, t)
                - p.S_EG0_to_R_transition_ode * S_EG0_to_R_transition_activity(u, p, t)
                + p.R_EG0_to_S_transition_ode * R_EG0_to_S_transition_activity(u, p, t))
            du.nR_EG0 = (p.R_EG0_symmetric_birth_ode * R_EG0_symmetric_birth_activity(u, p, t)
                - p.R_EG0_death_ode * R_EG0_death_activity(u, p, t)
                + p.S_EG0_to_R_transition_ode * S_EG0_to_R_transition_activity(u, p, t)
                - p.R_EG0_to_S_transition_ode * R_EG0_to_S_transition_activity(u, p, t)
                - p.R_EG0_to_E_transition_ode * R_EG0_to_E_transition_activity(u, p, t))
            du.nE_EG0 = (p.E_EG0_symmetric_birth_ode * E_EG0_symmetric_birth_activity(u, p, t)
                - p.E_EG0_death_ode * E_EG0_death_activity(u, p, t)
                + p.R_EG0_to_E_transition_ode * R_EG0_to_E_transition_activity(u, p, t))
            du.nS_EG1 = (p.S_EG1_symmetric_birth_ode * S_EG1_symmetric_birth_activity(u, p, t)
                - p.S_EG1_death_ode * S_EG1_death_activity(u, p, t)
                - p.S_EG1_to_R_transition_ode * S_EG1_to_R_transition_activity(u, p, t)
                + p.R_EG1_to_S_transition_ode * R_EG1_to_S_transition_activity(u, p, t))
            du.nR_EG1 = (p.R_EG1_symmetric_birth_ode * R_EG1_symmetric_birth_activity(u, p, t)
                - p.R_EG1_death_ode * R_EG1_death_activity(u, p, t)
                + p.S_EG1_to_R_transition_ode * S_EG1_to_R_transition_activity(u, p, t)
                - p.R_EG1_to_S_transition_ode * R_EG1_to_S_transition_activity(u, p, t)
                - p.R_EG1_to_E_transition_ode * R_EG1_to_E_transition_activity(u, p, t))
            du.nE_EG1 = (p.E_EG1_symmetric_birth_ode * E_EG1_symmetric_birth_activity(u, p, t)
                - p.E_EG1_death_ode * E_EG1_death_activity(u, p, t)
                + p.R_EG1_to_E_transition_ode * R_EG1_to_E_transition_activity(u, p, t))
        else
            du.nS_EG0 = (p.S_EG0_symmetric_birth_ode * S_EG0_symmetric_birth_activity(u, p, t)
                - p.S_EG0_death_ode * S_EG0_death_activity(u, p, t)
                + p.R_EG0_asymmetric_S_birth_ode * R_EG0_asymmetric_S_birth_activity(u, p, t))
            du.nR_EG0 = (p.S_EG0_asymmetric_R_birth_ode * S_EG0_asymmetric_R_birth_activity(u, p, t)
                + p.R_EG0_symmetric_birth_ode * R_EG0_symmetric_birth_activity(u, p, t)
                - p.R_EG0_death_ode * R_EG0_death_activity(u, p, t))
            du.nE_EG0 = (p.R_EG0_asymmetric_E_birth_ode * R_EG0_asymmetric_E_birth_activity(u, p, t)
                + p.E_EG0_symmetric_birth_ode * E_EG0_symmetric_birth_activity(u, p, t)
                - p.E_EG0_death_ode * E_EG0_death_activity(u, p, t))
            du.nS_EG1 = (p.S_EG1_symmetric_birth_ode * S_EG1_symmetric_birth_activity(u, p, t)
                - p.S_EG1_death_ode * S_EG1_death_activity(u, p, t)
                + p.R_EG1_asymmetric_S_birth_ode * R_EG1_asymmetric_S_birth_activity(u, p, t))
            du.nR_EG1 = (p.S_EG1_asymmetric_R_birth_ode * S_EG1_asymmetric_R_birth_activity(u, p, t)
                + p.R_EG1_symmetric_birth_ode * R_EG1_symmetric_birth_activity(u, p, t)
                - p.R_EG1_death_ode * R_EG1_death_activity(u, p, t))
            du.nE_EG1 = (p.R_EG1_asymmetric_E_birth_ode * R_EG1_asymmetric_E_birth_activity(u, p, t)
                + p.E_EG1_symmetric_birth_ode * E_EG1_symmetric_birth_activity(u, p, t)
                - p.E_EG1_death_ode * E_EG1_death_activity(u, p, t))
        end
        du.Pass_num = 0.0
        return nothing
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    ################################
    # Jump Problem: Division Outcomes, Death and Independent Conversion
    ################################

    function S_EG0_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG0_symmetric_birth_rate(u, p, t)
        p.S_EG0_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG0_INDEX] > 0.0 || return 0.0
        return S_EG0_symmetric_birth_activity(u, p, t)
    end
    S_EG0_symmetric_birth_jump = VariableRateJump(S_EG0_symmetric_birth_rate, S_EG0_symmetric_birth!)

    function S_EG0_death!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG0_death_rate(u, p, t)
        p.S_EG0_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG0_INDEX] >= 1.0 || return 0.0
        return S_EG0_death_activity(u, p, t)
    end
    S_EG0_death_jump = VariableRateJump(S_EG0_death_rate, S_EG0_death!)

    function R_EG0_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_symmetric_birth_rate(u, p, t)
        p.R_EG0_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return 0.0
        return R_EG0_symmetric_birth_activity(u, p, t)
    end
    R_EG0_symmetric_birth_jump = VariableRateJump(R_EG0_symmetric_birth_rate, R_EG0_symmetric_birth!)

    function R_EG0_death!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_death_rate(u, p, t)
        p.R_EG0_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return 0.0
        return R_EG0_death_activity(u, p, t)
    end
    R_EG0_death_jump = VariableRateJump(R_EG0_death_rate, R_EG0_death!)

    function E_EG0_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_EG0_symmetric_birth_rate(u, p, t)
        p.E_EG0_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NE_EG0_INDEX] > 0.0 || return 0.0
        return E_EG0_symmetric_birth_activity(u, p, t)
    end
    E_EG0_symmetric_birth_jump = VariableRateJump(E_EG0_symmetric_birth_rate, E_EG0_symmetric_birth!)

    function E_EG0_death!(integrator)
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_EG0_death_rate(u, p, t)
        p.E_EG0_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NE_EG0_INDEX] >= 1.0 || return 0.0
        return E_EG0_death_activity(u, p, t)
    end
    E_EG0_death_jump = VariableRateJump(E_EG0_death_rate, E_EG0_death!)

    # Asymmetric daughter births.

    function S_EG0_asymmetric_R_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG0_asymmetric_R_birth_rate(u, p, t)
        p.S_EG0_asymmetric_R_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG0_INDEX] > 0.0 || return 0.0
        return S_EG0_asymmetric_R_birth_activity(u, p, t)
    end
    function S_EG0_to_R_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG0_to_R_transition_rate(u, p, t)
        p.S_EG0_to_R_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG0_INDEX] >= 1.0 || return 0.0
        return S_EG0_to_R_transition_activity(u, p, t)
    end
    S_EG0_to_R_transition_jump = independent ? VariableRateJump(S_EG0_to_R_transition_rate, S_EG0_to_R_transition!) : nothing
    S_EG0_asymmetric_R_birth_jump = independent ? nothing : VariableRateJump(S_EG0_asymmetric_R_birth_rate, S_EG0_asymmetric_R_birth!)

    function R_EG0_asymmetric_S_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_asymmetric_S_birth_rate(u, p, t)
        p.R_EG0_asymmetric_S_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return 0.0
        return R_EG0_asymmetric_S_birth_activity(u, p, t)
    end
    function R_EG0_to_S_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NS_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_to_S_transition_rate(u, p, t)
        p.R_EG0_to_S_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return 0.0
        return R_EG0_to_S_transition_activity(u, p, t)
    end
    R_EG0_to_S_transition_jump = independent ? VariableRateJump(R_EG0_to_S_transition_rate, R_EG0_to_S_transition!) : nothing
    R_EG0_asymmetric_S_birth_jump = independent ? nothing : VariableRateJump(R_EG0_asymmetric_S_birth_rate, R_EG0_asymmetric_S_birth!)

    function R_EG0_asymmetric_E_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_asymmetric_E_birth_rate(u, p, t)
        p.R_EG0_asymmetric_E_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] > 0.0 || return 0.0
        return R_EG0_asymmetric_E_birth_activity(u, p, t)
    end
    function R_EG0_to_E_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG0_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NE_EG0_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG0_to_E_transition_rate(u, p, t)
        p.R_EG0_to_E_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG0_INDEX] >= 1.0 || return 0.0
        return R_EG0_to_E_transition_activity(u, p, t)
    end
    R_EG0_to_E_transition_jump = independent ? VariableRateJump(R_EG0_to_E_transition_rate, R_EG0_to_E_transition!) : nothing
    R_EG0_asymmetric_E_birth_jump = independent ? nothing : VariableRateJump(R_EG0_asymmetric_E_birth_rate, R_EG0_asymmetric_E_birth!)


    function S_EG1_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG1_symmetric_birth_rate(u, p, t)
        p.S_EG1_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG1_INDEX] > 0.0 || return 0.0
        return S_EG1_symmetric_birth_activity(u, p, t)
    end
    S_EG1_symmetric_birth_jump = VariableRateJump(S_EG1_symmetric_birth_rate, S_EG1_symmetric_birth!)

    function S_EG1_death!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG1_death_rate(u, p, t)
        p.S_EG1_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG1_INDEX] >= 1.0 || return 0.0
        return S_EG1_death_activity(u, p, t)
    end
    S_EG1_death_jump = VariableRateJump(S_EG1_death_rate, S_EG1_death!)

    function R_EG1_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_symmetric_birth_rate(u, p, t)
        p.R_EG1_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return 0.0
        return R_EG1_symmetric_birth_activity(u, p, t)
    end
    R_EG1_symmetric_birth_jump = VariableRateJump(R_EG1_symmetric_birth_rate, R_EG1_symmetric_birth!)

    function R_EG1_death!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_death_rate(u, p, t)
        p.R_EG1_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return 0.0
        return R_EG1_death_activity(u, p, t)
    end
    R_EG1_death_jump = VariableRateJump(R_EG1_death_rate, R_EG1_death!)

    function E_EG1_symmetric_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_EG1_symmetric_birth_rate(u, p, t)
        p.E_EG1_symmetric_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NE_EG1_INDEX] > 0.0 || return 0.0
        return E_EG1_symmetric_birth_activity(u, p, t)
    end
    E_EG1_symmetric_birth_jump = VariableRateJump(E_EG1_symmetric_birth_rate, E_EG1_symmetric_birth!)

    function E_EG1_death!(integrator)
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] -= 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function E_EG1_death_rate(u, p, t)
        p.E_EG1_death_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NE_EG1_INDEX] >= 1.0 || return 0.0
        return E_EG1_death_activity(u, p, t)
    end
    E_EG1_death_jump = VariableRateJump(E_EG1_death_rate, E_EG1_death!)

    # Asymmetric daughter births.

    function S_EG1_asymmetric_R_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG1_asymmetric_R_birth_rate(u, p, t)
        p.S_EG1_asymmetric_R_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG1_INDEX] > 0.0 || return 0.0
        return S_EG1_asymmetric_R_birth_activity(u, p, t)
    end
    function S_EG1_to_R_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function S_EG1_to_R_transition_rate(u, p, t)
        p.S_EG1_to_R_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NS_EG1_INDEX] >= 1.0 || return 0.0
        return S_EG1_to_R_transition_activity(u, p, t)
    end
    S_EG1_to_R_transition_jump = independent ? VariableRateJump(S_EG1_to_R_transition_rate, S_EG1_to_R_transition!) : nothing
    S_EG1_asymmetric_R_birth_jump = independent ? nothing : VariableRateJump(S_EG1_asymmetric_R_birth_rate, S_EG1_asymmetric_R_birth!)

    function R_EG1_asymmetric_S_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_asymmetric_S_birth_rate(u, p, t)
        p.R_EG1_asymmetric_S_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return 0.0
        return R_EG1_asymmetric_S_birth_activity(u, p, t)
    end
    function R_EG1_to_S_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NS_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_to_S_transition_rate(u, p, t)
        p.R_EG1_to_S_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return 0.0
        return R_EG1_to_S_transition_activity(u, p, t)
    end
    R_EG1_to_S_transition_jump = independent ? VariableRateJump(R_EG1_to_S_transition_rate, R_EG1_to_S_transition!) : nothing
    R_EG1_asymmetric_S_birth_jump = independent ? nothing : VariableRateJump(R_EG1_asymmetric_S_birth_rate, R_EG1_asymmetric_S_birth!)

    function R_EG1_asymmetric_E_birth!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return nothing
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_asymmetric_E_birth_rate(u, p, t)
        p.R_EG1_asymmetric_E_birth_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] > 0.0 || return 0.0
        return R_EG1_asymmetric_E_birth_activity(u, p, t)
    end
    function R_EG1_to_E_transition!(integrator)
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return nothing
        integrator.u[RESPOP_INVIVO_NR_EG1_INDEX] -= 1
        integrator.u[RESPOP_INVIVO_NE_EG1_INDEX] += 1
        refresh_integrator!(integrator; modified = true)
        return nothing
    end

    function R_EG1_to_E_transition_rate(u, p, t)
        p.R_EG1_to_E_transition_ode == 0.0 || return 0.0
        u[RESPOP_INVIVO_NR_EG1_INDEX] >= 1.0 || return 0.0
        return R_EG1_to_E_transition_activity(u, p, t)
    end
    R_EG1_to_E_transition_jump = independent ? VariableRateJump(R_EG1_to_E_transition_rate, R_EG1_to_E_transition!) : nothing
    R_EG1_asymmetric_E_birth_jump = independent ? nothing : VariableRateJump(R_EG1_asymmetric_E_birth_rate, R_EG1_asymmetric_E_birth!)

    ###########
    # Callbacks
    ###########

    n0 = sim_eff.n0

    # Decide each complete event from its source and full biological activity.
    event_flags = independent ? (
            :S_EG0_symmetric_birth_ode, :S_EG0_to_R_transition_ode, :S_EG0_death_ode,
            :R_EG0_symmetric_birth_ode, :R_EG0_to_S_transition_ode, :R_EG0_to_E_transition_ode,
            :R_EG0_death_ode, :E_EG0_symmetric_birth_ode, :E_EG0_death_ode,
            :S_EG1_symmetric_birth_ode, :S_EG1_to_R_transition_ode, :S_EG1_death_ode,
            :R_EG1_symmetric_birth_ode, :R_EG1_to_S_transition_ode, :R_EG1_to_E_transition_ode,
            :R_EG1_death_ode, :E_EG1_symmetric_birth_ode, :E_EG1_death_ode
    ) : (
            :S_EG0_symmetric_birth_ode, :S_EG0_asymmetric_R_birth_ode, :S_EG0_death_ode,
            :R_EG0_symmetric_birth_ode, :R_EG0_asymmetric_S_birth_ode, :R_EG0_asymmetric_E_birth_ode,
            :R_EG0_death_ode, :E_EG0_symmetric_birth_ode, :E_EG0_death_ode,
            :S_EG1_symmetric_birth_ode, :S_EG1_asymmetric_R_birth_ode, :S_EG1_death_ode,
            :R_EG1_symmetric_birth_ode, :R_EG1_asymmetric_S_birth_ode, :R_EG1_asymmetric_E_birth_ode,
            :R_EG1_death_ode, :E_EG1_symmetric_birth_ode, :E_EG1_death_ode
    )

    function classify_event_modes!(u, p, t)
        previous = Tuple(getproperty(p, name) for name in event_flags)
        p.S_EG0_symmetric_birth_ode = u[RESPOP_INVIVO_NS_EG0_INDEX] >= sim_eff.Nswitch &&
            S_EG0_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.S_EG0_to_R_transition_ode = u[RESPOP_INVIVO_NS_EG0_INDEX] >= sim_eff.Nswitch &&
                S_EG0_to_R_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.S_EG0_asymmetric_R_birth_ode = u[RESPOP_INVIVO_NS_EG0_INDEX] >= sim_eff.Nswitch &&
                S_EG0_asymmetric_R_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.S_EG0_death_ode = u[RESPOP_INVIVO_NS_EG0_INDEX] >= sim_eff.Nswitch &&
            S_EG0_death_activity(u, p, t) >= sim_eff.Aswitch
        p.R_EG0_symmetric_birth_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
            R_EG0_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.R_EG0_to_S_transition_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
                R_EG0_to_S_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_EG0_asymmetric_S_birth_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
                R_EG0_asymmetric_S_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        if independent
            p.R_EG0_to_E_transition_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
                R_EG0_to_E_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_EG0_asymmetric_E_birth_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
                R_EG0_asymmetric_E_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.R_EG0_death_ode = u[RESPOP_INVIVO_NR_EG0_INDEX] >= sim_eff.Nswitch &&
            R_EG0_death_activity(u, p, t) >= sim_eff.Aswitch
        p.E_EG0_symmetric_birth_ode = u[RESPOP_INVIVO_NE_EG0_INDEX] >= sim_eff.Nswitch &&
            E_EG0_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        p.E_EG0_death_ode = u[RESPOP_INVIVO_NE_EG0_INDEX] >= sim_eff.Nswitch &&
            E_EG0_death_activity(u, p, t) >= sim_eff.Aswitch
        p.S_EG1_symmetric_birth_ode = u[RESPOP_INVIVO_NS_EG1_INDEX] >= sim_eff.Nswitch &&
            S_EG1_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.S_EG1_to_R_transition_ode = u[RESPOP_INVIVO_NS_EG1_INDEX] >= sim_eff.Nswitch &&
                S_EG1_to_R_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.S_EG1_asymmetric_R_birth_ode = u[RESPOP_INVIVO_NS_EG1_INDEX] >= sim_eff.Nswitch &&
                S_EG1_asymmetric_R_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.S_EG1_death_ode = u[RESPOP_INVIVO_NS_EG1_INDEX] >= sim_eff.Nswitch &&
            S_EG1_death_activity(u, p, t) >= sim_eff.Aswitch
        p.R_EG1_symmetric_birth_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
            R_EG1_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        if independent
            p.R_EG1_to_S_transition_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
                R_EG1_to_S_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_EG1_asymmetric_S_birth_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
                R_EG1_asymmetric_S_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        if independent
            p.R_EG1_to_E_transition_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
                R_EG1_to_E_transition_activity(u, p, t) >= sim_eff.Aswitch
        else
            p.R_EG1_asymmetric_E_birth_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
                R_EG1_asymmetric_E_birth_activity(u, p, t) >= sim_eff.Aswitch
        end
        p.R_EG1_death_ode = u[RESPOP_INVIVO_NR_EG1_INDEX] >= sim_eff.Nswitch &&
            R_EG1_death_activity(u, p, t) >= sim_eff.Aswitch
        p.E_EG1_symmetric_birth_ode = u[RESPOP_INVIVO_NE_EG1_INDEX] >= sim_eff.Nswitch &&
            E_EG1_symmetric_birth_activity(u, p, t) >= sim_eff.Aswitch
        p.E_EG1_death_ode = u[RESPOP_INVIVO_NE_EG1_INDEX] >= sim_eff.Nswitch &&
            E_EG1_death_activity(u, p, t) >= sim_eff.Aswitch
        return previous != Tuple(getproperty(p, name) for name in event_flags)
    end

    function refresh_event_modes!(u, p, t)
        changed = classify_event_modes!(u, p, t)
        # Evaluate all eligibility before changing any count; net drift is not sufficient.
        NS_EG0_has_ode = (p.S_EG0_symmetric_birth_ode == 1.0 && S_EG0_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.S_EG0_death_ode == 1.0 && S_EG0_death_activity(u, p, t) > 0.0)
        NR_EG0_has_ode = (p.R_EG0_symmetric_birth_ode == 1.0 && R_EG0_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.R_EG0_death_ode == 1.0 && R_EG0_death_activity(u, p, t) > 0.0)
        NE_EG0_has_ode = (p.E_EG0_symmetric_birth_ode == 1.0 && E_EG0_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.E_EG0_death_ode == 1.0 && E_EG0_death_activity(u, p, t) > 0.0)
        if independent
            S_EG0_to_R_transition_active = p.S_EG0_to_R_transition_ode == 1.0 && S_EG0_to_R_transition_activity(u, p, t) > 0.0
            R_EG0_to_S_transition_active = p.R_EG0_to_S_transition_ode == 1.0 && R_EG0_to_S_transition_activity(u, p, t) > 0.0
            R_EG0_to_E_transition_active = p.R_EG0_to_E_transition_ode == 1.0 && R_EG0_to_E_transition_activity(u, p, t) > 0.0
            NS_EG0_has_ode |= S_EG0_to_R_transition_active || R_EG0_to_S_transition_active
            NR_EG0_has_ode |= S_EG0_to_R_transition_active || R_EG0_to_S_transition_active || R_EG0_to_E_transition_active
            NE_EG0_has_ode |= R_EG0_to_E_transition_active
        else
            NS_EG0_has_ode |= p.R_EG0_asymmetric_S_birth_ode == 1.0 && R_EG0_asymmetric_S_birth_activity(u, p, t) > 0.0
            NR_EG0_has_ode |= p.S_EG0_asymmetric_R_birth_ode == 1.0 && S_EG0_asymmetric_R_birth_activity(u, p, t) > 0.0
            NE_EG0_has_ode |= p.R_EG0_asymmetric_E_birth_ode == 1.0 && R_EG0_asymmetric_E_birth_activity(u, p, t) > 0.0
        end
        NS_EG1_has_ode = (p.S_EG1_symmetric_birth_ode == 1.0 && S_EG1_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.S_EG1_death_ode == 1.0 && S_EG1_death_activity(u, p, t) > 0.0)
        NR_EG1_has_ode = (p.R_EG1_symmetric_birth_ode == 1.0 && R_EG1_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.R_EG1_death_ode == 1.0 && R_EG1_death_activity(u, p, t) > 0.0)
        NE_EG1_has_ode = (p.E_EG1_symmetric_birth_ode == 1.0 && E_EG1_symmetric_birth_activity(u, p, t) > 0.0) ||
            (p.E_EG1_death_ode == 1.0 && E_EG1_death_activity(u, p, t) > 0.0)
        if independent
            S_EG1_to_R_transition_active = p.S_EG1_to_R_transition_ode == 1.0 && S_EG1_to_R_transition_activity(u, p, t) > 0.0
            R_EG1_to_S_transition_active = p.R_EG1_to_S_transition_ode == 1.0 && R_EG1_to_S_transition_activity(u, p, t) > 0.0
            R_EG1_to_E_transition_active = p.R_EG1_to_E_transition_ode == 1.0 && R_EG1_to_E_transition_activity(u, p, t) > 0.0
            NS_EG1_has_ode |= S_EG1_to_R_transition_active || R_EG1_to_S_transition_active
            NR_EG1_has_ode |= S_EG1_to_R_transition_active || R_EG1_to_S_transition_active || R_EG1_to_E_transition_active
            NE_EG1_has_ode |= R_EG1_to_E_transition_active
        else
            NS_EG1_has_ode |= p.R_EG1_asymmetric_S_birth_ode == 1.0 && R_EG1_asymmetric_S_birth_activity(u, p, t) > 0.0
            NR_EG1_has_ode |= p.S_EG1_asymmetric_R_birth_ode == 1.0 && S_EG1_asymmetric_R_birth_activity(u, p, t) > 0.0
            NE_EG1_has_ode |= p.R_EG1_asymmetric_E_birth_ode == 1.0 && R_EG1_asymmetric_E_birth_activity(u, p, t) > 0.0
        end
        projected = false
        if !NS_EG0_has_ode && 0.0 <= u[RESPOP_INVIVO_NS_EG0_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NS_EG0_INDEX])
            lower = floor(u[RESPOP_INVIVO_NS_EG0_INDEX])
            u[RESPOP_INVIVO_NS_EG0_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NS_EG0_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NR_EG0_has_ode && 0.0 <= u[RESPOP_INVIVO_NR_EG0_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NR_EG0_INDEX])
            lower = floor(u[RESPOP_INVIVO_NR_EG0_INDEX])
            u[RESPOP_INVIVO_NR_EG0_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NR_EG0_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NE_EG0_has_ode && 0.0 <= u[RESPOP_INVIVO_NE_EG0_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NE_EG0_INDEX])
            lower = floor(u[RESPOP_INVIVO_NE_EG0_INDEX])
            u[RESPOP_INVIVO_NE_EG0_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NE_EG0_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NS_EG1_has_ode && 0.0 <= u[RESPOP_INVIVO_NS_EG1_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NS_EG1_INDEX])
            lower = floor(u[RESPOP_INVIVO_NS_EG1_INDEX])
            u[RESPOP_INVIVO_NS_EG1_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NS_EG1_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NR_EG1_has_ode && 0.0 <= u[RESPOP_INVIVO_NR_EG1_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NR_EG1_INDEX])
            lower = floor(u[RESPOP_INVIVO_NR_EG1_INDEX])
            u[RESPOP_INVIVO_NR_EG1_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NR_EG1_INDEX] - lower ? 1.0 : 0.0)
            projected = true
        end
        if !NE_EG1_has_ode && 0.0 <= u[RESPOP_INVIVO_NE_EG1_INDEX] < sim_eff.Nswitch && !isinteger(u[RESPOP_INVIVO_NE_EG1_INDEX])
            lower = floor(u[RESPOP_INVIVO_NE_EG1_INDEX])
            u[RESPOP_INVIVO_NE_EG1_INDEX] = lower + (rand(rng) < u[RESPOP_INVIVO_NE_EG1_INDEX] - lower ? 1.0 : 0.0)
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

    # Sample first, then thin engraftment; both use the realized current population.
    function attempt_passage!(integrator, n0; require_enough::Bool)
        counts = round.(Int, integrator.u[2:7])
        require_enough && sum(counts) < n0 && return false
        draw = multivariate_hypergeometric_draw(counts, n0)
        eg0, eg1, stats = engraftment_selection(draw[[1, 3, 5]], draw[[2, 4, 6]], params.pEG, params.sEG)
        integrator.u[2:7] .= (eg0[1], eg1[1], eg0[2], eg1[2], eg0[3], eg1[3])
        integrator.u[RESPOP_INVIVO_PASS_INDEX] += 1
        on_engraft(round(Int, integrator.u[RESPOP_INVIVO_PASS_INDEX]), stats)
        refresh_integrator!(integrator; modified = true)
        respop_invivo_total_population(integrator.u) < 1.0 && terminate!(integrator)
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
    gam_above_max(u, t, integrator) = (u[RESPOP_INVIVO_GAM_INDEX] - 1.0)
    function clamp_gam_max!(integrator)
        integrator.u[RESPOP_INVIVO_GAM_INDEX] = 1.0
        integrator.p.kp = 0.0
        refresh_integrator!(integrator; modified = true)
    end

    gam_below_min(u, t, integrator) = (u[RESPOP_INVIVO_GAM_INDEX])
    function clamp_gam_min!(integrator)
        integrator.u[RESPOP_INVIVO_GAM_INDEX] = 0.0
        integrator.p.kp = 0.0
        refresh_integrator!(integrator; modified = true)
    end

    gam_max_cb = ContinuousCallback(gam_above_max, clamp_gam_max!,
                                    save_positions = (false, true))
    gam_min_cb = ContinuousCallback(gam_below_min, clamp_gam_min!,
                                    save_positions = (false, true))

    # Nmax, Extinction & Passage
    nmax_reached(u, t, integrator) = (respop_invivo_total_population(integrator.u) >= sim_eff.Nmax)

    function handle_nmax!(integrator)
        if integrator.u[RESPOP_INVIVO_PASS_INDEX] >= n_pass_eff
            terminate!(integrator)
        else
            attempt_passage!(integrator, n0; require_enough = false)
        end
    end

    nmax_cb = DiscreteCallback(nmax_reached, handle_nmax!,
                               save_positions = (false, true))

    extinction(u, t, integrator) = (respop_invivo_total_population(integrator.u) < 1.0)
    extinction_cb = DiscreteCallback(extinction, terminate!,
                                     save_positions = (false, false))

    function apply_scheduled_passage!(integrator)
        curr_pass = round(Int, integrator.u[RESPOP_INVIVO_PASS_INDEX])
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

    # Population clamps preserve the same event refresh contract in both EG strata.
    cb_list = Any[event_mode_cb, nmax_cb, extinction_cb,
        (build_clamp_callback(idx) for idx in 2:7)...]

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
                S_EG0_symmetric_birth_activity, (independent ? S_EG0_to_R_transition_activity : S_EG0_asymmetric_R_birth_activity), S_EG0_death_activity,
                R_EG0_symmetric_birth_activity, (independent ? R_EG0_to_S_transition_activity : R_EG0_asymmetric_S_birth_activity), (independent ? R_EG0_to_E_transition_activity : R_EG0_asymmetric_E_birth_activity),
                R_EG0_death_activity, E_EG0_symmetric_birth_activity, E_EG0_death_activity,
                S_EG1_symmetric_birth_activity, (independent ? S_EG1_to_R_transition_activity : S_EG1_asymmetric_R_birth_activity), S_EG1_death_activity,
                R_EG1_symmetric_birth_activity, (independent ? R_EG1_to_S_transition_activity : R_EG1_asymmetric_S_birth_activity), (independent ? R_EG1_to_E_transition_activity : R_EG1_asymmetric_E_birth_activity),
                R_EG1_death_activity, E_EG1_symmetric_birth_activity, E_EG1_death_activity
            ),
            jump_rates = (
                S_EG0_symmetric_birth_rate, (independent ? S_EG0_to_R_transition_rate : S_EG0_asymmetric_R_birth_rate), S_EG0_death_rate,
                R_EG0_symmetric_birth_rate, (independent ? R_EG0_to_S_transition_rate : R_EG0_asymmetric_S_birth_rate), (independent ? R_EG0_to_E_transition_rate : R_EG0_asymmetric_E_birth_rate),
                R_EG0_death_rate, E_EG0_symmetric_birth_rate, E_EG0_death_rate,
                S_EG1_symmetric_birth_rate, (independent ? S_EG1_to_R_transition_rate : S_EG1_asymmetric_R_birth_rate), S_EG1_death_rate,
                R_EG1_symmetric_birth_rate, (independent ? R_EG1_to_S_transition_rate : R_EG1_asymmetric_S_birth_rate), (independent ? R_EG1_to_E_transition_rate : R_EG1_asymmetric_E_birth_rate),
                R_EG1_death_rate, E_EG1_symmetric_birth_rate, E_EG1_death_rate
            ),
            flags = event_flags,
            sources = (
                RESPOP_INVIVO_NS_EG0_INDEX, RESPOP_INVIVO_NS_EG0_INDEX, RESPOP_INVIVO_NS_EG0_INDEX,
                RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NR_EG0_INDEX,
                RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NE_EG0_INDEX, RESPOP_INVIVO_NE_EG0_INDEX,
                RESPOP_INVIVO_NS_EG1_INDEX, RESPOP_INVIVO_NS_EG1_INDEX, RESPOP_INVIVO_NS_EG1_INDEX,
                RESPOP_INVIVO_NR_EG1_INDEX, RESPOP_INVIVO_NR_EG1_INDEX, RESPOP_INVIVO_NR_EG1_INDEX,
                RESPOP_INVIVO_NR_EG1_INDEX, RESPOP_INVIVO_NE_EG1_INDEX, RESPOP_INVIVO_NE_EG1_INDEX
            ),
            classify! = classify_event_modes!, refresh! = refresh_event_modes!)
    end

    cbs = CallbackSet(cb_list...)

    sjm_prob = JumpProblem(prob, Direct(),
        S_EG0_symmetric_birth_jump,
        (independent ? S_EG0_to_R_transition_jump : S_EG0_asymmetric_R_birth_jump),
        S_EG0_death_jump,
        R_EG0_symmetric_birth_jump,
        (independent ? R_EG0_to_S_transition_jump : R_EG0_asymmetric_S_birth_jump),
        (independent ? R_EG0_to_E_transition_jump : R_EG0_asymmetric_E_birth_jump),
        R_EG0_death_jump,
        E_EG0_symmetric_birth_jump,
        E_EG0_death_jump,
        S_EG1_symmetric_birth_jump,
        (independent ? S_EG1_to_R_transition_jump : S_EG1_asymmetric_R_birth_jump),
        S_EG1_death_jump,
        R_EG1_symmetric_birth_jump,
        (independent ? R_EG1_to_S_transition_jump : R_EG1_asymmetric_S_birth_jump),
        (independent ? R_EG1_to_E_transition_jump : R_EG1_asymmetric_E_birth_jump),
        R_EG1_death_jump,
        E_EG1_symmetric_birth_jump,
        E_EG1_death_jump; rng = rng)

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
