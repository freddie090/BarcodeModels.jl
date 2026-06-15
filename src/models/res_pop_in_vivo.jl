"""Resistance population in vivo hybrid model with engraftment phenotype parameters."""
struct ResPopInVivo <: HybridModel
    params::ResPopInVivoParams
    function ResPopInVivo(params::ResPopInVivoParams)
        validate_model_params(params)
        new(params)
    end
end

"""Construct a `ResPopInVivo` from keyword arguments forwarded to `ResPopInVivoParams`."""
ResPopInVivo(; kwargs...) = ResPopInVivo(ResPopInVivoParams(; kwargs...))

function _respop_params_from_invivo(params::ResPopInVivoParams; al = params.al, rho = params.rho, drug_effect = params.drug_effect)
    return ResPopParams(
        b = params.b,
        d = params.d,
        rho = rho,
        mu = params.mu,
        sig = params.sig,
        del = params.del,
        al = al,
        Dc = params.Dc,
        k = params.k,
        psi = params.psi,
        drug_effect = drug_effect
    )
end

"""Run ResPop dynamics on six S/R/E x EG compartments with static EG labels."""
function run_model_core_hybrid(model::ResPopInVivo, state::ResPopInVivoState, sim::SimParams; treat::Bool = sim.treat)
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
        N_trans_switch = sim.N_trans_switch
    )

    params = model.params
    de = params.drug_effect

    if de == :b
        bR_check = params.b * (1 - params.del)
        psi_scale = 1 - params.psi
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")
        if psi_scale > 0.0
            params.Dc <= (bR_check / psi_scale) || error("When drug_effect == :b, Dc*(1-psi) must be <= b*(1-del) to keep resistant birth non-negative. Use drug_effect == :c if stronger drug effects are intended.")
        end
    end

    base_params = _respop_params_from_invivo(params)
    p, bS, dS, bR, dR, bE, dE = build_component_params(base_params)

    u0 = to_componentarray(state)
    tspan = (sim_eff.t0, sim_eff.tmax)

    rate_modifier = select_drug_effect(de,
        apply_drug_effect_death,
        apply_drug_effect_birth,
        apply_drug_effect_combined
    )

    phenotype_total(u, idx0, idx1) = u[idx0] + u[idx1]
    logistic_scale(u) = max(logistic_factor(respop_invivo_total_population(u), sim_eff.Cc), 0.0)

    function ode_fxn!(du, u, p, t)
        @unpack kp, psi, mu_o, sig_o, al_o = p
        @unpack gam, nS_EG0, nS_EG1, nR_EG0, nR_EG1, nE_EG0, nE_EG1 = u

        N = nS_EG0 + nS_EG1 + nR_EG0 + nR_EG1 + nE_EG0 + nE_EG1

        bS_o_mod, dS_o_mod = rate_modifier(bS, dS, p.bS_o, p.dS_o, p.Dc, gam, N, sim_eff.Cc, 0.0)
        bR_o_mod, dR_o_mod = rate_modifier(bR, dR, p.bR_o, p.dR_o, p.Dc, gam, N, sim_eff.Cc, psi)
        bE_o_mod, dE_o_mod = rate_modifier(bE, dE, p.bE_o, p.dE_o, p.Dc, gam, N, sim_eff.Cc, psi)

        du.gam = kp
        du.nS_EG0 = (bS_o_mod - dS_o_mod) * nS_EG0 - mu_o * nS_EG0 * bS_o_mod + sig_o * nR_EG0 * bR_o_mod
        du.nS_EG1 = (bS_o_mod - dS_o_mod) * nS_EG1 - mu_o * nS_EG1 * bS_o_mod + sig_o * nR_EG1 * bR_o_mod
        du.nR_EG0 = (bR_o_mod - dR_o_mod) * nR_EG0 + mu_o * nS_EG0 * bS_o_mod - sig_o * nR_EG0 * bR_o_mod - al_o * gam * nR_EG0 * bR_o_mod
        du.nR_EG1 = (bR_o_mod - dR_o_mod) * nR_EG1 + mu_o * nS_EG1 * bS_o_mod - sig_o * nR_EG1 * bR_o_mod - al_o * gam * nR_EG1 * bR_o_mod
        du.nE_EG0 = (bE_o_mod - dE_o_mod) * nE_EG0 + al_o * gam * nR_EG0 * bR_o_mod
        du.nE_EG1 = (bE_o_mod - dE_o_mod) * nE_EG1 + al_o * gam * nR_EG1 * bR_o_mod
        du.Pass_num = 0.0
        nothing
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    function build_birth_death_rates(idx, b_ref, d_ref, b_sym::Symbol, d_sym::Symbol; psi_fn = p -> 0.0)
        function birth_d(u, p, t)
            return u[idx] * getproperty(p, b_sym) * logistic_scale(u)
        end

        function death_d(u, p, t)
            drug_interaction = safe_ratio(getproperty(p, d_sym), d_ref) * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            return u[idx] * (getproperty(p, d_sym) + drug_interaction) * logistic_scale(u)
        end

        function birth_b(u, p, t)
            drug_interaction = safe_ratio(getproperty(p, b_sym), b_ref) * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            return u[idx] * (getproperty(p, b_sym) - drug_interaction) * logistic_scale(u)
        end

        function death_b(u, p, t)
            return u[idx] * getproperty(p, d_sym) * logistic_scale(u)
        end

        function birth_c(u, p, t)
            drug_interaction = safe_ratio(getproperty(p, b_sym), b_ref) * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            return u[idx] * max(getproperty(p, b_sym) - drug_interaction, 0.0) * logistic_scale(u)
        end

        function death_c(u, p, t)
            drug_interaction = safe_ratio(getproperty(p, b_sym), b_ref) * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            b_drug = getproperty(p, b_sym) - drug_interaction
            death_rate = b_drug >= 0.0 ? getproperty(p, d_sym) : getproperty(p, d_sym) + abs(b_drug)
            return u[idx] * death_rate * logistic_scale(u)
        end

        return select_drug_effect(de, birth_d, birth_b, birth_c),
               select_drug_effect(de, death_d, death_b, death_c)
    end

    function build_switch_rate(idx, b_total_fn, rate_fn; gamma_factor = false, psi_fn = p -> 0.0)
        function b_effective(u, p)
            b_total = b_total_fn(p)
            drug_interaction = b_total * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            if de === :d
                return b_total
            elseif de === :b
                return b_total - drug_interaction
            else
                return max(b_total - drug_interaction, 0.0)
            end
        end

        function rate(u, p, t)
            out = u[idx] * rate_fn(p) * b_effective(u, p)
            if gamma_factor
                out *= u[RESPOP_INVIVO_GAM_INDEX]
            end
            return out * logistic_scale(u)
        end
        return rate
    end

    function build_switch_activity_proxy(idxs, b_total_fn, rate_value; gamma_factor = false, psi_fn = p -> 0.0)
        function proxy(integrator)
            u = integrator.u
            p = integrator.p
            b_total = b_total_fn(p)
            drug_interaction = b_total * p.Dc * u[RESPOP_INVIVO_GAM_INDEX] * (1 - psi_fn(p))
            b_drug = b_total - drug_interaction
            b_effective = if de === :d
                b_total
            elseif de === :b
                b_drug
            else
                max(b_drug, 0.0)
            end

            activity = sum(u[idx] for idx in idxs) * rate_value * b_effective
            if gamma_factor
                activity *= u[RESPOP_INVIVO_GAM_INDEX]
            end
            return activity * logistic_scale(u)
        end
        return proxy
    end

    function birth_affect!(idx)
        return integrator -> begin
            integrator.u[idx] += 1
            nothing
        end
    end

    function death_affect!(idx)
        return integrator -> begin
            integrator.u[idx] -= 1
            nothing
        end
    end

    function switch_affect!(from_idx, to_idx)
        return integrator -> begin
            integrator.u[from_idx] -= 1
            integrator.u[to_idx] += 1
            nothing
        end
    end

    S0b_rate, S0d_rate = build_birth_death_rates(RESPOP_INVIVO_NS_EG0_INDEX, bS, dS, :bS_j, :dS_j)
    S1b_rate, S1d_rate = build_birth_death_rates(RESPOP_INVIVO_NS_EG1_INDEX, bS, dS, :bS_j, :dS_j)
    R0b_rate, R0d_rate = build_birth_death_rates(RESPOP_INVIVO_NR_EG0_INDEX, bR, dR, :bR_j, :dR_j; psi_fn = p -> p.psi)
    R1b_rate, R1d_rate = build_birth_death_rates(RESPOP_INVIVO_NR_EG1_INDEX, bR, dR, :bR_j, :dR_j; psi_fn = p -> p.psi)
    E0b_rate, E0d_rate = build_birth_death_rates(RESPOP_INVIVO_NE_EG0_INDEX, bE, dE, :bE_j, :dE_j; psi_fn = p -> p.psi)
    E1b_rate, E1d_rate = build_birth_death_rates(RESPOP_INVIVO_NE_EG1_INDEX, bE, dE, :bE_j, :dE_j; psi_fn = p -> p.psi)

    jumps = Any[
        VariableRateJump(S0b_rate, birth_affect!(RESPOP_INVIVO_NS_EG0_INDEX)),
        VariableRateJump(S0d_rate, death_affect!(RESPOP_INVIVO_NS_EG0_INDEX)),
        VariableRateJump(S1b_rate, birth_affect!(RESPOP_INVIVO_NS_EG1_INDEX)),
        VariableRateJump(S1d_rate, death_affect!(RESPOP_INVIVO_NS_EG1_INDEX)),
        VariableRateJump(R0b_rate, birth_affect!(RESPOP_INVIVO_NR_EG0_INDEX)),
        VariableRateJump(R0d_rate, death_affect!(RESPOP_INVIVO_NR_EG0_INDEX)),
        VariableRateJump(R1b_rate, birth_affect!(RESPOP_INVIVO_NR_EG1_INDEX)),
        VariableRateJump(R1d_rate, death_affect!(RESPOP_INVIVO_NR_EG1_INDEX)),
        VariableRateJump(E0b_rate, birth_affect!(RESPOP_INVIVO_NE_EG0_INDEX)),
        VariableRateJump(E0d_rate, death_affect!(RESPOP_INVIVO_NE_EG0_INDEX)),
        VariableRateJump(E1b_rate, birth_affect!(RESPOP_INVIVO_NE_EG1_INDEX)),
        VariableRateJump(E1d_rate, death_affect!(RESPOP_INVIVO_NE_EG1_INDEX)),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NS_EG0_INDEX, p -> (p.bS_o + p.bS_j), p -> p.mu_j),
            switch_affect!(RESPOP_INVIVO_NS_EG0_INDEX, RESPOP_INVIVO_NR_EG0_INDEX)
        ),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NS_EG1_INDEX, p -> (p.bS_o + p.bS_j), p -> p.mu_j),
            switch_affect!(RESPOP_INVIVO_NS_EG1_INDEX, RESPOP_INVIVO_NR_EG1_INDEX)
        ),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NR_EG0_INDEX, p -> (p.bR_o + p.bR_j), p -> p.sig_j; psi_fn = p -> p.psi),
            switch_affect!(RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NS_EG0_INDEX)
        ),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NR_EG1_INDEX, p -> (p.bR_o + p.bR_j), p -> p.sig_j; psi_fn = p -> p.psi),
            switch_affect!(RESPOP_INVIVO_NR_EG1_INDEX, RESPOP_INVIVO_NS_EG1_INDEX)
        ),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NR_EG0_INDEX, p -> (p.bR_o + p.bR_j), p -> p.al_j; gamma_factor = true, psi_fn = p -> p.psi),
            switch_affect!(RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NE_EG0_INDEX)
        ),
        VariableRateJump(
            build_switch_rate(RESPOP_INVIVO_NR_EG1_INDEX, p -> (p.bR_o + p.bR_j), p -> p.al_j; gamma_factor = true, psi_fn = p -> p.psi),
            switch_affect!(RESPOP_INVIVO_NR_EG1_INDEX, RESPOP_INVIVO_NE_EG1_INDEX)
        )
    ]

    function round_idxs!(integrator, idxs)
        for idx in idxs
            integrator.u[idx] = round(integrator.u[idx])
        end
        nothing
    end

    function build_bd_toggle_callbacks(total_fn, idxs, b_o_sym, b_j_sym, d_o_sym, d_j_sym, Nswitch)
        cond_to_ode(u, t, integrator) = total_fn(integrator.u) >= Nswitch
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

        cond_to_jump(u, t, integrator) = total_fn(integrator.u) < Nswitch
        function switch_to_jump!(integrator)
            round_idxs!(integrator, idxs)
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

        return DiscreteCallback(cond_to_ode, switch_to_ode!; save_positions = (false, true)),
               DiscreteCallback(cond_to_jump, switch_to_jump!; save_positions = (false, true))
    end

    function build_rate_toggle_callbacks(activity_fn, rate_o_sym, rate_j_sym, N_trans_switch)
        cond_to_ode(u, t, integrator) = activity_fn(integrator) >= N_trans_switch
        function switch_to_ode!(integrator)
            if getproperty(integrator.p, rate_o_sym) == 0.0
                setproperty!(integrator.p, rate_o_sym, deepcopy(getproperty(integrator.p, rate_j_sym)))
                setproperty!(integrator.p, rate_j_sym, 0.0)
            end
            nothing
        end

        cond_to_jump(u, t, integrator) = activity_fn(integrator) < N_trans_switch
        function switch_to_jump!(integrator)
            if getproperty(integrator.p, rate_j_sym) == 0.0
                setproperty!(integrator.p, rate_j_sym, deepcopy(getproperty(integrator.p, rate_o_sym)))
                setproperty!(integrator.p, rate_o_sym, 0.0)
            end
            nothing
        end

        return DiscreteCallback(cond_to_ode, switch_to_ode!; save_positions = (false, true)),
               DiscreteCallback(cond_to_jump, switch_to_jump!; save_positions = (false, true))
    end

    s_idxs = (RESPOP_INVIVO_NS_EG0_INDEX, RESPOP_INVIVO_NS_EG1_INDEX)
    r_idxs = (RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NR_EG1_INDEX)
    e_idxs = (RESPOP_INVIVO_NE_EG0_INDEX, RESPOP_INVIVO_NE_EG1_INDEX)

    s_total(u) = phenotype_total(u, RESPOP_INVIVO_NS_EG0_INDEX, RESPOP_INVIVO_NS_EG1_INDEX)
    r_total(u) = phenotype_total(u, RESPOP_INVIVO_NR_EG0_INDEX, RESPOP_INVIVO_NR_EG1_INDEX)
    e_total(u) = phenotype_total(u, RESPOP_INVIVO_NE_EG0_INDEX, RESPOP_INVIVO_NE_EG1_INDEX)

    S_cb_switch1, S_cb_switch2 = build_bd_toggle_callbacks(s_total, s_idxs, :bS_o, :bS_j, :dS_o, :dS_j, sim_eff.Nswitch)
    R_cb_switch1, R_cb_switch2 = build_bd_toggle_callbacks(r_total, r_idxs, :bR_o, :bR_j, :dR_o, :dR_j, sim_eff.Nswitch)
    E_cb_switch1, E_cb_switch2 = build_bd_toggle_callbacks(e_total, e_idxs, :bE_o, :bE_j, :dE_o, :dE_j, sim_eff.Nswitch)

    StoR_activity_proxy = build_switch_activity_proxy(s_idxs, p -> (p.bS_o + p.bS_j), params.mu)
    RtoS_activity_proxy = build_switch_activity_proxy(r_idxs, p -> (p.bR_o + p.bR_j), params.sig; psi_fn = p -> p.psi)
    RtoE_activity_proxy = build_switch_activity_proxy(r_idxs, p -> (p.bR_o + p.bR_j), params.al; gamma_factor = true, psi_fn = p -> p.psi)

    StoR_cb_switch1, StoR_cb_switch2 = build_rate_toggle_callbacks(StoR_activity_proxy, :mu_o, :mu_j, sim_eff.N_trans_switch)
    RtoS_cb_switch1, RtoS_cb_switch2 = build_rate_toggle_callbacks(RtoS_activity_proxy, :sig_o, :sig_j, sim_eff.N_trans_switch)
    RtoE_cb_switch1, RtoE_cb_switch2 = build_rate_toggle_callbacks(RtoE_activity_proxy, :al_o, :al_j, sim_eff.N_trans_switch)

    function set_treatment_on!(integrator)
        integrator.p.kp = params.k
    end

    function set_treatment_off!(integrator)
        integrator.p.kp = -params.k
    end

    filter_times_in_span(times) = sort(unique(filter(t -> (t >= sim_eff.t0 && t <= sim_eff.tmax), times)))
    treat_on_times = filter_times_in_span(sim_eff.treat_ons)
    treat_off_times = filter_times_in_span(sim_eff.treat_offs)

    treat_on_cb = isempty(treat_on_times) ? nothing :
                  PresetTimeCallback(treat_on_times, set_treatment_on!; save_positions = (false, true))
    treat_off_cb = isempty(treat_off_times) ? nothing :
                   PresetTimeCallback(treat_off_times, set_treatment_off!; save_positions = (false, true))

    function clamp_component_callback(idx)
        condition(u, t, integrator) = u[idx] < 0.0
        function affect!(integrator)
            integrator.u[idx] = 0.0
        end
        return DiscreteCallback(condition, affect!; save_positions = (false, true))
    end

    nmax_reached(u, t, integrator) = respop_invivo_total_population(integrator.u) >= sim_eff.Nmax
    extinction(u, t, integrator) = respop_invivo_total_population(integrator.u) < 1.0

    gam_above_max(u, t, integrator) = u[RESPOP_INVIVO_GAM_INDEX] - 1.0
    function clamp_gam_max!(integrator)
        integrator.u[RESPOP_INVIVO_GAM_INDEX] = 1.0
        integrator.p.kp = 0.0
    end

    gam_below_min(u, t, integrator) = u[RESPOP_INVIVO_GAM_INDEX]
    function clamp_gam_min!(integrator)
        integrator.u[RESPOP_INVIVO_GAM_INDEX] = 0.0
        integrator.p.kp = 0.0
    end

    cb_list = Any[
        S_cb_switch1, S_cb_switch2,
        R_cb_switch1, R_cb_switch2,
        E_cb_switch1, E_cb_switch2,
        StoR_cb_switch1, StoR_cb_switch2,
        RtoS_cb_switch1, RtoS_cb_switch2,
        RtoE_cb_switch1, RtoE_cb_switch2,
        DiscreteCallback(nmax_reached, terminate!; save_positions = (false, true)),
        DiscreteCallback(extinction, terminate!; save_positions = (false, false)),
        clamp_component_callback(RESPOP_INVIVO_NS_EG0_INDEX),
        clamp_component_callback(RESPOP_INVIVO_NS_EG1_INDEX),
        clamp_component_callback(RESPOP_INVIVO_NR_EG0_INDEX),
        clamp_component_callback(RESPOP_INVIVO_NR_EG1_INDEX),
        clamp_component_callback(RESPOP_INVIVO_NE_EG0_INDEX),
        clamp_component_callback(RESPOP_INVIVO_NE_EG1_INDEX)
    ]

    if sim_eff.treat
        if treat_on_cb !== nothing
            push!(cb_list, treat_on_cb)
        end
        if treat_off_cb !== nothing
            push!(cb_list, treat_off_cb)
        end
        push!(cb_list,
              ContinuousCallback(gam_above_max, clamp_gam_max!; save_positions = (false, true)),
              ContinuousCallback(gam_below_min, clamp_gam_min!; save_positions = (false, true)))
    end

    base_tstops = collect(sim_eff.t0:sim_eff.save_at:sim_eff.tmax)
    event_times = sim_eff.treat ? vcat(treat_on_times, treat_off_times) : Float64[]
    tstops = isempty(event_times) ? base_tstops : sort(unique(vcat(base_tstops, event_times)))

    sjm_prob = JumpProblem(prob, Direct(), jumps...)
    sol = @maybe_suppress_solver solve(sjm_prob, Tsit5(), callback = CallbackSet(cb_list...), tstops = tstops)
    return sol
end

