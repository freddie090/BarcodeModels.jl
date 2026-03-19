"""Resistance + damage (res_dmg) hybrid model with an explicit damaged compartment."""
struct ResDmg <: HybridModel
    params::ResDmgParams
    function ResDmg(params::ResDmgParams)
        validate_model_params(params)
        new(params)
    end
end

ResDmg(; kwargs...) = ResDmg(ResDmgParams(; kwargs...))

"""Build component-rate parameters and return them with baseline rates."""
function build_component_params(params::ResDmgParams, ome::Float64, zet::Float64)
    bS = params.b
    dS = params.d

    bD = 0.0
    dD = params.d + params.Dc

    bR = params.b * (1 - params.del)
    dR = params.d * (1 - params.del)

    bE = params.b
    dE = params.d

    p = ComponentArray(
        bS_o = 0.0, dS_o = 0.0, bS_j = bS, dS_j = dS,
        bD_o = 0.0, dD_o = 0.0, bD_j = bD, dD_j = dD,
        bR_o = 0.0, dR_o = 0.0, bR_j = bR, dR_j = dR,
        bE_o = 0.0, dE_o = 0.0, bE_j = bE, dE_j = dE,
        Dc = params.Dc, kp = 0.0, psi = params.psi,
        mu_o = 0.0, sig_o = 0.0, al_o = 0.0,
        ome_S_o = 0.0, ome_R_o = 0.0, zet_o = 0.0,
        mu_j = params.mu, sig_j = params.sig, al_j = params.al,
        ome_S_j = ome, ome_R_j = ome, zet_j = zet
    )

    return p, bS, dS, bD, dD, bR, dR, bE, dE
end

"""
Run the hybrid growth/kill model with a model, state, and simulation settings.

Arguments: `model::ResDmg`, `state::ResDmgState`, `sim::SimParams`.
Returns: `ODESolution`.
"""
function run_model_core_hybrid(model::ResDmg, state::ResDmgState, sim::SimParams; treat::Bool = sim.treat)
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
        bR = params.b * (1 - params.del)
        psi_scale = 1 - params.psi
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")
        if psi_scale > 0.0
            params.Dc <= (bR / psi_scale) || error("When drug_effect == :b, Dc*(1-psi) must be <= b*(1-del) to keep resistant birth non-negative. Use drug_effect == :c if stronger drug effects are intended.")
        end
    end

    p, bS, dS, bD, dD, bR, dR, bE, dE = build_component_params(params, params.ome, params.zet)

    u0 = to_componentarray(state)
    tspan = (sim_eff.t0, sim_eff.tmax)

    rate_modifier = select_drug_effect(de,
        apply_drug_effect_death,
        apply_drug_effect_birth,
        apply_drug_effect_combined
    )

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

    function build_birth_death_rates(de, idx, b_ref, d_ref, b_sym::Symbol, d_sym::Symbol; psi_fn = p -> 0.0)
        safe_ratio(num, den) = abs(den) > eps(Float64) ? (num / den) : 0.0

        function logistic_scale(u)
            return max(resdmg_logistic_factor(u, sim_eff.Cc), 0.0)
        end

        function birth_drug_adjustment(u, p)
            return (safe_ratio(getproperty(p, b_sym), b_ref) * p.Dc * u[RESDMG_GAM_INDEX] * (1 - psi_fn(p)))
        end

        function death_drug_adjustment(u, p)
            return (safe_ratio(getproperty(p, d_sym), d_ref) * p.Dc * u[RESDMG_GAM_INDEX] * (1 - psi_fn(p)))
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

    function build_switch_rate(de, idx, b_total_fn, rate_fn; gamma_factor = false, psi_fn = p -> 0.0, logistic_scale = true)
        function logistic_term(u)
            return max(resdmg_logistic_factor(u, sim_eff.Cc), 0.0)
        end

        function rate_d(u, p, t)
            base = u[idx] * rate_fn(p) * b_total_fn(p)
            if gamma_factor
                base *= u[RESDMG_GAM_INDEX]
            end
            return logistic_scale ? (base * logistic_term(u)) : base
        end

        function rate_b(u, p, t)
            b_total = b_total_fn(p)
            drug_interaction = (b_total * p.Dc * u[RESDMG_GAM_INDEX] * (1 - psi_fn(p)))
            b_drug = b_total - drug_interaction
            base = u[idx] * b_drug * rate_fn(p)
            if gamma_factor
                base *= u[RESDMG_GAM_INDEX]
            end
            return logistic_scale ? (base * logistic_term(u)) : base
        end

        function rate_c(u, p, t)
            b_total = b_total_fn(p)
            drug_interaction = (b_total * p.Dc * u[RESDMG_GAM_INDEX] * (1 - psi_fn(p)))
            b_drug = b_total - drug_interaction
            if b_drug >= 0.0
                base = u[idx] * b_drug * rate_fn(p)
                if gamma_factor
                    base *= u[RESDMG_GAM_INDEX]
                end
                return logistic_scale ? (base * logistic_term(u)) : base
            else
                return 0.0
            end
        end

        return select_drug_effect(de, rate_d, rate_b, rate_c)
    end

    function ode_fxn!(du, u, p, t)
        @unpack kp, psi, mu_o, sig_o, al_o, ome_S_o, ome_R_o, zet_o = p
        @unpack gam, nS, nD, nR, nE = u

        N = nS + nD + nR + nE

        bS_o_mod, dS_o_mod = rate_modifier(bS, dS, p.bS_o, p.dS_o, p.Dc, gam, N, sim_eff.Cc, 0.0)
        bR_o_mod, dR_o_mod = rate_modifier(bR, dR, p.bR_o, p.dR_o, p.Dc, gam, N, sim_eff.Cc, psi)
        bE_o_mod, dE_o_mod = rate_modifier(bE, dE, p.bE_o, p.dE_o, p.Dc, gam, N, sim_eff.Cc, 0.0)

        bD_o_mod = p.bD_o * max(1 - (N / sim_eff.Cc), 0.0)
        dD_o_mod = p.dD_o * max(1 - (N / sim_eff.Cc), 0.0)

        du.gam = kp
        du.nS = (bS_o_mod - dS_o_mod) * nS -
            ome_S_o * gam * nS -
            mu_o * nS * bS_o_mod +
            sig_o * nR * bR_o_mod +
            zet_o * nD

        du.nR = (bR_o_mod - dR_o_mod) * nR +
            mu_o * nS * bS_o_mod -
            sig_o * nR * bR_o_mod -
            al_o * nR * bR_o_mod -
            ome_R_o * gam * (1 - psi) * nR

        du.nD = (bD_o_mod - dD_o_mod) * nD +
            ome_S_o * gam * nS +
            ome_R_o * gam * (1 - psi) * nR -
            zet_o * nD

        du.nE = (bE_o_mod - dE_o_mod) * nE +
            al_o * nR * bR_o_mod
    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    function S_birth!(integrator)
        integrator.u[RESDMG_NS_INDEX] += 1
        nothing
    end

    function S_death!(integrator)
        integrator.u[RESDMG_NS_INDEX] -= 1
        nothing
    end

    Sb_rate, Sd_rate = build_birth_death_rates(de, RESDMG_NS_INDEX, bS, dS, :bS_j, :dS_j)
    Sb_jump = VariableRateJump(Sb_rate, S_birth!)
    Sd_jump = VariableRateJump(Sd_rate, S_death!)

    function D_birth!(integrator)
        integrator.u[RESDMG_ND_INDEX] += 1
        nothing
    end

    function D_death!(integrator)
        integrator.u[RESDMG_ND_INDEX] -= 1
        nothing
    end

    Db_rate, Dd_rate = build_birth_death_rates(de, RESDMG_ND_INDEX, bD, dD, :bD_j, :dD_j)
    Db_jump = VariableRateJump(Db_rate, D_birth!)
    Dd_jump = VariableRateJump(Dd_rate, D_death!)

    function R_birth!(integrator)
        integrator.u[RESDMG_NR_INDEX] += 1
        nothing
    end

    function R_death!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        nothing
    end

    Rb_rate, Rd_rate = build_birth_death_rates(de, RESDMG_NR_INDEX, bR, dR, :bR_j, :dR_j; psi_fn = p -> p.psi)
    Rb_jump = VariableRateJump(Rb_rate, R_birth!)
    Rd_jump = VariableRateJump(Rd_rate, R_death!)

    function E_birth!(integrator)
        integrator.u[RESDMG_NE_INDEX] += 1
        nothing
    end

    function E_death!(integrator)
        integrator.u[RESDMG_NE_INDEX] -= 1
        nothing
    end

    Eb_rate, Ed_rate = build_birth_death_rates(de, RESDMG_NE_INDEX, bE, dE, :bE_j, :dE_j)
    Eb_jump = VariableRateJump(Eb_rate, E_birth!)
    Ed_jump = VariableRateJump(Ed_rate, E_death!)

    function SR_switch!(integrator)
        integrator.u[RESDMG_NS_INDEX] -= 1
        integrator.u[RESDMG_NR_INDEX] += 1
        nothing
    end

    SR_switch_rate = build_switch_rate(de, RESDMG_NS_INDEX, p -> (p.bS_o + p.bS_j), p -> p.mu_j)
    SR_switch_jump = VariableRateJump(SR_switch_rate, SR_switch!)

    function SD_switch!(integrator)
        integrator.u[RESDMG_NS_INDEX] -= 1
        integrator.u[RESDMG_ND_INDEX] += 1
        nothing
    end

    function SD_switch_rate(u, p, t)
        return u[RESDMG_NS_INDEX] * p.ome_S_j * u[RESDMG_GAM_INDEX]
    end
    SD_switch_jump = VariableRateJump(SD_switch_rate, SD_switch!)

    function RS_switch!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        integrator.u[RESDMG_NS_INDEX] += 1
        nothing
    end

    RS_switch_rate = build_switch_rate(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), p -> p.sig_j; psi_fn = p -> p.psi)
    RS_switch_jump = VariableRateJump(RS_switch_rate, RS_switch!)

    function RD_switch!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        integrator.u[RESDMG_ND_INDEX] += 1
        nothing
    end

    function RD_switch_rate(u, p, t)
        return u[RESDMG_NR_INDEX] * p.ome_R_j * u[RESDMG_GAM_INDEX] * (1 - p.psi)
    end
    RD_switch_jump = VariableRateJump(RD_switch_rate, RD_switch!)

    function RE_switch!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        integrator.u[RESDMG_NE_INDEX] += 1
        nothing
    end

    RE_switch_rate = build_switch_rate(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), p -> p.al_j;
                                       psi_fn = p -> p.psi)
    RE_switch_jump = VariableRateJump(RE_switch_rate, RE_switch!)

    function DS_switch!(integrator)
        integrator.u[RESDMG_ND_INDEX] -= 1
        integrator.u[RESDMG_NS_INDEX] += 1
        nothing
    end

    function DS_switch_rate(u, p, t)
        return u[RESDMG_ND_INDEX] * p.zet_j
    end
    DS_switch_jump = VariableRateJump(DS_switch_rate, DS_switch!)

    mu = params.mu
    sig = params.sig
    al = params.al
    ome = params.ome
    zet = params.zet
    n0 = sim_eff.n0

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
        cb1 = DiscreteCallback(cond_to_ode, switch_to_ode!, save_positions = (false, true))

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
        cb2 = DiscreteCallback(cond_to_jump, switch_to_jump!, save_positions = (false, true))
        return cb1, cb2
    end

    function build_rate_toggle_callbacks(idx, rate_value, rate_o_sym, rate_j_sym, N_trans_switch;
                                         activity_fn = integrator -> (integrator.u[idx] * rate_value))
        cond_to_ode(u, t, integrator) = activity_fn(integrator) >= N_trans_switch
        function switch_to_ode!(integrator)
            if getproperty(integrator.p, rate_o_sym) == 0.0
                setproperty!(integrator.p, rate_o_sym, deepcopy(getproperty(integrator.p, rate_j_sym)))
                setproperty!(integrator.p, rate_j_sym, 0.0)
            end
            nothing
        end
        cb1 = DiscreteCallback(cond_to_ode, switch_to_ode!, save_positions = (false, true))

        cond_to_jump(u, t, integrator) = activity_fn(integrator) < N_trans_switch
        function switch_to_jump!(integrator)
            if getproperty(integrator.p, rate_j_sym) == 0.0
                setproperty!(integrator.p, rate_j_sym, deepcopy(getproperty(integrator.p, rate_o_sym)))
                setproperty!(integrator.p, rate_o_sym, 0.0)
            end
            nothing
        end
        cb2 = DiscreteCallback(cond_to_jump, switch_to_jump!, save_positions = (false, true))

        return cb1, cb2
    end

    function build_switch_activity_proxy(de, idx, b_total_fn, rate_value;
                                         gamma_factor = false, psi_fn = p -> 0.0, logistic_scale = true)
        function proxy(integrator)
            u = integrator.u
            p = integrator.p

            b_total = b_total_fn(p)
            drug_interaction = b_total * p.Dc * u[RESDMG_GAM_INDEX] * (1 - psi_fn(p))
            b_drug = b_total - drug_interaction

            b_effective = if de === :d
                b_total
            elseif de === :b
                b_drug
            else
                max(b_drug, 0.0)
            end

            activity = u[idx] * rate_value * b_effective
            if gamma_factor
                activity *= u[RESDMG_GAM_INDEX]
            end
            if logistic_scale
                activity *= max(resdmg_logistic_factor(u, sim_eff.Cc), 0.0)
            end
            return activity
        end

        return proxy
    end

    function resample_passage_bottleneck!(integrator, n0)
        curr_nS = round(Int, integrator.u[RESDMG_NS_INDEX])
        curr_nD = round(Int, integrator.u[RESDMG_ND_INDEX])
        curr_nR = round(Int, integrator.u[RESDMG_NR_INDEX])
        curr_nE = round(Int, integrator.u[RESDMG_NE_INDEX])

        counts = multivariate_hypergeometric_draw([curr_nS, curr_nD, curr_nR, curr_nE], n0)

        integrator.u[RESDMG_NS_INDEX] = Float64(counts[1])
        integrator.u[RESDMG_ND_INDEX] = Float64(counts[2])
        integrator.u[RESDMG_NR_INDEX] = Float64(counts[3])
        integrator.u[RESDMG_NE_INDEX] = Float64(counts[4])
        nothing
    end

    function attempt_passage!(integrator, n0; require_enough::Bool)
        if require_enough
            curr_nS = round(Int, integrator.u[RESDMG_NS_INDEX])
            curr_nD = round(Int, integrator.u[RESDMG_ND_INDEX])
            curr_nR = round(Int, integrator.u[RESDMG_NR_INDEX])
            curr_nE = round(Int, integrator.u[RESDMG_NE_INDEX])
            if (curr_nS + curr_nD + curr_nR + curr_nE) < n0
                return false
            end
        end

        resample_passage_bottleneck!(integrator, n0)
        integrator.u[RESDMG_PASS_INDEX] += 1
        return true
    end

    function build_clamp_callback(idx)
        condition(u, t, integrator) = u[idx] < 0
        function affect!(integrator)
            integrator.u[idx] = 0.0
        end
        return DiscreteCallback(condition, affect!, save_positions = (false, true))
    end

    S_cb_switch1, S_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NS_INDEX, :bS_o, :bS_j, :dS_o, :dS_j, sim_eff.Nswitch)
    D_cb_switch1, D_cb_switch2 = build_bd_toggle_callbacks(RESDMG_ND_INDEX, :bD_o, :bD_j, :dD_o, :dD_j, sim_eff.Nswitch)
    R_cb_switch1, R_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NR_INDEX, :bR_o, :bR_j, :dR_o, :dR_j, sim_eff.Nswitch)
    E_cb_switch1, E_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NE_INDEX, :bE_o, :bE_j, :dE_o, :dE_j, sim_eff.Nswitch)

    SR_activity_proxy = build_switch_activity_proxy(de, RESDMG_NS_INDEX, p -> (p.bS_o + p.bS_j), mu)
    SR_cb_switch1, SR_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NS_INDEX, mu, :mu_o, :mu_j, sim_eff.N_trans_switch;
                                                                activity_fn = SR_activity_proxy)

    SD_activity_proxy = integrator -> (integrator.u[RESDMG_NS_INDEX] * ome * integrator.u[RESDMG_GAM_INDEX])
    SD_cb_switch1, SD_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NS_INDEX, ome, :ome_S_o, :ome_S_j, sim_eff.N_trans_switch;
                                                                activity_fn = SD_activity_proxy)

    RS_activity_proxy = build_switch_activity_proxy(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), sig;
                                                    psi_fn = p -> p.psi)
    RS_cb_switch1, RS_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NR_INDEX, sig, :sig_o, :sig_j, sim_eff.N_trans_switch;
                                                                activity_fn = RS_activity_proxy)

    RD_activity_proxy = integrator -> (integrator.u[RESDMG_NR_INDEX] * ome * (1 - params.psi) * integrator.u[RESDMG_GAM_INDEX])
    RD_cb_switch1, RD_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NR_INDEX, ome, :ome_R_o, :ome_R_j, sim_eff.N_trans_switch;
                                                                activity_fn = RD_activity_proxy)

    RE_activity_proxy = build_switch_activity_proxy(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), al;
                                                    psi_fn = p -> p.psi)
    RE_cb_switch1, RE_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NR_INDEX, al, :al_o, :al_j, sim_eff.N_trans_switch;
                                                                activity_fn = RE_activity_proxy)

    DS_activity_proxy = integrator -> (integrator.u[RESDMG_ND_INDEX] * zet)
    DS_cb_switch1, DS_cb_switch2 = build_rate_toggle_callbacks(RESDMG_ND_INDEX, zet, :zet_o, :zet_j, sim_eff.N_trans_switch;
                                                                activity_fn = DS_activity_proxy)

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

    gam_above_max(u, t, integrator) = (u[RESDMG_GAM_INDEX] - 1.0)
    function clamp_gam_max!(integrator)
        integrator.u[RESDMG_GAM_INDEX] = 1.0
        integrator.p.kp = 0.0
    end

    gam_below_min(u, t, integrator) = (u[RESDMG_GAM_INDEX])
    function clamp_gam_min!(integrator)
        integrator.u[RESDMG_GAM_INDEX] = 0.0
        integrator.p.kp = 0.0
    end

    gam_max_cb = ContinuousCallback(gam_above_max, clamp_gam_max!, save_positions = (false, true))
    gam_min_cb = ContinuousCallback(gam_below_min, clamp_gam_min!, save_positions = (false, true))

    nmax_reached(u, t, integrator) = (resdmg_total_population(integrator.u) >= sim_eff.Nmax)

    function handle_nmax!(integrator)
        if integrator.u[RESDMG_PASS_INDEX] >= n_pass_eff
            terminate!(integrator)
        else
            attempt_passage!(integrator, n0; require_enough = false)
        end
    end

    nmax_cb = DiscreteCallback(nmax_reached, handle_nmax!, save_positions = (false, true))

    extinction(u, t, integrator) = (resdmg_total_population(integrator.u) < 1.0)
    extinction_cb = DiscreteCallback(extinction, terminate!, save_positions = (false, false))

    function apply_scheduled_passage!(integrator)
        curr_pass = round(Int, integrator.u[RESDMG_PASS_INDEX])
        if curr_pass < n_pass_eff
            expected_t = t_pass_vec[curr_pass]
            if integrator.t == expected_t
                attempt_passage!(integrator, n0; require_enough = true)
            end
        end
    end

    passage_times = filter_times_in_span(t_pass_vec)
    passage_cb = isempty(passage_times) ? nothing :
                 PresetTimeCallback(passage_times, apply_scheduled_passage!; save_positions = (false, true))

    cb7a = build_clamp_callback(RESDMG_NS_INDEX)
    cb7b = build_clamp_callback(RESDMG_ND_INDEX)
    cb7c = build_clamp_callback(RESDMG_NR_INDEX)
    cb7d = build_clamp_callback(RESDMG_NE_INDEX)

    cb_list = Any[
        S_cb_switch1, S_cb_switch2,
        D_cb_switch1, D_cb_switch2,
        R_cb_switch1, R_cb_switch2,
        E_cb_switch1, E_cb_switch2,
        SR_cb_switch1, SR_cb_switch2,
        SD_cb_switch1, SD_cb_switch2,
        RS_cb_switch1, RS_cb_switch2,
        RD_cb_switch1, RD_cb_switch2,
        RE_cb_switch1, RE_cb_switch2,
        DS_cb_switch1, DS_cb_switch2,
        nmax_cb, extinction_cb,
        cb7a, cb7b, cb7c, cb7d
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

    cbs = CallbackSet(cb_list...)

    sjm_prob = JumpProblem(prob, Direct(),
                           Sb_jump, Sd_jump,
                           Db_jump, Dd_jump,
                           Rb_jump, Rd_jump,
                           Eb_jump, Ed_jump,
                           SR_switch_jump,
                           SD_switch_jump,
                           RS_switch_jump,
                           RD_switch_jump,
                           RE_switch_jump,
                           DS_switch_jump)

    base_tstops = collect(sim_eff.t0:sim_eff.save_at:sim_eff.tmax)
    event_times = passage_times
    if sim_eff.treat == true
        event_times = vcat(event_times, treat_on_times, treat_off_times)
    end
    tstops = isempty(event_times) ? base_tstops : sort(unique(vcat(base_tstops, event_times)))

    sol = @suppress solve(sjm_prob, Tsit5(), callback = cbs, tstops = tstops)

    return sol
end






