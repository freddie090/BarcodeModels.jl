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
function build_component_params(params::ResDmgParams)
    bS = params.b
    dS = params.d

    bDS = 0.0
    dDS = params.d + params.Dc

    bDR = 0.0
    dDR = params.d + params.Dc

    bR = params.b * (1 - params.del)
    dR = params.d * (1 - params.del)

    p = ComponentArray(
        bS_o = 0.0, dS_o = 0.0, bS_j = bS, dS_j = dS,
        bDS_o = 0.0, dDS_o = 0.0, bDS_j = bDS, dDS_j = dDS,
        bDR_o = 0.0, dDR_o = 0.0, bDR_j = bDR, dDR_j = dDR,
        bR_o = 0.0, dR_o = 0.0, bR_j = bR, dR_j = dR,
        Dc = params.Dc, kp = 0.0, psi = params.psi,
        mu_o = 0.0, sig_o = 0.0,
        ome_S_o = 0.0, ome_R_o = 0.0, zet_S_o = 0.0, zet_R_o = 0.0,
        mu_j = params.mu, sig_j = params.sig,
        ome_S_j = params.ome, ome_R_j = params.ome, zet_S_j = params.zet_S, zet_R_j = params.zet_R
    )

    return p, bS, dS, bDS, dDS, bDR, dDR, bR, dR
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

    p, bS, dS, bDS, dDS, bDR, dDR, bR, dR = build_component_params(params)

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
        @unpack kp, psi, mu_o, sig_o, ome_S_o, ome_R_o, zet_S_o, zet_R_o = p
        @unpack gam, nS, nDS, nDR, nR = u

        N = nS + nDS + nDR + nR

        bS_o_mod, dS_o_mod = rate_modifier(bS, dS, p.bS_o, p.dS_o, p.Dc, gam, N, sim_eff.Cc, 0.0)
        bR_o_mod, dR_o_mod = rate_modifier(bR, dR, p.bR_o, p.dR_o, p.Dc, gam, N, sim_eff.Cc, psi)

        bDS_o_mod = p.bDS_o * max(1 - (N / sim_eff.Cc), 0.0)
        dDS_o_mod = p.dDS_o * max(1 - (N / sim_eff.Cc), 0.0)
        bDR_o_mod = p.bDR_o * max(1 - (N / sim_eff.Cc), 0.0)
        dDR_o_mod = p.dDR_o * max(1 - (N / sim_eff.Cc), 0.0)

        du.gam = kp

        du.nS = (bS_o_mod - dS_o_mod) * nS -
            ome_S_o * gam * nS -
            mu_o * nS * bS_o_mod +
            sig_o * nR * bR_o_mod +
            zet_S_o * nDS

        du.nR = (bR_o_mod - dR_o_mod) * nR +
            mu_o * nS * bS_o_mod -
            sig_o * nR * bR_o_mod -
            ome_R_o * gam * (1 - psi) * nR +
            zet_R_o * nDR

        du.nDS = (bDS_o_mod - dDS_o_mod) * nDS +
            ome_S_o * gam * nS -
            zet_S_o * nDS

        du.nDR = (bDR_o_mod - dDR_o_mod) * nDR +
            ome_R_o * gam * (1 - psi) * nR -
            zet_R_o * nDR
            
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

    function DS_birth!(integrator)
        integrator.u[RESDMG_NDS_INDEX] += 1
        nothing
    end

    function DS_death!(integrator)
        integrator.u[RESDMG_NDS_INDEX] -= 1
        nothing
    end

    DSb_rate, DSd_rate = build_birth_death_rates(de, RESDMG_NDS_INDEX, bDS, dDS, :bDS_j, :dDS_j)
    DSb_jump = VariableRateJump(DSb_rate, DS_birth!)
    DSd_jump = VariableRateJump(DSd_rate, DS_death!)

    function DR_birth!(integrator)
        integrator.u[RESDMG_NDR_INDEX] += 1
        nothing
    end

    function DR_death!(integrator)
        integrator.u[RESDMG_NDR_INDEX] -= 1
        nothing
    end

    DRb_rate, DRd_rate = build_birth_death_rates(de, RESDMG_NDR_INDEX, bDR, dDR, :bDR_j, :dDR_j)
    DRb_jump = VariableRateJump(DRb_rate, DR_birth!)
    DRd_jump = VariableRateJump(DRd_rate, DR_death!)

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

    function StoR_switch!(integrator)
        integrator.u[RESDMG_NS_INDEX] -= 1
        integrator.u[RESDMG_NR_INDEX] += 1
        nothing
    end

    StoR_switch_rate = build_switch_rate(de, RESDMG_NS_INDEX, p -> (p.bS_o + p.bS_j), p -> p.mu_j)
    StoR_switch_jump = VariableRateJump(StoR_switch_rate, StoR_switch!)

    function StoDS_switch!(integrator)
        integrator.u[RESDMG_NS_INDEX] -= 1
        integrator.u[RESDMG_NDS_INDEX] += 1
        nothing
    end

    function StoDS_switch_rate(u, p, t)
        return u[RESDMG_NS_INDEX] * p.ome_S_j * u[RESDMG_GAM_INDEX]
    end
    StoDS_switch_jump = VariableRateJump(StoDS_switch_rate, StoDS_switch!)

    function RtoS_switch!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        integrator.u[RESDMG_NS_INDEX] += 1
        nothing
    end

    RtoS_switch_rate = build_switch_rate(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), p -> p.sig_j; psi_fn = p -> p.psi)
    RtoS_switch_jump = VariableRateJump(RtoS_switch_rate, RtoS_switch!)

    function RtoDR_switch!(integrator)
        integrator.u[RESDMG_NR_INDEX] -= 1
        integrator.u[RESDMG_NDR_INDEX] += 1
        nothing
    end

    function RtoDR_switch_rate(u, p, t)
        return u[RESDMG_NR_INDEX] * p.ome_R_j * u[RESDMG_GAM_INDEX] * (1 - p.psi)
    end
    RtoDR_switch_jump = VariableRateJump(RtoDR_switch_rate, RtoDR_switch!)

    function DStoS_switch!(integrator)
        integrator.u[RESDMG_NDS_INDEX] -= 1
        integrator.u[RESDMG_NS_INDEX] += 1
        nothing
    end

    function DStoS_switch_rate(u, p, t)
        return u[RESDMG_NDS_INDEX] * p.zet_S_j
    end
    DStoS_switch_jump = VariableRateJump(DStoS_switch_rate, DStoS_switch!)

    function DRtoR_switch!(integrator)
        integrator.u[RESDMG_NDR_INDEX] -= 1
        integrator.u[RESDMG_NR_INDEX] += 1
        nothing
    end

    function DRtoR_switch_rate(u, p, t)
        return u[RESDMG_NDR_INDEX] * p.zet_R_j
    end
    DRtoR_switch_jump = VariableRateJump(DRtoR_switch_rate, DRtoR_switch!)

    mu = params.mu
    sig = params.sig
    ome = params.ome
    zet_S = params.zet_S
    zet_R = params.zet_R
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
        curr_nDS = round(Int, integrator.u[RESDMG_NDS_INDEX])
        curr_nDR = round(Int, integrator.u[RESDMG_NDR_INDEX])
        curr_nR = round(Int, integrator.u[RESDMG_NR_INDEX])

        counts = multivariate_hypergeometric_draw([curr_nS, curr_nDS, curr_nDR, curr_nR], n0)

        integrator.u[RESDMG_NS_INDEX] = Float64(counts[1])
        integrator.u[RESDMG_NDS_INDEX] = Float64(counts[2])
        integrator.u[RESDMG_NDR_INDEX] = Float64(counts[3])
        integrator.u[RESDMG_NR_INDEX] = Float64(counts[4])
        nothing
    end

    function attempt_passage!(integrator, n0; require_enough::Bool)
        if require_enough
            curr_nS = round(Int, integrator.u[RESDMG_NS_INDEX])
            curr_nDS = round(Int, integrator.u[RESDMG_NDS_INDEX])
            curr_nDR = round(Int, integrator.u[RESDMG_NDR_INDEX])
            curr_nR = round(Int, integrator.u[RESDMG_NR_INDEX])
            if (curr_nS + curr_nDS + curr_nDR + curr_nR) < n0
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
    DS_cb_switch1, DS_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NDS_INDEX, :bDS_o, :bDS_j, :dDS_o, :dDS_j, sim_eff.Nswitch)
    DR_cb_switch1, DR_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NDR_INDEX, :bDR_o, :bDR_j, :dDR_o, :dDR_j, sim_eff.Nswitch)
    R_cb_switch1, R_cb_switch2 = build_bd_toggle_callbacks(RESDMG_NR_INDEX, :bR_o, :bR_j, :dR_o, :dR_j, sim_eff.Nswitch)

    StoR_activity_proxy = build_switch_activity_proxy(de, RESDMG_NS_INDEX, p -> (p.bS_o + p.bS_j), mu)
    StoR_cb_switch1, StoR_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NS_INDEX, mu, :mu_o, :mu_j, sim_eff.N_trans_switch;
                                                                    activity_fn = StoR_activity_proxy)

    StoDS_activity_proxy = integrator -> (integrator.u[RESDMG_NS_INDEX] * ome * integrator.u[RESDMG_GAM_INDEX])
    StoDS_cb_switch1, StoDS_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NS_INDEX, ome, :ome_S_o, :ome_S_j, sim_eff.N_trans_switch;
                                                                      activity_fn = StoDS_activity_proxy)

    RtoS_activity_proxy = build_switch_activity_proxy(de, RESDMG_NR_INDEX, p -> (p.bR_o + p.bR_j), sig;
                                                      psi_fn = p -> p.psi)
    RtoS_cb_switch1, RtoS_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NR_INDEX, sig, :sig_o, :sig_j, sim_eff.N_trans_switch;
                                                                    activity_fn = RtoS_activity_proxy)

    RtoDR_activity_proxy = integrator -> (integrator.u[RESDMG_NR_INDEX] * ome * (1 - params.psi) * integrator.u[RESDMG_GAM_INDEX])
    RtoDR_cb_switch1, RtoDR_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NR_INDEX, ome, :ome_R_o, :ome_R_j, sim_eff.N_trans_switch;
                                                                      activity_fn = RtoDR_activity_proxy)

    DStoS_activity_proxy = integrator -> (integrator.u[RESDMG_NDS_INDEX] * zet_S)
    DStoS_cb_switch1, DStoS_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NDS_INDEX, zet_S, :zet_S_o, :zet_S_j, sim_eff.N_trans_switch;
                                                                      activity_fn = DStoS_activity_proxy)

    DRtoR_activity_proxy = integrator -> (integrator.u[RESDMG_NDR_INDEX] * zet_R)
    DRtoR_cb_switch1, DRtoR_cb_switch2 = build_rate_toggle_callbacks(RESDMG_NDR_INDEX, zet_R, :zet_R_o, :zet_R_j, sim_eff.N_trans_switch;
                                                                      activity_fn = DRtoR_activity_proxy)

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
    cb7b = build_clamp_callback(RESDMG_NDS_INDEX)
    cb7c = build_clamp_callback(RESDMG_NDR_INDEX)
    cb7d = build_clamp_callback(RESDMG_NR_INDEX)

    cb_list = Any[
        S_cb_switch1, S_cb_switch2,
        DS_cb_switch1, DS_cb_switch2,
        DR_cb_switch1, DR_cb_switch2,
        R_cb_switch1, R_cb_switch2,
        StoR_cb_switch1, StoR_cb_switch2,
        StoDS_cb_switch1, StoDS_cb_switch2,
        RtoS_cb_switch1, RtoS_cb_switch2,
        RtoDR_cb_switch1, RtoDR_cb_switch2,
        DStoS_cb_switch1, DStoS_cb_switch2,
        DRtoR_cb_switch1, DRtoR_cb_switch2,
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
                           DSb_jump, DSd_jump,
                           DRb_jump, DRd_jump,
                           Rb_jump, Rd_jump,
                           StoR_switch_jump,
                           StoDS_switch_jump,
                           RtoS_switch_jump,
                           RtoDR_switch_jump,
                           DStoS_switch_jump,
                           DRtoR_switch_jump)

    base_tstops = collect(sim_eff.t0:sim_eff.save_at:sim_eff.tmax)
    event_times = passage_times
    if sim_eff.treat == true
        event_times = vcat(event_times, treat_on_times, treat_off_times)
    end
    tstops = isempty(event_times) ? base_tstops : sort(unique(vcat(base_tstops, event_times)))

    sol = @suppress solve(sjm_prob, Tsit5(), callback = cbs, tstops = tstops)

    return sol
end







