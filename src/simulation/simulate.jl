"""
    simulate_experiment(model, exp; kwargs...)

Run a full experiment. ABMs accept `lin_track`, `lin_track_rec_t` and `top_x`
keyword overrides of the corresponding `ExperimentParams` settings.

With tracking enabled, `out["lin_track_df"]` contains the largest `top_x`
barcodes at each fixed-interval snapshot, plus segment starts/ends and both
sides of instantaneous treatment. Columns identify `snapshot`, `t`, `cond`,
`rep`, `passage`, `event`, `bc`, barcode total `N`, full `population_N` and
phenotype counts. `snapshot` preserves ordering at equal times. Empty populations
have a missing `bc` and zero counts. Omitted barcodes may be below the top_x cutoff;
their absence does not establish extinction. `population_N - sum(N)` per snapshot
gives the population excluded by the cutoff.

POT uses `cond="POT"`, `rep=0`, `passage=0` and a continuous expansion clock,
including across expansion bottlenecks. Replicate clocks start at zero. POT
tracking is included whenever `lin_track=true`, independently of `inc_pot`.
"""
function simulate_experiment(model::Union{HybridModel, ABMModel}, exp::ExperimentParams; kwargs...)
    if model isa HybridModel
        return simulate_experiment_hybrid(model, exp; kwargs...)
    end
    tracker = _lin_track_recorder(model, exp, kwargs)
    out = simulate_experiment_abm(model, exp; kwargs..., lin_tracker = tracker)
    tracker === nothing || (out["lin_track_df"] = tracker.counts)
    return out
end

simulate_experiment_hybrid(model::HybridModel, exp::ExperimentParams; kwargs...) =
    error("simulate_experiment_hybrid is not implemented for $(typeof(model)).")

simulate_experiment_abm(model::ABMModel, exp::ExperimentParams; kwargs...) =
    error("simulate_experiment_abm is not implemented for $(typeof(model)).")

simulate_experiment_hybrid(model::ResPop, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_hybrid(model, exp; kwargs...)

simulate_experiment_hybrid(model::ResDmg, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_hybrid(model, exp; kwargs...)

simulate_experiment_abm(model::ResPop_ABM, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_abm(model, exp; kwargs...)

simulate_experiment_abm(model::ResDmg_ABM, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_abm(model, exp; kwargs...)

simulate_experiment_abm(model::ResPop_ABM_EvBC, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_abm(model, exp; kwargs...)

simulate_experiment_abm(model::ResDmg_ABM_EvBC, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_abm(model, exp; kwargs...)

"""
    simulate_simple(model, sim; kwargs...)

Run a single population. ABMs accept `lin_track=false`, `lin_track_rec_t=0.1`
and `top_x=100` through `SimpleSimParams` or keyword overrides. Tracking adds
`lin_track_df` with the schema described in [`simulate_experiment`](@ref), using
`cond="DT"`, `rep=1`, `passage=1`. Recording does not consume random draws.
"""
function simulate_simple(model::Union{HybridModel, ABMModel}, sim::SimpleSimParams; kwargs...)
    if model isa HybridModel
        return simulate_simple_hybrid(model, sim; kwargs...)
    end
    tracker = _lin_track_recorder(model, sim, kwargs)
    out = simulate_simple_abm(model, sim; kwargs..., lin_tracker = tracker)
    tracker === nothing || (out["lin_track_df"] = tracker.counts)
    return out
end

simulate_simple_hybrid(model::HybridModel, sim::SimpleSimParams; kwargs...) =
    error("simulate_simple_hybrid is not implemented for $(typeof(model)).")

simulate_simple_abm(model::ABMModel, sim::SimpleSimParams; kwargs...) =
    error("simulate_simple_abm is not implemented for $(typeof(model)).")

simulate_simple_hybrid(model::ResPop, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_hybrid(model, sim; kwargs...)

simulate_simple_hybrid(model::ResDmg, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_hybrid(model, sim; kwargs...)

simulate_simple_abm(model::ResPop_ABM, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_abm(model, sim; kwargs...)

simulate_simple_abm(model::ResDmg_ABM, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_abm(model, sim; kwargs...)

simulate_simple_abm(model::ResPop_ABM_EvBC, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_abm(model, sim; kwargs...)

simulate_simple_abm(model::ResDmg_ABM_EvBC, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_abm(model, sim; kwargs...)

simulate_experiment_hybrid(model::ResPopInVivo, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_hybrid(model, exp; kwargs...)
simulate_experiment_abm(model::ResPopInVivo_ABM, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_abm(model, exp; kwargs...)
simulate_simple_hybrid(model::ResPopInVivo, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_hybrid(model, sim; kwargs...)
simulate_simple_abm(model::ResPopInVivo_ABM, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_abm(model, sim; kwargs...)
