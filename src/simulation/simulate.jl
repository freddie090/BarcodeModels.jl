function simulate_experiment(model::Union{HybridModel, ABMModel}, exp::ExperimentParams; kwargs...)
    if model isa HybridModel
        return simulate_experiment_hybrid(model, exp; kwargs...)
    end
    return simulate_experiment_abm(model, exp; kwargs...)
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

