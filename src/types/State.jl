# State types

struct ModelState
    gam::Float64
    nS::Float64
    nR::Float64
    nE::Float64
    pass_num::Int64
end

function ModelState(nS::Real, nR::Real, nE::Real; gam=0.0, pass_num=1)
    return ModelState(Float64(gam), Float64(nS), Float64(nR), Float64(nE), Int64(pass_num))
end

function to_componentarray(state::ModelState)
    return ComponentArray(
        gam = state.gam,
        nS = state.nS,
        nR = state.nR,
        nE = state.nE,
        Pass_num = state.pass_num
    )
end
