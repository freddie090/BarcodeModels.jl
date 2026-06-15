# State types

struct ResPopState
    gam::Float64
    nS::Float64
    nR::Float64
    nE::Float64
    pass_num::Int64
end

function ResPopState(nS::Real, nR::Real, nE::Real; gam=0.0, pass_num=1)
    return ResPopState(Float64(gam), Float64(nS), Float64(nR), Float64(nE), Int64(pass_num))
end

struct ResPopInVivoState
    gam::Float64
    nS_EG0::Float64
    nS_EG1::Float64
    nR_EG0::Float64
    nR_EG1::Float64
    nE_EG0::Float64
    nE_EG1::Float64
    pass_num::Int64
end

function ResPopInVivoState(nS_EG0::Real, nS_EG1::Real, nR_EG0::Real, nR_EG1::Real,
    nE_EG0::Real, nE_EG1::Real; gam=0.0, pass_num=1)

    return ResPopInVivoState(
        Float64(gam),
        Float64(nS_EG0), Float64(nS_EG1),
        Float64(nR_EG0), Float64(nR_EG1),
        Float64(nE_EG0), Float64(nE_EG1),
        Int64(pass_num)
    )
end

struct ResDmgState
    gam::Float64
    nS::Float64
    nDS::Float64
    nDR::Float64
    nR::Float64
    pass_num::Int64
end

function ResDmgState(nS::Real, nDS::Real, nDR::Real, nR::Real; gam=0.0, pass_num=1)
    return ResDmgState(Float64(gam), Float64(nS), Float64(nDS), Float64(nDR), Float64(nR), Int64(pass_num))
end

"""Shared lineage record used by EvBC model variants and utilities."""
struct LineageRecord
    id::Int64
    parent_id::Int64
    birth_time::Float64
    parent_pheno::String
    child_pheno::String
    barcode::Float64
end

function to_componentarray(state::ResPopState)
    return ComponentArray(
        gam = state.gam,
        nS = state.nS,
        nR = state.nR,
        nE = state.nE,
        Pass_num = state.pass_num
    )
end

function to_componentarray(state::ResPopInVivoState)
    return ComponentArray(
        gam = state.gam,
        nS_EG0 = state.nS_EG0,
        nS_EG1 = state.nS_EG1,
        nR_EG0 = state.nR_EG0,
        nR_EG1 = state.nR_EG1,
        nE_EG0 = state.nE_EG0,
        nE_EG1 = state.nE_EG1,
        Pass_num = state.pass_num
    )
end

function to_componentarray(state::ResDmgState)
    return ComponentArray(
        gam = state.gam,
        nS = state.nS,
        nDS = state.nDS,
        nDR = state.nDR,
        nR = state.nR,
        Pass_num = state.pass_num
    )
end
