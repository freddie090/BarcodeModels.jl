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
