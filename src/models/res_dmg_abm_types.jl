# Standard ABM types and setup

mutable struct ResDmgCell
    barcode::Float64
    DS::Bool
    DR::Bool
    R::Bool
    alive::Bool
end

mutable struct ResDmgPhenoCounts
    Scount::Int64
    DScount::Int64
    DRcount::Int64
    Rcount::Int64
end

mutable struct ResDmgGrowOut
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    Svec::Vector{Int64}
    DSvec::Vector{Int64}
    DRvec::Vector{Int64}
    Rvec::Vector{Int64}
    Pvec::Vector{Int64}
    fin_t::Float64
end

mutable struct ResDmgABMState
    cells::Vector{ResDmgCell}
end

make_dead_resdmg_cell() = ResDmgCell(0.0, false, false, false, false)

function seed_resdmg_cells(N::Int64, rho::Float64, Nbuff::Int64;
    skew_lib::Bool = false, use_lib_probs::Bool = false, bc_unif::Float64 = 0.0,
    Nbc::Int64 = 0, bc_probs = Float64[])

    cells = Vector{ResDmgCell}(undef, max(N, Nbuff))

    samp_bcs = _sample_initial_barcodes(N;
        skew_lib = skew_lib, use_lib_probs = use_lib_probs, bc_unif = bc_unif,
        Nbc = Nbc, bc_probs = bc_probs)

    for i in 1:N
        cells[i] = ResDmgCell(samp_bcs[i], false, false, false, true)
    end

    if rho > 0.0
        nR = Int64(round(rho * N))
        R_cells = sample(1:N, nR, replace = false)
        for i in 1:nR
            cells[R_cells[i]].R = true
        end
    end

    if Nbuff > N
        for i in (N + 1):Nbuff
            cells[i] = make_dead_resdmg_cell()
        end
    end

    return cells
end

# EvBC types and setup

mutable struct ResDmgCellEvBC
    barcode::Float64
    DS::Bool
    DR::Bool
    R::Bool
    alive::Bool
    id::Int64
    parent_id::Int64
    birth_time::Float64
end

mutable struct ResDmgABMEvBCState
    cells::Vector{ResDmgCellEvBC}
    next_cell_id::Int64
    lineage_records::Vector{LineageRecord}
end

make_dead_resdmg_cell_evbc() = ResDmgCellEvBC(0.0, false, false, false, false, 0, 0, -1.0)

function _resdmg_pheno_label(cell::ResDmgCellEvBC)
    if cell.DS
        return "DS"
    elseif cell.DR
        return "DR"
    elseif cell.R
        return "R"
    end
    return "S"
end

"""Initialise root lineage records using the model phenotype labels."""
function _init_lineage_state!(cells::Vector{ResDmgCellEvBC}, t0::Float64 = 0.0;
    start_id::Int64 = 1)
    return _init_lineage_state!(cells, _resdmg_pheno_label, t0; start_id = start_id)
end

function _to_evbc_cells(cells::Vector{ResDmgCell}; t0::Float64 = 0.0, start_id::Int64 = 1)
    evbc_cells = Vector{ResDmgCellEvBC}(undef, length(cells))
    for i in eachindex(cells)
        c = cells[i]
        evbc_cells[i] = ResDmgCellEvBC(c.barcode, c.DS, c.DR, c.R, c.alive, 0, 0, c.alive ? t0 : -1.0)
    end
    next_cell_id, lineage_records = _init_lineage_state!(evbc_cells, t0; start_id = start_id)
    return evbc_cells, next_cell_id, lineage_records
end

