"""Apply a ResDmg ABM birth event with resistance mutation/reversion."""
function resdmg_birth_mutate_event!(cell_arr::Vector{ResDmgCell},
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    if cell_arr[cell_pos].DS || cell_arr[cell_pos].DR
        throw(ArgumentError("Cannot perform birth event attempt on a damaged cell."))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].DS = cell_arr[cell_pos].DS
    cell_arr[birth_pos].DR = cell_arr[cell_pos].DR
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].alive = true

    mut_p = rand()

    if cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        else
            phen_counts.Rcount += 1
        end
    else
        if mu > mut_p
            cell_arr[birth_pos].R = true
            phen_counts.Rcount += 1
        else
            phen_counts.Scount += 1
        end
    end
end

"""Apply a lineage-aware ResDmg ABM birth event and record the child lineage."""
function resdmg_birth_mutate_event_evbc!(state::ResDmgABMEvBCState,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64,
    phen_counts::ResDmgPhenoCounts,
    curr_t::Float64)

    cell_arr = state.cells
    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    if cell_arr[cell_pos].DS || cell_arr[cell_pos].DR
        throw(ArgumentError("Cannot perform birth event attempt on a damaged cell."))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].DS = cell_arr[cell_pos].DS
    cell_arr[birth_pos].DR = cell_arr[cell_pos].DR
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].alive = true
    child_id = state.next_cell_id
    cell_arr[birth_pos].id = child_id
    cell_arr[birth_pos].parent_id = cell_arr[cell_pos].id
    cell_arr[birth_pos].birth_time = curr_t
    parent_pheno = _resdmg_pheno_label(cell_arr[cell_pos])

    mut_p = rand()

    if cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        else
            phen_counts.Rcount += 1
        end
    else
        if mu > mut_p
            cell_arr[birth_pos].R = true
            phen_counts.Rcount += 1
        else
            phen_counts.Scount += 1
        end
    end

    child_pheno = _resdmg_pheno_label(cell_arr[birth_pos])
    push!(state.lineage_records, LineageRecord(child_id, cell_arr[cell_pos].id, curr_t, parent_pheno, child_pheno, cell_arr[birth_pos].barcode))
    state.next_cell_id += 1
end

"""Move a ResDmg ABM cell into its damaged phenotype and update counts."""
function resdmg_damage_event!(cell_arr::Vector{<:Union{ResDmgCell, ResDmgCellEvBC}},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    if !cell_arr[cell_pos].alive
        return
    end

    if cell_arr[cell_pos].DS || cell_arr[cell_pos].DR
        throw(ArgumentError("Cell is already damaged"))
    end

    if cell_arr[cell_pos].R
        cell_arr[cell_pos].R = false
        cell_arr[cell_pos].DR = true
        phen_counts.Rcount -= 1
        phen_counts.DRcount += 1
    else
        cell_arr[cell_pos].DS = true
        phen_counts.Scount -= 1
        phen_counts.DScount += 1
    end
end

"""Repair a damaged ResDmg ABM cell and update counts."""
function resdmg_repair_event!(cell_arr::Vector{<:Union{ResDmgCell, ResDmgCellEvBC}},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    if !cell_arr[cell_pos].alive
        return
    end

    if cell_arr[cell_pos].DS
        cell_arr[cell_pos].DS = false
        phen_counts.DScount -= 1
        phen_counts.Scount += 1
    elseif cell_arr[cell_pos].DR
        cell_arr[cell_pos].DR = false
        cell_arr[cell_pos].R = true
        phen_counts.DRcount -= 1
        phen_counts.Rcount += 1
    else
        return
    end
end

"""Mark a ResDmg ABM cell dead and decrement its phenotype count."""
function resdmg_death_event!(cell_arr::Vector{<:Union{ResDmgCell, ResDmgCellEvBC}},
    cell_pos::Int64,
    phen_counts::ResDmgPhenoCounts)

    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    cell_arr[cell_pos].alive = false

    if cell_arr[cell_pos].R
        phen_counts.Rcount -= 1
    elseif cell_arr[cell_pos].DR
        phen_counts.DRcount -= 1
    elseif cell_arr[cell_pos].DS
        phen_counts.DScount -= 1
    else
        phen_counts.Scount -= 1
    end
end
