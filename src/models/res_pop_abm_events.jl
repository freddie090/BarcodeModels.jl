# Standard ABM birth and shared events

"""Apply a ResPop ABM birth with a possible daughter phenotype change."""
function birth_mutate_event!(cell_arr,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts)

    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].E = cell_arr[cell_pos].E
    cell_arr[birth_pos].alive = true

    mut_p = rand()

    if cell_arr[birth_pos].E
        phen_counts.Ecount += 1
    elseif cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        elseif sig <= mut_p < (sig + al)
            cell_arr[birth_pos].R = false
            cell_arr[birth_pos].E = true
            phen_counts.Ecount += 1
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

# Independent conversion changes phenotype only, retaining barcode and cell identity.
# AbstractVector also serves the existing EvBC cell workflow.
function respop_S_to_R_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && !cell.R && !cell.E || return
    cell.R = true
    phen_counts.Scount -= 1
    phen_counts.Rcount += 1
end

function respop_R_to_S_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && cell.R && !cell.E || return
    cell.R = false
    phen_counts.Rcount -= 1
    phen_counts.Scount += 1
end

function respop_R_to_E_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && cell.R && !cell.E || return
    cell.R = false
    cell.E = true
    phen_counts.Rcount -= 1
    phen_counts.Ecount += 1
end

"""Mark a ResPop-family ABM cell dead and decrement its phenotype count."""
function death_event!(cell_arr, cell_pos::Int64, phen_counts)
    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    cell_arr[cell_pos].alive = false

    if cell_arr[cell_pos].E
        phen_counts.Ecount -= 1
    elseif cell_arr[cell_pos].R
        phen_counts.Rcount -= 1
    else
        phen_counts.Scount -= 1
    end
end

# EvBC birth event

function _birth_mutate_event_evbc!(state::ResPopABMEvBCState,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts::PhenoCounts,
    curr_t::Float64)

    cell_arr = state.cells
    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    mother_id = cell_arr[cell_pos].id
    mother_barcode = cell_arr[cell_pos].barcode
    mother_R = cell_arr[cell_pos].R
    mother_E = cell_arr[cell_pos].E
    mother_pheno = _respop_pheno_label(cell_arr[cell_pos])

    daughter_a_id = state.next_cell_id
    cell_arr[cell_pos].id = daughter_a_id
    cell_arr[cell_pos].parent_id = mother_id
    cell_arr[cell_pos].birth_time = curr_t
    push!(state.lineage_records, LineageRecord(daughter_a_id, mother_id, curr_t, mother_pheno, mother_pheno, mother_barcode))
    state.next_cell_id += 1

    cell_arr[birth_pos].barcode = mother_barcode
    cell_arr[birth_pos].R = mother_R
    cell_arr[birth_pos].E = mother_E
    cell_arr[birth_pos].alive = true
    daughter_b_id = state.next_cell_id
    cell_arr[birth_pos].id = daughter_b_id
    cell_arr[birth_pos].parent_id = mother_id
    cell_arr[birth_pos].birth_time = curr_t

    mut_p = rand()

    if cell_arr[birth_pos].E
        phen_counts.Ecount += 1
    elseif cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        elseif sig <= mut_p < (sig + al)
            cell_arr[birth_pos].R = false
            cell_arr[birth_pos].E = true
            phen_counts.Ecount += 1
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

    child_pheno = _respop_pheno_label(cell_arr[birth_pos])
    push!(state.lineage_records, LineageRecord(daughter_b_id, mother_id, curr_t, mother_pheno, child_pheno, cell_arr[birth_pos].barcode))
    state.next_cell_id += 1
end

