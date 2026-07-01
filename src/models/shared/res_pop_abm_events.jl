"""Apply a ResPop ABM birth event with mutation/switching into the target dead slot."""
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
