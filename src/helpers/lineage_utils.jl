"""Build the public lineage DataFrame from lineage records and extant cell ids."""
function _lineage_df(records::Vector{LineageRecord}, rep::Int64, alive_ids::Vector{Int64} = Int64[])
    if isempty(records)
        return DataFrame(
            id = Int64[],
            parent_id = Int64[],
            birth_time = Float64[],
            parent_pheno = String[],
            child_pheno = String[],
            barcode = Float64[],
            alive_at_end = Bool[],
            rep = Int64[]
        )
    end
    ids = [r.id for r in records]
    parent_ids = [r.parent_id for r in records]
    birth_times = [r.birth_time for r in records]
    parent_phenos = [r.parent_pheno for r in records]
    child_phenos = [r.child_pheno for r in records]
    barcodes = [r.barcode for r in records]
    alive_id_set = Set(alive_ids)
    alive_at_end = [id in alive_id_set for id in ids]
    return DataFrame(
        id = ids,
        parent_id = parent_ids,
        birth_time = birth_times,
        parent_pheno = parent_phenos,
        child_pheno = child_phenos,
        barcode = barcodes,
        alive_at_end = alive_at_end,
        rep = fill(rep, length(records))
    )
end

"""Initialize cell ids and root lineage records for live EvBC cells."""
function initialize_lineage_state!(cells, phenotype_label::Function, t0::Float64 = 0.0)
    lineage_records = LineageRecord[]
    next_cell_id = Int64(1)
    for i in eachindex(cells)
        if cells[i].alive
            cells[i].id = next_cell_id
            cells[i].parent_id = 0
            cells[i].birth_time = t0
            push!(lineage_records, LineageRecord(next_cell_id, 0, t0, "ROOT", phenotype_label(cells[i]), cells[i].barcode))
            next_cell_id += 1
        else
            cells[i].id = 0
            cells[i].parent_id = 0
            cells[i].birth_time = -1.0
        end
    end
    return next_cell_id, lineage_records
end

function _lineage_required_columns_present(lineage_df::DataFrame)
    required = ["id", "parent_id", "birth_time"]
    missing_cols = filter(c -> !(c in names(lineage_df)), required)
    isempty(missing_cols) || throw(ArgumentError("lineage_df is missing required columns: $(join(missing_cols, ", "))."))
    return nothing
end

function _lineage_extant_ancestor_closure(lineage_df::DataFrame)
    "alive_at_end" in names(lineage_df) || throw(ArgumentError("lineage_df is missing required column when extant_only=true: alive_at_end."))

    parent_lookup = Dict{Int64, Int64}()
    for row in eachrow(lineage_df)
        parent_lookup[Int64(row.id)] = Int64(row.parent_id)
    end

    alive_ids = Set{Int64}()
    for row in eachrow(lineage_df)
        if Bool(row.alive_at_end)
            push!(alive_ids, Int64(row.id))
        end
    end

    keep_ids = Set{Int64}(alive_ids)
    stack = collect(alive_ids)
    while !isempty(stack)
        node_id = pop!(stack)
        parent_id = get(parent_lookup, node_id, 0)
        if parent_id != 0 && !(parent_id in keep_ids)
            push!(keep_ids, parent_id)
            push!(stack, parent_id)
        end
    end

    return keep_ids
end

function _lineage_rows(lineage_df::DataFrame; extant_only::Bool = false)
    if !extant_only
        return lineage_df
    end

    keep_ids = _lineage_extant_ancestor_closure(lineage_df)
    return filter(row -> Int64(row.id) in keep_ids, lineage_df)
end

function build_phylogeny(lineage_df::DataFrame; extant_only::Bool = false)
    _lineage_required_columns_present(lineage_df)
    lineage_rows = _lineage_rows(lineage_df; extant_only = extant_only)
    edges = Tuple{Int64, Int64}[]
    for row in eachrow(lineage_rows)
        parent_id = Int64(row.parent_id)
        if parent_id != 0
            push!(edges, (parent_id, Int64(row.id)))
        end
    end
    return edges
end

function build_tree(lineage_df::DataFrame; extant_only::Bool = false)
    _lineage_required_columns_present(lineage_df)
    lineage_rows = _lineage_rows(lineage_df; extant_only = extant_only)
    children = Dict{Int64, Vector{Int64}}()

    for row in eachrow(lineage_rows)
        parent_id = Int64(row.parent_id)
        if parent_id != 0
            push!(get!(children, parent_id, Int64[]), Int64(row.id))
        end
    end

    for key in keys(children)
        sort!(children[key])
    end

    return children
end

function to_newick(node::Int64, children::Dict{Int64, Vector{Int64}})
    if !haskey(children, node)
        return string(node)
    end

    subtrees = [to_newick(child, children) for child in children[node]]
    return "(" * join(subtrees, ",") * ")" * string(node)
end

function lineage_to_newick(lineage_df::DataFrame, root_id::Int64; extant_only::Bool = false)
    children = build_tree(lineage_df; extant_only = extant_only)
    return to_newick(root_id, children) * ";"
end

function population_to_newick(lineage_df::DataFrame, root_id::Int64; extant_only::Bool = false)
    return lineage_to_newick(lineage_df, root_id; extant_only = extant_only)
end

function lineage_edge_barcodes(lineage_df::DataFrame)
    _lineage_required_columns_present(lineage_df)
    "barcode" in names(lineage_df) || throw(ArgumentError("lineage_df is missing required column: barcode."))

    nodes = select(lineage_df, :id, :parent_id, :barcode)
    child_edges = filter(:parent_id => !=(0), nodes)
    parent_nodes = select(nodes, :id, :barcode)
    rename!(parent_nodes, :id => :parent_id, :barcode => :parent_barcode)

    edge_df = leftjoin(child_edges, parent_nodes, on = :parent_id)
    rename!(edge_df, :barcode => :child_barcode)
    select!(edge_df, :parent_id, :id, :parent_barcode, :child_barcode)
    return edge_df
end

function lineage_node_metadata(lineage_df::DataFrame)
    _lineage_required_columns_present(lineage_df)

    cols = ["id", "parent_id", "birth_time"]
    if "parent_pheno" in names(lineage_df)
        push!(cols, "parent_pheno")
    end
    if "child_pheno" in names(lineage_df)
        push!(cols, "child_pheno")
    end
    if "barcode" in names(lineage_df)
        push!(cols, "barcode")
    end
    if "alive_at_end" in names(lineage_df)
        push!(cols, "alive_at_end")
    end
    if "rep" in names(lineage_df)
        push!(cols, "rep")
    end

    return select(lineage_df, Symbol.(cols))
end
