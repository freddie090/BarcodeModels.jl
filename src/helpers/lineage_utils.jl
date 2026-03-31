function _lineage_required_columns_present(lineage_df::DataFrame)
    required = ["id", "parent_id", "birth_time"]
    missing_cols = filter(c -> !(c in names(lineage_df)), required)
    isempty(missing_cols) || throw(ArgumentError("lineage_df is missing required columns: $(join(missing_cols, ", "))."))
    return nothing
end

function build_phylogeny(lineage_df::DataFrame)
    _lineage_required_columns_present(lineage_df)
    edges = Tuple{Int64, Int64}[]
    for row in eachrow(lineage_df)
        parent_id = Int64(row.parent_id)
        if parent_id != 0
            push!(edges, (parent_id, Int64(row.id)))
        end
    end
    return edges
end

function build_tree(lineage_df::DataFrame)
    _lineage_required_columns_present(lineage_df)
    children = Dict{Int64, Vector{Int64}}()

    for row in eachrow(lineage_df)
        parent_id = Int64(row.parent_id)
        if parent_id != 0
            push!(get!(children, parent_id, Int64[]), Int64(row.id))
        end
    end

    # Deterministic output helps testing and reproducible tree exports.
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

function lineage_to_newick(lineage_df::DataFrame, root_id::Int64)
    children = build_tree(lineage_df)
    return to_newick(root_id, children) * ";"
end

function population_to_newick(lineage_df::DataFrame, root_id::Int64)
    return lineage_to_newick(lineage_df, root_id)
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
    if "rep" in names(lineage_df)
        push!(cols, "rep")
    end

    return select(lineage_df, Symbol.(cols))
end
