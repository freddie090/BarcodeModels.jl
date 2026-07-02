using Test
using BarcodeModels
using DataFrames
using Random

function simulate_hybrid_experiment(model)
    exp = ExperimentParams(
        n0 = 10,
        t_exp = 8.0,
        tmax = 10.0,
        t_Pass = Float64[],
        Nseed = 10,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [10.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 4
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment(model::ResPop_ABM)
    exp = ExperimentParams(
        n0 = 5,
        t_exp = 1.0,
        tmax = 2.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 20,
        Cc = 20,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        n_rep = 1
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment(model::ResPop_ABM_EvBC)
    exp = ExperimentParams(
        n0 = 5,
        t_exp = 1.0,
        tmax = 2.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 20,
        Cc = 20,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        n_rep = 1
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment(model::ResDmg_ABM)
    exp = ExperimentParams(
        n0 = 5,
        t_exp = 1.0,
        tmax = 2.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 20,
        Cc = 20,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 1
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment(model::ResDmg_ABM_EvBC)
    exp = ExperimentParams(
        n0 = 5,
        t_exp = 1.0,
        tmax = 2.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 20,
        Cc = 20,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 1
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment(model::ResPopInVivo_ABM)
    exp = ExperimentParams(
        n0 = 8,
        t_exp = 1.0,
        tmax = 2.0,
        t_Pass = [1.0],
        Nseed = 4,
        Nmax = 20,
        Cc = 20,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 1
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_resdmg_experiment(model::ResDmg)
    exp = ExperimentParams(
        n0 = 10,
        t_exp = 8.0,
        tmax = 10.0,
        t_Pass = Float64[],
        Nseed = 10,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [10.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_hybrid_experiment_vector_tmax(model)
    exp = ExperimentParams(
        n0 = 10,
        t_exp = 4.0,
        tmax = [2.0, 4.0],
        t_Pass = Float64[],
        Nseed = 10,
        Nmax = 1000,
        Cc = 1000,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0, 4.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_abm_experiment_vector_tmax(model)
    exp = ExperimentParams(
        n0 = 20,
        t_exp = 1.0,
        tmax = [1.0, 2.0],
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 50,
        Cc = 50,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [1.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2
    )
    return BarcodeModels.simulate_experiment(model, exp)
end

function simulate_simple_run(model)
    sim = SimpleSimParams(
        n0 = 10,
        tmax = 2.0,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[]
    )
    return BarcodeModels.simulate_simple(model, sim)
end

function simulate_simple_run_treated(model)
    sim = SimpleSimParams(
        n0 = 10,
        tmax = 2.0,
        Nmax = 100,
        Cc = 100,
        treat_ons = [0.5],
        treat_offs = [1.5]
    )
    return BarcodeModels.simulate_simple(model, sim)
end

@testset "BarcodeModels integration" begin
    params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    model = ResPop(params)
    @test model isa ResPop

    result = simulate_hybrid_experiment(model)
    @test result !== nothing
end

@testset "BarcodeModels integration (ABM)" begin
    params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5)
    model = ResPop_ABM(params; abm = abm)
    @test model isa ResPop_ABM

    result = simulate_abm_experiment(model)
    @test result !== nothing
end

@testset "ABM barcode library probabilities" begin
    @test_throws ErrorException ABMParams(
        use_lib_probs = true,
        skew_lib = true,
        Nbc = 2,
        bc_probs = [0.5, 0.5]
    )

    @test_throws ErrorException ABMParams(
        use_lib_probs = true,
        Nbc = 2,
        bc_probs = [1.0]
    )

    @test_throws ErrorException ABMParams(
        use_lib_probs = true,
        Nbc = 2,
        bc_probs = [0.7, 0.4]
    )

    abm = ABMParams(Nbuff = 10, use_lib_probs = true, Nbc = 2, bc_probs = [1.0, 0.0])

    respop_cells = BarcodeModels.seed_cells(4, 0.0, 4;
                                            use_lib_probs = abm.use_lib_probs,
                                            Nbc = abm.Nbc,
                                            bc_probs = abm.bc_probs)
    resdmg_cells = BarcodeModels.seed_resdmg_cells(4, 0.0, 4;
                                                   use_lib_probs = abm.use_lib_probs,
                                                   Nbc = abm.Nbc,
                                                   bc_probs = abm.bc_probs)
    invivo_cells = BarcodeModels.seed_invivo_cells(4, 0.0, 0.5, 4;
                                                   use_lib_probs = abm.use_lib_probs,
                                                   Nbc = abm.Nbc,
                                                   bc_probs = abm.bc_probs)

    @test all(cell -> cell.barcode == 1.0, respop_cells)
    @test all(cell -> cell.barcode == 1.0, resdmg_cells)
    @test all(cell -> cell.barcode == 1.0, invivo_cells)
    @test all(cell -> 1.0 <= cell.barcode <= 2.0, respop_cells)
    @test all(cell -> 1.0 <= cell.barcode <= 2.0, resdmg_cells)
    @test all(cell -> 1.0 <= cell.barcode <= 2.0, invivo_cells)
end

@testset "ABM phenotype-by-barcode output" begin
    default_abm = ABMParams(Nbuff = 20)
    @test default_abm.full_pheno_bc == false

    full_abm = ABMParams(Nbuff = 20, full_pheno_bc = true)
    @test full_abm.full_pheno_bc == true

    respop_cells = [
        BarcodeModels.CancerCell(1.0, false, false, true),
        BarcodeModels.CancerCell(1.0, true, false, true),
        BarcodeModels.CancerCell(1.0, false, true, true),
        BarcodeModels.CancerCell(2.0, false, false, true),
        BarcodeModels.CancerCell(2.0, true, false, false)
    ]
    respop_counts = BarcodeModels.get_pheno_counts(BarcodeModels.alive_cells(respop_cells), "DT1_P1")
    @test all(in(names(respop_counts)).(["bc", "DT1_P1_S", "DT1_P1_R", "DT1_P1_E"]))
    bc1 = respop_counts[respop_counts.bc .== 1.0, :][1, :]
    @test bc1.DT1_P1_S == 1
    @test bc1.DT1_P1_R == 1
    @test bc1.DT1_P1_E == 1

    invivo_cells = [
        BarcodeModels.InVivoCancerCell(1.0, false, false, false, true),
        BarcodeModels.InVivoCancerCell(1.0, false, false, true, true),
        BarcodeModels.InVivoCancerCell(1.0, true, false, false, true),
        BarcodeModels.InVivoCancerCell(2.0, false, true, true, true)
    ]
    invivo_counts = BarcodeModels.get_pheno_counts(invivo_cells, "DT1_P1")
    @test all(in(names(invivo_counts)).(["bc", "DT1_P1_S", "DT1_P1_R", "DT1_P1_E"]))
    @test !any(contains("EG"), string.(names(invivo_counts)))
    invivo_bc1 = invivo_counts[invivo_counts.bc .== 1.0, :][1, :]
    @test invivo_bc1.DT1_P1_S == 2
    @test invivo_bc1.DT1_P1_R == 1
    @test invivo_bc1.DT1_P1_E == 0

    resdmg_cells = [
        BarcodeModels.ResDmgCell(1.0, false, false, false, true),
        BarcodeModels.ResDmgCell(1.0, true, false, false, true),
        BarcodeModels.ResDmgCell(1.0, false, true, false, true),
        BarcodeModels.ResDmgCell(1.0, false, false, true, true)
    ]
    resdmg_counts = BarcodeModels.get_pheno_counts(resdmg_cells, "DT1_P1")
    @test all(in(names(resdmg_counts)).(["bc", "DT1_P1_S", "DT1_P1_DS", "DT1_P1_DR", "DT1_P1_R"]))
    resdmg_bc1 = resdmg_counts[resdmg_counts.bc .== 1.0, :][1, :]
    @test resdmg_bc1.DT1_P1_S == 1
    @test resdmg_bc1.DT1_P1_DS == 1
    @test resdmg_bc1.DT1_P1_DR == 1
    @test resdmg_bc1.DT1_P1_R == 1
end

@testset "ABM split-after-barcoding mode" begin
    @test_throws ErrorException ABMParams(
        split_after_barcoding = true,
        use_lib_probs = false
    )

    abm = ABMParams(
        Nbuff = 20,
        use_lib_probs = true,
        split_after_barcoding = true,
        Nbc = 3,
        bc_probs = [0.7, 0.2, 0.1]
    )
    @test abm.split_after_barcoding

    params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.2,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )
    model = ResPop_ABM(params; abm = abm)
    exp = ExperimentParams(
        n0 = 9,
        t_exp = 1.0,
        tmax = 1.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = Float64[],
        Nswitch = 10,
        n_rep = 2
    )
    @test_throws ErrorException BarcodeModels._expand_split_cells_abm(model, exp, 2)

    # Uniform over unique barcodes: with a 99:1 clone-size imbalance across two barcodes,
    # selecting one barcode should still be close to 50/50 across repeated draws.
    trials = 200
    picked_minor = 0
    Random.seed!(1234)
    for _ in 1:trials
        cells = BarcodeModels.CancerCell[]
        for _ in 1:99
            push!(cells, BarcodeModels.CancerCell(1.0, false, false, true))
        end
        push!(cells, BarcodeModels.CancerCell(2.0, false, false, true))
        BarcodeModels._assign_resistance_by_barcode!(cells, 0.5)
        resistant_barcodes = unique([cell.barcode for cell in cells if cell.R])
        @test length(resistant_barcodes) == 1
        if resistant_barcodes[1] == 2.0
            picked_minor += 1
        end
    end
    @test picked_minor / trials > 0.35
    @test picked_minor / trials < 0.65

    invivo_cells = BarcodeModels.InVivoCancerCell[]
    for _ in 1:3
        push!(invivo_cells, BarcodeModels.InVivoCancerCell(1.0, false, false, false, true))
    end
    for _ in 1:2
        push!(invivo_cells, BarcodeModels.InVivoCancerCell(2.0, false, false, false, true))
    end
    for _ in 1:4
        push!(invivo_cells, BarcodeModels.InVivoCancerCell(3.0, false, false, false, true))
    end
    BarcodeModels._assign_resistance_by_barcode!(invivo_cells, 0.34)
    BarcodeModels._assign_engraftment_by_barcode!(invivo_cells, 0.67)
    resistant_bcs = unique([c.barcode for c in invivo_cells if c.R])
    eg1_bcs = unique([c.barcode for c in invivo_cells if c.EG])
    @test length(resistant_bcs) == 1
    @test length(eg1_bcs) == 2
    for bc in unique([c.barcode for c in invivo_cells])
        bc_cells = [c for c in invivo_cells if c.barcode == bc]
        @test length(unique([c.R for c in bc_cells])) == 1
        @test length(unique([c.EG for c in bc_cells])) == 1
    end
end

@testset "BarcodeModels integration (ResDmg hybrid)" begin
    params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.01,
        zet_R = 0.01,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    model = ResDmg(params)
    @test model isa ResDmg

    result = simulate_resdmg_experiment(model)
    @test result !== nothing
    @test haskey(result, "t")
    @test haskey(result, "u")
    @test haskey(result, "sol_df")
    @test "nDS" in names(result["sol_df"])
    @test "nDR" in names(result["sol_df"])
end

@testset "BarcodeModels integration (ResDmg ABM)" begin
    params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.01,
        zet_R = 0.01,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5)
    model = ResDmg_ABM(params; abm = abm)
    @test model isa ResDmg_ABM

    result = simulate_abm_experiment(model)
    @test result !== nothing
    @test haskey(result, "sol_df")
    @test "nDS" in names(result["sol_df"])
    @test "nDR" in names(result["sol_df"])
end

@testset "BarcodeModels integration (ResPop ABM EvBC)" begin
    params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5)
    model = ResPop_ABM_EvBC(params; abm = abm)
    @test model isa ResPop_ABM_EvBC

    result = simulate_abm_experiment(model)
    @test result !== nothing
    @test haskey(result, "lineage_df")
    @test all(in(names(result["lineage_df"])).(["id", "parent_id", "birth_time"]))
    @test all(in(names(result["lineage_df"])).(["parent_pheno", "child_pheno"]))
    @test "barcode" in names(result["lineage_df"])
    @test "alive_at_end" in names(result["lineage_df"])
end

@testset "BarcodeModels integration (ResDmg ABM EvBC)" begin
    params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.01,
        zet_R = 0.01,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5)
    model = ResDmg_ABM_EvBC(params; abm = abm)
    @test model isa ResDmg_ABM_EvBC

    result = simulate_abm_experiment(model)
    @test result !== nothing
    @test haskey(result, "lineage_df")
    @test all(in(names(result["lineage_df"])).(["id", "parent_id", "birth_time"]))
    @test all(in(names(result["lineage_df"])).(["parent_pheno", "child_pheno"]))
    @test "barcode" in names(result["lineage_df"])
    @test "alive_at_end" in names(result["lineage_df"])
end

@testset "Lineage utilities" begin
    lineage_df = DataFrame(
        id = Int64[1, 2, 3, 4, 5],
        parent_id = Int64[0, 1, 1, 3, 1],
        birth_time = Float64[0.0, 1.0, 1.2, 2.0, 2.4],
        parent_pheno = ["ROOT", "S", "S", "R", "S"],
        child_pheno = ["S", "S", "R", "R", "S"],
        barcode = Float64[10.0, 10.0, 10.0, 10.0, 10.0],
        alive_at_end = Bool[false, true, false, true, false],
        rep = Int64[1, 1, 1, 1, 1]
    )

    edges = build_phylogeny(lineage_df)
    @test length(edges) == 4
    @test (1, 2) in edges
    @test (1, 3) in edges
    @test (3, 4) in edges
    @test (1, 5) in edges

    extant_edges = build_phylogeny(lineage_df; extant_only = true)
    @test length(extant_edges) == 3
    @test (1, 2) in extant_edges
    @test (1, 3) in extant_edges
    @test (3, 4) in extant_edges
    @test !((1, 5) in extant_edges)

    children = build_tree(lineage_df)
    @test haskey(children, 1)
    @test children[1] == Int64[2, 3, 5]
    @test children[3] == Int64[4]

    extant_children = build_tree(lineage_df; extant_only = true)
    @test haskey(extant_children, 1)
    @test extant_children[1] == Int64[2, 3]
    @test extant_children[3] == Int64[4]
    @test !(5 in extant_children[1])

    newick = lineage_to_newick(lineage_df, 1)
    @test newick == "(2,(4)3,5)1;"
    @test population_to_newick(lineage_df, 1) == newick

    extant_newick = lineage_to_newick(lineage_df, 1; extant_only = true)
    @test extant_newick == "(2,(4)3)1;"
    @test population_to_newick(lineage_df, 1; extant_only = true) == extant_newick

    meta_df = lineage_node_metadata(lineage_df)
    @test all(in(names(meta_df)).(["id", "parent_id", "birth_time", "parent_pheno", "child_pheno", "barcode", "alive_at_end", "rep"]))
    @test nrow(meta_df) == nrow(lineage_df)

    edge_bc_df = lineage_edge_barcodes(lineage_df)
    @test nrow(edge_bc_df) == 4
    @test all(in(names(edge_bc_df)).(["parent_id", "id", "parent_barcode", "child_barcode"]))
    @test all(edge_bc_df.parent_barcode .== edge_bc_df.child_barcode)

    lineage_df_no_alive = select(lineage_df, Not(:alive_at_end))
    @test_throws ArgumentError build_tree(lineage_df_no_alive; extant_only = true)
    @test_throws ArgumentError build_phylogeny(lineage_df_no_alive; extant_only = true)
    @test_throws ArgumentError lineage_to_newick(lineage_df_no_alive, 1; extant_only = true)
end

@testset "simulate_simple API" begin
    respop_params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )
    resdmg_params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.01,
        zet_R = 0.01,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    respop_hybrid = ResPop(respop_params)
    resdmg_hybrid = ResDmg(resdmg_params)
    respop_abm = ResPop_ABM(respop_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))
    resdmg_abm = ResDmg_ABM(resdmg_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))
    respop_evbc = ResPop_ABM_EvBC(respop_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))
    resdmg_evbc = ResDmg_ABM_EvBC(resdmg_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))

    hybrid_out = simulate_simple_run(respop_hybrid)
    @test haskey(hybrid_out, "sol_df")
    @test !haskey(hybrid_out, "t")
    @test !haskey(hybrid_out, "u")

    hybrid_dmg_out = simulate_simple_run(resdmg_hybrid)
    @test haskey(hybrid_dmg_out, "sol_df")
    @test !haskey(hybrid_dmg_out, "t")
    @test !haskey(hybrid_dmg_out, "u")

    abm_out = simulate_simple_run(respop_abm)
    @test haskey(abm_out, "sol_df")
    @test haskey(abm_out, "lin_df")
    @test !haskey(abm_out, "t")
    @test !haskey(abm_out, "u")

    abm_dmg_out = simulate_simple_run(resdmg_abm)
    @test haskey(abm_dmg_out, "sol_df")
    @test haskey(abm_dmg_out, "lin_df")
    @test !haskey(abm_dmg_out, "t")
    @test !haskey(abm_dmg_out, "u")

    evbc_out = simulate_simple_run(respop_evbc)
    @test haskey(evbc_out, "sol_df")
    @test haskey(evbc_out, "lin_df")
    @test haskey(evbc_out, "lineage_df")
    @test all(in(names(evbc_out["lineage_df"])).(["id", "parent_id", "birth_time", "parent_pheno", "child_pheno", "barcode", "alive_at_end"]))
    @test any(evbc_out["lineage_df"].alive_at_end)
    @test !haskey(evbc_out, "t")
    @test !haskey(evbc_out, "u")

    evbc_dmg_out = simulate_simple_run(resdmg_evbc)
    @test haskey(evbc_dmg_out, "sol_df")
    @test haskey(evbc_dmg_out, "lin_df")
    @test haskey(evbc_dmg_out, "lineage_df")
    @test all(in(names(evbc_dmg_out["lineage_df"])).(["id", "parent_id", "birth_time", "parent_pheno", "child_pheno", "barcode", "alive_at_end"]))
    @test any(evbc_dmg_out["lineage_df"].alive_at_end)
    @test !haskey(evbc_dmg_out, "t")
    @test !haskey(evbc_dmg_out, "u")

    treated_out = simulate_simple_run_treated(respop_abm)
    @test haskey(treated_out, "sol_df")
    @test haskey(treated_out, "lin_df")

    respop_abm_full_pheno = ResPop_ABM(respop_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5, full_pheno_bc = true))
    full_pheno_out = simulate_simple_run(respop_abm_full_pheno)
    @test haskey(full_pheno_out, "pheno_bc_df")
    @test !any(contains("EG"), string.(names(full_pheno_out["pheno_bc_df"])))
end

@testset "ResPopInVivo implementation" begin
    @test_throws ErrorException ResPopInVivoParams(
        b = 1.0, d = 0.1, rho = 0.1, mu = 0.01, sig = 0.01, del = 0.0, al = 0.0,
        Dc = 0.0, k = 0.0, psi = 0.0, drug_effect = :d,
        fEG1 = 1.1, pEG = 0.7, sEG = 0.5
    ) |> ResPopInVivo

    params = ResPopInVivoParams(
        b = 1.0,
        d = 0.1,
        rho = 0.2,
        mu = 0.01,
        sig = 0.01,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d,
        fEG1 = 0.4,
        pEG = 0.7,
        sEG = 0.5
    )

    invivo_sol_cols = [
        "nS_EG0", "nS_EG1", "nR_EG0", "nR_EG1", "nE_EG0", "nE_EG1",
        "nS", "nR", "nE", "n_EG0", "n_EG1", "N"
    ]
    function test_invivo_sol_df_counts(df)
        @test all(in(names(df)).(invivo_sol_cols))
        @test all(isapprox.(df.nS, df.nS_EG0 .+ df.nS_EG1; atol = 1e-8))
        @test all(isapprox.(df.nR, df.nR_EG0 .+ df.nR_EG1; atol = 1e-8))
        @test all(isapprox.(df.nE, df.nE_EG0 .+ df.nE_EG1; atol = 1e-8))
        @test all(isapprox.(df.n_EG0, df.nS_EG0 .+ df.nR_EG0 .+ df.nE_EG0; atol = 1e-8))
        @test all(isapprox.(df.n_EG1, df.nS_EG1 .+ df.nR_EG1 .+ df.nE_EG1; atol = 1e-8))
        @test all(isapprox.(df.N, df.n_EG0 .+ df.n_EG1; atol = 1e-8))
    end

    eg0_new, eg1_new, stats = engraftment_selection([10, 0, 0], [0, 10, 0], 1.0, 1.0)
    @test eg0_new == [0, 0, 0]
    @test eg1_new == [0, 10, 0]
    @test stats["N_engraft"] == 10

    hybrid = ResPopInVivo(params)
    hybrid_simple = simulate_simple(
        hybrid,
        SimpleSimParams(
            n0 = 20,
            tmax = 2.0,
            Nmax = 100,
            Cc = 100,
            treat_ons = Float64[],
            treat_offs = Float64[]
        )
    )
    @test haskey(hybrid_simple, "sol_df")
    test_invivo_sol_df_counts(hybrid_simple["sol_df"])

    eg_params = ResPopInVivoParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.2,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d,
        fEG1 = 0.5,
        pEG = 1.0,
        sEG = 0.0
    )
    eg_hybrid = ResPopInVivo(eg_params)
    eg_state = ResPopInVivoState(10, 0, 0, 10, 0, 0)
    eg_core_sim = BarcodeModels.SimParams(
        n0 = 20,
        t0 = 0.0,
        tmax = 0.2,
        t_Pass = -1.0,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        Nswitch = 100
    )
    eg_core = BarcodeModels.run_model_core_hybrid(eg_hybrid, eg_state, eg_core_sim; treat = false)
    @test all(u -> u[BarcodeModels.RESPOP_INVIVO_NS_EG1_INDEX] == 0.0, eg_core.u)
    @test all(u -> u[BarcodeModels.RESPOP_INVIVO_NE_EG0_INDEX] == 0.0, eg_core.u)
    @test all(u -> u[BarcodeModels.RESPOP_INVIVO_NE_EG1_INDEX] == 0.0, eg_core.u)

    ode_switch_sim = BarcodeModels.SimParams(
        n0 = 40,
        t0 = 0.0,
        tmax = 0.2,
        t_Pass = -1.0,
        Nmax = 1000,
        Cc = 1000,
        treat_ons = Float64[],
        treat_offs = Float64[],
        Nswitch = 1,
        save_at = 0.05
    )
    ode_switch_state = ResPopInVivoState(20, 20, 0, 0, 0, 0)
    ode_switch_sol = BarcodeModels.run_model_core_hybrid(hybrid, ode_switch_state, ode_switch_sim; treat = false)
    @test any(u -> abs(u[BarcodeModels.RESPOP_INVIVO_NS_EG0_INDEX] - round(u[BarcodeModels.RESPOP_INVIVO_NS_EG0_INDEX])) > 1e-6, ode_switch_sol.u)

    abm = ResPopInVivo_ABM(params; abm = ABMParams(Nbuff = 300, t_frac = 0.2, dt_save_at = 0.2))
    abm_simple = simulate_simple(
        abm,
        SimpleSimParams(
            n0 = 20,
            tmax = 2.0,
            Nmax = 100,
            Cc = 100,
            treat_ons = Float64[],
            treat_offs = Float64[]
        )
    )
    @test haskey(abm_simple, "sol_df")
    @test haskey(abm_simple, "lin_df")
    test_invivo_sol_df_counts(abm_simple["sol_df"])

    exp = ExperimentParams(
        n0 = 30,
        t_exp = 1.0,
        tmax = 3.0,
        t_Pass = [1.5],
        Nseed = 10,
        Nmax = 200,
        Cc = 200,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [3.0],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2
    )

    hybrid_exp = simulate_experiment(hybrid, exp)
    @test haskey(hybrid_exp, "engraft_df")
    @test haskey(hybrid_exp, "sol_df")
    test_invivo_sol_df_counts(hybrid_exp["sol_df"])

    abm_exp = simulate_experiment(abm, exp)
    @test haskey(abm_exp, "engraft_df")
    @test haskey(abm_exp, "sol_df")
    @test haskey(abm_exp, "lin_df")
    test_invivo_sol_df_counts(abm_exp["sol_df"])
    @test "passage" in names(abm_exp["sol_df"])

    pot_exp = ExperimentParams(
        n0 = 40,
        t_exp = 0.5,
        tmax = 1.0,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 200,
        Cc = 200,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = Float64[],
        Nswitch = 100,
        full_sol = true,
        n_rep = 1,
        inc_pot = true
    )

    hybrid_pot_exp = simulate_experiment(hybrid, pot_exp)
    @test haskey(hybrid_pot_exp, "sol_df")
    @test first(hybrid_pot_exp["sol_df"].cond) == "POT"
    hybrid_pot_rows = hybrid_pot_exp["sol_df"][hybrid_pot_exp["sol_df"].cond .== "POT", :]
    @test !isempty(hybrid_pot_rows)
    @test all(hybrid_pot_rows.rep .== 0)
    @test all(hybrid_pot_rows.passage .== 0)
    test_invivo_sol_df_counts(hybrid_pot_exp["sol_df"])

    abm_pot_exp = simulate_experiment(abm, pot_exp)
    @test haskey(abm_pot_exp, "sol_df")
    @test first(abm_pot_exp["sol_df"].cond) == "POT"
    abm_pot_rows = abm_pot_exp["sol_df"][abm_pot_exp["sol_df"].cond .== "POT", :]
    @test !isempty(abm_pot_rows)
    @test all(abm_pot_rows.rep .== 0)
    @test all(abm_pot_rows.passage .== 0)
    @test string(names(abm_pot_exp["lin_df"])[2]) == "POT_P0"
    test_invivo_sol_df_counts(abm_pot_exp["sol_df"])

    early_stop_params = ResPopInVivoParams(
        b = 1.0,
        d = 0.0,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d,
        fEG1 = 0.0,
        pEG = 1.0,
        sEG = 0.0
    )
    early_stop_exp = ExperimentParams(
        n0 = 50,
        t_exp = 0.1,
        tmax = [5.0, 6.0],
        t_Pass = Float64[],
        Nseed = 10,
        Nmax = 20,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = Float64[],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2
    )
    early_stop_hybrid = simulate_experiment(ResPopInVivo(early_stop_params), early_stop_exp)
    @test all(early_stop_hybrid["t"] .< early_stop_exp.tmax)
    @test all(isapprox.(
        early_stop_hybrid["t"],
        [maximum(early_stop_hybrid["sol_df"][early_stop_hybrid["sol_df"].rep .== i, :t]) for i in 1:2];
        atol = 1e-8
    ))

    control_exp = ExperimentParams(
        n0 = 80,
        t_exp = 0.2,
        tmax = 0.4,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 200,
        Cc = 200,
        treat_ons = Float64[0.1],
        treat_offs = Float64[0.4],
        t_keep = Float64[],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2,
        drug_treatment = true,
        inc_control = true
    )

    hybrid_control_exp = simulate_experiment(hybrid, control_exp)
    @test hybrid_control_exp["cond"] == ["CO", "CO", "DT", "DT"]
    @test hybrid_control_exp["rep"] == [1, 2, 1, 2]
    @test all(in(names(hybrid_control_exp["sol_df"])).(["cond", "rep"]))
    @test all(in(names(hybrid_control_exp["engraft_df"])).(["cond", "rep"]))
    @test unique(hybrid_control_exp["engraft_df"].cond) == ["CO", "DT"]

    abm_control_exp = simulate_experiment(abm, control_exp)
    @test abm_control_exp["cond"] == ["CO", "CO", "DT", "DT"]
    @test abm_control_exp["rep"] == [1, 2, 1, 2]
    @test all(in(names(abm_control_exp["sol_df"])).(["cond", "rep"]))
    @test all(in(names(abm_control_exp["engraft_df"])).(["cond", "rep"]))
    @test unique(abm_control_exp["engraft_df"].cond) == ["CO", "DT"]
    abm_lin_names = string.(names(abm_control_exp["lin_df"]))
    @test any(name -> startswith(name, "CO"), abm_lin_names)
    @test any(name -> startswith(name, "DT"), abm_lin_names)

    control_only_exp = ExperimentParams(
        n0 = 40,
        t_exp = 0.2,
        tmax = 0.4,
        t_Pass = Float64[],
        Nseed = 5,
        Nmax = 200,
        Cc = 200,
        treat_ons = Float64[0.1],
        treat_offs = Float64[0.4],
        t_keep = Float64[],
        Nswitch = 100,
        full_sol = true,
        n_rep = 2,
        drug_treatment = false
    )
    hybrid_control_only = simulate_experiment(hybrid, control_only_exp)
    @test hybrid_control_only["cond"] == ["CO", "CO"]
    @test hybrid_control_only["rep"] == [1, 2]
    @test unique(hybrid_control_only["engraft_df"].cond) == ["CO"]

    abm_control_only = simulate_experiment(abm, control_only_exp)
    @test abm_control_only["cond"] == ["CO", "CO"]
    @test abm_control_only["rep"] == [1, 2]
    @test unique(abm_control_only["engraft_df"].cond) == ["CO"]
    @test all(name -> startswith(name, "CO"), string.(names(abm_control_only["lin_df"]))[2:end])

    no_intermediate_engraft_params = ResPopInVivoParams(
        b = 1.0,
        d = 0.0,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d,
        fEG1 = 0.0,
        pEG = 0.0,
        sEG = 0.0
    )
    no_intermediate_engraft_exp = ExperimentParams(
        n0 = 30,
        t_exp = [0.1, 0.1],
        tmax = 0.2,
        t_Pass = Float64[],
        Nseed = [10, 10],
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [0.2],
        Nswitch = 100,
        full_sol = true,
        n_rep = 1
    )
    no_intermediate_hybrid = simulate_experiment(ResPopInVivo(no_intermediate_engraft_params), no_intermediate_engraft_exp)
    @test no_intermediate_hybrid["t"] != [-1.0]
    @test no_intermediate_hybrid["engraft_df"].N_engraft == [0]

    no_intermediate_abm = ResPopInVivo_ABM(no_intermediate_engraft_params; abm = ABMParams(Nbuff = 100, t_frac = 0.5, dt_save_at = 0.1))
    _, no_intermediate_abm_engraft = BarcodeModels._expand_split_cells_abm(no_intermediate_abm, no_intermediate_engraft_exp, 1)
    @test no_intermediate_abm_engraft.N_engraft == [0]

    # Backwards compatibility sanity check on original model family.
    compat_out = simulate_hybrid_experiment(ResPop(ResPopParams(
        b = 1.0, d = 0.1, rho = 0.0, mu = 0.0, sig = 0.0, del = 0.0, al = 0.0,
        Dc = 0.0, k = 0.0, psi = 0.0, drug_effect = :d
    )))
    @test haskey(compat_out, "t")
    @test haskey(compat_out, "u")
end

@testset "Vector tmax support (single-passage only)" begin
    respop_params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )
    resdmg_params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.01,
        zet_R = 0.01,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    respop_hybrid = ResPop(respop_params)
    resdmg_hybrid = ResDmg(resdmg_params)
    respop_abm = ResPop_ABM(respop_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))
    resdmg_abm = ResDmg_ABM(resdmg_params; abm = ABMParams(Nbuff = 200, t_frac = 0.2, dt_save_at = 0.5))

    result_respop_hybrid = simulate_hybrid_experiment_vector_tmax(respop_hybrid)
    result_resdmg_hybrid = simulate_hybrid_experiment_vector_tmax(resdmg_hybrid)
    result_respop_abm = simulate_abm_experiment_vector_tmax(respop_abm)
    result_resdmg_abm = simulate_abm_experiment_vector_tmax(resdmg_abm)
    abm_tmax_tol = 0.5

    @test maximum(result_respop_hybrid["sol_df"][result_respop_hybrid["sol_df"].rep .== 1, :t]) <= 2.0
    @test maximum(result_respop_hybrid["sol_df"][result_respop_hybrid["sol_df"].rep .== 2, :t]) <= 4.0
    @test maximum(result_resdmg_hybrid["sol_df"][result_resdmg_hybrid["sol_df"].rep .== 1, :t]) <= 2.0
    @test maximum(result_resdmg_hybrid["sol_df"][result_resdmg_hybrid["sol_df"].rep .== 2, :t]) <= 4.0
    @test maximum(result_respop_abm["sol_df"][result_respop_abm["sol_df"].rep .== 1, :t]) <= 1.0 + abm_tmax_tol
    @test maximum(result_respop_abm["sol_df"][result_respop_abm["sol_df"].rep .== 2, :t]) <= 2.0 + abm_tmax_tol
    @test maximum(result_resdmg_abm["sol_df"][result_resdmg_abm["sol_df"].rep .== 1, :t]) <= 1.0 + abm_tmax_tol
    @test maximum(result_resdmg_abm["sol_df"][result_resdmg_abm["sol_df"].rep .== 2, :t]) <= 2.0 + abm_tmax_tol

    @test_throws ErrorException BarcodeModels.simulate_experiment(
        respop_hybrid,
        ExperimentParams(
            n0 = 10,
            t_exp = 4.0,
            tmax = [2.0, 4.0],
            t_Pass = Float64[],
            Nseed = 10,
            Nmax = 100,
            Cc = 100,
            treat_ons = Float64[],
            treat_offs = Float64[],
            t_keep = [2.0],
            Nswitch = 100,
            n_rep = 2
        );
        n_rep = 3
    )

    @test_throws ErrorException ExperimentParams(
        n0 = 10,
        t_exp = 4.0,
        tmax = [2.0, 4.0],
        t_Pass = [1.0],
        Nseed = 10,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        n_rep = 2
    )
end

@testset "Hybrid invalid-parameter handling for inference mode" begin
    exp_dummy = ExperimentParams(
        n0 = 10,
        t_exp = 2.0,
        tmax = 2.0,
        t_Pass = Float64[],
        Nseed = 10,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [2.0],
        Nswitch = 100,
        full_sol = false,
        n_rep = 2
    )

    respop_invalid_psi = ResPop(ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = -0.1,
        drug_effect = :d
    ))

    res_pop_dummy = BarcodeModels.simulate_experiment(respop_invalid_psi, exp_dummy)
    @test haskey(res_pop_dummy, "t")
    @test haskey(res_pop_dummy, "u")
    @test all(<(0.0), res_pop_dummy["t"])
    @test all(<(0.0), res_pop_dummy["u"])

    @test_throws ErrorException BarcodeModels.simulate_experiment(respop_invalid_psi, exp_dummy; full_sol = true)

    resdmg_invalid = ResDmg(ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.2,
        zet_R = 0.1,
        Dc = 0.0,
        k = 0.0,
        psi = -0.1,
        drug_effect = :d
    ))

    res_dmg_dummy = BarcodeModels.simulate_experiment(resdmg_invalid, exp_dummy)
    @test haskey(res_dmg_dummy, "t")
    @test haskey(res_dmg_dummy, "u")
    @test all(<(0.0), res_dmg_dummy["t"])
    @test all(<(0.0), res_dmg_dummy["u"])

    @test_throws ErrorException BarcodeModels.simulate_experiment(resdmg_invalid, exp_dummy; full_sol = true)

    resdmg_no_resistant = ResDmg(ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.2,
        zet_R = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    ))

    res_dmg_no_resistant_out = BarcodeModels.simulate_experiment(resdmg_no_resistant, exp_dummy)
    @test haskey(res_dmg_no_resistant_out, "t")
    @test haskey(res_dmg_no_resistant_out, "u")
    @test !all(<(0.0), res_dmg_no_resistant_out["t"])
    @test !all(<(0.0), res_dmg_no_resistant_out["u"])
end

@testset "ABM strict invalid-parameter validation" begin
    invalid_respop_params = ResPopParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = -0.1,
        drug_effect = :d
    )
    @test_throws ErrorException ResPop_ABM(invalid_respop_params)

    invalid_resdmg_params = ResDmgParams(
        b = 1.0,
        d = 0.1,
        rho = 0.0,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        ome = 0.01,
        zet_S = 0.2,
        zet_R = 0.1,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )
    @test_throws ErrorException ResDmg_ABM(invalid_resdmg_params)
end

