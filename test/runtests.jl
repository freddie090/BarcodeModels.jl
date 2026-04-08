using Test
using BarcodeModels

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
        n0 = 5,
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

    @test maximum(result_respop_hybrid["sol_df"][result_respop_hybrid["sol_df"].rep .== 1, :t]) <= 2.0
    @test maximum(result_respop_hybrid["sol_df"][result_respop_hybrid["sol_df"].rep .== 2, :t]) <= 4.0
    @test maximum(result_resdmg_hybrid["sol_df"][result_resdmg_hybrid["sol_df"].rep .== 1, :t]) <= 2.0
    @test maximum(result_resdmg_hybrid["sol_df"][result_resdmg_hybrid["sol_df"].rep .== 2, :t]) <= 4.0
    @test maximum(result_respop_abm["sol_df"][result_respop_abm["sol_df"].rep .== 1, :t]) <= 1.0
    @test maximum(result_respop_abm["sol_df"][result_respop_abm["sol_df"].rep .== 2, :t]) <= 2.0
    @test maximum(result_resdmg_abm["sol_df"][result_resdmg_abm["sol_df"].rep .== 1, :t]) <= 1.0
    @test maximum(result_resdmg_abm["sol_df"][result_resdmg_abm["sol_df"].rep .== 2, :t]) <= 2.0

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

