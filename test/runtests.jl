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

