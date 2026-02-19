using Test
using BarcodeModels

function simulate_hybrid_experiment(model)
    exp = ExperimentParams(
        n0 = 10,
        t_exp = 8.0,
        tmax = 10.0,
        t_Pass = -1.0,
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
        t_Pass = -1.0,
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

@testset "BarcodeModels integration" begin
    params = ModelParams(
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
    params = ModelParams(
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
