using Test
using BarcodeModels

function simulate(model::AbstractBarcodeModel)
    exp = ExperimentParams(
        n0 = 10,
        rho = 0.0,
        t_exp = 6.0,
        tmax = 10.0,
        t_Pass = -1.0,
        Nseed = 10,
        Nmax = 100,
        Cc = 100,
        treat_ons = Float64[],
        treat_offs = Float64[],
        t_keep = [10.0],
        Nswitch = 100
    )
    return BarcodeModels.simulate(model, exp)
end

@testset "BarcodeModels integration" begin
    params = ModelParams(
        b = 1.0,
        d = 0.1,
        mu = 0.0,
        sig = 0.0,
        del = 0.0,
        al = 0.0,
        Dc = 0.0,
        k = 0.0,
        psi = 0.0,
        drug_effect = :d
    )

    model = MonoModel(params)
    @test model isa MonoModel

    result = simulate(model)
    @test result !== nothing
end
