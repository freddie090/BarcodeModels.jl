@testset "ABM seeding retains barcode, phenotype and buffer contracts" begin
    for seeder in (BarcodeModels.seed_cells, BarcodeModels.seed_resdmg_cells)
        for settings in ((;),
            (; use_lib_probs=true, Nbc=3, bc_probs=[0.0, 1.0, 0.0]),
            (; skew_lib=true, Nbc=3, bc_unif=0.0),
            (; skew_lib=true, Nbc=3, bc_unif=0.6))
            Random.seed!(20261014)
            cells = seeder(12, 0.25, 16; settings...)
            @test length(cells) == 16
            @test all(c -> c.alive, cells[1:12])
            @test count(c -> c.R, cells) == 3
            @test all(c -> !c.alive && c.barcode == 0.0 && !c.R, cells[13:16])
            if isempty(settings)
                @test [c.barcode for c in cells[1:12]] == collect(1.0:12.0)
            elseif get(settings, :use_lib_probs, false)
                @test all(c -> c.barcode == 2.0, cells[1:12])
            elseif settings.bc_unif == 0.0
                @test all(c -> c.barcode == 1.0, cells[1:12])
            else
                @test all(c -> c.barcode in 1:3, cells[1:12])
            end
            if seeder === BarcodeModels.seed_cells
                @test all(c -> !c.E, cells)
            else
                @test all(c -> !c.DS && !c.DR, cells)
            end
        end
        @test length(seeder(4, 0.0, 2)) == 4
        @test isempty(seeder(0, 0.0, 0))
    end
end

@testset "EvBC root initialisation retains phenotype labels and birth ancestry" begin
    pop_cells = [BarcodeModels.make_dead_cell_evbc(),
        BarcodeModels.CancerCellEvBC(1.0, false, false, true, 99, 98, 97.0),
        BarcodeModels.CancerCellEvBC(2.0, true, false, true, 99, 98, 97.0),
        BarcodeModels.make_dead_cell_evbc(),
        BarcodeModels.CancerCellEvBC(3.0, false, true, true, 99, 98, 97.0)]
    dmg_cells = [BarcodeModels.make_dead_resdmg_cell_evbc(),
        BarcodeModels.ResDmgCellEvBC(1.0, false, false, false, true, 99, 98, 97.0),
        BarcodeModels.ResDmgCellEvBC(2.0, true, false, false, true, 99, 98, 97.0),
        BarcodeModels.make_dead_resdmg_cell_evbc(),
        BarcodeModels.ResDmgCellEvBC(3.0, false, true, false, true, 99, 98, 97.0),
        BarcodeModels.ResDmgCellEvBC(4.0, false, false, true, true, 99, 98, 97.0)]
    for (cells, labels) in ((pop_cells, ["S", "R", "E"]),
        (dmg_cells, ["S", "DS", "DR", "R"]))
        cells[1].id, cells[1].parent_id, cells[1].birth_time = 99, 98, 97.0
        n_live = length(labels)
        next_id, records = BarcodeModels._init_lineage_state!(cells, 2.5)
        @test next_id == n_live + 1
        @test [c.id for c in cells if c.alive] == collect(1:n_live)
        @test all(c -> c.parent_id == 0 && c.birth_time == 2.5, filter(c -> c.alive, cells))
        @test all(c -> c.id == 0 && c.parent_id == 0 && c.birth_time == -1.0,
            filter(c -> !c.alive, cells))
        @test [r.id for r in records] == collect(1:n_live)
        @test [r.child_pheno for r in records] == labels
        @test [r.barcode for r in records] == [c.barcode for c in cells if c.alive]
        @test all(r -> r.parent_id == 0 && r.birth_time == 2.5 && r.parent_pheno == "ROOT", records)
        roots = copy(records)
        if cells === pop_cells
            state = BarcodeModels.ResPopABMEvBCState(cells, next_id, records)
            counts = BarcodeModels.PhenoCounts(1, 1, 1)
            BarcodeModels._birth_mutate_event_evbc!(state, 2, 1, 0.0, 0.0, 0.0, counts, 3.0)
        else
            state = BarcodeModels.ResDmgABMEvBCState(cells, next_id, records)
            counts = BarcodeModels.ResDmgPhenoCounts(1, 1, 1, 1)
            BarcodeModels.resdmg_birth_mutate_event_evbc!(state, 2, 1, 0.0, 0.0, counts, 3.0)
        end
        @test records[1:n_live] == roots
        @test state.next_cell_id == n_live + 3
        @test cells[2].parent_id == cells[1].parent_id == 1
        @test cells[2].birth_time == cells[1].birth_time == 3.0
        @test [r.id for r in records[(n_live + 1):end]] == [n_live + 1, n_live + 2]
        @test all(r -> r.parent_id == 1 && r.parent_pheno == "S" && r.child_pheno == "S",
            records[(n_live + 1):end])

        alive_ids = Int64[c.id for c in cells if c.alive]
        df = BarcodeModels._lineage_df(records, 7, alive_ids)
        @test names(df) == ["id", "parent_id", "birth_time", "parent_pheno", "child_pheno",
            "barcode", "alive_at_end", "rep"]
        @test df.alive_at_end == [r.id in alive_ids for r in records]
        @test !df.alive_at_end[1]
        @test all(==(7), df.rep)

        dead_cells = cells[[4]]
        dead_next, dead_records = BarcodeModels._init_lineage_state!(dead_cells)
        @test dead_next == 1 && isempty(dead_records)
        empty_next, empty_records = BarcodeModels._init_lineage_state!(empty(cells))
        @test empty_next == 1 && isempty(empty_records)
    end
    empty_df = BarcodeModels._lineage_df(BarcodeModels.LineageRecord[], 7)
    @test size(empty_df) == (0, 8)
    @test eltype.(eachcol(empty_df)) == [Int64, Int64, Float64, String, String, Float64, Bool, Int64]
end

@testset "ABM public outputs retain passage and subsampling structure" begin
    abm = ABMParams(Nbuff=64, t_frac=0.2, dt_save_at=0.05, sub_sample_cells=true, K=3)
    pop = ResPopParams(b=1.0, d=0.0, mu=0.0, sig=0.0, del=0.0, al=0.0,
        Dc=0.0, k=0.0, psi=0.0)
    dmg = ResDmgParams(b=1.0, d=0.0, mu=0.0, sig=0.0, del=0.0, ome=0.0,
        zet_S=0.0, zet_R=0.0, Dc=0.0, k=0.0, psi=0.0)
    simple = SimpleSimParams(n0=8, tmax=0.3, Nmax=100, Cc=100,
        treat_ons=Float64[], treat_offs=Float64[])
    experiment = ExperimentParams(n0=8, t_exp=0.0, tmax=0.3, t_Pass=[0.15],
        Nseed=4, Nmax=100, Cc=100, treat_ons=Float64[], treat_offs=Float64[],
        t_keep=[0.3], Nswitch=100, n_rep=1, full_sol=true, drug_treatment=false)
    for model in (ResPop_ABM(pop; abm), ResPop_ABM_EvBC(pop; abm),
        ResDmg_ABM(dmg; abm), ResDmg_ABM_EvBC(dmg; abm))
        Random.seed!(20261014)
        simple_out = simulate_simple(model, simple)
        @test all(k -> haskey(simple_out, k), ("sol_df", "lin_df"))
        @test all(>=(8), simple_out["sol_df"].N)
        Random.seed!(20261014)
        out = simulate_experiment(model, experiment)
        @test all(k -> haskey(out, k), ("sol_df", "lin_df", "sub_lin_df", "t", "u"))
        @test "DT1_P1" in names(out["lin_df"])
        @test "DT1_P2" in names(out["lin_df"])
        for col in ("DT1_P1", "DT1_P2")
            @test sum(out["sub_lin_df"][!, col]) == 3
            @test sum(out["lin_df"][!, col]) >= 4
        end
        if model isa Union{ResPop_ABM_EvBC, ResDmg_ABM_EvBC}
            @test haskey(simple_out, "lineage_df") && haskey(out, "lineage_df")
            @test "alive_at_end" in names(out["lineage_df"])
        end
    end
end
