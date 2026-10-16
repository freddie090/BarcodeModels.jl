@testset "EvBC id allocation preserves ancestry through cell movement" begin
    for (seeder, state_type, dead_cell) in (
        (BarcodeModels.seed_cells, BarcodeModels.ResPopABMEvBCState, BarcodeModels.make_dead_cell_evbc),
        (BarcodeModels.seed_resdmg_cells, BarcodeModels.ResDmgABMEvBCState, BarcodeModels.make_dead_resdmg_cell_evbc))
        base_cells = seeder(2, 0.0, 4)
        cells, next_id, records = BarcodeModels._to_evbc_cells(base_cells; start_id = 41)
        @test [c.id for c in cells if c.alive] == [41, 42]
        @test next_id == 43
        @test [c.barcode for c in cells] == [c.barcode for c in base_cells]
        roots = copy(records)
        state = state_type(cells, next_id, records)
        if state_type === BarcodeModels.ResPopABMEvBCState
            counts = BarcodeModels.PhenoCounts(2, 0, 0)
            BarcodeModels._birth_mutate_event_evbc!(state, 2, 3, 0.0, 0.0, 0.0, counts, 0.5)
        else
            counts = BarcodeModels.ResDmgPhenoCounts(2, 0, 0, 0)
            BarcodeModels.resdmg_birth_mutate_event_evbc!(state, 2, 3, 0.0, 0.0, counts, 0.5)
        end
        @test state.next_cell_id == 45
        @test cells[2].id == 43 && cells[2].parent_id == 42
        @test cells[3].id == 44 && cells[3].parent_id == 42

        # A bottleneck or transfer retains the cell and its existing ancestry.
        retained = cells[[2]]
        identity_before = (retained[1].id, retained[1].parent_id,
            retained[1].birth_time, retained[1].barcode)
        BarcodeModels.extend_with_dead_cells!(retained, 3, dead_cell)
        @test (retained[1].id, retained[1].parent_id,
            retained[1].birth_time, retained[1].barcode) == identity_before
        moved_state = state_type(retained, state.next_cell_id, records)
        if state_type === BarcodeModels.ResPopABMEvBCState
            counts = BarcodeModels.PhenoCounts(1, 0, 0)
            BarcodeModels._birth_mutate_event_evbc!(moved_state, 1, 2, 0.0, 0.0, 0.0, counts, 1.0)
        else
            counts = BarcodeModels.ResDmgPhenoCounts(1, 0, 0, 0)
            BarcodeModels.resdmg_birth_mutate_event_evbc!(moved_state, 1, 2, 0.0, 0.0, counts, 1.0)
        end
        @test records[1:2] == roots
        @test [r.parent_id for r in records[3:end]] == [42, 42, 43, 43]
        @test all(r -> r.barcode == 2.0, records[3:end])
        @test moved_state.next_cell_id == 47

        later_cells, later_next, later_records = BarcodeModels._to_evbc_cells(base_cells;
            start_id = moved_state.next_cell_id)
        @test [r.id for r in later_records] == [47, 48]
        @test later_next == 49
        @test allunique([r.id for r in vcat(records, later_records)])
        default_cells, default_next, _ = BarcodeModels._to_evbc_cells(base_cells)
        @test [c.id for c in default_cells if c.alive] == [1, 2]
        @test default_next == 3
        for empty_cells in (empty(later_cells), [dead_cell()])
            empty_next, empty_records = BarcodeModels._init_lineage_state!(empty_cells; start_id = 49)
            @test empty_next == 49 && isempty(empty_records)
        end
        @test_throws ArgumentError BarcodeModels._to_evbc_cells(base_cells; start_id = 0)
    end
end

@testset "Lineage utilities reject ambiguous ids and allow cross-replicate ancestry" begin
    lineage = DataFrame(id = Int64[10, 11], parent_id = Int64[0, 10],
        birth_time = [0.0, 1.0], barcode = [7.0, 7.0],
        alive_at_end = [false, true], rep = Int64[1, 2])
    @test build_phylogeny(lineage) == [(10, 11)]
    @test build_tree(lineage; extant_only = true) == Dict(10 => [11])
    @test lineage_to_newick(lineage, 10; extant_only = true) == "(11)10;"
    edges = lineage_edge_barcodes(lineage)
    @test nrow(edges) == 1
    @test edges.parent_barcode == edges.child_barcode == [7.0]

    legacy = vcat(lineage, lineage)
    legacy.rep = [1, 1, 2, 2]
    for extant_only in (false, true)
        @test_throws ArgumentError build_phylogeny(legacy; extant_only = extant_only)
        @test_throws ArgumentError build_tree(legacy; extant_only = extant_only)
        @test_throws ArgumentError lineage_to_newick(legacy, 10; extant_only = extant_only)
        @test_throws ArgumentError population_to_newick(legacy, 10; extant_only = extant_only)
    end
    @test_throws ArgumentError lineage_edge_barcodes(legacy)
    @test build_tree(filter(:rep => ==(1), legacy)) == Dict(10 => [11])
    @test nrow(lineage_node_metadata(legacy)) == 4
    @test isempty(build_tree(lineage[Int[], :]))
    @test isempty(build_phylogeny(lineage[Int[], :]; extant_only = true))
    @test nrow(lineage_edge_barcodes(lineage[Int[], :])) == 0
end

@testset "EvBC experiment ids are unique across replicates and passages" begin
    abm = ABMParams(Nbuff = 128, t_frac = 0.2,
        use_lib_probs = true, Nbc = 1, bc_probs = [1.0])
    experiment = ExperimentParams(n0 = 8, t_exp = 0.0, tmax = 0.4, t_Pass = [0.2],
        Nseed = 4, Nmax = 100, Cc = 100, treat_ons = Float64[], treat_offs = Float64[],
        t_keep = [0.4], Nswitch = 100, n_rep = 2, drug_treatment = false)
    for coupling in (:birth_coupled, :independent_rate)
        pop = ResPopParams(b = 1.0, d = 0.0, rho = 0.25, mu = 0.0, sig = 0.0,
            del = 0.0, al = 0.0, Dc = 0.0, k = 0.0, psi = 0.0, transition_coupling = coupling)
        dmg = ResDmgParams(b = 1.0, d = 0.0, rho = 0.25, mu = 0.0, sig = 0.0,
            del = 0.0, ome = 0.0, zet_S = 0.0, zet_R = 0.0, Dc = 0.0, k = 0.0,
            psi = 0.0, transition_coupling = coupling)
        for (model, standard) in ((ResPop_ABM_EvBC(pop; abm), ResPop_ABM(pop; abm)),
            (ResDmg_ABM_EvBC(dmg; abm), ResDmg_ABM(dmg; abm)))
            Random.seed!(20261015)
            out = simulate_experiment(model, experiment)
            next_draw = rand()
            Random.seed!(20261015)
            standard_out = simulate_experiment(standard, experiment)
            @test rand() == next_draw
            for key in ("sol_df", "lin_df", "t", "u")
                @test out[key] == standard_out[key]
            end

            lineage = out["lineage_df"]
            reps = [filter(:rep => ==(i), lineage) for i in 1:2]
            @test allunique(lineage.id)
            @test all(==(1.0), lineage.barcode)
            @test sort(lineage.id) == collect(1:maximum(lineage.id))
            @test maximum(reps[1].id) < minimum(reps[2].id)
            @test all(rep -> count(==(0), rep.parent_id) == 4, reps)
            @test all(id -> id == 0 || id in lineage.id, lineage.parent_id)
            @test all(<=(2), length.(values(build_tree(lineage))))
            for extant_only in (false, true)
                expected = vcat(build_phylogeny(reps[1]; extant_only = extant_only),
                    build_phylogeny(reps[2]; extant_only = extant_only))
                @test sort(build_phylogeny(lineage; extant_only = extant_only)) == sort(expected)
                @test build_tree(lineage; extant_only = extant_only) ==
                    merge(build_tree(reps[1]; extant_only = extant_only),
                        build_tree(reps[2]; extant_only = extant_only))
            end
            edges = lineage_edge_barcodes(lineage)
            @test nrow(edges) == count(!=(0), lineage.parent_id)
            @test !any(ismissing, edges.parent_barcode)
            @test edges.parent_barcode == edges.child_barcode
            for rep in reps, root in rep.id[rep.parent_id .== 0]
                @test lineage_to_newick(lineage, root) == lineage_to_newick(rep, root)
            end
            for i in 1:2, passage in 1:2
                @test "DT$(i)_P$(passage)" in names(out["lin_df"])
            end
        end
    end
end
