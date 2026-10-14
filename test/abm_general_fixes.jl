@testset "Positive baseline growth is required in both coupling modes" begin
    for coupling in (:birth_coupled, :independent_rate),
        (param_type, constructors) in ((ResPopParams, (ResPop, ResPop_ABM, ResPop_ABM_EvBC)),
                                      (ResDmgParams, (ResDmg, ResDmg_ABM, ResDmg_ABM_EvBC)))
        extras = param_type === ResPopParams ? (al=0.0,) : (ome=0.0, zet_S=0.0, zet_R=0.0)
        for (b,d) in ((0.5,1.0), (1.0,1.0))
            params = param_type(;b,d,mu=0.0,sig=0.0,del=0.0,Dc=0.0,k=0.0,psi=0.0,
                transition_coupling=coupling,extras...)
            for constructor in constructors
                caught = try
                    constructor(params)
                    nothing
                catch err
                    err
                end
                @test caught isa ErrorException
                @test occursin("Baseline population growth rate", sprint(showerror,caught))
            end
        end
        params = param_type(;b=1.0,d=0.5,mu=0.0,sig=0.0,del=0.0,Dc=0.0,k=0.0,psi=0.0,
            transition_coupling=coupling,extras...)
        @test all(constructor -> constructor(params).params === params, constructors)
    end
end

@testset "ABM buffer and rate safeguards are common to both modes" begin
    constructors = (ResPop_ABM, ResPop_ABM_EvBC, ResDmg_ABM, ResDmg_ABM_EvBC)
    for constructor in constructors
        ispop = constructor in (ResPop_ABM, ResPop_ABM_EvBC)
        base = ispop ? coupling_pop_params : coupling_dmg_params
        states = Any[]
        outputs = Any[]
        for coupling in (:birth_coupled, :independent_rate)
            params = base(b=1.0,d=0.2,mu=0.0,sig=0.0,Dc=3.0,k=100.0,
                drug_effect=:c,transition_coupling=coupling)
            model = constructor(params)
            # No spare slots: treatment suppresses birth, but deaths must still run.
            state = coupling_abm_state(model,20;rho=1.0)
            Random.seed!(20261014)
            sim = BarcodeModels.ABMSimParams(t0=0.1,tmax=0.5,Nmax=1000,Cc=1000,
                treat_ons=[0.0],treat_offs=Float64[],dt_save_at=0.01,t_frac=0.1,drug_effect=:c)
            out = BarcodeModels.run_model_core_abm(model,state,sim;treat=true)
            @test count(c -> c.alive,state.cells) < 20
            @test all(n -> n <= 20,out.Nvec)
            push!(states,state)
            push!(outputs,out)
            # Above capacity, births/deaths have nonnegative zero density scaling.
            frozen = coupling_abm_state(constructor(base(mu=0.0,sig=0.0,
                transition_coupling=coupling)),4;spare=1)
            reverse!(frozen.cells)
            barcodes = sort([c.barcode for c in frozen.cells if c.alive])
            ids = hasproperty(frozen,:lineage_records) ? sort([c.id for c in frozen.cells if c.alive]) : nothing
            model = constructor(base(mu=0.0,sig=0.0,transition_coupling=coupling))
            out = BarcodeModels.run_model_core_abm(model,frozen,coupling_abm_sim(tmax=0.1,Cc=2))
            @test all(==(4),out.Nvec)
            @test sort([c.barcode for c in frozen.cells if c.alive]) == barcodes
            if ids !== nothing
                @test sort([c.id for c in frozen.cells if c.alive]) == ids
            end
        end
        @test outputs[1].Nvec == outputs[2].Nvec
        @test outputs[1].tvec == outputs[2].tvec
        @test [c.alive for c in states[1].cells] == [c.alive for c in states[2].cells]
        # Exercise the guarded ordering helper directly, including degenerate buffers.
        model = constructor(base(mu=0.0,sig=0.0))
        factory = ispop ? (constructor === ResPop_ABM ? BarcodeModels.make_dead_cell : BarcodeModels.make_dead_cell_evbc) :
            (constructor === ResDmg_ABM ? BarcodeModels.make_dead_resdmg_cell : BarcodeModels.make_dead_resdmg_cell_evbc)
        for live in (Bool[], [true], [false], [true,true], [false,false], [false,true], [true,false,true,false])
            state = coupling_abm_state(model,length(live))
            for (cell,alive) in zip(state.cells,live)
                cell.alive = alive
            end
            identities = sort([hasproperty(c,:id) ? (c.barcode,c.id,c.parent_id,c.birth_time) : (c.barcode,)
                for c in state.cells if c.alive])
            BarcodeModels.order_cells!(state.cells,factory)
            @test count(c -> c.alive,state.cells) == count(identity,live)
            @test sort([hasproperty(c,:id) ? (c.barcode,c.id,c.parent_id,c.birth_time) : (c.barcode,)
                for c in state.cells if c.alive]) == identities
            @test issorted([!c.alive for c in state.cells])
        end
    end
    for coupling in (:birth_coupled, :independent_rate), constructor in (ResDmg_ABM, ResDmg_ABM_EvBC)
        model = constructor(coupling_dmg_params(b=1.0,d=0.2,mu=0.0,sig=0.0,Dc=3.0,
            transition_coupling=coupling))
        state = coupling_abm_state(model,20)
        for cell in state.cells
            cell.DS = true
        end
        Random.seed!(20261014)
        BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.5,Cc=1000))
        @test count(c -> c.alive,state.cells) < 20 # d+Dc applies without treatment.
    end
    for coupling in (:birth_coupled, :independent_rate), constructor in (ResDmg_ABM, ResDmg_ABM_EvBC)
        model = constructor(coupling_dmg_params(mu=0.0,sig=0.0,ome=100.0,k=10000.0,
            transition_coupling=coupling))
        state = coupling_abm_state(model,1)
        Random.seed!(20261014)
        sim = BarcodeModels.ABMSimParams(tmax=0.2,Nmax=1000,Cc=1,
            treat_ons=[0.1],treat_offs=Float64[],dt_save_at=0.001,t_frac=0.1)
        out = BarcodeModels.run_model_core_abm(model,state,sim;treat=true)
        @test state.cells[1].DS # Treatment begins after the initial proposal time.
        @test all(==(1),out.Nvec)
    end
    # Defensive zero-envelope handling in the numerical core; public zero-growth
    # model construction is deliberately rejected by the new baseline check.
    for coupling in (:birth_coupled, :independent_rate)
        cells = BarcodeModels.seed_cells(4,0.0,4)
        params = coupling_pop_params(b=0.0,d=0.0,mu=0.0,sig=0.0,transition_coupling=coupling)
        out = BarcodeModels._core_grow_kill_abm!(cells,params,coupling_abm_sim(tmax=0.1))
        @test out.Nvec == [4,4] && out.fin_t == 0.1
        cells = BarcodeModels.seed_resdmg_cells(4,0.0,4)
        params = coupling_dmg_params(b=0.0,d=0.0,mu=0.0,sig=0.0,transition_coupling=coupling)
        out = BarcodeModels._core_grow_kill_abm!(cells,params,coupling_abm_sim(tmax=0.1))
        @test out.Nvec == [4,4] && out.fin_t == 0.1
    end
end
