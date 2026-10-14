using Random

# Conversion-only checks suppress division at capacity using a valid growing baseline.
function coupling_pop_params(; kwargs...)
    ResPopParams(; b=1.0, d=0.0, mu=0.7, sig=0.3, del=0.0,
        al=0.0, Dc=0.0, k=0.0, psi=0.0, transition_coupling=:independent_rate,
        kwargs...)
end

function coupling_dmg_params(; kwargs...)
    ResDmgParams(; b=1.0, d=0.0, mu=0.7, sig=0.3, del=0.0,
        ome=0.0, zet_S=0.0, zet_R=0.0, Dc=0.0, k=0.0, psi=0.0,
        transition_coupling=:independent_rate, kwargs...)
end

function coupling_abm_state(model, n=20; rho=0.0, spare=0)
    if model isa Union{ResPop_ABM,ResPop_ABM_EvBC}
        cells = BarcodeModels.seed_cells(Int64(n), Float64(rho), Int64(n+spare))
        if model isa ResPop_ABM_EvBC
            cells, next_id, records = BarcodeModels._to_evbc_cells(cells)
            return BarcodeModels.ResPopABMEvBCState(cells,next_id,records)
        end
        return BarcodeModels.ResPopABMState(cells)
    end
    cells = BarcodeModels.seed_resdmg_cells(Int64(n),Float64(rho),Int64(n+spare))
    if model isa ResDmg_ABM_EvBC
        cells,next_id,records = BarcodeModels._to_evbc_cells(cells)
        return BarcodeModels.ResDmgABMEvBCState(cells,next_id,records)
    end
    return BarcodeModels.ResDmgABMState(cells)
end

function coupling_abm_sim(; tmax=0.7, Cc=20, dt_save_at=0.01, treat=false, treat_on=0.02)
    BarcodeModels.ABMSimParams(; tmax,Nmax=1000,Cc,
        treat_ons=treat ? [treat_on] : Float64[], treat_offs=Float64[],
        dt_save_at,t_frac=0.1)
end

@testset "Transition coupling constructors, validation and copies" begin
    pop_kw = (b=1.0,d=0.1,mu=0.1,sig=0.1,del=0.0,al=0.1,Dc=0.0,k=0.0,psi=0.0)
    dmg_kw = (b=1.0,d=0.1,mu=0.1,sig=0.1,del=0.0,ome=0.0,zet_S=0.0,
        zet_R=0.0,Dc=0.0,k=0.0,psi=0.0)
    for (param_type,hybrid,kwargs,copy_params) in (
        (ResPopParams,ResPop,pop_kw,BarcodeModels._copy_respop_params),
        (ResDmgParams,ResDmg,dmg_kw,BarcodeModels._copy_resdmg_params))
        @test param_type(;kwargs...).transition_coupling === :birth_coupled
        for mode in (:birth_coupled,:independent_rate)
            for value in (mode,String(mode))
                p = param_type(;kwargs...,transition_coupling=value)
                @test p.transition_coupling === mode
                @test copy_params(p;rho=0.5,drug_effect=:c).transition_coupling === mode
            end
        end
        for mode in (:invalid,"INDEPENDENT_RATE","independent-rate",nothing,1)
            @test_throws ErrorException param_type(;kwargs...,transition_coupling=mode)
        end
        @test_throws ErrorException hybrid(;kwargs...,mu=1.1)
        @test_throws ErrorException hybrid(;kwargs...,sig=1.1)
        @test hybrid(;kwargs...,mu=2.0,sig=3.0,transition_coupling=:independent_rate).params.mu == 2.0
        for name in (:mu,:sig), value in (-0.01,NaN,Inf,-Inf)
            bad = NamedTuple{(name,)}((value,))
            @test_throws ErrorException hybrid(;kwargs...,bad...,transition_coupling=:independent_rate)
        end
    end
    @test_throws ErrorException ResPop(;pop_kw...,sig=0.8,al=0.3)
    @test ResPop(;pop_kw...,sig=2.0,al=3.0,transition_coupling=:independent_rate).params.al == 3.0
    for value in (-0.01,NaN,Inf,-Inf)
        @test_throws ErrorException ResPop(;pop_kw...,al=value,transition_coupling=:independent_rate)
    end
    @test ResPopParams(1.0,0.1,0.0,0.1,0.1,0.0,0.1,0.0,0.0,0.0,:d).transition_coupling === :birth_coupled
    @test ResDmgParams(1.0,0.1,0.0,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,:d).transition_coupling === :birth_coupled
    @test_throws ErrorException ResPop_ABM(coupling_pop_params(psi=-0.1))
    @test_throws ErrorException ResDmg_ABM(coupling_dmg_params(zet_S=0.4,zet_R=0.2))
    @test_throws ErrorException ResPop(coupling_pop_params(drug_effect=:b,Dc=1.1))
    for (model_type,kwargs) in ((ResPop,pop_kw),(ResDmg,dmg_kw))
        implicit = model_type(;kwargs...)
        explicit = model_type(;kwargs...,transition_coupling=:birth_coupled)
        state = model_type === ResPop ? ResPopState(20,10,0) : ResDmgState(20,0,0,10)
        a = hybrid_test_setup(implicit,state,hybrid_test_sim())
        b = hybrid_test_setup(explicit,state,hybrid_test_sim())
        @test a.flags == b.flags
        @test [f(a.u0,a.p,0.0) for f in a.activities] == [f(b.u0,b.p,0.0) for f in b.activities]
    end
end

@testset "Independent transition accounting for every representation mask" begin
    sim = hybrid_test_sim(Cc=1000)
    pop = ResPop(b=1.2,d=0.3,mu=1.3,sig=0.2,del=0.2,al=0.4,
        Dc=0.5,k=0.0,psi=0.6,transition_coupling=:independent_rate)
    dmg = ResDmg(b=1.2,d=0.3,mu=1.3,sig=0.2,del=0.2,
        ome=0.7,zet_S=0.2,zet_R=0.15,Dc=0.5,k=0.0,psi=0.6,
        transition_coupling=:independent_rate)
    # BS=144, BR=86.4, BE=36; deaths=(60,28.8,11.4); conversion=(260,30,24).
    pop_drift = [0.0,144-60-260+30,86.4-28.8+260-30-24,36-11.4+24,0.0]
    pop_changes = (((2,1),),((2,-1),(3,1)),((2,-1),),((3,1),),
        ((3,-1),(2,1)),((3,-1),(4,1)),((3,-1),),((4,1),),((4,-1),))
    # Damage=(56,16.8), repair=(6,3), damaged death=(14.4,9.6).
    dmg_drift = [0.0,144-60-260+30-56+6,56-6-14.4,
        16.8-3-9.6,86.4-28.8+260-30-16.8+3,0.0]
    dmg_changes = (((2,1),),((2,-1),(5,1)),((2,-1),),((5,1),),
        ((5,-1),(2,1)),((5,-1),),((3,-1),),((4,-1),),
        ((2,-1),(3,1)),((5,-1),(4,1)),((3,-1),(2,1)),((4,-1),(5,1)))
    for (model,state,drift,changes) in (
        (pop,ResPopState(200,150,50;gam=0.4),pop_drift,pop_changes),
        (dmg,ResDmgState(200,30,20,150;gam=0.4),dmg_drift,dmg_changes))
        setup = hybrid_test_setup(model,state,sim)
        u,p = setup.u0,setup.p
        @test length(setup.flags) == length(changes)
        @test :S_to_R_transition_ode in setup.flags
        @test :S_asymmetric_R_birth_ode ∉ setup.flags
        activities = [activity(u,p,0.0) for activity in setup.activities]
        for mask in 0:(2^length(setup.flags)-1)
            for (i,flag) in enumerate(setup.flags)
                setproperty!(p,flag,!iszero(mask & (1 << (i-1))))
            end
            du = similar(u)
            setup.ode!(du,u,p,0.0)
            combined = collect(du)
            for (i,activity) in enumerate(setup.activities)
                @test activity(u,p,0.0) == activities[i]
                jump_rate = setup.jump_rates[i](u,p,0.0)
                @test jump_rate == (getproperty(p,setup.flags[i]) == 1.0 ? 0.0 : activities[i])
                for (index,change) in changes[i]
                    combined[index] += change*jump_rate
                end
            end
            @test combined ≈ drift atol=1e-10 rtol=1e-12
        end
    end
end

@testset "Independent conversion thresholds, fraction retention and donor guards" begin
    for (model,state,Rindex) in (
        (ResPop(coupling_pop_params(mu=0.5,sig=2.0)),ResPopState(40,0,0),3),
        (ResDmg(coupling_dmg_params(mu=0.5,sig=2.0)),ResDmgState(40,0,0,0),5))
        setup = hybrid_test_setup(model,state,hybrid_test_sim(Nswitch=40,Aswitch=20.0,Cc=40))
        u,p = setup.u0,setup.p
        @test p.S_to_R_transition_ode == 1.0 # Inclusive bounds; recipient empty; at carrying capacity.
        @test setup.activities[1](u,p,0.0) == 0.0
        @test setup.activities[2](u,p,0.0) == 20.0
        u[Rindex] = 0.4
        setup.refresh!(u,p,0.0)
        @test u[Rindex] == 0.4
        @test setup.activities[5](u,p,0.0) == 0.8
        @test setup.jump_rates[5](u,p,0.0) == 0.0
        u[Rindex] = 1.4
        setup.refresh!(u,p,0.0)
        @test u[Rindex] == 1.4
        @test setup.jump_rates[5](u,p,0.0) ≈ 2.8
        p.mu = prevfloat(0.5)
        setup.classify!(u,p,0.0)
        @test p.S_to_R_transition_ode == 0.0
        p.mu = 0.0
        u[Rindex] = 0.4
        setup.refresh!(u,p,0.0)
        @test u[Rindex] in (0.0,1.0)
        snapshot = copy(u)
        for _ in 1:20
            setup.refresh!(u,p,0.0)
        end
        @test u == snapshot
        # A low-population S donor can remain fractional because ODE R->S also affects it.
        u.nS = 0.4
        u[Rindex] = 40.0
        p.mu = 2.0
        setup.refresh!(u,p,0.0)
        @test p.R_to_S_transition_ode == 1.0
        @test u.nS == 0.4
        @test setup.activities[2](u,p,0.0) == 0.8
        @test setup.jump_rates[2](u,p,0.0) == 0.0
        u.nS = 40.2
        setup.refresh!(u,p,0.0)
        @test u.nS == 40.2
    end
    pop = ResPop(coupling_pop_params(mu=0.0,sig=0.0,al=2.0))
    setup = hybrid_test_setup(pop,ResPopState(0,40,0;gam=0.25),
        hybrid_test_sim(Nswitch=40,Aswitch=20.0,Cc=40))
    @test setup.p.R_to_E_transition_ode == 1.0
    @test setup.activities[6](setup.u0,setup.p,0.0) == 20.0
    setup.u0.gam = 2.0
    @test setup.activities[6](setup.u0,setup.p,0.0) == 80.0
    setup.u0.gam = -1.0
    @test setup.activities[6](setup.u0,setup.p,0.0) == 0.0
    setup.classify!(setup.u0,setup.p,0.0)
    @test setup.p.R_to_E_transition_ode == 0.0
end

@testset "ABM conversions preserve cells, barcodes and EvBC ancestry" begin
    for model in (ResPop_ABM(coupling_pop_params(mu=10.0,sig=0.0)),
        ResPop_ABM_EvBC(coupling_pop_params(mu=10.0,sig=0.0)),
        ResDmg_ABM(coupling_dmg_params(mu=10.0,sig=0.0)),
        ResDmg_ABM_EvBC(coupling_dmg_params(mu=10.0,sig=0.0)))
        Random.seed!(182)
        state = coupling_abm_state(model,20)
        before_barcodes = [cell.barcode for cell in state.cells]
        isevbc = hasproperty(state,:lineage_records)
        if isevbc
            identities = [(c.id,c.parent_id,c.birth_time) for c in state.cells]
            records = copy(state.lineage_records)
            next_id = state.next_cell_id
        end
        out = BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim())
        @test all(==(20),out.Nvec)
        @test all(c->c.alive,state.cells)
        @test [c.barcode for c in state.cells] == before_barcodes
        @test count(c->c.R,state.cells) > 0
        @test sum(c->c.R,state.cells) + sum(c->!c.R,state.cells) == 20
        @test out.fin_t == 0.7
        if isevbc
            @test [(c.id,c.parent_id,c.birth_time) for c in state.cells] == identities
            @test state.lineage_records == records
            @test state.next_cell_id == next_id
        end
    end
    for (model_type,params) in (
        (ResPop_ABM,coupling_pop_params(mu=0.0,sig=0.0)),
        (ResPop_ABM_EvBC,coupling_pop_params(mu=0.0,sig=0.0)),
        (ResDmg_ABM,coupling_dmg_params(mu=0.0,sig=0.0)),
        (ResDmg_ABM_EvBC,coupling_dmg_params(mu=0.0,sig=0.0)))
        model = model_type(params)
        state = coupling_abm_state(model,4)
        out = BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.2,Cc=4))
        @test all(==(4),out.Nvec)
        @test out.fin_t == 0.2
        @test all(c->c.alive && !c.R,state.cells)
        # An unsorted spare slot must not overwrite the last living cell.
        state = coupling_abm_state(model,4;spare=1)
        reverse!(state.cells)
        @test !state.cells[1].alive && state.cells[end].alive
        barcodes = [c.barcode for c in state.cells if c.alive]
        if hasproperty(state,:lineage_records)
            ids = [c.id for c in state.cells if c.alive]
        end
        out = BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.2,Cc=4))
        @test all(==(4),out.Nvec)
        @test sort([c.barcode for c in state.cells if c.alive]) == sort(barcodes)
        if hasproperty(state,:lineage_records)
            @test sort([c.id for c in state.cells if c.alive]) == sort(ids)
        end
    end
    for model in (ResDmg_ABM(coupling_dmg_params(mu=100.0,sig=100.0)),
        ResDmg_ABM_EvBC(coupling_dmg_params(mu=100.0,sig=100.0)))
        state = coupling_abm_state(model,4)
        state.cells[1].DS = true
        state.cells[2].DR = true
        state.cells[2].R = false
        BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.1,Cc=4))
        @test state.cells[1].DS && !state.cells[1].R
        @test state.cells[2].DR && !state.cells[2].R
    end
    for (model_type,params) in (
        (ResPop_ABM,coupling_pop_params(b=100.0,mu=0.0,sig=0.0)),
        (ResPop_ABM_EvBC,coupling_pop_params(b=100.0,mu=0.0,sig=0.0)),
        (ResDmg_ABM,coupling_dmg_params(b=100.0,mu=0.0,sig=0.0)),
        (ResDmg_ABM_EvBC,coupling_dmg_params(b=100.0,mu=0.0,sig=0.0)))
        Random.seed!(55)
        model = model_type(params)
        state = coupling_abm_state(model,4)
        @test_throws ErrorException BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(Cc=1000))
    end

    # The birth records stay historical; subsequent division records current parent phenotype.
    pop = ResPop_ABM_EvBC(coupling_pop_params())
    state = coupling_abm_state(pop,1;spare=1)
    counts = BarcodeModels.PhenoCounts(1,0,0)
    record = state.lineage_records[1]
    BarcodeModels.respop_S_to_R_transition_event!(state.cells,1,counts)
    @test counts.Scount == 0 && counts.Rcount == 1
    @test state.cells[1].alive && state.cells[1].R && !state.cells[1].E
    @test state.lineage_records[1] == record
    BarcodeModels._birth_mutate_event_evbc!(state,1,2,0.0,0.0,0.0,counts,0.2)
    @test state.lineage_records[end].parent_pheno == "R"
    @test state.lineage_records[end].child_pheno == "R"
    @test state.cells[2].parent_id == state.cells[1].id
    BarcodeModels.respop_R_to_E_transition_event!(state.cells,1,counts)
    BarcodeModels.respop_R_to_S_transition_event!(state.cells,1,counts)
    @test state.cells[1].E && !state.cells[1].R # E is terminal.
    @test counts.Ecount == 1 && counts.Rcount == 1 && counts.Scount == 0

    dmg = ResDmg_ABM_EvBC(coupling_dmg_params())
    state = coupling_abm_state(dmg,1;spare=1)
    counts = BarcodeModels.ResDmgPhenoCounts(1,0,0,0)
    BarcodeModels.resdmg_S_to_R_transition_event!(state.cells,1,counts)
    @test counts.Scount == 0 && counts.Rcount == 1
    BarcodeModels.resdmg_birth_mutate_event_evbc!(state,1,2,0.0,0.0,counts,0.2)
    @test state.lineage_records[end].parent_pheno == "R"
    @test state.lineage_records[end].child_pheno == "R"
    BarcodeModels.resdmg_R_to_S_transition_event!(state.cells,1,counts)
    @test counts.Scount == 1 && counts.Rcount == 1
end

@testset "Independent transition public simulation APIs and primitive wrappers" begin
    simple = SimpleSimParams(n0=6,tmax=0.1,Nmax=100,Cc=100,Nswitch=1,
        Aswitch=0.01,save_at=0.05,drug_treatment=false)
    exp = ExperimentParams(n0=6,t_exp=0.1,tmax=0.1,t_Pass=Float64[],Nseed=6,
        Nmax=100,Cc=100,treat_ons=Float64[],treat_offs=Float64[],t_keep=[0.1],
        Nswitch=1,Aswitch=0.01,n_rep=1,save_at=0.05,full_sol=true,drug_treatment=false)
    abm = ABMParams(Nbuff=20)
    models = (ResPop(coupling_pop_params()),ResDmg(coupling_dmg_params()),
        ResPop_ABM(coupling_pop_params();abm),ResPop_ABM_EvBC(coupling_pop_params();abm),
        ResDmg_ABM(coupling_dmg_params();abm),ResDmg_ABM_EvBC(coupling_dmg_params();abm))
    for model in models
        @test BarcodeModels._with_drug_effect(model,:c).params.transition_coupling === :independent_rate
        Random.seed!(713)
        out = simulate_simple(model,simple;drug_effect=:c)
        @test haskey(out,"sol_df")
        @test all(x->x >= 6.0-1e-8,out["sol_df"].N)
        out = simulate_experiment(model,exp;drug_effect=:c)
        @test out isa Dict
        @test haskey(out,"sol_df") || haskey(out,"Nvec")
    end
    sol = BarcodeModels.simulate_grow_kill(6,6.0,0.0,0.0,
        1.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.1,Float64[],100,100,Float64[],Float64[],1;
        Aswitch=0.01,transition_coupling="independent_rate",save_at=0.05)
    @test sol.u[end][2] < 6.0
    @test sum(sol.u[end][2:4]) > 6.0
    popcells = BarcodeModels.seed_cells(4,0.0,4)
    out = BarcodeModels.grow_kill_lin_kmc!(popcells,1.0,0.0,2.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.1,100,4,Float64[],Float64[],0.01;
        transition_coupling="independent_rate")
    @test all(==(4),out.Nvec)
    dmgcells = BarcodeModels.seed_resdmg_cells(4,0.0,4)
    out = BarcodeModels.grow_kill_lin_kmc!(dmgcells,1.0,0.0,2.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,100,4,Float64[],Float64[],0.01;
        transition_coupling=:independent_rate)
    @test all(==(4),out.Nvec)
end

@testset "Omitted coupling matches explicit birth coupling for seeded simulations" begin
    pop_kw = (b=0.5,d=0.0,mu=0.2,sig=0.2,del=0.0,al=0.0,Dc=0.0,k=0.0,psi=0.0)
    dmg_kw = (b=0.5,d=0.0,mu=0.2,sig=0.2,del=0.0,ome=0.0,zet_S=0.0,
        zet_R=0.0,Dc=0.0,k=0.0,psi=0.0)
    for (model_type,param_type,kwargs) in (
        (ResPop,ResPopParams,pop_kw),(ResDmg,ResDmgParams,dmg_kw),
        (ResPop_ABM,ResPopParams,pop_kw),(ResPop_ABM_EvBC,ResPopParams,pop_kw),
        (ResDmg_ABM,ResDmgParams,dmg_kw),(ResDmg_ABM_EvBC,ResDmgParams,dmg_kw))
        models = (model_type(param_type(;kwargs...)),
            model_type(param_type(;kwargs...,transition_coupling=:birth_coupled)))
        if model_type in (ResPop,ResDmg)
            state = model_type === ResPop ? ResPopState(20,0,0) : ResDmgState(20,0,0,0)
            solutions = Any[]
            for model in models
                Random.seed!(572)
                push!(solutions,BarcodeModels.run_model_core_hybrid(model,state,
                    hybrid_test_sim(tmax=0.2);rng=MersenneTwister(572)))
            end
            @test solutions[1].t == solutions[2].t
            # Compare the biological state; ExtendedJumpArray also stores private clocks.
            nstate = model_type === ResPop ? 5 : 6
            observed = [[[u[i] for i in 1:nstate] for u in sol.u] for sol in solutions]
            @test observed[1] == observed[2]
        else
            states = [coupling_abm_state(model,20;spare=100) for model in models]
            outputs = Any[]
            for (model,state) in zip(models,states)
                Random.seed!(572)
                push!(outputs,BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.2,Cc=1000)))
            end
            @test outputs[1].tvec == outputs[2].tvec
            @test outputs[1].Nvec == outputs[2].Nvec
            contents = [map(c->Tuple(getfield(c,name) for name in fieldnames(typeof(c))),state.cells) for state in states]
            @test contents[1] == contents[2]
        end
    end
end

@testset "Independent ABM envelopes cover resistant costs and untreated damaged death" begin
    n,trials,tmax = 20,80,0.3
    pop = coupling_pop_params(b=2.0,d=0.4,mu=0.0,sig=0.0,del=0.5)
    dmg = coupling_dmg_params(b=2.0,d=0.4,mu=0.0,sig=0.0,del=0.5)
    for model in (ResPop_ABM(pop),ResPop_ABM_EvBC(pop),ResDmg_ABM(dmg),ResDmg_ABM_EvBC(dmg)),
        cost in ("b","d")
        birth,death = cost == "b" ? (1.2,0.4) : (2.0,1.2)
        factor = exp((birth-death)*tmax)
        expected = n*factor
        variance_per_cell = (birth+death)/(birth-death)*factor*(factor-1)
        tolerance = 6sqrt(n*variance_per_cell/trials)
        total = 0
        for trial in 1:trials
            Random.seed!(99300+trial)
            state = coupling_abm_state(model,n;rho=1.0,spare=100)
            sim = BarcodeModels.ABMSimParams(;tmax,Nmax=1000,Cc=10^9,
                treat_ons=Float64[],treat_offs=Float64[],R_real=cost,t_frac=0.1)
            BarcodeModels.run_model_core_abm(model,state,sim)
            total += count(c->c.alive,state.cells)
            @test all(c->!c.alive || c.R,state.cells)
        end
        observed = total/trials
        @info "Cost-envelope birth/death ensemble" model=typeof(model) cost observed expected tolerance
        @test abs(observed-expected) <= tolerance
    end
    expected = exp(-0.7)
    tolerance = 6sqrt(expected*(1-expected)/(n*trials))
    dmg = coupling_dmg_params(mu=0.0,sig=0.0,Dc=1.0)
    for model in (ResDmg_ABM(dmg),ResDmg_ABM_EvBC(dmg))
        survived = 0
        for trial in 1:trials
            Random.seed!(99100+trial)
            state = coupling_abm_state(model,n)
            for cell in state.cells
                cell.DS = true
            end
            BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(tmax=0.7,Cc=10^9))
            survived += count(c->c.alive,state.cells)
        end
        @info "Untreated damaged-death ensemble" model=typeof(model) frequency=survived/(n*trials) expected tolerance
        @test abs(survived/(n*trials)-expected) <= tolerance
    end
end

@testset "Seeded conversion ensembles agree with two-state analytical kinetics" begin
    n,trials,tmax = 20,80,0.7
    expected_R = 0.7*(-expm1(-tmax))
    tolerance = 6sqrt(expected_R*(1-expected_R)/(n*trials))
    pop,dmg = coupling_pop_params(),coupling_dmg_params()
    models = (ResPop(pop),ResDmg(dmg),ResPop_ABM(pop),ResPop_ABM_EvBC(pop),
        ResDmg_ABM(dmg),ResDmg_ABM_EvBC(dmg))
    for (model_index,model) in enumerate(models)
        resistant = 0.0
        for trial in 1:trials
            seed = 84000+1000model_index+trial
            Random.seed!(seed)
            if model isa Union{ResPop,ResDmg}
                state = model isa ResPop ? ResPopState(n,0,0) : ResDmgState(n,0,0,0)
                sol = BarcodeModels.run_model_core_hybrid(model,state,
                    hybrid_test_sim(Nswitch=1000,Aswitch=1000.0,tmax=tmax,Cc=n,Nmax=1000);
                    rng=MersenneTwister(seed))
                population_indices = model isa ResPop ? (2:4) : (2:5)
                @test sum(sol.u[end][population_indices]) == n
                resistant += sol.u[end][model isa ResPop ? 3 : 5]
            else
                state = coupling_abm_state(model,n)
                BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(;tmax,Cc=n))
                @test all(c->c.alive,state.cells)
                resistant += count(c->c.R,state.cells)
            end
        end
        frequency = resistant/(n*trials)
        @info "Two-state conversion ensemble" model=typeof(model) n trials frequency expected_R tolerance
        @test abs(frequency-expected_R) <= tolerance
    end
end

@testset "Deterministic conversions and numerical step sensitivity" begin
    n,tmax = 1000,0.02
    for model in (ResPop(coupling_pop_params()),ResDmg(coupling_dmg_params()))
        endpoints = Float64[]
        for (dtmax,tol) in ((0.001,1e-10),(0.01,1e-8))
            state = model isa ResPop ? ResPopState(n,100,0) : ResDmgState(n,0,0,100)
            sol = BarcodeModels.run_model_core_hybrid(model,state,
                hybrid_test_sim(Nswitch=1,Aswitch=0.01,Cc=1100,tmax=tmax,Nmax=10000);
                rng=MersenneTwister(614),dtmax=dtmax,abstol=tol,reltol=tol)
            # Start with 100 R so both directions stay deterministic throughout.
            expected = 1100*0.7 + (100-1100*0.7)*exp(-tmax)
            endpoint = sol.u[end][model isa ResPop ? 3 : 5]
            @test endpoint ≈ expected atol=1e-5
            population_indices = model isa ResPop ? (2:4) : (2:5)
            @test sum(sol.u[end][population_indices]) ≈ 1100.0 atol=1e-8
            push!(endpoints,endpoint)
        end
        @info "Deterministic conversion step sensitivity" model=typeof(model) endpoints maximum_difference=abs(endpoints[1]-endpoints[2])
        @test endpoints[1] ≈ endpoints[2] atol=1e-5
    end
    # R->E survival under gamma(t)=k*t until it reaches one.
    model = ResPop(coupling_pop_params(mu=0.0,sig=0.0,al=0.8,k=2.0,Dc=0.0))
    sim = BarcodeModels.SimParams(n0=n,tmax=0.8,Nmax=10000,Cc=n,Nswitch=1,
        Aswitch=0.01,treat=true,treat_ons=[0.0],treat_offs=Float64[],save_at=0.05)
    for dtmax in (0.001,0.005)
        sol = BarcodeModels.run_model_core_hybrid(model,ResPopState(0,n,0),sim;
            rng=MersenneTwister(72),dtmax,abstol=1e-10,reltol=1e-10)
        @test sol.u[end][3] ≈ n*exp(-0.8*(0.25+0.3)) atol=0.01
        @test sum(sol.u[end][2:4]) ≈ n atol=1e-8
        @info "Hybrid treatment-dependent conversion" dtmax observed_R=sol.u[end][3] expected_R=n*exp(-0.8*(0.25+0.3))
    end
end

@testset "Treatment-dependent R to E conversion and ABM lookup sensitivity" begin
    n,trials,tmax = 20,80,0.8
    # Retain the matched positive-onset comparison; onset zero is covered below.
    treat_on = 0.02
    expected_survival = exp(-0.8*(0.25+tmax-treat_on-0.5))
    tolerance = 6sqrt(expected_survival*(1-expected_survival)/(n*trials))
    params = coupling_pop_params(mu=0.0,sig=0.0,al=0.8,k=2.0,Dc=0.0)
    survivors = 0.0
    model = ResPop(params)
    sim = BarcodeModels.SimParams(;n0=n,tmax,Nmax=1000,Cc=n,Nswitch=1000,
        Aswitch=1000.0,treat=true,treat_ons=[treat_on],treat_offs=Float64[],save_at=0.1)
    for trial in 1:trials
        sol = BarcodeModels.run_model_core_hybrid(model,ResPopState(0,n,0),sim;
            rng=MersenneTwister(93200+trial),dtmax=0.01,abstol=1e-8,reltol=1e-8)
        @test sum(sol.u[end][2:4]) == n
        survivors += sol.u[end][3]
    end
    @info "Hybrid stochastic R to E ensemble" n trials frequency=survivors/(n*trials) expected_survival tolerance
    @test abs(survivors/(n*trials)-expected_survival) <= tolerance
    for model in (ResPop_ABM(params),ResPop_ABM_EvBC(params)), spacing in (0.01,0.002)
        survivors = 0
        for trial in 1:trials
            Random.seed!(91200+trial)
            state = coupling_abm_state(model,n;rho=1.0)
            sim = coupling_abm_sim(;tmax,Cc=n,dt_save_at=spacing,treat=true,treat_on=treat_on)
            BarcodeModels.run_model_core_abm(model,state,sim;treat=true)
            @test all(c->c.alive && xor(c.R,c.E),state.cells)
            survivors += count(c->c.R,state.cells)
        end
        @info "ABM treatment lookup ensemble" model=typeof(model) spacing n trials frequency=survivors/(n*trials) expected_survival tolerance
        @test abs(survivors/(n*trials)-expected_survival) <= tolerance
    end
    # No treatment means no R->E activity, even when alpha is large.
    model = ResPop_ABM(coupling_pop_params(mu=0.0,sig=0.0,al=10.0))
    state = coupling_abm_state(model,n;rho=1.0)
    BarcodeModels.run_model_core_abm(model,state,coupling_abm_sim(Cc=n))
    @test all(c->c.R && !c.E,state.cells)
end

@testset "ABM treatment starts at time zero" begin
    for onset in (0.0, 0.2)
        profile = BarcodeModels.drug_treat_concs(1.5, 2.0, 0.7, [onset], [0.9], 0.05)
        expected = [t <= 0.9 ? clamp(2*(t-onset), 0.0, 1.0) :
            max(1-2*(t-0.9), 0.0) for t in profile["t"]]
        @test profile["rconc"] ≈ expected atol=1e-8
        @test profile["dconc"] ≈ 0.7 .* expected atol=1e-8
        @test profile["rconc"][1] == 0.0
        @test profile["rconc"][end] ≈ 0.0 atol=1e-8
        @test maximum(profile["rconc"]) ≈ 1.0 atol=1e-8
    end
    for (ons, offs) in ((Float64[], Float64[]), ([0.0], [0.0]))
        profile = BarcodeModels.drug_treat_concs(0.2, 2.0, 0.7, ons, offs, 0.05)
        @test all(iszero, profile["rconc"])
        @test all(iszero, profile["dconc"])
    end
    # Both ABM variants consume this same profile for treatment-dependent conversion.
    params = coupling_pop_params(mu=0.0, sig=0.0, al=20.0, k=100.0)
    for model in (ResPop_ABM(params), ResPop_ABM_EvBC(params))
        Random.seed!(20261013)
        state = coupling_abm_state(model, 12; rho=1.0)
        sim = BarcodeModels.ABMSimParams(tmax=0.2, Nmax=1000, Cc=12,
            treat_ons=[0.0], treat_offs=Float64[], dt_save_at=0.01, t_frac=0.1)
        out = BarcodeModels.run_model_core_abm(model, state, sim; treat=true)
        @test any(c -> c.E, state.cells)
        @test all(c -> c.alive && xor(c.R, c.E), state.cells)
        @test all(==(12), out.Nvec)
    end
end
