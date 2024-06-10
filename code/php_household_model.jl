#demographic info on number of children
#income with and without cash transfer
using Random, Distributions
using StatFiles, DataFrames, GLM
using Optim
using Plots, StatsPlots, LaTeXTabulars, LaTeXStrings

includet("estimator_binary.jl")
includet("hte_analysis.jl")

#demand model parameters
struct DemandModel
    θ0::Vector
    θp::Float64
    θw::Float64
end

#supply model parameters
struct SupplyModel
    θ0::Float64
    θp::Float64
end

#outcome model parameters
struct OutcomeModel
    θ0::Vector
    θd::Float64
    θw::Float64
end

struct Model
    d::DemandModel
    s::SupplyModel
    y::OutcomeModel
    μw::Float64
    μ_elig::Float64
    dfh::DataFrame
end

function Agent(nh, m::Model, augment=false)
    #elig, baby, child ϵd, ϵs, ϵy
        hhi = rand(1:nrow(m.dfh), nh)
        dfh = m.dfh[hhi, :]
        dfh.nbabies = Int64.(dfh.nbabies)
        n = sum(dfh.hsize)
        #change this to ed1 and ed0
        ϵd1 = rand(Normal(0, 1/3), n)
        ϵd0 = rand(Normal(0, 1/3), n)
        ϵy1 = rand(Normal(0, 1.0), n)
        ϵy0 = rand(Normal(0, 1.0), n)
        ϵdh = rand(Normal(0, 1/3), nh)
        ϵyh = rand(Normal(0, 1/3), nh)
        eligh = rand(Bernoulli(m.μ_elig), nh)
        #hsize = zeros(n); elig = zeros(n)

        k = 10
        X = rand(Normal(0, 1.5), (nh, k))

        if augment
            #ζ_ϵ = rand(Normal(0, 0.01), (nh, 4))
            ζ_ϵ = rand(Gumbel(), (nh, 4))
            ζ_ϵ[:, 1] .+= 1.0 #0.00.5
            ζ_ϵ[:, 3] .+= -1.0#0.0-0.5

            β = zeros((k, 4))
            #used to be all divided by 4
            β[2, 1] = 1.0
            β[3, 1] = -0.5
            β[1, 2] = 1.0
            β[3, 3] = 1.0
            β[4, 4] = 1.0
            S = X*β .+ ζ_ϵ
            Sstar = argmax(S, dims=2)
            ζ = [Sstar[i][2] for i in 1:nh]
        else
            ζ = ones(nh)
        end

        #println("proportions: ")
        #println(mean(ζ.==1))
        #println(mean(ζ.==2))
        #println(mean(ζ.==3))
        #println(mean(ζ.==4))



        hmap = reduce(vcat, [repeat([v], dfh.hsize[v]) for v in 1:nh])
        child = [vcat(repeat([1], dfh.nbabies[v]), repeat([0], dfh.hsize[v] - dfh.nbabies[v])) for v in 1:nh]
        child = reduce(vcat, child)

        dfh_out = DataFrame(hmap = collect(1:nh), ϵd1 = ϵdh, ϵd0 = ϵdh, ϵy1 = ϵyh, ϵy0 = ϵyh, hsize = dfh.hsize, elig = eligh, type = ζ,
                            hchild = dfh.nbabies, haschild=dfh.nbabies .>0)
        dfh_out = hcat(dfh_out, DataFrame(X, :auto))
        dfi = DataFrame(hmap = hmap, child = child )
        dfi = leftjoin(dfi, dfh_out, on=:hmap)


        dfi.ϵy1 .+= ϵy1
        dfi.ϵy0 .+= ϵy0
        dfi.ϵd1 .+= ϵd1
        dfi.ϵd0 .+= ϵd0
        #get child, elig, hh_index from a dataframe?
    #    for i in 1:n
        #    dfi.ϵd[i] += ϵd[i] + ϵdh[hmap[i]]
        #    ϵy1[i] = ϵy1[i] + ϵyh[hmap[i]]
        #    ϵy0[i]  = ϵy0[i]+ ϵyh[hmap[i]]
        #    hsize[i] = dfh.hsize[hmap[i]]
    #        elig[i] = eligh[hmap[i]]
    #    end

        #return DataFrame(elig = elig, child = child, hsize=hsize, ϵd = ϵd, ϵy1 = ϵy1, ϵy0 = ϵy0, h = hmap)
        return dfi, dfh_out
end

### Functions for estimating the model

function get_moments_from_real_data()
    dpath = "../data/philippines/"
    df = DataFrame(load(dpath*"moduleA-all_1.dta"))

    gps = groupby(df, :bgycode)
    vpx = combine(gps, :treated => mean, :exposure => mean, :lneggprice => mean,
                                :lnriceprice => mean, :remote_centile => mean, nrow, renamecols=false)
    vpx_remote = vpx[ (vpx.remote_centile .> 70), :]

    bgys = unique(vpx_remote.bgycode)
    #average egg price in highly saturated villages that are treated/controlled
    p1 = mean(vpx_remote.lneggprice[vpx_remote.treated .==1 ])
    p0 = mean(vpx_remote.lneggprice[vpx_remote.treated .==0 ])
    #println(exp(p1))
    #println(exp(p0))

    dfc = DataFrame(load(dpath*"moduleD.dta"))
    #dfc.Dnumegg = dfc.Dnumegg
    #compute increase in children's egg consumption 18 months - 60 months ; in eggs per week

    dfc.remote = [x in bgys for x in dfc.bgycode]
    #don't just exclude remote villages for better point estimates
    dfc_elig = dfc[(dfc.sample .==1) .& dfc.remote, :]
    dfc_inelig = dfc[(dfc.sample .==2) .& dfc.remote, : ]


    μ_elig = mean(vpx.exposure[(vpx.remote_centile .> 70)] )

    μw = mean(dfc_elig.treated)

    adult_adj = 1.4
    μd = (elig = mean(skipmissing(dfc_elig.Dnumegg[dfc_elig.treated .==0])),
          inelig = mean(skipmissing(dfc_inelig.Dnumegg[dfc_inelig.treated .==0])))
    μy = (elig = mean(skipmissing(dfc_elig.Wz_heightforage5[dfc_elig.treated .==0])),
                inelig = mean(skipmissing(dfc_inelig.Wz_heightforage5[dfc_inelig.treated .==0])))
    μd_a = (elig = adult_adj * μd.elig, inelig = adult_adj * μd.inelig)

    τd_c_elig = mean(skipmissing(dfc_elig.Dnumegg[dfc_elig.treated .==1])) - μd.elig
    τd_c_inelig = mean(skipmissing(dfc_inelig.Dnumegg[dfc_inelig.treated .==1])) - μd.inelig

    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9373515/
    #Adolescent/ adult intake is ~40% greater than for younger for the bottom tertile

    τg_pct = τd_c_elig/μd.elig
    τi_pct = τd_c_inelig/μd.inelig
    τd_a_elig = 1.4*μd.elig*(1+τg_pct) -
                    1.4*μd.elig
    τd_a_inelig = 1.4*μd.inelig*(1+τi_pct) -
                    1.4*μd.inelig

    #for children 6 - 36 months old
    τy_elig = mean(skipmissing(dfc_elig.Wz_heightforage5[dfc_elig.treated .==1])) - μy.elig
    τy_inelig = mean(skipmissing(dfc_inelig.Wz_heightforage5[dfc_inelig.treated .==1])) - μy.inelig

    dfh =DataFrame(load(dpath*"balance_data_final.dta"))

    hsize = dfh.hhsize[dfh.exposure65.==1]
    nbabies = dfh.child05[dfh.exposure65 .==1 ]
    nchild = dfh.child614[dfh.exposure65 .==1 ] .+ dfh.child05[dfh.exposure65 .==1]
    dfd = DataFrame(hsize = hsize, nbabies = nbabies, nchild = nchild)
    dropmissing!(dfd)
    dfd = dfd[dfd.nchild .> 0, :]
    # what are sample = 3 and 4 ?
    egg_elasticity = -1.3
    return (elast = egg_elasticity, μ_elig = μ_elig, hh = dfd,
            μy = μy, μd = μd, μd_a = μd_a,
            τy = (elig = τy_elig, inelig = τy_inelig), τd_c = (elig = τd_c_elig, inelig = τd_c_inelig),
            τd_a = (elig = τd_a_elig, inelig = τd_a_inelig), p0 = exp(p0), p1 = exp(p1), μw = μw)

end

function solve_price_coefficient(elast, μd, p0)
    delta = 1.01*p0
    dp = elast*μd/delta
end

function estimate_model()
    Γ = get_moments_from_real_data()
    #println(Γ)
    # demand coefficients
    μd0 = Γ.μ_elig*Γ.μd.elig + (1 - Γ.μ_elig)*Γ.μd.inelig
    θdp = solve_price_coefficient(Γ.elast, μd0, Γ.p0)
    θd0_elig = Γ.μd.elig - θdp*Γ.p0
    θd0_inelig = Γ.μd.inelig - θdp*Γ.p0
    θdw = Γ.τd_c.elig - θdp*(Γ.p1 - Γ.p0)

    dm = DemandModel([θd0_elig, θd0_inelig], θdp, θdw)

    #check τd_c.inelig versus implied by model
    #println("check: ", demand(0, 0, 0, Γ.p1, dm, 0) - demand(0, 0, 0, Γ.p0, dm, 0))


    μd1 = Γ.μ_elig*demand(dm, Γ.p1, 1, (elig =1, ϵd1  = 0, ϵd0 =0, type=1)) + (1 - Γ.μ_elig)*demand(dm, Γ.p1, 0, (elig=0, ϵd1  = 0, ϵd0 =0, type=1))
    θsp = (μd1 - μd0) / (Γ.p1 - Γ.p0)
    θs0 = μd0  - θsp*Γ.p0
    sm = SupplyModel(θs0, θsp)

    #outcome coefficients
    θyd = Γ.τy.inelig/Γ.τd_c.inelig #it matches this moment ; but not the point estimates
    θyw = Γ.τy.elig - θyd*Γ.τd_c.elig
    θy0_elig = Γ.μy.elig - θyd*Γ.μd.elig #- θyw*Γ.μy.elig
    θy0_inelig = Γ.μy.inelig - θyd*Γ.μd.inelig #- θyw*Γ.μy.inelig

    om = OutcomeModel([θy0_elig, θy0_inelig], θyd, θyw)



    # structure for model
    return Model(dm, sm, om, Γ.μw, Γ.μ_elig, Γ.hh)

    #alternative for non-linear model:
    # draw a bunch of households
    #Optimizer minimizes difference in moments from these households using below functions compared to real data moments

end


function expected_demand(dm, p, W, μ_elig)
    dat1 = (elig = 1, ϵd1 = 0, ϵd0 =0,  type =1)
    dat0 = (elig= 0, ϵd1 = 0, ϵd0 =0, type=1)
    μd = μ_elig*demand(dm, p, W, dat1) + (1 - μ_elig)*demand(dm, p, 0, dat0)
    return μd
end

function expected_outcome(om, p, D, W, μ_elig)
    dat1 = (elig = 1, child =1, ϵy1 = 0, ϵy0=0, type=1)
    dat0 = (elig = 0, child =1, ϵy1 = 0, ϵy0 = 0, type=1)
    μy = μ_elig*outcome(om, D, W, dat1) + ( 1- μ_elig)*outcome(om, D, 0, dat0)
    return μy

end


### Functions for sampling from the model


function demand(dm::DemandModel, p, W, dfr)
    D = dm.θ0[1]* dfr.elig + dm.θ0[2]*(1 - dfr.elig) + dm.θw *W* dfr.elig + dm.θp*p + dfr.ϵd0*W + dfr.ϵd1*(1 -W)

    #if dfr.type == 1
    #    D = D + 0.5*dm.θ0[1]*W*dfr.elig
    #end
    if dfr.type == 2
        D = D + 0.75*dm.θ0[1]*W*dfr.elig - dm.θw *W* dfr.elig
        #D = D - dm.θw *W* dfr.elig
    elseif dfr.type ==3
        D = D - 0.75*dm.θ0[1]*W*dfr.elig*dfr.child - dm.θw *W* dfr.elig
    elseif dfr.type ==4
        D = D - dm.θw*W*dfr.elig + 1.5*dm.θw*W*dfr.elig*(1 - dfr.child)#0.05*dm.θ0[1]*W*dfr.elig*(1 - dfr.child)# + 0.5*dm.θ0[1]*W*dfr.elig*dfr.child
    end

    return D
end

function supply(sm::SupplyModel, p, ϵs)
    #can also do log-log supply here to make things a bit more interesting
    q = sm.θ0 + sm.θp * p + ϵs
    return q
end



#child-level outcomes, sum over number of children
function outcome(om::OutcomeModel, D, W, dfr)
    Y = om.θ0[1]*dfr.elig + om.θ0[2]*(1 - dfr.elig) + om.θd * D + om.θw * W *dfr.elig + W*dfr.ϵy1 + (1 - W)*dfr.ϵy0

    if dfr.type == 2
        Y = Y + 0.2*W *dfr.elig #- 3/4*om.θd*D
    elseif dfr.type ==3
        Y = Y #- 3/4*om.θd*D
    elseif dfr.type ==4
        Y = Y + 2*om.θw*W*dfr.elig
    end

    #can add heterogeneity easily using some other part of dataframe row
    return Y*dfr.child
end



#function map_outputs_hh(Y, D, S, Yw, Dw, hmap)
    #gfd = groupby(DataFrame(Y=Y, D=D, S=S, h=hmap), :h)
#    agg = combine(gfd, :Y => sum, :D => sum, :S => sum, renamecols=false)
#    return agg.Y .*(maximum(hmap)/Yw), agg.D .*(maximum(hmap)/Dw), agg.S.*(maximum(hmap)/Dw)

#end




function map_to_hh(Y, factor, hmap)
    nh = maximum(hmap)
    Yh = zeros(nh)
    for i in 1:length(hmap)
        Yh[hmap[i]] += Y[i] * nh/factor
    end
    return Yh
end

function map_treatments(W, U, hmap)
    return (map_to_individual(W, hmap), map_to_individual(U, hmap))
end

function map_to_individual(Vh, hmap)
    return [Vh[hmap[i]] for i in 1:length(hmap)]
end

#dfi is pre-treatment variables
function aggregate_outcome(om::OutcomeModel, D::Vector, W::Vector, df::DataFrame)
    Y = [outcome(om, D[i], W[i], df[i, :]) for i in 1:nrow(df)]
    return Y
end

function aggregate_demand(md::DemandModel, p::Float64, W::Vector, df::DataFrame)
    return aggregate_demand(md, p .+ zeros(nrow(df)), W, df)
end

function aggregate_demand(md::DemandModel, p::Vector, W::Vector, df::DataFrame)
    D = [demand(md, p[i], W[i], df[i, :]) for i in 1:nrow(df)]
    return D
end

function aggregate_supply(sm::SupplyModel, p::Float64, df::DataFrame)
    return aggregate_supply(sm, p .+ zeros(nrow(df)), df)
end

function aggregate_supply(sm::SupplyModel, p::Vector, df::DataFrame)
    S = [supply(sm, p[i], 0) for i in 1:nrow(df)]
    return S
end


#this is the main bottleneck, could be made much faster
function eqm_price(m::Model, df::DataFrame, W, U)
        noiseU = mean(m.d.θp*U) - mean(m.s.θp*U)
        dmμ = mean([demand(m.d, 6, W[i], df[i, :]) for i in 1:length(W)]) - 6*m.d.θp
        smμ = m.s.θ0
        function eqm_condition(p)
            zp =dmμ - smμ + noiseU + m.d.θp*p - m.s.θp*p
            return zp
        end
        r = fzero(eqm_condition, 6)
        return r
end


#augment adds additional heterogeneity beyond what is identifiable in the data
function sample_eqm_data(m::Model, df::DataFrame, dfh::DataFrame)
    nh = nrow(dfh)
    #types = UnobservedType(n)
    Wh = rand(Bernoulli(m.μw), nh)
    #W = append!(ones(5000), zeros(5000))
    Uh = 0.15 .* rand([-1,1], nh)



    W, U = map_treatments(Wh, Uh, df.hmap)
    pstar = eqm_price(m, df, W, U)

    #generate some additional covariates for heterogeneity simulation
    #xadd = rand(Normal(), (nh, 10))
    #df.x2 = map_to_individual(xadd[:, 2], df.h)
    #df.x3 = map_to_individual(xadd[:, 3], df.h)
    D = aggregate_demand(m.d, pstar .+ U, W, df)
    S = aggregate_supply(m.s, pstar .+ U, df)
    Y = aggregate_outcome(m.y, D, W, df)

    Yh = map_to_hh(Y, sum(df.child), df.hmap)
    Dh = map_to_hh(D, nrow(df), df.hmap)
    Sh = map_to_hh(S, nrow(df), df.hmap)

    k = 10
    cols = ["x"*string(i) for i in 1:k]
    cols = vcat(cols, ["hchild", "haschild"])
    Xh = Matrix(dfh[:, cols])
    #hajek has p hat rather than m.μw for Horvtiz-Thompson
    #println(m.μw)
    #println(mean(Wh))
    return EqmData(Yh, Dh, Sh, Wh, Uh, Xh, mean(Wh))
end


function compute_population_gt(m::Model)
    function objectiveπ(p)
        zp = expected_demand(m.d, p, m.μw, m.μ_elig) - supply(m.s, p, 0)
        return abs(zp)
    end

    function objective1(p)
        zp  = expected_demand(m.d, p, 1, m.μ_elig) - supply(m.s, p, 0)
        return abs(zp)
    end

    function objective0(p)
        zp  = expected_demand(m.d, p, 0, m.μ_elig) - supply(m.s, p, 0)
        return abs(zp)
    end

    p1 = optimize(objective1, 5.0, 8.0).minimizer
    p0 = optimize(objective0, 5.0, 8.0).minimizer
    pπ = optimize(objectiveπ, 5.0, 8.0).minimizer


    D1π = expected_demand(m.d, pπ, 1, m.μ_elig)
    D0π = expected_demand(m.d, pπ, 0, m.μ_elig)
    Y1π = expected_outcome(m.y, pπ, D1π, 1, m.μ_elig)
    Y0π = expected_outcome(m.y, pπ, D0π, 0, m.μ_elig)

    D1 = expected_demand(m.d, p1, 1, m.μ_elig)
    D0 = expected_demand(m.d, p0, 0, m.μ_elig)

    Y1 = expected_outcome(m.y, p1, D1, 1, m.μ_elig)
    Y0 = expected_outcome(m.y, p0, D0, 0, m.μ_elig)
    return (gte=Y1 - Y0, ade = Y1π - Y0π)

end


function compute_sample_ade(df::DataFrame, m::Model)
    n = nrow(df)
    Wh = rand(Bernoulli(m.μw), maximum(df.hmap))
    W = map_to_individual(Wh, df.hmap)
    #println(length(W))
    W1 = ones(n)
    W0 = zeros(n)
    p = eqm_price(m, df, W, zeros(n))

    #println(p)
    D0 = aggregate_demand(m.d, p, W0, df)

    D1 = aggregate_demand(m.d, p, W1, df)
    Y0 = aggregate_outcome(m.y, D0, W0, df)
    Y1 = aggregate_outcome(m.y, D1, W1, df)

    Y1h = map_to_hh(Y1, sum(df.child), df.hmap)
    Y0h = map_to_hh(Y0, sum(df.child), df.hmap)
    return mean(Y1h) - mean(Y0h)
end

function compute_sample_gt(df::DataFrame, m::Model)
    n = nrow(df)
    W1 = ones(n)
    W0 = zeros(n)
    p1 = eqm_price(m, df, W1, zeros(n))
    #println(p1)
    p0 = eqm_price(m, df, W0, zeros(n))
    #println(p0)
    D0 = aggregate_demand(m.d, p0, W0, df)
    D1 = aggregate_demand(m.d, p1, W1, df)
    Y0 = aggregate_outcome(m.y, D0, W0, df)
    Y1 = aggregate_outcome(m.y, D1, W1, df)

    Y1h = map_to_hh(Y1, sum(df.child), df.hmap)
    Y0h = map_to_hh(Y0, sum(df.child), df.hmap)
    return mean(Y1h) - mean(Y0h)
end

function compute_gts(m::Model, agents::DataFrame)
    gte = compute_sample_gt(agents, m)
    ade =  compute_sample_ade(agents, m)
    #gte faster than ade
    #aie = gte - ade
    return (gte = gte, ade = ade)
end

function run_c_check()
    n = 1000
    model = estimate_model()
    agts, hh = Agent(n, model, true)
    data =  sample_eqm_data(model, agts, hh)
    gt_ite(agts, data.W, model)
end


function gt_ite(df::DataFrame, Wh, m::Model)
    n = nrow(df)
    W = map_to_individual(Wh, df.hmap)
    pstar = eqm_price(m, df, W, zeros(n))
    W1 = ones(n)
    W0 = zeros(n)
    D0 = aggregate_demand(m.d, pstar, W0, df)
    D1 = aggregate_demand(m.d, pstar, W1, df)

    ited = D1 .- D0

    D1h = map_to_hh(D1, nrow(df), df.hmap)
    D0h = map_to_hh(D0, nrow(df), df.hmap)

#    println(mean(ited[df.type.==2]))
#    println(mean(ited[df.type.==3]))
    #println("ATED", mean(ited))

    itedh = map_to_hh(ited, nrow(df), df.hmap)

    Y0 = aggregate_outcome(m.y, D0, W0, df)
    Y1 = aggregate_outcome(m.y, D1, W1, df)

    Y0h = map_to_hh(Y0, sum(df.child), df.hmap)
    Y1h = map_to_hh(Y1, sum(df.child), df.hmap)

    itey = Y1h .- Y0h

    println(sum(itey .==0))
    println(mean(itedh[itey.==0]))

    #println(sum(itey .>0))
    #println(mean(itedh .* (itey .> 0)))
    #println(mean(itedh)*mean(Wh))

    #println("here")
    #plot(itedh[1:100], itey[1:100])

    #println("how many eq", sum(itey .== itedh))

    function demand_at_c(c)
        S = itey .> (itedh .* c)

        #println(sum(S))

        rule = mean(D1h .* S) + mean(D0h .*( 1 .- S))
        statusq = mean(D1h).* mean(Wh) + mean(D0h) .* ( 1 - mean(Wh))
        #println("rule", rule)
        #println("sq", statusq)
        return statusq - rule
    end

    #println("datc0: ", demand_at_c(0.0))
    #println("datc1: ", demand_at_c(0.0000001))
    #println("datc10: ", demand_at_c(10.0))
    cstar =  find_zero(demand_at_c,(0.0, 10.0))
    #println("CSTAR TRUTH", cstar)


    return (Y1h, Y0h)
end
####SIMULATION CODE###

function supply_plot()
    rangep = 5:0.01:7
    model = estimate_model()
    D0 = [expected_demand(model.d, p, 0, model.μ_elig) for p in rangep]
    D1 = [expected_demand(model.d, p, 1, model.μ_elig) for p in rangep]
    S = [supply(model.s, p, 0) for p in rangep]
    plot(Vector(rangep), [D0 D1 S], label=["Aggregate Demand (Control)" "Aggregate Demand (Treated)" "Aggregate Supply"], legend=:top, linewidth=2.0, xlabel="Price", ylabel="Quantity")
    savefig("../exhibits/fig1a.pdf")
end

function ade_plot(n, S)
    m = estimate_model()
    τ_truth = compute_population_gt(m).gte
    results = zeros((S, 2))
     for i in 1:S
        agents, hhs = Agent(n, m)
        estimates = compute_gts(m, agents)
        results[i, :] .= collect(estimates)
    end

    density(results[:, 1:2], label=["Sample Total Treatment Effect" "Sample Average Direct Effect"], linewidth=2.0, linestyle = [:solid :dash], legend=:outertop, xlabel="Treatment Effect", ylabel="Density")
    vline!([τ_truth], label= "Population Total Treatment Effect", linestyle = :dot, linecolor=:black, linewidth=2.0)
    savefig("../exhibits/fig1b.pdf")
end



function hte_plot(n)
    model = estimate_model()
    agts, hh = Agent(n, model, true)
    data = sample_eqm_data(model, agts, hh)

    test = gt_ite(agts, data.W, model)
    cstar, Yhat, Zhat, a, b = compute_rule(data.Y, data.D, data.X, data.W, mean(data.W), data.X)
    #println(sum(Zhat .< 0))
    #println("Type 3,", mean(Zhat[hh.type .== 3]))
    #println("Type 4,", mean(Zhat[hh.type .== 4]))
    #println("Type 2,", mean(Zhat[hh.type .== 2]))
    samplesub =  rand(1:length(Zhat), 250)
    W = Yhat[samplesub] .> Zhat[samplesub] .* cstar
    plot(Zhat[samplesub][W .== 1], Yhat[samplesub][W .==1], seriestype=:scatter, framestyle=:origin, xlabel="ZCADEᵢ",
            ylabel="CADEᵢ", label="Treated")
    plot!(Zhat[samplesub][W .== 0], Yhat[samplesub][W .==0], seriestype=:scatter, framestyle=:origin, label="Control")
    xs = Vector(minimum(Zhat):0.01:maximum(Zhat))
    ys = xs.*cstar
    println(cstar)
    plot!(xs, ys, label="Treatment Rule", legend=:topright, linewidth=2.0, xlims=(-1, 3), ylims =(-0.5, 1) )
    savefig("../exhibits/fig2.pdf")
end


#coverage_simulation(2000, 10)
function coverage_simulation(n, S)
    #AIE overcovers a bit; could be because of perturbation size; try adjusting?
        results = zeros((S, 5))
        model = estimate_model()
        gtsample = zeros((S, 2))
        for i in 1:S
                agents, hh = Agent(n, model)
                data = sample_eqm_data(model, agents, hh)
                estimates = estimate_τ(data)

                gts = compute_gts(model, agents)
                gtsample[i, 1] = gts[2]
                gtsample[i, 2] = gts[1] - gts[2]
                results[i, 1] = estimates.τD
                results[i, 2] = estimates.τI
                results[i, 3] = estimates.σD
                results[i, 4] = estimates.σI
                results[i, 5] = estimates.σRCT
        end

        popgt = compute_population_gt(model)
        #gte_ade =
        gtpop = [popgt.ade, popgt.gte - popgt.ade]
        #println(gt)
        header = ["", "Estimate", "Bias", "S.D.", L"Coverage for $\tau$", L"Coverage for $\tau^*$"]
        content = Array{Any}(zeros(2, 6))
        content[1, 1] = L"\hat\tau_{\text{ADE}}"
        content[2, 1] = L"\hat\tau_{\text{AIE}}"

        #fill in table here
        for j in 1:2
            content[j, 2] = round(mean(results[:, j]), digits = 3)
            content[j, 3] =  round(mean(results[:, j] .- gtpop[j]), digits=3)
            content[j, 4] = round(std(results[:, j]), digits=3)
            ciu = results[:, j] .+ results[:, j+2]*1.96
            cil = results[:, j] .- results[:, j+2]*1.96
            content[j, 5] = round(mean((gtsample[:, j] .< ciu) .& (gtsample[:, j] .> cil)), digits=3)
            content[j, 6] = round(mean((gtpop[j] .< ciu) .& (gtpop[j] .> cil)), digits=3)
        end

        #ciu = results[:, 1] .+ results[:, 5]*1.96
        #cil = results[:, 1] .- results[:, 5]*1.96
        #covg = (gt[1] .< ciu) .& (gt[1] .> cil)
        #println("Coverage RCT: ", mean(covg))
        #println("CI Width RCT:", mean(ciu .- cil))
        println(content)
        output_latex(header, content, "../exhibits/table1.tex")
        #return results
end


function hte_gain(S)
    model = estimate_model()

    outcomeNaive = zeros(S)
    outcomeRule = zeros(S)
    outcomeRandom = zeros(S)

    gains = zeros(S)
    gainsNaive = zeros(S)
    eqm = zeros(S)
    eqmNaive = zeros(S)

    cs = zeros(S)
    as = zeros(S)

    for s in 1:S

        agts, hh = Agent(5000, model, true)
        data = sample_eqm_data(model, agts, hh)
        agtsb, hh = Agent(5000, model, true)
        datab = sample_eqm_data(model, agtsb, hh)

        cstar, Yhat, Zhat, Yhatb, Zhatb, astar = compute_rule(data.Y, data.D, data.X, data.W, mean(data.W), datab.X)
        println(sum(Zhat .< 0))
        #println(astar)
        Y1b, Y0b = gt_ite(agtsb, datab.W, model)
        as[s] = astar
        cs[s] = cstar

        outcomeNaive[s] =rule_value(Y1b, Y0b, (Yhatb.> Zhatb).*astar)
        outcomeRule[s] = rule_value(Y1b, Y0b, Yhatb .> Zhatb .* cstar)
        outcomeRandom[s] = rule_value(Y1b, Y0b, mean(data.W))
        #println(outcomeNaive)
        #println(outcomeRule)
        #println(outcomeRandom)
        gains[s] = compute_gain(datab.Y, datab.W, Yhatb .> Zhatb .* cstar, mean(data.W))


        gainsNaive[s] = compute_gain(datab.Y, datab.W, (Yhatb .> 0) .* astar, mean(data.W))
        eqm[s] = compute_gain(datab.Z, datab.W, Yhatb .> Zhatb .* cstar, mean(data.W))
        eqmNaive[s] = compute_gain(datab.Z, datab.W, (Yhatb .> 0) .* astar, mean(data.W))
    end


    #c and equilibrium stability
    header = ["", "Optimal Equilibrium-Stable Rule", "Direct Rule", "Randomized Rule"]
    content = Array{Any}(zeros(2, 4))

    content[1, 1] = "Height-For-Age Z-Score"
    content[2, 1] = "Standard Deviation"
    content[1, 2] = round(mean(outcomeRule), digits=3)
    content[2, 2] = round(std(outcomeRule), digits=3)
    content[1, 3] = round(mean(outcomeNaive), digits=3)
    content[2, 3] = round(std(outcomeNaive), digits=3)
    content[1, 4] = round(mean(outcomeRandom), digits=3)
    content[2, 4] = round(std(outcomeRandom), digits=3)


    println(content)

    output_latex(header, content, "../exhibits/table2.tex")
    #println("Gains: ", mean(gains),", ", std(gains))
    println("Gains: ", mean(gains),", ", std(gains))
    #println("Eqm: ", mean(eqm),", ", std(eqm))
    println("Eqm: ", mean(eqm),", ", std(eqm))
    println("Cstar: ", mean(cs),", ", std(cs))
    println("EqmNaive :", mean(eqmNaive), ", ", std(eqmNaive))
    println("astar : ", mean(as), ", ", std(as))
end



function output_latex(header, content::Matrix, tname::String)
    latex_tabular(tname,
              Tabular("l"*repeat("r", length(header)-1)),
              [Rule(:top),
               header,
               Rule(:mid),
               content,
               Rule(:bottom)])
end

function generate_model_estimation_table()
    Γ = get_moments_from_real_data()
    m = estimate_model()

    header1 = MultiColumn(10, :c, "Estimated Moments from \\citet{filmer2023}")
    header2 = MultiColumn(10, :c, "Estimated Parameters")

    contentmom = [Γ.μd.elig, Γ.μd.inelig, Γ.μy.elig, Γ.μy.inelig, Γ.p0, Γ.p1, Γ.elast, Γ.τy.elig, Γ.τd_c.elig, Γ.τy.inelig/Γ.τd_c.inelig]
    contentmom = [round(x, digits = 2) for x in contentmom]
    contenttheta = [m.d.θ0[2], m.d.θ0[1], m.d.θw, m.d.θp, m.s.θ0, m.s.θp, m.y.θ0[2], m.y.θ0[1], m.y.θd, m.y.θw]
    contenttheta = [round(x, digits=2) for x in contenttheta]
    headermom = [L"\hat \mu^d_{elig}", L"\hat \mu^d_{inelig}", L"\hat \mu^y_{elig}",
                    L"\hat \mu^y_{inelig}", L"\hat p^0", L"\hat p^1", L"\eta_{egg}", L"\hat \tau^y_{elig}",
                     L"\hat \tau^d_{elig}", L"\frac{ \hat \tau^y_{inelig}} {\hat \tau^d_{inelig}}"]

    headertheta = [L"\hat \theta_{d00}",  L"\hat \theta_{d01}", L"\hat \theta_{dw}", L"\hat \theta_{dp}",
                   L"\hat \theta_{s0}",  L"\hat \theta_{sp}", L"\hat \theta_{y00}", L"\hat \theta_{y01}",
                   L"\hat \theta_{yd}", L"\hat \theta_{yw}"]

    latex_tabular("../exhibits/table3.tex",
                 Tabular("l"*repeat("l", length(headermom)-1)),
                 [ Rule(:top), 
                    [header1],
                  Rule(:mid),
                  headermom,
                  Rule(:mid),
                  contentmom,
                  Rule(:mid),
                  Rule(:mid),
                  [header2],
                  Rule(:mid),
                  headertheta,
                  Rule(:mid),
                  contenttheta,
                  Rule(:bottom)]
                 )
end



 #$\hat \mu^d_{elig} $ & $ \hat \mu^d_{inelig}$ & $\hat \mu^y_{elig}$ & $\hat \mu^y_{inelig}$ & $\hat p^0$ & $\hat p^1$ & $ \eta_{egg}$ &  $\hat \tau^y_{elig}$ & $\hat \tau^d_{elig} $& $\frac{ \hat \tau^y_{inelig}} {\hat \tau^d_{inelig}}$  \\

# 1.60 & 2.20 & -2.47 & -1.44 & 6.07 & 6.26 & -1.30 & 0.32 & 0.12 & 1.85 \\
#$\hat \theta_{d00}$ & $\hat \theta_{d01}$ & $\hat \theta_{dw}$ &$ \hat \theta_{dp}$ & $ \hat \theta_{s0}$ & $ \hat \theta_{sp}$ & $ \hat \theta_{y00} $ & $\hat \theta_{y01}$ & $\hat \theta_{yd} $&$ \hat \theta_{yw}$ \\
#4.49 & 3.89 & 0.19 & -0.38 & -0.21 & 0.32 & -5.36 & -5.18 & 1.85 & 0.10 \\
