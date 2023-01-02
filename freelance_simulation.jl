using Optim, Random, Distributions, Suppressor, Statistics, LinearAlgebra
using FixedEffectModels, DataFrames, Plots, StatsPlots
includet("estimator_binary.jl")

struct UnobservedType
        v::Array{Float64}
        c1::Array{Float64}
        c2::Array{Float64}
end

function UnobservedType(n)
        v = rand(Uniform(1, 2), n)
        c1 = rand(Uniform(1, 2), n)
        c2 = rand(Uniform(0.0, 1.0), n)
        return UnobservedType(v, c1, c2)
end

function aggregate_demand(V, p::Vector)
    return [good_demand(p[i], V[i]) for i in 1:length(V)]
end

function aggregate_supply(C1, C2, p::Vector, W)
    return [good_supply(p[i], C1[i], C2[i], W[i]) for i in 1:length(C1)]
end

function aggregate_demand(V, p::Float64)
    return [good_demand(p, V[i]) for i in 1:length(V)]
end

function aggregate_supply(C1, C2, p::Float64, W)
    return [good_supply(p, C1[i], C2[i], W[i]) for i in 1:length(C1)]
end

function good_demand(p, v)
        return  (p < v)
end

function good_supply(p, c1, c2, W)
        if Bool(W)
            threshold = 0.3
        else
            threshold = 0.7
        end
        return (p - c1 >0) + (p - c1>0)*(c2 >threshold)
end

function eqm_price(U, types::UnobservedType, W)
        function eqm_condition(p)
            zp = mean(aggregate_demand(types.v, p .+ U)) - mean(aggregate_supply(types.c1, types.c2, p .+ U, W))
            return abs(zp)
        end
        r = optimize(eqm_condition, 1.0, 5.0).minimizer
        return r
end

function compute_gts(types, prnt=false)
    gte = compute_global_effect(types)
    lte = compute_lte(types)
    dte = compute_direct_gt(types)
    eqmd = sample_eqm_data(types)

    if prnt
        println("gte: ", gte)
        println("lte-approx: ", lte)
        println("dte: ", dte)
        println("indirect: ", gte - dte)
    end
    #τ = estimate_τ(eqmd)
    #println(τ)
    #println("direct effect", τ.τD)
    #println("indirect effect", τ.τI)
    #println("lte", τ.τD + τ.τI)
    return (gte=gte, dte = dte)
end

function compute_lte(types)
    n = length(types.c1)
    W1 = rand(Bernoulli(0.55), n)
    W0 = rand(Bernoulli(0.45), n)
    p1 = eqm_price(zeros(n), types, W1)
    p0 = eqm_price(zeros(n), types, W0)
    S0 = aggregate_supply(types.c1, types.c2, p0, W0)
    S1 = aggregate_supply(types.c1, types.c2, p1, W1)
    Y0 = S0.*(p0)
    Y1 = S1.*(p1)
    return(mean(Y1) - mean(Y0))/0.1
end

function compute_direct_gt(types)
    n = length(types.c1)
    W = rand(Bernoulli(0.5), n)
    W1 = ones(n)
    W0 = zeros(n)
    p = eqm_price(zeros(n), types, W)
    S0 = aggregate_supply(types.c1, types.c2, p, W0)
    S1 = aggregate_supply(types.c1, types.c2, p, W1)
    Y0 = S0.*(p)
    Y1 = S1.*(p)
    return mean(Y1) - mean(Y0)
end

function compute_global_effect(types)
    n = length(types.c1)
    W1 = ones(n)
    W0 = zeros(n)
    p1 = eqm_price(zeros(n), types, W1)
    p0 = eqm_price(zeros(n), types, W0)
    S0 = aggregate_supply(types.c1, types.c2, p0, W0)
    S1 = aggregate_supply(types.c1, types.c2, p1, W1)
    #println("p0: ", p0)
    #println("p1:", p1)
    Y0 = S0.*(p0)
    Y1 = S1.*(p1)
    return mean(Y1) - mean(Y0)
end

#technology that increases production by 20% at the same cost
function sample_eqm_data(types; pi = 0.5)
    n = length(types.c1)
    #types = UnobservedType(n)
    W = rand(Bernoulli(pi), n)
    #W = append!(ones(5000), zeros(5000))
    U = 0.05 .* rand([-1,1], n)
    pstar = eqm_price(U, types, W)
    #println(pstar)
    #println("Ground truth direct: ", direct_ground(types, pstar))
    D = aggregate_demand(types.v, pstar .+ U)
    S = aggregate_supply(types.c1, types.c2, pstar .+ U, W)
    Y = S.*(pstar .+ U)
    return EqmData(Y, D, S, W, U)
end

function supply_plots(n)
    rangep = 1:0.01:2
    W1 = ones(n)
    W0 = zeros(n)
    types = UnobservedType(n)
    D = [mean(aggregate_demand(types.v, p)) for p in rangep]
    S0 = [ mean(aggregate_supply(types.c1, types.c2, p, W0)) for p in rangep]
    S1 = [mean(aggregate_supply(types.c1, types.c2, p, W1)) for p in rangep]

    plot(Vector(rangep), [D S0 S1], label=["Demand" "Aggregate Supply (Control)" "Aggregate Supply (Treated)"], legend=:top, linewidth=2.0, xlabel="Price", ylabel="Quantity")
    savefig("fig1b.pdf")
end

function ade_plot(n, S)
    τ_truth = 0.051786
    results = zeros((S, 2))
    for i in 1:S
        types = UnobservedType(n)
        estimates = compute_gts(types)
        results[i, :] .= collect(estimates)
    end

    density(results[:, 1:2], label=["Sample Global Treatment Effect" "Sample Average Direct Effect"], linewidth=2.0, linestyle = [:solid :dash])
    vline!([τ_truth], label= "Population GTE", linestyle = :dot, linecolor=:black, linewidth=2.0)
    savefig("fig1a.pdf")
end

function coverage_simulation(n, S)
        results = zeros((S, 5))
        for i in 1:S
                types = UnobservedType(n)
                data = sample_eqm_data(types)
                estimates = estimate_τ(data)
                results[i, 1] = estimates.τD
                results[i, 2] = estimates.τI
                results[i, 3] = estimates.σD
                results[i, 4] = estimates.σI
                results[i, 5] = estimates.σRCT
        end

        gt = [0.22375087334020516, -0.1719645337622137]
        for j in 1:2
            println("Effect Size: ", mean(results[:, j]))
            println("Bias: ", mean(results[:, j].- gt[j]))
            println("SD: ", std(results[:, j]))
            ciu = results[:, j] .+ results[:, j+2]*1.96
            cil = results[:, j] .- results[:, j+2]*1.96
            covg = (gt[j] .< ciu) .& (gt[j] .> cil)
            println("Coverage: ", mean(covg))
            println("CI Width: ", mean(ciu .- cil))
        end


        ciu = results[:, 1] .+ results[:, 5]*1.96
        cil = results[:, 1] .- results[:, 5]*1.96
        covg = (gt[1] .< ciu) .& (gt[1] .> cil)
        println("Coverage RCT: ", mean(covg))
        println("CI Width RCT:", mean(ciu .- cil))
        return results
end
