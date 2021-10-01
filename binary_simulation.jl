using Optim, Random, Distributions, Suppressor, Statistics, LinearAlgebra
using FixedEffectModels, DataFrames
includet("estimator_binary.jl")

struct UnobservedType
        v::Array{Float64}
        c::Array{Float64}
end

function UnobservedType(n)
        v = rand(Uniform(7, 12), n)
        c = rand(Uniform(5, 10), n)
        return UnobservedType(v, c)
end

function good_demand(p, v)
        return (p.<v).*1.0
end

function good_supply(p, c, tech, subsidy::Bool)
        if subsidy
            return (p .- tech .> c)
        else
            return tech.*(p .> c)
        end

end

function eqm_price(U, types::UnobservedType, tech, subsidy)
        function eqm_condition(p)
            zp = mean(good_demand(p .+ U, types.v)) - mean(good_supply(p .+ U, types.c, tech, subsidy))
            return abs(zp)
        end
        r = optimize(eqm_condition, 0.1, 50.0).minimizer
        return r
end

#technology that increases production by 20% at the same cost
function simulate_eqm_data(types, n; subsidy= false, p = 0.5)
    #types = UnobservedType(n)
    W = rand(Bernoulli(p), n)
    #W = append!(ones(5000), zeros(5000))
    U = 0.18 .* rand([-1,1], n)
    #binary subsidy
    if subsidy
        tech = W.*0.5 + (1 .- W).*0.0
    #production technology intervention
    else
        tech = W*1.2 + (1 .- W).*1.0
    end
    pstar = eqm_price(U, types, tech, subsidy)
    #println("Ground truth direct: ", direct_ground(types, pstar))
    D = good_demand(pstar .+ U, types.v)
    S = good_supply(pstar .+ U, types.c, tech, subsidy)
    Y = S.*(pstar .+ tech.*subsidy .+ U .- types.c)
    return EqmData(Y, D, S, W, U)
end

function run_simulation(n, S)
        results = zeros((S, 4))
        for i in 1:S
                types = UnobservedType(n)
                data = simulate_eqm_data(types, n)
                estimates = estimate_τ(data)
                results[i, 1] = estimates.τD
                results[i, 2] = estimates.τI
                results[i, 3] = estimates.σD
                results[i, 4] = estimates.σI
        end

        gt = [0.22222222222, 0.7333333333333*35/10.5^2*(-1.0)]
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
        return results
end
