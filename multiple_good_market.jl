#simulation of subsidizing
using Optim, LinearAlgebra, Random, Distributions, Suppressor, Statistics

struct UnobservedType
    v::Array{Float64}
    c::Array{Float64, 2}
end

function UnobservedType(n)
    v = rand(Uniform(7, 15), n)
    cW = rand(Uniform(5, 10), n)
    cM = rand(Uniform(2, 5), n)
    return UnobservedType(v, [cW cM])
end

function calculate_ground_truth()
    n = 1000000
    types = UnobservedType(n)
    p0 = eqm_price(zeros(n, 2), types)
    function direct_effect(η)
        y = get_supplier_welfare(transpose(p0) .+ zeros(n, 2), types; η = η)
        return mean(y)
    end

    function total_effect(η)
        p1 = eqm_price(zeros(n, 2), types; η = η)
        y1 = get_supplier_welfare(transpose(p1) .+ zeros(n, 2), types; η = η)
        return mean(y1)
    end
    #println(eqm_price(zeros(n, 2), types; η = 0.1) - eqm_price(zeros(n, 2), types; η = 0.0))
    println((total_effect(0.2) - total_effect(-0.2))/0.4)
    println((direct_effect(0.2) - direct_effect(-0.2))/0.4)
    #τ = central_fdm(5, 1; factor=1e6)(total_effect, 0.0)
    #τD = central_fdm(5, 1; factor=1e6)(direct_effect, 0.0)

    #return (τ = τ, τD = τD)
end

function mc_simulation(S, n)
    results = zeros((S, 7))
    for i in 1:S
        estimates = run_simulation(n)
        results[i, :] .= collect(estimates)
    end
    return results
end

function run_simulation(n)
    types = UnobservedType(n)
    return calculate_effects(types)
end

function get_supplier_welfare(P, types; η= 0.0)
    qs = aggregate_supply(P, types; η=η)
    return qs[:, 1].*(P[:, 1] .+ η .- P[:, 2] .- types.c[:, 1])
end

# calculate general eqm effect using local experimemnt
function calculate_effects(types)
    n = length(types.v)
    U = 0.6 .*rand([-1,1], (n, 2))
    pstar = eqm_price(U, types)
    Z = aggregate_excess_demand(transpose(pstar) .+ U, types)
    qs = aggregate_supply(transpose(pstar) .+ U, types)
    qd = aggregate_demand(transpose(pstar) .+ U, types)
    ys = get_supplier_welfare(transpose(pstar) .+ U, types)
    ∂Zπ = zeros(2, 1)
    ∂Zp = zeros(2, 2)
    ∂ys = inv(transpose(U)*U)*(transpose(U)*ys)

    for i in 1:2
        ∂Zp[:, i] .= inv(transpose(U)*U)*(transpose(U)*Z[:, i])
    end
    ∂Zπ[1] = -1.0 .* inv(transpose(U[:, 1])*U[:, 1])*(transpose(U[:, 1])*qs[:, 1])
    ∂Zπ[2] = inv(transpose(U[:, 1])*U[:, 1])*(transpose(U[:, 1])*Z[:, 2])
    ∂pstar = - inv(∂Zp)*(∂Zπ)
    τI = transpose(∂ys)*∂pstar
    dy = ∂ys[1] + τI[1]
    return (τ = dy, τD  = ∂ys[1], τI = τI[1], ∂pstarG = ∂pstar[1], ∂pstarH = ∂pstar[2],
            ∂ysG = ∂ys[1], ∂ysW = ∂ys[2])
end

function aggregate_demand(P, types)
    n = size(types.c)[1]
    qd = zeros((n, 2))
    for i in 1:n
        qd[i, :] .=  individual_excess_demand(P[i, :], types.c[i, :], types.v[i]).qd
    end
    return qd
end

function aggregate_supply(P, types; η=0.0)
    n = size(types.c)[1]
    qs = zeros((n, 2))
    for i in 1:n
        qs[i, :] .=  individual_excess_demand(P[i, :], types.c[i, :], types.v[i]; η= η).qs
    end
    return qs
end

function aggregate_excess_demand(P, types; η = 0.0)
    n = size(types.c)[1]
    Z = zeros((n, 2))
    for i in 1:n
        Z[i, :] .=  individual_excess_demand(P[i, :], types.c[i, :], types.v[i]; η = η).z
    end
    return Z
end

function individual_excess_demand(p, c, v; η = 0.0)
     #println(p)
     profits = [p[1] + η - p[2] - c[1], p[2] - c[2]]
     k = argmax(profits)
     l = argmin(profits)
     z = zeros(2)
     z[k] = -1.0 .* (profits[k] .>0)
     z[l] = 0.0
     supply = -1.0 .* z
     demand = [(p[1] < v).*1.0, supply[1]]
     z[2] = z[2] + demand[2]
     z[1] = z[1] + demand[1]
     return (z= z, qs = supply, qd = demand)
 end

function eqm_price(U, types::UnobservedType; η = 0.0)
    function eqm_condition(p)
        zp = aggregate_excess_demand(transpose(p) .+ U, types; η= η)
        return norm(mean(zp, dims=1))
    end
    r = optimize(eqm_condition, [8.0, 3.0]).minimizer
    #rintln(r)
    return r
end
