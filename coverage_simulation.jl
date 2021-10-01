using StatsPlots, Revise
include("binary_simulation.jl")
#coverage results printed here
results = run_simulation(2500, 10000)
#also can plot smoothed empirical distribution versus actual
density(results[:, 1:2], label = ["Direct Effect" "Indirect Effect"])
μ = mean(results, dims=1)
println("Estimated: ", μ[3], " Indirect: ", μ[4])
println("Simulated: ", std(results[:, 1]), " ", std(results[:, 2]))
plot!(Normal(μ[1], μ[3]), label="Approx Direct")
plot!(Normal(μ[2], μ[4]), label="Approx Indirect")
#plot(values, seriestype=:scatterhist, linestyle=:solid, size=(600,150))
