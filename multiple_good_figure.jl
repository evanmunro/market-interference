
using StatsPlots, Revise
includet("multiple_good_market.jl")

τ_truth = 0.172 #0global_estimate()
simulation = mc_simulation(10000, 1000)
println("Mean: ", mean(simulation; dims=1))
println("SD: ", std(simulation; dims = 1))
density(simulation[:, 1:2], label=["Estimated Marginal Policy Effect" "Estimated Direct Effect"], linestyle = [:solid :dash])
vline!([τ_truth], label= "True Marginal Policy Effect", linestyle = :dot, linecolor=:black, linewidth=2.0)
savefig("simulation.pdf")
