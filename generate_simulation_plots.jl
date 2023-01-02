includet("freelance_simulation.jl")
Random.seed!(1234)

#compute ground truths to verify line 174 of freelance_simulation.jl
compute_gts(UnobservedType(Int(1e7)), true)

#Figure 1b
supply_plots(10000)
#Figure 1a
ade_plot(2000, 10000)



#Components of Table 2 printed below
coverage_simulation(2000, 10000)
#Table 1
