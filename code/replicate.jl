
includet("php_household_model.jl")


#Figure 1a
Random.seed!(1)
supply_plot()

#Figure 1b, 1-2 mins
Random.seed!(1)
ade_plot(5000, 1000)

#Figure 2
Random.seed!(1)
hte_plot(5000)

#Table 1
Random.seed!(1)
coverage_simulation(2000, 1000)

#Table 2
Random.seed!(1)
hte_gain(50)

#Table 3
Random.seed(1)
generate_model_estimation_table()
