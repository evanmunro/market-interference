using Plots, Random
using CSV
using DataFrames
Random.seed!(125)
data = DataFrame(CSV.File("hattau.csv"))

samplesub = rand(1:length(data.Z), 125)

c = 3.726081
W = data.Y[samplesub] .> c .* data.Z[samplesub]
plot(data.Z[samplesub][W .== 1], data.Y[samplesub][W .==1], seriestype=:scatter, framestyle=:origin, xlabel="ZCADEᵢ",
        ylabel="CADEᵢ", label="Treated")
plot!(data.Z[samplesub][W .== 0], data.Y[samplesub][W .==0], seriestype=:scatter, framestyle=:origin, label="Control")
xs = Vector(-0.75:0.01:maximum(data.Z))
ys = xs.*c
plot!(xs, ys, label="Treatment Rule", legend=:topleft, linewidth=2.0, ylim = [-10, 10], xlim = [-0.75, 1.0])
savefig("fig2.pdf")
