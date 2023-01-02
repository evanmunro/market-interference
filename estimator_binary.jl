using Statistics, FixedEffectModels, DataFrames

struct EqmData
    Y::Array{Float64}
    D::Array{Float64}
    S::Array{Float64}
    Z::Array{Float64}
    W::Array{Bool, 1}
    U::Array{Float64}
    p::Float64
end

function EqmData(Y, D, S, W, U; p=0.5)
    return EqmData(Y, D, S, D .- S, W, U, p)
end

struct Partials
    yp1::Float64
    yp0::Float64
    dp1::Float64
    dp0::Float64
    sp1::Float64
    sp0::Float64
    μzp::Float64
end

function Partials(Y, D, S, W, U; p=0.5)
    ∂yp1, ∂yp0 = deriv_estimate(Y, W, U)
    ∂dp1, ∂dp0 = deriv_estimate(D, W, U)
    ∂sp1, ∂sp0 = deriv_estimate(S, W, U)
    ∂μzp = p*(∂dp1 - ∂sp1) + (1 - p)*(∂dp0 - ∂sp0)
    return Partials(∂yp1, ∂yp0, ∂dp1, ∂dp0, ∂sp1, ∂sp0, ∂μzp)
end

struct DiffMeans
    y::Float64
    d::Float64
    s::Float64
    z::Float64
end

function DiffMeans(Y, D, S, W; p=0.5)
    y = horvitz_thompson(Y, W)
    d = horvitz_thompson(D, W)
    s = horvitz_thompson(S, W)
    DiffMeans(y, d, s, d - s)
end

function horvitz_thompson(Y, W; p=0.5)
    data = (W./p - (1 .-W)./(1 - p)).*Y
    return mean(data)
end

function deriv_estimate(Y, W, U)
    ∂y1 = cov(Y[W.==1], U[W.==1])/var(U[W.==1])
    ∂y0 = cov(Y[W.==0], U[W.==0])/var(U[W.==0])
    return (∂y1, ∂y0)
end

function varianceDirect(Y, Z, W, ∂::Partials; p=0.5)
    a = (∂.yp1 - ∂.yp0)
    H = ∂.μzp
    A = -1 .*a./H.*Z
    #println("diff ", a)
    #println("H: ", H)
    #println("Cov", cov(Y, Z))
    data = (W./p - (1 .-W)./(1 - p)).*(Y) .+ A
    #println(std(A)/sqrt(length(Y)))
    #println()
    return std(data)/sqrt(length(Y)), std(data .- A)/sqrt(length(Y))
end

function estimate_τ(data::EqmData)
    ∂ = Partials(data.Y, data.D, data.S, data.W, data.U)
    Δ = DiffMeans(data.Y, data.D, data.S, data.W)
    τD = Δ.y
    #println((data.p*∂.yp1 + (1- data.p)*∂.yp0))
    τI = -1.0*(data.p*∂.yp1 + (1- data.p)*∂.yp0)*Δ.z/∂.μzp
    σD, σRCT = varianceDirect(data.Y, data.Z, data.W, ∂)
    σI = varianceIndirect(data.Y, data.Z, data.W, data.U, ∂, Δ)
    return (τD = τD, σD = σD, τI = τI, σI = σI, σRCT = σRCT)
end

function ivCOV(Y, X, Z)
    df = DataFrame(Y=Y, X=X, Z=Z)
    regmodel = reg(df, @formula(Y ~ (X ~ Z)), Vcov.robust())
    return (coef(regmodel)[2], stderror(regmodel)[2], coef(regmodel)[1])
end
function varianceIndirect(Y, Z, W, U, ∂::Partials, Δ::DiffMeans; p = 0.5)
    dpstar = -1.0*Δ.z
    coefIV, covIV, int = ivCOV(Y, Z, U)
    #a = (∂.yp1 - ∂.yp0)
    b = (∂.dp1 - ∂.sp1) - (∂.dp0 - ∂.sp0)
    H = ∂.μzp
    B = (W./p - (1 .-W)./(1 - p)).*Z .- b/H.*Z
    covZ = var(B)/length(B)
    #σ2 = dpstar*(covIV^2*dpstar) + coefIV^2*covZ
    #println("Cov direct:", sqrt(dpstar*(covIV^2*dpstar)))
    ehat = Y .- Z.* coefIV .- int
    estimatehat = mean(ehat.^2)*(Δ.z/∂.μzp)^2/(length(Y)*U[1]^2)
    σ2 = estimatehat + coefIV^2*covZ
    return sqrt(σ2)#/abs(U[1])
end
