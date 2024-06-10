using RCall, Roots

function compute_rule(Y, Z, X, W, μw, Xpred)
    #use 80\% to train and 20\% to choose c
    samp = rand(Bernoulli(0.8), length(Y))
    hattauY, tauYpred = cforest(Y[samp.==1], X[samp.==1, :], W[samp.==1], Xpred, X)
    hattauZ, tauZpred = cforest(Z[samp.==1], X[samp.==1, :], W[samp.==1], Xpred, X)


    function demand_at_c(c)
        S = hattauY .> hattauZ *c
        compute_gain(Z[samp.==0], W[samp.==0], S[samp.==0], μw)
    end

    if demand_at_c(0.0)*demand_at_c(10.0) > 0
        #println("here")
        cstar = 0
    else
        cstar =  find_zero(demand_at_c,(0.0, 10.0))
    end
    Snaive = hattauY .> 0
    αstar = get_naive_alpha(Z[samp.==0], W[samp.==0], Snaive[samp .==0] , μw)
    #cstar = 0.7
    #println("cstar: ", cstar)
    #println("demand: ", demand_at_c(cstar))
    Sstarpred = tauYpred .> tauZpred*cstar
    #println(mean(Sstar[X[:, 12].==0]))
    #println(mean(Sstar[X[:, 12] .> 0]))
    #println("lift:", compute_gain(Y, W, Sstarpred, μw))
    return cstar, hattauY, hattauZ, tauYpred, tauZpred, αstar
end

function rule_value(Y1, Y0, S)
    return mean(Y1.* S) + mean(Y0  .* (1 .- S))
end


function get_naive_alpha(Z, W, S, μw)
    Z1 = W .*Z ./ mean(W)
    Z0 = (1 .- W).*Z ./ (1 - mean(W))


    function solve_a(a)
        squo = mean(Z1)*μw + mean(Z0)*(1- μw)
        rule = mean(Z1 .* S)*a + mean(Z0 .* S)*(1 - a) + mean(Z0 .* (1 .- S))
        return squo - rule
    end

    if solve_a(0.0)* solve_a(1.0) > 0
        return 1
    else
        return find_zero(solve_a,(0.0, 1.0))
    end

end

function compute_gain(Y, W, S, μw)
    scoreLeft = (W .*Y ./ mean(W) .- (1 .- W).*Y ./ (1 - mean(W))) .* S
    scoreRight = (W .* Y ./ mean(W) .- (1 .- W) .* Y ./ (1 - mean(W)))
    #println(mean(scoreRight))
    #println(mean(scoreRight*μw))
    #println(mean(scoreLeft))
    return mean(scoreLeft) - mean(scoreRight)*μw
end

function cforest(Y, X, W, Xpred, Xpred2)
    tauhat = R"""
        set.seed(1)
        eval.forest <- grf::causal_forest(data.frame($X), as.numeric($Y), ci.group.size =1, honesty=F, num.trees=2000, as.numeric($W))
        list(predict(eval.forest, $Xpred2)[,1], predict(eval.forest, $Xpred)[,1])
    """
    return rcopy(tauhat)
end
