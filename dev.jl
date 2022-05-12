# Development, testing new functions

begin
    using MCMC
    # using .MetropolisUpdate
    # using .UsersGuide
    using Random
    using Plots
    using Statistics
    using Distributions
    using FFTW
    using StatsBase
    using Base.Threads
    using BenchmarkTools
    
    # using Printf
    # using DataFrames
    # using CSV
    
    const obsf = "results/measuredObs.csv"
    const meanf = "results/expfull.csv"
    
    # n_tau, m, ω
    rng = MersenneTwister(11111)
    gaussianD = Normal(0.16, 1.5)
    # n=rand(d,1000)
end




"""Does a metroswipe by testing a change for each element,  
and adds this new coord weighted by change in action.  
MultiThreaded"""
function MetroSwipeMT(n_tau, m, ω, h, idrate, rng, Path)
    accept = 0
    Threads.@threads for i = 1:2:n_tau
        x_new = Path[i] + h*2*(rand(rng)-1/2)
        s_old = HO_Action(n_tau, m, ω, Path, i, Path[i])
        s_new = HO_Action(n_tau, m, ω, Path, i, x_new)
        # printf("s_old: %f, s_new: %f\n", s_old, s_new)
        if rand(rng) < exp(s_old-s_new)
            Path[i] = x_new
            Threads.atomic_add!(accept,1/n_tau)
        end
    end
    Threads.@threads for i = 2:2:n_tau
        x_new = Path[i] + h*2*(rand(rng)-1/2)
        s_old = HO_Action(n_tau, m, ω, Path, i, Path[i])
        s_new = HO_Action(n_tau, m, ω, Path, i, x_new)
        # printf("s_old: %f, s_new: %f\n", s_old, s_new)
        if rand(rng) < exp(s_old-s_new)
            Path[i] = x_new
            Threads.atomic_add!(accept,1/n_tau)
        end
    end
    # rand!(gaussianD, Path)
    # randn!(rng, Path)
    h *= accept / idrate
    # println(h)
    return Path, accept, h
end
# Path = zeros(n_tau);
# rand!(gaussianD, Path)
# Path


#                               #
# Defining a main function      #
#                               #
function SimMetro(n_tau,meanfname,obsfname)
    Path = zeros(n_tau)
    sum1 = zeros(n_tau); sum2 = zeros(n_tau); sum3 = zeros(n_tau)
    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    # sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    n_burn = 200
    n_skip = 20
    n_total = n_burn + 1 + (n_skip+1)*10000#20000#200100
    accept = 0
    idrate = 0.8
    h = 1
    #                               #
    # Prepare writing to file       #
    #                               #
    rfold = "results/"
    touch(string(rfold,meanfname))
    touch(string(rfold,obsfname))

    # show(IOContext(stdout, :limit => true),"text/plain",Path)
    # println()
    itt = 0
    Randlist = []
    #                                                       #
    # Do n_total swipes, removing n_burn first swipes       #
    #                                                       #
    a = @timed if n_burn != 0
        for i = 1:n_burn
            Path, accept, h = MetroSwipe(n_tau, m, ω, h, idrate, rng, Path)
        end
    end
    println("Burn-in $(n_burn) complete!")
    b = @timed for i = 1:n_total-n_burn
        Path, accept, h = MetroSwipe(n_tau, m, ω, h, idrate, rng, Path)
        if n_skip == 0 || (i-1-n_burn)%n_skip == 0
            sum1, sum2, sum3 = MeasureObs(n_tau, sum1, sum2, sum3, Path)
            append!(Randlist,h)
            itt += 1
            # println("Measured ",itt)
            # itt = (i-1-n_burn)/(n_skip)+1
            exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
            writee123tofile(n_tau,rfold,meanfname, exp_x, exp_x2, exp_x0x1, itt)
            writec123tofile(rfold,obsfname, Path, itt)
        end
    end
    # exp_x /= (n_total-n_burn)/n_skip
    #exp_x, exp_x2, exp_x0x1 += MetroUpdate(n)
    # println("Mean h = ", Statistics.mean(Randlist))
    println("Time: ",a.time, " ", b.time)
    return (a.time, b.time)
end

