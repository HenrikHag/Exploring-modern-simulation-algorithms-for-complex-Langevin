begin
    using .UsersGuide
    using Plots
    using BenchmarkTools
    # using StatsBase

end

measf = "results/measuredObsa7.csv"
expf = "results/expfulla7.csv"

#                                               #
#               Plotting of data                #
#                                               #

# Probability Distribution Diagram #
if true
    PlotProbDD(measf,1)
    PlotProbDDe(0.1,0.1,1,20)
end


#                                   #
#   Two-Point Correlation Function  #
#                                   #
# begin # Jackknife estimate of error
@benchmark GetTwoPointData(measf)
twopointD = GetTwoPointData(measf)
@benchmark Jackknife1(twopointD)
jfd = Jackknife1(twopointD)
    plot(jfd[1:100,1],yerr=jfd[1:100,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G(Δτ)")
# end
begin # Naive estimate of error
    erd = Err1(twopointD)
    plot(erd[1:100,1],yerr=erd[1:100,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G(Δτ)")
end

# For just the ⟨x₁xᵢ⟩

begin # Jackknife estimate of error
    twopointD1 = GetTP1data(measf)
    jfd1 = Jackknife1(twopointD1)
    plot(jfd1[1:100,1],yerr=jfd1[1:100,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G₁(Δτ)")
end
begin # Naive estimate of error
    erd1 = Err1(twopointD1)
    plot(erd1[1:100,1],yerr=erd1[1:100,2].*300,yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G₁(Δτ)")
end









#   ⟨x̂⟩   #
expxData = transpose(GetExpXData(expf,1))
expxDatawErr = Jackknife1(expxData)
plot(expxDatawErr[:,1],yerr=expxDatawErr[:,2])
