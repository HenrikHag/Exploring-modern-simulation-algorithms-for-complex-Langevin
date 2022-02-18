begin
    using .UsersGuide
    using Plots
    using BenchmarkTools
    using StatsBase

end
# 
measf = "results/measuredObsHO_10_1.csv"#"results/measuredObsB100S0_7.csv"
expf = "results/expfullHO_10_1.csv"

#                                               #
#               Plotting of data                #
#                                               #

#############################
#       AutoCorrelation     #
#############################

# data = GetColumn(2,measf)
# LastRowFromFile("results/measuredObs.csv")
function PlotAC(filename,leng)
    data1 = GetData(filename,4,1)
    if leng > length(data1[:,1])
        leng = length(data1[:,1])
        println("PlotAC: Length specified to large, using length(data1[:,1]) = N_meas")
    end
    autocorrdata = transpose(StatsBase.autocor(data1,[i for i=0:leng-1]))
    jkf1 = Jackknife1(autocorrdata)
    plot(jkf1[:,1],yerr=jkf1[:,2],title="AutoCorr by StatsBase package")
end

PlotAC(measf,)
data1 = GetData(measf,4,1)
autocorrdata = Matrix{Float64}(undef,length(data1[1,:]),length(data1[:,1]))
for i = 1:length(data1[1,:])
    autocorrdata[i,:] = StatsBase.autocor(data1[:,i],[i for i=0:length(data1[:,1])-1];demean=true)
end

plot(autocorrdata[:,1],title="AutoCorr by StatsBase package")

data2 = append!(copy(data1),[0 for i=1:length(data1)])
autocorrdata2 = real.(AutoCorrR(data2))[1:length(data)]
plot(autocorrdata2,title="AutoCorr by padded data by fourier transform")




########################################
#   Probability Distribution Diagram   #
########################################

if true
    PlotProbDD(measf,0.1)
    PlotProbDDe(1,1,1,2)
end
begin
    PlotProbDD("results/measuredObsb.csv",0.1)
    PlotProbDDe(m,ω,1,2)
end




#####################################
#   Two-Point Correlation Function  #
#####################################
# begin # Jackknife estimate of error
# @benchmark GetTwoPointData(measf)
twopointD = GetTwoPointData(measf)
# @benchmark Jackknife1(twopointD)
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
    plot(jfd1[:,1],yerr=jfd1[:,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G₁(Δτ)")
end
begin # Naive estimate of error
    erd1 = Err1(twopointD1)
    plot(erd1[1:100,1],yerr=erd1[1:100,2].*300,yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G₁(Δτ)")
end

PlotTPCF(measf)
# PlotProbDD(measf,0.1)

begin
    plot_x("results/expfulla7.csv",1,[1,3,8,12])
    hline!([0])
end

jkxdat = Jackknife1(transpose(GetData("results/expfulla7.csv",3,1)))
plot(jkxdat[:,1],yerr=jkxdat[:,2])


#   ⟨x̂⟩   #
expxData = transpose(GetExpXData(expf,1))
expxDatawErr = Jackknife1(expxData)
plot(expxDatawErr[:,1],yerr=expxDatawErr[:,2])
