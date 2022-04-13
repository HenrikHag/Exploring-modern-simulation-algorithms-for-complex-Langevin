begin
    using .UsersGuide
    using .MetropolisUpdate
    using Plots
    using BenchmarkTools
    using StatsBase
    using DataFrames
    using GLM
end
# 
measf = "results/measuredObsHO_1_β8_16.csv"#"results/measuredObsB100S0_7.csv"
expf = "results/expfullHO_10_1.csv"

#                                               #
#               Plotting of data                #
#                                               #

#############################
#       AutoCorrelation     #
#############################

# data = GetColumn(2,measf)
# LastRowFromFile("results/measuredObs.csv")


PlotAC(measf)

"""Computes the integrated autocorrelationtime from "filename"  
"""
function ACI(filename)
    dat1 = GetData(filename,4,1)
    autocorrdata = AutoCorrR(dat1)
    auto1 = Err1(autocorrdata)[:,1]
    index1 = length(auto1)
    for i=1:length(auto1)
        if auto1[i] < 0
            index1 = i
            println("First negative value of autocorr at index: ",i)
            break
        end
    end
    auto1 = auto1[1:index1]
    acInt = sum(auto1)-auto1[1]
    # acInt *= (length(auto1)-1)/length(auto1)
    acInt += 0.5
    # acInt *= var(dat1)/length(auto1)
    return acInt#, length(dat1[:,1])
end

PlotAC("results/measuredObsB1.csv",100)
ACI("results/measuredObsB1.csv")

######
for i in ["results/measuredObsHO_1_β1_16.csv","results/measuredObsHO_1_β4_16.csv","results/measuredObsHO_1_β8_16.csv","results/measuredObsHO_1_β16_16.csv"]
    display(PlotAC(i,100000))
end
######

data1 = GetData(measf,4,1)
autocorrdata = Matrix{Float64}(undef,length(data1[1,:]),length(data1[:,1]))
for i = 1:length(data1[1,:])
    autocorrdata[i,:] = StatsBase.autocor(data1[:,i],[i for i=0:length(data1[:,1])-1];demean=true)
end
autocorrdata
plot(autocorrdata[1,:],title="AutoCorr by StatsBase package")
sbAjk = Jackknife1(autocorrdata)
plot(sbAjk[:,1],yerr=sbAjk[:,2])
sbAjk[1:5,1]

#                       #
#       AutoCorrR       #
#                       #
data2 = copy(data1)
autocorrdata2 = AutoCorrR(data2)#[1:length(data2[:,1])]
# plot(autocorrdata2[1,:])
autocorrdataJK2 = Jackknife1(autocorrdata2)
plot(autocorrdataJK2[:,1],yerr=autocorrdataJK2[:,2],title="AutoCorr by padded data by fourier transform")




########################################
#   Probability Distribution Diagram   #
########################################

begin
    PlotProbDD(measf,0.1)
    PlotProbDDe(1,1,1,2)
end
begin
    PlotProbDD("results/measuredObsb.csv",0.1)
    PlotProbDDe(m,ω,1,2)
end

######
for i in ["results/measuredObsHO_1_β1_16.csv","results/measuredObsHO_1_β4_16.csv","results/measuredObsHO_1_β8_16.csv","results/measuredObsHO_1_β16_16.csv"]
    PlotProbDD(i,0.1)
    display(PlotProbDDe(1,1,1,2))
end
######

for i in [i for i = -1:0.1:1]
    # println(i)
    println(cosh(i))
end


#####################################
#   Two-Point Correlation Function  #
#####################################
# Naive error
PlotTPCF(measf)
PlotTPCF(measf)

# Jackknife error
PlotTPCF(measf,true)
PlotTPCF(measf,true)

# Effective Mass(Δτ)
PlotEffM(measf)
PlotEffM(measf,true)



# begin # Jackknife estimate of error
# @benchmark GetTwoPointData(measf)
twopointD = GetTwoPointData(measf)
# @benchmark Jackknife1(twopointD)
jfd = Jackknife1(twopointD)
plot(jfd[:,1],yerr=jfd[:,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G(Δτ)")
# end
begin # Naive estimate of error
    erd = Err1(twopointD)
    plot(erd[:,1],yerr=erd[:,2],yrange=[1.4*10^-3,10^2],yaxis=:log,xlabel="Δτ",ylabel="G(Δτ)")
end

for i in ["results/measuredObsHO_1_β1_16.csv","results/measuredObsHO_1_β4_16.csv","results/measuredObsHO_1_β8_16.csv","results/measuredObsHO_1_β16_16.csv"]
    PlotTPCF(i)
end
PlotTPCF("results/measuredObsHO_1_β8_16.csv")


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

