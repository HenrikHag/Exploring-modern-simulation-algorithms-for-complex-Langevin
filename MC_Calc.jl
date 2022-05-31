begin
    using MCMC
    # using .UsersGuide
    # using .MetropolisUpdate
    using Plots
    using BenchmarkTools
    using StatsBase
    using DataFrames
    save_date = findDate()
    save_folder = "plots/"
    # using GLM
end
# 
measf = "results/measuredObsHO_1_β8_16.csv"#"results/measuredObsB100S0_7.csv"
measf = "results/22.05.24_M_β8_16_fullAC_measuredObs.csv"
expf = "results/expfullHO_10_1.csv"

#                                               #
#               Plotting of data                #
#                                               #

mean(GetData(measf,4,2)[:,1])
mean(GetData(measf,4,1)[:,1].^2)

#############################
#       AutoCorrelation     #
#############################

# data = GetColumn(2,measf)
# LastRowFromFile("results/measuredObs.csv")


PlotAC(measf)
TPCF(measf)[end,:]
TPCF(measf,true)[end,:]


PlotAC("results/measuredObsB1.csv",100)
dat1 = GetData("results/measuredObsB1.csv",4,1)
TPCF(dat1[1,:])
ACIntegrated(dat1)

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
savefig("plots/22.04.27_M_pdd10_b16.pdf")

arr1 = GetColumn(2+16:2*16+1,measf)
arr1 = reshape(arr1,:)
histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:0.01:(floor(maximum(arr1)*10)+1)/10],normed=true,xlabel="x²",ylabel="|ψ₀|²",legend=false)

arr1 = transpose(GetColumn(2+16:2*16+1,measf))
ACfuncTPCF = Jackknife1(real.(AutoCorrR(arr1,false,false)))
plot(ACfuncTPCF[:,1],yerr=ACfuncTPCF[:,2])
PlotTPCF(measf,true, false)
begin
    TPC1=[]
    Path=arr1
    pathl=length(arr1)
    for i = 0:pathl-2       # Two-Point Correlation
        twopointcorr=0
        for ii=1:pathl
            twopointcorr += Path[ii]*Path[(ii+i-1)%pathl+1]
        end
        append!(TPC1,twopointcorr/pathl)
    end
    twopointcorr=0
    for ii=1:pathl
        twopointcorr += Path[ii]*Path[(ii+pathl-2)%pathl+1]
    end
    append!(TPC1,twopointcorr/pathl)
    TPC1
end
plot(TPC1)

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





#   Action   #
begin   # Action(xᵢ)
    save_name = "$(save_folder)$(save_date)_M_badIC_action"
    a1=[]
    a = GetData(measf,4,1)
    for i = 1:100#length(a[:,1])
        append!(a1,HO_fullAction(a[i,:],0.5,1,1))
    end
    # scatter(a)
    plt = scatter(a1,legend=false)
    display(plt)                                   # Save as png manually
    savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end



#   ⟨x̂⟩   #
# expxData = transpose(GetExpXData(expf,1))
# expxDatawErr = Jackknife1(expxData)
# plot(expxDatawErr[:,1],yerr=expxDatawErr[:,2])

begin   # ⟨x₁⟩
    save_name = "$(save_folder)$(save_date)_M_badIC_expect_x_1"
    a1=[]
    a = GetData(measf,4,1)[1:400,1]
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)
    savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end

begin   # ⟨x₁²⟩
    save_name = "$(save_folder)$(save_date)_M_expect_x2_1"
    a1=[]
    a = GetData(measf,4,1)[1:4000,1].^2
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)
    savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end
Exp_x2(16,0.95,1,1)

begin   # ⟨xᵢ⟩, ⟨xᵢ²⟩
    save_name = "$(save_folder)$(save_date)_M_expect_x_x2"
    a1 = Jackknife1(GetData(measf,4,1))
    plot(a1[:,1],yerr=a1[:,2],legend=false)
    a2 = Jackknife1(GetData(measf,4,1).^2)
    plot!(a2[:,1],yerr=a2[:,2])
    plt = hline!([Exp_x2e(16,0.5,1,1),0])
    display(plt)
    # savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    # savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end