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
measf = "saved_results/22.06.07_M_shortSim_obs.csv"
measf = "saved_results/22.06.07_M_shortSim_HighAC_obs.csv"
measf = "results/22.06.07_M_midSim_HighAC_obs.csv"
measf = "results/22.06.08_M_midSim_obs.csv"
measf = "results/22.06.08_M_midSim_b120_obs.csv"
measf = "results/22.06.10_M_longSim_100k_1skip_obs.csv"
measf = "results/22.06.10_M_longSim_10k_1skip_obs.csv"
measf = "results/22.06.10_M_longSim_10k_20skip_obs.csv"
measf = "results/22.06.10_M_longSim_10k_40skip_obs.csv"
expf = "results/expfullHO_10_1.csv"

# ⟨x⟩ and ⟨x²⟩
measf = "results/22.06.12_M_longSim_1M_1skip_obs.csv"
measf = "results/22.06.12_M_longSim_10k_40skip_obs.csv"

# PDD
measf = "results/22.06.12_M_longSim_1M_160skip_b2_obs.csv"
measf = "results/22.06.12_M_longSim_1M_40skip_obs.csv"
measf = "results/22.06.12_M_longSim_1M_40skip_b32_obs.csv"

# Plot for enough samples in results and Jackknife
# PDD for β = 2 and 8



#                                               #
#               Plotting of data                #
#                                               #
PlotAC(measf,200)
# PlotACsb(measf)
data1 = GetData(measf,4,1)#[1:200,:]
data1 = [1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1]
@benchmark Autocorrelation_BySummation(data1)
@benchmark AutoCorrR(data1)
PlotAC_BySummation(data1,false)
plot(StatsBase.autocor(data1,[i for i=0:length(data1)-1]))

PlotAC_BySummation(data1,200)
PlotAC(data1,false)


mean(GetData(measf,4,2)[:,1])
mean(GetData(measf,4,1)[:,1].^2)

#############################
#       AutoCorrelation     #
#############################

# data = GetColumn(2,measf)
# LastRowFromFile("results/measuredObs.csv")


PlotAC(measf,300)
# title!("Autocorrelation")
savefig("$(save_folder)$(save_date)_M_1M_1skip_300AC.pdf")
savefig("$(save_folder)$(save_date)_M_1M_1skip_300AC.png")

#############################
#   Two-Point correlation   #
#############################

PlotTPCF(measf,true,true)
PlotTPCFe!(0.5,1,1,16)
# title!("Two-Point Correlation")
savefig("$(save_folder)$(save_date)_M_10k_40skip_TPCF.pdf")
savefig("$(save_folder)$(save_date)_M_10k_40skip_TPCF.png")

TPCF(measf,true)[end,:]

PlotEffM(measf,16)
# title!("Effective mass")
savefig("$(save_folder)$(save_date)_M_10k_40skip_EffM.pdf")
savefig("$(save_folder)$(save_date)_M_10k_40skip_EffM.png")

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
    save_name = "$(save_folder)$(save_date)_M_1M_40skip_b32"
    plt = PlotProbDD(measf)
    PlotProbDDe!(plt,1,1,1,4)
    display(plt)
    savefig(plt,"$(save_name)_PDD.pdf")
    savefig(plt,"$(save_name)_PDD.png")
end

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
    save_name = "$(save_folder)$(save_date)_M_longSim_10k_40skip"
    a1=[]
    a = GetData(measf,4,1)[1:1000,1]
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    plt = scatter(a,xlabel="t",ylabel="⟨x⟩",label="x₁")
    plot!(plt,a1,width=4,label="⟨x₁⟩")
    hline!(plt,[0],label="⟨x⟩ₜₕₑₒ")
    display(plt)
    savefig(plt,"$(save_name)_x_1.pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name)_x_1.png")   # Save as png in folder "plots"
end

begin   # ⟨x₁²⟩
    save_name = "$(save_folder)$(save_date)_M_longSim_10k_40skip"
    a1=[]
    a = GetData(measf,4,1)[1:1000,1].^2
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    plt = scatter(a,xlabel="t",ylabel="⟨x²⟩",label="x₁²")
    plot!(plt,a1,width=4,label="⟨x₁²⟩")
    hline!(plt,[Exp_x2e(16,0.5,1,1)],label="⟨x²⟩ₜₕₑₒ continuum")
    display(plt)
    savefig(plt,"$(save_name)_x2_1.pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name)_x2_1.png")   # Save as png in folder "plots"
end
Exp_x2(16,0.95,1,1)

begin   # ⟨xᵢ⟩, ⟨xᵢ²⟩
    save_name = "$(save_folder)$(save_date)_M_longSim_10k_40skip_er"
    data1 = GetData(measf,4,1)
    # a1 = Jackknife1(data1,true)
    a1 = Err1(data1)
    plt = plot(a1[:,1],yerr=a1[:,2],legend=:right,xlabel="τ",ylabel="⟨O⟩",label="⟨x(τ)⟩")
    # a2 = Jackknife1(data1.^2,true)
    a2 = Err1(data1.^2)
    plot!(plt,a2[:,1],yerr=a2[:,2],label="⟨x²(τ)⟩")
    hline!(plt,[[Exp_x2e(16,0.5,1,1)],[Exp_x2(16,0.5,1,1)],[0]],label=["⟨xᵢ²⟩ₜₕₑₒ continuum" "⟨xᵢ²⟩ₜₕₑₒ discretized" ""],color=["black" "green" "black"])
    display(plt)
    savefig(plt,"$(save_name)_x_x2.pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name)_x_x2.png")   # Save as png in folder "plots"
end






begin # TPCF
    save_name = "$(save_folder)$(save_date)_M_1M_40skip_b8"
    plt = PlotTPCF(measf)
    display(plt)
    # savefig(plt,"$(save_name)_TPCF.pdf")   # Save as pdf in folder "plots"
    # savefig(plt,"$(save_name)_TPCF.png")   # Save as p
end

begin # Effective mass
    save_name = "$(save_folder)$(save_date)_M_shortSim"
    plt = PlotEffM(measf)
    display(plt)
    # savefig(plt,"$(save_name)_EffM.pdf")   # Save as pdf in folder "plots"
    # savefig(plt,"$(save_name)_EffM.png")   # Save as p
end