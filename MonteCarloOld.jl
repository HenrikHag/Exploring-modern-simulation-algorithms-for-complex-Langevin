
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
    rng = MersenneTwister(11111)  # Reseed the generator
    gaussianD = Normal(0.16, 1.5)
    save_path = "results/"
    save_date = findDate()
    save_pre = "$(save_path)$(save_date)_M_"
    # n=rand(d,1000)
end


##############################################################
# 
# This program will run a Monte Carlo Simulation, simulating a system of different parameters.
# 
# The Action:   The action of the system simulated. Here the physics of the system is defined.
#   Implemented:    HO_Action(), AHO_Action()
# 
# n_tau:        The number of parts in the timelattice
# iδτ:          The index of Euclidean time, i ∈ [1,...,n_tau]
# δτ:           Lattice spacing
# n_burn:       The number of configurations to thrash at start of simulations
# n_skip - 1:   The number of configurations to thrash between each measurement
# m̃:            Dimentionless effective mass,    ̃m = m δτ
# ̃ω:            Dimentionless frequency,         ̃ω = ω δτ
# 
# 
##############################################################



#                   #
#       TODO        #
#                   #

# VIII. SUGGESTED PROBLEMS
#   - b, c, d, e, f, g
#


# MetroSwipe(16,1,1,0,0.5,1,0.8,rng,zeros(16))

#########################
# Metropolis Simulation #
#########################

phys_param = getAHO_M_param(16,8,1,1,0)
a = phys_param.a
sim_param = getSim_M_param(1000,200,20,0.8,rng)
MetropolisSim(phys_param, sim_param, save_pre, "trash")

PlotAC("$(save_pre)shortSim_obs.csv")

####################### main end ##########################################


n_tau=16;m=1;β=8;a=β/n_tau;λ=0;
# n_skip=0; n_burn=2500; n_total=n_burn + 1 + (1+n_skip)*15000;
main(n_tau,"$(save_date)_M_β$(β)_$(n_tau)_fullAC_expfull.csv","$(save_date)_M_β$(β)_$(n_tau)_fullAC_measuredObs.csv","$(save_date)_M_β$(β)_$(n_tau)_fullAC_expect.csv",m,m,a,λ)

main(n_tau,"22.05.12_M_β$(β)_$(n_tau)_fullAC_expfull.csv","22.05.12_M_β$(β)_$(n_tau)_fullAC_measuredObs.csv","22.05.12_M_β$(β)_$(n_tau)_fullAC_expect.csv",m,m,a,λ)

analysisFile = "results/22.05.12_M_β8_16_fullAC_measuredObs.csv"
PlotAC(analysisFile,300)
# title!("Autocorrelation")
savefig("plots/22.05.12_M_β8_16_fullAC.pdf")
savefig("plots/22.05.12_M_β8_16_fullAC.png")
PlotACsb(analysisFile,300)
# title!("Autocorrelation StatsBase")
savefig("plots/22.05.12_M_β8_16_fullACsb.pdf")
savefig("plots/22.05.12_M_β8_16_fullACsb.png")
PlotTPCF(analysisFile,true)
# title!("Two-Point Correlation")
savefig("plots/22.05.12_M_β8_16_fullAC_TPCF.pdf")
savefig("plots/22.05.12_M_β8_16_fullAC_TPCF.png")
PlotEffM(analysisFile)
# title!("Effective Mass")




configs1 = [1 120; 0.8 150; 0.6 200; 0.5 240; 0.3 400; 0.2 600; 0.1 1200; 0.08 1500; 0.06 2000; 0.05 2400; 0.03 4000; 0.02 6000; 0.01 12000]
configs1 = [1 1 16; 1 4 16; 1 8 16; 1 16 16]
for i=1:4
    m = configs1[i,1]
    β = configs1[i,2]
    n_tau = configs1[i,3]
    a = β/n_tau
    λ = 0
    main(n_tau,"expfullHO_1_β$(β)_16.csv","measuredObsHO_1_β$(β)_16.csv","expectationvalHO_1_β$(β)_16.csv",m,m,a,λ)
end
for i=1:7
    m = configs1[i,1]
    n_tau = configs1[i,2]
    main(n_tau,"expfullAHO$(i).csv","measuredObsAHO$(i).csv","expectationvalAHO$(i).csv",m,m,a,0)
end


# configs1[1,2]
main(16,"expfullB1.csv","measuredObsB1.csv","expectationvalB1.csv",1,1,0.5,0)
twopoint=[]
for i = (120*3+2):(120*4+1)
    append!(twopoint,mean(GetColumn(i,"results/measuredObsB1.csv")))
end
# twopoint
plot(twopoint)
# GetLastMean("results/measuredObsB1.csv",120*4)[120*3+1:120*4]
plot(GetLastMean("results/measuredObsB1.csv",120*4)[120*3+1:120*4])

a2=[]
for i=7:7#length(configs1[:,1])-5
    m = configs1[i,1]
    n_tau = configs1[i,2]
    push!(a2,main(n_tau,"expfulla$(i).csv","measuredObsa$(i).csv","expectationvala$(i).csv",m,m))
end
a2
twopointD = GetTwoPointData("results/measuredObsa7.csv")
jfd=Jackknife1(twopointD)
println(jfd[1:100])
jfd=a2[7]
jfd0 = Matrix{Float64}(undef,100,2)
for i=1:100
    if jfd[i,1]<=0.01
        jfd0[i,:]=jfd0[i-1,:]
    else
        jfd0[i,:]=jfd[i,:]
    end
end
plot(jfd0[1:100,1],yerr=jfd0[1:100,2],yaxis=:log,xlabel="Δτ",ylabel="G(Δτ)")
plot(jfd[:,1],yerr=jfd1[:,2],xlabel="Δτ")

function EstEffM(Gm1,Gp1)
    return log10(Gm1/Gp1)/2
end
EstEffM(jfd[30,1],jfd[32,1])

plot(GetColumn(4,"results/expectationvala$(num1).csv"))
begin
    num1 = 1
    plot(GetColumn(2,"results/expectationvala$(num1).csv"))    # ⟨x̂⟩
    hline!([0])
end
begin
    plot(GetColumn(3,"results/expectationvala$(num1).csv"))# ⟨x̂²⟩
    expected_x2 = Exp_x2(configs1[num1,2], configs1[num1,1], configs1[num1,1])
    hline!([expected_x2])
end
expected_x4 = 3*Exp_x2(120, m, ω)^2






###### Plot ⟨xᵢ⟩, ⟨xᵢ²⟩
plotexpx(1,4,meanf)         # corr1
plotexpx1(2,true)           # exp_x_1
plotexpx1(2,false)          # exp_x_1
plotexpx2(2,n_tau,false)    # exp_x_2
plotexpx(n_tau+1, 20)       # corr1


######### Plot final ⟨xᵢ⟩, and ⟨xᵢ²⟩ #########
begin
    lastRow = LastRowFromFile("results/expfullB100S0_1.csv")
    exp3 = [lastRow[i] for i=2:n_tau+1]
    PlotExp(exp3,3)         # ⟨x̂ᵢ⟩  (Sim)
    exp3 = [lastRow[i+n_tau] for i=2:n_tau+1]
    PlotExp(exp3,3)         # ⟨x̂ᵢ²⟩ (Sim)
    hline!([0])             # ⟨x̂⟩ₐ
    expected_x2 = Exp_x2(120, 1, 1)
    hline!([expected_x2/2]) # ⟨x̂²⟩ₐ
end




######### Plot ⟨x⟩ #########
n_tau
begin
    plot_x("results/expfulla7.csv",1,[20,40,60,80])#[5,12,22,32,42,52])#0)#[i for i=1:n_tau])
    plot!([0 for i=1:length(readlines("results/expfulla7.csv"))], title="Expectationvalue of x", label="⟨x⟩ₜₕₑ",legend=:bottom)
end




######### Big Simulation ###################
Path = zeros(n_tau)
sum1 = zeros(n_tau)
sum2 = zeros(n_tau)
sum3 = zeros(n_tau)
runtotal = 40
for i=1:runtotal
    main(n_tau,"expfull$(i).csv","measuredObs$(i).csv")
    Path = zeros(n_tau)
    println("Simulations run: (",i,"/",runtotal,")")
    # show(IOContext(stdout, :limit => true),"text/plain",Path)
    # println()
end

begin
    errScale = 2
    ter = 20
    a = []
    for i=1:runtotal
        # push!(a,[1,2,3,i])
        push!(a,AutoC(2,ter,"results/expfull$(i).csv","results/measuredObs$(i).csv")[:,1])
    end
    mean1 = zeros(Float64,ter-1,1)
    for ii = 1:ter-1
        for i=1:runtotal
            mean1[ii] += a[i][ii]/runtotal
        end
    end
    variance = zeros(Float64,ter-1,1)
    err = zeros(Float64,ter-1,1)
    for ii = 1:ter-1
        for i=1:runtotal
            variance[ii] += (a[i][ii]-mean1[ii])^2
        end
        variance[ii] /= runtotal-1
        err[ii] = errScale*√(variance[ii]/runtotal)
    end
end
a[1]
mean1
plot(mean1,yerr=err,title="Autocorrelation many runs")












# Create a "path of a gaussian random coordinate" to perform checks against
function GaussianPath(numb, rng)
    path = zeros(numb)
    randn!(rng, path)
    return path
end

path = GaussianPath(5000, rng)
begin
    for i=2:length(path)
        path[i] += path[i-1]
    end
    histogram(path)
    
end

plot(AutoCorrR(path))

path[1:4999-20]
mean(path[1:5000])
tau_MC = 5
N_meas = 4999
a = [path[i] for i=1:N_meas-tau_MC]
mean(a)



column = GetColumn(1,"results/measuredObs1.csv")
plot(AutoCorrR(column))


plot(AutoCorrR(GetColumn(2,"results/measuredObsG2.csv")))
AutoCorrR(rand(gaussianD,200))

plot(AutoCorrR(rand(gaussianD,2000)))



# printarray([1.1,-1.02,1.50,-0.4])
# printarray([-1.02,1.50,-0.4,1.1])