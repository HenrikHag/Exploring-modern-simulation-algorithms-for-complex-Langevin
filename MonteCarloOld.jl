# include("MetropolisUpdate.jl")
using Revise
begin
    using .MetropolisUpdate
    using .UsersGuide
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
    # expected_x2 = Exp_x2(n_tau, m, ω)
    rng = MersenneTwister(11111)
    gaussianD = Normal(0.16, 1.5)
    # n=rand(d,1000)
end


### Probability density diagram ###
function PlotProbDD(file)
    arr1 = Vector{Float64}(undef,0)
    for i=2:Int((length(LastRowFromFile(file))-1)/3)+1
        append!(arr1,GetColumn(i,file))
    end
    histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:0.1:(floor(maximum(arr1)*10)+1)/10],normed=true)#,weights=repeat([1/length(arr1)],length(arr1))
end

function PlotProbDDe(m,ω,ħ)
    println("m = ",m,", ω = ",ω,", ħ = ",ħ)
    plot!([x for x=-2:0.01:2],[((m*ω/(π*ħ))^(1/4)*exp(-m*ω*x^2/(2*ħ)))^2 for x=-2:0.01:2],linewidth=2)
end

begin
    PlotProbDD("results/measuredObsb.csv")
    PlotProbDDe(m,ω,1)
end


# TODO:
# VIII. SUGGESTED PROBLEMS
    # b, c, d, e, f, g



#                                   #
# Measure the current observables   #
#                                   #
"""Adds Path coords to the sums  
```julia
function MeasureObs()
    for i = 1:n_tau 
        sum1[i] += Path[i]
        sum2[i] += Path[i]^2
        sum3[i] += Path[1]*Path[i]
```"""
function MeasureObs(n_tau, sum1, sum2, sum3, Path)
    for i = 1:n_tau
        sum1[i] += Path[i]
        sum2[i] += Path[i]^2
        sum3[i] += Path[1]*Path[i]
    end
    return sum1, sum2, sum3
end

"""Append the mean of exp_x, exp_x2, exp_x0x1 on one line to a file  
"""
function writeeMean(path,filename,exp_x,exp_x2,exp_x0x1,itt)
    open(string(path,filename),"a") do file
        write(file,string(Int64(itt),","))
        write(file,string(mean(exp_x),","))
        write(file,string(mean(exp_x2),","))
        write(file,string(mean(exp_x0x1),"\n"))
    end
end



#                               #
# Defining a metropolis swipe   #
#                               #
"""Does a metroswipe by testing a change for each element,  
and adds this new coord weighted by change in action."""
function MetroSwipe(n_tau, m, ω, h, idrate, rng, Path)
    accept = 0
    for i = 1:n_tau
        x_new = Path[i] + h*2*(rand(rng)-1/2)
        s_old = HO_Action(n_tau, m, ω, Path, i, Path[i])
        s_new = HO_Action(n_tau, m, ω, Path, i, x_new)
        # printf("s_old: %f, s_new: %f\n", s_old, s_new)
        if rand(rng) < exp(s_old-s_new)
            Path[i] = x_new
            accept += 1/n_tau
        end
    end
    # rand!(gaussianD, Path)
    # randn!(rng, Path)
    h *= accept / idrate
    # println(h)
    return Path, accept, h
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


function main(n_tau,meanfname,obsfname,expfname,m,ω)
    # Logic of program should go here
    # We want to print final expectation values
    # rng = MersenneTwister(11111)
    n_tau = Int(n_tau)
    println("n_tau = ",n_tau,", m,ω = ", m,", ",ω)
    Path = zeros(n_tau)
    # sum1 = zeros(n_tau); sum2 = zeros(n_tau); sum3 = zeros(n_tau)
    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    n_burn = 100
    n_skip = n_tau/10#12
    n_total = n_burn + 1 + (n_skip)*10000#20000#200100
    accept = 0
    idrate = 0.8
    h = 1
    #                               #
    # Prepare writing to file       #
    #                               #
    rfold = "results/"
    touch(string(rfold,meanfname))
    touch(string(rfold,obsfname))
    touch(string(rfold,expfname))

    # show(IOContext(stdout, :limit => true),"text/plain",Path)
    # println()
    itt = 0
    Randlist = []
    #                                                       #
    # Do n_total swipes, removing n_burn first swipes       #
    #                                                       #
    a = @timed for i = 1:n_total
        Path, accept, h = MetroSwipe(n_tau, m, ω, h, idrate, rng, Path)
        if i > n_burn
            if n_skip == 0 || (i-1-n_burn)%n_skip == 0
                sum1, sum2, sum3 = MeasureObs(n_tau, sum1, sum2, sum3, Path)
                append!(Randlist,h)
                itt += 1
                # println("Measured ",itt)
                # itt = (i-1-n_burn)/(n_skip)+1
                exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
                # print("Path      "); printarray(Path); print("<s_i>     "); printarray(sum1)
                
                # print("<x_i>     "); printarray(exp_x)
                # print("<(x_i)^2> "); printarray(exp_x2)
                # print("<x_1*x_i> "); printarray(exp_x0x1)
                
                # exp1, exp2, exp3      // Append each itt
                # writee123tofile(n_tau,rfold,meanfname, exp_x, exp_x2, exp_x0x1, itt)
                # curr1, curr2, curr3   // Append each itt
                # writec123tofile(rfold,obsfname, Path, itt)
                # curr3, sum3           // New file
                #filenamec3s3 = string("curr3sum3n",Int64(itt),".csv")
                #writec3s3tofile(rfold,filenamec3s3, sum3)
                # sum3, exp3            // New file
                #filenames3e3 = string("sum3exp3n",Int64(itt),".csv")
                #writes3e3tofile(rfold,filenames3e3, sum3, exp_x0x1)
                # ⟨exp1⟩, ⟨exp2⟩, ⟨exp3⟩  // Append each itt
                writeeMean(rfold,expfname,exp_x,exp_x2,exp_x0x1,itt)
            end
        end
    end
    # exp_x /= (n_total-n_burn)/n_skip
    #exp_x, exp_x2, exp_x0x1 += MetroUpdate(n)
    # println("Mean h = ", Statistics.mean(Randlist))
    println("Time: ",a.time)
    return a.time
end

configs1 = [1 120; 0.8 150; 0.6 200; 0.5 240; 0.3 400; 0.2 600; 0.1 1200; 0.08 1500; 0.06 2000; 0.05 2400; 0.03 4000; 0.02 6000; 0.01 12000]
# configs1[1,2]

for i=1:length(configs1[:,1])
    m = configs1[i,1]
    n_tau = configs1[i,2]
    main(n_tau,"expfulla$(i).csv","measuredObsa$(i).csv","expectationvala$(i).csv",m,m)
end


begin
    plot(GetColumn(2,"results/expectationvala5.csv"))    # ⟨x̂⟩
    hline!([0])
end
begin
    plot(GetColumn(3,"results/expectationvala2.csv"))# ⟨x̂²⟩
    expected_x2 = Exp_x2(12000, 0.01, 0.01)
    hline!([expected_x2])
end
expected_x4 = 3*Exp_x2(120, m, ω)^2

begin
    PlotProbDD("results/measuredObsb.csv")
    PlotProbDDe(m,ω,1)
end

LastRowFromFile("results/measuredObsA14.csv")
LastRowFromFile("results/measuredObsA13.csv")

LastRowFromFile("results/expfullA14.csv")
LastRowFromFile("results/expfullA13.csv")
plot!(LastRowFromFile("results/expfullA13.csv")[2:120+1])

@time for i=1:10000
    writec123tofile("Benchmarks/","old.csv",[i for i=1:1000],i)
end






#############################
#       AutoCorrelation     #
#############################
data = GetColumn(2,"results/measuredObsA12.csv")
# LastRowFromFile("results/measuredObs.csv")
autocorrdata = StatsBase.autocor(data,[i for i=0:length(data)-1];demean=true)
plot(autocorrdata,title="AutoCorr by StatsBase package")

data2 = append!(copy(data),[0 for i=1:length(data)])
autocorrdata2 = real.(AutoCorrR(data2))[1:length(data)]
plot(autocorrdata2,title="AutoCorr by padded data by fourier transform")


###### Plot ⟨xᵢ⟩, ⟨xᵢ²⟩
plotexpx(1,4,meanf)         # corr1
plotexpx1(2,true)           # exp_x_1
plotexpx1(2,false)          # exp_x_1
plotexpx2(2,n_tau,false)    # exp_x_2
plotexpx(n_tau+1, 20)       # corr1


######### Plot final ⟨xᵢ⟩, and ⟨xᵢ²⟩ #########
begin
    plot([0 for i=1:n_tau])
    lastRow = LastRowFromFile("results/expfull40.csv")
    exp3 = [lastRow[i] for i=2:n_tau+1]
    PlotExp(exp3,3)
    plot!([expected_x2/2 for i=1:n_tau])
    exp3 = [lastRow[i+n_tau] for i=2:n_tau+1]
    PlotExp(exp3,3)
end


######### Plot final Two-Point Correlation in ⟨xᵢ⟩ #########
lastRow = LastRowFromFile("results/expfull.csv")
for i = 1:length(lastRow)
    if lastRow[i] <= 0
        lastRow[i] = NaN
    end
end
PlotTwoPointCorrelation(lastRow[2*n_tau+1:3*n_tau+1])
# exp3 = [lastRow[1:16],lastRow[17:32]]
# exp3 = [lastRow[i] for i=2:n_tau+1]
# exp3 = [lastRow[i+16] for i=2:17]
# PlotExp(exp3,3)


######### Plot ⟨x⟩ #########
n_tau
begin
    plot_x(n_tau,"results/expfullA12.csv",[2,4,6,8])#[5,12,22,32,42,52])#0)#[i for i=1:n_tau])
    plot!([0 for i=1:length(readlines("results/expfullA12.csv"))], title="Expectationvalue of x", label="⟨x⟩ₜₕₑ",legend=:bottom)
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







@time err = Err(n_tau,meanf,obsf)



function PlotExpE_x(num, meanf, obsf, err)
    if length(err) < num
        println("length(err) < num")
        err = Err(num, meanf, obsf)
    end
    plot([0 for i=1:num], label="⟨x⟩ₑₓₚₑ=0")
    lastRow = LastRowFromFile(meanf)
    exp3 = [lastRow[i] for i=2:num+1]
    title = "Expectationvalue x"
    label = "⟨xᵢ⟩"
    plot!(exp3, yerr=err[1:num], title=title, label=label)
    #plote_x(exp3,err)
    # exp3 = [lastRow[i+n_tau] for i=2:n_tau+1]
    # PlotExp(exp3,3)
end

PlotExpE_x(n_tau, meanf, obsf, err)






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

plot(real.(AutoCorrR(path)))

path[1:4999-20]
mean(path[1:5000])
tau_MC = 5
N_meas = 4999
a = [path[i] for i=1:N_meas-tau_MC]
mean(a)



column = GetColumn(1,"results/measuredObs1.csv")
plot(real.(AutoCorrR(column)))


plot(real.(AutoCorrR(GetColumn(2,"results/measuredObsG2.csv"))))
AutoCorrR(rand(gaussianD,200))

plot(real.(AutoCorrR(rand(gaussianD,2000))))


"""(DEVELOPMENT)  
AutoCorrelation function
"""
function AutoC(n_tau, reducedlength, meanf, obsf)
    if typeof(meanf) == String
        if reducedlength == 0
            mean3 = LastRowFromFile(meanf)[2:n_tau+1]      # Final mean
        else
            mean3 = parse.(Float64,split(readlines(meanf)[reducedlength],","))[2:n_tau+1]
        end
    else
        mean3 = meanf
    end
    #
    # Append all obs to rowData
    rowData = []
    if typeof(obsf) == String
        @time if reducedlength == 0
            for r = readlines(obsf)
                push!(rowData,parse.(Float64,split(r,",")[2:n_tau+1]))
            end
        else
            for r = readlines(obsf)[1:reducedlength]
                push!(rowData,parse.(Float64,split(r,",")[2:n_tau+1]))
            end
        end
    else
        rowData = obsf
    end
    println("Created rowData")
    # println(rowData)
    # Calculate autocorrelation
    N_meas = length(rowData)
    autocorr = zeros(Float64, N_meas-1, n_tau)
    if n_tau == 1
        for tau_MC=0:N_meas-2
            for ii = 1:N_meas-tau_MC
                # println(tau_MC+1," ",1," ",ii," ",ii+tau_MC," ",)
                # println([rowData[i] for i=N_meas-tau_MC:length(rowData)],length(rowData))
                meand1 = [rowData[i] for i=1:N_meas-tau_MC]
                meand2 = [rowData[i] for i=tau_MC+1:N_meas]
                if length(meand1) == 1
                    mean1 = meand1[1]
                else
                    mean1 = mean(meand1)
                    # println(length(rowData))
                    # println(length(meand1),typeof(meand1))
                    # println(mean([1,2,3]))
                    # println(success!)
                end
                if length(meand2) == 1
                    mean2 = meand2[1]
                else
                    mean2 = mean(meand2)
                end
                autocorr[tau_MC+1,1] += (rowData[ii]-mean1)*(rowData[ii+tau_MC]-mean2)/(N_meas-tau_MC-1)
                # println("Term1 ",rowData[ii]-mean([rowData[i] for i=1:N_meas-tau_MC]))
                # println("Term2 ",rowData[ii+tau_MC]-mean([rowData[i] for i=(N_meas-tau_MC):length(rowData)]))
                # println("Term3 ",(N_meas-tau_MC-1))
            end
            println("Finished tau_MC ",tau_MC," / ", N_meas-2)
        end
        for ii = 2:N_meas-1
            autocorr[ii,1] /= autocorr[1,1]
        end
        autocorr[1,1] /= autocorr[1,1]
        println("Calculated the autocorrelation for the (",1,"/",n_tau,")")
        return autocorr
    end
    # println("Autocorrelation: ")
    # show(IOContext(stdout, :limit => true),"text/plain",autocorr)
    println("\nRowData: ")
    show(IOContext(stdout, :limit => true),"text/plain",rowData)
    # println(rowData[100][100], " ",rowData[2][2])
    println("\nAutoCorrelation[2,2]: ",autocorr[2,2])
    for i=1:n_tau
        @time for tau_MC=0:N_meas-2
            mean1 = mean([rowData[iii][i] for iii=1:N_meas-tau_MC])
            mean2 = mean([rowData[iii][i] for iii=tau_MC+1:N_meas])
            for ii = 1:N_meas-tau_MC
                # println(tau_MC+1," ",i," ",ii," ",ii+tau_MC)
                autocorr[tau_MC+1,i] += (rowData[ii][i]-mean1)*(rowData[ii+tau_MC][i]-mean2)/(N_meas-tau_MC-1)
            end
            # if false#tau_MC > N_meas-5
            #     show(IOContext(stdout, :limit => true),"text/plain",rowData)
            #     println("\nmean1= ", mean1," of: ",[rowData[iii][i] for iii=1:N_meas-tau_MC])
            #     println("mean2= ", mean2," of: ",[rowData[iii][i] for iii=tau_MC+1:N_meas])
            #     for ii = 1:N_meas-tau_MC
            #         println("e1: ",rowData[ii][i]," e2: ",rowData[ii+tau_MC][i])
            #     end
            #     println()
            # end
        end
        for ii = 2:N_meas-1
            autocorr[ii,i] /= autocorr[1,i]
        end
        autocorr[1,i] /= autocorr[1,i]
        println("Calculated autocorrelation for (",i,"/",n_tau,")")
        println("Autocorrelation: ")
        show(IOContext(stdout, :limit => true),"text/plain",autocorr[i])
    end
    return autocorr
end

@time autoc = AutoC(2, 0,"results/expfullG.csv","results/measuredObsG.csv")
plot(autoc[:,1])

begin
    raw = []
    for r = readlines(obsf)
        push!(raw,parse.(Float64,split(r,",")[2]))
    end
    mean4 = LastRowFromFile(meanf)[2]
    # raw=[rowdata[i][1] for i=1:length(rowdata)]
    # mean4
    a_0=0
    a_5=0
    a_τ=0
    for i=1:length(raw)
        a_0+=(raw[i]-mean4)*(raw[i]-mean4)/(length(raw)-1)
    end
    m1=0
    m2=0
    τ_MC=length(raw)-5
    for i = 1:length(raw)-τ_MC
        m1 += raw[i]/(length(raw)-τ_MC)
        m2 += raw[length(raw)-i+1]/(length(raw)-τ_MC)
    end
    for i=1:length(raw)-τ_MC
        a_τ += (raw[i]-m1)*(raw[i+τ_MC]-m2)/(length(raw)-τ_MC-1)
    end
end
# Autocorrelation does not work for gaussian even

begin
    autocorrGaussian = AutoC(1,0,mean(path),path)
    autocorrGaussian
    plot(autocorrGaussian[:,1])
end
m1
m2
a_0
# a
a_τ
a_τ/a_0
# n_tau

# Ph = Pka + log10([a]/[b])
# ph-pka=log10([a]/[b])
# 10^(ph-pka)=[a]/[b]
# pka=2.2
# ph=4.2
# x=10^(ph-pka)
# 10^(ph-pka)/(10^(ph-pka)+1)*100


auto = sum(autoc, dims=2)
for i=1:1019
    auto[i] /= 1200
end
plot(auto)
auto
row(autoc)

@time begin
    plot([0 for i=1:length(autoc[:,1])],legend=:bottom)
    for i = 1:400:n_tau-1
        plot!(autoc[:,i], label="A(x_$(i)))")
    end
    plot!(autoc[:,n_tau], label="A(x_$(n_tau))")
end


"""
# n = 1; randinc = 0;
# Path, accept, sum1, sum2, sum3, (exp_x, exp_x2, exp_x0x1, randinc) = zeros(n_tau), zeros(n_tau), zeros(n_tau), zeros(n_tau), zeros(n_tau), MetropolisUpdate(n, randinc);
# sum1
# exp_x, exp_x2, exp_x0x1 = MetropolisUpdate(20);


# ab = DataFrame(Name = ["AKANKSHA", "TANYA", "PREETIKA", "VRINDA", "JAHNVI"],
#                Age = [42, 44, 22, 81, 93],
#                Salary = [540000, 650000, 900000, 770000, 850000],
#          RESIDENCE=["DELHI", "DELHI", "UP", "HARYANA", "UP"]
#                )

# CSV.write(string(rfold,"sum3exp3.csv"),ab)
# file = CSV.read(string(rfold,"sum3exp3n1.csv"),DataFrame)


# CSV.write()

# randlist=zeros(100)
# for i = 1:100
#     randlist[i]=rand(rng)
# end

# randlist
# M = zeros(10,10)
# M[1,1] = randn()
# M

# length(a)
# size(M)

# z = zeros(Integer, 10)
# o = ones(2,2)
# r = randn(Float64,2*n_tau*100*10000)


#global randinc = 0
# for i = 1:n_tau
#     Path[i] = 1
# end
# accept

# randinc
function printarrayold(array,sums)
    exval = zeros(length(array))
    for i = 1:length(array)
        if array[i] < 0
            @printf("%.2f ", (array[i]/sums))
        else
            @printf(" %.2f ", (array[i]/sums))
        end
        exval[i] = array[i]/sums[i]
    end
    print("\n")
    return exval
end

#                               #
# Defining a metropolis swipe   #
#                               #
function MetropolisUpdate(n_burn, n, randinc)
    for i = 1:n
        MetroSwipe()
    end
    
    for i = 1:n
        MetroSwipe()
    end
    exp_x, exp_x2, exp_x0x1 = E_Vals()
    acceptrate = 0
    for i=1:length(accept)
        acceptrate += accept[i]
    end
    # acceptrate = acceptrate/n
    @printf("accept:    %.3f,  %d\n", acceptrate/(n*length(accept)),acceptrate)
    @printf("randinc:   %.3f,  %d\n", randinc/(length(r)),randinc)
    return exp_x, exp_x2, exp_x0x1, randinc
end


# MetroSwipe()
# exp_x, exp_x2, exp_x0x1 = E_Vals();
# sum1


# Correlation(3,"results/expfull.csv")
# nrow(CSV.read("results/expfull.csv", DataFrame))
# println(Gnuplot.gpversion())
# test_terminal()
# for i = 1:length(Path)
#     print(i,"\n")
# end

# begin
#     auto = zeros(Float64,size(autoc,1))
#     @time for i=1:size(autoc,1)
#         for ii=1:size(autoc,2)
#             auto[i] += autoc[i,ii]
#         end
#         println("Calculated ",i)
#     end
#     plot(auto)
# end
    
# autocorrel[:,1]
# ray = readlines("results/measuredObs.csv")[1:48]
# rayow=split(ray, ",")


# begin
#     rowData = []
#     for r = readlines("results/measuredObs.csv")[1:50]
#         push!(rowData,parse.(Float64,split(r,",")[2:n_tau]))
#     end
# end
# rowData

# function r1(l)
#     rowData=Vector{Float64}(undef,0)
#     vec = []
#     for r = readlines("results/measuredObs.csv")[1:l]
#         # vec = Vector{Float64}(split(r,",")[2])
#         append!(rowData,parse.(Float64,split(r,",")[2]))
#         # println(split(r,",")[2])
#     end
#     # println(rowData)
#     plot(rowData)
# end

# r1(1000)

# autocorr = zeros(Float64,20)


# function Count()
#     number = 0
#     while true
#         print(number,"\n")
#         number += 1
#         sleep(1)
#     end
# end
# # Count()

# function testbc()
#     # Testing the boundary conditions of our chain
#     print("\nN_(n+1): ")
#     for i = 1:n_tau
#         coord = i
#         n_p1 = coord % n_tau + 1
#         print(n_p1," ")
#     end
#     print("\nN_(n-1): ")
#     for i = 1:n_tau
#         coord = i
#         n_m1 = (coord - 2 + n_tau) % n_tau + 1
#         print(n_m1," ")
#     end
# end
# # testbc()

# testMatrix = [1 2 3; 4 5 6; 7 8 9]
# for i = 1:len(testMatrix)
#     print(testMatrix[i][1])
# end
# sum(testMatrix, dims=1)


function plotgp()
    #picturename = string(path,figname,".png")
    #@gp "set datafile separator ','"
    x = 1:24
    @gp "set terminal png size 400,300"
    @gp x x.^2 "with lines"
    save("ex1.png"; term="pngcairo size 800,600", output="ex1.png")
    @gp "load(\"ex1.png\")"
    #@gp "plot x myProject/results/expfull.csv using 1 with lines"
    #save(term="pngcairo size 480,360",output=picturename)
end
plotgp()
# plotgp("myProject/results/","expfullpic")
# @gp (1:20).^2


###### Old AutoCorrelation ######
autocorrel = AutoCorrelation(n_tau, 500,"results/expfull.csv","results/measuredObs.csv")
plot(autocorrel[:,1])
begin
    plot([0 for i=1:n_tau-1],legend=:bottom)
    for i = 1:n_tau
        plot!(autocorrel[:,i])
    end
    plot!(autocorrel[:,n_tau])
end
#################################




"""