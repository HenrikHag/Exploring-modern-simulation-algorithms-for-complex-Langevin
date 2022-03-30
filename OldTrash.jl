# Dump old functions and calls that have become obsolete


n = 1; randinc = 0;
Path, accept, sum1, sum2, sum3, (exp_x, exp_x2, exp_x0x1, randinc) = zeros(n_tau), zeros(n_tau), zeros(n_tau), zeros(n_tau), zeros(n_tau), MetropolisUpdate(n, randinc);
sum1
exp_x, exp_x2, exp_x0x1 = MetropolisUpdate(20);






#                                               #
#               Dataframes                      #
#                                               #
ab = DataFrame(Name = ["AKANKSHA", "TANYA", "PREETIKA", "VRINDA", "JAHNVI"],
               Age = [42, 44, 22, 81, 93],
               Salary = [540000, 650000, 900000, 770000, 850000],
         RESIDENCE=["DELHI", "DELHI", "UP", "HARYANA", "UP"]
               )

CSV.write(string(rfold,"sum3exp3.csv"),ab)
file = CSV.read(string(rfold,"sum3exp3n1.csv"),DataFrame)


CSV.write()


###########################################################################






#                                               #
#               Error Old impl.                 #
#                                               #
"""Calculate the error  
for observables in obsf  
with mean in meanf
"""
function Err(n_tau, meanf, obsf)
    mean = GetLastMean(meanf, n_tau)
    variance = zeros(Float64,n_tau,1)
    err = zeros(Float64,n_tau,1)
    numb = 0
    @time for r = readlines(obsf)
        for i = 1:n_tau
            variance[i] += (parse.(Float64,split(r,",")[2:n_tau+1])[i] - mean[i])^2
        end
        numb += 1
        if numb%20==0
            println(numb)
        end
    end
    for i=1:n_tau
        variance[i] /= numb-1
    end
    for i = 1:n_tau
        err[i] = √(variance[i]/numb)
    end
    return err
end


###########################################################################






#                                               #
#               Random variables                #
#                                               #
randlist=zeros(100)
for i = 1:100
    randlist[i]=rand(rng)
end

randlist
M = zeros(10,10)
M[1,1] = randn()
M

length(a)
size(M)

z = zeros(Integer, 10)
o = ones(2,2)
r = randn(Float64,2*n_tau*100*10000)


global randinc = 0
for i = 1:n_tau
    Path[i] = 1
end
accept

randinc
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


MetroSwipe()
exp_x, exp_x2, exp_x0x1 = E_Vals();
sum1


Correlation(3,"results/expfull.csv")
nrow(CSV.read("results/expfull.csv", DataFrame))
println(Gnuplot.gpversion())
test_terminal()
for i = 1:length(Path)
    print(i,"\n")
end

begin
    auto = zeros(Float64,size(autoc,1))
    @time for i=1:size(autoc,1)
        for ii=1:size(autoc,2)
            auto[i] += autoc[i,ii]
        end
        println("Calculated ",i)
    end
    plot(auto)
end
    
autocorrel[:,1]
ray = readlines("results/measuredObs.csv")[1:48]
rayow=split(ray, ",")


begin
    rowData = []
    for r = readlines("results/measuredObs.csv")[1:50]
        push!(rowData,parse.(Float64,split(r,",")[2:n_tau]))
    end
end
rowData

function r1(l)
    rowData=Vector{Float64}(undef,0)
    vec = []
    for r = readlines("results/measuredObs.csv")[1:l]
        # vec = Vector{Float64}(split(r,",")[2])
        append!(rowData,parse.(Float64,split(r,",")[2]))
        # println(split(r,",")[2])
    end
    # println(rowData)
    plot(rowData)
end

r1(1000)

autocorr = zeros(Float64,20)


function Count()
    number = 0
    while true
        print(number,"\n")
        number += 1
        sleep(1)
    end
end
# Count()

function testbc()
    # Testing the boundary conditions of our chain
    print("\nN_(n+1): ")
    for i = 1:n_tau
        coord = i
        n_p1 = coord % n_tau + 1
        print(n_p1," ")
    end
    print("\nN_(n-1): ")
    for i = 1:n_tau
        coord = i
        n_m1 = (coord - 2 + n_tau) % n_tau + 1
        print(n_m1," ")
    end
end
# testbc()

testMatrix = [1 2 3; 4 5 6; 7 8 9]
for i = 1:len(testMatrix)
    print(testMatrix[i][1])
end
sum(testMatrix, dims=1)


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


###########################################################################






#                                               #
#               Two-Point Correlation           #
#                                               #
"""Plots two-point correlation from array exp_xjxi, normalized by first element
"""
function PlotTwoPointCorrelation(exp_xjxi)
    # println(exp_xjxi)
    plot([(i,exp_xjxi[i]/exp_xjxi[1]) for i=1:length(exp_xjxi)], title="Two-Point Correlation", label="⟨x₁xᵢ⟩", yaxis=:log10, size=[700,500])
end

PlotTwoPointCorrelation(LastRowFromFile("results/expfull.csv")[33:48])
plotattr(:Plot)

"""G(Δτ)=⟨x(τ)x(τ+Δτ)⟩
"""
function AutoCorrelation(array, n)
    autoCorrelation = 0
    for i = 1:n
        autoCorrelation += array[i]
    end
    return autoCorrelation
end

lastRow = LastRowFromFile("results/expfull.csv")
AutoCorrelation(lastRow, 4)
######### Plot final Two-Point Correlation in ⟨xᵢ⟩ #########
lastRow = LastRowFromFile("results/expfullB1.csv")
for i = 1:length(lastRow)
    if lastRow[i] <= 0
        lastRow[i] = NaN
    end
end

"""Plots two-point correlation from array exp_xjxi, normalized by first element
"""
function PlotTwoPointCorrelation1(exp_xjxi)
    # println(exp_xjxi)
    plot([(i,exp_xjxi[i]) for i=1:length(exp_xjxi)], title="Two-Point Correlation", label="⟨x₁xᵢ⟩")#, yaxis=:log10, size=[700,500])
end
PlotTwoPointCorrelation1(lastRow[2*120+2:3*120+1])

err1 = Err(120, "results/expfullB1.csv", "results/measuredObsB1.csv")

PlotTwoPointCorrelation(lastRow[2*n_tau+1:3*n_tau+1])
exp3 = [lastRow[1:16],lastRow[17:32]]
exp3 = [lastRow[i] for i=2:n_tau+1]
exp3 = [lastRow[i+16] for i=2:17]
PlotExp(exp3,3)


function TwoPointC(data1::AbstractArray)
    array1 = []
    for i in data1
        push!(array1,Err1(i))
    end
    return array1
end
function TwoPointC(data1::AbstractMatrix)
    array1 = Matrix{Float64}(undef,length(data1[1,:]),2)
    for i = 1:length(data1[1,:])
        array1[i,:] = Err1(data1[:,i])
    end
    return array1
end

begin
    tp1 = TwoPointC(twopointD1)
    plot(tp1[:,1],yerr=tp1[:,2])
end
begin
    t2 = TwoPointC(twopointD)
    plot(t2[:,1],yerr=t2[:,2])
end


###########################################################################





#                                               #
#               AutoCorrelation                 #
#                                               #
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


###### Old AutoCorrelation ######
# autocorrel = AutoCorrelation(n_tau, 500,"results/expfull.csv","results/measuredObs.csv")
plot(autocorrel[:,1])
begin
    plot([0 for i=1:n_tau-1],legend=:bottom)
    for i = 1:n_tau
        plot!(autocorrel[:,i])
    end
    plot!(autocorrel[:,n_tau])
end


###########################################################################





#                                               #
#       Convergence of expectationvalue         #
#                                               #
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


###########################################################################


# arr1 = [Langv1[i][j] for i = 1:length(Langv1) for j = 1:length(Langv1[1])]
# histogram(arr1,normed=true,xlabel="x",ylabel="|ψ_0|²")
# histogram(Langv1,normed=true,xlabel="x",ylabel="|ψ_0|²")