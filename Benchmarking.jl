# File for benchmarking functions in Metropolis Algorithm

begin
    using .MetropolisUpdate
    using .UsersGuide
    using Plots
    using BenchmarkTools
    using Base.Threads
end

# Define functions to test




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


#                               #
# Defining a main function      #
#                               #
function SimMetro(n_tau,meanfname,obsfname)
    Path = zeros(n_tau)
    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    n_burn = 200
    n_skip = 20
    n_total = n_burn + 1 + (n_skip+1)*10000#20000#200100
    accept = 0
    idrate = 0.8
    h = 1
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
    # println("Burn-in $(n_burn) complete!")
    b = @timed for i = 1:n_total-n_burn
        Path, accept, h = MetroSwipe(n_tau, m, ω, h, idrate, rng, Path)
        if n_skip == 0 || (i-1-n_burn)%n_skip == 0
            sum1, sum2, sum3 = MeasureObs(n_tau, sum1, sum2, sum3, Path)
            append!(Randlist,h)
            itt += 1
            exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
            # writee123tofile(n_tau,rfold,meanfname, exp_x, exp_x2, exp_x0x1, itt)
            # writec123tofile(rfold,obsfname, Path, itt)
        end
    end
    # println("Time: ",a.time, " ", b.time)
    return (a.time, b.time)
end


function main(n_tau,meanfname,obsfname)
    Path = zeros(n_tau)
    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    n_burn = 200
    n_skip = 20
    n_total = n_burn + 1 + (n_skip+1)*10000#20000#200100
    accept = 0
    idrate = 0.8
    h = 1
    rfold = "results/"
    touch(string(rfold,meanfname))
    touch(string(rfold,obsfname))
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
                exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
                # exp1, exp2, exp3      // Append each itt
                # writee123tofile(n_tau,rfold,meanfname, exp_x, exp_x2, exp_x0x1, itt)
                # curr1, curr2, curr3   // Append each itt
                # writec123tofile(rfold,obsfname, Path, itt)
            end
        end
    end
    # println("Time: ",a.time)
    return a.time
end







#                                                   #
#   Benchmark MT or single threaded simulations     #
#                                                   #
mainT = zeros(Float64,0)
SimT = zeros(Float64,0)

for n_tau = 2:20
    println(n_tau)
    # main()
    append!(mainT,main(n_tau,"expfullA14.csv","measuredObsA14.csv"))
    # $ rm -v *.csv
    append!(SimT,SimMetro(n_tau,"expfullA13.csv","measuredObsA13.csv"))
end

plot([(i,mainT[i-1]) for i=2:length(mainT)], ylabel="time (s)", xlabel="n_tau")
mainT
SimT
for i=1:Int(length(SimT)/2)
    SimT[i] = SimT[2*i] + SimT[2*i-1]
end
ab=[(SimT[Int(2*i)] + SimT[Int(2*i-1)]) for i=1:(length(SimT)/2)]
plot!([(i,ab[i-1]) for i=2:length(ab)])


@benchmark main(n_tau,"expfullA14.csv","measuredObsA14.csv")
@benchmark SimMetro(n_tau,"expfullA13.csv","measuredObsA13.csv")




#                                   #
#           Writing to file         #
#                                   #
function WritePathToFile(path, filename, Path, itt)
    row = [itt]
    # Path = repr.(Path)
    println(skipchars(isspace,strip(repr(Path),['[',']'])))
    for i=1:length(Path)
        append!(row,[","*repr(Path[i])])
    end
    for i=1:length(Path)
        append!(row,[","*repr(Path[i]^2)])
    end
    for i=1:length(Path)
        append!(row,[","*repr(Path[i]*Path[1])])
    end
    append!(row,["\n"])
    # row = row + Path[length(Path)]
    show(row)
    open(string(path, filename),"a") do file
        write(file,row)
    end
end

@time for i=1:10000
    WritePathToFile("Benchmarks/","new.csv",[i for i=1:16],i)
end


x = AbstractArray{}

#                                   #
#           Jackknife analysis      #
#                                   #
twopointD = GetTwoPointData("results/measuredObsHO_1_β8_16.csv")
@benchmark Jackknife1(twopointD[:,1])
@benchmark Err1(twopointD[:,1])
# @benchmark Jackknife2(twopointD[:,1])




function Jackknifeold(array1::AbstractArray)
    jf = [mean(array1[2:length(array1)])]
    for i=2:length(array1)-1
        append!(jf,mean(append!(array1[1:i-1],array1[(i+1):length(array1)])))
    end
    append!(jf,mean(array1[1:length(array1)-1]))
    jf = jf.-mean(array1)
    jf = jf.^2
    jfvm = mean(jf)
    # jfvm = 0
    # for i=1:length(jf)
    #     jfvm += (jf[i]-mean(jf))^2
    # end
    jfvm *= (length(array1)-1)#/length(jf)
    # println("Done")
    return [mean(array1),√(jfvm)]
end

Jackknifeold(twopointD[:,1])    # Now named Jackknifeold
Jackknife1(twopointD[:,1])      # New Jackknife1
@benchmark Jackknifeold(twopointD[:,1]) # Old Jackknife1
@benchmark Jackknife1(twopointD[:,1])   # New Jackknife1
