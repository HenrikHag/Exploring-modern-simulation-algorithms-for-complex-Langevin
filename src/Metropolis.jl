
using FFTW, Plots, Printf, Random, Statistics, DelimitedFiles

export FFTW, Plots, Printf, Random, Statistics
# export n_tau, idrate, h, m, ω, accept, sum1, sum2, sum3
export AHO_M_param, getAHO_M_param, Sim_M_param, getSim_M_param
export MetropolisSim
# export MetroSwipe, PreSim
export HO_fullAction, HO_Action, AHO_Action, difActionHO, difActionAHO
export MeasureObs, E_Vals#, MeasureObs
export printarray, printmatrix
export writee123tofile, writec123tofile, writec3s3tofile, writes3e3tofile, writeeMean



#                       #
# Defining variables    #
#                       #
"""Number of timesteps/coords
"""
# global n_tau = 12#1200      # Number of coords in Path
"""Optimal acceptrate for simulation, used to adjust h (random change scaler)
"""
# const idrate = 0.8



"""Scale for random change per coord per swipe
"""
# global h = 1           # Scale for random change in coords
# global Path = zeros(n_tau)
# global accept = 0
# global sum1 = zeros(n_tau)   # Path_n
# global sum2 = zeros(n_tau)  # Path_n^2
# global sum3 = zeros(n_tau)  # Path_1 * Path_n





"""
Physical parameters for Metropolis simulation
"""
struct AHO_M_param
    n_tau::Integer
    a::Real
    m::Real
    ω::Real
    λ::Real
end


"""
returns the AHO parameters for Metropolis simulation
"""
function getAHO_M_param(n_tau::Integer,β::Real,m::Real,ω::Real,λ::Real)
    a = β/n_tau
    return AHO_M_param(n_tau,a,m,ω,λ)
end

"""
Simulation parameters for Metropolis simulation
"""
struct Sim_M_param
    N::Integer
    n_burn::Integer
    n_skip::Integer
    idrate::Real
    rng
end

"""
returns the simulation parameters for Metropolis simulation
"""
function getSim_M_param(N::Integer,n_burn::Integer,n_skip::Integer,idrate::Real,rng)
    if n_skip < 1
        println("n_skip must be 1 or higher\nEvery n_skip-1 sample is discarded")
        return
    end
    if idrate < 0 || N < 1 || n_burn < 0
        println("Parameters must be 0 or higher (N 1 or higher)")
        return
    end
    return Sim_M_param(N,n_burn,n_skip,idrate,rng)
end




#                                #
# Defining metropolis functions  #
#                                #
"""
Does a metroswipe by testing a change for each element,  
and adds this new coord weighted by change in action.
"""
function MetroSwipe(n_tau::Int64, m, ω, λ, a, h, idrate, rng, Path)
    accept = 0
    for i = 1:n_tau
        x_new = Path[i] + h*2*(rand(rng)-1/2)
        if rand(rng) < exp(-difActionAHO(n_tau,a,m,ω,λ,Path,i,x_new))
            Path[i] = x_new
            accept += 1/n_tau
        end
    end
    if accept > 0
        h *= accept / idrate
    end
    return Path, accept, h
end

"""
Run small Simulation to estimate Autocorrelation at τ = 1  
Assumes n_burn is large enough to overcome burn-in time
"""
function PreSim(n_tau, m, ω, λ, a, h, idrate, rng, Path, n_burn, n_skip, accept, runt)
    h2 = copy(h)
    Path2 = copy(Path)
    accept2 = copy(accept)
    for i=1:n_burn
        Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
    end
    configs = Matrix{Float64}(undef,runt,n_tau)
    for i=1:runt
        for ii=1:n_skip
            Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
        end
        configs[i,:] = Path
    end
    return 10*mean(AutoCorrR(configs)[:,2]), accept2, h2, Path2
end

"""
Simulation with Metropolis algorithm
"""
function MetropolisSim(phys_param::AHO_M_param,sim_param::Sim_M_param,save_pre,filename)
    # Logic of program should go here
    # We want to print final expectation values
    # rng = MersenneTwister(11111)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.a
    m = phys_param.m
    ω = phys_param.ω
    λ = phys_param.λ
    println("\nn_tau = ",n_tau,", δτ = ", a," m,ω = ", m,", ",ω)
    Path = [0. for i=1:n_tau]

    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    
    # Simulation parameters
    N = sim_param.N
    n_burn = sim_param.n_burn
    n_skip = sim_param.n_skip
    idrate = sim_param.idrate
    rng = sim_param.rng
    accept = 0
    h = 1

    # Prepare writing to file       #
    obsfname = "$(save_pre)$(filename)_obs.csv"
    touch(obsfname)
    
    # meanfname = "$(filename)_expfull.csv"
    # touch("$(save_pre)$(meanfname)")
    # expfname = "$(filename)_expect.csv"
    # touch("$(save_pre)$(expfname)")
    # show(IOContext(stdout, :limit => true),"text/plain",Path); println();
    
    # Burn in, removing first n_burn values
    for i = 1:n_burn
        Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
    end
    
    # Simulation saving configurations
    itt = 0
    Randlist = []
    a = @timed for i = 1:N*n_skip
        Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
        if n_skip == 0 || (i-1)%n_skip == 0
            append!(Randlist,accept)
            itt += 1
            writec123tofile(obsfname, Path, itt)   # curr1, curr2, curr3  // Append each itt
            # println("Measured ",itt)
            # sum1, sum2, sum3 = MeasureObs(n_tau, sum1, sum2, sum3, Path)
            # exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
            # writee123tofile(n_tau,save_pre,meanfname, exp_x,exp_x2,exp_x0x1, itt) # exp1, exp2, exp3  // Append each itt
            # writeeMean(save_pre,expfname,exp_x,exp_x2,exp_x0x1,itt)   # ⟨exp1⟩, ⟨exp2⟩, ⟨exp3⟩         // Append each itt
        end
    end
    println("Mean acceptrate = ", Statistics.mean(Randlist))
    println("Main simulation time: ",a.time)
    return
end
function MetropolisSim(phys_param::AHO_M_param,sim_param::Sim_M_param,save_pre,filename,presim_AC::Bool)
    # Logic of program should go here
    # We want to print final expectation values
    # rng = MersenneTwister(11111)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.a
    m = phys_param.m
    ω = phys_param.ω
    λ = phys_param.λ
    println("\nn_tau = ",n_tau,", δτ = ", a," m,ω = ", m,", ",ω)
    Path = [0. for i=1:n_tau]

    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    sum1, sum2, sum3 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    
    # Simulation parameters
    N = sim_param.N
    n_burn = sim_param.n_burn
    n_skip = sim_param.n_skip
    accept = 0
    idrate = 0.8
    h = 1

    # Make autocorrelation negligible (optional)  #
    while presim_AC
        runt = 200
        ato2, accept, h, Path = PreSim(n_tau, m, ω, λ, a, h, idrate, rng, Path, n_burn, n_skip, accept, runt)
        if ato2 > 2
            println("Autocorrelation after $(runt) runs show Aₒ(1)= $(ato2/10)")
            println("Increasing n_skip from: $(n_skip) to $(round(n_skip*ato2))\n")
            n_skip *= round(ato2)
        else
            println("Autocorrelation after $(runt) runs show Aₒ(1)= $(ato2/10)")
            break
        end
    end

    # Prepare writing to file       #
    obsfname = "$(save_pre)$(filename)_obs.csv"
    touch(obsfname)
    
    meanfname = "$(filename)_expfull.csv"
    touch("$(save_pre)$(meanfname)")
    expfname = "$(filename)_expect.csv"
    touch("$(save_pre)$(expfname)")
    # show(IOContext(stdout, :limit => true),"text/plain",Path); println();
    
    # Burn in, removing first n_burn values
    for i = 1:n_burn
        Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
    end
    
    # Simulation saving configurations
    itt = 0
    Randlist = []
    a = @timed for i = 1:N*n_skip
        Path, accept, h = MetroSwipe(n_tau, m, ω, λ, a, h, idrate, rng, Path)
        if n_skip == 0 || (i-1)%n_skip == 0
            append!(Randlist,accept)
            itt += 1
            writec123tofile(obsfname, Path, itt)             # curr1, curr2, curr3            // Append each itt
            # println("Measured ",itt)
            sum1, sum2, sum3 = MeasureObs(n_tau, sum1, sum2, sum3, Path)
            exp_x, exp_x2, exp_x0x1 = E_Vals(n_tau,sum1,sum2,sum3,itt)
            writee123tofile(n_tau,save_pre,meanfname, exp_x,exp_x2,exp_x0x1, itt) # exp1, exp2, exp3  // Append each itt
            #filenamec3s3 = string("curr3sum3n",Int64(itt),".csv") # curr3, sum3                   // New file
            #writec3s3tofile(rfold,filenamec3s3, sum3)
            #filenames3e3 = string("sum3exp3n",Int64(itt),".csv")  # sum3, exp3                  // New file
            #writes3e3tofile(rfold,filenames3e3, sum3, exp_x0x1)
            writeeMean(save_pre,expfname,exp_x,exp_x2,exp_x0x1,itt)   # ⟨exp1⟩, ⟨exp2⟩, ⟨exp3⟩      // Append each itt
        end
    end
    println("Mean acceptrate = ", Statistics.mean(Randlist))
    println("Main simulation time: ",a.time)
    return
end


#                       #
# Defining the action   #
#                       #
"""
Computes the full Harmonic Oscillator action of array1
"""
function HO_fullAction(array1::AbstractArray,a,m,ω)
    return 0.5*m*sum([(array1[(i%length(array1))+1]+array1[i])^2/a + a*ω^2*array1[i]^2 for i=1:length(array1)])
end


"""Harmonic Oscillator change of Action  
Returns the part of the action dependent on coord at index i of Path  
Path:       The global Path  
Path[i]:    The index we want to propose a change to  
x:          Either the old coord, or the proposed new coord  
```julia
return part_of_action_value
```"""
function HO_Action(n_tau, m, ω, a, Path, coord, x)
    # n_p1 = coord % n_tau + 1
    # n_m1 = (coord-2) % n_tau + 1
    # print(n_p1)
    # print(n_m1)
    return 0.5*m*a*(((Path[((coord)%n_tau)+1]-x)/a)^2
    + ((x-Path[((coord-2+n_tau)%n_tau)+1])/a)^2
    + ω^2*x^2)
end
# Action(Path,2,Path[2])
# print(Path[((2-2)%16+1)])

function AHO_Action(n_tau, m, ω, a, λ, Path, coord, x)
    return 0.5*m*a*(((Path[((coord)%n_tau)+1]-x)/a)^2
    + ((x-Path[((coord-2+n_tau)%n_tau)+1])/a)^2
    + ω^2*x^2)# + 1/(4*a^4)*λ*x^4
end

"""
return change in HO action (s_new-s_old) by trial x_new
"""
function difActionHO(n_tau,a,m,ω,Path,index,x_new)
    return m*((x_new^2-Path[index]^2+(Path[index]-x_new)*(Path[(index %n_tau)+1]+Path[((index-2+n_tau)%n_tau)+1]))/a
              + 0.5*a*ω^2*(x_new^2-Path[index]^2))
end

"""
return change in AHO action (s_new-s_old) by trial x_new
"""
function difActionAHO(n_tau,a,m,ω,λ,Path,index,x_new)
    return (m*((x_new^2-Path[index]^2+(Path[index]-x_new)*(Path[(index %n_tau)+1]+Path[((index-2+n_tau)%n_tau)+1]))/a
               + 0.5*a*ω^2*(x_new^2-Path[index]^2))
            + λ*(x_new^4-Path[index]^4)/(24. *a^4))
end




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






#                                   #
# Defining the expectationvalues    #
#                                   #
"""Updates the expectationvalues based on the current value of sums and n_accept:
```julia
for i = 1:n_tau
    exp_x[i] = sum1[i]/n_accept
end
return exp_x, exp_x2, exp_x0x1
```"""
function E_Vals(n_tau,sum1,sum2,sum3,n_accept)
    # Print the expectation values
    # print("<x_i>     "); exp_x = zeros(n_tau)
    # print("<(x_i)^2> "); exp_x2 = zeros(n_tau)
    # print("<x_1*x_i> "); exp_x0x1 = zeros(n_tau)
    exp_x, exp_x2, exp_x0x1 = zeros(n_tau), zeros(n_tau), zeros(n_tau)
    for i = 1:n_tau
        exp_x[i] = sum1[i]/n_accept
        exp_x2[i] = sum2[i]/n_accept
        exp_x0x1[i] = sum3[i]/n_accept
    end
    return exp_x, exp_x2, exp_x0x1
end






#                                       #
# Beautify output of expectationvalues  #
#                                       #
"""Beautify printing of float-array elements.  
Prints array on one line with spaces between rounded .2f elements
"""
function printarray(array)
    for i = 1:length(array)
        if array[i] < 0
            @printf("%.2f ", array[i])
        else
            @printf(" %.2f ", array[i])
        end
    end
    print("\n")
    return
end

function printmatrix(matrix)
    show(IOContext(stdout, :limit => true),"text/plain",matrix);println();
end



#                       #
# Writing to files      #
#                       #
"""Append the current exp_x, exp_x2, exp_x0x1 on one line to a file"""
function writee123tofile(n_tau,path,filename,exp_x,exp_x2,exp_x0x1, itt)
    open(string(path,filename),"a") do file
        # write(file,"exp_x,exp_x2,exp_x0x1\n")
        write(file,string(Int64(itt),","))
        for i = 1:n_tau
            write(file,string(exp_x[i],","))
        end
        for i = 1:n_tau
            write(file,string(exp_x2[i]^2,","))
        end
        for i = 1:n_tau-1
            write(file,string(exp_x0x1[i],","))
        end
        write(file,string(exp_x0x1[n_tau],"\n"))
    end
end

"""
Append the current Path, Path^2 and Path1*x on one line to a file
"""
function writec123tofile(path, filename, Path, itt)
    open(string(path,filename),"a") do file
        # write(file,"c_x,c_x2,c_x0x1\n")
        write(file,string(Int64(itt),","))
        pathl=length(Path)
        for i = 1:pathl         # Path
            write(file,string(Path[i],","))
        end
        for i = 1:pathl         # Path.^2
            write(file,string(Path[i]^2,","))
        end
        for i = 1:pathl         # Path_1*Path_i
            write(file,string(Path[1]*Path[i],","))
        end
        for i = 0:pathl-2       # Two-Point Correlation
            twopointcorr=0
            for ii=1:pathl
                twopointcorr += Path[ii]*Path[(ii+i-1)%pathl+1]
            end
            write(file,string(twopointcorr/pathl,","))
        end
        twopointcorr=0
        for ii=1:pathl
            twopointcorr += Path[ii]*Path[(ii+pathl-2)%pathl+1]
        end
        write(file,string(twopointcorr/pathl,"\n"))
    end
end
function writec123tofile(save_name::AbstractString, Path, itt)
    open(save_name,"a") do file
        # write(file,"c_x,c_x2,c_x0x1\n")
        write(file,string(Int64(itt),","))
        pathl=length(Path)
        for i = 1:pathl         # Path
            write(file,string(Path[i],","))
        end
        for i = 1:pathl         # Path.^2
            write(file,string(Path[i]^2,","))
        end
        for i = 1:pathl         # Path_1*Path_i
            write(file,string(Path[1]*Path[i],","))
        end
        for i = 0:pathl-2       # Two-Point Correlation
            twopointcorr=0
            for ii=1:pathl
                twopointcorr += Path[ii]*Path[(ii+i-1)%pathl+1]
            end
            write(file,string(twopointcorr/pathl,","))
        end
        twopointcorr=0
        for ii=1:pathl
            twopointcorr += Path[ii]*Path[(ii+pathl-2)%pathl+1]
        end
        write(file,string(twopointcorr/pathl,"\n"))
    end
end
# function writec123tofile(save_name::AbstractString,Path,itt::Integer)
#     open(save_name,"a") do file
#         # write(file,"$(itt),")
#         writedlm(file,transpose([itt,Path...,Path.^2...,
#                         [Path[1]*Path[i] for i=1:length(Path)]...#[sum(Path[ii].*circshift(Path,i)[ii] for i=0:length(Path)-1]./length(Path)]
#                         ]),',')
#         # write(file,",")
#     end
# end

function writec3s3tofile(path,filename,sum3)
    open(string(path,filename),"a") do file
        write(file,"current3,sum3\n")
        for i = 1:n_tau
            write(file,string(Path[1]*Path[i],",",sum3[i],"\n"))
        end
    end
end

function writes3e3tofile(path,filename,sum3,exp_x0x1)
    open(string(path,filename),"a") do file
        write(file,"sum3,exp_x0x1\n")
        for i = 1:n_tau
            write(file,string(sum3[i],",",exp_x0x1[i],"\n"))
        end
    end
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

