
using FFTW, Plots, Printf, Random, Statistics

export FFTW, Plots, Printf, Random, Statistics
# export n_tau, idrate, h, m, ω, accept, sum1, sum2, sum3
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




"""Dimensionless mass
"""
# global m = 1
"""Dimensionless natural frequency
"""
# global ω = 1


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

"""Append the current Path, Path^2 and Path1*x on one line to a file"""
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

