module MetropolisUpdate

using Random, Printf, Plots, Gnuplot, Statistics, FFTW

export Random, Printf, Plots, Gnuplot
# export n_tau, idrate, h, m, ω, accept, sum1, sum2, sum3
export Exp_x2, AutoCorrR, HO_Action, AHO_Action, difActionAHO, MeasureObs, E_Vals#, MeasureObs
export printarray, writee123tofile, writec123tofile, writec3s3tofile, writes3e3tofile, writeeMean
export plotexpx, plotexpx1, plotexpx2

#, AutoCorrelation



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

#                                                              #
# Defining the analytic formulas for values to be calculated   #
#                                                              #
function Exp_x2(n_tau, a, m, ω)
    R = 1+((a*ω)^2)/2-a*ω*sqrt(1+((a*ω)^2)/4)
    return 1/(2*a^2*m*ω*sqrt(1+(a*ω/2)^2))*(1+R^n_tau)/(1-R^n_tau)
end


# nrow()
"""Returns AutoCorrelation of arrayC.  
Returns a matrix of correlations → if matrix of configs ↓ is passed.  
"""
function AutoCorrR(arrayC)
    mean1 = mean(arrayC)
    arrayCm = arrayC .- mean1
    autoCorr = fft(arrayCm)
    arrayCm = (abs.(autoCorr)).^2
    autoCorr = ifft(arrayCm)
    e1 = autoCorr[1]
    return (autoCorr)./e1
end
function AutoCorrR(arrayC,norm::Bool)
    mean1 = mean(arrayC)
    arrayCm = arrayC .- mean1
    autoCorr = fft(arrayCm)
    arrayCm = (abs.(autoCorr)).^2
    autoCorr = ifft(arrayCm)
    e1 = autoCorr[1]
    if !norm
        return autoCorr
    end
    return (autoCorr)./e1
end
function AutoCorrR(matrixC::AbstractMatrix)
    CorrD=Matrix{Float64}(undef,length(matrixC[1,:]),length(matrixC[:,1]))
    for i=1:length(matrixC[1,:])
        CorrD[i,:]=real.(AutoCorrR(append!(matrixC[:,i],[0 for i=0:length(matrixC[:,1])])))[1:length(matrixC[:,1])]
    end
    return CorrD
end
function AutoCorrR(matrixC::AbstractMatrix,norm::Bool)
    CorrD=Matrix{Float64}(undef,length(matrixC[1,:]),length(matrixC[:,1]))
    for i=1:length(matrixC[1,:])
        CorrD[i,:]=real.(AutoCorrR(append!(matrixC[:,i],[0 for i=0:length(matrixC[:,1])]),norm::Bool))[1:length(matrixC[:,1])]
    end
    return CorrD
end
function AutoCorrR(matrixC::AbstractMatrix,norm::Bool,padded::Bool)
    CorrD=Matrix{Float64}(undef,length(matrixC[1,:]),length(matrixC[:,1]))
    if !padded
        for i=1:length(matrixC[1,:])
            CorrD[i,:]=real.(AutoCorrR(matrixC[:,i],norm::Bool))[1:length(matrixC[:,1])]
        end
        return CorrD
    end
    for i=1:length(matrixC[1,:])
        CorrD[i,:]=real.(AutoCorrR(append!(matrixC[:,i],[0 for i=0:length(matrixC[:,1])]),norm::Bool))[1:length(matrixC[:,1])]
    end
    return CorrD
end

#                       #
# Defining the action   #
#                       #
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
    + ω^2*(x)^2)
end
# Action(Path,2,Path[2])
# print(Path[((2-2)%16+1)])

function AHO_Action(n_tau, m, ω, a, λ, Path, coord, x)
    return 0.5*m*a*(((Path[((coord)%n_tau)+1]-x)/a)^2
    + ((x-Path[((coord-2+n_tau)%n_tau)+1])/a)^2
    + ω^2*(x)^2) + 1/(4*a^4)*λ*x^4
end

function difActionAHO(n_tau,a,m,ω,λ,Path,index,proposedx)
    return m*((proposedx^2-Path[index]^2+(Path[index]-proposedx)*(Path[((index)%n_tau)+1]+Path[((index-2+n_tau)%n_tau)+1]))/a
    + 0.5*a*ω^2*(proposedx^2-Path[index]^2))
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



#                       #
# Creating plots        #
#                       #
"""Used to plot the n first values of exp_x in expfull.csv"""
function plotexpx(b, n, meanf)
    # x = nrow(CSV.read("results/expfull.csv", DataFrame))
    @gp "set datafile separator ','" :-
    @gp :- "set terminal png size 800,600" tit="exp_x" :-
    for i = b:n
        @gp :- "plot \"$(meanf)\" using 1:$(i+1) w l tit 'exp(x_$i)'" :-
    end
    @gp :- "set xzeroaxis; set yrange [-2:2]"
    #@gp :- """plot "results/expfull.csv" using 1:3 w l tit 'x_2'"""
    save(term="pngcairo size 800,600", output="corr1.png")
    #@gp """load("corr1.png")"""
end

"""Used to plot the nth value of exp_x in expfull.csv  
individualm true/false: Plot the individual measurements / don't"""
function plotexpx1(n,individualm)
    # x = nrow(CSV.read("results/expfull.csv", DataFrame))
    @gp "set datafile separator ','" :-
    @gp :- "set terminal png size 400,300" tit="exp_x" :-
    for i = [n]
        @gp :- "plot \"results/expfull.csv\" using 1:$(i+1) w l tit 'exp(x_$i)'" :-
        if individualm
            @gp :- "plot \"results/measuredObs.csv\" using 1:$(i+1) tit 'obs(x_$i)'" :-
        end
    end
    @gp :- "set xzeroaxis"
    #@gp :- """plot "results/expfull.csv" using 1:3 w l tit 'x_2'"""
    save(term="pngcairo size 800,600", output="exp_x_1.png")
    #@gp """load("corr1.png")"""
end


"""Used to plot the nth value of exp_x2 in expfull.csv  
individualm true/false: Plot the individual measurements / don't"""
function plotexpx2(n, n_tau, individualm)
    # x = nrow(CSV.read("results/expfull.csv", DataFrame))
    @gp "set datafile separator ','" :-
    @gp :- "set terminal png size 400,300" tit="exp_x2" :-
    for i = [n+n_tau]
        @gp :- "plot \"results/expfull.csv\" using 1:$(i+1) w l tit 'exp(x_$i)'" :-
        if individualm
            @gp :- "plot \"results/measuredObs.csv\" using 1:$(i+1) tit 'obs(x_$i)'" :-
        end
    end
    @gp :- "set xzeroaxis"
    #@gp :- """plot "results/expfull.csv" using 1:3 w l tit 'x_2'"""
    save(term="pngcairo size 800,600", output="exp_x_2.png")
    #@gp """load("corr1.png")"""
end










#                                   #
#   Old AutoCorrelation (broken)    #
#                                   #
# """Returns a matrix of autocorrelation for observables with respect to MC time
# """
# function AutoCorrelation(n_tau, reducedlength, meanf, obsf)
#     if reducedlength == 0
#         mean = LastRowFromFile(meanf)[2:n_tau+1]      # Final mean
#     else
#         mean = parse.(Float64,split(readlines(meanf)[reducedlength],","))[2:n_tau+1]
#     end
#     #
#     # Append all obs to rowData
#     rowData = []
#     if reducedlength == 0
#         for r = readlines(obsf)
#             push!(rowData,parse.(Float64,split(r,",")[2:n_tau+1]))
#         end
#     else
#         for r = readlines(obsf)[1:reducedlength]
#             push!(rowData,parse.(Float64,split(r,",")[2:n_tau+1]))
#         end
#     end
#     println("Created rowData")
#     #
#     # Calculate first term for each tau, tau_MC
#     N_meas = length(rowData)
#     autocorr = zeros(Float64, N_meas, n_tau)    # autocorr saves for each tau, tau_MC
#     for i = 1:n_tau
#         # Sum over N_meas - tau_MC
#         # autoc = zeros(Float64, n_tau)               # autoc saves for each tau
#         for tau_MC=0:N_meas-1
#             for ii = 1:N_meas-tau_MC
#                 # println(tau_MC+1," ",i," ",ii," ",ii+tau_MC)
#                 autocorr[tau_MC+1,i] += rowData[ii][i]*rowData[ii+tau_MC][i]/(N_meas-tau_MC)
#             end
#             autocorr[tau_MC+1,i] -= mean[i]^2
#         end
#         println("Calculated first terms (",i,"/",n_tau,")")
#     end
#     #
#     # Normalize by standard deviation
#     for i = 1:n_tau
#         # autoc = zeros(Float64, n_tau)               # autoc saves for each tau
#         autoc=0
#         for ii = 1:N_meas
#             # autoc[i] += rowData[ii][i]^2
#             autoc += rowData[ii][i]^2
#         end
#         for ii = 1:N_meas
#             # println(i," ",ii," ",autoc)
#             # autocorr[ii,i] /= autoc[i]/N_meas - mean[i]^2
#             autocorr[ii,i] /= autoc/N_meas - mean[i]^2
#         end
#         println("Calculated last terms (",i,"/",n_tau,")")
#     end
    

#     # row = split(r,",")[1:17]        # How many observables to calc
#     # rowData = Vector{Float64}(undef,length(row))
#     # for i = 2:length(row)
#     #     rowData[i-1] = parse.(Float64,row[i])-mean[i]
#     # end
#     # println(rowData)
#     # for c = 1:(length(row)-1)
#     #     append!(autoc,rowData[c]*rowData[(c+tau-1)%length(row)+1])
#     # end
#     # tau += 1
#     # append!(autocorr,[autoc])
#     # end

#     # println(autocorr)
    
#     # Normalize by first element
#     # for i=1:length(autocorr)
#     #     for ii=1:length(autocorr[1])
#     #         autocorr[i][ii] = autocorr[i][ii]/autocorr[i][1]
#     #     end
#     # end
#     return autocorr
# end

# #                               #
# # Defining the autocorrelation  #
# #                               #
# """Calculates the autocorrelation of the values in a column from filename for the MC_tstep first elements
# """
# function Correlation(MC_tstep, filename, mean_filename)
#     correlation = 0
#     data = ReadFromFile(filename)
#     mean = ReadFromFile(mean_filename)[length(data)]    # Assume length of data corresponds to elements in mean_file
#     for i = 2:length(data)
#         correlation += (data[i]-mean)*data[i+MC_tstep] - data[]
#     end
# end

# """Reads column from file"""
# function ReadFromFile(filename)
#     open(filename, "r") do file         # Get column as array
#         #for i = 1:length(file)
#         #println(read(file,String))
#         println(length(file))
#     end
# end

# module
end