using Distributions
using Plots
using .UsersGuide
gaussianD = Normal(0,1)

#       μ/2*ϕ²+λ/4!*ϕ⁴
#       -> (μϕ+λ/6*ϕ³)*Δt
function Action(mu,la,F,dt)
   return -(mu*F+0*la/6*F^3)*dt
end

function Langevin(N,mu,la,gaussianD)
    F = [20. for i = 1:16]
    Flist = []
    dt = 0.01
    for i=1:N
        push!(Flist,F)
        for ii = 1:length(F)
            F[ii] += Action(mu,la,F[ii],dt) + sqrt(2*dt)*rand(gaussianD)
        end
    end
    return Flist
end
function Langevin(N,mu,la,gaussianD,filename)
    path = "results/"
    F = [20. for i = 1:16]
    Flist = []
    dt = 0.01
    n_burn = 200
    for i=1:n_burn
        for ii = 1:length(F)
            F[ii] += Action(mu,la,F[ii],dt) + sqrt(2*dt)*rand(gaussianD)
        end
    end
    for i=1:N
        push!(Flist,F)
        open(string(path,filename),"a") do file
            # write(file,"exp_x,exp_x2,exp_x0x1\n")
            write(file,string(Int64(i),","))
            for ii = 1:length(F)-1
                write(file,string(F[ii],","))
            end
            write(file,string(F[end],"\n"))
        end
        for ii = 1:length(F)
            F[ii] += Action(mu,la,F[ii],dt) + sqrt(2*dt)*rand(gaussianD)
        end
    end
    return Flist
end
println(Langevin(20,1,0.4,gaussianD,"CL_1.csv"))
Langv1=Langevin(100000,1,0.4,gaussianD,"CL_1.csv")
# Probability density diagram #
arr1 = [Langv1[i][j] for i = 1:length(Langv1) for j = 1:length(Langv1[1])]
histogram(arr1,normed=true,xlabel="x",ylabel="|ψ_0|²")
histogram(Langv1,normed=true,xlabel="x",ylabel="|ψ_0|²")
PlotProbDD("results/CL_1.csv",0.1)
PlotProbDDe(1/2,1,1,3)
# sampling
plot(reshape(GetColumn(2,"results/CL_1.csv"),:))#:Int((length(LastRowFromFile(file))-1)/4)+1
# Autocorrelation
PlotAC("results/CL_1.csv",1000)
PlotAC("results/CL_1.csv",false)
