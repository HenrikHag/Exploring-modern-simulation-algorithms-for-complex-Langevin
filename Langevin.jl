using Distributions
using Plots
using .UsersGuide
using .MetropolisUpdate
gaussianD = Normal(0,1)

#       μ/2*ϕ² + λ/4!*ϕ⁴ + 1/2*m*((P[i+1]-P[i])²+(P[i]-P[i-1])²)/a²+ 
#       -> (μϕ+λ/6*ϕ³)*Δt
function ActionDer(a,m,mu,la,F,f₋₁,f₊₁)
   return m/a*(2*F-(f₋₁+f₊₁)) + a*m*mu*F# + a*m*la/6*F^3
end
# Add cupling terms 2f(i)-f(i+1)-f(i-i)
# Understand the discretizing integral and meeting mat. from 16.03

function Langevin(N,a,m,mu,la,gaussianD)
    F = [20. for i = 1:16]
    Flist = []
    push!(Flist,F)
    dt = 0.01
    for i=1:N
        for ii = 1:length(F)
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt)*rand(gaussianD)
        end
    end
    return Flist
end
function Langevin(N,a,m,mu,la,gaussianD,filename)
    path = "results/"
    F = [20. for i = 1:16]
    dt = 0.01
    n_burn = 3/dt
    n_skip = 3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F)
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt)*rand(gaussianD)
        end
    end
    show(F);println()

    # Simulate
    t1 = @timed for i=1:N
        writec123tofile(path,filename,F,i)
        for iii = 1:n_skip+1
            for ii = 1:length(F)
                F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt)*rand(gaussianD)
            end
        end
        if i%2000==0
            println(i)
        end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return
end
println()
Langevin(20,0.5,1,1,0,gaussianD)
begin
    n_tau=16
    β=8
    a=β/n_tau
    Langv1=Langevin(15000,a,1,1,0,gaussianD,"CL_1.csv")
end

    # Probability density diagram #
# arr1 = [Langv1[i][j] for i = 1:length(Langv1) for j = 1:length(Langv1[1])]
# histogram(arr1,normed=true,xlabel="x",ylabel="|ψ_0|²")
# histogram(Langv1,normed=true,xlabel="x",ylabel="|ψ_0|²")
PlotProbDD("results/CL_1.csv",0.1)
PlotProbDDe(1,1,1,3)
    # sampling
plot(reshape(GetColumn(2,"results/CL_1.csv"),:))#:Int((length(LastRowFromFile(file))-1)/4)+1
    # Autocorrelation
PlotAC("results/CL_1.csv",1000)
# PlotACsb("results/CL_1.csv",1000)
# PlotAC("results/CL_1.csv",false)
# PlotAC("results/CL_1.csv",true)
    # Twopoint Correlation
PlotTPCF("results/CL_1.csv")                # Naive error
PlotTPCF("results/CL_1.csv",true)           # For autocorrelated data
a = [0.990894    0.00115783;
 -0.0086682   0.000855048;
  0.00979547  0.000858429;
 -0.0109495   0.00088357;
 -0.0096282   0.000881336;
 -9.15654e-6  0.000868927;
  0.00361841  0.000861336;
  0.00736314  0.000860249;
  0.0157388   0.00113745;
  0.00736314  0.000860249;
  0.00361841  0.000861336;
 -9.15654e-6  0.000868927;
 -0.0096282   0.000881336;
 -0.0109495   0.00088357;
  0.00979547  0.000858429;
 -0.0086682   0.000855048;]
 plot(a[:,1],yerr=a[:,2])
 PlotEffM("results/CL_1.csv")