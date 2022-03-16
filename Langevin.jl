using Distributions
using Plots
using .UsersGuide
gaussianD = Normal(0,1)

#       μ/2*ϕ²+λ/4!*ϕ⁴
#       -> (μϕ+λ/6*ϕ³)*Δt
function ActionDer(a,m,mu,la,F,f₋₁,f₊₁)
   return m/a*(F-1/2*(f₋₁+f₊₁)) + mu*F + la/6*F^3
end
# Add cupling terms 2f(i)-f(i+1)-f(i-i)
# Understand the discretizing integral and meeting mat. from 16.03

function Langevin(N,mu,la,gaussianD)
    F = [20. for i = 1:16]
    Flist = []
    dt = 0.01
    for i=1:N
        push!(Flist,F)
        for ii = 1:length(F)
            F[ii] -= Action(mu,la,F[ii],dt) - sqrt(2*dt)*rand(gaussianD)
        end
    end
    return Flist
end
function Langevin(N,a,m,mu,la,gaussianD,filename)
    path = "results/"
    F = [20. for i = 1:16]
    Flist = []
    dt = 0.01
    n_burn = 200
    n_skip = 200
    for i=1:n_burn
        for ii = 1:length(F)
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt)*rand(gaussianD)
        end
    end
    show(F);println()
    for i=1:N
        push!(Flist,F)
        open(string(path,filename),"a") do file
            # write(file,"exp_x,exp_x2,exp_x0x1\n")
            write(file,string(Int64(i),","))
            pathl=length(F)
            for ii = 1:pathl
                write(file,string(F[ii],","))
            end
            # write(file,string(F[end],"\n"))
            for ii = 1:pathl         # Path.^2
                write(file,string(F[ii]^2,","))
            end
            for ii = 1:pathl         # Path_1*Path_i
                write(file,string(F[1]*F[ii],","))
            end
            for ii = 0:pathl-2       # Two-Point Correlation
                twopointcorr=0
                for iii=1:pathl
                    twopointcorr += F[iii]*F[(iii+ii-1)%pathl+1]
                end
                write(file,string(twopointcorr/pathl,","))
            end
            twopointcorr=0
            for ii=1:pathl
                twopointcorr += F[ii]*F[(ii+pathl-2)%pathl+1]
            end
            write(file,string(twopointcorr/pathl,"\n"))
        end
        for iii = 1:n_skip+1
            for ii = 1:length(F)
                F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt)*rand(gaussianD)
            end
        end
    end
    return Flist
end
# println(Langevin(20,1,0.4,gaussianD,"CL_1.csv"))

begin
    n_tau=16
    β=8
    a=β/n_tau
    Langv1=Langevin(15000,a,1,1,0,gaussianD,"CL_1.csv")
end

# Probability density diagram #
arr1 = [Langv1[i][j] for i = 1:length(Langv1) for j = 1:length(Langv1[1])]
histogram(arr1,normed=true,xlabel="x",ylabel="|ψ_0|²")
histogram(Langv1,normed=true,xlabel="x",ylabel="|ψ_0|²")
PlotProbDD("results/CL_1.csv",0.1)
PlotProbDDe(1,1,1,3)
# sampling
plot(reshape(GetColumn(2,"results/CL_1.csv"),:))#:Int((length(LastRowFromFile(file))-1)/4)+1
# Autocorrelation
PlotAC("results/CL_1.csv",1000)
PlotACsb("results/CL_1.csv",1000)
PlotAC("results/CL_1.csv",false)
PlotAC("results/CL_1.csv",true)
# Twopoint Correlation
PlotTPCF("results/CL_1.csv")
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