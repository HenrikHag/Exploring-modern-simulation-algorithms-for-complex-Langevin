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
println(Langevin(20,1,0.4,gaussianD))
histogram(Langevin(100000,1,0.4,gaussianD)[200:end],normed=true,xlabel="x",ylabel="|ψ_0|²")
PlotProbDDe(1/2,1,1,3)
