begin
    using MCMC
    using StochasticDiffEq, DifferentialEquations
    using LabelledArrays
    using BenchmarkTools
    using Plots
    using Base
end

struct AHO_CL_Param
    a::Real
    m::Real
    mu::Complex
    la::Real
end

begin
    a=0.5; m=1; mu = exp(im*π/3.); la = 0;
    n_tau = 16;
    # function CLangevinSchem(N,a,m,mu,la)
    F0 = append!([20. for i = 1:n_tau],[0. for i = 1:n_tau])
    # Flist = Matrix{Float64}(undef,N+1,n_tau)
    # Flist[1,:] = F0
    dt = 1e-4;
    timespan = (0.0,10)
    # params = Lvector(p=struct w fields m μ λ)
    params = LVector(p=AHO_CL_Param(a,m,mu,la))
    # Function to calculate change in action for whole path
    sdeprob1 = SDEProblem(ActionCLDerSchem,RandScale,F0,timespan,params)

    @time sol = solve(sdeprob1, Euler(), saveat=0.01, save_start=false,
                dtmax=1e-3, dt=dt, abstol=5e-2, reltol=5e-2)
end

hcat(sol.u[1:end]...)
plot(sol,label="")
savefig("plots/$(findDate())_CL_SDESolver_thermalization.pdf")
savefig("plots/$(findDate())_CL_SDESolver_thermalization.png")
cov(sol)
# Euler() 1st time from clear 64.7 s
# Euler() 2nd time from clear 0.56 s
