begin
    using MCMC
    using StochasticDiffEq, DifferentialEquations
    using LabelledArrays
    using BenchmarkTools
    using Plots
    using Base
    save_date = findDate() # Save date for plot names
    save_folder = "plots/" # Where to store plots
end

struct AHO_CL_Params
    n_tau::Integer
    a::Real
    m::Real
    mu::Complex  # (m ω²)
    la::Real
end
function getAHO_CL_Params(n_tau,β,m,mu,la)
    a = β/n_tau
    AHO_CL_Params(n_tau,a,m,mu,la)
end

begin # Prepare SDEProblem
    n_tau=16; β=8; m=1; μ=1. +0*im +0*exp(im*π/3.); λ=0;
    # params = Lvector(p=struct w fields m μ λ)
    params = LVector(p=getAHO_CL_Params(n_tau,β,m,μ,λ))
    # function CLangevinSchem(N,a,m,mu,la)
    F0 = append!([20. for i = 1:n_tau],[0. for i = 1:n_tau])
    # Flist = Matrix{Float64}(undef,N+1,n_tau)
    # Flist[1,:] = F0
    dt = 1e-4;
    timespan = (0.0,10.)
    # Generate the SDE to be solved
    sdeprob1 = SDEProblem(ActionCLDerSchem,RandScale,F0,timespan,params)
end

begin # Simulations
    # @time sol = solve(sdeprob1, Euler(), saveat=0.01, save_start=false,
    #             dtmax=1e-3, dt=dt, abstol=5e-2, reltol=5e-2);
    @time sol = solve(sdeprob1, EM(), saveat=0.01, save_start=false,
                dtmax=1e-3, dt=dt, abstol=5e-2, reltol=5e-2);
    @time sol = solve(sdeprob1, ImplicitEM(theta=0, symplectic=false),
                saveat=0.01, save_start=false, adaptive=true,
                dtmin=1e-5, dtmax=1e-3, dt=dt, abstol=5e-2, reltol=5e-2);
    @time sol = solve(sdeprob1, ImplicitEM(theta=1, symplectic=false),
                saveat=0.01, save_start=false, adaptive=true,
                dtmin=1e-5, dtmax=1e-3, dt=dt, abstol=5e-2, reltol=5e-2);
end;
# Euler() 1st time from clear; tspan(0,10), dt=1e-4; 64.7 s
# Euler() 2nd time from clear; tspan(0,10), dt=1e-4; 0.56 s

# timespan = (0.,1000.)
# 71.8 s
# 191. s
# 34.3 s # Error
# 1010 s

# plot(sol,label="")
# savefig("$(save_folder)$(save_date)_CL_SDESolver_thermalization.pdf")
# savefig("$(save_folder)$(save_date)_CL_SDESolver_thermalization.png")
sol(1:2)
Therm = 1000 # Throw away unthermalized configurations
# length(sol.u)
n_skip = 100  # Throw away every n_skip to reduce correlation
# Real part of ϕ without Therm first elements and every n_skip intermediate elements
sol1 = transpose(hcat(sol.u[Therm:n_skip:end]...))[:,1:n_tau]
plot(sol1)
# Autocorrelation
PlotAC(sol1,1000)
# savefig("$(save_folder)$(save_date)_CL_SDESolver_AC.pdf")
# savefig("$(save_folder)$(save_date)_CL_SDESolver_AC.png")
# Twopoint Correlation
PlotTPCF(sol1)
PlotTPCFe!(params.p.a,m,sqrt(real(μ)/m),n_tau)
# PlotTPCFe!(params.p)
# savefig("$(save_folder)$(save_date)_CL_SDESolver_TPCF.pdf")
# savefig("$(save_folder)$(save_date)_CL_SDESolver_TPCF.png")
PlotProbDD(sol1)
PlotProbDDe(m,sqrt(real(μ)/m),1,2)

# Taking the last configuration from last simulation as initial condition
sol.u[end]