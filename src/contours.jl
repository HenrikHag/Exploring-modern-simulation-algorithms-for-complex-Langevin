# Defining contours onto which the simulation of systems occur

using BenchmarkTools
using Plots
using LabelledArrays
import Plots.scatter

struct EuclideanContour
    β::Real # Inverse temperature of system
end

struct HO_Param
    m::Float64
    σ::Float64
end
struct D_HO_Param
    a::Float64
    m::Float64
    σ::Float64
end

struct AHO_Param
    m::Float64  # mass
    μ::Float64  # harmonic coupling
    λ::Float64  # anharmonic coupling
end
struct D_AHO_Param
    a::Float64  # Lattice spacing
    m::Float64  # mass
    μ::Float64  # harmonic coupling
    λ::Float64  # anharmonic coupling
end

struct AHO_CL_Param
    m::Float64
    μ::Complex
    λ::Float64
end
function Param(p::AHO_CL_Param,a::AbstractArray)
    return LVector(p=p, a=a, vals=fill(zeros(length(a)),6),
                   as=(circshift(a,1),circshift(a,-1),circshift(a,2)))
end



function DiscretizeContour(C::EuclideanContour,N)
    return [i for i=0:1/N:1]
end

# C1 = EuclideanContour(8)
# N_tau = 16
# DC1 = DiscretizeContour(C1,N_tau)
# a = im*(DC1[1]-DC1[2])
# CC = [-im*DC1[i] for i=1:length(DC1)]
# scatter(CC,legend=false)

# @benchmark DiscretizeContour($C1,$N_tau) # 188ns

# struct EucledianContourDiscretization
#     NR_Points::Integer
# end

# CD1 = EucledianContourDiscretization(N_tau)

function ContourDistances(CD,C::EuclideanContour)
    return [CD[i+1]-CD[i] for i=1:length(CD)-1]
end

function Contour(CD,C::EuclideanContour)
    return [-im * u for u in CD] .* C.β
end

# function ContourDerivative()
    
# end

# function spread_timepoints(CD::EucledianContourDiscretization,
#     C::EuclideanContour)
#     return collect(range(0,stop=1,length=CD.NR_Points+1))
# end
# @benchmark spread_timepoints($CD1,$C1) # 378 ns

function scatter(discretizedContour,C::EuclideanContour)
    scatter(imag.(discretizedContour),real.(discretizedContour))
end

function a_AHO_r(du,u,param,t)
    p = param.p
    a = param.a
    a_m1, a_p1, _ = param.as
    xR = u[1:div(end,2)]
    xI = u[div(end,2)+1:end]
    ###### Copying copying copying 😰
end

function b_AHO_r(du,u,param,t)
    a = param.a
    a_m1, a_p1, _ = param.as
    du[1:div(end,2)] .= sqrt.(2 ./(0.5 .*(abs.(a)+abs.(a_m1))))
end

function SDEfunc(::AHO_Param)
    return a_AHO_r, b_AHO_r, false
end

begin # main
    p = AHO_CL_Param(1,1,24)
    contour = EuclideanContour(8)
    t_steps = 16
    discretizedContour = DiscretizeContour(contour,t_steps)
    a = ContourDistances(discretizedContour,contour)
    tp = Contour(discretizedContour,contour)
    scatter(discretizedContour,contour)
end

## simulation
args = (N_tr=10, θ=0.6, dt=1e-4, tol= 5e-2, dtmax=1e-3, tspan=(0.,10.), adaptive=true)
begin
    params = Param(p,a)

    sde_a, sde_b, _ = SDEfuncs(p)

    sdeProblem = SDEProblem(sde_a,sde_b,y0,args.tspan,params)
    ensembleProblem = EnsembleProblem(sdeProblem)

    @time sol = solve(ensembleProblem, ImplicitEM(theta=args.θ, symplectic=false),
                EnsembleThreads(), trajectories=args.N_tr,
                progress=true, saveat=0.01, savestart=false,
                dtmax=args.adaptive, abstol=args.tol, reltol=args.tol)
end