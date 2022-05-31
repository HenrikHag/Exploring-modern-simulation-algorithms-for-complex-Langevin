# Defining contours onto which the simulation of systems occur

using BenchmarkTools
using Plots 

struct EuclideanContour
    β::Real # Inverse temperature of system
end

struct HO_Param
    m::Float64
    σ::Float64
end

struct AHO_Param
    m::Float64  # mass
    μ::Float64  # harmonic coupling
    λ::Float64  # anharmonic coupling
end

struct AHO_CL_Param
    m::Float64
    μ::Complex
    λ::Float64
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

begin # main
    p = AHO_CL_Param(1,1,24)
    contour = EuclideanContour(8)
    t_steps = 16
    discretizedContour = DiscretizeContour(contour,t_steps)
    a = ContourDistances(discretizedContour,contour)
    tp = Contour(discretizedContour,contour)
    scatter(discretizedContour)
end