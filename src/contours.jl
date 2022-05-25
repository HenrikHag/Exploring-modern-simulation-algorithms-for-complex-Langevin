# Defining contours onto which the simulation of systems occur

using BenchmarkTools
using Plots

struct EuclideanContour
    β::Real # Inverse temperature of system
end

function DiscretizeContour(C::EuclideanContour,N)
    return [i for i=0:1/N:1]
end

C1 = EuclideanContour(8)
N_tau = 16
DC1 = DiscretizeContour(C1,N_tau)
a = im*(DC1[1]-DC1[2])
CC = [-im*DC1[i] for i=1:length(DC1)]
scatter(CC,legend=false)

@benchmark DiscretizeContour($C1,$N_tau) # 188ns

struct EucledianContourDiscretization
    NR_Points::Integer
end

CD1 = EucledianContourDiscretization(N_tau)

function spread_timepoints(CD::EucledianContourDiscretization,
    C::EuclideanContour)
return collect(range(0,stop=1,length=CD.NR_Points+1))
end

@benchmark spread_timepoints($CD1,$C1) # 378 ns
