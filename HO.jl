## Calculate expectationvalues of the Harmonic Oscillator ##

import math

Nt = 12
λ = 24
β = 1
m = 1
dt = β/Nt
Basis = 32

HAMOC = zeros(Float64, Basis, Basis)

for i = 1:Basis
    for j = 1:Basis
        if i == j
            HAMOC[i,j] = m*(i+1/2)
        end
    end
end

HAMOC

function HTEvol(t)
    return exp(-im*t*HAMOC)
end

ex1 = HTEvol(34)

HBoltzmann = exp(-β*HAMOC)

PartFunc = tr(HBoltzmann)