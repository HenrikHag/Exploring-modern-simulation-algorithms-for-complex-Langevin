
begin
    # using Distributions
end

export AHO_L_param_old, Simulation_param_old, Langevin_old#, AHO_L_Der_old




# Langevin functions

struct AHO_L_param_old
    β::Float64      # Inverse temperature
    n_tau::Integer  # Lattice points
    m::Float64      # mass
    μ::Float64      # harmonic coupling
    λ::Float64      # anharmonic coupling
end

struct Simulation_param_old
    N::Integer
    n_burn::Integer
    n_skip::Integer
    dt::Real
end

"""
Derivative of discretized AHO (a,m,μ,λ)∈R action  
`S = a(1/2*m*(P[i+1]-P[i])²/a² + μ/2*ϕ²) + λ*ϕ⁴/(4!*a⁴)`  
`∂S/∂ϕⱼ = m*(2*P[i]-P[i+1]-P[i-1])/a² + μϕ + λ*ϕ³/(6*a³)`  
`μ = m*ω²`
"""
function AHO_L_Der_old(a,m,mu,la,F,f₋₁,f₊₁)
    # return m/a*(2*F-(f₋₁+f₊₁)) + a*mu*F + la*F^3/(6*a^4)
    return m*(2*F-(f₋₁+f₊₁))/a^2 + mu*F + la*F^3/(6*a^3)
end

"""
Simulates an Euclidean time Anharmonic Oscillator  
If provided a filename, stores the configurations there instead
"""
function Langevin_old(phys_param::AHO_L_param_old,sim_param::Simulation_param_old,gaussianD)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.β/n_tau
    m =  phys_param.m
    mu = phys_param.μ
    la = phys_param.λ
    # Initial condition
    F = [20. for i = 1:n_tau]
    # Simulation parameters
    N = sim_param.N
    n_burn = sim_param.n_burn
    n_skip = sim_param.n_skip
    dt = sim_param.dt
    
    # Burn in
    randoms1 = rand(gaussianD,(N*n_skip+n_burn)*n_tau)
    for i=1:n_burn
        for ii = 1:n_tau
            F[ii] -= AHO_L_Der_old(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*randoms1[n_tau*(i-1)+ii]
        end
    end
    
    # Simulate
    Flist = Matrix{Float64}(undef,N,n_tau)#fld(N,n_skip),n_tau)
    # Flist[1,:] = F#; println(Flist)
    for i=0:N*n_skip-1
        # println(F)
        for ii = 1:n_tau
            F[ii] -= AHO_L_Der_old(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*randoms1[n_tau*(i+1)+ii]
        end
        if i%n_skip == 0 # Save once every n_skip-1
            Flist[div(i,n_skip)+1,:] = F
        end
    end
    return Flist
end
function Langevin_old(a,m,mu,la,sim_param,gaussianD,filename::AbstractString)
    # Initial condition
    F = [20. for i = 1:n_tau]
    dt = sim_param.dt
    # timespan = (0.0,dt)
    n_burn = 3/dt
    n_skip = 3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F)
            # Other solvers
            # ϕ₋₁ = F[(ii-2+n_tau)%n_tau+1]; ϕ₊₁ = F[(ii)%n_tau+1]; ϕ0 = F[ii]
            # f(ϕ,t,p) = (m/a*(ϕ^2-ϕ*(ϕ₊₁+ϕ₋₁)) + 0.5*m*mu*a*ϕ^2)*dt
            # prob = ODEProblem(f,ϕ0,timespan)
            # sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
            # F[ii] -= sol(dt)*dt - sqrt(2*dt/a)*rand(gaussianD)
            # Self implemented explicit Euler solver
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
        end
    end
    show(F);println()

    # Simulate
    t1 = @timed for i=1:N
        writec123tofile(save_path,filename,F,i)
        for iii = 1:n_skip+1
            for ii = 1:length(F)
                # Other solvers
                # ϕ₋₁ = F[(ii-2+n_tau)%n_tau+1]; ϕ₊₁ = F[(ii)%n_tau+1]; ϕ0 = F[ii]
                # f(ϕ,t,p) = (m/a*(ϕ^2-ϕ*(ϕ₊₁+ϕ₋₁)) + 0.5*m*mu*a*ϕ^2)*dt
                # prob = ODEProblem(f,ϕ0,timespan)
                # sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
                # F[ii] -= sol(dt)*dt - sqrt(2*dt/a)*rand(gaussianD)
                # Self implemented explicit Euler solver
                F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
            end
        end
        if i%2000==0
            println(i)
        end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return
end

# res1 = Langevin_old(100,0.5,1,1,0,gaussianD)








