
# begin
#     # using Distributions
# end

export AHO_L_param, getAHO_L_param, Sim_L_param, getSim_L_param, Langevin_AHO #, AHO_L_ActionDer

export Gaussian_CL_param, getGaussian_CL_param, Sim_CL_param, getSim_CL_param, CLangevin_Gauss

export AHO_CL_param, getAHO_CL_param, CLangevin_AHO #, AHO_CL_ActionDer

# export AHO_CL_param
# export ActionCLDerSchem, RandScale, CLangevinSchem



#                       #
#  Langevin functions   #
#                       #
"""
Physical parameters for Langevin simulation
"""
struct AHO_L_param
    n_tau::Integer
    a::Real
    m::Real
    μ::Real
    λ::Real
end

"""
returns the AHO parameters for Langevin simulation
"""
function getAHO_L_param(n_tau::Integer,β::Real,m::Real,μ::Real,λ::Real)
    a = β/n_tau
    return AHO_L_param(n_tau,a,m,μ,λ)
end


"""
Simulation parameters for Langevin simulation
"""
struct Sim_L_param
    N::Integer
    n_burn::Integer
    n_skip::Integer
    dt::Real
end

"""
returns the simulation parameters for Langevin simulation
"""
function getSim_L_param(N::Integer,n_burn::Integer,n_skip::Integer,dt::Real)
    if n_skip < 1
        println("n_skip must be 1 or higher\nEvery n_skip-1 sample is discarded")
        return
    end
    if dt < 0 || N < 1 || n_burn < 0
        println("Parameters must be 0 or higher (N 1 or higher)")
        return
    end
    return Sim_L_param(N,n_burn,n_skip,dt)
end

# struct AHO_L_param_old
#     β::Float64      # Inverse temperature
#     n_tau::Integer  # Lattice points
#     m::Float64      # mass
#     μ::Float64      # harmonic coupling
#     λ::Float64      # anharmonic coupling
# end

# struct Simulation_param_old
#     N::Integer
#     n_burn::Integer
#     n_skip::Integer
#     dt::Real
# end

"""
Derivative of discretized AHO (a,m,μ,λ)∈R action  
`S = a(1/2*m*(P[i+1]-P[i])²/a² + μ/2*ϕ²) + λ*ϕ⁴/(4!*a⁴)`  
`∂S/∂ϕⱼ = m*(2*P[i]-P[i+1]-P[i-1])/a² + μϕ + λ*ϕ³/(6*a³)`  
`μ = m*ω²`
"""
function AHO_L_ActionDer(a,m,mu,la,F,f₋₁,f₊₁)
    # return m/a*(2*F-(f₋₁+f₊₁)) + a*mu*F + la*F^3/(6*a^4)
    return m*(2*F-(f₋₁+f₊₁))/a^2 + mu*F + la*F^3/(6*a^3)
end

"""
Simulates an Euclidean time Anharmonic Oscillator  
If provided a filename, stores the configurations there instead
"""
function Langevin_AHO(phys_param::AHO_L_param,sim_param::Sim_L_param,gaussianD)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.a
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
    # randoms1 = rand(gaussianD,n_burn*n_tau)
    for i=1:n_burn
        for ii = 1:n_tau
            F[ii] -= AHO_L_ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)#randoms1[n_tau*(i-1)+ii]
        end
    end
    
    # Simulate
    Flist = Matrix{Float64}(undef,N,n_tau)#fld(N,n_skip),n_tau)
    # Flist[1,:] = F#; println(Flist)
    for i=1:N
        # println(F)
        Flist[i,:] = F
        for iii=1:n_skip
            for ii = 1:n_tau
                F[ii] -= AHO_L_ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)#randoms1[n_tau*(i+1)+ii]
            end
        end
    end
    return Flist
end
function Langevin_AHO(phys_param::AHO_L_param,sim_param::Sim_L_param,gaussianD,filename::AbstractString)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.a
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
    for i=1:n_burn
        for ii = 1:n_tau
            # Self implemented explicit Euler solver
            F[ii] -= AHO_L_Der_old(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
        end
    end
    # show(F);println()

    # Simulate
    t1 = @timed for i=1:N
        writec123tofile(filename,F,i)
        for iii = 1:n_skip
            for ii = 1:n_tau
                # Self implemented explicit Euler solver
                F[ii] -= AHO_L_Der_old(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
            end
        end
        # if i%2000==0
        #     println(i)
        # end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return
end









#                               #
#  Complex Langevin functions   #
#                               #

# For the Gaussian system: S = 1/2 μ ϕ²
# Then the imaginary part drifts as:
#       ϕᵢⁱ⁺¹ = ϕᵢⁱ - m μ ϕᵢⁱ dt
"""
Physical parameters for Langevin simulation
"""
struct Gaussian_CL_param
    μ::Real
end


"""
returns the gaussian system parameters for Langevin simulation
"""
function getGaussian_CL_param(μ::Number)
    if isa(μ,Type{Complex})
        return Gaussian_CL_param(μ)
    end
    μ = μ + im*0
    return Gaussian_CL_param(μ)
end

"""
Simulation parameters for Langevin simulation
"""
struct Sim_CL_param
    N::Integer
    n_burn::Integer
    n_skip::Integer
    dt::Real
end

"""
returns the simulation parameters for complex Langevin simulation
"""
function getSim_CL_param(N::Integer,n_burn::Integer,n_skip::Integer,dt::Real)
    if n_skip < 1
        println("n_skip must be 1 or higher\nEvery n_skip-1 sample is discarded")
        return
    end
    if dt < 0 || N < 1 || n_burn < 0
        println("Parameters must be 0 or higher (N 1 or higher)")
        return
    end
    return Sim_CL_param(N,n_burn,n_skip,dt)
end

"""
Complex Langevin for S = - μϕ²
returns list of real and list of imaginary parts of ϕ
"""
function CLangevin_Gauss(phys_param::Gaussian_CL_param,sim_param::Sim_CL_param,gaussianD)
    # Physical parameters
    mu_r = real(phys_param.μ)
    mu_i = imag(phys_param.μ)
    # Simulation parameters
    N = sim_param.N
    n_burn = sim_param.n_burn
    n_skip = sim_param.n_skip
    dt = sim_param.dt
    # Initial conditions
    # [x₁ʳ,x₂ʳ,...]
    F_r = 0.1 #[20. for i = 1:n_tau]
    # [x₁ᶜ,x₂ᶜ,...]
    F_i = 0. #[1. for i = 1:n_tau]

    # Thermalize
    for i=1:n_burn
        # derAction = CActionDer(a, m, mu, la, F_r, F_i)
        # Forward Euler: F_{n+1} = F_{n} + f(t_n,F_{n})
        F_r = F_r - (mu_r*F_r + mu_i*F_i)*dt + sqrt(2*dt)*rand(gaussianD) # N_R = 1 => N_I = 0
        F_i = F_i - (mu_i*F_r + mu_r*F_i)*dt
        # F_r[ii] -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD)
        # F_i[ii] -= derAction[2]*dt
    end
    # show(F_r);println()

    Flist_r = Array{Float64}(undef,N)
    Flist_i = Array{Float64}(undef,N)
    # push!(Flist_r,F_r)
    # push!(Flist_i,F_i)

    # Simulate
    t1 = @timed for i=1:N
        Flist_r[i] = F_r
        Flist_i[i] = F_i
        # writeto(filename,[F_r,F_i],i)
        for iii = 1:n_skip
            # derAction = CActionDer(a,m,mu,la,F_r[ii],F_r[(ii-2+n_tau)%n_tau+1],F_r[(ii)%n_tau+1,F_c[ii],F_c[(ii-2+n_tau)%n_tau+1],F_c[(ii)%n_tau+1]])
            F_r = F_r - (mu_r*F_r + mu_i*F_i)*dt + sqrt(2*dt)*rand(gaussianD) # N_R = 1 => N_I = 0
            F_i = F_i - (mu_i*F_r + mu_r*F_i)*dt
            # F_r[ii] -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD)
            # F_i[ii] -= derAction[2]*dt
        end
        # if i%2000==0
        #     println(i)
        # end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return Flist_r, Flist_i#WeightP_r, WeightP_i, WeightN_r, WeightN_i, 
end









# HS transform gives an equation where fermions are integrated out by an auxiliary field
# function AHO_CL_HS_ActionDer(a,m,μ,λ,F_r,F_i,ii)
#     # z = m*μ*(Fᵣ+im*Fᵢ) + m*λ/6*(Fᵣ+im*Fᵢ)^3 #+ m/a^2*(2*F-(f₋₁+f₊₁))
#     # z = λ/12*(Fᵣ + im*Fᵢ) + im*λ/(12μ+2*im*λ*(Fᵣ+im*Fᵢ))
#     return real.(z),imag.(z)
# end
# """
# Calculated real and complex part of the weight term in Eq.(110)
# exp( - S_B(σ) )
# """
# function Weight_func(μ,λ,Fᵣ,Fᵢ)
#     z = λ/24*(Fᵣ+im*Fᵢ)^2-0.5*log(λ/(12*μ+2*im*λ*(Fᵣ+im*Fᵢ)))
#     return real(exp(-z)), imag(exp(-z)) # e⁻ˢ
# end
# Weight_func(1,0.4,1,0)
    # WeightP_r = []
    # WeightP_i = []
    # WeightN_r = []
    # WeightN_i = []
    # weight_c = 0# Weight_func(mu,la,F_r,F_i)
    # push!(Weight_r,weight_c[1])
    # push!(Weight_i,weight_c[2])
    # weight_c = Weight_func(mu,la,F_r,F_i)
    # if F_r > 0
    #     push!(WeightP_r,weight_c[1])
    #     push!(WeightP_i,weight_c[2])
    # else
    #     push!(WeightN_r,weight_c[1])
    #     push!(WeightN_i,weight_c[2])
    # end
    

# Compute ⟨ϕ²⟩ = ⟨(ϕ_r + Im*ϕ_i)^2⟩ᵩ
# Instead of doing a Hubbard-Stratonovich transformation described in paper [2], 
#   we just do the mentioned choice of complex μ (or λ).

"""
Physical parameters for Langevin simulation
"""
struct AHO_CL_param
    n_tau::Integer
    a::Real
    m::Real
    μ::Complex
    λ::Real
end

"""
returns the AHO parameters for complex Langevin simulation
"""
function getAHO_CL_param(n_tau::Integer,β::Real,m::Real,μ::Complex,λ::Real)
    a = β/n_tau
    if isa(μ,Type{Complex})
        return AHO_CL_param(n_tau,a,m,μ,λ)
    end
    μ = μ + im*0
    return AHO_CL_param(n_tau,a,m,μ,λ)
end

"""
Complex and real part of the derivative of the action wrt. ϕ
"""
function AHO_CL_ActionDer(a,m,μ,λ,F_r,F_i,ii)
    n_tau = length(F_r)
    # return m/a*(2*F-(f₋₁+f₊₁)) + a*mu*F + la*F^3/(6*a^4)
    # m*(2*F-(f₋₁+f₊₁))/a^2 + mu*F + la*F^3/(6*a^3)
    z_r = (m*(2*F_r[ii]-(F_r[(ii-2+n_tau)%n_tau+1]+F_r[(ii)%n_tau+1]))/a^2 + real(μ)*F_r[ii] + imag(μ)*F_i[ii] + λ*F_r^3/(6*a^3))
    z_i = (m*(2*F_i[ii]-(F_i[(ii-2+n_tau)%n_tau+1]+F_i[(ii)%n_tau+1]))/a^2 + real(μ)*F_i[ii] + imag(μ)*F_r[ii] + λ*F_i^3/(6*a^3))
    return z_r, z_i
end


function CLangevin_AHO(phys_param::AHO_CL_param,sim_param::Sim_CL_param,gaussianD)
    # Physical parameters
    n_tau = phys_param.n_tau
    a = phys_param.a
    m =  phys_param.m
    mu = phys_param.μ
    la = phys_param.λ
    # Initial condition
    F_r = [1. for i = 1:n_tau]
    F_i = [0.1 for i = 1:n_tau]
    # Simulation parameters
    N = sim_param.N
    n_burn = sim_param.n_burn
    n_skip = sim_param.n_skip
    dt = sim_param.dt
    
    # Burn in
    # randoms1 = rand(gaussianD,n_burn*n_tau)
    for i=1:n_burn
        for ii = 1:n_tau
            dS_r, dS_i = AHO_L_ActionDer(a,m,mu,la,F_r,F_i,ii)
            F_r[ii] = F_r[ii] - dS_r * dt - sqrt(2*dt/a)*rand(gaussianD)# randoms1[n_tau*(i-1)+ii]
            F_i[ii] = F_i[ii] - dS_i * dt
        end
    end
    
    # Simulate
    F_r_list = Matrix{Float64}(undef,N,n_tau)#fld(N,n_skip),n_tau)
    F_i_list = Matrix{Float64}(undef,N,n_tau)#fld(N,n_skip),n_tau)
    # Flist[1,:] = F#; println(Flist)
    # randoms1 = rand(gaussianD,N*n_tau)
    for i=1:N
        # println(F_r)
        F_r_list[div(i,n_skip)+1,:] = F_r
        F_i_list[div(i,n_skip)+1,:] = F_i
        for iii = 1:n_skip
            for ii = 1:n_tau
                dS_r, dS_i = AHO_L_ActionDer(a,m,mu,la,F_r,F_i,ii)
                F_r[ii] = F_r[ii] - dS_r * dt + sqrt(2*dt/a)*rand(gaussianD)# randoms1[n_tau*(i-1)+ii]
                F_i[ii] = F_i[ii] - dS_i * dt
                # F[ii] -= AHO_L_ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*randoms1[n_tau*(i+1)+ii]
            end
        end
    end
    return F_r_list, F_i_list
end




















######################################
## Complex Langevin solver package ###
######################################

# """
# returns the array [dx_r...,dx_i...], the derivative of the action for StochasticDiffEq  
# `ϕ = ∑ᵢ 1/2 a[m ((ϕᵢ₊₁-ϕᵢ)/a)² + μϕᵢ²]`  
# `ϕⱼ = 1/2 a[m ((ϕⱼ₊₁-ϕⱼ)/a)² + ((ϕⱼ-ϕⱼ₋₁)/a)² + μϕⱼ²]`  
# `∂ϕ/∂ϕⱼ = ∂/∂ϕⱼ a [m (ϕⱼ²-(ϕⱼ)(ϕⱼ₊₁+ϕⱼ₋₁))/2a² + μ/2 ϕⱼ²]`  
# `       = [m (2ϕⱼ - ϕⱼ₊₁ - ϕⱼ₋₁)/a² + μϕⱼ]`
# """
# function ActionCLDerSchem(du, u, params, t)
#     p = params.p
#     xR = @view u[1:div(end,2)]
#     xI = @view u[div(end,2)+1:end]
#     Fr_diff_m1 = xR .- xR[vcat(end,1:end-1)]   # dx_j - dx_{j-1}
#     Fr_diff_p1 = xR[vcat(2:end,1)] .- xR       # dx_{j+1} - dx_j
#     Fi_diff_m1 = xI .- xI[vcat(end,1:end-1)]   # dx_j - dx_{j-1}
#     Fi_diff_p1 = xI[vcat(2:end,1)] .- xI       # dx_{j+1} - dx_j
#     # dx_j - dx_{j-1} - (dx_{j+1} - dx_j) = 2dx_j - dx_{j+1} - dx_{j-1}
#     du[1:div(end,2)] .= (p.m .* real.(Fr_diff_p1 .- Fr_diff_m1 .+ im .* (Fi_diff_p1 .- Fi_diff_m1)) ./ p.a^2) .- real.(p.mu .* xR .+ im .* (p.mu .* xI))
#     du[div(end,2)+1:end] .= (p.m .* imag.(im .* (Fi_diff_p1 .- Fi_diff_m1) .+ (Fr_diff_p1 .- Fr_diff_m1)) ./ p.a^2) .- imag.(p.mu .* xR .+ im .* (p.mu .* xI))
# end

# function RandScale(du, u, param, t)
#     a = param.p.a
#     du[1:div(end,2)] .= sqrt.(2. ./ a)
# end

# function CLangevinSchem(N,a,m,mu,la)
#     n_tau = 16
#     F0 = append!([20. for i = 1:n_tau],[0. for i = 1:n_tau])
#     # Flist = Matrix{Float64}(undef,N+1,n_tau)
#     # Flist[1,:] = F0
#     dt = 0.01
#     timespan = (0.0,3*N)
#     # params = Lvector(p=struct w fields m μ λ)
#     params = LVector(p=AHO_CL_Param(a,m,mu,la))
#     # Function to calculate change in action for whole path
#     sdeprob1 = SDEProblem(ActionCLDerSchem,RandScale,F0,timespan,params)

#     @time sol = solve(sdeprob1, Euler(), progress=true, saveat=0.01/dt, save_start=false,
#                 dtmax=1e-3, dt=dt, abstol=5e-2,reltol=5e-2)
# end





# function writeto(filename, Path, itt::Integer)
#     open(filename,"a") do file
#         # write(file,"c_x,c_x2,c_x0x1\n")
#         write(file,string(itt,","))
#         pathl=length(Path)
#         for i = 1:pathl-1         # Path
#             write(file,string(Path[i],","))
#         end
#         write(file,string(Path[end],"\n"))
#     end
# end