begin
    using Distributions
    using Plots
    using .UsersGuide
    using .MetropolisUpdate
    using DifferentialEquations, StochasticDiffEq#, SimpleDiffEq
    using LabelledArrays
    save_path = "results/"
    gaussianD = Normal(0,1)
end
#       a(1/2*m*(P[i+1]-P[i])²/a² + m*μ/2*ϕ² + m*λ/4!*ϕ⁴)
#       -> (m/a²*(2*P[i]-P[i+1]-P[i-1]) + mμϕ + mλ/6*ϕ³)*Δt
function ActionDer(a,m,mu,la,F,f₋₁,f₊₁)
    # return m/a*(2*F-(f₋₁+f₊₁)) + a*m*mu*F# + a*m*la/6*F^3
    return m/a^2*(2*F-(f₋₁+f₊₁)) + m*mu*F# + m*la/6*F^3
end
# DONE: Add coupling terms 2f(i)-f(i+1)-f(i-i)
# Understand the discretizing integral and meeting mat. from 16.03


# Use the StochasticDiffEq package to achive the correct result for different solvers
# SDEFunction()
# SDEProblem()
# Use the SimpleDiffEq package to get the SimpleEM for fixed stepsize Euler-Maruyama




# The equation to solve using different schemes:
# dS/dϕ     , S = 1/2 m δτ ∑[(ϕ_i+1 - ϕ_i)^2 / δτ^2 + ω^2ϕ_i^2]
# This means feeding S to a solver
# But here ϕ is an array, so can make a system of equations:
# dS/dϕ_i   , S = 1/2 m δτ [((\phi_i+1 - \phi_i)^2 +(ϕ_i - ϕ_i-1)^2) / δτ^2 + ω^2ϕ_i^2]
#               = 1/2 m δτ [((\phi_i+1 - \phi_i)^2 +(ϕ_i - ϕ_i-1)^2) / δτ^2 + ω^2ϕ_i^2]

using DifferentialEquations
dt = 0.01
ϕ0 = 20.; a=0.5; m=1; mu=1; ϕ₊₁=20.; ϕ₋₁=20.
f(ϕ,t,p) = (m/a*(ϕ^2-ϕ*(ϕ₊₁+ϕ₋₁)) + 0.5*m*mu*a*ϕ)*dt
timespan = (0.0,0.01)
prob = ODEProblem(f,ϕ0,timespan)
sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
# New path after drift from forward Euler:
ϕ0 - sol(0.01)*dt   # 19.800795
plot(sol)
sol(0.01)
sol = solve(prob,ImplicitEuler(),dt=dt,abstol=1e-8,reltol=1e-8)
# New path after drift from forward Euler:
ϕ0 - sol(0.01)*dt   # 19.800795007234047
plot(sol)
sol(0.01)
# 19.800795007234047
F = [20. for i = 1:16]
F[2] -= ActionDer(a,m,mu,0,F[2],F[1],F[3])*dt
F[2]    # 19.8

# function Langevin(N,a,m,mu,la,gaussianD)
#     F = [20. for i = 1:16]
#     Flist = []
#     push!(Flist,F)
#     dt = 0.01
#     prob = ODEProblem(f,ϕ0,timespan)
#     for i=1:N
#         for ii = 1:length(F)
#             sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
#             F[ii] -= sol(0.01)*dt - sqrt(2*dt/a)*rand(gaussianD)
#             # push!(Flist,F)
#         end
#     end
#     return Flist
# end
# Langevin(20,0.5,1,1,0,gaussianD)


# TODO
    # Use the newly found method to run Langevin simulations, and see if systems are more stable for
        # ImplicitEuler
    # And if the system is the same for
        # Forward Euler




# Using DifferentialEquations.jl to solve with different solvers
using DifferentialEquations
using Plots
u0 = 1.5
a=0.5; m=1; mu=1; la=0; f₋₁=1; f₊₁=1
f(t,p,F) = m/a^2*(2*F-(f₋₁+f₊₁)) + m*mu*F
# g(gaussianD,p,P) =
timespan = (0.0,0.01) # Solve from time = 0 to
# time = 1
prob = ODEProblem(f,u0,timespan)
dt = 0.001
sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8) # Solves the ODE
sol(0.001)
plot(sol)

using DifferentialEquations
α=1
u0=1/2
f2(t,p,u) = α*u
tspan = (0.0,1.0)
prob = ODEProblem(f2,u0,tspan)
sol = solve(prob,Euler(),dt=1/2^4) # Solves the ODE
sol[5]      # 0.5234375 not .637
sol.t[8]    # 0.4375    not .438
sol = solve(prob,Vern7())
plot(sol)

ActionDer(a,m,mu,la,u0,1,1)














struct param
    a::Float64
    m::Float64
    mu::Float64
end

function Langevin(N,a,m,mu,la,gaussianD)
    n_tau = 16
    F = [20. for i = 1:n_tau]
    F2 = [20. for i = 1:n_tau]
    F3 = [20. for i = 1:n_tau]
    Flist = Matrix{Float64}(undef,N+1,n_tau)
    F2list = Matrix{Float64}(undef,N+1,n_tau)
    F3list = Matrix{Float64}(undef,N+1,n_tau)
    Flist[1,:] = F
    F2list[1,:] = F2
    F3list[1,:] = F3
    # println(Flist)
    dt = 0.001
    timespan = (0.0,dt)
    randoms1 = rand(gaussianD,N*n_tau)
    for i=1:N
        # println(F)
        for ii = 1:n_tau
            ϕ₋₁ = F[(ii-2+n_tau)%n_tau+1]; ϕ₊₁ = F[(ii)%n_tau+1]; ϕ0 = F[ii]
            f(ϕ,t,p) = (m/a*(2ϕ-(ϕ₊₁+ϕ₋₁)) + m*mu*a*ϕ)*dt
            prob = ODEProblem(f,ϕ0,timespan)
            sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
            sol3 = solve(prob,ImplicitEuler(),dt=dt,abstol=1e-8,reltol=1e-8)
            # println(sol(0.01))
            F[ii] -= sol(dt)*dt - sqrt(2*dt/a)*randoms1[n_tau*(i-1)+ii]
            F2[ii] -= ActionDer(a,m,mu,la,F2[ii],F2[(ii-2+n_tau)%n_tau+1],F2[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*randoms1[n_tau*(i-1)+ii]
            F3[ii] -= sol3(dt)*dt - sqrt(2*dt/a)*randoms1[n_tau*(i-1)+ii]
        end
        Flist[i+1,:] = F
        F2list[i+1,:] = F2
        F3list[i+1,:] = F3
    end
    return Flist, F2list, F3list
end
res1, res2, res3 = Langevin(10,0.5,1,1,0,gaussianD);

res1
res2
res3
begin
    plot(res1[:,1])
    plot!(res2[:,1])
    plot!(res3[:,1])
end



struct AHO_Param
    a::Float64  # Lattice spacing
    m::Float64  # Mass
    mu::Float64 # μ
    λ::Float64
end

function ActionDerSchem(du, u, params, t)
    p = params.p
    xR = @view u[:]
    F_diff_m1 = xR .- xR[vcat(end,1:end-1)]   # dx_j - dx_{j-1}
    F_diff_p1 = xR[vcat(2:end,1)] .- xR       # dx_{j+1} - dx_j

    du .= p.m .* (F_diff_p1 .- F_diff_m1) ./ p.a^2 .- (p.mu .* xR)
end
# ActionDerSchem([1,2,1,1,2,1],params)

function RandScale(du, u, param, t)
    a = param.p.a
    du .= sqrt.(2. ./ a)
end

function LangevinSchem(N,a,m,mu,la,gaussianD)
    n_tau = 16
    F0 = [20. for i = 1:n_tau]
    Flist = Matrix{Float64}(undef,N+1,n_tau)
    Flist[1,:] = F0
    dt = 0.01
    timespan = (0.0,3*N)
    # params = Lvector(p=struct w fields m μ λ)
    params = LVector(p=AHO_Param(a,m,mu,la))
    # Function to calculate change in action for whole path
    sdeprob1 = SDEProblem(ActionDerSchem,RandScale,F0,timespan,params)

    @time sol = solve(sdeprob1, Euler(), progress=true, saveat=0.1/dt, savestart=false,
                dtmax=1e-3, dt=dt, abstol=5e-2,reltol=5e-2)
end

Solution1 = LangevinSchem(80000,0.5,1,1,0,gaussianD)
plot(Solution1)
Solution1.u

begin
    n_burn = 2
    Set1 = Matrix{Float64}(undef,length(Solution1.u[n_burn:end]),length(Solution1.u[1]))
    for i = n_burn:length(Solution1.u)
        Set1[i-n_burn+1,:] = Solution1.u[i]
    end

    # Plot AutoCorrelation
    autocorrdata = AutoCorrR(Set1)
    jkf1 = Jackknife1(autocorrdata)
    # jkf1[:,1]
    display(plot(jkf1[:,1],yerr=jkf1[:,2],title="AutoCorrelation",xlabel="τ",ylabel="Aₒ(τ)"))

    # Plot TPCF
    arr1 = reshape(Set1,:)
    histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:0.1:(floor(maximum(arr1)*10)+1)/10],normed=true,xlabel="x",ylabel="|ψ₀|²",legend=false)
    display(PlotProbDDe(1,1,1,2))

    println("⟨x⟩ = ",Err1(arr1)[1]," with err: ",Err1(arr1)[2])         # - 4.9*10^-4   ± 0.002048
    println("⟨x²⟩ = ",Err1(arr1.^2)[1]," with err: ",Err1(arr1.^2)[2])  # 0.5102        ± 0.002060

    println("⟨x⟩ = ",Jackknife1(arr1)[1]," with err: ",Jackknife1(arr1)[2])         # - 4.9*10^-4   ± 0.002048
    println("⟨x²⟩ = ",Jackknife1(arr1.^2)[1]," with err: ",Jackknife1(arr1.^2)[2])  # 0.5102        ± 0.002060
end





function Langevin(N,a,m,mu,la,gaussianD,filename)
    path = "results/"
    F = [20. for i = 1:n_tau]
    dt = 0.01
    timespan = (0.0,dt)
    n_burn = 3/dt
    n_skip = 3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F)
            # Other solvers
            ϕ₋₁ = F[(ii-2+n_tau)%n_tau+1]; ϕ₊₁ = F[(ii)%n_tau+1]; ϕ0 = F[ii]
            f(ϕ,t,p) = (m/a*(ϕ^2-ϕ*(ϕ₊₁+ϕ₋₁)) + 0.5*m*mu*a*ϕ^2)*dt
            prob = ODEProblem(f,ϕ0,timespan)
            sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
            F[ii] -= sol(dt)*dt - sqrt(2*dt/a)*rand(gaussianD)
            # Own solver
            # F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
        end
    end
    show(F);println()

    # Simulate
    t1 = @timed for i=1:N
        writec123tofile(path,filename,F,i)
        for iii = 1:n_skip+1
            for ii = 1:length(F)
                # Other solvers
                ϕ₋₁ = F[(ii-2+n_tau)%n_tau+1]; ϕ₊₁ = F[(ii)%n_tau+1]; ϕ0 = F[ii]
                f(ϕ,t,p) = (m/a*(ϕ^2-ϕ*(ϕ₊₁+ϕ₋₁)) + 0.5*m*mu*a*ϕ^2)*dt
                prob = ODEProblem(f,ϕ0,timespan)
                sol = solve(prob,Euler(),dt=dt,abstol=1e-8,reltol=1e-8)
                F[ii] -= sol(dt)*dt - sqrt(2*dt/a)*rand(gaussianD)
                # Own solver
                # F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
            end
        end
        if i%2000==0
            println(i)
        end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return
end
println()
Langevin(20,0.5,1,1,0,gaussianD)

begin
    n_tau=16
    β=8
    a=β/n_tau
    Langv1=Langevin(10000,a,1,1,0,gaussianD,"L_dt0.01_Euler_b8.csv")
end

# Langevin expectationvalues
x1 = Err1(GetData("results/CL_4.csv",4,1))
x2 = Err1(GetData("results/CL_4.csv",4,2))
plot(x1[:,1],yerr=x1[:,2])
plot(x2[:,1],yerr=x2[:,2])

# Metropolis expectationvalues
x2 = Err1(GetData("results/measuredObsHO_1_β8_16.csv",4,2))
plot(x2[:,1],yerr=x2[:,2])





#Computed by a Hubbard-Stratonovich transformation.
# Instead of doing the HS transformation described in paper [2], 
#   we just do the mentioned choice of complex μ (or λ).
# Then the imaginary part drifts as:
#       ϕᵢⁱ⁺¹ = ϕᵢⁱ - m μ ϕᵢⁱ dt
# The real and complex parts are decoupled for just the HO, coupling comes from AHO, λ
"""Complex and real part of the derivative of the action wrt. ϕ  
"""
function CActionDer(a,m,μ,λ,Fᵣ,Fᵢ)
    z = m*μ*(Fᵣ+im*Fᵢ) + m*λ/6*(Fᵣ+im*Fᵢ)^3 #+ m/a^2*(2*F-(f₋₁+f₊₁))
    # z = λ/12*(Fᵣ + im*Fᵢ) + im*λ/(12μ+2*im*λ*(Fᵣ+im*Fᵢ))
    return real(z),imag(z)# m/a*(2*Fᵣ-(f₋₁+f₊₁)) + a*m*mu*Fᵣ# + a*m*la/6*Fᵣ^3
end

"""Calculated real and complex part of the weight term in Eq.(110)
exp( - S_B(σ) )"""
function Weight_func(μ,λ,Fᵣ,Fᵢ)
    z = λ/24*(Fᵣ+im*Fᵢ)^2-0.5*log(λ/(12*μ+2*im*λ*(Fᵣ+im*Fᵢ)))
    return real(exp(-z)), imag(exp(-z)) # e⁻ˢ
end


Weight_func(1,0.4,1,0)

# Compute ⟨ϕ²⟩ = ⟨(ϕ_r + Im*ϕ_i)^2⟩ᵩ
function CLangevin(N,a,m,mu,la,gaussianD,filename)
    save_path = "results/"
    intScheme = "fEuler"
    # [x₁ʳ,x₂ʳ,...]
    # [x₁ᶜ,x₂ᶜ,...]
    F_r = 0.1 #[20. for i = 1:n_tau]     # n_tau global
    F_i = 0. #[1. for i = 1:n_tau]
    # println("The complex weight starts at: ",Weight_func(mu,la,F_r,F_i))

    dt = 0.01
    n_burn = 0#3/dt
    n_skip = 0#3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F_r)
            if intScheme == "fEuler"
                derAction = CActionDer(a, m, mu, la, F_r, F_i)
                # Forward Euler: F_{n+1} = F_{n} + f(t_n,F_{n})
                F_r -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD) # N_R = 1 => N_I = 0
                F_i -= derAction[2]*dt
                # F_r[ii] -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD)
                # F_i[ii] -= derAction[2]*dt
            else
                # Backward Euler: F_{n+1} = F_{n} + f(t_{n+1},F_{n+1})
                F_r0,F_i0 = copy(F_r),copy(F_i)
                F_r, F_i = (F_r0,F_i0).+CActionDer(a, m, mu, la, F_r, F_i).*dt
                # F_r, F_i = (1,2).+(1,2).*0.1 = (1.1, 2.2)
                
                # Instead of two variables, can go to one complex
                # Matricies / vectors?
                # Then f(F,t,F_0) = F_0 + CActionDer(a,m,mu,la,F)*dt = F_0 + m*μ*F + m*λ/6*F^3
                    # + sqrt(2*dt)*rand(gaussianD)

                # Then using the DifferentialEquations package:
                # prob = ODEProblem(f,F_0,)
                # sol = solve(prob,ImplicitEuler())
                F_r += sqrt(2*dt)*rand(gaussianD)
            end
        end
    end
    show(F_r);println()

    Flist_r = []
    Flist_i = []
    # push!(Flist_r,F_r)
    # push!(Flist_i,F_i)
    # WeightP_r = []
    # WeightP_i = []
    # WeightN_r = []
    # WeightN_i = []
    weight_c = 0# Weight_func(mu,la,F_r,F_i)
    # push!(Weight_r,weight_c[1])
    # push!(Weight_i,weight_c[2])
    
    # Simulate
    t1 = @timed for i=1:N
        # writec123tofile(path,filename,F_r,i)
        push!(Flist_r,F_r)
        push!(Flist_i,F_i)
        writec123tofile(save_path,filename,F_r,i)
        # weight_c = Weight_func(mu,la,F_r,F_i)
        # if F_r > 0
        #     push!(WeightP_r,weight_c[1])
        #     push!(WeightP_i,weight_c[2])
        # else
        #     push!(WeightN_r,weight_c[1])
        #     push!(WeightN_i,weight_c[2])
        # end
        for iii = 1:n_skip+1
            for ii = 1:length(F_r)
                # derAction = CActionDer(a,m,mu,la,F_r[ii],F_r[(ii-2+n_tau)%n_tau+1],F_r[(ii)%n_tau+1,F_c[ii],F_c[(ii-2+n_tau)%n_tau+1],F_c[(ii)%n_tau+1]])
                # Possible for hashtable ii -> index in F_r / F_c
                derAction = CActionDer(a, m, mu, la, F_r, F_i)
                F_r -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD) # N_R = 1 => N_I = 0
                F_i -= derAction[2]*dt
                # F_r[ii] -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD)
                # F_i[ii] -= derAction[2]*dt
            end 
        end
        if i%2000==0
            println(i)
        end
    end
    println("t: ",t1.time, " t2:", t1.gctime)
    return Flist_r, Flist_i#WeightP_r, WeightP_i, WeightN_r, WeightN_i, 
end

# μ = e^iϕ, ϕ = (0,2π) (+n*2π)
ComplexSys = CLangevin(20000,0.5,1,0.05*im+1,0,gaussianD,"CL_2")
# ComplexSys = CLangevin(20000,0.5,1,1,0,gaussianD,"CL_2")
for i = 0:11
    if i==0
        # Calculate the analytical result 1/μ = ⟨z²⟩, where μ = exp(nπi/6), n = (0,11)
        # ⟹ 1/μ = exp(-nπi/6), n = (0,11)
        scatter([cos(ii*π/6) for ii=0:11],[sin(ii*π/6) for ii=0:11],color="red",marker=:x,legend=false)#:inside)
        # scatter([real(exp(-im*ii*π/6)) for ii=0:11],[imag(exp(-im*ii*π/6)) for ii=0:11],color="red",legend=:inside,marker=:x)
    end
    arr2=[]
    for runs = 1:64
        ComplexSys = CLangevin(2000,0.5,1,exp(i*im*π/6),0,gaussianD,"CL_2")
        append!(arr2,getExp2(ComplexSys[1],ComplexSys[2])[1])     # ⟨x²⟩
    end
    # display(scatter(ComplexSys[1],ComplexSys[2]))
    # arr1 = float.(ComplexSys[1])
    println("i: ",i,"e^z:",exp(i*im))
    # display(histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:incsize1:(floor(maximum(arr1)*10)+1)/10],normed=true,xlabel="x",ylabel="|ψ_0|²"))
    if in(i,[0,1,2,3,9,10,11])
        arr3 = [mean(arr2),Err1(real.(arr2))[2],Err1(imag.(arr2))[2]]
        fig1 = scatter!([real(arr3[1])],[imag(arr3[1])],xerr=arr3[2],yerr=arr3[3],color="blue",marker=:cross)
        # fig1 = scatter!([real(arr2[1])],[imag(arr2[1])],xerr=arr2[2],yerr=arr2[3],color="blue",marker=:cross)
        display(fig1)
        if true
            if i == 11
                savefig(fig1,"plots/22.04.22_CL_gauss_mod2.pdf") # This is how to save a Julia plot as pdf !!!
            end
        end
    end
end
# Diverges to NaN coordinates late in τ time when Re(z) → 0, Im(z) → 1
# Scatterplot showes values collected moves from only real part to uniform real/complex parts
# Find out for which values eᶻ should diverge


# scatter(ComplexSys[3],ComplexSys[4],yrange=[-0.004,0.003],xlabel="Re[ρ]",ylabel="Im[ρ]")
ComplexSys[1]
scatter(ComplexSys[1],ComplexSys[2])

# NIntegrate[ x^2*Exp[-(1/2)*m*\[Mu]*x^2 - m*(\[Lambda]/24)*x^4], {x, -\[Infinity], \[Infinity]}]
function getExp2(field_r,field_c)
    z = []
    for i = 1:length(field_r)
        append!(z,(field_r[i]+im*field_c[i])^2)
    end
    return append!([mean(z)], Err1(real.(z))[2], Err1(imag.(z))[2])
end


# CL expectationvalues = 1 ???
Err1(ComplexSys[1])                         # ⟨x_r⟩
Err1(ComplexSys[2])                         # ⟨x_i⟩
getExp2(ComplexSys[1],ComplexSys[2])[1]     # ⟨x²⟩
1/(1+0.05*im) #1/μ
arr1 = float.(ComplexSys[1])
incsize1= 0.1
histogram(arr1,bins=[i for i=floor(minimum(arr1)*10)/10:incsize1:(floor(maximum(arr1)*10)+1)/10],normed=true,xlabel="x",ylabel="|ψ_0|²")




# mean(ComplexSys[5])
# mean(ComplexSys[6])


# for i = 1:length(ComplexSys)

# end
# mean(6/(6*(1+1im)+im*(0.4+1im)*ComplexSys))

TwopointE = []
n_tau = 16
for i=0:n_tau-1
    m=1
    mu=1
    dt = 0.5
    R = 1+mu^2*dt^2/2-mu*dt*sqrt(1+mu^2*dt^2/4)
    res= (R^i+R^(n_tau-i))/(1-R^n_tau)
    res/=2*m*mu
    push!(TwopointE,res)
end


Exp_x2e(16, 0.5, 1, 1)
plot!(TwopointE)

    # Probability density diagram #
PlotProbDD("results/CL_4.csv",0.1)
PlotProbDDe(1,1,1,3)
    # sampling
scatter(reshape(GetColumn(2,"results/CL_4.csv"),:))#:Int((length(LastRowFromFile(file))-1)/4)+1
    # Autocorrelation
PlotAC("results/CL_1.csv",1000)
# PlotACsb("results/CL_1.csv",1000)
# PlotAC("results/CL_1.csv",false)
# PlotAC("results/CL_1.csv",true)
    # Twopoint Correlation
PlotTPCF("results/CL_3.csv")                # Naive error
PlotTPCF("results/CL_1.csv",true)           # For autocorrelated data
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
 PlotEffM("results/CL_1.csv")