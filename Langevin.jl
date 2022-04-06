begin
    using Distributions
    using Plots
    using .UsersGuide
    using .MetropolisUpdate
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

function Langevin(N,a,m,mu,la,gaussianD)
    F = [20. for i = 1:16]
    Flist = []
    push!(Flist,F)
    dt = 0.01
    for i=1:N
        for ii = 1:length(F)
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
            # push!(Flist,F)
        end
    end
    return Flist
end
function Langevin(N,a,m,mu,la,gaussianD,filename)
    path = "results/"
    F = [20. for i = 1:n_tau]
    dt = 0.01
    n_burn = 3/dt
    n_skip = 3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F)
            F[ii] -= ActionDer(a,m,mu,la,F[ii],F[(ii-2+n_tau)%n_tau+1],F[(ii)%n_tau+1])*dt - sqrt(2*dt/a)*rand(gaussianD)
        end
    end
    show(F);println()

    # Simulate
    t1 = @timed for i=1:N
        writec123tofile(path,filename,F,i)
        for iii = 1:n_skip+1
            for ii = 1:length(F)
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
println()
Langevin(20,0.5,1,1,0,gaussianD)

begin
    n_tau=16
    β=8
    a=β/n_tau
    Langv1=Langevin(1500,a,1,1,0,gaussianD,"CL_3.csv")
end


"""Complex and real part of the derivative of the action wrt. ϕ  
Computed by a Hubbard-Stratonovich transformation."""
function CActionDer(μ,λ,Fᵣ,Fᵢ)
    z = λ/12*(Fᵣ + im*Fᵢ) + im*λ/(12μ+2*im*λ*(Fᵣ+im*Fᵢ))
    return real(z),imag(z)# m/a*(2*Fᵣ-(f₋₁+f₊₁)) + a*m*mu*Fᵣ# + a*m*la/6*Fᵣ^3
end

"""Calculated real and complex part of the weight term in Eq.(110)
exp( - S_B(σ) )"""
function Weight_func(μ,λ,Fᵣ,Fᵢ)
    z = λ/24*(Fᵣ+im*Fᵢ)^2-0.5*log(λ/(12*μ+2*im*λ*(Fᵣ+im*Fᵢ)))
    return real(exp(-z)), imag(exp(-z)) # e⁻ˢ
end


Weight_func(1,0.4,1,0)

function CLangevin(N,a,m,mu,la,gaussianD,filename)
    path = "results/"
    # [x₁ʳ,x₂ʳ,...]
    # [x₁ᶜ,x₂ᶜ,...]
    F_r = 0.1 #[20. for i = 1:n_tau]     # n_tau global
    F_i = 0. #[1. for i = 1:n_tau]
    println("The complex weight starts at: ",Weight_func(mu,la,F_r,F_i))

    dt = 0.001
    n_burn = 3/dt
    n_skip = 3/dt

    # Thermalize
    for i=1:n_burn
        for ii = 1:length(F_r)
            derAction = CActionDer(mu,la,F_r,F_i)
            F_r -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD) # N_R = 1 => N_I = 0
            F_i -= derAction[2]*dt
            # F_r[ii] -= derAction[1]*dt - sqrt(2*dt)*rand(gaussianD)
            # F_i[ii] -= derAction[2]*dt
        end
    end
    show(F_r);println()

    Flist_r = []
    Flist_i = []
    push!(Flist_r,F_r)
    push!(Flist_i,F_i)
    WeightP_r = []
    WeightP_i = []
    WeightN_r = []
    WeightN_i = []
    weight_c = 0# Weight_func(mu,la,F_r,F_i)
    # push!(Weight_r,weight_c[1])
    # push!(Weight_i,weight_c[2])
    
    # Simulate
    t1 = @timed for i=1:N
        # writec123tofile(path,filename,F_r,i)
        push!(Flist_r,F_r)
        push!(Flist_i,F_i)
        weight_c = Weight_func(mu,la,F_r,F_i)
        if F_r > 0
            push!(WeightP_r,weight_c[1])
            push!(WeightP_i,weight_c[2])
        else
            push!(WeightN_r,weight_c[1])
            push!(WeightN_i,weight_c[2])
        end
        for iii = 1:n_skip+1
            for ii = 1:length(F_r)
                # derAction = CActionDer(a,m,mu,la,F_r[ii],F_r[(ii-2+n_tau)%n_tau+1],F_r[(ii)%n_tau+1,F_c[ii],F_c[(ii-2+n_tau)%n_tau+1],F_c[(ii)%n_tau+1]])
                # Possible for hashtable ii -> index in F_r / F_c
                derAction = CActionDer(mu,la,F_r,F_i)
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
    return WeightP_r, WeightP_i, WeightN_r, WeightN_i#Flist_r, Flist_i
end

ComplexSys = CLangevin(2000,0.5,1,1,0.4,gaussianD,"CL_2")

scatter(ComplexSys[3],ComplexSys[4],yrange=[-0.004,0.003],xlabel="Re[ρ]",ylabel="Im[ρ]")
scatter!(ComplexSys[1],ComplexSys[2])

for i = 1:length(ComplexSys)

end
mean(6/(6*(1+1im)+im*(0.4+1im)*ComplexSys))

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
plot!(TwopointE)

    # Probability density diagram #
PlotProbDD("results/CL_3.csv",0.1)
PlotProbDDe(1,1,1,3)
    # sampling
plot(reshape(GetColumn(2,"results/CL_1.csv"),:))#:Int((length(LastRowFromFile(file))-1)/4)+1
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