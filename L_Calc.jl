# Plot results of Langevin simulation

# julia> ] activate .
begin
    using MCMC
    # using .UsersGuide
    # using .MetropolisUpdate
    using Plots
    using BenchmarkTools
    using StatsBase
    using DataFrames
    save_date = findDate()
    save_folder = "plots/"
    # using GLM
end

# Filename of datafile to be analysed
measf = "results/CL_1.csv"#"results/measuredObsB100S0_7.csv"
measf = "results/CL_4.csv"
measf = "results/L_dt0.01_Euler_b8.csv"
measf = "results/L_dt0.1_b8.csv"
measf = "results/22.05.21_L_dt0.001_b8.csv"





begin   # ⟨x₁⟩
    # Plot of x₁ with running mean
    save_name = "$(save_folder)$(save_date)_L_expect_x_1"
    a1=[]
    a = GetData(measf,4,1)[1:400,1]
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)
    # savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    # savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end

begin   # ⟨x₁²⟩
    # Plot of x₁² with running mean
    save_name = "$(save_folder)$(save_date)_L_expect_x2_1"
    a1=[]
    a = GetData(measf,4,1)[1:400,1].^2
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    plt = scatter(a)
    plot!(plt,a1,width=4)
    display(plt)
    # savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    # savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
    Jackknife1(a,true),Exp_x2e(16,0.5,1,1) # Estimated ⟨x₁²⟩ with error, and analytical
end

begin   # ⟨xᵢ⟩, ⟨xᵢ²⟩
    # Plot of xᵢ's and xᵢ²'s with expectation values
    save_name = "$(save_folder)$(save_date)_L_dt0.001_b8_x_x2"
    
    a1 = Jackknife1(GetData(measf,4,1)[:,:],true)
    a2 = Jackknife1(GetData(measf,4,1)[:,:].^2,true)
    
    plt = plot(a1[:,1],yerr=a1[:,2],label="⟨xᵢ⟩",legend=:right,xlabel="O(τ)",ylabel="⟨O⟩")
    plot!(plt,a2[:,1],yerr=a2[:,2],label="⟨xᵢ²⟩")
    hline!(plt,[[Exp_x2e(16,0.5,1,1)],[Exp_x2(16,0.5,1,1)],[0]],label=["⟨xᵢ²⟩ₜₕₑₒ continuum" "⟨xᵢ²⟩ₜₕₑₒ dicretized" ""],color=["black" "green" "black"]) # Analytical expectation values for system
    # hline!(plt,[],color="green")

    display(plt)
    savefig(plt,"$(save_name).pdf")   # Save as pdf in folder "plots"
    savefig(plt,"$(save_name).png")   # Save as png in folder "plots"
end