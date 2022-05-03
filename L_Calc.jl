begin
    using .UsersGuide
    using .MetropolisUpdate
    using Plots
    using BenchmarkTools
    using StatsBase
    using DataFrames
    using GLM
end
# 
measf = "results/CL_1.csv"#"results/measuredObsB100S0_7.csv"
measf = "results/CL_4.csv"
measf = "results/L_dt0.01_Euler_b8.csv"
measf = "results/L_dt0.1_b8.csv"




begin   # ⟨x₁⟩
    a1=[]
    a = GetData(measf,4,1)[1:400,1]
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)                                   # Save as png manually
    # savefig(plt,"plots/22.05.03_L_expect_x1.pdf")   # Save as pdf in folder "plots"
end

begin   # ⟨x₁²⟩
    a1=[]
    a = GetData(measf,4,1)[1:400,1].^2
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)                                   # Save as png manually
    # savefig(plt,"plots/22.05.03_L_expect_x2.pdf")   # Save as pdf in folder "plots"
end
Exp_x2e(16,0.5,1,1)

begin   # ⟨xᵢ⟩, ⟨xᵢ²⟩
    a1 = Jackknife1(GetData(measf,4,1)[:,:])
    plot(a1[:,1],yerr=a1[:,2],legend=false)
    a2 = Jackknife1(GetData(measf,4,1)[:,:].^2)
    plot!(a2[:,1],yerr=a2[:,2])
    plt = hline!([Exp_x2e(16,0.5,1,1),0])
    display(plt)                                   # Save as png manually
    # savefig(plt,"plots/22.05.03_L_expect_x_i.pdf")   # Save as pdf in folder "plots"
end