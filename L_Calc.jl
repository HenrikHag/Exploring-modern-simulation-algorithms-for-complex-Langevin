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
measf = "results/CL_4.csv"#"results/measuredObsB100S0_7.csv"



begin   # ⟨x₁⟩
    a1=[]
    a = GetData(measf,4,1)[1:400,1]
    for i = 1:length(a)
        append!(a1,mean(a[1:i]))
    end
    scatter(a)
    plt = plot!(a1,width=4)
    display(plt)                                   # Save as png manually
    savefig(plt,"plots/22.05.01_L_expect_x1.pdf")   # Save as pdf in folder "plots"
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
    savefig(plt,"plots/22.05.01_L_expect_x2.pdf")   # Save as pdf in folder "plots"
end

begin   # ⟨xᵢ⟩, ⟨xᵢ²⟩
    a1 = Jackknife1(GetData(measf,4,1)[1:400,:])
    plot(a1[:,1],yerr=a1[:,2],legend=false)
    a2 = Jackknife1(GetData(measf,4,1)[1:400,:].^2)
    plt = plot!(a2[:,1],yerr=a2[:,2])
    display(plt)                                   # Save as png manually
    savefig(plt,"plots/22.05.01_L_expect_x_i.pdf")   # Save as pdf in folder "plots"
end