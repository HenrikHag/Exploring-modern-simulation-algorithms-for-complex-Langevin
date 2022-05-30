# Plot different values
using MCMC
using Plots
using Random


# Stiff vs. order Solvers diagram
LambaEM = (0.5,0.2)
DRI1 = (1.5,0)
ImplicitEM = (0.5,0.9)
SKenCarp = (1.5,1)

Solverlist = [LambaEM, DRI1, ImplicitEM, SKenCarp]
Namelist = ["LambaEM","DRI1","ImplicitEM","SKenCarp"]

begin
    # fig_Solvers = hline([0],label=false,color="black")
    # fig_Solvers = vline!([0],label=false,color="black")
    fig_Solvers = scatter([Solverlist[1][2]],[Solverlist[1][1]],label=Namelist[1],xlim=[0,1],ylim=[0,2],xlabel="Stiffness",ylabel="Order")
    for i=2:length(Solverlist)
        fig_Solvers = scatter!([Solverlist[i][2]],[Solverlist[i][1]],label=Namelist[i])
    end
    display(fig_Solvers)
    savefig(fig_Solvers,"plots/$(findDate())_Solvers_orderStiff.pdf")
    savefig(fig_Solvers,"plots/$(findDate())_Solvers_orderStiff.png")
end

# Euclidean time path
rng = MersenneTwister(11111)

"""Does a metroswipe by testing a change for each element,  
and adds this new coord weighted by change in action."""
function MetroSwipe(n_tau::Int64, m, ω, λ, a, h, idrate, rng, Path)
    accept = 0
    for i = 1:n_tau
        x_new = Path[i] + h*2*(rand(rng) - 0.5)
        if rand(rng) < exp(-difActionAHO(n_tau,a,m,ω,λ,Path,i,x_new))
            Path[i] = x_new
            accept += 1/n_tau
        end
    end
    # rand!(gaussianD, Path)
    # randn!(rng, Path)
    h *= accept / idrate
    # println(h)
    return Path, accept, h
end
N_tau = 5
a = 1
m,ω,λ,h,idrate=1,1,0,1,0.8
path1 = [0. for i=1:N_tau]
n_burn = 20
begin
    for i=1:n_burn
        path1, accept, h = MetroSwipe(N_tau, m, ω, λ, a, h, idrate, rng, path1)
    end
    vline([0],label=false,xlabel)
    scatter!([0 for i=1:N_tau],[i*a for i=1:N_tau],label=false,color="black")
    plot!(path1,[i*a for i=1:N_tau],label="X₀")
    path1, accept, h = MetroSwipe(N_tau, m, ω, λ, a, h, idrate, rng, path1)
    plot!(path1,[i*a for i=1:N_tau],label="X₁")
    for i=1:19
        path1, accept, h = MetroSwipe(N_tau, m, ω, λ, a, h, idrate, rng, path1)
        println(i)
    end
    plot!(path1,[i*a for i=1:N_tau],label="X₂₀")
end
xlabel!("Real time t")
ylabel!("Euclidean time τ")
savefig("plots/$(findDate())_EuclideanTimeContourPath.pdf")
savefig("plots/$(findDate())_EuclideanTimeContourPath.png")