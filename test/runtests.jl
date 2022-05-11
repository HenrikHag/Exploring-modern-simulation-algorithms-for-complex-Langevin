using MetroLangCLsimulation
using Random
using Test








# difActionAHO vs AHO_Action
testrng = MersenneTwister(11111)
n_tau=16; a=0.5; m=1; ω=1; λ=0;
Path = [2*(rand(testrng)-1/2) for i=1:n_tau]
i = 3
ip1 = ((i)%n_tau)+1
im1 = ((i-2+n_tau)%n_tau)+1
x_new = 1*(Path[i]+rand(testrng))

println("Old: ",AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i]))
println("New: ",difActionAHO(n_tau,a,m,ω,λ,Path,i,x_new))

@testset "difActionAHO" begin
    @test difActionAHO(n_tau,a,m,ω,λ,Path,i,Path[i]) ≈ 0. atol=1e-10
    # @test 1*(0.5*m*a*(((Path[ip1]-x_new)/a)^2 + ((x_new-Path[im1])/a)^2 + ω^2*x_new^2)-0.5*m*a*(((Path[ip1]-Path[i])/a)^2 + ((Path[i]-Path[im1])/a)^2 + ω^2*Path[i]^2)) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
    # @test 1*(0.5*m*a*((x_new^2-2*x_new*Path[ip1])/a^2 + ((x_new-Path[im1])/a)^2 + ω^2*x_new^2)-0.5*m*a*((Path[i]^2-2*Path[i]*Path[ip1])/a^2 + ((Path[i]-Path[im1])/a)^2 + ω^2*Path[i]^2)) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
    # @test 1*0.5*m*a*(((x_new^2-2*x_new*Path[ip1])/a^2 + (x_new^2-2*x_new*Path[im1])/a^2 + ω^2*(x_new^2- Path[i]^2)) - ((Path[i]^2-2*Path[i]*Path[ip1])/a^2 + (Path[i]^2-2*Path[i]*Path[im1])/a^2)) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
    # @test m*a*(((x_new^2-Path[i]^2+(Path[i]-x_new)*(Path[ip1]+Path[im1]))/a^2 + 0.5*ω^2*(x_new^2 - Path[i]^2))) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
    # @test m*(x_new^2-Path[i]^2+(Path[i]-x_new)*(Path[ip1]+Path[im1]))/a + 0.5*m*a*ω^2*(x_new^2 - Path[i]^2) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
    @test difActionAHO(n_tau,a,m,ω,λ,Path,i,x_new) ≈ (AHO_Action(n_tau,m,ω,a,λ,Path,i,x_new)-AHO_Action(n_tau,m,ω,a,λ,Path,i,Path[i])) atol=1e-8
end
