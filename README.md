# Exploring modern simulation algorithms for complex Langevin

This is the Julia code used in the developement of my Masters thesis "Exploring modern simulation algorithms for complex Langevin".

---

# To run

Download to your local enviroment.  
Install a stable version of Julia.  
Open a Julia REPL in the folder `myProject` ( `> julia` in ubuntu )  
Run `] activate .`  
Run `] instantiate`

All dependencies will be downloaded and the two scripts `MonteCarloOld.jl` and `Langevin.jl` can be used to run simulations.  
The scripts `MC_Calc.jl` and `L_Calc.jl` can be run to do analysis on simulations run in the previous.

By using VScode, select or open the working enviroment to `myProject`  
Do the same steps in the VSCode terminal.  
All plots will be shown in its own plotting window.

---

# Figures
Figures to visualize the theoretical and simulated expectation values obtained from an example simulation.  

## Expectation values

<figure>
  <image src="saved_plots/22.06.12_M_longSim_10k_40skip_jk_x_x2.png" alt="" width="450" />
  <figcaption>Expectation values for the different Euclidean time points</figcaption>
</figure>

## Autocorrelation

<figure>
  <image src="saved_plots/22.06.10_M_10k_40skip_20AC.png" alt="" width="450" />
  <figcaption>Autocorrelation as a function of simulation time inteval</figcaption>
</figure>

## Two-point correlation

<figure>
  <image src="saved_plots/22.06.10_M_10k_40skip_TPCF.png" alt="" width="450" />
  <figcaption>Two-point correlation as a function of euclidean time interval</figcaption>
</figure>
