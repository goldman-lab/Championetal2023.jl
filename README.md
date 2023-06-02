# Championetal2023

## Code to accompany the paper 'An Oscillatory Mechanism for Multi-level Storage in Short-term Memory'

## The best way to use this is to install it as a package directly from Julia.
### from the Julia REPL:

   
    ] add https://github.com/goldman-lab/Championetal2023.jl.git

 
### If you plan to make extensions or develop from it, then:
 
    ] dev https://github.com/goldman-lab/Championetal2023.jl.git




### The best way to get going here is to use the notebooks included in the package. I've tied up all of the code you'll need to run the simulations in each paper figure into individual notebooks. The first time you use the notebooks, you may be prompted to add some packages that you'll need.

### Bear in mind that most figures in the paper were plotted using DataGraph, so expect minor differences in interpolation, smoothing, etc. 
### The code has been re-written here for easier readability and extension, at the expense of some speed. I've plotted the simulation durations for most of the runs that will take some time.

## General Flow of the simulations
### Most of the figures (Figs. 2,3,5,6,7,8) use calls to the simulation backend for most of the computations. The main simulation functions ingest a dictionary of parameter values that are wrapped and sent to the primary simulation functions. Results are typically returned as a data structure. 
