# Championetal2023

## Code to accompany the paper 'An Oscillatory Mechanism for Multi-level Storage in Short-term Memory'

## To use the code, install it as a package directly from Julia.
### from the Julia REPL:

   
    ] add https://github.com/goldman-lab/Championetal2023.jl.git

 
### If you plan to make extensions or develop from it, then:
 
    ] dev https://github.com/goldman-lab/Championetal2023.jl.git


### All code in this respository was written and tested in Julia 1.7.1 on Mac OS; head to https://julialang.org/downloads/oldreleases/ for old releases. 

### The best way to get going here is to use the notebooks included in the package. I've tied up all of the code you'll need to run the simulations in each paper figure into individual notebooks. The first time you use the notebooks, you'll likely be prompted to add some packages that you'll need.

### Bear in mind that most figures in the paper were plotted using DataGraph, so expect minor differences in interpolation, smoothing, etc. 
### I've plotted the simulation durations for most of the runs that will take some time.

## General Flow of the simulations
### Most of the figures (Figs. 2,3,5,6,7,8) use calls to the simulation backend for the majority of the computations. The user-facing simulation functions ingest a dictionary of parameter values that are wrapped and then sent to the primary simulation functions. Results are typically returned as a data structure. 
