module Championetal2023

using SparseArrays, LinearAlgebra, Distributions


include("types.jl")
include("functions.jl")
include("construction.jl")
include("runsingleneuron.jl")
include("runsingleneuron_stdystate.jl")
include("runalltoallnetwork.jl")
include("runalltoallnetwork_stdystate.jl")
include("runseungautapse.jl")
include("seungmodelconstruction.jl")
include("runwbsynapsephaseshift.jl")


export oscintegrator_sim
export oscintegrator_stdystate
export oscint_ata_stdystate
export getedges
export seungautapse_sim
export getsequenceweightmatrix


end
