struct Wts
    w_EE::Float64      #excitatory auto-recurrence
end

struct NetParams
    g_Na::Float64
    g_K::Float64
    g_L::Float64
    E_Na::Float64
    phi::Float64
    E_K::Float64
    C_m::Float64
    E_L::Float64
    tau_syn::Float64
    g_syn::Float64
    spike_threshold::Float64
end

struct SeungNetParams
    g_Na::Float64
    g_K::Float64
    g_L::Float64
    E_Na::Float64
    phi::Float64
    E_K::Float64
    C_m::Float64
    E_L::Float64
    E_E::Float64
    E_I::Float64
    spike_threshold::Float64
    theta_s::Float64
    sigma_s::Float64
    alpha_s::Float64
    tau_S::Float64
    tau_O::Float64
    tau_E::Float64
    tau_I::Float64
    tau_A::Float64
    g_A::Float64
    w_S::Float64
    w_O::Float64
    w_E::Float64
    w_I::Float64
    I_app_O::Float64
    I_app_E::Float64
    I_app_I::Float64
end

struct RunParams
    I_O::Float64
    gamma_I_t::Float64
    omega_I_t
    dt::Float64
    tmax::Float64
    downsamplespacing::Int64
end

struct RunParamsStdyState
    I_O::Float64
    I_ext::Float64
    gamma_I_t::Float64
    omega_I_t #::Type inferred at runtime
    dt::Float64
    tmax::Float64
    downsamplespacing::Int64
end

struct SimVals #pushing most of these into local storage in nodes
    t_vect::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    plot_t_vect::Array{Float64,1}
    input_t_vect::Array{Float64,1}
end

mutable struct SingleIntegratorNode
    EE::Float64 #placeholder for code-sharing, will become vector of connectivity weights in network
    v_vect::Array{Float64,1}
    s_vect::Array{Float64,1}
    spikes_vect::SparseVector
    fn::Float64
    n::Float64
    fh::Float64
    h::Float64
    fs::Float64
    s::Float64
    fv::Float64
    v::Float64
    vstore::Array{Float64,1}
    spikeyesno::Int64
    noise::Float64
end

mutable struct SeungNode
    v_vect::Array{Float64,1}
    s_vect::Array{Float64,1}
    spikes_vect::SparseVector
    fn::Float64
    n::Float64
    fh::Float64
    h::Float64
    fs::Float64
    s::Float64
    fv::Float64
    v::Float64
    vstore::Array{Float64,1}
    spikeyesno::Int64
    noise::Float64
    fb::Float64
    b::Float64
end

mutable struct SingleIntegratorNet
    EE::Array{Float64,2} #placeholder for code-sharing, will become vector of connectivity weights in network
    s_mat::Array{Float64,2}
    spikes_mat::SparseMatrixCSC
    fn::Array{Float64,1}
    n::Array{Float64,1}
    fh::Array{Float64,1}
    h::Array{Float64,1}
    fs::Array{Float64,1}
    s::Array{Float64,1}
    fv::Array{Float64,1}
    v::Array{Float64,1}
    vstore::Array{Float64,2}
    spikeyesno::Array{Int64,1}
    noise::Array{Float64,1}
    phase::Array{Float64,1}
end

mutable struct SequenceIntegratorNode
    EE::Array{Float64,1} #Float64 #placeholder for code-sharing, will become vector of connectivity weights in network
    s_vect::Array{Float64,1}
    spikes_vect::SparseVector
    fn::Float64
    n::Float64
    fh::Float64
    h::Float64
    fs::Float64
    s::Float64
    fv::Float64
    v::Float64
    vstore::Array{Float64,1}
    spikeyesno::Int64
    noise::Float64
    phase::Float64
end

struct SingleNeuronResult
    timevect::Array{Float64,1}
    v_vect::Array{Float64,1}
    spikes::Array{Float64,1}
end

struct NoiseParams
    sigma::Float64 # input noise level
    tau_n::Float64     #input noise time constant
    gammasigma::Float64 # input noise level
    tau_n_gamma::Float64 # input noise level
    omegasigma::Float64     #input noise time constant
    tau_n_omega::Float64     #input noise time constant
end

abstract type InputContainer end

mutable struct InputParams <: InputContainer
    pulsestarts::Array{Float64,1}
    pulsestops::Array{Float64,1}
    pulseheights::Array{Float64,1}
    pulsecounter::Int64
    pulsewindow::Int64
    input::Float64
end

mutable struct SequenceInputParams <: InputContainer
    pulsestarts::Array{Float64,1}
    pulsestops::Array{Float64,1}
    pulseheights::Array{Float64}
    pulsecounter::Int64
    pulsewindow::Int64
    input::Array{Float64,1}
    inputnull::Array{Float64,1}
end
