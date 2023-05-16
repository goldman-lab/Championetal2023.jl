function seungautapse_sim(params)

    _netparams = SeungNetParams(
        params["g_Na"],#::Float64
        params["g_K"],#::Float64
        params["g_L"],#::Float64
        params["E_Na"],#::Float64
        params["phi"],#::Float64     
        params["E_K"],#::Float64     
        params["C_m"],#::Float64     
        params["E_L"],#::Float64     
        params["E_E"],#::Float64     
        params["E_I"],#::Float64     
        params["spike_threshold"],  #::Float64
        params["theta_s"],  #::Float64
        params["sigma_s"],  #::Float64
        params["alpha_s"],  #::Float64
        params["tau_S"],  #::Float64
        params["tau_O"],  #::Float64
        params["tau_E"],  #::Float64
        params["tau_I"],  #::Float64
        params["tau_A"],  #::Float64
        params["g_A"],  #::Float64
        params["w_S"],  #::Float64
        params["w_O"],  #::Float64
        params["w_E"],  #::Float64
        params["w_I"],  #::Float64
        params["I_app_O"],  #::Float64
        params["I_app_E"],  #::Float64
        params["I_app_I"]  #::Float64
        )

    _runparams = RunParams(
        params["I_O"],#::Float64      
        params["gamma_I_t"],#::Float64      
        params["omega_I_t"],#::Float64      
        params["dt"],#::Float64        #time step in ms
        params["tmax"],#::Float64     #time in ms to run simulation
        params["downsamplespacing"]#::Int64 #downsample for plotting
        )

    _simvals = SimVals(
        0:_runparams.dt:_runparams.tmax, #t_vect
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #plot_t_vect::Array{Float64,1}
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))) #input_t_vect::Array{Float64,1}
        )

    _noiseparams = NoiseParams(
        params["noisesigma"],#::Float64
        params["tau_n"],#::Float64     noise time constant
        params["gammasigma"],#::Float64
        params["tau_n_gamma"],#::Float64     noise time constant
        params["omegasigma"],#::Float64
        params["tau_n_omega"],#::Float64     noise time constant
        )

    _inputparams = InputParams(
        vcat(params["pulsestarts"], 1e23), #pulse starts plus a realmax pad to make things easy
        vcat(params["pulsestops"], 1e23), #pulse stops plus a realmax pad to make things easy
        vcat(params["pulseheights"], 1e23),
        1, #pulsecounter::Int64
        0, #pulsewindow::Int64
        0.0 #input::Float64
        )


    _snode = SeungNode(
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), 
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), 
        spzeros(Float64, length(0:_runparams.dt:_runparams.tmax)), 
        0.1, # fn
        0.1, # n
        0.9, #  fh
        0.9, # h
        0.0, # fs
        0.0, # s
        -55.0, # fv
        -55.0, # v
        zeros(3), # vstore
        0, # spikeyesno
        0.0, # noise
        0.02, # fb
        0.02 # b
        );

    if params["simtype"] == :stdystate
        # Call primary implementation
        return _seungautapse_stdystatesim(_netparams, _runparams, _simvals, _snode, _noiseparams, _inputparams, params["synapseactivity"])
    else
        # Call primary implementation
        return _seungautapse_sim(_netparams, _runparams, _simvals, _snode,  _noiseparams, _inputparams)
    end

end
