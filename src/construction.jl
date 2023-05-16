function oscintegrator_sim(params)

    _netparams = NetParams(
        params["g_Na"],#::Float64
        params["g_K"],#::Float64
        params["g_L"],#::Float64
        params["E_Na"],#::Float64
        params["phi"],#::Float64     
        params["E_K"],#::Float64     
        params["C_m"],#::Float64     
        params["E_L"],#::Float64     
        params["tau_syn"],#::Float64     
        params["g_syn"],#::Float64     
        params["spike_threshold"] #::Float64
        )

    _runparams = RunParams(
        params["I_O"],#::Float64      #offset
        params["gamma_I_t"],#::Float64      #amplitude
        params["omega_I_t"],#::Type inferred at runtime      #frequency
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
        params["tau_n_gamma"],#::Float64     osc amplitude noise time constant
        params["omegasigma"],#::Float64
        params["tau_n_omega"],#::Float64     osc frequency noise time constant
        )

    if params["nnodes"] == 1

        _wts = Wts(
                    params["w_EE"]#::Float64      excitatory auto-recurrence
        );
        
        if haskey(params, "inputalignment") && params["inputalignment"] == false
            realstarts = params["pulsestarts"] # no locking of input pulses to phase:
        else
            realstarts = (2*pi/params["omega_I_t"])*round.((params["pulsestarts"]/(2*pi/params["omega_I_t"]))) .- (pi/params["omega_I_t"]);
        end


        # typical params:
        _inputparams = InputParams(
            vcat(realstarts, 1e23), #pulse starts plus a realmax pad to make things easy
            vcat(realstarts + params["pulsedurations"], 1e23),
            vcat(params["pulseheights"], 1e23),
            1, #pulsecounter::Int64
            0, #pulsewindow::Int64
            0.0 #input::Float64, initialize at zero, gets updated during sim
        )

        if haskey(params, "perturbationstarts")
            _perturbationparams = InputParams(
                vcat(params["perturbationstarts"], 1e23), #pulse starts plus a realmax pad to make things easy
                vcat(params["perturbationstarts"] + params["perturbationdurations"], 1e23),
                vcat(params["perturbationheights"], 1e23),
                1, #pulsecounter::Int64
                0, #pulsewindow::Int64
                0.0 #input::Float64, initialize at zero, gets updated during sim
            )
        else
            _perturbationparams = InputParams(
                [1e23], #pulse starts plus a realmax pad to make things easy
                [1e23],
                [1e23],
                1, #pulsecounter::Int64
                0, #pulsewindow::Int64
                0.0 #input::Float64, initialize at zero, gets updated during sim
            )
        end


        _integratornetwork = SingleIntegratorNode(
            1.0,
            zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), 
            zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), 
            spzeros(Float64, length(0:_runparams.dt:_runparams.tmax)), 
            0.09,
            0.09,
            0.78,
            0.78,
            0.0,
            0.0,
            -64.0,
            -64.0,
            zeros(3),
            0,
            0.0
        );

        # Call primary implementation
        if params["simtype"] == :synapsephase
            return _oscintegrator_synapsephaseshift(_wts, _netparams, _runparams, _simvals, 
                                                    _integratornetwork, _noiseparams, 
                                                    _inputparams, params["synapsephase"],
                                                    params["synapseamplitude"]
                                                    )

        else
            return _oscintegrator_sim(_wts, _netparams, _runparams, _simvals, _integratornetwork, _noiseparams, _inputparams, _perturbationparams) 
        end
    elseif params["nnodes"] > 1 && params["simtype"] == :alltoall

        if haskey(params, "inputalignment") && params["inputalignment"] == false
            realstarts = params["pulsestarts"] # no locking of input pulses to phase:
        else
            realstarts = (2*pi/params["omega_I_t"])*round.((params["pulsestarts"]/(2*pi/params["omega_I_t"]))) .- (pi/params["omega_I_t"]);
        end

        _inputparams = InputParams(
            vcat(realstarts, 1e23), #pulse starts plus a realmax pad to make things easy
            vcat(realstarts + params["pulsedurations"], 1e23),
            vcat(params["pulseheights"], 1e23),
            1, #pulsecounter::Int64
            0, #pulsewindow::Int64
            0.0 #input::Float64, initialize at zero, gets updated during sim
        )

        if haskey(params, "phasedistribution")
            phaseexpr = Expr(:call, params["phasedistribution"], params["phaseparam1"], params["phaseparam2"])
            phasevect = rand(eval(phaseexpr), params["nnodes"])
        else
            phasevect = zeros(params["nnodes"])
        end

        
        if haskey(params, "sparsenet_N")
            _integratornetwork = SingleIntegratorNet(
                getsparseweights(params["sparsenet_N"], params["w_EE"], params["nnodes"]),
                zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
                spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
                0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
                0.09*ones(params["nnodes"]), #n::Array{Float64,1}
                0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
                0.78*ones(params["nnodes"]), #h::Array{Float64,1}
                0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
                0.0*ones(params["nnodes"]), #s::Array{Float64,1}
                -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
                -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
                zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
                zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
                zeros(params["nnodes"]), #noise::Array{Float64,1}
                phasevect #phase::Array{Float64,1}
            )
        else
            _integratornetwork = SingleIntegratorNet(
                params["w_EE"] * 1/params["nnodes"] * ones(params["nnodes"], params["nnodes"]) + rand( Normal(0.0,params["weightssigma"]), params["nnodes"], params["nnodes"]), #EE::Array{Float64,2}
                zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
                spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
                0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
                0.09*ones(params["nnodes"]), #n::Array{Float64,1}
                0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
                0.78*ones(params["nnodes"]), #h::Array{Float64,1}
                0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
                0.0*ones(params["nnodes"]), #s::Array{Float64,1}
                -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
                -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
                zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
                zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
                zeros(params["nnodes"]), #noise::Array{Float64,1}
                phasevect #phase::Array{Float64,1}
            )
        end


        # Call primary implementation
        return _oscintegrator_sim_alltoall(_netparams, _runparams, _simvals, _integratornetwork, _noiseparams, _inputparams) 

    elseif params["nnodes"] > 1 && params["simtype"] == :sequence

        if haskey(params, "inputalignment") && params["inputalignment"] == false
            realstarts = params["pulsestarts"] # no locking of input pulses to phase
        else
            realstarts = (2*pi/params["omega_I_t"])*round.((params["pulsestarts"]/(2*pi/params["omega_I_t"]))) .- (pi/params["omega_I_t"]);
        end

        if haskey(params, "phasedistribution")
            phaseexpr = Expr(:call, params["phasedistribution"], params["phaseparam1"], params["phaseparam2"])
            phasevect = rand(eval(phaseexpr), params["nnodes"])
        else
            phasevect = zeros(params["nnodes"])
        end


        _inputparams = SequenceInputParams(
            vcat(realstarts, 1e23), #pulse starts plus a realmax pad to make things easy
            vcat(realstarts + params["pulsedurations"], 1e23),
            hcat(generateinputbumps(params["nnodes"], params["pulselocs"], params["pulseheights"], params["pulsewidths"]), ones(params["nnodes"], size(params["pulseheights"],2)) * 1e23),
            1, #pulsecounter::Int64
            0, #pulsewindow::Int64
            zeros(params["nnodes"]), #input::Array{Float64,1}
            zeros(params["nnodes"]) #inputnull::Array{Float64,1}
        )

        
        _integratornetwork = SingleIntegratorNet(
            getsequenceweightmatrix(params["weightoffset_A"], params["weightmult_B"], params["weight_pimult"], params["nnodes"], params["weight_posoffset"], params["heavysideoffset"]) +  rand( Normal(0.0,params["weightssigma"]), params["nnodes"], params["nnodes"]), #EE::Array{Float64,2}
            zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
            spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
            0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
            0.09*ones(params["nnodes"]), #n::Array{Float64,1}
            0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
            0.78*ones(params["nnodes"]), #h::Array{Float64,1}
            0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
            0.0*ones(params["nnodes"]), #s::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
            zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
            zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
            zeros(params["nnodes"]), #noise::Array{Float64,1}
            phasevect #phase::Array{Float64,1}
        )


        # Call primary implementation
        return _oscintegrator_sim_alltoall(_netparams, _runparams, _simvals, _integratornetwork, _noiseparams, _inputparams)

    elseif params["nnodes"] > 1 && params["simtype"] == :ring


        if haskey(params, "inputalignment") && params["inputalignment"] == false
            realstarts = params["pulsestarts"] # no locking of input pulses to phase:
        else
            realstarts = (2*pi/params["omega_I_t"])*round.((params["pulsestarts"]/(2*pi/params["omega_I_t"]))) .- (pi/params["omega_I_t"]);
        end

        if haskey(params, "phasedistribution")
            phaseexpr = Expr(:call, params["phasedistribution"], params["phaseparam1"], params["phaseparam2"])
            phasevect = rand(eval(phaseexpr), params["nnodes"])
        else
            phasevect = zeros(params["nnodes"])
        end

        _inputparams = SequenceInputParams(
            vcat(realstarts, 1e23), #pulse starts plus a realmax pad to make things easy
            vcat(realstarts + params["pulsedurations"], 1e23),
            hcat(generateinputbumps(params["nnodes"], params["pulselocs"], params["pulseheights"], params["pulsewidths"]), ones(params["nnodes"], size(params["pulseheights"],2)) * 1e23),
            1, #pulsecounter::Int64
            0, #pulsewindow::Int64
            zeros(params["nnodes"]), #input::Array{Float64,1}
            zeros(params["nnodes"]) #inputnull::Array{Float64,1}
        )

        
        if haskey(params, "heavysideoffset")
            hvyoff = params["heavysideoffset"]
            println("yep")
        else
            hvyoff = params["nnodes"]*2
        end

        _integratornetwork = SingleIntegratorNet(
            getringweightmatrix(params["weightoffset_A"], params["weightmult_B"], params["weight_pimult"], params["nnodes"], hvyoff) +  rand( Normal(0.0,params["weightssigma"]), params["nnodes"], params["nnodes"]), #EE::Array{Float64,2}
            zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
            spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
            0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
            0.09*ones(params["nnodes"]), #n::Array{Float64,1}
            0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
            0.78*ones(params["nnodes"]), #h::Array{Float64,1}
            0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
            0.0*ones(params["nnodes"]), #s::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
            zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
            zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
            zeros(params["nnodes"]), #noise::Array{Float64,1}
            phasevect #phase::Array{Float64,1}
        )

        # Call primary implementation
        return _oscintegrator_sim_alltoall(_netparams, _runparams, _simvals, _integratornetwork, _noiseparams, _inputparams) #_noiseparams,
    end

end


function oscintegrator_stdystate(params)

    _netparams = NetParams(
        params["g_Na"],#::Float64
        params["g_K"],#::Float64
        params["g_L"],#::Float64
        params["E_Na"],#::Float64
        params["phi"],#::Float64     
        params["E_K"],#::Float64     
        params["C_m"],#::Float64     
        params["E_L"],#::Float64     
        params["tau_syn"],#::Float64     
        params["g_syn"],#::Float64     
        params["spike_threshold"]  #::Float64
        )

    _runparams = RunParamsStdyState(
        params["I_O"],#::Float64      #offset
        params["I_syn"],#::Float64      #synaptic input
        params["gamma_I_t"],#::Float64      #osc amplitude
        params["omega_I_t"],#::Float64      #osc freq
        params["dt"],#::Float64        #time step in ms
        params["tmax"],#::Float64     #time in ms to run simulation
        params["downsamplespacing"]#::Int64 #downsample for plotting
        )

    _simvals = SimVals(
        0:_runparams.dt:_runparams.tmax, #t_vect
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #plot_t_vect::Array{Float64,1}
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))) #input_t_vect::Array{Float64,1}
        )

    _wts = Wts(
                params["w_EE"]#::Float64      excitatory auto-recurrence
                );

    _integratornetwork = SingleIntegratorNode(
        1.0,
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))),
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), 
        spzeros(Float64, length(0:_runparams.dt:_runparams.tmax)), 
        0.09,
        0.09,
        0.78,
        0.78,
        0.0,
        0.0,
        -64.0,
        -64.0,
        zeros(3),
        0,
        0.0
        );

    # Call primary implementation
    _oscintegrator_stdystate(_wts, _netparams, _runparams, _simvals, _integratornetwork)



end



function oscint_ata_stdystate(params)

    _netparams = NetParams(
        params["g_Na"],#::Float64
        params["g_K"],#::Float64
        params["g_L"],#::Float64
        params["E_Na"],#::Float64
        params["phi"],#::Float64     
        params["E_K"],#::Float64     
        params["C_m"],#::Float64     
        params["E_L"],#::Float64     
        params["tau_syn"],#::Float64     
        params["g_syn"],#::Float64     
        params["spike_threshold"]  #::Float64
    )

    _runparams = RunParamsStdyState(
        params["I_O"],#::Float64      #offset
        params["I_syn"],#::Float64      #synaptic input
        params["gamma_I_t"],#::Float64      #osc amplitude
        params["omega_I_t"],#::Type inferred at runtime      #frequency
        params["dt"],#::Float64        #time step in ms
        params["tmax"],#::Float64     #time in ms to run simulation
        params["downsamplespacing"]#::Int64 #downsample for plotting
    )

    _simvals = SimVals(
        0:_runparams.dt:_runparams.tmax, #t_vect
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #plot_t_vect::Array{Float64,1}
        zeros(Float64, convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))) #input_t_vect::Array{Float64,1}
    )

    if haskey(params, "phasedistribution")
        phaseexpr = Expr(:call, params["phasedistribution"], params["phaseparam1"], params["phaseparam2"])
        phasevect = rand(eval(phaseexpr), params["nnodes"])
    else
        phasevect = zeros(params["nnodes"])
    end

    if haskey(params, "sparsenet_N")
        _integratornetwork = SingleIntegratorNet(
            getsparseweights(params["sparsenet_N"], params["w_EE"], params["nnodes"]),
            zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
            spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
            0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
            0.09*ones(params["nnodes"]), #n::Array{Float64,1}
            0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
            0.78*ones(params["nnodes"]), #h::Array{Float64,1}
            0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
            0.0*ones(params["nnodes"]), #s::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
            zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
            zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
            zeros(params["nnodes"]), #noise::Array{Float64,1}
            phasevect #phase::Array{Float64,1}
        )
    else
        _integratornetwork = SingleIntegratorNet(
            params["w_EE"] * 1/params["nnodes"] * ones(params["nnodes"], params["nnodes"]),
            zeros(Float64, params["nnodes"], convert(Int64, round(length(0:_runparams.dt:_runparams.tmax)/_runparams.downsamplespacing))), #s_mat::Array{Float64,1}
            spzeros(Float64, params["nnodes"], length(0:_runparams.dt:_runparams.tmax)), #spikes_mat::SparseMatrix
            0.09*ones(params["nnodes"]), #fn::Array{Float64,1}
            0.09*ones(params["nnodes"]), #n::Array{Float64,1}
            0.78*ones(params["nnodes"]), #fh::Array{Float64,1}
            0.78*ones(params["nnodes"]), #h::Array{Float64,1}
            0.0*ones(params["nnodes"]), #fs::Array{Float64,1}
            0.0*ones(params["nnodes"]), #s::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #fv::Array{Float64,1}
            -64.0*ones(params["nnodes"]), #v::Array{Float64,1}
            zeros(params["nnodes"], 3), #vstore::Array{Float64,2}
            zeros(Int64, params["nnodes"]), #spikeyesno::Array{Int64,1}
            zeros(params["nnodes"]), #noise::Array{Float64,1}
            phasevect #phase::Array{Float64,1}
        )
    end


    # Call primary implementation
    return _oscint_ata_stdystate(_netparams, _runparams, _simvals, _integratornetwork) 

end