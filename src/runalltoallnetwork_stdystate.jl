
function _oscint_ata_stdystate(
                                netparams::NetParams,
                                runparams::RunParamsStdyState,
                                simvals::SimVals,
                                integratornet::SingleIntegratorNet,
                                )

    a=0
    b=0
    dt = runparams.dt
    oscinput = 0.0
    netsize = length(integratornet.v)

    kn1 = zeros(length(integratornet.n))
    kh1 = zeros(length(integratornet.h))
    kv1 = zeros(length(integratornet.v))

    kn2 = zeros(length(integratornet.n))
    kh2 = zeros(length(integratornet.h))
    kv2 = zeros(length(integratornet.v))

    kn3 = zeros(length(integratornet.n))
    kh3 = zeros(length(integratornet.h))
    kv3 = zeros(length(integratornet.v))

    kn4 = zeros(length(integratornet.n))
    kh4 = zeros(length(integratornet.h))
    kv4 = zeros(length(integratornet.v))

    svect = runparams.I_ext*ones(netsize)
    
    for t = simvals.t_vect
        a+=1

        #Calculate whether there has been a spike
        integratornet.spikeyesno, integratornet.vstore = netspiketest(  integratornet.spikeyesno,
                                                                        integratornet.vstore,
                                                                        integratornet.v,
                                                                        netparams.spike_threshold
        );

        #Store spike time if there was a spike
        spikeidxs =  integratornet.spikeyesno .== 1
        integratornet.spikes_mat[spikeidxs, a] .= t
        

        oscinput = I(t, runparams.gamma_I_t, runparams.omega_I_t, integratornet.phase) 

        
        ### Runge steps ###
        #step 1
        kn1 = dt * dndt(integratornet.n, netparams.phi, integratornet.v)
        kh1 = dt * dhdt(integratornet.h, netparams.phi, integratornet.v)
        
        kv1 = dt * dvdt(integratornet.v, 
                        netparams.g_L, netparams.E_L, 
                        integratornet.h, netparams.g_Na, netparams.E_Na, 
                        integratornet.n, netparams.g_K, netparams.E_K,
                        svect, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, 0.0, 0.0, #swapped 0.0 into integratornet.noise and inputparams.input
                        netparams.C_m
                        )
        #step 2
        kn2 = dt * dndt(integratornet.n + kn1/2, netparams.phi, integratornet.v + kv1/2)
        kh2 = dt * dhdt(integratornet.h + kh1/2, netparams.phi, integratornet.v + kv1/2)
        
        kv2 = dt * dvdt(integratornet.v + kv1/2, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh1/2, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn1/2, netparams.g_K, netparams.E_K,
                        svect, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, 0.0, 0.0, #swapped 0.0 into integratornet.noise and inputparams.input
                        netparams.C_m
                        )
        #step 3
        kn3 = dt * dndt(integratornet.n + kn2/2, netparams.phi, integratornet.v + kv2/2)
        kh3 = dt * dhdt(integratornet.h + kh2/2, netparams.phi, integratornet.v + kv2/2)
        
        kv3 = dt * dvdt(integratornet.v + kv2/2, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh2/2, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn2/2, netparams.g_K, netparams.E_K,
                        svect, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, 0.0, 0.0, #swapped 0.0 into integratornet.noise and inputparams.input
                        netparams.C_m
                        )
        #step 4
        kn4 = dt * dndt(integratornet.n + kn3, netparams.phi, integratornet.v + kv3)
        kh4 = dt * dhdt(integratornet.h + kh3, netparams.phi, integratornet.v + kv3)
        
        kv4 = dt * dvdt(integratornet.v + kv3, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh3, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn3, netparams.g_K, netparams.E_K,
                        svect, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, 0.0, 0.0, #swapped 0.0 into integratornet.noise and inputparams.input
                        netparams.C_m
                        )
        #reduction
        integratornet.n = integratornet.n + (kn1 + 2*kn2 + 2*kn3 + kn4)/6
        integratornet.h = integratornet.h + (kh1 + 2*kh2 + 2*kh3 + kh4)/6
        integratornet.s = svect
        integratornet.v = integratornet.v + (kv1 + 2*kv2 + 2*kv3 + kv4)/6



        # Downsample storage for plotting
        if mod(a,runparams.downsamplespacing) == 0
            b += 1;
            # integratornet.v_vect[b] = integratornet.v; vmat taken out to save space
            integratornet.s_mat[:,b] = integratornet.s; 
            simvals.plot_t_vect[b] = simvals.t_vect[a]; 
        end

    end

    return simvals, integratornet

end
