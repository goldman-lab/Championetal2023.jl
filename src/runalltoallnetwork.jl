
function _oscintegrator_sim_alltoall(
                            netparams::NetParams,
                            runparams::RunParams,
                            simvals::SimVals,
                            integratornet::SingleIntegratorNet,
                            noiseparams::NoiseParams,
                            inputparams::InputContainer
                            )

    a=0
    b=0
    omeganoisevar = ones(length(integratornet.v))
    gammanoisevar = ones(length(integratornet.v))
    dt = runparams.dt
    oscinput = 0.0
    netsize = length(integratornet.v)

    kn1 = zeros(length(integratornet.n))
    kh1 = zeros(length(integratornet.h))
    ks1 = zeros(length(integratornet.s))
    kv1 = zeros(length(integratornet.v))

    kn2 = zeros(length(integratornet.n))
    kh2 = zeros(length(integratornet.h))
    ks2 = zeros(length(integratornet.s))
    kv2 = zeros(length(integratornet.v))

    kn3 = zeros(length(integratornet.n))
    kh3 = zeros(length(integratornet.h))
    ks3 = zeros(length(integratornet.s))
    kv3 = zeros(length(integratornet.v))

    kn4 = zeros(length(integratornet.n))
    kh4 = zeros(length(integratornet.h))
    ks4 = zeros(length(integratornet.s))
    kv4 = zeros(length(integratornet.v))
    
    for t = simvals.t_vect
        a+=1

        #Get external input from user-supplied parameters
        updateinput!(inputparams, t)

        #Calculate whether there has been a spike
        integratornet.spikeyesno, integratornet.vstore = netspiketest(  integratornet.spikeyesno,
                                                                        integratornet.vstore,
                                                                        integratornet.v,
                                                                        netparams.spike_threshold
        );

        #Store spike time if there was a spike
        spikeidxs =  integratornet.spikeyesno .== 1
        integratornet.spikes_mat[spikeidxs, a] .= t
        

        #Update input noise term
        if noiseparams.sigma > 0
            integratornet.noise =  integratornet.noise + (-integratornet.noise/noiseparams.tau_n)*runparams.dt +
                                    noiseparams.sigma*randn(netsize)*sqrt(runparams.dt/noiseparams.tau_n);
        end

        # Update frequency noising
        if t>1 
            if noiseparams.gammasigma > 0
                gammanoisevar = gammanoisevar + ((1 .- gammanoisevar)/noiseparams.tau_n_gamma)*runparams.dt + (noiseparams.gammasigma)*randn(netsize)*sqrt(runparams.dt/noiseparams.tau_n_gamma)
            elseif noiseparams.omegasigma > 0
                omeganoisevar = omeganoisevar + ((1 .- omeganoisevar)/noiseparams.tau_n_omega)*runparams.dt + (noiseparams.omegasigma/(t/1000))*randn(netsize)*sqrt(runparams.dt/noiseparams.tau_n_omega)
            end
        end
        oscinput = I(t, runparams.gamma_I_t*gammanoisevar, runparams.omega_I_t.*omeganoisevar, integratornet.phase) 


        
        ### Runge steps ###
        #step 1
        kn1 = dt * dndt(integratornet.n, netparams.phi, integratornet.v)
        kh1 = dt * dhdt(integratornet.h, netparams.phi, integratornet.v)
        ks1 = dt * dsdt(integratornet.s, netparams.tau_syn)
        kv1 = dt * dvdt(integratornet.v, 
                        netparams.g_L, netparams.E_L, 
                        integratornet.h, netparams.g_Na, netparams.E_Na, 
                        integratornet.n, netparams.g_K, netparams.E_K,
                        integratornet.s, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, integratornet.noise, inputparams.input,
                        netparams.C_m
                        )
        #step 2
        kn2 = dt * dndt(integratornet.n + kn1/2, netparams.phi, integratornet.v + kv1/2)
        kh2 = dt * dhdt(integratornet.h + kh1/2, netparams.phi, integratornet.v + kv1/2)
        ks2 = dt * dsdt(integratornet.s + ks1/2, netparams.tau_syn)
        kv2 = dt * dvdt(integratornet.v + kv1/2, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh1/2, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn1/2, netparams.g_K, netparams.E_K,
                        integratornet.s + ks1/2, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, integratornet.noise, inputparams.input,
                        netparams.C_m
                        )
        #step 3
        kn3 = dt * dndt(integratornet.n + kn2/2, netparams.phi, integratornet.v + kv2/2)
        kh3 = dt * dhdt(integratornet.h + kh2/2, netparams.phi, integratornet.v + kv2/2)
        ks3 = dt * dsdt(integratornet.s + ks2/2, netparams.tau_syn)
        kv3 = dt * dvdt(integratornet.v + kv2/2, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh2/2, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn2/2, netparams.g_K, netparams.E_K,
                        integratornet.s + ks2/2, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, integratornet.noise, inputparams.input,
                        netparams.C_m
                        )
        #step 4
        kn4 = dt * dndt(integratornet.n + kn3, netparams.phi, integratornet.v + kv3)
        kh4 = dt * dhdt(integratornet.h + kh3, netparams.phi, integratornet.v + kv3)
        ks4 = dt * dsdt(integratornet.s + ks3, netparams.tau_syn)
        kv4 = dt * dvdt(integratornet.v + kv3, netparams.g_L, netparams.E_L, 
                        integratornet.h + kh3, netparams.g_Na, netparams.E_Na, 
                        integratornet.n + kn3, netparams.g_K, netparams.E_K,
                        integratornet.s + ks3, netparams.g_syn, integratornet.EE,
                        runparams.I_O, oscinput, integratornet.noise, inputparams.input,
                        netparams.C_m
                        )
        #reduction
        integratornet.n = integratornet.n + (kn1 + 2*kn2 + 2*kn3 + kn4)/6
        integratornet.h = integratornet.h + (kh1 + 2*kh2 + 2*kh3 + kh4)/6
        integratornet.s = integratornet.s + (ks1 + 2*ks2 + 2*ks3 + ks4)/6 + integratornet.spikeyesno/netparams.tau_syn;
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
