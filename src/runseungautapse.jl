function sam_dndt(n, phi, v)
    phi*( sam_alpha_n(v)*(1-n) - sam_beta_n(v)*n)
end

function sam_dhdt(h, phi, v)
    phi*(sam_alpha_h(v)*(1-h) - sam_beta_h(v)*h)
end

function sam_dbdt(b, tau, v)
    (sam_b_inf(v) - b)/tau
end

function sam_dsdt(s, tau, alpha_s, theta_s, sigma_s, v)
    ( -s + alpha_s*(1-s)*sam_sigma_k(v, theta_s, sigma_s) ) / tau
end

function sam_dvmemdt(   v, g_L, E_L,
                    h, g_Na, E_Na,
                    n, g_K, E_K,
                    b, g_A,
                    w_S, s, w_O, O_s, E_E,
                    currentinjection,
                    C_membrane
                )

    return (
            -sam_I_L(v, g_L, E_L) -
            sam_I_na(v, h, g_Na, E_Na) -
            sam_I_k(v, n, g_K, E_K) -
            sam_I_A(v, b, g_A, E_K) + currentinjection +
            sam_I_app_S(w_S, s, w_O, O_s, v, E_E)
            ) / C_membrane;

end

# Leak current
function sam_I_L(v, g_L, E_L)
    g_L*(v - E_L)
end

# Sodium current I_na
function sam_I_na(v, h, g_Na, E_Na)
    g_Na*(sam_m_inf(v)^3)*h*(v-E_Na)
end

function sam_m_inf(v)
    sam_alpha_m(v)/( sam_alpha_m(v) + sam_beta_m(v) )
end

function sam_alpha_m(v)
     ((v+30)/10)/(1 - exp(-(v+30)/10))
end

function sam_beta_m(v)
   4*exp( -(v+55)/18 )
end

function sam_alpha_h(v)
    0.07*exp( -(v+44)/20 )
end

function sam_beta_h(v)
    1/(exp( -(v+14)/10 )+1)
end

# Delayed rectifier potassium current I_k
function sam_I_k(v, n, g_k, E_K)
    g_k*(n^4)*(v-E_K)
end

function sam_alpha_n(v)
    ((v+34)/100)/(1 - exp(-(v+34)/10))
end

function sam_beta_n(v)
    0.125*exp( -(v+44)/80 )
end

#A-type potassium current I_A
function sam_I_A(v, b, g_A, E_K)
    g_A*(sam_a_inf(v)^3)*b*(v-E_K)
end

function sam_a_inf(v)
    1/( exp( -(v+50)/20 ) + 1 )
end

function sam_b_inf(v)
    1/( exp( (v+80)/6 ) + 1 )
end

# other terms
function sam_sigma_k(v, theta_s, sigma_s)
    1 / ( 1 + exp( -(v-theta_s)/sigma_s ) )
end

function sam_I_app_S(w_S, s_S, w_O, s_O, v_S, E_E)
    -(w_S*s_S + w_O*s_O)*(v_S - E_E)
end


function _seungautapse_stdystatesim(
    netparams::SeungNetParams,
    runparams::RunParams,
    simvals::SimVals,
    seungnode_M::SeungNode,
    noiseparams::NoiseParams,
    inputparams::InputParams,
    synapseactivity
    )

    kn1 = 0.0
    kh1 = 0.0
    kb1 = 0.0
    kv1 = 0.0

    kn2 = 0.0
    kh2 = 0.0
    kb2 = 0.0
    kv2 = 0.0

    kn3 = 0.0
    kh3 = 0.0
    kb3 = 0.0
    kv3 = 0.0

    kn4 = 0.0
    kh4 = 0.0
    kb4 = 0.0
    kv4 = 0.0

    O_s = 0.0093

    a=0
    b=0
    for t = simvals.t_vect
        a+=1

        seungnode_M.spikeyesno, seungnode_M.vstore = spiketest(   seungnode_M.vstore,
                                                        seungnode_M.v,
                                                        netparams.spike_threshold
                                                        );

        if seungnode_M.spikeyesno == 1
            seungnode_M.spikes_vect[a] = t
        end

        oscinput = -runparams.I_O + I(t, runparams.gamma_I_t, runparams.omega_I_t)

        # Sodium currents:
        # function calls: dhdt(h, phi, v)
        #Delayed rectifier potassium currents:
        # function calls: dndt(n, phi, v)
        #A-type potassium currents:
        # function calls: dbdt(b, tau, v)
        #Synaptic dynamics:
        #function calls: dsdt(s, tau_s, alpha_s, theta_s, sigma_s, v)
        kh1 = runparams.dt * sam_dhdt(seungnode_M.h, netparams.phi, seungnode_M.v)
        kn1 = runparams.dt * sam_dndt(seungnode_M.n, netparams.phi, seungnode_M.v)
        kb1 = runparams.dt * sam_dbdt(seungnode_M.b, netparams.tau_A, seungnode_M.v)
        ks1 = runparams.dt * sam_dsdt(seungnode_M.s, netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v)
        kv1 = runparams.dt * sam_dvmemdt(seungnode_M.v, netparams.g_L, netparams.E_L,
                            seungnode_M.h, netparams.g_Na, netparams.E_Na,
                            seungnode_M.n, netparams.g_K, netparams.E_K,
                            seungnode_M.b, netparams.g_A,
                            netparams.w_S, synapseactivity, netparams.w_O, O_s, netparams.E_E,
                            oscinput,
                            netparams.C_m)


        kh2 = runparams.dt * sam_dhdt(seungnode_M.h + (0.5 * kh1 ), netparams.phi, seungnode_M.v + (0.5*kv1))
        kn2 = runparams.dt * sam_dndt(seungnode_M.n + (0.5 * kn1 ), netparams.phi, seungnode_M.v + (0.5*kv1))
        kb2 = runparams.dt * sam_dbdt(seungnode_M.b + (0.5 * kb1 ), netparams.tau_A, seungnode_M.v + (0.5*kv1))
        ks2 = runparams.dt * sam_dsdt(seungnode_M.s + (0.5 * ks1 ), netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + (0.5*kv1))
        kv2 = runparams.dt * sam_dvmemdt(seungnode_M.v + (0.5 * kv1 ), netparams.g_L, netparams.E_L,
                            seungnode_M.h + (0.5 * kh1 ), netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + (0.5 * kn1 ), netparams.g_K, netparams.E_K,
                            seungnode_M.b + (0.5 * kb1 ), netparams.g_A,
                            netparams.w_S, synapseactivity, netparams.w_O, O_s, netparams.E_E,
                            oscinput,
                            netparams.C_m)


        kh3 = runparams.dt * sam_dhdt(seungnode_M.h + (0.5 * kh2 ), netparams.phi, seungnode_M.v + (0.5*kv2))
        kn3 = runparams.dt * sam_dndt(seungnode_M.n + (0.5 * kn2 ), netparams.phi, seungnode_M.v + (0.5*kv2))
        kb3 = runparams.dt * sam_dbdt(seungnode_M.b + (0.5 * kb2 ), netparams.tau_A, seungnode_M.v + (0.5*kv2))
        ks3 = runparams.dt * sam_dsdt(seungnode_M.s + (0.5 * ks2 ), netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + (0.5*kv2))
        kv3 = runparams.dt * sam_dvmemdt(seungnode_M.v + (0.5 * kv2 ), netparams.g_L, netparams.E_L,
                            seungnode_M.h + (0.5 * kh2 ), netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + (0.5 * kn2 ), netparams.g_K, netparams.E_K,
                            seungnode_M.b + (0.5 * kb2 ), netparams.g_A,
                            netparams.w_S, synapseactivity, netparams.w_O, O_s, netparams.E_E,
                            oscinput,
                            netparams.C_m)


        kh4 = runparams.dt * sam_dhdt(seungnode_M.h + kh3, netparams.phi, seungnode_M.v + kv3)
        kn4 = runparams.dt * sam_dndt(seungnode_M.n + kn3, netparams.phi, seungnode_M.v + kv3)
        kb4 = runparams.dt * sam_dbdt(seungnode_M.b + kb3, netparams.tau_A, seungnode_M.v + kv3)
        ks4 = runparams.dt * sam_dsdt(seungnode_M.s + ks3, netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + kv3)
        kv4 = runparams.dt * sam_dvmemdt(seungnode_M.v + kv3, netparams.g_L, netparams.E_L,
                            seungnode_M.h + kh3, netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + kn3, netparams.g_K, netparams.E_K,
                            seungnode_M.b + kb3, netparams.g_A,
                            netparams.w_S, synapseactivity, netparams.w_O, O_s, netparams.E_E,
                            oscinput,
                            netparams.C_m)

        seungnode_M.h =  seungnode_M.h + (kh1 + 2 * kh2 + 2 * kh3 + kh4) / 6
        seungnode_M.n = seungnode_M.n + (kn1 + 2 * kn2 + 2 * kn3 + kn4) / 6  
        seungnode_M.b =  seungnode_M.b + (kb1 + 2 * kb2 + 2 * kb3 + kb4) / 6
        seungnode_M.s =  seungnode_M.s + (ks1 + 2 * ks2 + 2 * ks3 + ks4) / 6
        seungnode_M.v =  seungnode_M.v + (kv1 + 2 * kv2 + 2 * kv3 + kv4) / 6



        # Downsample storage for plotting
        if mod(a,runparams.downsamplespacing) == 0
            b += 1;
            seungnode_M.v_vect[b] = seungnode_M.v; 
            seungnode_M.s_vect[b] = seungnode_M.s; 
            simvals.plot_t_vect[b] = simvals.t_vect[a]; 
        end

    end

    return simvals, seungnode_M

end

function _seungautapse_sim(
    netparams::SeungNetParams,
    runparams::RunParams,
    simvals::SimVals,
    seungnode_M::SeungNode,
    noiseparams::NoiseParams,
    inputparams::InputParams,
    )

    kn1 = 0.0
    kh1 = 0.0
    kb1 = 0.0
    kv1 = 0.0

    kn2 = 0.0
    kh2 = 0.0
    kb2 = 0.0
    kv2 = 0.0

    kn3 = 0.0
    kh3 = 0.0
    kb3 = 0.0
    kv3 = 0.0

    kn4 = 0.0
    kh4 = 0.0
    kb4 = 0.0
    kv4 = 0.0

    O_s = 0.0093

    srate = 0.0
    sratevect = zeros(length(simvals.plot_t_vect))

    a=0
    b=0
    for t = simvals.t_vect
        a+=1

        seungnode_M.spikeyesno, seungnode_M.vstore = spiketest(   seungnode_M.vstore,
                                                        seungnode_M.v,
                                                        netparams.spike_threshold
                                                        );

        if seungnode_M.spikeyesno == 1
            seungnode_M.spikes_vect[a] = t
        end

        #Get external input from user-supplied parameters
        updateinput!(inputparams, t)
        #Update oscillatory input
        oscinput = -runparams.I_O + I(t, runparams.gamma_I_t, runparams.omega_I_t)
        #Update input noise
        if noiseparams.sigma > 0
            seungnode_M.noise =  seungnode_M.noise + (-seungnode_M.noise/noiseparams.tau_n)*runparams.dt +
                        noiseparams.sigma*randn()*sqrt(runparams.dt/noiseparams.tau_n);
        end

        # Sodium currents:
        # function calls: dhdt(h, phi, v)
        #Delayed rectifier potassium currents:
        # function calls: dndt(n, phi, v)
        #A-type potassium currents:
        # function calls: dbdt(b, tau, v)
        #Synaptic dynamics:
        #function calls: dsdt(s, tau_s, alpha_s, theta_s, sigma_s, v)
        kh1 = runparams.dt * sam_dhdt(seungnode_M.h, netparams.phi, seungnode_M.v)
        kn1 = runparams.dt * sam_dndt(seungnode_M.n, netparams.phi, seungnode_M.v)
        kb1 = runparams.dt * sam_dbdt(seungnode_M.b, netparams.tau_A, seungnode_M.v)
        ks1 = runparams.dt * sam_dsdt(seungnode_M.s, netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v)
        kv1 = runparams.dt * sam_dvmemdt(seungnode_M.v, netparams.g_L, netparams.E_L,
                            seungnode_M.h, netparams.g_Na, netparams.E_Na,
                            seungnode_M.n, netparams.g_K, netparams.E_K,
                            seungnode_M.b, netparams.g_A,
                            netparams.w_S, seungnode_M.s, netparams.w_O, O_s, netparams.E_E,
                            oscinput + seungnode_M.noise + inputparams.input,
                            netparams.C_m)


        kh2 = runparams.dt * sam_dhdt(seungnode_M.h + (0.5 * kh1 ), netparams.phi, seungnode_M.v + (0.5*kv1))
        kn2 = runparams.dt * sam_dndt(seungnode_M.n + (0.5 * kn1 ), netparams.phi, seungnode_M.v + (0.5*kv1))
        kb2 = runparams.dt * sam_dbdt(seungnode_M.b + (0.5 * kb1 ), netparams.tau_A, seungnode_M.v + (0.5*kv1))
        ks2 = runparams.dt * sam_dsdt(seungnode_M.s + (0.5 * ks1 ), netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + (0.5*kv1))
        kv2 = runparams.dt * sam_dvmemdt(seungnode_M.v + (0.5 * kv1 ), netparams.g_L, netparams.E_L,
                            seungnode_M.h + (0.5 * kh1 ), netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + (0.5 * kn1 ), netparams.g_K, netparams.E_K,
                            seungnode_M.b + (0.5 * kb1 ), netparams.g_A,
                            netparams.w_S, seungnode_M.s, netparams.w_O, O_s, netparams.E_E,
                            oscinput + seungnode_M.noise + inputparams.input,
                            netparams.C_m)


        kh3 = runparams.dt * sam_dhdt(seungnode_M.h + (0.5 * kh2 ), netparams.phi, seungnode_M.v + (0.5*kv2))
        kn3 = runparams.dt * sam_dndt(seungnode_M.n + (0.5 * kn2 ), netparams.phi, seungnode_M.v + (0.5*kv2))
        kb3 = runparams.dt * sam_dbdt(seungnode_M.b + (0.5 * kb2 ), netparams.tau_A, seungnode_M.v + (0.5*kv2))
        ks3 = runparams.dt * sam_dsdt(seungnode_M.s + (0.5 * ks2 ), netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + (0.5*kv2))
        kv3 = runparams.dt * sam_dvmemdt(seungnode_M.v + (0.5 * kv2 ), netparams.g_L, netparams.E_L,
                            seungnode_M.h + (0.5 * kh2 ), netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + (0.5 * kn2 ), netparams.g_K, netparams.E_K,
                            seungnode_M.b + (0.5 * kb2 ), netparams.g_A,
                            netparams.w_S, seungnode_M.s, netparams.w_O, O_s, netparams.E_E,
                            oscinput + seungnode_M.noise + inputparams.input,
                            netparams.C_m)


        kh4 = runparams.dt * sam_dhdt(seungnode_M.h + kh3, netparams.phi, seungnode_M.v + kv3)
        kn4 = runparams.dt * sam_dndt(seungnode_M.n + kn3, netparams.phi, seungnode_M.v + kv3)
        kb4 = runparams.dt * sam_dbdt(seungnode_M.b + kb3, netparams.tau_A, seungnode_M.v + kv3)
        ks4 = runparams.dt * sam_dsdt(seungnode_M.s + ks3, netparams.tau_S, netparams.alpha_s, netparams.theta_s, netparams.sigma_s, seungnode_M.v + kv3)
        kv4 = runparams.dt * sam_dvmemdt(seungnode_M.v + kv3, netparams.g_L, netparams.E_L,
                            seungnode_M.h + kh3, netparams.g_Na, netparams.E_Na,
                            seungnode_M.n + kn3, netparams.g_K, netparams.E_K,
                            seungnode_M.b + kb3, netparams.g_A,
                            netparams.w_S, seungnode_M.s, netparams.w_O, O_s, netparams.E_E,
                            oscinput + seungnode_M.noise + inputparams.input,
                            netparams.C_m)

        seungnode_M.h =  seungnode_M.h + (kh1 + 2 * kh2 + 2 * kh3 + kh4) / 6
        seungnode_M.n = seungnode_M.n + (kn1 + 2 * kn2 + 2 * kn3 + kn4) / 6  
        seungnode_M.b =  seungnode_M.b + (kb1 + 2 * kb2 + 2 * kb3 + kb4) / 6
        seungnode_M.s =  seungnode_M.s + (ks1 + 2 * ks2 + 2 * ks3 + ks4) / 6
        seungnode_M.v =  seungnode_M.v + (kv1 + 2 * kv2 + 2 * kv3 + kv4) / 6

        srate += runparams.dt*(-0.01*srate + 0.01*seungnode_M.spikeyesno);




        # Downsample storage for plotting
        if mod(a,runparams.downsamplespacing) == 0
            b += 1;
            seungnode_M.v_vect[b] = seungnode_M.v; 
            seungnode_M.s_vect[b] = seungnode_M.s; 
            simvals.plot_t_vect[b] = simvals.t_vect[a]; 
            sratevect[b] = srate*(1000/runparams.dt)
        end

    end

    return simvals, seungnode_M, sratevect

end
