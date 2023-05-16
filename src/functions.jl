beta_n(v)=0.125*exp( -((v+44)/80) )

alpha_n(v)=-0.01*( (v+34)/(exp(-0.1*(v+34))-1) )

beta_h(v)=1/(exp( -0.1*(v+28) )+1)

alpha_h(v)=0.07*exp( -((v+58)/20) )

beta_m(v)=4*exp( -((v+60)/18) )

alpha_m(v)=-0.1*( (v+35)/(exp(-0.1*(v+35))-1) )

tau_n(v)=1 ./ ( alpha_n.(v) + beta_n.(v) )

tau_h(v)=1 ./ ( alpha_h.(v) + beta_h.(v) )

n_inf(v)=alpha_n.(v) ./ ( alpha_n.(v) + beta_n.(v) )

h_inf(v)=alpha_h.(v) ./ ( alpha_h.(v) + beta_h.(v) )

m_inf(v)=alpha_m.(v) ./ ( alpha_m.(v) + beta_m.(v) )

I(t, gamma, omega)=gamma*cos(omega*t)
I(t, gamma, omega, phase)=gamma .* cos.(omega*t - phase)
    
I_syn(s, g_syn, w_ee) = w_ee*s * g_syn

I_L(v, g_L, E_L)=g_L*(v .- E_L)

I_k(v, n, g_k, E_K)=g_k*(n.^4) .* (v .- E_K)

I_na(v, h, g_Na, E_Na)=g_Na*(m_inf(v).^3) .* h .* (v .- E_Na)

function spiketest(vstore, v, spike_thresh)

    vstore[1] = vstore[2]
    vstore[2] = vstore[3]
    vstore[3] = max(0, sign(v-spike_thresh))*v

    spiketest = 0
    if vstore[1] > 0
        if vstore[2]>vstore[1] && vstore[3]<vstore[2]
            spiketest = 1
        end
    end

    return spiketest, vstore
end

function netspiketest(spiketestvect, vstore, v, spike_thresh)

    vstore = circshift(vstore, (0,-1))
    vstore[:,3] = max.(0, sign.(v .- spike_thresh)) .* v

    for i in 1:length(v)
        if vstore[i,1] > 0
            if vstore[i,2] > vstore[i,1] && vstore[i,3] < vstore[i,2]
                spiketestvect[i] = 1
            else
                spiketestvect[i] = 0
            end
        else
            spiketestvect[i] = 0
        end
    end

    return spiketestvect, vstore
end

function getedges(omega, tmax)
    edgeholder = spzeros(100000)
    t = 0
    n=1
    a=0

    while t < tmax
        if a>0
            edgeholder[a] = t
        end

        a+=1
        t = n*(pi/2)/omega
        n += 4
    end

    return nonzeros(edgeholder)
end


function updateinput!(inputstructure::InputParams, currenttime)
    # inputstructure.pulsecounter initialized at 1
    if currenttime >= inputstructure.pulsestarts[inputstructure.pulsecounter] &&
        currenttime <= inputstructure.pulsestops[inputstructure.pulsecounter]

        inputstructure.input = inputstructure.pulseheights[inputstructure.pulsecounter]
        inputstructure.pulsewindow = 1
    end

    if inputstructure.pulsewindow == 1 &&
        currenttime > inputstructure.pulsestops[inputstructure.pulsecounter]

        inputstructure.pulsewindow = 0
        inputstructure.pulsecounter += 1
        inputstructure.input = 0.0
    end

    return nothing

end

function updateinput!(inputstructure::SequenceInputParams, currenttime)
    # inputstructure.pulsecounter initialized at 1
    if currenttime >= inputstructure.pulsestarts[inputstructure.pulsecounter] &&
        currenttime <= inputstructure.pulsestops[inputstructure.pulsecounter]

        inputstructure.input = inputstructure.pulseheights[:,inputstructure.pulsecounter]
        inputstructure.pulsewindow = 1
    end

    if inputstructure.pulsewindow == 1 &&
        currenttime > inputstructure.pulsestops[inputstructure.pulsecounter]

        inputstructure.pulsewindow = 0
        inputstructure.pulsecounter += 1
        inputstructure.input = inputstructure.inputnull
    end

    return nothing

end

function generateinputbumps(nnodes, bumplocs, bumpheights, bumpwidths)

    bumpstorage = zeros(nnodes, size(bumpheights,2))

    if maximum(bumplocs + round.(bumpwidths ./ 2)) <= nnodes && minimum(bumplocs - round.(bumpwidths ./ 2)) > 0
        for j in 1:size(bumpheights,2)

            holder = zeros(nnodes)
            for i in 1:size(bumplocs,1)
                holder[bumplocs[i,j]-round(Int64,bumpwidths[i,j]/2):bumplocs[i,j]+round(Int64,bumpwidths[i,j]/2)] .= bumpheights[i,j]
            end

            bumpstorage[:,j] = holder

        end

        return bumpstorage

    else
        error("bump runs off the edge of the network. this isn't actually a ring.")
    end

end


function getsequenceweights(i, weightoffset_A, weightmult_B, weight_pimult, nnodes, weight_posoffset, heavysideoffset)

    return  (
            weightoffset_A .+
            weightmult_B *
            cos.(
                    (
                    weight_pimult*pi*( i .- collect(1:nnodes) )
                    )/nnodes .+ 0.1
                )
            ) .*
            max.(
                0,
                sign.( heavysideoffset .- abs.((i+1) .- collect(1:nnodes)) )
                )

end

function getsequenceweightmatrix(weightoffset_A, weightmult_B, weight_pimult, nnodes, weight_posoffset, heavysideoffset)

    w_mat = zeros(nnodes, nnodes);

    for i in 1:nnodes

        weightvect_i =  (  weightoffset_A .+ weightmult_B *
                            cos.((weight_pimult*pi*( i .- collect(1:nnodes) ))/nnodes .+ 0.1)
                        ) .*
                        max.( 0,
                              sign.( heavysideoffset .- abs.((i+1) .- collect(1:nnodes)) )
                        )
        w_mat[i,:] = weightvect_i'

    end

    return w_mat

end

function getringweights(i, weightoffset_A, weightmult_B, weight_pimult, nnodes)

    return  weightoffset_A .+
            weightmult_B *
            cos.(
                    (
                    weight_pimult*pi*( i .- collect(1:nnodes) )
                    )/nnodes
                )

end


function getringweightmatrix(weightoffset_A, weightmult_B, weight_pimult, nnodes, heavysideoffset)

    w_mat = zeros(nnodes, nnodes);
    for i in 1:nnodes

        weightvect_i =   (weightoffset_A .+ 
                        weightmult_B * 
                        cos.((weight_pimult*pi*( i .- collect(1:nnodes) ))/nnodes)) .* 
                        max.( 0,
                            sign.( heavysideoffset .- abs.((i+1) .- collect(1:nnodes)) )
                            )

        w_mat[i,:] = weightvect_i'

    end

    return w_mat

end

function getsparseweights(nafferents, w_EE, nnodes)
    W = zeros(nnodes, nnodes)

    for i in 1:nnodes
        w_vect = zeros(nnodes)
        sparsepicks = rand(1:nnodes, nafferents)
        w_vect[sparsepicks] .= w_EE/nafferents
        W[i,:] = w_vect'
    end
    
    return W
end

dndt(xn, phi, v) = phi*(n_inf(v) - xn) ./ tau_n(v);
        
dhdt(xh, phi, v) = phi*(h_inf(v) - xh) ./ tau_h(v)

dsdt(xs, tau_syn) = -xs ./ tau_syn; 

function dvdt(xv, g_L, E_L, h, g_Na, E_Na, n, g_K, E_K, s, g_syn, w_EE, I_O, oscinput, noise, input, C_m)

    return (
            -I_L(xv, g_L, E_L) -
            I_na(xv, h, g_Na, E_Na) -
            I_k(xv, n, g_K, E_K) +
            I_syn(s, g_syn, w_EE) .+
            I_O + #if net: constant, Float64
            oscinput .+ #if net: vector
            noise .+ #if net: vector
            input #can be Float or vector
            ) / C_m;

end