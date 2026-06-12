# File with auxiliar functions

# Function to calculate the AC power flow
# Auxiliar function to calculate the AC power flow at each branch
function Calculate_AC_Power_Flow(DCIR::DataFrame, nCIR::Int64, V::Vector{Float64}, θ::Vector{Float64}, base_MVA::Float64)
    # In MATPOWER and POWERMODELS, the TAP and SHIFT of the transformers are treated as "from" to "to"
    # This is the reason why they divide by TAP and the signal of the shift is switched

    Pik = zeros(Float64, nCIR)           # Initializing vector Pik
    Qik = zeros(Float64, nCIR)           # Initializing vector Qik
    Sik = zeros(Float64, nCIR)           # Initializing vector Sik
    Pki = zeros(Float64, nCIR)           # Initializing vector Pki
    Qki = zeros(Float64, nCIR)           # Initializing vector Qki
    Ski = zeros(Float64, nCIR)           # Initializing vector Ski
    Plosses = zeros(Float64, nCIR)       # Initializing vector Plosses
    Qlosses = zeros(Float64, nCIR)       # Initializing vector Qlosses
    circ_loading = zeros(Float64, nCIR)  # Initializing vector of percentage loading

    # Loop to calculate the power flow in the lines and transformers

    for lin = 1:nCIR
        if DCIR.l_status[lin] == true # Check if the branch is ON

            i       = DCIR.from_bus[lin]                           # Bus i (from)
            k       = DCIR.to_bus[lin]                             # Bus k (to)

            yik     = 1 / (DCIR.l_res[lin] + 1im*DCIR.l_reac[lin]) # Series admittance
            bik_sh  = DCIR.l_sh_susp[lin] / 2                      # Shunt suscpetance
            g       = real(yik)                                    # Branch series conductance
            b       = imag(yik)                                    # Branch series susceptance

            t_tap   = DCIR.t_tap[lin]                              # Transformer tap ratio (tap:1)
            t_shift = deg2rad(DCIR.t_shift[lin])                   # Transformer shift angle

            t_r     = t_tap * cos(t_shift)                         # Real number for transformers
            t_i     = t_tap * sin(t_shift)                         # Imaginary number for transformers

            angik = θ[i] - θ[k]                                    # Angular difference between bus i and k
            angki = θ[k] - θ[i]                                    # Angular difference between bus k and i

            # Active power flow from i to k
            # Pik[lin] = g*((1/t_tap) * V[i])^2 - ((1/t_tap) * V[i] * V[k]) * (g*cos(angik - t_shift) + b*sin(angik - t_shift))
            Pik[lin] = g*((1/t_tap) * V[i])^2 + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*sin(angik))

            # Reactive power flow from i to k
            # Qik[lin] = -(b + bik_sh)*((1/t_tap) * V[i])^2 + ((1/t_tap) * V[i] * V[k]) * (b*cos(angik - t_shift) - g*sin(angik - t_shift))
            Qik[lin] = -(b + bik_sh)*((1/t_tap) * V[i])^2 - (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*sin(angik))

            # Active power flow from k to i
            # Pki[lin] = g*(V[k]^2) - ((1/t_tap) * V[i] * V[k]) * (g*cos(angki + t_shift) + b*sin(angki + t_shift))
            Pki[lin] = g*(V[k]^2) + (-g*t_r-b*t_i)/t_tap^2*(V[i] * V[k] * cos(angki)) + (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * sin(angki))

            # Reactive power flow from k to i
            # Qki[lin] = -(b + bik_sh)*(V[k]^2) + ((1/t_tap) * V[i] * V[k]) * (b*cos(angki + t_shift) - g*sin(angki + t_shift))
            Qki[lin] = -(b + bik_sh)*(V[k]^2)  - (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * cos(angki)) + (-g*t_r-b*t_i)/t_tap^2*( V[i] * V[k] * sin(angki))

            # Apparent power flow from i to k
            Sik[lin] = abs(Pik[lin] + 1im*Qik[lin])

            # Apparent power flow from k to i
            Ski[lin] = abs(Pki[lin] + 1im*Qki[lin])

            # Active power losses in the branch
            Plosses[lin] = Pik[lin] + Pki[lin]

            # Reactive power losses in the branch
            Qlosses[lin] = Qik[lin] + Qki[lin]

            all_cap = [DCIR.l_cap_1[lin]; DCIR.l_cap_2[lin]; DCIR.l_cap_3[lin]]
            circ_cap = 0.0
            if any(!iszero, all_cap)
                index_cap = findfirst(!iszero, all_cap)
                circ_cap = all_cap[index_cap]
            else
                circ_cap = Inf
            end
            circ_loading[lin] = DCIR.l_status[lin] * abs(max(Sik[lin], Ski[lin])) / (circ_cap / base_MVA)

        end
        
    end

    return Pik, Qik, Sik, Pki, Qki, Ski, Plosses, Qlosses, circ_loading
end

# Function used to map the generators, circuits connected and adjacent buses for each bus
function Manage_Bus_Gen_Circ(DBUS::DataFrame, DGEN::DataFrame, DCIR::DataFrame)

    # Initialize the dictionary with empty vectors for each bus
    bus_gen_circ_dict = OrderedDict(b => Dict(:gen_ids => Int[], :gen_status => Bool[], :circ => Int[], :adj_buses => Int[], :vg_spe => 1.0) for b in DBUS.bus)

    bus_gen_circ_dict_ON = OrderedDict(b => Dict(:gen_ids => Int[], :gen_status => Bool[], :circ => Int[], :adj_buses => Int[], :vg_spe => 1.0) for b in DBUS.bus)

    # Loop through the generators and fill the dictionary
    # Associate the number of the bus with a vector containing
    # the ids of the generators
    for (gen_idx, bus_id) in enumerate(DGEN.bus)
        push!(bus_gen_circ_dict[bus_id][:gen_ids], gen_idx)
        push!(bus_gen_circ_dict[bus_id][:gen_status], DGEN.g_status[gen_idx])

        local a = unique(DGEN.vg_spe[gen_idx])
        if isempty(a)
            bus_gen_circ_dict[bus_id][:vg_spe] = 0.0
        elseif length(a) > 1
            throw(ArgumentError("The voltages specified for the generators connected to bus $bus_id must be equal."))
        else
            bus_gen_circ_dict[bus_id][:vg_spe] = a[1]
        end

        if DGEN.g_status[gen_idx] == 1
            push!(bus_gen_circ_dict_ON[bus_id][:gen_ids], gen_idx)
            push!(bus_gen_circ_dict_ON[bus_id][:gen_status], DGEN.g_status[gen_idx])
            local a = unique(DGEN.vg_spe[gen_idx])
            if isempty(a)
                bus_gen_circ_dict_ON[bus_id][:vg_spe] = 0.0
            elseif length(a) > 1
                throw(ArgumentError("The voltages specified for the generators connected to bus $bus_id must be equal."))
            else
                bus_gen_circ_dict_ON[bus_id][:vg_spe] = a[1]
            end
        end

    end

    # Loop through circuits and populate ids of circuits connectec
    # and adjacent buses
    for (i, from) in enumerate(DCIR.from_bus)
        to = DCIR.to_bus[i]
        
        # Add the number of the circuit
        push!(bus_gen_circ_dict[from][:circ], i)
        push!(bus_gen_circ_dict[to][:circ], i)

        # Add each bus as adjacent to the other
        push!(bus_gen_circ_dict[from][:adj_buses], to)
        push!(bus_gen_circ_dict[to][:adj_buses], from)

        if DCIR.l_status[i] == 1
            # Add the number of the circuit
            push!(bus_gen_circ_dict_ON[from][:circ], i)
            push!(bus_gen_circ_dict_ON[to][:circ], i)

            # Add each bus as adjacent to the other
            push!(bus_gen_circ_dict_ON[from][:adj_buses], to)
            push!(bus_gen_circ_dict_ON[to][:adj_buses], from)
        end

    end
    
    return bus_gen_circ_dict, bus_gen_circ_dict_ON
end

