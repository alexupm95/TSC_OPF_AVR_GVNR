# ===========================================
# Function to calculate the admittance matrix
# ===========================================
function Calculate_Ybus(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64, base_MVA::Float64)
    # Use Coordinate vectors (I, J, V) to build the sparse matrix directly
    I = Int[]
    J = Int[]
    V = ComplexF64[]
    
    # Pre-allocate memory to make push! operations
    sizehint!(I, 4 * nCIR + nBUS)
    sizehint!(J, 4 * nCIR + nBUS)
    sizehint!(V, 4 * nCIR + nBUS)

    # Calculate the admittance matrix including the line data
    for i = 1:nCIR
        if DCIR.l_status[i] == 1
            k = DCIR.from_bus[i] # Index from bus
            m = DCIR.to_bus[i]   # Index to bus
            
            ykm = 1.0 / (DCIR.l_res[i] + 1im * DCIR.l_reac[i]) # Series admittance
            bkm_sh = DCIR.l_sh_susp[i] / 2.0                 # Shunt admittance
            t_shift = deg2rad(DCIR.t_shift[i])
            tap = DCIR.t_tap[i]

            # k, k
            push!(I, k)
            push!(J, k)
            push!(V, ((1.0 / tap)^2 * ykm) + (1im * bkm_sh))
            
            # k, m
            push!(I, k)
            push!(J, m)
            push!(V, -((1.0 / tap) * ykm * exp(1im * t_shift)))
            
            # m, k
            push!(I, m)
            push!(J, k)
            push!(V, -((1.0 / tap) * ykm * exp(-1im * t_shift)))
            
            # m, m
            push!(I, m)
            push!(J, m)
            push!(V, ykm + (1im * bkm_sh))
        end
    end

    # Include the shunt components of the nodes
    for bus = 1:nBUS
        push!(I, bus)
        push!(J, bus)
        push!(V, (DBUS.g_sh[bus] + 1im * DBUS.b_sh[bus]) / base_MVA)
    end

    # sparse() automatically sums duplicate (i,j) entries (e.g., parallel lines or shunts)
    return sparse(I, J, V, nBUS, nBUS)
end

# =========================================================================
# Function to calculate the admittance matrix for the "during fault" period
# =========================================================================
function Calculate_Ybus_fault(Ybus_pref::SparseMatrixCSC, bus_fault::Int64)
    Ybus_fault = copy(Ybus_pref)

    Ybus_fault[bus_fault, bus_fault] += 1e9

    return Ybus_fault
end

# =======================================================================
# Function to calculate the admittance matrix for the "post-fault" period
# =======================================================================
function Calculate_Ybus_postf(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64, base_MVA::Float64, l_status::Vector{Int64})
    
    I = Int[]
    J = Int[]
    V = ComplexF64[]
    
    sizehint!(I, 4 * nCIR + nBUS)
    sizehint!(J, 4 * nCIR + nBUS)
    sizehint!(V, 4 * nCIR + nBUS)

    # Calculate the admittance matrix including the line data
    for i = 1:nCIR
        if l_status[i] == 1  #reflects the tripped line
            k = DCIR.from_bus[i] 
            m = DCIR.to_bus[i]   
            
            ykm = 1.0 / (DCIR.l_res[i] + 1im * DCIR.l_reac[i]) 
            bkm_sh = DCIR.l_sh_susp[i] / 2.0                 
            t_shift = deg2rad(DCIR.t_shift[i])
            tap = DCIR.t_tap[i]

            push!(I, k); push!(J, k); push!(V, ((1.0 / tap)^2 * ykm) + (1im * bkm_sh))
            push!(I, k); push!(J, m); push!(V, -((1.0 / tap) * ykm * exp(1im * t_shift)))
            push!(I, m); push!(J, k); push!(V, -((1.0 / tap) * ykm * exp(-1im * t_shift)))
            push!(I, m); push!(J, m); push!(V, ykm + (1im * bkm_sh))
        end
    end

    # Include the shunt components of the nodes
    for bus = 1:nBUS
        push!(I, bus)
        push!(J, bus)
        push!(V, (DBUS.g_sh[bus] + 1im * DBUS.b_sh[bus]) / base_MVA)
    end

    return sparse(I, J, V, nBUS, nBUS)
end