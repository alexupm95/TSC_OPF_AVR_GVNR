# Function to build the AC OPF model
function Make_ACOPF_Model!(model::Model, 
    Ybus::SparseMatrixCSC,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    bus_gen_circ_dict::OrderedDict,
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64)

    """
    This function builds a general AC OPF Model with the objective function, variables and constraints
    """
    #---------------------------------------------
    # Get some variables associated with the buses
    #---------------------------------------------
    P_d  = deepcopy(DBUS.p_d) ./ base_MVA  # Active power demanded by the loads     [p.u.]
    Q_d  = deepcopy(DBUS.q_d) ./ base_MVA  # Reactive power demanded by the loads   [p.u.]

    #------------------------------------------
    # Check if there is at least one swing bus
    #------------------------------------------
    SW = findall(x -> x == 3, DBUS.type)
    if isempty(SW) 
        throw(ArgumentError("You must define one bus as the SLACK BUS (type 3).")) 
    elseif length(SW) > 1 
        throw(ArgumentError("This code still does not support more than one SLACK BUS (type 3).")) 
    end

    #--------------------------------------------------------------------------
    # Defining dictionaries to save everything related to the optmization model
    #--------------------------------------------------------------------------
    model.ext[:objective] = OrderedDict{Symbol, Any}() # Objective

    # ========================================================================
    #              DEFINE THE VARIABLES OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                                Buses
    #-------------------------------------------------------------------------

    V = OrderedDict{Int, JuMP.VariableRef}()
    θ = OrderedDict{Int, JuMP.VariableRef}()

    for i in eachindex(DBUS.bus)
        V[i] = JuMP.@variable(model,              
        lower_bound = DBUS.v_min[i], # Lower Bounds
        upper_bound = DBUS.v_max[i], # Upper Bounds
        base_name = "V[$i]",          # Variable name -> Voltage magnitude of buses [p.u.]
        start = 1.0
        )

        θ[i] = JuMP.@variable(model,
        lower_bound = -π, # Lower Bounds
        upper_bound = π,  # Upper Bounds
        base_name = "θ[$i]",          # Variable name -> Voltage magnitude of buses [p.u.]
        start = 0.0
        )

    end

    #-------------------------------------------------------------------------
    #                            Generators
    #-------------------------------------------------------------------------

    P_g = OrderedDict{Int, JuMP.VariableRef}()
    Q_g = OrderedDict{Int, JuMP.VariableRef}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            P_g[i] = @variable(model,
                lower_bound = DGEN.pg_min[i] / base_MVA, # Lower Bounds
                upper_bound = DGEN.pg_max[i] / base_MVA, # Upper Bounds
                base_name = "P_g[$i]",
                start = max(0.0, DGEN.pg_min[i] / base_MVA)

            )

            Q_g[i] = @variable(model,
                lower_bound = DGEN.qg_min[i] / base_MVA, # Lower Bounds
                upper_bound = DGEN.qg_max[i] / base_MVA, # Upper Bounds
                base_name = "Q_g[$i]",
                start = max(0.0, DGEN.qg_min[i] / base_MVA)

            )
        end
    end
  
    #-------------------------------------------------------------------------
    #                                Branches
    #-------------------------------------------------------------------------

    # Check if the circuits have capacity/thermal limits defined in the input data
    all_cap           = [DCIR.l_cap_1 DCIR.l_cap_2 DCIR.l_cap_3] # Matrix containing all three capacity limits defined in the input data file
    lower_bounds_circ = OrderedDict{Int, Float64}()                     # Dictionary to save lower bounds of capacity limits for circuits if it is defined in the input data file
    upper_bounds_circ = OrderedDict{Int, Float64}()                     # Dictionary to save upper bounds of capacity limits for circuits if it is defined in the input data file

    # Loop to check if the branch has at least one capacity limit defined in the input data file
    # If not defined, the power flow through the branch has no bounds
    for i in 1:nCIR
        if any(!iszero, all_cap[i, :])  # Line ON or OFF and has capacity
            index_cap = findfirst(!iszero, all_cap[i, :])
            cap = all_cap[i, index_cap]
            lower_bounds_circ[i] = -cap / base_MVA # Lower Bounds
            upper_bounds_circ[i] =  cap / base_MVA # Upper Bounds
        end
    end
    
    # ========================================================================
    #            DEFINE THE OBJECTIVE FUNCTION OF THE MODEL
    # ========================================================================
    # Minimize total fuel cost

    total_cost = 0.0  # Start from zero

    for gen in 1:nGEN
        if DGEN.g_status[gen] == 1 && haskey(P_g, gen)
            if DGEN.g_cost_0[gen] != 0.0
                total_cost += DGEN.g_cost_0[gen]
            end
            if DGEN.g_cost_1[gen] != 0.0
                total_cost += DGEN.g_cost_1[gen] * P_g[gen] * base_MVA
            end
            if DGEN.g_cost_2[gen] != 0.0
                total_cost += DGEN.g_cost_2[gen] * (P_g[gen] * base_MVA)^2
            end
        end
    end

    model.ext[:objective] = @objective(model, Min, total_cost)
    
    # ========================================================================
    #            DEFINE THE CONSTRAINTS OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                  Equality and Inequality Constraints 
    #-------------------------------------------------------------------------

    # *********************************
    # Constraint for angle of swing bus
    # *********************************
    eq_const_angle_sw = OrderedDict{Int, JuMP.ConstraintRef}()

    if any(bus_gen_circ_dict[SW[1]][:gen_status] .== 1)                          # First check if the swing bus has at least one generator connected to it
        eq_const_angle_sw = JuMP.@constraint(model, θ[DBUS.bus[SW[1]]] == 0.0)   # Set the constraint -> Angle == 0
    else
        throw(ArgumentError("Swing bus $(SW[1]) must have at least one connected generator with status ON.")) # Throw an error if the swing bus has no generator connected to it
    end

    # ******************************************************************************************************
    # Create expressions for the net inflow/outflow of power in all buses according to the admittance matrix
    # ******************************************************************************************************
    p_net_dict = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of net active power at the bus
    q_net_dict = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of net active power at the bus

    for i in 1:nBUS # Loop in all buses

        indices_circ_connected = bus_gen_circ_dict[i][:circ]    # Circuits connected to bus i
        if isempty(indices_circ_connected)
            throw(ArgumentError("The bus $i is islanded, i.e., there is no line or transformer connected to it."))
        end

        # terms_p_net = JuMP.NonlinearExpr[] # List of expressions
        # terms_q_net = JuMP.NonlinearExpr[] # List of expressions
    #     for k in 1:nBUS
    #         angik   = θ[i] - θ[k]         # Angular difference between bus i and k
    #         Gik = real(Ybus[i,k])
    #         Bik = imag(Ybus[i,k])

    #         if !iszero(Gik) && !iszero(Bik)
    #             push!(terms_p_net, V[k] * (Gik * cos(angik) + Bik * sin(angik))) 
    #             push!(terms_q_net, V[k] * (Gik * sin(angik) - Bik * cos(angik)))

    #         elseif iszero(Gik) && !iszero(Bik)
    #             push!(terms_p_net, V[k] * (Bik * sin(angik)))   
    #             push!(terms_q_net, V[k] * (- Bik * cos(angik)))

    #         elseif !iszero(Gik) && iszero(Bik)
    #             push!(terms_p_net, V[k] * (Gik * cos(angik)))
    #             push!(terms_q_net, V[k] * (Gik * sin(angik)))
    #         end
    #     end
    #     p_net_dict[i] = V[i] * sum(terms_p_net)
    #     q_net_dict[i] = V[i] * sum(terms_q_net)
    # end
        # Get the non-zero elements for row 'i' (connections to bus 'i')
        row_Y = Ybus[i, :] 
        connected_buses = row_Y.nzind # Indices of connected buses (the 'k's)
        y_ik_vals = row_Y.nzval       # The complex admittance values

        p_net_dict[i] = @expression(model, 
            V[i] * sum(
                V[k] * (real(y_ik) * cos(θ[i] - θ[k]) + imag(y_ik) * sin(θ[i] - θ[k]))
                for (k, y_ik) in zip(connected_buses, y_ik_vals)
            )
        )

        q_net_dict[i] = @expression(model, 
            V[i] * sum(
                V[k] * (real(y_ik) * sin(θ[i] - θ[k]) - imag(y_ik) * cos(θ[i] - θ[k]))
                for (k, y_ik) in zip(connected_buses, y_ik_vals)
            )
        )
    end

   
    # **************************************************
    # Constraints for Active and Reactive Power Balance
    # *************************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_q_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in 1:nBUS
        indices_bus_gen        = bus_gen_circ_dict[i][:gen_ids] # Generators at bus i
        indices_circ_connected = bus_gen_circ_dict[i][:circ]    # Circuits connected to bus i

        if isempty(indices_circ_connected) # Check if the bus has at least one branch connected to it
            throw(ArgumentError("The bus $i is islanded, i.e., there is no line or transformer connected to it."))
        end

        # Get generator variables from the dictionary for those ON at this bus
        Pg_terms = [P_g[g] for g in indices_bus_gen if haskey(P_g, g)]
        Qg_terms = [Q_g[g] for g in indices_bus_gen if haskey(Q_g, g)]

        Pg_sum = isempty(Pg_terms) ? 0.0 : sum(Pg_terms) # Sum P_g terms; zero otherwise
        Qg_sum = isempty(Qg_terms) ? 0.0 : sum(Qg_terms) # Sum Q_g terms; zero otherwise

        eq_const_p_balance[i] = @constraint(model, p_net_dict[i] == Pg_sum - P_d[i])
        eq_const_q_balance[i] = @constraint(model, q_net_dict[i] == Qg_sum - Q_d[i])

    end

    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    # Loop in all branches
    for lin = 1:nCIR
        if DCIR.l_status[lin] == 1 # Check if the branch is ON

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

            angik   = θ[i] - θ[k]                                  # Angular difference between bus i and k
            angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i
            
            inv_tap2    = 1.0 / (t_tap^2)
            c_Pik_cos   = (-g * t_r + b * t_i) * inv_tap2
            c_Pik_sin   = (-b * t_r - g * t_i) * inv_tap2
            c_Qik_cos   = -(-b * t_r - g * t_i) * inv_tap2
            c_Qik_sin   = (-g * t_r + b * t_i) * inv_tap2
            
            c_Pki_cos   = (-g * t_r - b * t_i) * inv_tap2
            c_Pki_sin   = (-b * t_r + g * t_i) * inv_tap2
            c_Qki_cos   = -(-b * t_r + g * t_i) * inv_tap2
            c_Qki_sin   = (-g * t_r - b * t_i) * inv_tap2

            # ***************************************************
            # Equality Constraints for power flow in the branches
            # ***************************************************

            # Line flow from i to k
            # Same equations used in POWERMODELS
            P_ik_expr = g*((1/t_tap) * V[i])^2 + c_Pik_cos*(V[i]*V[k]*cos(angik)) + c_Pik_sin*(V[i]*V[k]*sin(angik))
            Q_ik_expr = -(b + bik_sh)*((1/t_tap) * V[i])^2 + c_Qik_cos*(V[i]*V[k]*cos(angik)) + c_Qik_sin*(V[i]*V[k]*sin(angik))


            # Line flow from k to i
            # Same equations used in POWERMODELS
            P_ki_expr = g*(V[k]^2) + c_Pki_cos*(V[i] * V[k] * cos(angki)) + c_Pki_sin*( V[i] * V[k] * sin(angki))
            Q_ki_expr = -(b + bik_sh)*(V[k]^2)  + c_Qki_cos*( V[i] * V[k] * cos(angki)) + c_Qki_sin*( V[i] * V[k] * sin(angki))

            # ****************************************************************
            # Inequality Constraint for capacity/thermal limits of power flows
            # ****************************************************************
            ineq_const_s_ik[lin] = JuMP.@constraint(model, P_ik_expr^2 + Q_ik_expr^2 <= (get(upper_bounds_circ, lin, Inf))^2)
            ineq_const_s_ki[lin] = JuMP.@constraint(model, P_ki_expr^2 + Q_ki_expr^2 <= (get(upper_bounds_circ, lin, Inf))^2)
         
        end
    end

    # Initialize dictionary with count, min_ang, max_ang
    pair_info = Dict{Tuple{Int, Int}, NamedTuple{(:count, :min_ang, :max_ang), Tuple{Int, Float64, Float64}}}()
    pair_circ_map = Dict{Tuple{Int, Int}, Int}()

    # This loop maps the circuits that have parallel branches, count the number of parallel lines/transformers,
    # and save the minimum and maximum angular difference for the buses related to these circuits
    for i in 1:nCIR
        if DCIR.l_status[i] == 1                     # Only consider active branches
            a, b = DCIR.from_bus[i], DCIR.to_bus[i]  # Bus from and bus to
            pair = (min(a, b), max(a, b))            # Sort the buses id
            circ_num = DCIR.circ[i]                  # Get the number of the branch

            if !haskey(pair_circ_map, pair)          # Check if these buses were already addded
                pair_circ_map[pair] = circ_num       # Get only the number of the first branch ON connecting these buses
            else
                pair_circ_map[pair] = min(pair_circ_map[pair], circ_num) # Get the number of the first branch ON connecting these buses
            end
            
            if haskey(pair_info, pair) # Check if this pair of buses were already added
                # Update count
                old       = pair_info[pair]
                new_count = old.count + 1                       # Count the number of branches ON connecting these buses
                new_min   = min(old.min_ang, DCIR.ang_min[i])   # Get the minimum angle defined for the branches connecting these buses
                new_max   = max(old.max_ang, DCIR.ang_max[i])   # Get the maximum angle defined for the branches connecting these buses
                pair_info[pair] = (new_count, new_min, new_max)
            else
                pair_info[pair] = (1, DCIR.ang_min[i], DCIR.ang_max[i]) # Add this pair of buses for the first time
            end
        end
    end
    sorted_pair_info  = sort(collect(pair_info), by = x -> pair_circ_map[x[1]]) # Sort the data inside pair_info Dict

    ineq_const_diff_ang  = OrderedDict{Int, JuMP.ConstraintRef}() # Vector to save inequality constraints of angle difference between adjacent buses

    for (pair, pair_data) in sorted_pair_info
        bus_from, bus_to = pair
        ang_ik = θ[bus_from] - θ[bus_to]

        circ_id = pair_circ_map[pair]

        min_ang = pair_data.min_ang
        max_ang = pair_data.max_ang

        # ***************************************************
        # Inequality Constraint for voltage angle differences
        # ***************************************************
        if min_ang >= -60 && max_ang <= 60 
            ineq_const_diff_ang[circ_id] = JuMP.@constraint(model, deg2rad(min_ang) <= ang_ik <= deg2rad(max_ang))

        elseif min_ang < -60 && max_ang <= 60
            println("Correcting angle constraints between adjacent buses ($bus_from, $bus_to): setting ang_min to -60°.")
            ineq_const_diff_ang[circ_id] = JuMP.@constraint(model, deg2rad(-60) <= ang_ik <= deg2rad(max_ang))

        elseif min_ang >= -60 && max_ang > 60
            println("Correcting angle constraints between adjacent buses ($bus_from, $bus_to): setting ang_max to +60°.")
            ineq_const_diff_ang[circ_id] = JuMP.@constraint(model, deg2rad(min_ang) <= ang_ik <= deg2rad(60))

        else # min_ang < -60 && max_ang > 60
            println("Correcting angle constraints between adjacent buses ($bus_from, $bus_to): setting ang_min to -60° and ang_max to +60°.")
            ineq_const_diff_ang[circ_id] = JuMP.@constraint(model, deg2rad(-60) <= ang_ik <= deg2rad(60))
        end
    end

    return model, V, θ, P_g, Q_g, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, ineq_const_s_ik, ineq_const_s_ki, ineq_const_diff_ang # Return the model to the main function
end
