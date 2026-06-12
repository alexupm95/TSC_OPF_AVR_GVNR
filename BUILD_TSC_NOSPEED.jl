# Functions created to Build the Transient Stability Constraints into the Model

# ===================================================================================
#                          DEFINE INITIAL VARIABLES PRE-FAULT
# ===================================================================================
# Function used to define E, δ and Pm variables of each active generator
function Define_Dyn_Var_EδPm!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef},
    val_V::Dict, val_θ::Dict, val_Pg::Dict, val_Qg::Dict
    )

    # ======================================================================================
    # Setting the variables of initial internal voltage magnitude and angle of the generator
    # ======================================================================================

    E_fd = OrderedDict{Int, JuMP.VariableRef}()  # Initialize the dictionary to save the internal voltage magnitudes
    δ = OrderedDict{Int, JuMP.VariableRef}()  # Initialize the dictionary to save the rotor angles
    δCOI_aux = []
    Pm = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the mechanical power

    Ed = OrderedDict{Int, JuMP.VariableRef}()
    Eq = OrderedDict{Int, JuMP.VariableRef}()
    Id   = OrderedDict{Int, JuMP.VariableRef}()
    Iq   = OrderedDict{Int, JuMP.VariableRef}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]       
            Xd = DGEN_DYN.Xd[i] 
            Xq = DGEN_DYN.Xq[i]
            Xd_tr = DGEN_DYN.Xd_tr[i] 
            Xq_tr = DGEN_DYN.Xq_tr[i] 
            Ra = DGEN_DYN.Ra[i] 

            # --- CALCULATE EXACT WARM START W/ SOLVED ACOPF VALUES ---
            v_val  = val_V[bus]
            th_val = val_θ[bus]
            pg_val = val_Pg[i]
            qg_val = val_Qg[i]

            # Convert to complex domain
            V_c = v_val * exp(1im * th_val)
            I_c = (pg_val - 1im * qg_val) / conj(V_c)

            Eq_c = V_c + I_c * (Ra + 1im * Xq) # Calculate exact internal voltage E'
            start_δ = angle(Eq_c)

            # Project I and V onto d-q axes
            I_mag = abs(I_c)
            I_ang = angle(I_c)
            start_Id = I_mag * sin(start_δ - I_ang)
            start_Iq = I_mag * cos(start_δ - I_ang)
            
            start_Vd = v_val * sin(start_δ - th_val)
            start_Vq = v_val * cos(start_δ - th_val)
            
            start_Ed = start_Vd + Ra * start_Id - Xq_tr * start_Iq
            start_Eq = start_Vq + Ra * start_Iq + Xd_tr * start_Id
            start_Efd  = start_Eq + (Xd - Xd_tr) * start_Id

            E_fd[i] = JuMP.@variable(model,  # Define the variable related to the internal voltage magnitude
                lower_bound = 0.0,   # Set Lower Bounds
                upper_bound = 4.0,   # Set Upper Bounds
                base_name = "E_fd[$i]", # Set a name
                start = start_Efd          # Set a flat start
            )
            δ[i] = JuMP.@variable(model,  # Define the variable related to the initial internal angle
                lower_bound = -π,    # Set Lower Bounds
                upper_bound = π,     # Set Upper Bounds
                base_name = "δ[$i]", # Set a name
                start = start_δ          # Set a flat start
            )

            Pm[i] = JuMP.@variable(model,                     # Define the variable related to the initial mechanical power
                lower_bound = DGEN.pg_min[i] / base_MVA, # Set Lower Bounds
                upper_bound = DGEN.pg_max[i] / base_MVA, # Set Upper Bounds
                base_name = "Pm[$i]",                    # Set a name
                start = pg_val        # Set a flat start
            )

            Ed[i] = JuMP.@variable(model,
                lower_bound = -Inf, upper_bound = Inf,
                base_name   = "Ed_p[$i]",
                start       = start_Ed
            )
            Eq[i] = JuMP.@variable(model,
                lower_bound = -Inf,  upper_bound = Inf,
                base_name   = "Eq_p[$i]",
                start       = start_Eq
            )
            Id[i] = JuMP.@variable(model,
                lower_bound = -Inf,
                upper_bound =  Inf,
                base_name   = "Id[$i]",
                start       = start_Id
            )
            Iq[i] = JuMP.@variable(model,
                lower_bound = -Inf,
                upper_bound =  Inf,
                base_name   = "Iq[$i]",
                start       = start_Iq
            )

            push!(δCOI_aux, δ[i]*DGEN_DYN.H[i])
        end
    end

    # Initialize the dictionary to save the rotor angle of the COI
    δCOI = JuMP.@variable(model, base_name = "δCOI")
    JuMP.@constraint(model, δCOI == sum(δCOI_aux) / sum(DGEN_DYN.H .* DGEN.g_status))

    # ======================================================================================
    #           Setting the constraints of initial active and reactive
    #               power to define the initial value of E and δ
    # ======================================================================================

    eq_const_Ed_init = OrderedDict{Int, JuMP.ConstraintRef}() 
    eq_const_Eq_init = OrderedDict{Int, JuMP.ConstraintRef}() 
    eq_const_Vd_init = OrderedDict{Int, JuMP.ConstraintRef}() 
    eq_const_Vq_init = OrderedDict{Int, JuMP.ConstraintRef}() 
    eq_const_P_init = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the constraints of active power generated in the pre-fault stage
    eq_const_Q_init = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the constraints of reactive power generated in the pre-fault stage
    eq_const_Pm     = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the equality constraint of mechanical power in the pre-fault stage

    # Loop to create the equality constraints
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]         # ID of the bus in which the generator is connected

            bus = DGEN.bus[i]       
            Xd = DGEN_DYN.Xd[i] 
            Xq = DGEN_DYN.Xq[i]
            Xd_tr = DGEN_DYN.Xd_tr[i] 
            Xq_tr = DGEN_DYN.Xq_tr[i] 
            Ra = DGEN_DYN.Ra[i]

            # Pre-fault steady state constraints (derivatives = 0)
            eq_const_Ed_init[i]=JuMP.@constraint(model, Ed[i] == (Xq - Xq_tr)*Iq[i])
            eq_const_Eq_init[i]=JuMP.@constraint(model, Eq[i] + (Xd - Xd_tr)*Id[i] == E_fd[i])
            
            eq_const_Vd_init[i]=JuMP.@constraint(model, V[bus] * sin(δ[i] - θ[bus]) - Ed[i] + Ra * Id[i] - Xq_tr * Iq[i] == 0)
            eq_const_Vq_init[i]=JuMP.@constraint(model, V[bus] * cos(δ[i] - θ[bus]) - Eq[i] + Ra * Iq[i] + Xd_tr * Id[i] == 0)
            
            # Match grid steady-state injection
            eq_const_P_init[i] = JuMP.@constraint(model, P_g[i] == V[bus] * sin(δ[i] - θ[bus]) * Id[i] + V[bus] * cos(δ[i] - θ[bus]) * Iq[i])
            eq_const_Q_init[i] = JuMP.@constraint(model, Q_g[i] == V[bus] * cos(δ[i] - θ[bus]) * Id[i] - V[bus] * sin(δ[i] - θ[bus]) * Iq[i])

            # Constraint for mechanical power
            eq_const_Pm[i] = JuMP.@constraint(model,
            Pm[i] == P_g[i]
            )
        end
    end

    return model, E_fd, δ, Pm, δCOI, Ed, Eq, Id, Iq, eq_const_P_init, eq_const_Q_init, eq_const_Ed_init, eq_const_Eq_init, eq_const_Vd_init, eq_const_Vq_init,  eq_const_Pm
end

# ===================================================================================
#                                  FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the fault period
function Def_Dyn_Fault_All!(
    model::Model,
    DBUS::DataFrame,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nBUS::Int64,
    nGEN::Int64,
    base_MVA::Float64,
    ZIP::Vector{Float64},
    t_window_fault::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64},
    Ybus_fault::SparseMatrixCSC,
    V::OrderedDict{Int64, VariableRef},
    θ::OrderedDict{Int64, VariableRef},
    E_fd::OrderedDict{Int64, VariableRef},
    V_ref::OrderedDict{Int64, VariableRef},
    P_ref::OrderedDict, P_valve::OrderedDict,
    Pm::OrderedDict{Int64, VariableRef},
    δ::OrderedDict{Int64, VariableRef},
    Δω_0::Float64,
    P_g::OrderedDict{Int64, VariableRef},
    ω_syn::Float64,
    Δt::Float64,
    val_V::Dict, val_θ::Dict, val_Pg::Dict, val_Qg::Dict,
    Ed::OrderedDict, Eq::OrderedDict, Id::OrderedDict, Iq::OrderedDict, 
	)

    # ======================================================================================
    #               Setting the variables used in the fault period
    # ======================================================================================

    Pe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables electrical active power for each period of time (fault period)
    Qe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables electrical reactive power for each period of time (fault period)
    δ_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables internal angle for each period of time (fault period)
    Δω_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables speed deviation for each period of time (fault period)
    V_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables bus voltage magnitude for each period of time (fault period)
    θ_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables bus voltage angle for each period of time (fault period)
    Ed_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Eq_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Id_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Iq_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()

    # Terms that will be summed to calculate the dynamics of the COI
    terms_δCOI_tf_dict  = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary
    terms_ΔωCOI_tf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary

    # Loop to create the variables
    for i in 1:nBUS
        # Initialize inner dict for this bus
        V_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
        θ_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
        
        v_start = val_V[i]
        th_start = val_θ[i]

        for t in eachindex(t_window_fault)
            V_tf[i][t] = JuMP.@variable(model, # Define the variable related to the voltage magnitude at the buses during fault
            lower_bound = 0.0,
            base_name = "V_tf[$i, $t]",         # Set the name
            start = v_start
            )

            θ_tf[i][t] = JuMP.@variable(model, # Define the variable related to the voltage angle at the buses during fault
            base_name = "θ_tf[$i, $t]",         # Set the name
            start = th_start
            )
        end
    end

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            Pe_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Qe_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            δ_tf[i]      = OrderedDict{Int, JuMP.VariableRef}()
            Δω_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Ed_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Eq_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Id_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Iq_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()

            terms_δCOI_tf_dict[i]  = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()
            terms_ΔωCOI_tf_dict[i] = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()

            pg_start = val_Pg[i]
            qg_start = val_Qg[i]
            delta_start = JuMP.start_value(δ[i]) 
            Ed_start=JuMP.start_value(Ed[i])
            Eq_start=JuMP.start_value(Eq[i]) 
            Id_start=JuMP.start_value(Id[i]) 
            Iq_start=JuMP.start_value(Iq[i]) 

            for t in eachindex(t_window_fault)
                # Define the variable related to the electrical active power during fault
                Pe_tf[i][t] = JuMP.@variable(model, 
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "Pe_tf[$i, $t]", start=pg_start)

                # Define the variable related to the electrical reactive power during fault
                Qe_tf[i][t] = JuMP.@variable(model, 
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "Qe_tf[$i, $t]", start=qg_start)

                # Define the variable related to the internal angle during fault
                δ_tf[i][t] = JuMP.@variable(model, 
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "δ_tf[$i, $t]", start=delta_start)

                # Define the variable related to the speed deviation during fault
                Δω_tf[i][t] = JuMP.@variable(model, 
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "Δω_tf[$i, $t]", start=0.0)

                Ed_tf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Ed_tf[$i, $t]", start=Ed_start)
                Eq_tf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Eq_tf[$i, $t]", start=Eq_start)
                Id_tf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Id_tf[$i, $t]", start=Id_start)
                Iq_tf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Iq_tf[$i, $t]", start=Iq_start)

                terms_δCOI_tf_dict[i][t]  = DGEN_DYN.H[i] * δ_tf[i][t]
                terms_ΔωCOI_tf_dict[i][t] = DGEN_DYN.H[i] * Δω_tf[i][t]

            end
        end
    end


    # Gives the expressions to calculate the angle and speed deviation of the COI
    expr_δCOI_per_time  = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_δCOI_tf_dict))  for t in eachindex(t_window_fault))
    expr_ΔωCOI_per_time = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_ΔωCOI_tf_dict)) for t in eachindex(t_window_fault))

    #Exciter
    model, E_fd_tf, eq_const_E_fd_tf = Def_Dyn_Fault_Exciter!(model, DGEN, DGEN_DYN, nGEN, E_fd, V_ref, V_tf, t_window_fault, Δt)

    #Governor
    model, P_valve_tf, P_mech_tf, eq_const_P_valve_tf, eq_const_P_mech_tf = Def_Dyn_Fault_Governor!(model, DGEN, DGEN_DYN, nGEN, P_ref, P_valve, Pm, Δω_tf, t_window_fault, Δt, base_MVA)

    # Define the variables and equality constraints of the COI
    model, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf = Def_Dyn_Fault_COI!(model, DGEN, DGEN_DYN, nGEN, t_window_fault, expr_δCOI_per_time, expr_ΔωCOI_per_time)

    # Define the inequality constraints of angle and speed deviation related to their limits
    model, ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf = Def_Dyn_Fault_rel_COI!(model, DGEN, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, nGEN, t_window_fault, δ_tol, Δω_tol)

    # Define the equality constraints for Pe in the fault period
    model, eq_const_Pe_tf, eq_const_Qe_tf, eq_const_Pbalance_tf, eq_const_Qbalance_tf, eq_const_Vd_tf, eq_const_Vq_tf, Te_tf, eq_const_Te_tf = Def_Dyn_Fault_eqconst_PeQe!(model, DBUS, nBUS, bus_gen_circ_dict, DGEN, DGEN_DYN, nGEN, base_MVA, V, θ, E_fd, δ, Pe_tf, Qe_tf, δ_tf, V_tf, θ_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf,Δω_tf, t_window_fault, Ybus_fault, ZIP)

    # Define the equality constraints of the swing equation for δ and Δω during the fault period
    model, eq_const_δ_tf, eq_const_Δω_tf = Def_Dyn_Fault_Swing!(model, DGEN, DGEN_DYN, nGEN, P_mech_tf, Pe_tf, δ_tf, Δω_tf, δ, Δω_0, P_g, t_window_fault, ω_syn, Δt)
    
    model, eq_const_Ed_tf, eq_const_Eq_tf = Def_Dyn_Fault_EMF!(model, DGEN, DGEN_DYN, nGEN, E_fd_tf, Ed, Eq, Id, Iq, Ed_tf, Eq_tf, Id_tf, Iq_tf, t_window_fault, Δt)

    return model, Pe_tf, Qe_tf, V_tf, θ_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf, E_fd_tf, P_valve_tf, P_mech_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf, eq_const_Pe_tf, eq_const_Qe_tf, eq_const_Pbalance_tf, eq_const_Qbalance_tf, eq_const_Vd_tf, eq_const_Vq_tf, eq_const_δ_tf, eq_const_Δω_tf, eq_const_Ed_tf, eq_const_Eq_tf, eq_const_E_fd_tf, eq_const_P_valve_tf, eq_const_P_mech_tf
end

# Function that define the variables of the COI as well as its equality constraints across the fault period
function Def_Dyn_Fault_COI!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_fault::Vector{Float64},
    aux_expr_δCOI_per_time::OrderedDict{Int64, AffExpr},
    aux_expr_ΔωCOI_per_time::OrderedDict{Int64, AffExpr}
    )

    H_COI = sum(DGEN_DYN.H[i] for i in 1:nGEN if DGEN.g_status[i] == 1) # Inertia of the COI

    δCOI_tf  = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables angle of the COI for each period of time (fault period)
    ΔωCOI_tf = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables speed deviation of the COI for each period of time (fault period)

    eq_const_δCOI_tf  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_ΔωCOI_tf = OrderedDict{Int, JuMP.ConstraintRef}()

    for t in eachindex(t_window_fault)
        δCOI_tf[t]  = JuMP.@variable(model, base_name = "δCOI_tf[$t]")
        ΔωCOI_tf[t] = JuMP.@variable(model, base_name = "ΔωCOI_tf[$t]")
        
        eq_const_δCOI_tf[t]  = JuMP.@constraint(model, δCOI_tf[t]  == aux_expr_δCOI_per_time[t]  / H_COI)
        eq_const_ΔωCOI_tf[t] = JuMP.@constraint(model, ΔωCOI_tf[t] == aux_expr_ΔωCOI_per_time[t] / H_COI)

    end

    return model, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf

end

# Function that define the constraints of the angle of generators in relation to the COI variables across the fault period
function Def_Dyn_Fault_rel_COI!(model::Model,
    DGEN::DataFrame,
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δCOI_tf::OrderedDict{Int, JuMP.VariableRef},
    ΔωCOI_tf::OrderedDict{Int, JuMP.VariableRef},
    nGEN::Int64,
    t_window_fault::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64}
    )

    ineq_const_δ_COI_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    ineq_const_Δω_COI_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            ineq_const_δ_COI_tf[i]  = OrderedDict{Int, JuMP.ConstraintRef}()
            ineq_const_Δω_COI_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_fault)
                ineq_const_δ_COI_tf[i][t] = JuMP.@constraint(model,
                δ_tol[1] <= δ_tf[i][t] - δCOI_tf[t] <= δ_tol[2]
                )

                ineq_const_Δω_COI_tf[i][t] = JuMP.@constraint(model,
                Δω_tol[1] <= Δω_tf[i][t] - ΔωCOI_tf[t] <= Δω_tol[2]
                )
            end
        end
    end

    return model, ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf


end

# Function that create the constraints of the elecrical power across the fault period
function Def_Dyn_Fault_eqconst_PeQe!(model::Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    base_MVA::Float64,
    V::OrderedDict{Int64, VariableRef},
    θ::OrderedDict{Int64, VariableRef},
    E_fd::OrderedDict{Int64, VariableRef},
    δ ::OrderedDict{Int64, VariableRef},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Qe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    V_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    θ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Ed_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Eq_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Id_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Iq_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_fault::Vector{Float64},
    Ybus_fault::SparseMatrixCSC,
    ZIP::Vector{Float64}
    )
    
    eq_const_Pe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Qe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Vd_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Vq_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    eq_const_Pbalance_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Qbalance_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    terms_Pb_tf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() # Dictionary of Terms that will be summed to calculate Pb
    terms_Qb_tf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() # Dictionary of Terms that will be summed to calculate Qb

    Te_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq_const_Te_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    active_gen = findall(isone, DGEN.g_status) # Vector with the indices of the active generators

    # Generate the expression of power generated by each generator
    for g in active_gen
        bus_g = DGEN.bus[g]
        Xd_tr = DGEN_DYN.Xd_tr[g]; Xq_tr = DGEN_DYN.Xq_tr[g]; Ra = DGEN_DYN.Ra[g]

        eq_const_Pe_tf[g] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints
        eq_const_Qe_tf[g] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints
        eq_const_Vd_tf[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_const_Vq_tf[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        Te_tf[g]          = OrderedDict{Int, JuMP.VariableRef}()
        eq_const_Te_tf[g] = OrderedDict{Int, JuMP.ConstraintRef}()       

        for t in eachindex(t_window_fault) # Loop over the time
            Te_tf[g][t] = JuMP.@variable(model,
            lower_bound = -Inf, upper_bound = Inf,
            base_name = "Te_tf[$g,$t]",
            start = val_Pg[g] / base_MVA)

            eq_const_Te_tf[g][t] = JuMP.@constraint(model,
            Te_tf[g][t] == Ed_tf[g][t]*Id_tf[g][t] + Eq_tf[g][t]*Iq_tf[g][t])

            Vd_expr = @expression(model, V_tf[bus_g][t] * sin(δ_tf[g][t] - θ_tf[bus_g][t]))
            Vq_expr = @expression(model, V_tf[bus_g][t] * cos(δ_tf[g][t] - θ_tf[bus_g][t]))

            # Stator algebraic equations linking internal voltage to terminal voltage
            eq_const_Vd_tf[g][t] = JuMP.@constraint(model, Vd_expr - (1 + Δω_tf[g][t])*Ed_tf[g][t] + Ra * Id_tf[g][t] - Xq_tr * Iq_tf[g][t] == 0)
            eq_const_Vq_tf[g][t] = JuMP.@constraint(model, Vq_expr - (1 + Δω_tf[g][t])*Eq_tf[g][t] + Ra * Iq_tf[g][t] + Xd_tr * Id_tf[g][t] == 0)

            eq_const_Pe_tf[g][t] = JuMP.@constraint(model, Pe_tf[g][t] == (1 + Δω_tf[g][t])*Te_tf[g][t]) # Equality constraint for the active power generated
            #eq_const_Pe_tf[g][t] = JuMP.@constraint(model, Pe_tf[g][t] == Vq_expr*Iq_tf[g][t] + Vd_expr*Id_tf[g][t]) # Equality constraint for the active power generated
            eq_const_Qe_tf[g][t] = JuMP.@constraint(model, Qe_tf[g][t] == Vq_expr*Id_tf[g][t] - Vd_expr*Iq_tf[g][t]) # Equality constraint for the reactive power generated
        end
    end

    # Generate the expressions for the power balance during the fault
    for i in 1:nBUS

        indices_circ_connected = bus_gen_circ_dict[i][:circ]    # Circuits connected to bus i

        if isempty(indices_circ_connected) # Check if the bus has at least one branch connected to it
            throw(ArgumentError("The bus $i is islanded, i.e., there is no line or transformer connected to it."))
        end

        terms_Pb_tf_dict[i] = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of terms to form a nonlinear expression
        terms_Qb_tf_dict[i] = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of terms to form a nonlinear expression

        indices_bus_gen        = bus_gen_circ_dict[i][:gen_ids] # Generators at bus i

        eq_const_Pbalance_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints
        eq_const_Qbalance_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints


        # Load as constant admittance Y_d = (S_d)* / |V|² = (P_d - jQ_d) / |V|²
        p_d = DBUS.p_d[i] / base_MVA
        q_d = (1.0 .* DBUS.q_d[i]) / base_MVA # It must be multiplied by (-1) here to avoid confusion hereafter

        for t in eachindex(t_window_fault) # Loop over the time

            row_Y = Ybus_fault[i, :] # Grab only the non-zero connections for Bus i
            
            #avoid zero checks with row_Y.nzind
            terms_Pb_tf_dict[i][t] = @expression(model, 
                V_tf[i][t] * sum(
                    V_tf[j][t] * (real(y_ij) * cos(θ_tf[i][t] - θ_tf[j][t]) + imag(y_ij) * sin(θ_tf[i][t] - θ_tf[j][t]))
                    for (j, y_ij) in zip(row_Y.nzind, row_Y.nzval)
                )
            )

            terms_Qb_tf_dict[i][t] = @expression(model, 
                V_tf[i][t] * sum(
                    V_tf[j][t] * (real(y_ij) * sin(θ_tf[i][t] - θ_tf[j][t]) - imag(y_ij) * cos(θ_tf[i][t] - θ_tf[j][t]))
                    for (j, y_ij) in zip(row_Y.nzind, row_Y.nzval)
                )
            )

            Total_P_inj = isempty(indices_bus_gen) ? 0.0 : @expression(model, sum(V_tf[i][t] * Id_tf[g][t] * sin(δ_tf[g][t] - θ_tf[i][t]) + V_tf[i][t] * Iq_tf[g][t] * cos(δ_tf[g][t] - θ_tf[i][t]) for g in indices_bus_gen))
            Total_Q_inj = isempty(indices_bus_gen) ? 0.0 : @expression(model, sum(V_tf[i][t] * Id_tf[g][t] * cos(δ_tf[g][t] - θ_tf[i][t]) - V_tf[i][t] * Iq_tf[g][t] * sin(δ_tf[g][t] - θ_tf[i][t]) for g in indices_bus_gen))
                    # if p_d == 0.0
            #     eq_const_Pbalance_tf[i][t] = JuMP.@constraint(model, terms_Pb_tf_dict[i][t] == 0.0) # Equality constraint for balance of the active power
            # else
            #     eq_const_Pbalance_tf[i][t] = JuMP.@constraint(model, terms_Pb_tf_dict[i][t] + p_d*ZIP[1] + p_d*ZIP[2]*(V_tf[i][t] / V[i]) + p_d*ZIP[3]*((V_tf[i][t] / V[i])^2) == 0.0) # Equality constraint for balance of the active power

            # end
            # if q_d == 0.0
            #     eq_const_Qbalance_tf[i][t] = JuMP.@constraint(model, terms_Qb_tf_dict[i][t] == 0.0) # Equality constraint for balance of the reactive power
            # else
            #     eq_const_Qbalance_tf[i][t] = JuMP.@constraint(model, terms_Qb_tf_dict[i][t] + q_d*ZIP[1] + q_d*ZIP[2]*(V_tf[i][t] / V[i]) + q_d*ZIP[3]*((V_tf[i][t] / V[i])^2) == 0.0) # Equality constraint for balance of the reactive power
            # end

            if p_d == 0.0
                eq_const_Pbalance_tf[i][t] = JuMP.@constraint(model, terms_Pb_tf_dict[i][t] == Total_P_inj) # Equality constraint for balance of the active power
            else
                eq_const_Pbalance_tf[i][t] = JuMP.@constraint(model, terms_Pb_tf_dict[i][t]*V[i]^2 == Total_P_inj*V[i]^2 - p_d*(ZIP[1]*V[i]^2 + ZIP[2]*(V_tf[i][t] * V[i]) + ZIP[3]*(V_tf[i][t]^2)) ) # Equality constraint for balance of the active power
                #eq_const_Pbalance_tf[i][t] = JuMP.@constraint(model, terms_Pb_tf_dict[i][t] == Total_P_inj - p_d*(ZIP[1] + ZIP[2]*(V_tf[i][t]) + ZIP[3]*(V_tf[i][t]^2)) ) # Equality constraint for balance of the active power
            end

            if q_d == 0.0
                eq_const_Qbalance_tf[i][t] = JuMP.@constraint(model, terms_Qb_tf_dict[i][t] == Total_Q_inj) # Equality constraint for balance of the reactive power
            else
                eq_const_Qbalance_tf[i][t] = JuMP.@constraint(model, terms_Qb_tf_dict[i][t]*V[i]^2 == Total_Q_inj*V[i]^2 - q_d*(ZIP[1]*V[i]^2 + ZIP[2]*(V_tf[i][t] * V[i]) + ZIP[3]*(V_tf[i][t]^2))) # Equality constraint for balance of the reactive power
                #eq_const_Qbalance_tf[i][t] = JuMP.@constraint(model, terms_Qb_tf_dict[i][t] == Total_Q_inj - q_d*(ZIP[1] + ZIP[2]*(V_tf[i][t]) + ZIP[3]*(V_tf[i][t]^2)) ) # Equality constraint for balance of the active power

            end
        end
    end

    return model, eq_const_Pe_tf, eq_const_Qe_tf, eq_const_Pbalance_tf, eq_const_Qbalance_tf, eq_const_Vd_tf, eq_const_Vq_tf, Te_tf, eq_const_Te_tf

end

# Function that model the swing equation using Trapezoidal Approximation
function Def_Dyn_Fault_Swing!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    P_mech_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ::OrderedDict{Int64, VariableRef},
    Δω_0::Float64,
    P_g::OrderedDict{Int64, VariableRef},
    t_window_fault::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Δω_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 

    for i in 1:nGEN
        if DGEN.g_status[i] == 1

            eq_const_δ_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Δω_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            H = DGEN_DYN.H[i]
            D = DGEN_DYN.D[i]

            for t in eachindex(t_window_fault)                
                if t == 1
                    Te_prev = P_g[i]
                    # --- BACKWARD EULER FOR STEP 1 ---
                    # Ignores pre-fault Pe completely (avoids integrating the jump)
                    eq_const_δ_tf[i][t] = JuMP.@constraint(model, 
                        δ_tf[i][t] - δ[i] - (ω_syn * Δt * Δω_tf[i][t]) == 0.0)

                    # eq_const_Δω_tf[i][t] = JuMP.@constraint(model, 
                    #     (1.0 + ((D*Δt) / (2*H))) * Δω_tf[i][t] - Δω_0 - 
                    #     (Δt / (2*H))*(P_mech_tf[i][t] - Te_tf[i][t]-Te_prev) == 0.0)
                    eq_const_Δω_tf[i][t] = JuMP.@constraint(model, 
                        (1.0 + ((D*Δt) / (2*H))) * Δω_tf[i][t] - Δω_0 - 
                        (Δt / (2*H))*(P_mech_tf[i][t] - Pe_tf[i][t]) == 0.0)
                else
                    # --- TRAPEZOIDAL FOR REMAINDER ---
                    eq_const_δ_tf[i][t] = JuMP.@constraint(model, 
                        δ_tf[i][t] - δ_tf[i][t-1] - (ω_syn * (Δt/2) * (Δω_tf[i][t] + Δω_tf[i][t-1])) == 0.0)

                    eq_const_Δω_tf[i][t] = JuMP.@constraint(model, 
                        (1.0 + ((D*Δt) / (4*H))) * Δω_tf[i][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_tf[i][t-1] - 
                        (Δt / (4*H))*(P_mech_tf[i][t]+P_mech_tf[i][t-1] - Pe_tf[i][t] - Pe_tf[i][t-1])  == 0.0)
                end
            end
        end
    end
    return model, eq_const_δ_tf, eq_const_Δω_tf

end

function Def_Dyn_Fault_EMF!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    E_fd_tf::OrderedDict,
    Ed::OrderedDict,   # initial values (t=0)
    Eq::OrderedDict,
    Id::OrderedDict,
    Iq::OrderedDict,
    Ed_tf::OrderedDict, Eq_tf::OrderedDict,
    Id_tf::OrderedDict,   Iq_tf::OrderedDict,
    t_window_fault::Vector{Float64}, Δt::Float64
    )

    eq_const_Ed_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_Eq_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            Xd  = DGEN_DYN.Xd[i];  Xd_tr = DGEN_DYN.Xd_tr[i]
            Xq  = DGEN_DYN.Xq[i];  Xq_tr = DGEN_DYN.Xq_tr[i]
            Td  = DGEN_DYN.Td[i];  Tq    = DGEN_DYN.Tq[i]

            eq_const_Ed_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Eq_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_fault)
                if t == 1
                    # --- BACKWARD EULER FOR STEP 1 ---
                    # Ignores pre-fault discontinuous currents (Id, Iq)
                    Ed_prev = Ed[i]
                    Eq_prev = Eq[i]

                    eq_const_Ed_tf[i][t] = JuMP.@constraint(model,
                        Ed_tf[i][t] * (1 + Δt/Tq) - Ed_prev -
                        (Δt/Tq) * (Xq - Xq_tr) * Iq_tf[i][t] == 0.0)

                    eq_const_Eq_tf[i][t] = JuMP.@constraint(model,
                        Eq_tf[i][t] * (1 + Δt/Td) - Eq_prev -
                        (Δt/Td) * (E_fd_tf[i][t] - (Xd - Xd_tr) * Id_tf[i][t]) == 0.0)
                else
                    # --- TRAPEZOIDAL FOR REMAINDER ---
                    Ed_prev = Ed_tf[i][t-1]
                    Eq_prev = Eq_tf[i][t-1]
                    Iq_prev = Iq_tf[i][t-1]
                    Id_prev = Id_tf[i][t-1]

                    eq_const_Ed_tf[i][t] = JuMP.@constraint(model,
                        Ed_tf[i][t] * (1 + Δt/(2*Tq)) - Ed_prev * (1 - Δt/(2*Tq)) -
                        (Δt/(2*Tq)) * (Xq - Xq_tr) * (Iq_tf[i][t] + Iq_prev) == 0.0)

                    eq_const_Eq_tf[i][t] = JuMP.@constraint(model,
                        Eq_tf[i][t] * (1 + Δt/(2*Td)) - Eq_prev * (1 - Δt/(2*Td)) -
                        (Δt/(2*Td)) * (E_fd_tf[i][t] + E_fd_tf[i][t-1] - (Xd - Xd_tr)*(Id_tf[i][t] + Id_prev)) == 0.0)
                end
            end
        end
    end
    return model, eq_const_Ed_tf, eq_const_Eq_tf
end

# ===================================================================================
#                               POST-FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the post-fault period
function Def_Dyn_PostF_All!(
    model::Model,
    DBUS::DataFrame,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nBUS::Int64,
    nGEN::Int64,
    base_MVA::Float64,
    ZIP::Vector{Float64},
    t_window_postf::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64},
    Ybus_postf::SparseMatrixCSC,
    V::OrderedDict{Int64, VariableRef},
    θ::OrderedDict{Int64, VariableRef},
    E_fd::OrderedDict{Int64, VariableRef},
    V_ref::OrderedDict{Int64, VariableRef},
    P_ref::OrderedDict,
    Pm::OrderedDict,
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    V_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    E_fd_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Ed_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}}, 
    Eq_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}}, 
    Id_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}}, 
    Iq_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}}, 
    P_valve_tf::OrderedDict, P_mech_tf::OrderedDict,
    ω_syn::Float64,
    Δt::Float64,
    val_V::Dict, val_θ::Dict, val_Pg::Dict, val_Qg::Dict,
    δ::OrderedDict, Ed::OrderedDict, Eq::OrderedDict, Id::OrderedDict, Iq::OrderedDict 
    )
    # ======================================================================================
    #               Setting the variables used in the post-fault period
    # ======================================================================================

    Pe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 
    Qe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 
    δ_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 
    Δω_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 
    V_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 
    θ_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() 

    Ed_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Eq_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Id_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Iq_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()

    terms_δCOI_tpf_dict  = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() 
    terms_ΔωCOI_tpf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() 

    # Loop to create the network variables (V and θ)
    for i in 1:nBUS
        # Initialize inner dict for this bus
        V_tpf[i] = OrderedDict{Int, JuMP.VariableRef}()
        θ_tpf[i] = OrderedDict{Int, JuMP.VariableRef}()
        v_start = val_V[i]
        th_start = val_θ[i]

        for t in eachindex(t_window_postf)
            V_tpf[i][t] = JuMP.@variable(model, lower_bound = 0.0, base_name = "V_tpf[$i, $t]", start = v_start)
            θ_tpf[i][t] = JuMP.@variable(model, base_name = "θ_tpf[$i, $t]", start = th_start)
        end
    end

    # Loop to create the generator variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            Pe_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Qe_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            δ_tpf[i]      = OrderedDict{Int, JuMP.VariableRef}()
            Δω_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            
            Ed_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Eq_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Id_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            Iq_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()

            terms_δCOI_tpf_dict[i]  = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()
            terms_ΔωCOI_tpf_dict[i] = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()

            pg_start = val_Pg[i]
            qg_start = val_Qg[i]
            
            delta_start = JuMP.start_value(δ[i])
            Ed_start = JuMP.start_value(Ed[i])
            Eq_start = JuMP.start_value(Eq[i]) 
            Id_start = JuMP.start_value(Id[i]) 
            Iq_start = JuMP.start_value(Iq[i]) 

            for t in eachindex(t_window_postf)
                Pe_tpf[i][t] = JuMP.@variable(model, lower_bound = -Inf, upper_bound = Inf, base_name = "Pe_tpf[$i, $t]", start=pg_start)
                Qe_tpf[i][t] = JuMP.@variable(model, lower_bound = -Inf, upper_bound = Inf, base_name = "Qe_tpf[$i, $t]", start=qg_start)
                δ_tpf[i][t]  = JuMP.@variable(model, lower_bound = -Inf, upper_bound = Inf, base_name = "δ_tpf[$i, $t]", start=delta_start)
                Δω_tpf[i][t] = JuMP.@variable(model, lower_bound = -Inf, upper_bound = Inf, base_name = "Δω_tpf[$i, $t]", start=0.0)

                Ed_tpf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Ed_tpf[$i, $t]", start=Ed_start)
                Eq_tpf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Eq_tpf[$i, $t]", start=Eq_start)
                Id_tpf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Id_tpf[$i, $t]", start=Id_start)
                Iq_tpf[i][t] = JuMP.@variable(model,lower_bound = -Inf, upper_bound = Inf, base_name="Iq_tpf[$i, $t]", start=Iq_start)

                terms_δCOI_tpf_dict[i][t]  = DGEN_DYN.H[i] * δ_tpf[i][t]
                terms_ΔωCOI_tpf_dict[i][t] = DGEN_DYN.H[i] * Δω_tpf[i][t]
            end
        end
    end

    # Gives the expressions to calculate the angle and speed deviation of the COI
    expr_δCOI_per_time  = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_δCOI_tpf_dict))  for t in eachindex(t_window_postf))
    expr_ΔωCOI_per_time = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_ΔωCOI_tpf_dict)) for t in eachindex(t_window_postf))

    model, E_fd_tpf, eq_const_E_fd_tpf = Def_Dyn_PostF_Exciter!(model, DGEN, DGEN_DYN, nGEN, E_fd_tf, V_ref, V_tf, V_tpf, t_window_postf, Δt)

    model, P_valve_tpf, P_mech_tpf, eq_const_P_valve_tpf, eq_const_P_mech_tpf = Def_Dyn_PostF_Governor!(model, DGEN, DGEN_DYN, nGEN, P_ref, P_valve_tf, P_mech_tf, Δω_tf, Δω_tpf, t_window_postf, Δt, base_MVA)

    # Define COI variables and constraints
    model, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf = Def_Dyn_PostF_COI!(model, DGEN_DYN, nGEN, t_window_postf, expr_δCOI_per_time, expr_ΔωCOI_per_time)

    # Define COI limits
    model, ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf = Def_Dyn_PostF_rel_COI!(model, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, nGEN, t_window_postf, δ_tol, Δω_tol)

    # Define the equality constraints for Pe, Qe, and Nodal Balances in the post-fault period 
    model, eq_const_Pe_tpf, eq_const_Qe_tpf, eq_const_Pbalance_tpf, eq_const_Qbalance_tpf, eq_const_Vd_tpf, eq_const_Vq_tpf, Te_tpf, eq_const_Te_tpf = Def_Dyn_PostF_eqconst_PeQe!(model, DBUS, nBUS, bus_gen_circ_dict, DGEN, DGEN_DYN, nGEN, base_MVA, V, θ, E_fd, Pe_tpf, Qe_tpf, δ_tpf, V_tpf, θ_tpf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf,Δω_tpf, t_window_postf, Ybus_postf, ZIP)

    # Define the swing equations
    model, eq_const_δ_tpf, eq_const_Δω_tpf = Def_Dyn_PostF_Swing!(model, DGEN, DGEN_DYN, nGEN, P_mech_tpf, Pe_tpf, δ_tpf, Δω_tpf, δ_tf, Δω_tf, Pe_tpf, Pe_tf, t_window_postf, ω_syn, Δt)    # Define 4th order EMF dynamics
    
    model, eq_const_Ed_tpf, eq_const_Eq_tpf = Def_Dyn_PostF_EMF!(model, DGEN, DGEN_DYN, nGEN, E_fd_tpf, Ed_tf, Eq_tf, Id_tf, Iq_tf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, t_window_postf, Δt)

    return model, Pe_tpf, Qe_tpf, V_tpf, θ_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf, E_fd_tpf, P_valve_tpf, P_mech_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf, eq_const_Pe_tpf, eq_const_Qe_tpf, eq_const_Pbalance_tpf, eq_const_Qbalance_tpf, eq_const_Vd_tpf, eq_const_Vq_tpf, eq_const_δ_tpf, eq_const_Δω_tpf, eq_const_Ed_tpf, eq_const_Eq_tpf, eq_const_E_fd_tpf, eq_const_P_valve_tpf, eq_const_P_mech_tpf
end

# Function that define the variables of the COI as well as its equality constraints across the post-fault period
function Def_Dyn_PostF_COI!(model::Model,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_postf::Vector{Float64},
    aux_expr_δCOI_per_time::OrderedDict{Int64, AffExpr},
    aux_expr_ΔωCOI_per_time::OrderedDict{Int64, AffExpr}
    )

    H_COI = sum(DGEN_DYN.H[i] for i in 1:nGEN if DGEN.g_status[i] == 1) 

    δCOI_tpf  = OrderedDict{Int, JuMP.VariableRef}() 
    ΔωCOI_tpf = OrderedDict{Int, JuMP.VariableRef}() 

    eq_const_δCOI_tpf  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_ΔωCOI_tpf = OrderedDict{Int, JuMP.ConstraintRef}()

    for t in eachindex(t_window_postf)
        δCOI_tpf[t]  = JuMP.@variable(model, base_name = "δCOI_tpf[$t]",  start = 0.0)
        ΔωCOI_tpf[t] = JuMP.@variable(model, base_name = "ΔωCOI_tpf[$t]", start = 0.0)
        
        eq_const_δCOI_tpf[t]  = JuMP.@constraint(model, δCOI_tpf[t]  == aux_expr_δCOI_per_time[t]  / H_COI)
        eq_const_ΔωCOI_tpf[t] = JuMP.@constraint(model, ΔωCOI_tpf[t] == aux_expr_ΔωCOI_per_time[t] / H_COI)
    end

    return model, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf
end

# Function that define the constraints of the angle of generators in relation to the COI variables across the post-fault period
function Def_Dyn_PostF_rel_COI!(model::Model,
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δCOI_tpf::OrderedDict{Int, JuMP.VariableRef},
    ΔωCOI_tpf::OrderedDict{Int, JuMP.VariableRef},
    nGEN::Int64,
    t_window_postf::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64}
    )

    ineq_const_δ_COI_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    ineq_const_Δω_COI_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            ineq_const_δ_COI_tpf[i]  = OrderedDict{Int, JuMP.ConstraintRef}()
            ineq_const_Δω_COI_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_postf)
                ineq_const_δ_COI_tpf[i][t] = JuMP.@constraint(model, δ_tol[1] <= δ_tpf[i][t] - δCOI_tpf[t] <= δ_tol[2])
                ineq_const_Δω_COI_tpf[i][t] = JuMP.@constraint(model, Δω_tol[1]<= Δω_tpf[i][t] - ΔωCOI_tpf[t] <= Δω_tol[2])
            end
        end
    end

    return model, ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf
end

# Function that creates the constraints of the electrical power across the post-fault period (4th order)
function Def_Dyn_PostF_eqconst_PeQe!(model::Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    base_MVA::Float64,
    V::OrderedDict{Int64, VariableRef},
    θ::OrderedDict{Int64, VariableRef},
    E_fd::OrderedDict{Int64, VariableRef},
    Pe_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Qe_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    V_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    θ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Ed_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Eq_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Id_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Iq_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_postf::Vector{Float64},
    Ybus_postf::SparseMatrixCSC,
    ZIP::Vector{Float64}
    )
    
    eq_const_Pe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Qe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Vd_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Vq_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 

    eq_const_Pbalance_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Qbalance_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 

    terms_Pb_tpf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() 
    terms_Qb_tpf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() 

    Te_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq_const_Te_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    active_gen = findall(isone, DGEN.g_status) 

    # Generate the expression of power generated by each generator
    for g in active_gen
        bus_g = DGEN.bus[g]
        Xd_tr = DGEN_DYN.Xd_tr[g]; Xq_tr = DGEN_DYN.Xq_tr[g]; Ra = DGEN_DYN.Ra[g]

        eq_const_Pe_tpf[g] = OrderedDict{Int, JuMP.ConstraintRef}() 
        eq_const_Qe_tpf[g] = OrderedDict{Int, JuMP.ConstraintRef}() 
        eq_const_Vd_tpf[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_const_Vq_tpf[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        
        Te_tpf[g]          = OrderedDict{Int, JuMP.VariableRef}()
        eq_const_Te_tpf[g] = OrderedDict{Int, JuMP.ConstraintRef}()

        for t in eachindex(t_window_postf) 

            Te_tpf[g][t] = JuMP.@variable(model,
                lower_bound = -Inf, upper_bound = Inf,
                base_name = "Te_tpf[$g,$t]",
                start = val_Pg[g] / base_MVA)

            eq_const_Te_tpf[g][t] = JuMP.@constraint(model,
                Te_tpf[g][t] == Ed_tpf[g][t]*Id_tpf[g][t] + Eq_tpf[g][t]*Iq_tpf[g][t])


            Vd_expr = @expression(model, V_tpf[bus_g][t] * sin(δ_tpf[g][t] - θ_tpf[bus_g][t]))
            Vq_expr = @expression(model, V_tpf[bus_g][t] * cos(δ_tpf[g][t] - θ_tpf[bus_g][t]))

            # Stator algebraic equations linking internal voltage to terminal voltage
            eq_const_Vd_tpf[g][t] = JuMP.@constraint(model, Vd_expr - (1 + Δω_tpf[g][t])*Ed_tpf[g][t] + Ra * Id_tpf[g][t] - Xq_tr * Iq_tpf[g][t] == 0)
            eq_const_Vq_tpf[g][t] = JuMP.@constraint(model, Vq_expr - (1 + Δω_tpf[g][t])*Eq_tpf[g][t] + Ra * Iq_tpf[g][t] + Xd_tr * Id_tpf[g][t] == 0)
            
            eq_const_Pe_tpf[g][t] = JuMP.@constraint(model, Pe_tpf[g][t] == (1 + Δω_tpf[g][t])*Te_tpf[g][t]) # Equality constraint for the active power generated
            #eq_const_Pe_tpf[g][t] = JuMP.@constraint(model, Pe_tpf[g][t] == Vq_expr*Iq_tpf[g][t] + Vd_expr*Id_tpf[g][t]) 
            eq_const_Qe_tpf[g][t] = JuMP.@constraint(model, Qe_tpf[g][t] == Vq_expr*Id_tpf[g][t] - Vd_expr*Iq_tpf[g][t]) 
        end
    end

    # Generate the expressions for the power balance during the post-fault
    for i in 1:nBUS
        indices_circ_connected = bus_gen_circ_dict[i][:circ]    

        if isempty(indices_circ_connected) 
            throw(ArgumentError("The bus $i is islanded, i.e., there is no line or transformer connected to it."))
        end

        terms_Pb_tpf_dict[i] = OrderedDict{Int, JuMP.NonlinearExpr}() 
        terms_Qb_tpf_dict[i] = OrderedDict{Int, JuMP.NonlinearExpr}() 

        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids] 

        eq_const_Pbalance_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}() 
        eq_const_Qbalance_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}() 

        p_d = DBUS.p_d[i] / base_MVA
        q_d = (1.0 .* DBUS.q_d[i]) / base_MVA 

        for t in eachindex(t_window_postf) 
            row_Y = Ybus_postf[i, :] 

            terms_Pb_tpf_dict[i][t] = @expression(model, 
                V_tpf[i][t] * sum(
                    V_tpf[j][t] * (real(y_ij) * cos(θ_tpf[i][t] - θ_tpf[j][t]) + imag(y_ij) * sin(θ_tpf[i][t] - θ_tpf[j][t]))
                    for (j, y_ij) in zip(row_Y.nzind, row_Y.nzval)
                )
            )

            terms_Qb_tpf_dict[i][t] = @expression(model, 
                V_tpf[i][t] * sum(
                    V_tpf[j][t] * (real(y_ij) * sin(θ_tpf[i][t] - θ_tpf[j][t]) - imag(y_ij) * cos(θ_tpf[i][t] - θ_tpf[j][t]))
                    for (j, y_ij) in zip(row_Y.nzind, row_Y.nzval)
                )
            )
            
            Total_P_inj = isempty(indices_bus_gen) ? 0.0 : @expression(model, sum(V_tpf[i][t] * Id_tpf[g][t] * sin(δ_tpf[g][t] - θ_tpf[i][t]) + V_tpf[i][t] * Iq_tpf[g][t] * cos(δ_tpf[g][t] - θ_tpf[i][t]) for g in indices_bus_gen))
            Total_Q_inj = isempty(indices_bus_gen) ? 0.0 : @expression(model, sum(V_tpf[i][t] * Id_tpf[g][t] * cos(δ_tpf[g][t] - θ_tpf[i][t]) - V_tpf[i][t] * Iq_tpf[g][t] * sin(δ_tpf[g][t] - θ_tpf[i][t]) for g in indices_bus_gen))

            if p_d == 0.0
                eq_const_Pbalance_tpf[i][t] = JuMP.@constraint(model, terms_Pb_tpf_dict[i][t] == Total_P_inj) 
            else
                eq_const_Pbalance_tpf[i][t] = JuMP.@constraint(model, terms_Pb_tpf_dict[i][t]*V[i]^2 == Total_P_inj*V[i]^2 - p_d*(ZIP[1]*V[i]^2 + ZIP[2]*(V_tpf[i][t] * V[i]) + ZIP[3]*(V_tpf[i][t]^2)) ) 
                #eq_const_Pbalance_tpf[i][t] = JuMP.@constraint(model, terms_Pb_tpf_dict[i][t] == Total_P_inj - p_d*(ZIP[1] + ZIP[2]*(V_tpf[i][t]) + ZIP[3]*(V_tpf[i][t]^2)) ) 
            end
            
            if q_d == 0.0
                eq_const_Qbalance_tpf[i][t] = JuMP.@constraint(model, terms_Qb_tpf_dict[i][t] == Total_Q_inj) 
            else
                eq_const_Qbalance_tpf[i][t] = JuMP.@constraint(model, terms_Qb_tpf_dict[i][t]*V[i]^2 == Total_Q_inj*V[i]^2 - q_d*(ZIP[1]*V[i]^2 + ZIP[2]*(V_tpf[i][t] * V[i]) + ZIP[3]*(V_tpf[i][t]^2))) 
                #eq_const_Qbalance_tpf[i][t] = JuMP.@constraint(model, terms_Qb_tpf_dict[i][t] == Total_Q_inj - q_d*(ZIP[1] + ZIP[2]*(V_tpf[i][t]) + ZIP[3]*(V_tpf[i][t]^2)) )             
            end
        end
    end

    return model, eq_const_Pe_tpf, eq_const_Qe_tpf, eq_const_Pbalance_tpf, eq_const_Qbalance_tpf, eq_const_Vd_tpf, eq_const_Vq_tpf, Te_tpf, eq_const_Te_tpf
end

# Function that model the swing equation using Trapezoidal Approximation (Post-Fault period)
function Def_Dyn_PostF_Swing!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    P_mech_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Pe_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Te_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Te_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_postf::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    eq_const_Δω_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() 
    
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            eq_const_δ_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Δω_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            H = DGEN_DYN.H[i]
            D = DGEN_DYN.D[i]
            
            for t in eachindex(t_window_postf)
                if t == 1
                    # --- BACKWARD EULER FOR STEP 1 ---
                    # Ignores discontinuous last_var_Pe_tf
                    last_var_δ_tf = last(δ_tf[i])[2]   
                    last_var_Δω_tf = last(Δω_tf[i])[2] 
                    last_var_Pe_tf = last(Pe_tf[i])[2]

                    eq_const_δ_tpf[i][t] = JuMP.@constraint(model, 
                        δ_tpf[i][t] - last_var_δ_tf - (ω_syn * Δt * Δω_tpf[i][t]) == 0.0)
                        
                    # eq_const_Δω_tpf[i][t] = JuMP.@constraint(model, 
                    #     (1.0 + ((D*Δt) / (2*H))) * Δω_tpf[i][t] - last_var_Δω_tf - 
                    #     (Δt / (2*H))*(P_mech_tpf[i][t] - Te_tpf[i][t]- last_var_Te_tf) == 0.0)
                    eq_const_Δω_tpf[i][t] = JuMP.@constraint(model, 
                        (1.0 + ((D*Δt) / (2*H))) * Δω_tpf[i][t] - last_var_Δω_tf - 
                        (Δt / (2*H))*(P_mech_tpf[i][t] - Pe_tpf[i][t]) == 0.0)
                else
                    # --- TRAPEZOIDAL FOR REMAINDER ---
                    eq_const_δ_tpf[i][t] = JuMP.@constraint(model, 
                        δ_tpf[i][t] - δ_tpf[i][t-1] - (ω_syn * (Δt/2) * (Δω_tpf[i][t] + Δω_tpf[i][t-1])) == 0.0)
                        
                    eq_const_Δω_tpf[i][t] = JuMP.@constraint(model, 
                        (1.0 + ((D*Δt) / (4*H))) * Δω_tpf[i][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_tpf[i][t-1] - 
                        (Δt / (4*H))*(P_mech_tpf[i][t]+P_mech_tpf[i][t-1] - Pe_tpf[i][t] - Pe_tpf[i][t-1]) == 0.0)
                end
            end
        end
    end
    return model, eq_const_δ_tpf, eq_const_Δω_tpf
end

# Function that applies Trapezoidal Approximation to Generator Electromotive Forces (EMF)
function Def_Dyn_PostF_EMF!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    E_fd_tpf::OrderedDict,
    Ed_tf::OrderedDict, Eq_tf::OrderedDict, # last states (t=end) from fault
    Id_tf::OrderedDict, Iq_tf::OrderedDict, 
    Ed_tpf::OrderedDict, Eq_tpf::OrderedDict,
    Id_tpf::OrderedDict, Iq_tpf::OrderedDict,
    t_window_postf::Vector{Float64}, Δt::Float64
    )

    eq_const_Ed_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_Eq_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            Xd = DGEN_DYN.Xd[i]; Xd_tr = DGEN_DYN.Xd_tr[i]
            Xq = DGEN_DYN.Xq[i]; Xq_tr = DGEN_DYN.Xq_tr[i]
            Td = DGEN_DYN.Td[i]; Tq  = DGEN_DYN.Tq[i]

            eq_const_Ed_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Eq_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_postf)
                if t == 1
                    # --- BACKWARD EULER FOR STEP 1 ---
                    # Ignores discontinuous last_var_Iq_tf
                    Ed_prev = last(Ed_tf[i])[2]
                    Eq_prev = last(Eq_tf[i])[2]

                    eq_const_Ed_tpf[i][t] = JuMP.@constraint(model,
                        Ed_tpf[i][t] * (1 + Δt/Tq) - Ed_prev -
                        (Δt/Tq) * (Xq - Xq_tr) * Iq_tpf[i][t] == 0.0)

                    eq_const_Eq_tpf[i][t] = JuMP.@constraint(model,
                        Eq_tpf[i][t] * (1 + Δt/Td) - Eq_prev -
                        (Δt/Td) * (E_fd_tpf[i][t] - (Xd - Xd_tr) * Id_tpf[i][t]) == 0.0)
                else
                    # --- TRAPEZOIDAL FOR REMAINDER ---
                    Ed_prev = Ed_tpf[i][t-1]
                    Eq_prev = Eq_tpf[i][t-1]
                    Iq_prev = Iq_tpf[i][t-1]
                    Id_prev = Id_tpf[i][t-1]

                    eq_const_Ed_tpf[i][t] = JuMP.@constraint(model,
                        Ed_tpf[i][t] * (1 + Δt/(2*Tq)) - Ed_prev * (1 - Δt/(2*Tq)) -
                        (Δt/(2*Tq)) * (Xq - Xq_tr) * (Iq_tpf[i][t] + Iq_prev) == 0.0)

                    eq_const_Eq_tpf[i][t] = JuMP.@constraint(model,
                        Eq_tpf[i][t] * (1 + Δt/(2*Td)) - Eq_prev * (1 - Δt/(2*Td)) -
                        (Δt/(2*Td)) * (E_fd_tpf[i][t] + E_fd_tpf[i][t-1] - (Xd - Xd_tr)*(Id_tpf[i][t] + Id_prev)) == 0.0)
                end
            end
        end
    end
    return model, eq_const_Ed_tpf, eq_const_Eq_tpf
end

# ===================================================================
#             EXCITER 
# ================================================================

# ===================================================================================
#                      PRE-FAULT EXCITER / AVR INITIALIZATION
# ===================================================================================
# Completely decoupled function to initialize Exciter/AVR models pre-fault
function Def_Dyn_Init_Exciter!(model::Model, DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64, V::OrderedDict{Int, JuMP.VariableRef}, E_fd::OrderedDict{Int, JuMP.VariableRef}, val_V::Dict)
    
    V_ref = OrderedDict{Int, JuMP.VariableRef}()
    eq_const_Efd_init = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]
            K_exc = DGEN_DYN.K_exc[i]
            
            # Use the E_fd start value automatically guessed by the generator's internal calculation
            start_Efd = JuMP.start_value(E_fd[i])
            v_val = val_V[bus]

            V_ref[i] = JuMP.@variable(model,  
                lower_bound = 0.0,   
                base_name = "V_ref[$i]", 
                start = v_val + start_Efd / K_exc          
            )

            # Link the generator's field voltage to the exciter equation
            eq_const_Efd_init[i] = JuMP.@constraint(model, E_fd[i] == K_exc * (V_ref[i] - V[bus]))
        end
    end

    return model, V_ref, eq_const_Efd_init
end

# Function defining dynamic Exciter variables and constraints during Fault
function Def_Dyn_Fault_Exciter!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    E_fd::OrderedDict, V_ref::OrderedDict, 
    V_tf::OrderedDict, t_window_fault::Vector{Float64}, Δt::Float64)

    E_fd_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq_const_E_fd_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    E_fd_unlim_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq_const_E_fd_unlim_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    # Smooth Min/Max parameters
    rho = 0.0001  # Smoothing factor (keeps gradients stable for IPOPT)
    E_max = 4.5 # Exciter Ceiling
    E_min = 0.0 # Exciter Floor

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]
            K_exc = DGEN_DYN.K_exc[i]
            T_exc = DGEN_DYN.T_exc[i]
            
            E_fd_tf[i] = OrderedDict{Int, JuMP.VariableRef}()
            eq_const_E_fd_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            E_fd_unlim_tf[i] = OrderedDict{Int, JuMP.VariableRef}()
            eq_const_E_fd_unlim_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            
            Efd_start = JuMP.start_value(E_fd[i])
            
            for t in eachindex(t_window_fault)
                E_fd_tf[i][t] = JuMP.@variable(model, 
                    # lower_bound = -9999, 
                    # upper_bound = 9999, 
                    base_name = "E_fd_tf[$i, $t]", 
                    start = Efd_start)

                E_fd_unlim_tf[i][t] = JuMP.@variable(model, base_name = "E_fd_unlim_tf[$i, $t]", start = Efd_start)
                
                if t == 1
                    # --- BACKWARD EULER FOR STEP 1 ---
                    eq_const_E_fd_tf[i][t] = JuMP.@constraint(model,
                        (E_fd_unlim_tf[i][t] * (1 + Δt / T_exc))/K_exc - E_fd[i]/K_exc - 
                        (Δt / T_exc) * (V_ref[i] - V_tf[bus][t]) == 0.0)
                else
                    # --- TRAPEZOIDAL FOR REMAINDER ---
                    E_fd_prev = E_fd_tf[i][t-1]
                    V_prev = V_tf[bus][t-1]
                    V_curr = V_tf[bus][t]
                    
                    eq_const_E_fd_tf[i][t] = JuMP.@constraint(model,
                        (E_fd_unlim_tf[i][t] * (1 + Δt / (2 * T_exc)))/K_exc - 
                        (E_fd_prev * (1 - Δt / (2 * T_exc)))/K_exc - 
                        (Δt / (2 * T_exc)) * (2 * V_ref[i] - V_curr - V_prev) == 0.0)
                end
                # --- MATHEMATICAL ANTI-WINDUP LIMITER (SMOOTH MIN/MAX) ---
                # 1. Smooth Minimum (Limits to E_max)
                E_aux = @expression(model, (E_fd_unlim_tf[i][t] + E_max - sqrt((E_fd_unlim_tf[i][t] - E_max)^2 + rho)) / 2)
                # 2. Smooth Maximum (Limits to E_min)
                eq_const_E_fd_tf[i][t] = JuMP.@constraint(model, 
                    E_fd_tf[i][t] == (E_aux + E_min + sqrt((E_aux - E_min)^2 + rho)) / 2)
            end
        end
    end
    
    return model, E_fd_tf, eq_const_E_fd_tf
end

# Function defining dynamic Exciter variables and constraints during Post-Fault
function Def_Dyn_PostF_Exciter!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    E_fd_tf::OrderedDict, V_ref::OrderedDict, V_tf::OrderedDict,
    V_tpf::OrderedDict, t_window_postf::Vector{Float64}, Δt::Float64)

    E_fd_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    E_fd_unlim_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    
    eq_const_E_fd_unlim_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_E_fd_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    # Smooth Min/Max parameters
    rho = 0.0001
    E_max = 4.5
    E_min = 0.0

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]
            K_exc = DGEN_DYN.K_exc[i]
            T_exc = DGEN_DYN.T_exc[i]
            
            E_fd_tpf[i] = OrderedDict{Int, JuMP.VariableRef}()
            E_fd_unlim_tpf[i] = OrderedDict{Int, JuMP.VariableRef}()
            
            eq_const_E_fd_unlim_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_E_fd_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            
            Efd_start = JuMP.start_value(last(E_fd_tf[i])[2])
            
            for t in eachindex(t_window_postf)
                # NO HARD BOUNDS
                E_fd_tpf[i][t] = JuMP.@variable(model, base_name = "E_fd_tpf[$i, $t]", start = Efd_start)
                E_fd_unlim_tpf[i][t] = JuMP.@variable(model, base_name = "E_fd_unlim_tpf[$i, $t]", start = Efd_start)
                
                if t == 1
                    # --- SCALED BACKWARD EULER ---
                    E_fd_prev = last(E_fd_tf[i])[2]
                    V_prev = last(V_tf[bus])[2]
                    V_curr = V_tpf[bus][t]
                    
                    eq_const_E_fd_unlim_tpf[i][t] = JuMP.@constraint(model,
                        (E_fd_unlim_tpf[i][t] / K_exc) * (1 + Δt / T_exc) - (E_fd_prev / K_exc) - 
                        (Δt / T_exc) * (V_ref[i] - V_curr) == 0.0)
                else
                    # --- SCALED TRAPEZOIDAL ---
                    E_fd_prev = E_fd_tpf[i][t-1]
                    V_prev = V_tpf[bus][t-1]
                    V_curr = V_tpf[bus][t]
                    
                    eq_const_E_fd_unlim_tpf[i][t] = JuMP.@constraint(model,
                        (E_fd_unlim_tpf[i][t] / K_exc) * (1 + Δt / (2 * T_exc)) - 
                        (E_fd_prev / K_exc) * (1 - Δt / (2 * T_exc)) - 
                        (Δt / (2 * T_exc)) * (2 * V_ref[i] - V_curr - V_prev) == 0.0)
                end
                
                # --- MATHEMATICAL ANTI-WINDUP LIMITER (SMOOTH MIN/MAX) ---
                E_aux = @expression(model, (E_fd_unlim_tpf[i][t] + E_max - sqrt((E_fd_unlim_tpf[i][t] - E_max)^2 + rho)) / 2)
                eq_const_E_fd_tpf[i][t] = JuMP.@constraint(model, 
                    E_fd_tpf[i][t] == (E_aux + E_min + sqrt((E_aux - E_min)^2 + rho)) / 2)
            end
        end
    end
    return model, E_fd_tpf, eq_const_E_fd_tpf
end

# ===================================================================================
#                      GOVERNOR MODEL IMPLEMENTATION
# ===================================================================================

# Function to initialize Governor variables pre-fault
function Def_Dyn_Init_Governor!(model::Model, DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64, Pm::OrderedDict)
    P_ref   = OrderedDict{Int, JuMP.VariableRef}()
    P_valve = OrderedDict{Int, JuMP.VariableRef}()
    
    eq_const_Pref_init   = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_Pvalve_init = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            R = DGEN_DYN.R[i]
            # Since Pm is the steady state initial mechanical power, we use it to calculate Pref
            start_Pm = JuMP.start_value(Pm[i])

            # Pref is a constant setpoint during the transient
            P_ref[i] = JuMP.@variable(model,  
                base_name = "P_ref[$i]", 
                start = start_Pm * R          
            )
            P_valve[i] = JuMP.@variable(model,  
                base_name = "P_valve[$i]", 
                start = start_Pm          
            )

            # Guarantees initialization: pmech0 = pvalve0 = (pref)/R
            eq_const_Pref_init[i]   = JuMP.@constraint(model, P_ref[i] == Pm[i] * R)
            eq_const_Pvalve_init[i] = JuMP.@constraint(model, P_valve[i] == Pm[i])
        end
    end

    return model, P_ref, P_valve, eq_const_Pref_init, eq_const_Pvalve_init
end

# Function defining dynamic Governor variables and constraints during Fault
function Def_Dyn_Fault_Governor!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    P_ref::OrderedDict, P_valve::OrderedDict, Pm::OrderedDict, 
    Δω_tf::OrderedDict, t_window_fault::Vector{Float64}, Δt::Float64, base_MVA::Float64)

    P_valve_tf       = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    P_valve_unlim_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    P_mech_tf        = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    
    eq_const_P_valve_unlim_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_P_valve_tf       = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_P_mech_tf        = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    rho = 0.0001  # Smoothing factor

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            T1 = DGEN_DYN.T1[i]
            T2 = DGEN_DYN.T2[i]
            T3 = DGEN_DYN.T3[i]
            R  = DGEN_DYN.R[i]
            
            # Dynamic Limits adjusted to System Base
            P_max = DGEN.pg_max[i] / base_MVA
            P_min = 0.0 
            
            P_valve_tf[i]       = OrderedDict{Int, JuMP.VariableRef}()
            P_valve_unlim_tf[i] = OrderedDict{Int, JuMP.VariableRef}()
            P_mech_tf[i]        = OrderedDict{Int, JuMP.VariableRef}()
            
            eq_const_P_valve_unlim_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_P_valve_tf[i]       = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_P_mech_tf[i]        = OrderedDict{Int, JuMP.ConstraintRef}()
            
            P_mech_start  = JuMP.start_value(Pm[i])
            P_valve_start = JuMP.start_value(P_valve[i])
            
            for t in eachindex(t_window_fault)
                P_valve_tf[i][t]       = JuMP.@variable(model, base_name = "P_valve_tf[$i, $t]", start = P_valve_start)
                P_valve_unlim_tf[i][t] = JuMP.@variable(model, base_name = "P_valve_unlim_tf[$i, $t]", start = P_valve_start)
                P_mech_tf[i][t]        = JuMP.@variable(model, base_name = "P_mech_tf[$i, $t]",  start = P_mech_start)
                
                if t == 1
                    P_v_prev = P_valve[i] 
                    P_m_prev = Pm[i]
                    Δω_curr  = Δω_tf[i][t]
                    
                    eq_const_P_valve_unlim_tf[i][t] = JuMP.@constraint(model,
                        P_valve_unlim_tf[i][t] * (1 + Δt / T1) - P_v_prev - (Δt / (R * T1)) * (P_ref[i] - Δω_curr) == 0.0)

                    eq_const_P_mech_tf[i][t] = JuMP.@constraint(model,
                        P_mech_tf[i][t] * (1 + Δt / T3) - P_m_prev -
                        (Δt / T3) * ((1 - T2/T1) * P_valve_tf[i][t] + (T2 / T1) * (P_ref[i] - Δω_curr) / R) == 0.0)
                else
                    P_v_prev = P_valve_tf[i][t-1] 
                    P_m_prev = P_mech_tf[i][t-1]
                    Δω_prev  = Δω_tf[i][t-1]
                    Δω_curr  = Δω_tf[i][t]
                    
                    eq_const_P_valve_unlim_tf[i][t] = JuMP.@constraint(model,
                        P_valve_unlim_tf[i][t] * (1 + Δt / (2*T1)) - P_v_prev * (1 - Δt / (2*T1)) -
                        (Δt / (2*R*T1)) * (2*P_ref[i] - Δω_prev - Δω_curr) == 0.0)

                    eq_const_P_mech_tf[i][t] = JuMP.@constraint(model,
                        P_mech_tf[i][t] * (1 + Δt / (2*T3)) - P_m_prev * (1 - Δt / (2*T3)) -
                        (Δt / (2*T3)) * ((1 - T2/T1) * (P_valve_tf[i][t] + P_v_prev) + (T2/T1) * (2*P_ref[i] - Δω_prev - Δω_curr)/R) == 0.0)
                end
                
                P_aux = @expression(model, (P_valve_unlim_tf[i][t] + P_max - sqrt((P_valve_unlim_tf[i][t] - P_max)^2 + rho)) / 2)
                eq_const_P_valve_tf[i][t] = JuMP.@constraint(model, 
                    P_valve_tf[i][t] == (P_aux + P_min + sqrt((P_aux - P_min)^2 + rho)) / 2)
            end
        end
    end
    return model, P_valve_tf, P_mech_tf, eq_const_P_valve_tf, eq_const_P_mech_tf
end

# Function defining dynamic Governor variables and constraints during Post-Fault
function Def_Dyn_PostF_Governor!(model::Model,
    DGEN::DataFrame, DGEN_DYN::DataFrame, nGEN::Int64,
    P_ref::OrderedDict, P_valve_tf::OrderedDict, P_mech_tf::OrderedDict, 
    Δω_tf::OrderedDict, Δω_tpf::OrderedDict, t_window_postf::Vector{Float64}, Δt::Float64, base_MVA::Float64)

    P_valve_tpf       = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    P_valve_unlim_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    P_mech_tpf        = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    
    eq_const_P_valve_unlim_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_P_valve_tpf       = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_const_P_mech_tpf        = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    rho = 0.0001

    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            T1 = DGEN_DYN.T1[i]
            T2 = DGEN_DYN.T2[i]
            T3 = DGEN_DYN.T3[i]
            R  = DGEN_DYN.R[i]
            
            # Dynamic Limits adjusted to System Base
            P_max = DGEN.pg_max[i] / base_MVA
            P_min = 0.0 
            
            P_valve_tpf[i]       = OrderedDict{Int, JuMP.VariableRef}()
            P_valve_unlim_tpf[i] = OrderedDict{Int, JuMP.VariableRef}()
            P_mech_tpf[i]        = OrderedDict{Int, JuMP.VariableRef}()
            
            eq_const_P_valve_unlim_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_P_valve_tpf[i]       = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_P_mech_tpf[i]        = OrderedDict{Int, JuMP.ConstraintRef}()
            
            P_mech_start  = JuMP.start_value(last(P_mech_tf[i])[2])
            P_valve_start = JuMP.start_value(last(P_valve_tf[i])[2])
            
            for t in eachindex(t_window_postf)
                P_valve_tpf[i][t]       = JuMP.@variable(model, base_name = "P_valve_tpf[$i, $t]", start = P_valve_start)
                P_valve_unlim_tpf[i][t] = JuMP.@variable(model, base_name = "P_valve_unlim_tpf[$i, $t]", start = P_valve_start)
                P_mech_tpf[i][t]        = JuMP.@variable(model, base_name = "P_mech_tpf[$i, $t]",  start = P_mech_start)
                
                if t == 1
                    P_v_prev = last(P_valve_tf[i])[2]
                    P_m_prev = last(P_mech_tf[i])[2]
                    Δω_curr  = Δω_tpf[i][t]
                    
                    eq_const_P_valve_unlim_tpf[i][t] = JuMP.@constraint(model,
                        P_valve_unlim_tpf[i][t] * (1 + Δt / T1) - P_v_prev - (Δt / (R * T1)) * (P_ref[i] - Δω_curr) == 0.0)

                    eq_const_P_mech_tpf[i][t] = JuMP.@constraint(model,
                        P_mech_tpf[i][t] * (1 + Δt / T3) - P_m_prev -
                        (Δt / T3) * ((1 - T2/T1) * P_valve_tpf[i][t] + (T2 / T1) * (P_ref[i] - Δω_curr) / R) == 0.0)
                else
                    P_v_prev = P_valve_tpf[i][t-1]
                    P_m_prev = P_mech_tpf[i][t-1]
                    Δω_prev  = Δω_tpf[i][t-1]
                    Δω_curr  = Δω_tpf[i][t]
                    
                    eq_const_P_valve_unlim_tpf[i][t] = JuMP.@constraint(model,
                        P_valve_unlim_tpf[i][t] * (1 + Δt / (2*T1)) - P_v_prev * (1 - Δt / (2*T1)) -
                        (Δt / (2*R*T1)) * (2*P_ref[i] - Δω_prev - Δω_curr) == 0.0)

                    eq_const_P_mech_tpf[i][t] = JuMP.@constraint(model,
                        P_mech_tpf[i][t] * (1 + Δt / (2*T3)) - P_m_prev * (1 - Δt / (2*T3)) -
                        (Δt / (2*T3)) * ((1 - T2/T1) * (P_valve_tpf[i][t] + P_v_prev) + (T2/T1) * (2*P_ref[i] - Δω_prev - Δω_curr)/R) == 0.0)
                end
                
                P_aux = @expression(model, (P_valve_unlim_tpf[i][t] + P_max - sqrt((P_valve_unlim_tpf[i][t] - P_max)^2 + rho)) / 2)
                eq_const_P_valve_tpf[i][t] = JuMP.@constraint(model, 
                    P_valve_tpf[i][t] == (P_aux + P_min + sqrt((P_aux - P_min)^2 + rho)) / 2)
            end
        end
    end
    return model, P_valve_tpf, P_mech_tpf, eq_const_P_valve_tpf, eq_const_P_mech_tpf
end