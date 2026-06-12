# ===================================================================================
#                  PRINT THE INPUT PARAMETERS IN TXT FILE
# ===================================================================================
function Print_Input_Parameters(case::String, 
    base_MVA::Float64, 
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64},
    f_syn::Float64,
    bus_fault::Int,
    circ_trip::Vector,
    t_start_sim::Float64,
    t_end_sim::Float64,
    t_step::Float64,
    t_start_fault::Float64,
    clearing_time::Float64,
    t_clear_fault::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(path_folder_results)

    filename = "input_parameters.txt"
    open(filename, "w") do io
        println(io, "******** Simulation Input Parameters ********")
        println(io, "==============================================")
        println(io, "Case:                         $case")
        println(io, "Base MVA:                     $base_MVA")
        println(io, "δ tolerance (deg):[$(rad2deg(δ_tol[1])); +$(rad2deg(δ_tol[2]))]")
        println(io, "Δω tolerance (p.u.):          [$(Δω_tol[1]); +$(Δω_tol[2])]")
        println(io, "Synchronous frequency (Hz):   $f_syn")
        println(io, "Fault bus:                    $bus_fault")
        println(io, "Tripped circuits:             $(join(string.(circ_trip), ", "))")
        println(io, "Simulation start time (s):    $t_start_sim")
        println(io, "Simulation end time (s):      $t_end_sim")
        println(io, "Time step (s):                $t_step")
        println(io, "Fault start time (s):         $t_start_fault")
        println(io, "Fault lasts (s):              $clearing_time")
        println(io, "Fault is cleared at time (s): $t_clear_fault")
        println(io, "==============================================")
    end

    println("Input parameters successfully saved in: ", path_folder_results)

    cd(current_path_folder)
    
end

# ===================================================================================
#                   PRINT THE OPTIMIZATION MODEL IN TXT FILE
# ===================================================================================

# Function to write the AC-OPF model in a txt file
function Export_ACOPF_Model(model::Model, 
    V::OrderedDict{Int64, VariableRef}, 
    θ::OrderedDict{Int64, VariableRef}, 
    P_g::OrderedDict{Int64, VariableRef}, 
    Q_g::OrderedDict{Int64, VariableRef}, 
    eq_const_angle_sw::ConstraintRef, 
    eq_const_p_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_balance::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_s_ik::OrderedDict{Int64, ConstraintRef},
    ineq_const_s_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_diff_ang::OrderedDict{Int64, ConstraintRef},
    current_path_folder::String, 
    path_folder_results::String
    )

    cd(joinpath(path_folder_results,"ACOPF")) # Load the path folder for results

    # Open the file for writing
    open("model_summary.txt", "w") do io
        # Print the model to the file
        show(io, model)
    end

    # Desired key order to print the variables
    vector_dict_var = [V, θ, P_g, Q_g]

    open("ACOPF_model_details.txt", "w") do io

        # ------------------
        # Objective fuction
        # ------------------
        println(io, "=========")
        println(io, "Objective ")
        println(io, "=========")
        println(io, model.ext[:objective])
        println(io, "\n")

        # ---------------------------
        # Variables used in the model
        # ---------------------------
        println(io, "=========")
        println(io, "Variables")
        println(io, "=========")
        for i in eachindex(vector_dict_var)
            for (j, info) in vector_dict_var[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # --------------------
        # Equality constraint
        # --------------------
        println(io, "===============================")
        println(io, "Equality Constraint Angle Swing ")
        println(io, "===============================")
        println(io, "1: ", eq_const_angle_sw)
        println(io, "\n")

        println(io, "===================================================")
        println(io, "Equality Constraints Active Power Balance for Buses ")
        println(io, "===================================================")
        for (i, info) in eq_const_p_balance 
            println(io, "$i: ", info) 
            println(io, "")
        end
        println(io, "\n")

        println(io, "=====================================================")
        println(io, "Equality Constraints Reactive Power Balance for Buses ")
        println(io, "=====================================================")
        for (i, info) in eq_const_q_balance 
            println(io, "$i: ", info) 
            println(io, "")
        end
        println(io, "\n")

        # ----------------------
        # Inequality constraint
        # ----------------------
        println(io, "===========================================================")
        println(io, "Inequality Constraints Capacity Power Flow from Line i to k ")
        println(io, "===========================================================")
        for (i, info) in ineq_const_s_ik
            println(io, "$i: ", info)
            println(io, "")
        end
        println(io, "\n")

        println(io, "===========================================================")
        println(io, "Inequality Constraints Capacity Power Flow from Line k to i ")
        println(io, "===========================================================")
        for (i, info) in ineq_const_s_ki
            println(io, "$i: ", info)
            println(io, "")
        end
        println(io, "\n")

        println(io, "==============================================================")
        println(io, "Inequality Constraints Voltage Angle Differences between Buses ")
        println(io, "==============================================================")
        for (i, info) in ineq_const_diff_ang
            println(io, "$i: ", info)
        end
        println(io, "\n")

        # Dicts with variables
        dicts_of_vars =[V, θ, P_g, Q_g]  
        all_vars = Set(v for d in dicts_of_vars for v in values(d))

        println(io, "=========================================================")
        println(io, "Inequality Constraints Inferior Limits Decision Variables ")
        println(io, "=========================================================")
        const_lim_inf_decision_var = JuMP.all_constraints(model, VariableRef, MOI.LessThan{Float64})

        for (i, cref) in enumerate(const_lim_inf_decision_var)
            c_obj = JuMP.constraint_object(cref)  
            var = c_obj.func                      
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end
        println(io, "\n")

        println(io, "=========================================================")
        println(io, "Inequality Constraints Superior Limits Decision Variables ")
        println(io, "=========================================================")
        const_lim_sup_decision_var = JuMP.all_constraints(model, VariableRef, MOI.GreaterThan{Float64})

        for (i, cref) in enumerate(const_lim_sup_decision_var)
            c_obj = JuMP.constraint_object(cref)  
            var = c_obj.func                      
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end   
        println(io, "\n")

    end
    cd(current_path_folder)

    println("AC-OPF Model successfully saved as TXT file in: ", joinpath(path_folder_results,"ACOPF"))

end


# Function to print and save variables according to the solution of the model
function Save_Solution_Model_KCL(model::Model, 
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame,
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict,
    current_path_folder::String,
    path_folder_results::String
    )

    println("===================================================")
    println("Objective Function: € "*string(round(JuMP.value.(model.ext[:objective]), digits=2))*"")
    println("===================================================")

    P_g_optim =[JuMP.value(v) for (i, v) in P_g]   
    Q_g_optim =[JuMP.value(v) for (i, v) in Q_g]   
    S_g_optim = abs.(P_g_optim .+ 1im .* Q_g_optim) 
    V_optim   =[JuMP.value(v) for (i, v) in V]     
    θ_optim   =[JuMP.value(v) for (i, v) in θ]     

    P_ik, Q_ik, S_ik, P_ki, Q_ki, S_ki, Plosses, Qlosses, circ_loading = Calculate_AC_Power_Flow(DCIR, nCIR, V_optim, θ_optim, base_MVA)

    # =============================================================================
    #                                   Generators
    # =============================================================================
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()
    Q_g_dict = Dict{Int, Float64}()
    S_g_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
            Q_g_dict[id] = Float64(Q_g_optim[aux_count])
            S_g_dict[id] = Float64(S_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  
            Q_g_all[i] = Q_g_dict[id]  
            S_g_all[i] = S_g_dict[id]
        end
    end

    # Calculate loading
    gen_loading_p =[DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] 
    gen_loading_q = zeros(Float64, nGEN) 
    for i in 1:nGEN
        if Q_g_all[i] >= 0.0 
            gen_loading_q[i] = DGEN.g_status[i] * (Q_g_all[i] / (DGEN.qg_max[i] / base_MVA))
        else 
            gen_loading_q[i] = -DGEN.g_status[i] * (abs(Q_g_all[i]) / (abs(DGEN.qg_min[i]) / base_MVA))
        end
    end

    # =============================================================================
    #                                   Buses
    # =============================================================================
    P_g_bus = zeros(Float64, nBUS) 
    Q_g_bus = zeros(Float64, nBUS) 
    S_g_bus = zeros(Float64, nBUS) # Added for Bus apparent power
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
            Q_g_bus[i] = sum(Q_g_all[indices_bus_gen])
            S_g_bus[i] = sqrt((P_g_bus[i]^2) + (Q_g_bus[i]^2)) 
        end
    end

    all_cap =[DCIR.l_cap_1 DCIR.l_cap_2 DCIR.l_cap_3]
    circ_cap = zeros(Float64, nCIR)
    for i in 1:nCIR
        if any(!iszero, all_cap[i,:])
            index_cap = findfirst(!iszero, all_cap[i,:])
            circ_cap[i] = all_cap[i, index_cap]
        else
            circ_cap[i] = Inf
        end
    end

    # Correcting the buses labels
    bus_bus      =[reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      =[reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     =[reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       =[reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             
        v    = V_optim,                                                 
        θ    = round.(rad2deg.(θ_optim), digits=3),                   
        p    = round.((P_g_bus .* base_MVA) .- DBUS.p_d, digits=3), 
        q    = round.((Q_g_bus .* base_MVA) .- DBUS.q_d, digits=3),  
        p_g  = round.(P_g_bus .* base_MVA, digits=3),                
        q_g  = round.(Q_g_bus .* base_MVA, digits=3),
        s_g  = round.(S_g_bus .* base_MVA, digits=3), 
        p_d  = round.(DBUS.p_d, digits=3),                            
        q_d  = round.(DBUS.q_d, digits=3),                            
        p_sh = round.(DBUS.g_sh .* (V_optim.^2), digits=3),           
        q_sh = round.(DBUS.b_sh .* (V_optim.^2), digits=3)            
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        circ      = DCIR.circ,                              
        from_bus  = from_bus,                               
        to_bus    = to_bus,                                 
        p_ik      = round.(P_ik .* base_MVA, digits=3), 
        q_ik      = round.(Q_ik .* base_MVA, digits=3), 
        s_ik      = round.(S_ik .* base_MVA, digits=3), 
        p_ki      = round.(P_ki .* base_MVA, digits=3), 
        q_ki      = round.(Q_ki .* base_MVA, digits=3), 
        s_ki      = round.(S_ki .* base_MVA, digits=3), 
        p_losses  = round.(Plosses .* base_MVA, digits=3),  
        q_losses  = round.(Qlosses .* base_MVA, digits=3),  
        s_cap     = circ_cap,                               
        loading   = circ_loading                            
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id_gen    = DGEN.id,                                 
        id_bus    = gen_bus,                                 
        p_g       = round.((P_g_all .* base_MVA), digits=3), 
        q_g       = round.((Q_g_all .* base_MVA), digits=3), 
        s_g       = round.((S_g_all .* base_MVA), digits=3), 
        loading_p = round.(gen_loading_p,         digits=3), 
        loading_q = round.(gen_loading_q,         digits=3)  
    )

    Save_ResultsTXT_ACOPF(path_folder_results, RBUS, nBUS, RGEN, nGEN, RCIR, nCIR, Float64(JuMP.value.(model.ext[:objective]))) 
    Save_ResultsCSV_ACOPF(path_folder_results, RBUS, nBUS, RGEN, nGEN, RCIR, nCIR, Float64(JuMP.value.(model.ext[:objective]))) 

    cd(current_path_folder)

    return RBUS, RGEN, RCIR
end

# Save the reports of the power flow in CSV files
function Save_ResultsCSV_ACOPF(path_folder_results::String, RBUS::DataFrame, nBUS::Int64, RGEN::DataFrame, nGEN::Int64, RCIR::DataFrame, nCIR::Int64, objective::Float64)
    cd(joinpath(path_folder_results, "ACOPF\\CSV"))
  
    # Save Bus Report as CSV 
    df_buses = DataFrame(
        BUS = RBUS.bus,
        V_pu = round.(RBUS.v, digits = 4),
        Theta_deg = round.(RBUS.θ, digits = 4),
        P_MW = round.(RBUS.p, digits = 4),
        Q_MVAr = round.(RBUS.q, digits = 4),
        PG_MW = round.(RBUS.p_g, digits = 4),
        QG_MVAr = round.(RBUS.q_g, digits = 4),
        SG_MVA = round.(RBUS.s_g, digits = 4), 
        PD_MW = round.(RBUS.p_d, digits = 4),
        QD_MVAr = round.(RBUS.q_d, digits = 4),
        Psh_MW =  round.(RBUS.p_sh, digits = 4),
        Qsh_MVAr =  round.(RBUS.q_sh, digits = 4)
    )
    CSV.write("buses_report.csv", df_buses; delim=';', writeheader=true)

    # Save Generators Report as CSV
    df_generators = DataFrame(
        ID = RGEN.id_gen,
        BUS = RGEN.id_bus,
        P_MW = round.(RGEN.p_g, digits = 4),
        Q_MVAr = round.(RGEN.q_g, digits = 4),
        S_MVA = round.(RGEN.s_g, digits = 4),
        Loading_P = round.(RGEN.loading_p, digits = 4),
        Loading_Q = round.(RGEN.loading_q, digits = 4)
    )
    CSV.write("generators_report.csv", df_generators; delim=';', writeheader=true)

    # Save Circuit Report as CSV
    df_circuits = DataFrame(
        ID_CIRC = RCIR.circ,
        FROM_BUS = RCIR.from_bus,
        TO_BUS = RCIR.to_bus,
        Pik_MW = round.(RCIR.p_ik, digits = 4),
        Qik_MVAr = round.(RCIR.q_ik, digits = 4),
        Sik_MVA = round.(RCIR.s_ik, digits = 4),
        Pki_MW = round.(RCIR.p_ki, digits = 4),
        Qki_MVAr = round.(RCIR.q_ki, digits = 4),
        Ski_MVA = round.(RCIR.s_ki, digits = 4),
        Cap_MVA = round.(RCIR.s_cap, digits = 4),
        Loading = round.(RCIR.loading, digits = 4),
        Ploss_MW = round.(RCIR.p_losses, digits = 4),
        Qloss_MVAr = round.(RCIR.q_losses, digits = 4)
    )
    CSV.write("circuits_report.csv", df_circuits; delim=';', writeheader=true)

    # Save Optimization Report as CSV
    df_optimization = DataFrame(
        Metric =["Total Cost (Euros)"],
        Value = [round(objective, digits=2)]
    )
    CSV.write("optimization_report.csv", df_optimization; delim=';', writeheader=true)

    println("AC-OPF results successfully saved as CSV files in: ", joinpath(path_folder_results, "ACOPF\\CSV"))
end

# Save the reports of the power flow in TXT files
function Save_ResultsTXT_ACOPF(path_folder_results::String, RBUS::DataFrame, nBUS::Int64, RGEN::DataFrame, nGEN::Int64, RCIR::DataFrame, nCIR::Int64, objective::Float64)
    cd(joinpath(path_folder_results, "ACOPF"))

    sum_Pg = 0.0
    sum_Qg = 0.0
    sum_Sg = 0.0
    sum_Pd = 0.0
    sum_Qd = 0.0
    sum_Psh = 0.0
    sum_Qsh = 0.0
    sum_losses_P = 0.0
    sum_losses_Q = 0.0
    sum_gen_Pg = 0.0
    sum_gen_Qg = 0.0
    sum_gen_Sg = 0.0

    io = open("buses_report.txt", "w")
    @printf(io, "BUSES REPORT\n")
    @printf(io, "======================================================================================================================================== \n")
    @printf(io, "   BUS    V (pu)     O (º)     P (MW)     Q (MVAr)      PG (MW)    QG (MVAr)    SG (MVA)    PD (MW)   QD (MVAr)     Psh (MW)   Qsh (MVAr) \n")
    @printf(io, "---------------------------------------------------------------------------------------------------------------------------------------- \n")
    for i = 1:nBUS
        @printf(io, " %4d    %6.4f    %6.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n", RBUS.bus[i], RBUS.v[i], RBUS.θ[i], RBUS.p[i], RBUS.q[i], RBUS.p_g[i], RBUS.q_g[i], RBUS.s_g[i], RBUS.p_d[i], RBUS.q_d[i], RBUS.p_sh[i], RBUS.q_sh[i])
        sum_Pg += RBUS.p_g[i]
        sum_Qg += RBUS.q_g[i]
        sum_Sg += RBUS.s_g[i]
        sum_Pd += RBUS.p_d[i]
        sum_Qd += RBUS.q_d[i]
        sum_Psh += RBUS.p_sh[i]
        sum_Qsh += RBUS.q_sh[i]
    end
    @printf(io, "---------------------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                               %8.2f    %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f\n", sum_Pg, sum_Qg, sum_Sg, sum_Pd , sum_Qd, sum_Psh , sum_Qsh)
    @printf(io, "======================================================================================================================================== \n")
    @printf(io, "\n")
    close(io)  

    io = open("generators_report.txt", "w")
    @printf(io, "GENERATORS REPORT\n")
    @printf(io, "========================================================================== \n")
    @printf(io, "   ID     BUS     P_g (MW)  Q_g (MVAr)   S_g (MVA)   Loading_P   Loading_Q \n")
    @printf(io, "-------------------------------------------------------------------------- \n")
    for i = 1:nGEN
        @printf(io, " %4d   %4d    %8.2f    %8.2f    %8.2f    %8.4f    %8.4f   \n", RGEN.id_gen[i], RGEN.id_bus[i], RGEN.p_g[i], RGEN.q_g[i], RGEN.s_g[i], RGEN.loading_p[i], RGEN.loading_q[i])
        sum_gen_Pg += RGEN.p_g[i]
        sum_gen_Qg += RGEN.q_g[i]
        sum_gen_Sg += RGEN.s_g[i]
    end
    @printf(io, "-------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:         %8.4f   %8.4f   %8.4f\n", round(sum_gen_Pg, digits=3), round(sum_gen_Qg, digits=3), round(sum_gen_Sg, digits=3))
    @printf(io, "========================================================================== \n")
    close(io)
    
    io = open("circuits_report.txt", "w")
    @printf(io, "CIRCUITS REPORT\n")
    @printf(io, "=============================================================================================================================================== \n")
    @printf(io, "  CIRC    FROM    TO      Pik (MW)  Qik (MVAr)   Sik (MVA)   Pki (MW)  Qki (MVAr)    Ski (MVA)  Cap (MVA)    Loading    Ploss (MW)  Qloss (MVAr)\n")
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    for i = 1:nCIR
        @printf(io, " %4d   %4d   %4d    %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f   %8.4f     %8.4f     %8.3f     %8.3f\n", RCIR.circ[i], RCIR.from_bus[i], RCIR.to_bus[i], RCIR.p_ik[i], RCIR.q_ik[i], RCIR.s_ik[i], RCIR.p_ki[i], RCIR.q_ki[i], RCIR.s_ki[i], RCIR.s_cap[i], RCIR.loading[i], RCIR.p_losses[i], RCIR.q_losses[i])
        sum_losses_P += RCIR.p_losses[i]
        sum_losses_Q += RCIR.q_losses[i]
    end
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                                                                                                  %8.3f     %8.3f\n", sum_losses_P , sum_losses_Q)
    @printf(io, "=============================================================================================================================================== \n")
    close(io)
    
    io = open("optimization_report.txt", "w")
    @printf(io, "OBJECTIVE\n")
    @printf(io, "================================ \n")
    @printf(io, "Total cost: (Euros) %8.2f \n", round(objective, digits = 2))
    @printf(io, "================================ \n")
    close(io)

    println("AC-OPF results successfully saved as TXT files in: ", joinpath(path_folder_results, "ACOPF"))
end

# Function to save results across the whole time window (fault and post-fault)
function Manage_Dyn_Results(model::Model,
    DGEN_DYN::DataFrame,
    E_fd::OrderedDict{Int, VariableRef}, 
    δ::OrderedDict{Int, VariableRef}, 
    Pm::OrderedDict{Int, VariableRef},
    Ed::OrderedDict{Int, VariableRef},
    Eq::OrderedDict{Int, VariableRef},
    Id::OrderedDict{Int, VariableRef},
    Iq::OrderedDict{Int, VariableRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Qe_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}}, 
    V_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},  
    θ_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},  
    δ_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Δω_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Ed_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Eq_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Id_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Iq_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Te_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    E_fd_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    P_valve_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}}, # Added Governor
    P_mech_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},  # Added Governor
    δCOI_tf::OrderedDict{Int, VariableRef},
    ΔωCOI_tf::OrderedDict{Int, VariableRef},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Qe_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    V_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}}, 
    θ_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}}, 
    δ_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Δω_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Ed_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Eq_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Id_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Iq_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Te_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    E_fd_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    P_valve_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}}, # Added Governor
    P_mech_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},  # Added Governor
    δCOI_tpf::OrderedDict{Int, VariableRef},
    ΔωCOI_tpf::OrderedDict{Int, VariableRef},
    t_window_total::Vector{Float64},
    t_clear_fault::Float64,
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64},
    base_MVA::Float64,
    f_syn::Float64,
    ω_syn::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability"))

    E_fd_values =[JuMP.value(v) for (i, v) in E_fd]
    δ_values  =[JuMP.value(v) for (i, v) in δ]
    Pm_values =[JuMP.value(v) for (i, v) in Pm] .* base_MVA
    Ed_values =[JuMP.value(v) for (i, v) in Ed]
    Eq_values =[JuMP.value(v) for (i, v) in Eq]
    Id_values =[JuMP.value(v) for (i, v) in Id]
    Iq_values =[JuMP.value(v) for (i, v) in Iq]

    # ====================================
    # Save short results in TXT
    # ====================================
    open("dynamic_model_short_results.txt", "w") do io
        println(io, "==========================")
        println(io, "         E_fd[p.u.]")
        println(io, "==========================")
        for (i, info) in E_fd
            println(io, "$info: ", E_fd_values[i])
        end
        println(io, "\n")

        println(io, "==========================")
        println(io, "         δ [deg]")
        println(io, "==========================")
        for (i, info) in δ
            println(io, "$info: ", rad2deg(δ_values[i]))
        end
        println(io, "\n")

        println(io, "==========================")
        println(io, "          Pm [MW]")
        println(io, "==========================")
        for (i, info) in δ
            println(io, "$info: ", Pm_values[i])
        end
        println(io, "\n")
        
        println(io, "==========================")
        println(io, "   Initial Ed, Eq, Id, Iq")
        println(io, "==========================")
        for (i, _) in Ed
            println(io, "Gen $i -> Ed: $(Ed_values[i]), Eq: $(Eq_values[i]), Id: $(Id_values[i]), Iq: $(Iq_values[i])")
        end
        println(io, "\n")
    end

    # ====================================
    #      Pe, Qe, Pmech, Pvalve
    # ====================================
    Pe_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Pe_tf)
    Pe_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Pe_tpf)
    Pet = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Pe_tf_values)
        Pet[k] = vcat(Pe_tf_values[k], Pe_tpf_values[k])
    end

    Qe_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Qe_tf)
    Qe_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Qe_tpf)
    Qet = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Qe_tf_values)
        Qet[k] = vcat(Qe_tf_values[k], Qe_tpf_values[k])
    end

    P_mech_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in P_mech_tf)
    P_mech_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in P_mech_tpf)
    P_mecht = OrderedDict{Int, Vector{Float64}}()
    for k in keys(P_mech_tf_values)
        P_mecht[k] = vcat(P_mech_tf_values[k], P_mech_tpf_values[k])
    end

    P_valve_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in P_valve_tf)
    P_valve_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in P_valve_tpf)
    P_valvet = OrderedDict{Int, Vector{Float64}}()
    for k in keys(P_valve_tf_values)
        P_valvet[k] = vcat(P_valve_tf_values[k], P_valve_tpf_values[k])
    end

    # ====================================
    #      Ed, Eq, Id, Iq (Dynamics)
    # ====================================
    Edt = OrderedDict{Int, Vector{Float64}}()
    Eqt = OrderedDict{Int, Vector{Float64}}()
    Idt = OrderedDict{Int, Vector{Float64}}()
    Iqt = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Ed_tf)
        Edt[k] = vcat([JuMP.value(v) for (_, v) in Ed_tf[k]],[JuMP.value(v) for (_, v) in Ed_tpf[k]])
        Eqt[k] = vcat([JuMP.value(v) for (_, v) in Eq_tf[k]],[JuMP.value(v) for (_, v) in Eq_tpf[k]])
        Idt[k] = vcat([JuMP.value(v) for (_, v) in Id_tf[k]],[JuMP.value(v) for (_, v) in Id_tpf[k]])
        Iqt[k] = vcat([JuMP.value(v) for (_, v) in Iq_tf[k]],[JuMP.value(v) for (_, v) in Iq_tpf[k]])
    end

    # ====================================
    #              V and θ (Nodes)
    # ====================================
    V_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in V_tf)
    V_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in V_tpf)
    Vt = OrderedDict{Int, Vector{Float64}}()
    for k in keys(V_tf_values)
        Vt[k] = vcat(V_tf_values[k], V_tpf_values[k])
    end

    θ_tf_values  = OrderedDict(i =>[rad2deg(JuMP.value(v)) for (_, v) in inner] for (i, inner) in θ_tf)
    θ_tpf_values = OrderedDict(i =>[rad2deg(JuMP.value(v)) for (_, v) in inner] for (i, inner) in θ_tpf)
    θt = OrderedDict{Int, Vector{Float64}}()
    for k in keys(θ_tf_values)
        θt[k] = vcat(θ_tf_values[k], θ_tpf_values[k])
    end

    # ====================================
    #                 δ
    # ====================================
    δ_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in δ_tf)   
    δ_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in δ_tpf) 

    δCOI_tf_values  =[JuMP.value(v) for (i, v) in δCOI_tf]  
    δCOI_tpf_values =[JuMP.value(v) for (i, v) in δCOI_tpf] 

    δt = OrderedDict{Int, Vector{Float64}}()      
    δ_COIt = OrderedDict{Int, Vector{Float64}}()  

    δCOIt = vcat(δCOI_tf_values, δCOI_tpf_values) 

    for k in keys(δ_tf_values)
        δt[k]    = vcat(δ_tf_values[k], δ_tpf_values[k])
        δ_COIt[k] = δt[k] .- δCOIt
    end

    # ====================================
    #                 Δω
    # ====================================
    Δω_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in Δω_tf)  
    Δω_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in Δω_tpf) 

    ΔωCOI_tf_values  =[JuMP.value(v) for (i, v) in ΔωCOI_tf]  
    ΔωCOI_tpf_values =[JuMP.value(v) for (i, v) in ΔωCOI_tpf] 

    Δωt = OrderedDict{Int, Vector{Float64}}()     
    Δω_COIt = OrderedDict{Int, Vector{Float64}}() 

    ΔωCOIt = vcat(ΔωCOI_tf_values, ΔωCOI_tpf_values) 

    for k in keys(Δω_tf_values)
        Δωt[k]    = vcat(Δω_tf_values[k], Δω_tpf_values[k])
        Δω_COIt[k] = Δωt[k] .- ΔωCOIt
    end

    # ====================================
    #              RoCoF
    # ====================================
    time_RoCoF =[]
    RoCoF = OrderedDict{Int, Vector{Float64}}()  
    for (i, info) in Δωt
        time_RoCoF, RoCoF[i] = Calculate_Derivative_CFD(t_window_total, f_syn .* (1.0 .+ info))
    end
    time_RoCoF, RoCoFCOI = Calculate_Derivative_CFD(t_window_total, f_syn .* (1.0 .+ ΔωCOIt))

    # ====================================
    #          Kinetic Energy
    # ====================================
    Vke = OrderedDict{Int, Vector{Float64}}()  
    for (i, info) in Δω_COIt
        Vke[i] = DGEN_DYN.H[i] .* ω_syn .* (info .^2)
    end

    # ====================================
    #        Accelerating Power
    # ====================================
    H_total   = sum(DGEN_DYN.H[Int.(keys(Pet))])

    Pacc      = OrderedDict{Int, Vector{Float64}}() 
    PaccCOI   = OrderedDict{Int, Vector{Float64}}() 
    Pacc_COI  = OrderedDict{Int, Vector{Float64}}() 
    aux_count = 0
    for (i, info) in Pet
        Pacc[i] = (P_mecht[i] .- info) # Calculated dynamically via Governor Mechanical Power

        aux_count += 1
        if aux_count == 1
            PaccCOI[1] = (P_mecht[i] .- info)
        else
            PaccCOI[1] = PaccCOI[1] .+ (P_mecht[i] .- info)
        end
    end
    for (i, info) in Pacc
        Pacc_COI[i] = (info .- ((DGEN_DYN.H[i] .* PaccCOI[1]) ./ H_total)) ./ base_MVA
    end
    
    # ====================================
    #          Potential Energy
    # ====================================
    Vpe = OrderedDict{Int, Vector{Float64}}() 
    for k in keys(Pacc_COI)
        integral =[]
        for (i, _) in enumerate(Pacc_COI[k])
            if i == 1
                push!(integral, 0.0)
            else
                push!(integral, -trapz(δ_COIt[k][1:i], Pacc_COI[k][1:i]))
            end
        end
        Vpe[k] = Float64.(integral) 
    end
    for (i, info) in Pacc_COI
        Pacc_COI[i] = info .* base_MVA
    end

    # Field voltage time series (AVR dynamics)
    E_fd_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in E_fd_tf)
    E_fd_tpf_values = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] for (i, inner) in E_fd_tpf)
    E_fdt = OrderedDict{Int, Vector{Float64}}()
    for k in keys(E_fd_tf_values)
        E_fdt[k] = vcat(E_fd_tf_values[k], E_fd_tpf_values[k])
    end

    # Electrical torque time series
    Te_tf_values  = OrderedDict(i =>[JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Te_tf)
    Te_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Te_tpf)
    Tet = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Te_tf_values)
        Tet[k] = vcat(Te_tf_values[k], Te_tpf_values[k])
    end

    # ====================================
    #           SAVE FIGURES
    # ====================================
    Save_Dyn_Results_Plots(t_window_total, Pm_values,
    Pet, Qet, P_mecht, P_valvet, Vt, θt, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt,
    Edt, Eqt, Idt, Iqt, Tet, E_fdt,
    time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe,
    base_MVA, f_syn, δ_tol, Δω_tol, t_clear_fault,
    current_path_folder, path_folder_results)

    # ====================================
    #           SAVE CSV FILES
    # ====================================
    Save_Dyn_Results_CSV(t_window_total, Pm_values, E_fd_values,
        Pet, Qet, P_mecht, P_valvet, Vt, θt, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt,
        Edt, Eqt, Idt, Iqt, Tet, E_fdt,
        time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe,
        base_MVA, f_syn, current_path_folder, path_folder_results
    )
    cd(current_path_folder)
end

# Function to save the transient stability results in CSV files
function Save_Dyn_Results_CSV(t_window_total::Vector{Float64},
    Pm_values::Vector{Float64},
    E_fd_values::Vector{Float64},
    Pet::OrderedDict{Int, Vector{Float64}},
    Qet::OrderedDict{Int, Vector{Float64}}, 
    P_mecht::OrderedDict{Int, Vector{Float64}}, 
    P_valvet::OrderedDict{Int, Vector{Float64}}, 
    Vt::OrderedDict{Int, Vector{Float64}},  
    θt::OrderedDict{Int, Vector{Float64}},  
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    Edt::OrderedDict{Int, Vector{Float64}},
    Eqt::OrderedDict{Int, Vector{Float64}},
    Idt::OrderedDict{Int, Vector{Float64}},
    Iqt::OrderedDict{Int, Vector{Float64}},
    Tet::OrderedDict{Int, Vector{Float64}},
    E_fdt::OrderedDict{Int, Vector{Float64}},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability\\CSV"))
    
    # Names of the variables
    gen_names =["G$(i)" for (i, inner) in Pet]
    bus_names =["B$(i)" for (i, inner) in Vt]
    delta_names =["G$(i)" for (i, inner) in δt]
    push!(delta_names, "COI")
    omega_names =["G$(i)" for (i, inner) in Δωt]
    push!(omega_names, "COI")

    # Constant Generator Variables
    gen_ids = collect(keys(Pet))
    df_gen_consts = DataFrame(
        GEN_ID = gen_ids,
        E_fd_pu = round.(E_fd_values, digits=6),
        Pm_MW = round.(Pm_values, digits=4)
    )

    # Electrical Active Power (Pe)
    Pe_matrix = zeros(Float64, length(t_window_total), length(Pet))
    aux_count = 0; for (gen_id, values) in Pet aux_count += 1; Pe_matrix[:,aux_count] = values end
    df_Pe = DataFrame(hcat(t_window_total, Pe_matrix), vcat("t", gen_names))

    # Electrical Reactive Power (Qe)
    Qe_matrix = zeros(Float64, length(t_window_total), length(Qet))
    aux_count = 0; for (gen_id, values) in Qet aux_count += 1; Qe_matrix[:,aux_count] = values end
    df_Qe = DataFrame(hcat(t_window_total, Qe_matrix), vcat("t", gen_names))

    # Mechanical Power (Pmech)
    P_mech_matrix = zeros(Float64, length(t_window_total), length(P_mecht))
    aux_count = 0; for (gen_id, values) in P_mecht aux_count += 1; P_mech_matrix[:,aux_count] = values end
    df_P_mech = DataFrame(hcat(t_window_total, P_mech_matrix), vcat("t", gen_names))

    # Valve Position (Pvalve)
    P_valve_matrix = zeros(Float64, length(t_window_total), length(P_valvet))
    aux_count = 0; for (gen_id, values) in P_valvet aux_count += 1; P_valve_matrix[:,aux_count] = values end
    df_P_valve = DataFrame(hcat(t_window_total, P_valve_matrix), vcat("t", gen_names))

    # 4th Order DQ variables
    Ed_matrix = zeros(Float64, length(t_window_total), length(Edt))
    Eq_matrix = zeros(Float64, length(t_window_total), length(Eqt))
    Id_matrix = zeros(Float64, length(t_window_total), length(Idt))
    Iq_matrix = zeros(Float64, length(t_window_total), length(Iqt))
    aux_count = 0
    for (gen_id, _) in Edt
        aux_count += 1
        Ed_matrix[:,aux_count] = Edt[gen_id]; Eq_matrix[:,aux_count] = Eqt[gen_id]
        Id_matrix[:,aux_count] = Idt[gen_id]; Iq_matrix[:,aux_count] = Iqt[gen_id]
    end
    df_Ed = DataFrame(hcat(t_window_total, Ed_matrix), vcat("t", gen_names))
    df_Eq = DataFrame(hcat(t_window_total, Eq_matrix), vcat("t", gen_names))
    df_Id = DataFrame(hcat(t_window_total, Id_matrix), vcat("t", gen_names))
    df_Iq = DataFrame(hcat(t_window_total, Iq_matrix), vcat("t", gen_names))

    # Bus Voltage Magnitude & Angle
    V_matrix = zeros(Float64, length(t_window_total), length(Vt))
    aux_count = 0; for (bus_id, values) in Vt aux_count += 1; V_matrix[:,aux_count] = values end
    df_V = DataFrame(hcat(t_window_total, V_matrix), vcat("t", bus_names))

    θ_matrix = zeros(Float64, length(t_window_total), length(θt))
    aux_count = 0; for (bus_id, values) in θt aux_count += 1; θ_matrix[:,aux_count] = values end
    df_θ = DataFrame(hcat(t_window_total, θ_matrix), vcat("t", bus_names))

    # Accelerating Power
    Pacc_matrix = zeros(Float64, length(t_window_total), length(Pacc)+1)
    aux_count = 0; for (gen_id, values) in Pacc aux_count += 1; Pacc_matrix[:,aux_count] = values end
    Pacc_matrix[:,end] = PaccCOI[1]
    df_Pacc = DataFrame(hcat(t_window_total, Pacc_matrix), vcat("t", delta_names))

    Pacc_COI_matrix = zeros(Float64, length(t_window_total), length(Pacc_COI))
    aux_count = 0; for (gen_id, values) in Pacc_COI aux_count += 1; Pacc_COI_matrix[:,aux_count] = values end
    df_Pacc_COI = DataFrame(hcat(t_window_total, Pacc_COI_matrix), vcat("t", gen_names))

    # δ and Δω
    δ_matrix = zeros(Float64, length(t_window_total), length(δt)+1)
    aux_count = 0; for (gen_id, values) in δt aux_count += 1; δ_matrix[:,aux_count] = rad2deg.(values) end
    δ_matrix[:,end] = rad2deg.(δCOIt)
    df_δ = DataFrame(hcat(t_window_total, δ_matrix), vcat("t", delta_names))

    δ_COI_matrix = zeros(Float64, length(t_window_total), length(δ_COIt))
    aux_count = 0; for (gen_id, values) in δ_COIt aux_count += 1; δ_COI_matrix[:,aux_count] = rad2deg.(values) end
    df_δ_COI = DataFrame(hcat(t_window_total, δ_COI_matrix), vcat("t", gen_names))
    
    Δω_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0; for (gen_id, values) in Δωt aux_count += 1; Δω_matrix[:,aux_count] = values end
    Δω_matrix[:,end] = ΔωCOIt
    df_Δω = DataFrame(hcat(t_window_total, Δω_matrix), vcat("t", omega_names))

    Δω_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0; for (gen_id, values) in Δω_COIt aux_count += 1; Δω_COI_matrix[:,aux_count] = values end
    df_Δω_COI = DataFrame(hcat(t_window_total, Δω_COI_matrix), vcat("t", gen_names))

    # Frequency
    f_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0; for (gen_id, values) in Δωt aux_count += 1; f_matrix[:,aux_count] = f_syn .* (1.0 .+ values) end
    f_matrix[:,end] = f_syn .* (1.0 .+ ΔωCOIt)
    df_f = DataFrame(hcat(t_window_total, f_matrix), vcat("t", omega_names))

    f_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0; for (gen_id, values) in Δω_COIt aux_count += 1; f_COI_matrix[:,aux_count] = f_syn .* (1.0 .+ values) end
    df_f_COI = DataFrame(hcat(t_window_total, f_COI_matrix), vcat("t", gen_names))
    
    # RoCoF, Vke, Vpe, E_fd, Te
    RoCoF_matrix = zeros(Float64, length(time_RoCoF), length(RoCoF)+1)
    aux_count = 0; for (gen_id, values) in RoCoF aux_count += 1; RoCoF_matrix[:,aux_count] = values end
    RoCoF_matrix[:,end] = RoCoFCOI
    df_RoCoF = DataFrame(hcat(time_RoCoF, RoCoF_matrix), vcat("t", omega_names))

    Vke_matrix = zeros(Float64, length(t_window_total), length(Vke))
    aux_count = 0; for (gen_id, values) in Vke aux_count += 1; Vke_matrix[:,aux_count] = values end
    df_Vke = DataFrame(hcat(t_window_total, Vke_matrix), vcat("t", gen_names))

    Vpe_matrix = zeros(Float64, length(t_window_total), length(Vpe))
    aux_count = 0; for (gen_id, values) in Vpe aux_count += 1; Vpe_matrix[:,aux_count] = values end
    df_Vpe = DataFrame(hcat(t_window_total, Vpe_matrix), vcat("t", gen_names))

    E_fd_matrix = zeros(Float64, length(t_window_total), length(E_fdt))
    aux_count = 0; for (gen_id, values) in E_fdt aux_count += 1; E_fd_matrix[:, aux_count] = values end
    df_E_fd = DataFrame(hcat(t_window_total, E_fd_matrix), vcat("t", gen_names))

    Te_matrix = zeros(Float64, length(t_window_total), length(Tet))
    aux_count = 0; for (gen_id, values) in Tet aux_count += 1; Te_matrix[:, aux_count] = values end
    df_Te = DataFrame(hcat(t_window_total, Te_matrix), vcat("t", gen_names))   

    # CSV saving
    CSV.write("generator_constants.csv", df_gen_consts; delim=';')
    CSV.write("electrical_active_power.csv", df_Pe; delim=';')
    CSV.write("electrical_reactive_power.csv", df_Qe; delim=';') 
    CSV.write("mechanical_power.csv", df_P_mech; delim=';') 
    CSV.write("valve_position.csv", df_P_valve; delim=';') 
    
    CSV.write("bus_voltage_mag.csv", df_V; delim=';')             
    CSV.write("bus_voltage_angle.csv", df_θ; delim=';')           
    
    CSV.write("transient_voltage_d.csv", df_Ed; delim=';')
    CSV.write("transient_voltage_q.csv", df_Eq; delim=';')
    CSV.write("stator_current_d.csv", df_Id; delim=';')
    CSV.write("stator_current_q.csv", df_Iq; delim=';')
    
    CSV.write("accelerating_power.csv", df_Pacc; delim=';')
    CSV.write("accelerating_power_COI.csv", df_Pacc_COI; delim=';')
    CSV.write("angle_abs.csv", df_δ; delim=';')
    CSV.write("angle_rel_COI.csv", df_δ_COI; delim=';')
    CSV.write("speed_dev.csv", df_Δω; delim=';')
    CSV.write("speed_dev_rel_COI.csv", df_Δω_COI; delim=';')
    CSV.write("frequency.csv", df_f; delim=';')
    CSV.write("frequency_rel_COI.csv", df_f_COI; delim=';')
    CSV.write("RoCoF.csv", df_RoCoF; delim=';')
    CSV.write("kinetic_energy.csv", df_Vke; delim=';')
    CSV.write("potential_energy.csv", df_Vpe; delim=';')
    CSV.write("field_voltage.csv", df_E_fd; delim=';')
    CSV.write("electrical_torque.csv", df_Te; delim=';')

    cd(current_path_folder)
    println("Results of the dynamic model successfully saved as CSV files in: ", joinpath(path_folder_results, "Transient_Stability\\CSV"))
end

# Function to save the transient stability results in figures
function Save_Dyn_Results_Plots(t_window_total::Vector{Float64},
    Pm::Vector{Float64},
    Pet::OrderedDict{Int, Vector{Float64}},
    Qet::OrderedDict{Int, Vector{Float64}},
    P_mecht::OrderedDict{Int, Vector{Float64}}, 
    P_valvet::OrderedDict{Int, Vector{Float64}}, 
    Vt::OrderedDict{Int, Vector{Float64}},
    θt::OrderedDict{Int, Vector{Float64}},
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    Edt::OrderedDict{Int, Vector{Float64}},
    Eqt::OrderedDict{Int, Vector{Float64}},
    Idt::OrderedDict{Int, Vector{Float64}},
    Iqt::OrderedDict{Int, Vector{Float64}},
    Tet::OrderedDict{Int, Vector{Float64}},
    E_fdt::OrderedDict{Int, Vector{Float64}},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    δ_tol::Tuple{Float64, Float64},
    Δω_tol::Tuple{Float64, Float64},
    t_clear_fault::Float64,
    current_path_folder::String,
    path_folder_results::String
    )
    cd(joinpath(path_folder_results, "Transient_Stability\\Figures"))

    index_t_ct = findfirst(x -> x >= t_clear_fault,  t_window_total)

    plot_Pe = plot()
    for (gen_id, values) in Pet plot!(plot_Pe, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Pe = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"P_e \,\, \left(MW\right)", title="Electrical Active Power", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Qe = plot()
    for (gen_id, values) in Qet plot!(plot_Qe, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Qe = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"Q_e \,\, \left(MVAr\right)", title="Electrical Reactive Power", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Qe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Pmech = plot()
    for (gen_id, values) in P_mecht plot!(plot_Pmech, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Pmech = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"P_m \,\, \left(MW\right)", title="Mechanical Power", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pmech = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Pvalve = plot()
    for (gen_id, values) in P_valvet plot!(plot_Pvalve, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Pvalve = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"P_{valve} \,\, \left(p.u.\right)", title="Valve Position", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pvalve = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # 4th Order State Variables
    plot_Ed = plot()
    for (gen_id, values) in Edt plot!(plot_Ed, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Ed = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"E'_d \,\, \left(p.u.\right)", title="Transient d-axis Voltage", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Ed = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Eq = plot()
    for (gen_id, values) in Eqt plot!(plot_Eq, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Eq = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"E'_q \,\, \left(p.u.\right)", title="Transient q-axis Voltage", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Eq = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Id = plot()
    for (gen_id, values) in Idt plot!(plot_Id, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Id = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"I_d \,\, \left(p.u.\right)", title="Stator d-axis Current", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Id = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Iq = plot()
    for (gen_id, values) in Iqt plot!(plot_Iq, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Iq = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"I_q \,\, \left(p.u.\right)", title="Stator q-axis Current", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Iq = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_V = plot()
    for (bus_id, values) in Vt plot!(plot_V, t_window_total, values, lw=2, label="Bus $bus_id", ls=:solid) end
    plot_V = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"V \,\, \left(p.u.\right)", title="Bus Voltage Magnitude", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(10), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_V = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_θ = plot()
    for (bus_id, values) in θt plot!(plot_θ, t_window_total, values, lw=2, label="Bus $bus_id", ls=:solid) end
    plot_θ = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\theta \,\, \left(deg\right)", title="Bus Voltage Angle", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(10), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_θ = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # Rotor angles
    plot_δ = plot()
    for (gen_id, values) in δt plot!(plot_δ, t_window_total, rad2deg.(values), lw = 3, label = "G$gen_id", ls=:solid) end
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt), label="COI", lw=2, ls=:solid, lc=:black)
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt .+ δ_tol[1]) , label="- tol", lw=2, ls=:dash,  lc=:red)
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt .+ δ_tol[2]) , label="+ tol", lw=2, ls=:dash,  lc=:red)
    plot_δ = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\delta \,\, \left(^o\right)", title="Rotor Angles", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_δ = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_δ_COI = plot()
    for (gen_id, values) in δ_COIt plot!(plot_δ_COI, t_window_total, rad2deg.(values), lw = 3, label = "G$gen_id", ls=:solid) end
    plot_δ_COI = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\delta _{COI} \,\, \left(^o\right)", title="Rotor Angles (ref. COI)", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_δ_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # Speed deviation
    plot_Δω = plot()
    for (gen_id, values) in Δωt plot!(plot_Δω, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Δω = plot!(t_window_total, ΔωCOIt, label="COI", lw=2, ls=:solid, lc=:black)
    plot_Δω = plot!(t_window_total, ΔωCOIt .+ Δω_tol[1], label="- tol", lw=2, ls=:dash,  lc=:red)
    plot_Δω = plot!(t_window_total, ΔωCOIt .+ Δω_tol[2], label="+ tol", lw=2, ls=:dash,  lc=:red)
    plot_Δω = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\Delta \omega \,\, \left(p.u.\right)", title="Speed Deviation", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Δω_COI = plot()
    for (gen_id, values) in Δω_COIt plot!(plot_Δω_COI, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Δω_COI = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\Delta \omega _{COI} \,\, \left(p.u.\right)", title="Speed Deviation (ref. COI)", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # Rotor speed
    plot_Δω_1 = plot()
    for (gen_id, values) in Δωt plot!(plot_Δω_1, t_window_total, 1 .+ values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Δω_1 = plot!(t_window_total, 1 .+ (ΔωCOIt), label="COI", lw=2, ls=:solid, lc=:black)
    plot_Δω_1 = plot!(t_window_total, 1 .+ (ΔωCOIt .+ Δω_tol[1]), label="- tol", lw=2, ls=:dash,  lc=:red)
    plot_Δω_1 = plot!(t_window_total, 1 .+ (ΔωCOIt .+ Δω_tol[2]), label="+ tol", lw=2, ls=:dash,  lc=:red)
    plot_Δω_1 = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\omega \,\, \left(p.u.\right)", title="Rotor Speed", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Δω_COI_1 = plot()
    for (gen_id, values) in Δω_COIt plot!(plot_Δω_COI_1, t_window_total, 1 .+ values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Δω_COI_1 = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"\omega _{COI} \,\, \left(p.u.\right)", title="Rotor Speed (ref. COI)", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_COI_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # Frequency
    plot_f = plot()
    for (gen_id, values) in Δωt plot!(plot_f, t_window_total, f_syn .* (1 .+ values), lw = 3, label = "G$gen_id", ls=:solid) end
    plot_f = plot!(t_window_total, f_syn .* (1 .+ (ΔωCOIt)), label="COI", lw=2, ls=:solid, lc=:black)
    plot_f = plot!(t_window_total, f_syn .* (1 .+ (ΔωCOIt .+ Δω_tol[1])), label="- tol", lw=2, ls=:dash,  lc=:red)
    plot_f = plot!(t_window_total, f_syn .* (1 .+ (ΔωCOIt .+ Δω_tol[2])), label="+ tol", lw=2, ls=:dash,  lc=:red)
    plot_f = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"f \,\, \left(Hz\right)", title="Frequency", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_f = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_f_COI = plot()
    for (gen_id, values) in Δω_COIt plot!(plot_f_COI, t_window_total, f_syn .* (1 .+ values), lw = 3, label = "G$gen_id", ls=:solid) end
    plot_f_COI = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"f \,\, \left(Hz\right)", title="Frequency (ref. COI)", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_f_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # RoCoF
    plot_RoCoF = plot()
    for (gen_id, values) in RoCoF plot!(plot_RoCoF, time_RoCoF, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_RoCoF = plot!(time_RoCoF, RoCoFCOI, label="COI", lw=2, ls=:solid, lc=:black)
    plot_RoCoF = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"f \,\, \left(Hz / sec\right)", title="RoCoF", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_RoCoF = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # δ vs ω
    plot_δω = plot()
    for (gen_id, values) in δ_COIt plot!(plot_δω, rad2deg.(values), 1.0 .+ Δω_COIt[gen_id], lw = 3, label = "G$gen_id", ls=:solid) end
    plot_δω = plot!(xlabel=L"\delta \,\, \left(^o\right)", ylabel=L"\Delta \omega \,\, \left(p.u.\right)", title="Rotor speed vs Angle", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    for (gen_id, values) in δ_COIt scatter!(plot_δω, [rad2deg.(values[index_t_ct])],[1.0 .+ Δω_COIt[gen_id][index_t_ct]], shape=:circle, msize= 7, color=:black, label="") end

    # δ vs Pe and Pmech (Updated with Gov Dynamic power)
    plot_δ_COIPePm = plot()
    for (gen_id, values) in δ_COIt
        plot!(plot_δ_COIPePm, rad2deg.(values), Pet[gen_id], lw = 3, label = "Pe G$gen_id", ls=:solid)
        plot!(plot_δ_COIPePm, rad2deg.(values), P_mecht[gen_id], lw = 3, label = "Pmech G$gen_id", ls=:dash)
    end
    plot_δ_COIPePm = plot!(xlabel=L"\delta \,\, \left(^o\right)", ylabel=L"P \,\, \left(MW\right)", title="P vs Angle", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)

    # Accelerating Power
    plot_Pacc = plot()
    for (gen_id, values) in Pacc plot!(plot_Pacc, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Pacc = plot!(t_window_total, PaccCOI[1], label="COI", lw=2, ls=:solid, lc=:black)
    plot_Pacc = plot!(xlabel=L"t \,\, \left(s\right)", ylabel=L"P \_{acc} \,\, \left(MW\right)", title="Accelerating Power", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pacc = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Pacc_COI = plot()
    for (gen_id, values) in Pacc_COI plot!(plot_Pacc_COI, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Pacc_COI = plot!(xlabel=L"t \,\, \left(s\right)", ylabel=L"P \_{acc\_COI} \,\, \left(MW\right)", title="Accelerating Power (ref. COI)", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pacc_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # Energies
    plot_Vke = plot()
    for (gen_id, values) in Vke plot!(plot_Vke, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Vke = plot!(xlabel=L"t \,\, \left(s\right)", ylabel=L"V \_{ke} \,\, \left(p.u. \cdot rad\right)", title="Kinetic Energy vs Time", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Vke = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Vpe = plot()
    for (gen_id, values) in Vpe plot!(plot_Vpe, t_window_total, values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Vpe = plot!(xlabel=L"t \,\, \left(s\right)", ylabel=L"V \_{pe} \,\, \left(p.u. \cdot rad\right)", title="Potential Energy vs Time", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Vpe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Ve = plot()
    for (gen_id, values) in Vpe plot!(plot_Ve, t_window_total, Vke[gen_id] .+ values, lw = 3, label = "G$gen_id", ls=:solid) end
    plot_Ve = plot!(xlabel=L"t \,\, \left(s\right)", ylabel=L"V \_{e} \,\, \left(p.u. \cdot rad\right)", title="Total Energy vs Time", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Ve = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Efd = plot()
    for (gen_id, values) in E_fdt plot!(plot_Efd, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Efd = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"E_{fd} \,\, \left(p.u.\right)", title="Field Voltage", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Efd = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    plot_Te = plot()
    for (gen_id, values) in Tet plot!(plot_Te, t_window_total, values, lw=3, label="G$gen_id", ls=:solid) end
    plot_Te = plot!(xlabel=L"t \, \left(s\right)", ylabel=L"T_e \,\, \left(MW\right)", title="Electrical Torque", size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm, gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Te = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    # ====================================
    #            Save Plots
    # ====================================
    savefig(plot_Pe,        "Electrical Active Power vs time.png")
    savefig(plot_Qe,        "Electrical Reactive Power vs time.png")
    savefig(plot_Pmech,     "Mechanical Power vs time.png")
    savefig(plot_Pvalve,    "Valve Position vs time.png")
    savefig(plot_V,         "Bus Voltage Mag vs time.png")
    savefig(plot_θ,         "Bus Voltage Angle vs time.png")
    
    savefig(plot_Ed,        "Ed vs time.png")
    savefig(plot_Eq,        "Eq vs time.png")
    savefig(plot_Id,        "Id vs time.png")
    savefig(plot_Iq,        "Iq vs time.png")

    savefig(plot_δ,         "Delta vs time.png")
    savefig(plot_δ_COI,     "Delta (ref COI) vs time.png")
    savefig(plot_Δω,        "Speed Deviation vs time.png")
    savefig(plot_Δω_COI,    "Speed Deviation (ref COI) vs time.png")
    savefig(plot_Δω_1,      "Rotor speed vs time.png")
    savefig(plot_Δω_COI_1,  "Rotor speed (ref COI) vs time.png")
    savefig(plot_f,         "Frequency vs time.png")
    savefig(plot_f_COI,     "Frequency (ref COI) vs time.png")
    savefig(plot_RoCoF,     "RoCoF vs time.png")
    savefig(plot_δω,        "Rotor speed vs delta.png")
    savefig(plot_δ_COIPePm, "Power vs delta.png")
    savefig(plot_Pacc,      "Accelerating power vs time.png")
    savefig(plot_Pacc_COI,  "Accelerating power (ref COI) vs time.png")
    savefig(plot_Vke,       "Kinetic Energy vs time.png")
    savefig(plot_Vpe,       "Potential Energy vs time.png")
    savefig(plot_Ve,        "Total Energy vs time.png")
    savefig(plot_Efd,       "Field Voltage vs time.png")
    savefig(plot_Te,        "Electrical Torque vs time.png")

    cd(current_path_folder)
    println("Figures of the dynamic model successfully saved as PNG files in: ", joinpath(path_folder_results, "Transient_Stability\\Figures"))
end

# Function to write the Dynamic model to a txt file
function Export_Dyn_Model(model::Model, 
    E_fd::OrderedDict{Int, VariableRef}, δ::OrderedDict{Int, VariableRef}, Pm::OrderedDict{Int, VariableRef}, Ed::OrderedDict{Int, VariableRef}, Eq::OrderedDict{Int, VariableRef}, Id::OrderedDict{Int, VariableRef}, Iq::OrderedDict{Int, VariableRef}, V_ref::OrderedDict{Int, VariableRef}, δCOI::VariableRef,
    eq_const_Ed_init::OrderedDict{Int, ConstraintRef}, eq_const_Eq_init::OrderedDict{Int, ConstraintRef}, eq_const_Vd_init::OrderedDict{Int, ConstraintRef}, eq_const_Vq_init::OrderedDict{Int, ConstraintRef}, eq_const_Pm::OrderedDict{Int, ConstraintRef}, eq_const_Efd_init::OrderedDict{Int, ConstraintRef}, eq_const_P_init::OrderedDict{Int, ConstraintRef}, eq_const_Q_init::OrderedDict{Int, ConstraintRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Qe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, V_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, θ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Δω_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Ed_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Eq_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Id_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Iq_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Te_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, E_fd_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, P_valve_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, P_mech_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δCOI_tf::OrderedDict{Int, VariableRef}, ΔωCOI_tf::OrderedDict{Int, VariableRef},
    eq_const_δCOI_tf::OrderedDict{Int, ConstraintRef}, eq_const_ΔωCOI_tf::OrderedDict{Int, ConstraintRef}, eq_const_Pe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Pbalance_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qbalance_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vd_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vq_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_δ_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Δω_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Ed_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Eq_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_E_fd_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_valve_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_mech_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_δ_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_Δω_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Qe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, V_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, θ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Δω_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Ed_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Eq_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Id_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Iq_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Te_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, E_fd_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, P_valve_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, P_mech_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δCOI_tpf::OrderedDict{Int64, VariableRef}, ΔωCOI_tpf::OrderedDict{Int64, VariableRef},
    eq_const_δCOI_tpf::OrderedDict{Int, ConstraintRef}, eq_const_ΔωCOI_tpf::OrderedDict{Int, ConstraintRef}, eq_const_Pe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Pbalance_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qbalance_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vd_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vq_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_δ_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Δω_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Ed_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Eq_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_E_fd_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_valve_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_mech_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_δ_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_Δω_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    current_path_folder::String, path_folder_results::String
    )

    cd(joinpath(path_folder_results,"Transient_Stability"))

    vector_dict_var_pref      =[E_fd, δ, Pm, Ed, Eq, Id, Iq, V_ref]
    vector_dict_var_fault     =[Pe_tf, Qe_tf, V_tf, θ_tf, δ_tf, Δω_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf, E_fd_tf, P_valve_tf, P_mech_tf]
    vector_dict_var_fault_COI =[δCOI_tf, ΔωCOI_tf]
    vector_dict_var_postf     =[Pe_tpf, Qe_tpf, V_tpf, θ_tpf, δ_tpf, Δω_tpf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf, E_fd_tpf, P_valve_tpf, P_mech_tpf]
    vector_dict_var_postf_COI =[δCOI_tpf, ΔωCOI_tpf]

    open("dynamic_model_details.txt", "w") do io

        println(io, "=================================")
        println(io, "Variables in the Pre-Fault Period")
        println(io, "=================================")
        for i in eachindex(vector_dict_var_pref)
            for (j, info) in vector_dict_var_pref[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "=============================")
        println(io, "Variables in the Fault Period")
        println(io, "=============================")
        for i in eachindex(vector_dict_var_fault)
            for (j, infoj) in vector_dict_var_fault[i]
                for (k, infok) in infoj
                    println(io, "$k: ", infok)
                end
            end
        end
        for i in eachindex(vector_dict_var_fault_COI)
            for (j, info) in vector_dict_var_fault_COI[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "==================================")
        println(io, "Variables in the Post-Fault Period")
        println(io, "==================================")
        for i in eachindex(vector_dict_var_postf)
            for (j, infoj) in vector_dict_var_postf[i]
                for (k, infok) in infoj
                    println(io, "$k: ", infok)
                end
            end
        end
        for i in eachindex(vector_dict_var_postf_COI)
            for (j, info) in vector_dict_var_postf_COI[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "==========================================================")
        println(io, "Equality Constraints Pre-Fault Stator & EMF Initialization")
        println(io, "==========================================================")
        for (i, info) in eq_const_Vd_init println(io, "Gen $i (Vd_init): ", info) end
        for (i, info) in eq_const_Vq_init println(io, "Gen $i (Vq_init): ", info) end
        for (i, info) in eq_const_Ed_init println(io, "Gen $i (Ed_init): ", info) end
        for (i, info) in eq_const_Eq_init println(io, "Gen $i (Eq_init): ", info) end
        for (i, info) in eq_const_Pm println(io, "Gen $i (Pm_init): ", info) end
        println(io, "\n")

        println(io, "=====================================")
        println(io, "Equality Constraints Angle of the COI")
        println(io, "=====================================")
        for (i, info) in eq_const_δCOI_tf println(io, "$i (Fault): ", info) end
        for (i, info) in eq_const_δCOI_tpf println(io, "$i (Post-F): ", info) end
        println(io, "\n")

        println(io, "================================================")
        println(io, "Equality Constraints Speed Deviation of the COI ")
        println(io, "================================================")
        for (i, info) in eq_const_ΔωCOI_tf println(io, "$i (Fault): ", info) end
        for (i, info) in eq_const_ΔωCOI_tpf println(io, "$i (Post-F): ", info) end
        println(io, "\n")

        println(io, "=========================================================================")
        println(io, "Equality Constraints Electrical Power, Stator Alg. and Nodal Balances (P, Q)")
        println(io, "=========================================================================")
        for i in keys(eq_const_Pe_tf) for (j, info) in eq_const_Pe_tf[i] println(io, "Gen $i (Fault Pe): ", info) end end
        for i in keys(eq_const_Qe_tf) for (j, info) in eq_const_Qe_tf[i] println(io, "Gen $i (Fault Qe): ", info) end end
        for i in keys(eq_const_Vd_tf) for (j, info) in eq_const_Vd_tf[i] println(io, "Gen $i (Fault Vd): ", info) end end
        for i in keys(eq_const_Vq_tf) for (j, info) in eq_const_Vq_tf[i] println(io, "Gen $i (Fault Vq): ", info) end end
        for i in keys(eq_const_Pbalance_tf) for (j, info) in eq_const_Pbalance_tf[i] println(io, "Bus $i (Fault Pb): ", info) end end
        for i in keys(eq_const_Qbalance_tf) for (j, info) in eq_const_Qbalance_tf[i] println(io, "Bus $i (Fault Qb): ", info) end end
        
        for i in keys(eq_const_Pe_tpf) for (j, info) in eq_const_Pe_tpf[i] println(io, "Gen $i (Post-F Pe): ", info) end end
        for i in keys(eq_const_Qe_tpf) for (j, info) in eq_const_Qe_tpf[i] println(io, "Gen $i (Post-F Qe): ", info) end end
        for i in keys(eq_const_Vd_tpf) for (j, info) in eq_const_Vd_tpf[i] println(io, "Gen $i (Post-F Vd): ", info) end end
        for i in keys(eq_const_Vq_tpf) for (j, info) in eq_const_Vq_tpf[i] println(io, "Gen $i (Post-F Vq): ", info) end end
        for i in keys(eq_const_Pbalance_tpf) for (j, info) in eq_const_Pbalance_tpf[i] println(io, "Bus $i (Post-F Pb): ", info) end end
        for i in keys(eq_const_Qbalance_tpf) for (j, info) in eq_const_Qbalance_tpf[i] println(io, "Bus $i (Post-F Qb): ", info) end end
        println(io, "\n")

        println(io, "==========================================")
        println(io, "Equality Constraints Angle Swing Equation")
        println(io, "==========================================")
        for i in keys(eq_const_δ_tf) for (j, info) in eq_const_δ_tf[i] println(io, "Gen $i (Fault): ", info) end end
        for i in keys(eq_const_δ_tpf) for (j, info) in eq_const_δ_tpf[i] println(io, "Gen $i (Post-F): ", info) end end
        println(io, "\n")

        println(io, "===================================================")
        println(io, "Equality Constraints Speed Deviation Swing Equation")
        println(io, "===================================================")
        for i in keys(eq_const_Δω_tf) for (j, info) in eq_const_Δω_tf[i] println(io, "Gen $i (Fault): ", info) end end
        for i in keys(eq_const_Δω_tpf) for (j, info) in eq_const_Δω_tpf[i] println(io, "Gen $i (Post-F): ", info) end end
        println(io, "\n")
        
        println(io, "===================================================")
        println(io, "Equality Constraints EMF Transient Voltage Dynamics")
        println(io, "===================================================")
        for i in keys(eq_const_Ed_tf) for (j, info) in eq_const_Ed_tf[i] println(io, "Gen $i (Fault Ed): ", info) end end
        for i in keys(eq_const_Ed_tpf) for (j, info) in eq_const_Ed_tpf[i] println(io, "Gen $i (Post-F Ed): ", info) end end
        for i in keys(eq_const_Eq_tf) for (j, info) in eq_const_Eq_tf[i] println(io, "Gen $i (Fault Eq): ", info) end end
        for i in keys(eq_const_Eq_tpf) for (j, info) in eq_const_Eq_tpf[i] println(io, "Gen $i (Post-F Eq): ", info) end end
        println(io, "\n")

        println(io, "===================================================")
        println(io, "Inequality Constraints Angle in Relation to the COI")
        println(io, "===================================================")
        for i in keys(ineq_const_δ_COI_tf) for (j, info) in ineq_const_δ_COI_tf[i] println(io, "Gen $i (Fault): ", info) end end
        for i in keys(ineq_const_δ_COI_tpf) for (j, info) in ineq_const_δ_COI_tpf[i] println(io, "Gen $i (Post-F): ", info) end end
        println(io, "\n")

        println(io, "=============================================================")
        println(io, "Inequality Constraints Speed Deviation in Relation to the COI")
        println(io, "=============================================================")
        for i in keys(ineq_const_Δω_COI_tf) for (j, info) in ineq_const_Δω_COI_tf[i] println(io, "Gen $i (Fault): ", info) end end
        for i in keys(ineq_const_Δω_COI_tpf) for (j, info) in ineq_const_Δω_COI_tpf[i] println(io, "Gen $i (Post-F): ", info) end end
        println(io, "\n")

        # New pre-fault constraint sections:
        println(io, "==============================================")
        println(io, "Equality Constraints Exciter Initialization  ")
        println(io, "==============================================")
        for (i, info) in eq_const_Efd_init println(io, "Gen $i (Efd_init): ", info) end
        println(io, "\n")

        println(io, "==============================================")
        println(io, "Equality Constraints Power Initialization    ")
        println(io, "==============================================")
        for (i, info) in eq_const_P_init println(io, "Gen $i (P_init): ", info) end
        for (i, info) in eq_const_Q_init println(io, "Gen $i (Q_init): ", info) end
        println(io, "\n")

        # New time-series variable blocks (fault and post-fault):
        println(io, "============================")
        println(io, "Variables Te")
        println(io, "============================")
        for (j, infoj) in Te_tf
            for (k, infok) in infoj
                println(io, "$k: ", infok)
            end
        end
        for (j, infoj) in Te_tpf
            for (k, infok) in infoj
                println(io, "$k: ", infok)
            end
        end
        # New time-series variable blocks (fault and post-fault):
        println(io, "============================")
        println(io, "Variables Pm")
        println(io, "============================")
        for (j, infoj) in P_mech_tf
            for (k, infok) in infoj
                println(io, "$k: ", infok)
            end
        end
        for (j, infoj) in P_mech_tpf
            for (k, infok) in infoj
                println(io, "$k: ", infok)
            end
        end
    end

    cd(current_path_folder)
    println("Dynamic Model successfully saved as TXT file in: ", joinpath(path_folder_results,"Transient_Stability"))
end

function Save_Duals_Dynamic_Model(model::Model,
    E_fd::OrderedDict{Int, VariableRef}, δ::OrderedDict{Int, VariableRef}, Pm::OrderedDict{Int, VariableRef}, eq_const_Pm::OrderedDict{Int, ConstraintRef}, eq_const_Efd_init::OrderedDict{Int, ConstraintRef}, eq_const_P_init::OrderedDict{Int, ConstraintRef}, eq_const_Q_init::OrderedDict{Int, ConstraintRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Qe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, V_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, θ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Δω_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δCOI_tf::OrderedDict{Int, VariableRef}, ΔωCOI_tf::OrderedDict{Int, VariableRef},
    eq_const_δCOI_tf::OrderedDict{Int, ConstraintRef}, eq_const_ΔωCOI_tf::OrderedDict{Int, ConstraintRef}, eq_const_Pe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Pbalance_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qbalance_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vd_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vq_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_δ_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Δω_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Ed_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Eq_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_E_fd_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_valve_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_mech_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_δ_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_Δω_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Qe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, V_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, θ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, Δω_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, δCOI_tpf::OrderedDict{Int64, VariableRef}, ΔωCOI_tpf::OrderedDict{Int64, VariableRef},
    eq_const_δCOI_tpf::OrderedDict{Int, ConstraintRef}, eq_const_ΔωCOI_tpf::OrderedDict{Int, ConstraintRef}, eq_const_Pe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Pbalance_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Qbalance_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vd_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Vq_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_δ_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Δω_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Ed_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_Eq_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_E_fd_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_valve_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, eq_const_P_mech_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_δ_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}}, ineq_const_Δω_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    base_MVA::Float64, current_path_folder::String, path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability"))

    dual_Pm =[JuMP.dual(info) for (i, info) in eq_const_Pm] ./ base_MVA

    dual_δCOI_tf  =[JuMP.dual(info) for (i, info) in eq_const_δCOI_tf ] 
    dual_δCOI_tpf =[JuMP.dual(info) for (i, info) in eq_const_δCOI_tpf]
    dual_δCOI = vcat(dual_δCOI_tf, dual_δCOI_tpf)

    dual_ΔωCOI_tf  =[JuMP.dual(info) for (i, info) in eq_const_ΔωCOI_tf ] 
    dual_ΔωCOI_tpf =[JuMP.dual(info) for (i, info) in eq_const_ΔωCOI_tpf]
    dual_ΔωCOI = vcat(dual_ΔωCOI_tf, dual_ΔωCOI_tpf)

    # Generators Power Equations
    dual_Pe = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Pe_tf)
        dual_Pe[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Pe_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Pe_tpf[k]]) ./ base_MVA
    end
    dual_Qe = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Qe_tf)
        dual_Qe[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Qe_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Qe_tpf[k]]) ./ base_MVA
    end
    
    # 4th Order Stator Algebraic Equations
    dual_Vd = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Vd_tf)
        dual_Vd[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Vd_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Vd_tpf[k]]) 
    end
    dual_Vq = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Vq_tf)
        dual_Vq[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Vq_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Vq_tpf[k]]) 
    end

    # Nodal Balance Equations
    dual_Pbalance = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Pbalance_tf)
        dual_Pbalance[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Pbalance_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Pbalance_tpf[k]]) ./ base_MVA
    end
    dual_Qbalance = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Qbalance_tf)
        dual_Qbalance[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Qbalance_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Qbalance_tpf[k]]) ./ base_MVA
    end

    # Dynamics Equations
    dual_δ = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_δ_tf)
        dual_δ[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_δ_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_δ_tpf[k]])
    end
    dual_Δω = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Δω_tf)
        dual_Δω[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Δω_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Δω_tpf[k]])
    end
    
    # 4th Order EMF Dynamics Equations
    dual_Ed = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Ed_tf)
        dual_Ed[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Ed_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Ed_tpf[k]])
    end
    dual_Eq = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_Eq_tf)
        dual_Eq[k] = vcat([JuMP.dual(v) for (_, v) in eq_const_Eq_tf[k]],[JuMP.dual(v) for (_, v) in eq_const_Eq_tpf[k]])
    end

    dual_δ_COI = OrderedDict{Int, Vector{Float64}}()
    for k in keys(ineq_const_δ_COI_tf)
        dual_δ_COI[k] = vcat([JuMP.dual(v) for (_, v) in ineq_const_δ_COI_tf[k]],[JuMP.dual(v) for (_, v) in ineq_const_δ_COI_tpf[k]])
    end
    dual_Δω_COI = OrderedDict{Int, Vector{Float64}}()
    for k in keys(ineq_const_Δω_COI_tf)
        dual_Δω_COI[k] = vcat([JuMP.dual(v) for (_, v) in ineq_const_Δω_COI_tf[k]],[JuMP.dual(v) for (_, v) in ineq_const_Δω_COI_tpf[k]])
    end

    dual_LB_E_fd =[JuMP.dual(LowerBoundRef(info)) for (i, info) in E_fd]
    dual_UB_E_fd =[JuMP.dual(UpperBoundRef(info)) for (i, info) in E_fd]
    dual_LB_δ =[JuMP.dual(LowerBoundRef(info)) for (i, info) in δ]
    dual_UB_δ =[JuMP.dual(UpperBoundRef(info)) for (i, info) in δ]
    dual_LB_Pm =[JuMP.dual(LowerBoundRef(info)) for (i, info) in Pm] ./ base_MVA
    dual_UB_Pm =[JuMP.dual(UpperBoundRef(info)) for (i, info) in Pm] ./ base_MVA

    # E_fd dynamics (AVR)
    dual_E_fd = OrderedDict{Int, Vector{Float64}}()
    for k in keys(eq_const_E_fd_tf)
        dual_E_fd[k] = vcat(
            [JuMP.dual(v) for (_, v) in eq_const_E_fd_tf[k]],
            [JuMP.dual(v) for (_, v) in eq_const_E_fd_tpf[k]]
        )
    end

    # Write to file
    open("dynamic_model_duals.txt", "w") do io
        function write_dual_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec) println(io, "[$i] =\t €/MW $val") end
            println(io)
        end
        function write_dual_others(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec) println(io, "[$i] =\t $val") end
            println(io)
        end

        write_dual_others(io, "dual_eq_Pm",     dual_Pm)
        write_dual_others(io, "dual_eq_δCOI",  dual_δCOI)
        write_dual_others(io, "dual_eq_ΔωCOI", dual_ΔωCOI)

        for (i, info) in dual_Pe write_dual_power(io, "dual_eq_Pe[G$i]", info) end
        for (i, info) in dual_Qe write_dual_power(io, "dual_eq_Qe[G$i]", info) end
        for (i, info) in dual_Pbalance write_dual_power(io, "dual_eq_Pbalance[Bus $i]", info) end
        for (i, info) in dual_Qbalance write_dual_power(io, "dual_eq_Qbalance[Bus $i]", info) end
        
        for (i, info) in dual_Vd write_dual_others(io, "dual_eq_Vd[G$i]", info) end
        for (i, info) in dual_Vq write_dual_others(io, "dual_eq_Vq[G$i]", info) end

        for (i, info) in dual_δ write_dual_others(io, "dual_eq_δ[G$i]", info) end
        for (i, info) in dual_Δω write_dual_others(io, "dual_eq_Δω[G$i]", info) end
        for (i, info) in dual_Ed write_dual_others(io, "dual_eq_Ed[G$i]", info) end
        for (i, info) in dual_Eq write_dual_others(io, "dual_eq_Eq[G$i]", info) end
        for (i, info) in dual_E_fd write_dual_others(io, "dual_E_fd[G$i]", info) end

        for (i, info) in dual_δ_COI write_dual_others(io, "dual_ineq_δ_COI[G$i]", info) end
        for (i, info) in dual_Δω_COI write_dual_others(io, "dual_ineq_Δω_COI[G$i]", info) end

        write_dual_others(io, "dual_LB_E_fd",  dual_LB_E_fd)
        write_dual_others(io, "dual_UB_E_fd",  dual_UB_E_fd)
        write_dual_others(io, "dual_LB_δ",  dual_LB_δ)
        write_dual_others(io, "dual_UB_δ",  dual_UB_δ)
        write_dual_power(io,  "dual_LB_Pm", dual_LB_Pm)
        write_dual_power(io,  "dual_UB_Pm", dual_UB_Pm)
    end

    cd(current_path_folder)
    println("Duals of the dynamic model successfully saved as TXT file in: ", joinpath(path_folder_results, "Transient_Stability"))
end


# ===================================================================================
#                  PRINT THE DUALS OF THE OPTIMIZATION PROBLEM
# ===================================================================================
# Function to obtain the duals of the AC-OPF
function Save_Duals_ACOPF_Model(model::Model,
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    eq_const_angle_sw::ConstraintRef, 
    eq_const_p_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_balance::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_s_ik::OrderedDict{Int64, ConstraintRef},
    ineq_const_s_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_diff_ang::OrderedDict{Int64, ConstraintRef},
    base_MVA::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "ACOPF")) # Load the results path folder

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint of angle at the swing bus
    dual_θ_SW = JuMP.dual(eq_const_angle_sw) 

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint active power balance (Nodal KCL)
    dual_P_balance =[JuMP.dual(info) for (i, info) in eq_const_p_balance] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint reactive power balance (Nodal KCL)
    dual_Q_balance =[JuMP.dual(info) for (i, info) in eq_const_q_balance] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint apparent power flow capacity from i to k
    dual_Sik =[JuMP.dual(info) for (i, info) in ineq_const_s_ik] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint apparent power flow capacity from k to i
    dual_Ski =[JuMP.dual(info) for (i, info) in ineq_const_s_ki] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint voltage angle differences between buses
    dual_diff_ang =[JuMP.dual(info) for (i, info) in ineq_const_diff_ang]

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the active power of each generator
    dual_LB_Pg =[JuMP.has_lower_bound(info) ? JuMP.dual(LowerBoundRef(info)) : 0.0 for (i, info) in P_g] ./ base_MVA

    # Dual of the UPPER bound of the active power of each generator
    dual_UB_Pg =[JuMP.has_upper_bound(info) ? JuMP.dual(UpperBoundRef(info)) : 0.0 for (i, info) in P_g] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the reactive power of each generator
    dual_LB_Qg =[JuMP.has_lower_bound(info) ? JuMP.dual(LowerBoundRef(info)) : 0.0 for (i, info) in Q_g] ./ base_MVA

    # Dual of the UPPER bound of the reactive power of each generator
    dual_UB_Qg =[JuMP.has_upper_bound(info) ? JuMP.dual(UpperBoundRef(info)) : 0.0 for (i, info) in Q_g] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of each bus voltage
    dual_LB_V =[JuMP.has_lower_bound(info) ? JuMP.dual(LowerBoundRef(info)) : 0.0 for (i, info) in V]

    # Dual of the UPPER bound of each bus voltage
    dual_UB_V =[JuMP.has_upper_bound(info) ? JuMP.dual(UpperBoundRef(info)) : 0.0 for (i, info) in V]


    # ========== WRITE TO TXT FILE ==========
    open("ACOPF_duals.txt", "w") do io
        function write_dual_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MW $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_others(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t $val")
            end
            println(io)  # empty line between sections
        end

        println(io, "======================================")
        println(io, "          dual_θ_SW:")
        println(io, "======================================")
        println(io, "[1] =\t $dual_θ_SW")
        println(io)

        write_dual_power(io, "dual_P_balance", dual_P_balance)
        write_dual_power(io, "dual_Q_balance", dual_Q_balance)
        write_dual_power(io, "dual_Sik",       dual_Sik)
        write_dual_power(io, "dual_Ski",       dual_Ski)
        write_dual_others(io, "dual_diff_ang", dual_diff_ang)

        write_dual_power(io,  "dual_LB_Pg", dual_LB_Pg)
        write_dual_power(io,  "dual_UB_Pg", dual_UB_Pg)
        write_dual_power(io,  "dual_LB_Qg", dual_LB_Qg)
        write_dual_power(io,  "dual_UB_Qg", dual_UB_Qg)
        write_dual_others(io, "dual_LB_V",  dual_LB_V)
        write_dual_others(io, "dual_UB_V",  dual_UB_V)

    end

    println("Duals of the AC-OPF model successfully saved as TXT file in: ", joinpath(path_folder_results, "ACOPF"))

    cd(current_path_folder)

end