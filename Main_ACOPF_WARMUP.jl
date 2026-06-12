cd(dirname(@__FILE__));

#=
CODE FOR SOLVING THE AC POWER FLOW WITHOUT RELAXATIONS
AND ADDING TRANSIENT STABILITY CONSTRAINTS

Author:      Alex Junior da Cunha Coelho
Supervisors: Luis Badesa Bernardo and Araceli Hernandez Bayo
Affiliation: Technical University of Madrid
August 2025

===================================================================
                        IMPORTANT NOTES 
===================================================================
Taps of transformers are considered tap:1 (i.e., from tap -> to 1)
Phase shift of transformers must be in degrees in the input file
=#

#---------------------------
# INCLUDE THE PACKAGES USED
#---------------------------
# Packages related to linear algebra
using LinearAlgebra, SparseArrays

# Packages related to treatement of data
using Dates, NumericIO, DataFrames, Printf, CSV, DataStructures

# Packages related to the optimization
using JuMP, Ipopt, MadNLP
using HSL_jll

# Packages for plotting
using Plots, LaTeXStrings, Measures

# Packages for numerical integration
using Trapz

#---------------------------------
# INCLUDE AUXILIAR FUNCTION FILES
#--------------------------------
include("AF_CLEAN_TERMINAL.jl")         # Auxiliar function to clean the terminal
include("AF_DERIVATIVE_TECHINIQUES.jl") # Auxiliar function to calculate numerical derivatives
include("AF_READ_DATA.jl")              # Auxiliar functions used to read input data
include("AF_YBUS_SPARSE.jl")                   # Auxiliar function to create Ybus
include("BUILD_ACOPF_MODEL_KCL.jl")     # Auxiliar function to create the AC OPF model for optimization based on the KCL (this creates a reduced version neglecting OFF components)
include("AF_SAVE_OUTPUT.jl")            # Auxiliar function to save the output results
include("AF_MANAGEMENT.jl")             # Auxiliar function that calculates AC power flow and manage some data
include("BUILD_TSC_NOSPEED.jl") 


Clean_Terminal() # Clean the terminal

#-----------------------------------------
# Generate a folder to export the results
#-----------------------------------------
current_path_folder = pwd()                                            # Directory of the current folder
name_path_results   = "Results"                                        # Name of the folder to save the results (it must be created in advance)
path_folder_results = joinpath(current_path_folder, name_path_results) # Results directory
cd(current_path_folder)                                                # Load the current folder

#=
----------------------------------------------
 Relevant input variables to solve the AC-OPF
----------------------------------------------
** Select a system from the options below: **
9bus
9bus_Conejo_Paper
39bus
39bus_Conejo_Paper
118bus
=#

case     = "9bus" # Case under study (folder name)
base_MVA = 100.0  # Base Power [MVA]
load_factor = 1.5 # Coefficient factor that multiply the loads

δ_tol  = 100 # Module of maximum tolerance for internal angle [degrees]
Δω_tol = 0.5 # Module of maximum tolerance for angular speed deviation [p.u]

δ_tol  = (deg2rad(-δ_tol), deg2rad(δ_tol)) # Tupple with the min and max tolerance for internal angle
Δω_tol = (-Δω_tol, Δω_tol)                 # Tupple with the min and max tolerance for angular speed deviation

# System Parameters
f_syn = 50.0        # Synchronous frequency in [Hz]
ω_syn = 2π * f_syn  # Synchronous speed in     [rad/s]
Δω_0  = 0.0         # Initial speed deviation  [p.u.]

# ----------- Contingency Number 2 -------------- #
bus_fault = 7   # Faulted bus (3Φ short circuit)
circ_trip = [6] # Circuit(s) that trip(s) and clear(s) the fault


# Coefficients of the ZIP model
a = 0.0; b = 0.0; c = 1.0;
if (a + b + c) != 1.0
    throw(ArgumentError("Summation of ZIP load coefficients must be equal 1."))
end
ZIP = [a, b, c]

# Variables related to the simulation time
t_start_sim    = 0.0    # Time when the simulation is started           [seconds] 
t_end_sim      = 5.0    # Time to stop the simulation                   [seconds]
t_step_fault   = 0.01   # Time step to solve the differential equations [seconds]
t_step_postf   = 0.01
t_start_fault  = 0.01   # Time when the fault is started                [seconds]
clearing_time  = 0.30   # Time the fault lasts until clearing           [seconds] 
t_clear_fault  = clearing_time + t_start_fault   # Time to clear the fault [seconds] 
t_window_fault = collect(t_start_fault:t_step_fault:t_clear_fault)
t_window_postf = collect(t_clear_fault:t_step_postf:t_end_sim)
t_window_total = vcat(t_window_fault, t_window_postf)

if t_end_sim <= t_clear_fault
    throw(ArgumentError("t_end_sim is lower than t_clear_fault; PLEASE CORRECT."))
end

println("--------------------------------------------------------------------------------------------------------------------------------------")
Print_Input_Parameters(case, base_MVA, δ_tol, Δω_tol, f_syn, bus_fault, circ_trip, t_start_sim, t_end_sim, t_step_fault, t_start_fault, clearing_time, t_clear_fault, current_path_folder, path_folder_results)

#------------------------------------------
# Call the function to read the input data
#------------------------------------------
input_data_path_folder = joinpath(current_path_folder, "Input_Data", case) # Folder name where the input data is located
# Get the structs with data related to buses, generators and circuits
DBUS, DGEN, DGEN_DYN, DCIR, bus_mapping, reverse_bus_mapping = Read_Input_Data(input_data_path_folder) 

# Variable to multiply the power demanded by the loads
DBUS.p_d = load_factor .* (DBUS.p_d)
DBUS.q_d = load_factor .* (DBUS.q_d)
# ---------------------------------------------------

nBUS = length(DBUS.bus)      # Number of buses in the system
nGEN = length(DGEN.id)       # Number of generators in the system
if nGEN != length(DGEN_DYN.id) throw(ArgumentError("Dynamic and static data of generators do not match in length.")) end
nCIR = length(DCIR.from_bus) # Number of circuits in the system
cd(current_path_folder)      # Load the current folder

# ------------------
# Some sanity checks
# ------------------
if !(haskey(bus_mapping, bus_fault))
    throw(ArgumentError("The ID of the faulted bus does not exist in the CSV of input data."))
end

if maximum(circ_trip) > nCIR
    throw(ArgumentError("The IDs of the faulted circuits are greater than the ids in DCIR."))
end

#-------------------------------------------------------------------------------------------------
# Associates the buses with the generators and circuits connected to it, as well as adjacent buses
#-------------------------------------------------------------------------------------------------
bus_gen_circ_dict, bus_gen_circ_dict_ON = Manage_Bus_Gen_Circ(DBUS, DGEN, DCIR) 

# -----------------------------------
# Calculate the admittance matrices
# -----------------------------------

# Calculate the admittance matrix for the pre-fault period
Ybus_pref = Calculate_Ybus(DBUS, DCIR, nBUS, nCIR, base_MVA) # Admittance matrix for the pre-fault stage

# Calculate the new admittance matrix for the fault period and its reduced version 
Ybus_fault = Calculate_Ybus_fault(Ybus_pref, bus_fault)  # Admittance matrix augmented for the fault

# Calculate the new admittance matrix for the post-fault period and its reduced version 
l_status_postf = deepcopy(DCIR.l_status)
l_status_postf[circ_trip] .= 0
Ybus_postf = Calculate_Ybus_postf(DBUS, DCIR, nBUS, nCIR, base_MVA, l_status_postf)      # Admittance matrix augmented for the post-fault

cd(joinpath(path_folder_results,"Admittance_Matrices"))
df_Ybus_pref = DataFrame(Matrix(Ybus_pref), :auto)        # Convert the admittance matrix of the pre-fault into a DataFrame to save it
CSV.write("df_Ybus_pref.csv", df_Ybus_pref; delim=';')    # Save the admittance matrix of the pre-fault in a CSV file

df_Ybus_fault = DataFrame(Matrix(Ybus_fault), :auto)      # Convert the admittance matrix of the fault into a DataFrame to save it
CSV.write("df_Ybus_fault.csv", df_Ybus_fault; delim=';')  # Save the admittance matrix of the fault in a CSV file

df_Ybus_postf = DataFrame(Matrix(Ybus_postf), :auto)      # Convert the admittance matrix of the post-fault into a DataFrame to save it
CSV.write("df_Ybus_postf.csv", df_Ybus_postf; delim=';')  # Save the admittance matrix of the post-fault in a CSV file

println("Admittance matrices successfully saved in: ", joinpath(path_folder_results,"Admittance_Matrices"))
println("--------------------------------------------------------------------------------------------------------------------------------------")
cd(current_path_folder)

# ########################################################################################
#                                 STARTS OPTIMIZATION PROCESS 

#-----------------------------------
# Optimization model -> Setup
#-----------------------------------
optimizer = Ipopt.Optimizer
#optimizer = MadNLP.Optimizer
#solver = "ma57"
#solver = "ma97"
solver="ma57"

if optimizer == Ipopt.Optimizer
    optimizer = Ipopt.Optimizer
    model = JuMP.Model(optimizer)
    JuMP.set_optimizer_attribute(model, "tol", 1e-9)
    JuMP.set_optimizer_attribute(model, "acceptable_tol", 1e-9)
    JuMP.set_optimizer_attribute(model, "print_level", 5) # Verbosity (0–12, default = 5)
    JuMP.set_optimizer_attribute(model, "output_file", "ipopt_log.txt")     # Prints to console instead of a file
    JuMP.set_optimizer_attribute(model, "max_iter", 5_000)
    JuMP.set_silent(model)
    if solver=="ma57" || solver=="ma97"
        set_attribute(model, "hsllib", HSL_jll.libhsl_path)
        set_attribute(model, "linear_solver", solver)
    end

elseif optimizer == MadNLP.Optimizer
    optimizer = MadNLP.Optimizer
    model = JuMP.Model(optimizer)
    JuMP.set_optimizer_attribute(model, "tol", 1e-9)
    # JuMP.set_optimizer_attribute(model, "acceptable_tol", 1e-9)
    JuMP.set_optimizer_attribute(model, "print_level", MadNLP.INFO) # Verbosity
    JuMP.set_optimizer_attribute(model, "output_file", "ipopt_log.txt")     # Prints to console instead of a file
    JuMP.set_optimizer_attribute(model, "max_iter", 5_000)
    # JuMP.set_silent(model)
end

#------------------------------
# Build the Optimization Model
#------------------------------
time_to_build_model = time() # Start the timer to build the Optimization Model

# model, V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, 
# eq_const_p_ik, eq_const_q_ik, eq_const_p_ki, eq_const_q_ki, ineq_const_s_ik, ineq_const_s_ki, 
# ineq_const_diff_ang = Make_ACOPF_Model!(model, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR)

model, V, θ, P_g, Q_g, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, ineq_const_s_ik, ineq_const_s_ki, 
ineq_const_diff_ang = Make_ACOPF_Model!(model, Ybus_pref, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR)
#model, V, θ, P_g, Q_g, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, ineq_const_s_ik, ineq_const_s_ki, 
#ineq_const_diff_ang = Make_ACOPF_Model!(model, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR)

# =====================================================================================
# WARM START: SOLVE THE STEADY-STATE ACOPF FIRST
# =====================================================================================
println("\n--- Solving Steady-State ACOPF for Warm Start ---")
JuMP.optimize!(model)

if termination_status(model) ∈ (JuMP.OPTIMAL, JuMP.LOCALLY_SOLVED)
    println("ACOPF Warm Start Solved Successfully! Injecting exact roots into dynamic variables...\n")
else
    @warn "ACOPF Warm Start failed! The dynamic model will likely struggle."
end
val_V  = Dict(i => JuMP.value(v) for (i, v) in V)
val_θ  = Dict(i => JuMP.value(v) for (i, v) in θ)
val_Pg = Dict(i => JuMP.value(v) for (i, v) in P_g)
val_Qg = Dict(i => JuMP.value(v) for (i, v) in Q_g)
# =====================================================================================


# =====================================================================================

# Adding transient stability constraints
model, E_fd, δ, Pm, δCOI, Ed, Eq, Id, Iq, eq_const_P_init, eq_const_Q_init, eq_const_Ed_init, eq_const_Eq_init, eq_const_Vd_init, eq_const_Vq_init, eq_const_Pm = 
Define_Dyn_Var_EδPm!(
    model, DGEN, DGEN_DYN, nGEN, base_MVA, V, θ, P_g, Q_g, val_V, val_θ, val_Pg, val_Qg
)

model, V_ref, eq_const_Efd_init = Def_Dyn_Init_Exciter!(
model, DGEN, DGEN_DYN, nGEN, V, E_fd, val_V
)

# Initialize Governor
model, P_ref, P_valve, eq_const_Pref_init, eq_const_Pvalve_init = Def_Dyn_Init_Governor!(
model, DGEN, DGEN_DYN, nGEN, Pm
)

# Add variables and constraints of the fault period
model, Pe_tf, Qe_tf, V_tf, θ_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf, E_fd_tf, P_valve_tf, P_mech_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf, eq_const_Pe_tf, eq_const_Qe_tf, eq_const_Pbalance_tf, eq_const_Qbalance_tf, eq_const_Vd_tf, eq_const_Vq_tf, eq_const_δ_tf, eq_const_Δω_tf, eq_const_Ed_tf, eq_const_Eq_tf, eq_const_E_fd_tf, eq_const_P_valve_tf, eq_const_P_mech_tf = Def_Dyn_Fault_All!(
model, DBUS, bus_gen_circ_dict_ON, DGEN, DGEN_DYN, nBUS, nGEN, base_MVA, ZIP, t_window_fault, δ_tol, Δω_tol, Ybus_fault, V, θ, E_fd, V_ref, P_ref, P_valve, Pm, δ, Δω_0, P_g, ω_syn, t_step_fault, val_V, val_θ, val_Pg, val_Qg, Ed, Eq, Id, Iq
)

# Add variables and constraints of the post-fault period
model, Pe_tpf, Qe_tpf, V_tpf, θ_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf, E_fd_tpf, P_valve_tpf, P_mech_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf, eq_const_Pe_tpf, eq_const_Qe_tpf, eq_const_Pbalance_tpf, eq_const_Qbalance_tpf, eq_const_Vd_tpf, eq_const_Vq_tpf, eq_const_δ_tpf, eq_const_Δω_tpf, eq_const_Ed_tpf, eq_const_Eq_tpf, eq_const_E_fd_tpf, eq_const_P_valve_tpf, eq_const_P_mech_tpf = Def_Dyn_PostF_All!(
model, DBUS, bus_gen_circ_dict_ON, DGEN, DGEN_DYN, nBUS, nGEN, base_MVA, ZIP, t_window_postf, δ_tol, Δω_tol, Ybus_postf, V, θ, E_fd, V_ref, P_ref, Pm, δ_tf, Δω_tf, Pe_tf, V_tf, E_fd_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf, P_valve_tf, P_mech_tf, ω_syn, t_step_postf, val_V, val_θ, val_Pg, val_Qg, δ, Ed, Eq, Id, Iq
)

# # =====================================================================================

#-------------------------------------------------------------------------------------
#                         SAVE MODEL SUMMARY AND DETAILS
#-------------------------------------------------------------------------------------
println("--------------------------------------------------------------------------------------------------------------------------------------")
#Export_ACOPF_Model(model, V, θ, P_g, Q_g, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, 
#ineq_const_s_ik, ineq_const_s_ki, ineq_const_diff_ang, current_path_folder, path_folder_results)

# # Export_Dyn_Model(model, E, δ, Pm, eq_const_Pm, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, eq_const_Pe_tf, eq_const_δ_tf,
# # eq_const_Δω_tf, ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf, Pe_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, eq_const_Pe_tpf,
# # eq_const_δ_tpf, eq_const_Δω_tpf, ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf, current_path_folder, path_folder_results)
# # println("--------------------------------------------------------------------------------------------------------------------------------------")

# ---------------------------------
#  Solve the optmization problem
# ---------------------------------
cd(path_folder_results)
time_to_solve_model = time()                       # Start the timer to solve the Optimization Model
JuMP.optimize!(model)                              # Optimize model
time_to_solve_model = time() - time_to_solve_model # End the timer to build the Optimization Model
println("\nTime to solve the model: $time_to_solve_model sec")
status_model = JuMP.termination_status(model)
println("Termination Status: $status_model \n")
println("--------------------------------------------------------------------------------------------------------------------------------------")
cd(current_path_folder)

#                               ENDS OPTIMIZATION PROCESS 
# ########################################################################################

RBUS::Union{Nothing, DataFrame} = nothing
RGEN::Union{Nothing, DataFrame} = nothing
RCIR::Union{Nothing, DataFrame} = nothing

if status_model == OPTIMAL || status_model == LOCALLY_SOLVED || status_model == ITERATION_LIMIT
    #-------------------------------------------------------------------------------------
    #                             SAVE RESULTS 
    #-------------------------------------------------------------------------------------
    RBUS, RGEN, RCIR = Save_Solution_Model_KCL(model, V, θ, P_g, Q_g, bus_gen_circ_dict_ON, DBUS, DGEN, DCIR,
    base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping, current_path_folder, path_folder_results)

    #-------------------------------------------------------------------------------------
    #                          SAVE DYNAMIC RESULTS 
    #-------------------------------------------------------------------------------------
     Manage_Dyn_Results(model,
        DGEN_DYN,
        E_fd, δ, Pm, Ed, Eq, Id, Iq,
        Pe_tf, Qe_tf, V_tf, θ_tf, δ_tf, Δω_tf, Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf, E_fd_tf, P_valve_tf, P_mech_tf,
        δCOI_tf, ΔωCOI_tf,
        Pe_tpf, Qe_tpf, V_tpf, θ_tpf, δ_tpf, Δω_tpf, Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf, E_fd_tpf, P_valve_tpf, P_mech_tpf,
        δCOI_tpf, ΔωCOI_tpf,
        t_window_total, t_clear_fault, δ_tol, Δω_tol, base_MVA, f_syn, ω_syn,
        current_path_folder, path_folder_results
    )

    #-------------------------------------------------------------------------------------
    #                             SAVE DUALS 
    #-------------------------------------------------------------------------------------

    Save_Duals_ACOPF_Model(model, V, θ, P_g, Q_g, eq_const_angle_sw, 
                       eq_const_p_balance, eq_const_q_balance, 
                       ineq_const_s_ik, ineq_const_s_ki, 
                       ineq_const_diff_ang, base_MVA, 
                       current_path_folder, path_folder_results)


    Save_Duals_Dynamic_Model(
    model,
    # Initial Variables
    E_fd, δ, Pm, eq_const_Pm,
    eq_const_Efd_init, eq_const_P_init, eq_const_Q_init,  # Missing these 3 arguments
    
    # Fault Variables
    Pe_tf, Qe_tf, V_tf, θ_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf,
    
    # Fault Constraints
    eq_const_δCOI_tf, eq_const_ΔωCOI_tf,
    eq_const_Pe_tf, eq_const_Qe_tf,
    eq_const_Pbalance_tf, eq_const_Qbalance_tf,
    eq_const_Vd_tf, eq_const_Vq_tf,
    eq_const_δ_tf, eq_const_Δω_tf,
    eq_const_Ed_tf, eq_const_Eq_tf,
    eq_const_E_fd_tf, eq_const_P_valve_tf, eq_const_P_mech_tf,  # Missing these 3 constraints
    ineq_const_δ_COI_tf, ineq_const_Δω_COI_tf,
    
    # Post-Fault Variables
    Pe_tpf, Qe_tpf, V_tpf, θ_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf,
    
    # Post-Fault Constraints
    eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf,
    eq_const_Pe_tpf, eq_const_Qe_tpf,
    eq_const_Pbalance_tpf, eq_const_Qbalance_tpf,
    eq_const_Vd_tpf, eq_const_Vq_tpf,
    eq_const_δ_tpf, eq_const_Δω_tpf,
    eq_const_Ed_tpf, eq_const_Eq_tpf,
    eq_const_E_fd_tpf, eq_const_P_valve_tpf, eq_const_P_mech_tpf,  # Missing these 3 constraints
    ineq_const_δ_COI_tpf, ineq_const_Δω_COI_tpf,
    
    # Final Parameters
    base_MVA, current_path_folder, path_folder_results
)

    
#     V_tf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in V_tf)
#     θ_tf_values = OrderedDict(i => [rad2deg.(JuMP.value(v)) for (_, v) in inner] for (i, inner) in θ_tf)

else
    JuMP.@warn "Optmization process failed. No feasible solution found."
end
println("--------------------------------------------------------------------------------------------------------------------------------------")


# E_optim = [JuMP.value(v) for (i, v) in E]
# δ_optim = [rad2deg(JuMP.value(v)) for (i, v) in δ]
# Pm_optim = [JuMP.value(v) for (i, v) in Pm]
# δCOI_optim = rad2deg(JuMP.value(δCOI))