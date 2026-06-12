# Function used to Read the Input data from the CSV files and store them into Structs
function Read_Input_Data(folder_path::String)

    # ====================================================================================================
    # If some modification is done in the name of the variables in the CSV files, it must be modified here
    # ==================================================================================================== 

    # Function to read buses data and store in DataFrame
    function read_bus_data()
        df = CSV.read("bus_data.csv", DataFrame; delim=';')  # Read CSV file
        
        df_bus = deepcopy(df)
        rename!(df_bus, Dict(
            names(df_bus)[1] => :bus,       # Bus number
            names(df_bus)[2] => :type,      # Bus type: (3 = SW), (2 = PV) or (1 = PQ) 
            names(df_bus)[3] => :p_d,       # Active power demanded (MW)
            names(df_bus)[4] => :q_d,       # Reactive power demanded (MVAr) 
            names(df_bus)[5] => :g_sh,      # Shunt conductance (MW) 
            names(df_bus)[6] => :b_sh,      # Shunt susceptance (MVAr) 
            names(df_bus)[7] => :area,      # Area where the bus is located
            names(df_bus)[8] => :v_spe,     # Voltage specified (p.u)
            names(df_bus)[9] => :v_a,       # Voltage angles (rad)
            names(df_bus)[10] => :base_kV,  # Base voltage of the system (kV)
            names(df_bus)[11] => :zone,     # Zone where the bus is located
            names(df_bus)[12] => :v_max,    # Maximum bus voltage (p.u)
            names(df_bus)[13] => :v_min,    # Minimum bus voltage (p.u)
        ))

        df_bus.bus = Int64.(df_bus.bus)            # Bus number   
        df_bus.type = Int64.(df_bus.type)          # Bus type: (3 = SW), (2 = PV) or (1 = PQ)   
        df_bus.p_d = Float64.(df_bus.p_d)          # Active power demanded (MW) 
        df_bus.q_d = Float64.(df_bus.q_d)          # Reactive power demanded (MVAr)    
        df_bus.g_sh = Float64.(df_bus.g_sh)        # Shunt conductance (MW)      
        df_bus.b_sh = Float64.(df_bus.b_sh)        # Shunt susceptance (MVAr)    
        df_bus.area = Int64.(df_bus.area)          # Area where the bus is located   
        df_bus.v_spe = Float64.(df_bus.v_spe)      # Voltage specified (p.u)      
        df_bus.v_a = Float64.(df_bus.v_a)          # Voltage angles (rad)      
        df_bus.base_kV = Float64.(df_bus.base_kV)  # Base voltage of the system (kV)
        df_bus.zone = Int64.(df_bus.zone)          # Zone where the bus is located 
        df_bus.v_max = Float64.(df_bus.v_max)      # Maximum bus voltage (p.u)
        df_bus.v_min = Float64.(df_bus.v_min)      # Minimum bus voltage (p.u)   

        return df_bus
    end

    # Function to read generators data and store in DataFrame
    function read_gen_data()
        df = CSV.read("generators_data.csv", DataFrame; delim=';') # Read CSV file
        num_gen = length(df.bus)
        id = collect(1:num_gen)

        df_gen = hcat(DataFrame(id = id), df; makeunique=true)

        rename!(df_gen, Dict(
            names(df_gen)[2] => :bus,        # Bus number
            names(df_gen)[3] => :pg_spe,     # Active power for generation specified (MW)
            names(df_gen)[4] => :qg_spe,     # Reactive power for generation specified (MVAr)
            names(df_gen)[5] => :qg_max,     # Maximum reactive power capacity (MVAr)
            names(df_gen)[6] => :qg_min,     # Minimum reactive power capacity (MVAr)
            names(df_gen)[7] => :vg_spe,     # Voltage specified for the generator (p.u)
            names(df_gen)[8] => :base_MVA,   # Base power of the system (MVA)
            names(df_gen)[9] => :g_status,   # Generator status -> ON or OFF
            names(df_gen)[10] => :pg_max,    # Maximum active power capacity (MW)
            names(df_gen)[11] => :pg_min,    # Minimum active power capacity (MW)  
            names(df_gen)[12] => :g_cost_2,  # Generation cost (quadratic) (€/MW)
            names(df_gen)[13] => :g_cost_1,  # Generation cost (linear) (€/MW)
            names(df_gen)[14] => :g_cost_0,  # Generation cost (fixed) (€/MW)
        ))

        df_gen.id = Int64.(df_gen.id)               # Generator identifier 
        df_gen.bus = Int64.(df_gen.bus)             # Bus number 
        df_gen.pg_spe = Float64.(df_gen.pg_spe)     # Active power for generation specified (MW)
        df_gen.qg_spe = Float64.(df_gen.qg_spe)     # Reactive power for generation specified (MVAr)
        df_gen.qg_max = Float64.(df_gen.qg_max)     # Maximum reactive power capacity (MVAr)
        df_gen.qg_min = Float64.(df_gen.qg_min)     # Minimum reactive power capacity (MVAr)
        df_gen.vg_spe = Float64.(df_gen.vg_spe)     # Voltage specified for the generator (p.u)
        df_gen.base_MVA = Float64.(df_gen.base_MVA) # Base power of the system (MVA)
        df_gen.g_status = Int64.(df_gen.g_status)   # Generator status   
        df_gen.pg_max = Float64.(df_gen.pg_max)     # Maximum active power capacity (MW)  
        df_gen.pg_min = Float64.(df_gen.pg_min)     # Minimum active power capacity (MW)  
        df_gen.g_cost_2 = Float64.(df_gen.g_cost_2) # Generation cost (quadratic) (€/MW)
        df_gen.g_cost_1 = Float64.(df_gen.g_cost_1) # Generation cost (linear) (€/MW)
        df_gen.g_cost_0 = Float64.(df_gen.g_cost_0) # Generation cost (fixed) (€/MW)
        
        return df_gen
    end

    # Function to read generators dynamic data and store in DGEN_DYNAMIC_Struct
    function read_gen_dynamic_data()
        df = CSV.read("gen_dynamic_data.csv", DataFrame; delim=';') # Read CSV file
        num_gen = length(df.bus)
        id = collect(1:num_gen)

        df_gen_dyn = hcat(DataFrame(id = id), df; makeunique=true)

        rename!(df_gen_dyn, Dict(
            names(df_gen_dyn)[1] => :id,       # Generator ID
            names(df_gen_dyn)[2] => :bus,      # Bus ID
            names(df_gen_dyn)[3] => :Eg,       # Internal Voltage Magnitude
            names(df_gen_dyn)[4] => :Xd_tr,    # Direct-Axis Transient Reactance (original "Xd" values)
            names(df_gen_dyn)[5] => :Xq_tr,    # Quadrature-Axis Transient Reactance
            names(df_gen_dyn)[6] => :Xd,       # Direct-Axis Synchronous Reactance
            names(df_gen_dyn)[7] => :Xq,       # Quadrature-Axis Synchronous Reactance
            names(df_gen_dyn)[8] => :Td,       # Direct-Axis Transient Time Constant
            names(df_gen_dyn)[9] => :Tq,       # Quadrature-Axis Transient Time Constant
            names(df_gen_dyn)[10] => :H,       # Inertia Constant
            names(df_gen_dyn)[11] => :D,       # Damping Coefficient
            names(df_gen_dyn)[12] => :Ra,      # Armature Resistance
            names(df_gen_dyn)[13] => :T_exc,    # Exciter Time Constant
            names(df_gen_dyn)[14] => :K_exc,    # Exciter Gain
            names(df_gen_dyn)[15] => :R,
            names(df_gen_dyn)[16] => :T1,
            names(df_gen_dyn)[17] => :T2,
            names(df_gen_dyn)[18] => :T3,
        ))

        df_gen_dyn.id = Int64.(df_gen_dyn.id)           # Generator ID
        df_gen_dyn.bus = Int64.(df_gen_dyn.bus)         # Bus ID
        df_gen_dyn.Eg = Float64.(df_gen_dyn.Eg)         # Internal Voltage Magnitude
        df_gen_dyn.Xd_tr = Float64.(df_gen_dyn.Xd_tr)   # Direct-Axis Transient Reactance
        df_gen_dyn.Xq_tr = Float64.(df_gen_dyn.Xq_tr)   # Quadrature-Axis Transient Reactance
        df_gen_dyn.Xd = Float64.(df_gen_dyn.Xd)         # Direct-Axis Synchronous Reactance
        df_gen_dyn.Xq = Float64.(df_gen_dyn.Xq)         # Quadrature-Axis Synchronous Reactance
        df_gen_dyn.Td = Float64.(df_gen_dyn.Td)         # Direct-Axis Transient Time Constant
        df_gen_dyn.Tq = Float64.(df_gen_dyn.Tq)         # Quadrature-Axis Transient Time Constant
        df_gen_dyn.H = Float64.(df_gen_dyn.H)           # Inertia Constant
        df_gen_dyn.D = Float64.(df_gen_dyn.D)           # Damping Coefficient
        df_gen_dyn.Ra = Float64.(df_gen_dyn.Ra)         # Armature Resistance
        df_gen_dyn.T_exc = Float64.(df_gen_dyn.T_exc) 
        df_gen_dyn.K_exc = Float64.(df_gen_dyn.K_exc)
        df_gen_dyn.R = Float64.(df_gen_dyn.R)
        df_gen_dyn.T1 = Float64.(df_gen_dyn.T1)
        df_gen_dyn.T2 = Float64.(df_gen_dyn.T2)
        df_gen_dyn.T3 = Float64.(df_gen_dyn.T3)
        return df_gen_dyn
    end

    # Function to read circuits data and store in DataFrame
    function read_circuit_data()
        df = CSV.read("line_data.csv", DataFrame; delim=';')
        num_circ = length(df.fbus)
        id = collect(1:num_circ)

        df_cir = hcat(DataFrame(circ = id), df; makeunique=true)

        rename!(df_cir, Dict(
            names(df_cir)[2] => :from_bus,  # "From" bus  
            names(df_cir)[3] => :to_bus,    # "To" bus
            names(df_cir)[4] => :l_res,     # Branch resistance (p.u)  
            names(df_cir)[5] => :l_reac,    # Branch reactance (p.u) 
            names(df_cir)[6] => :l_sh_susp, # Branch shunt susceptance (Π model) (p.u)
            names(df_cir)[7] => :l_cap_1,   # Branch maximum capacity 1 (MW or MVA) 
            names(df_cir)[8] => :l_cap_2,   # Branch maximum capacity 2 (MW or MVA) 
            names(df_cir)[9] => :l_cap_3,   # Branch maximum capacity 3 (MW or MVA)
            names(df_cir)[10] => :t_tap,    # Transformer tap (p.u) 
            names(df_cir)[11] => :t_shift,  # Shift angle of the transformer (degrees)   
            names(df_cir)[12] => :l_status, # Branch ON or Branch OFF
            names(df_cir)[13] => :ang_min,  # Minimum angle (degrees)
            names(df_cir)[14] => :ang_max,  # Maximum angle (degrees)
        ))

        df_cir.circ = Int64.(df_cir.circ)             # Circuit identifier
        df_cir.from_bus = Int64.(df_cir.from_bus)     # "From" bus 
        df_cir.to_bus = Int64.(df_cir.to_bus)         # "To" bus
        df_cir.l_res = Float64.(df_cir.l_res)         # Branch resistance (p.u)
        df_cir.l_reac = Float64.(df_cir.l_reac)       # Branch reactance (p.u)
        df_cir.l_sh_susp = Float64.(df_cir.l_sh_susp) # Branch shunt susceptance (Π model) (p.u)
        df_cir.l_cap_1 = Float64.(df_cir.l_cap_1)     # Branch maximum capacity 1 (MW or MVA)
        df_cir.l_cap_2 = Float64.(df_cir.l_cap_2)     # Branch maximum capacity 2 (MW or MVA)
        df_cir.l_cap_3 = Float64.(df_cir.l_cap_3)     # Branch maximum capacity 3 (MW or MVA)
        df_cir.t_tap = Float64.(df_cir.t_tap)         # Transformer tap (p.u)
        df_cir.t_shift = Float64.(df_cir.t_shift)     # Shift angle of the transformer (degrees) 
        df_cir.l_status = Int64.(df_cir.l_status)     # Branch ON or Branch OFF 
        df_cir.ang_min = Float64.(df_cir.ang_min)     # Minimum angle (degrees) 
        df_cir.ang_max = Float64.(df_cir.ang_max)     # Maximum angle (degrees) 

        return df_cir
    end

    cd(folder_path) # Load the folder were the input data files are stored

    DBUS     = read_bus_data()          # Generate the DataFrame with Buses data
    DGEN     = read_gen_data()          # Generate the DataFrame with Generators data
    DGEN_DYN = read_gen_dynamic_data()  # Generate the DataFrame with Generators Dynamic data
    DCIR     = read_circuit_data()      # Generate the DataFrame with Circuits data

    # For the code to work properly, the bus indices must be set in ascending order from 1 to nBUS
    bus_mapping, reverse_bus_mapping = Mapping_Buses_Labels(DBUS) # Map the buses labels from old to new nomeclature
    DBUS.bus = [bus_mapping[b] for b in DBUS.bus]                 # Rename the buses labels from 1 to nBUS

    # Map the buses labels to be in ascending order from 1 to nBUS
    DGEN.bus      = [bus_mapping[b] for b in DGEN.bus]
    DGEN_DYN.bus  = [bus_mapping[b] for b in DGEN_DYN.bus]
    DCIR.from_bus = [bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [bus_mapping[b] for b in DCIR.to_bus]

    return DBUS, DGEN, DGEN_DYN, DCIR, bus_mapping, reverse_bus_mapping # Return the data
end

# Function used to map the from old to new nomeclature
function Mapping_Buses_Labels(DBUS::DataFrame)

    # Given bus numbers
    original_buses = DBUS.bus

    # Create a dictionary that maps original bus labels to new indices
    bus_mapping = OrderedDict(original_buses[i] => i for i in eachindex(original_buses))

    # Reverse mapping (for converting back later)
    reverse_bus_mapping = OrderedDict(i => original_buses[i] for i in eachindex(original_buses))

    return bus_mapping, reverse_bus_mapping
end

# Function that can change the buses labels according to the new nomenclature
function Change_Buses_Labels(DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, DCIR::DataFrame, bus_mapping::OrderedDict)
    
    # Convert using the reverse mapping
    DBUS.bus      = [bus_mapping[b] for b in DBUS.bus]
    DGEN.bus      = [bus_mapping[b] for b in DGEN.bus]
    DGEN_DYN.bus  = [bus_mapping[b] for b in DGEN_DYN.bus]
    DCIR.from_bus = [bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [bus_mapping[b] for b in DCIR.to_bus]

    return DBUS, DGEN, DGEN_DYN, DCIR
end

# Function that can return the buses labels according to the original nomenclature
function Reverse_Buses_Labels(DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, DCIR::DataFrame, reverse_bus_mapping::OrderedDict)
        
    # Convert using the reverse mapping
    DBUS.bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    DGEN.bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    DGEN_DYN.bus  = [reverse_bus_mapping[b] for b in DGEN_DYN.bus]
    DCIR.from_bus = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    return DBUS, DGEN, DGEN_DYN, DCIR
end
