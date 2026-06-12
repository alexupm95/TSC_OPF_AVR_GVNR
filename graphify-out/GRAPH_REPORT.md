# Graph Report - .  (2026-06-12)

## Corpus Check
- Corpus is ~22,945 words - fits in a single context window. You may not need a graph.

## Summary
- 85 nodes · 127 edges · 11 communities (10 shown, 1 thin omitted)
- Extraction: 74% EXTRACTED · 26% INFERRED · 0% AMBIGUOUS · INFERRED: 33 edges (avg confidence: 0.9)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Julia Package Dependencies|Julia Package Dependencies]]
- [[_COMMUNITY_AC Power Flow and Network Setup|AC Power Flow and Network Setup]]
- [[_COMMUNITY_Results Export and Output|Results Export and Output]]
- [[_COMMUNITY_Post-Fault Dynamics Constraints|Post-Fault Dynamics Constraints]]
- [[_COMMUNITY_Fault-Period Dynamics|Fault-Period Dynamics]]
- [[_COMMUNITY_Numerical Derivative Methods|Numerical Derivative Methods]]
- [[_COMMUNITY_Bus Data and Input Reading|Bus Data and Input Reading]]
- [[_COMMUNITY_Center of Inertia Dynamics|Center of Inertia Dynamics]]
- [[_COMMUNITY_ZIP Load and Power Equations|ZIP Load and Power Equations]]
- [[_COMMUNITY_Generator and AVR State Variables|Generator and AVR State Variables]]
- [[_COMMUNITY_Contingency and Stability Validation|Contingency and Stability Validation]]

## God Nodes (most connected - your core abstractions)
1. `Main_ACOPF_WARMUP (entry point)` - 13 edges
2. `Def_Dyn_Fault_All!()` - 12 edges
3. `Def_Dyn_PostF_All!()` - 12 edges
4. `Admittance Matrix (Ybus)` - 6 edges
5. `Save_Solution_Model_KCL()` - 5 edges
6. `Manage_Dyn_Results()` - 5 edges
7. `Make_ACOPF_Model!()` - 5 edges
8. `Calculate_Derivatives_Lagrange_Interpolator()` - 4 edges
9. `Def_Dyn_Fault_eqconst_PeQe!()` - 4 edges
10. `Def_Dyn_PostF_eqconst_PeQe!()` - 4 edges

## Surprising Connections (you probably didn't know these)
- `9-Bus Generator Data (UC3M Paper)` --references--> `Make_ACOPF_Model!()`  [INFERRED]
  Input_Data/9bus/generators_data_paperUC3M.txt → BUILD_ACOPF_MODEL_KCL.jl
- `Contingency List (9-bus, cross-validated Anatem/Julia)` --references--> `Transient Stability Constraints (TSC)`  [INFERRED]
  graphify-out/converted/Contingency List_4bcb6975.md → BUILD_TSC_NOSPEED.jl
- `Contingency List with Circuit Mapping (9-bus)` --references--> `Transient Stability Constraints (TSC)`  [INFERRED]
  graphify-out/converted/Contingency List_562779c0.md → BUILD_TSC_NOSPEED.jl
- `Main_ACOPF_WARMUP (entry point)` --calls--> `Print_Input_Parameters()`  [EXTRACTED]
  Main_ACOPF_WARMUP.jl → AF_SAVE_OUTPUT.jl
- `Def_Dyn_Fault_eqconst_PeQe!()` --shares_data_with--> `Admittance Matrix (Ybus)`  [INFERRED]
  BUILD_TSC_NOSPEED.jl → AF_YBUS_SPARSE.jl

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **TSC-OPF Full Solution Pipeline (warm-start ACOPF then TSC)** — main_acopf_warmup_main, build_acopf_model_kcl_make_acopf_model, build_tsc_nospeed_def_dyn_fault_all, build_tsc_nospeed_def_dyn_postf_all, af_save_output_manage_dyn_results [EXTRACTED 1.00]
- **Admittance Matrix for Three Operating Periods (pre-fault, fault, post-fault)** — af_ybus_sparse_calculate_ybus, af_ybus_sparse_calculate_ybus_fault, af_ybus_sparse_calculate_ybus_postf [INFERRED 0.95]
- **Dynamic Generator Model Components (AVR, Governor, Swing, 4th-Order EMF, COI)** — concept_avr, concept_governor, concept_swing_equation, concept_4th_order_generator, concept_coi [INFERRED 0.95]

## Communities (11 total, 1 thin omitted)

### Community 0 - "Julia Package Dependencies"
Cohesion: 0.12
Nodes (16): CSV, DataFrames, DataStructures, Dates, HSL_jll, Ipopt, JuMP, LaTeXStrings (+8 more)

### Community 1 - "AC Power Flow and Network Setup"
Cohesion: 0.19
Nodes (11): 9-Bus Generator Data (UC3M Paper), Clean_Terminal(), Calculate_AC_Power_Flow(), Manage_Bus_Gen_Circ(), Calculate_Ybus(), Calculate_Ybus_fault(), Calculate_Ybus_postf(), Make_ACOPF_Model!() (+3 more)

### Community 2 - "Results Export and Output"
Cohesion: 0.23
Nodes (7): Manage_Dyn_Results(), Print_Input_Parameters(), Save_Dyn_Results_CSV(), Save_Dyn_Results_Plots(), Save_ResultsCSV_ACOPF(), Save_ResultsTXT_ACOPF(), Save_Solution_Model_KCL()

### Community 3 - "Post-Fault Dynamics Constraints"
Cohesion: 0.31
Nodes (6): Def_Dyn_PostF_All!(), Def_Dyn_PostF_EMF!(), Def_Dyn_PostF_Exciter!(), Def_Dyn_PostF_Governor!(), Def_Dyn_PostF_rel_COI!(), Def_Dyn_PostF_Swing!()

### Community 4 - "Fault-Period Dynamics"
Cohesion: 0.29
Nodes (8): Def_Dyn_Fault_All!(), Def_Dyn_Fault_EMF!(), Def_Dyn_Fault_Exciter!(), Def_Dyn_Fault_Governor!(), Def_Dyn_Fault_rel_COI!(), Def_Dyn_Fault_Swing!(), Speed Governor (GVNR), Swing Equation (Trapezoidal Discretization)

### Community 5 - "Numerical Derivative Methods"
Cohesion: 0.53
Nodes (5): Calculate_Derivative_CFD(), Calculate_Derivatives_CFDD(), Calculate_Derivatives_Lagrange_Interpolator(), Lagrange Polynomial Interpolation for Derivatives, Numerical Derivative Techniques

### Community 7 - "Center of Inertia Dynamics"
Cohesion: 0.67
Nodes (3): Def_Dyn_Fault_COI!(), Def_Dyn_PostF_COI!(), Center of Inertia (COI) Reference Frame

### Community 8 - "ZIP Load and Power Equations"
Cohesion: 0.67
Nodes (3): Def_Dyn_Fault_eqconst_PeQe!(), Def_Dyn_PostF_eqconst_PeQe!(), ZIP Load Model

### Community 9 - "Generator and AVR State Variables"
Cohesion: 1.00
Nodes (3): Define_Dyn_Var_EδPm!, 4th-Order Generator Model (Ed, Eq, Id, Iq), Automatic Voltage Regulator (AVR / Exciter)

### Community 10 - "Contingency and Stability Validation"
Cohesion: 1.00
Nodes (3): Transient Stability Constraints (TSC), Contingency List (9-bus, cross-validated Anatem/Julia), Contingency List with Circuit Mapping (9-bus)

## Knowledge Gaps
- **17 isolated node(s):** `LinearAlgebra`, `SparseArrays`, `Dates`, `NumericIO`, `DataFrames` (+12 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Main_ACOPF_WARMUP (entry point)` connect `AC Power Flow and Network Setup` to `Results Export and Output`, `Post-Fault Dynamics Constraints`, `Fault-Period Dynamics`, `Bus Data and Input Reading`, `Generator and AVR State Variables`?**
  _High betweenness centrality (0.428) - this node is a cross-community bridge._
- **Why does `Def_Dyn_Fault_All!()` connect `Fault-Period Dynamics` to `AC Power Flow and Network Setup`, `Post-Fault Dynamics Constraints`, `Center of Inertia Dynamics`, `ZIP Load and Power Equations`, `Generator and AVR State Variables`, `Contingency and Stability Validation`?**
  _High betweenness centrality (0.161) - this node is a cross-community bridge._
- **Why does `Def_Dyn_PostF_All!()` connect `Post-Fault Dynamics Constraints` to `AC Power Flow and Network Setup`, `Fault-Period Dynamics`, `Center of Inertia Dynamics`, `ZIP Load and Power Equations`, `Generator and AVR State Variables`, `Contingency and Stability Validation`?**
  _High betweenness centrality (0.156) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `Def_Dyn_Fault_All!()` (e.g. with `Automatic Voltage Regulator (AVR / Exciter)` and `Speed Governor (GVNR)`) actually correct?**
  _`Def_Dyn_Fault_All!()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `Def_Dyn_PostF_All!()` (e.g. with `Automatic Voltage Regulator (AVR / Exciter)` and `Speed Governor (GVNR)`) actually correct?**
  _`Def_Dyn_PostF_All!()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 6 inferred relationships involving `Admittance Matrix (Ybus)` (e.g. with `Calculate_Ybus()` and `Calculate_Ybus_fault()`) actually correct?**
  _`Admittance Matrix (Ybus)` has 6 INFERRED edges - model-reasoned connections that need verification._
- **What connects `LinearAlgebra`, `SparseArrays`, `Dates` to the rest of the system?**
  _18 weakly-connected nodes found - possible documentation gaps or missing edges._