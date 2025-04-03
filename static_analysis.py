"""
Static analysis module using Davies' method.
"""
from sage.all import vector, QQ, RDF, N, matrix, SR  # Added SR
from common_utils import Mechanism, get_constraints_for_joint_type, format_vector
from graph_utils import (
    build_coupling_graph, get_incidence_matrix,
    get_cutset_matrix, get_circuit_matrix,
    expand_cutset_matrix_for_action_graph
)
from matrix_utils import (
    select_independent_rows, partition_system_matrix,
    combine_magnitudes, assemble_network_matrix
)

def print_step(step_name, outputs=None):
    """Print a step with its outputs nicely formatted"""
    print(f"\nStep {step_name}...")
    if outputs:
        for key, value in outputs.items():
            print(f"  → {key}: {value}")

def compute_unit_wrench(joint_id, constraint_index, mechanism, coord_system, external_actions, lambda_dim=3, ring=QQ):
    """Step 4a (Static): Compute Unit Wrench ˆ$^A$ for a specific constraint"""
    joint_type = mechanism.joint_types[joint_id]
    geom = mechanism.geometry.get(joint_id, {})
    point_coords = geom.get('point', [0, 0])  # Planar point Px, Py
    point = vector(ring, point_coords[:2])  # Use specified ring for point
    Px, Py = point[0], point[1]

    if lambda_dim != 3:
        raise NotImplementedError("Only planar (lambda=3) wrenches implemented here")

    # Determine passive constraint count
    cp = get_constraints_for_joint_type(joint_type, lambda_dim)
    ca = external_actions.get(joint_id, 0)

    unit_wrench = vector(ring, [0] * lambda_dim)  # Initialize with specified ring

    # Determine which constraint this index refers to (passive first, then active)
    if constraint_index < cp:  # Passive constraint
        if joint_type == 'revolute':
            if constraint_index == 0:  # Rx constraint
                unit_wrench = vector(ring, [-Py, 1, 0])
            elif constraint_index == 1:  # Ry constraint
                unit_wrench = vector(ring, [Px, 0, 1])
        elif joint_type == 'prismatic':
            dir_coords = geom.get('direction', [1, 0])  # ux, uy
            direction = vector(ring, dir_coords[:2])
            perp_dir = vector(ring, [-direction[1], direction[0]])  # vx, vy = -uy, ux
            if constraint_index == 0:  # Force constraint Rv along perp_dir
                unit_wrench = vector(ring, [Px * perp_dir[1] - Py * perp_dir[0], perp_dir[0], perp_dir[1]])
            elif constraint_index == 1:  # Torque constraint Tz
                unit_wrench = vector(ring, [1, 0, 0])
        else:
            print(f"Warning: Passive unit wrench for joint type '{joint_type}' index {constraint_index} not implemented.")
    else:  # Active constraint (external action treated as constraint)
        active_constraint_idx = constraint_index - cp
        if active_constraint_idx == 0:  # Assume first active is Tz
            unit_wrench = vector(ring, [1, 0, 0])
        else:
            print(f"Warning: Active unit wrench index {active_constraint_idx} for joint '{joint_id}' not implemented.")

    return unit_wrench

def assemble_unit_wrench_matrix(mechanism, gc_edges, ga_edge_map, total_C, external_actions, lambda_dim=3, ring=QQ):
    """Step 4b (Static): Assemble Unit Wrench Matrix [ÂD]"""
    A_hat_D = matrix(ring, lambda_dim, total_C)  # Use specified ring
    unit_wrenches_list = []

    coord_system = {'origin': vector(ring, [0] * lambda_dim)}  # Use specified ring

    # Ensure column order matches ga_edge_map
    ga_col_to_joint_constraint = {}
    for joint_id, indices in ga_edge_map.items():
        for i, col_idx in enumerate(indices):
            ga_col_to_joint_constraint[col_idx] = (joint_id, i)  # Map col to joint and constraint index

    for col_idx in range(total_C):
        joint_id, constraint_index = ga_col_to_joint_constraint[col_idx]
        unit_wrench = compute_unit_wrench(joint_id, constraint_index, mechanism, coord_system, external_actions, lambda_dim, ring)
        A_hat_D[:, col_idx] = unit_wrench
        unit_wrenches_list.append(unit_wrench)

    return A_hat_D, unit_wrenches_list

def format_static_report(report, detailed=False, precision=4):
    """Format the static analysis report dictionary into a readable string."""
    if not report:
        return "No report available."

    lines = ["Davies' Static Analysis Report", "=" * 40]

    # Mechanism information
    mechanism = report.get('mechanism')
    if mechanism:
        lines.append("\n## Mechanism Information")
        lines.append(f"- Bodies: {sorted(mechanism.bodies)}")
        lines.append(f"- Joints: {[j[0] for j in mechanism.joints]}")

        lines.append("\n### Joint Details:")
        for j_id, b1, b2 in mechanism.joints:
            j_type = mechanism.joint_types.get(j_id, 'N/A')
            j_dof = mechanism.joint_dof.get(j_id, 'N/A')
            geom_str = ", ".join([f"{k}: {format_vector(v, precision)}" if hasattr(v,'__iter__') else f"{k}: {v}"
                                for k, v in mechanism.geometry.get(j_id, {}).items()])
            lines.append(f"- {j_id} ({j_type}, f={j_dof}): {b1} <-> {b2} [{geom_str}]")

    # External actions
    external_actions = report.get('external_actions')
    if external_actions:
        lines.append("\n## External Actions")
        for joint_id, num_actions in external_actions.items():
            lines.append(f"- Joint {joint_id}: {num_actions} active constraint(s)")

    # Graph information
    gc = report.get('coupling_graph')
    ga_edges_ordered = report.get('ga_edges_ordered')
    if gc:
        lines.append("\n## Graph Information")
        lines.append(f"- GC: {gc.order()} nodes, {gc.size()} edges")
    if ga_edges_ordered:
        lines.append(f"- GA: {len(ga_edges_ordered)} constraints (edges)")

    # Matrix Dimensions
    lines.append("\n## Matrix Dimensions")
    if 'incidence_matrix' in report: lines.append(f"- IC: {report['incidence_matrix'].dimensions()}")
    if 'cutset_matrix_gc' in report: lines.append(f"- QC (GC): {report['cutset_matrix_gc'].dimensions()}")
    if 'cutset_matrix_ga' in report: lines.append(f"- QA (GA): {report['cutset_matrix_ga'].dimensions()}")
    if 'unit_wrench_matrix' in report: lines.append(f"- ÂD: {report['unit_wrench_matrix'].dimensions()}")
    if 'network_action_matrix' in report: lines.append(f"- ÂN: {report['network_action_matrix'].dimensions()}")
    if 'network_action_matrix_reduced' in report: lines.append(f"- ÂN_reduced: {report['network_action_matrix_reduced'].dimensions()}")
    if 'A_N_S' in report: lines.append(f"- ÂNS: {report['A_N_S'].dimensions()}")
    if 'A_N_P' in report: lines.append(f"- ÂNP: {report['A_N_P'].dimensions()}")

    # Unit Wrenches
    unit_wrenches = report.get('unit_wrenches_ordered')
    if unit_wrenches:
        lines.append("\n## Unit Wrenches (ÂD Columns)")
        ga_edges = report.get('ga_edges_ordered', [])
        for i, wrench in enumerate(unit_wrenches):
            edge_label = ga_edges[i] if i < len(ga_edges) else f"Constraint_{i}"
            lines.append(f"- {edge_label:<12}: {format_vector(wrench, precision)}")

    # Network Action Matrix
    a_hat_n = report.get('network_action_matrix')
    if a_hat_n is not None and detailed:
        lines.append("\n## Network Action Matrix (ÂN)")
        lines.append(f"{a_hat_n}")

    # Results
    psi_total = report.get('Psi_total')
    ga_edges = report.get('ga_edges_ordered')
    if psi_total is not None and ga_edges is not None:
        lines.append("\n## Final Wrench Magnitudes {Ψ_total}")
        for i, magnitude in enumerate(psi_total):
            edge_label = ga_edges[i] if i < len(ga_edges) else f"Constraint_{i}"
            unit = "N" if ("Rx" in edge_label or "Ry" in edge_label or "Rv" in edge_label) else "Nm"
            lines.append(f"- {edge_label:<12}: {format_vector(vector([magnitude]), precision+1)} {unit}")

    return "\n".join(lines)

def save_static_report_to_file(report, filename, detailed=True, precision=6):
    """Save the formatted static report to a file."""
    with open(filename, 'w') as f:
        f.write(format_static_report(report, detailed, precision))
    print(f"Static report saved to {filename}")

def davies_static_analysis(mechanism, primary_constraint_ids, primary_values, external_actions=None, lambda_dim=3, generate_report=False, ring=QQ):
    """
    Main function implementing Davies' method for static analysis.

    Args:
        mechanism: A Mechanism object.
        primary_constraint_ids: List of IDs for the primary known actions/reactions.
                                These need a clear mapping to specific constraints
                                (e.g., 'a_Tz', 'b_Rx', 'd_Tout').
        primary_values: List/Vector of corresponding primary force/torque magnitudes.
        external_actions: Dict mapping joint_id -> number of active constraints ('ca').
                          Defaults to 0 if None. Used for expanding QA.
        lambda_dim: Dimension of the workspace (e.g., 3 for planar, 6 for spatial).
        generate_report: Whether to generate a detailed report dictionary.
        ring: The SageMath ring (QQ, RDF, SR) for calculations. Defaults to QQ.

    Returns:
        Psi_total: Vector of all constraint variable magnitudes (forces/torques).
                   Order corresponds to the columns defined by GA_Edges_Ordered.
        GA_Edges_Ordered: List of GA constraint identifiers corresponding to Psi_total.
        report: Dictionary containing detailed information about the analysis (if generate_report=True).
    """
    print(f"\n=== Davies Static Analysis (Ring: {ring}) ===\n")
    report = {} if generate_report else None
    
    if external_actions is None:
        external_actions = {}

    if generate_report:
        report['mechanism'] = mechanism
        report['external_actions'] = external_actions
        report['calculation_ring'] = ring  # Store the ring used

    # Step 1 & 2 (GC, IC, QC): Same as kinematics
    GC = build_coupling_graph(mechanism)
    IC, gc_nodes, gc_edges = get_incidence_matrix(GC)
    num_gc_nodes = len(gc_nodes)
    num_gc_edges = len(gc_edges)
    QC = get_cutset_matrix(IC, ring=ring)
    print_step("1 & 2: Building GC, IC, QC", {
        "Graph": f"GC: {num_gc_nodes} nodes, {num_gc_edges} edges",
        "Matrix": f"QC dimensions: {QC.dimensions()}"
    })
    
    if generate_report:
        report['coupling_graph'] = GC
        report['incidence_matrix'] = IC
        report['cutset_matrix_gc'] = QC
        report['gc_edges_ordered'] = gc_edges

    # Step 3: Expand [QC] to Action Graph Cut-set Matrix [QA]
    QA, ga_edge_map, total_C = expand_cutset_matrix_for_action_graph(QC, mechanism, gc_edges, external_actions, ring=ring)
    num_cuts = QA.nrows()  # k
    print_step("3: Expanding [QC] to [QA] for Action Graph", {
        "Dimensions": f"QA dimensions: {num_cuts} cuts, {total_C} total constraints (C)"
    })
    
    if generate_report:
        report['cutset_matrix_ga'] = QA
        report['ga_edge_map'] = ga_edge_map

    # Create ordered list of GA 'edges' (constraints)
    GA_Edges_Ordered = [None] * total_C
    coord_system = {'origin': vector(ring, [0] * lambda_dim)}  # Use specified ring
    for joint_id, indices in ga_edge_map.items():
         cp = get_constraints_for_joint_type(mechanism.joint_types[joint_id], lambda_dim)
         ca = external_actions.get(joint_id, 0)
         for i, col_idx in enumerate(indices):
              constraint_name = f"{joint_id}_constraint_{i}"  # Placeholder name
              if i < cp:  # Passive
                  if mechanism.joint_types[joint_id] == 'revolute':
                       if i == 0: constraint_name = f"{joint_id}_Rx"
                       elif i == 1: constraint_name = f"{joint_id}_Ry"
                  elif mechanism.joint_types[joint_id] == 'prismatic':
                       if i == 0: constraint_name = f"{joint_id}_Rv"  # Perpendicular force
                       elif i == 1: constraint_name = f"{joint_id}_Tz"
              else:  # Active
                  constraint_name = f"{joint_id}_Tz_active"  # Assume Tz for active constraints
                  
              GA_Edges_Ordered[col_idx] = constraint_name
    
    if generate_report:
        report['ga_edges_ordered'] = GA_Edges_Ordered

    # Step 4: Assemble Unit Wrench Matrix [ÂD]
    A_hat_D, unit_wrenches_ordered = assemble_unit_wrench_matrix(
        mechanism, gc_edges, ga_edge_map, total_C, external_actions, lambda_dim, ring=ring
    )
    print_step("4: Assembling Unit Wrench Matrix [ÂD]", {
        "Dimensions": f"ÂD dimensions: {A_hat_D.dimensions()}"
    })
    
    if generate_report:
        report['unit_wrench_matrix'] = A_hat_D
        report['unit_wrenches_ordered'] = unit_wrenches_ordered

    # Step 5: Assemble Network Action Matrix [ÂN]
    A_hat_N = assemble_network_matrix(A_hat_D, QA, lambda_dim, ring=ring)
    print_step("5: Assembling Network Action Matrix [ÂN]", {
        "Dimensions": f"ÂN dimensions: {A_hat_N.dimensions()}"
    })
    
    if generate_report:
        report['network_action_matrix'] = A_hat_N

    # Step 6: Handle Sub-Restriction (Rank Reduction)
    A_hat_N_reduced, rank_a = select_independent_rows(A_hat_N, ring=ring)
    print_step("6: Reducing ÂN to independent rows", {
        "Rank": f"Rank a = {rank_a}",
        "Reduced Dimensions": f"Reduced ÂN dimensions: {A_hat_N_reduced.dimensions()}"
    })
    
    if generate_report:
        report['network_action_matrix_reduced'] = A_hat_N_reduced
        report['rank_a'] = rank_a
    
    # Check for static indeterminacy
    FN_static = lambda_dim * num_cuts - rank_a
    if FN_static > 0:
        print(f"Warning: System appears statically indeterminate (FN = {FN_static}). May need more constraints or primary inputs.")
    elif FN_static < 0:
         print(f"Warning: System appears overconstrained based on rank calculation (FN = {FN_static}).")

    # Step 7: Partition System
    primary_indices = []
    primary_value_map = {}
    if len(primary_constraint_ids) != len(primary_values):
         raise ValueError("Mismatch between number of primary constraint IDs and primary values")

    # Create mapping from descriptive ID to column index
    constraint_id_to_col = {name: idx for idx, name in enumerate(GA_Edges_Ordered)}

    primary_values_ring = [ring(v) for v in primary_values]

    for constraint_id, value in zip(primary_constraint_ids, primary_values_ring):
         if constraint_id not in constraint_id_to_col:
              print(f"Warning: Primary constraint ID '{constraint_id}' not found in generated constraints: {GA_Edges_Ordered}")
              continue
         col_idx = constraint_id_to_col[constraint_id]
         primary_indices.append(col_idx)
         primary_value_map[col_idx] = value

    CN_provided = len(primary_indices)
    print_step("7: System Partitioning", {
        "Primary": f"Identified {CN_provided} primary constraint(s) (columns): {primary_indices}"
    })

    try:
        A_hat_N_S, A_hat_N_P, secondary_indices = partition_system_matrix(
            A_hat_N_reduced, total_C, primary_indices
        )
        print(f"  → Secondary: Secondary constraint indices (columns): {secondary_indices}")
        print(f"  → Matrix: ÂNS {A_hat_N_S.dimensions()}, ÂNP {A_hat_N_P.dimensions()}")
        
        if generate_report:
            report['primary_indices'] = primary_indices 
            report['secondary_indices'] = secondary_indices
            report['A_N_S'] = A_hat_N_S
            report['A_N_P'] = A_hat_N_P
            
    except ValueError as e:
        print(f"Error during partitioning: {e}")
        print(f"Expected {A_hat_N_reduced.nrows()} secondary variables, but got {total_C - len(primary_indices)}.")
        print("Check mechanism definition, primary variable selection, and matrix ranks.")
        return None, None, report if generate_report else None
    except Exception as e:
        print(f"An unexpected error occurred during partitioning: {e}")
        return None, None, report if generate_report else None

    # Step 8: Solve for Secondary Variables {ΨS}
    Psi_P = vector(ring, [primary_value_map[i] for i in primary_indices])

    try:
        rhs = -A_hat_N_P * Psi_P
        Psi_S = A_hat_N_S.solve_right(rhs)
        
        if generate_report:
            report['Psi_P'] = Psi_P
            report['Psi_S'] = Psi_S
            
    except Exception as e:
        print(f"Error solving system: {e}")
        print("ÂNS might be singular. Check mechanism constraints or primary variable choice.")
        return None, None, report if generate_report else None

    # Step 9: Combine Results {Ψ}
    Psi_total = combine_magnitudes(Psi_P, Psi_S, primary_indices, secondary_indices, total_C, ring=ring)
    
    if generate_report:
        report['Psi_total'] = Psi_total

    print("\n=== Davies Static Analysis Complete ===\n")
    return (Psi_total, GA_Edges_Ordered, report) if generate_report else (Psi_total, GA_Edges_Ordered)