"""
Kinematic analysis module using Davies' method.
"""
from sage.all import vector, QQ, N, matrix, SR  # Added SR
from common_utils import Mechanism, format_vector
from graph_utils import (
    build_coupling_graph, get_incidence_matrix, 
    get_cutset_matrix, get_circuit_matrix,
    expand_circuit_matrix_for_motion_graph
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

def compute_unit_twist(joint_id, dof_index, mechanism, coord_system, lambda_dim=3, ring=QQ):
    """Step 6: Compute Unit Twist ˆ$^M$ for a specific DoF"""
    joint_type = mechanism.joint_types[joint_id]
    geom = mechanism.geometry.get(joint_id, {})
    point_coords = geom.get('point', [0, 0, 0])
    point = vector(ring, point_coords[:lambda_dim])  # Use specified ring

    if lambda_dim != 3:
        raise NotImplementedError("Only planar (lambda=3) twists implemented")

    unit_twist = vector(ring, [0] * lambda_dim)  # Initialize with specified ring

    if joint_type == 'revolute' and dof_index == 0:
        # ωz=1; Vp = [-Py*1; Px*1]
        unit_twist = vector(ring, [1, -point[1], point[0]])
    elif joint_type == 'prismatic' and dof_index == 0:
        dir_coords = geom.get('direction', [1, 0, 0])
        direction = vector(ring, dir_coords[:lambda_dim])
        # ωz=0; Vp = direction
        unit_twist = vector(ring, [0, direction[0], direction[1]])
    elif joint_type == 'planar':
        if dof_index == 0:  # Rotation about P
            unit_twist = vector(ring, [1, -point[1], point[0]])
        elif dof_index == 1:  # Translation X
            unit_twist = vector(ring, [0, 1, 0])
        elif dof_index == 2:  # Translation Y
            unit_twist = vector(ring, [0, 0, 1])
    
    return unit_twist

def assemble_unit_twist_matrix(mechanism, gc_edges, gm_edge_map, lambda_dim=3, ring=QQ):
    """Step 7: Assemble Unit Twist Matrix [M̂D]"""
    F = sum(dof for dof in mechanism.joint_dof.values())  # Total DoF = columns
    M_hat_D = matrix(ring, lambda_dim, F)  # Use specified ring
    unit_twists_list = []  # Keep track of the order

    coord_system = {'origin': vector(ring, [0] * lambda_dim)}  # Use specified ring

    gm_col_to_joint_dof = {}
    for joint_id, indices in gm_edge_map.items():
        for i, col_idx in enumerate(indices):
            gm_col_to_joint_dof[col_idx] = (joint_id, i)  # Map col to joint and DoF index

    for col_idx in range(F):
         joint_id, dof_index = gm_col_to_joint_dof[col_idx]
         unit_twist = compute_unit_twist(joint_id, dof_index, mechanism, coord_system, lambda_dim, ring)
         M_hat_D[:, col_idx] = unit_twist
         unit_twists_list.append(unit_twist)  # Store in correct order

    return M_hat_D, unit_twists_list

def map_primary_inputs_to_indices(primary_joint_ids, primary_values, gm_edge_map, mechanism):
    """Helper to map primary joint IDs/values to column indices and values"""
    primary_indices = []
    primary_value_map = {}  # Store value for each primary *column index*
    processed_base_ids = set()

    if len(primary_joint_ids) != len(primary_values):
        raise ValueError("Mismatch between number of primary joint IDs and primary values")

    for base_joint_id, value in zip(primary_joint_ids, primary_values):
        if base_joint_id not in gm_edge_map:
            raise ValueError(f"Primary joint ID '{base_joint_id}' not found in mechanism.")
        if base_joint_id in processed_base_ids:
            raise ValueError(f"Primary joint ID '{base_joint_id}' specified multiple times.")

        gm_cols = gm_edge_map[base_joint_id]
        primary_col_index = gm_cols[0]
        primary_indices.append(primary_col_index)
        primary_value_map[primary_col_index] = value
        processed_base_ids.add(base_joint_id)
        for i in range(1, len(gm_cols)):
            print(f"Warning: Joint {base_joint_id} has multiple DoFs. Assuming DoF {i} is implicitly primary with value 0.")
            implicit_primary_col_index = gm_cols[i]
            primary_indices.append(implicit_primary_col_index)
            primary_value_map[implicit_primary_col_index] = 0  # Assign 0 value

    return sorted(list(set(primary_indices))), primary_value_map  # Ensure unique indices

def format_kinematic_report(report, detailed=False, precision=4):
    """Format the report dictionary into a readable string."""
    if not report:
        return "No report available."

    lines = ["Davies' Kinematic Analysis Report", "=" * 40]

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

    gc = report.get('coupling_graph')
    gm_edges_ordered = report.get('gm_edges_ordered')
    if gc:
        lines.append("\n## Graph Information")
        lines.append(f"- GC: {gc.order()} nodes, {gc.size()} edges")
    if gm_edges_ordered:
        lines.append(f"- GM: {len(gm_edges_ordered)} DoFs (edges)")

    lines.append("\n## Matrix Dimensions")
    if 'incidence_matrix' in report: lines.append(f"- IC: {report['incidence_matrix'].dimensions()}")
    if 'cutset_matrix_gc' in report: lines.append(f"- QC (GC): {report['cutset_matrix_gc'].dimensions()}")
    if 'circuit_matrix_gc' in report: lines.append(f"- BC (GC): {report['circuit_matrix_gc'].dimensions()}")
    if 'circuit_matrix_gm' in report: lines.append(f"- BM (GM): {report['circuit_matrix_gm'].dimensions()}")
    if 'unit_twist_matrix' in report: lines.append(f"- M̂D: {report['unit_twist_matrix'].dimensions()}")
    if 'network_motion_matrix' in report: lines.append(f"- M̂N: {report['network_motion_matrix'].dimensions()}")
    if 'network_motion_matrix_reduced' in report: lines.append(f"- M̂N_reduced: {report['network_motion_matrix_reduced'].dimensions()}")
    if 'M_N_S' in report: lines.append(f"- M̂NS: {report['M_N_S'].dimensions()}")
    if 'M_N_P' in report: lines.append(f"- M̂NP: {report['M_N_P'].dimensions()}")

    circuits = report.get('circuit_matrix_gm')
    if circuits is not None:
        lines.append(f"\n## Circuit Information (BM - {circuits.nrows()} Circuits)")
        if detailed:
            lines.append(f"{circuits}")

    unit_twists = report.get('unit_twists_ordered')
    if unit_twists:
        lines.append("\n## Unit Twists (M̂D Columns)")
        gm_edges = report.get('gm_edges_ordered', [])
        for i, twist in enumerate(unit_twists):
            edge_label = gm_edges[i][2] if i < len(gm_edges) else f"DoF_{i}"
            lines.append(f"- {edge_label:<8}: {format_vector(twist, precision)}")

    m_hat_n = report.get('network_motion_matrix')
    if m_hat_n is not None and detailed:
        lines.append("\n## Network Motion Matrix (M̂N)")
        lines.append(f"{m_hat_n}")

    phi_total = report.get('Phi_total')
    gm_edges = report.get('gm_edges_ordered')
    if phi_total is not None and gm_edges is not None:
        lines.append("\n## Final Twist Magnitudes {Φ_total}")
        for i, magnitude in enumerate(phi_total):
            edge_label = gm_edges[i][2] if i < len(gm_edges) else f"DoF_{i}"
            lines.append(f"- {edge_label:<10}: {format_vector(vector([magnitude]), precision+1)}")

    return "\n".join(lines)

def save_report_to_file(report, filename, detailed=True, precision=6):
    """Save the formatted report to a file."""
    with open(filename, 'w') as f:
        f.write(format_kinematic_report(report, detailed, precision))
    print(f"Report saved to {filename}")

def davies_kinematic_analysis(mechanism, primary_joint_ids, primary_values, lambda_dim=3, generate_report=False, ring=QQ):
    """
    Main function implementing Davies' method for kinematic analysis.
    
    Args:
        mechanism: A Mechanism object.
        primary_joint_ids: List of IDs for the primary known joint velocities.
        primary_values: List/Vector of corresponding primary velocity magnitudes.
        lambda_dim: Dimension of the workspace (e.g., 3 for planar, 6 for spatial).
        generate_report: Whether to generate a detailed report dictionary.
        ring: The SageMath ring (QQ, RDF, SR) for calculations. Defaults to QQ.
        
    Returns:
        Phi_total: Vector of all twist magnitudes (velocities) over the specified ring.
        GM_Edges_Ordered: List of GM edge identifiers corresponding to Phi_total.
        report: Dictionary containing detailed information about the analysis (if generate_report=True).
    """
    print(f"\n=== Davies Kinematic Analysis (Ring: {ring}) ===\n")
    report = {} if generate_report else None

    if generate_report:
        report['mechanism'] = mechanism
        report['calculation_ring'] = ring  # Store the ring used

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

    BC = get_circuit_matrix(QC, num_gc_edges, num_gc_nodes, ring=ring)
    print_step("3: Calculating Circuit Matrix [BC]", {
        "Dimensions": f"BC dimensions: {BC.dimensions()}"
    })
    if generate_report:
        report['circuit_matrix_gc'] = BC

    BM, gm_edge_map = expand_circuit_matrix_for_motion_graph(BC, mechanism, gc_edges, ring=ring)
    F = BM.ncols()  # Total DoF
    num_circuits = BM.nrows()
    print_step("4: Expanding [BC] to [BM] for Motion Graph", {
        "Dimensions": f"BM dimensions: {num_circuits} circuits, {F} total DoF (F)"
    })
    if generate_report:
        report['circuit_matrix_gm'] = BM
        report['gm_edge_map'] = gm_edge_map

    GM_Edges_Ordered = [None] * F
    for joint_id, indices in gm_edge_map.items():
        dof = mechanism.joint_dof[joint_id]
        gc_edge = next((edge for edge in gc_edges if edge[2] == joint_id), None)
        if gc_edge is None: 
            raise Exception(f"Cannot find GC edge for joint {joint_id}")

        for i, col_idx in enumerate(indices):
             gm_edge_name = f"{joint_id}" if dof == 1 else f"{joint_id}_{i}"
             GM_Edges_Ordered[col_idx] = (gc_edge[0], gc_edge[1], gm_edge_name)
    if generate_report:
        report['gm_edges_ordered'] = GM_Edges_Ordered

    coord_system = {'origin': vector(ring, [0] * lambda_dim)}

    M_hat_D, unit_twists_ordered = assemble_unit_twist_matrix(mechanism, gc_edges, gm_edge_map, lambda_dim, ring=ring)
    print_step("6 & 7: Assembling Unit Twist Matrix [M̂D]", {
        "Dimensions": f"M̂D dimensions: {M_hat_D.dimensions()}"
    })
    if generate_report:
        report['unit_twist_matrix'] = M_hat_D
        report['unit_twists_ordered'] = unit_twists_ordered

    M_hat_N = assemble_network_matrix(M_hat_D, BM, lambda_dim, ring=ring)
    print_step("8: Assembling Network Motion Matrix [M̂N]", {
        "Dimensions": f"M̂N dimensions: {M_hat_N.dimensions()}"
    })
    if generate_report:
        report['network_motion_matrix'] = M_hat_N

    M_hat_N_reduced, rank_m = select_independent_rows(M_hat_N, ring=ring)
    print_step("9: Reducing M̂N to independent rows", {
        "Rank": f"m = {rank_m}",
        "Matrix": f"Reduced M̂N dimensions: {M_hat_N_reduced.dimensions()}"
    })
    if generate_report:
        report['network_motion_matrix_reduced'] = M_hat_N_reduced
        report['rank_m'] = rank_m

    try:
        primary_indices, primary_value_map = map_primary_inputs_to_indices(
            primary_joint_ids, primary_values, gm_edge_map, mechanism
        )
        primary_value_map_ring = {k: ring(v) for k, v in primary_value_map.items()}

        FN_actual = len(primary_indices)
        FN_expected = F - rank_m
        if FN_actual != FN_expected:
             print(f"Warning: Number of specified primary DoFs ({FN_actual}) does not match expected free DoFs ({FN_expected}).")

        M_hat_N_S, M_hat_N_P, secondary_indices = partition_system_matrix(
            M_hat_N_reduced, F, primary_indices
        )
        print_step("10: Partitioning the system", {
            "Primary": f"{len(primary_indices)} primary variable(s): {primary_indices}",
            "Secondary": f"{len(secondary_indices)} secondary variable(s): {secondary_indices}"
        })
        if generate_report:
            report['primary_indices'] = primary_indices
            report['secondary_indices'] = secondary_indices
            report['M_N_S'] = M_hat_N_S
            report['M_N_P'] = M_hat_N_P
    except ValueError as e:
        print(f"Error during partitioning: {e}")
        return None, None, report if generate_report else None
    except Exception as e:
        print(f"An unexpected error occurred during partitioning: {e}")
        return None, None, report if generate_report else None

    Phi_P = vector(ring, [primary_value_map_ring[i] for i in primary_indices])
    try:
        rhs = -M_hat_N_P * Phi_P
        Phi_S = M_hat_N_S.solve_right(rhs)
    except Exception as e:
        print(f"Error solving system: {e}")
        print("M̂NS might be singular. Check mechanism constraints or primary variable choice.")
        return None, None, report if generate_report else None
    print_step("11: Solving for secondary variables {ΦS}")
    if generate_report:
        report['Phi_P'] = Phi_P
        report['Phi_S'] = Phi_S

    Phi_total = combine_magnitudes(Phi_P, Phi_S, primary_indices, secondary_indices, F, ring=ring)
    print_step("12: Combining primary and secondary variables")
    if generate_report:
        report['Phi_total'] = Phi_total

    return (Phi_total, GM_Edges_Ordered, report) if generate_report else (Phi_total, GM_Edges_Ordered)
