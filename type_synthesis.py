"""
Type synthesis module using Davies' method and matroid theory.

This module implements automated type synthesis algorithms that use Davies' method
to systematically enumerate and select self-aligning mechanism types.

Joint Types in Planar Mechanisms (lambda=3):
-----------------------------------------
- revolute:       1 DOF, rotation about Z-axis
- prismatic_x:    1 DOF, translation along X-axis
- prismatic_y:    1 DOF, translation along Y-axis
- pin_in_slot_x:  2 DOF, rotation + translation along X-axis
- pin_in_slot_y:  2 DOF, rotation + translation along Y-axis
- planar_trans:   2 DOF, translation along X and Y axes
- planar:         3 DOF, rotation + translation along X and Y axes
- rigid:          0 DOF, fully constrained

Joint Types in Spatial Mechanisms (lambda=6):
-----------------------------------------
- revolute_x:        1 DOF, rotation about X-axis
- revolute_y:        1 DOF, rotation about Y-axis
- revolute_z:        1 DOF, rotation about Z-axis
- prismatic_x:       1 DOF, translation along X-axis
- prismatic_y:       1 DOF, translation along Y-axis
- prismatic_z:       1 DOF, translation along Z-axis
- planar_translation_xy: 2 DOF, translation in XY plane
- planar_translation_xz: 2 DOF, translation in XZ plane
- planar_translation_yz: 2 DOF, translation in YZ plane
- spherical:         3 DOF, rotation about X, Y, and Z axes
- spatial_translation: 3 DOF, translation along X, Y, and Z axes
- rigid:             0 DOF, fully constrained

Additional custom joint types may be identified during type synthesis,
labeled with their degrees of freedom as "spatial_dofN_custom".
"""
from sage.all import matrix, vector, QQ, SR, MatrixSpace
from sage.matroids.constructor import Matroid
from sage.combinat.subset import Subsets
import itertools

from common_utils import Mechanism, get_constraints_for_joint_type
from static_analysis import davies_static_analysis
from graph_utils import expand_cutset_matrix_for_action_graph, build_coupling_graph, get_incidence_matrix, get_cutset_matrix
from matrix_utils import select_independent_rows, assemble_network_matrix

def print_step(step_name, outputs=None):
    """Print a step with its outputs nicely formatted."""
    print(f"\n=== {step_name} ===")
    if outputs:
        for k, v in outputs.items():
            print(f"  {k}: {v}")

# Step 2: Create Seed Mechanism Model using Reshetov Virtual Joints
def create_seed_mechanism_using_reshetov_virtual_joints(mechanism, lambda_dim):
    """
    Create a seed mechanism where each joint initially imposes all λ constraints.
    
    Args:
        mechanism: A Mechanism object with the connectivity topology.
        lambda_dim: Dimension of the motion space (e.g., 3 for planar, 6 for spatial).
        
    Returns:
        A new Mechanism object with Reshetov virtual joints.
    """
    seed_mechanism = Mechanism()
    
    # Copy all bodies
    for body_id in mechanism.bodies:
        seed_mechanism.add_body(body_id)
    
    # Replace each conceptual joint with a Reshetov virtual joint (all constraints active)
    for joint_id, body1, body2 in mechanism.joints:
        geometry = mechanism.geometry.get(joint_id, {})
        
        # For a Reshetov virtual joint, DOF is 0 (fully constrained) and joint_type is 'virtual'
        seed_mechanism.add_joint(joint_id, body1, body2, "virtual", 0, geometry)
    
    print(f"Created seed mechanism with {len(seed_mechanism.bodies)} bodies and {len(seed_mechanism.joints)} virtual joints.")
    return seed_mechanism

# Step 3: Compute Network Unit Action Matrix
def compute_network_unit_action_matrix(seed_mechanism, lambda_dim, ring=QQ):
    """
    Use the static analysis machinery to compute the Network Unit Action Matrix [AN].
    
    Args:
        seed_mechanism: A Mechanism object with Reshetov virtual joints.
        lambda_dim: Dimension of the motion space.
        ring: The SageMath ring for calculations.
        
    Returns:
        The Network Unit Action Matrix [AN] and related data.
    """
    # Build the coupling graph and its matrices
    GC = build_coupling_graph(seed_mechanism)
    IC, gc_nodes, gc_edges = get_incidence_matrix(GC)
    QC = get_cutset_matrix(IC, ring=ring)
    
    # External actions are not needed for type synthesis
    external_actions = {}
    
    # Expand QC to Action Graph cutset matrix QA
    QA, ga_edge_map, total_C = expand_cutset_matrix_for_action_graph(
        QC, seed_mechanism, gc_edges, external_actions, ring=ring
    )
    
    # Create ordered list of constraints for each joint
    GA_Edges_Ordered = [None] * total_C
    for joint_id, indices in ga_edge_map.items():
        for i, col_idx in enumerate(indices):
            # For virtual joints, we'll name constraints based on their freedom type
            # This helps with interpretation later
            if lambda_dim == 3:  # Planar
                if i == 0: constraint_name = f"{joint_id}_Tx"
                elif i == 1: constraint_name = f"{joint_id}_Ty"
                elif i == 2: constraint_name = f"{joint_id}_Rz"
            elif lambda_dim == 6:  # Spatial
                if i == 0: constraint_name = f"{joint_id}_Tx"
                elif i == 1: constraint_name = f"{joint_id}_Ty"
                elif i == 2: constraint_name = f"{joint_id}_Tz"
                elif i == 3: constraint_name = f"{joint_id}_Rx"
                elif i == 4: constraint_name = f"{joint_id}_Ry"
                elif i == 5: constraint_name = f"{joint_id}_Rz"
            GA_Edges_Ordered[col_idx] = constraint_name
    
    # Set up coordinate system and compute unit wrenches
    coord_system = {'origin': vector(ring, [0] * lambda_dim)}
    
    # Build the unit wrench matrix manually (we can't use the static analysis directly
    # because our virtual joints are a special case)
    A_hat_D = matrix(ring, lambda_dim, total_C)
    
    # For each constraint, set up the appropriate unit wrench
    for col_idx in range(total_C):
        constraint_name = GA_Edges_Ordered[col_idx]
        joint_id = constraint_name.split('_')[0]
        constraint_type = constraint_name.split('_')[1]
        
        # Get joint geometry
        geom = seed_mechanism.geometry.get(joint_id, {})
        point_coords = geom.get('point', [0, 0, 0])
        
        # Create unit wrench based on constraint type
        # For planar (lambda=3)
        if lambda_dim == 3:
            # Point coordinates
            Px, Py = point_coords[0], point_coords[1]
            
            if constraint_type == "Tx":  # Constraint on x translation
                A_hat_D[:, col_idx] = vector(ring, [0, 1, -Py])
            elif constraint_type == "Ty":  # Constraint on y translation
                A_hat_D[:, col_idx] = vector(ring, [0, 0, Px])
            elif constraint_type == "Rz":  # Constraint on z rotation
                A_hat_D[:, col_idx] = vector(ring, [1, 0, 0])
        
        # For spatial (lambda=6), we'd need to add more cases here
        # This is a simplification for demonstration
        elif lambda_dim == 6:
            Px, Py, Pz = point_coords[0], point_coords[1], point_coords[2]
            
            # This is a simplified placeholder - actual spatial wrenches would be more complex
            if constraint_type == "Tx":
                A_hat_D[:, col_idx] = vector(ring, [0, 1, 0, 0, 0, -Pz])
            # ... additional cases for other constraint types ...
    
    # Assemble the network unit action matrix
    A_hat_N = assemble_network_matrix(A_hat_D, QA, lambda_dim, ring=ring)
    
    print_step("Network Unit Action Matrix Computed", {
        "Dimensions": f"AN: {A_hat_N.dimensions()}",
        "Constraint Count": total_C
    })
    
    return A_hat_N, GA_Edges_Ordered, gc_edges, ga_edge_map

# Step 4: Define the Initial Matroid
def define_linear_matroid(A_hat_N):
    """
    Define a linear matroid M = (E, I) on the columns of [AN].
    
    Args:
        A_hat_N: The Network Unit Action Matrix.
        
    Returns:
        A SageMath Matroid object.
    """
    # Create a linear matroid using the columns of A_hat_N
    matroid = Matroid(matrix=A_hat_N)
    
    print_step("Linear Matroid Defined", {
        "Ground Set Size": len(matroid.groundset()),
        "Rank": matroid.rank()
    })
    
    return matroid

# Step 5: Translate Design Requirements into Matroid Operations
def determine_deletions_from_design_requirements(design_requirements, mechanism, ground_set, constraint_names):
    """
    Identify constraints to be removed based on design requirements.
    
    Args:
        design_requirements: Dictionary of design requirements.
        mechanism: The original Mechanism object.
        ground_set: The ground set of the matroid.
        constraint_names: Ordered list of constraint names.
        
    Returns:
        Set of elements to be deleted from the matroid.
    """
    required_deletions = set()
    
    # If no design requirements, return empty set
    if not design_requirements or 'required_freedoms' not in design_requirements:
        return required_deletions
    
    # Map constraint names to their index in the ground set
    constraint_to_idx = {name: idx for idx, name in enumerate(constraint_names) if name is not None}
    
    # Process required freedoms
    for joint_id, freedoms in design_requirements.get('required_freedoms', {}).items():
        for freedom in freedoms:
            # For each freedom, identify the corresponding constraint to remove
            constraint_to_remove = f"{joint_id}_{freedom}"
            if constraint_to_remove in constraint_to_idx:
                required_deletions.add(constraint_to_idx[constraint_to_remove])
    
    print_step("Design Requirements - Deletions", {
        "Number of Required Deletions": len(required_deletions)
    })
    
    return required_deletions

def determine_contractions_from_design_requirements(design_requirements, mechanism, ground_set, constraint_names):
    """
    Identify constraints to be preserved based on design requirements.
    
    Args:
        design_requirements: Dictionary of design requirements.
        mechanism: The original Mechanism object.
        ground_set: The ground set of the matroid.
        constraint_names: Ordered list of constraint names.
        
    Returns:
        Set of elements to be contracted in the matroid.
    """
    required_contractions = set()
    
    # If no design requirements, return empty set
    if not design_requirements or 'required_constraints' not in design_requirements:
        return required_contractions
    
    # Map constraint names to their index in the ground set
    constraint_to_idx = {name: idx for idx, name in enumerate(constraint_names) if name is not None}
    
    # Process required constraints
    for joint_id, constraints in design_requirements.get('required_constraints', {}).items():
        for constraint in constraints:
            # For each required constraint, identify its index
            constraint_to_keep = f"{joint_id}_{constraint}"
            if constraint_to_keep in constraint_to_idx:
                required_contractions.add(constraint_to_idx[constraint_to_keep])
    
    print_step("Design Requirements - Contractions", {
        "Number of Required Contractions": len(required_contractions)
    })
    
    return required_contractions

# Step 6: Apply Matroid Contraction and Deletion
def apply_matroid_contraction_and_deletion(matroid, contractions, deletions):
    """
    Apply contraction and deletion operations to the matroid.
    
    Args:
        matroid: The original Matroid object.
        contractions: Set of elements to contract.
        deletions: Set of elements to delete.
        
    Returns:
        A new filtered Matroid object.
    """
    # Start with the original matroid
    filtered_matroid = matroid
    
    # Apply contractions first (if any)
    if contractions:
        filtered_matroid = filtered_matroid.contract(contractions)
    
    # Then apply deletions (if any)
    if deletions:
        filtered_matroid = filtered_matroid.delete(deletions)
    
    print_step("Matroid Operations Applied", {
        "Original Ground Set Size": len(matroid.groundset()),
        "Filtered Ground Set Size": len(filtered_matroid.groundset()),
        "Filtered Matroid Rank": filtered_matroid.rank()
    })
    
    return filtered_matroid

# Step 7: Enumerate Bases of the Filtered Matroid
def enumerate_all_bases(matroid):
    """
    Enumerate all bases of the matroid.
    
    Args:
        matroid: A Matroid object.
        
    Returns:
        List of all bases (each basis is a set of elements).
    """
    # Use SageMath's built-in basis enumeration
    bases = list(matroid.bases())
    
    print_step("Bases Enumeration", {
        "Number of Bases Found": len(bases)
    })
    
    return bases

# Step 8: Reconstruct and Interpret Candidate Mechanisms
def reconstruct_constraint_set(basis, contractions, ground_set):
    """
    Reconstruct the complete set of constraints for a candidate.
    
    Args:
        basis: A basis of the filtered matroid.
        contractions: Set of contracted elements.
        ground_set: The original ground set.
        
    Returns:
        The complete set of constraints for this candidate.
    """
    # The constraints in the mechanism are those in the basis plus the contracted elements
    constraint_set = set(basis).union(contractions)
    
    print_step("Constraint Set Reconstruction", {
        "Constraints in Basis": len(basis),
        "Contracted Constraints": len(contractions),
        "Total Constraints": len(constraint_set)
    })
    
    return constraint_set

def determine_joint_types_from_constraint_set(constraint_set, mechanism, constraint_names, lambda_dim):
    """
    Determine the joint types based on the constraint set.
    
    Args:
        constraint_set: Set of constraints in the mechanism.
        mechanism: The original Mechanism object.
        constraint_names: Ordered list of constraint names.
        lambda_dim: Dimension of the motion space.
        
    Returns:
        Dictionary mapping joint IDs to their determined types.
    """
    joint_types = {}
    joint_constraints = {}
    
    # Group constraints by joint
    for idx in constraint_set:
        if idx < len(constraint_names) and constraint_names[idx]:
            parts = constraint_names[idx].split('_')
            joint_id = parts[0]
            constraint_type = parts[1]
            
            if joint_id not in joint_constraints:
                joint_constraints[joint_id] = set()
            joint_constraints[joint_id].add(constraint_type)
    
    # Determine joint type based on constraints
    for joint_id, constraints in joint_constraints.items():
        dof = lambda_dim - len(constraints)
        
        # For planar mechanisms (lambda=3)
        if lambda_dim == 3:
            if dof == 1:  # 1 DOF joints
                if 'Tx' in constraints and 'Ty' in constraints:  # Only rotation allowed
                    joint_types[joint_id] = "revolute"
                elif 'Tx' in constraints and 'Rz' in constraints:  # Only Y translation allowed
                    joint_types[joint_id] = "prismatic_y"
                elif 'Ty' in constraints and 'Rz' in constraints:  # Only X translation allowed
                    joint_types[joint_id] = "prismatic_x"
                else:
                    joint_types[joint_id] = f"unknown_dof{dof}"
            elif dof == 2:  # 2 DOF joints
                if 'Tx' in constraints:  # Y translation and rotation allowed
                    # This is a pin-in-slot joint with Y-direction slot
                    joint_types[joint_id] = "pin_in_slot_y"
                elif 'Ty' in constraints:  # X translation and rotation allowed
                    # This is a pin-in-slot joint with X-direction slot
                    joint_types[joint_id] = "pin_in_slot_x"
                elif 'Rz' in constraints:  # X and Y translation allowed
                    joint_types[joint_id] = "planar_trans"
                else:
                    joint_types[joint_id] = f"unknown_dof{dof}"
            elif dof == 3:
                joint_types[joint_id] = "planar"
            elif dof == 0:
                joint_types[joint_id] = "rigid"
            else:
                joint_types[joint_id] = f"unknown_dof{dof}"
        
        # For spatial mechanisms (lambda=6)
        elif lambda_dim == 6:
            # Joint type determination for spatial mechanisms
            if dof == 1:  # 1-DOF joints in spatial mechanisms
                # Determine specific 1-DOF joint types based on constraints
                # (revolute, prismatic, etc.)
                if len({'Tx', 'Ty', 'Tz', 'Rx', 'Ry'} & constraints) == 5:
                    joint_types[joint_id] = "revolute_z"
                elif len({'Tx', 'Ty', 'Tz', 'Rx', 'Rz'} & constraints) == 5:
                    joint_types[joint_id] = "revolute_y"
                elif len({'Tx', 'Ty', 'Tz', 'Ry', 'Rz'} & constraints) == 5:
                    joint_types[joint_id] = "revolute_x"
                elif len({'Rx', 'Ry', 'Rz', 'Ty', 'Tz'} & constraints) == 5:
                    joint_types[joint_id] = "prismatic_x"
                elif len({'Rx', 'Ry', 'Rz', 'Tx', 'Tz'} & constraints) == 5:
                    joint_types[joint_id] = "prismatic_y"
                elif len({'Rx', 'Ry', 'Rz', 'Tx', 'Ty'} & constraints) == 5:
                    joint_types[joint_id] = "prismatic_z"
                else:
                    joint_types[joint_id] = f"spatial_dof1_custom"
            elif dof == 2:
                # Determine specific 2-DOF joint types
                if 'Rz' in constraints and 'Ry' in constraints and 'Rx' in constraints:
                    # All rotations constrained - planar translation
                    if 'Tz' in constraints:
                        joint_types[joint_id] = "planar_translation_xy"
                    elif 'Ty' in constraints:
                        joint_types[joint_id] = "planar_translation_xz"
                    elif 'Tx' in constraints:
                        joint_types[joint_id] = "planar_translation_yz"
                    else:
                        joint_types[joint_id] = "spatial_dof2_custom"
                else:
                    # Other 2-DOF combinations
                    joint_types[joint_id] = "spatial_dof2_custom"
            elif dof == 3:
                # 3-DOF joints like planar, spherical
                if 'Rx' in constraints and 'Ry' in constraints and 'Rz' in constraints:
                    # All rotations constrained
                    joint_types[joint_id] = "spatial_translation"
                elif 'Tx' in constraints and 'Ty' in constraints and 'Tz' in constraints:
                    # All translations constrained
                    joint_types[joint_id] = "spherical"
                else:
                    joint_types[joint_id] = "spatial_dof3_custom"
            elif dof == 0:
                joint_types[joint_id] = "rigid"
            else:
                joint_types[joint_id] = f"spatial_dof{dof}_custom"
    
    print_step("Joint Type Determination", {
        "Number of Joints": len(joint_types)
    })
    
    return joint_types

# Step 9: Refine Selection
def refine_candidate_mechanisms(feasible_mechanisms, design_requirements):
    """
    Apply further filtering to remove isomorphic or undesired solutions.
    
    Args:
        feasible_mechanisms: List of feasible mechanism candidates.
        design_requirements: Dictionary of design requirements.
        
    Returns:
        Filtered list of feasible mechanisms.
    """
    # This is a placeholder for more sophisticated filtering
    # For now, we'll just remove potential duplicates based on joint types
    
    # Create a unique representation for each mechanism based on sorted joint types
    unique_mechanisms = {}
    for i, (joint_types, constraint_set) in enumerate(feasible_mechanisms):
        # Create a sortable representation of the joint types
        joint_type_repr = tuple(sorted((j, t) for j, t in joint_types.items()))
        
        if joint_type_repr not in unique_mechanisms:
            unique_mechanisms[joint_type_repr] = (joint_types, constraint_set)
    
    refined_mechanisms = list(unique_mechanisms.values())
    
    print_step("Mechanism Refinement", {
        "Original Count": len(feasible_mechanisms),
        "Refined Count": len(refined_mechanisms)
    })
    
    return refined_mechanisms

# Main function
def automate_type_synthesis(mechanism, design_requirements=None, lambda_dim=3, ring=QQ, max_bases=1000, is_gripper=False):
    """
    Systematically enumerate and select self-aligning mechanism types.
    
    Args:
        mechanism: A Mechanism object describing the seed mechanism topology.
        design_requirements: Dict of design criteria (can be None).
            Format: {
                'required_freedoms': {joint_id: [freedom_types]},
                'required_constraints': {joint_id: [constraint_types]}
            }
            Example: {
                'required_freedoms': {'a': ['Rz'], 'b': ['Tx', 'Ty']}
                'required_constraints': {'c': ['Tx', 'Ty', 'Rz']}
            }
        lambda_dim: Dimension of the motion space (e.g., 3 for planar, 6 for spatial).
        ring: The SageMath ring for calculations (QQ, RDF, SR).
        max_bases: Maximum number of bases to enumerate (for computational efficiency).
        is_gripper: Flag to indicate if this is a gripper mechanism (for reporting).
        
    Returns:
        List of synthesized self-aligning mechanism types, each as a tuple:
            (joint_types_dict, constraint_set)
    """
    print("\n=== AUTOMATED TYPE SYNTHESIS ===\n")
    
    # Step 1: Extract information from the mechanism
    n = len(mechanism.bodies)
    j = len(mechanism.joints)
    
    print_step("Initial Mechanism", {
        "Number of Links": n,
        "Number of Joints": j,
        "Motion Space Dimension": lambda_dim
    })
    
    # Step 2: Create Seed Mechanism with Reshetov Virtual Joints
    seed_mechanism = create_seed_mechanism_using_reshetov_virtual_joints(mechanism, lambda_dim)
    
    # Step 3: Compute Network Unit Action Matrix via Davies' Method
    A_hat_N, constraint_names, gc_edges, ga_edge_map = compute_network_unit_action_matrix(
        seed_mechanism, lambda_dim, ring
    )
    
    # Step 4: Define the Initial Matroid
    matroid = define_linear_matroid(A_hat_N)
    
    # Step 5: Translate Design Requirements into Matroid Operations
    ground_set = list(matroid.groundset())
    if design_requirements:
        required_deletions = determine_deletions_from_design_requirements(
            design_requirements, mechanism, ground_set, constraint_names
        )
        required_contractions = determine_contractions_from_design_requirements(
            design_requirements, mechanism, ground_set, constraint_names
        )
    else:
        required_deletions = set()
        required_contractions = set()
    
    # Step 6: Apply Matroid Contraction and Deletion
    filtered_matroid = apply_matroid_contraction_and_deletion(
        matroid, required_contractions, required_deletions
    )
    
    # Step 7: Enumerate Bases of the Filtered Matroid (limit to max_bases)
    all_bases = enumerate_all_bases(filtered_matroid)
    if len(all_bases) > max_bases:
        print(f"Warning: Limiting to {max_bases} bases out of {len(all_bases)} total.")
        bases = all_bases[:max_bases]
    else:
        bases = all_bases
    
    # Step 8: Reconstruct and Interpret Candidate Mechanisms
    feasible_mechanisms = []
    for basis in bases:
        constraint_set = reconstruct_constraint_set(basis, required_contractions, ground_set)
        joint_types = determine_joint_types_from_constraint_set(
            constraint_set, mechanism, constraint_names, lambda_dim
        )
        feasible_mechanisms.append((joint_types, constraint_set))
    
    # Step 9: Refine Selection (Optional)
    refined_mechanisms = refine_candidate_mechanisms(feasible_mechanisms, design_requirements)
    
    print("\n=== TYPE SYNTHESIS COMPLETE ===\n")
    print(f"Found {len(refined_mechanisms)} unique mechanism types.")
    
    # Print a summary of the results
    for i, (joint_types, _) in enumerate(refined_mechanisms[:10]):  # Show at most 10 results
        print(f"\nMechanism Type #{i+1}:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    
    if len(refined_mechanisms) > 10:
        print(f"\n... and {len(refined_mechanisms) - 10} more mechanism types.")
    
    return refined_mechanisms

# Function to save results to file
def save_type_synthesis_results(results, filename, mechanism=None, design_requirements=None, lambda_dim=3, is_gripper=False):
    """
    Save type synthesis results to a file with improved formatting and application guidance.
    
    Args:
        results: List of feasible mechanism types.
        filename: Path to output file.
        mechanism: Original Mechanism object (optional).
        design_requirements: Design requirements used in synthesis (optional).
        lambda_dim: Dimension of motion space.
        is_gripper: Whether to include gripper-specific recommendations.
    """
    with open(filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("                      TYPE SYNTHESIS RESULTS REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # Write summary information
        f.write("SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total unique mechanism types found: {len(results)}\n")
        f.write(f"Motion space dimension: {lambda_dim} ({'Planar' if lambda_dim == 3 else 'Spatial'})\n")
        
        # Include information about the original mechanism if provided
        if mechanism:
            f.write(f"Original mechanism: {len(mechanism.bodies)} links, {len(mechanism.joints)} joints\n")
            joint_list = ", ".join([j[0] for j in mechanism.joints])
            f.write(f"Joints: {joint_list}\n")
        
        # Include design requirements if provided
        if design_requirements:
            f.write("\nDESIGN REQUIREMENTS\n")
            f.write("-" * 40 + "\n")
            
            if 'required_freedoms' in design_requirements:
                f.write("Required freedoms:\n")
                for joint_id, freedoms in design_requirements['required_freedoms'].items():
                    f.write(f"  Joint {joint_id}: {', '.join(freedoms)}\n")
            
            if 'required_constraints' in design_requirements:
                f.write("Required constraints:\n")
                for joint_id, constraints in design_requirements['required_constraints'].items():
                    f.write(f"  Joint {joint_id}: {', '.join(constraints)}\n")
        
        # Write joint type reference table
        f.write("\nJOINT TYPE REFERENCE\n")
        f.write("-" * 40 + "\n")
        if lambda_dim == 3:  # Planar
            f.write("revolute:       1 DOF, rotation about Z-axis\n")
            f.write("prismatic_x:    1 DOF, translation along X-axis\n")
            f.write("prismatic_y:    1 DOF, translation along Y-axis\n")
            f.write("pin_in_slot_x:  2 DOF, rotation + translation along X-axis\n")
            f.write("pin_in_slot_y:  2 DOF, rotation + translation along Y-axis\n")
            f.write("planar_trans:   2 DOF, translation along X and Y axes\n")
            f.write("planar:         3 DOF, rotation + translation along X and Y axes\n")
            f.write("rigid:          0 DOF, fully constrained\n")
        else:  # Spatial
            f.write("revolute_x/y/z:    1 DOF, rotation about specified axis\n")
            f.write("prismatic_x/y/z:   1 DOF, translation along specified axis\n")
            f.write("spherical:         3 DOF, rotation about X, Y, and Z axes\n")
            f.write("planar:            3 DOF, translation in plane + rotation about normal\n")
            f.write("spatial_dofN:      N DOF, various combinations\n")
        
        # Write detailed results for each mechanism
        f.write("\nDETAILED MECHANISM TYPES\n")
        f.write("=" * 80 + "\n\n")
        
        for i, (joint_types, constraint_set) in enumerate(results):
            f.write(f"MECHANISM TYPE #{i+1}\n")
            f.write("-" * 40 + "\n")
            
            # Create a table of joint types
            f.write("Joint Types:\n")
            for joint_id, joint_type in sorted(joint_types.items()):
                dof = 0
                if "revolute" in joint_type or "prismatic" in joint_type and "_" not in joint_type:
                    dof = 1
                elif "pin_in_slot" in joint_type:
                    dof = 2
                elif "planar" in joint_type:
                    dof = 3
                
                f.write(f"  {joint_id:<5} : {joint_type:<15} ({dof} DOF)\n")
            
            # Compute mobility
            if mechanism:
                mobility = sum(dof for joint_id, joint_type in joint_types.items() 
                              if "rigid" not in joint_type)
                mobility -= (lambda_dim * (len(mechanism.bodies) - 1))
                f.write(f"\nMobility: {mobility}\n")
            
            # Gripper-specific recommendations if requested
            if is_gripper:
                f.write("\nGRIPPER RECOMMENDATIONS:\n")
                
                # Identify potential contact points (revolute joints are good candidates)
                contact_candidates = [j for j, t in joint_types.items() if "revolute" in t]
                f.write(f"Potential contact points: {', '.join(contact_candidates) if contact_candidates else 'None'}\n")
                
                # Identify potential actuation points (prismatic joints are good candidates)
                actuation_candidates = [j for j, t in joint_types.items() if "prismatic" in t]
                f.write(f"Potential actuation joints: {', '.join(actuation_candidates) if actuation_candidates else 'None'}\n")
                
                # Identify potential fixed joints (rigid connections)
                fixed_joints = [j for j, t in joint_types.items() if "rigid" in t]
                f.write(f"Fixed connections: {', '.join(fixed_joints) if fixed_joints else 'None'}\n")
                
                # Provide general guidance
                f.write("\nGuidance:\n")
                f.write("- For grippers, contact points should typically be revolute joints\n")
                f.write("- Linear actuation often uses prismatic joints\n")
                f.write("- Pin-in-slot joints can provide helpful degrees of freedom for self-alignment\n")
            
            f.write("\n" + "=" * 80 + "\n\n")
        
        f.write("END OF REPORT\n")
    
    print(f"Detailed results saved to {filename}")
