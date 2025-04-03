"""
Graph utilities for Davies' method implementation.
"""
from sage.all import Graph, matrix, QQ, identity_matrix, block_matrix

def build_coupling_graph(mechanism):
    """Step 1a: Build the Coupling Graph (GC)"""
    G = Graph(multiedges=False)  # Use simple graph for GC
    # Ensure consistent vertex order for matrix construction
    sorted_bodies = sorted(mechanism.bodies)
    G.add_vertices(sorted_bodies)
    # Assign canonical direction (lower index to higher index body)
    for joint_id, body1, body2 in mechanism.joints:
        u, v = sorted((body1, body2))
        G.add_edge(u, v, label=joint_id)  # Store joint_id as label
    return G

def get_incidence_matrix(graph):
    """Step 1d: Get Incidence Matrix [IC] from a graph"""
    # SageMath's incidence_matrix uses different conventions (+1/-1 reversed, nodes x edges)
    # We need bodies(rows) x couplings(cols) with +1 leaving, -1 entering
    nodes = sorted(graph.vertices())
    edges = sorted(graph.edges(), key=lambda x: x[2])  # Sort by label for consistency
    node_map = {node: i for i, node in enumerate(nodes)}

    n = len(nodes)
    e = len(edges)
    IC = matrix(QQ, n, e)  # Use QQ for exact rational arithmetic

    for j, edge in enumerate(edges):
        u, v, label = edge
        # Use original direction based on how edge was added (u -> v)
        row_u = node_map[u]
        row_v = node_map[v]
        IC[row_u, j] = 1   # Leaving u
        IC[row_v, j] = -1  # Entering v
    return IC, nodes, edges

def get_cutset_matrix(IC):
    """Step 1e: Get Cut-set Matrix [QC] from [IC] via echelon form"""
    # Perform Gaussian elimination (row reduction)
    echelon_form = IC.echelon_form()
    # Remove the last row (guaranteed to be zero for a connected graph)
    QC = echelon_form[:-1, :]
    return QC

def find_pivot_columns(matrix_in_ref):
    """Find pivot columns in a matrix that is in row echelon form."""
    pivots = []
    for row in range(matrix_in_ref.nrows()):
        # Find the first non-zero element in this row
        for col in range(matrix_in_ref.ncols()):
            if matrix_in_ref[row, col] != 0:
                pivots.append(col)
                break
    return pivots

def get_circuit_matrix(QC, num_edges, num_nodes):
    """Step 1f & 1g: Get Circuit Matrix [BC] from [QC] using orthogonality"""
    k = num_nodes - 1  # Number of cuts = number of branches
    l = num_edges - k  # Number of circuits = number of chords

    if QC.ncols() != num_edges or QC.nrows() != k:
        raise ValueError("QC dimensions do not match expected values")

    # Find pivot columns to identify branches and chords implicitly
    pivots = find_pivot_columns(QC)  # Use our new function instead of pivot_columns()
    if len(pivots) != k:
        # This might happen if the graph isn't connected, handle appropriately
        raise ValueError(f"Could not find {k} pivot columns in QC (found {len(pivots)}). Graph might not be connected?")

    non_pivots = sorted(list(set(range(num_edges)) - set(pivots)))
    col_order = list(pivots) + non_pivots  # Branches first, then chords

    # Reorder QC columns to get [[Ub] | [Qc]] form
    QC_reordered = QC.matrix_from_columns(col_order)
    Ub = QC_reordered[:, :k]
    Qc = QC_reordered[:, k:]

    # Check if Ub is identity (it should be if QC is row-reduced correctly)
    if Ub != identity_matrix(QQ, k):
        print("Warning: QC matrix pivots did not form Identity. Check echelon form.")

    # Orthogonality: [Bb] = -[Qc]^T
    Bb = -Qc.transpose()  # Dimensions: l x k

    # Form [BC] = [[Bb] | [Uc]]
    Uc = identity_matrix(QQ, l)  # Dimensions: l x l
    BC_reordered = block_matrix([[Bb, Uc]], subdivide=False)  # Dimensions: l x e

    # Apply the *inverse* column permutation to get BC in the original edge order
    inv_col_order = [0] * num_edges
    for i, original_index in enumerate(col_order):
        inv_col_order[original_index] = i
    BC = BC_reordered.matrix_from_columns(inv_col_order)  # Apply inverse permutation

    return BC

def expand_circuit_matrix_for_motion_graph(BC, mechanism, gc_edges):
    """Step 4: Expand [BC] columns for Motion Graph [BM]"""
    new_cols = []
    gm_edge_map = {}  # Map original joint_id to list of new GM edge indices
    current_gm_col = 0

    for gc_col_idx, gc_edge in enumerate(gc_edges):
        joint_id = gc_edge[2]
        dof = mechanism.joint_dof.get(joint_id, 1)  # Get DOF with default=1
        
        original_col = BC.column(gc_col_idx)
        gm_edge_indices_for_joint = []
        for _ in range(dof):
            new_cols.append(original_col)
            gm_edge_indices_for_joint.append(current_gm_col)
            current_gm_col += 1
        gm_edge_map[joint_id] = gm_edge_indices_for_joint

    # Create BM by concatenating the columns
    if not new_cols:  # Handle case with no edges
        return matrix(QQ, BC.nrows(), 0), {}

    # Construct matrix from columns
    BM = matrix(QQ, BC.nrows(), len(new_cols))
    for i, col in enumerate(new_cols):
        BM.set_column(i, col)

    return BM, gm_edge_map

def expand_cutset_matrix_for_action_graph(QC, mechanism, gc_edges, external_actions=None):
    """Step 3b (Static): Expand [QC] columns for Action Graph [QA]"""
    from common_utils import get_constraints_for_joint_type
    
    if external_actions is None:
        external_actions = {}  # Dict: joint_id -> number of active constraints 'ca'

    num_rows = QC.nrows()
    new_cols = []
    ga_edge_map = {}  # Map original joint_id to list of new GA edge indices (constraints)
    current_ga_col = 0
    total_C = 0  # Total number of constraints (columns in QA)

    for gc_col_idx, gc_edge in enumerate(gc_edges):
        joint_id = gc_edge[2]
        joint_type = mechanism.joint_types[joint_id]
        # Get PASSIVE constraint count
        lambda_dim = 3  # Hardcoding planar for now based on previous example
        cp = get_constraints_for_joint_type(joint_type, lambda_dim)

        # Get ACTIVE constraint count (applied external actions treated as constraints)
        ca = external_actions.get(joint_id, 0)
        c_total_joint = cp + ca

        original_col_vector = QC.column(gc_col_idx)
        ga_edge_indices_for_joint = []
        for _ in range(c_total_joint):
            new_cols.append(original_col_vector)
            ga_edge_indices_for_joint.append(current_ga_col)
            current_ga_col += 1
        ga_edge_map[joint_id] = ga_edge_indices_for_joint  # Store mapping
        total_C += c_total_joint

    if not new_cols:  # Handle case with no edges/constraints
        return matrix(QQ, num_rows, 0), {}, 0

    # Create QA directly using column vectors
    QA = matrix(QQ, num_rows, total_C)
    for col_idx, col_vector in enumerate(new_cols):
        QA.set_column(col_idx, col_vector)

    return QA, ga_edge_map, total_C
