"""
Matrix utilities for Davies' method implementation.
"""
from sage.all import matrix, vector, QQ, RDF, block_matrix, diagonal_matrix

def select_independent_rows(A):
    """Step 9: Select linearly independent rows"""
    # Ensure matrix is over QQ or RDF for rank/pivot methods
    A_field = A.change_ring(RDF) if A.base_ring() != RDF else A
    try:
        rank = A_field.rank()
        pivot_rows = A_field.pivot_rows()
        # Convert back to original ring if needed, selecting original rows
        A_reduced = A.matrix_from_rows(pivot_rows)
        if A_reduced.nrows() != rank:
            print(f"Warning: Number of pivot rows ({len(pivot_rows)}) != rank ({rank}). Using first {rank} rows.")
            A_reduced = A.matrix_from_rows(list(range(rank)))
    except Exception as e:
        print(f"Error during rank/pivot calculation: {e}. Attempting rank via RDF.")
        rank = A.change_ring(RDF).rank()
        print(f"Using rank {rank}. Selecting first {rank} rows as fallback.")
        A_reduced = A.matrix_from_rows(list(range(rank)))

    return A_reduced, rank

def partition_system_matrix(A_reduced, num_total_vars, primary_indices):
    """Step 10: Partition the reduced system matrix"""
    m, F = A_reduced.dimensions()
    if F != num_total_vars:
        raise ValueError("Number of columns in reduced matrix doesn't match total variables")

    secondary_indices = sorted(list(set(range(F)) - set(primary_indices)))

    # Check if the number of secondary variables matches the rank (number of equations)
    if len(secondary_indices) != m:
        raise ValueError(f"Number of secondary variables ({len(secondary_indices)}) does not match system rank ({m}). "
                         f"FN ({F-m}) != specified primary ({len(primary_indices)}). Check primary variable selection.")

    A_S = A_reduced.matrix_from_columns(secondary_indices)
    A_P = A_reduced.matrix_from_columns(primary_indices)

    return A_S, A_P, secondary_indices

def combine_magnitudes(P_primary, P_secondary, primary_indices, secondary_indices, total_size):
    """Step 12: Combine results into the full vector"""
    P_total = vector(QQ, total_size)  # Use QQ for consistency
    for i, idx in enumerate(primary_indices):
        P_total[idx] = P_primary[i]
    for i, idx in enumerate(secondary_indices):
        P_total[idx] = P_secondary[i]
    return P_total

def assemble_network_matrix(D_hat, B_matrix, lambda_dim):
    """Assemble Network Matrix (Motion or Action)"""
    rows, cols = B_matrix.dimensions()
    if cols == 0:  
        return matrix(QQ, lambda_dim * rows, 0)  # Handle case with no edges

    blocks = []
    for i in range(rows):
        row_i = B_matrix.row(i)
        # Create diagonal matrix from the row's elements
        diag_i = diagonal_matrix(QQ, row_i.list())  # Size cols x cols
        if diag_i.dimensions() != (cols, cols):
            raise ValueError(f"Diagonal matrix creation failed for row {i}. "
                            f"Expected size {(cols, cols)}, got {diag_i.dimensions()}.")
        # Check dimensions before multiplication
        if D_hat.ncols() != diag_i.nrows():
            raise ValueError(f"Dimension mismatch: D_hat ({D_hat.dimensions()}) and diag_i ({diag_i.dimensions()})")
        block_i = D_hat * diag_i  # Result lambda x cols
        blocks.append(block_i)

    # Stack the blocks vertically
    if not blocks:  # Handle case with no rows
        return matrix(QQ, 0, cols)
    N_hat = block_matrix(blocks, nrows=rows, subdivide=False)
    return N_hat
