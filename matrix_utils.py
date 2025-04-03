"""
Matrix utilities for Davies' method implementation.
"""
from sage.all import matrix, vector, QQ, RDF, SR, block_matrix, diagonal_matrix

def select_independent_rows(A, ring=QQ):
    """Step 9: Select linearly independent rows"""
    # Ensure matrix is over the appropriate ring for rank/pivot methods
    # For symbolic, rank/pivot work directly on SR
    A_field = A if ring is SR else A.change_ring(RDF)

    try:
        rank = A_field.rank()
        # Pivot rows might not be reliable/meaningful for SR in the same way
        # For SR, we often want to keep all rows unless rank deficiency is proven
        if ring is SR:
            # For symbolic, assume full rank unless proven otherwise, or handle rank deficiency later
            # Alternatively, one could try symbolic simplification and row reduction,
            # but let's keep it simple first and return the original matrix rows if rank matches.
            if rank == A.nrows():
                pivot_rows = list(range(A.nrows()))
                A_reduced = A  # Keep original symbolic matrix if rank matches rows
            else:
                # If rank deficient symbolically, more complex handling is needed.
                # For now, let's try echelon form, but this can be complex symbolically.
                print(f"Warning: Symbolic matrix rank ({rank}) < nrows ({A.nrows()}). Attempting echelon form.")
                try:
                    A_reduced = A.echelon_form()
                    # Recalculate rank based on non-zero rows of echelon form
                    rank = sum(1 for r in range(A_reduced.nrows()) if not A_reduced.row(r).is_zero())
                    A_reduced = A_reduced.matrix_from_rows(list(range(rank)))
                    print(f"Using {rank} rows from symbolic echelon form.")
                except Exception as e_ech:
                    print(f"Symbolic echelon form failed: {e_ech}. Falling back to first {rank} rows.")
                    A_reduced = A.matrix_from_rows(list(range(rank)))

        else:  # Numerical case (QQ or RDF)
            pivot_rows = A_field.pivot_rows()
            # Convert back to original ring if needed, selecting original rows
            A_reduced = A.matrix_from_rows(pivot_rows)
            if A_reduced.nrows() != rank:
                print(f"Warning: Number of pivot rows ({len(pivot_rows)}) != rank ({rank}). Using first {rank} rows.")
                A_reduced = A.matrix_from_rows(list(range(rank)))

    except Exception as e:
        print(f"Error during rank/pivot calculation: {e}. Attempting rank via RDF.")
        # Fallback for numerical might still be needed
        rank_rdf = A.change_ring(RDF).rank()
        print(f"Using rank {rank_rdf}. Selecting first {rank_rdf} rows as fallback.")
        A_reduced = A.matrix_from_rows(list(range(rank_rdf)))
        rank = rank_rdf  # Update rank based on fallback

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

def combine_magnitudes(P_primary, P_secondary, primary_indices, secondary_indices, total_size, ring=QQ):
    """Step 12: Combine results into the full vector"""
    P_total = vector(ring, total_size)  # Use specified ring
    for i, idx in enumerate(primary_indices):
        P_total[idx] = P_primary[i]
    for i, idx in enumerate(secondary_indices):
        P_total[idx] = P_secondary[i]
    return P_total

def assemble_network_matrix(D_hat, B_matrix, lambda_dim, ring=QQ):
    """Assemble Network Matrix (Motion or Action)"""
    rows, cols = B_matrix.dimensions()
    if cols == 0:
        return matrix(ring, lambda_dim * rows, 0)  # Use specified ring

    # Ensure B_matrix is over the correct ring before creating diagonal matrix
    B_matrix_ring = B_matrix.change_ring(ring)

    blocks = []
    for i in range(rows):
        row_i = B_matrix_ring.row(i)
        # Create diagonal matrix from the row's elements over the specified ring
        diag_i = diagonal_matrix(ring, row_i.list())  # Size cols x cols
        if diag_i.dimensions() != (cols, cols):
            raise ValueError(f"Diagonal matrix creation failed for row {i}. "
                             f"Expected size {(cols, cols)}, got {diag_i.dimensions()}.")
        # Ensure D_hat is also over the correct ring before multiplication
        D_hat_ring = D_hat.change_ring(ring)
        if D_hat_ring.ncols() != diag_i.nrows():
            raise ValueError(f"Dimension mismatch: D_hat ({D_hat_ring.dimensions()}) and diag_i ({diag_i.dimensions()})")
        block_i = D_hat_ring * diag_i  # Result lambda x cols
        blocks.append(block_i)

    # Stack the blocks vertically
    if not blocks:  # Handle case with no rows
        return matrix(ring, 0, cols)
    N_hat = block_matrix(ring, blocks, nrows=rows, subdivide=False)  # Use specified ring
    return N_hat
