"""Example demonstrating both kinematic and static analysis using Davies' method."""
import os
from sage.all import vector, QQ, RDF, pi, cos, sin, N, SR, var
from davies_method import (
    Mechanism, DaviesKinematicAnalysis, DaviesStaticAnalysis,
    save_report_to_file, save_static_report_to_file, format_report, format_static_report
)

# Create reports directory if it doesn't exist
REPORTS_DIR = "reports"
if not os.path.exists(REPORTS_DIR):
    os.makedirs(REPORTS_DIR)

def print_section_header(title, char='='):
    """Print a section header with consistent formatting"""
    print(f"\n{title}")
    print(char * len(title))
    print()

def print_step_header(title):
    """Print a step header with consistent formatting"""
    print(f"\n--- {title} ---")

def user_verified_fourbar_example():
    """Run verified four-bar linkage example"""
    print_section_header("Verified Four-Bar Linkage Example")
    
    # Create mechanism for kinematic analysis
    print("Creating mechanism...")
    kin_mechanism = Mechanism()
    
    # Add bodies with verified geometry for kinematics
    kin_mechanism.add_body("ground")
    kin_mechanism.add_body("link_ab")
    kin_mechanism.add_body("link_bc")
    kin_mechanism.add_body("link_cd")
    
    # Add joints with original verified geometry
    kin_mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0, 0, 0])
    })
    
    kin_mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0, 0.35, 0])
    })
    
    kin_mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0.45, 0.15, 0])
    })
    
    kin_mechanism.add_joint("d", "link_cd", "ground", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0.40, 0, 0])
    })
    
    print(f"Kinematic mechanism created with {len(kin_mechanism.bodies)} bodies and {len(kin_mechanism.joints)} joints.")
    
    # ===== Kinematic Analysis =====
    print_step_header("Running Kinematic Analysis")
    
    # Use verified input values
    primary_joint_ids = ['d']  # Input at joint d
    input_speed = QQ(15)/100  # 0.15 rad/s verified input
    primary_vals = [input_speed]
    
    # Run kinematic analysis
    Phi_solution, ordered_gm_edges, kin_report = DaviesKinematicAnalysis(
        kin_mechanism, primary_joint_ids, primary_vals, lambda_dim=3, generate_report=True
    )
    
    # Print kinematic results with expected values
    if Phi_solution is not None:
        expected = {
            "a": QQ(31)/420,    # ≈ 0.07380952
            "b": -QQ(2)/35,     # ≈ -0.05714286
            "c": QQ(2)/15,      # ≈ 0.13333333
            "d": QQ(3)/20       # = 0.15 (input)
        }
        
        print("\nKinematic Analysis Results (with verification):")
        for i, val in enumerate(Phi_solution):
            edge_label = ordered_gm_edges[i][2]
            val_rdf = N(val, digits=6)
            expected_val = N(expected[edge_label], digits=6)
            diff = abs(val_rdf - expected_val)
            print(f"  ω_{edge_label:<5} = {val_rdf:>8.4f} rad/s (expected: {expected_val:>8.4f}, diff: {diff:.6f})")
        
        save_report_to_file(kin_report, os.path.join(REPORTS_DIR, "verified_fourbar_kinematic_report.txt"))
    else:
        print("Kinematic analysis failed.")
    
    # ===== Static Analysis =====
    print_step_header("Running Static Analysis")
    
    # Create a new mechanism with test case geometry
    print("Creating mechanism for static analysis...")
    static_mechanism = Mechanism()
    
    # Add joints with test_static_analysis.py verified geometry
    static_mechanism.add_joint("a", 0, 1, "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0, 0, 0])
    })
    
    static_mechanism.add_joint("b", 1, 2, "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([2, 0, 0])
    })
    
    static_mechanism.add_joint("c", 2, 3, "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([4, 2, 0])
    })
    
    static_mechanism.add_joint("d", 3, 0, "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([0, 2, 0])
    })
    
    print(f"Static mechanism created with {len(static_mechanism.bodies)} bodies and {len(static_mechanism.joints)} joints.")
    
    # Use verified static test values
    input_torque = QQ(10)  # 10 Nm input torque
    external_actions = {'a': 1, 'd': 1}  # Input torque at a, reaction at d
    primary_constraints = ['a_Tz_active']  # Specify input torque as primary
    primary_constraint_values = [input_torque]
    
    # Run static analysis with proper geometry
    Psi_solution, ordered_ga_constraints, static_report = DaviesStaticAnalysis(
        static_mechanism,
        primary_constraints,
        primary_constraint_values,
        external_actions,
        lambda_dim=3,
        generate_report=True
    )
    
    # Print static results (without expected values)
    if Psi_solution is not None:
        print("\nStatic Analysis Results:")
        for i, val in enumerate(Psi_solution):
            constraint_id = ordered_ga_constraints[i]
            unit = "N" if ("Rx" in constraint_id or "Ry" in constraint_id or "Rv" in constraint_id) else "Nm"
            val_rdf = N(val, digits=6)
            print(f"  {constraint_id:<12} = {val_rdf:>10.4f} {unit}")
        
        save_static_report_to_file(static_report, os.path.join(REPORTS_DIR, "verified_fourbar_static_report.txt"))
    else:
        print("Static analysis failed.")
    
    return kin_mechanism, Phi_solution, Psi_solution

def example_slider_crank():
    """Run verified slider-crank example"""
    print_section_header("Slider-Crank Mechanism Example")
    
    print("Creating mechanism...")
    mechanism = Mechanism()
    
    # Add bodies
    mechanism.add_body(0)  # Ground
    mechanism.add_body(1)  # Crank
    mechanism.add_body(2)  # Connecting Rod
    mechanism.add_body(3)  # Slider
    
    # Define geometry
    crank_length = QQ(2)
    slider_pos_x = QQ(5)
    
    # Add joints with verified geometry
    mechanism.add_joint("A", 0, 1, "revolute", 1, {
        'axis': vector(QQ, [0, 0, 1]),
        'point': vector(QQ, [0, 0, 0])
    })
    
    mechanism.add_joint("B", 1, 2, "revolute", 1, {
        'axis': vector(QQ, [0, 0, 1]),
        'point': vector(QQ, [crank_length, 0, 0])
    })
    
    mechanism.add_joint("C", 2, 3, "revolute", 1, {
        'axis': vector(QQ, [0, 0, 1]),
        'point': vector(QQ, [slider_pos_x, 0, 0])
    })
    
    mechanism.add_joint("D", 3, 0, "prismatic", 1, {
        'direction': vector(QQ, [1, 0, 0])
    })
    
    print(f"Mechanism created with {len(mechanism.bodies)} bodies and {len(mechanism.joints)} joints.")
    
    # ===== Kinematic Analysis =====
    print_step_header("Running Kinematic Analysis")
    primary_joint_ids = ['A']
    input_speed = QQ(1)
    primary_vals = [input_speed]
    
    Phi_solution, ordered_gm_edges, kin_report = DaviesKinematicAnalysis(
        mechanism, primary_joint_ids, primary_vals, lambda_dim=3, generate_report=True
    )
    
    if Phi_solution is not None:
        print("\nKinematic Analysis Results:")
        for i, val in enumerate(Phi_solution):
            edge_label = ordered_gm_edges[i][2]
            val_rdf = N(val, digits=6)
            is_revolute = edge_label in ['A', 'B', 'C']
            prefix = 'ω_' if is_revolute else 'v_'
            unit = 'rad/s' if is_revolute else 'm/s'
            print(f"  {prefix}{edge_label:<5} = {val_rdf:>8.4f} {unit}")
        
        save_report_to_file(kin_report, os.path.join(REPORTS_DIR, "slider_crank_kinematic_report.txt"))
    else:
        print("Kinematic analysis failed.")
    
    return mechanism, Phi_solution, None  # Static analysis not verified for slider-crank

def example_symbolic_fourbar():
    """Run four-bar linkage example with symbolic dimensions/input."""
    print_section_header("Symbolic Four-Bar Linkage Example")

    # Define symbolic variables
    var('L_ab, L_bc, L_cd, L_ad_x, L_ad_y, theta_a, omega_a')

    # Create mechanism with symbolic geometry
    print("Creating symbolic mechanism...")
    mechanism = Mechanism()
    mechanism.add_body("ground")
    mechanism.add_body("link_ab")
    mechanism.add_body("link_bc")
    mechanism.add_body("link_cd")

    # Joint positions defined symbolically
    # Assume joint 'a' is at origin for simplicity here
    # Joint 'd' is fixed relative to 'a'
    # Joint 'b' depends on L_ab and theta_a (implicitly via analysis)
    # Joint 'c' depends on other links (implicitly via analysis)
    # NOTE: Davies method uses fixed geometry. For symbolic *positions*,
    # a different approach (like loop closure equations) is usually needed.
    # Here, we'll make the *lengths* symbolic, but define an initial numeric pose
    # for calculating twists/wrenches, and use a symbolic input velocity.

    # Let's use the verified four-bar geometry but make input velocity symbolic
    # Define numeric parts using integers/fractions
    ax_num, ay_num = 0, 0
    # Create fractions using division rather than QQ(num, denom)
    by_frac = QQ(35)/QQ(100)  # 0.35
    cx_frac = QQ(45)/QQ(100)  # 0.45
    cy_frac = QQ(15)/QQ(100)  # 0.15
    dx_frac = QQ(40)/QQ(100)  # 0.40
    
    # Create vectors directly over SR, letting it coerce numeric types
    mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'point': vector(SR, [0, 0, 0])
    })
    mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'point': vector(SR, [0, by_frac, 0])
    })
    mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'point': vector(SR, [cx_frac, cy_frac, 0])
    })
    mechanism.add_joint("d", "link_cd", "ground", "revolute", 1, {
        'point': vector(SR, [dx_frac, 0, 0])
    })

    print(f"Symbolic mechanism created with {len(mechanism.bodies)} bodies and {len(mechanism.joints)} joints.")

    # ===== Symbolic Kinematic Analysis =====
    print_step_header("Running Symbolic Kinematic Analysis")

    primary_joint_ids = ['d']
    primary_vals = [omega_a] # Symbolic input velocity

    # Run kinematic analysis using the Symbolic Ring (SR)
    Phi_solution, ordered_gm_edges, kin_report = DaviesKinematicAnalysis(
        mechanism, primary_joint_ids, primary_vals,
        lambda_dim=3, generate_report=True, ring=SR # Specify SR
    )

    if Phi_solution is not None:
        print("\nSymbolic Kinematic Analysis Results (Velocities):")
        results = {}
        for i, val in enumerate(Phi_solution):
            edge_label = ordered_gm_edges[i][2]
            # Check if val is symbolic before simplifying
            if hasattr(val, 'simplify_full'):
                 results[edge_label] = val.simplify_full() # Simplify symbolic results
            else:
                 results[edge_label] = val # Keep as is if not symbolic (e.g., numeric result from SR)
            print(f"  ω_{edge_label:<5} = {results[edge_label]}")

        # Optionally substitute numerical values into symbolic results
        try:
            # Ensure substitution happens only on symbolic elements
            subs_results = {}
            # Use a floating point number directly instead of QQ for substitution
            omega_a_val = 0.15  # This avoids the QQ(num,den) syntax
            for k, v in results.items():
                if hasattr(v, 'subs'):
                    subs_results[k] = v.subs({omega_a: omega_a_val})
                else:
                    subs_results[k] = v # Keep numeric results as they are

            print("\nResults with omega_a = 0.15 substituted:")
            for joint_id in ['a', 'b', 'c', 'd']:
                 # Handle potential errors during numerical evaluation
                 try:
                     print(f"  ω_{joint_id}: {N(subs_results[joint_id], digits=8)} rad/s")
                 except TypeError as e_num:
                     print(f"  ω_{joint_id}: Could not evaluate numerically - {subs_results[joint_id]}")
                     
        except Exception as e:
            print(f"\nCould not substitute numerical values: {e}")

        # Save symbolic report
        save_report_to_file(kin_report, os.path.join(REPORTS_DIR, "symbolic_fourbar_kinematic_report.txt"), detailed=True)
        # Also print report to console
        print("\n--- Symbolic Kinematic Report ---")
        # Use detailed=True for symbolic report formatting
        print(format_report(kin_report, detailed=True, precision=8))
        print("--- End Report ---")

    else:
        print("Symbolic kinematic analysis failed.")

    # Static symbolic analysis could be added similarly if needed
    # ...

    return mechanism, Phi_solution, None

def example_symbolic_fourbar_full():
    """Run four-bar linkage example with fully symbolic joint positions and input."""
    print_section_header("Fully Symbolic Four-Bar Linkage Example")

    # Define symbolic variables for joint positions and input velocity
    var('a_x, a_y, b_x, b_y, c_x, c_y, d_x, d_y, omega_d')

    # Create mechanism with symbolic geometry
    print("Creating fully symbolic mechanism...")
    mechanism = Mechanism()
    mechanism.add_body("ground")
    mechanism.add_body("link_ab")
    mechanism.add_body("link_bc")
    mechanism.add_body("link_cd")
    
    # Define all joints with symbolic positions
    mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'point': vector(SR, [a_x, a_y, 0])
    })
    mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'point': vector(SR, [b_x, b_y, 0])
    })
    mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'point': vector(SR, [c_x, c_y, 0])
    })
    mechanism.add_joint("d", "link_cd", "ground", "revolute", 1, {
        'point': vector(SR, [d_x, d_y, 0])
    })

    print(f"Fully symbolic mechanism created with {len(mechanism.bodies)} bodies and {len(mechanism.joints)} joints.")

    # ===== Symbolic Kinematic Analysis =====
    print_step_header("Running Fully Symbolic Kinematic Analysis")

    primary_joint_ids = ['d']
    primary_vals = [omega_d] # Symbolic input velocity
    
    # Run kinematic analysis using the Symbolic Ring (SR)
    Phi_solution, ordered_gm_edges, kin_report = DaviesKinematicAnalysis(
        mechanism, primary_joint_ids, primary_vals,
        lambda_dim=3, generate_report=True, ring=SR
    )

    if Phi_solution is not None:
        print("\nFully Symbolic Kinematic Analysis Results (Velocities):")
        results = {}
        for i, val in enumerate(Phi_solution):
            edge_label = ordered_gm_edges[i][2]
            # Check if val is symbolic before simplifying
            if hasattr(val, 'simplify_full'):
                 results[edge_label] = val.simplify_full() # Simplify symbolic results
            else:
                 results[edge_label] = val # Keep as is if not symbolic
            print(f"  ω_{edge_label:<5} = {results[edge_label]}")

        # Substitute sample values to verify the solution
        try:
            # Define sample values that match the verified four-bar
            subs_dict = {
                a_x: 0, a_y: 0,
                b_x: 0, b_y: 0.35,
                c_x: 0.45, c_y: 0.15,
                d_x: 0.40, d_y: 0,
                omega_d: 0.15
            }
            
            subs_results = {}
            for k, v in results.items():
                if hasattr(v, 'subs'):
                    subs_results[k] = v.subs(subs_dict)
                else:
                    subs_results[k] = v
            
            print("\nResults with sample geometry and omega_d = 0.15 substituted:")
            for joint_id in ['a', 'b', 'c', 'd']:
                try:
                    print(f"  ω_{joint_id}: {N(subs_results[joint_id], digits=8)} rad/s")
                except Exception as e_num:
                    print(f"  ω_{joint_id}: Could not evaluate numerically - {subs_results[joint_id]}")
        except Exception as e:
            print(f"\nCould not substitute sample values: {e}")

        # Save symbolic report
        report_path = os.path.join(REPORTS_DIR, "fully_symbolic_fourbar_report.txt")
        save_report_to_file(kin_report, report_path, detailed=True)
        print(f"\nFully symbolic report saved to {report_path}")

    else:
        print("Fully symbolic kinematic analysis failed.")

    return mechanism, Phi_solution, None

def example_fully_symbolic_static_fourbar():
    """Run four-bar linkage static analysis with fully symbolic geometry and input torque."""
    print_section_header("Fully Symbolic Static Four-Bar Linkage Example")

    # Define symbolic variables for geometry and input torque
    var('a_x, a_y, b_x, b_y, c_x, c_y, d_x, d_y, T_a_in')

    # Create mechanism using symbolic coordinates and SR ring
    print("Creating fully symbolic mechanism for static analysis...")
    mechanism = Mechanism()

    # Add bodies
    mechanism.add_body(0) # Ground
    mechanism.add_body(1)
    mechanism.add_body(2)
    mechanism.add_body(3)

    # Define all joints with symbolic positions
    mechanism.add_joint("a", 0, 1       , "revolute", 1, {
        'point': vector(SR, [a_x, a_y, 0])
    })
    mechanism.add_joint("b", 1, 2, "revolute", 1, {
        'point': vector(SR, [b_x, b_y, 0])
    })
    mechanism.add_joint("c", 2, 3, "revolute", 1, {
        'point': vector(SR, [c_x, c_y, 0])
    })
    mechanism.add_joint("d", 3, 0, "revolute", 1, {
        'point': vector(SR, [d_x, d_y, 0])
    })

    print(f"Fully symbolic static mechanism created with {len(mechanism.bodies)} bodies and {len(mechanism.joints)} joints.")

    # ===== Fully Symbolic Static Analysis =====
    print_step_header("Running Fully Symbolic Static Analysis")

    external_actions = {'a': 1, 'd': 1}  # Active torques at a and d
    primary_constraints = ['a_Tz_active'] # Input torque at 'a' is primary
    primary_constraint_values = [T_a_in]  # Use symbolic variable T_a_in

    # Run static analysis using the Symbolic Ring (SR)
    Psi_solution, ordered_ga_constraints, static_report = DaviesStaticAnalysis(
        mechanism,
        primary_constraints,
        primary_constraint_values,
        external_actions,
        lambda_dim=3,
        generate_report=True,
        ring=SR # Specify SR
    )

    if Psi_solution is not None:
        print("\nFully Symbolic Static Analysis Results (Forces/Torques):")
        symbolic_results = {}
        for i, val in enumerate(Psi_solution):
            constraint_id = ordered_ga_constraints[i]
            # Simplify symbolic results
            if hasattr(val, 'simplify_full'):
                simplified_val = val.simplify_full()
                symbolic_results[constraint_id] = simplified_val
            else:
                symbolic_results[constraint_id] = val # Keep as is if not symbolic
            unit = "N" if ("Rx" in constraint_id or "Ry" in constraint_id or "Rv" in constraint_id) else "Nm"
            print(f"  {constraint_id:<12} = {symbolic_results[constraint_id]} {unit}")

        # Verify by substituting numerical values matching the verified static case
        try:
            subs_dict = {
                a_x: 0, a_y: 0,
                b_x: 2, b_y: 0,
                c_x: 4, c_y: 2,
                d_x: 0, d_y: 2,
                T_a_in: 10 # Use numeric value from verified example
            }
            subs_results = {}
            print(f"\nVerifying with geometry a({subs_dict[a_x]},{subs_dict[a_y]}), b({subs_dict[b_x]},{subs_dict[b_y]}), c({subs_dict[c_x]},{subs_dict[c_y]}), d({subs_dict[d_x]},{subs_dict[d_y]}) and T_a_in = {subs_dict[T_a_in]} Nm:")

            for k, v in symbolic_results.items():
                if hasattr(v, 'subs'):
                    subs_results[k] = v.subs(subs_dict)
                else:
                    subs_results[k] = v # Keep numeric results as they are

            for constraint_id in ordered_ga_constraints:
                unit = "N" if ("Rx" in constraint_id or "Ry" in constraint_id or "Rv" in constraint_id) else "Nm"
                try:
                    # Use N() for numerical evaluation after substitution
                    print(f"  {constraint_id:<12}: {N(subs_results[constraint_id], digits=8)} {unit}")
                except TypeError:
                    print(f"  {constraint_id:<12}: Could not evaluate numerically - {subs_results[constraint_id]} {unit}")

        except Exception as e:
            print(f"\nCould not substitute numerical values for verification: {e}")

        # Save fully symbolic static report
        report_path = os.path.join(REPORTS_DIR, "fully_symbolic_fourbar_static_report.txt")
        save_static_report_to_file(static_report, report_path, detailed=True)
        print(f"\nFully symbolic static report saved to {report_path}")
        # Also print report to console
        print("\n--- Fully Symbolic Static Report ---")
        print(format_static_report(static_report, detailed=True, precision=8))
        print("--- End Report ---")

    else:
        print("Fully symbolic static analysis failed.")

    return mechanism, Psi_solution, None


if __name__ == "__main__":
    print_section_header("Davies' Method Examples", char='=')
    
    # Run examples
    fourbar_mech, fourbar_kin, fourbar_static = user_verified_fourbar_example()
    slider_mech, slider_kin, _ = example_slider_crank()
    symbolic_mech, symbolic_kin, _ = example_symbolic_fourbar()
    
    # Run the fully symbolic four-bar example
    fully_symbolic_mech, fully_symbolic_kin, _ = example_symbolic_fourbar_full()

    # Run the fully symbolic static four-bar example
    fully_symbolic_static_mech, fully_symbolic_static_psi, _ = example_fully_symbolic_static_fourbar()

    print("\nExamples completed. Reports saved to 'reports' directory.")
