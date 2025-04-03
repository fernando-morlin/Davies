"""Example demonstrating both kinematic and static analysis using Davies' method."""
import os
from sage.all import vector, QQ, RDF, pi, cos, sin, N
from davies_method import (
    Mechanism, DaviesKinematicAnalysis, DaviesStaticAnalysis,
    save_report_to_file, save_static_report_to_file
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

if __name__ == "__main__":
    print_section_header("Davies' Method Examples", char='=')
    
    # Run examples
    fourbar_mech, fourbar_kin, fourbar_static = user_verified_fourbar_example()
    slider_mech, slider_kin, _ = example_slider_crank()
    
    print("\nExamples completed. Reports saved to 'reports' directory.")
