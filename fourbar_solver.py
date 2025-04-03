#!/usr/bin/env python3
"""
Four-Bar Linkage Solver

This script solves the kinematic analysis of a four-bar linkage mechanism
using Davies' method. The user can specify the coordinates of the four joints,
which joint is the input, and the angular velocity of the input joint.
"""
import sys
from sage.all import vector, QQ, N, pi, cos, sin
from davies_method import Mechanism, DaviesKinematicAnalysis, format_report

def solve_fourbar(ax, ay, bx, by, cx, cy, dx, dy, input_joint='d', input_velocity=1.0, 
                  generate_report=True, report_path=None):
    """
    Solve the kinematic analysis of a four-bar linkage mechanism.
    
    Parameters:
    -----------
    ax, ay : float or exact fraction
        Coordinates of joint a
    bx, by : float or exact fraction
        Coordinates of joint b
    cx, cy : float or exact fraction
        Coordinates of joint c
    dx, dy : float or exact fraction
        Coordinates of joint d
    input_joint : str, optional
        Specify which joint is the input ('a', 'b', 'c', or 'd'), default is 'd'
    input_velocity : float or exact fraction, optional
        Angular velocity at the input joint, default is 1.0
    generate_report : bool, optional
        Whether to generate a detailed report, default is True
    report_path : str, optional
        Path to save the report if not None, otherwise print to terminal
        
    Returns:
    --------
    dict
        Dictionary containing the angular velocities of all joints
    """
    print(f"Solving four-bar linkage with joints:")
    print(f"  a: ({ax}, {ay})")
    print(f"  b: ({bx}, {by})")
    print(f"  c: ({cx}, {cy})")
    print(f"  d: ({dx}, {dy})")
    print(f"\nInput Configuration:")
    print(f"  Input joint: {input_joint}")
    print(f"  Input angular velocity: {input_velocity} rad/s")
    
    # Create the mechanism
    mechanism = Mechanism()
    
    # Add bodies (links)
    mechanism.add_body("ground")     # Fixed link
    mechanism.add_body("link_ab")    # Link connecting a and b 
    mechanism.add_body("link_bc")    # Link connecting b and c
    mechanism.add_body("link_cd")    # Link connecting c and d
    
    # Add joints with provided geometry
    mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'axis': vector([0, 0, 1]),  # z-axis rotation
        'point': vector([ax, ay, 0])
    })
    
    mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([bx, by, 0])
    })
    
    mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([cx, cy, 0])
    })
    
    mechanism.add_joint("d", "link_cd", "ground", "revolute", 1, {
        'axis': vector([0, 0, 1]),
        'point': vector([dx, dy, 0])
    })
    
    # Convert input angular velocity to QQ for exact arithmetic
    try:
        input_velocity_exact = QQ(input_velocity)
    except (ValueError, TypeError):
        # Fall back to floating point if conversion fails
        input_velocity_exact = float(input_velocity)
        print("Warning: Using floating point arithmetic instead of exact fractions.")
    
    # Validate input joint
    if input_joint not in ['a', 'b', 'c', 'd']:
        print(f"Error: Invalid input joint '{input_joint}'. Must be one of 'a', 'b', 'c', or 'd'.")
        return None
    
    # Perform the kinematic analysis
    primary_joint_ids = [input_joint]
    primary_vals = [input_velocity_exact]
    
    print(f"\nRunning Davies' method with {input_joint} as primary joint...")
    
    Phi_solution, ordered_gm_edges, kin_report = DaviesKinematicAnalysis(
        mechanism, primary_joint_ids, primary_vals, 
        lambda_dim=3, generate_report=generate_report
    )
    
    # Extract results and map joint IDs to angular velocities
    if Phi_solution is None:
        print("Error: Kinematic analysis failed!")
        return None
    
    results = {}
    for i, val in enumerate(Phi_solution):
        joint_id = ordered_gm_edges[i][2]
        results[joint_id] = val
    
    # Display results
    print("\nResults (Angular Velocities):")
    for joint_id in ['a', 'b', 'c', 'd']:
        print(f"  ω_{joint_id}: {N(results[joint_id], digits=8)} rad/s")
    
    # If report was generated, either save to file or print to terminal
    if generate_report and kin_report:
        if report_path:
            # Save to file if path specified
            with open(report_path, 'w') as f:
                f.write(format_report(kin_report))
            print(f"Report saved to {report_path}")
        else:
            # Print detailed report to terminal
            print("\n" + "=" * 50)
            print("DETAILED ANALYSIS REPORT")
            print("=" * 50)
            print(format_report(kin_report))
    
    return results

def parse_fraction(value):
    """Parse a string into a fraction if possible, otherwise return float."""
    try:
        if '/' in value:
            num, denom = value.split('/')
            return QQ(int(num), int(denom))
        else:
            return QQ(value)
    except (ValueError, TypeError):
        return float(value)

def main():
    """Command line interface for the four-bar solver."""
    if len(sys.argv) < 9:
        print("Usage: python fourbar_solver.py ax ay bx by cx cy dx dy [input_joint] [input_velocity] [report_path]")
        print("\nParameters:")
        print("  ax, ay: coordinates of joint a")
        print("  bx, by: coordinates of joint b")
        print("  cx, cy: coordinates of joint c")
        print("  dx, dy: coordinates of joint d")
        print("  input_joint: which joint is the input ('a', 'b', 'c', or 'd', default: 'd')")
        print("  input_velocity: angular velocity of input joint (default: 1.0)")
        print("  report_path: path to save report to file (optional, omit to print to terminal)")
        
        print("\nExamples:")
        print("  python fourbar_solver.py 0 0 0 0.35 0.45 0.15 0.4 0")
        print("  python fourbar_solver.py 0 0 0 0.35 0.45 0.15 0.4 0 a 2.5")
        print("  python fourbar_solver.py 0 0 0 0.35 0.45 0.15 0.4 0 c 1.5 report.txt")
        sys.exit(1)
    
    # Parse coordinates
    coords = [parse_fraction(arg) for arg in sys.argv[1:9]]
    ax, ay, bx, by, cx, cy, dx, dy = coords
    
    # Parse optional arguments with improved defaults
    input_joint = 'd'  # Default is joint d
    input_velocity = QQ(1)  # Default is 1 rad/s
    report_path = None
    
    # Process command line arguments based on how many were provided
    remaining_args = sys.argv[9:]
    if len(remaining_args) >= 1:
        # First optional arg is the input joint
        if remaining_args[0] in ['a', 'b', 'c', 'd']:
            input_joint = remaining_args[0]
        else:
            # If not a valid joint, assume it's the velocity
            try:
                input_velocity = parse_fraction(remaining_args[0])
                # Keep default input_joint = 'd'
            except (ValueError, TypeError):
                print(f"Warning: Invalid input joint or velocity '{remaining_args[0]}', using defaults")
    
    if len(remaining_args) >= 2:
        # Second optional arg is input velocity (only if joint was specified)
        if input_joint in ['a', 'b', 'c', 'd'] and input_joint == remaining_args[0]:
            try:
                input_velocity = parse_fraction(remaining_args[1])
            except (ValueError, TypeError):
                print(f"Warning: Invalid input velocity '{remaining_args[1]}', using default 1.0 rad/s")
    
    if len(remaining_args) >= 3:
        # Third optional arg is report path
        report_path = remaining_args[2]
    
    # Output the configuration being used
    print(f"Configuration:")
    print(f"  Joint coordinates: a({ax},{ay}), b({bx},{by}), c({cx},{cy}), d({dx},{dy})")
    print(f"  Input joint: {input_joint}")
    print(f"  Input velocity: {input_velocity} rad/s")
    print(f"  Report path: {report_path if report_path else 'console output'}")
    
    # Run the analysis
    solve_fourbar(ax, ay, bx, by, cx, cy, dx, dy,
                  input_joint=input_joint, 
                  input_velocity=input_velocity,
                  generate_report=True, 
                  report_path=report_path)

if __name__ == "__main__":
    main()
