"""
Example demonstrating exhaustive type synthesis for a five-bar gripper mechanism.
This example explores joint type combinations while ensuring revolute joints
at finger-object contact points to model friction.
"""
import os
import sys

# Add parent directory to path so we can import the davies_method module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results
from sage.all import vector, QQ
from itertools import product

def filter_out_rigid_joints(results):
    """Filter out mechanism types that contain rigid joints."""
    filtered_results = []
    
    for result in results:
        # Handle both old (2-element) and new (3-element) formats
        if len(result) >= 3:
            joint_types, constraint_set, mobility = result
        else:
            joint_types, constraint_set = result
            mobility = None
        
        # Check if any joints are rigid
        has_rigid = any("rigid" in jt for jt in joint_types.values())
        
        if not has_rigid:
            # Maintain the same format as input
            if mobility is not None:
                filtered_results.append((joint_types, constraint_set, mobility))
            else:
                filtered_results.append((joint_types, constraint_set))
    
    print(f"Filtered out {len(results) - len(filtered_results)} mechanism types with rigid joints.")
    return filtered_results

def fivebar_gripper_synthesis():
    """
    Type synthesis for a five-bar gripper mechanism.
    
    Key requirements:
    - Revolute joints at finger-object contact points to model friction
    - Analysis for both binary object (2 connections) and ternary object (3 connections)
    """
    print("=== Five-Bar Gripper Mechanism Type Synthesis ===")
    
    # Create directory for results
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reports")
    os.makedirs(results_dir, exist_ok=True)
    
    # ========================================================================
    # CASE A: Binary object (degree 2) - object connects through two joints
    # ========================================================================
    print("\n\n===== CASE A: Binary Object Configuration =====")
    print("Implementing type synthesis for gripper with binary object...")
    
    # Create the mechanism topology for binary object case
    binary_mechanism = Mechanism()
    
    # Add bodies (links)
    binary_mechanism.add_body("link_1")      # Base/ground link
    binary_mechanism.add_body("link_2")      # Left finger
    binary_mechanism.add_body("link_3")      # Object (binary - connects to 2 points)
    binary_mechanism.add_body("link_4")      # Right finger
    binary_mechanism.add_body("link_5")      # Connection link
    
    # Add joints with geometry for binary object case
    # In this configuration, the object connects through joints b and c
    binary_mechanism.add_joint("a", "link_1", "link_2", "virtual", 0, {
        'point': vector(QQ, [0, 0, 0])
    })
    binary_mechanism.add_joint("b", "link_2", "link_3", "virtual", 0, {  # Object contact point (MUST be revolute)
        'point': vector(QQ, [2, 0, 0])
    })
    binary_mechanism.add_joint("c", "link_3", "link_4", "virtual", 0, {  # Object contact point (MUST be revolute)
        'point': vector(QQ, [3, 1, 0])
    })
    binary_mechanism.add_joint("d", "link_4", "link_1", "virtual", 0, {
        'point': vector(QQ, [4, 0, 0])
    })
    binary_mechanism.add_joint("e", "link_2", "link_5", "virtual", 0, {
        'point': vector(QQ, [1, -1, 0])
    })
    binary_mechanism.add_joint("f", "link_5", "link_1", "virtual", 0, {
        'point': vector(QQ, [0, -2, 0])
    })
    
    print(f"Created binary object gripper with {len(binary_mechanism.bodies)} links and {len(binary_mechanism.joints)} joints")
    print("Joints b and c are designated as object connection points (must be revolute)")
    
    # Define design requirements for binary object
    binary_requirements = {
        'required_freedoms': {
            # Object contact points MUST be revolute (for friction)
            'b': ['Rz'],
            'c': ['Rz'],
            # Other joints are unconstrained for now
        },
        'required_constraints': {
            # Constrain translations ONLY at contact points to ensure revolute joints with friction
            'b': ['Tx', 'Ty'],
            'c': ['Tx', 'Ty'],
            # Remove translational constraints from other joints to allow for more design possibilities
        }
    }
    
    # Run the synthesis for binary object case
    print("\nRunning type synthesis for binary object case...")
    binary_results = DaviesTypeSynthesis(
        binary_mechanism,
        design_requirements=binary_requirements,
        lambda_dim=3,     # Planar mechanism
        ring=QQ,          # Use exact arithmetic
        max_bases=5000,   # Explore many possibilities
        is_gripper=True   # This is a gripper mechanism
    )
    
    # Filter out solutions with rigid joints
    print("\nFiltering binary results to remove solutions with rigid joints...")
    binary_results = filter_out_rigid_joints(binary_results)
    
    # Save binary results
    binary_results_file = os.path.join(results_dir, "fivebar_binary_object_synthesis_results.txt")
    save_type_synthesis_results(
        binary_results, 
        binary_results_file,
        mechanism=binary_mechanism,
        design_requirements=binary_requirements,
        lambda_dim=3
    )
    
    # ========================================================================
    # CASE B: Ternary object (degree 3) - object connects to three joints
    # ========================================================================
    print("\n\n===== CASE B: Ternary Object Configuration =====")
    print("Implementing type synthesis for gripper with ternary object...")
    
    # Create the mechanism topology for ternary object case
    ternary_mechanism = Mechanism()
    
    # Add bodies (links)
    ternary_mechanism.add_body("link_1")      # Base/ground link
    ternary_mechanism.add_body("link_2")      # First connected link
    ternary_mechanism.add_body("link_3")      # Object (ternary - connects to 3 joints: a, b, e)
    ternary_mechanism.add_body("link_4")      # Second connected link
    ternary_mechanism.add_body("link_5")      # Third connected link
    
    # Add joints with geometry for ternary object case
    ternary_mechanism.add_joint("a", "link_1", "link_3", "virtual", 0, {  # Object connection (MUST be revolute)
        'point': vector(QQ, [0, 0, 0])
    })
    ternary_mechanism.add_joint("b", "link_3", "link_2", "virtual", 0, {  # Object connection (MUST be revolute)
        'point': vector(QQ, [2, 0, 0])
    })
    ternary_mechanism.add_joint("c", "link_2", "link_4", "virtual", 0, {
        'point': vector(QQ, [3, 1, 0])
    })
    ternary_mechanism.add_joint("d", "link_4", "link_1", "virtual", 0, {
        'point': vector(QQ, [4, 0, 0])
    })
    ternary_mechanism.add_joint("e", "link_3", "link_5", "virtual", 0, {  # Object connection (MUST be revolute)
        'point': vector(QQ, [1, 2, 0])
    })
    ternary_mechanism.add_joint("f", "link_5", "link_1", "virtual", 0, {
        'point': vector(QQ, [0, 2, 0])
    })
    
    print(f"Created ternary object gripper with {len(ternary_mechanism.bodies)} links and {len(ternary_mechanism.joints)} joints")
    print("Joints a, b and e are designated as object connection points (must be revolute)")
    
    # Define design requirements for ternary object
    ternary_requirements = {
        'required_freedoms': {
            # Object connections MUST be revolute (for friction)
            'a': ['Rz'],
            'b': ['Rz'],
            'e': ['Rz'],
            # Other joints are unconstrained for now
        },
        'required_constraints': {
            # Constrain translations ONLY at object connection points for friction modeling
            'a': ['Tx', 'Ty'],
            'b': ['Tx', 'Ty'],
            'e': ['Tx', 'Ty'],
            # Remove constraints from non-contact joints to explore more design possibilities
        }
    }
    
    # Run the synthesis for ternary object case
    print("\nRunning type synthesis for ternary object case...")
    ternary_results = DaviesTypeSynthesis(
        ternary_mechanism,
        design_requirements=ternary_requirements,
        lambda_dim=3,     # Planar mechanism
        ring=QQ,          # Use exact arithmetic
        max_bases=5000,   # Explore many possibilities
        is_gripper=True   # This is a gripper mechanism
    )
    
    # Filter out solutions with rigid joints
    print("\nFiltering ternary results to remove solutions with rigid joints...")
    ternary_results = filter_out_rigid_joints(ternary_results)
    
    # Save ternary results
    ternary_results_file = os.path.join(results_dir, "fivebar_ternary_object_synthesis_results.txt")
    save_type_synthesis_results(
        ternary_results, 
        ternary_results_file,
        mechanism=ternary_mechanism,
        design_requirements=ternary_requirements,
        lambda_dim=3
    )
    
    # Process and display results for both cases
    print("\n\n===== SYNTHESIS RESULTS SUMMARY =====")
    
    # Process binary object results
    print("\nBinary Object Results:")
    if binary_results and len(binary_results) > 0:
        print(f"Found {len(binary_results)} valid mechanism types")
        
        # Categorize the binary object results
        binary_categories = categorize_results(binary_results)
        print(f"Unique joint combinations: {len(binary_categories)}")
        
        # Display top configurations
        print("\nTop 5 binary object configurations:")
        for i, (sig, indices) in enumerate(sorted(binary_categories.items(), 
                                               key=lambda x: len(x[1]), 
                                               reverse=True)[:5]):
            print(f" {i+1}. {sig}: {len(indices)} variations")
    else:
        print("No valid binary object mechanism types found")
    
    # Process ternary object results
    print("\nTernary Object Results:")
    if ternary_results and len(ternary_results) > 0:
        print(f"Found {len(ternary_results)} valid mechanism types")
        
        # Categorize the ternary object results
        ternary_categories = categorize_results(ternary_results)
        print(f"Unique joint combinations: {len(ternary_categories)}")
        
        # Display top configurations
        print("\nTop 5 ternary object configurations:")
        for i, (sig, indices) in enumerate(sorted(ternary_categories.items(), 
                                               key=lambda x: len(x[1]), 
                                               reverse=True)[:5]):
            print(f" {i+1}. {sig}: {len(indices)} variations")
    else:
        print("No valid ternary object mechanism types found")
    
    print("\nType synthesis complete.")
    print(f"Binary object results saved to: {binary_results_file}")
    print(f"Ternary object results saved to: {ternary_results_file}")
    
    return binary_results, ternary_results

def categorize_results(results):
    """Helper function to categorize joint types in results"""
    joint_type_categories = {}
    
    for i, result in enumerate(results):
        # Handle both old (2-element) and new (3-element) formats
        if len(result) >= 3:
            joint_types, _, _ = result
        else:
            joint_types, _ = result
        
        # Create a signature for this configuration
        signature = []
        for joint_id in sorted(joint_types.keys()):
            joint_type = joint_types.get(joint_id, 'unknown')
            # Simplify the joint type name for categorization
            if 'revolute' in joint_type:
                simple_type = 'R'  # Revolute
            elif 'prismatic_x' in joint_type:
                simple_type = 'Px' # Prismatic in x
            elif 'prismatic_y' in joint_type:
                simple_type = 'Py' # Prismatic in y
            elif 'pin_in_slot_x' in joint_type:
                simple_type = 'PSx' # Pin-in-slot in x
            elif 'pin_in_slot_y' in joint_type:
                simple_type = 'PSy' # Pin-in-slot in y
            else:
                simple_type = '?'
            signature.append(f"{joint_id}:{simple_type}")
        
        sig_key = ', '.join(signature)
        if sig_key not in joint_type_categories:
            joint_type_categories[sig_key] = []
        joint_type_categories[sig_key].append(i)
    
    return joint_type_categories

if __name__ == "__main__":
    try:
        binary_results, ternary_results = fivebar_gripper_synthesis()
        print("\nAnalysis complete. Review the results files for detailed mechanism specifications.")
    except Exception as e:
        import traceback
        print(f"An error occurred: {e}")
        traceback.print_exc()
