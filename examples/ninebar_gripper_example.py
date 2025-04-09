"""
Example demonstrating type synthesis for a gripper mechanism using a Baranov chain.
This example focuses on a nine-bar mechanism where link 3 serves as a ternary object
connected via three revolute joints (c, d, j) to model friction at contact points.
"""
import os
import sys

# Add parent directory to path so we can import the davies_method module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results
from sage.all import vector, QQ

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

def baranov_gripper_type_synthesis():
    """
    Type synthesis for a nine-bar Baranov chain gripper mechanism.
    
    Key requirements:
    - Link 3 is a ternary object connecting to three joints (c, d, j)
    - Revolute joints at object-contact points to model friction
    """
    print("=== Nine-Bar Baranov Chain Gripper Mechanism Type Synthesis ===")
    
    # Create a seed mechanism with the Baranov chain topology from the example
    mechanism = Mechanism()
    
    # Add links/bodies - using named links for better readability
    mechanism.add_body("link_1")  # Base/ground link
    mechanism.add_body("link_2")  # Connected to multiple joints
    mechanism.add_body("link_3")  # This is the ternary object being grasped
    mechanism.add_body("link_4")  # Part of the finger structure
    mechanism.add_body("link_5")  # Part of the finger structure
    mechanism.add_body("link_6")  # Part of the finger structure
    mechanism.add_body("link_7")  # Part of the finger structure
    mechanism.add_body("link_8")  # Part of the finger structure
    mechanism.add_body("link_9")  # Part of the finger structure
    
    # Create incidence matrix for the Baranov chain as specified
    incidence_matrix = [
        [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],  # link 1
        [0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0],  # link 2
        [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],  # link 3 (object)
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1],  # link 4
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],  # link 5
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],  # link 6
        [0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # link 7
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # link 8
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],  # link 9
    ]
    
    # Joint names
    joint_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
    
    # Add joints based on the incidence matrix
    for joint_idx, joint_id in enumerate(joint_names):
        # Find the two links connected by this joint
        connected_links = []
        for link_idx, row in enumerate(incidence_matrix):
            if row[joint_idx] == 1:
                connected_links.append(f"link_{link_idx+1}")
        
        if len(connected_links) == 2:
            # Place joints at arbitrary but distinct positions
            position = [joint_idx*2, joint_idx*1.5, 0]
            mechanism.add_joint(joint_id, connected_links[0], connected_links[1], "virtual", 0, {
                'point': vector(QQ, position)
            })
    
    print(f"Created Baranov chain mechanism with {len(mechanism.bodies)} links and {len(mechanism.joints)} joints.")
    print("Link 3 is designated as the ternary object being grasped via joints c, d, and j.")
    
    # Create directory for results
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reports")
    os.makedirs(results_dir, exist_ok=True)
    
    # ========================================================================
    # CASE: Ternary object (link 3) connected at three points (c, d, j)
    # ========================================================================
    print("\n===== Ternary Object Configuration =====")
    print("Implementing type synthesis for gripper with ternary object (link 3)...")
    print("- Joints c, d, j: contact points with object (MUST be revolute for friction)")
    print("- Other joints: unconstrained for exploration of design possibilities")
    
    # Define design requirements for ternary object
    design_requirements = {
        'required_freedoms': {
            # Object connections MUST be revolute (for friction)
            'c': ['Rz'],
            'd': ['Rz'],
            'j': ['Rz'],
            # Other joints are unconstrained for now
        },
        'required_constraints': {
            # Constrain translations ONLY at object connection points for friction modeling
            'c': ['Tx', 'Ty'],
            'd': ['Tx', 'Ty'],
            'j': ['Tx', 'Ty'],
            # Remove constraints from non-contact joints to explore more design possibilities
        }
    }
    
    # Run the synthesis for ternary object case
    print("\nRunning type synthesis for ternary object case...")
    results = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements,
        lambda_dim=3,     # Planar mechanism
        ring=QQ,          # Use exact arithmetic
        max_bases=5000,   # Explore many possibilities
        is_gripper=True   # This is a gripper mechanism
    )
    
    # Filter out solutions with rigid joints
    print("\nFiltering results to remove solutions with rigid joints...")
    filtered_results = filter_out_rigid_joints(results)
    
    # Save results
    results_file = os.path.join(results_dir, "ninebar_ternary_object_synthesis_results.txt")
    save_type_synthesis_results(
        filtered_results, 
        results_file,
        mechanism=mechanism,
        design_requirements=design_requirements,
        lambda_dim=3,
        is_gripper=True
    )
    
    # Process and display results
    print("\n===== SYNTHESIS RESULTS SUMMARY =====")
    
    if filtered_results and len(filtered_results) > 0:
        print(f"Found {len(filtered_results)} valid mechanism types")
        
        # Categorize the results
        categories = categorize_results(filtered_results)
        print(f"Unique joint combinations: {len(categories)}")
        
        # Display top configurations
        print("\nTop 5 configurations:")
        for i, (sig, indices) in enumerate(sorted(categories.items(), 
                                               key=lambda x: len(x[1]), 
                                               reverse=True)[:5]):
            print(f" {i+1}. {sig}: {len(indices)} variations")
        
        # Display first solution's joint types
        print("\nExample joint type configuration (first solution):")
        joint_types = filtered_results[0][0]
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
        
        # Verify that object contact joints are revolute
        contact_joints_revolute = all("revolute" in joint_types.get(j, '') for j in ['c', 'd', 'j'])
        if contact_joints_revolute:
            print("\n✓ Verified: All object contact joints (c, d, j) are revolute as required.")
        else:
            print("\n✗ Warning: Not all object contact joints are revolute!")
            for j in ['c', 'd', 'j']:
                print(f"  Joint {j}: {joint_types.get(j, 'not found')}")
    else:
        print("No valid mechanism types found")
    
    print("\nType synthesis complete.")
    print(f"Results saved to: {results_file}")
    
    return filtered_results

if __name__ == "__main__":
    try:
        results = baranov_gripper_type_synthesis()
        print("\nAnalysis complete. Review the results file for detailed mechanism specifications.")
    except Exception as e:
        import traceback
        print(f"An error occurred: {e}")
        traceback.print_exc()
