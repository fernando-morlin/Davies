"""
Example demonstrating type synthesis for a five-bar mechanism.
This example follows the same structure as the nine-bar Baranov chain example.
"""
import os
import sys

# Add parent directory to path so we can import the davies_method module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results
from sage.all import vector, QQ

def fivebar_type_synthesis_example():
    """
    Example of automated type synthesis for a five-bar linkage.
    """
    print("=== Five-Bar Mechanism Type Synthesis Example ===")
    print("\nImplementing type synthesis for a five-bar mechanism...")
    
    # Create a seed mechanism with the desired topology
    # The seed just defines the connectivity, not the joint types
    mechanism = Mechanism()
    
    # Add bodies (links)
    mechanism.add_body("ground")    # Fixed base/ground
    mechanism.add_body("link_ab")   # Link connecting joints a and b
    mechanism.add_body("link_bc")   # Link connecting joints b and c
    mechanism.add_body("link_cd")   # Link connecting joints c and d
    mechanism.add_body("link_de")   # Link connecting joints d and e
    
    # Add joints with geometry (points are arbitrary at this stage)
    mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'point': vector(QQ, [0, 0, 0])
    })
    mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'point': vector(QQ, [2, 0, 0])
    })
    mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'point': vector(QQ, [4, 1, 0])
    })
    mechanism.add_joint("d", "link_cd", "link_de", "revolute", 1, {
        'point': vector(QQ, [3, 3, 0])
    })
    mechanism.add_joint("e", "link_de", "ground", "revolute", 1, {
        'point': vector(QQ, [0, 2, 0])
    })
    
    print(f"Created five-bar mechanism with {len(mechanism.bodies)} links and {len(mechanism.joints)} joints.")
    
    # Create directory for results
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reports")
    os.makedirs(results_dir, exist_ok=True)
    
    # =====================================================================
    # CASE 1: Basic five-bar with all revolute joints
    # =====================================================================
    print("\nCASE 1: Basic five-bar with all revolute joints")
    print("- All joints (a, b, c, d, e) are revolute")
    print("- Expected mechanism: Typical five-bar linkage with 1 DOF")
    
    # Define specific requirements for all revolute joints
    design_requirements_revolute = {
        'required_freedoms': {
            'a': ['Rz'], 'b': ['Rz'], 'c': ['Rz'], 'd': ['Rz'], 'e': ['Rz']
        },
        'required_constraints': {
            'a': ['Tx', 'Ty'], 'b': ['Tx', 'Ty'], 'c': ['Tx', 'Ty'], 
            'd': ['Tx', 'Ty'], 'e': ['Tx', 'Ty']
        }
    }
    
    # Run the type synthesis for Case 1
    print("\nRunning analysis for basic five-bar with revolute joints...")
    results_revolute = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_revolute,
        lambda_dim=3,  # Planar mechanism
        ring=QQ,       # Use exact arithmetic
        max_bases=10   # Limit the number of results
    )
    
    # Save results to file
    results_file_revolute = os.path.join(results_dir, "fivebar_revolute_synthesis_results.txt")
    save_type_synthesis_results(
        results_revolute, 
        results_file_revolute,
        mechanism=mechanism,
        design_requirements=design_requirements_revolute,
        lambda_dim=3
    )
    
    # Display and analyze the results
    if results_revolute and len(results_revolute) > 0:
        joint_types = results_revolute[0][0]  # First solution's joint types
        print("\nVerifying Case 1 results:")
        
        # Check if all joints are revolute
        all_revolute = all("revolute" in joint_types.get(j, '') for j in ['a', 'b', 'c', 'd', 'e'])
        
        if all_revolute:
            print("✓ MATCH: All joints are revolute as specified.")
        else:
            print("✗ MISMATCH: Not all joints are revolute.")
            for j in ['a', 'b', 'c', 'd', 'e']:
                print(f"  Joint {j}: {joint_types.get(j, 'not found')}")
        
        # Show detailed results
        print("\nJoint types for the all-revolute five-bar mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 1")
    
    # =====================================================================
    # CASE 2: Five-bar with pin-in-slot joints for increased mobility
    # =====================================================================
    print("\nCASE 2: Five-bar with pin-in-slot joints")
    print("- Joints a, e are fixed revolute joints to ground")
    print("- Joints b, d are pin-in-slot joints")
    print("- Joint c is revolute")
    print("- Expected mechanism: Five-bar with 3 DOF")
    
    # Define requirements for pin-in-slot design
    design_requirements_pin_in_slot = {
        'required_freedoms': {
            # Fixed revolute ground joints
            'a': ['Rz'], 'e': ['Rz'],
            
            # Pin-in-slot joints - allow rotation and translation in one direction
            'b': ['Rz', 'Ty'],  # Pin-in-slot allowing y translation
            'd': ['Rz', 'Tx'],  # Pin-in-slot allowing x translation
            
            # Regular revolute joint
            'c': ['Rz']
        },
        'required_constraints': {
            # For revolute joints, constrain all translations
            'a': ['Tx', 'Ty'], 'e': ['Tx', 'Ty'], 'c': ['Tx', 'Ty'],
            
            # For pin-in-slot joints, constrain only one translation direction
            'b': ['Tx'],  # Constrain x translation
            'd': ['Ty']   # Constrain y translation
        }
    }
    
    # Run the type synthesis for Case 2
    print("\nRunning analysis for five-bar with pin-in-slot joints...")
    results_pin_in_slot = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_pin_in_slot,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 2 results
    results_file_pin_in_slot = os.path.join(results_dir, "fivebar_pin_in_slot_synthesis_results.txt")
    save_type_synthesis_results(
        results_pin_in_slot, 
        results_file_pin_in_slot,
        mechanism=mechanism,
        design_requirements=design_requirements_pin_in_slot,
        lambda_dim=3
    )
    
    # Display Case 2 results
    if results_pin_in_slot and len(results_pin_in_slot) > 0:
        joint_types = results_pin_in_slot[0][0]  # First solution's joint types
        print("\nVerifying Case 2 results:")
        
        # Check pin-in-slot joints
        pin_slot_b = "pin_in_slot" in joint_types.get('b', '')
        pin_slot_d = "pin_in_slot" in joint_types.get('d', '')
        
        if pin_slot_b and pin_slot_d:
            print("✓ MATCH: Joints b and d are pin-in-slot joints as specified.")
        else:
            print("✗ MISMATCH: Expected pin-in-slot joints at b and d.")
            print(f"  Joint b: {joint_types.get('b', 'not found')}")
            print(f"  Joint d: {joint_types.get('d', 'not found')}")
        
        # Show detailed results
        print("\nJoint types for the pin-in-slot five-bar mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 2")
    
    # =====================================================================
    # CASE 3: Five-bar with prismatic joints for linear actuation
    # =====================================================================
    print("\nCASE 3: Five-bar with prismatic joints for linear actuation")
    print("- Joints a, e are fixed revolute joints to ground")
    print("- Joints c, d are prismatic for linear actuation")
    print("- Joint b is revolute")
    print("- Expected mechanism: Modified five-bar for linear motion")
    
    # Define requirements for prismatic joint design
    design_requirements_prismatic = {
        'required_freedoms': {
            # Fixed revolute ground joints
            'a': ['Rz'], 'e': ['Rz'],
            
            # Regular revolute joint
            'b': ['Rz'],
            
            # Prismatic joints
            'c': ['Tx'],  # Prismatic in x-direction
            'd': ['Ty']   # Prismatic in y-direction
        },
        'required_constraints': {
            # For revolute joints, constrain translations
            'a': ['Tx', 'Ty'], 'b': ['Tx', 'Ty'], 'e': ['Tx', 'Ty'],
            
            # For prismatic joints, constrain rotation and one translation
            'c': ['Rz', 'Ty'],  # Constrain rotation and y-translation
            'd': ['Rz', 'Tx']   # Constrain rotation and x-translation
        }
    }
    
    # Run the type synthesis for Case 3
    print("\nRunning analysis for five-bar with prismatic joints...")
    results_prismatic = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_prismatic,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 3 results
    results_file_prismatic = os.path.join(results_dir, "fivebar_prismatic_synthesis_results.txt")
    save_type_synthesis_results(
        results_prismatic, 
        results_file_prismatic,
        mechanism=mechanism,
        design_requirements=design_requirements_prismatic,
        lambda_dim=3
    )
    
    # Display Case 3 results
    if results_prismatic and len(results_prismatic) > 0:
        joint_types = results_prismatic[0][0]  # First solution's joint types
        print("\nVerifying Case 3 results:")
        
        # Check prismatic joints
        prismatic_c = "prismatic_x" in joint_types.get('c', '')
        prismatic_d = "prismatic_y" in joint_types.get('d', '')
        
        if prismatic_c and prismatic_d:
            print("✓ MATCH: Joint c is prismatic_x and d is prismatic_y as specified.")
        else:
            print("✗ MISMATCH: Expected specific prismatic joints.")
            print(f"  Joint c: {joint_types.get('c', 'not found')}")
            print(f"  Joint d: {joint_types.get('d', 'not found')}")
        
        # Show detailed results
        print("\nJoint types for the prismatic five-bar mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 3")
    
    # =====================================================================
    # CASE 4: Five-bar with mixed joint types for specialized motion
    # =====================================================================
    print("\nCASE 4: Five-bar with mixed joint types for specialized motion")
    print("- Joint a is prismatic (actuator)")
    print("- Joint e is a pin-in-slot joint")
    print("- Joints b, c, d are revolute")
    print("- Expected mechanism: Specialized five-bar for complex motion")
    
    # Define requirements for mixed joint type design
    design_requirements_mixed = {
        'required_freedoms': {
            # Prismatic actuator
            'a': ['Tx'],
            
            # Pin-in-slot joint
            'e': ['Rz', 'Tx'],
            
            # Regular revolute joints
            'b': ['Rz'], 'c': ['Rz'], 'd': ['Rz']
        },
        'required_constraints': {
            # For prismatic actuator, constrain rotation and y-translation
            'a': ['Rz', 'Ty'],
            
            # For pin-in-slot joint, constrain only y-translation
            'e': ['Ty'],
            
            # For revolute joints, constrain translations
            'b': ['Tx', 'Ty'], 'c': ['Tx', 'Ty'], 'd': ['Tx', 'Ty']
        }
    }
    
    # Run the type synthesis for Case 4
    print("\nRunning analysis for five-bar with mixed joint types...")
    results_mixed = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_mixed,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 4 results
    results_file_mixed = os.path.join(results_dir, "fivebar_mixed_synthesis_results.txt")
    save_type_synthesis_results(
        results_mixed, 
        results_file_mixed,
        mechanism=mechanism,
        design_requirements=design_requirements_mixed,
        lambda_dim=3
    )
    
    # Display Case 4 results
    if results_mixed and len(results_mixed) > 0:
        joint_types = results_mixed[0][0]  # First solution's joint types
        print("\nJoint types for the mixed five-bar mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 4")
    
    print("\nType synthesis complete. Results saved to:")
    print(f"- Case 1 (all revolute): {results_file_revolute}")
    print(f"- Case 2 (pin-in-slot): {results_file_pin_in_slot}")
    print(f"- Case 3 (prismatic): {results_file_prismatic}")
    print(f"- Case 4 (mixed): {results_file_mixed}")
    
    return {
        "case1": results_revolute,
        "case2": results_pin_in_slot,
        "case3": results_prismatic,
        "case4": results_mixed
    }

if __name__ == "__main__":
    try:
        results = fivebar_type_synthesis_example()
        print("\nAll analyses complete.")
    except Exception as e:
        import traceback
        print(f"An error occurred: {e}")
        traceback.print_exc()
