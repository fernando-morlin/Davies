"""
Example demonstrating type synthesis for Watt's mechanism with various joint configurations.

Watt's mechanism is a six-bar linkage known for generating approximate straight-line motion.
This example explores various configurations including hydraulic actuation and self-aligning joints.
"""
import os
import sys

# Add parent directory to path so we can import the davies_method module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results
from sage.all import vector, QQ

def watt_type_synthesis_example():
    """
    Type synthesis analysis for various configurations of a Watt six-bar linkage.
    
    Watt's mechanism: A six-bar linkage that can generate approximate straight-line motion
    at a point on the coupler link. The mechanism consists of two separate four-bar linkages
    that share a common link, allowing for interesting motion characteristics.
    """
    print("=== Watt Six-Bar Mechanism Type Synthesis Example ===")
    print("\nImplementing type synthesis for a Watt six-bar mechanism...")
    
    # Create a seed mechanism with Watt's six-bar linkage topology
    mechanism = Mechanism()
    
    # Add bodies (links)
    mechanism.add_body("ground")     # Fixed base/ground
    mechanism.add_body("crank_1")    # First input crank
    mechanism.add_body("coupler")    # Coupler link with approximate straight-line motion
    mechanism.add_body("rocker")     # Rocker link
    mechanism.add_body("crank_2")    # Second crank
    mechanism.add_body("output")     # Output link
    
    # Add joints with geometry
    # First four-bar chain:
    mechanism.add_joint("a", "ground", "crank_1", "revolute", 1, {
        'point': vector(QQ, [0, 0, 0])
    })
    mechanism.add_joint("b", "crank_1", "coupler", "revolute", 1, {
        'point': vector(QQ, [2, 0, 0])
    })
    mechanism.add_joint("c", "coupler", "rocker", "revolute", 1, {
        'point': vector(QQ, [4, 1, 0])
    })
    mechanism.add_joint("d", "rocker", "ground", "revolute", 1, {
        'point': vector(QQ, [6, 0, 0])
    })
    
    # Second four-bar chain (shares the coupler link with the first chain):
    mechanism.add_joint("e", "coupler", "crank_2", "revolute", 1, {
        'point': vector(QQ, [3, 3, 0])
    })
    mechanism.add_joint("f", "crank_2", "output", "revolute", 1, {
        'point': vector(QQ, [5, 4, 0])
    })
    mechanism.add_joint("g", "output", "ground", "revolute", 1, {
        'point': vector(QQ, [7, 2, 0])
    })
    
    print(f"Created Watt six-bar mechanism with {len(mechanism.bodies)} links and {len(mechanism.joints)} joints.")
    print("This mechanism can generate approximate straight-line motion at certain points on the coupler.")
    
    # Create directory for results
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reports")
    os.makedirs(results_dir, exist_ok=True)
    
    # =====================================================================
    # CASE 1: Classic Watt mechanism with all revolute joints
    # =====================================================================
    print("\nCASE 1: Classic Watt mechanism with all revolute joints")
    print("- All joints (a-g) are revolute")
    print("- This is the traditional Watt six-bar linkage")
    print("- Expected mobility: 1 DOF")
    
    # Define specific requirements for the classic design (all revolute)
    design_requirements_classic = {
        'required_freedoms': {
            'a': ['Rz'], 'b': ['Rz'], 'c': ['Rz'], 'd': ['Rz'], 
            'e': ['Rz'], 'f': ['Rz'], 'g': ['Rz']
        },
        'required_constraints': {
            'a': ['Tx', 'Ty'], 'b': ['Tx', 'Ty'], 'c': ['Tx', 'Ty'], 'd': ['Tx', 'Ty'],
            'e': ['Tx', 'Ty'], 'f': ['Tx', 'Ty'], 'g': ['Tx', 'Ty']
        }
    }
    
    # Run the type synthesis for Case 1
    print("\nRunning analysis for classic Watt mechanism with all revolute joints...")
    results_classic = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_classic,
        lambda_dim=3,  # Planar mechanism
        ring=QQ,       # Use exact arithmetic
        max_bases=10   # Limit the number of results
    )
    
    # Save results to file
    results_file_classic = os.path.join(results_dir, "watt_classic_results.txt")
    save_type_synthesis_results(
        results_classic, 
        results_file_classic,
        mechanism=mechanism,
        design_requirements=design_requirements_classic,
        lambda_dim=3
    )
    
    # Display results
    if results_classic and len(results_classic) > 0:
        joint_types = results_classic[0][0]  # First solution's joint types
        print("\nVerifying Case 1 results:")
        
        # Check if all joints are revolute
        all_revolute = all("revolute" in joint_types.get(j, '') for j in ['a', 'b', 'c', 'd', 'e', 'f', 'g'])
        
        if all_revolute:
            print("✓ MATCH: All joints are revolute as specified.")
        else:
            print("✗ MISMATCH: Not all joints are revolute.")
            for j in ['a', 'b', 'c', 'd', 'e', 'f', 'g']:
                print(f"  Joint {j}: {joint_types.get(j, 'not found')}")
        
        # Show detailed results
        print("\nJoint types for the classic Watt mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 1")
    
    # =====================================================================
    # CASE 2: Watt mechanism with hydraulic actuator
    # =====================================================================
    print("\nCASE 2: Watt mechanism with hydraulic actuator")
    print("- Joint a is a hydraulic actuator (prismatic)")
    print("- The remaining joints are revolute")
    print("- This design is used in heavy machinery where hydraulic power is preferred")
    print("- Expected mobility: 1 DOF")
    
    # Define requirements for the hydraulic actuator design
    design_requirements_hydraulic = {
        'required_freedoms': {
            # Hydraulic actuator (prismatic joint)
            'a': ['Tx'],  # Linear actuation along x-axis
            
            # All other joints remain revolute
            'b': ['Rz'], 'c': ['Rz'], 'd': ['Rz'], 
            'e': ['Rz'], 'f': ['Rz'], 'g': ['Rz']
        },
        'required_constraints': {
            # Hydraulic cylinder constraints
            'a': ['Ty', 'Rz'],  # Constrain y-translation and rotation
            
            # Revolute joint constraints
            'b': ['Tx', 'Ty'], 'c': ['Tx', 'Ty'], 'd': ['Tx', 'Ty'],
            'e': ['Tx', 'Ty'], 'f': ['Tx', 'Ty'], 'g': ['Tx', 'Ty']
        }
    }
    
    # Run the type synthesis for Case 2
    print("\nRunning analysis for Watt mechanism with hydraulic actuator...")
    results_hydraulic = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_hydraulic,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 2 results
    results_file_hydraulic = os.path.join(results_dir, "watt_hydraulic_results.txt")
    save_type_synthesis_results(
        results_hydraulic, 
        results_file_hydraulic,
        mechanism=mechanism,
        design_requirements=design_requirements_hydraulic,
        lambda_dim=3
    )
    
    # Display Case 2 results
    if results_hydraulic and len(results_hydraulic) > 0:
        joint_types = results_hydraulic[0][0]  # First solution's joint types
        print("\nVerifying Case 2 results:")
        
        # Check hydraulic actuator joint
        hydraulic_a = "prismatic_x" in joint_types.get('a', '')
        
        if hydraulic_a:
            print("✓ MATCH: Joint a is a prismatic joint (hydraulic actuator) as specified.")
        else:
            print(f"✗ MISMATCH: Joint a is {joint_types.get('a', 'not found')} instead of prismatic_x.")
        
        # Show detailed results
        print("\nJoint types for the hydraulically actuated Watt mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 2")
    
    # =====================================================================
    # CASE 3: Self-adjusting Watt mechanism with pin-in-slot joints
    # =====================================================================
    print("\nCASE 3: Self-adjusting Watt mechanism with pin-in-slot joints")
    print("- Joints b and e are pin-in-slot joints for self-adjustment")
    print("- Joint a is a prismatic actuator")
    print("- Remaining joints are revolute")
    print("- This design allows for self-alignment and adaptability")
    print("- Expected mobility: 3 DOF")
    
    # Define requirements for self-adjusting mechanism
    design_requirements_self_adjusting = {
        'required_freedoms': {
            # Prismatic actuator
            'a': ['Tx'],
            
            # Pin-in-slot joints for self-adjustment
            'b': ['Rz', 'Ty'],  # Rotation + vertical adjustment
            'e': ['Rz', 'Tx'],  # Rotation + horizontal adjustment
            
            # Standard revolute joints
            'c': ['Rz'], 'd': ['Rz'], 'f': ['Rz'], 'g': ['Rz']
        },
        'required_constraints': {
            # Actuator constraints
            'a': ['Ty', 'Rz'],
            
            # Pin-in-slot constraints (each constrains one translation direction)
            'b': ['Tx'],  # Constrains only x-translation
            'e': ['Ty'],  # Constrains only y-translation
            
            # Revolute joint constraints
            'c': ['Tx', 'Ty'], 'd': ['Tx', 'Ty'], 
            'f': ['Tx', 'Ty'], 'g': ['Tx', 'Ty']
        }
    }
    
    # Run the type synthesis for Case 3
    print("\nRunning analysis for self-adjusting Watt mechanism...")
    results_self_adjusting = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_self_adjusting,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 3 results
    results_file_self_adjusting = os.path.join(results_dir, "watt_self_adjusting_results.txt")
    save_type_synthesis_results(
        results_self_adjusting, 
        results_file_self_adjusting,
        mechanism=mechanism,
        design_requirements=design_requirements_self_adjusting,
        lambda_dim=3
    )
    
    # Display Case 3 results
    if results_self_adjusting and len(results_self_adjusting) > 0:
        joint_types = results_self_adjusting[0][0]  # First solution's joint types
        print("\nVerifying Case 3 results:")
        
        # Check pin-in-slot joints
        pin_in_slot_b = "pin_in_slot" in joint_types.get('b', '')
        pin_in_slot_e = "pin_in_slot" in joint_types.get('e', '')
        
        if pin_in_slot_b and pin_in_slot_e:
            print("✓ MATCH: Joints b and e are pin-in-slot joints as specified.")
        else:
            print("✗ MISMATCH with pin-in-slot requirements:")
            print(f"  Joint b: {joint_types.get('b', 'not found')}")
            print(f"  Joint e: {joint_types.get('e', 'not found')}")
        
        # Show detailed results
        print("\nJoint types for the self-adjusting Watt mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
    else:
        print("\nNo solutions found for Case 3")
    
    # =====================================================================
    # CASE 4: Specialized Watt mechanism for precise path generation
    # =====================================================================
    print("\nCASE 4: Specialized Watt mechanism for precise path generation")
    print("- Ground joints (a, d, g) are precise ball-joint revolutes")
    print("- Joints b, f are prismatic for adjustable link lengths")
    print("- Joints c, e are pin-in-slot for path correction")
    print("- This specialized design allows for precise path generation with adjustment capability")
    
    # Define requirements for specialized path generation mechanism
    design_requirements_specialized = {
        'required_freedoms': {
            # Ground joints (revolute)
            'a': ['Rz'], 'd': ['Rz'], 'g': ['Rz'],
            
            # Adjustable length joints (prismatic)
            'b': ['Tx'], 'f': ['Ty'],
            
            # Path correction joints (pin-in-slot)
            'c': ['Rz', 'Tx'], 'e': ['Rz', 'Ty']
        },
        'required_constraints': {
            # Ground revolute joints
            'a': ['Tx', 'Ty'], 'd': ['Tx', 'Ty'], 'g': ['Tx', 'Ty'],
            
            # Prismatic joints
            'b': ['Ty', 'Rz'], 'f': ['Tx', 'Rz'],
            
            # Pin-in-slot joints
            'c': ['Ty'], 'e': ['Tx']
        }
    }
    
    # Run the type synthesis for Case 4
    print("\nRunning analysis for specialized path generation Watt mechanism...")
    results_specialized = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements_specialized,
        lambda_dim=3,
        ring=QQ,
        max_bases=10
    )
    
    # Save Case 4 results
    results_file_specialized = os.path.join(results_dir, "watt_specialized_results.txt")
    save_type_synthesis_results(
        results_specialized, 
        results_file_specialized,
        mechanism=mechanism,
        design_requirements=design_requirements_specialized,
        lambda_dim=3
    )
    
    # Display Case 4 results
    if results_specialized and len(results_specialized) > 0:
        joint_types = results_specialized[0][0]  # First solution's joint types
        
        print("\nJoint types for the specialized path generation Watt mechanism:")
        for joint_id, joint_type in sorted(joint_types.items()):
            print(f"  Joint {joint_id}: {joint_type}")
            
        # Calculate theoretical mobility
        revolute_count = sum(1 for j, t in joint_types.items() if "revolute" in t)
        prismatic_count = sum(1 for j, t in joint_types.items() if "prismatic" in t and "pin_in_slot" not in t)
        pin_in_slot_count = sum(1 for j, t in joint_types.items() if "pin_in_slot" in t)
        
        total_dof = revolute_count + prismatic_count + 2 * pin_in_slot_count
        mobility = total_dof - 3 * (len(mechanism.bodies) - 1)
        
        print(f"\nTheoretical mobility analysis:")
        print(f"  - Revolute joints: {revolute_count} × 1 DOF = {revolute_count}")
        print(f"  - Prismatic joints: {prismatic_count} × 1 DOF = {prismatic_count}")
        print(f"  - Pin-in-slot joints: {pin_in_slot_count} × 2 DOF = {pin_in_slot_count * 2}")
        print(f"  - Total DOF: {total_dof}")
        print(f"  - Mobility: M = {total_dof} - 3({len(mechanism.bodies) - 1}) = {mobility}")
    else:
        print("\nNo solutions found for Case 4")
    
    # =====================================================================
    # Print summary of all cases
    # =====================================================================
    print("\n=== Type Synthesis Summary for Watt's Mechanism ===")
    print("1. Classic Watt mechanism: All revolute joints")
    print("2. Hydraulic Watt mechanism: Prismatic actuator at joint a")
    print("3. Self-adjusting Watt mechanism: Pin-in-slot joints at b and e")
    print("4. Specialized path generation: Mixed joint types for precision")
    
    print("\nType synthesis complete. Results saved to:")
    print(f"- Case 1 (Classic): {results_file_classic}")
    print(f"- Case 2 (Hydraulic): {results_file_hydraulic}")
    print(f"- Case 3 (Self-adjusting): {results_file_self_adjusting}")
    print(f"- Case 4 (Specialized): {results_file_specialized}")
    
    return {
        "case1": results_classic,
        "case2": results_hydraulic,
        "case3": results_self_adjusting,
        "case4": results_specialized
    }

if __name__ == "__main__":
    try:
        results = watt_type_synthesis_example()
        print("\nAll analyses complete. Watt mechanism type synthesis successful.")
        print("\nApplication notes:")
        print("- The hydraulic configuration is ideal for heavy machinery applications")
        print("- Self-adjusting configuration helps with manufacturing tolerances")
        print("- Specialized path generation is useful for precision mechanisms")
        print("- All solutions maintain the fundamental characteristics of Watt's mechanism")
    except Exception as e:
        import traceback
        print(f"An error occurred: {e}")
        traceback.print_exc()
