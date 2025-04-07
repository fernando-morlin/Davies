"""
Example demonstrating type synthesis for a four-bar mechanism.
"""
import os
import sys

# Add parent directory to path so we can import the davies_method module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results
from sage.all import vector, QQ

def fourbar_type_synthesis_example():
    """
    Example of automated type synthesis for a four-bar linkage.
    """
    print("=== Four-Bar Mechanism Type Synthesis Example ===")
    
    # Create a seed mechanism with the desired topology
    # The seed just defines the connectivity, not the joint types
    mechanism = Mechanism()
    
    # Add bodies (links)
    mechanism.add_body("ground")
    mechanism.add_body("link_ab")
    mechanism.add_body("link_bc")
    mechanism.add_body("link_cd")
    
    # Add joints with geometry (points are arbitrary at this stage)
    mechanism.add_joint("a", "ground", "link_ab", "revolute", 1, {
        'point': vector(QQ, [0, 0, 0])
    })
    mechanism.add_joint("b", "link_ab", "link_bc", "revolute", 1, {
        'point': vector(QQ, [2, 0, 0])
    })
    mechanism.add_joint("c", "link_bc", "link_cd", "revolute", 1, {
        'point': vector(QQ, [4, 2, 0])
    })
    mechanism.add_joint("d", "link_cd", "ground", "revolute", 1, {
        'point': vector(QQ, [0, 2, 0])
    })
    
    # Define design requirements (optional)
    # For example, require joint 'a' to be a revolute joint
    design_requirements = {
        'required_freedoms': {
            'a': ['Rz']  # 'a' should allow rotation about z-axis
        },
        'required_constraints': {
            'a': ['Tx', 'Ty']  # 'a' should constrain x and y translations
        }
    }
    
    # Run the type synthesis
    results = DaviesTypeSynthesis(
        mechanism,
        design_requirements=design_requirements,
        lambda_dim=3,  # Planar mechanism
        ring=QQ,       # Use exact arithmetic
        max_bases=50   # Limit the number of results for this example
    )
    
    # Save results to file
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reports")
    os.makedirs(results_dir, exist_ok=True)
    results_file = os.path.join(results_dir, "fourbar_type_synthesis_results.txt")
    save_type_synthesis_results(results, results_file)
    
    return results

if __name__ == "__main__":
    results = fourbar_type_synthesis_example()
    print("\nType synthesis complete. Check the results file for details.")
