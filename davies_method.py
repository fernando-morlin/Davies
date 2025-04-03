"""
Davies' method implementation for mechanism kinematics and statics.
This is the main module that re-exports all functionality.
"""

# Re-export Mechanism class from common_utils
from common_utils import Mechanism

# Re-export kinematic analysis functionality
from kinematic_analysis import (
    davies_kinematic_analysis as DaviesKinematicAnalysis,
    save_report_to_file,
    format_kinematic_report as format_report
)

# Re-export static analysis functionality
from static_analysis import (
    davies_static_analysis as DaviesStaticAnalysis,
    save_static_report_to_file,
    format_static_report
)

# Define a version
__version__ = '1.0.0'

# If run directly, show basic usage information
if __name__ == '__main__':
    print("Davies' Method Implementation for Mechanism Analysis")
    print("====================================================")
    print("\nImport the following classes and functions to use the package:")
    print("  from davies_method import Mechanism")
    print("  from davies_method import DaviesKinematicAnalysis")
    print("  from davies_method import DaviesStaticAnalysis")
    print("\nSee example.py for detailed usage examples.")