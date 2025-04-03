"""
Common utilities shared between kinematic and static analyses for Davies' method.
"""
from sage.all import vector, QQ, N, SR

class Mechanism:
    """Class representing a mechanical system with bodies and joints."""
    def __init__(self):
        self.bodies = [] # List of body IDs (e.g., 0, 1, 2...)
        self.joints = [] # List of tuples: (joint_id, body1, body2)
        self.joint_types = {} # Dict: joint_id -> 'revolute', 'prismatic', etc.
        self.joint_dof = {} # Dict: joint_id -> degrees of freedom (f)
        self.geometry = {} # Dict: joint_id -> {'point': vector(...), 'axis': vector(...), ...}

    def add_body(self, body_id):
        """Add a body to the mechanism."""
        if body_id not in self.bodies:
            self.bodies.append(body_id)
        return self

    def add_joint(self, joint_id, body1, body2, joint_type, joint_dof=1, geometry=None):
        """Add a joint between two bodies."""
        # Ensure bodies exist
        self.add_body(body1)
        self.add_body(body2)
        # Store joint info with explicit DOF
        if not any(j[0] == joint_id for j in self.joints):
            self.joints.append((joint_id, body1, body2))
            self.joint_types[joint_id] = joint_type
            self.joint_dof[joint_id] = joint_dof  # Store DOF explicitly
            if geometry:
                self.geometry[joint_id] = geometry
            else:
                self.geometry[joint_id] = {}
        return self

def format_vector(v, precision=4):
    """Format a vector with clean representation of near-zero values.
    Handle both numeric and symbolic vectors."""
    try:
        components = []
        for x in v:
            # Check if x is symbolic using the parent ring
            if hasattr(x, 'parent') and x.parent() is SR:
                components.append(str(x))
            else:
                # Attempt numerical conversion and formatting
                try:
                    x_num = N(x, prec=53) # Use standard precision for check
                    if abs(x_num) < 1e-10:
                        components.append("0")
                    else:
                        # Format with requested precision
                        components.append(f"{N(x, digits=precision+2):.{precision}f}")
                except TypeError: # Handle cases where N() fails (like QQ)
                    if abs(float(x)) < 1e-10:
                         components.append("0")
                    else:
                         # Try formatting QQ directly if possible
                         try:
                              components.append(f"{x:.{precision}f}")
                         except TypeError:
                              components.append(str(x)) # Fallback to string
        return f"({', '.join(components)})"
    except TypeError: # Fallback for non-iterable v
        return str(v)

def get_dof_for_joint_type(joint_type, lambda_dim=3):
    """Return the number of degrees of freedom (f) for each joint type."""
    # For lambda=3 (planar)
    dof_map_planar = {'revolute': 1, 'prismatic': 1, 'planar': 3, 'fixed': 0}
    return dof_map_planar.get(joint_type, 0)  # Default to 0 DoF if unknown

def get_constraints_for_joint_type(joint_type, lambda_dim=3):
    """Return the number of constraint components (c) for each joint type."""
    # For lambda=3 (planar): f + c = 3
    f = get_dof_for_joint_type(joint_type, lambda_dim)
    c = lambda_dim - f
    return c  # Returns passive constraints 'cp'
