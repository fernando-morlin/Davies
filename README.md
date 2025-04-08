# Davies' Method Implementation for Mechanism Analysis

An implementation of Davies' method for analyzing planar and spatial mechanisms using SageMath. This library enables both kinematic and static analysis through a systematic graph-theory approach, with support for both numeric and symbolic calculations, and specialized features for gripper analysis.

## Quick Start

1. Install SageMath 9.0+ from [sagemath.org](https://www.sagemath.org)

2. Analyze a four-bar mechanism in just a few lines:

```python
from davies_method import Mechanism, DaviesKinematicAnalysis
from sage.all import vector, QQ

# Create the mechanism
m = Mechanism()
m.add_joint('a', 'ground', 'crank',  'revolute', 1, {'point': vector(QQ,[0,0,0])})
m.add_joint('b', 'crank',  'coupler', 'revolute', 1, {'point': vector(QQ,[2,0,0])})
m.add_joint('c', 'coupler', 'rocker', 'revolute', 1, {'point': vector(QQ,[4,2,0])})
m.add_joint('d', 'rocker', 'ground',  'revolute', 1, {'point': vector(QQ,[0,2,0])})

# Run kinematic analysis with input at joint 'a' rotating at 1 rad/s
Phi, edges, report = DaviesKinematicAnalysis(m, ['a'], [QQ(1)], generate_report=True)

# Print angular velocities
print("\nJoint Angular Velocities:")
for i, edge in enumerate(edges):
    print(f"{edge[2]}: {float(Phi[i])} rad/s")
```

## Core Concepts

### 1. Mechanism Definition

A mechanism is defined by:
- Bodies (links)
- Joints connecting the bodies
- For each joint:
  - Type (revolute, prismatic, etc.)
  - Degrees of Freedom (f)
  - Geometry (position, axis)

### 2. Kinematic Analysis

Calculates velocities throughout the mechanism:
```python
# Example: Slider-crank with input at joint 'a'
mechanism = Mechanism()
mechanism.add_joint('a', 'ground', 'crank', 'revolute', 1, 
                   {'point': vector(QQ,[0,0,0])})
mechanism.add_joint('b', 'crank', 'rod', 'revolute', 1,
                   {'point': vector(QQ,[2,0,0])})
mechanism.add_joint('c', 'rod', 'slider', 'revolute', 1,
                   {'point': vector(QQ,[5,0,0])})
mechanism.add_joint('d', 'slider', 'ground', 'prismatic', 1,
                   {'direction': vector(QQ,[1,0,0])})

# Analyze with crank rotating at 0.15 rad/s
Phi, edges, report = DaviesKinematicAnalysis(
    mechanism, ['a'], [QQ(15)/100], generate_report=True
)
```

### 3. Static Analysis 

Calculates forces and moments for static equilibrium:
```python
# Example: Finding reaction forces with input torque
mechanism = Mechanism()
# ... add joints as before ...

# Specify active torque at joint 'a' and reaction at 'd'
external_actions = {'a': 1, 'd': 1}  

# Apply 10 Nm input torque at joint 'a'
Psi, constraints, report = DaviesStaticAnalysis(
    mechanism,
    ['a_Tz_active'],   # Known constraint
    [QQ(10)],          # 10 Nm torque
    external_actions,
    generate_report=True
)
```

#### External Forces and Constraints

The static analysis supports defining external forces and constraints:

```python
# Define external forces at specific joints
external_actions = {
    'a': 1,  # 1 external action at joint 'a' (e.g., input torque)
    'c': 2   # 2 external actions at joint 'c' (e.g., contact forces)
}

# Specify known constraint magnitudes
primary_constraint_ids = ['a_Tz_active', 'c_Rx']  # Torque at 'a', force at 'c'
primary_values = [QQ(10), QQ(5)]  # 10 Nm torque, 5 N force

Psi, constraints, report = DaviesStaticAnalysis(
    mechanism,
    primary_constraint_ids,
    primary_values,
    external_actions,
    generate_report=True
)
```

### 4. Symbolic Analysis

Performs analysis using symbolic variables to produce exact parametric solutions:

```python
from sage.all import vector, SR, var
from davies_method import Mechanism, DaviesKinematicAnalysis

# Define symbolic variables
var('a_x, a_y, b_x, b_y, c_x, c_y, d_x, d_y, omega_d')

# Create mechanism with symbolic coordinates
mechanism = Mechanism()
mechanism.add_body("ground")
mechanism.add_body("link_ab")
mechanism.add_body("link_bc")
mechanism.add_body("link_cd")

# Define joints with symbolic positions
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

# Run symbolic analysis (specify SR as the ring for calculations)
Phi, edges, report = DaviesKinematicAnalysis(
    mechanism, ['d'], [omega_d], lambda_dim=3, 
    generate_report=True, ring=SR
)

# Results will be symbolic expressions in terms of defined variables
for i, val in enumerate(Phi):
    joint_name = edges[i][2]
    print(f"ω_{joint_name} = {val.simplify_full()}")
```

## Supported Joint Types

### Planar Mechanism Joints (lambda_dim=3)

#### Revolute Joint (R)
- 1 degree of freedom (rotation)
- 2 constraint components (forces)
```python
mechanism.add_joint('a', 'body1', 'body2', 'revolute', 1,
                   {'point': vector(QQ,[x,y,0])})
```

#### Prismatic Joint (P)
- 1 degree of freedom (translation)
- 2 constraint components (perpendicular force + moment)
```python
mechanism.add_joint('p', 'body1', 'body2', 'prismatic', 1,
                   {'direction': vector(QQ,[ux,uy,0])})
```

#### Pin-in-Slot Joint
- 2 degrees of freedom (rotation + translation along slot)
- 1 constraint component
- Available as 'pin_in_slot_x' and 'pin_in_slot_y'
```python
mechanism.add_joint('p', 'body1', 'body2', 'pin_in_slot_x', 2,
                   {'point': vector(QQ,[x,y,0])})
```

#### Planar Joint (E)
- 3 degrees of freedom (rotation + translation)
- No constraints
```python
mechanism.add_joint('e', 'body1', 'body2', 'planar', 3,
                   {'point': vector(QQ,[x,y,0])})
```

### Spatial Mechanism Joints (lambda_dim=6)

#### Revolute Joint
- 1 degree of freedom (rotation about specified axis)
- Available as 'revolute_x', 'revolute_y', 'revolute_z'

#### Prismatic Joint
- 1 degree of freedom (translation along specified axis)
- Available as 'prismatic_x', 'prismatic_y', 'prismatic_z'

#### Spherical Joint
- 3 degrees of freedom (rotation about X, Y, and Z axes)

#### Rigid Connection
- 0 degrees of freedom (fully constrained)

## Gripper Analysis

This library includes specialized features for gripper analysis, including:

### Gripper Type Synthesis

```python
from davies_method import Mechanism, DaviesTypeSynthesis, save_type_synthesis_results

# Create a gripper seed mechanism
mechanism = Mechanism()
mechanism.add_body('ground')
mechanism.add_body('left_finger')
mechanism.add_body('right_finger')

# Add virtual joints (will be determined by type synthesis)
mechanism.add_joint('a', 'ground', 'left_finger', 'planar', 3, {'point': vector(QQ,[0,0,0])})
mechanism.add_joint('b', 'ground', 'right_finger', 'planar', 3, {'point': vector(QQ,[0,0,0])})
mechanism.add_joint('c', 'left_finger', 'right_finger', 'planar', 3, {'point': vector(QQ,[1,1,0])})

# Design requirements for gripper mechanism
design_requirements = {
    'required_freedoms': {
        'a': ['Rz'],      # Left finger should rotate
        'b': ['Rz']       # Right finger should rotate
    },
    'required_constraints': {
        'c': ['Tx', 'Ty'] # Contact point constraints
    }
}

# Run type synthesis specifically for grippers
results = DaviesTypeSynthesis(
    mechanism,
    design_requirements,
    lambda_dim=3
)

# Save results with gripper-specific recommendations
save_type_synthesis_results(
    results, 
    'gripper_types.txt', 
    mechanism, 
    design_requirements,
    lambda_dim=3,
    is_gripper=True  # Enable gripper-specific recommendations
)
```

### Gripper Static Analysis

For gripper mechanisms, static analysis can be used to calculate:
- Gripping forces at contact points
- Required actuation forces/torques
- Reaction forces at joints

```python
# Define external actions for gripper analysis
external_actions = {
    'c': 1,  # Contact point external force
    'a': 1   # Actuation force/torque
}

# Apply a known actuation force and solve for the resulting grip force
Psi, constraints, report = DaviesStaticAnalysis(
    mechanism,
    ['a_Tz_active'],  # Known actuation torque
    [QQ(10)],         # 10 Nm actuation
    external_actions,
    generate_report=True
)

# Find the grip force at contact point 'c'
for i, constraint_name in enumerate(constraints):
    if constraint_name.startswith('c_'):
        print(f"Grip force at {constraint_name}: {float(Psi[i])} N")
```

## Reporting and Debugging 

1. Generate detailed analysis report:
```python
Phi, edges, report = DaviesKinematicAnalysis(..., generate_report=True)
save_report_to_file(report, 'analysis_report.txt')
```

2. Report contents include:
- Mechanism topology
- Matrix calculations
- System equations
- Final results

3. Common issues:
- Singular matrices: Check mechanism constraints
- Rank mismatch: Verify number of primary inputs
- NaN results: Look for numerical precision issues

4. Static analysis reports:
```python
Psi, constraints, report = DaviesStaticAnalysis(..., generate_report=True)
save_static_report_to_file(report, 'static_analysis.txt')
```

## Type Synthesis

The library includes tools for automated type synthesis using Davies' method and matroid theory:

```python
from davies_method import Mechanism, DaviesTypeSynthesis

# Create a seed mechanism (topology only)
mechanism = Mechanism()
# ...add links and joints...

# Define design requirements (optional)
design_requirements = {
    'required_freedoms': {
        'a': ['Rz'],  # Joint 'a' should be revolute
        'b': ['Tx']   # Joint 'b' should allow x translation
    },
    'required_constraints': {
        'c': ['Tx', 'Ty']  # Joint 'c' should constrain x,y translation
    }
}

# Run type synthesis
results = DaviesTypeSynthesis(
    mechanism,
    design_requirements,
    lambda_dim=3  # Planar mechanism
)

# Each result contains a dictionary of joint types
for i, (joint_types, _) in enumerate(results[:5]):
    print(f"\nMechanism Type #{i+1}:")
    for joint_id, joint_type in sorted(joint_types.items()):
        print(f"  Joint {joint_id}: {joint_type}")
```

The type synthesis algorithm:
1. Builds a seed mechanism with virtual joints
2. Computes the Network Unit Action Matrix
3. Uses matroid theory to systematically identify all valid constraint patterns
4. Determines appropriate joint types for each pattern

For grippers, it provides additional recommendations:
- Potential contact points (typically revolute joints)
- Potential actuation joints (often prismatic joints)
- Fixed connections

See examples in the `examples/` directory for more usage details.

## Best Practices

1. Use exact arithmetic with `QQ` when possible:
```python
from sage.all import QQ
omega = QQ(1)    # Exact 1
omega = QQ(1,2)  # Exact 1/2
```

2. Use symbolic calculation with `SR` for parametric analysis:
```python
from sage.all import SR, var
var('omega_a')  # Define symbolic variable
# Use ring=SR in analysis function
Phi, edges, report = DaviesKinematicAnalysis(..., ring=SR)
```

3. Consistent coordinate systems:
- Use right-handed coordinate system
- Z-axis points out of plane
- Define joint positions relative to global origin

4. Mechanism topology:
- Start with ground as body 0
- Number bodies consistently
- Verify joint connectivity

5. Gripper-specific considerations:
- Define potential contact points between fingers
- Include actuation mechanisms (typically prismatic joints)
- Consider using pin-in-slot joints for self-alignment

## Further Reading

- Davies, T. H. "Mechanical Networks" (original paper)
- Repository examples: See `example.py` and `fourbar_solver.py`
- Generated reports in `reports/` directory after analysis
- Symbolic examples: See `example_symbolic_fourbar()` and `example_symbolic_fourbar_full()` in `example.py`