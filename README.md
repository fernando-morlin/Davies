# Davies' Method Implementation for Mechanism Analysis

An implementation of Davies' method for analyzing planar mechanisms using SageMath. This library enables both kinematic and static analysis through a systematic graph-theory approach.

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

## Supported Joint Types

### Revolute Joint (R)
- 1 degree of freedom (rotation)
- 2 constraint components (forces)
```python
mechanism.add_joint('a', 'body1', 'body2', 'revolute', 1,
                   {'point': vector(QQ,[x,y,0])})
```

### Prismatic Joint (P)
- 1 degree of freedom (translation)
- 2 constraint components (perpendicular force + moment)
```python
mechanism.add_joint('p', 'body1', 'body2', 'prismatic', 1,
                   {'direction': vector(QQ,[ux,uy,0])})
```

### Planar Joint (E)
- 3 degrees of freedom (rotation + translation)
- No constraints
```python
mechanism.add_joint('e', 'body1', 'body2', 'planar', 3,
                   {'point': vector(QQ,[x,y,0])})
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

## Best Practices

1. Use exact arithmetic with `QQ` when possible:
```python
from sage.all import QQ
omega = QQ(1)    # Exact 1
omega = QQ(1,2)  # Exact 1/2
```

2. Consistent coordinate systems:
- Use right-handed coordinate system
- Z-axis points out of plane
- Define joint positions relative to global origin

3. Mechanism topology:
- Start with ground as body 0
- Number bodies consistently
- Verify joint connectivity

## Further Reading

- Davies, T. H. "Mechanical Networks" (original paper)
- Repository examples: See `example.py` and `fourbar_solver.py`
- Generated reports in `reports/` directory after analysis