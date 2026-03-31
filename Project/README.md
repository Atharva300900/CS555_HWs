# 2D Viscous Burgers' Equation Solver

This code solves the 2D viscous Burgers' equations on the domain [0,1] x [0,1] using finite differences. It supports multiple time integration schemes and allows comparison with an exact solution.

## Usage

Run the solver using:

```python
u, v, X, Y, info = solve(N, nu, T, dt, method='rk4', adv_type='central')
```

### Parameters

| Parameter  | Description                  |
|------------|------------------------------|
| `N`        | Grid resolution              |
| `nu`       | Viscosity                    |
| `T`        | Final time                   |
| `dt`       | Time step                    |
| `method`   | `'euler'`, `'rk4'`, `'be'`, `'cn'` |
| `adv_type` | `'central'` or `'upwind'`    |

## Reproducing Results

Reference case:

```python
N  = 80
nu = 0.01
T  = 0.5
dt = 1e-4

u, v, X, Y, info = solve(N, nu, T, dt, method='rk4')

# Compare with exact solution
ue, ve = exact_solution(X, Y, T, nu)
```

## Experiments

### Convergence Study

```python
for N in [20, 40, 80]:
    solve(N, nu, T, dt=1e-4)
```

### Stability Analysis

```python
probe_explicit_stability(N, nu, T, method='euler')
probe_explicit_stability(N, nu, T, method='rk4')
```

### Method Comparison

```python
solve(N, nu, T, dt, method='euler')
solve(N, nu, T, dt, method='rk4')
solve(N, nu, T, dt, method='be')
solve(N, nu, T, dt, method='cn')
```

## Output

The solver returns:

- `u`, `v` -- numerical solution
- `X`, `Y` -- grid coordinates
- `info` -- dictionary containing:
  - L2 error
  - Linf error
  - Number of time steps
  - Runtime
  - CFL information (for RK4)

## Visualization

```python
import matplotlib.pyplot as plt

plt.contourf(X, Y, u)
plt.colorbar()
plt.title("u field")
plt.show()
```

## Notes

- Explicit methods require CFL stability conditions
- RK4 includes a CFL check
- Implicit methods are more stable for larger time steps
- Upwind scheme is more stable but less accurate than central
