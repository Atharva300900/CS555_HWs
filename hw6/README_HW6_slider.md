# CS555 HW6 Slider Bearing MATLAB Files

This folder contains separate MATLAB scripts for the different parts of the slider-bearing assignment. The scripts use the helper functions provided with the starter code, including the mesh generation, basis functions, quadrature routines, height profile, boundary detection, and restriction matrix setup.

## Required helper files

Keep the following helper files in the same MATLAB working directory:

- `abqfem.m`
- `basis_deriv_12.m`
- `basis_tri_12.m`
- `boundedges.m`
- `box_elem.m`
- `hdr.m`
- `profile_taper_flat.m`
- `restriction.m`
- `trigausspoints.m`
- `triquad.m`
- `zwgll.m`
- `zwuni.m`

These are used by the main scripts to build the triangular mesh, assemble FEM matrices, evaluate the taper-flat height profile, apply boundary restrictions, and perform numerical quadrature.

## Main scripts

### `part1a_load_center_pressure.m`

This script solves the incompressible Reynolds equation for the original taper-flat slider configuration.

Settings used:

- `Ta = pi/180`, corresponding to the 1 degree leading taper
- Dirichlet pressure condition on the entire boundary
- `Ex = 120`, `Ey = 40`

The script computes:

- Load `F`
- Center-of-pressure `xp`
- Nondimensional center-of-pressure `xp/L`
- Average pressure `F/(L W)` in Pa and atm
- A mesh plot of the scaled pressure distribution

The load and center-of-pressure are computed using the assembled FEM mass matrix:

```matlab
F  = ones(1,nb)*Bb*Pb;
xp = (x')*Bb*Pb/F;
```

### `part1b_wedge_convergence.m`

This script verifies the code using the analytical solution for the 1D wedge-bearing problem.

Settings used:

- `Ta = 0`, so the leading 1 degree taper is removed
- `Ey = 40`
- `Ex` is varied over several mesh sizes
- Dirichlet pressure condition is applied only at `x = 0` and `x = L`
- Top and bottom boundaries are left as natural Neumann boundaries

The boundary condition change is implemented using:

```matlab
tol = 1e-12;
boundary_nodes = find(abs(x) < tol | abs(x-L) < tol);
R = restriction(nb,boundary_nodes);
```

The script computes:

- Numerical FEM pressure
- Analytical wedge pressure
- Maximum pointwise pressure error
- Numerical load `F_num`
- Analytical load `F_exact`
- Relative load error
- `xp/L`

It also creates the log-log convergence plot of maximum pointwise pressure error versus `Ex`, which is used to demonstrate second-order convergence.

### `part2a_compressible_taper_flat.m`

This script solves the compressible Reynolds equation for the original taper-flat slider configuration.

Settings used:

- `Ta = pi/180`
- Dirichlet pressure condition on the full boundary
- `Ex = 120`, `Ey = 40`

The compressible pressure is written as:

```matlab
pa = p + patm;
```

The nonlinear dependence of the diffusion coefficient on pressure is handled using fixed-point iteration. At each iteration:

1. Update `pa = p + patm`
2. Rebuild the pressure-dependent diffusion matrix
3. Solve the linearized compressible system
4. Repeat until the relative pressure change is below the tolerance

The script reports:

- Trailing height `h2`
- Pitch angle `gamma`
- Number of fixed-point iterations
- Final relative pressure change
- Load `F`
- Center-of-pressure `xp`
- `xp/L`
- Average relative pressure `F/(L W)`
- Maximum relative pressure
- Maximum absolute pressure

The load is computed using the relative pressure `p`, not the absolute pressure `pa`.

### `part2bc_compare_cases.m`

This script runs the three requested configurations and compares the pressure distributions, loads, and center-of-pressure values.

The three cases are:

1. `1a: incompressible taper-flat`
2. `1b: incompressible wedge`
3. `2a: compressible taper-flat`

For each case, the script computes:

- Load `F`
- Nondimensional center-of-pressure `xp/L`
- Minimum relative pressure
- Maximum relative pressure
- Number of iterations

The script also generates the pressure distribution comparison figure. In the current version, the three plots are vertically stacked.

The plotted pressure height is scaled as:

```matlab
W*Pb/max(abs(Pb))
```

This means the vertical coordinate is a normalized pressure field multiplied by the slider width `W`. The scaling is used only for visualization, so the pressure surface has a height comparable to the bearing width.

## Recommended run order

Run the scripts in this order:

1. `part1a_load_center_pressure.m`
2. `part1b_wedge_convergence.m`
3. `part2a_compressible_taper_flat.m`
4. `part2bc_compare_cases.m`

This order matches the assignment structure and makes it easier to transfer results into the report.

## Figures to export for the report

Export the following MATLAB figures as PNG files and use the same names in the LaTeX report, or update the LaTeX file names accordingly.

Recommended file names:

- `part1b_convergence.png`
  - Log-log plot of maximum pressure error versus `Ex`
- `pressure_comparison_stacked.png`
  - Vertically stacked pressure distributions for the three configurations

Optional figure:

- `part1b_pressure_comparison.png`
  - Finest-grid FEM wedge pressure compared against the analytical wedge pressure

## Notes on pressure quantities

For incompressible cases, the pressure variable is the relative pressure `p`.

For the compressible case, the absolute pressure is:

```matlab
pa = p + patm;
```

However, the load and center-of-pressure should be computed using the relative pressure `p`, not `pa`. Integrating `pa` would add a uniform ambient-pressure contribution over the slider area, which is not the bearing-generated load.

## Notes on the wedge verification case

The wedge verification case is not meant to match the shape or magnitude of the taper-flat pressure field. It is a verification problem.

Differences from the taper-flat cases:

- The leading 1 degree taper is removed by setting `Ta = 0`
- Pressure is prescribed only at `x = 0` and `x = L`
- Top and bottom boundaries are Neumann boundaries
- The solution is nearly one-dimensional
- The pressure and load are much larger than in the taper-flat case

This is expected and is used to verify convergence against the analytical wedge-bearing solution.
