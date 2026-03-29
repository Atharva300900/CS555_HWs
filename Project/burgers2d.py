#!/usr/bin/env python3
"""
================================================================================
Mini-App: 2D Viscous Burgers' Equation
================================================================================
CS555: Numerical Methods for PDEs — Spring 2026 Final Project

Solves the coupled 2D viscous Burgers' equations:

    u_t + u*u_x + v*u_y = nu*(u_xx + u_yy)
    v_t + u*v_x + v*v_y = nu*(v_xx + v_yy)

on [0,1]^2 with Dirichlet BCs from the Cole-Hopf exact solution.

Features
--------
- Spatial:  Central differences (2nd order), first-order upwind
- Temporal: Forward Euler, RK4, Backward Euler, Crank-Nicolson
- Verification via Cole-Hopf analytical solution
- Convergence studies (space and time)
- Stability analysis
- Cost vs. accuracy comparison
- Viscosity sweep (nu -> 0)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.sparse import kron, eye, diags, csc_matrix
from scipy.sparse.linalg import spsolve, factorized
import time as timer


# ================================================================
# 1. EXACT SOLUTION  (Cole-Hopf transformation)
# ================================================================

def exact_solution(X, Y, t, nu):
    """
    Cole-Hopf exact solution for the 2D Burgers' system.

        u(x,y,t) = 3/4 - 1 / (4*(1 + exp((-4x + 4y - t) / (32*nu))))
        v(x,y,t) = 3/4 + 1 / (4*(1 + exp((-4x + 4y - t) / (32*nu))))

    Reference: Fletcher (1983), Bahadir (2003).
    """
    eta = (-4.0 * X + 4.0 * Y - t) / (32.0 * nu)
    eta = np.clip(eta, -500, 500)
    e = np.exp(eta)
    u = 0.75 - 0.25 / (1.0 + e)
    v = 0.75 + 0.25 / (1.0 + e)
    return u, v


# ================================================================
# 2. GRID
# ================================================================

def make_grid(N):
    """Uniform grid on [0,1]^2 with (N+1)^2 points.  Returns X, Y, h."""
    h = 1.0 / N
    x = np.linspace(0, 1, N + 1)
    X, Y = np.meshgrid(x, x, indexing='ij')
    return X, Y, h


# ================================================================
# 3. SPATIAL OPERATORS  (operate on full arrays, return interior)
# ================================================================

def laplacian_fd(phi, h):
    """Central-difference Laplacian at interior points."""
    return ((phi[2:, 1:-1] - 2 * phi[1:-1, 1:-1] + phi[:-2, 1:-1])
          + (phi[1:-1, 2:] - 2 * phi[1:-1, 1:-1] + phi[1:-1, :-2])) / h**2


def advection_central(u, v, phi, h):
    """Central-difference advection  u*phi_x + v*phi_y  (2nd order)."""
    px = (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2 * h)
    py = (phi[1:-1, 2:] - phi[1:-1, :-2]) / (2 * h)
    return u[1:-1, 1:-1] * px + v[1:-1, 1:-1] * py


def advection_upwind(u, v, phi, h):
    """First-order upwind advection  u*phi_x + v*phi_y."""
    ui = u[1:-1, 1:-1]
    vi = v[1:-1, 1:-1]

    px_f = (phi[2:, 1:-1] - phi[1:-1, 1:-1]) / h      # forward
    px_b = (phi[1:-1, 1:-1] - phi[:-2, 1:-1]) / h      # backward
    px = np.where(ui >= 0, px_b, px_f)

    py_f = (phi[1:-1, 2:] - phi[1:-1, 1:-1]) / h
    py_b = (phi[1:-1, 1:-1] - phi[1:-1, :-2]) / h
    py = np.where(vi >= 0, py_b, py_f)

    return ui * px + vi * py


def full_rhs(u, v, nu, h, adv_type='central'):
    """Full explicit RHS = -advection + nu*Laplacian  for both components."""
    adv = advection_central if adv_type == 'central' else advection_upwind
    ru = -adv(u, v, u, h) + nu * laplacian_fd(u, h)
    rv = -adv(u, v, v, h) + nu * laplacian_fd(v, h)
    return ru, rv


# ================================================================
# 4. LAPLACIAN MATRIX  (for implicit methods)
# ================================================================

def build_laplacian(m, h):
    """
    2D Laplacian for m x m interior points via Kronecker products.
    Returns sparse CSC matrix of size m^2 x m^2.
    Ordering: column-major (Fortran) flattening of the m x m interior block.
    """
    e = np.ones(m)
    D = diags([e, -2 * e, e], [-1, 0, 1], shape=(m, m)) / h**2
    I = eye(m)
    return csc_matrix(kron(I, D) + kron(D, I))


def bc_contribution(u_bc, v_bc, m, h):
    """Boundary-value contributions to the Laplacian RHS for interior points."""
    h2 = h**2
    bu = np.zeros((m, m))
    bv = np.zeros((m, m))

    bu[0, :]  += u_bc[0, 1:-1]  / h2      # left   x = 0
    bu[-1, :] += u_bc[-1, 1:-1] / h2      # right  x = 1
    bu[:, 0]  += u_bc[1:-1, 0]  / h2      # bottom y = 0
    bu[:, -1] += u_bc[1:-1, -1] / h2      # top    y = 1

    bv[0, :]  += v_bc[0, 1:-1]  / h2
    bv[-1, :] += v_bc[-1, 1:-1] / h2
    bv[:, 0]  += v_bc[1:-1, 0]  / h2
    bv[:, -1] += v_bc[1:-1, -1] / h2

    return bu, bv


# ================================================================
# 5. BOUNDARY CONDITIONS
# ================================================================

def apply_bc(u, v, X, Y, t, nu):
    """Set Dirichlet BCs from the exact solution on all four edges."""
    ue, ve = exact_solution(X, Y, t, nu)
    u[:, 0]  = ue[:, 0];   v[:, 0]  = ve[:, 0]     # bottom
    u[:, -1] = ue[:, -1];  v[:, -1] = ve[:, -1]     # top
    u[0, :]  = ue[0, :];   v[0, :]  = ve[0, :]      # left
    u[-1, :] = ue[-1, :];  v[-1, :] = ve[-1, :]     # right


# ================================================================
# 6. TIME INTEGRATORS
# ================================================================

def step_euler(u, v, nu, h, dt, X, Y, t, adv_type='central', **kw):
    """Forward Euler — fully explicit, O(dt) in time."""
    ru, rv = full_rhs(u, v, nu, h, adv_type)
    un = u.copy();  vn = v.copy()
    un[1:-1, 1:-1] += dt * ru
    vn[1:-1, 1:-1] += dt * rv
    apply_bc(un, vn, X, Y, t + dt, nu)
    return un, vn


def step_rk4(u, v, nu, h, dt, X, Y, t, adv_type='central', **kw):
    """Classical RK4 — fully explicit, O(dt^4) in time."""
    def rhs(uu, vv):
        return full_rhs(uu, vv, nu, h, adv_type)

    def stage(u0, v0, ku, kv, a):
        us = u0.copy();  vs = v0.copy()
        us[1:-1, 1:-1] += a * dt * ku
        vs[1:-1, 1:-1] += a * dt * kv
        apply_bc(us, vs, X, Y, t + a * dt, nu)
        return us, vs

    k1u, k1v = rhs(u, v)
    k2u, k2v = rhs(*stage(u, v, k1u, k1v, 0.5))
    k3u, k3v = rhs(*stage(u, v, k2u, k2v, 0.5))
    k4u, k4v = rhs(*stage(u, v, k3u, k3v, 1.0))

    un = u.copy();  vn = v.copy()
    un[1:-1, 1:-1] += (dt / 6) * (k1u + 2*k2u + 2*k3u + k4u)
    vn[1:-1, 1:-1] += (dt / 6) * (k1v + 2*k2v + 2*k3v + k4v)
    apply_bc(un, vn, X, Y, t + dt, nu)
    return un, vn


def step_be(u, v, nu, h, dt, X, Y, t, adv_type='central',
            L=None, solve_fn=None, **kw):
    """
    Backward Euler — semi-implicit.
    Advection explicit at t^n, diffusion implicit at t^{n+1}.

        (I - dt*nu*L) u^{n+1} = u^n - dt*adv^n + dt*nu*bc^{n+1}
    """
    m = u.shape[0] - 2
    adv = advection_central if adv_type == 'central' else advection_upwind

    au = adv(u, v, u, h)
    av = adv(u, v, v, h)

    # Boundary values at t+dt
    ubc = np.zeros_like(u);  vbc = np.zeros_like(v)
    apply_bc(ubc, vbc, X, Y, t + dt, nu)
    bu, bv = bc_contribution(ubc, vbc, m, h)

    # RHS vectors (column-major flatten to match Kronecker ordering)
    ru = (u[1:-1, 1:-1] - dt * au + dt * nu * bu).flatten('F')
    rv = (v[1:-1, 1:-1] - dt * av + dt * nu * bv).flatten('F')

    # Solve
    ui = solve_fn(ru).reshape((m, m), order='F')
    vi = solve_fn(rv).reshape((m, m), order='F')

    un = ubc.copy();  vn = vbc.copy()
    un[1:-1, 1:-1] = ui
    vn[1:-1, 1:-1] = vi
    return un, vn


def step_cn(u, v, nu, h, dt, X, Y, t, adv_type='central',
            L=None, solve_fn=None, **kw):
    """
    Crank-Nicolson with predictor-corrector — O(dt^2) overall.

    Predictor  (Forward Euler):
        u* = u^n + dt * (-adv^n + nu*lap(u^n))

    Corrector  (trapezoidal advection + CN diffusion):
        (I - dt/2*nu*L) u^{n+1} = u^n - dt/2*(adv^n + adv(u*))
                                   + dt/2*nu*lap(u^n)
                                   + dt/2*nu*bc^{n+1}

    The predictor gives a second-order estimate of advection at t^{n+1},
    so averaging adv^n and adv(u*) yields a trapezoidal (2nd-order)
    treatment of the nonlinear advection term.
    """
    m = u.shape[0] - 2
    adv = advection_central if adv_type == 'central' else advection_upwind

    # ── stage 1: quantities at t^n ──
    au_n = adv(u, v, u, h)
    av_n = adv(u, v, v, h)
    du_n = laplacian_fd(u, h)
    dv_n = laplacian_fd(v, h)

    # ── stage 2: Forward-Euler predictor ──
    u_star = u.copy();  v_star = v.copy()
    u_star[1:-1, 1:-1] += dt * (-au_n + nu * du_n)
    v_star[1:-1, 1:-1] += dt * (-av_n + nu * dv_n)
    apply_bc(u_star, v_star, X, Y, t + dt, nu)

    # Advection evaluated at predicted state
    au_s = adv(u_star, v_star, u_star, h)
    av_s = adv(u_star, v_star, v_star, h)

    # ── stage 3: corrector solve ──
    ubc = np.zeros_like(u);  vbc = np.zeros_like(v)
    apply_bc(ubc, vbc, X, Y, t + dt, nu)
    bu, bv = bc_contribution(ubc, vbc, m, h)

    c = 0.5 * dt * nu
    ru = (u[1:-1, 1:-1]
          - 0.5 * dt * (au_n + au_s)
          + c * du_n + c * bu).flatten('F')
    rv = (v[1:-1, 1:-1]
          - 0.5 * dt * (av_n + av_s)
          + c * dv_n + c * bv).flatten('F')

    ui = solve_fn(ru).reshape((m, m), order='F')
    vi = solve_fn(rv).reshape((m, m), order='F')

    un = ubc.copy();  vn = vbc.copy()
    un[1:-1, 1:-1] = ui
    vn[1:-1, 1:-1] = vi
    return un, vn


METHODS = {
    'euler': step_euler,
    'rk4':   step_rk4,
    'be':    step_be,
    'cn':    step_cn,
}

METHOD_NAMES = {
    'euler': 'Forward Euler',
    'rk4':   'RK4',
    'be':    'Backward Euler',
    'cn':    'Crank-Nicolson',
}


# ================================================================
# 7. SOLVER
# ================================================================

def solve(N, nu, T, dt, method='rk4', adv_type='central'):
    """
    Solve 2D Burgers' on [0,1]^2 from t=0 to t=T.

    Parameters
    ----------
    N        : int    — grid points per side (N+1 total, N-1 interior)
    nu       : float  — viscosity
    T        : float  — final time
    dt       : float  — time-step size
    method   : str    — 'euler', 'rk4', 'be', 'cn'
    adv_type : str    — 'central' or 'upwind'

    Returns
    -------
    u, v     : final numerical solution  (N+1) x (N+1)
    X, Y     : meshgrid arrays
    info     : dict with err, wall time, nsteps, etc.
    """
    X, Y, h = make_grid(N)
    m = N - 1                         # interior points per direction

    # Initial condition from exact solution
    u, v = exact_solution(X, Y, 0.0, nu)

    # Number of steps (adjust dt so we land exactly on T)
    nsteps = max(1, int(np.round(T / dt)))
    dt_actual = T / nsteps

    # Pre-build & pre-factor the implicit matrix (constant across steps)
    L = None
    solve_fn = None
    if method == 'be':
        L = build_laplacian(m, h)
        A = eye(m * m, format='csc') - dt_actual * nu * L
        solve_fn = factorized(A)
    elif method == 'cn':
        L = build_laplacian(m, h)
        A = eye(m * m, format='csc') - 0.5 * dt_actual * nu * L
        solve_fn = factorized(A)

    step = METHODS[method]

    # Time loop
    t = 0.0
    tic = timer.perf_counter()
    for _ in range(nsteps):
        u, v = step(u, v, nu, h, dt_actual, X, Y, t,
                    adv_type=adv_type, L=L, solve_fn=solve_fn)
        t += dt_actual
    wall = timer.perf_counter() - tic

    # Error vs exact solution
    ue, ve = exact_solution(X, Y, T, nu)
    err_u = np.sqrt(h**2 * np.sum((u - ue)**2))       # discrete L2
    err_v = np.sqrt(h**2 * np.sum((v - ve)**2))
    linf_u = np.max(np.abs(u - ue))                    # L-infinity
    linf_v = np.max(np.abs(v - ve))

    info = dict(
        err_u=err_u, err_v=err_v,
        err=max(err_u, err_v),
        linf=max(linf_u, linf_v),
        wall=wall, nsteps=nsteps,
        dt=dt_actual, h=h, N=N,
    )
    return u, v, X, Y, info


# ================================================================
# 8. CONVERGENCE STUDIES
# ================================================================

def convergence_space(Ns, nu=0.01, T=0.5, method='rk4', adv_type='central'):
    """
    Spatial convergence: vary N with dt small enough that time error
    is negligible.
    """
    results = []
    for N in Ns:
        h = 1.0 / N
        # dt ~ h^3 keeps O(h^3) temporal error below O(h^2) spatial error
        if method == 'euler':
            dt = 0.25 * h**2 / (4 * nu + 1.0)   # within stability limit
        elif method == 'rk4':
            dt = 0.5 * h**2 / (4 * nu + 1.0)
        else:
            dt = h**2                              # implicit — no stability limit
        _, _, _, _, info = solve(N, nu, T, dt, method, adv_type)
        results.append(info)
        print(f"    N={N:4d}  h={h:.4e}  dt={info['dt']:.2e}  "
              f"err={info['err']:.4e}  steps={info['nsteps']:6d}  "
              f"wall={info['wall']:.2f}s")
    return results


def convergence_time(dts, N=80, nu=0.01, T=0.5, method='rk4', adv_type='central'):
    """Temporal convergence: fix N (large), vary dt."""
    results = []
    for dt in dts:
        try:
            _, _, _, _, info = solve(N, nu, T, dt, method, adv_type)
            results.append(info)
            print(f"    dt={dt:.2e}  err={info['err']:.4e}  "
                  f"steps={info['nsteps']:6d}  wall={info['wall']:.2f}s")
        except Exception as e:
            print(f"    dt={dt:.2e}  FAILED: {e}")
    return results


# ================================================================
# 9. STABILITY STUDY
# ================================================================

def stability_boundary(N, nu, T, method='euler', adv_type='central', n_test=20):
    """
    Probe stability of Forward Euler by testing a range of dt values
    around the theoretical diffusion stability limit  dt <= h^2 / (4*nu).

    Returns list of (dt, ratio, stable_bool)  and  dt_limit.
    """
    h = 1.0 / N
    dt_limit = h**2 / (4 * nu)        # theoretical diffusion limit (2D)

    ratios = np.linspace(0.1, 2.0, n_test)
    results = []
    for r in ratios:
        dt_test = r * dt_limit
        X, Y, h_ = make_grid(N)
        u, v = exact_solution(X, Y, 0, nu)
        step = METHODS[method]

        stable = True
        t = 0.0
        nsteps = max(1, int(T / dt_test))
        for _ in range(nsteps):
            u, v = step(u, v, nu, h_, dt_test, X, Y, t, adv_type=adv_type)
            t += dt_test
            if np.any(np.isnan(u)) or np.max(np.abs(u)) > 10:
                stable = False
                break

        results.append((dt_test, r, stable))
    return results, dt_limit


# ================================================================
# 10. VISUALIZATION
# ================================================================

def plot_solution_3d(u, v, X, Y, title='', save=None):
    """Side-by-side 3-D surface plots for u and v."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             subplot_kw={'projection': '3d'})
    for ax, field, name in zip(axes, [u, v], ['u', 'v']):
        ax.plot_surface(X, Y, field, cmap=cm.viridis, alpha=0.9,
                        rstride=max(1, X.shape[0]//40),
                        cstride=max(1, X.shape[1]//40), linewidth=0)
        ax.set_xlabel('x');  ax.set_ylabel('y');  ax.set_zlabel(name)
        ax.set_title(f'{name}(x, y)')
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_contours(u, v, X, Y, title='', save=None):
    """Side-by-side filled contour plots."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    c1 = ax1.contourf(X, Y, u, levels=30, cmap=cm.viridis)
    ax1.set_xlabel('x');  ax1.set_ylabel('y');  ax1.set_title('u(x, y)')
    plt.colorbar(c1, ax=ax1)

    c2 = ax2.contourf(X, Y, v, levels=30, cmap=cm.viridis)
    ax2.set_xlabel('x');  ax2.set_ylabel('y');  ax2.set_title('v(x, y)')
    plt.colorbar(c2, ax=ax2)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_convergence(results_dict, xkey='h', title='', save=None):
    """Log-log convergence plot for multiple method/advection combos."""
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    for i, (label, results) in enumerate(results_dict.items()):
        if not results:
            continue
        x = np.array([r[xkey] for r in results])
        y = np.array([r['err'] for r in results])
        ax.loglog(x, y, f'-{markers[i % len(markers)]}',
                  label=label, markersize=7, linewidth=1.5)

    # Reference slopes
    first_results = list(results_dict.values())[0]
    if first_results:
        x0 = np.array([r[xkey] for r in first_results])
        if xkey == 'h':
            ax.loglog(x0, 0.5 * (x0 / x0[0])**1 * first_results[0]['err'],
                      'k--', alpha=0.3, label='O(h)')
            ax.loglog(x0, 0.5 * (x0 / x0[0])**2 * first_results[0]['err'],
                      'k-.', alpha=0.3, label=r'O(h$^2$)')
        elif xkey == 'dt':
            ax.loglog(x0, (x0 / x0[0])**1 * first_results[0]['err'],
                      'k--', alpha=0.3, label=r'O($\Delta t$)')
            ax.loglog(x0, (x0 / x0[0])**2 * first_results[0]['err'],
                      'k-.', alpha=0.3, label=r'O($\Delta t^2$)')
            ax.loglog(x0, (x0 / x0[0])**4 * first_results[0]['err'],
                      'k:', alpha=0.3, label=r'O($\Delta t^4$)')

    ax.set_xlabel(xkey, fontsize=12)
    ax.set_ylabel('L2 Error', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_cost_vs_accuracy(results_dict, title='Cost vs. Accuracy', save=None):
    """Wall time vs. error for multiple methods."""
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    for i, (label, results) in enumerate(results_dict.items()):
        if not results:
            continue
        walls = [r['wall'] for r in results]
        errs  = [r['err'] for r in results]
        ax.loglog(walls, errs, f'-{markers[i % len(markers)]}',
                  label=label, markersize=7, linewidth=1.5)

    ax.set_xlabel('Wall time (s)', fontsize=12)
    ax.set_ylabel('L2 Error', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_viscosity_sweep(nus, N=80, T=0.5, save=None):
    """Contour plots of u at different viscosities."""
    ncols = len(nus)
    fig, axes = plt.subplots(2, ncols, figsize=(4 * ncols, 7))
    if ncols == 1:
        axes = axes.reshape(2, 1)

    for j, nu_val in enumerate(nus):
        h = 1.0 / N
        dt = min(0.001, 0.4 * h**2 / (4 * nu_val + 1.0))
        u, v, X, Y, info = solve(N, nu_val, T, dt, 'rk4')
        print(f"  nu={nu_val:.4f}: err={info['err']:.4e}, wall={info['wall']:.2f}s")

        c1 = axes[0, j].contourf(X, Y, u, levels=30, cmap=cm.viridis)
        axes[0, j].set_title(f'u,  nu={nu_val}')
        axes[0, j].set_xlabel('x');  axes[0, j].set_ylabel('y')
        plt.colorbar(c1, ax=axes[0, j])

        c2 = axes[1, j].contourf(X, Y, v, levels=30, cmap=cm.viridis)
        axes[1, j].set_title(f'v,  nu={nu_val}')
        axes[1, j].set_xlabel('x');  axes[1, j].set_ylabel('y')
        plt.colorbar(c2, ax=axes[1, j])

    fig.suptitle(f'Effect of Viscosity at T={T}', fontsize=14)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_stability(stab_results, dt_limit, title='', save=None):
    """Bar chart showing stable / unstable dt values."""
    fig, ax = plt.subplots(figsize=(10, 4))
    for dt_val, ratio, stable in stab_results:
        color = 'green' if stable else 'red'
        ax.bar(ratio, 1, width=0.08, color=color, alpha=0.7)
    ax.axvline(x=1.0, color='black', linestyle='--', label=f'dt_limit={dt_limit:.2e}')
    ax.set_xlabel('dt / dt_limit', fontsize=12)
    ax.set_yticks([])
    ax.set_title(title or 'Forward Euler Stability', fontsize=14)
    ax.legend(fontsize=10)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


# ================================================================
# 11. MAIN DRIVER
# ================================================================

def main():
    nu = 0.01
    T  = 0.5

    # ── A. Reference solution ────────────────────────────────────
    print("=" * 65)
    print("A. Reference solution  (RK4, N=80, central)")
    print("=" * 65)
    u, v, X, Y, info = solve(80, nu, T, dt=1e-4, method='rk4')
    print(f"   L2 err = {info['err']:.4e},  Linf = {info['linf']:.4e},  "
          f"wall = {info['wall']:.2f}s")

    ue, ve = exact_solution(X, Y, T, nu)
    plot_solution_3d(u, v, X, Y,
                     title=f'Numerical solution  T={T}, nu={nu}',
                     save='fig_solution_3d.png')
    plot_contours(u, v, X, Y,
                  title=f'Numerical solution  T={T}, nu={nu}',
                  save='fig_solution_contour.png')
    plot_contours(u - ue, v - ve, X, Y,
                  title='Pointwise error  (numerical - exact)',
                  save='fig_error_contour.png')

    # ── B. Spatial convergence ───────────────────────────────────
    print("\n" + "=" * 65)
    print("B. Spatial convergence")
    print("=" * 65)
    Ns = [10, 20, 40, 80, 160]
    space = {}
    for meth in ['rk4', 'cn']:
        for adv in ['central', 'upwind']:
            key = f'{METHOD_NAMES[meth]} + {adv}'
            print(f"\n  {key}:")
            space[key] = convergence_space(Ns, nu, T, meth, adv)

    plot_convergence(space, xkey='h',
                     title='Spatial Convergence (L2 error vs h)',
                     save='fig_convergence_space.png')

    # ── C. Temporal convergence ──────────────────────────────────
    print("\n" + "=" * 65)
    print("C. Temporal convergence  (N=80 fixed)")
    print("=" * 65)
    N_ref = 80
    dts_ex  = [0.005, 0.002, 0.001, 0.0005, 0.0002]
    dts_im  = [0.05, 0.02, 0.01, 0.005, 0.002, 0.001]

    time_res = {}
    for meth, dts in [('euler', dts_ex), ('rk4', dts_ex),
                      ('be', dts_im), ('cn', dts_im)]:
        key = METHOD_NAMES[meth]
        print(f"\n  {key}:")
        time_res[key] = convergence_time(dts, N_ref, nu, T, meth)

    plot_convergence(time_res, xkey='dt',
                     title='Temporal Convergence (L2 error vs dt)',
                     save='fig_convergence_time.png')

    # ── D. Cost vs accuracy ──────────────────────────────────────
    print("\n" + "=" * 65)
    print("D. Cost vs. accuracy")
    print("=" * 65)
    plot_cost_vs_accuracy(time_res, save='fig_cost_vs_accuracy.png')

    # ── E. Viscosity sweep ───────────────────────────────────────
    print("\n" + "=" * 65)
    print("E. Viscosity sweep")
    print("=" * 65)
    plot_viscosity_sweep([0.1, 0.01, 0.005, 0.001], N=80, T=T,
                         save='fig_viscosity_sweep.png')

    # ── F. Stability boundary (Forward Euler) ────────────────────
    print("\n" + "=" * 65)
    print("F. Stability boundary  (Forward Euler, N=80)")
    print("=" * 65)
    stab, dt_lim = stability_boundary(80, nu, T, 'euler')
    for dt_val, ratio, ok in stab:
        tag = "STABLE" if ok else "UNSTABLE"
        print(f"    dt/dt_limit = {ratio:.2f}  ->  {tag}")
    plot_stability(stab, dt_lim, save='fig_stability.png')

    print("\n" + "=" * 65)
    print("All done.  Figures saved to current directory.")
    print("=" * 65)


if __name__ == '__main__':
    main()
