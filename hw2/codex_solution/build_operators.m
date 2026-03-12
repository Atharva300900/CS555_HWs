function [A, D, M, x_int, xb] = build_operators(N, nu, grid_type)
% BUILD_OPERATORS  Build FD operators for Burgers equation on [0,1]
%
%   [A, D, M, x_int, xb] = build_operators(N, nu, grid_type)
%
%   Inputs:
%     N         - number of grid intervals
%     nu        - viscosity
%     grid_type - 'uniform' or 'chebyshev'
%
%   Outputs:
%     A     - viscous operator = -nu * Laplacian, (N-1)x(N-1)
%     D     - first derivative operator, (N-1)x(N-1)
%     M     - mass matrix (identity for FD), (N-1)x(N-1)
%     x_int - interior grid points, (N-1)x1
%     xb    - all grid points including boundaries, (N+1)x1

    % -----------------------------
    % Build grid
    % -----------------------------
    if strcmpi(grid_type, 'uniform')
        xb = (0:N)'/N;
    elseif strcmpi(grid_type, 'chebyshev')
        j  = (0:N)';
        xb = (1 - cos(pi*j/N))/2;
    else
        error('grid_type must be ''uniform'' or ''chebyshev''.');
    end

    n = N - 1;              % number of interior points
    x_int = xb(2:end-1);    % interior grid points
    M = speye(n);

    % Preallocate sparse triplets
    iA = zeros(3*n,1); jA = zeros(3*n,1); vA = zeros(3*n,1); cA = 0;
    iD = zeros(3*n,1); jD = zeros(3*n,1); vD = zeros(3*n,1); cD = 0;

    % -----------------------------
    % Build A and D row by row
    % -----------------------------
    for i = 1:n
        % Interior point x_i corresponds to xb(i+1)
        ip = i + 1;

        dm = xb(ip)   - xb(ip-1);   % x_i   - x_{i-1}
        dp = xb(ip+1) - xb(ip);     % x_{i+1} - x_i

        % =========================================
        % Second derivative operator:
        % u_xx(x_i) ≈ cL*u_{i-1} + c0*u_i + cR*u_{i+1}
        %
        % cL =  2 / (dm*(dm+dp))
        % c0 = -2 / (dm*dp)
        % cR =  2 / (dp*(dm+dp))
        %
        % A = -nu * u_xx
        % =========================================
        cL =  2 / (dm*(dm+dp));
        c0 = -2 / (dm*dp);
        cR =  2 / (dp*(dm+dp));

        if i > 1
            cA = cA + 1;
            iA(cA) = i;
            jA(cA) = i-1;
            vA(cA) = -nu*cL;
        end

        cA = cA + 1;
        iA(cA) = i;
        jA(cA) = i;
        vA(cA) = -nu*c0;

        if i < n
            cA = cA + 1;
            iA(cA) = i;
            jA(cA) = i+1;
            vA(cA) = -nu*cR;
        end

        % =========================================
        % First derivative operator from the HW statement:
        %
        % u_x(x_i) ≈ (u_{i+1} - u_{i-1}) / (x_{i+1} - x_{i-1})
        %
        % This yields a tridiagonal operator with zero diagonal.
        % =========================================
        di = 1 / (dm + dp);
        dL = -di;
        d0 = 0;
        dR =  di;

        if i > 1
            cD = cD + 1;
            iD(cD) = i;
            jD(cD) = i-1;
            vD(cD) = dL;
        end

        cD = cD + 1;
        iD(cD) = i;
        jD(cD) = i;
        vD(cD) = d0;

        if i < n
            cD = cD + 1;
            iD(cD) = i;
            jD(cD) = i+1;
            vD(cD) = dR;
        end
    end

    % Assemble sparse matrices
    A = sparse(iA(1:cA), jA(1:cA), vA(1:cA), n, n);
    D = sparse(iD(1:cD), jD(1:cD), vD(1:cD), n, n);

end
