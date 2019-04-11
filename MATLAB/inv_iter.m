function [ z, P, iterations_actual ] = inv_iter( K, M, x, k, P)
    z0 = ones(size(M,1),1);         %Starting vector. 
    z0 = z0 / norm(z0);     %normalizing starting vector. 
    tol = 0.0001;             % accepted convergence accuracy tolerance. 
    epsilon = 1;            % error measurement
    
    
    if k == 1       % if first eigenmode.
        P = eye(size(M,1));         % Projection operator set to identity matrix
        
    else
        P = P * (eye(size(M,1)) - (x(:, k - 1) * x(:, k - 1).' * M ) / (x(:, k - 1).' * M * x(:, k - 1)));
    end 
        
    [K_L, K_U] = lu(K);         % Factorizing K to a lower and upper matrix.
    z = z0;     %z is initialized with the normalized, random vector. 
    iterations_actual = 0;
    %for p = 1:1:iterations
        while epsilon >= tol
            omega2_est_0 = (z.' * K * z) / (z.' * M * z);       % estimating Rayleigh quotient (3.4.15)
            y = M * z;             % eq. (3.4.23)
            b = P.' * y;                % want to solve Kz = P^T*y (3.4.24) with forward/backward substitution. 
            b_hat = K_L \ b;            % forward substitution (2.1.9)
            z_hat = K_U \ b_hat;        % backward substitution (2.1.10)
            z = P * z_hat;         % Projecting z on to the correct dimension.
            %z_hat = K \ P.' * y;
            %z = P * z_hat;
            
            z = z / norm(z);     % Normalizing z
            omega2_est_1 = (z.' * K * z) / (z.' * M * z);       % estimating Rayleigh quotient at p + 1 (3.4.15)
            epsilon = norm(omega2_est_1 - omega2_est_0);        % updating errror measurement.
            iterations_actual = iterations_actual + 1;
        end
    %end
    
    fprintf('----Results of inv_iter-------\n')
    fprintf('\nEigenmode no. %d\nEpsilon: %d\nIterations: %d\nz(1,1): %d\n', k, epsilon, iterations_actual, z(1,1));
    fprintf('-------------------------------\n')
                
end

