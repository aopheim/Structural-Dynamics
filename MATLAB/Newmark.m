function [q, q_dot, q_dot_dot, h, num_t_entries, t] = Newmark(gamma, beta, K, M, C, h, q0, q0_dot, omega)
% Creating vector p, containing the triangular force subjected to the
% structure and time vector. Collected in function for easy reuse.

deltaT = 10.5;            %Time after the force is applied [s]     


[p, num_t_entries, t] = createForceAndTime(h, omega, 100, deltaT );
p_new = zeros(size(M,1), num_t_entries);            % The force must be a coloumn vector, with force only in the node subjected to force (node 51)
p_new(51,:) = p';
p = -p_new;


% Defining all needed matrices
q_star = zeros(size(K,1), num_t_entries);        % q of size (1 x n)
q = zeros(size(K,1), num_t_entries);        % q of size (1 x n)
q(1,1) = q0;

q_dot_star = zeros(size(K,1), num_t_entries);        % q_dot of size (1 x n)
q_dot = zeros(size(K,1), num_t_entries);        % q_dot of size (1 x n)
q_dot(1,1) = q0_dot;

q_dot_dot = zeros(size(K,1), num_t_entries);        % q_dot_dot of size (1 x n)
%q_dot_dot_imp(1,1) = q0_dot_dot_imp(1,1);

S = M + h * gamma * C + h^2 * beta * K;
[S_L, S_U] = lu(S);


for i = 2 : 1 : num_t_entries       % TODO: change to half the period of the fourth eigenfrequency.
    % Time incrementation (done at the start of every for loop)
    
    % Prediction
    q_dot_star(:, i) = q_dot(:, i - 1) + ((1 - gamma) * h * q_dot_dot(:, i - 1));
    q_star(:,i) = q(:, i - 1) + (h * q_dot(:, i - 1)) + (0.5 - beta) * h^2 * q_dot_dot(:, i - 1);
    
    % Acceleration computing
    % S, S_L and S_K is computed before the for loop.
    
                                    % Solving q_dot_dot(:,i) = S \ ( p(i,1)
                                    % - K * q_star(:,i)) with LU
                                    % factorization
    b = p(:,i) - (C * q_dot_star(:,i)) - K * q_star(:,i);
    b_hat = S_L \ b;        %Forward substitution
    q_dot_dot(:,i) = S_U \ b_hat;           %Backwards substitution.
        
    % Correction
    q_dot(:, i) = q_dot_star(:,i) + h * gamma * q_dot_dot(:,i);
    q(:,i) = q_star(:,i) + h^2 * beta * q_dot_dot(:, i);
end

fprintf('----Results of Newmark-------\n')
fprintf('\nno_t_entries: %d\nh: %d\n', num_t_entries, h);
fprintf('-------------------------------\n')
       

end

