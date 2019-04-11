function [ x ] = LDL_solv( A, b )
% Solves the system Ax = b using LDL^T factorization and forwards and
% backwards substitution. 

% A: Non-singular, symmetric, square matrix (n x n)
% b: Size (n x 1)
% x: Resulting vector, size (n x 1)

x = A\b;


%{
x = zeros(size(A),1);
x(:,1) = 1;                 % x is originally set to all zeroes. 


[L, D] = ldl(A);            % LDL' factorization. A = LDL'
U = D * L';                 % Using notation as in LU factorization (pg. 43)



b_hat = U * x;              % Eq. (2.1.8)

% Solving L * b_hat = b: (2.1.9)

for i = 1:1:size(A,1)
   
    for j = 1:1:size(A,2)
        k = i - 1;
        
        l_ik = U(i, k) / U(k,k);                        % eq. (2.1.3)
        b_hat(i) = b_hat(i) - l_ik * b_hat(k);          % eq. (2.1.4) 
    end
end

x = zeros(size(A, 1), 1);


%}

end

