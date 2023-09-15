% For a polynomial p of degree at most d with m monomials in n variables, Exponents:
%an m*n matrix, Coeffs: an m*1 vector
% Returns lambda, U such that p = \sum_{i=1}^r lambda(i) * ith column of U^d.
function [lambda,U] = sym_decom(Exponents,Coeffs, Max_error,d) 

[m,n] = size(Exponents); % m: number of monomials, n: number of variables

% Define the coefficient tensor
T = symtensor(@zeros, d, n+1);
for i = 1:m
    index = zeros(1,d);
    for j = 1:n
        index(1, sum(Exponents(i,1:j-1))+1: sum(Exponents(i,1:j))) = repmat(j, 1, Exponents(i,j));
    end
    index(1, sum(Exponents(i,:))+1: d) = repmat(n+1, 1, d-sum(Exponents(i,:)));
    T(index) = Coeffs(i)*prod(factorial(Exponents(i,:)))*factorial(d-sum(Exponents(i,:)))/factorial(d);
end
T = full(T);

Max_T = max(abs(T(:)));
error = Max_error+1;
r = 0; %Initial assumed symmetric rank of T

optparams = lbfgs('defaults'); % Get the optimization parameters
optparams.RelFuncTol = 1e-10; % Tighten the stopping tolerance
optparams.StopTol = 1e-6; % Tighten the stopping tolerance
rng(5); % Set random number generator state for consistent results


while error >= Max_error
    r = r+1;
    [S,~] = cp_sym(T,r,'unique',false,'l1param',0,'alg_options',optparams); % Find the symmetric decomposition of T with rank r
    T_obtained = double(symktensor(S.lambda, S.U, d)); % Compute the actual tensor euqal to the symmetric CP decomposition obtained
    E = T - T_obtained;
    % error = max(abs(E(:)))/Max_T; 
    error = max(abs(E(:)));
end

lambda = S.lambda;
U = S.U;
end