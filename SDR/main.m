clc
clear all

addpath('tensor_toolbox')
savepath
addpath("poblano_toolbox-1.2")
savepath
addpath("tensorlab_2016-03-28")
savepath

syms x0 x1 x2 x3 x4 x5 x6 x y z % list the variables here
syms posv positive


vars = [x0 x1 x2 x3 x4 x5 x6]; % make a vector consisting of all the variables
% vars = [x y z]
% p = x*y*z - 2*x*y + 3*x^2*z+5*y*z^2;
  p =  x0^{3}*x3^{2}+x0^{3}*x4^{2}+x0^{2}*x3*x4^{2}+x0^{2}*x4^{3}-2*x0^{3}*x3*x5-2*x0^{2}*x3^{2}*x5-2*x0^{2}*x3*x4*x5+x0^{2}*x4^{2}*x5+x0^{3}*x5^{2}+x0*x3^{2}*x5^{2}+2*x0^{2}*x4*x5^{2}+2*x0*x3*x4*x5^{2}+x0*x4^{2}*x5^{2} ...
 +2*x0^{2}*x5^{3}+2*x0*x3*x5^{3}+2*x0*x4*x5^{3}+x0*x5^{4}-2*x0^{3}*x4*x6-4*x0^{2}*x3*x4*x6-2*x0*x3^{2}*x4*x6-3*x0^{2}*x4^{2}*x6-4*x0*x3*x4^{2}*x6-2*x0*x4^{3}*x6-2*x0^{2}*x3*x5*x6 ...
 -4*x0^{2}*x4*x5*x6-4*x0*x3*x4*x5*x6-4*x0*x4^{2}*x5*x6+2*x0^{2}*x5^{2}*x6+2*x0*x3*x5^{2}*x6+2*x0*x5^{3}*x6+x0^{3}*x6^{2}+3*x0^{2}*x3*x6^{2}+3*x0*x3^{2}*x6^{2}+x3^{3}*x6^{2}- x0^{2}*x4*x6^{2} ...
     +2*x0*x3*x4*x6^{2}+3*x3^{2}*x4*x6^{2}- x0*x4^{2}*x6^{2}+3*x3*x4^{2}*x6^{2}+x4^{3}*x6^{2}+3*x0^{2}*x5*x6^{2}+6*x0*x3*x5*x6^{2}+3*x3^{2}*x5*x6^{2}+2*x0*x4*x5*x6^{2}+6*x3*x4*x5*x6^{2}+3*x4^{2}*x5*x6^{2} ...
      +4*x0*x5^{2}*x6^{2}+3*x3*x5^{2}*x6^{2}+3*x4*x5^{2}*x6^{2}+x5^{3}*x6^{2}+3*x0^{2}*x6^{3}+6*x0*x3*x6^{3}+3*x3^{2}*x6^{3}+4*x0*x4*x6^{3}+6*x3*x4*x6^{3}+3*x4^{2}*x6^{3}+6*x0*x5*x6^{3}+6*x3*x5*x6^{3} ...
     +6*x4*x5*x6^{3}+3*x5^{2}*x6^{3}+3*x0*x6^{4}+3*x3*x6^{4}+3*x4*x6^{4}+3*x5*x6^{4}+x6^{5};
% p = 7/25*x^2-y^2-48/25*x*z-7/25*z^2+1; % write the polynomial here
Max_error = 0.2; % to be chosen
d = 5; % to be chosen


n = length(vars);
% Exponents = [2 0 0; 0 2 0; 1 0 1; 0 0 2; 0 0 0];
% Coeffs = [7/25; -1; -48/25; -7/25; 1];
[c,t] = coeffs(p);
Coeffs = transpose(c); %Coefficient vector
[~,m] = size(c);
% Exponent matrix
Exponents = zeros(m,n); 
for i = 1:n
    vars_new = vars;
    vars_new(i) = [];
    only_one_var = subs(t, vars_new, ones(1,length(vars)-1));
    Exponents(:,i) = simplify(subs(log(only_one_var)./log(vars(i)),vars(i),posv));
end

[lambda,U] = tensor_decomposition(Exponents,Coeffs, Max_error,d);
X = transpose(U)*[transpose(vars);1]; % affine polynomials f_i.
q = simplify(sum(lambda.*(X.^d)));
error = eval(coeffs(p-q))
[r,~] = size(lambda) % This is equal to r (symmetric rank).

Pencil = SDR(lambda, d);
[~,s,~] = size(Pencil);
s
% A1 = zeros(s,s);
% A2 = zeros(s,s);
% A3 = zeros(s,s);
% A4 = zeros(s,s);
% A5 = zeros(s,s);
% A6 = zeros(s,s);
% A1(:,:) = Pencil(1,:,:);
% A2(:,:) = Pencil(2,:,:);
% A3(:,:) = Pencil(3,:,:);
% A4(:,:) = Pencil(4,:,:);
% A5(:,:) = Pencil(5,:,:);
% A6(:,:) = Pencil(6,:,:);
A = zeros(s,s);
for i = 1:r
    B = zeros(s,s);
    B(:,:) = Pencil(i,:,:);
    A = A + X(i)*B;
end
B = zeros(s,s);
B(:,:) = Pencil(r+1,:,:);
A = A+B;

% for i=1:s
%     for j=1:s
%         A(i,j) = X(1)*Pencil(1,i,j) + X(2)*Pencil(2,i,j) + X(3)*Pencil(3,i,j) + Pencil(4,i,j);
%     end
% end
% A = X(1)*Pencil(1,:,:) + X(2)*Pencil(2,:,:) + X(3)*Pencil(3,:,:) + Pencil(4,:,:);
% qq = det( X(1)*A1 + X(2)*A2 + X(3)*A3 + X(4)* A4 + X(5)*A5 + A6)
% qq = det(A);
% error2 = eval(coeffs(p-qq))

test = zeros(1,100);
 for i=1:100
     x = randn(1,7);
     a = subs(p,vars,x);
     b = det(subs(A,vars,x));
     test(i) = eval(abs(a-b));
 end

 test






