% returns SDR for the polynomial q=\sum_{k=1}^r lambda_k y_k^d.
% q = \sum_{k=1}^r y_k Pencil(k,:,:) + Pencil(r+1,:,:).
function Pencil = SDR(lambda, d)
[r,~] = size(lambda); % r: number of monomials 
m = floor(log2(d))+1; % m: number of rounds of simple substitutions
base2 = dec2base(d,2);

% d = sum_{l=1}^M 2^{p_l}.
p = [];
for i = 1:length(base2)
    if base2(i) =='1'
        p = [p, length(base2)-i];
    end
end
M = length(p);

% T is a 4-way tensor whose (k,j,:,:)th element is the coefficient of the
% variable u_j^(k) in the SDR of p, k \in [r], j \in [M].
max_SDR = 2^(m-1)*(r^m)*factorial(m+2);
T = zeros(r,M,max_SDR,max_SDR);

A0 = zeros(max_SDR); % matrix A_0 in the SDR (coefficient of 1).

% SDR of \sum_{k=1}^r lambda_k u_1^(k):
for k = 1:r
    T(k,1,1,1) = lambda(k);
end

s = 1; % records the SDR size

for i = 1:m
    s_new = s; % stores the new size of SDR
    deter = 1; % stores the determinant of the (2,2)-block of the matrix that's being constructed
    for j = 1:min(i,M)
        if (j==i && M>=j+1)
            for k = 1:r
                C = zeros(s,s);
                C(:,:) = T(k,j,1:s,1:s);
                inverse = C\eye(s);
                if (sum(sum(inverse==Inf | inverse==-Inf))>=1) %if T(k,j,1:s,1:s) is not invertible
                    [~,E] = eig(C);
                    E = max(abs(E),[],"all")+1; % E is non-zero and not an eigenvalue of C.
                    T(k,j,:,:) = zeros(max_SDR);
                    T(k,j+1,:,:) = zeros(max_SDR);

                    % 1
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    Z = zeros(max_SDR);
                    D = (C-E*eye(s))\eye(s);
                    deter_D = det(D);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = -D;
                    A0 = A0 + Z;
                    deter = deter * (-deter_D);

                    s_new = s_new + s;

                    % 2
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = -0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = -0.5 * eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = D;
                    A0 = A0 + Z;
                    deter = deter*deter_D;

                    s_new = s_new + s;

                    % 3
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = -eye(s)/E;
                    A0 = A0 + Z;
                    deter = deter*((-1)^s)/(E^s);

                    s_new = s_new + s;

                    % 4
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = -0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = -0.5 * eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = eye(s)/E;
                    A0 = A0 + Z;
                    deter = deter/(E^s);

                    s_new = s_new + s;

                else
                    T(k,j,:,:) = zeros(max_SDR);
                    T(k,j+1,:,:) = zeros(max_SDR);
                    % 1
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = -inverse;
                    A0 = A0 + Z;
                    deter = deter*(-det(inverse));

                    s_new = s_new + s;

                    % 2
                    T(k, j, 1:s, s_new+1:s_new+s) = 0.5 * eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = 0.5 * eye(s);

                    T(k, j+1, 1:s, s_new+1:s_new+s) = -0.5 * eye(s);
                    T(k, j+1, s_new+1:s_new+s, 1:s) = -0.5 * eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = inverse;
                    A0 = A0 + Z;
                    deter = deter*det(inverse);

                    s_new = s_new + s;
                end
            end
        elseif ((j==M && p(M)>=i-M+1) || (j<min(i,M) && p(j)>=i-j))
            for k = 1:r
                C = zeros(s,s);
                C(:,:) = T(k,j,1:s,1:s);
                inverse = C \eye(s);
                if (sum(sum(inverse==Inf | inverse==-Inf))>=1) % if T(k,j,1:s,1:s) is not invertible
                    [~,E] = eig(C);
                    E = max(abs(E),[],"all")+1; % E is non-zero and not an eigenvalue of C.
                    T(k,j,:,:) = zeros(max_SDR);

                    % 1
                    T(k, j, 1:s, s_new+1:s_new+s) = eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = eye(s);

                    Z = zeros(max_SDR);
                    D = -(C-E*eye(s))\eye(s);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = D;
                    A0 = A0 + Z;
                    deter = deter*det(D);

                    s_new = s_new + s;

                    % 2
                    T(k, j, 1:s, s_new+1:s_new+s) = eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = -eye(s)/E;
                    A0 = A0 + Z;
                    deter = deter*((-1)^s)/(E^s);

                    s_new = s_new + s;

                else
                    T(k,j,:,:) = zeros(max_SDR);
                    % 1
                    T(k, j, 1:s, s_new+1:s_new+s) = eye(s);
                    T(k, j, s_new+1:s_new+s, 1:s) = eye(s);

                    Z = zeros(max_SDR);
                    Z(s_new+1:s_new+s, s_new+1:s_new+s) = -inverse;
                    A0 = A0 + Z;
                    deter = deter*det(-inverse);

                    s_new = s_new + s;
                end
            end
        end
    end
    if (deter>0 || mod(s_new,2)==1)
        T = deter^(-1/s_new)*T;
        A0 = deter^(-1/s_new)*A0;
        s = s_new;
    else
        A0(s_new+1,s_new+1) = 1/deter;
        s = s_new+1;
    end
end
Pencil = zeros(r+1,s,s);
for k = 1:r
    for j = 1:M
        Aux = zeros(1,s,s);
        Aux(1,:,:) = T(k,j,1:s,1:s);
        Pencil(k,:,:) = Pencil(k,:,:)+Aux(1,:,:);
    end
end
Pencil(r+1,:,:) = A0(1:s,1:s);
end