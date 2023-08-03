SMatrix = method();
SMatrix = (a, b, x) -> (
   n := a + b;
   matrix for i from 1 to 2 list (
    	for j from 1 to n list (
	    if i == 1 then (
	    	if j <= a then x_(j-1) else x_j
	    	)
	    else if j <= a then x_j else x_(j+1)
    	    )
    	)
    )

a = 1
b = 3
R = QQ[x_0..x_(a+b+1)]
X = vars R

M = SMatrix(a,b,R_*)
I = minors(2, M)
netList I_*

degree I -- check that the degree is a + b
dim I -- check that the projective dimension is 2 (so affine dimension is 3)

JacI = jacobian I -- since S(a,b) is smooth we expect the Jacobian to have rank #vars - dim = (a + b + 2) - 2 = a + b, except possibly at the origin.

J = saturate(I, ideal (x_0..x_(a+b)))
I == J -- saturation doesn't change the ideal in this case

-- with random subspaces
MM = random(QQ^(a + b - 1), QQ^(a + b + 2))
MM = sub(MM, R)
KMat = MM*transpose(X)
K = ideal KMat
dim K
L = saturate(I + K, ideal(R_*))

MMker = ker MM
MMkerMat = transpose matrix {{576, -7418, -29256}, {-4164, 7230, -28965}, {-337, -5682, 11406}, {3244, 0, 0}, {0, 20275, 0}, {0, 0, 20275}} -- need to do this step systematically
MMkerMat = sub(MMkerMat, R)
Kperp = ideal (MMkerMat*transpose(X))
netList Kperp_*
dim Kperp

loadPackage "Resultants"
Sstar = dualVariety I
Xideal = saturate(Sstar + Kperp, ideal(R_*))
dim Xideal
degree Xideal

-- coming up with a linear space that works
K = ideal (x_0 - x_1, x_5 + x_6, x_0 + x_2, x_3 + x_5)
L = saturate(I + K, ideal(x_0..x_6))

dim(K)
