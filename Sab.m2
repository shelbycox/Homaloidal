SMat = method();
SMat = (a, b, x) -> (
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
b = 2
n = a + b
R = QQ[x_0..x_(n+1)]

M = SMat(a,b,R_*)
I = minors(2, M)
netList I_*

degree I -- check that the degree is a + b
dim I -- check that the projective dimension is 2 (so affine dimension is 3)

JacI = jacobian I -- since S(a,b) is smooth we expect the Jacobian to have rank #vars - dim = (n + 2) - 2 = n, except possibly at the origin. This seems hard to check...

J = saturate(I, ideal (x_0..x_n))
I == J -- saturation doesn't change the ideal in this case

