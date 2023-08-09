--=================================================
-- construct S(a, b)
--=================================================
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

S = method();
S = (a, b, x) -> (
    M = SMatrix(a, b, x);
    minors(2, M)
    )

a = 1;
b = 3;
R = QQ[x_0..x_(a+b+1)];
X = transpose vars R;

I = S(a, b, R_*);
netList I_* -- see the generators of I

degree I == a + b -- check that the degree is a + b
dim I == 3 -- check that the projective dimension is 2 (so affine dimension is 3)

JacI = jacobian I -- since S(a,b) is smooth we expect the Jacobian to have rank #vars - dim = (a + b + 2) - 2 = a + b, except possibly at the origin.

I == saturate(I, ideal (x_0..x_(a+b))) -- saturation doesn't change the ideal in this case
netList I_*

--=================================================
-- compute E (and its orthogonal complement)
--=================================================

-- compute S(a,b)*
loadPackage "Resultants"
Istar = dualVariety I
dim Istar
degree Istar

-- compute the singular locus (up to degree b)
JacStar = jacobian Istar
Jstar = ideal jacobian Istar -- This variety really is linear (see the decomposition), but it doesn't look linear.
D = decompose Jstar -- But the primary decomposition shows that it is a linear space of dimension 5 - 2 = 3 = b.

-- According to the primary decomposition, the singularities of the dual lie in this linear space.
Eperp = D_0

-- We only need Eperp for the computation of XabStar, so I don't think we need the code below.
EperpMat = matrix {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}} -- TODO: generalize this.
Eperp == ideal (EperpMat*X)
E = ideal ((transpose gens ker EperpMat)*X)

--==================================================
-- project S(a, b) to get X(a, b)

-- NOTE: What we will actually compute is X(a, b)*
-- as the intersection S(a, b)* \cap Psi, for some
-- Psi \supset E^Perp a linear space.
--==================================================

-- to find a subspace that contains Eperp, we just need to remove some generators from Eperp.
Psi = ideal x_0

XabStar = I + Psi
Xab = dualVariety XabStar

dim Xab
degree Xab

--==================================================
-- find the line directrix Lambda
--==================================================

Lambda = E + Psi -- I am confused. 
dim Lambda
degree Lambda

--==================================================
-- find rulings on X
--==================================================



--==================================================
-- compute Phi
--==================================================



--==================================================
-- compute Y(a, b)*
--==================================================



--==================================================
--- row spans and their orthogonal complements
--==================================================
T = QQ[x,y,z]
M = random(QQ^2, QQ^3)
L = ideal (M*(transpose vars T)) -- linear space is the kernel of M i.e. orthogonal complement of the row span, so the row span is the orthogonal complement of L
rowSpan = ideal ((transpose gens ker M)*(transpose vars T)) -- this is the orthogonal complement of L, the row span of M

dim L == 1 -- orthogonal complement of row span has dim: n - row rank = 1
dim rowSpan == 2 -- row span has dim: row rank = 2
saturate(L + rowSpan, ideal T_*) == ideal 1_T -- the row span and its orthogonal complement intersect only at the origin

P = ideal ((transpose vars T) - (matrix (transpose M)_1)) -- a point in the row span
-- rowSpan contains P
degree localize(rowSpan, P)
dim (rowSpan + P)
degree (rowSpan + P)

Q = ideal ((vars T) - (matrix {{1,1,1}})) -- a point not in the row span
-- rowSpan does not contain Q
degree localize(rowSpan, Q)
dim (rowSpan + Q)
degree (rowSpan + Q)


--=================================================
--- column spans and their orthogonal complement
--=================================================
R = QQ[u,v,w]
N = random(QQ^3, QQ^2)
colSpan = ideal ((transpose gens ker transpose N)*(transpose vars R))
O = ideal ((transpose N)*(transpose vars R))

dim(colSpan) == 2 -- column span had dim: column rank = 2
dim(O) == 1 -- orthogonal complement has dim: kernel rank = 1
saturate(colSpan + O, ideal R_*) == ideal 1_R -- colSpan and its orthogonal complement intersect only at the origin

Q = ideal ((transpose vars R) - (matrix N_1)) -- get a point in the column span
-- colSpan contains Q
degree localize(colSpan, Q)
dim (colSpan + Q)
degree (colSpan + Q)

P = ideal ((vars R) - (matrix {{1,1,1}})) -- a point not in the column span
-- colSpan does not contain P
degree localize(colSpan, P)
dim (colSpan + P)
degree (colSpan + P)
