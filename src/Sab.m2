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

a = 2;
b = 3;
R = QQ[x_0..x_(a+b+1)];
X = transpose vars R;

I = S(a, b, R_*);
netList I_* -- see the generators of I

degree I == a + b -- check that the degree is a + b
dim I == 3 -- check that the projective dimension is 2 (so affine dimension is 3)
codim I == 4 -- projective dimension 2 in a PP6, so codimension is 4.

JacI = jacobian I -- since S(a,b) is smooth we expect the Jacobian to have rank #vars - dim = (a + b + 2) - 2 = a + b, except possibly at the origin.

I == saturate(I, ideal (x_0..x_(a+b))) -- saturation doesn't change the ideal in this case
netList I_*

--=================================================
-- compute E (and its orthogonal complement)
--=================================================

-- compute S(a,b)*
loadPackage "Resultants"
Istar = dualVariety I
codim Istar == 1 -- the dual is a hypersurface
degree Istar == a + b -- the dual has degree a + b

-- compute the singular locus (up to degree b)
JacStar = jacobian Istar
Jstar = ideal jacobian Istar -- This variety is not linear! (see the decomposition), but it doesn't look linear.
D = decompose Jstar -- Jstar is actually a prime ideal of degree 6.
netList (D_0)_*

P = ideal (x_0..x_2, x_3 + random(QQ), x_4 - random(QQ), x_5 - random(QQ), x_6 - random(QQ))
IstarP = localize (Istar, P)
degree IstarP -- this point is a singularity of degree 5 in Istar

JJstar = ideal jacobian Jstar -- singularities of degree at least 3 (= b in this example)
decompose JJstar

-- According to the primary decomposition, the singularities of the dual lie in this linear space.
Eperp = ideal(x_0..x_2) -- This is a linear space with sigularities containing some singularities of degree 5.
codim Eperp == b -- codim b --> proj dimension E = b - 1 -- not sure if this holds in general

-- We only need Eperp for the computation of XabStar, so I don't think we need the code below.
EperpMat = matrix {{1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}} -- TODO: generalize this.
Eperp == ideal (EperpMat*X)
E = ideal ((transpose gens ker EperpMat)*X)
codim E -- so E is 2 dimensional

--==================================================
-- project S(a, b) to get X(a, b)

-- NOTE: What we will actually compute is X(a, b)*
-- as the intersection S(a, b)* \cap Psi, for some
-- Psi \supset E^Perp a linear space.
--==================================================

-- to find a dim a - 2 = 0 subspace contained in E
PsiMat = matrix {(entries random(QQ^1, QQ^7))_0, -- Psi generators in matrix form
				(entries random(QQ^1, QQ^7))_0, 
				{0, 0, 0, 1, 0, 0, 0}, 
				{0, 0, 0, 0, 1, 0, 0}, 
				{0, 0, 0, 0, 0, 1, 0}, 
				{0, 0, 0, 0, 0, 0, 1}}
Psi = ideal (PsiMat*X) -- double check the matrix form
isSubset(E, Psi) -- check that V(Psi) contains V(E)
Psiperp = ideal ((transpose gens ker PsiMat)*X) -- get the orthogonal complement


EminusSab = saturate(E, I) -- see what generators we need in Psi
saturate(Psi + I, ideal (x_0..x_6)) == ideal 1_R -- test that Psi is not a point on S(a,b)


XabStar = Istar + Psiperp -- get the ideal of XabStar
codim XabStar

Xab = dualVariety XabStar -- This computation now takes a long time...
decompose Xab
codim Xab == 4 -- I think Xab should be a surface, but it has dimension 3.
dim Xab -- 
degree Xab -- 

--==================================================
-- try computing the dual by hand -- IT WORKS!
--==================================================

T = QQ[y_0..y_(a+b+1)][x_0..x_(a+b+1)];
IT = sub(I, T);
N = transpose jacobian IT

S = QQ[y_0..y_(a+b+1),x_0..x_(a+b+1)];
N = sub(N, S);
NN = (vars S)_{0..(a+b+1)} || N
IT = sub(IT, S);

J = IT + minors(5, NN);
J = saturate(J, entries (vars S)_{(a+b+2)..(2*a+2*b+3)})
Xab = eliminate(J, (entries (vars S)_{(a+b+2)..(2*a+2*b+3)})_0)

codim Xab

--==================================================
-- find the line directrix Lambda
--==================================================

LambdaStar = Eperp + Psiperp
Lambda = dualVariety LambdaStar
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
