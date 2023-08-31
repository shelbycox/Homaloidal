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
n = a + b + 1;
d = a + b;
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

EminusSab = saturate(E, I) -- see what generators we need in Psi

-- to find a dim a - 2 = 0 subspace (i.e., a point) contained in E
PsiMat = matrix {{1, 0, 0, 0, 0, 0, 0}, -- Psi generators in matrix form
				{0, 0, 1, 0, 0, 0, 0}, 
				{0, 0, 0, 1, 0, 0, 0}, 
				{0, 0, 0, 0, 1, 0, 0}, 
				{0, 0, 0, 0, 0, 1, 0}, 
				{0, 0, 0, 0, 0, 0, 1}}
Psi = ideal (PsiMat*X) -- define the ideal

-- checks
isSubset(E, Psi) -- check that V(E) contains V(Psi)
codim Psi == 6 -- check that Psi is a point
saturate(Psi + I, ideal (x_0..x_6)) == ideal 1_R -- test that Psi is not a point on S(a,b)

Psiperp = ideal ((transpose gens ker PsiMat)*X) -- get the orthogonal complement
dim Psiperp == n

XabStar = Istar + Psiperp -- get the ideal of XabStar
codim XabStar

--==================================================
-- computing X(a,b) = (X(a,b)*)*
--==================================================

Xab = dualVariety XabStar -- WARNING: times out for some Psi!
decompose Xab
codim Xab -- I think Xab should be a surface, but it has dimension 3.
dim Xab -- 
degree Xab == a + b -- check the degree

--==================================================
-- computing the dual by hand --
--==================================================

T = QQ[y_0..y_(a+b+1)][x_0..x_(a+b+1)];
IT = sub(XabStar, T);
N = transpose jacobian IT

S = QQ[y_0..y_(a+b+1),x_0..x_(a+b+1)];
N = sub(N, S);
NN = (vars S)_{0..(a+b+1)} || N
IT = sub(IT, S);

J = IT + minors(3, NN);
singJ = ideal singularLocus J
J = saturate(J, singJ)
Xab = eliminate(J, (entries (vars S)_{(a+b+2)..(2*a+2*b+3)})_0)

codim Xab

--==================================================
-- find the line directrix Lambda
--==================================================

ESab = E + I
mingens ESab
ESabStar = dualVariety ESab
LambdaStar = ESabStar + Psiperp
Lambda = dualVariety LambdaStar
dim Lambda -- Lambda should be a line
degree Lambda
E

--==================================================
-- find b - a rulings on Xab
--==================================================

-- first find some lines on S(a,b)
L1mat = matrix {{1, -1, 0, 0, 0, 0, 0},
				{1, 0, -1, 0, 0, 0, 0},
				{0, 0, 0, 1, -1, 0, 0},
				{0, 0, 0, 1, 0, -1, 0},
				{0, 0, 0, 1, 0, 0, -1}
}
L1 = ideal (L1mat*X)
codim L1 == n - 1 -- double check that L1 is a line
-- b - a = 3 - 2 = 1, so we just need one ruling

-- now project those lines to X(a,b)
L1perp = ideal ((transpose gens ker L1mat)*X)
F1perp = L1perp + Psiperp
F1perpmat = matrix {{1, 1, 1, 0, 0, 0, 0},
					{0, 0, 0, 1, 1, 1, 1},
					{0, 1, 0, 0, 0, 0, 0}}
F1perp == ideal (F1perpmat*X)
F1 = ideal ((transpose gens ker F1perpmat)*X)

--==================================================
-- compute Phi (spanned by line directrix + rulings)
--==================================================

Phiperp = Lambda + F1 -- not sure about this...
mingens Phiperp

--==================================================
-- compute Y(a, b)*
--==================================================

Yabstar = XabStar + Phiperp
decompose Yabstar
