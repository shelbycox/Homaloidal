--=================================================
-- (0) Set-up
--=================================================
a = 2;
b = 3;
n = a + b + 1;
d = a + b;
R = QQ[x_0..x_(a+b+1)];
X = transpose vars R;

--=================================================
-- (1) construct S(a, b)
--=================================================

-- (1a) function to compute Sab matrix ============
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

-- (1b) compute the Sab ideal =====================
I = S(a, b, R_*);
netList I_* -- see the generators of I

-- checks
degree I == a + b -- check that the degree is a + b
dim I == 3 -- check that the projective dimension is 2 (so affine dimension is 3)
codim I == 4 -- projective dimension 2 in a PP6, so codimension is 4.

--=================================================
-- (2) compute E (and its orthogonal complement)
--=================================================

-- (2a) compute S(a,b)* ===========================
loadPackage "Resultants"
Istar = dualVariety I
codim Istar == 1 -- the dual is a hypersurface (defined by one equation)
degree Istar == a + b -- the dual has degree a + b

-- (2b) compute the singular locus of deg >= b ====
JacStar = jacobian Istar
Jstar = ideal jacobian Istar -- This variety is not linear! (see the decomposition), but it doesn't look linear.
JJstar = ideal jacobian Jstar -- singularities of degree at least 3 (= b in this example)
JJJstar = ideal jacobian JJstar
decompose JJJstar
netList Jstar_*

Eperp = ideal(x_0..x_2) -- This is a linear space with sigularities containing some singularities of degree 5.
codim Eperp == b -- codim b --> proj dimension E = b - 1 -- not sure if this holds in general

-- (2c) compute the orthogonal complement of Eperp=
EperpMat = matrix {{1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}} -- TODO: generalize this.
Eperp == ideal (EperpMat*X)
E = ideal ((transpose gens ker EperpMat)*X)

--==================================================
-- (3) project S(a, b) to get X(a, b)

-- NOTE: What we will actually compute is X(a, b)*
-- as the intersection S(a, b)* \cap Psi, for some
-- Psi \supset E^Perp a linear space.
--==================================================

-- (3a) compute a suitable Psi ====================
EminusSab = saturate(E, I) -- see what generators we need in Psi
-- to find a dim a - 2 = 0 subspace (i.e., a point) contained in E
PsiMat = matrix {{1, 0, 0, 0, 0, 0, 0}, -- Psi generators in matrix form
				{0, 0, 1, 0, 0, 0, 0}, 
				{0, 0, 0, 1, 0, 0, 0}, 
				{0, 0, 0, 0, 1, 0, 0}, 
				{0, 0, 0, 0, 0, 1, 0}, 
				{0, 0, 0, 0, 0, 0, 1}}
Psi = ideal (PsiMat*X) -- define the ideal

-- (3b) check that this Psi is good ===============
isSubset(E, Psi) -- check that V(E) contains V(Psi)
codim Psi == 6 -- check that Psi is a point
saturate(Psi + I, ideal (x_0..x_6)) == ideal 1_R -- test that Psi is not a point on S(a,b)

-- (3c) compute the orthogonal complement of Psi ==
Psiperp = ideal ((transpose gens ker PsiMat)*X) -- get the orthogonal complement
dim Psiperp == n

-- (3d) compute X(a, b)* ==========================
XabStar = Istar + Psiperp -- get the ideal of XabStar
codim XabStar == 2 -- X(a,b) is a hypersurface in P^(a + b)
-- Note: here we have the extra variable x_1, which should be eliminated
degree XabStar == a + b

--==================================================
-- (4) computing X(a,b) = (X(a,b)*)*
--==================================================

Xab = dualVariety XabStar -- WARNING: times out for some Psi!
decompose Xab
codim Xab
dim Xab == 4 -- 
degree Xab == a + b -- check the degree

T = QQ[x_0,x_2..x_6]
RtoT = map(T, R, {x_0, 0, x_2..x_6})
TXab = RtoT(Xab)

--==================================================
-- (5) find the line directrix Lambda
--==================================================

ESab = E + I
mingens ESab
ESabStar = dualVariety ESab
LambdaStar = ESabStar + Psiperp
Lambda = dualVariety LambdaStar
dim Lambda -- Lambda should be a line
degree Lambda
E
eliminate(x_1, E)

Lambda = RtoT(E)

--==================================================
-- (6) find b - a rulings on Xab
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
TF1 = RtoT(F1)

--==================================================
-- (7) compute Phi (spanned by line directrix + rulings)
--==================================================

-- Phi has dim b - a - 1
Phi = Lambda + TF1 -- not sure about this...
Phimat = matrix {{1, 0, 1, 0, 0, 0, 0},
				 {1, 0, 0, 0, 0, 0, -1},
				 {1, 0, 0, 0, 0, -1, 0},
				 {1, 0, 0, 0, -1, 0, 0},
				 {1, 0, 0, -1, 0, 0, 0}}
inverse Phimat
Phi = RtoT(ideal(Phimat*X))
Phiperp = RtoT(ideal ((transpose gens ker Phimat)*X))

--==================================================
-- (8a) compute Y(a, b) with elimination
--==================================================

newT = QQ[x_0,x_2..x_6,z_0..z_4]
zvars = ((entries vars newT)_0)_{6..10}
Phi = sub(Phi, newT)
elimMat = gens Phi || (vars newT)_{6..10}
Xab = sub(Xab, newT)
projXab = Xab + minors(2, elimMat)
satXab = saturate(projXab, product Phi_*)
use newT
x_0
Yab = eliminate(satXab, {x_0,x_2,x_3,x_4,x_5,x_6})
isPrime Yab == true

--==================================================
-- (8b) compute Y(a, b) with a linear change of coordinates
--==================================================
inverse Phimat
YabStar = RtoT(XabStar) + Phiperp
YabStar = sub(YabStar, T)
U = QQ[y_0..y_4]
coordChangeMat = matrix {{1, 1, 1, 1, 1, 1},
						 {0, -1, -1, -1, -1, -1},
						 {1, 1, 0, 1, 1, 1},
						 {1, 1, 1, 0, 1, 1},
						 {1, 1, 1, 1, 0, 1},
						 {1, 1, 1, 1, 1, 0}}
coordChange = map(T, T, entries (coordChangeMat*(transpose vars T))_0)
coordPhi = coordChange(Phi)
coordTXab = coordChange(TXab)
Yab = eliminate(x_2, coordTXab)
degree Yab == a + b

--==================================================
-- (9) compute Y(a, b)*
--==================================================

zT = QQ[z_0..z_4]
Yab = sub(Yab, zT) -- don't need UYabstar after doing this
Yabstar = dualVariety Yab
codim Yabstar == 1 -- Yabstar should be a hypersurface
degree Yabstar == a + b

U = QQ[x_0, x_3..x_6]
TtoU = map(U, T, {x_0, 0, x_3..x_6})
UYabstar = TtoU(Yabstar)

codim UYabstar == 1
degree UYabstar == 5
jacobian (UYabstar)_1
#(terms Yabstar_0)
sum for i in (flatten entries (coefficients Yabstar_0)_1) list abs sub(i, QQ)
texMath Yabstar_0

--==================================================
-- (10) check to see if Yabstar is homaloidal
--==================================================
loadPackage "RationalMaps"
polarYab = map(zT, zT, entries (jacobian Yabstar_0)_0) -- find the polar map
isBirationalMap(polarYab) -- if true, then Yabstar is homaloidal!
inverseOfMap(polarYab); -- computes an explicit inverse to polarYab

--==================================================
-- BONUS: compute the Newton polytope of Yabstar
--==================================================
loadPackage("Polyhedra")
N = newtonPolytope(Yabstar_0)
fVector(N)
vertices N
volume N
