--=====================================
-- Definitions
--=====================================
R = QQ[x_0..x_2] -- working in the projective plane
f = x_0^2 + x_1*x_2 -- here is a conic
C = ideal f

--=====================================
-- Sanity checks
--=====================================
dim C == 2 -- C has affine dimension 2 = projective dimension 1
degree C == 2 -- C has degree 2, so it really is a conic

--=====================================
-- Jacobian of C tells us...
-- 1) The conic is smooth because the tangent space has dimension 1 except at the origin.
-- 2) f is homaloidal because the polar map is birational.
--=====================================
jacobian C

loadPackage "DeterminantalRepresentations"
detRep f -- the package doesn't compute a determinantal representation for f :(

M = matrix {{x_0, 0}, {0, x_0}}

--=======================================================
-- Using the process in Stefan and Welter's paper (2021)
--=======================================================
-- STEP 1: Find an affine linear seed.
S = QQ[x_0..x_2,y,z]
f0 = y + z -- start with an affine polynomial in one variable

M0 = detRep f0 -- the package still doesn't compute the determinantal representation, not sure why

A1 = matrix {{1}}
A2 = matrix {{1}}
-- NOTE: A_1 = A_2 = [1] are not linearly independent.

M0 = y*A1 + z*A2
det M0 == f0 -- this matrix works

-- STEP 2: y -> x_0^2, our first simple product substitution.
f1 = x_0^2 + z
f1 = det matrix {{ x_0^2 + z }} -- This is not a symmetric determinantal representation.

B0 = matrix {{ 0, 0}, {0, -1}}
B1 = matrix {{ 0, 1 }, { 1, 0 }}
B2 = matrix {{ 1, 0 }, {0, 0}}
M1 = (B0 + x_0*B1 + z*B2) ++ matrix {{ -1 }}
det M1 == f1 -- this matrix works

-- STEP 3: z -> x_1*x_2
f2 = sub(f, S)
B0' = matrix {{ 0, 0, 0}, {0, -1, 0}, {0, 0, -1}}
B1' = matrix {{ 0, 1, 0}, {1, 0, 0}, {0, 0, 0}}
B2' = matrix {{ 1, 0, 0}, {0, 0, 0}, {0, 0, 0}}
det (B0' + x_0*B1' + x_1*x_2*B2') == f2 -- This is not a symmetric determinantal representation.

C = matrix {{ 0, (x_1 + x_2)/2, (x_1 - x_2)/2}, {0, -1/2, 0}, {0, 0, 1/2}}
C = C + transpose C

block2 = matrix {{ 0, 0}, {0, 0}}
M11 = (B0' + x_0*B1') ++ block2
M2 = M11 + C
