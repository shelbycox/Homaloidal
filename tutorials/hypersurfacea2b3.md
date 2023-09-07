# Homaloidal Hypersurfaces from Scroll Surfaces (Example)
## Scroll Surfaces

```
restart
needsPackage "Homaloidal"
--viewHelp "Homaloidal"
```

### What is a scroll surface?

\( math mode \)

### Let's define a scroll surface

From now on, let \( a = 2, b = 3 \). 

'''
a = 2;
b = 3;
R = QQ[x_0..x_(a+b+1)];
X = transpose vars R;
'''

We will use the ideal of minors definition to produce \( S(2,3) \).

'''
M = SMatrix(a, b, R_*)
Sab = minors(2, M);
netList Sab
'''

'''
degree Sab == a + b -- check that degree is a + b.
codim Sab == 4 -- so dim(Sab) = 6 - 4 = 2.
'''

## Dual Varieties

'''
needsPackage "Resultants"
--viewHelp "Resultants"
'''

'''
SabStar = dualVariety Sab;
codim SabStar == 1 -- the dual is a hypersurface
degree SabStar == a + b -- the dual has degree a + b
'''

The hypersurface \(S(a,b)^\ast\) has a special property: the subvariety of singularities of degree at least \(b\) is a linear space, \(E^\perp\). The singularities of degree 2 and higher are described by the jacobian ideal. The singularities of degree 3 and higher are described by the jacobian of the jacobian ideal, and so on.

'''
Jstar = ideal jacobian SabStar -- singularities of degree at least 2.
JJstar = ideal jacobian Jstar -- singularities of degree at least 3.
decompose JJstar
'''

'''
EperpMat = matrix {{1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}}
Eperp = ideal (EperpMat*X)
E = ideal ((transpose gens ker EperpMat)*X)
'''

## Projections

A projection from a linear space \( \Pi \) is ...

There are two ways to compute this projection in Macaulay2.

### Method 1: elimination

If \( \Pi = \langle x_0, \ldots, x_k \), \( k < n \), the projection is a coordinate projection. On the ring side, this corresponds to eliminating the variables \( x_{k + 1}, \ldots, x_n \) from the ideal.

```
code block
```

### Method 2: dual intersection dual

In many cases, the projection is the intersection of \( X^\ast \) with \( \Pi^\perp \), the orthogonal complement of \( \Pi \).

## \( \Phi \), the first linear space

The first projection we make from \( S(a,b) \) is from \( \Phi \), which is any linear space satisfying the following:
* \( \dim \Phi = a - 2 \). In our case, \( a = 2 \), so \( \Phi \) is a point.
* \( \Phi \) \subseteq \langle E \rangle.

### \( E^\perp \)

### \( E \)

### \( \Phi \)
    
## \(X(a,b)\)

### subblock
* unordered list

## \( \Psi \), the second linear space

## \(Y(a,b)\)

1. numbered list
