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
a = 2
b = 3
R = QQ[x_0..x_(a+b+1)]
'''

We will use the ideal of minors definition to produce \( S(2,3) \).

'''
M = getSabMat(2,3)
'''

## Dual Varieties

'''
needsPackage "Resultants"
--viewHelp "Resultants"
'''

## Projections

```
code block
```

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
