# and example document to calculate the homology of the three-torus
# T^3 := S^1 x S^1 x S^1
# for this example see also [fenn1983techniques]

# implements the definition [x,y] := xyx^-1y^-1
commutator = lambda x, y: x*y*x^-1*y^-1

G.<a,b,c> = FreeGroup()
r = commutator(a,b)
s = commutator(b,c)
t = commutator(c,a)

# K2 is the two-skeleton of the 3-torus, here we are using the correspondence
# {finite 2-dim CW-complexes} <--> {group presentations}
K2 = G.quotient([r,s,t])

# from the documentation at [http://doc.sagemath.org/html/en/reference/groups/sage/groups/finitely_presented.html]:
# This matrix is given by the fox derivatives of the relations with respect to the generators.
print "Alexander Matrix of K2 = "
print K2.alexander_matrix()

# constuct the group ring of K2 over the Integers
ZK2 = GroupAlgebra(K2, ZZ)
