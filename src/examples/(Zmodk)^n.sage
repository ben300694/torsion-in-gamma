"""
The n-fold product of the cyclic group with k elements
G = ( Z/(k) ) ^ n = Z/(k) x Z/(k) x ... x Z/(k) (n times)
  = < t0, t1, ..., tn | t0^k = t1^k = ... = tn^k = e, 
                        [t0, t1], [t0, t2], ..., [t0, tn],
                        [t1, t2], ..., [t1,tn],
                        ...,
                        [t(n-1), tn] >
where [a, b] := aba^(-1)b^(-1) is the commutator of a and b 
"""

def product_of_cyclic_groups(k, n):
    Fn = FreeGroup(n)
    words_order = [Fn([i+1])^k for i in range(0,n)]
    words_commutators = []
    for i in range(0, n-1):
        for j in range(i+1, n):
            words_commutators.append(commutator(Fn([i+1]), Fn([j+1])))
    words = words_order + words_commutators
    Zk_n = Fn.quotient(words)
    return Zk_n

# Some definition of special cases

#----
# Product of the cyclic group of order 2 with itself
# H = Z/(2) x Z/(2) = < a, b | a^2, b^2, aba^{-1}b^{-1}>
#----
Z2xZ2 = product_of_cyclic_groups(2,2)
Z_Z2xZ2 = GroupAlgebra(Z2xZ2, ZZ)
