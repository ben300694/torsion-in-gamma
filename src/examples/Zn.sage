"""
Cyclic group of order n
G = Z/(n) = < a | a^n > 
"""

def cyclic_group(n):
    F1.<a> = FreeGroup(1)
    relations_Zn = [a^n]
    Zn = F1.quotient(relations_Zn)
    return Zn

Z2 = cyclic_group(2)
Z3 = cyclic_group(3)
Z4 = cyclic_group(4)

# Shotcut to get the group algebra over 
# cyclic groups represented as subgroup of permutation group
def Z_cyclic_group(n):
    Zn = cyclic_group(n)    
    return GroupAlgebra(Zn.as_permutation_group(), ZZ)

Z_Z2 = Z_cyclic_group(2)
Z_Z3 = Z_cyclic_group(3)
Z_Z4 = Z_cyclic_group(4)
