"""
Dihedral group of order 2*n, symmetry group of the regular n-gon
D_n = < s, d | s^2, d^n, (s*d)^2 >
"""

def dihedral_group(n):
    F2.<s, d> = FreeGroup(2)
    relations_Dn = [s^2, d^n, (s*d)^2]
    Dn = F2.quotient(relations_Dn)
    return Dn

D3 = dihedral_group(3)
D4 = dihedral_group(4)

D8 = dihedral_group(8)
