"""
Module for searching torsion in Gamma(pi_2 K)/(pi_1 K) 
for K a presentation of a finite finitely presented group
"""

cyclics = lambda x : [cyclic_group(n) for n in range(2,x)]
dihedrals = lambda x : [dihedral_group(n) for n in range(3,x)]

list_of_groups_example = cyclics(12) + dihedrals(5) + [product_of_cyclic_groups(2,2)]


def test_list_of_groups(list_of_groups):
    """
    Checks Gamma(pi_2 G)/(pi_1 G) for torsion for the groups
    in the list 'list_of_groups'.
    Prints out status messages in between
    """    
    for G in list_of_groups:
        print G
        print "Checking for torsion in Gamma(pi_2 G)/(pi_1 G)"
        print "Rank of pi_2(G): %d" % pi_2(G).rank()        
        factors = struc(G)
        print "Invariant factors:"        
        print factors
        result = contains_element_not_zero_or_one(factors)
        if result == False:
            print "Gamma(pi_2(G)/(pi_1 G): No torsion found"
            print "Gamma(pi_2(G)/(pi_1 G)) isomorphic to Z^%d" % count(0, factors)
        else:
            print "Gamma(pi_2(G)/(pi_1 G): Has torsion"
        print "-----------"

# Some specializations of the function 'test_list_of_groups(list_of_groups)'

# Presentations of cyclic groups of the form
# G = < a | a^n >
#   = Z/(n)

def test_cyclic_groups(limit=12):
    """
    Check Gamma(pi_2 G)/(pi_1 G) for torsion for the groups G = Z/(2), Z/(3), ..., Z/(limit)
    Prints out status messages in between
    """
    test_list_of_groups(cyclics(limit+1))

# Presentations of dihedral groups of the form
# G = < s, d | s^2, d^n, (s*d)^2 >
#   = D_n

def test_dihedral_groups(limit=12):
    """
    Check Gamma(pi_2 G)/(pi_1 G) for torsion for the groups G = D_2, D_3, ..., D_(limit)
    Prints out status messages in between
    """
    test_list_of_groups(dihedrals(limit+1))

def run_test():
    test_list_of_groups(list_of_groups_example)
