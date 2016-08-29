"""
A module with small helper functions
"""

def matrix_map(f,AA):
    """
    The usual map function implemented for the matrix type
    
    INPUT:
        :AA: A matrix with entries of type t   

        :f: a function from type t to type t'

    OUTPUT:
        A matrix of with entries of type t',
        where the function f was applied to any entry in AA
    """
    apply_to_rows = lambda x: map(f, x)    
    list_of_new_rows = map(apply_to_rows, AA.rows())
    return matrix(AA.base_ring(), list_of_new_rows)

def matrix_to_list(AA):
    """
    INPUT:
        :AA: A matrix

    OUTPUT:
        The matrix in form of a list, where the rows are
        consecutively inserted
    """
    return [entry for (i, row) in enumerate(AA) for entry in row]    

def diagonal(AA):
    """
    INPUT:
        :AA: matrix    
    
    OUTPUT:
        The main diagonal of AA in list form
    """
    n = min(AA.nrows(), AA.ncols())
    return [AA[i,i] for i in range(0,n)]

def contains_element_not_zero_or_one(l):
    """
    INPUT:
        :l: list of integers

    OUTPUT:
        True if l contains an entry not equal to zero or one,
        otherwise False
    """
    for x in l:
        if x != 0 and x != 1:
            return True
    return False

def smith_form_with_gap(AA):
    """
    Wrapper to use the algorithm in GAP for the smith normal form of a matrix
    """
    gAA = libgap(AA)
    gAA_smith = gAA.SmithNormalFormIntegerMat()
    return matrix(gAA_smith)

def count(i, lis):
    """
    Count how often 'i' occurs in the list 'lis'
    """
    occurences = 0
    for x in lis:
        if i == x:
            occurences=occurences+1
    return occurences

# lambda function to generate the commutator of two group elements
commutator = lambda x, y: x*y*x^(-1)*y^(-1)

# For testing the functions on the group algebra
# S3 is the symmetric group
# S3 = SymmetricGroup(3)
# ZS3 is the Integer Group Algebra of the group S3
# Z_S3 = GroupAlgebra(S3, ZZ)
# s1, s2 = Z_S3.gens()

