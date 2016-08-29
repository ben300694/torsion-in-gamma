"""
Module with all the functions needed to calculate the second homotopy
group of a cell complex associated to a finitely presented finite group.

Author: Benjamin Ruppik (University of Bonn)
Date: June 2016

Credit:
This uses the packages
- sage.groups.free group written by Miguel Angel Marco Buzunariz and Volker Braun
- sage.groups.finitely presented written by Miguel Angel Marco Buzunariz
- sage.algebras.group algebra.GroupAlgebra written by David Loeffler, Martin Raum, John Palmieri.
Some of the functions use algorithms which are implemented in GAP.
"""

attach("util.sage")

def second_boundary_matrix(grouppresentation, im_gens=None):
    """
    Calculates the linear second boundary map
    in the cellular chain complex belonging to the group presentation.

    The matrix will act by multiplication from the right, more precisely:
    Let F be the matrix returned by this function, v in C_2 (the second
    chain group).
    Then the image is calculated as v*F.

    Remark: There is a library function on finitely presented groups
    called '.alexander_matrix()' which also returns the the matrix calculated here.
    We could use this function, but it doesn't hurt to implement this ourselves.

    INPUT:
        :grouppresentation: a finitely presented group of type
        'sage.groups.finitely_presented.FinitelyPresentedGroup_with_category'
        
        :im_gens: Here one can give images of the generators of the free group

    OUTPUT: 
        A matrix with entries the Fox-derivatives of the relations
        Specifically this matrix has in the (i,j)-entry the value
        del_{j}(rel[i])
    """
    rel = grouppresentation.relations()
    gen = grouppresentation._free_group.gens()
    # The matrix contructor takes the arguments 'ring', 'nrows', 'ncolumns', 'generating function'
    return matrix(len(rel), len(gen), lambda i,j: fox(rel[i], gen[j], im_gens))

def groupalgebra_element_to_matrix(g, side):
    """
    Assumption: G is a finite group.
    As R-modules there is an isomorphism between 
    the group ring R[G] and R^ord(G). 
    Multiplication by an element of the group algebra gives a R-linear
    map R[G] --> R[G] which can be represented by a
    ord(G) x ord(G) matrix with entries in R.
    This function takes an element in the groupalgebra R[G] and returns this
    matrix with entries in R.
    
    INPUT:
        :g: an element in the group algebra
    
    OUTPUT:
        A matrix which represents the linear map 
    """
    return g.to_matrix(base_ring=ZZ,side=side)

def to_integer_matrix(AA):
    """
    This uses the fact that for a ring R and a group G
    we have an isomorphism of R-modules between
    R[G] and R^ord(G)
    Thus every linear map between free R[G]-modules
    can also be seen as a map between free R modules
    
    INPUT:
        :AA: Matrix with entries in the group algebra over the integer
             ring

    RETURNS:
        The linear map given by the matrix as a matrix with
        integer coefficients
    """
    apply_to_rows = lambda x: map(lambda y: y.to_matrix(side='right').transpose(), x)
    list_of_matrices = map(apply_to_rows, AA.rows())
    if len(list_of_matrices) == 1:
        return list_of_matrices[0][0]
    else:
        return block_matrix(list_of_matrices)
    
def fox(word, gen, im_gens=None, ring=None):
    """
    Disclaimer: This is almost identical to the function already present
    in the sage library
    Here just a small thing was changed to make it work smoothly with the group ring    
        
    INPUT:
        :word: The word in the free group we want to derive
        
        :gen: A generator of the free group with respect to which
        the Fox-derivative is taken
    
        :im_gens: The user can give images of the generators of the free group
        The function derives the word, then sends each generator to its image
        (if those are specified).

    OUTPUT:
        The Fox-derivative of word with respect to the generator gen
    """
    if not gen in word.parent().generators():
            raise ValueError("Fox derivative can only be computed with respect to generators of the group")
    l = list(word.Tietze())
    if im_gens is None:
        F = word.parent()
        R = F.algebra(IntegerRing())
        R_basis = R.basis()
        symb = [R_basis[a] for a in F.gens()]
        symb += reversed([R_basis[a.inverse()] for a in F.gens()])
        if ring is not None:
            R = ring
            symb = [R(i) for i in symb]
    else:
        if ring is None:
            R = Sequence(im_gens).universe()
        else:
            R = ring
        symb = list(im_gens)
        symb += reversed([R((a.trailing_support())^(-1)) for a in im_gens]) # We had to change this line
    i = gen.Tietze()[0]  # So ``gen`` is the `i`-th
                         # generator of the free group.
    a = R.zero()
    coef = R.one()
    while len(l) > 0:
        b = l.pop(0)
        if b == i:
            a += coef * R.one()
            coef *= symb[b-1]
        elif b == -i:
            a -= coef * symb[b]
            coef *= symb[b]
        elif b > 0:
            coef *= symb[b-1]
        else:
            coef *= symb[b]
    return a

def second_boundary_matrix_with_permutation(grouppresentation):
    """
    Get a subgroup of some permutation group which is isomorphic
    to the group 'grouppresentation'.
    Then write the matrix with elements in the
    group algebra over this permutation group.
    """
    perm = grouppresentation.as_permutation_group()
    Z_perm = GroupAlgebra(perm, ZZ)    
    images_in_groupalgebra = map(Z_perm, perm.gens())
    return second_boundary_matrix(grouppresentation, images_in_groupalgebra)

def second_boundary_matrix_with_integers(grouppresentation):
    """
    Calculated the second boundary matrix but with coefficients
    in the Integers
    """
    # d_2 is the second boundary matrix with coefficients in
    # the group ring
    d_2 = second_boundary_matrix_with_permutation(grouppresentation)
    return to_integer_matrix(copy(d_2))

def pi_2(grouppresentation):
    """
    Calculates the second homotopy group of the 2-dim CW-complex K
    associated to the group presentation given
    
    Uses the fact that the second homotopy group is isomorphic to
    the second homology group of the universal cover of K,
    H_2(K~, ZZ) = Ker(d_2)/Im(d_3)

    This in turn can be calculated as the kernel of the second
    boundary map, because Im(d_3) = 0 (no 3-cells)
    """
    d_2_integer = second_boundary_matrix_with_integers(grouppresentation)
    return d_2_integer.left_kernel()

# Functions for representing the kernel as Z[G]-module
# i.e. this is the direction Z^n ---> Z[G] of the isomorphism
def to_vector_with_permutation(vec, Zperm):
    order = len(Zperm.basis())
    # A vector can only be converted to the group algebra format if its length
    # is a multiple of the rank of the group algebra
    assert (vec.length() % order) == 0, "Vector length does not match up"
    no_components = vec.length() // order
    list = []
    for i in range(0, no_components):
        vector_part = vec[i*order : (i+1)*order]
        list.append(Zperm.from_vector(vector_part))
    return list

def to_groupalgebra_matrix(AA, Zperm):
    list_of_rows = map(lambda i: to_vector_with_permutation(i, Zperm), AA.rows())
    return matrix(list_of_rows)
    
def pi_2_with_permutation(grouppresentation):
    Zperm = GroupAlgebra(grouppresentation.as_permutation_group(),ZZ)    
    AA = pi_2(grouppresentation).matrix()
    return to_groupalgebra_matrix(AA, Zperm)
