"""
Module to implement the Gamma-Functor for a finitely generated free module M.
It also contains methods to calculate the diagonal action of a group G on a tensor product.

Suppose we are given an action G x M ---> M
                               (g,m) |--> g*m

If we view Gamma(M) as subgroup in M \otimes M we get the diagonal action
g * (m \otimes n) := (g*m) \otimes (g*n)


Author: Benjamin Ruppik (University of Bonn)
Date: June 2016

Credit:
This uses the packages
- sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule
- sage.tensor.modules.tensor_free_module.TensorFreeModule
both written by Eric Gourgoulhon and Michal Bejger
"""

def action(g, vect):
    """
    INPUT
        :g: is an element of the group algebra Z[G] of the group G
        
        :vect: A vector with integer coefficients 
            is interpreted as an element of a Z[G]-module where the action of Z[G]
            is given by multiplication from the left
            
            The output basis is the same as the input basis (i.e. the canonical
            one for the group ring)
            
    """
    order = len(g.parent().basis())
    # g can act on a vector only if its length
    # is a multiple of the rank of the group algebra
    assert (vect.length() % order) == 0, "Vector length does not match up"
    no_components = vect.length() // order
    # get the matrix representation of the action of g on the group algebra
    g_matrix = g.to_matrix(side='left')
    result = []
    for i in range(0, no_components):
        vector_part = vect[i*order : (i+1)*order]
        new_vector_part = g_matrix * vector_part
        result += new_vector_part.list()
    return result

def action_with_g_given_as_matrix(g_matrix, vect):
    """
    INPUT
        :g: is matrix representing left multiplication with g
        
        :vect: A vector with integer coefficients 
            is interpreted as an element of a Z[G]-module where the action of Z[G]
            is given by multiplication from the left
            
            The output basis is the same as the input basis (i.e. the canonical
            one for the group ring)
            
    """
    order = len(g_matrix.rows())
    # g can act on a vector only if its length
    # is a multiple of the rank of the group algebra
    assert (vect.length() % order) == 0, "Vector length does not match up"
    no_components = vect.length() // order
    # get the matrix representation of the action of g on the group algebra
    result = []
    for i in range(0, no_components):
        vector_part = vect[i*order : (i+1)*order]
        new_vector_part = g_matrix * vector_part
        result += new_vector_part.list()
    return result

def action_in_coordinates(g, vect, module):
    """    
    The .coordinates() method returns the coordinates of the
    vector in the default given basis of the module
    It is these coordinates which we will later use in the application
    of the gamma functor
    """
    return module.coordinates(action_with_g_given_as_matrix(g, vect))

def diagonal_tensor(m, module, M):
    """
    Return m \otimes m
    """    
    m_coord = module.coordinates(m)
    return M(m_coord)*M(m_coord)

def image_of_diagonal_tensor(g, m, module, M):
    """
    Return g*(m \otimes m) := (g*m) \otimes (g*m)
    """
    gm_coord = action_in_coordinates(g, m, module)
    return M(gm_coord)*M(gm_coord)

def other_tensor(m, n, module, M):
    """
    Return m \otimes n + n \otimes m

    Observe that other_tensor(m, n, module) == other_tensor(n, m, module)
    """
    m_coord = module.coordinates(m)
    n_coord = module.coordinates(n)
    return M(m_coord)*M(n_coord) + M(n_coord)*M(m_coord)

def image_of_other_tensor(g, m, n, module, M):
    """
    Return g*(m \otimes n + n \otimes m) := (g*m) \otimes (g*n) + (g*n) \otimes (g*m)

    Observe that image_of_other_tensor(m, n, module) == image_of_other_tensor(n, m, module)
    """
    gm_coord = action_in_coordinates(g, m, module)
    gn_coord = action_in_coordinates(g, n, module)
    return M(gm_coord)*M(gn_coord) + M(gn_coord)*M(gm_coord)


def tensor_to_list(t):
    """
    Only possible for type (2,0)-tensors
    Converts the matrix representation of 't'
    into a (long) list
    """
    assert (t.tensor_type() == (2, 0)), "Tensor of this type cannot be converted to list"
    return matrix_to_list(matrix(t.components()))
