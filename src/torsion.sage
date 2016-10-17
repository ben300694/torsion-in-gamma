"""
Module for combining all the previous functions
and defining the function which finds the torsion in Gamma(pi_2 K)/(pi_1 K)

Author: Benjamin Ruppik (University of Bonn)
Date: June 2016

Credit:
This uses pi_2.sage, gamma.sage and all their dependencies (see the credit there)
Some of the functions use algorithms which are implemented in GAP.
"""

def basis_for_symmetric_matrices(n):
    """    
    Gives a specific basis for the vector space of symmetric n x n - matrices
    These matrices give the components for the basis of the subspace Gamma(M) 
    in M \otimes M we will use
    """
    diagonal = [matrix(n, n, {(i,i): 1}) for i in range(0,n)]
    others = [matrix(n, n, {(i,j): 1, (j,i): 1}) for i in range(0,n) for j in range(i+1, n)]
    return diagonal + others

def subspace_of_symmetric_matrices(n):
    """
    INPUT:
        :dim: integer specifying the dimension of the matrices involved

    OUTPUT:
        The submodule spanned by the symmetric (n x n)-matrices
    """
    basis_vectors = map(matrix_to_list, basis_for_symmetric_matrices(n))
    return span(basis_vectors, ZZ)

def tensors_identified_with_zero(grouppresentation):
    """    
    To caculate the quotient of the pi_1-action on pi_2
    we need to find the map
    ZZ^(rk(Gamma(M))*ord(G)) ---> Gamma(M) = {symmetric rk(M) x rk(M) - matrices}
    """
    H_2 = pi_2(grouppresentation)

    perm = grouppresentation.as_permutation_group()
    Z_perm = GroupAlgebra(perm, ZZ)
    groupelements_in_algebra = map(Z_perm, perm.list())
    groupelements_in_algebra_as_matrix = map(lambda x: x.to_matrix(side='left'), groupelements_in_algebra)
    
    free_module = FiniteRankFreeModule(ZZ, H_2.rank())
    f = free_module.basis('f')

    tensors = []
    
    # the tensors of the form m \otimes m
    for m in H_2.basis():
        for g in groupelements_in_algebra_as_matrix[1:]:
            # We need to calculate m \otimes m - g*(m \otimes m)
                t = diagonal_tensor(m, H_2, free_module) - image_of_diagonal_tensor(g, m, H_2, free_module)
                tensors.append(t)
    
    # the tensors of the form m \otimes n + n \otimes m
    for (i, m) in enumerate(H_2.basis()):
        for j in range(i+1, len(H_2.basis())):
            n = H_2.basis()[j]
            for g in groupelements_in_algebra_as_matrix[1:]:
                # We need to calculate m \otimes n + n \otimes m - g*(m \otimes n + n \otimes m)                
                t = other_tensor(m, n, H_2, free_module) - image_of_other_tensor(g, m, n, H_2, free_module)
                tensors.append(t)

    return tensors

def tensors_identified_with_zero_use_just_generators(grouppresentation):
    """
    The same as tensors_identified_with_zero(grouppresentation),
    but only iterates over the generators of the group
    """    
    H_2 = pi_2(grouppresentation)
    perm = grouppresentation.as_permutation_group()
    Z_perm = GroupAlgebra(perm, ZZ)
    groupgens_in_algebra = map(Z_perm, perm.gens())
    groupgens_in_algebra_as_matrix = map(lambda x: x.to_matrix(side='left'), groupgens_in_algebra)
    # print groupgens_in_algebra_as_matrix    

    free_module = FiniteRankFreeModule(ZZ, H_2.rank())
    f = free_module.basis('f')

    tensors = []
    
    # the tensors of the form m \otimes m
    for m in H_2.basis():
        for g in groupgens_in_algebra_as_matrix:
            # We need to calculate m \otimes m - g*(m \otimes m)
                t = diagonal_tensor(m, H_2, free_module) - image_of_diagonal_tensor(g, m, H_2, free_module)
                tensors.append(t)
    
    # the tensors of the form m \otimes n + n \otimes m
    for (i, m) in enumerate(H_2.basis()):
        for j in range(i+1, len(H_2.basis())):
            n = H_2.basis()[j]
            for g in groupgens_in_algebra_as_matrix:
                # We need to calculate m \otimes n + n \otimes m - g*(m \otimes n + n \otimes m)                
                t = other_tensor(m, n, H_2, free_module) - image_of_other_tensor(g, m, n, H_2, free_module)
                tensors.append(t)

    return tensors

def tensors_identified_with_zero_use_just_generators_with_progress(grouppresentation, file_progress):
    """
    The same as tensors_identified_with_zero(grouppresentation),
    but only iterates over the generators of the group
    """
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero preparations\n')
        
    H_2 = pi_2(grouppresentation)
    perm = grouppresentation.as_permutation_group()
    Z_perm = GroupAlgebra(perm, ZZ)
    groupgens_in_algebra = map(Z_perm, perm.gens())
    groupgens_in_algebra_as_matrix = map(lambda x: x.to_matrix(side='left'), groupgens_in_algebra)
    # print groupgens_in_algebra_as_matrix    

    free_module = FiniteRankFreeModule(ZZ, H_2.rank())
    f = free_module.basis('f')

    tensors = []
    
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero starting with tensors of the form m \otimes m\n')
    # the tensors of the form m \otimes m
    for m in H_2.basis():
        for g in groupgens_in_algebra_as_matrix:
            # We need to calculate m \otimes m - g*(m \otimes m)
                t = diagonal_tensor(m, H_2, free_module) - image_of_diagonal_tensor(g, m, H_2, free_module)
                tensors.append(t)
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero finished with tensors of the form m \otimes m\n')
    
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero starting with tensors of the form m \otimes n + n \otimes m\n')
    # the tensors of the form m \otimes n + n \otimes m
    for (i, m) in enumerate(H_2.basis()):
        file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' in outer loop ' + str(i+1) + ' of ' + str(H_2.rank()) + '\n')
        for j in range(i+1, len(H_2.basis())):
            n = H_2.basis()[j]
            for g in groupgens_in_algebra_as_matrix:
                # We need to calculate m \otimes n + n \otimes m - g*(m \otimes n + n \otimes m)                
                t = other_tensor(m, n, H_2, free_module) - image_of_other_tensor(g, m, n, H_2, free_module)
                tensors.append(t)
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero finished with tensors of the form m \otimes n + n \otimes m\n')

    return tensors

def matrix_of_surjection(grouppresentation):
    """
    Gives a matrix representation of the map from
    a free module to Gamma(pi_2)
    that sends each basis element to one of the tensors identified
    with zero
    """
    M = pi_2(grouppresentation)   
    Sym = subspace_of_symmetric_matrices(M.rank())
    list_of_tensors = tensors_identified_with_zero(grouppresentation)

    # Just write the coefficients of the linear combination
    # as the columns
    columns_coordinates = map(lambda x: Sym.coordinates(tensor_to_list(x)), list_of_tensors)
    return matrix(ZZ, columns_coordinates).transpose()

def matrix_of_surjection_use_just_generators(grouppresentation):
    M = pi_2(grouppresentation)   
    Sym = subspace_of_symmetric_matrices(M.rank())
    list_of_tensors = tensors_identified_with_zero_use_just_generators(grouppresentation)

    columns_coordinates = map(lambda x: Sym.coordinates(tensor_to_list(x)), list_of_tensors)
    return matrix(ZZ, columns_coordinates).transpose()

def matrix_of_surjection_use_just_generators_with_progress(grouppresentation, file_progress):
    M = pi_2(grouppresentation)   
    
    file_progress.write('\n' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' generating subspace of symmetric matrices (this might take a long time, no progress shown in between)\n')    
    Sym = subspace_of_symmetric_matrices(M.rank())
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished with subspace of symmetric matrices \n\n')

    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' calling tensors_identified_with_zero \n\n')
    list_of_tensors = tensors_identified_with_zero_use_just_generators_with_progress(grouppresentation, file_progress)
    file_progress.write('\n' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' tensors_identified_with_zero returned \n\n')

    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' calculating coordinates from list of tensors \n')

    list_of_tensors_len = len(list_of_tensors)
    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' the list_of_tensors has length ' + str(list_of_tensors_len) + ' \n')
    
    # Split the list_of_tensors into chunks of length n = list_of_tensors_len/planned_no_of_chunks or if this is less than 200
    # into chunks of length 200
    planned_no_of_chunks = 30
    list_of_tensors_len = len(list_of_tensors)    
    n = max(list_of_tensors_len // planned_no_of_chunks, 200)
    chunks = [list_of_tensors[i:i + n] for i in xrange(0, list_of_tensors_len, n)]
    total_no_of_chunks = len(chunks)    

    # Call the Sym_coordinates function on each chunk (this executes in parallel)
    r = Sym_coordinates([(chunk, Sym, file_progress, j+1, total_no_of_chunks) for (j, chunk) in enumerate(chunks)])    
    list_columns_coordinates = []
    for x in r:
        list_columns_coordinates.append(x[1])
    columns_coordinates = [item for sublist in list_columns_coordinates for item in sublist] 

    file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished with coordinates from list of tensors \n')
    
    return matrix(ZZ, columns_coordinates).transpose()


@parallel
def Sym_coordinates(list_of_tensors, subspace, file_progress, no_of_chunk, total_no_of_chunks):
    list_of_tensors_len = len(list_of_tensors)    
    result = []
    for (i, x) in enumerate(list_of_tensors):
        file_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' entry ' + str(i+1) + ' of ' + str(list_of_tensors_len) + ' in chunk ' + str(no_of_chunk) + ' of ' + str(total_no_of_chunks) + '\n')
        y = subspace.coordinates(tensor_to_list(x))
        result.append(y)
    return result


#
# Test to see if removing the abstraction of the tensor product
# increases the performance
#

def tensors_identified_with_zero_as_matrices(grouppresentation):
    H_2 = pi_2(grouppresentation)
    # To caculate the quotient of the pi_1-action on pi_2
    # we need to find the map
    # ZZ^(rk(Gamma(M))*#gens.) ---> Gamma(M) = {symmetric rk(M) x rk(M) - matrices}

    perm = grouppresentation.as_permutation_group()
    Z_perm = GroupAlgebra(perm, ZZ)
    groupgens_in_algebra = map(Z_perm, perm.gens())
    groupgens_in_algebra_as_matrix = map(lambda x: x.to_matrix(side='left'), groupgens_in_algebra)
    # print groupelements_in_algebra_as_matrix

    tensors_as_matrix = []
    
    # the tensors of the form m \otimes m
    for m in H_2.basis():
        for g in groupgens_in_algebra_as_matrix:
            # We need to calculate m \otimes m - g*(m \otimes m)
                m_coord_mat=matrix(H_2.coordinates(m))
                gm_coord_mat=matrix(action_in_coordinates(g, m, H_2))          
                t = (m_coord_mat.transpose())*m_coord_mat - (gm_coord_mat.transpose())*gm_coord_mat
                tensors_as_matrix.append(t)
    
    # the tensors of the form m \otimes n + n \otimes m
    for (i, m) in enumerate(H_2.basis()):
        for j in range(i+1, len(H_2.basis())):
            n = H_2.basis()[j]
            for g in groupgens_in_algebra_as_matrix:
                # We need to calculate m \otimes n + n \otimes m - g*(m \otimes n + n \otimes m)                
                m_coord_mat=matrix(H_2.coordinates(m))
                gm_coord_mat=matrix(action_in_coordinates(g, m, H_2))          
                n_coord_mat=matrix(H_2.coordinates(n))
                gn_coord_mat=matrix(action_in_coordinates(g, n, H_2))          
                                
                other_tensor = (m_coord_mat.transpose())*n_coord_mat + (n_coord_mat.transpose())*m_coord_mat
                image_of_other_tensor = (gm_coord_mat.transpose())*gn_coord_mat + (gn_coord_mat.transpose())*gm_coord_mat
                t = other_tensor - image_of_other_tensor                
                tensors_as_matrix.append(t)

    return tensors_as_matrix

def matrix_of_surjection_as_matrices(grouppresentation):
    M = pi_2(grouppresentation)   
    Sym = subspace_of_symmetric_matrices(M.rank())
    list_of_tensors_as_matrices = tensors_identified_with_zero_as_matrices(grouppresentation)

    columns_coordinates = map(lambda x: Sym.coordinates(matrix_to_list(x)), list_of_tensors_as_matrices)
    return matrix(columns_coordinates).transpose()

# Sparse
#
# Converts the matrices to sparse matrices before applying the
# .elementary_divisors() method
def structure_of_gamma_pi_2_quotient_sparse(grouppresentation):
    AA = matrix_of_surjection(grouppresentation)
    BB = AA.sparse_matrix() 
    return BB.elementary_divisors()

def has_torsion_in_gamma_pi_2_quotient(grouppresentation):
    factors = structure_of_gamma_pi_2_quotient_sparse(grouppresentation)
    return contains_element_not_zero_or_one(factors)

# GAP
# 
# The algorithms that use a GAP-implementation
def structure_from_matrix_gap(AA):
    AA_smith = smith_form_with_gap(AA)
    return diagonal(AA_smith)

def structure_of_gamma_pi_2_quotient_gap(grouppresentation):
    AA = matrix_of_surjection(grouppresentation)
    return structure_from_matrix_gap(AA)

def structure_of_gamma_pi_2_quotient_gap_use_just_generators(grouppresentation):
    AA = matrix_of_surjection_use_just_generators(grouppresentation)
    return structure_from_matrix_gap(AA)

def structure_of_gamma_pi_2_quotient_gap_as_matrices(grouppresentation):
    AA = matrix_of_surjection_as_matrices(grouppresentation)
    return structure_from_matrix_gap(AA)

def has_torsion_in_gamma_pi_2_quotient_gap_use_just_generators(grouppresentation):
    factors = structure_of_gamma_pi_2_quotient_gap_use_just_generators(grouppresentation)
    return contains_element_not_zero_or_one(factors)

# Alias for the most important functions
# that the user will use

def struc(grouppresentation):
    return structure_of_gamma_pi_2_quotient_gap_use_just_generators(grouppresentation)

def has_tors(grouppresentation):
    return has_torsion_in_gamma_pi_2_quotient_gap_use_just_generators(grouppresentation)
    
