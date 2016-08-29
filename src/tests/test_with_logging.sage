"""
Function to allow logging of some intermediate results
in the test of a group presentation G
"""
from datetime import datetime

# For logging Sage call
# logstart -o -t Z2xZ2xZ4_log
# before running the test

def test_with_logging(grouppresentation, grouppresentation_string):    
    """
    INPUT:
        :grouppresentation: The presentation of a finite, finitely presented group
        
        :grouppresentation_string: A text description of the grouppresentation,
        this will be used as the prefix for the log-files

    OUTPUT:
        Calculates the invariant factors of the grouppresentation
        and saves those in a Sage object file.
        Saves the string representations of some intermediate results
        in a text file.
        
        After starting the function Sage can be left unsupervised,
        all the relevant results are saved. 
    """
    gp_string=grouppresentation_string    
    f = file(gp_string + '_output.txt','w')
    f.write(str(grouppresentation) + '\n\n')
    
    f.write(str(second_boundary_matrix(grouppresentation)) + '\n\n')
    
    H_2=pi_2(grouppresentation)
    f.write(str(H_2) + '\n\n')
    f.write(H_2.matrix().str() + '\n\n')

    f.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' beginning to calculate the matrix \n\n')
    AA = matrix_of_surjection_use_just_generators(grouppresentation)
    f.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished to calculate the matrix \n\n')    
    save(AA, gp_string + '_matrix_of_surjection.sobj')
    f.write(str(AA) + '\n\n')
    
    f.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' beginning with invariant factors from matrix \n\n')
    inv_factors = structure_from_matrix_gap(AA)    
    f.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished with invariant factors from matrix \n\n')
    save(inv_factors, gp_string + '_structure.sobj')
    f.write(str(inv_factors) + '\n\n')
    
    f.close()

def test_with_logging_and_progress(grouppresentation, grouppresentation_string):
    gp_string=grouppresentation_string    
    f_output = file(gp_string + '_output.txt','w')
    # The argument '1' for the open method sets the file to be line buffered
    f_progress = file(gp_string + '_progress.txt', 'w', 1)
    
    f_output.write(str(grouppresentation) + '\n\n')
    f_output.write(str(second_boundary_matrix(grouppresentation)) + '\n\n')
    
    f_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' starting pi_2 calculation \n\n')
    H_2=pi_2(grouppresentation)
    f_output.write(str(H_2) + '\n\n')
    f_output.write(H_2.matrix().str() + '\n\n')    
    
    H_2_rk = H_2.rank()
    Gamma_rk = H_2_rk*(H_2_rk+1)/2
    f_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished pi_2 calculation \n\n')
    f_progress.write('rank_Z (pi_2 K) = ' + str(H_2_rk) + ' and thus rank_Z (Gamma pi_2 K) = ' + str(Gamma_rk) + '\n\n')

    f_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' calling matrix_of_surjection function \n')
    AA = matrix_of_surjection_use_just_generators_with_progress(grouppresentation, f_progress)
    f_progress.write('\n' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' matrix_of_surjection function returned \n\n')    
    save(AA, gp_string + '_matrix_of_surjection.sobj')
    #f_output.write(str(AA) + '\n\n')
    
    f_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' beginning with invariant factors from matrix (this might take a long time, no progress shown in between) \n')
    inv_factors = structure_from_matrix_gap(AA)    
    f_progress.write(datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ' finished with invariant factors from matrix \n\n')
    save(inv_factors, gp_string + '_structure.sobj')
    f_output.write(str(inv_factors) + '\n\n')
    
    f_output.close()
    f_progress.close()


# Specializations for some example groups

# G = < a, b, c | a^2, b^2, c^3, a*b*a^-1*b^-1, a*c*a^-1*c^-1, b*c*b^-1*c^-1 >
#   = Z/(2) x Z/(2) x Z/(3)
def test_Z2xZ2xZ3():
    # Warning: The functions exits Sage after it is done    
    F.<a,b,c>=FreeGroup(3)
    relations_Z2xZ2xZ3=[a^2, b^2, c^3, commutator(a,b), commutator(a,c), commutator(b,c)]
    G=F.quotient(relations_Z2xZ2xZ3)
    test_with_logging_and_progress(G, "Z2xZ2xZ3")
    exit()


# G = < a, b, c | a^2, b^2, c^4, a*b*a^-1*b^-1, a*c*a^-1*c^-1, b*c*b^-1*c^-1 >
#   = Z/(2) x Z/(2) x Z/(4)
def test_Z2xZ2xZ4():
    # Warning: The functions exits Sage after it is done    
    F.<a,b,c>=FreeGroup(3)
    relations_Z2xZ2xZ4=[a^2, b^2, c^4, commutator(a,b), commutator(a,c), commutator(b,c)]
    G=F.quotient(relations_Z2xZ2xZ4)
    test_with_logging_and_progress(G, "Z2xZ2xZ4")
    exit()

def test_Z2xZ2xZ2xZ2():
    # Warning: The functions exits Sage after it is done    
    G=product_of_cyclic_groups(2,4)
    test_with_logging_and_progress(G, "Z2^4")
    exit()

# G = < a, b, c | a^2, b^3, (a*b)^2, c^2, a^-1*c^-1*a*c, b^-1*c^-1*b*c >
#   = D_3 x Z/(2)
def test_D3xZ2():
    # Warning: The functions exits Sage after it is done
    G=dihedral_group(3).direct_product(cyclic_group(2))
    test_with_logging(G, "D3xZ2")
    exit()

# G = < a, b | a^2, b^8, a*b*a^-1*b^-1 >
#   = Z/(2) x Z/(8)
def test_Z2xZ8():
    # Warning: The functions exits Sage after it is done
    F.<a,b>=FreeGroup(2)
    relations_Z2xZ8=[a^2, b^8, a*b*a^-1*b^-1]
    G=F.quotient(relations_Z2xZ8)
    test_with_logging(G, "Z2xZ8")
    exit()
