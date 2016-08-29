import collections

# testMatrix from exercise 2.9/4 in [bosch2013algebra]
# The build in method .smith_form() returns (
# [ 1  0  0]  [ 0  0  1]  [  1   1   5]
# [ 0  2  0]  [-1  1 -3]  [ -4  -5 -29]
# [ 0  0 36], [ 3 -2  4], [  3   4  25]
# )
testMatrix = matrix([[2,6,8],[3,1,2],[9,5,4]])
testMatrix2 = matrix(ZZ,3,5, lambda i, j: - i - j)
print "testMatrix = " ; print testMatrix

# Function for calculating the Smith Normal Form of a matrix [https://en.wikipedia.org/wiki/Smith_normal_form]
# We want both the normal form and the base changes that transform the matrix
# Maybe later we will use the library function already present for this, but at this point implement it on our own
def snf(AA, Left, Right):
    """:AA: input matrix  
       :Left: the left base change matrix
       :Right: the right base change matrix
    """
    if Left == None:
        Left = matrix.identity(len(AA.rows()))
    if Right == None:
        Right = matrix.identity(len(AA.columns()))

    MatrixBaseChange = collections.namedtuple('MatrixBaseChange', ['A', 'Left', 'Right'])
    
    # If we have the zero-matrix then we are done        
    if AA == 0:
        result = MatrixBaseChange(AA, Left, Right)
        return result
    
    pos, val = minimal_delta_position_and_value(AA)
    i, j = pos
    assert val != 0
    
    if pos != [0,0]:
        map(lambda x: x.swap_rows(0, i), [AA, Left])
        map(lambda x: x.swap_columns(0, j), [AA, Right])
    
    # !!! TODO
    for i, row in enumerate(AA.rows()):
        quot, rem = divmod(AA[i,0], AA[0,0])
        if rem != 0:
            AA.add_multiple_of_row(i, 0, -quot)
            Left.add_multiple_of_row(i, 0, quot)
            #return snf(AA, Left, Right)

    return MatrixBaseChange(AA, Left, Right)
        

def minimal_delta_position(AA):
        """:AA:
        
           :return:
        """       
        #Just return the first element of the tuple from the maximal_delta_position_and_value function
        return minimal_delta_position_and_value(AA)[0]

def minimal_delta_value(AA):
        """:AA: input matrix
           
           :return: the maximal absolute value of the entries of the matrix         
        """
        #Just return the second element of the tuple from the maximal_delta_position_and_value function
        return minimal_delta_position_and_value(AA)[1]

def minimal_delta_position_and_value(AA):
    """:AA: input matrix with values in in Euclidean ring
    
       :return: a named tuple, the first entry is the position of the element not equal to 0 with minmal euclidean degree,
                the second entry is the value of this delta
                returns 0 iff the matrix is equal to zero
    """
    PositionValue = collections.namedtuple('PositionValue', ['position', 'value'])

    if AA == 0:
        return PositionValue([0,0],0)
       
        #The degree function is always an function delta: R-{0} --> NN
        #so minimal_delta takes values in the natural numbers
        #Observe that for R = ZZ the ring of integers, the euclidean degree is just the absolute value
    position = nonzero_entry_position(AA)
    i, j = position            
    minimal_delta = AA[i,j].euclidean_degree()  

    for i, row in enumerate(AA.rows()):
        for j, column in enumerate(AA.columns()):
            if (AA[i,j] != 0) and (AA[i,j].euclidean_degree() < minimal_delta):
                minimal_delta = AA[i,j].euclidean_degree()
                position = [i,j]                        
            else:
                continue

    return PositionValue(position, minimal_delta)

def nonzero_entry_position(AA):
    """:AA: input matrix which is nonzero
       :return: the position of the first entry which is nonzero, when going over the rows from top to bottom    
    """
    for i, row in enumerate(AA.rows()):
        for j, column in enumerate(AA.columns()):
            if AA[i,j] != 0:
                return [i,j]
    raise ValueError("nonzero_entry_position called with zero matrix")

def apply_base_change((AA, Left, Right)):
    return Left * AA * Right
    
