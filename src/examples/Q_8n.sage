"""
Generalized quaternion group of order 8*n,
Q_8n = < x, y | x^(2*n)*y^(-2), x*y*x*y^(-1) >
"""

def quaternion_group(n):
    F2.<x, y> = FreeGroup(2)
    relations_Q8n = [x^(2*n)*y^(-2), x*y*x*y^(-1)]
    Q8n = F2.quotient(relations_Q8n)
    return Q8n

Q8 = quaternion_group(1)
Q16 = quaternion_group(2)
Q24 = quaternion_group(3)
Q32 = quaternion_group(4)
