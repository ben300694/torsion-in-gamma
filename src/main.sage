"""
Main document to combine all the other modules and for testing purposes
Attach this to sage via the command "attach('main.sage')" to use all the functions

Author: Benjamin Ruppik (University of Bonn)
Date: June 2016
"""

# load the parts of the algorithm
attach("pi_2.sage")
attach("gamma.sage")
attach("torsion.sage")
attach("util.sage")

# load all the files containing examples of groups
attach("examples/Zn.sage")
attach("examples/Zk^n.sage")
attach("examples/Dn.sage")

# load the tests
attach("tests/test_suite.sage")
attach("tests/test_with_logging.sage")
