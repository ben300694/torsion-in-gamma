"""
Main document to combine all the other modules and for testing purposes
Attach this to sage via the command "attach('main.sage')" to use all the functions

"""

# load the parts of the algorithm
attach("pi_2.sage")
attach("gamma.sage")
attach("torsion.sage")
attach("util.sage")

# load all the files containing examples of groups
attach("examples/Zmodn.sage")
attach("examples/(Zmodk)^n.sage")
attach("examples/D_n.sage")
attach("examples/Q_8n.sage")

# load the tests
attach("tests/test_suite.sage")
attach("tests/test_with_logging.sage")
