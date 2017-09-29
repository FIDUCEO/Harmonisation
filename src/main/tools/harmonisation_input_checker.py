"""
Check harmonisation input files for internal consistancy

Usage:
python add_W_to_matchup_file.py path/to/matchup/file.nc
"""

'''___Python Modules___'''

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "28/09/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

# Ideas for tests

# > Variable content:
# - Must be variables: H, Ur, Us, K, Kr, Kr
# - If any of w matrix variables all must be there

# > Dimension of variables
# - H(M, m), Ur(M, m), Us(M, m)
# - K(m), Kr(m), Ks(m)
# - M must be even

# > Consistancy of W matrix variables
# - Either Ur, Ur & Us or W
# - Number of Ws linked to maps length of w_matrix_nnz entries
# - 
