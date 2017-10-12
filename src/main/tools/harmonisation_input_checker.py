"""
Check harmonisation input files match specification and are internally consistant

Usage:
python harmonisation_input_checker.py path/to/matchup/file.nc
"""

'''___Python Modules___'''
import sys

'''___Third Party Modules___'''
from netCDF4 import Dataset
from numpy import sum, where, array_equal, isnan, zeros

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "28/09/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

'''___Constants___'''
# VARIABLE DATA

# Required variables and dimensions
VARIABLE_DATA = {"H": {"dim": ["M", "m"], "dtype": "float64"},
                 "Ur": {"dim": ["M", "m"], "dtype": "float64"},
                 "Us": {"dim": ["M", "m"], "dtype": "float64"},
                 "K": {"dim": ["M"], "dtype": "float64"},
                 "Kr": {"dim": ["M"], "dtype": "float64"},
                 "Ks": {"dim": ["M"], "dtype": "float64"}}

# Optional W matrix variable dimensions (if included complete set required)
W_VARIABLE_DATA = {"w_matrix_nnz": {"dim": ["w_matrix_count"], "dtype": "int32"},
                   "w_matrix_row": {"dim": ['w_matrix_count', 'w_matrix_num_row'], "dtype": "int32"},
                   "w_matrix_col": {"dim": ["w_matrix_sum_nnz"], "dtype": "int32"},
                   "w_matrix_val": {"dim": ["w_matrix_sum_nnz"], "dtype": "float64"},
                   "w_matrix_use": {"dim": ["m"], "dtype": "int32"},
                   "uncertainty_vector_row_count": {"dim": ["uncertainty_vector_count"], "dtype": "int32"},
                   "uncertainty_vector": {"dim": ["uncertainty_vector_sum_row"], "dtype": "float64"},
                   "uncertainty_vector_use": {"dim": ["m"], "dtype": "int32"}}

# Attributes
ATTRS = ["sensor_i_name", "sensor_j_name"]


def check_variable_included(rootgrp):
    """
    Return errors found with matchup dataset included variables

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    # TEST: Assert basic required variables included -------------------------------------------------------------------

    # Find if any variables from the required list not contained in file
    required_variables = VARIABLE_DATA.keys()
    missing_variables = []
    for required_variable in required_variables:
        if required_variable not in rootgrp.variables.keys():
            missing_variables.append(required_variable)

    # If variables missing return error message with missing variable names
    if missing_variables != []:
        errors.append("Missing Variable[s]: '" + str(missing_variables) + "' not in file")
    # ------------------------------------------------------------------------------------------------------------------

    # TEST: Assert if any w variables included then all are included ---------------------------------------------------

    # Find if any w variables from the list are not contained in file
    w_variables = W_VARIABLE_DATA.keys()
    missing_w_variables = []
    for w_variable in w_variables:
        if w_variable not in rootgrp.variables.keys():
            missing_w_variables.append(w_variable)

    # Return error message with file w variables is partial set of full list of w variables
    if (missing_w_variables != []) and (missing_w_variables != w_variables):
        errors.append("Incomplete W Variable[s]: '" + str(missing_w_variables) + "' not in file")
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_variable_dimension(rootgrp):
    """
    Return errors found with matchup dataset variable dimensions

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    # TEST: Assert variables in file are of specified required dimensions ----------------------------------------------

    # Read file variable names and dimensions
    file_dims = {}
    for variable in rootgrp.variables.keys():
        file_dims[str(variable)] = [str(i) for i in rootgrp.variables[variable]._getdims()]

    # Check dimensions of variables in file against required values
    for variable in file_dims.keys():
        test_dims = file_dims[variable]

        expected_dims = []
        if variable in VARIABLE_DATA.keys():
            expected_dims = VARIABLE_DATA[variable]["dim"]
        elif variable in W_VARIABLE_DATA.keys():
            expected_dims = W_VARIABLE_DATA[variable]["dim"]

        if (expected_dims != []) and (test_dims != expected_dims):
            errors.append("Dimension Error: Dimension of '" + variable + "' must be " + str(expected_dims) + ", not "
                          + str(test_dims))
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_variable_dtype(rootgrp):
    """
    Return errors found with matchup dataset variable dtypes

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    # TEST: Assert variables in file are of specified required dimensions ----------------------------------------------

    # Read file variable names and dimensions
    file_dtype = {}
    for variable in rootgrp.variables.keys():
        file_dtype[str(variable)] = str(rootgrp.variables[variable].dtype)

    # Check dimensions of variables in file against required values
    for variable in file_dtype.keys():
        test_dtype = file_dtype[variable]

        expected_dtype = []
        if variable in VARIABLE_DATA.keys():
            expected_dtype = VARIABLE_DATA[variable]["dtype"]
        elif variable in W_VARIABLE_DATA.keys():
            expected_dtype = W_VARIABLE_DATA[variable]["dtype"]

        if (expected_dtype != []) and (test_dtype != expected_dtype):
            errors.append("Data Type Error: dtype of '" + variable + "' must be " + str(expected_dtype) + ", not "
                          + str(test_dtype))
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_attributes(rootgrp):
    """
    Return errors found with matchup dataset attributes

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    # TEST: Assert required attributes in file  ------------------------------------------------------------------------

    test_attrs = [str(attr) for attr in rootgrp.ncattrs()]

    for attr in ATTRS:
        if attr not in test_attrs:
            errors.append("Attribute Error: Attribute '"+str(attr)+"' missing")
    # ------------------------------------------------------------------------------------------------------------------

    # TEST: Assert sensor_i_name and sensor_j_name different -----------------------------------------------------------
    if ("sensor_i_name" in test_attrs) and ("sensor_j_name" in test_attrs):
        if rootgrp.sensor_i_name == rootgrp.sensor_j_name:
            errors.append("Attribute Error: Attributes sensor_i_name & sensor_j_name must have different values")
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_dimension_sizes(rootgrp):
    """
    Return errors found with matchup dataset dimension sizes

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    # TEST: Assert that dimension 'm' is even --------------------------------------------------------------------------
    if "m" in rootgrp.dimensions:
        if rootgrp.dimensions["m"].size % 2 != 0:
            errors.append("Dimension Error: Dimension 'm' must be even -> "
                          "sensor 1 and sensor 2 must share same measurement equation")
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_w_matrix_variable_dimensions(rootgrp):
    """
    Return errors found with matchup dataset w matrix variable dimensions

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    if set(W_VARIABLE_DATA.keys()).issubset(rootgrp.variables.keys()):

        # TEST: Assert the size of dimension 'w_matrix_count' is equal to the number of w matrices indexed in variable -
        #       'w_matrix_use' -----------------------------------------------------------------------------------------
        dim_num_ws = rootgrp.dimensions["w_matrix_count"].size
        use_num_ws = max(rootgrp.variables["w_matrix_use"][:])
        if dim_num_ws != use_num_ws:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_count' ("+str(dim_num_ws) +
                          ") must match number of labelled W matrices in variable 'w_matrix_use'(" + str(use_num_ws) +
                          ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'uncertainty_vector_count' is equal to the number of uncertainty vectors --
        #       indexed in variable 'uncertainty_vector_use' -----------------------------------------------------------
        dim_num_u_vecs = rootgrp.dimensions["uncertainty_vector_count"].size
        use_num_u_vecs = max(rootgrp.variables["uncertainty_vector_use"][:])
        if dim_num_u_vecs != use_num_u_vecs:
            errors.append("Uncertainty Vector Dimension Error: Size of dimension 'uncertainty_vector_count' ("
                          + str(dim_num_u_vecs) + ") must match number of labelled uncertainty vectors in variable"
                          " 'uncertainty_vector_use'(" + str(use_num_u_vecs) + ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert that the row size of w matrices is equal to the number of matchups ------------------------------
        w_matrix_num_row_size = rootgrp.dimensions["w_matrix_num_row"].size
        M_size = rootgrp.dimensions["M"].size
        if w_matrix_num_row_size != M_size+1:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_row_num' (" +
                          str(w_matrix_num_row_size) + ") must equal number of matchups + 1 (" + str(M_size) +
                          " + 1 = "+str(M_size+1)+")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'w_matrix_sum_nnz' is equal to the sum of each w matrix nnz combined, as --
        #       given in the variable 'w_matrix_nnz' -------------------------------------------------------------------
        w_matrix_sum_nnz_size = rootgrp.dimensions["w_matrix_sum_nnz"].size
        w_matrix_nnz_combined = sum(rootgrp.variables["w_matrix_nnz"][:])
        if w_matrix_sum_nnz_size != w_matrix_nnz_combined:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_sum_nnz' (" +
                          str(w_matrix_sum_nnz_size) + ") must equal combined per w matrix nnz's contained in variable"
                          "'w_matrix_nnz' (" + str(w_matrix_nnz_combined) + ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'uncertainty_vector_sum_row' is equal to the sum of each uncertainty vector
        #       number of rows combined, as given in the variable 'uncertainty_vector_row_count' -----------------------
        uncertainty_vector_sum_row_size = rootgrp.dimensions["uncertainty_vector_sum_row"].size
        uncertainty_vector_row_count_combined = sum(rootgrp.variables["uncertainty_vector_row_count"][:])
        if uncertainty_vector_sum_row_size != uncertainty_vector_row_count_combined:
            errors.append("Uncertainty Vector Dimension Error: Size of dimension 'uncertainty_vector_sum_row' ("
                          + str(uncertainty_vector_sum_row_size) + ") must equal combined per uncertainty vector"
                          " row counts contained in variable 'uncertainty_vector_row_count' (" +
                          str(uncertainty_vector_row_count_combined) + ")")
        # --------------------------------------------------------------------------------------------------------------

    return errors


def check_variable_values(rootgrp):
    """
    Return errors found with matchup dataset required variable values

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []  # Initialise list to store error messages

    # TEST: Assert uncertainties in Kr and Ks must positive and above zero in combination ------------------------------
    if set(["Kr", "Ks"]).issubset(rootgrp.variables.keys()):
        if not (rootgrp.variables['Kr'][:] >= 0).all():
            errors.append("Value Error: Negative value[s] found in 'Kr' variable")
        if not (rootgrp.variables['Ks'][:] >= 0).all():
            errors.append("Value Error: Negative value[s] found in 'Ks' variable")
        if not (rootgrp.variables['Ks'][:]**2+rootgrp.variables['Kr'][:]**2 > 0).all():
            errors.append("Value Error: Combined uncertainties of 'Kr' and 'Ks' values must be greater than zero, "
                          "exception[s] found")
    # ------------------------------------------------------------------------------------------------------------------

    # TEST: Assert Required Values of H, Ur , Us NaN for reference sensor pair -----------------------------------------
    test_variables = ["H", "Ur", "Us"]
    if (set(test_variables).issubset(rootgrp.variables.keys())) and ("sensor_i_name" in rootgrp.ncattrs()):
        if rootgrp.sensor_i_name == -1:
            for test_variable in test_variables:
                if not isnan(rootgrp.variables[test_variable][:, 1:rootgrp.dimensions['m'].size/2]).all():
                    errors.append("Value Error: All values "+test_variable+"[:, 1:"
                                  + str(rootgrp.dimensions['m'].size/2) + "] for reference-sensor pair must be NaN")

        for test_variable in test_variables:
            if isnan(rootgrp.variables[test_variable][:, 1:rootgrp.dimensions['m'].size / 2]).all():
                if rootgrp.sensor_i_name != -1:
                    errors.append("Value Error: All values " + test_variable + "[:, 1:" +
                                  str(rootgrp.dimensions['m'].size / 2) +
                                  "] are NaN, only allowed for reference-sensor pair - i.e. attribute 'sensor_i_name'"
                                  "must be -1 (not +" + str(rootgrp.sensor_i_name) + ")")

    # ------------------------------------------------------------------------------------------------------------------

    # TEST: Assert no negative uncertainties in Ur and Us --------------------------------------------------------------
    test_variables = ["Ur", "Us"]
    if set(test_variables).issubset(rootgrp.variables.keys()):
        for test_variable in test_variables:
            test_variable_val = rootgrp.variables[test_variable][:]
            if not (test_variable_val.flatten()[~isnan(test_variable_val.flatten())] >= 0).all():
                errors.append("Value Error: All values in " + test_variable + " must be greater than 0")
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def check_w_variable_values(rootgrp):
    """
    Return errors found with matchup dataset w matrix variable values

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []     # Initialise list to store error messages

    if set(W_VARIABLE_DATA.keys()).issubset(rootgrp.variables.keys()):

        # TEST: Assert W Matrix indices in w_matrix_use are in numerical order from 1 ----------------------------------
        w_matrix_use = rootgrp.variables["w_matrix_use"][:]
        if set(w_matrix_use[w_matrix_use!=0]) != set(range(1, max(w_matrix_use+1))):
            errors.append("W Matrix Value Error: Variable 'w_matrix_use' must index w matrix "
                          "use in H in numerical order from 1")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert Uncertainty Vector indices in uncertainty_vector_use are in numerical order from 1 --------------
        uncertainty_vector_use = rootgrp.variables["uncertainty_vector_use"][:]
        if set(uncertainty_vector_use[uncertainty_vector_use!=0]) != set(range(1, max(uncertainty_vector_use + 1))):
            errors.append("Uncertainty Vector Value Error: Variable 'uncertainty_vector_use' must index uncertainty "
                          "vector use in H in numerical order from 1")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert H columns assigned with W matrix also assign with Uncertainty Vector ----------------------------
        if not array_equal(where(w_matrix_use!=0)[0], where(uncertainty_vector_use!=0)[0]):
            errors.append("W Matrix Value Error: Mismatch between H columns with w matrix given by 'w_matrix_use' ("
                          + str(where(w_matrix_use!=0)[0]) + ") and H columns with an uncertainty vector given by "
                          "'uncertainty_vector_use' (" + str(where(uncertainty_vector_use!=0)[0])
                          + ") - should be the same")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert W matrix column widths equal to length of uncertainty vectors they are paired with --------------
        w_matrix_nnz = rootgrp.variables["w_matrix_nnz"][:]
        w_matrix_col = rootgrp.variables["w_matrix_col"][:]
        uncertainty_vector_row_count = rootgrp.variables["uncertainty_vector_row_count"][:]
        istart = 0
        iend = 0
        for i, w_matrix_nnz_i in enumerate(w_matrix_nnz):
            # Find highest column number in w matrix
            iend += w_matrix_nnz_i
            max_col = max(w_matrix_col[istart:iend]) + 1

            # Compare w matrix column width to size of uncertainty vectors paired with it
            i_w_u_vectors = uncertainty_vector_use[where(w_matrix_use == i+1)[0]] - 1
            for i_w_u_vector in i_w_u_vectors:
                max_row = uncertainty_vector_row_count[i_w_u_vector]
                if max_row != max_col:
                    errors.append("W Matrix Value Error: Number of columns of w matrix " + str(i+1) + " (" +
                                  str(max_col) + ") do no correspond to number of rows of uncertainty vector " +
                                  str(i_w_u_vector+1) + " (" + str(max_row) + ")")
            istart = iend
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert all uncertainty vector values > 0 ---------------------------------------------------------------
        uncertainty_vector = rootgrp.variables['uncertainty_vector'][:]
        if not (uncertainty_vector > 0).all():
            errors.append("W Matrix Value Error: Not all values of variable 'uncertainty_vector' > 0 - uncertainties "
                          "must be greater than zero")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert all index variables values >= 0 -----------------------------------------------------------------
        index_variables = ["w_matrix_nnz", "w_matrix_row", "w_matrix_col", "w_matrix_use",
                           "uncertainty_vector_row_count", "uncertainty_vector_use"]
        for variable in index_variables:
            if not (rootgrp.variables[variable][:] >= 0).all():
                errors.append("W Matrix Value Error: Not all values of variable '"+variable+"' >= 0 - indices "
                              "must be greater than or equal to zero")
        # --------------------------------------------------------------------------------------------------------------

    return errors


def check_uncertainty_assignment(rootgrp):
    """
    Return errors found with matchup dataset uncertainty type assignment

    :type rootgrp: netCDF4.Dataset
    :param rootgrp: In memory representation of matchup dataset

    :return:
        :error: *list:str:

        List of errors found during processing
    """

    errors = []  # Initialise list to store error messages

    # TEST: Assert assert 1 and only 1 error correlation form is attributed to each sensor state variable --------------
    if set(["H", "Us", "Ur"]).issubset(rootgrp.variables.keys()):

        m = rootgrp.dimensions['m'].size

        w_matrix_use = zeros(m)
        if "w_matrix_use" in rootgrp.variables.keys():
            w_matrix_use = rootgrp.variables['w_matrix_use']

        Ur = rootgrp.variables["Ur"][:]
        Us = rootgrp.variables["Us"][:]

        for col in range(m):
            rand = False
            randsys = False
            w = False
            if not isnan(rootgrp.variables['H'][0,col]):
                if (Ur[:, col] > 0).all() and (Us[:, col] == 0).all():
                    rand = True
                if (Ur[:, col] > 0).all() and (Us[:, col] > 0).all():
                    randsys = True
                if w_matrix_use[col] > 0:
                    w = True

                if (rand and (not randsys and not w)) \
                        or (randsys and (not rand and not w)) or (w and (not rand and not randsys)):
                    pass
                else:
                    errors.append("Variable Correlation Form Assignement Error: Variable in H[:, " + str(col) +
                                  "] correlation assigned as:"
                                  "\n - random correlation: " + str(rand) +
                                  "\n - random+systematic correlation: " + str(randsys) +
                                  "\n - w-matrix correlation: " + str(w) +
                                  "\n Must have one and only one form")
    # ------------------------------------------------------------------------------------------------------------------

    return errors


def main(path):
    """
    Run routine to check harmonisation input files match specification and are internally consistant

    :type path: str
    :param path: path of harmonisation input file to
    """

    print "Testing File:", path

    # 1. Find any errors!
    rootgrp = Dataset(path)

    errors = []
    errors += check_variable_included(rootgrp)
    errors += check_variable_dtype(rootgrp)
    errors += check_variable_dimension(rootgrp)
    errors += check_attributes(rootgrp)
    errors += check_dimension_sizes(rootgrp)
    errors += check_w_matrix_variable_dimensions(rootgrp)
    errors += check_variable_values(rootgrp)
    errors += check_w_variable_values(rootgrp)
    errors += check_uncertainty_assignment(rootgrp)

    rootgrp.close()

    # 2. Report errors

    if errors == []:
        print "Test Passed"
    else:
        print "The following errors have been detected:"
        for error in errors:
            print ">", error
    print ""

    return 0

if __name__ == "__main__":
    main(sys.argv[1])