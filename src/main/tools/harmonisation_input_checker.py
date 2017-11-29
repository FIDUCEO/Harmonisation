"""
Check harmonisation input files match specification and are internally consistant

Usage:
python harmonisation_input_checker.py path/to/matchup/file.nc
"""

'''___Python Modules___'''
import sys

'''___Third Party Modules___'''
from netCDF4 import Dataset
from numpy import sum, where, array_equal, zeros, append

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
VARIABLE_DATA = {"X1": {"dim": ["M", "m1"], "dtype": ["float32", "float64"]},
                 "X2": {"dim": ["M", "m2"], "dtype": ["float32", "float64"]},
                 "Ur1": {"dim": ["M", "m1"], "dtype": ["float32", "float64"]},
                 "Ur2": {"dim": ["M", "m2"], "dtype": ["float32", "float64"]},
                 "Us1": {"dim": ["M", "m1"], "dtype": ["float32", "float64"]},
                 "Us2": {"dim": ["M", "m2"], "dtype": ["float32", "float64"]},
                 "K": {"dim": ["M"], "dtype": ["float32", "float64"]},
                 "Kr": {"dim": ["M"], "dtype": ["float32", "float64"]},
                 "Ks": {"dim": ["M"], "dtype": ["float32", "float64"]},
                 "time1": {"dim": ["M"], "dtype": ["float32", "float64"]},
                 "time2": {"dim": ["M"], "dtype": ["float32", "float64"]},
                 "uncertainty_type1": {"dim": ["m1"], "dtype": ["int32"]},
                 "uncertainty_type2": {"dim": ["m2"], "dtype": ["int32"]}}

# Optional W matrix variable dimensions (if included complete set required)
W_VARIABLE_DATA = {"w_matrix_nnz": {"dim": ["w_matrix_count"], "dtype": ["int32"]},
                   "w_matrix_row": {"dim": ['w_matrix_count', 'w_matrix_row_count'], "dtype": ["int32"]},
                   "w_matrix_col": {"dim": ["w_matrix_nnz_sum"], "dtype": ["int32"]},
                   "w_matrix_val": {"dim": ["w_matrix_nnz_sum"], "dtype": ["float32", "float64"]},
                   "w_matrix_use1": {"dim": ["m1"], "dtype": ["int32"]},
                   "w_matrix_use2": {"dim": ["m2"], "dtype": ["int32"]},
                   "u_matrix_row_count": {"dim": ["u_matrix_count"], "dtype": ["int32"]},
                   "u_matrix_val": {"dim": ["u_matrix_row_count_sum"], "dtype": ["float32", "float64"]},
                   "u_matrix_use1": {"dim": ["m1"], "dtype": ["int32"]},
                   "u_matrix_use2": {"dim": ["m2"], "dtype": ["int32"]}}

# Attributes
ATTRS = ["sensor_1_name", "sensor_2_name"]


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

        if (expected_dims != []) and (test_dims not in expected_dims):
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

    # TEST: Assert sensor_1_name and sensor_2_name different -----------------------------------------------------------
    if ("sensor_1_name" in test_attrs) and ("sensor_2_name" in test_attrs):
        if rootgrp.sensor_1_name == rootgrp.sensor_2_name:
            errors.append("Attribute Error: Attributes sensor_1_name & sensor_2_name must have different values")
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
        use_num_ws = max(max(rootgrp.variables["w_matrix_use1"][:]),max(rootgrp.variables["w_matrix_use2"][:]))
        if dim_num_ws != use_num_ws:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_count' ("+str(dim_num_ws) +
                          ") must match number of labelled W matrices in variable 'w_matrix_use'(" + str(use_num_ws) +
                          ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'u_matrix_count' is equal to the number of u matrices ---------------------
        #       indexed in variable 'u_matrix_use1' and 'u_matrix_use2' ------------------------------------------------
        dim_num_u_vecs = rootgrp.dimensions["u_matrix_count"].size
        use_num_u_vecs = max(max(rootgrp.variables["u_matrix_use1"][:]),
                             max(rootgrp.variables["u_matrix_use2"][:]))
        if dim_num_u_vecs != use_num_u_vecs:
            errors.append("U Matrix Dimension Error: Size of dimension 'u_matrix_count' ("
                          + str(dim_num_u_vecs) + ") must match number of labelled u matrices in variable"
                          " 'u_matrix_use'(" + str(use_num_u_vecs) + ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert that the row size of w matrices is equal to the number of matchups ------------------------------
        w_matrix_row_count_size = rootgrp.dimensions["w_matrix_row_count"].size
        M_size = rootgrp.dimensions["M"].size
        if w_matrix_row_count_size != M_size+1:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_row_count_num' (" +
                          str(w_matrix_row_count_size) + ") must equal number of matchups + 1 (" + str(M_size) +
                          " + 1 = "+str(M_size+1)+")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'w_matrix_nnz_sum' is equal to the sum of each w matrix nnz combined, as --
        #       given in the variable 'w_matrix_nnz' -------------------------------------------------------------------
        w_matrix_nnz_sum_size = rootgrp.dimensions["w_matrix_nnz_sum"].size
        w_matrix_nnz_combined = sum(rootgrp.variables["w_matrix_nnz"][:])
        if w_matrix_nnz_sum_size != w_matrix_nnz_combined:
            errors.append("W Matrix Dimension Error: Size of dimension 'w_matrix_nnz_sum' (" +
                          str(w_matrix_nnz_sum_size) + ") must equal combined per w matrix nnz's contained in variable"
                          "'w_matrix_nnz' (" + str(w_matrix_nnz_combined) + ")")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert the size of dimension 'u_matrix_row_count_sum' is equal to the sum of each u matrix number of ---
        #       rows combined, as given in the variable 'u_matrix_row_count' -------------------------------------------
        u_matrix_row_count_sum_size = rootgrp.dimensions["u_matrix_row_count_sum"].size
        u_matrix_row_count_combined = sum(rootgrp.variables["u_matrix_row_count"][:])
        if u_matrix_row_count_sum_size != u_matrix_row_count_combined:
            errors.append("U Matrix Dimension Error: Size of dimension 'u_matrix_row_count_sum' ("
                          + str(u_matrix_row_count_sum_size) + ") must equal combined per u matrix"
                          " row counts contained in variable 'u_matrix_row_count' (" +
                          str(u_matrix_row_count_combined) + ")")
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

    # TEST: Assert no negative uncertainties in Ur1, Ur2, Us1 and Us2 --------------------------------------------------
    test_variables = ["Ur1", "Ur2", "Us1", "Ur2"]
    if set(test_variables).issubset(rootgrp.variables.keys()):
        for test_variable in test_variables:
            test_variable_val = rootgrp.variables[test_variable][:]
            if not (test_variable_val.flatten() >= 0).all():
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
        w_matrix_use = append(rootgrp.variables["w_matrix_use1"][:], rootgrp.variables["w_matrix_use2"][:])
        u_matrix_use = append(rootgrp.variables["u_matrix_use1"][:], rootgrp.variables["u_matrix_use2"][:])

        # TEST: Assert W Matrix indices in w_matrix_use are in numerical order from 1 ----------------------------------
        if set(w_matrix_use[w_matrix_use!=0]) != set(range(1, max(w_matrix_use+1))):
            errors.append("W Matrix Value Error: Variable 'w_matrix_use' must index w matrix "
                          "use in X1 then X2 in numerical order from 1")
        # --------------------------------------------------------------------------------------------------------------
        # TEST: Assert U Matrix indices in u_matrix_use are in numerical order from 1 ----------------------------------

        if set(u_matrix_use[u_matrix_use!=0]) != set(range(1, max(u_matrix_use + 1))):
            errors.append("U Matrix Value Error: Variable 'u_matrix_use' must index u matrix "
                          "use in X1 then X2 in numerical order from 1")
        # --------------------------------------------------------------------------------------------------------------
        # TEST: Assert X1 and X2 columns assigned with W matrix also assign with U Matrix ------------------------------
        if not array_equal(where(w_matrix_use!=0)[0], where(u_matrix_use!=0)[0]):
            errors.append("W Matrix Value Error: Mismatch between X1/X2 columns with w matrix given by 'w_matrix_use' ("
                          + str(where(w_matrix_use!=0)[0]) + ") and X1/X2 columns with an u matrix given by "
                          "'u_matrix_use' (" + str(where(u_matrix_use!=0)[0])
                          + ") - should be the same")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert W matrix column widths equal to rows of u matrix they are paired with ---------------------------
        w_matrix_nnz = rootgrp.variables["w_matrix_nnz"][:]
        w_matrix_col = rootgrp.variables["w_matrix_col"][:]
        u_matrix_row_count = rootgrp.variables["u_matrix_row_count"][:]
        istart = 0
        iend = 0
        for i, w_matrix_nnz_i in enumerate(w_matrix_nnz):
            # Find highest column number in w matrix
            iend += w_matrix_nnz_i
            max_col = max(w_matrix_col[istart:iend]) + 1

            # Compare w matrix column width number of rows in u matrix paired with it
            i_w_u_vectors = u_matrix_use[where(w_matrix_use == i+1)[0]] - 1
            for i_w_u_vector in i_w_u_vectors:
                max_row = u_matrix_row_count[i_w_u_vector]
                if max_row != max_col:
                    errors.append("W Matrix Value Error: Number of columns of w matrix " + str(i+1) + " (" +
                                  str(max_col) + ") do not correspond to number of rows of u matrix " +
                                  str(i_w_u_vector+1) + " (" + str(max_row) + ")")
            istart = iend
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert all U Matrix values > 0 ---------------------------------------------------------------
        u_matrix_val = rootgrp.variables['u_matrix_val'][:]
        if not (u_matrix_val > 0).all():
            errors.append("U Matrix Value Error: Not all values of variable 'u_matrix' > 0 - uncertainties "
                          "must be greater than zero")
        # --------------------------------------------------------------------------------------------------------------

        # TEST: Assert all index variables values >= 0 -----------------------------------------------------------------
        index_variables = ["w_matrix_nnz", "w_matrix_row", "w_matrix_col", "u_matrix_row_count"]
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
    if set(["X1", "X2", "Ur1", "Ur2", "Us1", "Us2", "uncertainty_type1", "uncertainty_type2"])\
            .issubset(rootgrp.variables.keys()):

        m1 = rootgrp.dimensions['m1'].size
        m2 = rootgrp.dimensions['m2'].size
        uncertainty_type1 = rootgrp.variables["uncertainty_type1"][:]
        uncertainty_type2 = rootgrp.variables["uncertainty_type2"][:]

        w_matrix_use1 = zeros(m1)
        if "w_matrix_use1" in rootgrp.variables.keys():
            w_matrix_use1 = rootgrp.variables['w_matrix_use1']

        w_matrix_use2 = zeros(m2)
        if "w_matrix_use2" in rootgrp.variables.keys():
            w_matrix_use2 = rootgrp.variables['w_matrix_use2']

        # 1. Sensor 1
        Ur1 = rootgrp.variables["Ur1"][:]
        Us1 = rootgrp.variables["Us1"][:]

        for col in range(m1):

            # (a) Find uncertainty_type asserted in file variable
            if uncertainty_type1[col] == 1:
                uncertainty_str = "independent"
            if uncertainty_type1[col] == 2:
                uncertainty_str = "independent+systematic"
            if uncertainty_type1[col] == 3:
                uncertainty_str = "structured"
            if uncertainty_type1[col] not in set([1, 2, 3]):
                uncertainty_str = "[NOT ASSIGNED!]"
                "Value Error: All values of 'uncertainty_type' variable must be equal to either 1, 2 or 3"

            # (b) Find uncertainty_type implied by the data
            rand = False
            randsys = False
            w = False
            if (Ur1[:, col] > 0).all() and (Us1[:, col] == 0).all():
                rand = True
            if (Ur1[:, col] > 0).all() and (Us1[:, col] > 0).all():
                randsys = True
            if w_matrix_use1[col] > 0:
                w = True

            # Compare how (a) and (b) match:
            if uncertainty_type1[col] == 1 and rand is True:
                pass
            elif uncertainty_type1[col] == 2 and randsys is True:
                pass
            elif uncertainty_type1[col] == 3 and w is True:
                pass
            else:
                errors.append("Variable Correlation Form Assignement Error: Variable in X1[:, " + str(col) +
                              "] correlation assigned as "+uncertainty_str+" but data for: "
                              "\n - independent error correlation: " + str(rand) +
                              "\n - independent+systematic correlation: " + str(randsys) +
                              "\n - structured correlation: " + str(w))

            if (rand and (not randsys and not w)) \
                    or (randsys and (not rand and not w)) or (w and (not rand and not randsys)):
                pass
            else:
                errors.append("Variable Correlation Form Assignement Error: Variable in X1[:, " + str(col) +
                              "] correlation assigned as "+uncertainty_str+" but data for: "
                              "\n - random correlation: " + str(rand) +
                              "\n - random+systematic correlation: " + str(randsys) +
                              "\n - w-matrix correlation: " + str(w) +
                              "\n Must have one and only one form!")

        # 2. Sensor 2
        Ur2 = rootgrp.variables["Ur2"][:]
        Us2 = rootgrp.variables["Us2"][:]

        for col in range(m2):

            # (a) Find uncertainty_type asserted in file variable
            if uncertainty_type2[col] == 1:
                uncertainty_str = "independent"
            if uncertainty_type2[col] == 2:
                uncertainty_str = "independent+systematic"
            if uncertainty_type2[col] == 3:
                uncertainty_str = "structured"
            if uncertainty_type2[col] not in set([1, 2, 3]):
                "Value Error: All values of 'uncertainty_type' variable must be equal to either 1, 2 or 3"

            # (b) Find uncertainty_type implied by the data
            rand = False
            randsys = False
            w = False
            if (Ur2[:, col] > 0).all() and (Us2[:, col] == 0).all():
                rand = True
            if (Ur2[:, col] > 0).all() and (Us2[:, col] > 0).all():
                randsys = True
            if w_matrix_use2[col] > 0:
                w = True

            # Compare how (a) and (b) match:
            if uncertainty_type2[col] == 1 and rand is True:
                pass
            elif uncertainty_type2[col] == 2 and randsys is True:
                pass
            elif uncertainty_type2[col] == 3 and w is True:
                pass
            else:
                errors.append("Variable Correlation Form Assignement Error: Variable in X2[:, " + str(col) +
                              "] correlation assigned as " + uncertainty_str + " but data for: "
                              "\n - independent error correlation: " + str(rand) +
                              "\n - independent+systematic correlation: " + str(randsys) +
                              "\n - structured correlation: " + str(w))

            if (rand and (not randsys and not w)) \
                    or (randsys and (not rand and not w)) or (w and (not rand and not randsys)):
                pass
            else:
                errors.append("Variable Correlation Form Assignement Error: Variable in X2[:, " + str(col) +
                              "] correlation assigned as " + uncertainty_str + " but data for: "
                              "\n - random correlation: " + str(rand) +
                              "\n - random+systematic correlation: " + str(randsys) +
                              "\n - w-matrix correlation: " + str(w) +
                              "\n Must have one and only one form!")

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
