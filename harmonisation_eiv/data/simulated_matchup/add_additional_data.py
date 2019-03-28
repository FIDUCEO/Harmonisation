from netCDF4 import Dataset
from numpy import ones


def add_additional_data(path, variable_name):
    dataset = Dataset(path, "a")

    M = dataset.dimensions["M"].size

    additional_data = ones(M)

    additional_var = dataset.createVariable(variable_name, 'f8', ('M',), zlib=True, complevel=9)
    additional_var.description = ""
    additional_var[:] = additional_data[:]

    return 0


def main():
    add_additional_data("0_1.nc", "across_track_index1")
    add_additional_data("2_3.nc", "across_track_index1")


if __name__ == "__main__":
    main()
