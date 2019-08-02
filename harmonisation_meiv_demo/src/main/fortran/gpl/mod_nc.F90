!> @file mod_nc.F90
!! Reading and writing of netCDF datasets.
!! Copyright (C) 2019 FastOpt GmbH, Hamburg, Germany (info@fastopt.de)
!!
!! This code was developed for the EC project "Fidelity and Uncertainty in
!! Climate Data Records from Earth Observations (FIDUCEO)".
!! Grant Agreement: 638822
!!
!! This program is free software; you can redistribute it and/or modify it
!! under the terms of the GNU General Public License as published by the Free
!! Software Foundation; either version 3 of the License, or (at your option)
!! any later version.
!! This program is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
!! more details.
!!
!! A copy of the GNU General Public License should have been supplied along
!! with this program; if not, see http://www.gnu.org/licenses/
!!

!> @brief Module for reading and writing netCDF datasets.
!! @authors Ralf Quast, Ralf Giering
!! @copyright GNU Public License.
module mod_nc
  use mod_base, only: dp, sp, long

  implicit none

  private

  public add_attribute
  public add_global_attribute
  public add_dimension
  public add_variable
  public close_dataset
  public create_dataset
  public finish_dataset
  public has_global_attribute
  public has_dimension
  public has_variable
  public open_dataset
  public read_dimension
  public read_global_attribute
  public read_variable
  public write_variable

  interface add_attribute
    module procedure add_attribute_scalar_int
    module procedure add_attribute_scalar_real
    module procedure add_attribute_scalar_text
    module procedure add_attribute_vector_int
    module procedure add_attribute_vector_real
  end interface

  interface add_global_attribute
    module procedure add_global_attribute_scalar_int
    module procedure add_global_attribute_scalar_real
    module procedure add_global_attribute_scalar_text
    module procedure add_global_attribute_vector_int
    module procedure add_global_attribute_vector_real
  end interface

  interface add_dimension
    module procedure add_dimension_int
  end interface

  interface add_variable
    module procedure add_variable_matrix_int
    module procedure add_variable_matrix_real__dp
    module procedure add_variable_matrix_real__sp
    module procedure add_variable_scalar_int
    module procedure add_variable_scalar_real__dp
    module procedure add_variable_scalar_real__sp
    module procedure add_variable_scalar_text
    module procedure add_variable_vector_int
    module procedure add_variable_vector_real__dp
    module procedure add_variable_vector_real__sp
    module procedure add_variable_vector_text
  end interface

  interface read_dimension
    module procedure read_dimension_int
    module procedure read_dimension_int__long
    module procedure read_dimension_with_default_int
    module procedure read_dimension_with_default_int__long
  end interface

  interface read_global_attribute
    module procedure read_global_attribute_scalar_text
  end interface

  interface read_variable
    module procedure read_variable_matrix_int
    module procedure read_variable_matrix_real__dp
    module procedure read_variable_matrix_real__sp
    module procedure read_variable_matrix_with_default_real__dp
    module procedure read_variable_matrix_with_default_real__sp
    module procedure read_variable_scalar_int
    module procedure read_variable_scalar_real__dp
    module procedure read_variable_scalar_real__sp
    module procedure read_variable_vector_int
    module procedure read_variable_vector_real__dp
    module procedure read_variable_vector_real__sp
    module procedure read_variable_vector_with_default_int
    module procedure read_variable_vector_with_default_real__dp
    module procedure read_variable_vector_with_default_real__sp
  end interface

  interface write_variable
    module procedure write_variable_matrix_int
    module procedure write_variable_matrix_real__dp
    module procedure write_variable_matrix_real__sp
    module procedure write_variable_scalar_int
    module procedure write_variable_scalar_real__dp
    module procedure write_variable_scalar_real__sp
    module procedure write_variable_scalar_text
    module procedure write_variable_vector_int
    module procedure write_variable_vector_real__dp
    module procedure write_variable_vector_real__sp
    module procedure write_variable_vector_text
  end interface

  interface assert_dimension
    module procedure assert_dimension_int
    module procedure assert_dimension_int__long
  end interface

  integer, parameter :: SHUFFLE = 0
  integer, parameter :: DEFLATE = 1
  integer, parameter :: DEFLATE_LEVEL = 9

contains

  !> @brief Opens a netCDF dataset (for reading only).
  !! @param[in] path The path of the netCDF dataset.
  !! @param[out] ncid The ID of the netCDF dataset.
  subroutine open_dataset( path, ncid )
    use netcdf4_f03, only: nf_nowrite
    use netcdf4_f03, only: nf_open

    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid

    call assert_status( nf_open( trim( path ), nf_nowrite, ncid ), path )
  end subroutine open_dataset

  !> @brief Closes a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  subroutine close_dataset( ncid )
    use netcdf4_f03, only: nf_close

    integer, intent(in) :: ncid

    call assert_status( nf_close( ncid ), "close" )
  end subroutine close_dataset

  !> @brief Creates a new netCDF dataset.
  !! @param[in] path The path of the netCDF dataset.
  !! @param[out] ncid The ID of the netCDF dataset.
  subroutine create_dataset( path, ncid )
    use netcdf4_f03, only: nf_classic_model
    use netcdf4_f03, only: nf_create
    use netcdf4_f03, only: nf_netcdf4

    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer, parameter :: mode = ior( nf_classic_model, nf_netcdf4 )

    call assert_status( nf_create( trim( path ), mode, ncid ), path )
  end subroutine create_dataset

  !> @brief Finishes the definition of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  subroutine finish_dataset( ncid )
    use netcdf4_f03, only: nf_enddef

    integer, intent(in) :: ncid

    call assert_status( nf_enddef( ncid ), "finish" )
  end subroutine finish_dataset

  !> @brief Adds a dimension to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The dimension name.
  !! @param[in] name The length of the dimension.
  !! @param[out] dimid The ID of the dimension.
  subroutine add_dimension_int( ncid, name, nlen, dimid )
    use netcdf4_f03, only: nf_def_dim

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlen
    integer, intent(out) :: dimid

    call assert_status( nf_def_dim( ncid, trim( name ), nlen, dimid ), name )
  end subroutine add_dimension_int

  integer function my_inq_dimlen( ncid, dimid, llen )
    use netcdf4_nc_interfaces, only: nc_inq_dimlen

    integer, intent(in) :: ncid
    integer, intent(in) :: dimid
    integer(kind=long), intent(out) :: llen

    ! subtract 1 to obtain the C dimension ID
    my_inq_dimlen = nc_inq_dimlen( ncid, dimid - 1, llen )
  end function my_inq_dimlen

  !> @brief Adds a scalar variable (of double-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_scalar_real__dp( ncid, name, template, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_double

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: template
    integer, intent(out) :: varid

    integer, parameter :: xtype = nf_double
    integer, parameter :: ndims = 0
    integer, parameter :: dimids(0) = 0

    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, dimids, varid ), name )
  end subroutine add_variable_scalar_real__dp

  !> @brief Adds a scalar variable (of single-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_scalar_real__sp( ncid, name, template, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_float

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(in) :: template
    integer, intent(out) :: varid

    integer, parameter :: xtype = nf_float
    integer, parameter :: ndims = 0
    integer, parameter :: dimids(0) = 0

    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, dimids, varid ), name )
  end subroutine add_variable_scalar_real__sp

  !> @brief Adds a scalar variable (of integer type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_scalar_int( ncid, name, template, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_int

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: template
    integer, intent(out) :: varid

    integer, parameter :: xtype = nf_int
    integer, parameter :: ndims = 0
    integer, parameter :: dimids(0) = 0

    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, dimids, varid ), name )
  end subroutine add_variable_scalar_int

  !> @brief Adds a scalar variable (of character type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] lenid The ID of the dimension which defines the length of the character type.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_scalar_text( ncid, name, template, lenid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_char

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: template
    integer, intent(in) :: lenid
    integer, intent(out) :: varid

    integer, parameter :: xtype = nf_char
    integer, parameter :: ndims = 1

    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ lenid /), varid ), name )
  end subroutine add_variable_scalar_text

  !> @brief Adds a vector variable (of double-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] dimid The ID of the dimension which defines the length of the vector.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_vector_real__dp( ncid, name, template, dimid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_double

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: template(:)
    integer, intent(in ) :: dimid
    integer, intent(out) :: varid

    integer(kind=long) :: llen
    integer, parameter :: xtype = nf_double
    integer, parameter :: ndims = 1

    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )
    call assert_dimension( llen, size( template, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ dimid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_vector_real__dp

  !> @brief Adds a vector variable (of single-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] dimid The ID of the dimension which defines the length of the vector.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_vector_real__sp( ncid, name, template, dimid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_float

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(in) :: template(:)
    integer, intent(in ) :: dimid
    integer, intent(out) :: varid

    integer(kind=long) :: llen
    integer, parameter :: xtype = nf_float
    integer, parameter :: ndims = 1

    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )
    call assert_dimension( llen, size( template, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ dimid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_vector_real__sp

  !> @brief Adds a vector variable (of integer type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] dimid The ID of the dimension which defines the length of the vector.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_vector_int( ncid, name, template, dimid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_int

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: template(:)
    integer, intent(in ) :: dimid
    integer, intent(out) :: varid

    integer(kind=long) :: llen
    integer, parameter :: xtype = nf_int
    integer, parameter :: ndims = 1

    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )
    call assert_dimension( llen, size( template, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ dimid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_vector_int

  !> @brief Adds a vector variable (of character type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] lenid The ID of the dimension which defines the length of the character type.
  !! @param[in] dimid The ID of the dimension which defines the length of the vector.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_vector_text( ncid, name, template, lenid, dimid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_char

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: template(:)
    integer, intent(in ) :: lenid
    integer, intent(in ) :: dimid
    integer, intent(out) :: varid

    integer(kind=long) :: llen
    integer, parameter :: xtype = nf_char
    integer, parameter :: ndims = 2

    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )
    call assert_dimension( llen, size( template, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ lenid, dimid /), varid ), name )
  end subroutine add_variable_vector_text

  !> @brief Adds a matrix variable (of double-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] rowid The ID of the dimension which defines the number of matrix rows.
  !! @param[in] colid The ID of the dimension which defines the number of matrix columns.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_matrix_real__dp( ncid, name, template, rowid, colid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_double

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in ) :: template(:,:)
    integer, intent(in) :: rowid
    integer, intent(in) :: colid
    integer, intent(out) :: varid

    integer(kind=long) :: llen(2)
    integer, parameter :: xtype = nf_double
    integer, parameter :: ndims = 2

    call assert_status( my_inq_dimlen( ncid, rowid, llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, colid, llen(2) ), name )
    call assert_dimension( llen(1), size( template, 1, kind=long ), name )
    call assert_dimension( llen(2), size( template, 2, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ rowid, colid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_matrix_real__dp

  !> @brief Adds a matrix variable (of single-precision real type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] rowid The ID of the dimension which defines the number of matrix rows.
  !! @param[in] colid The ID of the dimension which defines the number of matrix columns.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_matrix_real__sp( ncid, name, template, rowid, colid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_float

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(in ) :: template(:,:)
    integer, intent(in) :: rowid
    integer, intent(in) :: colid
    integer, intent(out) :: varid

    integer(kind=long) :: llen(2)
    integer, parameter :: xtype = nf_float
    integer, parameter :: ndims = 2

    call assert_status( my_inq_dimlen( ncid, rowid, llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, colid, llen(2) ), name )
    call assert_dimension( llen(1), size( template, 1, kind=long ), name )
    call assert_dimension( llen(2), size( template, 2, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ rowid, colid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_matrix_real__sp

  !> @brief Adds a matrix variable (of integer type) to a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[in] template A template of the variable.
  !! @param[in] rowid The ID of the dimension which defines the number of matrix rows.
  !! @param[in] colid The ID of the dimension which defines the number of matrix columns.
  !! @param[out] varid The ID of the variable.
  subroutine add_variable_matrix_int( ncid, name, template, rowid, colid, varid )
    use netcdf4_f03, only: nf_def_var
    use netcdf4_f03, only: nf_def_var_deflate
    use netcdf4_f03, only: nf_int

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in ) :: template(:,:)
    integer, intent(in) :: rowid
    integer, intent(in) :: colid
    integer, intent(out) :: varid

    integer(kind=long) :: llen(2)
    integer, parameter :: xtype = nf_int
    integer, parameter :: ndims = 2

    call assert_status( my_inq_dimlen( ncid, rowid, llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, colid, llen(2) ), name )
    call assert_dimension( llen(1), size( template, 1, kind=long ), name )
    call assert_dimension( llen(2), size( template, 2, kind=long ), name )
    call assert_status( nf_def_var( ncid, trim( name ), xtype, ndims, (/ rowid, colid /), varid ), name )
    call assert_status( nf_def_var_deflate( ncid, varid, SHUFFLE, DEFLATE, DEFLATE_LEVEL ), name )
  end subroutine add_variable_matrix_int

  !> @brief Adds a scalar attribute (of real type) to a variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The ID of the variable.
  !! @param[in] name The attribute name.
  !! @param[in] value The attribute value.
  subroutine add_attribute_scalar_real( ncid, varid, name, value )
    use netcdf, only: nf90_put_att

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: value

    call assert_status( nf90_put_att( ncid, varid, trim( name ), value ), name )
  end subroutine add_attribute_scalar_real

  !> @brief Adds a scalar attribute (of integer type) to a variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The ID of the variable.
  !! @param[in] name The attribute name.
  !! @param[in] value The attribute value.
  subroutine add_attribute_scalar_int( ncid, varid, name, value )
    use netcdf, only: nf90_put_att

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in)  :: name
    integer, intent(in) :: value

    call assert_status( nf90_put_att( ncid, varid, trim( name ), value ), name )
  end subroutine add_attribute_scalar_int

  !> @brief Adds a scalar attribute (of character type) to a variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The ID of the variable.
  !! @param[in] name The attribute name.
  !! @param[in] text The attribute text.
  subroutine add_attribute_scalar_text( ncid, varid, name, text )
    use netcdf, only: nf90_put_att

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text

    call assert_status( nf90_put_att( ncid, varid, trim( name ), text ), name )
  end subroutine add_attribute_scalar_text

  !> @brief Adds a vector attribute (of real type) to a variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The ID of the variable.
  !! @param[in] name The attribute name.
  !! @param[in] values The attribute values.
  subroutine add_attribute_vector_real( ncid, varid, name, values )
    use netcdf, only: nf90_put_att

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: values(:)

    call assert_status( nf90_put_att( ncid, varid, trim( name ), values ), name )
  end subroutine add_attribute_vector_real

  !> @brief Adds a vector attribute (of integer) to a variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The ID of the variable.
  !! @param[in] name The attribute name.
  !! @param[in] values The attribute values.
  subroutine add_attribute_vector_int( ncid, varid, name, values )
    use netcdf, only: nf90_put_att

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: name
    integer, intent(in) :: values(:)

    call assert_status( nf90_put_att( ncid, varid, trim( name ), values ), name )
  end subroutine add_attribute_vector_int

  !> @brief Adds a scalar attribute (of real type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[in] value The attribute value.
  subroutine add_global_attribute_scalar_real( ncid, name, value )
    use netcdf, only: nf90_put_att, nf90_global

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: value

    call assert_status( nf90_put_att( ncid, nf90_global, trim( name ), value ), name )
  end subroutine add_global_attribute_scalar_real

  !> @brief Adds a scalar attribute (of integer type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[in] value The attribute value.
  subroutine add_global_attribute_scalar_int( ncid, name, value )
    use netcdf, only: nf90_put_att, nf90_global

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: value

    call assert_status( nf90_put_att( ncid, nf90_global, trim( name ), value ), name )
  end subroutine add_global_attribute_scalar_int

  !> @brief Adds a scalar attribute (of character type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[in] text The attribute text.
  subroutine add_global_attribute_scalar_text( ncid, name, text )
    use netcdf, only: nf90_put_att, nf90_global

    integer, intent(in)  :: ncid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text

    call assert_status( nf90_put_att( ncid, nf90_global, trim( name ), text ), name )
  end subroutine add_global_attribute_scalar_text

  !> @brief Adds a vector attribute (of real type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[in] values The attribute values.
  subroutine add_global_attribute_vector_real( ncid, name, values )
    use netcdf, only: nf90_put_att, nf90_global

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: values(:)

    call assert_status( nf90_put_att( ncid, nf90_global, trim( name ), values ), name )
  end subroutine add_global_attribute_vector_real

  !> @brief Adds a vector attribute (of integer type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[in] values The attribute values.
  subroutine add_global_attribute_vector_int( ncid, name, values )
    use netcdf, only: nf90_put_att, nf90_global

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: values(:)

    call assert_status( nf90_put_att( ncid, nf90_global, trim( name ), values ), name )
  end subroutine add_global_attribute_vector_int

  !> @brief Tests if a netCDF dataset has a certain global attribute.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @return @c .true. if the dataset has the global attribute, @c .false. otherwise.
  logical function has_global_attribute( ncid, name )
    use netcdf, only: nf90_global
    use netcdf, only: nf90_inquire_attribute
    use netcdf4_f03, only: nf_noerr

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: dimid
    integer :: status

    status = nf90_inquire_attribute( ncid, nf90_global, trim( name ) )
    has_global_attribute = (status == nf_noerr)
  end function has_global_attribute

  !> @brief Tests if a netCDF dataset has a certain dimension.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The dimension name.
  !! @return @c .true. if the dataset has the dimension, @c .false. otherwise.
  logical function has_dimension( ncid, name )
    use netcdf4_f03, only: nf_inq_dimid
    use netcdf4_f03, only: nf_noerr

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: dimid
    integer :: status

    status = nf_inq_dimid( ncid, trim( name ), dimid )
    has_dimension = (status == nf_noerr)
  end function has_dimension

  !> @brief Tests if a netCDF dataset has a certain variable.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @return @c .true. if the dataset has the variable, @c .false. otherwise.
  logical function has_variable( ncid, name )
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_noerr

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: varid
    integer :: status

    status = nf_inq_varid( ncid, trim( name ), varid )
    has_variable = (status == nf_noerr)
  end function has_variable

  !> @brief Reads the length of a dimension (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The dimension name.
  !! @param[out] nlen The length of the dimension.
  subroutine read_dimension_int( ncid, name, nlen )
    use netcdf4_f03, only: nf_inq_dimid

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: nlen
    integer :: dimid
    integer(kind=long) :: llen

    call assert_status( nf_inq_dimid( ncid, trim( name ), dimid ), name )
    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )

    nlen = int( llen )
  end subroutine read_dimension_int

  !> @brief Reads the length of a dimension (of long integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The dimension name.
  !! @param[out] llen The length of the dimension.
  subroutine read_dimension_int__long( ncid, name, llen )
    use netcdf4_f03, only: nf_inq_dimid

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(kind=long), intent(out) :: llen
    integer :: dimid

    call assert_status( nf_inq_dimid( ncid, trim( name ), dimid ), name )
    call assert_status( my_inq_dimlen( ncid, dimid, llen ), name )
  end subroutine read_dimension_int__long

  subroutine read_dimension_with_default_int( ncid, name, nlen, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: nlen
    integer, intent(in) :: default_value

    if (has_dimension( ncid, name )) then
      call read_dimension_int( ncid, name, nlen )
    else
      nlen = default_value
    end if
  end subroutine read_dimension_with_default_int

  subroutine read_dimension_with_default_int__long( ncid, name, llen, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(kind=long), intent(out) :: llen
    integer(kind=long), intent(in) :: default_value

    if (has_dimension( ncid, name )) then
      call read_dimension_int__long( ncid, name, llen )
    else
      llen = default_value
    end if
  end subroutine read_dimension_with_default_int__long

  !> @brief Reads a scalar attribute (of character type) to the global attributes of a netCDF dataset.
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The attribute name.
  !! @param[out] value The attribute value.
  !! @param[in] default_value The default attribute value to return, if the attribute is not present.
  subroutine read_global_attribute_scalar_text( ncid, name, value, default_value )
    use netcdf, only: nf90_get_att, nf90_global

    integer, intent(in)  :: ncid
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: value
    character(len=*), intent(in), optional :: default_value

    if (.not. has_global_attribute( ncid, name ) .and. present(default_value)) then
      value = default_value
    else
      call assert_status( nf90_get_att( ncid, nf90_global, trim( name ), value ), name )
    end if
  end subroutine read_global_attribute_scalar_text

  !> @brief Reads the datum of a scalar variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] datum The variable datum.
  subroutine read_variable_scalar_real__dp( ncid, name, datum )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: datum
    integer :: varid
    integer :: ndims

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_get_var( ncid, varid, datum ), name )
  end subroutine read_variable_scalar_real__dp

  !> @brief Reads the datum of a scalar variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] datum The variable datum.
  subroutine read_variable_scalar_real__sp( ncid, name, datum )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: datum
    integer :: varid
    integer :: ndims

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_get_var( ncid, varid, datum ), name )
  end subroutine read_variable_scalar_real__sp

  !> @brief Reads the datum of a scalar variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] datum The variable datum.
  subroutine read_variable_scalar_int( ncid, name, datum )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: datum
    integer :: varid
    integer :: ndims

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_get_var( ncid, varid, datum ), name )
  end subroutine read_variable_scalar_int

  !> @brief Reads the data of a vector variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_vector_real__dp( ncid, name, ddata )
    use netcdf4_f03, only: nf_float
    use netcdf4_f03, only: nf_get_var_double
    use netcdf4_f03, only: nf_get_var_real
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_inq_vartype

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: ddata(:)
    integer :: varid
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen
    integer :: xtype
    real(kind=sp), allocatable :: rdata(:)

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( ddata, kind=long ), name )
    call assert_status( nf_inq_vartype( ncid, varid, xtype ), name )

    if (xtype == nf_float) then !! convert data
      allocate (rdata(llen))
      call assert_status( nf_get_var_real( ncid, varid, rdata ), name )
      ddata = real( rdata, dp )
      deallocate (rdata)
    else
      call assert_status( nf_get_var_double( ncid, varid, ddata ), name )
    end if
  end subroutine read_variable_vector_real__dp

  !> @brief Reads the data of a vector variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_vector_real__sp( ncid, name, rdata )
    use netcdf4_f03, only: nf_double
    use netcdf4_f03, only: nf_get_var_double
    use netcdf4_f03, only: nf_get_var_real
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_inq_vartype

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: rdata(:)
    integer :: varid
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen
    integer :: xtype
    real(kind=dp), allocatable :: ddata(:)

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( rdata, kind=long ), name )
    call assert_status( nf_inq_vartype( ncid, varid, xtype ), name )

    if (xtype == nf_double) then !! convert data
      allocate (ddata(llen))
      call assert_status( nf_get_var_double( ncid, varid, ddata ), name )
      rdata = real( ddata, sp )
      deallocate (ddata)
    else
      call assert_status( nf_get_var_real( ncid, varid, rdata ), name )
    end if
  end subroutine read_variable_vector_real__sp

  !> @brief Reads the data of a vector variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_vector_int( ncid, name, data )
    use netcdf4_f03, only: nf_get_var_int
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: data(:)
    integer :: varid
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen

    call assert_status( nf_inq_varid( ncid, trim( name ), varid ), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( data, kind=long ), name )
    call assert_status( nf_get_var_int( ncid, varid, data ), name )
  end subroutine read_variable_vector_int

  subroutine read_variable_vector_with_default_real__dp( ncid, name, data, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: data(:)
    real(kind=dp), intent(in) :: default_value

    if (has_variable( ncid, name )) then
      call read_variable_vector_real__dp( ncid, name, data )
    else
      data = default_value
    end if
  end subroutine read_variable_vector_with_default_real__dp

  subroutine read_variable_vector_with_default_real__sp( ncid, name, data, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: data(:)
    real(kind=sp), intent(in) :: default_value

    if (has_variable( ncid, name )) then
      call read_variable_vector_real__sp( ncid, name, data )
    else
      data = default_value
    end if
  end subroutine read_variable_vector_with_default_real__sp

  subroutine read_variable_vector_with_default_int( ncid, name, data, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: data(:)
    integer, intent(in) :: default_value

    if (has_variable( ncid, name )) then
      call read_variable_vector_int( ncid, name, data )
    else
      data = default_value
    end if
  end subroutine read_variable_vector_with_default_int

  !> @brief Reads the data of a matrix variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_matrix_real__dp( ncid, name, data )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: data(:,:)
    integer :: varid
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varid( ncid, trim( name ), varid), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_get_var( ncid, varid, data ), name )
  end subroutine read_variable_matrix_real__dp

  !> @brief Reads the data of a matrix variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_matrix_real__sp( ncid, name, data )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: data(:,:)
    integer :: varid
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varid( ncid, trim( name ), varid), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_get_var( ncid, varid, data ), name )
  end subroutine read_variable_matrix_real__sp

  !> @brief Reads the data of a matrix variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] name The variable name.
  !! @param[out] data The variable data.
  subroutine read_variable_matrix_int( ncid, name, data )
    use netcdf, only: nf90_get_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varid
    use netcdf4_f03, only: nf_inq_varndims

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: data(:,:)
    integer :: varid
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varid( ncid, trim( name ), varid), name )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_get_var( ncid, varid, data ), name )
  end subroutine read_variable_matrix_int

  subroutine read_variable_matrix_with_default_real__dp( ncid, name, data, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: data(:,:)
    real(kind=dp), intent(in) :: default_value

    if (has_variable( ncid, name )) then
      call read_variable_matrix_real__dp( ncid, name, data )
    else
      data = default_value
    end if
  end subroutine read_variable_matrix_with_default_real__dp

  subroutine read_variable_matrix_with_default_real__sp( ncid, name, data, default_value )
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: data(:,:)
    real(kind=sp), intent(in) :: default_value

    if (has_variable( ncid, name )) then
      call read_variable_matrix_real__sp( ncid, name, data )
    else
      data = default_value
    end if
  end subroutine read_variable_matrix_with_default_real__sp

  !> @brief Writes the datum of a scalar variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] datum The variable datum.
  subroutine write_variable_scalar_real__dp( ncid, varid, datum )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=dp), intent(in) :: datum
    character(len=nf_max_name) :: name
    integer :: ndims

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_put_var( ncid, varid, datum ), name )
  end subroutine write_variable_scalar_real__dp

  !> @brief Writes the datum of a scalar variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] datum The variable datum.
  subroutine write_variable_scalar_real__sp( ncid, varid, datum )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=sp), intent(in) :: datum
    character(len=nf_max_name) :: name
    integer :: ndims

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_put_var( ncid, varid, datum ), name )
  end subroutine write_variable_scalar_real__sp

  !> @brief Writes the datum of a scalar variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] datum The variable datum.
  subroutine write_variable_scalar_int( ncid, varid, datum )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: datum
    character(len=nf_max_name) :: name
    integer :: ndims

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 0, ndims, name )
    call assert_status( nf90_put_var( ncid, varid, datum ), name )
  end subroutine write_variable_scalar_int

  !> @brief Writes the data of a scalar variable (of character type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable datum.
  subroutine write_variable_scalar_text( ncid, varid, data )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: data
    character(len=nf_max_name) :: name
    integer :: ndims

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf90_put_var( ncid, varid, data ), name )
  end subroutine write_variable_scalar_text

  !> @brief Writes the data of a vector variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_vector_real__dp( ncid, varid, data )
    use netcdf4_f03, only: nf_put_var_double
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=dp), intent(in) :: data(:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( data, kind=long ), name )
    call assert_status( nf_put_var_double( ncid, varid, data ), name )
  end subroutine write_variable_vector_real__dp

  !> @brief Writes the data of a vector variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_vector_real__sp( ncid, varid, data )
    use netcdf4_f03, only: nf_put_var_real
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=sp), intent(in) :: data(:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( data, kind=long ), name )
    call assert_status( nf_put_var_real( ncid, varid, data ), name )
  end subroutine write_variable_vector_real__sp

  !> @brief Writes the data of a vector variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_vector_int( ncid, varid, data )
    use netcdf4_f03, only: nf_put_var_int
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: data(:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(1)
    integer(kind=long) :: llen

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 1, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen ), name )
    call assert_dimension( llen, size( data, kind=long ), name )
    call assert_status( nf_put_var_int( ncid, varid, data ), name )
  end subroutine write_variable_vector_int

  !> @brief Writes the data of a vector variable (of character real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_vector_text( ncid, varid, data )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character(len=*), intent(in) :: data(:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen ), name )
    call assert_dimension( llen, size( data, kind=long ), name )
    call assert_status( nf90_put_var( ncid, varid, data ), name )
  end subroutine write_variable_vector_text

  !> @brief Writes the data of a matrix variable (of double-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_matrix_real__dp( ncid, varid, data )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=dp), intent(in) :: data(:,:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_put_var( ncid, varid, data ), name )
  end subroutine write_variable_matrix_real__dp

  !> @brief Writes the data of a vector variable (of single-precision real type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_matrix_real__sp( ncid, varid, data )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(kind=sp), intent(in) :: data(:,:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_put_var( ncid, varid, data ), name )
  end subroutine write_variable_matrix_real__sp

  !> @brief Writes the data of a vector variable (of integer type).
  !! @param[in] ncid The ID of the netCDF dataset.
  !! @param[in] varid The variable ID.
  !! @param[in] data The variable data.
  subroutine write_variable_matrix_int( ncid, varid, data )
    use netcdf, only: nf90_put_var
    use netcdf4_f03, only: nf_inq_vardimid
    use netcdf4_f03, only: nf_inq_varname
    use netcdf4_f03, only: nf_inq_varndims
    use netcdf4_f03, only: nf_max_name

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: data(:,:)
    character(len=nf_max_name) :: name
    integer :: ndims
    integer :: dimids(2)
    integer(kind=long) :: llen(2)

    call assert_status( nf_inq_varname( ncid, varid, name ), "unknown variable" )
    call assert_status( nf_inq_varndims( ncid, varid, ndims ), name )
    call assert_dimension( 2, ndims, name )
    call assert_status( nf_inq_vardimid( ncid, varid, dimids ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(1), llen(1) ), name )
    call assert_status( my_inq_dimlen( ncid, dimids(2), llen(2) ), name )
    call assert_dimension( llen(1), size( data, 1, kind=long ), name )
    call assert_dimension( llen(2), size( data, 2, kind=long ), name )
    call assert_status( nf90_put_var( ncid, varid, data ), name )
  end subroutine write_variable_matrix_int

  subroutine assert_dimension_int( expected, actual, what )
    integer, intent (in) :: expected
    integer, intent (in) :: actual
    character(len=*), intent (in) :: what

    if (expected /= actual) then
      call handle_error( message="expected and actual dimensions do not match", what=what )
    end if
  end subroutine assert_dimension_int

  subroutine assert_dimension_int__long( expected, actual, what )
    integer(kind=long), intent (in) :: expected
    integer(kind=long), intent (in) :: actual
    character(len=*), intent (in) :: what

    if (expected /= actual) then
      call handle_error( message="expected and actual dimensions do not match", what=what )
    end if
  end subroutine assert_dimension_int__long

  subroutine assert_status( status, what )
    use netcdf4_f03, only: nf_noerr

    integer, intent (in) :: status
    character(len=*), intent (in) :: what

    if (status /= nf_noerr) then
      call handle_error( status=status, what=what )
    end if
  end subroutine assert_status

  !> @brief Handles an error.
  !! @param[in] status The status returned by a netCDF API function (optional).
  !! @param[in] message The error message (optional).
  !! @param[in] what The cause of the error.
  !! @bug The affected file is not closed explicitly.
  subroutine handle_error( status, message, what )
    use netcdf4_f03, only: nf_strerror
    use mod_error_handler, only: handle => handle_error

    integer, intent (in), optional :: status
    character(len=*), intent (in), optional :: message
    character(len=*), intent (in), optional :: what

    if (present( status )) then
      if (present( what )) then
        call handle( status, nf_strerror( status ), what )
      else
        call handle( status, nf_strerror( status ) )
      end if
    else if (present( message )) then
      if (present( what )) then
        call handle( message=message, what=what )
      else
        call handle( message=message )
      end if
    end if
  end subroutine handle_error

end module mod_nc
