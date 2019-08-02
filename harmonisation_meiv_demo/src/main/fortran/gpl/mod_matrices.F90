!> @file mod_matrices.F90
!! Error covariance matrices data structure.
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

!> @brief Module to provide the error covariance matrices data structure.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_matrices
  use mod_base

  implicit none

  private

  public matrices_initialize
  public matrices_finalize
  public add_b_matrix
  public add_d_matrix
  public add_v_matrix
  public add_v_matrices
  public add_w_matrix
  public add_w_matrices
  public add_z_matrix
  public get_capacity
  public get_b_matrix_count
  public get_d_matrix_count
  public get_v_matrix_count
  public get_w_matrix_count
  public get_z_matrix_count
  public has_b_matrix
  public has_d_matrix
  public has_v_matrix
  public has_w_matrix
  public has_z_matrix
  public remove_b_matrices
  public remove_d_matrices
  public remove_v_matrices
  public remove_w_matrices
  public remove_z_matrices
  !! The methods below are public for testing only
  public delete_matrix
  public new_banded_matrix
  public new_diagonal_matrix
  public new_perturbed_diagonal_matrix
  public new_rectangular_matrix
  public new_sparse_matrix

  interface add_b_matrix
    module procedure add_b_matrix_from_arrays
    module procedure add_b_matrix_from_perturbed_diagonal_matrix
  end interface

  interface add_d_matrix
    module procedure add_d_matrix_from_array
    module procedure add_d_matrix_from_diagonal_matrix
  end interface

  interface add_v_matrix
    module procedure add_v_matrix_from_array
    module procedure add_v_matrix_from_diagonal_matrix
  end interface

  interface add_v_matrices
    module procedure add_v_matrices_from_storage
  end interface

  interface add_w_matrix
    module procedure add_w_matrix_from_sparse_arrays
    module procedure add_w_matrix_from_banded_matrix
  end interface

  interface add_w_matrices
    module procedure add_w_matrices_from_storage
  end interface

  interface add_z_matrix
    module procedure add_z_matrix_from_array
    module procedure add_z_matrix_from_rectangular_matrix
  end interface

  interface delete_matrix
    module procedure delete_banded_matrix
    module procedure delete_diagonal_matrix
    module procedure delete_perturbed_diagonal_matrix
    module procedure delete_rectangular_matrix
    module procedure delete_sparse_matrix
  end interface

  interface new_banded_matrix
    module procedure new_banded_matrix_from_sparse_matrix
  end interface

  interface new_diagonal_matrix
    module procedure new_diagonal_matrix_from_value
    module procedure new_diagonal_matrix_from_array
  end interface

  interface new_perturbed_diagonal_matrix
    module procedure new_perturbed_diagonal_matrix_from_arrays
  end interface

  interface new_rectangular_matrix
    module procedure new_rectangular_matrix_from_array
  end interface

  interface new_sparse_matrix
    module procedure new_sparse_matrix_from_sparse_arrays
  end interface

  !> @brief The type to represent a banded matrix.
  type, public :: banded_matrix_t
    !> @brief The number of rows.
    integer :: m = 0
    !> @brief The number of columns.
    integer :: n = 0
    !> @brief The row indices.
    integer, allocatable :: ia(:)
    !> @brief The column indices.
    integer, allocatable :: ja(:)
    !> @brief The elements.
    real(kind=cp), allocatable :: a(:)
  end type banded_matrix_t

  !> @brief The type to represent a diagonal matrix.
  type, public :: diagonal_matrix_t
    !> @brief The number of rows.
    integer :: m = 0
    !> @brief The diagonal elements.
    real(kind=cp), allocatable :: a(:)
  end type diagonal_matrix_t

  !> @brief The type to represent a diagonal matrix that is perturbed by an outer vector product.
  type, public :: perturbed_diagonal_matrix_t
    !> @brief The number of rows.
    integer :: m = 0
    !> @brief The diagonal elements.
    real(kind=cp), allocatable :: a(:)
    !> @brief The elements of the perturbation vector.
    real(kind=cp), allocatable :: u(:)
  end type perturbed_diagonal_matrix_t

  !> @brief The type to represent a rectangular matrix.
  type, public :: rectangular_matrix_t
    !> @brief The number of rows.
    integer :: m = 0
    !> @brief The number of columns.
    integer :: n = 0
    !> @brief The matrix elements.
    real(kind=cp), allocatable :: a(:,:)
  end type rectangular_matrix_t

  !> @brief The type to represent a sparse matrix in compressed row storage (CRS) scheme.
  type, public :: sparse_matrix_t
    !> @brief The number of rows.
    integer :: m = 0
    !> @brief The number of columns.
    integer :: n = 0
    !> @brief The row indices of the sparse matrix.
    integer, allocatable :: ia(:)
    !> @brief The column indices associated with the non-zero elements.
    integer, allocatable :: ja(:)
    !> @brief The non-zero elements.
    real(kind=cp), allocatable :: a(:)
  end type sparse_matrix_t

  !> @brief The constant to designate independent random errors.
  integer, public, parameter :: C_INDEPENDENT = 1
  !> @brief The constant to designate independent random errors plus a common random perturbation.
  integer, public, parameter :: C_INDEPENDENT_PERTURBED = 2
  !> @brief The constant to designate structured random errors.
  integer, public, parameter :: C_STRUCTURED = 3
  !> @brief The constant to designate structured random errors plus a common random perturbation.
  !! @remark Not supported yet.
  integer, public, parameter :: C_STRUCTURED_PERTURBED = 4

  !> @brief The maximum number of matrices of each type.
  integer :: capacity = 0
  !> @brief The actual number of B matrices.
  integer :: nb = 0
  !> @brief The actual number of D matrices.
  integer :: nd = 0
  !> @brief The actual number of V matrices.
  integer :: nv = 0
  !> @brief The actual number of W matrices.
  integer :: nw = 0
  !> @brief The actual number of Z matrices.
  integer :: nz = 0
  !> @brief The B matrices (perturbed diagonal matrices).
  !! Each B matrix represents the error covariance matrix of a sensor telemetry data column of
  !! error correlation type @c C_INDEPENDENT_PERTURBED.
  type(perturbed_diagonal_matrix_t), allocatable, public, protected :: b_matrices(:)
  !> @brief The D matrices (diagonal matrices).
  !! Each D matrix represents the error covariance matrix of a sensor telemetry data column of
  !! error correlation type @c C_INDEPENDENT.
  type(diagonal_matrix_t), allocatable, public, protected :: d_matrices(:)
  !> @brief The V matrices (variance matrices).
  !! Each V matrix represents the error variance matrix associated with a sensor telemetry data
  !! column of error correlation type @c C_STRUCTURED.
  type(diagonal_matrix_t), allocatable, public, protected :: v_matrices(:)
  !> @brief The W matrices (weight matrices).
  !! Each W matrix represents the weight matrix associated with a sensor telemetry data
  !! column of error correlation type @c C_STRUCTURED or @c C_STRUCTURED_PERTURBED.
  type(banded_matrix_t), allocatable, public, protected :: w_matrices(:)
  !> @brief The Z matrices (variance matrices).
  !! Each column of each Z matrix represents the diagonal part of the error covariance matrix
  !! of a sensor telemetry data column.
  type(rectangular_matrix_t), allocatable, public, protected :: z_matrices(:)

contains

  !> @brief Initializes the module.
  !! @param[in] capacity_in The maximum number of matrices (of each type) that can be added.
  subroutine matrices_initialize( capacity_in )
    integer, intent(in) :: capacity_in

    call matrices_finalize
    capacity = capacity_in
    nb = 0
    nd = 0
    nv = 0
    nw = 0
    nz = 0
    allocate (b_matrices(capacity))
    allocate (d_matrices(capacity))
    allocate (v_matrices(capacity))
    allocate (w_matrices(capacity))
    allocate (z_matrices(capacity))
  end subroutine matrices_initialize

  !> @brief Finalizes the module and deallocates all memory used.
  subroutine matrices_finalize
    integer :: i

    call remove_z_matrices
    call remove_w_matrices
    call remove_v_matrices
    call remove_d_matrices
    call remove_b_matrices

    if (allocated( z_matrices )) then
      deallocate (z_matrices)
    end if
    if (allocated( w_matrices )) then
      deallocate (w_matrices)
    end if
    if (allocated( v_matrices )) then
      deallocate (v_matrices)
    end if
    if (allocated( d_matrices )) then
      deallocate (d_matrices)
    end if
    if (allocated( b_matrices )) then
      deallocate (b_matrices)
    end if

    capacity = 0
  end subroutine matrices_finalize

  subroutine add_b_matrix_from_perturbed_diagonal_matrix( matrix, matrix_index )
    type(perturbed_diagonal_matrix_t), intent(in) :: matrix
    integer, intent(out) :: matrix_index

    if (nb < capacity) then
      nb = nb + 1
      b_matrices(nb) = matrix
      matrix_index = nb
    else
      matrix_index = 0
    end if
  end subroutine add_b_matrix_from_perturbed_diagonal_matrix

  !> @brief Adds a new B matrix.
  !! @param[in] d The uncertainty data (independent random).
  !! @param[in] u The uncertainty data (common random).
  !! @param[out] matrix_indices The indices of the added B matrix.
  subroutine add_b_matrix_from_arrays( d, u, matrix_index )
    real(kind=cp), intent(in) :: d(:)
    real(kind=cp), intent(in) :: u(:)
    integer, intent(out) :: matrix_index

    call add_b_matrix_from_perturbed_diagonal_matrix( &
      new_perturbed_diagonal_matrix_from_arrays( size( d ), d**2, u ), matrix_index )
  end subroutine add_b_matrix_from_arrays

  subroutine add_d_matrix_from_diagonal_matrix( matrix, matrix_index )
    type(diagonal_matrix_t), intent(in)  :: matrix
    integer,                 intent(out) :: matrix_index

    if (nd < capacity) then
      nd = nd + 1
      d_matrices(nd) = matrix
      matrix_index = nd
    else
      matrix_index = 0
    end if
  end subroutine add_d_matrix_from_diagonal_matrix

  !> @brief Adds a new D matrix.
  !! @param[in] u The uncertainty data (independent random).
  !! @param[out] matrix_indices The indices of the added D matrix.
  subroutine add_d_matrix_from_array( u, matrix_index )
    real(kind=cp), intent(in)  :: u(:)
    integer, intent(out) :: matrix_index

    call add_d_matrix_from_diagonal_matrix( &
      new_diagonal_matrix_from_array( size( u ), u**2 ), matrix_index )
  end subroutine add_d_matrix_from_array

  subroutine add_v_matrix_from_diagonal_matrix( matrix, matrix_index )
    type(diagonal_matrix_t), intent(in) :: matrix
    integer, intent(out) :: matrix_index

    if (nv < capacity) then
      nv = nv + 1
      v_matrices(nv) = matrix
      matrix_index = nv
    else
      matrix_index = 0
    end if
  end subroutine add_v_matrix_from_diagonal_matrix

  !> @brief Adds a new V matrix.
  !! @param[in] u The uncertainty data (independent random).
  !! @param[out] matrix_indices The indices of the added V matrix.
  subroutine add_v_matrix_from_array( u, matrix_index )
    real(kind=cp), intent(in) :: u(:)
    integer, intent(out) :: matrix_index

    call add_v_matrix_from_diagonal_matrix( &
      new_diagonal_matrix_from_array( size( u ), u**2 ), matrix_index )
  end subroutine add_v_matrix_from_array

  !> @brief Adds new V matrices from raw storage.
  !! @param[in] matrix_count The number of V matrices included with the raw storage.
  !! @param[in] m The number of rows in each V matrix.
  !! @param[in] u The uncertainty data (independent random).
  !! @param[out] matrix_indices The indices of the added V matrices.
  subroutine add_v_matrices_from_storage( matrix_count, m, u, matrix_indices )
    integer, intent(in) :: matrix_count
    integer, intent(in) :: m(matrix_count)
    real(kind=cp), intent(in) :: u(:)
    integer, intent(out) :: matrix_indices(matrix_count)
    integer :: k
    integer :: m_begin
    integer :: m_final

    m_final = 0
    do k = 1, matrix_count
      m_begin = m_final + 1
      m_final = m_final + m(k)
      call add_v_matrix_from_array( u(m_begin:m_final), matrix_indices(k) )
    end do
  end subroutine add_v_matrices_from_storage

  subroutine add_w_matrix_from_banded_matrix( matrix, matrix_index )
    type(banded_matrix_t), intent(in) :: matrix
    integer, intent(out) :: matrix_index

    if (nw < capacity) then
      nw = nw + 1
      w_matrices(nw) = matrix
      matrix_index = nw
    else
      matrix_index = 0
    end if
  end subroutine add_w_matrix_from_banded_matrix

  subroutine add_w_matrix_from_sparse_arrays( m, s, ia, ja, a, matrix_index )
    integer, intent(in) :: m
    integer, intent(in) :: s
    integer, intent(in) :: ia(m + 1)
    integer, intent(in) :: ja(:)
    real(kind=cp), intent(in) :: a(:)
    integer, intent(out) :: matrix_index

    call add_w_matrix_from_banded_matrix( &
      new_banded_matrix_from_sparse_arrays( m, s, maxval( ja ), ia, ja, a ), matrix_index )
  end subroutine add_w_matrix_from_sparse_arrays

  !> @brief Adds new W matrices from raw storage.
  !! @param[in] matrix_count The number of W matrices included with the raw storage.
  !! @param[in] m The number of rows in each W matrix, which is the same for all W matrices.
  !! @param[in] s The subsetting step.
  !! @param[in] nnz The number of non-zero elements in each W matrix.
  !! @param[in,out] ia The row pointers.
  !! @param[in,out] ja The column indices of non-zero elements.
  !! @param[in] a The non-zero elements.
  !! @param[out] matrix_indices The indices of the added W matrices.
  !! @remark Row pointers @c ia and column indices @c ja are converted from 0-based array indices (C or
  !! Python format) into 1-based array indices (Fortran format), if necessary.
  subroutine add_w_matrices_from_storage( matrix_count, m, s, nnz, ia, ja, a, matrix_indices )
    integer, intent(in) :: matrix_count
    integer, intent(in) :: m
    integer, intent(in) :: s
    integer, intent(in) :: nnz(matrix_count)
    integer, intent(inout) :: ia(m + 1, matrix_count)
    integer, intent(inout) :: ja(:)
    real(kind=cp), intent(in) :: a(:)
    integer, intent(out) :: matrix_indices(matrix_count)

    integer :: k
    integer(kind=long) :: n_begin
    integer(kind=long) :: n_final

    call assert_fortran_crs_format( matrix_count, m, nnz, ia, ja )

    n_final = 0
    do k = 1, matrix_count
      n_begin = n_final + 1
      n_final = n_final + nnz(k)
      call add_w_matrix_from_sparse_arrays( m, s, ia(:, k), ja(n_begin:n_final), a(n_begin:n_final), matrix_indices(k) )
    end do
  end subroutine add_w_matrices_from_storage

  subroutine add_z_matrix_from_array( m, n, z, matrix_index )
    integer, intent(in) :: m
    integer, intent(in) :: n
    real(kind=cp), intent(in) :: z(m, n)
    integer, intent(out) :: matrix_index

    call add_z_matrix_from_rectangular_matrix( new_rectangular_matrix_from_array( m, n, z ), matrix_index )
  end subroutine add_z_matrix_from_array

  subroutine add_z_matrix_from_rectangular_matrix( matrix, matrix_index )
    type(rectangular_matrix_t), intent(in) :: matrix
    integer, intent(out) :: matrix_index

    if (nz < capacity) then
      nz = nz + 1
      z_matrices(nz) = matrix
      matrix_index = nz
    else
      matrix_index = 0
    end if
  end subroutine add_z_matrix_from_rectangular_matrix

  subroutine delete_banded_matrix( matrix )
    type(banded_matrix_t), intent(inout) :: matrix

    if (allocated( matrix%a )) then
      matrix%a = 0.0_cp
      deallocate (matrix%a)
    end if
    if (allocated( matrix%ia )) then
      matrix%ia = 0
      deallocate (matrix%ia)
    end if
    if (allocated( matrix%ja )) then
      matrix%ja = 0
      deallocate (matrix%ja)
    end if
    matrix%m = 0
    matrix%n = 0
  end subroutine delete_banded_matrix

  subroutine delete_diagonal_matrix( matrix )
    type(diagonal_matrix_t), intent(inout) :: matrix

    matrix%m = 0

    if (allocated( matrix%a )) then
      matrix%a = 0.0_cp
      deallocate (matrix%a)
    end if
  end subroutine delete_diagonal_matrix

  subroutine delete_perturbed_diagonal_matrix( matrix )
    type(perturbed_diagonal_matrix_t), intent(inout) :: matrix

    matrix%m = 0

    if (allocated( matrix%u )) then
      matrix%u = 0.0_cp
      deallocate (matrix%u)
    end if
    if (allocated( matrix%a )) then
      matrix%a = 0.0_cp
      deallocate (matrix%a)
    end if
  end subroutine delete_perturbed_diagonal_matrix

  subroutine delete_rectangular_matrix( matrix )
    type(rectangular_matrix_t), intent(inout) :: matrix

    matrix%m = 0
    matrix%n = 0

    if (allocated( matrix%a )) then
      matrix%a = 0.0_cp
      deallocate (matrix%a)
    end if
  end subroutine delete_rectangular_matrix

  subroutine delete_sparse_matrix( matrix )
    type(sparse_matrix_t), intent(inout) :: matrix

    if (allocated( matrix%a )) then
      matrix%a = 0.0_cp
      deallocate (matrix%a)
    end if
    if (allocated( matrix%ia )) then
      matrix%ia = 0
      deallocate (matrix%ia)
    end if
    if (allocated( matrix%ja )) then
      matrix%ja = 0
      deallocate (matrix%ja)
    end if
    matrix%m = 0
    matrix%n = 0
  end subroutine delete_sparse_matrix

  pure integer function get_capacity()
    get_capacity = capacity
  end function get_capacity

  pure integer function get_b_matrix_count()
    get_b_matrix_count = nb
  end function get_b_matrix_count

  pure integer function get_d_matrix_count()
    get_d_matrix_count = nd
  end function get_d_matrix_count

  pure integer function get_v_matrix_count()
    get_v_matrix_count = nv
  end function get_v_matrix_count

  pure integer function get_z_matrix_count()
    get_z_matrix_count = nz
  end function get_z_matrix_count

  pure integer function get_w_matrix_count()
    get_w_matrix_count = nw
  end function get_w_matrix_count

  pure logical function has_b_matrix( b_index )
    integer, intent(in) :: b_index

    has_b_matrix = b_index >= 1 .and. b_index <= nb
  end function has_b_matrix

  pure logical function has_d_matrix( d_index )
    integer, intent(in) :: d_index

    has_d_matrix = d_index >= 1 .and. d_index <= nd
  end function has_d_matrix

  pure logical function has_v_matrix( v_index )
    integer, intent(in) :: v_index

    has_v_matrix = v_index >= 1 .and. v_index <= nv
  end function has_v_matrix

  pure logical function has_w_matrix( w_index )
    integer, intent(in) :: w_index

    has_w_matrix = w_index >= 1 .and. w_index <= nw
  end function has_w_matrix

  pure logical function has_z_matrix( z_index )
    integer, intent(in) :: z_index

    has_z_matrix = z_index >= 1 .and. z_index <= nz
  end function has_z_matrix

  pure function new_banded_matrix_from_sparse_matrix( w ) result (banded_matrix)
    type(sparse_matrix_t), intent(in) :: w
    type(banded_matrix_t) :: banded_matrix

    banded_matrix = new_banded_matrix_from_sparse_arrays( w%m, 1, w%n, w%ia, w%ja, w%a )
  end function new_banded_matrix_from_sparse_matrix

  !> @brief Creates a new banded matrix from the arrays of a sparse matrix.
  !! @param[in] m The number of matrix rows.
  !! @param[in] s The subsetting step.
  !! @param[in] n The number of matrix columns.
  !! @param[in] ia The row pointers.
  !! @param[in] ja The column indices of non-zero elements.
  !! @param[in] a The non-zero elements.
  !! @return a new banded matrix.
  pure function new_banded_matrix_from_sparse_arrays( m, s, n, ia, ja, a ) result (banded_matrix)
    integer, intent(in) :: m
    integer, intent(in) :: s
    integer, intent(in) :: n
    integer, intent(in) :: ia(m + 1)
    integer, intent(in) :: ja(:)
    real(kind=cp), intent(in) :: a(:)
    type(banded_matrix_t) :: banded_matrix
    integer, allocatable :: ib(:)
    integer, allocatable :: jb(:)
    real(kind=cp), allocatable :: b(:)
    integer :: i
    integer :: k
    integer :: l

    !! Count number of banded elements
    k = 0
    l = 0
    do i = 1, m, s
      if (ia(i + 1) > ia(i)) then
        !! Invariant with respect to ascending or descending ordering of column indices
        k = k + abs( ja(ia(i + 1) - 1) - ja(ia(i)) ) + 1
      end if
      l = l + 1
    end do

    allocate (ib(l + 1))
    allocate (jb(l))
    allocate (b(k))
    jb = 0
    b = 0.0_cp

    !! Traverse banded elements
    k = 1
    l = 0
    do i = 1, m, s
      ib(l + 1) = k
      if (ia(i + 1) > ia(i)) then
        !! Invariance with respect to ascending or descending ordering of column indices
        jb(l + 1) = min( ja(ia(i)), ja(ia(i + 1) - 1) )
        b(ja(ia(i):ia(i + 1) - 1) - jb(l + 1) + k) = a(ia(i):ia(i + 1) - 1)
        b(ja(ia(i):ia(i + 1) - 1) - jb(l + 1) + k) = a(ia(i):ia(i + 1) - 1)
        k = k + abs( ja(ia(i + 1) - 1) - ja(ia(i)) ) + 1
      end if
      l = l + 1
    end do
    ib(l + 1) = k

    banded_matrix = banded_matrix_t( l, n, ib, jb, b )
  end function new_banded_matrix_from_sparse_arrays

  pure function new_diagonal_matrix_from_value( m, a ) result (diagonal_matrix)
    integer, intent(in) :: m
    real(kind=cp), intent(in) :: a
    type(diagonal_matrix_t) :: diagonal_matrix

    allocate (diagonal_matrix%a(m))
    diagonal_matrix%m = m
    diagonal_matrix%a = a
  end function new_diagonal_matrix_from_value

  pure function new_diagonal_matrix_from_array( m, a ) result (diagonal_matrix)
    integer, intent(in) :: m
    real(kind=cp), intent(in) :: a(m)
    type(diagonal_matrix_t) :: diagonal_matrix

    diagonal_matrix = diagonal_matrix_t( m, a )
  end function new_diagonal_matrix_from_array

  pure function new_perturbed_diagonal_matrix_from_arrays( m, a, v ) result (perturbed_diagonal_matrix)
    integer, intent(in) :: m
    real(kind=cp), intent(in) :: a(m)
    real(kind=cp), intent(in) :: v(m)
    type(perturbed_diagonal_matrix_t) :: perturbed_diagonal_matrix

    perturbed_diagonal_matrix = perturbed_diagonal_matrix_t( m, a, v )
  end function new_perturbed_diagonal_matrix_from_arrays

  pure function new_rectangular_matrix_from_array( m, n, a ) result (rectangular_matrix)
    integer, intent(in) :: m
    integer, intent(in) :: n
    real(kind=cp), intent(in) :: a(m, n)
    type(rectangular_matrix_t) :: rectangular_matrix

    rectangular_matrix = rectangular_matrix_t( m, n, a )
  end function new_rectangular_matrix_from_array

  pure function new_sparse_matrix_from_sparse_arrays( m, n, ia, ja, a ) result (sparse_matrix)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: ia(m + 1)
    integer, intent(in) :: ja(:)
    real(kind=cp), intent(in) :: a(:)
    type(sparse_matrix_t) :: sparse_matrix

    sparse_matrix = sparse_matrix_t( m, n, ia, ja, a )
  end function new_sparse_matrix_from_sparse_arrays

  subroutine remove_b_matrices
    integer :: i

    do i = nb, 1, -1
      call delete_matrix( b_matrices(i) )
    end do
    nb = 0
  end subroutine remove_b_matrices

  subroutine remove_d_matrices
    integer :: i

    do i = nd, 1, -1
      call delete_matrix( d_matrices(i) )
    end do
    nd = 0
  end subroutine remove_d_matrices

  subroutine remove_v_matrices
    integer :: i

    do i = nv, 1, -1
      call delete_matrix( v_matrices(i) )
    end do
    nv = 0
  end subroutine remove_v_matrices

  subroutine remove_w_matrices
    integer :: i

    do i = nw, 1, -1
      call delete_matrix( w_matrices(i) )
    end do
    nw = 0
  end subroutine remove_w_matrices

  subroutine remove_z_matrices
    integer :: i

    do i = nz, 1, -1
      call delete_matrix( z_matrices(i) )
    end do
    nz = 0
  end subroutine remove_z_matrices

  !> @brief Asserts that the input data comply with Fortran compressed row storage (CRS) format.
  !! @details Refer to [Barret et al. (1994, p. 57)](http://www.netlib.org/linalg/html_templates/Templates.html) for details.
  !!
  !! @param[in] matrix_count The number of matrices stored in the input data.
  !! @param[in] m The number of rows in each matrix stored, which is the same for all matrices.
  !! @param[in] nnz The number of non-zero elements in each matrix stored.
  !! @param[in,out] ia The row pointers.
  !! @param[in,out] ja The column indices of non-zero elements.
  !! @remark Row pointers @c ia and column indices @c ja are converted from 0-based array indices (C or
  !! Python format) into 1-based array indices (Fortran format), if necessary.
  subroutine assert_fortran_crs_format( matrix_count, m, nnz, ia, ja )
    integer, intent(in) :: matrix_count
    integer, intent(in) :: m
    integer, intent(in) :: nnz(matrix_count)
    integer, intent(inout) :: ia(m + 1, matrix_count)
    integer, intent(inout) :: ja(:)
    integer :: k
    integer :: n_begin
    integer :: n_final
    integer start_index

    select case (ia(m + 1, 1) - nnz(1))
      case (0)
        start_index = 0 !! we have 0-based array indices (C format)
      case (1)
        start_index = 1 !! we have 1-based array indices (Fortran format)
      case default
        error stop "Invalid CRS format" !! @todo handle error
    end select

    do k = 2, matrix_count
      if (ia(m + 1, k) - nnz(k) /= start_index) then
        error stop "Invalid CRS format" !! @todo handle error
      end if
    end do

    if (start_index == 0) then !! convert from C format into Fortan format
      ia = ia + 1
      ja = ja + 1
    end if
  end subroutine assert_fortran_crs_format

end module mod_matrices
