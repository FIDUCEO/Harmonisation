!> @file mod_reader.F90
!! Reading of matchup datasets.
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

!> @brief Module for reading matchup datasets.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_reader
  use mod_base
  use mod_assert
  use mod_config
  use mod_logger

  implicit none

  private

  public read_matchup_dataset

  interface read_matchup_dataset
    module procedure read_matchup_dataset_nc
  end interface

contains

  !> @brief Reads a netCDF matchup dataset.
  !! @param[in] path The path to the matchup dataset.
  !! @param[in] name The name of the matchup dataset.
  !! @param[in] i The index of the first sensor.
  !! @param[in] j The index of the second sensor.
  !! @param[in] s The subsetting step.
  subroutine read_matchup_dataset_nc( path, name, i, j, s )
    use mod_matrices, only: add_v_matrices
    use mod_matrices, only: add_w_matrices
    use mod_matchup
    use mod_sensors
    use mod_nc

    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: name
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: s
    !> @brief The telemetry data of the 1st sensor.
    real(kind=cp), allocatable :: x1(:,:)
    !> @brief The telemetry data of the 2nd sensor.
    real(kind=cp), allocatable :: x2(:,:)
    !> @brief The independent random uncertainties associated with @c x1.
    real(kind=cp), allocatable :: ur1(:,:)
    !> @brief The independent random uncertainties associated with @c x2.
    real(kind=cp), allocatable :: ur2(:,:)
    !> @brief The common random uncertainties associated with @c x1.
    real(kind=cp), allocatable :: us1(:,:)
    !> @brief The common random uncertainties associated with @c x2.
    real(kind=cp), allocatable :: us2(:,:)
    !> @brief The expected radiance differences.
    real(kind=cp), allocatable :: k(:)
    !> @brief The independent random uncertainties associated with @c k_.
    real(kind=cp), allocatable :: kr(:)
    !> @brief The common random uncertainties associated with @c k_.
    real(kind=cp), allocatable :: ks(:)
    !> @brief The time stamps for the 1st sensor.
    real(kind=tp), allocatable :: t1(:)
    !> @brief The time stamps for the 2nd sensor.
    real(kind=tp), allocatable :: t2(:)
    !> @brief The pixel positions for the 1st sensor.
    integer, allocatable :: p1(:)
    !> @brief The pixel positions for the 2nd sensor.
    integer, allocatable :: p2(:)
    !> @brief The uncertainty types associated with the telemetry data columns of the 1st sensor.
    integer, allocatable :: u_types1(:)
    !> @brief The uncertainty types associated with the telemetry data columns of the 2nd sensor.
    integer, allocatable :: u_types2(:)
    !> @brief The U-matrix associations for the telemetry data columns of the 1st sensor.
    integer, allocatable :: u_matrix_use1(:)
    !> @brief The U-matrix associations for the telemetry data columns of the 2nd sensor.
    integer, allocatable :: u_matrix_use2(:)
    !> @brief The W-matrix associations for the telemetry data columns of the 1st sensor.
    integer, allocatable :: w_matrix_use1(:)
    !> @brief The W-matrix associations for the telemetry data columns of the 2nd sensor.
    integer, allocatable :: w_matrix_use2(:)
    !> @brief The number of U matrices.
    integer :: u_matrix_count
    !> @brief The total number of rows in all U matrices.
    integer :: u_matrix_row_count_sum
    !> @brief The number of rows for each U matrix.
    integer, allocatable :: u_matrix_row_count(:)
    !> @brief The diagonal elements of all U matrices.
    real(kind=cp), allocatable :: u_matrix_val(:)
    !> @brief The number of W matrices.
    integer :: w_matrix_count
    !> @brief The total number of rows in all W matrices.
    integer :: w_matrix_row_count
    !> @brief The total number of non-zero elements in all W matrices.
    integer(kind=long) :: w_matrix_nnz_sum
    !> @brief The total number of non-zero elements in each W matrix.
    integer, allocatable :: w_matrix_nnz(:)
    !> @brief The row pointers for each W matrix.
    integer, allocatable :: w_matrix_row(:,:)
    !> @brief The non-zero elements of all W matrices.
    real(kind=cp), allocatable :: w_matrix_val(:)
    !> @brief The column indices of non-zero elements.
    integer, allocatable :: w_matrix_col(:)
    !> @brief The V matrix indices.
    integer, allocatable :: v_matrix_indices(:)
    !> @brief The W matrix indices.
    integer, allocatable :: w_matrix_indices(:)
    !> @brief Used internally.
    integer, allocatable :: ic(:)
    !> @brief Used internally.
    integer, allocatable :: jc(:)
    !> @brief Used internally.
    integer, allocatable :: iv(:)
    !> @brief Used internally.
    integer, allocatable :: jv(:)
    !> @brief Used internally.
    integer, allocatable :: iw(:)
    !> @brief Used internally.
    integer, allocatable :: jw(:)
    !> @brief The dataset ID.
    integer :: ncid
    !> @brief The matchup ID.
    integer :: m_id
    !> @brief The number of columns in @c x1.
    integer :: m1
    !> @brief The number of columns in @c x2.
    integer :: m2
    !> @brief The number of telemetry data columns used in the measurement equation of the 1st sensor.
    integer :: ri
    !> @brief The number of telemetry data columns used in the measurement equation of the 2nd sensor.
    integer :: rj
    !> @brief The number of matchup data records.
    integer :: m
    !> @brief The number of matchup data records after subsetting.
    integer :: m_act
    !> @brief The UUID of the matchup dataset.
    character(len=36) :: uuid

    !! 1) Open the dataset
    call log_info( "mod_reader", "reading dataset ", path )
    call open_dataset( path, ncid )
    call read_global_attribute( ncid, config_name_att_matchup_dataset_uuid, uuid, "" )

    !! 2) Read and check dimensions
    call read_dimension( ncid, "M", m )
    call read_dimension( ncid, "m1", m1 )
    call read_dimension( ncid, "m2", m2 )
    ri = get_sensor_r( i )
    rj = get_sensor_r( j )
    if (.not. config_job_enable_covariance_common) then
      call assert_true( m1 + m1 >= ri, "invalid dimension m1" )
      call assert_true( m2 + m2 >= rj, "invalid dimension m2" )
    else
      call assert_true( m1 >= ri, "invalid dimension m1" )
      call assert_true( m2 >= rj, "invalid dimension m2" )
    end if

    !! 3) Create a new matchup dataset
    m_act = m / s + min( 1, mod( m, s ) )
    call log_debug( "mod_reader", "adding matchup dataset: ", m_act )
    call add_matchup( name, i, j, m_act, ri, rj, m_id, uuid )
    if (m_id == 0) then
      call log_error( "mod_reader", "the maximum number of matchup datasets is exceeded: ", CONFIG_CAPACITY_MATCHUP )
      error stop "The maximum number of matchup datasets is exceeded"
    end if
    call log_debug( "mod_reader", "adding matchup dataset completed: ", m_id )

    !! 4.a) Allocate storage for sensor telemetry data
    if (.not. config_job_enable_covariance_common) then
      allocate (x1(m1 + m1, m))
      allocate (x2(m2 + m2, m))
    else
      allocate (x1(m1, m))
      allocate (x2(m2, m))
    end if
    !! 4.b) Read and check all sensor telemetry data
    call read_variable( ncid, "X1", x1(1:m1,:) )
    call read_variable( ncid, "X2", x2(1:m2,:) )
    call assert_not_nan( x1(1:m1,:), "X1" )
    call assert_not_nan( x2(1:m2,:), "X2" )
    call assert_not_infinite( x1(1:m1,:), "X1" )
    call assert_not_infinite( x2(1:m2,:), "X2" )
    call assert_not_huge( x1(1:m1,:), "X1" )
    call assert_not_huge( x2(1:m2,:), "X2" )
    if (.not. config_job_enable_covariance_common) then
      call read_variable( ncid, "Us1", x1(m1 + 1:m1 + m1,:), 0.0_cp )
      call read_variable( ncid, "Us2", x2(m2 + 1:m2 + m2,:), 0.0_cp )
      call assert_not_nan( x1(m1 + 1:m1 + m1,:), "Us1" )
      call assert_not_nan( x2(m2 + 1:m2 + m2,:), "Us2" )
      call assert_not_infinite( x1(m1 + 1:m1 + m1,:), "Us1" )
      call assert_not_infinite( x2(m2 + 1:m2 + m2,:), "Us2" )
      call assert_not_huge( x1(m1 + 1:m1 + m1,:), "Us1" )
      call assert_not_huge( x2(m2 + 1:m2 + m2,:), "Us2" )
    end if
    if (config_job_sensor_columns_select(i)) then
      call set_qi( m_id, transpose( x1(config_job_sensor_columns(i, 1:ri), 1:m:s) ) )
    else
      call set_qi( m_id, transpose( x1(1:ri, 1:m:s) ) )
    end if
    if (config_job_sensor_columns_select(j)) then
      call set_qj( m_id, transpose( x2(config_job_sensor_columns(j, 1:rj), 1:m:s) ) )
    else
      call set_qj( m_id, transpose( x2(1:rj, 1:m:s) ) )
    end if
    !! 4.c) Deallocate storage for sensor telemetry data
    deallocate (x2)
    deallocate (x1)

    !! 5.a) Allocate storage for sensor telemetry data uncertainties (independent and common random)
    allocate (ur1(m1, m))
    allocate (ur2(m2, m))
    allocate (us1(m1, m))
    allocate (us2(m2, m))
    !! 5.b) Read and check sensor telemetry uncertainties (independent and common random)
    call read_variable( ncid, "Ur1", ur1, 0.0_cp )
    call read_variable( ncid, "Ur2", ur2, 0.0_cp )
    call assert_not_nan( ur1, "Ur1" )
    call assert_not_nan( ur2, "Ur2" )
    call assert_not_infinite( ur1, "Ur1" )
    call assert_not_infinite( ur2, "Ur2" )
    call assert_not_huge( ur1, "Ur1" )
    call assert_not_huge( ur1, "Ur2" )
    call read_variable( ncid, "Us1", us1, 0.0_cp )
    call read_variable( ncid, "Us2", us2, 0.0_cp )
    call assert_not_nan( us1, "Us1" )
    call assert_not_nan( us2, "Us2" )
    call assert_not_infinite( us1, "Us1" )
    call assert_not_infinite( us2, "Us2" )
    call assert_not_huge( us1, "Us1" )
    call assert_not_huge( us1, "Us2" )
    !! 5.c) Allocate storage for matrix metadata
    allocate (u_types1(m1 + m1))
    allocate (u_types2(m2 + m2))
    u_types1 = 0
    u_types2 = 0
    allocate (u_matrix_use1(m1))
    allocate (u_matrix_use2(m2))
    allocate (w_matrix_use1(m1))
    allocate (w_matrix_use2(m2))
    !! 5.d) Read matrix metadata
    call read_variable( ncid, "uncertainty_type1", u_types1(1:m1), 1 )
    call read_variable( ncid, "uncertainty_type2", u_types2(1:m2), 1 )
    call read_dimension( ncid, "u_matrix_count", u_matrix_count, 0 )
    call read_dimension( ncid, "w_matrix_count", w_matrix_count, 0 )
    if (u_matrix_count > 0 .and. w_matrix_count > 0) then
      call read_dimension( ncid, "u_matrix_row_count_sum", u_matrix_row_count_sum )
      call read_dimension( ncid, "w_matrix_nnz_sum", w_matrix_nnz_sum )
      call read_dimension( ncid, "w_matrix_row_count", w_matrix_row_count )
      call read_variable( ncid, "u_matrix_use1", u_matrix_use1 )
      call read_variable( ncid, "u_matrix_use2", u_matrix_use2 )
      call read_variable( ncid, "w_matrix_use1", w_matrix_use1 )
      call read_variable( ncid, "w_matrix_use2", w_matrix_use2 )
      !! 5.e) Allocate storage for matrix indices
      allocate (v_matrix_indices(u_matrix_count))
      allocate (w_matrix_indices(w_matrix_count))
      v_matrix_indices = 0
      w_matrix_indices = 0
      !! 5.f) Read U matrices
      call log_debug( "mod_reader", "reading U matrices: ", w_matrix_count )
      allocate (u_matrix_row_count(u_matrix_count))
      allocate (u_matrix_val(u_matrix_row_count_sum))
      call read_variable( ncid, "u_matrix_row_count", u_matrix_row_count )
      call read_variable( ncid, "u_matrix_val", u_matrix_val )
      call assert_not_nan( u_matrix_val, "u_matrix_val" )
      call assert_not_infinite( u_matrix_val, "u_matrix_val" )
      call assert_not_huge( u_matrix_val, "u_matrix_val" )
      call assert_not_zero( u_matrix_val, "u_matrix_val" )
      call add_v_matrices( u_matrix_count, u_matrix_row_count, u_matrix_val, v_matrix_indices )
      if (any(v_matrix_indices == 0)) then
        call log_error( "mod_reader", "the maximum number of V matrices is exceeded ", CONFIG_CAPACITY_MATRICES )
        error stop "The maximum number of V matrices is exceeded"
      end if
      deallocate (u_matrix_val)
      deallocate (u_matrix_row_count)
      call log_debug( "mod_reader", "reading U matrices completed" )
      !! 5.g) Read W matrices
      call log_debug( "mod_reader", "reading W matrices: ", w_matrix_count )
      allocate (w_matrix_nnz(w_matrix_count))
      allocate (w_matrix_row(w_matrix_row_count, w_matrix_count))
      allocate (w_matrix_col(w_matrix_nnz_sum))
      allocate (w_matrix_val(w_matrix_nnz_sum))
      call read_variable( ncid, "w_matrix_nnz", w_matrix_nnz )
      call read_variable( ncid, "w_matrix_row", w_matrix_row )
      call read_variable( ncid, "w_matrix_col", w_matrix_col )
      call read_variable( ncid, "w_matrix_val", w_matrix_val )
      call assert_not_nan( w_matrix_val, "w_matrix_val" )
      call assert_not_infinite( w_matrix_val, "w_matrix_val" )
      call assert_not_huge( w_matrix_val, "w_matrix_val" )
      call assert_positive( w_matrix_val, "w_matrix_val" )
      call add_w_matrices( w_matrix_count, m, s, w_matrix_nnz, w_matrix_row, w_matrix_col, w_matrix_val, w_matrix_indices )
      if (any(w_matrix_indices == 0)) then
        call log_error( "mod_reader", "the maximum number of W matrices is exceeded ", CONFIG_CAPACITY_MATRICES )
        error stop "The maximum number of W matrices is exceeded"
      end if
      deallocate (w_matrix_val)
      deallocate (w_matrix_col)
      deallocate (w_matrix_row)
      deallocate (w_matrix_nnz)
      call log_debug( "mod_reader", "reading W matrices completed" )
    end if
    !! 5.h) Allocate storage for mapping of sensor telemetry data columns to matrices
    allocate (ic(ri))
    allocate (jc(rj))
    allocate (iv(ri))
    allocate (jv(rj))
    allocate (iw(ri))
    allocate (jw(rj))
    ic = 0
    jc = 0
    iv = 0
    jv = 0
    iw = 0
    jw = 0
    !! 5.i) Establish the mapping
    call setup( i, ri, m, s, ur1, us1, u_matrix_use1, w_matrix_use1, v_matrix_indices, w_matrix_indices, u_types1, ic, iv, iw )
    call setup( j, rj, m, s, ur2, us2, u_matrix_use2, w_matrix_use2, v_matrix_indices, w_matrix_indices, u_types2, jc, jv, jw )
    if (config_job_sensor_columns_select(i)) then
      call set_ui( m_id, u_types1(config_job_sensor_columns(i, 1:ri)) )
    else
      call set_ui( m_id, u_types1(1:ri) )
    end if
    if (config_job_sensor_columns_select(j)) then
      call set_uj( m_id, u_types2(config_job_sensor_columns(j, 1:rj)) )
    else
      call set_uj( m_id, u_types2(1:rj) )
    end if
    call set_ic( m_id, ic )
    call set_jc( m_id, jc )
    call set_iv( m_id, iv )
    call set_jv( m_id, jv )
    call set_iw( m_id, iw )
    call set_jw( m_id, jw )
    !! 5.j) Deallocate storage for sensor telemetry data uncertainties (independent and common random)
    deallocate (us2)
    deallocate (us1)
    deallocate (ur2)
    deallocate (ur1)

    !! 6) Read and check expected radiance differences and uncertainties
    allocate (k(m))
    allocate (kr(m))
    allocate (ks(m))
    call read_variable( ncid, "K", k, 0.0_cp )
    call read_variable( ncid, "Kr", kr, 1.0_cp )
    call read_variable( ncid, "Ks", ks, 0.0_cp )
    call assert_not_nan( k, "K" )
    call assert_not_nan( kr, "Kr" )
    call assert_not_nan( ks, "Ks" )
    call assert_not_infinite( k, "K" )
    call assert_not_infinite( kr, "Kr" )
    call assert_not_infinite( ks, "Ks" )
    call assert_not_huge( k, "K" )
    call assert_not_huge( kr, "Kr" )
    call assert_not_huge( ks, "Ks" )
    call assert_not_zero( Kr, "Kr" )
    call set_hj( m_id, k(1:m:s) )
    call set_vhj( m_id, addsq( kr(1:m:s), ks(1:m:s) ) )
    deallocate (ks)
    deallocate (kr)
    deallocate (k)

    !! 7) Read and check matchup times
    allocate (t1(m))
    allocate (t2(m))
    call read_variable( ncid, "time1", t1 )
    call read_variable( ncid, "time2", t2 )
    call assert_not_infinite( t1, "time1" )
    call assert_not_infinite( t2, "time2" )
    call assert_not_nan( t1, "time1" )
    call assert_not_nan( t2, "time2" )
    call set_ti( m_id, t1(1:m:s) )
    call set_tj( m_id, t2(1:m:s) )
    deallocate (t2)
    deallocate (t1)

    !! 8) Read and check pixel positions, if there are any
    if (has_variable( ncid, "across_track_index1")) then
      allocate (p1(m))
      call read_variable( ncid, "across_track_index1", p1 )
      call set_pi( m_id, p1(1:m:s) )
      deallocate (p1)
    end if
    if (has_variable( ncid, "across_track_index2")) then
      allocate (p2(m))
      call read_variable( ncid, "across_track_index2", p2 )
      call set_pj( m_id, p2(1:m:s) )
      deallocate (p2)
    end if

    !! 9) Close the dataset
    call close_dataset( ncid )
    call log_info( "mod_reader", "reading dataset completed" )
  end subroutine read_matchup_dataset_nc

  subroutine setup( i, r, m, s, ur, us, u_matrix_use, w_matrix_use, v_matrix_indices, w_matrix_indices, u_types, ic, iv, iw )
    use mod_matrices, only: C_INDEPENDENT
    use mod_matrices, only: C_INDEPENDENT_PERTURBED
    use mod_matrices, only: C_STRUCTURED
    use mod_matrices, only: C_STRUCTURED_PERTURBED
    use mod_matrices, only: add_b_matrix
    use mod_matrices, only: add_d_matrix

    integer, intent(in) :: i
    integer, intent(in) :: r
    integer, intent(in) :: m
    integer, intent(in) :: s
    real(kind=cp), intent(in) :: ur(:,:)
    real(kind=cp), intent(in) :: us(:,:)
    integer, intent(in) :: u_matrix_use(:)
    integer, intent(in) :: w_matrix_use(:)
    integer, intent(in) :: v_matrix_indices(:)
    integer, intent(in) :: w_matrix_indices(:)
    integer, intent(inout) :: u_types(:)
    integer, intent(out) :: ic(:)
    integer, intent(out) :: iv(:)
    integer, intent(out) :: iw(:)
    integer :: l

    do l = 1, r
      select case (u_types(l))
        case (C_INDEPENDENT)
          if (config_job_sensor_columns_select(i)) then
            call add_d_matrix( ur(config_job_sensor_columns(i, l), 1:m:s), ic(l) )
          else
            call add_d_matrix( ur(l, 1:m:s), ic(l) )
          end if
          if (ic(l) == 0) then
            call log_error( "mod_reader", "the maximum number of C matrices is exceeded ", CONFIG_CAPACITY_MATRICES )
            error stop "The maximum number of C matrices is exceeded"
          end if
        case (C_INDEPENDENT_PERTURBED)
          if (.not. config_job_enable_covariance_common) then
            u_types(l) = C_INDEPENDENT
            if (config_job_sensor_columns_select(i)) then
              call add_d_matrix( ur(config_job_sensor_columns(i, l), 1:m:s), ic(l) )
            else
              call add_d_matrix( ur(l, 1:m:s), ic(l) )
            end if
          else
            if (config_job_sensor_columns_select(i)) then
              call add_b_matrix( ur(config_job_sensor_columns(i, l), 1:m:s), us(config_job_sensor_columns(i, l), 1:m:s), ic(l) )
            else
              call add_b_matrix( ur(l, 1:m:s), us(l, 1:m:s), ic(l) )
            end if
          end if
          if (ic(l) == 0) then
            call log_error( "mod_reader", "the maximum number of C matrices is exceeded ", CONFIG_CAPACITY_MATRICES )
            error stop "The maximum number of C matrices is exceeded"
          end if
        case (C_STRUCTURED)
          if (config_job_sensor_columns_select(i)) then
            iv(l) = v_matrix_indices(u_matrix_use(config_job_sensor_columns(i, l)))
            iw(l) = w_matrix_indices(w_matrix_use(config_job_sensor_columns(i, l)))
          else
            iv(l) = v_matrix_indices(u_matrix_use(l))
            iw(l) = w_matrix_indices(w_matrix_use(l))
          end if
        case (C_STRUCTURED_PERTURBED) !! not yet supported, treated like structured random correlation
          u_types(l) = C_STRUCTURED
          if (config_job_sensor_columns_select(i)) then
            iv(l) = v_matrix_indices(u_matrix_use(config_job_sensor_columns(i, l)))
            iw(l) = w_matrix_indices(w_matrix_use(config_job_sensor_columns(i, l)))
          else
            iv(l) = v_matrix_indices(u_matrix_use(l))
            iw(l) = w_matrix_indices(w_matrix_use(l))
          end if
      end select
    end do
  end subroutine

  elemental function addsq( x, y )
    real(kind=cp), intent(in) :: x
    real(kind=cp), intent(in) :: y
    real(kind=cp) :: addsq

    addsq = real( x, dp )**2 + real( y, dp )**2 !! to compensate roundoff error
  end function addsq

end module mod_reader
