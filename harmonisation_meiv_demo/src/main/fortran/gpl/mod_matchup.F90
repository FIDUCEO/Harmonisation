!> @file mod_matchup.F90
!! Matchup data structure.
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

!> @brief Module providing the matchup data structure.
!! @author Ralf Quast
!! @copyright GNU Public License.
!! @pre Requires proper initialization and set up of the @c mod_sensors module.
module mod_matchup
  use mod_type

  implicit none

  private

  public matchup_initialize
  public matchup_finalize
  public add_matchup
  public get_matchup_capacity
  public get_matchup_count
  public get_matchup_name
  public get_i
  public get_j
  public get_m
  public get_ri
  public get_rj
  public get_ti
  public get_tj
  public set_tj
  public set_ti
  public set_qi
  public set_qj
  public set_hi
  public set_hj
  public set_vhi
  public set_vhj
  public set_pi
  public set_pj
  public set_ui
  public set_uj
  public set_ic
  public set_jc
  public set_iv
  public set_jv
  public set_iw
  public set_jw
  public set_iz
  public set_jz

  !> @brief The maximum number of matchup datasets.
  integer, public, protected :: capacity = 0
  !> @brief The actual number of matchup datasets.
  integer, public, protected :: nmatchup = 0
  !> @brief The matchup datasets.
  type(matchup_t), allocatable, public, protected :: matchups(:)

  !> @brief The start indices into the total array of matchups.
  integer, allocatable, public, protected :: mstart(:)
  !> @brief The final indices into the total array of matchups.
  integer, allocatable, public, protected :: mfinal(:)

contains

  subroutine matchup_initialize( capacity_in )
    integer, intent(in) :: capacity_in

    call matchup_finalize
    capacity = capacity_in
    nmatchup = 0
    allocate (matchups(capacity))
    allocate (mstart(capacity))
    allocate (mfinal(capacity))
  end subroutine matchup_initialize

  subroutine matchup_finalize
    integer :: i

    do i = nmatchup, 1, -1
      call delete_matchup( matchups(i) )
    end do
    if (allocated( mfinal )) then
      deallocate (mfinal)
    end if
    if (allocated( mstart )) then
      deallocate (mstart)
    end if
    if (allocated( matchups )) then
      deallocate (matchups)
    end if

    capacity = 0
    nmatchup = 0
  end subroutine matchup_finalize

  subroutine delete_matchup( matchup )
    type(matchup_t), intent(inout) :: matchup

    matchup%name = ""
    matchup%i = 0
    matchup%j = 0
    matchup%m = 0
    matchup%ri = 0
    matchup%rj = 0

    if (allocated( matchup%qi )) then
      matchup%qi = 0.0_cp
      deallocate (matchup%qi)
    end if
    if (allocated( matchup%qj )) then
      matchup%qj = 0.0_cp
      deallocate (matchup%qj)
    end if
    if (allocated( matchup%hi )) then
      matchup%hi = 0.0_cp
      deallocate (matchup%hi)
    end if
    if (allocated( matchup%hj )) then
      matchup%hj = 0.0_cp
      deallocate (matchup%hj)
    end if
    if (allocated( matchup%vhi )) then
      matchup%vhi = 0.0_cp
      deallocate (matchup%vhi)
    end if
    if (allocated( matchup%vhj )) then
      matchup%vhj = 0.0_cp
      deallocate (matchup%vhj)
    end if
    if (allocated( matchup%ti )) then
      matchup%ti = 0.0_tp
      deallocate (matchup%ti)
    end if
    if (allocated( matchup%tj )) then
      matchup%tj = 0.0_tp
      deallocate (matchup%tj)
    end if
    if (allocated( matchup%pi )) then
      matchup%pi = 0
      deallocate (matchup%pi)
    end if
    if (allocated( matchup%pj )) then
      matchup%pj = 0
      deallocate (matchup%pj)
    end if
    if (allocated( matchup%ui )) then
      matchup%ui = 0
      deallocate (matchup%ui)
    end if
    if (allocated( matchup%uj )) then
      matchup%uj = 0
      deallocate (matchup%uj)
    end if
    if (allocated( matchup%ic )) then
      matchup%ic = 0
      deallocate (matchup%ic)
    end if
    if (allocated( matchup%jc )) then
      matchup%jc = 0
      deallocate (matchup%jc)
    end if
    if (allocated( matchup%iv )) then
      matchup%iv = 0
      deallocate (matchup%iv)
    end if
    if (allocated( matchup%jv )) then
      matchup%jv = 0
      deallocate (matchup%jv)
    end if
    if (allocated( matchup%iw )) then
      matchup%iw = 0
      deallocate (matchup%iw)
    end if
    if (allocated( matchup%jw )) then
      matchup%jw = 0
      deallocate (matchup%jw)
    end if

    matchup%iz = 0
    matchup%jz = 0
  end subroutine delete_matchup

  subroutine add_matchup( name, i, j, m, ri, rj, matchup_index, uuid )
    character(len=*), intent(in) :: name
    integer, intent(in)  :: i
    integer, intent(in)  :: j
    integer, intent(in)  :: m
    integer, intent(in)  :: ri
    integer, intent(in)  :: rj
    integer, intent(out) :: matchup_index
    character(len=36), intent(in), optional :: uuid

    associate (k => nmatchup)
      if (k < capacity .and. is_pair( i, j ) .and. .not. has_pair( i, j )) then
        k = k + 1

        matchups(k)%name = trim( name )
        if (present( uuid )) then
          matchups(k)%uuid = uuid
        end if
        matchups(k)%i = i
        matchups(k)%j = j
        matchups(k)%m = m
        matchups(k)%ri = ri
        matchups(k)%rj = rj

        if (k > 1) then
          mstart(k) = mfinal(k - 1) + 1
          mfinal(k) = mfinal(k - 1) + m
        else
          mstart(k) = 1
          mfinal(k) = m
        end if

        allocate (matchups(k)%qi(m, ri))
        matchups(k)%qi = 0.0_cp

        allocate (matchups(k)%qj(m, rj))
        matchups(k)%qj = 0.0_cp

        allocate (matchups(k)%hi(m))
        matchups(k)%hi = 0.0_cp

        allocate (matchups(k)%hj(m))
        matchups(k)%hj = 0.0_cp

        allocate (matchups(k)%vhi(m))
        matchups(k)%vhi = 0.0_cp

        allocate (matchups(k)%vhj(m))
        matchups(k)%vhj = 0.0_cp

        allocate (matchups(k)%ti(m))
        matchups(k)%ti = 0.0_tp

        allocate (matchups(k)%tj(m))
        matchups(k)%tj = 0.0_tp

        allocate (matchups(k)%pi(m))
        matchups(k)%pi = 0

        allocate (matchups(k)%pj(m))
        matchups(k)%pj = 0

        allocate (matchups(k)%ui(ri))
        matchups(k)%ui = 0

        allocate (matchups(k)%uj(rj))
        matchups(k)%uj = 0

        allocate (matchups(k)%ic(ri))
        matchups(k)%ic = 0

        allocate (matchups(k)%jc(rj))
        matchups(k)%jc = 0

        allocate (matchups(k)%iv(ri))
        matchups(k)%iv = 0

        allocate (matchups(k)%jv(rj))
        matchups(k)%jv = 0

        allocate (matchups(k)%iw(ri))
        matchups(k)%iw = 0

        allocate (matchups(k)%jw(rj))
        matchups(k)%jw = 0

        matchups(k)%iz = 0
        matchups(k)%jz = 0

        matchup_index = k
      else
        matchup_index = 0
      end if
    end associate
  end subroutine add_matchup

  pure integer function get_matchup_capacity()
    get_matchup_capacity = capacity
  end function get_matchup_capacity

  pure integer function get_matchup_count()
    get_matchup_count = nmatchup
  end function get_matchup_count

  function get_matchup_name( k ) result (name)
    integer, intent(in) :: k
    character(len=TYPE_MAX_LEN_NAME) :: name

    call assert_matchup( k )
    name = matchups(k)%name
  end function get_matchup_name

  integer function get_i( k )
    integer, intent(in) :: k

    call assert_matchup( k )
    get_i = matchups(k)%i
  end function get_i

  integer function get_j( k )
    integer, intent(in) :: k

    call assert_matchup( k )
    get_j = matchups(k)%j
  end function get_j

  integer function get_m( k )
    integer, intent(in) :: k

    call assert_matchup( k )
    get_m = matchups(k)%m
  end function get_m

  integer function get_ri( k )
    integer, intent(in) :: k

    call assert_matchup( k )
    get_ri = matchups(k)%ri
  end function get_ri

  integer function get_rj( k )
    integer, intent(in) :: k

    call assert_matchup( k )
    get_rj = matchups(k)%rj
  end function get_rj

  subroutine get_ti( k, t )
    integer, intent(in) :: k
    real(kind=tp), intent(out) :: t(:)

    call assert_matchup( k )
    t = matchups(k)%ti
  end subroutine get_ti

  subroutine set_ti( k, t )
    integer,       intent(in) :: k
    real(kind=tp), intent(in) :: t(:)

    call assert_matchup( k )
    matchups(k)%ti = t
  end subroutine set_ti

  subroutine get_tj( k, t )
    integer,  intent(in) :: k
    real(kind=tp), intent(out) :: t(:)

    call assert_matchup( k )
    t = matchups(k)%tj
  end subroutine get_tj

  subroutine set_tj( k, t )
    integer,       intent(in) :: k
    real(kind=tp), intent(in) :: t(:)

    call assert_matchup( k )
    matchups(k)%tj = t
  end subroutine set_tj

  subroutine set_qi( k, q )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: q(:,:)

    call assert_matchup( k )
    matchups(k)%qi = q
  end subroutine set_qi

  subroutine set_qj( k, q )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: q(:,:)

    call assert_matchup( k )
    matchups(k)%qj = q
  end subroutine set_qj

  subroutine set_hi( k, h )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: h(:)

    call assert_matchup( k )
    matchups(k)%hi = h
  end subroutine set_hi

  subroutine set_hj( k, h )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: h(:)

    call assert_matchup( k )
    matchups(k)%hj = h
  end subroutine set_hj

  subroutine set_vhi( k, v )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: v(:)

    call assert_matchup( k )
    matchups(k)%vhi = v
  end subroutine set_vhi

  subroutine set_vhj( k, v )
    integer,       intent(in) :: k
    real(kind=cp), intent(in) :: v(:)

    call assert_matchup( k )
    matchups(k)%vhj = v
  end subroutine set_vhj

  subroutine set_pi( k, p )
    integer, intent(in) :: k
    integer, intent(in) :: p(:)

    call assert_matchup( k )
    matchups(k)%pi = p
  end subroutine set_pi

  subroutine set_pj( k, p )
    integer, intent(in) :: k
    integer, intent(in) :: p(:)

    call assert_matchup( k )
    matchups(k)%pj = p
  end subroutine set_pj

  subroutine set_ui( k, uncertainty_types )
    integer, intent(in) :: k
    integer, intent(in) :: uncertainty_types(:)

    call assert_matchup( k )
    matchups(k)%ui = uncertainty_types
  end subroutine set_ui

  subroutine set_uj( k, uncertainty_types )
    integer, intent(in) :: k
    integer, intent(in) :: uncertainty_types(:)

    call assert_matchup( k )
    matchups(k)%uj = uncertainty_types
  end subroutine set_uj

  subroutine set_ic( k, c_indices )
    integer, intent(in) :: k
    integer, intent(in) :: c_indices(:)

    call assert_matchup( k )
    matchups(k)%ic = c_indices
  end subroutine set_ic

  subroutine set_jc( k, c_indices )
    integer, intent(in) :: k
    integer, intent(in) :: c_indices(:)

    call assert_matchup( k )
    matchups(k)%jc = c_indices
  end subroutine set_jc

  subroutine set_iv( k, v_indices )
    integer, intent(in) :: k
    integer, intent(in) :: v_indices(:)

    call assert_matchup( k )
    matchups(k)%iv = v_indices
  end subroutine set_iv

  subroutine set_jv( k, v_indices )
    integer, intent(in) :: k
    integer, intent(in) :: v_indices(:)

    call assert_matchup( k )
    matchups(k)%jv = v_indices
  end subroutine set_jv

  subroutine set_iw( k, w_indices )
    integer, intent(in) :: k
    integer, intent(in) :: w_indices(:)

    call assert_matchup( k )
    matchups(k)%iw = w_indices
  end subroutine set_iw

  subroutine set_jw( k, w_indices )
    integer, intent(in) :: k
    integer, intent(in) :: w_indices(:)

    call assert_matchup( k )
    matchups(k)%jw = w_indices
  end subroutine set_jw

  subroutine set_iz( k, z_index )
    integer, intent(in) :: k
    integer, intent(in) :: z_index

    call assert_matchup( k )
    matchups(k)%iz = z_index
  end subroutine set_iz

  subroutine set_jz( k, z_index )
    integer, intent(in) :: k
    integer, intent(in) :: z_index

    call assert_matchup( k )
    matchups(k)%jz = z_index
  end subroutine set_jz

  pure logical function has_pair( i, j )
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer :: k

    has_pair = .false.
    do k = 1, nmatchup
      if (matchups(k)%i == i .and. matchups(k)%j == j) then
        has_pair = .true.
        exit
      end if
    end do
  end function has_pair

  pure logical function is_pair( i, j )
    integer, intent(in) :: i
    integer, intent(in) :: j

    is_pair = i < j .and. i > 0
  end function is_pair

  subroutine assert_matchup( k )
    integer, intent(in) :: k

    if (k < 1 .or. k > nmatchup) then
      !! @todo handle error
      error stop "Invalid matchup index"
    end if
  end subroutine assert_matchup

end module mod_matchup
