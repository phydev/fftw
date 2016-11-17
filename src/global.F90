!! Copyright (C) 2016 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!!

module global_m

  implicit none

  private

  public :: mesh_t

  type mesh_t
     real :: phi, lapl_phi
  end type mesh_t
    !> some mathematical constants
    REAL(8), public, parameter :: M_Pi        = 3.1415926535897932384626433832795029
    REAL(8), public, parameter :: M_E         = 2.7182818284590452353602874713526625
    REAL(8), public, parameter :: M_TWO       = 2.0
    REAL(8), public, parameter :: M_THREE     = 3.0
    REAL(8), public, parameter :: M_FOUR      = 4.0
    REAL(8), public, parameter :: M_FIVE      = 5.0
    REAL(8), public, parameter :: M_SIX       = 6.0
    REAL(8), public, parameter :: M_SEVEN     = 7.0
    REAL(8), public, parameter :: M_EIGHT     = 8.0
    REAL(8), public, parameter :: M_NINE      = 9.0
    REAL(8), public, parameter :: M_TEN       = 10.0
    REAL(8), public, parameter :: M_HALF      = 0.5

      ! begin parameters
      real, allocatable, public :: kvec(:,:)
      integer, public :: i, j, k
      real, public ::  diffusion_const, interface_width, dt
      integer, public :: tstep, iseed, dr(3)
      character(len=6), public :: file_name, file_id
      character(len=17), public :: file_phi
      character(len=3),public :: sim_id
      ! mesh variables
      integer, public, allocatable :: lxyz(:,:), lxyz_inv(:,:,:)
      integer(4), public :: Lsize(1:3), Lx, Ly, Lz
      integer, public :: np, np_part, np_complex, ip, boundary_points


      ! misc
      integer, public :: nstep, counter = 0, ndim, output_period, extra_steps
      real, public :: hs, time_init, time_end, ctime, phi_max
      logical, public :: periodic


end module global_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
