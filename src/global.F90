!! Copyright (C) 2015 M. Moreira
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

      type(C_PTR) :: plan
      real(C_DOUBLE), dimension(L,M,N) :: in
      complex(C_DOUBLE_COMPLEX), dimension(L/2+1,M,N) :: out
      ! begin parameters
      integer, public :: i, j, k
      real, public ::  diffusion_const, interface_width, dt
      integer, public :: tstep, iseed, dr(3)
      character(len=6), public :: file_name, file_id
      character(len=17), public :: file_phi
      character(len=3),public :: sim_id
      ! mesh variables
      integer, public, allocatable :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, public :: Lsize(1:3), boundary_points
      integer, public :: np, np_part, ip

      type(mesh_t), public, allocatable :: rho(:)

      ! misc
      integer, public :: nstep, counter = 0, ndim, output_period, extra_steps
      real, public :: hs, time_init, time_end, ctime, phi_max
      logical, public :: periodic


end module global_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
