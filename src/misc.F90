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

module misc_m

  use global_m

  implicit none

  private

  public :: convert_arrays

  contains


    subroutine convert_arrays(cell, lxyz, lxyz_inv, Lsize, np)

      implicit none

      ! input variables
      type(mesh_t), allocatable, intent(inout) :: cell(:)



    end subroutine simul_box_init



  end module init_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
