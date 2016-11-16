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

module init_m

  use global_m

  implicit none

  private

  public :: simul_box_init, space_init, print_header, heaviside

  contains


    subroutine simul_box_init(cell, lxyz, lxyz_inv, Lsize, np)

      implicit none

      ! input variables
      type(mesh_t), allocatable, intent(inout) :: cell(:)

      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) ::np, Lsize(3)
      ! internal only
      integer ::  ip, i, j, k

      cell(:)%phi = 0.d8
      do ip = 1, np

         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)

         if(sqrt(real(i)**2 + real(j)**2 + real(k)**2) .le. 5.d0) then
           cell(ip)%phi = 1.d0
         end if


      end do


    end subroutine simul_box_init


    subroutine space_init(Lsize, lxyz, lxyz_inv, boundary_points, np)

      implicit none

      ! input/output variables
      integer, intent(in) ::  Lsize(1:3), boundary_points, np
      integer, allocatable, intent(inout) :: lxyz(:,:), lxyz_inv(:,:,:)
      ! internal variables
      integer :: i, j, k, l, m, n, ip, ip_part
      real :: hs(1:3)
      logical :: boundary

      ip = 0
      ip_part = np

      ! bulk points
      ! allocated from 1 to np

      do i=-Lsize(1), Lsize(1)-1
         do j=-Lsize(2), Lsize(2)-1
            do k=-Lsize(3), Lsize(3)-1

               ip = ip + 1

               lxyz(ip,1) = i
               lxyz(ip,2) = j
               lxyz(ip,3) = k

               lxyz_inv(i,j,k) = ip

            end do
         end do
      end do

      ! boundary points
      ! allocated from np to np_part

      do i=-Lsize(1)-boundary_points, Lsize(1)-1+boundary_points
         do j=-Lsize(2)-boundary_points, Lsize(2)-1+boundary_points
            do k=-Lsize(3)-boundary_points, Lsize(3)-1+boundary_points

               l = i
               m = j
               n = k

               boundary = .false.

               hs(1) = heaviside(real(i))
               hs(2) = heaviside(real(j))
               hs(3) = heaviside(real(k))


                if( abs(i)>Lsize(1)-hs(1)) then
                  boundary = .true.
                  l = i - SIGN(1,i)*(2*Lsize(1))!- SIGN(1,i)*heaviside(-real(i))
               end if

               if( abs(j)>Lsize(2)-hs(2)) then
                  boundary = .true.
                  m = j  - SIGN(1,j)*(2*Lsize(2))!- SIGN(1,j)*heaviside(-real(j))
               end if

               if( abs(k)>Lsize(3)-hs(3)) then
                  boundary = .true.
                  n = SIGN(1,k)*Lsize(3) - hs(3)*SIGN(1,k)
               end if


               if(boundary) then
                  ip_part = ip_part + 1

                  lxyz(ip_part,1) = l
                  lxyz(ip_part,2) = m
                  lxyz(ip_part,3) = n

                  lxyz_inv(i,j,k) = lxyz_inv(l, m, n)
               end if


               ! the updating of boundaries should be
               ! cell(ip_part) = cell( lxyz_inv( lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ) )
               ! because:
               ! lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ->(l, m, n) in the bulk
               ! then lxyz_inv(m,n,l) will give the updated value in the bulk
               ! for the respective boundary position cell(ip_part)

            end do
         end do
      end do

    end subroutine space_init




    subroutine print_header(Lsize, sim_id)

      implicit none
      integer, intent(in) :: Lsize(1:3)
      character(3),intent(in)  :: sim_id
      character(len=255) :: cwd, hostname
      character(len=32) :: username
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer, dimension(8) :: values


      call date_and_time(date,time,zone,values)
      call date_and_time(DATE=date,ZONE=zone)
      call date_and_time(TIME=time)

      call hostnm(hostname)
      call getcwd(cwd)
      call getlog(username)
      write(*,'(A)') "                                Running Poison Solver with FFTW"

      write(*,'(A)') "       "
      write(*,'(A)') "Version        :       0.0.1 (November 16, 2016)"
      write(*,'(A,A)') "Locate         :       ", trim(cwd)
      write(*,'(A,A)') "User           :       ", trim(username)
      write(*,'(A)') "Developer      :       Moreira, M."
      write(*,'(A)') "       "
      write(*,'(A,A)') "                      The code is running in ", trim(hostname)
      write(*,'(A)') "       "
      write(*,'(A)') "       "
      write(*,'(A,2X,A,2X,A)') "             Calculation started on", date(7:8)//"/"//date(5:6)//"/"//date(1:4),&
           time(1:2)//":"//time(3:4)//":"//time(5:6)
      write(*,'(A)') "       "
      write(*,'(A)') "************************************ Grid *************************************"
      write(*,'(A)') "Simulation Box:"
      write(*,'(A,I3,A,I3,A,I3,A)') "  Lengths = (",Lsize(1),",", Lsize(2),",", Lsize(3), ")"
      write(*,'(A)') "  the code will run in 3 dimension(s)."
      write(*,'(A)') "  the code will treat the system as periodic in 3 dimension(s)."

      write(*,'(A)') "*******************************************************************************"
      write(*,'(A,A)') "Simulation ID: ", sim_id
    end subroutine print_header


    function heaviside(x)

      implicit none

      real :: x, heaviside

      if ( x<0.d0) then
         heaviside = 0.d0
      else
         heaviside = 1.d0
      end if
    end function heaviside

  end module init_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
