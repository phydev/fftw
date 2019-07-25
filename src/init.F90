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

module init_m

  use global_m

  implicit none

  private

  public :: wave_function_init, simul_box_init, space_init, print_header, heaviside, ran2

  contains

    subroutine wave_function_init(kvec, Lx, Ly, Lz)

      implicit none
      real, allocatable, intent(inout) :: kvec(:,:)
      integer(4), intent(in) ::  Lx, Ly, Lz
      integer :: i,j,k
      ALLOCATE(kvec(0:max(Lx,Ly,Lz),1:3))
      do i=0, Lx/2
        kvec(i,1) = 2*M_PI*i/Lx
        kvec(Lx-i,1) = -2*M_PI*i/Lx
      end do

      do j=0, Ly/2
        kvec(j,2) = 2*M_PI*j/Ly
        kvec(Ly-j,2) = -2*M_PI*j/Ly
      end do

      do k=0, Lz/2
        kvec(k,3) = 2*M_PI*k/Lz
        kvec(Lz-k,3) = -2*M_PI*k/Lz
      end do


    end subroutine wave_function_init

    subroutine simul_box_init(rho, lxyz, lxyz_inv, Lsize, np)

      implicit none

      ! input variables
      !type(C_DOUBLE), allocatable, intent(inout) :: rho(:)
      real, allocatable, intent(inout) :: rho(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) ::np
      integer(4), intent(in) :: Lsize(3)
      ! internal only
      integer ::  ip, i, j, k
      iseed = -92999384
      rho(:) = 0.d8
      do ip = 0, np-1

         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)

         if(sqrt(real(i)**2 + real(j)**2 + real(k)**2) .le. 5.d0) then
           rho(ip) = 1.d0
         else
           rho(ip) = ran2(iseed)
         end if


      end do


    end subroutine simul_box_init


    subroutine space_init(Lsize, lxyz, lxyz_inv, boundary_points, np)

      implicit none

      ! input/output variables
      integer, intent(in) ::  boundary_points, np
      integer(4), intent(in) :: Lsize(3)
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



               lxyz(ip,1) = i
               lxyz(ip,2) = j
               lxyz(ip,3) = k

               lxyz_inv(i,j,k) = ip

               ip = ip + 1

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
                  n = k  - SIGN(1,k)*(2*Lsize(3))
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
      integer(4), intent(in) :: Lsize(3)
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

    function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer :: idum2,j,k,iv(NTAB),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1

            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11          continue
            iy=iv(1)
         endif
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         k=idum2/IQ2
         idum2=IA2*(idum2-k*IQ2)-k*IR2
         if (idum2.lt.0) idum2=idum2+IM2
         j=1+iy/NDIV
         iy=iv(j)-idum2
         iv(j)=idum
         if(iy.lt.1)iy=iy+IMM1
         ran2=min(AM*iy,RNMX)
         return
       end function ran2

  end module init_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
