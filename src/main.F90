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

program main

  use global_m
  use init_m
  use, intrinsic :: iso_c_binding


  implicit none
  include 'fftw3.f03'
  real(C_DOUBLE), allocatable :: rho(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: rho_complex(:)
  type(C_PTR) :: plan, data

  call get_command_argument(1,sim_id)
  Lsize(:) = 20
  Lx = 2*Lsize(1)
  Ly = 2*Lsize(2)
  Lz = 2*Lsize(3)
  boundary_points = 2.0
  np = 8*Lsize(1)*Lsize(2)*Lsize(3) ! number of points
  np_part = 8*(Lsize(1)+2*boundary_points)*(Lsize(2)+2*boundary_points)*(Lsize(3)+2*boundary_points) ! number of points plus boundary points
  np_complex = Lx*Ly*(Lz/2 +1)
  ! allocating matrices and vectors
  ALLOCATE(lxyz(0:np_part-1,1:3))
  ALLOCATE(lxyz_inv(-Lsize(1)-boundary_points:Lsize(1)+boundary_points, &
  -Lsize(2)-boundary_points:Lsize(2)+boundary_points, &
  -Lsize(3)-boundary_points:Lsize(3)+boundary_points))
  ALLOCATE(rho(0:np-1))
  ALLOCATE(rho_complex(0:np_complex))

!! f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
!! f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd

  !data = fftw_alloc_complex(int((L/2+1) * M * N, C_SIZE_T))
  !call c_f_pointer(data, in, [2*(L/2+1),M,N])
  !call c_f_pointer(data, out, [L/2+1,M,N])

  !call fftw_free(data)
   call print_header(Lsize, sim_id)
   call space_init(Lsize, lxyz, lxyz_inv, boundary_points, np)
   call simul_box_init(rho, lxyz, lxyz_inv, Lsize, np)
   call wave_function_init(kvec, Lx, Ly, Lz)

   open(unit=1000,file='rhoi.xyz')
   do ip=0,np-1
     write(1000,'(I10,I10,I10,F10.4)') lxyz(ip,1:3), rho(ip)
   end do
   close(1000)

  ! ! integer L, N, M, double in, complex out, unsigned flags
  ! plan = fftw_plan_dft_r2c_3d(Lsize(1),Lsize(2), Lsize(3),rho,rho_complex,FFTW_ESTIMATE)
  ! call fftw_execute_dft_r2c(plan, rho, rho_complex)
  ! call fftw_destroy_plan(plan)
    plan = fftw_plan_dft_r2c_3d(Lx,Ly,Lz, rho,rho_complex, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, rho, rho_complex)
    call fftw_destroy_plan(plan)
  !rho_complex = rho_complex(ip) ! normalizing
  !do ip=0, np_complex
  !  write(*,*) ip, rho_complex(ip)
    ! if(abs(rho_complex(ip))>0) i = i + 1
  !end do
!write(*,*)i, Lx*Ly*(Lz/2 +1), Lx*Ly*Lz
  plan = fftw_plan_dft_c2r_3d(Lx,Ly,Lz, rho_complex, rho, FFTW_ESTIMATE)
  call fftw_execute_dft_c2r(plan, rho_complex,  rho)
  call fftw_destroy_plan(plan)
  rho(:) = rho(:)/(Lx*Ly*Lz)
  open(unit=1000,file='rhof.xyz')
  do ip=0,np-1
    write(1000,'(I10,I10,I10,F10.4)') lxyz(ip,1:3), rho(ip)
  end do
  close(1000)
  ! FFTW computes an unnormalized transform, that is, the equation IFFT(FFT(X)) = n X holds. In other words, applying the forward and then the backward transform will multiply the input by n.



  !call fftw_plan_dft_c2r_3d(int n0, int n1, int n2, fftw_complex *in, double *out, unsigned flags);

end program main
