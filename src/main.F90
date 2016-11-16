program main

  use global_m
  use init_m

  implicit none

  call get_command_argument(1,sim_id)
  Lsize(:) = 20
  boundary_points = 2.0
  np = 8*Lsize(1)*Lsize(2)*Lsize(3) ! number of points
  np_part = 8*(Lsize(1)+2*boundary_points)*(Lsize(2)+2*boundary_points)*(Lsize(3)+2*boundary_points) ! number of points plus boundary points
  ! allocating matrices and vectors
  ALLOCATE(lxyz(np_part,1:3))
  ALLOCATE(lxyz_inv(-Lsize(1)-boundary_points:Lsize(1)+boundary_points, &
  -Lsize(2)-boundary_points:Lsize(2)+boundary_points, &
  -Lsize(3)-boundary_points:Lsize(3)+boundary_points))
  ALLOCATE(rho(np))

  call print_header(Lsize, sim_id)
  call space_init(Lsize, lxyz, lxyz_inv, boundary_points, np)
  call simul_box_init(rho, lxyz, lxyz_inv, Lsize, np)

  


end program main
