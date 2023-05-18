program main
  use precision
  use common
  implicit none

  !######################################################################
  !
  !             dq_i/dt + df_i/dx = 0, for x \in [a,b]
  !
  ! t=0                                 t=tEnd
  ! Density                             Density
  !   ****************|                 *********\
  !                   |                           \
  !                   |                            \
  !                   |                             ****|
  !                   |                                 |
  !                   |                                 ****|
  !                   ***************                       ***********
  !
  ! Domain cells (I{i}) reference:
  !
  !                |           |   u(i)    |           |
  !                |  u(i-1)   |___________|           |
  !                |___________|           |   u(i+1)  |
  !                |           |           |___________|
  !             ...|-----0-----|-----0-----|-----0-----|...
  !                |    i-1    |     i     |    i+1    |
  !                |-         +|-         +|-         +|
  !              i-3/2       i-1/2       i+1/2       i+3/2
  !
  !######################################################################


  !######################################################################
  
  ! Variables 
  integer :: halfcells
  integer :: it                                                                                    ! Iteration counter
  real(dp) :: t                                                                                    ! Time
  real(dp) :: x(0:ncells)                                                                          ! Meshs  
  real(dp) :: u0(0:ncells), u(0:ncells)                                                            ! Velocity
  real(dp) :: p0(0:ncells), p(0:ncells)                                                            ! Pressure
  real(dp) :: rho0(0:ncells), rho(0:ncells)                                                        ! Density
  real(dp) :: E0(0:ncells), E(0:ncells)                                                            ! Total energy density
  real(dp) :: a0(0:ncells), a(0:ncells)                                                            ! Speed of sound
  real(dp) :: q0(0:ncells, 0:2), q(0:ncells, 0:2)                                                  ! Vector of conserved variables (U_{i} for all i)
  real(dp) :: dF(0:ncells-2, 0:2)                                                                  ! U_{i}^{n+1} = U_{i}^{n} - dt/dx*dF_{i}^{n}
                                                                                                   ! dF is defined for each cell, except the boundaries
  character(1) :: M 
  character(20) :: arg1
  character(3) :: N 
  character(20) :: arg2
  
  !######################################################################

  if (iargc()<1) then
     write(*,*) 'Please, insert L for Lax-Wendroff, R for Roe or M for MUSCL'
     stop
  end if

  call getarg(1,arg1)
  read(arg1,*) M
  
  dx = (x_fin-x_ini)/ncells                    
  halfcells = int(ncells/2.0_dp)
  do i = 0, ncells
     x(i) = x_ini + i*dx
     !x(i) = x_ini+0.5_dp*dx+i*dx*(1.0_dp-1.0_dp/(2.0_dp*ncells))
  end do

  ! Initial conditions (Sod shock tube problem)
  u0(0:halfcells) = 0.0_dp
  u0(halfcells:ncells) = 0.0_dp
  p0(0:halfcells) = 1.0_dp
  p0(halfcells:ncells) = 0.1_dp
  rho0(0:halfcells) = 1.0_dp
  rho0(halfcells:ncells) = 0.125_dp

  E0 = p0/((gamma-1.0_dp)*rho0)+0.5_dp*u0**2
  a0 =  sqrt(gamma*p0/rho0)
  q0 = array_to_matrix(rho0, rho0*u0, rho0*E0, ncells)
  q = q0
  
  ! Solver loop
  t = 0
  it = 0
  a = a0
  dt = CFL*dx/maxval(abs(u0)+a0)                                                                  ! abs(u0) + a0 = module of the largest eigenvalue

  select case(M)
     
  case('L')
     print*, '###############################################################################'
     print*, 'Plotting Position, Denisty, Velocity, Pressure and Energy at various time steps'
     print*, '*******************************   Lax-Wendroff   *******************************'
     open(101, file = 'sol.dat')
     
  case('R')
     print*, '###############################################################################'
     print*, 'Plotting Position, Denisty, Velocity, Pressure and Energy at various time steps'
     print*, '***********************************   Roe   ***********************************'
     open(101, file = 'sol.dat')

  case('M')
     if (iargc()<2) then
        write(*,*) 'Please, insert MIN for MINMOD or SUP for SUPERBEE'
        stop
     end if
     call getarg(2,arg2)
     read(arg2,*) N
     select case(N)
     case('MIN')
        beta = 1.0_dp
     case('SUP')
        beta = 2.0_dp
     case default
        print*, 'Please, insert MIN (MINMOD) or SUP (SUPERBEE), no other cases are allowed'
        stop
     end select
     print*, '###############################################################################'
     print*, 'Plotting Position, Denisty, Velocity, Pressure and Energy at various time steps'
     print*, '**********************************   MUSCL   **********************************'
     open(101, file = 'sol.dat')

  case default
     print*, 'Please, insert a letter between L (Lax-Wendroff), R (Roe) and M (MUSCL)'

  end select

  do while(t < tEnd)

     select case(M)

     case('L')
        dF = flux_LW(q)

     case('R')        
        dF = flux_R(q)

     case('M')
        dF = flux_M(q)

     case default
        stop

     end select

     q0 = q
     q(1:ncells-1,:) = q0(1:ncells-1,:)-dt/dx*dF
     q(0,:) = q0(0,:)                                                                            ! Neumann boundary conditions
     q(ncells,:) = q0(ncells,:)

     ! Compute primary variables
     rho = q(:,0)
     u = q(:,1)/rho
     E = q(:,2)/rho
     p = (gamma-1.0_dp)*rho*(E-0.5_dp*u**2)
     a = sqrt(gamma*p/rho)   
     if (minval(p)<0) then
        print*, 'negative pressure found!'
     end if

     ! Writing the results
     if (MOD(it,40) == 0) then
        write(101,*) ""
        write(101,*) ""
        print*, '-------------------------------------------------------------------------------'
        print*, 'Iteration  =', it,  '        ----->        ', 't  =', t       
        do i = 0, ncells
           write(101,*) x(i), rho(i), u(i), p(i), E(i)
        end do
     end if
     
     ! Update/correct time step
     dt = CFL*dx/maxval(abs(u)+a)
    
     ! Update time and iteration counter
     t = t+dt
     it = it+1

  end do
     
 print*, '###############################################################################'
  
end program main
  
