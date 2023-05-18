module common
  use precision
  implicit none
  private

  ! Parameters
  real(dp), public, parameter :: CFL = 0.5_dp                         ! Courant number
  real(dp), public, parameter :: gamma = 1.4_dp                       ! Ratio of specific heats
  integer, public, parameter :: ncells = 400                          ! Number of cells
  real(dp), public, parameter :: x_ini = 0.0_dp                       ! Step size
  real(dp), public, parameter :: x_fin = 1.0_dp                       ! Number of points
  real(dp), public, parameter :: tEnd = 0.2_dp                        ! Final time
  character(10), public, parameter :: minus = 'minus'                 ! Logical variables
  character(10), public, parameter :: plus = 'plus'

  ! Variables
  real(dp), public :: dt                                              ! Time step
  real(dp), public :: dx                                              ! Spatial step
  real(dp), public :: beta                                            ! Limiter function selecter (beta=1-->MINMOD; beta=2-->SUPERBEE)
  integer, public :: i, j, k, l                                       ! Loop indices
  
  
  ! Functions and subroutines
  public :: array_to_matrix
  public :: func_flux
  public :: flux_LW
  public :: flux_R
  public :: flux_M


contains

    
  function array_to_matrix(M0, M1, M2, up) result(M)
    integer :: up                                                     ! Upper bound (it is ncells for volume variables and ncells-1 for surface ones)
    real(dp) :: M(0:up,0:2)
    real(dp) :: M0(0:up)
    real(dp) :: M1(0:up)
    real(dp) :: M2(0:up)
    
    do i = 0, up
       M(i, 0:2) = [M0(i), M1(i), M2(i)]
    end do
    
  end function array_to_matrix

  
  function func_flux(q, up) result(flux)
    integer :: up
    real(dp) :: rho(0:up)
    real(dp) :: u(0:up)
    real(dp) :: E(0:up)
    real(dp) :: p(0:up)
    real(dp) :: F0(0:up)                                              ! First component of the flux
    real(dp) :: F1(0:up)                                              ! Second component of the flux
    real(dp) :: F2(0:up)                                              ! Third component of the flux
    real(dp) :: q(0:up, 0:2)
    real(dp) :: flux(0:up, 0:2)                                       ! Vector of the flux (F_{i} for all i)
    
    ! Primitive variables
    rho = q(:,0)
    u = q(:,1)/rho
    E = q(:,2)/rho
    p = (gamma-1.0_dp)*rho*(E-0.5_dp*u**2)
    
    ! Flux vector of conserved variables
    F0 = rho*u
    F1 = rho*u**2+p
    F2 = u*(rho*E+p)
    flux = array_to_matrix(F0, F1, F2, up)
    
  end function func_flux

  
  function flux_LW(q) result(dF)
    
    real(dp) :: rho(0:ncells)
    real(dp) :: u(0:ncells)
    real(dp) :: E(0:ncells)
    real(dp) :: p(0:ncells)
    real(dp) :: F0(0:ncells)
    real(dp) :: F1(0:ncells)
    real(dp) :: F2(0:ncells)
    real(dp) :: q(0:ncells, 0:2)
    real(dp) :: qpHalf(0:ncells-1, 0:2)                               ! q_{i+1/2} (surface variables have one less component)
    real(dp) :: FqpHalf(0:ncells-1, 0:2)                              ! flux_{i+1/2}
    real(dp) :: flux(0:ncells, 0:2)                
    real(dp) :: dF(0:ncells-2, 0:2)

    ! Predictor step                                                  ! U_{i+1/2}^{tilde} = 1/2*(U_{i+1}^{n}+U_{i}^{n})-dt/(2*dx)*(F_{i+1}^{n}-F{i}^{n}) 

    flux = func_flux(q,ncells)

    qpHalf = (q(1:ncells,:)+q(0:ncells-1,:))/2.0_dp - dt/(2.0_dp*dx)*(flux(1:ncells,:)-flux(0:ncells-1,:))
           
    ! Corrector step                                                  ! U_{i}^{n+1} = U_{i}^{n}-dt/dx*(F_{i+1/2}^{tilde}-F_{i-1/2}^{tilde})  
    
    rho(0:ncells-1) = qpHalf(:,0)
    u(0:ncells-1) = qpHalf(:,1)/rho(0:ncells-1)
    E(0:ncells-1) = qpHalf(:,2)/rho(0:ncells-1)
    p = (gamma-1.0_dp)*rho*(E-0.5_dp*u**2)
    F0 = rho*u
    F1 = rho*u**2+p
    F2 = u*(rho*E+p)
    FqpHalf = array_to_matrix(F0(0:ncells-1), F1(0:ncells-1), F2(0:ncells-1), ncells-1)

    dF = FqpHalf(1:ncells-1,:) - FqpHalf(0:ncells-2,:)

  end function flux_LW

  
  function flux_R(q) result(dF)                                       ! Phi_{i+1/2}^{n} = 1/2*(F_{i+1}^{n}+F_{i}^{n})-1/2*|M_{1/2}|*(U_{i+1}^{n}-U_{i}^{n})
                                                                      ! U_{i}^{n+1} = U_{i}^{n} -dt/dx*(Phi_{i+1/2}^{n}-Phi_{i-1/2}^{n})
    real(dp) :: rho(0:ncells)
    real(dp) :: u(0:ncells)
    real(dp) :: E(0:ncells)
    real(dp) :: p(0:ncells)
    real(dp) :: H(0:ncells)                                           ! Enthalpy
    real(dp) :: rho_hat                                               ! _hat --> Roe average
    real(dp) :: u_hat
    real(dp) :: H_hat
    real(dp) :: a_hat                                                 ! Speed of sound
    real(dp) :: alph1                                                 ! Auxiliary variables
    real(dp) :: alph2
    real(dp) :: R                                                     ! Roe factor
    real(dp) :: M_hat(0:2,0:2)                                        ! Roe matrix
    real(dp) :: P_hat_inv(0:2,0:2)                                    ! M = P Lambda Pinv
    real(dp) :: P_hat(0:2,0:2)
    real(dp) :: Lambda_hat(0:2,0:2)
    real(dp) :: Phip(0:ncells-1, 0:2)                                 ! Roe flux at i+1/2
    real(dp) :: q(0:ncells, 0:2)
    real(dp) :: delta_q(0:2)
    real(dp) :: corr(0:ncells-1, 0:2)                                 ! Auxiliary quantity
    real(dp) :: flux(0:ncells, 0:2)
    real(dp) :: dF(0:ncells-2, 0:2)

    ! Primitive variables
    rho = q(:,0)
    u = q(:,1)/rho
    E = q(:,2)/rho
    p = (gamma-1.0_dp)*rho*(E-0.5_dp*u**2)
    H = gamma/(gamma-1)*p/rho+0.5_dp*u**2

    ! Compute Roe averages at positions i+1/2 for all i
    do i = 0, ncells-1                                                ! max_{i} = ncells-1 because Roe averages are surface variables
                                    
       R = sqrt(rho(i+1)/rho(i))
       u_hat = (R*u(i+1)+u(i))/(R+1.0_dp)
       H_hat = (R*H(i+1)+H(i))/(R+1.0_dp)     
       rho_hat = R*rho(i)
       a_hat = sqrt((gamma-1.0_dp)*(H_hat-u_hat**2/2.0_dp))

       !Compute matrices
       alph1 = (gamma-1.0_dp)*u_hat**2/(2.0_dp*a_hat**2)
       alph2 = (gamma-1.0_dp)/(a_hat**2)
       
       P_hat_inv = reshape((/ &
            0.5_dp*(alph1+u_hat/a_hat), -0.5_dp*(alph2*u_hat+1.0_dp/a_hat), alph2/2.0_dp, &
            1-alph1, alph2*u_hat, -alph2, &
            0.5_dp*(alph1-u_hat/a_hat), -0.5_dp*(alph2*u_hat-1/a_hat), alph2/2.0_dp &
            /), [3,3])

       P_hat = reshape((/ &
            1.0_dp, 1.0_dp, 1.0_dp, &
            u_hat-a_hat, u_hat, u_hat+a_hat, &
            H_hat-a_hat*u_hat, 0.5_dp*u_hat**2, H_hat+a_hat*u_hat &
            /), [3,3])

       Lambda_hat = reshape((/ &
            abs(u_hat-a_hat), 0.0_dp, 0.0_dp, &
            0.0_dp, abs(u_hat), 0.0_dp, &
            0.0_dp, 0.0_dp, abs(u_hat+a_hat) &
            /), [3,3])
 
       P_hat_inv = transpose(P_hat_inv)                               ! In Fortran matrices are defined with rows and columns reversed
       P_hat = transpose(P_hat)
       Lambda_hat = transpose(Lambda_hat)

       M_hat = matmul(P_hat,Lambda_hat)
       M_hat = matmul(M_hat, P_hat_inv)   

       ! Compute vector (U_{i+1}-U_{i})
       delta_q(:) = q(i+1,:)-q(i,:)

       ! corr_{i}=1/2*(A_{i+1/2})*(U_{i+1}-U_{i}) 
       corr(i,:) = 0.5_dp*matmul(M_hat,delta_q)     

    end do

    ! Compute Phi=1/2*(F(U_{i+1})+F(U_{i}))-1/2*|A_{i+1/2}|*(U_{i+1}-U_{i})
    flux = func_flux(q,ncells)
    Phip = 0.5_dp*(flux(1:ncells,:)+flux(0:ncells-1,:))-corr

    dF = Phip(1:ncells-1,:)-Phip(0:ncells-2,:)  
 
  end function flux_R
  

  function Psi(r) result(psi_r)                                       ! Flux limiter function: MINMOD   -->  Psi(r) = max(0,min(1,r))
                                                                      !                        SUPERBEE -->  Psi(r) = max(0,min(1,2r),min(2,r))
    real(dp) :: r
    real(dp) :: psi_r

    psi_r = max(max(0.0_dp, min(beta*r,1.0_dp)), min(r, beta))
  
  end function Psi
  

  function flux_M(q) result(dF)                                       ! U_{i+1/2}^{L} = U_{i}+1/2*Psi_{i}*(U_{i+1}-U_{i})
                                                                      ! U_{i+1/2}^{R} = U_{i+1}-1/2*Psi_{i+1}*(U_{i+2}-U_{i+1})
    real(dp) :: q(0:ncells, 0:2)
    real(dp) :: dqL(1:ncells-1, 0:2)                                  ! Matrix to store the values of Psi(r_{i})                                 
    real(dp) :: dqR(0:ncells-2, 0:2)                                  ! Matrix to store the values of Psi(r_{i+1})
    real(dp) :: qL(0:ncells-1, 0:2)                                   ! U_{i+1/2}^{L}
    real(dp) :: qR(0:ncells-1, 0:2)                                   ! U_{i+1/2}^{R}
    real(dp) :: num                                                   ! Flux limiter argument = r = num/den
    real(dp) :: den
    ! #############################
    ! Roe variables                                                   ! In MUSCL scheme these are surface variables (index i+1/2 goes up to i=ncells-1)
    real(dp) :: rhoL(0:ncells-1), rhoR(0:ncells-1)                                  
    real(dp) :: uL(0:ncells-1), uR(0:ncells-1)
    real(dp) :: EL(0:ncells-1), ER(0:ncells-1)
    real(dp) :: pL(0:ncells-1), pR(0:ncells-1)
    real(dp) :: HL(0:ncells-1), HR(0:ncells-1)
    real(dp) :: fluxL(0:ncells-1, 0:2), fluxR(0:ncells-1, 0:2)
    real(dp) :: rho_hat                                               
    real(dp) :: u_hat
    real(dp) :: H_hat
    real(dp) :: a_hat                                                 
    real(dp) :: alph1                                                 
    real(dp) :: alph2
    real(dp) :: R                                                     
    real(dp) :: M_hat(0:2,0:2)                                        
    real(dp) :: P_hat_inv(0:2,0:2)                                    
    real(dp) :: P_hat(0:2,0:2)
    real(dp) :: Lambda_hat(0:2,0:2)
    real(dp) :: Phip(0:ncells-1, 0:2)                                 
    real(dp) :: delta_q(0:2)
    real(dp) :: corr(0:ncells-1, 0:2)                                 
    real(dp) :: dF(0:ncells-2, 0:2)
    real(dp), parameter :: eps = 1e-8

    ! Compute limiter functions                                        ! Most internal cycle must be on the rows for cache memory optimization
    do k = 0, 2                                                        ! Since r_{i} = (U_{i}/U_{i-1})/(U_{i+1}-U_{i}): Psi(r_{i}) is not defined on the cells i=0
       do l = 1, ncells-1                                              
          num = q(l,k)-q(l-1,k)                                        ! Flux limiter argument = (U_{i}/U_{i-1})/(U_{i+1}-U_{i})
          den = q(l+1,k)-q(l,k)
          if (abs(num) .lt. eps) then                                  ! Conditions to take into account computational errors
             num = 0
             den = 1
          else if (num .gt. eps .and. den .lt. eps) then
             num = 1
             den = 1
          else if (num .lt. -eps .and. den .lt. eps) then
             num = -1
             den = 1
          end if
          dqL(l,k) = Psi(num/den)
       end do
    end do
    
    do k = 0, 2                                                        
       do l = 0, ncells-2                                              ! Psi(r_{i+1}) it is not defined for the cell i=ncells-1
          num = q(l+1,k)-q(l,k)
          den = q(l+2,k)-q(l+1,k)
          if (abs(num) .lt. eps) then
             num = 0
             den = 1
          else if (num .gt. eps .and. den .lt. eps) then
             num = 1
             den = 1
          else if (num .lt. -eps .and. den .lt. eps) then
             num = -1
             den = 1
          end if
          dqR(l,k) = Psi(num/den)       
       end do
    end do
    
    ! Compute qL and qR
    qL(1:ncells-1,:) = q(1:ncells-1,:)+0.5_dp*dqL(1:ncells-1,:)*(q(2:ncells,:)-q(1:ncells-1,:))
    qR(0:ncells-2,:) = q(1:ncells-1,:)-0.5_dp*dqR(0:ncells-2,:)*(q(2:ncells,:)-q(1:ncells-1,:))
    ! At the boundaries (i.e. i=1/2 and i=ncells-1/2) we require qL_{1/2}=qR_{1/2}, qR_{ncells-1/2}=qL{ncells-1/2})
    qL(0,:) = qR(0,:)
    qR(ncells-1,:) = qL(ncells-1,:)

    ! ###############################################################################
    
    ! We now have to choose a scheme to campute dF; in this case we use Roe's scheme
    ! To compute surface variables at i+1/2 we can use qL and qR 

    ! Primitive right variables
    rhoR = qR(:,0)
    uR = qR(:,1)/rhoR
    ER = qR(:,2)/rhoR
    pR = (gamma-1.0_dp)*rhoR*(ER-0.5_dp*uR**2)
    HR = gamma/(gamma-1)*pR/rhoR+0.5_dp*uR**2
    ! Primitive left variables
    rhoL = qL(:,0)
    uL = qL(:,1)/rhoL
    EL = qL(:,2)/rhoL
    pL = (gamma-1.0_dp)*rhoL*(EL-0.5_dp*uL**2)
    HL = gamma/(gamma-1)*pL/rhoL+0.5_dp*uL**2

    ! Compute Roe averages
    do i = 0, ncells-1                                                
                                    
       R = sqrt(rhoR(i)/rhoL(i))
       u_hat = (R*uR(i)+uL(i))/(R+1.0_dp)
       H_hat = (R*HR(i)+HL(i))/(R+1.0_dp)     
       rho_hat = R*rhoL(i)
       a_hat = sqrt((gamma-1.0_dp)*(H_hat-u_hat**2/2.0_dp))

       !Compute matrices
       alph1 = (gamma-1.0_dp)*u_hat**2/(2.0_dp*a_hat**2)
       alph2 = (gamma-1.0_dp)/(a_hat**2)
       
       P_hat_inv = reshape((/ &
            0.5_dp*(alph1+u_hat/a_hat), -0.5_dp*(alph2*u_hat+1.0_dp/a_hat), alph2/2.0_dp, &
            1-alph1, alph2*u_hat, -alph2, &
            0.5_dp*(alph1-u_hat/a_hat), -0.5_dp*(alph2*u_hat-1/a_hat), alph2/2.0_dp &
            /), [3,3])

       P_hat = reshape((/ &
            1.0_dp, 1.0_dp, 1.0_dp, &
            u_hat-a_hat, u_hat, u_hat+a_hat, &
            H_hat-a_hat*u_hat, 0.5_dp*u_hat**2, H_hat+a_hat*u_hat &
            /), [3,3])

       Lambda_hat = reshape((/ &
            abs(u_hat-a_hat), 0.0_dp, 0.0_dp, &
            0.0_dp, abs(u_hat), 0.0_dp, &
            0.0_dp, 0.0_dp, abs(u_hat+a_hat) &
            /), [3,3])

       P_hat_inv = transpose(P_hat_inv)                                ! In Fortran matrices are defined with rows and columns reversed
       P_hat = transpose(P_hat)
       Lambda_hat = transpose(Lambda_hat)

       M_hat = matmul(P_hat,Lambda_hat)
       M_hat = matmul(M_hat, P_hat_inv)   

       ! Compute vector (U_{i+1/2}^{R}-U_{i+1/2}^{L})
       delta_q(:) = qR(i,:)-qL(i,:)

       ! corr_{i}=1/2*(A_{i+1/2})*(U_{i+1/2}^{R}-U_{i+1/2}^{L}) 
       corr(i,:) = 0.5_dp*matmul(M_hat,delta_q)     

    end do

    ! Compute Phi=1/2*(F(U_{i+1/2}^{R})+F(U_{i+1/2}^{L}))-1/2*|A_{i+1/2}|*(U_{i+1/2}^{R}-U_{i+1/2}^{L})
    fluxR = func_flux(qR,ncells-1)
    fluxL = func_flux(qL,ncells-1)
    Phip = 0.5_dp*(fluxR+fluxL)-corr

    dF = Phip(1:ncells-1,:)-Phip(0:ncells-2,:)

  end function flux_M
 
  
end module common
