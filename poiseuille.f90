!==========================================!
! Gravity driven Poiseuille flow with LBM  !
!          Apr. 2020, Slovenia             !
!==========================================!

!=======================================================================!
! Author: Anze Hubman                                                   !
!         FCCT Ljubljana & National Institute of Chemistry/Theory dept. !
!=======================================================================!

!===================================================================!
!Code description:                                                  !
! This code produces pipe flow driven by gravity. D2Q9 velocity     !
! set is used in standard notation. At top and bottom (inlet and    !
! outlet) periodic boundary conditions are applied, while at        !
! the side walls half-wa bounce back is used. Gravity is a constant !
! external body force and is included via a Shan-Chen scheme. The   !
! simulation is initialsed as a fluid with u(x,y,t) = 0.0           !
!===================================================================!

!===============================================================================!
! DISCLAIMER: This code can be used and modified freely.                        !
!             I am not taking responsibility for any mistakes that might occur. !
!===============================================================================!

!=============================================================================!
! WARNINGS:                                                                   !
! 1.) Input parameters are given in common LB lattice units.                  !
!     User should make conversion to physical units on its own.               !
!                                                                             !
! 2.) Origin of coordinate system is placed at the top left corner.           !
!                                                                             !
! 3.) Gravity term should not be too high. Be aware of intrinsic              !
!     LBM limitations.                                                        !
!                                                                             !
! 4.) When increasing Reynolds number grid resolution has to be increased.    !
!                                                                             !
! 5.) When compiling the code I used -O3 optimisation level to speed up the   !
!     calculations. -Ofast option was not tested.                             !
!=============================================================================!

MODULE params
  REAL*8,  PARAMETER :: tau   = 1.0D0   !LBM relaxation parameter
  REAL*8,  PARAMETER :: dt    = 1.0D0   !Time step
  REAL*8,  PARAMETER :: g     = 1.0D-05 !Body force density (gravity)
  INTEGER, PARAMETER :: nstep = 20000   !Number of time steps
  INTEGER, PARAMETER :: N_x   = 50      !Number of nodes along x-axis          
  INTEGER, PARAMETER :: N_z   = 250     !Number of nodes along z-axis
END MODULE params


PROGRAM gravity
  USE params
  IMPLICIT NONE

  REAL*8, ALLOCATABLE :: rho(:,:), ux(:,:), uz(:,:), usq(:,:)
  REAL*8, ALLOCATABLE :: feq_k(:,:,:), f(:,:,:)
  INTEGER :: i,j,k

  ALLOCATE(rho(1:N_z,1:N_x), ux(1:N_z,1:N_x), uz(1:N_z,1:N_x), usq(1:N_z,1:N_x))
  ALLOCATE(feq_k(1:N_z,1:N_x,0:8), f(1:N_z,1:N_x,0:8))
  
  !initialise density and velocity field
  rho(:,:) = 1.0D0     !uniform density
  ux(:,:)  = 0.0D0     !fluid not moving
  uz(:,:)  = 0.0D0
  usq(:,:) = 0.0D0

  !main LBM run
  CALL feq(feq_k,rho,ux,uz,usq)
  f = feq_k
  CALL update(f,rho,ux,uz,usq)

  DO k = 1, nstep
     CALL feq(feq_k,rho,ux,uz,usq)
     CALL collision(feq_k,f)
     CALL stream(f)
     CALL update(f,rho,ux,uz,usq)
  END DO

  !produce output
  OPEN(UNIT=11, FILE="output_1")
  !output structure: x y pressure ux uy u
  DO i = 1, N_z
     DO j = 1, N_x
        WRITE(11,*) j, N_z-i, rho(i,j)/3.0D0, ux(i,j), uz(i,j), SQRT(usq(i,j))
     END DO
  END DO
    
  CLOSE(UNIT=11)
  
END PROGRAM gravity

!============================================!
! Compute equilibrium distribution function  !
!============================================!
SUBROUTINE feq(feq_k,rho,ux,uz,usq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: rho(1:N_z,1:N_x), ux(1:N_z,1:N_x), uz(1:N_z,1:N_x), usq(1:N_z,1:N_x)
  REAL*8, INTENT(INOUT) :: feq_k(1:N_z,1:N_x,0:8)
  INTEGER :: i,j

  DO i = 1, N_z
     DO j = 1, N_x

        feq_k(i,j,0) = (2.0D0*rho(i,j)/9.0D0)*(2.0D0-3.0D0*usq(i,j))
        
        feq_k(i,j,1) = (rho(i,j)/18.0D0)*(2.0D0 + 6.0D0*ux(i,j) + (9*ux(i,j)**2) - 3*usq(i,j))
        feq_k(i,j,2) = (rho(i,j)/18.0D0)*(2.0D0 + 6.0D0*uz(i,j) + (9*uz(i,j)**2) - 3*usq(i,j))
        feq_k(i,j,3) = (rho(i,j)/18.0D0)*(2.0D0 - 6.0D0*ux(i,j) + (9*ux(i,j)**2) - 3*usq(i,j))
        feq_k(i,j,4) = (rho(i,j)/18.0D0)*(2.0D0 - 6.0D0*uz(i,j) + (9*uz(i,j)**2) - 3*usq(i,j))

        feq_k(i,j,5) = (rho(i,j)/36.0D0)*(1.0D0 + 3.0D0*(ux(i,j) + uz(i,j)) + 9.0D0*ux(i,j)*uz(i,j) + 3*usq(i,j))
        feq_k(i,j,6) = (rho(i,j)/36.0D0)*(1.0D0 - 3.0D0*(ux(i,j) - uz(i,j)) - 9.0D0*ux(i,j)*uz(i,j) + 3*usq(i,j))
        feq_k(i,j,7) = (rho(i,j)/36.0D0)*(1.0D0 - 3.0D0*(ux(i,j) + uz(i,j)) + 9.0D0*ux(i,j)*uz(i,j) + 3*usq(i,j))
        feq_k(i,j,8) = (rho(i,j)/36.0D0)*(1.0D0 + 3.0D0*(ux(i,j) - uz(i,j)) - 9.0D0*ux(i,j)*uz(i,j) + 3*usq(i,j))

     END DO
  END DO
END SUBROUTINE feq

!==============!
!BGK collision !
!==============!
SUBROUTINE collision(feq_k,f)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: feq_k(1:N_z,1:N_x,0:8)
  REAL*8, INTENT(INOUT) :: f(1:N_z,1:N_x,0:8)
  INTEGER :: i,j,k

  DO i = 1,N_z
     DO j = 1,N_x
        DO k = 0,8
           f(i,j,k) = (1.0D0-(dt/tau))*f(i,j,k) + (dt/tau)*feq_k(i,j,k)
        END DO
     END DO
  END DO
END SUBROUTINE collision

!====================!
! Stream populations !
!====================!
SUBROUTINE stream(f)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: f(1:N_z,1:N_x,0:8)
  REAL*8 :: ftemp(1:N_z,1:N_x,0:8)

  !0:
  ftemp(:,:,0) = f(:,:,0)

  !1:
  ftemp(:,2:N_x,1) = f(:,1:N_x-1,1)
  ftemp(:,N_x,3) = f(:,N_x,1)

  !2: 
  ftemp(1:N_z-1,:,2) = f(2:N_z,:,2)
  ftemp(N_z,:,2) = f(1,:,2)

  !3:
  ftemp(:,1:N_x-1,3) = f(:,2:N_x,3)
  ftemp(:,1,1) = f(:,1,3)

  !4: 
  ftemp(2:N_z,:,4) = f(1:N_z-1,:,4)
  ftemp(1,:,4) = f(N_z,:,4)

  !5:
  ftemp(1:N_z-1,2:N_x,5) = f(2:N_z,1:N_x-1,5)
  ftemp(N_z,2:N_x,5) = f(1,1:N_x-1,5)
  ftemp(:,N_x,7) = f(:,N_x,5)
  
  !6:
  ftemp(1:N_z-1,1:N_x-1,6) = f(2:N_z,2:N_x,6)
  ftemp(N_z,1:N_x-1,6) = f(1,2:N_x,6)
  ftemp(:,1,8) = f(:,1,6)

  !7: 
  ftemp(2:N_z,1:N_x-1,7) = f(1:N_z-1,2:N_x,7)
  ftemp(1,1:N_x-1,7) = f(N_z,2:N_x,7)
  ftemp(:,1,5) = f(:,1,7)


  !8:
  ftemp(2:N_z,2:N_x,8) = f(1:N_z-1,1:N_x-1,8)
  ftemp(1,2:N_x,8) = f(N_z,1:N_x-1,8)
  ftemp(:,N_x,6) = f(:,N_x,8)

  f = ftemp

END SUBROUTINE stream

!==========================!
! Update macroscopic field !
!==========================!
SUBROUTINE update(f,rho,ux,uz,usq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: f(1:N_z,1:N_x,0:8)
  REAL*8, INTENT(INOUT) :: rho(1:N_z,1:N_x), ux(1:N_z,1:N_x), uz(1:N_z,1:N_x), usq(1:N_z,1:N_x)
  INTEGER :: i,j

  DO i = 1, N_z
     DO j = 1, N_x
        rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
        ux(i,j) = ((f(i,j,1)+f(i,j,5)+f(i,j,8))-(f(i,j,3)+f(i,j,6)+f(i,j,7)))/rho(i,j)
        uz(i,j) = (((f(i,j,2)+f(i,j,5)+f(i,j,6))-(f(i,j,4)+f(i,j,7)+f(i,j,8)))/rho(i,j))-(g*tau)
        usq(i,j) = (ux(i,j)**2)+(uz(i,j)**2)
     END DO
  END DO
END SUBROUTINE update
        
       
