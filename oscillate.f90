  !==================================================!
  ! Flow induced by vertically oscillating cylinder  !
  !        LBM/Immersed boundary condition           !
  !              Multi direct forcing                !
  !              April 2020, Slovenia                !
  !==================================================!

!=======================================================================!
! Author: Anze Hubman                                                   !
!         FCCT Ljubljana & National Institute of Chemistry/Theory dept. !
!=======================================================================!

!===================================================================!
!Code description:                                                  !
! This code is one of the standard benchmark problems for moving    !
! boundary problems. A circular cylinder is forced to oscillate     !
! in transverse (y) direction. Fluid is initially at rest i.e.      !
! u(x,y,t) = 0.0. Immersed boundary condition with multi direct     !
! forcing is applied at the cylinder to ensure no-slip condition.   !
! Inlet and outlet are subjected to periodic boundary conditions,   !
! while at the top and bottom half-way bounce back is used. Motion  !
! of a cylinder is described as: y(t) =y0 + A*sin(omega*t).         !
!===================================================================!

!================================================================================!
! DISCLAIMER: This code can be used and modified freely.                         !
!             I am not taking responsibility for any mistakes that might occur.  !
!================================================================================!

!===============================================================================!
! WARNINGS:                                                                     !
! 1.) Input parameters are given in common LB lattice units.                    !
!     User should make conversion to physical units on its own.                 !
!                                                                               ! 
! 2.) Choose LB frequency wisely. A frequency that is too high                  !
!     might yield numerically unstable solutions. Be aware of                   !
!     intrinsic LBM limitations.                                                !
!                                                                               !
! 3.) Spacing between Lagrangian markers is important. For 2D problems          !
!     it is recommended that d = 0.33*dh (approx.). Do your own testing!        !
!                                                                               !
! 4.) Obstacle should be well enough resolved in IBC. It is usually enough      !
!     to have at least 20 l.u. accorss diameter.                                !
!                                                                               !
! 5.) Internal mass effect is neglected here. For more accurate calculations    !
!     of forces this should be considered. For example Feng's rigid body        !
!     approximation would be an appropriate choice. It is trivial to change     !
!     this code to include that.                                                !
!                                                                               !
! 6.) When compiling the code I used -O3 optimisation level to speed up the     !
!     calculations. -Ofast option was not tested.                               !
!===============================================================================!

MODULE params
  INTEGER, PARAMETER :: xmax = 1000        !Number of Eulerian nodes along x-axis
  INTEGER, PARAMETER :: ymax = 201         !Number of Eulerian nodes along y-axis
  INTEGER, PARAMETER :: nstp = 500000      !Number of time steps
  INTEGER, PARAMETER :: rc   = 25.0D0      !Cylinder radius
  REAL*8,  PARAMETER :: tau  = 0.757D0     !LBM relaxation parameter
  REAL*8,  PARAMETER :: dh   = 1.0D0       !Lattice spacing
  REAL*8,  PARAMETER :: dt   = 1.0D0       !Time step
  INTEGER, PARAMETER :: nmrk = 420         !Number of Lagrangian marker points
  REAL*8,  PARAMETER :: drat = 1.0D0       !density(solid)/density(liquid)
  INTEGER, PARAMETER :: mmax = 5           !Number of MDF iterations per time step
  REAL*8,  PARAMETER :: xc_0 = 200.5D0     !Initial radius c.o.m. (x)
  REAL*8,  PARAMETER :: yc_0 = 100.5D0     !Initial radius c.o.m. (y)
END MODULE params

PROGRAM lbm
  USE params
  IMPLICIT NONE

  REAL*8, ALLOCATABLE :: lgrx(:), lgry(:), rho(:,:), ux(:,:), uy(:,:), usq(:,:)
  REAL*8, ALLOCATABLE :: ubx(:), uby(:), f(:,:,:), umx(:), umy(:), fmx(:), fmy(:)
  REAL*8, ALLOCATABLE :: fdm(:,:,:), fd(:,:,:), feq(:,:,:)
  REAL*8 :: pi,xc,yc,ampl,freq
  INTEGER :: i,j,k,m

  ALLOCATE(lgrx(1:nmrk),lgry(1:nmrk),rho(0:xmax,0:ymax),ux(0:xmax,0:ymax),uy(0:xmax,0:ymax),usq(0:xmax,0:ymax))
  ALLOCATE(ubx(1:nmrk),uby(1:nmrk),f(0:xmax,0:ymax,0:8),umx(1:nmrk),umy(1:nmrk),fmx(1:nmrk),fmy(1:nmrk))
  ALLOCATE(fdm(0:xmax,0:ymax,1:2),fd(0:xmax,0:ymax,1:2),feq(0:xmax,0:ymax,0:8))

  !define pi
  pi = 4.0D0*ATAN(1.0D0)

  !initial conditions
  xc     = xc_0
  yc     = yc_0
  freq   = 2.0D0*pi/10000.0D0    !angular frequency
  ampl   = 50.0D0                !oscillation amplitude

  OPEN(UNIT=13, FILE="output")

  !main LBM run
  CALL lagrange(lgrx,lgry,pi,xc,yc)
  CALL initialise(rho,ux,uy,usq,ubx,uby,ampl,freq)
  CALL equilibrium(rho,ux,uy,usq,feq)
  f = feq

  DO k = 0, nstp
     WRITE(13,*) k, xc, yc
     PRINT *, k, xc, yc
     CALL move(yc,uby,lgry,k,ampl,freq)
     fd(:,:,:) = 0.0D0
     CALL collision(f,feq)
     CALL stream(f)
     CALL update(rho,ux,uy,usq,f)
     !MDF-LBM
     DO m = 1, mmax
        umx(:) = 0.0D0
        umy(:) = 0.0D0
        fmx(:) = 0.0D0
        fmy(:) = 0.0D0
        fdm(:,:,:) = 0.0D0
        CALL interpolate(lgrx,lgry,ux,uy,umx,umy)
        CALL lagr_force(fmx,fmy,ubx,uby,umx,umy)
        CALL spread(fmx,fmy,fdm,fd,lgrx,lgry,pi)
        CALL correct(ux,uy,fdm,rho)
     END DO
     CALL fluid_solid(f,fd)
     CALL update(rho,ux,uy,usq,f)
     CALL equilibrium(rho,ux,uy,usq,feq)
  END DO

  OPEN(UNIT=12,FILE="vfield")
  DO i = 0, xmax
     DO j = 0, ymax
        WRITE(12,*) i, j, SQRT(usq(i,j)), rho(i,j)/3.0D0
     END DO
  END DO
  CLOSE(UNIT=12)
  CLOSE(UNIT=13)
END PROGRAM lbm

!==========================!
! Create Lagrangian points !
!==========================!
SUBROUTINE lagrange(lgrx,lgry,pi,xc,yc)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(OUT) :: lgrx(1:nmrk), lgry(1:nmrk)
  REAL*8, INTENT(IN)  :: pi, xc, yc
  REAL*8 :: phi,d
  INTEGER :: j

  !create Lagrangian points
  phi = 2.0D0*(pi/nmrk)
  d = 2.0D0*rc*SIN(phi*0.5D0)
  PRINT *, 'd = ', d

  DO j = 0, nmrk-1
     lgrx(j+1) = (COS(j*phi)*rc) + xc
     lgry(j+1) = (SIN(j*phi)*rc) + yc
  END DO
END SUBROUTINE lagrange


!=========================!
! Initialise a simulation !
!=========================!
SUBROUTINE initialise(rho,ux,uy,usq,ubx,uby,ampl,freq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(OUT) :: rho(0:xmax,0:ymax)
  REAL*8, INTENT(OUT) :: ux(0:xmax,0:ymax), uy(0:xmax,0:ymax), usq(0:xmax,0:ymax)
  REAL*8, INTENT(OUT) :: ubx(1:nmrk), uby(1:nmrk)
  REAL*8, INTENT(IN)  :: freq, ampl

  rho(:,:) = 1.0D0
  ux(:,:)  = 0.0D0
  uy(:,:)  = 0.0D0
  usq(:,:) = 0.0D0
  ubx(:)   = 0.0D0
  uby(:)   = ampl*freq
END SUBROUTINE initialise

!==========================!
! Update macroscopic field !
!==========================!
SUBROUTINE update(rho,ux,uy,usq,f)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: rho(0:xmax,0:ymax), ux(0:xmax,0:ymax), uy(0:xmax,0:ymax), usq(0:xmax,0:ymax)
  REAL*8, INTENT(IN) :: f(0:xmax,0:ymax,0:8)
  INTEGER :: i,j

  DO i = 0, xmax
     DO j = 0, ymax
        rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
        ux(i,j)  = ((f(i,j,1)+f(i,j,5)+f(i,j,8)) - (f(i,j,3)+f(i,j,6)+f(i,j,7)))/rho(i,j)
        uy(i,j)  = ((f(i,j,2)+f(i,j,5)+f(i,j,6)) - (f(i,j,4)+f(i,j,7)+f(i,j,8)))/rho(i,j)
        usq(i,j) = ux(i,j)**2 + uy(i,j)**2
     END DO
  END DO
END SUBROUTINE update

!===================================!
! Interpolate u(x,t) ==> u^(m)(r,t) !
!===================================!
SUBROUTINE interpolate(lgrx,lgry,ux,uy,umx,umy)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: lgrx(1:nmrk), lgry(1:nmrk), ux(0:xmax,0:ymax), uy(0:xmax,0:ymax)
  REAL*8, INTENT(INOUT) :: umx(1:nmrk), umy(1:nmrk)
  INTEGER :: i,j,k
  REAL*8 :: xk,yk,x,y,kx,ky

  DO k = 1, nmrk
     xk = lgrx(k)
     yk = lgry(k)

     DO i = 0, xmax
        DO j = 0, ymax

           x = ABS(dble(i)-xk)
           y = ABS(dble(j)-yk)

           IF ((x>=0.0D0) .AND. (x<=dh) .AND. (y>=0.0D0) .AND. (y<=dh)) THEN
              kx = 1.0D0-x
              ky = 1.0D0-y
              umx(k) = umx(k) + ux(i,j)*kx*ky
              umy(k) = umy(k) + uy(i,j)*kx*ky
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE interpolate

!=======================================!
! Evaluate Lagrangian force on boundary !
!=======================================!
SUBROUTINE lagr_force(fmx,fmy,ubx,uby,umx,umy)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: ubx(1:nmrk), uby(1:nmrk), umx(1:nmrk), umy(1:nmrk)
  REAL*8, INTENT(INOUT) :: fmx(1:nmrk), fmy(1:nmrk)
  INTEGER :: k

  DO k = 1, nmrk
     fmx(k) = (1.0D0/dt)*(ubx(k)-umx(k))  !modified from 2
     fmy(k) = (1.0D0/dt)*(uby(k)-umy(k))
  END DO
END SUBROUTINE lagr_force

!==================================!
! Spread force to Eulerian lattice !
!==================================!
SUBROUTINE spread(fmx,fmy,fdm,fd,lgrx,lgry,pi)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: fmx(1:nmrk), fmy(1:nmrk), lgrx(1:nmrk), lgry(1:nmrk)
  REAL*8, INTENT(INOUT) :: fdm(0:xmax,0:ymax,1:2), fd(0:xmax,0:ymax,1:2),pi
  INTEGER :: i,j,k
  REAL*8 :: xk,yk,x,y,kx,ky,ds

  ds = 2.0D0*pi*rc/nmrk

  DO k = 1, nmrk
     xk = lgrx(k)
     yk = lgry(k)

     DO i = 0, xmax
        DO j = 0, ymax
           x = ABS(dble(i)-xk)
           y = ABS(dble(j)-yk)
           IF ((x>=0.0D0) .AND. (x<=dh) .AND. (y>=0.0D0) .AND. (y<=dh)) THEN
              kx = 1.0D0 - x
              ky = 1.0D0 - y
              fd(i,j,1)  = fd(i,j,1) + fmx(k)*kx*ky*ds
              fd(i,j,2)  = fd(i,j,2) + fmy(k)*kx*ky*ds
              fdm(i,j,1) = fdm(i,j,1) + fmx(k)*kx*ky*ds
              fdm(i,j,2) = fdm(i,j,2) + fmy(k)*kx*ky*ds
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE spread

!========================!
! Correct fluid velocity !
!========================!
SUBROUTINE correct(ux,uy,fdm,rho)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: ux(0:xmax,0:ymax), uy(0:xmax,0:ymax)
  REAL*8, INTENT(INOUT) :: fdm(0:xmax,0:ymax,1:2), rho(0:xmax,0:ymax)
  INTEGER :: i,j

  DO i = 0, xmax
     DO j = 0, ymax
        ux(i,j) = ux(i,j) + (fdm(i,j,1)*dt)/(2.0D0*rho(i,j))
        uy(i,j) = uy(i,j) + (fdm(i,j,2)*dt)/(2.0D0*rho(i,j))
     END DO
  END DO
END SUBROUTINE correct

!===========================!
! Relax towards equilibrium !
!===========================!
SUBROUTINE equilibrium(rho,ux,uy,usq,feq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: rho(0:xmax,0:ymax), ux(0:xmax,0:ymax), usq(0:xmax,0:ymax), uy(0:xmax,0:ymax)
  REAL*8, INTENT(INOUT) :: feq(0:xmax,0:ymax,0:8)
  INTEGER :: i,j
  REAL*8 :: uxy,A,B,C

  DO i = 0, xmax
     DO j = 0, ymax
        uxy = ux(i,j)*uy(i,j)
        A = 2.0D0*rho(i,j)/9.0D0
        B = rho(i,j)/18.0D0
        C = rho(i,j)/36.0D0

        feq(i,j,0) = A*(2.0D0 - 3.0D0*usq(i,j))

        feq(i,j,1) = B*(2.0D0 + 6.0D0*ux(i,j) + 9.0D0*ux(i,j)**2 - 3.0D0*usq(i,j))
        feq(i,j,2) = B*(2.0D0 + 6.0D0*uy(i,j) + 9.0D0*uy(i,j)**2 - 3.0D0*usq(i,j))
        feq(i,j,3) = B*(2.0D0 - 6.0D0*ux(i,j) + 9.0D0*ux(i,j)**2 - 3.0D0*usq(i,j))
        feq(i,j,4) = B*(2.0D0 - 6.0D0*uy(i,j) + 9.0D0*uy(i,j)**2 - 3.0D0*usq(i,j))

        feq(i,j,5) = C*(1.0D0 + 3.0D0*(ux(i,j)+uy(i,j)) + 9.0D0*uxy + 3.0D0*usq(i,j))
        feq(i,j,6) = C*(1.0D0 - 3.0D0*(ux(i,j)-uy(i,j)) - 9.0D0*uxy + 3.0D0*usq(i,j))
        feq(i,j,7) = C*(1.0D0 - 3.0D0*(ux(i,j)+uy(i,j)) + 9.0D0*uxy + 3.0D0*usq(i,j))
        feq(i,j,8) = C*(1.0D0 + 3.0D0*(ux(i,j)-uy(i,j)) - 9.0D0*uxy + 3.0D0*usq(i,j))
     END DO
  END DO
END SUBROUTINE equilibrium

!==============================!
! Add fluid-solid interactions !
!==============================!
SUBROUTINE fluid_solid(f,fd)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: f(0:xmax,0:ymax,0:8), fd(0:xmax,0:ymax,1:2)
  REAL*8 :: fx,fy
  INTEGER :: i,j

  DO i = 0, xmax
     DO j = 0, ymax
        fx = fd(i,j,1)
        fy = fd(i,j,2)

        f(i,j,1) = f(i,j,1) + (1.0D0/3.0D0)*fx
        f(i,j,2) = f(i,j,2) + (1.0D0/3.0D0)*fy
        f(i,j,3) = f(i,j,3) - (1.0D0/3.0D0)*fx
        f(i,j,4) = f(i,j,4) - (1.0D0/3.0D0)*fy
        f(i,j,5) = f(i,j,5) + (1.0D0/12.0D0)*( fx + fy)
        f(i,j,6) = f(i,j,6) + (1.0D0/12.0D0)*(-fx + fy)
        f(i,j,7) = f(i,j,7) + (1.0D0/12.0D0)*(-fx - fy)
        f(i,j,8) = f(i,j,8) + (1.0D0/12.0D0)*( fx - fy)
     END DO
  END DO
END SUBROUTINE fluid_solid

!===================!
! Perform collision !
!===================!
SUBROUTINE collision(f,feq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: f(0:xmax,0:ymax,0:8), feq(0:xmax,0:ymax,0:8)
  INTEGER :: i,j,k
  REAL*8 :: omeg
  
  omeg = -1.0D0/tau

  DO i = 0, xmax
     DO j = 0, ymax
        DO k = 0, 8
           f(i,j,k) = f(i,j,k) + omeg*(f(i,j,k) - feq(i,j,k))
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

  REAL*8, INTENT(INOUT) :: f(0:xmax,0:ymax,0:8)
  REAL*8 :: ftemp(0:xmax,0:ymax,0:8)

  !0:
  ftemp(:,:,0) = f(:,:,0)   

  !1:
  ftemp(1:xmax,:,1) = f(0:xmax-1,:,1)  
  ftemp(0,:,1) = f(xmax,:,1)

  !2:
  ftemp(:,1:ymax,2) = f(:,0:ymax-1,2) 
  ftemp(:,ymax,4) = f(:,ymax,2)

  !3:
  ftemp(0:xmax-1,:,3) = f(1:xmax,:,3) 
  ftemp(xmax,:,3) = f(0,:,3)

  !4:
  ftemp(:,0:ymax-1,4) = f(:,1:ymax,4) 
  ftemp(:,0,2) = f(:,0,4)

  !5:
  ftemp(1:xmax,1:ymax,5) = f(0:xmax-1,0:ymax-1,5) 
  ftemp(0,1:ymax,5) = f(xmax,0:ymax-1,5) 
  ftemp(:,ymax,7) = f(:,ymax,5)  

  !6:
  ftemp(0:xmax-1,1:ymax,6) = f(1:xmax,0:ymax-1,6)  
  ftemp(xmax,1:ymax,6) = f(0,0:ymax-1,6)   
  ftemp(:,ymax,8) = f(:,ymax,6)   

  !7:
  ftemp(0:xmax-1,0:ymax-1,7) = f(1:xmax,1:ymax,7) 
  ftemp(xmax,0:ymax-1,7) = f(0,1:ymax,7) 
  ftemp(:,0,5) = f(:,0,7) 

  !8:
  ftemp(1:xmax,0:ymax-1,8) = f(0:xmax-1,1:ymax,8)
  ftemp(0,0:ymax-1,8) = f(xmax,1:ymax,8)  
  ftemp(:,0,6) = f(:,0,8)

  f = ftemp

END SUBROUTINE stream

!===============!
! Move particle !
!===============!
SUBROUTINE move(yc,uby,lgry,k,ampl,freq)
  USE params
  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: uby(1:nmrk), lgry(1:nmrk)
  REAL*8, INTENT(INOUT) :: ampl, freq, yc
  INTEGER, INTENT(IN)   :: k
  INTEGER :: j
  REAL*8 :: dy, yc_n

  !update center of mass
  yc_n = yc_0 + ampl*SIN(freq*k)
  dy   = yc_n - yc

  !update Lagrangian points
  DO j = 1, nmrk
     lgry(j) = lgry(j) + dy
  END DO

  !update velocity on boundary markers
  DO j = 1, nmrk
     uby(j) = ampl*freq*COS(freq*k)
  END DO

  !reverse
  yc = yc_n
END SUBROUTINE move
