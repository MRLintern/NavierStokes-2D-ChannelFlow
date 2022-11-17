!Solves the Navier-Stokes Equations for 2D Open Channel Flow
!Discretization: Finite Volume Method
!Method for solving the resulting System of Algebraic Equations: Successive Over-Relaxation (SOR)
!************************************************************************************************
      module VARIABLES

      implicit none

      save

      integer, parameter :: nx=31, ny=31

      integer :: iteration, number_iterations, iter, iprob

      real :: dx,dy,xmax,ymax,rho,mu,omegau,omegav,omegap,error,total, &
              xx,yy,velmag,min,mout

      real, dimension(0:nx+1,0:ny+1) :: u,v,p,uold,vold,pp, &
                     uu,vv,pcorner,apu,apv,app,source,ae,aw,an,as    

      end module VARIABLES
!**********************************************************************
      program NAVIERSTOKES

      use VARIABLES

      implicit none

!     set fluid properties
      mu=0.01    !dynamic viscosity
      rho=1.     !density
      number_iterations=400  !number of outer iterations

!     set relaxation factors
      omegau=0.7    !u-momentum underrelaxation
      omegav=0.7    !v-momentum underrelaxation
      omegap=0.3    !pressure update underrelaxation

!*********************************
!     CHANNEL FLOW SETUP

      xmax=10.   !set domain size for channel
      ymax=1.
!     initialize variables for channel flow
      u=0.0; v=0.; p=0.
!     set inlet and wall boundary conditions for channel flow
      u(1,:)=1.0    !one at inlet boundary
      u(nx+1,:)=u(1,:) !set outlet to inlet; conserves mass overall
      u(:,0)=0.0      !lower wall
      u(:,ny+1)=0.0   !upper wall

!********************************

!     compute mesh spacing (keep aspect ratio close to 1)
      dx=xmax/float(nx)
      dy=ymax/float(ny)

      uold=u; vold=v    !initialize uold, vold
!     initialize solver coefficients
      app=1.; apu=1.; apv=1.
      ae=0.; aw=0.; an=0.; as=0.

!     begin solver outer loop
      OUTERLOOP : DO iteration=1,number_iterations

          print*,'iteration = ',iteration

          call umom
          call vmom
          call pressure

          print*

      end do OUTERLOOP

      call output

      end program NAVIERSTOKES

!**********************************************************************************
      subroutine umom 
      
      use VARIABLES
      
      implicit none

      integer :: i,j
      real :: mdote,mdotw,mdotn,mdots,residual

!     set x-momentum equation coefficients
      do i=2,nx
            do j=1,ny

                  mdote=0.5*rho*dy*(uold(i+1,j)+uold(i,j))
                  mdotw=0.5*rho*dy*(uold(i-1,j)+uold(i,j))
                  mdotn=0.5*rho*dx*(vold(i,j+1)+vold(i-1,j+1))
                  mdots=0.5*rho*dx*(vold(i,j)+vold(i-1,j))
                  ae(i,j)=MAX(-mdote,0.) + mu*dy/dx
                  aw(i,j)=MAX( mdotw,0.) + mu*dy/dx
                  an(i,j)=MAX(-mdotn,0.) + mu*dx/dy
                  as(i,j)=MAX( mdots,0.) + mu*dx/dy

            end do
      end do
      

!     overwrite boundary coefficients along north/south walls with half cell (dy) size
      do i=2,nx

        mdotn=0.5*rho*dx*(vold(i,ny+1) + vold(i-1,ny+1))
        mdots=0.5*rho*dx*(vold(i,1) + vold(i-1,1))
        an(i,ny)=MAX(-mdotn,0.0) + mu*dx/(dy/2.)
        as(i,1)=MAX( mdots,0.0) + mu*dx/(dy/2.)

      end do

      apu=ae+aw+an+as
      apu=apu/omegau   !build underrelaxation into apu

!     iterate x-momentum equations
      do iter=1,10
            do i=2,nx
                  do j=1,ny

                        u(i,j)=(1.-omegau)*uold(i,j)+1./apu(i,j)*( &
                        ae(i,j)*u(i+1,j)+ &
                        aw(i,j)*u(i-1,j)+ &
                        an(i,j)*u(i,j+1)+ &
                        as(i,j)*u(i,j-1)+ &
                        dy*(p(i-1,j)-p(i,j)))

                  end do
            end do
      end do

      u(nx+1,:)=u(nx,:)  !for channel problem only; du/dx=0 at outlet

!     ensure overall mass conservation
      min=1.0   !analytical value for channel flow rate
      mout=0.0

         do j=1,ny

           mout=mout+dy*rho*u(nx+1,j)

         end do

      u(nx+1,:)=u(nx+1,:)*min/mout !outflow mass correction for channel

!     compute x-momentum residual
      residual=0
      do i=2,nx
            do j=1,ny

                  residual=residual+(u(i,j)-uold(i,j))**2
            end do
      end do
      
      residual=sqrt(residual)

      print*,'x-momentum residual = ',residual

      return

      end subroutine umom

!*************************************************************
      subroutine vmom

      use VARIABLES

      implicit none

      integer :: i,j

      real :: mdote,mdotw,mdotn,mdots,residual

      DO i=1,nx
            DO j=2,ny

                  mdote=0.5*rho*dy*(uold(i+1,j) + uold(i+1,j-1))
                  mdotw=0.5*rho*dy*(uold(i,j-1) + uold(i,j))
                  mdotn=0.5*rho*dx*(vold(i,j+1) + vold(i,j))
                  mdots=0.5*rho*dx*(vold(i,j-1) + vold(i,j))
                  ae(i,j)=max(-mdote,0.) + mu*dy/dx
                  aw(i,j)=max( mdotw,0.) + mu*dy/dx
                  an(i,j)=max(-mdotn,0.) + mu*dx/dy
                  as(i,j)=max( mdots,0.) + mu*dx/dy

            end do
      end do
      

!     overwrite boundary coefficients along east/west walls due to half cell (dx) size
      do j=2,ny

        mdote=0.5*rho*dy*(uold(nx+1,j) + uold(nx+1,j-1))
        mdotw=0.5*rho*dy*(uold(1,j-1) +  uold(1,j))
        ae(nx,j)=max(-mdote,0.)+mu*dy/(dx/2.0)
        aw(1,j)=max( mdotw,0.)+mu*dy/(dx/2.0)

      end do
      
      apv=ae+aw+an+as
      apv=apv/omegav  !build underrelaxation into apv

      do iter=1,10
            do i=1,nx
                  do j=2,ny

                        v(i,j)=(1.-omegav)*vold(i,j)+1./apv(i,j)*( &
                        ae(i,j)*v(i+1,j)+ &
                        aw(i,j)*v(i-1,j)+ &
                        an(i,j)*v(i,j+1)+ &
                        as(i,j)*v(i,j-1)+ &
                        dx*(p(i,j-1)-p(i,j)))

                  end do
            end do
      end do
      
      
!     compute y-momentum residual
      residual=0

      do i=1,nx
            do j=2,ny

                  residual=residual+(v(i,j)-vold(i,j))**2

            end do
      end do
                  
      residual=sqrt(residual)

      print*,'y-momentum residual = ',residual

      return

      end subroutine vmom

!*************************************************************
      subroutine pressure

      use VARIABLES

      implicit none

      integer :: i,j

!     PRESSURE CORRECTION ROUTINE
!     set coefficients
      do i=1,nx
            do j=1,ny

                  ae(i,j)=rho*dy**2/apu(i+1,j)
                  aw(i,j)=rho*dy**2/apu(i,j)
                  an(i,j)=rho*dx**2/apv(i,j+1)
                  as(i,j)=rho*dx**2/apv(i,j)

            end do
      end do
                  
!     set boundary values for coefficients
      ae(nx,:)=0.0
      aw(1,:)=0.0
      an(:,ny)=0.0
      as(:,1)=0.0

      app=ae+aw+an+as
!     app(1,1)=1.E30   !optional; designates a cell for a reference pressure
      pp=0.0   !initialize pressure corrections to zero

      source=0.0
!     compute the mass source term before corrections
      do i=1,nx
            do j=1,ny

                  source(i,j)=rho*dy*(u(i+1,j)-u(i,j)) + &
                        rho*dx*(v(i,j+1)-v(i,j))
            end do
      end do

!     compute square root of sum of squares of mass imbalance
      total=sum(source**2)
      total=sqrt(total)

      print*,'mass imbalance = ',total
 
!     SOR iterations to solve for pressure correction pp
      do iter=1,100
            do j=1,ny
                  do i=1,nx

                        pp(i,j)=pp(i,j)+1.7/app(i,j)*( &
                        ae(i,j)*pp(i+1,j)+ &
                        aw(i,j)*pp(i-1,j)+ &
                        an(i,j)*pp(i,j+1)+ &
                        as(i,j)*pp(i,j-1)- &
                        source(i,j)-pp(i,j)*app(i,j))

                  end do
            end do
      end do
      
!     Apply corrections to pressure
      do i=1,nx
            do j=1,ny

                  p(i,j)=p(i,j)+omegap*pp(i,j)

            end do
      end do
      

!     Apply corrections to u velocity
      do i=2,nx
            do j=1,ny

                  u(i,j)=u(i,j)+dy/apu(i,j)* &
                        (pp(i-1,j)-pp(i,j))

            end do
      end do
      
!     Apply corrections to v velocity
      do i=1,nx
            do j=2,ny

                  v(i,j)=v(i,j)+dx/apv(i,j)* &
                        (pp(i,j-1)-pp(i,j))

            end do
      end do
      
!     update velocity variables
      uold=u   !note uold and vold are now mass conserving
      vold=v

!     recompute source term to check that mass is being
!     conserved after corrections; for debugging purposes

!     do i=1,nx
!           do j=1,ny
!                 source(i,j)=rho*dy*(u(i+1,j)-u(i,j)) + &
!                             rho*dx*(v(i,j+1)-v(i,j))
!           end do
!      end do
!      
!      total=sum(source**2)
!      total=sqrt(total)

!      print*,'outer iteration = ',iteration,'  mass imbalance = ',total
!      print*
      end subroutine pressure

!*************************************************************
      subroutine output

      use variables

      implicit none 

      integer :: i,j

!     open output file for plotting data using Paraview
      open(UNIT=50, FILE='RESULTS.dat')

!     interpolate velocities to scalar cell corners
      do i=2,nx
            do j=2,ny

                  uu(i,j)=0.5*(u(i,j-1)+u(i,j))
                  vv(i,j)=0.5*(v(i-1,j)+v(i,j))
                  pcorner(i,j)=0.25*(p(i-1,j-1)+p(i-1,j)+p(i,j-1)+p(i,j))

            end do
      end do
      
      
!     fill in east and west boundaries
      uu(1,2:ny)=0.5*(u(1,1:ny-1)+u(1,2:ny))
      uu(nx+1,2:ny)=0.5*(u(nx+1,1:ny-1)+u(nx+1,2:ny))
      vv(1,2:ny)=v(0,2:ny)
      vv(nx+1,2:ny)=v(nx+1,2:ny)
      pcorner(1,2:ny)=0.5*(p(1,1:ny-1)+p(1,2:ny))
      pcorner(nx+1,2:ny)=0.5*(p(nx,1:ny-1)+p(nx,2:ny))

!     fill in north and south boundaries
      uu(2:nx,1)=u(2:nx,0)
      uu(2:nx,ny+1)=u(2:nx,ny+1)
      vv(2:nx,1)=0.5*(v(1:nx-1,1)+v(2:nx,1))
      vv(2:nx,ny+1)=0.5*(v(1:nx-1,ny+1)+v(2:nx,ny+1))
      pcorner(2:nx,1)=0.5*(p(1:nx-1,1)+p(2:nx,1))
      pcorner(2:nx,ny+1)=0.5*(p(1:nx-1,ny)+p(2:nx,ny))

!     fill in southwest corner
      uu(1,1)=0
      vv(1,1)=0
      pcorner(1,1)=p(1,1)

!     fill in southeast corner
      uu(nx+1,1)=0
      vv(nx+1,1)=0
      pcorner(nx+1,1)=p(nx,1)

!     fill in northeast corner
      uu(nx+1,ny+1)=0 
      vv(nx+1,ny+1)=0
      pcorner(nx+1,ny+1)=p(nx,ny)

!    fill in northwest corner
      uu(1,ny+1)=0
      vv(1,ny+1)=0
      pcorner(1,ny+1)=p(1,ny)

!     write .csv file for ParaView
      write(50,*)'x, y, z, pressure, velmag, vx, vy, vz'

      do i=1,nx+1

        xx=float(i-1)*dx

      do j=1,ny+1

        yy=float(j-1)*dy
        velmag=sqrt(uu(i,j)**2+vv(i,j)**2)

        write(50,100)xx,yy,0.0,pcorner(i,j),velmag,uu(i,j),vv(i,j),0.0

      end do
      end do

100   format(f8.4,', ',f8.4,', ',f8.4,', ',f8.4,', ',f8.5,', ',f8.5,', ',f8.5,', ',f8.5)

      close(UNIT=50)

      end subroutine output
!**********************************************************************
