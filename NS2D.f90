! Navier Stokes Equation Solver for 2D Channel flow
! Finite Volume Method is used
! The System of Equations is solved using the Successive Over-Relaxation (SOR) Method
! A good method for faster convergence

! Need to comment/uncomment certain statements to change problem type
      
      PROGRAM NAVIERSTOKES

        implicit none
        integer, parameter :: nx=31, ny=31
        integer :: i,j,iter,outeriter,outeriterations
        real :: dx,dy,xmax,ymax,rho,mu,omegau,omegav,omegap,error,total, &
                xx,yy,velmag,min,mout
        real, dimension(0:nx+1,0:ny+1) :: u,v,p,uold,vold,pp, &
                       uu,vv,pressure,apu,apv,app,source,ae,aw,an,as
  
  !     set fluid properties
        mu=0.01 !viscosity
        rho=1. !density
  
        outeriter=1000   !number of outer iterations
  
  !     set relaxation factors
        omegau=0.7    !u-momentum underrelaxation
        omegav=0.7    !v-momentum underrelaxation
        omegap=0.3    !pressure update underrelaxation
  
  !     initialize variables for channel flow
        u=0.0; v=0.; p=0. ! velocities and pressure respectively

  !     set inlet and wall boundary conditions for channel flow
        u(:,0)=0.0    !zero on top wall
        u(:,ny+1)=0.0 !zero on lower wall
        u(1,:)=1.0    !one at inlet boundary
        u(nx+1,:)=u(1,:) !set outlet to inlet; conserve mass overall
  
        uold=u; vold=v    !initialize uold, vold
  
        xmax=1.           !set domain size  (increase length if channel)
        ymax=1.
        dx=xmax/float(nx) !compute mesh spacing
        dy=ymax/float(ny)
  
  !     initialize solver coefficients
        app=1.; apu=1.; apv=1.
        ae=0.; aw=0.; an=0.; as=0.
  
  !     begin outer iteration loop for sequential solution procedure
        OUTERLOOP: DO OUTERITERATIONS=1,OUTERITER
  
  !     X-MOMENTUM EQUATION SOLUTION
  !     set x-momentum equation coefficients
        do i=2,nx
            do j=1,ny
                ae(i,j)=max(-0.5*rho*dy*(uold(i+1,j)+ &
                           uold(i,j)),0.)+mu*dy/dx
                aw(i,j)=max( 0.5*rho*dy*(uold(i-1,j)+ &
                          uold(i,j)),0.)+mu*dy/dx
                an(i,j)=max(-0.5*rho*dx*(vold(i,j+1)+ &
                           vold(i-1,j+1)),0.0)+mu*dx/dy
                as(i,j)=max( 0.5*rho*dx*(vold(i,j)+ &
                           vold(i-1,j)),0.0)+mu*dx/dy
            end do
        end do
  
  !     overwrite boundary coefficients along north/south walls with half cell (dy) size
        do i=2,nx
            an(i,ny)=max(-0.5*rho*dx*(vold(i,ny+1)+ &
                           vold(i-1,ny+1)),0.0)+mu*dx/(dy/2.)
       
            as(i,1)=max( 0.5*rho*dx*(vold(i,1)+ &
                           vold(i-1,1)),0.0)+mu*dx/(dy/2.)
        end do
  
        apu=ae+aw+an+as
  
  !     Optional Coding to Block Out Cells
        apu(20:21,1:12)=1.e30
  
  
        apu=apu/omegau
  
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
        
        u(nx+1,:)=u(nx,:)  !for channel problem du/dx=0 at outlet
        !end do
  
  !     for channel problem to ensure overall mass conservation
        min=1.0   !analytical value
        mout=0.0

           do j=1,ny
             mout=mout+dy*rho*u(nx+1,j)
           end do
        u(nx+1,:)=u(nx+1,:)*min/mout !outflow mass correction for channel
  
  !      to generate fully developed flow in channel
        u(1,:)=u(nx+1,:)
  
  !*************************************************************
  
  !     Y-MOMENTUM EQUATION SOLUTION
  
        do i=1,nx
            do j=2,ny
                ae(i,j)=max(-0.5*rho*dy*(uold(i+1,j)+ &
                           uold(i+1,j-1)),0.)+mu*dy/dx
                aw(i,j)=max( 0.5*rho*dy*(uold(i,j-1)+ &
                           uold(i,j)),0.)+mu*dy/dx
  
                an(i,j)=max(-0.5*rho*dx*(vold(i,j+1)+ &
                           vold(i,j)),0.)+mu*dx/dy
                as(i,j)=max( 0.5*rho*dx*(vold(i,j-1)+ &
                           vold(i,j)),0.)+mu*dx/dy
            end do
        end do
  
  !     overwrite boundary coefficients along east/west walls due to half cell (dx) size
        do j=2,ny
          ae(nx,j)=max(-0.5*rho*dy*(uold(nx+1,j)+ &
                           uold(nx+1,j-1)),0.)+mu*dy/(dx/2.0)
          aw(1,j)=max( 0.5*rho*dy*(uold(1,j-1)+ &
                           uold(1,j)),0.)+mu*dy/(dx/2.0)
        end do
        
        apv=ae+aw+an+as
  
  !Optional Coding to Block out cells
        apv(21:21,1:13)=1.e30
  
        apv=apv/omegav
  
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
        
  
  !*************************************************************
  
  !      PRESSURE CORRECTION EQUATION
  
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
        app(1,1)=1.E30   !set reference cell value for pressure
        pp=0.0   !initialize corrections to zero
  
        source=0.0
  !     compute the mass source term
        do i=1,nx
            do j=1,ny
                source(i,j)=rho*dy*(u(i+1,j)-u(i,j)) + &
                          rho*dx*(v(i,j+1)-v(i,j))
            end do
        end do

  !     compute square root of sum of squares of mass imbalance
        total=sum(source**2)
        total=sqrt(total)
        print*,'outer iteration = ',outeriterations,'  mass imbalance = ',total
   
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
  !     conserved after corrections
        do i=1,nx
            do j=1,ny
                source(i,j)=rho*dy*(u(i+1,j)-u(i,j)) + &
                          rho*dx*(v(i,j+1)-v(i,j))
            end do
        end do
  
        total=sum(source**2)
        total=sqrt(total)
        print*,'outer iteration = ',outeriterations,'  mass imbalance = ',total
        print*
  
  !*************************************************************
  
        END DO OUTERLOOP
  
  
  !*************************************************************
  
  !     output results to .csv file for plotting with ParaView
        open(unit=50, file='NS_RESULTS.csv')
  
  !     interpolate velocities to scalar cell corners
        do i=2,nx
            do j=2,ny
                uu(i,j)=0.5*(u(i,j-1)+u(i,j))
                vv(i,j)=0.5*(v(i-1,j)+v(i,j))
                pressure(i,j)=0.25*(p(i-1,j-1)+p(i-1,j)+p(i,j-1)+p(i,j))
            end do
        end do
        
  !     east and west boundaries
        uu(1,2:ny)=0.5*(u(1,1:ny-1)+u(1,2:ny))
        uu(nx+1,2:ny)=0.5*(u(nx+1,1:ny-1)+u(nx+1,2:ny))
        vv(1,2:ny)=v(0,2:ny)
        vv(nx+1,2:ny)=v(nx+1,2:ny)
        pressure(1,2:ny)=0.5*(p(1,1:ny-1)+p(1,2:ny))
        pressure(nx+1,2:ny)=0.5*(p(nx,1:ny-1)+p(nx,2:ny))
  
  !     north and south boundaries
        uu(2:nx,1)=u(2:nx,0)
        uu(2:nx,ny+1)=u(2:nx,ny+1)
        vv(2:nx,1)=0.5*(v(1:nx-1,1)+v(2:nx,1))
        vv(2:nx,ny+1)=0.5*(v(1:nx-1,ny+1)+v(2:nx,ny+1))
        pressure(2:nx,1)=0.5*(p(1:nx-1,1)+p(2:nx,1))
        pressure(2:nx,ny+1)=0.5*(p(1:nx-1,ny)+p(2:nx,ny))
  
  !     southwest corner
        uu(1,1)=0
        vv(1,1)=0
        pressure(1,1)=pressure(2,2)
  
  !     southeast corner
        uu(nx+1,1)=0
        vv(nx+1,1)=0
        pressure(nx+1,1)=p(nx,1)
  
  !     northeast corner
        uu(nx+1,ny+1)=0 
        vv(nx+1,ny+1)=0
        pressure(nx+1,ny+1)=p(nx,ny)
  
  !     northwest corner
        uu(1,ny+1)=0
        vv(1,ny+1)=0
        pressure(1,ny+1)=p(1,ny)
  
  !      for blocked regions only
  !      uu(20:21,1:13)=0.0
  !      vv(20:21,1:13)=0.0
  
  !     write .csv file for ParaView
        write(50,*)'x, y, z, pressure, velmag, vx, vy, vz'
        do i=1,nx+1
          xx=float(i-1)*dx
        do j=1,ny+1
          yy=float(j-1)*dy
          velmag=sqrt(uu(i,j)**2+vv(i,j)**2)
          write(50,100)xx,yy,0.0,pressure(i,j),velmag,uu(i,j),vv(i,j),0.0
        end do
        end do
        
  100   format(f8.4,', ',f8.4,', ',f8.4,', ',f8.4,', ',f8.5,', ',f8.5,', ',f8.5,', ',f8.5)
  
        END PROGRAM NAVIERSTOKES
  
  
  
  
  
  
  
  