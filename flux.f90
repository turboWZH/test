module flux_module
  use mpi
  use mpifunc
  use mpished
  use datastructure_module
contains
  
  subroutine calc_centrifugal_source(grid)
    implicit none
    type(grid_t)          :: grid
    real(8),dimension(3)  :: flux
    integer               :: i
    do i = 1, grid%ninternal_nodes
       call calc_centrifugal_source_point(grid%omega,grid%x(1:3,i),&
            grid%dv(1,i),grid%dv(2:4,i),grid%vol(i),flux(1:3))
       grid%r(2:4,i)=grid%r(2:4,i)+flux(1:3)
    end do
  end subroutine calc_centrifugal_source

  subroutine calc_centrifugal_source_point(omega,x,rho,vel,vol,flux)
    implicit none
    real(8), dimension(3), intent(in)  :: omega,x,vel
    real(8),               intent(in)  :: rho,vol
    real(8), dimension(3), intent(out) :: flux
    real(8)                            :: volRho
    volRho=vol*rho
    flux(1)=volRho*(+omega(2)*vel(3)-omega(3)*vel(2))
    flux(2)=volRho*(-omega(1)*vel(3)+omega(3)*vel(1))
    flux(3)=volRho*(+omega(1)*vel(2)-omega(2)*vel(1))
  end subroutine calc_centrifugal_source_point

  subroutine calc_centrifugal_source_d(grid)
    implicit none
    type(grid_t)          :: grid
    integer:: i,j
    real*8:: flux_d(6,6),omega(3),vel(3),rho,vol
    do i = 1, grid%ninternal_nodes
       do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
          if(grid%JA_1sto(j).eq.i) then
             vol=grid%vol(i)
             rho=grid%dv(1,i)
             vel(1:3)=grid%dv(2:4,i)
             omega(1:3)=grid%omega(1:3)
             flux_d=0
             flux_d(1,2)=vol*( omega(2)*vel(3)-omega(3)*vel(2))
             flux_d(1,3)=vol*(-omega(1)*vel(3)+omega(3)*vel(1))
             flux_d(1,4)=vol*( omega(1)*vel(2)-omega(2)*vel(1))
             flux_d(2,3)=+omega(3)*vol*rho
             flux_d(2,4)=-omega(2)*vol*rho
             flux_d(3,2)=-omega(3)*vol*rho
             flux_d(3,4)=+omega(1)*vol*rho
             flux_d(4,2)=+omega(2)*vol*rho
             flux_d(4,3)=-omega(1)*vol*rho
             grid%jacobian_1sto(1:5,1:5,j)=grid%jacobian_1sto(1:5,1:5,j)+flux_d(1:5,1:5)
          end if
       end do
    end do
  end subroutine calc_centrifugal_source_d

  subroutine calc_interior_flux(grid)
    implicit none
    type(grid_t)                  :: grid
    integer                       :: nedge, i1, i2, i
    real(8)                       :: ds(3),dvL(6),dvR(6)
    real(8)                       :: xm(3),dx1(3),dx2(3)
    real(8)                       :: flux(6)
    real(8),dimension(grid%mnode) :: lim
    real*8                        :: URot(3),vn1(3),vn2(3),vn_SA,area,norm(3)
    real*8:: vel1(3),vel2(3),res1,res2

    if(grid%limiter) then
       do i=1,grid%mnode
          if (grid%TwoDim) then
             grid%lim(4,i)=1
          end if
          ! warning: simplified implementation, need to check
          lim(i) = min(grid%lim(1, i), grid%lim(2, i))
          lim(i) = min(lim(i),         grid%lim(3, i))
          lim(i) = min(lim(i),         grid%lim(4, i))
          lim(i) = min(lim(i),         grid%lim(5, i))
       end do

       do i=1,grid%mNode
          if(grid%dist(i).lt.0.1) then
             lim(i)=0
          elseif(grid%dist(i).lt.1) then
             lim(i)=log10(grid%dist(i)/0.1)
          else
             !lim(i)=1
          end if
       end do
    else
       do i=1,grid%mnode
          lim(i)=1.d0
       end do
    end if

    ! even when limiter is not used, full limiting is applied to vis wall nodes
    ! --------------------------------------------------------------------------
    if(grid%viscous) then
       do i=1,grid%mNode
          if(grid%dist(i).eq.0) then
              lim(i)=0.d0
          end if
       end do
    end if
    
    ! loop over edges
    ! ---------------
    do nedge=1,grid%medge
       i1=grid%edgeList(1,nedge)
       i2=grid%edgeList(2,nedge)

       if(i1.gt.grid%ninternal_nodes.and.&
            i2.gt.grid%ninternal_nodes) cycle
       ! this can be avoided if the edge loop
       ! exclude L2 edges
       ds=grid%edgeweight(1:3,nedge)
       area=sqrt(sum(ds**2))
       norm=ds(1:3)/area

       ! left and right before extrapolation
       ! -----------------------------------
       dvL(1:6)=grid%dv(1:6,i1)
       dvR(1:6)=grid%dv(1:6,i2)

       ! MUSCL reconstruction
       ! --------------------
       dx1=grid%x(1:3,i1)
       dx2=grid%x(1:3,i2)
       xm=0.5*(grid%x(1:3,i2)+grid%x(1:3,i1))          

       if(grid%order.eq.2) then
          dx1=xm-grid%x(1:3,i1)
          dx2=xm-grid%x(1:3,i2)          
          do i=1,5
             dvL(i)=dvL(i)+sum(dx1*grid%pgrad(1:3,i,i1))*lim(i1)
             dvR(i)=dvR(i)+sum(dx2*grid%pgrad(1:3,i,i2))*lim(i2)
          end do
          if(dvL(1).lt.0.0.or.dvR(1).lt.0.0.or.dvL(5).lt.0.0.or.dvR(5).lt.0.0) then
             print*,'Negative density or pressure in MUSCL. '
             write(*,'(A, I, A,  6ES11.4)'),'rank=',rank,' location=',grid%x(1:3,i1),grid%x(1:3,i2)

          end if
       end if

       ! compute convective flux
       ! ======================================================
       ! mass, momentum and energy
       ! ------------------------------------------------------
       URot=0.d0
       if(grid%RPM.gt.0) call calc_URot(grid%omega,xm,URot)
       call iflux_roe(dvL,dvR,URot,ds,grid%RoeEntropyFixCoeff,flux)

       ! sa convective flux
       ! Replace SA convective flux with the strictly 1stO one.
       ! This is not based on the linearly reconstructed value.
       ! ------------------------------------------------------
       if(grid%turbulent) then
          vn1=(grid%dv(2:4,i1)-Urot)*norm
          vn2=(grid%dv(2:4,i2)-Urot)*norm
          vn_SA=abs(sum(vn1+vn2))
          flux(6)=sum(vn1)*grid%dv(6,i1)+sum(vn2)*grid%dv(6,i2)-&
               vn_SA*(grid%dv(6,i2)-grid%dv(6,i1))
          flux(6)=flux(6)*0.5d0*area
       else
          flux(6)=0.d0
       end if
       ! ------------------------------------------------------
       grid%r(1:6,i1)=grid%r(1:6,i1)+flux
       grid%r(1:6,i2)=grid%r(1:6,i2)-flux

       ! compute convective flux
       ! ======================================================
       if(grid%viscous) then
          ! for viscous nodes, its 'physical' value is used for
          ! computing the viscous flux on an internal edge:
          ! the velocity is forced to be zero and then reverted
          ! after the vflux subroutine call.
          !----------------------------------------------------

          if(grid%dist(i1).eq.0.and.grid%dist(i2).ne.0) then
             vel1=grid%dv(2:4,i1)
             grid%dv(2:4,i1)=0
          end if
          if(grid%dist(i1).ne.0.and.grid%dist(i2).eq.0) then
             vel2=grid%dv(2:4,i2)
             grid%dv(2:4,i2)=0
          end if

          call vflux_FT_sa_neg(grid%dv(1:6,i1),grid%x(1:3,i1),grid%pgrad(1:3,1:7,i1),&
               grid%dv(1:6,i2),grid%x(1:3,i2),grid%pgrad(1:3,1:7,i2),&
               ds,grid%turbulent,grid%sutherland,flux,res1,res2)

          if(grid%dist(i1).eq.0.and.grid%dist(i2).ne.0) then
             grid%dv(2:4,i1)=vel1
          end if
          if(grid%dist(i1).ne.0.and.grid%dist(i2).eq.0) then
             grid%dv(2:4,i2)=vel2
          end if

          grid%r(1:6,i1)=grid%r(1:6,i1)+flux
          grid%r(1:6,i2)=grid%r(1:6,i2)-flux

       end if
    end do

  end subroutine calc_interior_flux

  subroutine calc_interior_flux_d(grid)
    implicit none
    type(grid_t)                 :: grid
    integer                      :: nedge,i1,i2,i,j
    real(8)                      :: ds(3),dvL(6),dvR(6)
    real(8)                      :: flux(6)

    !derivative seed
    real(8)::dvL_d(6,6),dvR_d(6,6),flux_d(6,6)
    real(8) URot(3),xm(3),vel1(3),vel2(3)

    ! loop over volume edges
    ! ----------------------
    do nedge=1,grid%medge
       i1=grid%edgeList(1,nedge)
       i2=grid%edgeList(2,nedge) 

       ds=grid%edgeweight(1:3,nedge)
       dvL(1:6)=grid%dv(1:6,i1)
       dvR(1:6)=grid%dv(1:6,i2)
       xm=0.5*(grid%x(1:3,i2)+grid%x(1:3,i1))
       URot=0.d0
       if(grid%RPM.gt.0) then
          call calc_URot(grid%omega,xm,URot)
       end if
       ! inviscid flux
       ! ---------------
       dvL_d=0
       dvR_d=0
       ! inviscid flux derivative wrt the left state
       ! -------------------------------------------
       do i=1,6
          dvL_d(i,i)=1.0
          dvR_d(i,i)=0.0
       end do

       flux_d=0.0
       call iflux_roe_dv(grid%dv(1:6,i1),dvL_d,&
            grid%dv(1:6,i2),dvR_d,URot,ds,&
            grid%RoeEntropyFixCoeff,flux,flux_d,6)

       do i=grid%IA_1sto(i1),grid%IA_1sto(i1+1)-1
          j=grid%JA_1sto(i)
          if(j.eq.i1) then
             grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)+flux_d(1:5,1:5)
          end if
       end do

       do i=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1
          j=grid%JA_1sto(i)
          if(j.eq.i1) then
             grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)-flux_d(1:5,1:5)
          end if
       end do

       if(grid%turbulent) grid%jaca(i1)=grid%jaca(i1)+flux_d(6,6)

       ! inviscid flux derivative wrt the right state
       ! --------------------------------------------
       dvL_d=0
       dvR_d=0
       do i=1,6
          dvL_d(i,i)=0.0
          dvR_d(i,i)=1.0
       end do

       flux_d=0.0
       call iflux_roe_dv(grid%dv(1:6,i1),dvL_d,&
            grid%dv(1:6,i2),dvR_d,URot,ds,&
            grid%RoeEntropyFixCoeff,flux,flux_d,6)

       do i=grid%IA_1sto(i1),grid%IA_1sto(i1+1)-1
          j=grid%JA_1sto(i)
          if(j.eq.i2) then
             grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)+flux_d(1:5,1:5)
          end if
       end do

       do i=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1
          j=grid%JA_1sto(i)
          if(j.eq.i2) then
             grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)-flux_d(1:5,1:5)
          end if
       end do

       if(grid%turbulent) grid%jaca(i2)=grid%jaca(i2)-flux_d(6,6)

       if(grid%viscous) then
          ! viscous flux derivative wrt the left state
          ! ------------------------------------------
          dvL_d=0
          dvR_d=0
          do i=1,6
             dvL_d(i,i)=1.0
             dvR_d(i,i)=0.0
          end do
          flux_d=0.0

          if(grid%dist(i1).eq.0) then
             vel1=grid%dv(2:4,i1)
             grid%dv(2:4,i1)=0
             dvL_d(2,2)=0
             dvL_d(3,3)=0
             dvL_d(4,4)=0
          end if

          if(grid%dist(i2).eq.0) then
             vel2=grid%dv(2:4,i2)
             grid%dv(2:4,i2)=0
             dvR_d(2,2)=0
             dvR_d(3,3)=0
             dvR_d(4,4)=0
          end if

          call vflux_FT_sa_neg_dv(grid%dv(1:6,i1),dvL_d,grid%x(1:3,i1),grid%pgrad(1:3,1:7,i1),&
               grid%dv(1:6,i2),dvR_d,grid%x(1:3,i2),grid%pgrad(1:3,1:7,i2),ds,&
               grid%turbulent,grid%sutherland,&
               flux,flux_d,6)

          if(grid%dist(i1).eq.0) then
             grid%dv(2:4,i1)=vel1
          end if

          if(grid%dist(i2).eq.0) then
             grid%dv(2:4,i2)=vel2
          end if

          do i=grid%IA_1sto(i1),grid%IA_1sto(i1+1)-1
             j=grid%JA_1sto(i)
             if(j.eq.i1) then
                grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)+flux_d(1:5,1:5)
             end if
          end do

          do i=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1
             j=grid%JA_1sto(i)
             if(j.eq.i1) then
                grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)-flux_d(1:5,1:5)
             end if
          end do

          if(grid%turbulent) grid%jaca(i1)=grid%jaca(i1)+flux_d(6,6)

          ! viscous flux derivative wrt the right state
          ! -------------------------------------------
          dvL_d=0
          dvR_d=0
          do i=1,6
             dvL_d(i,i)=0.0
             dvR_d(i,i)=1.0
          end do

          flux_d=0.0

          call vflux_ft_sa_neg_dv(grid%dv(1:6,i1),dvL_d,grid%x(1:3,i1),grid%pgrad(1:3,1:7,i1),&
               grid%dv(1:6,i2),dvR_d,grid%x(1:3,i2),grid%pgrad(1:3,1:7,i2),&
               ds,grid%turbulent,grid%sutherland,&
               flux,flux_d,6)

          do i=grid%IA_1sto(i1),grid%IA_1sto(i1+1)-1
             j=grid%JA_1sto(i)
             if(j.eq.i2) then
                grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)+flux_d(1:5,1:5)
             end if
          end do

          do i=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1
             j=grid%JA_1sto(i)
             if(j.eq.i2) then
                grid%jacobian_1sto(1:5,1:5,i)=grid%jacobian_1sto(1:5,1:5,i)-flux_d(1:5,1:5)
             end if
          end do

          if(grid%turbulent) grid%jaca(i2)=grid%jaca(i2)-flux_d(6,6)

       end if
    end do
  end subroutine calc_interior_flux_d

  subroutine calc_gradient_greenGauss(grid)
    implicit none
    type(grid_t)  grid
    integer :: i1,i2,nedge,k,i,j,nnode
    real(8) :: UR,UL,vds(3),vel1(3),vel2(3)

    grid%pgrad=0.0
    ! loop over interior edge
    ! -----------------------
    do nedge=1, grid%medge
       i1=grid%edgeList(1,nedge)
       i2=grid%edgeList(2,nedge)
       
       if(grid%viscous) then
          if(grid%dist(i1).eq.0) then
             vel1=grid%dv(2:4,i1)
             grid%dv(2:4,i1)=0
          end if
          if(grid%dist(i2).eq.0) then
             vel2=grid%dv(2:4,i2)
             grid%dv(2:4,i2)=0
          end if
       end if

       do k=1,6
          vds=0.5*(grid%dv(k,i2)-grid%dv(k,i1))*grid%edgeweight(1:3,nedge)
          grid%pgrad(1:3,k,i1)=grid%pgrad(1:3,k,i1)+vds
          grid%pgrad(1:3,k,i2)=grid%pgrad(1:3,k,i2)+vds
       end do

       if(grid%viscous) then
          if(grid%dist(i1).eq.0) then
             grid%dv(2:4,i1)=vel1
          end if
          if(grid%dist(i2).eq.0) then
             grid%dv(2:4,i2)=vel2
          end if
       end if

       UR=grid%dv(5,i2)/grid%dv(1,i2)/gas_constant
       UL=grid%dv(5,i1)/grid%dv(1,i1)/gas_constant
       vds=0.5*(UR-UL)*grid%edgeweight(1:3,nedge)
       grid%pgrad(1:3,7,i1)=grid%pgrad(1:3,7,i1)+vds
       grid%pgrad(1:3,7,i2)=grid%pgrad(1:3,7,i2)+vds
    end do

    ! periodic boundary
    ! -----------------
    if(grid%mPeriodicEdge >0) then
       do i=1,grid%mPeriodicEdge
          i1=grid%periodicEdgeList(1,i)
          i2=grid%periodicEdgeList(2,i)
          do j=1,7
             call rotateL2U(grid%pgrad(1:3,j,i1),grid%theta)
          end do
          do j=1,3
             call rotateL2U(grid%pgrad(j,2:4,i1),grid%theta)
          end do
          do j=1,7
             grid%pgrad(1:3,j,i2)=grid%pgrad(1:3,j,i2)+grid%pgrad(1:3,j,i1)
             grid%pgrad(1:3,j,i1)=grid%pgrad(1:3,j,i2)
          end do
          do j=1,3
             call rotateU2L(grid%pgrad(j,2:4,i1),grid%theta)
          end do
          do j=1,7
             call rotateU2L(grid%pgrad(1:3,j,i1),grid%theta)
          end do
       end do
    end if

    ! divided by cell volumes to get gradient
    ! ---------------------------------------
    do nNode=1,grid%mNode
       grid%pgrad(1:3,1:7,nNode)=grid%pgrad(1:3,1:7,nNode)/grid%vol(nNode)
    end do

    ! quasi2D case: zero the z-component of the gradient
    ! --------------------------------------------------
    if(grid%twoDim) then
       grid%pgrad(3,1:7,1:grid%mNode)=0
    end if

  end subroutine calc_gradient_greenGauss

  subroutine calc_gradient_leastSquare(grid)
    implicit none
    type(grid_t)  grid
    integer :: k,i,j,nnode,m,n,i1,i2,n1,n2,n3
    real(8) :: UR,UL,dvi(6),dvj(6)
    real*8 :: a11,a12,a21,a22    

    ! compute least square related coefficient
    ! ----------------------------------------
    if(grid%mPeriodicEdge>0) then
       a11= cos(grid%theta)
       a12=-sin(grid%theta)
       a21= sin(grid%theta)
       a22= cos(grid%theta)
    end if
    

    grid%pgrad=0.0
    ! loop over interior edge
    ! -----------------------
    ! compute the gradient of (rho, u, v, w, p, sa, T) using least squares

    do i=1,grid%mNode  
       n1=grid%NumberOfNeighbours(i)
       n2=grid%NumberOfNeighboursFromUpper(i)
       n3=grid%NumberOfNeighboursFromLower(i)  
       dvi=grid%dv(1:6,i)
       do m=1,n1+n2+n3 !grid%NumberOfNeighbours(i)
          if(m.le.n1) then
             j=grid%ListOfNeighbours(m,i)
             dvj=grid%dv(1:6,j)             
          elseif(m.le.n1+n2) then
             j=grid%ListOfNeighboursFromUpper(m-n1,i)
             dvj=grid%dv(1:6,j)             
             dvj(2)=a11*grid%dv(2,j)+a12*grid%dv(3,j)
             dvj(3)=a21*grid%dv(2,j)+a22*grid%dv(3,j)
             dvj(4)=    grid%dv(4,j)
          else
             j=grid%ListOfNeighboursFromLower(m-n1-n2,i)
             dvj=grid%dv(1:6,j)             
             dvj(2)=a11*grid%dv(2,j)+a21*grid%dv(3,j)
             dvj(3)=a12*grid%dv(2,j)+a22*grid%dv(3,j)
             dvj(4)=    grid%dv(4,j)
          end if


          !testing 20180813
          if(grid%viscous) then
             if(grid%dist(j).eq.0) then
                dvj(2:4)=0
             end if
             if(grid%dist(i).eq.0) then
                dvi(2:4)=0
             end if
          end if


          do k=1,6
             grid%pgrad(1,k,i)=grid%pgrad(1,k,i)+&
                  grid%wij1(m,i)*(dvj(k)-dvi(k))
             grid%pgrad(2,k,i)=grid%pgrad(2,k,i)+&
                  grid%wij2(m,i)*(dvj(k)-dvi(k))
             grid%pgrad(3,k,i)=grid%pgrad(3,k,i)+&
                  grid%wij3(m,i)*(dvj(k)-dvi(k))
          end do

          UR=dvj(5)/dvj(1)/gas_constant
          UL=dvi(5)/dvi(1)/gas_constant


          grid%pgrad(1,7,i)=grid%pgrad(1,7,i)+&
               grid%wij1(m,i)*(UR-UL)
          grid%pgrad(2,7,i)=grid%pgrad(2,7,i)+&
               grid%wij2(m,i)*(UR-UL)
          grid%pgrad(3,7,i)=grid%pgrad(3,7,i)+&
               grid%wij3(m,i)*(UR-UL)
       end do
    end do


  end subroutine calc_gradient_leastSquare

  subroutine calc_gradient_weightedLeastSquare(grid)
    implicit none
    type(grid_t)  grid
    integer :: k,i,j,nnode,m,n,i1,i2,n1,n2,n3
    real(8) :: UR,UL,dvi(6),dvj(6)
    real*8 :: a11,a12,a21,a22    
    real*8 :: weight
    ! compute least square related coefficient
    ! ----------------------------------------
    if(grid%mPeriodicEdge>0) then
       a11= cos(grid%theta)
       a12=-sin(grid%theta)
       a21= sin(grid%theta)
       a22= cos(grid%theta)
    end if
    

    grid%pgrad=0.0
    ! loop over interior edge
    ! -----------------------
    ! compute the gradient of (rho, u, v, w, p, sa, T) using least squares

    do i=1,grid%mNode  
       n1=grid%NumberOfNeighbours(i)
       n2=grid%NumberOfNeighboursFromUpper(i)
       n3=grid%NumberOfNeighboursFromLower(i)  
       dvi=grid%dv(1:6,i)
       do m=1,n1+n2+n3 !grid%NumberOfNeighbours(i)
          if(m.le.n1) then
             j=grid%ListOfNeighbours(m,i)
             dvj=grid%dv(1:6,j)             
          elseif(m.le.n1+n2) then
             j=grid%ListOfNeighboursFromUpper(m-n1,i)
             dvj=grid%dv(1:6,j)             
             dvj(2)=a11*grid%dv(2,j)+a12*grid%dv(3,j)
             dvj(3)=a21*grid%dv(2,j)+a22*grid%dv(3,j)
             dvj(4)=    grid%dv(4,j)
          else
             j=grid%ListOfNeighboursFromLower(m-n1-n2,i)
             dvj=grid%dv(1:6,j)             
             dvj(2)=a11*grid%dv(2,j)+a21*grid%dv(3,j)
             dvj(3)=a12*grid%dv(2,j)+a22*grid%dv(3,j)
             dvj(4)=    grid%dv(4,j)
          end if

          weight=1.d0/sqrt(sum((grid%x(1:3,j)-grid%x(1:3,i))**2))
          !!! for periodic nodes, need to be modified
          ! to reflect the information of rotation.

          !testing 20180813
          !blanking is used only for viscous wall nodes
          if(grid%viscous) then
             if(grid%dist(j).eq.0) then
                dvj(2:4)=0
             end if
             if(grid%dist(i).eq.0) then
                dvi(2:4)=0
             end if
          end if

          do k=1,6
             grid%pgrad(1,k,i)=grid%pgrad(1,k,i)+&
                  grid%wij1(m,i)*(dvj(k)-dvi(k))*weight
             grid%pgrad(2,k,i)=grid%pgrad(2,k,i)+&
                  grid%wij2(m,i)*(dvj(k)-dvi(k))*weight
             grid%pgrad(3,k,i)=grid%pgrad(3,k,i)+&
                  grid%wij3(m,i)*(dvj(k)-dvi(k))*weight
          end do

          UR=dvj(5)/dvj(1)/gas_constant
          UL=dvi(5)/dvi(1)/gas_constant


          grid%pgrad(1,7,i)=grid%pgrad(1,7,i)+&
               grid%wij1(m,i)*(UR-UL)*weight
          grid%pgrad(2,7,i)=grid%pgrad(2,7,i)+&
               grid%wij2(m,i)*(UR-UL)*weight
          grid%pgrad(3,7,i)=grid%pgrad(3,7,i)+&
               grid%wij3(m,i)*(UR-UL)*weight
       end do
    end do


  end subroutine calc_gradient_weightedLeastSquare



  subroutine rotateL2U(vec,theta)
    implicit none
    real*8::theta, vec(3)
    real*8::a(3,3),vec0(3)

    a=0
    a(1,1)= cos(theta)
    a(1,2)=-sin(theta)
    a(2,1)= sin(theta)
    a(2,2)= cos(theta)
    a(3,3)= 1.d0
    vec0  = vec
    vec(1)= sum(a(1,1:3)*vec0(1:3))
    vec(2)= sum(a(2,1:3)*vec0(1:3))
    vec(3)= sum(a(3,1:3)*vec0(1:3))
  end subroutine rotateL2U

  subroutine rotateU2L(vec,theta)
    implicit none
    real*8::theta, vec(3)
    real*8::a(3,3),vec0(3)

    a=0
    a(1,1)= cos(theta)
    a(1,2)= sin(theta)
    a(2,1)=-sin(theta)
    a(2,2)= cos(theta)
    a(3,3)= 1.d0
    vec0  = vec
    vec(1)= sum(a(1,1:3)*vec0(1:3))
    vec(2)= sum(a(2,1:3)*vec0(1:3))
    vec(3)= sum(a(3,1:3)*vec0(1:3))    
  end subroutine rotateU2L




  subroutine calc_sourceterm_uns(grid)
    implicit none
    type(grid_t)          :: grid
    integer               :: i    
    real*8:: cv(6),cv0(6)
    do i = 1, grid%ninternal_nodes
       cv(1)=grid%dv(1,i)
       cv(2:4)=grid%dv(2:4,i)*cv(1)
       cv(5)=grid%dv(5,i)/0.4+0.5*sum(cv(2:4)**2)/cv(1)
       cv(6)=grid%dv(6,i)

       cv0(1)=grid%dv_old(1,i)
       cv0(2:4)=grid%dv_old(2:4,i)*cv0(1)
       cv0(5)=grid%dv_old(5,i)/0.4+0.5*sum(cv0(2:4)**2)/cv0(1)
       cv0(6)=grid%dv_old(6,i)
       grid%r(1:6,i)=grid%r(1:6,i)+&
            (cv(1:6)-cv0(1:6))*grid%vol(i)/grid%timestep
    end do
  end subroutine calc_sourceterm_uns

  subroutine volume_scale_residual(grid)
    implicit none
    type(grid_t) :: grid
    integer::i,j
    do i=1,grid%ninternal_nodes
       grid%r(1:6,i)=grid%r(1:6,i)/grid%vol(i)
    end do
  end subroutine volume_scale_residual

  subroutine flow_scale_residual(grid)
    implicit none
    type(grid_t) :: grid
    integer::i,j  
    real*8:: velMag
    velMag=sqrt(sum(grid%dv_inf(2:4)**2))
    do i=2,4
       grid%r(i,1:grid%ninternal_nodes)=grid%r(i,1:grid%ninternal_nodes)/velMag
    end do
    grid%r(5,1:grid%ninternal_nodes)=grid%r(5,1:grid%ninternal_nodes)/grid%dv_inf(5)
    if(grid%turbulent) then
       grid%r(6,1:grid%ninternal_nodes)=grid%r(6,1:grid%ninternal_nodes)/grid%dv_inf(6)
    end if
  end subroutine flow_scale_residual



  pure function turbulentDynamicViscosity(rho,mul,nutSA) result(mut)
    implicit none
    real(8), intent(in)    :: rho,mul,nutSA
    real(8)                :: mut
    real(8)                :: chi,chi3,fv1
    chi=rho*nutSA/muL
    chi3=chi**3
    fv1=chi3/(chi3+c_v1**3)
    mut=rho*nutSA*fv1
  end function turbulentDynamicViscosity

  subroutine iflux_roe(UL,UR,URot,ds,entropyCorrectionCoeff,iflux)
    implicit none
    real*8, dimension(3), intent(in) :: ds
    real*8, dimension(3), intent(in) :: entropyCorrectionCoeff
    real*8, dimension(6), intent(in) :: UL,UR
    real*8, dimension(3), intent(in) :: URot
    real*8, dimension(6), intent(out):: iflux

    real*8 :: norm(3), area, &
         HL,HR,rhoVL,rhoVR,rtemp,rtemp1,rRoe,uRoe(3),&
         HRoe,VRoe,qRoe,cRoe,c2Roe,c2Roe2,dr,du(3),dp,&
         deltaV,lambda1,lambda2,lambda3,delta,&
         h1,h2,h3,h4,h5,fluxRoe(6),fluxC(6)

    iflux=0.0
    area=sqrt(sum(ds**2))
    norm=ds(1:3)/area

    ! total enthalpy 
    HL = gamma/gamm1*UL(5)/UL(1) + 0.5*sum(UL(2:4)**2)
    HR = gamma/gamm1*UR(5)/UR(1) + 0.5*sum(UR(2:4)**2)

    ! Roe-averaged variables
    rRoe    = sqrt(UL(1)*UR(1)) 
    rtemp   = rRoe / UL(1)
    rtemp1  = 1.0 / (1.0+rtemp)
    uRoe    = (UL(2:4)+UR(2:4)*rtemp)*rtemp1    
    HRoe    = (HL+HR*rtemp)*rtemp1 
    !   VRoe  = sum(uRoe(1:3)*norm(1:3))
    VRoe  = sum((uRoe(1:3)-URot)*norm(1:3))

    qRoe    = 0.5*sum(uRoe(1:3)**2)
    c2Roe   = gamm1*(HRoe-qRoe)
    cRoe    = sqrt(c2Roe)
    c2Roe2  = 2.0*c2Roe

    ! additional variables
    dr      = UR(1)-UL(1)
    du      = UR(2:4)-UL(2:4)
    dp      = UR(5)-UL(5)
    deltaV  = sum(du*norm)

    ! eigenvalues
    lambda1 = abs(VRoe-cRoe)   ! 1
    lambda2 = abs(VRoe)        ! 2,3,4
    lambda3 = abs(VRoe+cRoe)   ! 5
    ! entropy correction
    delta   = entropyCorrectionCoeff(1)*cRoe 
    if(lambda1.lt.delta) lambda1=0.5*(lambda1**2/delta+delta)
    delta   = entropyCorrectionCoeff(2)*cRoe 
    if(lambda2.lt.delta) lambda2=0.5*(lambda2**2/delta+delta)
    delta   = entropyCorrectionCoeff(3)*cRoe 
    if(lambda3.lt.delta) lambda3=0.5*(lambda3**2/delta+delta)

    ! Roe flux
    h1 = rRoe*cRoe*deltaV
    h2 = lambda1*(dp-h1)/c2Roe2
    h3 = lambda2*(dr-dp/c2Roe)
    h4 = lambda2*rRoe
    h5 = lambda3*(dp+h1)/c2Roe2

    fluxRoe     = 0.0
    fluxRoe(1)  = h2+h3+h5
    fluxRoe(2:4)= h2*(uRoe-cRoe*norm)+h3*uRoe+h4*(du-deltaV*norm)+h5*(uRoe+cRoe*norm)
    fluxRoe(5)  = h2*(HRoe-cRoe*VRoe)+h3*qRoe+h4*(sum(uRoe*du)-VRoe*deltaV)+h5*(HRoe+cRoe*VRoe)
    fluxRoe(6)  = abs(VRoe)*(UR(6)-UL(6))

    ! --- convective flux ---
    rhoVL = UL(1)*sum((UL(2:4)-UROT)*norm(1:3)) ! rho V - left
    rhoVR = UR(1)*sum((UR(2:4)-UROT)*norm(1:3)) ! rho V - right

    fluxC(1) = rhoVL          +rhoVR
    fluxC(2) = rhoVL*UL(2)    +rhoVR*UR(2)+(UL(5)+UR(5))*norm(1)
    fluxC(3) = rhoVL*UL(3)    +rhoVR*UR(3)+(UL(5)+UR(5))*norm(2)
    fluxC(4) = rhoVL*UL(4)    +rhoVR*UR(4)+(UL(5)+UR(5))*norm(3)
    fluxC(5) = rhoVL*HL       +rhoVR*HR+(UL(5)+UR(5))*sum(URot*norm)
    fluxC(6) = rhoVL*UL(6)/UL(1)+rhoVR*UR(6)/UR(1)
    iflux    = 0.5*area*(fluxC-fluxRoe)

  end subroutine iflux_roe
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.12 (r6213) -  9 May 2017 14:18
  !
  !  Differentiation of iflux_roe in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
  !   variations   of useful results: iflux
  !   with respect to varying inputs: ul ur
  !   RW status of diff variables: iflux:out ul:in ur:in
  SUBROUTINE IFLUX_ROE_DV(ul, uld, ur, urd, urot, ds, &
       & entropycorrectioncoeff, iflux, ifluxd, nbdirs)
    !  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    integer,parameter::nbdirsmax=6
    REAL*8, DIMENSION(3), INTENT(IN) :: ds
    REAL*8, DIMENSION(3), INTENT(IN) :: entropycorrectioncoeff
    REAL*8, DIMENSION(6), INTENT(IN) :: ul, ur
    REAL*8, DIMENSION(nbdirsmax, 6), INTENT(IN) :: uld, urd
    REAL*8, DIMENSION(3), INTENT(IN) :: urot
    REAL*8, DIMENSION(6), INTENT(OUT) :: iflux
    REAL*8, DIMENSION(nbdirsmax, 6), INTENT(OUT) :: ifluxd
    REAL*8 :: norm(3), area, hl, hr, rhovl, rhovr, rtemp, rtemp1, rroe, &
         & uroe(3), hroe, vroe, qroe, croe, c2roe, c2roe2, dr, du(3), dp, deltav&
         & , lambda1, lambda2, lambda3, delta, h1, h2, h3, h4, h5, fluxroe(6), &
         & fluxc(6)
    REAL*8 :: hld(nbdirsmax), hrd(nbdirsmax), rhovld(nbdirsmax), rhovrd(&
         & nbdirsmax), rtempd(nbdirsmax), rtemp1d(nbdirsmax), rroed(nbdirsmax), &
         & uroed(nbdirsmax, 3), hroed(nbdirsmax), vroed(nbdirsmax), qroed(&
         & nbdirsmax), croed(nbdirsmax), c2roed(nbdirsmax), c2roe2d(nbdirsmax), &
         & drd(nbdirsmax), dud(nbdirsmax, 3), dpd(nbdirsmax), deltavd(nbdirsmax)&
         & , lambda1d(nbdirsmax), lambda2d(nbdirsmax), lambda3d(nbdirsmax), &
         & deltad(nbdirsmax), h1d(nbdirsmax), h2d(nbdirsmax), h3d(nbdirsmax), h4d&
         & (nbdirsmax), h5d(nbdirsmax), fluxroed(nbdirsmax, 6), fluxcd(nbdirsmax&
         & , 6)
    INTRINSIC SUM
    INTRINSIC SQRT
    INTRINSIC ABS
    REAL*8 :: abs0
    REAL*8, DIMENSION(nbdirsmax) :: abs0d
    REAL*8, DIMENSION(3) :: arg1
    REAL*8, DIMENSION(nbdirsmax, 3) :: arg1d
    REAL*8 :: arg2
    REAL*8 :: arg10
    REAL*8, DIMENSION(nbdirsmax) :: arg10d
    INTEGER :: nd
    INTEGER :: nbdirs
    iflux = 0.0
    arg1(:) = ds**2
    arg2 = SUM(arg1(:))
    area = SQRT(arg2)
    norm = ds(1:3)/area
    arg1(:) = ul(2:4)**2
    hl = gamma/gamm1*ul(5)/ul(1) + 0.5*SUM(arg1(:))
    arg1(:) = ur(2:4)**2
    hr = gamma/gamm1*ur(5)/ur(1) + 0.5*SUM(arg1(:))
    arg10 = ul(1)*ur(1)
    rroe = SQRT(arg10)
    rtemp = rroe/ul(1)
    rtemp1 = 1.0/(1.0+rtemp)
    uroe = (ul(2:4)+ur(2:4)*rtemp)*rtemp1
    hroe = (hl+hr*rtemp)*rtemp1
    arg1(:) = (uroe(1:3)-urot)*norm(1:3)
    vroe = SUM(arg1(:))
    arg1(:) = uroe(1:3)**2
    qroe = 0.5*SUM(arg1(:))
    c2roe = gamm1*(hroe-qroe)
    DO nd=1,nbdirs
       ! total enthalpy 
       arg1d(nd, :) = 2*ul(2:4)*uld(nd, 2:4)
       hld(nd) = (gamma*uld(nd, 5)*ul(1)/gamm1-gamma*ul(5)*uld(nd, 1)/gamm1&
            &     )/ul(1)**2 + 0.5*SUM(arg1d(nd, :))
       arg1d(nd, :) = 2*ur(2:4)*urd(nd, 2:4)
       hrd(nd) = (gamma*urd(nd, 5)*ur(1)/gamm1-gamma*ur(5)*urd(nd, 1)/gamm1&
            &     )/ur(1)**2 + 0.5*SUM(arg1d(nd, :))
       ! Roe-averaged variables
       arg10d(nd) = uld(nd, 1)*ur(1) + ul(1)*urd(nd, 1)
       IF (arg10 .EQ. 0.0) THEN
          rroed(nd) = 0.0_8
       ELSE
          rroed(nd) = arg10d(nd)/(2.0*SQRT(arg10))
       END IF
       rtempd(nd) = (rroed(nd)*ul(1)-rroe*uld(nd, 1))/ul(1)**2
       rtemp1d(nd) = -(rtempd(nd)/(1.0+rtemp)**2)
       uroed(nd, :) = (uld(nd, 2:4)+urd(nd, 2:4)*rtemp+ur(2:4)*rtempd(nd))*&
            &     rtemp1 + (ul(2:4)+ur(2:4)*rtemp)*rtemp1d(nd)
       hroed(nd) = (hld(nd)+hrd(nd)*rtemp+hr*rtempd(nd))*rtemp1 + (hl+hr*&
            &     rtemp)*rtemp1d(nd)
       !   VRoe  = sum(uRoe(1:3)*norm(1:3))
       arg1d(nd, :) = norm(1:3)*uroed(nd, 1:3)
       vroed(nd) = SUM(arg1d(nd, :))
       arg1d(nd, :) = 2*uroe(1:3)*uroed(nd, 1:3)
       qroed(nd) = 0.5*SUM(arg1d(nd, :))
       c2roed(nd) = gamm1*(hroed(nd)-qroed(nd))
       IF (c2roe .EQ. 0.0) THEN
          croed(nd) = 0.0_8
       ELSE
          croed(nd) = c2roed(nd)/(2.0*SQRT(c2roe))
       END IF
       c2roe2d(nd) = 2.0*c2roed(nd)
       ! additional variables
       drd(nd) = urd(nd, 1) - uld(nd, 1)
       dud(nd, :) = urd(nd, 2:4) - uld(nd, 2:4)
       dpd(nd) = urd(nd, 5) - uld(nd, 5)
       arg1d(nd, :) = norm*dud(nd, :)
       deltavd(nd) = SUM(arg1d(nd, :))
    END DO
    croe = SQRT(c2roe)
    c2roe2 = 2.0*c2roe
    dr = ur(1) - ul(1)
    du = ur(2:4) - ul(2:4)
    dp = ur(5) - ul(5)
    arg1(:) = du*norm
    deltav = SUM(arg1(:))
    IF (vroe - croe .GE. 0.) THEN
       DO nd=1,nbdirs
          lambda1d(nd) = vroed(nd) - croed(nd)
       END DO
       lambda1 = vroe - croe
    ELSE
       DO nd=1,nbdirs
          lambda1d(nd) = -(vroed(nd)-croed(nd))
       END DO
       lambda1 = -(vroe-croe)
    END IF
    IF (vroe .GE. 0.) THEN
       DO nd=1,nbdirs
          lambda2d(nd) = vroed(nd)
       END DO
       lambda2 = vroe
    ELSE
       DO nd=1,nbdirs
          lambda2d(nd) = -vroed(nd)
       END DO
       lambda2 = -vroe
    END IF
    IF (vroe + croe .GE. 0.) THEN
       DO nd=1,nbdirs
          lambda3d(nd) = vroed(nd) + croed(nd)
       END DO
       lambda3 = vroe + croe
    ELSE
       DO nd=1,nbdirs
          lambda3d(nd) = -(vroed(nd)+croed(nd))
       END DO
       lambda3 = -(vroe+croe)
    END IF
    DO nd=1,nbdirs
       ! entropy correction
       deltad(nd) = entropycorrectioncoeff(1)*croed(nd)
    END DO
    delta = entropycorrectioncoeff(1)*croe
    IF (lambda1 .LT. delta) THEN
       DO nd=1,nbdirs
          lambda1d(nd) = 0.5*((2*lambda1*lambda1d(nd)*delta-lambda1**2*&
               &       deltad(nd))/delta**2+deltad(nd))
       END DO
       lambda1 = 0.5*(lambda1**2/delta+delta)
    END IF
    DO nd=1,nbdirs
       deltad(nd) = entropycorrectioncoeff(2)*croed(nd)
    END DO
    delta = entropycorrectioncoeff(2)*croe
    IF (lambda2 .LT. delta) THEN
       DO nd=1,nbdirs
          lambda2d(nd) = 0.5*((2*lambda2*lambda2d(nd)*delta-lambda2**2*&
               &       deltad(nd))/delta**2+deltad(nd))
       END DO
       lambda2 = 0.5*(lambda2**2/delta+delta)
    END IF
    DO nd=1,nbdirs
       deltad(nd) = entropycorrectioncoeff(3)*croed(nd)
    END DO
    delta = entropycorrectioncoeff(3)*croe
    IF (lambda3 .LT. delta) THEN
       DO nd=1,nbdirs
          lambda3d(nd) = 0.5*((2*lambda3*lambda3d(nd)*delta-lambda3**2*&
               &       deltad(nd))/delta**2+deltad(nd))
       END DO
       lambda3 = 0.5*(lambda3**2/delta+delta)
    END IF
    h1 = rroe*croe*deltav
    h2 = lambda1*(dp-h1)/c2roe2
    h3 = lambda2*(dr-dp/c2roe)
    h4 = lambda2*rroe
    h5 = lambda3*(dp+h1)/c2roe2
    arg1(:) = uroe*du
    DO nd=1,nbdirs
       ! Roe flux
       h1d(nd) = (rroed(nd)*croe+rroe*croed(nd))*deltav + rroe*croe*deltavd&
            &     (nd)
       h2d(nd) = ((lambda1d(nd)*(dp-h1)+lambda1*(dpd(nd)-h1d(nd)))*c2roe2-&
            &     lambda1*(dp-h1)*c2roe2d(nd))/c2roe2**2
       h3d(nd) = lambda2d(nd)*(dr-dp/c2roe) + lambda2*(drd(nd)-(dpd(nd)*&
            &     c2roe-dp*c2roed(nd))/c2roe**2)
       h4d(nd) = lambda2d(nd)*rroe + lambda2*rroed(nd)
       h5d(nd) = ((lambda3d(nd)*(dp+h1)+lambda3*(dpd(nd)+h1d(nd)))*c2roe2-&
            &     lambda3*(dp+h1)*c2roe2d(nd))/c2roe2**2
       fluxroed(nd, :) = 0.0_8
       fluxroed(nd, 1) = h2d(nd) + h3d(nd) + h5d(nd)
       fluxroed(nd, 2:4) = h2d(nd)*(uroe-croe*norm) + h2*(uroed(nd, :)-norm&
            &     *croed(nd)) + h3d(nd)*uroe + h3*uroed(nd, :) + h4d(nd)*(du-deltav*&
            &     norm) + h4*(dud(nd, :)-norm*deltavd(nd)) + h5d(nd)*(uroe+croe*norm&
            &     ) + h5*(uroed(nd, :)+norm*croed(nd))
       arg1d(nd, :) = uroed(nd, :)*du + uroe*dud(nd, :)
       fluxroed(nd, 5) = h2d(nd)*(hroe-croe*vroe) + h2*(hroed(nd)-croed(nd)&
            &     *vroe-croe*vroed(nd)) + h3d(nd)*qroe + h3*qroed(nd) + h4d(nd)*(SUM&
            &     (arg1(:))-vroe*deltav) + h4*(SUM(arg1d(nd, :))-vroed(nd)*deltav-&
            &     vroe*deltavd(nd)) + h5d(nd)*(hroe+croe*vroe) + h5*(hroed(nd)+croed&
            &     (nd)*vroe+croe*vroed(nd))
    END DO
    fluxroe = 0.0
    fluxroe(1) = h2 + h3 + h5
    fluxroe(2:4) = h2*(uroe-croe*norm) + h3*uroe + h4*(du-deltav*norm) + &
         &   h5*(uroe+croe*norm)
    fluxroe(5) = h2*(hroe-croe*vroe) + h3*qroe + h4*(SUM(arg1(:))-vroe*&
         &   deltav) + h5*(hroe+croe*vroe)
    IF (vroe .GE. 0.) THEN
       DO nd=1,nbdirs
          abs0d(nd) = vroed(nd)
       END DO
       abs0 = vroe
    ELSE
       DO nd=1,nbdirs
          abs0d(nd) = -vroed(nd)
       END DO
       abs0 = -vroe
    END IF
    arg1(:) = (ul(2:4)-urot)*norm(1:3)
    DO nd=1,nbdirs
       fluxroed(nd, 6) = abs0d(nd)*(ur(6)-ul(6)) + abs0*(urd(nd, 6)-uld(nd&
            &     , 6))
       ! --- convective flux ---
       ! rho V - left
       arg1d(nd, :) = norm(1:3)*uld(nd, 2:4)
       rhovld(nd) = uld(nd, 1)*SUM(arg1(:)) + ul(1)*SUM(arg1d(nd, :))
       ! rho V - right
       arg1d(nd, :) = norm(1:3)*urd(nd, 2:4)
       fluxcd(nd, :) = 0.0_8
    END DO
    fluxroe(6) = abs0*(ur(6)-ul(6))
    rhovl = ul(1)*SUM(arg1(:))
    arg1(:) = (ur(2:4)-urot)*norm(1:3)
    rhovr = ur(1)*SUM(arg1(:))
    DO nd=1,nbdirs
       rhovrd(nd) = urd(nd, 1)*SUM(arg1(:)) + ur(1)*SUM(arg1d(nd, :))
       fluxcd(nd, 1) = rhovld(nd) + rhovrd(nd)
       fluxcd(nd, 2) = rhovld(nd)*ul(2) + rhovl*uld(nd, 2) + rhovrd(nd)*ur(&
            &     2) + rhovr*urd(nd, 2) + norm(1)*(uld(nd, 5)+urd(nd, 5))
       fluxcd(nd, 3) = rhovld(nd)*ul(3) + rhovl*uld(nd, 3) + rhovrd(nd)*ur(&
            &     3) + rhovr*urd(nd, 3) + norm(2)*(uld(nd, 5)+urd(nd, 5))
       fluxcd(nd, 4) = rhovld(nd)*ul(4) + rhovl*uld(nd, 4) + rhovrd(nd)*ur(&
            &     4) + rhovr*urd(nd, 4) + norm(3)*(uld(nd, 5)+urd(nd, 5))
    END DO
    fluxc(1) = rhovl + rhovr
    fluxc(2) = rhovl*ul(2) + rhovr*ur(2) + (ul(5)+ur(5))*norm(1)
    fluxc(3) = rhovl*ul(3) + rhovr*ur(3) + (ul(5)+ur(5))*norm(2)
    fluxc(4) = rhovl*ul(4) + rhovr*ur(4) + (ul(5)+ur(5))*norm(3)
    arg1(:) = urot*norm
    DO nd=1,nbdirs
       fluxcd(nd, 5) = rhovld(nd)*hl + rhovl*hld(nd) + rhovrd(nd)*hr + &
            &     rhovr*hrd(nd) + SUM(arg1(:))*(uld(nd, 5)+urd(nd, 5))
       fluxcd(nd, 6) = ((rhovld(nd)*ul(6)+rhovl*uld(nd, 6))*ul(1)-rhovl*ul(&
            &     6)*uld(nd, 1))/ul(1)**2 + ((rhovrd(nd)*ur(6)+rhovr*urd(nd, 6))*ur(&
            &     1)-rhovr*ur(6)*urd(nd, 1))/ur(1)**2
       ifluxd(nd, :) = 0.5*area*(fluxcd(nd, :)-fluxroed(nd, :))
    END DO
    fluxc(5) = rhovl*hl + rhovr*hr + (ul(5)+ur(5))*SUM(arg1(:))
    fluxc(6) = rhovl*ul(6)/ul(1) + rhovr*ur(6)/ur(1)
    iflux = 0.5*area*(fluxc-fluxroe)
  END SUBROUTINE IFLUX_ROE_DV


  subroutine calc_URot(omega,x,URot)
    implicit none
    real*8, dimension(3), intent(in)    :: omega,x
    real*8, dimension(3), intent(out)   :: URot
    URot(1)=+omega(2)*x(3)-omega(3)*x(2)
    URot(2)=-omega(1)*x(3)+omega(3)*x(1)
    URot(3)=+omega(1)*x(2)-omega(2)*x(1)
  end subroutine calc_URot


  subroutine iflux_wall(U,x,ds,omega,flux)
    implicit none
    real*8, dimension(3), intent(in)    :: omega
    real*8, dimension(3), intent(in)    :: ds
    real*8, dimension(6), intent(in)    :: U
    real*8, dimension(3), intent(in)    :: x
    real*8, dimension(6), intent(out)   :: flux

    real(8), dimension(3)               :: URot

    URot(1)=+omega(2)*x(3)-omega(3)*x(2)
    URot(2)=-omega(1)*x(3)+omega(3)*x(1)
    URot(3)=+omega(1)*x(2)-omega(2)*x(1)

    flux(2:4)=ds*U(5)
    flux(5)=u(5)*sum(URot*ds)    
  end subroutine iflux_wall

  subroutine calc_boundary_flux(grid)
    implicit none
    type(grid_t) :: grid
    integer      :: nBdryGrp,nBnd,i,ii
    real(8)      :: flux(6),dvB(6),UROT(3),area,norm(3)
    real(8)      :: r_ref,p_ref,t_ref
    real(8)      :: height,dvel(3),radius,veltemp(3)
    logical      :: writeCF
    real(8)      :: qe_ref,temperature,mul,mut,mulL,mulR,mutL,mutR,chi,chi3,fv1,muf
    real(8)      :: force_x, force_y, force_x1, force_y1,  lift_force, drag_force, lift_force1, drag_force1, lift_force2, drag_force2

    !force_x = 0.0; force_x1 = 0.0
    !force_y = 0.0; force_y1 = 0.0
    !lift_force = 0.0; lift_force1 = 0.0; lift_force2 = 0.0
    !drag_force = 0.0; drag_force1 = 0.0; drag_force2 = 0.0

    T_ref=grid%totalTemperature_inlet
    P_ref=grid%totalPressure_inlet
    r_ref=P_ref/T_ref/Gas_Constant
    qe_ref=0.5*r_ref*grid%ref_u**2
    
    writeCf=.false.
    
    if(writeCf) then
       open(213,file='cf.dat')
    end if

    grid%reverseflow_clipping=0
    do nBdryGrp=1,grid%mBdryGrp
       ! skip periodic boundary groups
       if(grid%bdryGrpType(nBdryGrp).lt.0) cycle
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          area=sqrt(sum(grid%bdry2Weight(1:3,nBnd)**2))
          norm=grid%bdry2Weight(1:3,nBnd)/area

          call calc_dv2dvB(grid%bdryGrpType(nBdryGrp),&
               grid%x(1:3,i),&
               grid%omega(1:3),&
               grid%dv_inf(1:6),&
               grid%totalPressure_inlet,&
               grid%totalTemperature_inlet,&
               grid%backpressure_outlet,&
               grid%uDir(1:3),&
               grid%bdry2Weight(1:3,nBnd),&
               grid%dv(1:6,i),&
               dvB(1:6),&
               grid%reverseflow_clipping,grid%radius_c)

          URot=0
          if(grid%RPM.gt.0) call calc_URot(grid%omega,grid%x(1:3,i),URot)
          ! --- inviscid flux
          call iflux_roe(grid%dv(1:6,i),dvB(1:6),URot,&
               grid%bdry2Weight(1:3,nBnd),&
               grid%RoeEntropyFixCoeff,flux)
          
          grid%r(1:6,i)=grid%r(1:6,i)+flux
          ! --- viscous flux
          if(grid%viscous) then
             if(grid%bdryGrpType(nBdryGrp).eq.2 .or. &
                  grid%bdryGrpType(nBdryGrp).eq.3) then                
                ! testing proper weak no slip wall bc
                if(grid%bdryGrpType(nBdryGrp).eq.2) then
                   URot=0
                elseif(grid%bdryGrpType(nBdryGrp).eq.3) then
                   radius=grid%x(1,i)**2+grid%x(2,i)**2
                   radius=sqrt(radius)
                   ! --- set part of the bdry patch to be non-rotating
                   if(radius.gt.grid%Radius_c)then
                      URot=0
                   end if
                end if
                dvel=grid%dv(2:4,i)-URot
                ii=grid%nearwallnode(nBnd)
                ! if (ii.eq.0) goto 555
                
                height=grid%frictionHeight(nBnd)/4.0
                if(grid%turbulent) then
                   temperature=grid%dv(5,ii)/grid%dv(1,ii)/gas_constant          
                   if(grid%sutherland) then
                      mul = 0.00000145d0*temperature**1.5d0/(110.d0+temperature)
                   else
                      mul =airViscosity
                   end if
                   if(grid%dv(6,ii).le.0.d0) then
                      mut=0
                   else
                      chi   = grid%dv(1,ii)*grid%dv(6,ii)/mul
                      chi3  = chi**3
                      fv1   = chi3/(chi3+c_v1**3)
                      mut  = grid%dv(1,ii)*grid%dv(6,ii)*fv1       
                   end if
                   mulL=mul
                   mutL=mut

                   temperature=grid%dv(5,i)/grid%dv(1,i)/gas_constant          
                   if(grid%sutherland) then
                      mul = 0.00000145d0*temperature**1.5d0/(110.d0+temperature)
                   else
                      mul =airViscosity
                   end if
                   if(grid%dv(6,i).le.0.d0) then
                      mut=0
                   else
                      chi   = grid%dv(1,i)*grid%dv(6,i)/mul
                      chi3  = chi**3
                      fv1   = chi3/(chi3+c_v1**3)
                      mut  = grid%dv(1,i)*grid%dv(6,i)*fv1       
                   end if
                   mulR=mul
                   mutR=mut

                   muF=0.5*(mulL+mutL+mulR+mutR)
                   dvel=dvel/height*area*muF
                else
                   dvel=dvel/height*area*airviscosity
                end if


                ! warning: what viscosity should be used here?
                flux(2:4)=dvel
                if(writeCf) then
                   write(213,'(3ES14.5)') grid%x(1,i),(grid%dv(5,i)-101325)/qe_ref,&
                        sqrt(sum(flux(2:4)**2))/area/qe_ref
                   force_x = force_x + (grid%dv(5,i) - 101325)/qe_ref * grid%bdry2Weight(1,nBnd)
                   force_y = force_y + (grid%dv(5,i) - 101325)/qe_ref * grid%bdry2Weight(2,nBnd)
                   force_x1 = force_x1 + sqrt(sum(flux(2:4)**2))/area/qe_ref * grid%bdry2Weight(1,nBnd)
                   force_y1 = force_y1 + sqrt(sum(flux(2:4)**2))/area/qe_ref * grid%bdry2Weight(2,nBnd)
                end if
                flux(5)=0
                if(grid%bdryGrpType(nBdryGrp).eq.3) then
                   flux(5)=sum(dvel*URot)
                end if
                grid%r(2:5,i)=grid%r(2:5,i)+flux(2:5)
! 555             continue ! for corner nodes, dont compute shear stress
             end if
          end if
       end do
    end do
    if(writeCf) then       
       close(213)
       lift_force1 = -force_x * grid%UDir(2) + force_y * grid%UDir(1)
       lift_force2 = -force_x1 * grid%UDir(2) + force_y1 * grid%UDir(1)
       drag_force1 = force_x * grid%UDir(1) + force_y * grid%UDir(2)
       drag_force2 = force_x1 * grid%UDir(1) + force_y1 * grid%UDir(2)
       lift_force  = lift_force1 + lift_force2
       drag_force  = drag_force1 + drag_force2 
       write(*,*) 'Cl and Cd from cp distribution can be shown as follows'
       write(*,*) lift_force1, drag_force1
       write(*,*) 'Cl and Cd from cf distribution can be shown as follows'
       write(*,*) lift_force2, drag_force2
       write(*,*) 'Cl and Cd can be shown as follows'
       write(*,*) lift_force, drag_force
    end if
  end subroutine calc_boundary_flux

  subroutine calc_boundary_flux_d(grid)
    use mpished
    implicit none
    type(grid_t) :: grid
    integer      :: nBdryGrp,nBnd,i,j,i1,dummy_int,ii
    real(8)      :: flux(6),dvB(6)
    real*8:: flux_d(6,6),dv_d(6,6),dvB_d(6,6),URot(3),coeff,area,norm(3)


    real(8)      :: height,temperature,mul,mut,mulL,mulR,mutL,mutR,chi,chi3,fv1,muf


    do nBdryGrp=1,grid%mBdryGrp
       if(grid%bdryGrpType(nBdryGrp).lt.0) cycle
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          area=sqrt(sum(grid%bdry2Weight(1:3,nBnd)**2))
          norm=grid%bdry2Weight(1:3,nBnd)/area
          flux=0.0
          dvB=0.0         
          call calc_dv2dvB(grid%bdryGrpType(nBdryGrp),&
               grid%x(1:3,i),&
               grid%omega(1:3),&
               grid%dv_inf(1:6),&
               grid%totalPressure_inlet,grid%totalTemperature_inlet,grid%backpressure_outlet,&
               grid%uDir(1:3),grid%bdry2Weight(1:3,nBnd),grid%dv(1:6,i),dvB(1:6),&
               dummy_int,grid%radius_c)
          ! dBFlux/dv
          dv_d=0.0
          dvB_d=0.0
          do j=1,6
             dv_d(j,j)=1.0
          end do
          call calc_dv2dvB_dv(grid%bdryGrpType(nBdryGrp),&
               grid%x(1:3,i),&
               grid%omega(1:3),&
               grid%dv_inf(1:6),&
               grid%totalPressure_inlet,grid%totalTemperature_inlet,&
               grid%backpressure_outlet,grid%uDir(1:3),grid%bdry2Weight(1:3,nBnd),&
               grid%dv(1:6,i),dv_d,dvB(1:6),dvB_d,dummy_int,grid%radius_c,6)
          flux_d=0.0
          URot=0
          if(grid%RPM.gt.0) then
             call calc_URot(grid%omega,grid%x(1:3,i),URot)
          end if

          call iflux_roe_dv(grid%dv(1:6,i),dv_d,dvB(1:6),dvB_d,&
               URot,grid%bdry2Weight(1:3,nBnd),&
               grid%RoeEntropyFixCoeff,flux,flux_d,6)         


          coeff=0.d0
          if(grid%viscous) then
             if(grid%bdryGrpType(nBdryGrp).eq.2 .or. &
                  grid%bdryGrpType(nBdryGrp).eq.3) then
                ii=grid%nearwallnode(nBnd)
                height=grid%frictionHeight(nBnd)/4.0
                if(grid%turbulent) then
                   temperature=grid%dv(5,ii)/grid%dv(1,ii)/gas_constant          
                   if(grid%sutherland) then
                      mul = 0.00000145d0*temperature**1.5d0/(110.d0+temperature)
                   else
                      mul =airViscosity
                   end if
                   if(grid%dv(6,ii).le.0.d0) then
                      mut=0
                   else
                      chi   = grid%dv(1,ii)*grid%dv(6,ii)/mul
                      chi3  = chi**3
                      fv1   = chi3/(chi3+c_v1**3)
                      mut  = grid%dv(1,ii)*grid%dv(6,ii)*fv1       
                   end if
                   mulL=mul
                   mutL=mut
                   temperature=grid%dv(5,i)/grid%dv(1,i)/gas_constant          
                   if(grid%sutherland) then
                      mul = 0.00000145d0*temperature**1.5d0/(110.d0+temperature)
                   else
                      mul =airViscosity
                   end if
                   if(grid%dv(6,i).le.0.d0) then
                      mut=0
                   else
                      chi   = grid%dv(1,i)*grid%dv(6,i)/mul
                      chi3  = chi**3
                      fv1   = chi3/(chi3+c_v1**3)
                      mut  = grid%dv(1,i)*grid%dv(6,i)*fv1       
                   end if
                   mulR=mul
                   mutR=mut
                   muF=0.5*(mulL+mutL+mulR+mutR)
                   coeff=area*muF/height
                else
                   coeff=area*airviscosity/height
                end if
             end if
          end if

          do i1=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
             j=grid%JA_1sto(i1)
             if(j.eq.i) then
                grid%jacobian_1sto(1:5,1:5,i1)=grid%jacobian_1sto(1:5,1:5,i1)+flux_d(1:5,1:5)

                grid%jacobian_1sto(2,2,i1)=grid%jacobian_1sto(2,2,i1)+coeff*(1-norm(1)*norm(1))
                grid%jacobian_1sto(2,3,i1)=grid%jacobian_1sto(2,3,i1)+coeff*(  norm(1)*norm(2))
                grid%jacobian_1sto(2,4,i1)=grid%jacobian_1sto(2,4,i1)+coeff*(  norm(1)*norm(3))
                grid%jacobian_1sto(3,2,i1)=grid%jacobian_1sto(3,2,i1)+coeff*(  norm(2)*norm(1))
                grid%jacobian_1sto(3,3,i1)=grid%jacobian_1sto(3,3,i1)+coeff*(1-norm(2)*norm(2))
                grid%jacobian_1sto(3,4,i1)=grid%jacobian_1sto(3,4,i1)+coeff*(  norm(2)*norm(3))
                grid%jacobian_1sto(4,2,i1)=grid%jacobian_1sto(4,2,i1)+coeff*(  norm(3)*norm(1))
                grid%jacobian_1sto(4,3,i1)=grid%jacobian_1sto(4,3,i1)+coeff*(  norm(3)*norm(2))
                grid%jacobian_1sto(4,4,i1)=grid%jacobian_1sto(4,4,i1)+coeff*(1-norm(3)*norm(3))

             end if
          end do

          if(grid%turbulent) grid%jaca(i)=grid%jaca(i)+flux_d(6,6)

       end do
    end do
  end subroutine calc_boundary_flux_d


  subroutine calc_dv2dvB(bdryGrpType,x,omega,dv_inf,&
       totalpressure,totaltemperature,backpressure,&
       UDir,ds,dv,dvB,reverseFlow_clipping,radius_c)
    implicit none
    integer,             intent(in) :: bdryGrpType
    real(8),dimension(6),intent(in) :: dv,dv_inf
    real*8              ,intent(in) :: UDir(3),backpressure,totalPressure,totalTemperature,radius_c
    real(8),dimension(3),intent(in) :: ds,x,omega

    real(8),dimension(6),intent(out):: dvB
    real(8),dimension(3)            :: norm
    real*8                           ::radius
    real(8)                         :: area,vn

    real*8::a,b,c,d,Hti,sspdi,Rplus,cbplus,cbminus,cb,Ub,Mb,pb,Tb
    integer:: reverseFlow_clipping   
    area=sqrt(sum(ds**2))
    norm=ds/area

    if(bdryGrpType.eq.0) then                         ! farfield
       dvB(1:6)=dv_inf(1:6)
    elseif(bdryGrpType.eq.1) then                     ! inviscid slip wall
       dvB(1:6)=dv(1:6)
       dvB(2:4)=dv(2:4)-2*norm*sum(dv(2:4)*norm)
    elseif(bdryGrpType.eq.2.or.bdryGrpType.eq.3) then ! viscous no-slip wall
       dvB(1:6)=dv(1:6)
       dvB(2:4)=dv(2:4)-2*norm*sum(dv(2:4)*norm)
       if(bdryGrpType.eq.3) then                      ! rotating visc wall
          dvB(2)=+omega(2)*x(3)-omega(3)*x(2)
          dvB(3)=-omega(1)*x(3)+omega(3)*x(1)
          dvB(4)=+omega(1)*x(2)-omega(2)*x(1)
          radius=sqrt(x(1)**2+x(2)**2)
          if(radius.gt.radius_c) dvB(2:4)=0
       end if
    elseif(bdryGrpType.eq.4) then                     ! pressure outlet (back pressure)
       dvB(1:6)=dv(1:6)
       if(sum(dv(2:4)*ds)<0.0) then
          reverseFlow_clipping=reverseFlow_clipping+1
          dvB(2:4)=-dv(2:4)
          dvB(6)=dv_inf(6)
       end if
       dvB(5)=backpressure     
    elseif(bdryGrpType.eq.5) then                     ! pressure inlet (Blazek)
       Hti=dv(5)/dv(1)*gamma/gamm1+0.5*sum(dv(2:4)**2)
       sspdi=sqrt(gamma*dv(5)/dv(1))
       Rplus=-sum(dv(2:4)*norm)-2*sspdi/gamm1
       a=1+2/gamm1
       b=2*RPlus
       c=0.2d0*(Rplus**2-2.d0*Hti)     
       cbplus =-b/(2.d0*a)+sqrt(b**2-4.d0*a*c)/(2.d0*a)
       cbminus=-b/(2.d0*a)-sqrt(b**2-4.d0*a*c)/(2.d0*a)
       cb=cbplus
       Ub=+2.d0*cb/0.4d0+Rplus
       if(Ub.lt.0) Ub=0
       Mb=Ub/cb
       pb=totalPressure/(1.d0+0.2d0*Mb**2)**(gamma/gamm1)
       Tb=totalTemperature/(totalpressure/pb)**(gamm1/gamma)       
       dvB(1)=pb/Gas_Constant/Tb
       dvB(2:4)=-Ub*norm
       dvB(5)=pb       
       dvB(6)=dv_inf(6)
    end if

  end subroutine calc_dv2dvB

!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666) - 27 Nov 2017 17:10
!
!  Differentiation of calc_dv2dvb in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: dvb
!   with respect to varying inputs: dv
!   RW status of diff variables: dv:in dvb:out
SUBROUTINE CALC_DV2DVB_DV(bdrygrptype, x, omega, dv_inf, totalpressure, &
& totaltemperature, backpressure, udir, ds, dv, dvd, dvb, dvbd, &
& reverseflow_clipping, radius_c, nbdirs)
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
  integer,parameter::nbdirsmax=6
  INTEGER, INTENT(IN) :: bdrygrptype
  REAL*8, DIMENSION(6), INTENT(IN) :: dv, dv_inf
  REAL*8, DIMENSION(nbdirsmax, 6), INTENT(IN) :: dvd
  REAL*8, INTENT(IN) :: udir(3), backpressure, totalpressure, &
& totaltemperature, radius_c
  REAL*8, DIMENSION(3), INTENT(IN) :: ds, x, omega
  REAL*8, DIMENSION(6), INTENT(OUT) :: dvb
  REAL*8, DIMENSION(nbdirsmax, 6), INTENT(OUT) :: dvbd
  REAL*8, DIMENSION(3) :: norm
  REAL*8 :: radius
  REAL*8 :: area, vn
  REAL*8 :: a, b, c, d, hti, sspdi, rplus, cbplus, cbminus, cb, ub, mb, &
& pb, tb
  REAL*8, DIMENSION(nbdirsmax) :: bd, cd, htid, sspdid, rplusd, cbplusd&
& , cbd, ubd, mbd, pbd, tbd
! real(8)                :: cStagSq, cspeed, riemann, Vimag
! real(8)                :: alphaSq, cGhost, vGhost, MGhost
! real(8)                :: pGhost, TGhost, rhoGhost, fact
! real(8), dimension(3)  ::  ei
! real(8), parameter     :: GAMMA_R = gamma / gamm1
  REAL*8, PARAMETER :: small=1.0d-20
  INTEGER :: reverseflow_clipping
  INTRINSIC SUM
  INTRINSIC SQRT
  REAL*8, DIMENSION(3) :: arg1
  REAL*8, DIMENSION(nbdirsmax, 3) :: arg1d
  REAL*8 :: arg2
  REAL*8 :: arg10
  REAL*8, DIMENSION(nbdirsmax) :: arg10d
  REAL*8 :: result1
  REAL*8, DIMENSION(nbdirsmax) :: result1d
  REAL*8 :: pwx1
  REAL*8, DIMENSION(nbdirsmax) :: pwx1d
  REAL :: pwy1
  REAL*8 :: pwr1
  REAL*8, DIMENSION(nbdirsmax) :: pwr1d
  INTEGER :: nd
  INTEGER :: nbdirs
  arg1(:) = ds**2
  arg2 = SUM(arg1(:))
  area = SQRT(arg2)
  norm = ds/area
  IF (bdrygrptype .EQ. 0) THEN
! farfield
    dvb(1:6) = dv_inf(1:6)
    DO nd=1,nbdirs
      dvbd(nd, :) = 0.0_8
    END DO
  ELSE IF (bdrygrptype .EQ. 1) THEN
    DO nd=1,nbdirs
! inviscid slip wall
      dvbd(nd, 1:6) = dvd(nd, 1:6)
      arg1d(nd, :) = norm*dvd(nd, 2:4)
      dvbd(nd, 2:4) = dvd(nd, 2:4) - 2*norm*SUM(arg1d(nd, :))
    END DO
    dvb(1:6) = dv(1:6)
    arg1(:) = dv(2:4)*norm
    dvb(2:4) = dv(2:4) - 2*norm*SUM(arg1(:))
  ELSE IF (bdrygrptype .EQ. 2 .OR. bdrygrptype .EQ. 3) THEN
! viscous no-slip wall
    dvb(6) = 0
    dvb(2:4) = 0
    DO nd=1,nbdirs
      dvbd(nd, :) = 0.0_8
      dvbd(nd, 1) = dvd(nd, 1)
      dvbd(nd, 5) = dvd(nd, 5)
    END DO
    dvb(1) = dv(1)
    dvb(5) = dv(5)
    IF (bdrygrptype .EQ. 3) THEN
      DO nd=1,nbdirs
! rotating visc wall
        dvbd(nd, 2) = 0.0_8
        dvbd(nd, 3) = 0.0_8
        dvbd(nd, 4) = 0.0_8
      END DO
      dvb(2) = omega(2)*x(3) - omega(3)*x(2)
      dvb(3) = -(omega(1)*x(3)) + omega(3)*x(1)
      dvb(4) = omega(1)*x(2) - omega(2)*x(1)
!! sxu20180509 
      radius = x(1)**2 + x(2)**2
      radius = SQRT(radius)
      IF (radius .GT. radius_c) THEN
        DO nd=1,nbdirs
          dvbd(nd, 2:4) = 0.0_8
        END DO
        dvb(2:4) = 0
      END IF
    END IF
  ELSE IF (bdrygrptype .EQ. 4) THEN
    DO nd=1,nbdirs
! pressure outlet (back pressure)
      dvbd(nd, 1:6) = dvd(nd, 1:6)
    END DO
    dvb(1:6) = dv(1:6)
    arg1(:) = dv(2:4)*ds
    IF (SUM(arg1(:)) .LT. 0.0) THEN
      reverseflow_clipping = reverseflow_clipping + 1
      DO nd=1,nbdirs
        dvbd(nd, 2:4) = 0.0_8
      END DO
      dvb(2:4) = 0
    END IF
    DO nd=1,nbdirs
      dvbd(nd, 5) = 0.0_8
    END DO
    dvb(5) = backpressure
  ELSE IF (bdrygrptype .EQ. 5) THEN
    arg1(:) = dv(2:4)**2
    hti = dv(5)/dv(1)*gamma/gamm1 + 0.5*SUM(arg1(:))
    arg10 = gamma*dv(5)/dv(1)
    sspdi = SQRT(arg10)
    arg1(:) = dv(2:4)*norm
    rplus = -SUM(arg1(:)) - 2*sspdi/gamm1
    a = 1 + 2/gamm1
    b = 2*rplus
    DO nd=1,nbdirs
! pressure inlet (Blazek)
      arg1d(nd, :) = 2*dv(2:4)*dvd(nd, 2:4)
      htid(nd) = gamma*(dvd(nd, 5)*dv(1)-dv(5)*dvd(nd, 1))/dv(1)**2/&
&       gamm1 + 0.5*SUM(arg1d(nd, :))
      arg10d(nd) = (gamma*dvd(nd, 5)*dv(1)-gamma*dv(5)*dvd(nd, 1))/dv(1)&
&       **2
      IF (arg10 .EQ. 0.0) THEN
        sspdid(nd) = 0.0_8
      ELSE
        sspdid(nd) = arg10d(nd)/(2.0*SQRT(arg10))
      END IF
      arg1d(nd, :) = norm*dvd(nd, 2:4)
      rplusd(nd) = -SUM(arg1d(nd, :)) - 2*sspdid(nd)/gamm1
      bd(nd) = 2*rplusd(nd)
      cd(nd) = 0.2d0*(2*rplus*rplusd(nd)-2.d0*htid(nd))
      arg10d(nd) = 2*b*bd(nd) - 4.d0*a*cd(nd)
    END DO
    c = 0.2d0*(rplus**2-2.d0*hti)
    arg10 = b**2 - 4.d0*a*c
    DO nd=1,nbdirs
      IF (arg10 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = arg10d(nd)/(2.0*SQRT(arg10))
      END IF
      cbplusd(nd) = result1d(nd)/(2.d0*a) - bd(nd)/(2.d0*a)
      cbd(nd) = cbplusd(nd)
      ubd(nd) = 2.d0*cbd(nd)/0.4d0 + rplusd(nd)
    END DO
    result1 = SQRT(arg10)
    cbplus = -(b/(2.d0*a)) + result1/(2.d0*a)
    arg10 = b**2 - 4.d0*a*c
    result1 = SQRT(arg10)
    cbminus = -(b/(2.d0*a)) - result1/(2.d0*a)
    cb = cbplus
    ub = 2.d0*cb/0.4d0 + rplus
    IF (ub .LT. 0) THEN
      ub = 0
      DO nd=1,nbdirs
        ubd(nd) = 0.0_8
      END DO
    END IF
    mb = ub/cb
    pwx1 = 1.d0 + 0.2d0*mb**2
    pwy1 = gamma/gamm1
    pwr1 = pwx1**pwy1
    pb = totalpressure/pwr1
    DO nd=1,nbdirs
      mbd(nd) = (ubd(nd)*cb-ub*cbd(nd))/cb**2
      pwx1d(nd) = 0.2d0*2*mb*mbd(nd)
      IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. pwy1 .EQ. INT(pwy1))) &
&     THEN
        pwr1d(nd) = pwy1*pwx1**(pwy1-1)*pwx1d(nd)
      ELSE IF (pwx1 .EQ. 0.0 .AND. pwy1 .EQ. 1.0) THEN
        pwr1d(nd) = pwx1d(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
      pbd(nd) = -(totalpressure*pwr1d(nd)/pwr1**2)
      pwx1d(nd) = -(totalpressure*pbd(nd)/pb**2)
      dvbd(nd, :) = 0.0_8
    END DO
    pwx1 = totalpressure/pb
    pwy1 = gamm1/gamma
    pwr1 = pwx1**pwy1
    tb = totaltemperature/pwr1
    DO nd=1,nbdirs
      IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. pwy1 .EQ. INT(pwy1))) &
&     THEN
        pwr1d(nd) = pwy1*pwx1**(pwy1-1)*pwx1d(nd)
      ELSE IF (pwx1 .EQ. 0.0 .AND. pwy1 .EQ. 1.0) THEN
        pwr1d(nd) = pwx1d(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
      tbd(nd) = -(totaltemperature*pwr1d(nd)/pwr1**2)
      dvbd(nd, 1) = (pbd(nd)*tb/gas_constant-pb*tbd(nd)/gas_constant)/tb&
&       **2
      dvbd(nd, 2:4) = -(norm*ubd(nd))
      dvbd(nd, 5) = pbd(nd)
      dvbd(nd, 6) = 0.0_8
    END DO
    dvb(1) = pb/gas_constant/tb
    dvb(2:4) = -(ub*norm)
    dvb(5) = pb
    dvb(6) = dv_inf(6)
  ELSE
    DO nd=1,nbdirs
      dvbd(nd, :) = 0.0_8
    END DO
  END IF
END SUBROUTINE CALC_DV2DVB_DV


  subroutine bdry_setflow_periodic(grid)
    implicit none    
    type(grid_t) :: grid    
    real(8)      :: dvRot(3)
    integer      :: i,i2,i1
    real*8       :: a11,a12,a21,a22

    if(grid%mPeriodicEdge>0) then
       a11= cos(grid%theta)
       a12=-sin(grid%theta)
       a21= sin(grid%theta)
       a22= cos(grid%theta)
       do i=1,grid%mPeriodicEdge
          i1=grid%periodicEdgeList(1,i)
          i2=grid%periodicEdgeList(2,i)
          dvRot(1)=a11*grid%dv(2,i1)+a12*grid%dv(3,i1)
          dvRot(2)=a21*grid%dv(2,i1)+a22*grid%dv(3,i1)
          dvRot(3)=                      grid%dv(4,i1)
          grid%dv(1,i2)=grid%dv(1,i1)
          grid%dv(2,i2)=dvRot(1)
          grid%dv(3,i2)=dvRot(2)
          grid%dv(4,i2)=dvRot(3)
          grid%dv(5,i2)=grid%dv(5,i1)
          grid%dv(6,i2)=grid%dv(6,i1)
       end do
    end if
  end subroutine bdry_setflow_periodic



  subroutine calc_objective(grid)
    use mpifunc
    implicit none
    type(grid_t)          :: grid
    integer::nBdryGrp,nBnd,i
    real(8)::area,norm(3),tt,temperature,pt,sspd,spd,mach,massflux,pt_sum(6)

    grid%pt_inlet=0.
    grid%pt_outlet=0.
    grid%tt_inlet=0.
    grid%tt_outlet=0.
    grid%mass_inlet=0.
    grid%mass_outlet=0.
    grid%eta_t2t=0.

    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)

          if(i.gt. grid%ninternal_nodes) cycle

          area=sqrt(sum(grid%bdry2weight(1:3,nBnd)**2))

          norm=grid%bdry2weight(1:3,nBnd)/area                   
          if(grid%bdryGrpType(nBdryGrp).eq.5 .or. grid%bdryGrpType(nBdryGrp).eq.6 .or.grid%bdryGrpType(nBdryGrp).eq.0) then ! inlet
             sspd=sqrt(gamma*grid%dv(5,i)/grid%dv(1,i))
             spd=sqrt(sum(grid%dv(2:4,i)**2))
             mach=spd/sspd             
             pt=grid%dv(5,i)*(1+0.2d0*mach**2)**(gamma/gamm1)
             temperature=grid%dv(5,i)/grid%dv(1,i)/gas_constant             
             tt=temperature*(1+0.2d0*mach**2)             
             massflux=-grid%dv(1,i)*sum(grid%dv(2:4,i)*norm)*area
             grid%pt_inlet=grid%pt_inlet+pt*massflux
             grid%tt_inlet=grid%tt_inlet+tt*massflux
             grid%mass_inlet=grid%mass_inlet+massflux


          elseif(grid%bdryGrpType(nBdryGrp).eq.4) then ! outflow
             sspd=sqrt(gamma*grid%dv(5,i)/grid%dv(1,i))
             spd=sqrt(sum(grid%dv(2:4,i)**2))
             mach=spd/sspd
             pt=grid%dv(5,i)*(1+0.2d0*mach**2)**(gamma/gamm1)
             temperature=grid%dv(5,i)/grid%dv(1,i)/gas_constant             
             tt=temperature*(1+0.2d0*mach**2)             
             massflux=grid%dv(1,i)*sum(grid%dv(2:4,i)*norm)*area
             grid%pt_outlet=grid%pt_outlet+pt*massflux
             grid%tt_outlet=grid%tt_outlet+tt*massflux
             grid%mass_outlet=grid%mass_outlet+massflux



          end if
       end do
    end do

    pt_sum(1)        = grid%pt_inlet
    pt_sum(2)        = grid%mass_inlet
    pt_sum(3)        = grid%pt_outlet
    pt_sum(4)        = grid%mass_outlet
    pt_sum(5)        = grid%tt_inlet
    pt_sum(6)        = grid%tt_outlet

    call reduce_sum(6, pt_sum)

    grid%pt_inlet    = pt_sum(1)
    grid%mass_inlet  = pt_sum(2)
    grid%pt_outlet   = pt_sum(3)
    grid%mass_outlet = pt_sum(4)
    grid%tt_inlet    = pt_sum(5)
    grid%tt_outlet   = pt_sum(6)
    grid%pt_inlet    = grid%pt_inlet / grid%mass_inlet
    grid%pt_outlet   = grid%pt_outlet / grid%mass_outlet

    grid%tt_inlet    = grid%tt_inlet / grid%mass_inlet
    grid%tt_outlet   = grid%tt_outlet / grid%mass_outlet

    grid%eta_t2t=(grid%pt_outlet/grid%pt_inlet)**(gamm1/gamma)
    grid%eta_t2t=(grid%eta_t2t-1.d0)/&
         (grid%tt_outlet/grid%tt_inlet-1.d0)

  end subroutine calc_objective

  subroutine calc_backpressure_method1(grid)
    use mpifunc
    implicit none
    type(grid_t)          :: grid
    integer::nBdryGrp,nBnd,i
    real(8)::area,norm(3),massflux,pt_sum(1),relaxFactor
    grid%mass_outlet=0.d0
    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          if(i.gt. grid%ninternal_nodes) cycle
          area=sqrt(sum(grid%bdry2weight(1:3,nBnd)**2))
          norm=grid%bdry2weight(1:3,nBnd)/area                   
          if(grid%bdryGrpType(nBdryGrp).eq.4) then ! outflow
             massflux=grid%dv(1,i)*sum(grid%dv(2:4,i)*norm)*area
             grid%mass_outlet=grid%mass_outlet+massflux
          end if
       end do
    end do
    pt_sum(1) = grid%mass_outlet
    call reduce_sum(1, pt_sum)   
    if(grid%massflow.ne.0) then
       relaxFactor=1-0.02d0*(grid%massflow-pt_sum(1))/grid%massflow
       grid%backPressure_outlet=grid%backPressure_outlet*relaxFactor
    end if
  end subroutine calc_backpressure_method1

  subroutine calc_backpressure_method2(grid)
    use mpifunc
    implicit none
    type(grid_t)          :: grid
    integer::nBdryGrp,nBnd,i
    real(8)::area,norm(3),massflux,pt_sum(2),&
         mr,spd,sspd,mach,coeff,rho0,p0,c0,un,ut,&
         cc,p_sum,p,lambda,vel(3),totalArea

    grid%mass_outlet=0.d0
    totalArea=0.d0

    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          if(i.gt. grid%ninternal_nodes) cycle
          if(grid%bdryGrpType(nBdryGrp).eq.4) then ! outflow
             area=sqrt(sum(grid%bdry2weight(1:3,nBnd)**2))
             norm=grid%bdry2weight(1:3,nBnd)/area                   
             massflux=grid%dv(1,i)*sum(grid%dv(2:4,i)*norm)*area
             grid%mass_outlet=grid%mass_outlet+massflux
             totalArea=totalArea+area
          end if
       end do
    end do
    pt_sum(1) =grid%mass_outlet
    pt_sum(2)=totalArea
    call reduce_sum(2, pt_sum)   
    grid%mass_outlet=pt_sum(1)
    totalArea=pt_sum(2)

    p_sum=0

    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          if(i.gt. grid%ninternal_nodes) cycle
          if(grid%bdryGrpType(nBdryGrp).eq.4) then ! outflow
             area=sqrt(sum(grid%bdry2weight(1:3,nBnd)**2))
             norm=grid%bdry2weight(1:3,nBnd)/area                   
             massflux=grid%dv(1,i)*sum(grid%dv(2:4,i)*norm)*area        


             spd=sqrt(sum(grid%dv(2:4,i)**2))
             sspd=sqrt(gamma*grid%dv(5,i)/grid%dv(1,i))
             mach=spd/sspd
             coeff=1+(gamma-1)/2*mach**2
             rho0=grid%dv(1,i)*coeff**(    1/(gamma-1))
             p0  =grid%dv(5,i)*coeff**(gamma/(gamma-1))
             c0=sqrt(gamma*p0/rho0)
             un=sum(grid%dv(2:4,i)*norm)


             mr=massflux/area*grid%massflow/grid%mass_outlet

             vel=grid%dv(2:4,i)
             vel=vel-norm*sum(vel*norm)
             ut=sqrt(sum(vel**2))

             ! initial guess for un
             call calc_un(un,ut,rho0,c0,mr)

             cc=c0*sqrt(2/(gamma+1))

             lambda=sqrt(un**2+ut**2)/cc

             p=p0*(1-(gamma-1)/(gamma+1)*lambda**2)**(gamma/(gamma-1))
             p_sum=p_sum+p*area

          end if
       end do
    end do
    pt_sum(1)=p_sum
    call reduce_sum(1,pt_sum)
    p_sum=pt_sum(1)
    grid%backpressure_outlet=p_sum/totalArea

  end subroutine calc_backpressure_method2


  subroutine  calc_un(un,ut,r,c,mr)
    implicit none
    real*8::un,ut,u,r,c,mr,coeff,F,dF,coeff1
    integer::i

    do i=1,10
       u=sqrt(un**2+ut**2)
       coeff=1-(gamma-1)/2*(u/c)**2
       F=r*un*coeff**(1/(gamma-1))-mr
       dF=r*coeff**(1/(gamma-1))
       coeff1=1-(gamma-1)/2*(u/c)**2-(un/c)**2

       dF=dF*coeff1/coeff
       !       print*,i,F,dF,un
       un=un-F/dF
    end do
    !    stop

  end subroutine calc_un



  subroutine calc_total_pressue_loss(dv, grid, pt_loss)
    use mpifunc
    implicit none
    type(grid_t) :: grid
    real(8)      :: dv(6, grid%mNode), pt_loss
    integer      :: nBdryGrp, nBnd, i
    real(8)      :: area, norm(3), pt, sspd, spd, mach, massflux, pt_sum(4)

    grid%pt_inlet=0.
    grid%pt_outlet=0.
    grid%mass_inlet=0.
    grid%mass_outlet=0.
    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          i=grid%bdry2node(nBnd)
          if(i.gt. grid%ninternal_nodes) cycle
          area=sqrt(sum(grid%bdry2weight(1:3,nBnd)**2))

          norm=grid%bdry2weight(1:3,nBnd)/area                   
          if(grid%bdryGrpType(nBdryGrp).eq.5.or.grid%bdryGrpType(nBdryGrp).eq.0) then ! inlet
             sspd=sqrt(gamma*dv(5,i)/dv(1,i))
             spd=sqrt(sum(dv(2:4,i)**2))
             mach=spd/sspd             
             pt=dv(5,i)*(1+0.2d0*mach**2)**(gamma/gamm1)
             massflux=-dv(1,i)*sum(dv(2:4,i)*norm)*area
             grid%pt_inlet=grid%pt_inlet+pt*massflux
             grid%mass_inlet=grid%mass_inlet+massflux
          elseif(grid%bdryGrpType(nBdryGrp).eq.4) then ! outflow
             sspd=sqrt(gamma*dv(5,i)/dv(1,i))
             spd=sqrt(sum(dv(2:4,i)**2))
             mach=spd/sspd             
             pt=dv(5,i)*(1+0.2d0*mach**2)**(gamma/gamm1)
             massflux=dv(1,i)*sum(dv(2:4,i)*norm)*area
             grid%pt_outlet=grid%pt_outlet+pt*massflux
             grid%mass_outlet=grid%mass_outlet+massflux
          end if
       end do
    end do
    pt_sum(1)        = grid%pt_inlet
    pt_sum(2)        = grid%mass_inlet
    pt_sum(3)        = grid%pt_outlet
    pt_sum(4)        = grid%mass_outlet
    call reduce_sum(4, pt_sum)
    grid%pt_inlet    = pt_sum(1)
    grid%mass_inlet  = pt_sum(2)
    grid%pt_outlet   = pt_sum(3)
    grid%mass_outlet = pt_sum(4)
    grid%pt_inlet    = grid%pt_inlet / grid%mass_inlet
    grid%pt_outlet   = grid%pt_outlet / grid%mass_outlet
    pt_loss          = grid%pt_outlet - grid%pt_inlet
  end subroutine calc_total_pressue_loss


  subroutine bdry_setres(grid)
    implicit none
    type(grid_t) :: grid
    integer      :: nBdryGrp,nBnd,i    ,i1,i2,k,j

    ! viscous no slip wall, SA no update
    do nBdryGrp=1,grid%mBdryGrp
       if(grid%bdryGrpType(nBdryGrp).eq.2.or.grid%bdryGrpType(nBdryGrp).eq.3) then ! 
          do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
             i=grid%bdry2node(nBnd)
             grid%r(6,i)=0
          end do
       end if
    end do

    ! grid%r(6,1:grid%mNode)=0

    ! grid%r(1:5,1:grid%mNode)=0

    ! periodic bc, discard shadow periodic bdry
    if(grid%mPeriodicEdge>0) then
       do i=1,grid%mPeriodicEdge
          i1=grid%periodicEdgeList(1,i)
          i2=grid%periodicEdgeList(2,i)
          grid%r(1:6,i2)=0.0
       end do
    end if

  end subroutine bdry_setres



  subroutine bdry_setres_jacobian(grid)
    implicit none
    type(grid_t) :: grid
    integer      :: nBdryGrp,nBnd,i,j,k,i1,i2
    real(8)      :: blockJac(grid%mEqn,grid%mEqn),dt_coeff

    ! transpose each blockk in Jacobian (this is due to AD)    
    ! and add local time step to diagonal blocks.
    ! ----------------------------------------------------    





    do i=1,grid%mNode
       do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
          blockJac(1:grid%mEqn,1:grid%mEqn)=grid%jacobian_1sto(1:grid%mEqn,1:grid%mEqn,j)
          do k=1,5
             grid%jacobian_1sto(k,1:5,j)=blockJac(1:5,k)
          end do
          if(i.eq.grid%JA_1sto(j)) then
             do k=1,5
                grid%jacobian_1sto(k,k,j)=grid%jacobian_1sto(k,k,j)+grid%dt(i)/grid%cfl_global_flow
             end do
          end if
       end do
    end do



    ! turbulence time stepping augmented by local times step    
    ! ------------------------------------------------------
    do i=1,grid%mNode
       grid%dt_turb(i)=grid%jaca(i)+abs(grid%jacs(i))+grid%dt(i)/grid%cfl_global_sa
    end do

    do nBdryGrp=1,grid%mBdryGrp
       do nBnd=grid%bdry2nodeStart(nBdryGrp),grid%bdry2nodeStart(nBdryGrp+1)-1
          ! viscous wall correction
          ! -----------------------
          if(grid%bdryGrpType(nBdryGrp).eq.2.or.grid%bdryGrpType(nBdryGrp).eq.3) then
             i=grid%bdry2node(nBnd)         
             grid%dt_turb(i)=1.0
          end if
       end do
    end do

    do i=1,grid%mNode
       do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
          if(i.eq.grid%JA_1sto(j)) then
             grid%jacobian_1sto(6,6,j)=grid%dt_turb(i)
          end if
       end do
    end do


    if(grid%volumeScaling)then
       do i=1,grid%mNode
          do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
             grid%jacobian_1sto(1:6,1:6,j)=grid%jacobian_1sto(1:6,1:6,j)/grid%vol(i)
          end do
       end do
    end if



    if(grid%flowScaling)then
       do i=1,grid%mNode
          do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
             do k=2,4
                grid%jacobian_1sto(k,1:6,j)=grid%jacobian_1sto(k,1:6,j)/grid%dv_inf(2)
             end do
             k=5
             grid%jacobian_1sto(k,1:6,j)=grid%jacobian_1sto(k,1:6,j)/grid%dv_inf(5)
             if(grid%turbulent)then
                k=6
                grid%jacobian_1sto(k,1:6,j)=grid%jacobian_1sto(k,1:6,j)/grid%dv_inf(6)
             end if
          end do
       end do
    end if


    ! for laminar and euler, de-activate S-A update    
    if(.not. grid%turbulent) then
       do i=1,grid%mNode
          do j=grid%IA_1sto(i),grid%IA_1sto(i+1)-1
             grid%jacobian_1sto(1:6,6,j)=0
             grid%jacobian_1sto(6,1:6,j)=0
             if(i.eq.grid%JA_1sto(j)) then
                grid%jacobian_1sto(6,6,j)=1.0
             end if
          end do
       end do
    end if

    ! suppress the lower periodic boundary
    ! ------------------------------------
    do i=1,grid%mPeriodicEdge
       i1=grid%periodicEdgeList(1,i)
       i2=grid%periodicEdgeList(2,i)
       do j=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1
          if(grid%JA_1sto(j).eq.i2) then
             grid%jacobian_1sto(1:6,1:6,j)=0
             grid%jacobian_1sto(1,1,j)=1.d0
             grid%jacobian_1sto(2,2,j)=1.d0
             grid%jacobian_1sto(3,3,j)=1.d0
             grid%jacobian_1sto(4,4,j)=1.d0
             grid%jacobian_1sto(5,5,j)=1.d0
             grid%jacobian_1sto(6,6,j)=1.d0
          else
             grid%jacobian_1sto(1:6,1:6,j)=0
          end if
       end do
    end do

  end subroutine bdry_setres_jacobian


  subroutine per_accumulate(grid)
    implicit none
    type(grid_t) :: grid
    integer      :: i,i1,i2
    real*8 :: a11,a12,a21,a22    

    a11= cos(grid%theta)
    a12=-sin(grid%theta)
    a21= sin(grid%theta)
    a22= cos(grid%theta)


    do i=1,grid%mPeriodicEdge
       i1=grid%periodicEdgeList(1,i)
       i2=grid%periodicEdgeList(2,i)

       grid%r(1,i2)=grid%r(1,i2)+grid%r(1,i1)
       grid%r(2,i2)=grid%r(2,i2)+a11*grid%r(2,i1)+a12*grid%r(3,i1)
       grid%r(3,i2)=grid%r(3,i2)+a21*grid%r(2,i1)+a22*grid%r(3,i1)
       grid%r(4,i2)=grid%r(4,i2)+grid%r(4,i1)
       grid%r(5,i2)=grid%r(5,i2)+grid%r(5,i1)
       grid%r(6,i2)=grid%r(6,i2)+grid%r(6,i1)       

       grid%r(1,i1)=grid%r(1,i2)
       grid%r(2,i1)=a11*grid%r(2,i2)+a21*grid%r(3,i2)
       grid%r(3,i1)=a12*grid%r(2,i2)+a22*grid%r(3,i2)
       grid%r(4,i1)=grid%r(4,i2)
       grid%r(5,i1)=grid%r(5,i2)
       grid%r(6,i1)=grid%r(6,i2)

    end do
  end subroutine per_accumulate

  subroutine per_accumulate_jacobian(grid)
    implicit none
    type(grid_t) :: grid
    integer      :: i,i1,i2,j1,j2,m,n,n1,n2

    real*8::B(5,5),BJAC(5,5),BJACBInv(5,5),oldBJ(5,5),newBJ(5,5)

    BJAC=0
    BJACBINV=0
    B(1:5,1:5)=0
    B(1,1)=1
    B(2,2)= cos(grid%theta) 
    B(2,3)=-sin(grid%theta)
    B(3,2)= sin(grid%theta)
    B(3,3)= cos(grid%theta)
    B(4,4)=1
    B(5,5)=1


    do i=1,grid%mPeriodicEdge
       i1=grid%periodicEdgeList(1,i)
       i2=grid%periodicEdgeList(2,i)


       grid%jaca(i1)=grid%jaca(i1)+grid%jaca(i2)
       grid%jaca(i2)=grid%jaca(i1)       

       do j1=grid%IA_1sto(i1),grid%IA_1sto(i1+1)-1
          do j2=grid%IA_1sto(i2),grid%IA_1sto(i2+1)-1             
             if(grid%JA_1sto(j1).eq.i1.and.grid%JA_1sto(j2).eq.i2) then
                BJAC=0
                do n1=1,5
                   do n2=1,5
                      oldBJ(n1,n2)=grid%jacobian_1sto(n2,n1,j1)
                   end do
                end do

                do m=1,5
                   do n=1,5
                      BJAC(m,n)=sum(B(m,1:5)*oldBJ(1:5,n))
                   end do
                end do

                do m=1,5
                   do n=1,5
                      newBJ(n,m)=sum(BJAC(m,1:5)*B(n,1:5))
                   end do
                end do
                grid%jacobian_1sto(1:5,1:5,j2)=grid%jacobian_1sto(1:5,1:5,j2)+newBJ(1:5,1:5)

                BJAC=0

                do n1=1,5
                   do n2=1,5
                      oldBJ(n1,n2)=grid%jacobian_1sto(n2,n1,j2)
                   end do
                end do

                do m=1,5
                   do n=1,5
                      BJAC(m,n)=sum(B(1:5,m)*oldBJ(1:5,n))
                   end do
                end do

                do m=1,5
                   do n=1,5
                      newBJ(n,m)=sum(BJAC(m,1:5)*B(1:5,n))                      
                   end do
                end do

                grid%jacobian_1sto(1:5,1:5,j1)=newBJ(1:5,1:5)

             end if
          end do
       end do
    end do

  end subroutine per_accumulate_jacobian

  subroutine calc_limiter(grid)
    implicit none
    type(grid_t) :: grid
    real(8)      :: Umin(grid%mNode),Umax(grid%mNode)!local extrema
    real(8)      :: UminDomain,UmaxDomain            !global extrema    
    integer      :: k

    ! venkat limiter
    ! --------------
    do k=1,5
       call calc_limiter_extrema(grid,Umin,Umax,k)


       if(k.eq.1) then
          UminDomain=0.1
          UmaxDomain=2.0
       else if (k.eq.2 .or. k.eq.3.or. k.eq.4) then
          UminDomain=0
          UmaxDomain=500
       else 
          UminDomain=50000
          UmaxDomain=200000
       end if

       call calc_limiter_venkatDifferentiable(grid%mNode,grid%mEdge,grid%edgeList,&
            grid%x,grid%vol,grid%dv(k,1:grid%mNode),&
            grid%pgrad(1:3,k,1:grid%mNode),grid%lim(k,1:grid%mNode),&
            Umin,Umax,UminDomain,UmaxDomain,grid%limiterCoeff)
    end do

    ! ! Barth-Jesperpson limiter
    ! do k=1,y6
    !    call calc_limiter_extrema(grid,Umin,Umax,UminDomain,UmaxDomain,k)
    !    call calc_limiter_BJ(grid%mNode,grid%mEdge,grid%edgeList,&
    !         grid%x,grid%dv(k,1:grid%mNode),&
    !         grid%pgrad(1:3,k,1:grid%mNode),grid%lim(k,1:grid%mNode),&
    !         Umin,Umax)

    ! end do



  end subroutine calc_limiter

  subroutine calc_limiter_extrema(grid,Umin,Umax,k)
    ! Umin Umax: stencil-local min/max
    ! UminDomain UmaxDomain: global min/max
    implicit none
    type(grid_t)  :: grid
    real(8)       :: Umin(grid%mNode),Umax(grid%mNode)
!    real(8)       :: UminDomain,UmaxDomain
    integer       :: k,i,i1,i2


    ! compute node/neighbour min/max
    ! ------------------------------    
    call localExtrema(grid%mNode,grid%mEdge,grid%edgeList,&
         grid%dv(k,1:grid%mNode),Umin,Umax)

    ! sync extrema on periodic boundaries
    ! -----------------------------------
    do i=1,grid%mPeriodicEdge
       i1=grid%periodicEdgeList(1,i)
       i2=grid%periodicEdgeList(2,i)
       if(Umin(i1).gt.Umin(i2)) then
          Umin(i1)=Umin(i2)
       else
          Umin(i2)=Umin(i1)
       end if
       if(Umax(i1).lt.Umax(i2)) then
          Umax(i1)=Umax(i2)
       else
          Umax(i2)=Umax(i1)
       end if
    end do

    ! ! global min/max
    ! ! ---------------
    ! call globalExtrema( grid%mNode, grid%dv(k,1:grid%mNode), &
    !      UminDomain, UmaxDomain)

  end subroutine calc_limiter_extrema



  subroutine calc_limiter_BJ( mNode,mEdge,edgeList,x,&
       U,gradU,limiter,&    
       Umin,Umax)
    implicit none
    integer, intent(in)    :: mNode, mEdge
    integer, intent(in)    :: edgeList(2, mEdge)
    real(8), intent(in)    :: x(3, mNode)
    real(8), intent(in)    :: U(mNode)
    real(8), intent(in)    :: gradU(3,mNode)
    real(8), intent(out)   :: limiter(mNode)
  
    ! Local variables
    integer :: i, j, k, ij(3), iEdge
    real(8) :: Umin(mNode), Umax(mNode)
    real(8) :: dx(3)

    real*8::delta2

    limiter=1.d0

    do iEdge = 1, mEdge
       ij(1) = edgeList(1, iEdge)
       ij(2) = edgeList(2, iEdge)
       ij(3) = ij(1)

       do k=1,2
          i=ij(k)
          j=ij(k+1)
          dx=0.5*(x(1:3,j)-x(1:3,i))
          delta2 = sum(gradU(:, i) * dx)
          if(delta2.gt.1e-12) then
             limiter(i)=min(limiter(i),min(1.d0,(umax(i)-u(i))/delta2))
          else if(delta2.le.-1e-12) then
             limiter(i)=min(limiter(i),min(1.d0,(umin(i)-u(i))/delta2))
          end if
       end do
    end do



  end subroutine calc_limiter_BJ



  subroutine calc_limiter_venkatDifferentiable( mNode,mEdge,edgeList,x,volume,&
       U,gradU,limiter,&    
       Umin,Umax,UminDomain,UmaxDomain,limiterCoeff)
    implicit none
    integer, intent(in)    :: mNode, mEdge
    integer, intent(in)    :: edgeList(2, mEdge)
    real(8), intent(in)    :: x(3, mNode)
    real(8), intent(in)    :: volume(mNode), U(mNode)
    real(8), intent(in)    :: gradU(3,mNode)
    real(8), intent(out)   :: limiter(mNode)
    real(8), intent(in)    :: limiterCoeff
    ! Local variables
    integer :: i, j, k, ij(3), iEdge
    real(8) :: Umin(mNode), Umax(mNode)
    real(8) :: UminDomain, UmaxDomain
    real(8) :: dMaxMinSq, epSq, dPlus, dMinus,dMinus0
    real(8) :: dx(3), limiterHat,sigma 

    ! limiterCoeff is received from the input file
    epSq=(limiterCoeff*(UmaxDomain-UminDomain))**2.0




    ! Initialise limiter to 1.0
    limiter(1:mNode) = 1.0
    ! Loop over edges to get the minimum limiter value
    do iEdge = 1, mEdge
       ij(1) = edgeList(1, iEdge)
       ij(2) = edgeList(2, iEdge)
       ij(3) = ij(1)
       ! compute venkat limiter value
       do k=1,2
          i=ij(k)
          j=ij(k+1)
          dx=0.5*(x(1:3,j)-x(1:3,i))
          dMaxMinSq=(Umax(i)-Umin(i))**2        
          limiterHat = 1.0
          if (dMaxMinSq .gt. epSq) then
             dMinus0 = sum(gradU(:, i) * dx)
             if(dMinus0.lt.0) then
                dMinus=dMinus0 - 1.e-12
             else
                dMinus=dMinus0 + 1.e-12
             end if
             if (dMinus .ge. 1.e-12) then
                dPlus=Umax(i)-U(i)
                limiterHat=venkatLimiter(dPlus, dMinus, epSq)
             else if (dMinus .lt. -1.e-12) then
                dPlus=Umin(i)-U(i)
                limiterHat=venkatLimiter(dPlus, dMinus, epSq)
             end if
             ! additional bit proposed by Krzysztof Michalak and Carl Ollivier-Gooch
             ! in 'Limiters for unstructured higher-order accurate solutions of the
             ! Euler equations' AIAA paper 2008
             if (dMaxMinSq .lt. 2.*epSq ) then
                sigma  = (dMaxMinSq - epSq) / epSq
                sigma  = 2.0*sigma**3-3.0*sigma**2+1.0
                limiterHat = sigma+(1.0-sigma)*limiterHat
             end if
          end if
          limiter(i) = min(limiter(i), limiterHat)
       end do
    end do
  end subroutine calc_limiter_venkatDifferentiable

  subroutine localExtrema( mNode, mEdge, edgeList, U, Umin, Umax )
    implicit none
    integer, intent(in)    :: mNode, mEdge
    integer, intent(in)    :: edgeList(2,mEdge)
    real(8), intent(in)    :: U(mNode)
    real(8), intent(inout) :: Umin(mNode), Umax(mNode)
    integer :: i, j, ie
    Umin = U
    Umax = U
    do ie=1,mEdge
       i=edgeList(1, ie)
       j=edgeList(2, ie)
       Umin(i)=min(Umin(i),U(j))
       Umin(j)=min(Umin(j),U(i))
       Umax(i)=max(Umax(i),U(j))
       Umax(j)=max(Umax(j),U(i))
    end do
  end subroutine localExtrema


  subroutine globalExtrema(nNodes, U, UminDomain, UmaxDomain)
    use mpifunc
    implicit none
    integer,                    intent(in)    :: nNodes
    real(8), dimension(nNodes), intent(in)    :: U
    real(8),                    intent(inout) :: UmaxDomain, UminDomain
    integer :: i, ierr



    UmaxDomain = U(1)
    UminDomain = U(1)
    do i = 1, nNodes
       UmaxDomain = max(UmaxDomain, U(i))
       UminDomain = min(UminDomain, U(i))
    end do
    call reduce_min_max(UminDomain, UmaxDomain)
  end subroutine globalExtrema


  pure function venkatLimiter( dPlus, dMinus, epSq ) result(rvalue)
    implicit none
    real(8), intent(in)    :: dPlus, dMinus, epSq
    real(8) :: rvalue, dPlusSq, dMinusSq, Num, Den   
    dPlusSq  = dPlus**2
    dMinusSq = dMinus**2
    Num = (dPlusSq+epSq)*dMinus+2*dMinusSq*dPlus 
    Den = dPlusSq+2*dMinusSq+dMinus*dPlus+epSq
    rvalue = Num/(Den*dMinus)
  end function venkatLimiter


  subroutine distwall(grid)
    use reader
    use mpished
    use iorestart
    use mpifunc
    implicit none
    type(grid_t)                 :: grid
    integer                      :: nBdryGrp,nBnd,i
    real(8), allocatable         :: bxyz(:, :), temp(:)
    integer, allocatable         :: offset(:), tqoffset(:), bmesh(:,:)
    integer(hid_t)               :: hfile
    logical                      :: is_wall_patch(grid%mBdryGrp)

    ! allocate(grid%nearWallNode(grid%mNode))
    ! grid%nearWallNode=0

    call init_reader ()
    call open_mesh (grid%meshfile, hfile)
    ! Read node values from HDF5 mesh file
    ! xyz of boundary nodes
    call read_vector_real8 (get_link(BOUNDARY_MESH_X_COORDINATE), hfile, temp)
    allocate(bxyz(3, size(temp))); bxyz(1, :) = temp(:); deallocate(temp)
    call read_vector_real8 (get_link(BOUNDARY_MESH_Y_COORDINATE), hfile, temp)
    bxyz(2, :) = temp(:); deallocate(temp)
    call read_vector_real8 (get_link(BOUNDARY_MESH_Z_COORDINATE), hfile, temp)
    bxyz(3, :) = temp(:); deallocate(temp)
    call read_vector_int (get_link(SURFACE_ELEMENT_OFFSET), hfile, offset)
    call read_vector_int (get_link(SURFACE_ELEMENT_TRI_QUAD_OFFSET), hfile, tqoffset)
    call read_matrix_int (get_link(SURFACE_ELEMENT_NODES), hfile, bmesh)
    offset = offset + 1; tqoffset = tqoffset + 1
    call close_mesh (hfile)
    call close_reader ()

    ! determine the wall patches
    is_wall_patch = .false.
    do nBdryGrp = 1, grid%mBdryGrp
       if ( (grid%bdrygrptype(nBdryGrp).eq.2).or.  &
            (grid%bdrygrptype(nBdryGrp).eq.3)    ) &
            is_wall_patch(nBdryGrp) = .true. 
    end do

    grid%dist = HUGE(0.d0)
    do nBdryGrp = 1, grid%mBdryGrp
       if ( (grid%bdrygrptype(nBdryGrp).eq.2).or.  &
            (grid%bdrygrptype(nBdryGrp).eq.3)    ) then
          do nBnd = grid%bdry2nodeStart(nBdryGrp), &
               grid%bdry2nodeStart(nBdryGrp+1)-1
             i = grid%bdry2node(nBnd)
             grid%dist(i) = 0.d0
          end do
       end if
    end do

    ! wall distance by ray tracing
    call calc_wall_distance &
         (grid%mNode, size(bxyz,2), size(bmesh, 2), grid%mBdryGrp, &
         grid%mPeriodicEdge, is_wall_patch, offset, tqoffset,     &
         grid%periodicEdgeList, grid%z_trans_periodic_value,      &
         grid%theta, (/0.d0, 0.d0, 1.d0 /), bmesh, bxyz, grid%x,  &
         grid%dist )

    deallocate ( bxyz )
    deallocate ( offset, tqoffset )

    call blocking_exchange(1, grid%dist)

    ! write wall distance and its error to file
    call write_real_var ('WDIST', grid%mNode,  &
         grid%ninternal_nodes, &
         1, 1, grid%dist )


  end subroutine distwall


  subroutine distwall_restart (grid)
    use iorestart
    use mpifunc
    implicit none
    type(grid_t) :: grid

    call read_real_var ('WDIST', grid%mNode,  &
         grid%ninternal_nodes, &
         1, 1, grid%dist )
    call blocking_exchange(1, grid%dist)
  end subroutine distwall_restart





  subroutine isNaN(x,n1,n2)
    use mpished
    implicit none
    integer:: n1,n2
    real(8):: x(n1,n2)
    integer:: i,j
    do i=1,n1
       do j=1,n2
          if(x(i,j).ne.x(i,j)) then
             print*,'Residual NAN, rank ', rank,' node',j,' variable', i
             stop
          end if
       end do
    end do
  end subroutine isNaN

  subroutine calc_SA(grid)
    use mpished
    implicit none
    type(grid_t)                :: grid
    integer :: i
    real(8) :: source_sa0

    do i = 1, grid%ninternal_nodes
       if(grid%dist(i).ne.0.0) then     
          call SA_neg(grid%vol(i),grid%dv(1:6,i),grid%pgrad(1:3,2:4,i),&
               grid%pgrad(1:3,6,i),grid%dist(i),grid%sutherland,source_sa0)
          if(source_sa0.ne.source_sa0) then
             !...
          end if
          
          grid%r(6, i) = grid%r(6, i) + source_sa0

       end if
    end do

  end subroutine calc_SA


  subroutine calc_SA_d(grid)
    implicit none
    type(grid_t)                :: grid

    integer :: i,j
    real(8) ::source_sa
    real*8:: dv_d(6,6),flux_d(6)

    dv_d      = 0
    source_sa = 0

    do i = 1, grid%ninternal_nodes
       if(grid%dist(i).ne.0.0) then
          dv_d=0
          do j = 1, 6
             dv_d(j,j) = 1.0
          end do
          flux_d = 0
          call SA_neg_dv(grid%vol(i),grid%dv(1:6,i),dv_d,grid%pgrad(1:3,2:4,i),&
               grid%pgrad(1:3,6,i),grid%dist(i),grid%sutherland,source_sa,flux_d,6)

          ! call sourceTermSA_dv(grid%vol(i),grid%dv(1:6,i),dv_d,grid%pgrad(1:3,2:4,i),&
          !      grid%pgrad(1:3,6,i),grid%dist(i),grid%sutherland,source_sa,flux_d,6)
          grid%jacs(i)=flux_d(6)
       end if
    end do

  end subroutine calc_SA_d


  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.12 (r6213) -  9 May 2017 14:18
  !
  !  Differentiation of vflux_ft_sa_neg in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
  !   variations   of useful results: flux
  !   with respect to varying inputs: ul ur
  !   RW status of diff variables: flux:out ul:in ur:in
  SUBROUTINE VFLUX_FT_SA_NEG_DV(ul, uld, xl, gdl, ur, urd, xr, gdr, ds, &
       & turbulent, sutherland, flux, fluxd, nbdirs)
    !  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    integer,parameter:: nbdirsmax=6
    LOGICAL, INTENT(IN) :: turbulent, sutherland
    REAL*8, DIMENSION(3), INTENT(IN) :: xl, xr, ds
    REAL*8, DIMENSION(6), INTENT(IN) :: ul, ur
    REAL*8, DIMENSION(nbdirsmax, 6), INTENT(IN) :: uld, urd
    REAL*8, DIMENSION(3, 7), INTENT(IN) :: gdl, gdr
    REAL*8, DIMENSION(6), INTENT(OUT) :: flux
    REAL*8, DIMENSION(nbdirsmax, 6), INTENT(OUT) :: fluxd
    REAL*8 :: norm(3), area
    REAL*8 :: rhof, muf, nuf, dtdn, dnusadn
    REAL*8, DIMENSION(nbdirsmax) :: rhofd, mufd, nufd, dtdnd, dnusadnd
    REAL*8 :: theta(3), edgetangent(3), vf(3), dudn(3), gdtf(3), gdkf(3)
    REAL*8 :: thetad(nbdirsmax, 3), vfd(nbdirsmax, 3), dudnd(nbdirsmax, 3)&
         & , gdtfd(nbdirsmax, 3), gdkfd(nbdirsmax, 3)
    REAL*8 :: gduf(3, 3), tau(3, 3)
    REAL*8 :: gdufd(nbdirsmax, 3, 3), taud(nbdirsmax, 3, 3)
    ! lam/turb dynamic viscosity
    REAL*8 :: mul, mull, mulr, mut, mutl, mutr
    REAL*8, DIMENSION(nbdirsmax) :: muld, mulld, mulrd, mutd, mutld, mutrd
    REAL*8 :: kappal, kappat, kappa, dl, tl, tr
    REAL*8, DIMENSION(nbdirsmax) :: kappald, kappatd, kappad, tld, trd
    REAL*8 :: chi, chi3, fv1, ndne(3)
    REAL*8, DIMENSION(nbdirsmax) :: chid, chi3d, fv1d
    REAL*8 :: f_n_l, f_n_r
    REAL*8, DIMENSION(nbdirsmax) :: f_n_ld, f_n_rd
    INTRINSIC SUM
    INTRINSIC SQRT
    REAL*8, DIMENSION(3) :: arg1
    REAL*8, DIMENSION(nbdirsmax, 3) :: arg1d
    REAL*8 :: arg2
    INTEGER :: nd
    INTEGER :: nbdirs
    flux = 0
    arg1(:) = ds**2
    arg2 = SUM(arg1(:))
    area = SQRT(arg2)
    norm = ds/area
    DO nd=1,nbdirs
       tld(nd) = (uld(nd, 5)*ul(1)-ul(5)*uld(nd, 1))/ul(1)**2/gas_constant
       trd(nd) = (urd(nd, 5)*ur(1)-ur(5)*urd(nd, 1))/ur(1)**2/gas_constant
    END DO
    tl = ul(5)/ul(1)/gas_constant
    tr = ur(5)/ur(1)/gas_constant
    IF (sutherland) THEN
       DO nd=1,nbdirs
          mulld(nd) = (0.00000145d0*1.5d0*tl**0.5D0*tld(nd)*(110.d0+tl)-&
               &       0.00000145d0*tl**1.5d0*tld(nd))/(110.d0+tl)**2
          mulrd(nd) = (0.00000145d0*1.5d0*tr**0.5D0*trd(nd)*(110.d0+tr)-&
               &       0.00000145d0*tr**1.5d0*trd(nd))/(110.d0+tr)**2
       END DO
       mull = 0.00000145d0*tl**1.5d0/(110.d0+tl)
       mulr = 0.00000145d0*tr**1.5d0/(110.d0+tr)
    ELSE
       mull = airviscosity
       mulr = airviscosity
       DO nd=1,nbdirs
          mulld(nd) = 0.0_8
          mulrd(nd) = 0.0_8
       END DO
    END IF
    DO nd=1,nbdirs
       ! value on flux face
       muld(nd) = 0.5d0*(mulld(nd)+mulrd(nd))
       ! thermal conductivity coeff 
       kappald(nd) = cp*muld(nd)/prl
    END DO
    mul = 0.5d0*(mull+mulr)
    kappal = cp*mul/prl
    IF (turbulent) THEN
       chi = ul(1)*ul(6)/mull
       chi3 = chi**3
       fv1 = chi3/(chi3+c_v1**3)
       DO nd=1,nbdirs
          chid(nd) = ((uld(nd, 1)*ul(6)+ul(1)*uld(nd, 6))*mull-ul(1)*ul(6)*&
               &       mulld(nd))/mull**2
          chi3d(nd) = 3*chi**2*chid(nd)
          fv1d(nd) = (chi3d(nd)*(chi3+c_v1**3)-chi3*chi3d(nd))/(chi3+c_v1**3&
               &       )**2
          mutld(nd) = (uld(nd, 1)*fv1+ul(1)*fv1d(nd))*ul(6) + ul(1)*fv1*uld(&
               &       nd, 6)
       END DO
       mutl = ul(1)*ul(6)*fv1
       IF (ul(6) .LT. 0) THEN
          mutl = 0.d0
          DO nd=1,nbdirs
             mutld(nd) = 0.0_8
          END DO
       END IF
       chi = ur(1)*ur(6)/mulr
       chi3 = chi**3
       fv1 = chi3/(chi3+c_v1**3)
       DO nd=1,nbdirs
          chid(nd) = ((urd(nd, 1)*ur(6)+ur(1)*urd(nd, 6))*mulr-ur(1)*ur(6)*&
               &       mulrd(nd))/mulr**2
          chi3d(nd) = 3*chi**2*chid(nd)
          fv1d(nd) = (chi3d(nd)*(chi3+c_v1**3)-chi3*chi3d(nd))/(chi3+c_v1**3&
               &       )**2
          mutrd(nd) = (urd(nd, 1)*fv1+ur(1)*fv1d(nd))*ur(6) + ur(1)*fv1*urd(&
               &       nd, 6)
       END DO
       mutr = ur(1)*ur(6)*fv1
       IF (ur(6) .LT. 0) THEN
          mutr = 0.d0
          DO nd=1,nbdirs
             mutrd(nd) = 0.0_8
          END DO
       END IF
       DO nd=1,nbdirs
          mutd(nd) = 0.5*(mutld(nd)+mutrd(nd))
          kappatd(nd) = cp*mutd(nd)/prt
       END DO
       mut = 0.5*(mutl+mutr)
       kappat = cp*mut/prt
    ELSE
       DO nd=1,nbdirs
          kappatd(nd) = 0.0_8
          mutd(nd) = 0.0_8
       END DO
    END IF
    DO nd=1,nbdirs
       ! --- Face values
       vfd(nd, :) = 0.5*(uld(nd, 2:4)+urd(nd, 2:4))
       rhofd(nd) = 0.5*(uld(nd, 1)+urd(nd, 1))
    END DO
    vf = 0.5*(ul(2:4)+ur(2:4))
    rhof = 0.5*(ul(1)+ur(1))
    IF (turbulent) THEN
       DO nd=1,nbdirs
          mufd(nd) = muld(nd) + mutd(nd)
          kappad(nd) = kappald(nd) + kappatd(nd)
       END DO
       muf = mul + mut
       kappa = kappal + kappat
       f_n_l = 1.d0
       f_n_r = 1.d0
       IF (ul(6) .LT. 0) THEN
          chi = ul(6)/mull*ul(1)
          DO nd=1,nbdirs
             chid(nd) = (uld(nd, 6)*mull-ul(6)*mulld(nd))*ul(1)/mull**2 + ul(&
                  &         6)*uld(nd, 1)/mull
             f_n_ld(nd) = (3*chi**2*chid(nd)*(c_n1-chi**3)+(c_n1+chi**3)*3*&
                  &         chi**2*chid(nd))/(c_n1-chi**3)**2
          END DO
          f_n_l = (c_n1+chi**3)/(c_n1-chi**3)
       ELSE
          DO nd=1,nbdirs
             f_n_ld(nd) = 0.0_8
          END DO
       END IF
       IF (ur(6) .LT. 0) THEN
          chi = ur(6)/mulr*ur(1)
          DO nd=1,nbdirs
             chid(nd) = (urd(nd, 6)*mulr-ur(6)*mulrd(nd))*ur(1)/mulr**2 + ur(&
                  &         6)*urd(nd, 1)/mulr
             f_n_rd(nd) = (3*chi**2*chid(nd)*(c_n1-chi**3)+(c_n1+chi**3)*3*&
                  &         chi**2*chid(nd))/(c_n1-chi**3)**2
          END DO
          f_n_r = (c_n1+chi**3)/(c_n1-chi**3)
       ELSE
          DO nd=1,nbdirs
             f_n_rd(nd) = 0.0_8
          END DO
       END IF
       DO nd=1,nbdirs
          nufd(nd) = (muld(nd)*rhof-mul*rhofd(nd))/rhof**2 + 0.5d0*(uld(nd, &
               &       6)*f_n_l+ul(6)*f_n_ld(nd)+urd(nd, 6)*f_n_r+ur(6)*f_n_rd(nd))
       END DO
       nuf = mul/rhof + 0.5d0*(ul(6)*f_n_l+ur(6)*f_n_r)
    ELSE
       DO nd=1,nbdirs
          mufd(nd) = muld(nd)
          kappad(nd) = kappald(nd)
       END DO
       muf = mul
       kappa = kappal
       DO nd=1,nbdirs
          nufd(nd) = 0.0_8
       END DO
    END IF
    ! --- gradients at face
    ! velocity
    gduf = 0.5*(gdl(1:3, 2:4)+gdr(1:3, 2:4))
    ! temperature
    gdtf = 0.5*(gdl(1:3, 7)+gdr(1:3, 7))
    IF (turbulent) gdkf = 0.5*(gdl(1:3, 6)+gdr(1:3, 6))
    ! sa variable
    ! directional derivative
    arg1(:) = (xr-xl)**2
    arg2 = SUM(arg1(:))
    dl = SQRT(arg2)
    edgetangent = (xr-xl)/dl
    DO nd=1,nbdirs
       ! vel
       dudnd(nd, :) = (urd(nd, 2:4)-uld(nd, 2:4))/dl
       ! temperature
       dtdnd(nd) = (trd(nd)-tld(nd))/dl
    END DO
    dudn = (ur(2:4)-ul(2:4))/dl
    dtdn = (tr-tl)/dl
    IF (turbulent) THEN
       DO nd=1,nbdirs
          ! sa
          dnusadnd(nd) = (urd(nd, 6)-uld(nd, 6))/dl
       END DO
       dnusadn = (ur(6)-ul(6))/dl
    ELSE
       DO nd=1,nbdirs
          dnusadnd(nd) = 0.0_8
       END DO
    END IF
    arg1(:) = norm*edgetangent
    ndne = norm/SUM(arg1(:))
    ! velocity u
    arg1(:) = gduf(1:3, 1)*edgetangent
    DO nd=1,nbdirs
       gdufd(nd, :, :) = 0.0_8
       gdufd(nd, 1:3, 1) = ndne*dudnd(nd, 1)
       ! velocity v
       arg1d(nd, :) = edgetangent*gdufd(nd, 1:3, 2)
       gdufd(nd, 1:3, 2) = gdufd(nd, 1:3, 2) - ndne*(SUM(arg1d(nd, :))-&
            &     dudnd(nd, 2))
       ! velocity w
       arg1d(nd, :) = edgetangent*gdufd(nd, 1:3, 3)
       gdufd(nd, 1:3, 3) = gdufd(nd, 1:3, 3) - ndne*(SUM(arg1d(nd, :))-&
            &     dudnd(nd, 3))
       gdtfd(nd, 1:3) = ndne*dtdnd(nd)
    END DO
    gduf(1:3, 1) = gduf(1:3, 1) - ndne*(SUM(arg1(:))-dudn(1))
    arg1(:) = gduf(1:3, 2)*edgetangent
    gduf(1:3, 2) = gduf(1:3, 2) - ndne*(SUM(arg1(:))-dudn(2))
    arg1(:) = gduf(1:3, 3)*edgetangent
    gduf(1:3, 3) = gduf(1:3, 3) - ndne*(SUM(arg1(:))-dudn(3))
    ! temperature
    arg1(:) = gdtf(1:3)*edgetangent
    gdtf(1:3) = gdtf(1:3) - ndne*(SUM(arg1(:))-dtdn)
    IF (turbulent) THEN
       ! SA variable
       arg1(:) = gdkf*edgetangent
       DO nd=1,nbdirs
          gdkfd(nd, 1:3) = ndne*dnusadnd(nd)
       END DO
       gdkf(1:3) = gdkf(1:3) - ndne*(SUM(arg1(:))-dnusadn)
    ELSE
       DO nd=1,nbdirs
          gdkfd(nd, :) = 0.0_8
       END DO
    END IF
    tau(1, 1) = 2.d0/3.d0*muf*(2.0*gduf(1, 1)-gduf(2, 2)-gduf(3, 3))
    tau(2, 2) = 2.d0/3.d0*muf*(2.0*gduf(2, 2)-gduf(1, 1)-gduf(3, 3))
    tau(3, 3) = 2.d0/3.d0*muf*(2.0*gduf(3, 3)-gduf(1, 1)-gduf(2, 2))
    tau(1, 2) = muf*(gduf(1, 2)+gduf(2, 1))
    tau(1, 3) = muf*(gduf(1, 3)+gduf(3, 1))
    tau(2, 3) = muf*(gduf(2, 3)+gduf(3, 2))
    tau(2, 1) = tau(1, 2)
    tau(3, 1) = tau(1, 3)
    tau(3, 2) = tau(2, 3)
    DO nd=1,nbdirs
       ! -- viscous stresses tau
       taud(nd, :, :) = 0.0_8
       taud(nd, 1, 1) = 2.d0*(mufd(nd)*(2.0*gduf(1, 1)-gduf(2, 2)-gduf(3, 3&
            &     ))+muf*(2.0*gdufd(nd, 1, 1)-gdufd(nd, 2, 2)-gdufd(nd, 3, 3)))/3.d0
       taud(nd, 2, 2) = 2.d0*(mufd(nd)*(2.0*gduf(2, 2)-gduf(1, 1)-gduf(3, 3&
            &     ))+muf*(2.0*gdufd(nd, 2, 2)-gdufd(nd, 1, 1)-gdufd(nd, 3, 3)))/3.d0
       taud(nd, 3, 3) = 2.d0*(mufd(nd)*(2.0*gduf(3, 3)-gduf(1, 1)-gduf(2, 2&
            &     ))+muf*(2.0*gdufd(nd, 3, 3)-gdufd(nd, 1, 1)-gdufd(nd, 2, 2)))/3.d0
       taud(nd, 1, 2) = mufd(nd)*(gduf(1, 2)+gduf(2, 1)) + muf*(gdufd(nd, 1&
            &     , 2)+gdufd(nd, 2, 1))
       taud(nd, 1, 3) = mufd(nd)*(gduf(1, 3)+gduf(3, 1)) + muf*(gdufd(nd, 1&
            &     , 3)+gdufd(nd, 3, 1))
       taud(nd, 2, 3) = mufd(nd)*(gduf(2, 3)+gduf(3, 2)) + muf*(gdufd(nd, 2&
            &     , 3)+gdufd(nd, 3, 2))
       taud(nd, 2, 1) = taud(nd, 1, 2)
       taud(nd, 3, 1) = taud(nd, 1, 3)
       taud(nd, 3, 2) = taud(nd, 2, 3)
       ! --- phi components for energy equation
       arg1d(nd, :) = vfd(nd, 1:3)*tau(1, 1:3) + vf(1:3)*taud(nd, 1, 1:3)
       thetad(nd, :) = 0.0_8
       thetad(nd, 1) = SUM(arg1d(nd, :)) + kappad(nd)*gdtf(1) + kappa*gdtfd&
            &     (nd, 1)
       arg1d(nd, :) = vfd(nd, 1:3)*tau(2, 1:3) + vf(1:3)*taud(nd, 2, 1:3)
       thetad(nd, 2) = SUM(arg1d(nd, :)) + kappad(nd)*gdtf(2) + kappa*gdtfd&
            &     (nd, 2)
       arg1d(nd, :) = vfd(nd, 1:3)*tau(3, 1:3) + vf(1:3)*taud(nd, 3, 1:3)
       thetad(nd, 3) = SUM(arg1d(nd, :)) + kappad(nd)*gdtf(3) + kappa*gdtfd&
            &     (nd, 3)
       arg1d(nd, :) = norm*taud(nd, 1, 1:3)
       fluxd(nd, :) = 0.0_8
       fluxd(nd, 2) = SUM(arg1d(nd, :))
       arg1d(nd, :) = norm*taud(nd, 2, 1:3)
       fluxd(nd, 3) = SUM(arg1d(nd, :))
       arg1d(nd, :) = norm*taud(nd, 3, 1:3)
       fluxd(nd, 4) = SUM(arg1d(nd, :))
       arg1d(nd, :) = norm*thetad(nd, :)
       fluxd(nd, 5) = SUM(arg1d(nd, :))
    END DO
    arg1(:) = vf(1:3)*tau(1, 1:3)
    theta(1) = SUM(arg1(:)) + kappa*gdtf(1)
    arg1(:) = vf(1:3)*tau(2, 1:3)
    theta(2) = SUM(arg1(:)) + kappa*gdtf(2)
    arg1(:) = vf(1:3)*tau(3, 1:3)
    theta(3) = SUM(arg1(:)) + kappa*gdtf(3)
    ! --- viscous fluxes
    flux(1) = 0.0
    arg1(:) = tau(1, 1:3)*norm
    flux(2) = SUM(arg1(:))
    arg1(:) = tau(2, 1:3)*norm
    flux(3) = SUM(arg1(:))
    arg1(:) = tau(3, 1:3)*norm
    flux(4) = SUM(arg1(:))
    arg1(:) = theta*norm
    flux(5) = SUM(arg1(:))
    IF (turbulent) THEN
       arg1(:) = gdkf*norm
       DO nd=1,nbdirs
          ! SA variable dissipation
          arg1d(nd, :) = norm*gdkfd(nd, :)
          fluxd(nd, 6) = nufd(nd)*SUM(arg1(:))/sigma_sa + nuf*SUM(arg1d(nd, &
               &       :))/sigma_sa
       END DO
       flux(6) = nuf/sigma_sa*SUM(arg1(:))
    END IF
    DO nd=1,nbdirs
       fluxd(nd, :) = -(area*fluxd(nd, :))
    END DO
    flux = -(flux*area)
  END SUBROUTINE VFLUX_FT_SA_NEG_DV


  subroutine vflux_FT_sa_neg(UL,xL,gdL,UR,xR,gdR,ds,&
       turbulent,sutherland,flux,res1,res2)
    implicit none
    logical,                   intent(in)    :: turbulent,sutherland
    real(8), dimension(3),     intent(in)    :: xL,xR,ds
    real(8), dimension(6),     intent(in)    :: UL,UR
    real(8), dimension(3,7),   intent(in)    :: gdL,gdR
    real(8), dimension(6),     intent(out)   :: flux
    real(8),  intent(out)   :: res1,res2
    real(8) :: norm(3), area
    real(8) :: rhoF,muF,nuF,dtdn,dnuSAdn  
    real(8) :: theta(3),edgeTangent(3),vF(3),dudn(3),gdtF(3),gdkF(3)
    real(8) :: gduF(3,3),tau(3,3)

    real(8) :: mul,mulL,mulR,mut,mutL,mutR ! lam/turb dynamic viscosity
    real(8) :: kappal,kappat,kappa,dl,tL,tR
    real*8  :: chi,chi3,fv1,ndne(3)

    real*8  :: f_n_L,f_n_R
    res1=0
    res2=0

    flux=0
    area=sqrt(sum(ds**2))
    norm=ds/area

    tL   = UL(5)/UL(1)/gas_constant
    tR   = UR(5)/UR(1)/gas_constant
    if(sutherland) then
       mulL = 0.00000145d0*tL**1.5d0/(110.d0+tL)
       mulR = 0.00000145d0*tR**1.5d0/(110.d0+tR)
    else
       mulL=airViscosity
       mulR=airViscosity
    end if
    mul  = 0.5d0*(mulL+mulR)
    ! value on flux face
    kappal= Cp*mul/PrL 

    ! compute turbulent dynamic viscosity and thermal conductivity coeff 
    if(turbulent) then
       if(UL(6).le.0.d0) then
          mutL=0.d0
          ! if SA-neg, if SA working variable is negative, eddy viscosity is zero
       else
          chi   = UL(1)*UL(6)/mulL                   
          chi3  = chi**3
          fv1   = chi3/(chi3+c_v1**3)
          mutL  = UL(1)*UL(6)*fv1       
       end if
       if(UR(6).le.0.d0) then
          mutR=0.d0
          ! if SA-neg, if SA working variable is negative, eddy viscosity is zero
       else
          chi   = UR(1)*UR(6)/mulR          
          chi3  = chi**3
          fv1   = chi3/(chi3+c_v1**3)
          mutR  = UR(1)*UR(6)*fv1
       end if
       mut   = 0.5*(mutL+mutR)       
       kappat= Cp*mut/Prt
    endif

    ! --- Face values
    vF   = 0.5*(UL(2:4)+UR(2:4))
    rhoF = 0.5*(UL(1)+UR(1))
    if(turbulent) then
       muF  = mul + mut   
       ! muF is used for meanflow diffusion
       kappa   = kappal  + kappat
       f_n_L=1.d0
       f_n_R=1.d0       

       if(UL(6).lt.0) then
          chi=UL(6)/mulL*UL(1)
          f_n_L=(c_n1+chi**3)/(c_n1-chi**3)
       end if

       if(UR(6).lt.0) then
          chi=UR(6)/mulR*UR(1)
          f_n_R=(c_n1+chi**3)/(c_n1-chi**3)
       end if

       nuF = mul/rhoF+0.5d0*(UL(6)*f_n_L+UR(6)*f_n_R)
       ! nuF is used for SA working variable diffusion
    else
       muF  = mul
       kappa   = kappal
    end if

    ! --- gradients at face
    gduF = 0.5*(gdL(1:3,2:4) + gdR(1:3,2:4)) ! velocity
    gdtF = 0.5*(gdL(1:3,  7) + gdR(1:3,  7)) ! temperature
    if(turbulent) then
       gdkF = 0.5*(gdL(1:3,6) + gdR(1:3,6))  ! sa variable
    end if

    ! --- directional derivative
    dl=sqrt(sum((xR-xL)**2))
    edgeTangent=(xR-xL)/dl
    dudn=(UR(2:4)-UL(2:4))/dl     ! vel
    dtdn=(tR-tL)/dl               ! temperature
    if(turbulent) then    
       dnuSAdn=(UR(6)-UL(6))/dl   ! sa
    endif

    ! averaged gradient
    ndne=norm

    ! face tangent corrected gradient
    !ndne=norm/sum(norm*edgeTangent)

    gduF(1:3,1) = gduF(1:3,1)-ndne*(sum(gduF(1:3,1)*edgeTangent)-dudn(1)) ! velocity u
    gduF(1:3,2) = gduF(1:3,2)-ndne*(sum(gduF(1:3,2)*edgeTangent)-dudn(2)) ! velocity v
    gduF(1:3,3) = gduF(1:3,3)-ndne*(sum(gduF(1:3,3)*edgeTangent)-dudn(3)) ! velocity w
    gdtF(1:3  ) = gdtF(1:3  )-ndne*(sum(gdtF(1:3)*edgeTangent)-dtdn) ! temperature
    if(turbulent) then    
       gdkF(1:3)=gdkF(1:3)-ndne*(sum(gdkF*edgeTangent)-dnuSAdn) ! SA variable
    endif

    ! -- viscous stresses tau
    tau(1,1) = 2.d0/3.d0*muF*(2.0*gduF(1,1)-gduF(2,2)-gduF(3,3))
    tau(2,2) = 2.d0/3.d0*muF*(2.0*gduF(2,2)-gduF(1,1)-gduF(3,3))
    tau(3,3) = 2.d0/3.d0*muF*(2.0*gduF(3,3)-gduF(1,1)-gduF(2,2))
    tau(1,2) = muF*(gduF(1,2)+gduF(2,1))
    tau(1,3) = muF*(gduF(1,3)+gduF(3,1))
    tau(2,3) = muF*(gduF(2,3)+gduF(3,2))
    tau(2,1) = tau(1,2)
    tau(3,1) = tau(1,3)
    tau(3,2) = tau(2,3)

    ! --- phi components for energy equation
    theta(1)    = sum(vF(1:3)*tau(1,1:3))+kappa*gdtF(1)
    theta(2)    = sum(vF(1:3)*tau(2,1:3))+kappa*gdtF(2)
    theta(3)    = sum(vF(1:3)*tau(3,1:3))+kappa*gdtF(3)

    ! --- viscous fluxes
    flux(1) = 0.0
    flux(2) = sum(tau(1,1:3)*norm)
    flux(3) = sum(tau(2,1:3)*norm)
    flux(4) = sum(tau(3,1:3)*norm)
    flux(5) = sum(theta*norm)

    ! --- diffusion of SA variable
    if(turbulent) then
       ! flux(6) = (1+c_b2)* nuF/sigma_sa * sum(gdkF*norm)  ! SA variable dissipation
       flux(6) = nuF/sigma_sa * sum(gdkF*norm)  ! SA variable dissipation
       ! res1= c_b2/sigma_sa*(UL(6)+mulL/UL(1))*sum(gdkF*norm)*area
       ! res2=-c_b2/sigma_sa*(UR(6)+mulR/UR(1))*sum(gdkF*norm)*area
    else
       flux(6)=0
       ! res1=0
       ! res2=0
    end if
    flux = -flux * area



  end subroutine vflux_FT_sa_neg


  subroutine SA_neg( volume, dv, gradu, gradk, dist, sutherland,sourceSA)
    implicit none
    !---------------------------------------------------------------
    logical,                 intent(in)  :: sutherland
    real(8),                 intent(in)  :: volume,dist 
    real(8), dimension(6),   intent(in)  :: dv 
    real(8), dimension(3),   intent(in)  :: gradk  
    real(8), dimension(3,3), intent(in)  :: gradu  
    real(8),                 intent(out) :: sourceSA 
    ! Local variables

    real(8) :: SA_omega,S_tilde,chi,chi3,r,g 
    real(8) :: f_w,f_v1,f_v2,f_t2
    real(8) :: G_v,Y_v,Diffn     
    real(8) :: T, muL, nuT,nuL
    integer :: i,j     
    real(8) :: div
    ! Sutherland formula
    T = dv(5)/dv(1)/gas_constant
    if(sutherland) then
       muL = 0.00000145d0*T**1.5d0/(110.d0+T)
    else
       muL=airViscosity
    end if
    nuT = dv(6)

    SA_Omega = 0.d0
    Diffn = 0.d0
    do i=1,3
       do j=1,3
          SA_Omega  = SA_Omega+(gradu(j,i)-gradu(i,j))**2
       end do
       Diffn = Diffn+gradk(i)**2
    end do
    SA_Omega = sqrt(0.5d0*SA_Omega)

    ! calculate production term coefficients
    nuL = muL/dv(1)
    chi = dv(6)/nuL 
   
    chi3=chi**3
    f_v1 = chi3/(chi3+c_v1**3)
    f_v2 = 1-chi/(1+chi*f_v1)     

    if(nut.lt.0) then
       G_v=c_b1*(1-c_t3)*SA_Omega*nuT
       Y_v=-c_w1*(nuT/dist)**2
    else
       ! SA-standard Note 1(c)
       S_tilde=dv(6)*f_v2/kappa_sa**2/dist**2
       if(S_tilde .ge. -c_2*SA_Omega) then
          S_tilde=S_tilde+SA_omega   
       else
          S_tilde=SA_Omega*(c_2**2*SA_Omega+c_3*S_tilde)/((c_3-2*c_2)*SA_Omega-S_tilde)
          S_tilde=S_tilde+SA_Omega
       end if
       ! calculate destruction term coefficients
       if(S_tilde.ne.0) then
          r = nuT/(S_tilde*kappa_sa**2*dist**2)
          r=min(r,10.d0)
       else
          r=10.d0
       end if
       g = r+c_w2*(r**6-r)
       f_w = g*((1+c_w3**6)/(g**6+c_w3**6))**(1.d0/6.d0)
       f_t2=c_t3*exp(-c_t4*chi**2)

       !calculate SA source term
       G_v      = c_b1 *(1-f_t2)* S_tilde * nuT
       Y_v      = (c_w1*f_w-c_b1*f_t2/kappa_sa**2)*(nuT/dist)**2
    end if

   Diffn    = c_b2/sigma_sa*Diffn    
! !    Diffn    = (1+c_b2)/sigma_sa*Diffn    
!      Diffn    = 0

     div=0

     ! div=dv(6)*(gradu(1,1)+gradu(2,2)+gradu(3,3))

    sourceSA = -(G_v-Y_v+Diffn+div)*volume

  end subroutine SA_neg


  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.12 (r6213) -  9 May 2017 14:18
  !
  !  Differentiation of sa_neg in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
  !   variations   of useful results: sourcesa
  !   with respect to varying inputs: dv
  !   RW status of diff variables: dv:in sourcesa:out
  SUBROUTINE SA_NEG_DV(volume, dv, dvd, gradu, gradk, dist, sutherland, &
       & sourcesa, sourcesad, nbdirs)
    !  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    !---------------------------------------------------------------
    integer,parameter::nbdirsmax=6
    LOGICAL, INTENT(IN) :: sutherland
    REAL*8, INTENT(IN) :: volume, dist
    REAL*8, DIMENSION(6), INTENT(IN) :: dv
    REAL*8, DIMENSION(nbdirsmax, 6), INTENT(IN) :: dvd
    REAL*8, DIMENSION(3), INTENT(IN) :: gradk
    REAL*8, DIMENSION(3, 3), INTENT(IN) :: gradu
    REAL*8, INTENT(OUT) :: sourcesa
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sourcesad
    ! Local variables
    REAL*8 :: sa_omega, s_tilde, chi, chi3, r, g
    REAL*8, DIMENSION(nbdirsmax) :: s_tilded, chid, chi3d, rd, gd
    REAL*8 :: f_w, f_v1, f_v2, f_t2
    REAL*8, DIMENSION(nbdirsmax) :: f_wd, f_v1d, f_v2d, f_t2d
    REAL*8 :: g_v, y_v, diffn
    REAL*8, DIMENSION(nbdirsmax) :: g_vd, y_vd
    REAL*8 :: t, mul, nut, nul
    REAL*8, DIMENSION(nbdirsmax) :: td, muld, nutd, nuld
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC MIN
    INTRINSIC EXP
    REAL*8 :: pwx1
    REAL*8, DIMENSION(nbdirsmax) :: pwx1d
    REAL*8 :: pwr1
    REAL*8, DIMENSION(nbdirsmax) :: pwr1d
    REAL*8 :: arg1
    REAL*8, DIMENSION(nbdirsmax) :: arg1d
    INTEGER :: nd
    INTEGER :: nbdirs
    DO nd=1,nbdirs
       ! Sutherland formula
       td(nd) = (dvd(nd, 5)*dv(1)-dv(5)*dvd(nd, 1))/dv(1)**2/gas_constant
    END DO
    t = dv(5)/dv(1)/gas_constant
    IF (sutherland) THEN
       DO nd=1,nbdirs
          muld(nd) = (0.00000145d0*1.5d0*t**0.5D0*td(nd)*(110.d0+t)-&
               &       0.00000145d0*t**1.5d0*td(nd))/(110.d0+t)**2
       END DO
       mul = 0.00000145d0*t**1.5d0/(110.d0+t)
    ELSE
       mul = airviscosity
       DO nd=1,nbdirs
          muld(nd) = 0.0_8
       END DO
    END IF
    DO nd=1,nbdirs
       nutd(nd) = dvd(nd, 6)
    END DO
    nut = dv(6)
    sa_omega = 0.d0
    diffn = 0.d0
    DO i=1,3
       DO j=1,3
          sa_omega = sa_omega + (gradu(j, i)-gradu(i, j))**2
       END DO
       diffn = diffn + gradk(i)**2
    END DO
    sa_omega = SQRT(0.5d0*sa_omega)
    nul = mul/dv(1)
    chi = dv(6)/nul
    chi3 = chi**3
    f_v1 = chi3/(chi3+c_v1**3)
    DO nd=1,nbdirs
       ! calculate production term coefficients
       nuld(nd) = (muld(nd)*dv(1)-mul*dvd(nd, 1))/dv(1)**2
       chid(nd) = (dvd(nd, 6)*nul-dv(6)*nuld(nd))/nul**2
       chi3d(nd) = 3*chi**2*chid(nd)
       f_v1d(nd) = (chi3d(nd)*(chi3+c_v1**3)-chi3*chi3d(nd))/(chi3+c_v1**3)&
            &     **2
       f_v2d(nd) = -((chid(nd)*(1+chi*f_v1)-chi*(chid(nd)*f_v1+chi*f_v1d(nd&
            &     )))/(1+chi*f_v1)**2)
    END DO
    f_v2 = 1 - chi/(1+chi*f_v1)
    IF (nut .LT. 0) THEN
       DO nd=1,nbdirs
          g_vd(nd) = c_b1*(1-c_t3)*sa_omega*nutd(nd)
          y_vd(nd) = -(c_w1*2*nut*nutd(nd)/dist**2)
       END DO
       g_v = c_b1*(1-c_t3)*sa_omega*nut
       y_v = -(c_w1*(nut/dist)**2)
    ELSE
       DO nd=1,nbdirs
          ! SA-standard Note 1(c)
          s_tilded(nd) = (dvd(nd, 6)*f_v2+dv(6)*f_v2d(nd))/kappa_sa**2/dist&
               &       **2
       END DO
       s_tilde = dv(6)*f_v2/kappa_sa**2/dist**2
       IF (s_tilde .GE. -(c_2*sa_omega)) THEN
          s_tilde = s_tilde + sa_omega
       ELSE
          DO nd=1,nbdirs
             s_tilded(nd) = (sa_omega*c_3*s_tilded(nd)*((c_3-2*c_2)*sa_omega-&
                  &         s_tilde)+sa_omega*(c_2**2*sa_omega+c_3*s_tilde)*s_tilded(nd))/&
                  &         ((c_3-2*c_2)*sa_omega-s_tilde)**2
          END DO
          s_tilde = sa_omega*(c_2**2*sa_omega+c_3*s_tilde)/((c_3-2*c_2)*&
               &       sa_omega-s_tilde)
          s_tilde = s_tilde + sa_omega
       END IF
       ! calculate destruction term coefficients
       IF (s_tilde .NE. 0) THEN
          DO nd=1,nbdirs
             rd(nd) = (nutd(nd)*s_tilde*kappa_sa**2*dist**2-nut*kappa_sa**2*&
                  &         dist**2*s_tilded(nd))/(s_tilde*kappa_sa**2*dist**2)**2
          END DO
          r = nut/(s_tilde*kappa_sa**2*dist**2)
          IF (r .GT. 10.d0) THEN
             r = 10.d0
             DO nd=1,nbdirs
                rd(nd) = 0.0_8
             END DO
          ELSE
             r = r
          END IF
       ELSE
          r = 10.d0
          DO nd=1,nbdirs
             rd(nd) = 0.0_8
          END DO
       END IF
       g = r + c_w2*(r**6-r)
       pwx1 = (1+c_w3**6)/(g**6+c_w3**6)
       pwr1 = pwx1**(1.d0/6.d0)
       f_w = g*pwr1
       arg1 = -(c_t4*chi**2)
       f_t2 = c_t3*EXP(arg1)
       DO nd=1,nbdirs
          gd(nd) = rd(nd) + c_w2*(6*r**5*rd(nd)-rd(nd))
          pwx1d(nd) = -((1+c_w3**6)*6*g**5*gd(nd)/(g**6+c_w3**6)**2)
          IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. 1.d0/6.d0 .EQ. INT(&
               &         1.d0/6.d0))) THEN
             pwr1d(nd) = pwx1**(1.d0/6.d0-1)*pwx1d(nd)/6.d0
          ELSE IF (pwx1 .EQ. 0.0 .AND. 1.d0/6.d0 .EQ. 1.0) THEN
             pwr1d(nd) = pwx1d(nd)
          ELSE
             pwr1d(nd) = 0.0
          END IF
          f_wd(nd) = gd(nd)*pwr1 + g*pwr1d(nd)
          arg1d(nd) = -(c_t4*2*chi*chid(nd))
          f_t2d(nd) = c_t3*arg1d(nd)*EXP(arg1)
          !calculate SA source term
          g_vd(nd) = c_b1*((1-f_t2)*(s_tilded(nd)*nut+s_tilde*nutd(nd))-&
               &       f_t2d(nd)*s_tilde*nut)
          y_vd(nd) = (c_w1*f_wd(nd)-c_b1*f_t2d(nd)/kappa_sa**2)*nut**2/dist&
               &       **2 + (c_w1*f_w-c_b1*f_t2/kappa_sa**2)*2*nut*nutd(nd)/dist**2
       END DO
       g_v = c_b1*(1-f_t2)*s_tilde*nut
       y_v = (c_w1*f_w-c_b1*f_t2/kappa_sa**2)*(nut/dist)**2
    END IF
    diffn = c_b2/sigma_sa*diffn
    DO nd=1,nbdirs
       sourcesad(nd) = -(volume*(g_vd(nd)-y_vd(nd)))
    END DO
    sourcesa = -((g_v-y_v+diffn)*volume)
  END SUBROUTINE SA_NEG_DV




  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.12 (r6213) -  9 May 2017 14:18
  !
  !  Differentiation of turbulentdynamicviscosity in forward (tangent) mode (with options multiDirectional):
  !   variations   of useful results: mut
  !   with respect to varying inputs: mul rho nutsa
  SUBROUTINE TURBULENTDYNAMICVISCOSITY_DV(rho, rhod, mul, muld, nutsa, &
       & nutsad, mut, mutd, nbdirs)

    !  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    integer,parameter::nbdirsmax=6
    REAL*8, INTENT(IN) :: rho, mul, nutsa
    REAL*8, DIMENSION(nbdirsmax), INTENT(IN) :: rhod, muld, nutsad
    REAL*8, INTENT(OUT) :: mut
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: mutd
    REAL*8 :: chi, chi3, fv1
    REAL*8, DIMENSION(nbdirsmax) :: chid, chi3d, fv1d
    INTEGER :: nd
    INTEGER :: nbdirs

    chi = rho*nutsa/mul
    chi3 = chi**3
    fv1 = chi3/(chi3+c_v1**3)
    DO nd=1,nbdirs
       chid(nd) = ((rhod(nd)*nutsa+rho*nutsad(nd))*mul-rho*nutsa*muld(nd))/&
            &     mul**2
       chi3d(nd) = 3*chi**2*chid(nd)
       fv1d(nd) = (chi3d(nd)*(chi3+c_v1**3)-chi3*chi3d(nd))/(chi3+c_v1**3)&
            &     **2
       mutd(nd) = (rhod(nd)*nutsa+rho*nutsad(nd))*fv1 + rho*nutsa*fv1d(nd)
    END DO
    mut = rho*nutsa*fv1
  END SUBROUTINE TURBULENTDYNAMICVISCOSITY_DV



  subroutine solutionSanityCheck(grid)
    implicit none
    type(grid_t)          :: grid
    integer::i,j,ierr

    do i = 1, grid%ninternal_nodes
       do j = 1, 6
          if(grid%dv(j,i).ne.grid%dv(j,i)) then
             print*,'after update NAN detected for variable ',j,' node ',i,' rank=',rank
             write(*,'(A,3F14.6)')'location=',grid%x(1:3,i)
             stop
             !             call MPI_Abort( MPI_COMM_WORLD, -1, ierr )
          end if
       end do
    end do

    do i = 1, grid%ninternal_nodes
       do j = 1, 6
          if(grid%dv(1,i).lt.0.0) then
             print*,'negative density after update, node',i,' rank=',rank
             print*,'density=',grid%dv(1,i)
             write(*,'(A,3F14.6)')'location=',grid%x(1:3,i)

             print*,'walldist',grid%dist(i)
             print*,'velocity',grid%dv(2:4,i)
             stop
             !call MPI_Abort( MPI_COMM_WORLD, -1, ierr)
          end if
          if(grid%dv(5,i).lt.0.0) then
             print*,'negative pressure after update, node',i,' rank=',rank
             print*,'pressure=',grid%dv(5,i)
             write(*,'(A,3F14.6)')'location=',grid%x(1:3,i)

             print*,'walldist',grid%dist(i)
             print*,'velocity',grid%dv(2:4,i)
             stop
             !call MPI_Abort( MPI_COMM_WORLD, -1, ierr)
          end if
       end do
    end do

  end subroutine solutionSanityCheck


  subroutine residualNanCheck(grid)
    implicit none
    type(grid_t)          :: grid
    integer::i,j,ierr

    do i = 1, grid%ninternal_nodes
       do j = 1, 6
          if(grid%r(j,i).ne.grid%r(j,i)) then
             print*,'res NAN detected for variable ',j,' node ',i
             write(*,'(A,3F14.6)')'location=',grid%x(1:3,i)
             stop
             !             call MPI_Abort( MPI_COMM_WORLD, -1, ierr )
          end if
       end do
    end do
  end subroutine residualNanCheck

end module flux_module
