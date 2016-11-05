      subroutine viscous_cmt(e,eq)
      include  'SIZE'
      include  'CMTDATA'
      include  'DG'
      include  'INPUT'
!-----------------------------------------------------------------------
! diagnostic
!-----------------------------------------------------------------------
      include 'GEOM' ! diagnostic
      include 'SOLN' ! diagnostic
!-----------------------------------------------------------------------
! diagnostic
!-----------------------------------------------------------------------

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelcmt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*3*lfq) ! might not need ldim
      common /CMTSURFLX/ fatface(heresize),graduf(hdsize)
      real fatface,graduf

      integer e,eq

      if (eq .lt. toteq) then ! not energy
         if (eq .gt. ndim+1) return ! not if3d
      endif

      nxyz=nx1*ny1*nz1
      nfq=nx1*nz1*2*ndim*nelt
      nstate = nqq
! where different things live
      iqm =1
      iqp =iqm+nstate*nfq
      iuj =iqp+nstate*nfq

! apply viscous flux jacobian A.
      call fluxj_ns(diffh,gradu,e,eq)
! monolithic regularization never
!     call fluxj(diffh,gradu,e,eq)

!-----------------------------------------------------------------------
! diagnostic
!-----------------------------------------------------------------------
      if (eq.eq.2) then
         if (e.eq.1) then
            write(100,*) '#xmom:x,y,tau11,Adu,diff'
            write(200,*) '#xmom:x,y,tau12,Adu,diff'
            write(101,*) '#x,y,du1dx,gradu,diff,du1dy,gradu,diff'
            write(102,*) '#x,y,du2dx,gradu,diff,du2dy,gradu,diff'
            write(103,*) '#x,y,du3dx,gradu,diff,du3dy,gradu,diff'
            write(105,*) '#x,y,du5dx,gradu,diff,du5dy,gradu,diff'
         endif
         pi=4.0*atan(1.0)
         do i=1,nxyz
            beta=5.0
            x=xm1(i,1,1,e)-5.0
            y=ym1(i,1,1,e)
            onemr2=1.0-(x**2+y**2)
            tau11=vdiff(i,1,1,e,imu)*(2.0*beta*x*y/pi)*exp(onemr2)
            tau12=vdiff(i,1,1,e,imu)*(beta/pi)*(y**2-x**2)*exp(onemr2)
            write(100,'(5e17.8)') x,y,tau11,diffh(i,1),diffh(i,1)-tau11
            du1=du1dx(beta,x,y,1.0,1.4)
            du2=du1dy(beta,x,y,1.0,1.4)
            write(101,'(8e17.8)') x,y,gradu(i,1,1),du1,gradu(i,1,1)-du1,
     >                                gradu(i,1,2),du2,gradu(i,1,2)-du2
            du1=du2dx(beta,x,y,1.0,1.4)
            du2=du2dy(beta,x,y,1.0,1.4)
            write(102,'(8e17.8)') x,y,gradu(i,2,1),du1,gradu(i,2,1)-du1,
     >                                gradu(i,2,2),du2,gradu(i,2,2)-du2
            du1=du3dx(beta,x,y,1.0,1.4)
            du2=du3dy(beta,x,y,1.0,1.4)
            write(103,'(8e17.8)') x,y,gradu(i,3,1),du1,gradu(i,3,1)-du1,
     >                                gradu(i,3,2),du2,gradu(i,3,2)-du2
            du1=du5dx(beta,x,y,1.0,1.4)
            du2=du5dy(beta,x,y,1.0,1.4)
            write(105,'(8e17.8)') x,y,gradu(i,5,1),du1,gradu(i,5,1)-du1,
     >                                gradu(i,5,2),du2,gradu(i,5,2)-du2
            write(200,'(5e17.8)') x,y,tau12,diffh(i,2),diffh(i,2)-tau12
         enddo
      elseif(eq.eq.3) then
         if (e.eq.1) then
            write(300,*) '#ymom:x,y,tau21,Adu,diff'
            write(400,*) '#ymom:x,y,tau22,Adu,diff'
         endif
         pi=4.0*atan(1.0)
         do i=1,nxyz
            beta=5.0
            x=xm1(i,1,1,e)-5.0
            y=ym1(i,1,1,e)
            onemr2=1.0-(x**2+y**2)
            tau22=-vdiff(i,1,1,e,imu)*(2.0*beta*x*y/pi)*exp(onemr2)
            tau12=vdiff(i,1,1,e,imu)*(beta/pi)*(y**2-x**2)*exp(onemr2)
            write(300,'(5e17.8)') x,y,tau12,diffh(i,1),diffh(i,1)-tau12
            write(400,'(5e17.8)') x,y,tau22,diffh(i,2),diffh(i,2)-tau22
         enddo
      elseif(eq.eq.5) then
         if (e.eq.1) then
            write(500,*) '#energy:x,y,h1,diff'
            write(600,*) '#energy:x,y,h2,diff'
         endif
         pi=4.0*atan(1.0)
         do i=1,nxyz
            beta=5.0
            x=xm1(i,1,1,e)-5.0
            y=ym1(i,1,1,e)
            onemr2=1.0-(x**2+y**2)
            uxinf=0.0
            uyinf=0.0 ! betta match useric
            zeu=uxinf-beta*exp(onemr2)*y/(2.0*pi)
            zev=uyinf+beta*exp(onemr2)*x/(2.0*pi)
            tau11=vdiff(i,1,1,e,imu)*(2.0*beta*x*y/pi)*exp(onemr2)
            tau12=vdiff(i,1,1,e,imu)*(beta/pi)*(y**2-x**2)*exp(onemr2)
            tau22=-vdiff(i,1,1,e,imu)*(2.0*beta*x*y/pi)*exp(onemr2)
            du1=(tau11*zeu+zev*tau12)+vdiff(i,1,1,e,iknd)*
     >          beta**2*x*0.25/1.4*0.4/pi/pi*exp(2.0*onemr2)
            du2=(tau12*zeu+zev*tau22)+vdiff(i,1,1,e,iknd)*
     >          beta**2*y*0.25/1.4*0.4/pi/pi*exp(2.0*onemr2)
            write(500,'(5e17.8)') x,y,diffh(i,1),du1,diffh(i,1)-du1
            write(600,'(5e17.8)') x,y,diffh(i,2),du2,diffh(i,2)-du2
         enddo
      endif
! diagnostic
!-----------------------------------------------------------------------

      call diffh2graduf(e,eq,graduf) ! on faces for QQ^T and igu_cmt

! volume integral involving "DG-ish" stiffness operator K
      call half_iku_cmt(res1(1,1,1,e,eq),diffh,e)

      return
      end

!-----------------------------------------------------------------------

      subroutine igtu_cmt(qminus,ummcu,hface)

!     Vol integral [[u]].{{gradv}}. adapted from Lu's dgf3.f;
! interior penalty stuff could go
! here too. Debug subroutine ihu in heat.usr and try to fold it in here.

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'DG'      ! iface
      include 'CMTDATA'

! arguments
      real qminus(nx1*nz1,2*ndim,nelt,*)    ! intent(in)
      real ummcu(nx1*nz1*2*ndim,nelt,toteq) ! intent(in)
      real hface(nx1*nz1*2*ndim*nelt,toteq,3) ! intent(out) scratch

! commons and scratch
      common /scrns/ superhugeh(lx1*ly1*lz1*lelt,3) ! like totalh, but super-huge
      common /scruz/ gradm1_t_overwrites(lx1*ly1*lz1*lelt) ! sigh
      real superhugeh,gradm1_t_overwrites
      common /ctmp0/ viscscr(lx1,ly1,lz1)
      real viscscr
!     parameter (lfq=lx1*lz1*2*ldim*lelt)
!     common /ctmp0/ ftmp1(lfq),ftmp2(lfq)
!     real ftmp1,ftmp2

      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/

      integer e, eq, n, npl, nf, f, i, k, eq2

      nvol = nx1*ny1*nz1*nelt

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nf     = nxz*nfaces ! 1 element's face points
      nfq    = nf*nelt ! all points in a pile of faces
      if (ifsip) then
         const=-1.0 ! SIP
      else
         const=1.0 ! Baumann-Oden
      endif

! compute (U-{{U}})_i * n_k
      l=1
      do e=1,nelt
         do eq=1,toteq
            call col3(hface(l,eq,3),ummcu(1,e,eq),area(1,1,1,e),nf)
            call col3(hface(l,eq,1),hface(l,eq,3),unx(1,1,1,e), nf)
            call col3(hface(l,eq,2),hface(l,eq,3),uny(1,1,1,e), nf)
            if(if3d) call col2(hface(l,eq,3),unz(1,1,1,e),nf)
         enddo
         l=l+nf
      enddo

      call rzero(superhugeh,3*lx1*ly1*lz1*lelt)
      nxyz=nx1*ny1*nz1
      ngradu=nxyz*toteq*3
      do eq=2,toteq ! MASS DIFFUSION GETS ITS OWN BLOCK. somewhere
         call rzero(superhugeh,3*lx1*ly1*lz1*lelt)
         if (eq .eq. 4 .and. .not. if3d) goto 133
         l=1
         m=1
         do e=1,nelt
            call rzero(gradu,ngradu) ! this too goes away when gradu is global
            do j=1,ndim
               do eq2=1,toteq ! sigh
                  call add_face2full_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                                gradu(1,eq2,j),hface(l,eq2,j))
               enddo
            enddo

            l=l+nf ! index for hface, which is global. this all goes away
                ! once you get rid of that execrable "element loop" in
                ! compute_rhs_and_dt
!           call fluxj_ns(superhugeh,... THIS will be correctly strided as well
            do j=1,ndim    ! flux direction
               do k=1,ndim ! dU   direction
                  ieijk=0
                  if (eq .lt. toteq) ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
                  if (ieijk .eq. 0) then
                     call agradu_ns(superhugeh(m,j),gradu(1,1,k),viscscr
     >                             ,e,eq,j,k) ! the worst stride ever
                  endif
               enddo
            enddo

            m=m+nxyz

         enddo ! element loop

         call gradm1_t(gradm1_t_overwrites,superhugeh(1,1),
     >                        superhugeh(1,2),superhugeh(1,3))
         call cmult(gradm1_t_overwrites,const,nvol)
         call add2(res1(1,1,1,1,eq),gradm1_t_overwrites,nvol)
133      continue
      enddo

!!! apply flux jacobian to get Ajac (U-{{U}})_i * n_k
!!! JH110116 DEPRECATED. Always apply A to volume, not surface points.
!!!                      igtu_cmt will reflect this change in philosophy.
!!!                      Less to debug that way
!!         do j=1,ndim
!!            call rzero(ftmp1,nfq)
!!!           call agradu_sfc(ftmp1,qminus,hface(1,1,j),eq)
!!! yes I know I should wrap this. Bite me.
!!            do k=1,ndim
!!               call agradu_ns_sfc(ftmp1,qminus,hface(1,1,k),
!!     >                               ftmp2,eq,j,k)
!!            enddo
!!            call add_face2full_cmt(nelt,nx1,ny1,nz1,iface_flux,
!!     >                      superhugeh(1,j),ftmp1)
!!         enddo
!

      return
      end

!-----------------------------------------------------------------------

      subroutine convective_cmt(e,eq)
! JH081916 convective flux divergence integrated in weak form and
!          placed in res1.
      include 'SIZE'
      include 'CMTDATA'
      integer e,eq

      n=3*lxd*lyd*lzd

      call rzero(convh,n)
      if (nxd.gt.nx1) then
         call evaluate_dealiased_conv_h(e,eq)
         call copy(totalh,convh,n)
         call flux_div_integral_dealiased(e,eq)
      else
         call evaluate_aliased_conv_h(e,eq)
         call copy(totalh,convh,n)
         call flux_div_integral_aliased(e,eq)
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine evaluate_dealiased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'
     
      integer  e,eq

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2

      n=nxd*nyd*nzd

      if (eq .eq. 1) then ! convective flux of mass=rho u_j=U_{j+1}

         do j=1,ndim
            call intp_rstd(convh(1,j),u(1,1,1,eq+j,e),nx1,nxd,if3d,0)
         enddo

      else

c To be consistent with momentum equation, for mass balance flux vector is 
c computed by multiplying rho by u_j
         call intp_rstd(ju1,phig(1,1,1,e),nx1,nxd,if3d,0)
         call intp_rstd(ju2,pr(1,1,1,e),nx1,nxd,if3d,0)

         if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
            call add2col2(convh(1,eq-1),ju1,ju2,n)

         elseif (eq .eq. 5) then

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            call add2col2(convh(1,1),ju1,ju2,n)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            call col2(convh(1,3),vzd(1,1,1,e),n)

         else
            if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
            if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
            call exitt
         endif

      endif
     
      return
      end

!-----------------------------------------------------------------------

      subroutine evaluate_aliased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2
      integer  e,eq

      n=nxd*nyd*nzd

      call copy(ju1,phig(1,1,1,e),n)
      call copy(ju2,pr(1,1,1,e),n)

      if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
         if(eq. gt. 1) call add2col2(convh(1,eq-1),ju1,ju2,n)

      elseif (eq .eq. 5) then

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         call add2col2(convh(1,1),ju1,ju2,n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         call col2(convh(1,3),vzd(1,1,1,e),n)

      else
         if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
         if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
         call exitt
      endif

      return
      end

!----------------------

      subroutine flux_div_integral_dealiased(e,eq)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      integer  dir
      parameter (ldd=lxd*lyd*lzd)
      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call get_dgl_ptr(ip,nxd,nxd) ! fills dg, dgt
      mdm1=nxd-1

      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ndim
         j0=j0+1
         call add2col2(ur,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      do j=1,ndim
         j0=j0+1
         call add2col2(us,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ndim
            j0=j0+1
            call add2col2(ut,totalh(1,j),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,dg(ip),dgt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,dg(ip),dgt(ip),wkd)
      endif

      call intp_rstd(tu,ud,nx1,nxd,if3d,1)

! multiply the residue by mass matrix. Non diagonal B should still be
! one block per element
!     call col2(ud,bm1(1,1,1,e),nxyz)  ! res = B*res  --- B=mass matrix
!     call add2(res1(1,1,1,e,eq),tu,nxyz)
! weak?
      call sub2(res1(1,1,1,e,eq),tu,nxyz)

      return
      end

!----------------------

      subroutine flux_div_integral_aliased(e,eq)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'CMTDATA'

      integer  e,eq
      integer  dir
      parameter (ldd=lxd*lyd*lzd)
      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      character*32 cname

      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call get_dgll_ptr(ip,nxd,nxd) ! fills dg, dgt
      mdm1=nxd-1

      call rzero(ur,nrstd)
      call rzero(us,nrstd)
      call rzero(ut,nrstd)
      call rzero(ud,nrstd)
      call rzero(tu,nrstd)

      j0=0
      do j=1,ndim
         j0=j0+1
         call add2col2(ur,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      do j=1,ndim
         j0=j0+1
         call add2col2(us,totalh(1,j),rx(1,j0,e),nrstd)
      enddo
      if (if3d) then
         do j=1,ndim
            j0=j0+1
            call add2col2(ut,totalh(1,j),rx(1,j0,e),nrstd)
         enddo
         call local_grad3_t(ud,ur,us,ut,mdm1,1,d(ip),dt(ip),wkd)
      else
         call local_grad2_t(ud,ur,us,   mdm1,1,d(ip),dt(ip),wkd)
      endif

      call copy(tu,ud,nxyz)

! needs fleg or removal altogether. not good modularity
      call sub2(res1(1,1,1,e,eq),tu,nxyz)

      return
      end

!-----------------------------------------------------------------------
      subroutine compute_forcing(e,eq_num)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'SOLN'
      include  'CMTDATA'
      include  'DEALIAS'
      
      integer e,eq_num
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      nrstd=ldd
      nxyz=nx1*ny1*nz1
      call rzero(ud,nxyz)
      if(eq_num.ne.1.and.eq_num.ne.5)then
        if(eq_num.eq.2)then
           j=1
        elseif(eq_num.eq.3)then
           j=2
        elseif(eq_num.eq.4)then
           j=2
           if(ldim.eq.3) j=3
        endif
c       write(6,*)'enter  compute_forcing ', j
        call gradl_rst(ur,us,ut,phig(1,1,1,e),lx1,if3d) ! navier1
        if (if3d) then
           j0=j+0
           j3=j+3
           j6=j+6
           do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
              ud(i)=rx(i,j0,e)*ur(i)+rx(i,j3,e)*us(i)+rx(i,j6,e)*ut(i)
           enddo
        else
           j0=j+0
           j2=j+2
           do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
              ud(i)=rx(i,j0,e)*ur(i)+rx(i,j2,e)*us(i)
           enddo
        endif
        if (eq_num.eq.4.and.ldim.eq.2)then

        else
           call col2(ud,pr(1,1,1,e),nxyz)
           call copy(convh(1,1),ud,nxyz)
           call col2(convh(1,1),jacmi(1,e),nxyz)
           call col2(convh(1,1),bm1(1,1,1,e),nxyz)  ! res = B*res
           call sub2(res1(1,1,1,e,eq_num),convh(1,1),nxyz)
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
        endif
      elseif(eq_num.eq.5)then
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cmtusrf(e)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'TSTEP'
      include 'PARALLEL'

      integer e,eg

      if(istep.eq.1)then
        n = nx1*ny1*nz1*5
        call rzero(usrf,n)
      endif
      eg = lglel(e)
      do k=1,nz1
         do j=1,ny1
            do i=1,nx1
               call NEKASGN(i,j,k,e)
               call userf(i,j,k,eg)
               usrf(i,j,k,2) = FFX
               usrf(i,j,k,3) = FFY
               usrf(i,j,k,4) = FFZ
               usrf(i,j,k,5) = (U(i,j,k,2,e)*FFX + U(i,j,k,3,e)*FFY
     &                       +  U(i,j,k,4,e)*FFZ)/ U(i,j,k,1,e)
            enddo
         enddo
      enddo

      return
      end 
