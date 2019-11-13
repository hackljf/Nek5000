! JH100919
! Now with some code introduced by Jacob's work on Tait-equation water
! mixture modeling. Commented out until shocks are working in essplit
C> @file driver3_cmt.f routines for primitive variables, usr-file interfaces
C> and properties. Also initializes flow field.

C> \ingroup state
C> @{
C> \brief Compute primitive variables (velocity, thermodynamic state) from 
C>        conserved unknowns U and store them in SOLN and CMTDATA
      subroutine compute_primitive_vars(ilim)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'GEOM'
      include 'CMTDATA'
      include 'SOLN'
      include 'DEALIAS' ! until we are comfortable with setup_convect

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ energy(lx1,ly1,lz1),scr(lx1,ly1,lz1)
      integer e, eq
      common /posflags/ ifailr,ifaile,ifailt,ilimflag
      integer ifailr,ifaile,ifailt,ilimflag
C> ilim is a flag for positivity checks. ilim==0 means ``do not
C> perform positivity checks.'' ilim !=0 means ``perform positivity
C> checks and exit with a diagnostic dump if density, energy or
C> temperature fall to zero or below at any GLL node.''
      integer ilim

      nxyz= lx1*ly1*lz1
      ntot=nxyz*nelt
C> Flags for density, energy and temperature positivity
      ifailr=-1 ! density
      ifaile=-1 ! energy
      ifailt=-1 ! temperature
      ilimflag=ilim

      do e=1,nelt
C> Density positivity check.
         dmin=vlmin(u(1,1,1,irg,e),nxyz)
         if (dmin .lt. 0.0 .and. ilim .ne. 0) then
            ifailr=lglel(e)
            write(6,*) nid,'***NEGATIVE DENSITY***',dmin,lglel(e)
         endif
C> Divide momentum by density to get velocity
         call invcol3(vx(1,1,1,e),u(1,1,1,irpu,e),u(1,1,1,irg,e),nxyz)
         call invcol3(vy(1,1,1,e),u(1,1,1,irpv,e),u(1,1,1,irg,e),nxyz)
!        if (if3d)
         call invcol3(vz(1,1,1,e),u(1,1,1,irpw,e),u(1,1,1,irg,e),nxyz)
C> Compute kinetic energy using vdot2/3
         if (if3d) then
            call vdot3(scr,
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
         else
            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
         endif
         call invcol2(scr,u(1,1,1,irg,e),nxyz)
         call cmult(scr,0.5,nxyz)
C> Compute internal energy. First, subtract volume-fraction-weighted kinetic energy from total
C> energy.
         call sub3(energy,u(1,1,1,iret,e),scr,nxyz)
C> Then, divide internal energy by density.
         call invcol2(energy,u(1,1,1,irg,e),nxyz)
C> Compute density by dividing U1 by gas volume fraction. store in vtrans(:,jrho)
         call invcol3(vtrans(1,1,1,e,jrho),u(1,1,1,irg,e),phig(1,1,1,e),
     >                nxyz)
C> Check positivity of internal energy.
         emin=vlmin(energy,nxyz)
         if (emin .lt. 0.0 .and. ilim .ne. 0) then
            ifaile=lglel(e)
            write(6,*) stage,nid, ' HAS NEGATIVE ENERGY ',emin,lglel(e)
         endif
!! JH070219 Tait mixture model mass fractions. just one for now
!c JB080119 go throug hmultiple species
!         call invcol3(t(1,1,1,e,2),u(1,1,1,imfrac,e),
!     >                u(1,1,1,irg,e),
!     >                nxyz)
!c        do iscal = 1,NPSCAL
!c        call invcol3(t(1,1,1,e,1+iscal),u(1,1,1,imfrac+iscal-1,e),
!c    >                u(1,1,1,irg,e),
!c    >                nxyz)
!c        enddo
C> Compute thermodynamic state variables cv, T and p. Check temperature
C> positivity using ifailt in /posflags/
         call tdstate(e,energy)
      enddo

! Avoid during EBDG testing
! JH070219 Tait mixture model: man up and test T(:,2) for positivity
!          someday.
C> call poscheck for each of the posflags and exit if ilim!=0 and
C> any posflag>0. Nonzero posflags are set to the global element
C> number where the first positivity failure on each MPI task was
C> encountered.
      call poscheck(ifailr,'density    ')
      call poscheck(ifaile,'energy     ')
      call poscheck(ifailt,'temperature')
C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup state
C> @{
C> \brief calls cmt_userEOS in the usr file.
C> Compute thermodynamic state for element e from internal energy and density.
      subroutine tdstate(e,energy)!,energy)
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'
      include 'PARALLEL'
      include 'NEKUSE'
      integer   e,eg
      real energy(lx1,ly1,lz1)

      common /posflags/ ifailr,ifaile,ifailt,ilimflag
      integer ifailr,ifaile,ifailt,ilimflag

      eg = lglel(e)
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
C> loop over GLL nodes in element e. Fill /NEKUSE/ and /nekuscmt/
C> by nekasgn and cmtasgn calls, respectively.
         call nekasgn(i,j,k,e)
         call cmtasgn(i,j,k,e)
         e_internal=energy(i,j,k) !cmtasgn should do this, but can't
C> Compute thermodynamic state from scalars declared in /NEKUSE/ and /nekuscmt/
C> store state variables in temp,cv,cp,pres and asnd
         call cmt_userEOS(i,j,k,eg)
! JH020718 long-overdue sanity checks
C> Check temperature positivity
         if (temp .lt. 0.0 .and. ilimflag .ne. 0) then
            ifailt=eg
            write(6,'(i6,a26,e12.4,3i2,i8,3e15.6)') ! might want to be less verbose
     >      nid,' HAS NEGATIVE TEMPERATURE ', x,i,j,k,eg,temp,rho,pres
         endif
C> Fill SOLN and CMTDATA arrays one GLL node at a time from scalars
C> in /NEKUSE/ and /nekuscmt/
         vtrans(i,j,k,e,jen)= e_internal
         vtrans(i,j,k,e,jcv)= cv*rho
         vtrans(i,j,k,e,jcp)= cp*rho
         t(i,j,k,e,1)       = temp
         pr(i,j,k,e)        = pres
         csound(i,j,k,e)    = asnd
      enddo
      enddo
      enddo
C> @}
      return
      end

c-----------------------------------------------------------------------

C> \ingroup state bcond initialconds
C> @{
C> \brief Fill /NEKUSE/ and /NEKUSCMT/ common blocks from a single GLL node
C>        extends nekasgn to CMT-nek without affecting core routines
      subroutine cmtasgn (ix,iy,iz,e)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
      include 'NEKUSE'

      integer e,eqnum
!     do eqnum=1,toteq
!        varsic(eqnum)=u(ix,iy,iz,eqnum,e)  
!     enddo
      phi  = phig  (ix,iy,iz,e)
      rho  = vtrans(ix,iy,iz,e,jrho)
      pres = pr    (ix,iy,iz,e)
      if (rho.ne.0) then
         cv   = vtrans(ix,iy,iz,e,jcv)/rho
         cp   = vtrans(ix,iy,iz,e,jcp)/rho
         e_internal = vtrans(ix,iy,iz,e,jen)
      endif
      asnd = csound(ix,iy,iz,e)
      mu     = vdiff(ix,iy,iz,e,jmu)
      udiff  = vdiff(ix,iy,iz,e,jknd)
! MAKE SURE WE''RE NOT USING UTRANS FOR ANYTHING IN pre-v16 code!!
      lambda = vdiff(ix,iy,iz,e,jlam)
C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup initialconds
C> @{
C> \brief set initial values of conserved variables in U for CMT-nek
C>
C> Over-engineered duplicate of setics in core nek5000.
C> Calls cmtuic for a fresh start or my_full_restart for restart.
C> cmtuic actually initializes the flow field through cmt-nek's
C> own dedicated calls to useric.
C> logs min and max of primitive variables as a sanity check in
C> diagnostic I/O labeled "Cuvwpt," etc.
      subroutine cmt_ics
! overlaps with setics. -DCMT will require IFDG as well
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      nxyz2=lx2*ly2*lz2       ! Initialize all fields:
      ntot2=nxyz2*nelv
      nxyz1=lx1*ly1*lz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
      ltott=lelt*nxyz1
      ntotcv=lelt*nxyz1*toteq
      call rone(phig,ltott)
      call rzero(csound,ltott)
      call rzero(vtrans,ltott*ldimt1)
      call rzero(vdiff ,ltott*ldimt1)
      call rzero(u,ntotcv)
! JH100919 where does particle stuff live these days?
!     call usr_particles_init
      call cmtuic
      if(ifrestart) call my_full_restart !  Check restart files. soon...

C print min values
      xxmax = glmin(xm1,ntott)
      yymax = glmin(ym1,ntott)
      zzmax = glmin(zm1,ntott)

      vxmax = glmin(vx,ntotv)
      vymax = glmin(vy,ntotv)
      vzmax = glmin(vz,ntotv)
      prmax = glmin(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmin(t ,ntott)

      if (nio.eq.0) then
         write(6,19) xxmax,yymax,zzmax
   19    format('Cxyz min  ',5g25.18)
      endif
      if (nio.eq.0) then
         write(6,20) vxmax,vymax,vzmax,prmax,ttmax
   20    format('Cuvwpt min',5g25.18)
      endif

c print max values
      xxmax = glmax(xm1,ntott)
      yymax = glmax(ym1,ntott)
      zzmax = glmax(zm1,ntott)

      vxmax = glmax(vx,ntotv)
      vymax = glmax(vy,ntotv)
      vzmax = glmax(vz,ntotv)
      prmax = glmax(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmax(t ,ntott)

      if (nio.eq.0) then
         write(6,16) xxmax,yymax,zzmax
   16    format('Cxyz max  ',5g25.18)
      endif

      if (nio.eq.0) then
         write(6,17) vxmax,vymax,vzmax,prmax,ttmax
   17    format('Cuvwpt max',5g25.18)
      endif

c     ! save velocity on fine mesh for dealiasing
!     call setup_convect(2) ! check what this does again. might be a good
!                           ! idea, or it might be counterproductive
      if(nio.eq.0) then
        write(6,*) 'done :: set initial conditions, CMT-nek'
        write(6,*) ' '
      endif
C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup initialconds
C> @{
C> \brief Fresh initialization of conserved variables from useric. 
C>
C> Calls cmtasgn
C> to interface with userbc, forms conserved variables from
C> scalar primitive variables one grid point at a time, and fills
C> U completely.
      subroutine cmtuic
! overlaps with setics. -DCMT will require IFDG as well
! need to make sure setics has no effect.
! JH070219 cmtuic now sets U and U alone. EVERYTHING else should come
!          from compute_primitive_variables
      include 'SIZE'
      include 'SOLN'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      integer e,eg
      do e=1,nelt
         eg = lglel(e)
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1           
            call nekasgn (i,j,k,e)
            call cmtasgn (i,j,k,e)
            call useric  (i,j,k,eg)
            phig(i,j,k,e)  = phi ! only sane way to run CMT-nek without
                                 ! particles is to have useric set phi=1
            u(i,j,k,irg,e) = phi*rho
            u(i,j,k,irpu,e)= phi*rho*ux
            u(i,j,k,irpv,e)= phi*rho*uy
            u(i,j,k,irpw,e)= phi*rho*uz
            u(i,j,k,iret,e)=phi*rho*(e_internal+0.5*(ux**2+uy**2+uz**2))
!            u(i,j,k,imfrac,e)=phi*rho*ps(1)
!c JB080119 multiple species
!               t(i,j,k,e,2) = ps(l)
!c           do l = 2,NPSCAL
!c               t(i,j,k,e,l) = ps(l-1)
!c           enddo
         enddo
         enddo
         enddo
      enddo
C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup state
C> @{
C> \brief if positive posflags, write failure message and exit
      subroutine poscheck(ifail,what)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
      include 'PARALLEL'
      include 'INPUT'
!JH020918 handles reporting, I/O and exit from failed positivity checks
!         in compute_primitive_variables
      character*11 what

      ifail0=iglmax(ifail,1)
      if(ifail0 .ne. -1) then
         if (nio .eq. 0)
     >   write(6,*) 'dumping solution after negative ',what,'@ eg=',
     >             ifail0
         ifxyo=.true.
!        call out_fld_nek
         call outpost2(vx,vy,vz,pr,t,ldimt,'EBL')
         call exitt
      endif
C> @}
      return
      end
