      subroutine bcflux(flux,agradu)
! Need proper indexing and nekasgn & cmtasgn calls
      include 'SIZE'
      include 'INPUT'
      include 'DG'
!     include 'NEKUSE'
      include 'TSTEP' ! wait how do we know what ifield is?
      integer e,eq,f
      real flux(nx1*nz1,2*ndim,nelt,toteq)
      real agradu(nx1*nz1,2*ndim,nelt,toteq)
      common /nekcb/ cb
      character*3 cb

      nfaces=2*ndim
      nxz=nx1*nz1
      ifield=1

      do e=1,nelt
         do f=1,nfaces
            if (cbc(f, e, ifield).ne.'E  '.and.
     >          cbc(f, e, ifield).ne.'P  ') then ! cbc bndy
               cb=cbc(f,e,ifield)
               if (cb .eq. 'I  ') then ! NEUMANN CONDITIONS GO HERE
!-------------------------------------------------------------
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
! METHOD "B", ADIABATIC NO-SLIP
                  call rzero(flux(1,f,e,1),nxz)
                  do eq=2,ndim+1
                     call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
                  enddo
                  call rzero(flux(1,f,e,toteq),nxz)
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
!                 cbu=cb
!                 do eq=1,toteq
!                    call userflux(flux(1,f,e,eq)) ! replace this with userbc
!                 enddo
               elseif (cb .eq. 'SYM') then
                  do eq=1,toteq
                     call rzero(flux(1,f,e,eq),nxz)
                  enddo
               endif
            endif
         enddo
      enddo

      return
      end
