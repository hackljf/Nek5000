      subroutine bcflux(flux)
! Need proper indexing and nekasgn & cmtasgn calls
      include 'SIZE'
      include 'INPUT'
      include 'DG'
!     include 'NEKUSE'
      include 'TSTEP' ! wait how do we know what ifield is?
      integer e,eq,f
      real flux(nx1*nz1,2*ndim,nelt,toteq)
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
               if (cb .eq. 'I  ') then
!                 cbu=cb
                  do eq=1,toteq
                     call userflux(flux(1,f,e,eq)) ! replace this with userbc
                  enddo
               endif
            endif
         enddo
      enddo

      return
      end
