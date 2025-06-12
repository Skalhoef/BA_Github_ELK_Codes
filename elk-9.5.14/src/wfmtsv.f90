
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfmtsv(tsh,lrstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,ld,wfmt)
use modmain
use modomp
use modsebbe

implicit none
! arguments

! SK: tsh = logical flag (whether to transform wavefunctions 
! SK: into spherical harmonics or not)
logical, intent(in) :: tsh
! SK: lrstp: determines the radial mesh type 
! SK: (full-potential vs. muffin-tin region)
integer, intent(in) :: lrstp,is,ias,nst,idx(*),ngp(nspnfv)
! SK: APW coefficients
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
! SK: First and second variational eigenvectors
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
! SK: Leading dimension of output wavefunction array.
integer, intent(in) :: ld
! SK: wfmt = output wavefunction matrix
complex(8), intent(out) :: wfmt(ld,nspinor,nst)

! local variables

! SK: tasv = logical flag that determines whether second variational states 
! SK: should be computed or not.
logical tasv




! SK: iro, nr, nri, nro, np and npi are variables that determine the ...
! SK: ... radial grid.

integer io,ilo,ispn,jspn
integer nr,nri,nro,iro
integer l,lm,np,npi
integer n,p,i,j,k,nthd
complex(8) zq(2),z1
! automatic arrays
complex(8) x(nstfv,nspnfv),y(nlmwf(is),nspinor,nst)
! external functions
complex(8), external :: zdotu

! >>>> Add this variable to control print statements <<<<
logical :: print99
print99 = .false.  ! Set to .false. to deactivate all write(99,...) statements


if (print99) then
  write(99, '()')
  write(99,'("       We are now computing the wavefunctions inside the muffin-tins.")')
  write(99, '()')
end if

iro=nrmti(is)+lrstp
if (lrstp == 1) then
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  npi=npmti(is)
else
  nr=nrcmt(is)
  nri=nrcmti(is)
  np=npcmt(is)
  npi=npcmti(is)
end if
nro=nr-nri

! SK Begin
if (print99) then
  write(99, '("       iro = nrmti(is) + lrstp = ", I10)') iro
  write(99, '("       nr = ", I10)') nr
  write(99, '("       nrmti(is) = nri = ", I10)') nri
  write(99, '("       np = ", I10)') np
  write(99, '("       npi = ", I10)') npi
  write(99, '("       nro = nr - nri = ", I10)') nro
  write(99, '()')
  write(99, '()')
end if
! SK End


! de-phasing factor for spin-spirals
if (ssdph) then
  zq(1)=zqss(ias)
  zq(2)=conjg(zq(1))
end if

!
! SK: Is it possible, that this statement, i.e. idx(1) == 0, is N E V E R fulfilled?
!
if (print99) then
  write(99,'("       We are now checking, whether idx(1)==0, which should N E V E R be fulfilled... (Weird?!)")')
  write(99,'("       If it is fulfilled, the code will apparently compute all the second-variational wavefunctions...")')
  write(99, '()')
end if

! check if all the second-variational wavefunctions should be calculated
if (idx(1) == 0) then
  tasv=.true.
else
  tasv=.false.
end if

if (print99) then
  write(99,'("       lmaxo = ", I10)') lmaxo
  write(99,'("       We now loop over local orbitals for APW-functions.")')
  write(99, '()')

  write(99,'("       For is = ", I3, " one has for apword(l, is) the values")') is

  do iseb=0, lmaxo
    write(99,'("       l = ", I3, " ==> apword(l, is) = ", I3)') iseb, apword(iseb, is)
  end do
  write(99, '()')
end if


call holdthd(nst,nthd)
!-----------------------!
!     APW functions     !
!-----------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(p,l,lm,io,ispn,jspn) &
!$OMP PRIVATE(n,i,j,k,z1,ilo) &
!$OMP NUM_THREADS(nthd)
p=0
do l=0,lmaxo
  
  if (print99) then
    iseb = l**2 + 1
    jseb = (l + 1)**2

    write(99,'("       For l = ", I3, " we lm between l**2+1 = ", I10, " and (l+1)**2 = ", I10)') l, iseb, jseb
    write(99, '()')
    write(99,'("       Computing the scalar-product for the first variation.")')
    write(99, '()')
  end if 

  do lm=l**2+1,(l+1)**2
    do io=1,apword(l,is)
      p=p+1
      ! tevecsv = True for spin-polarized calculations.
      if (tevecsv) then
        do jspn=1,nspnfv
          n=ngp(jspn)
!$OMP DO
          do j=1,nstfv
            x(j,jspn)=zdotu(n,evecfv(:,j,jspn),1,apwalm(:,io,lm,ias,jspn),1)
          end do
!$OMP END DO
        end do
! loop only over required states
!$OMP DO
        do j=1,nst
! index to state in evecsv
          if (tasv) then; k=j; else; k=idx(j); end if
          ! write(99,'("       Computing another scalar-product for p = ", I10)') p
          y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)

          ! SK Begin
          ! if (not_printed_zfzrf) then
            ! write(99, '("       Kays Version: y(p,1,j) = ", F12.6, " + i", F12.6)') real(y(p,1,j)), aimag(y(p,1,j))
            ! y(p,1,j)=zdotu(nstfv,evecsv(1:nstfv,k),1,x(:,1),1)
            ! write(99, '("       My Version: y(p,1,j) = ", F12.6, " + i", F12.6)') real(y(p,1,j)), aimag(y(p,1,j))
            ! write(99, '()')
          ! end if
          ! SK End

          if (spinpol) then
            jspn=jspnfv(2)
            ! write(99,'("       jspn = jspnfv(2) = ", I10)') jspn ! is constant to one, since we have no spirals
            y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)

            ! SK Begin
            ! if (not_printed_zfzrf) then
              ! write(99, '("       Kays Version: y(p,2,j) = ", F12.6, " + i", F12.6)') real(y(p,2,j)), aimag(y(p,2,j))
              ! y(p,1,2)=zdotu(nstfv,evecsv(nstfv+1:nstsv,k),1,x(:,1),1)
              ! write(99, '("       My Version: y(p,2,j) = ", F12.6, " + i", F12.6)') real(y(p,2,j)), aimag(y(p,2,j))
              ! not_printed_zfzrf = .false.  ! Prevent future prints
              ! write(99, '()')
            ! end if
            ! SK End


          end if
        end do
!$OMP END DO
      else
!$OMP DO
        do j=1,nst
          if (tasv) then; k=j; else; k=idx(j); end if
          y(p,1,j)=zdotu(ngp(1),evecfv(:,k,1),1,apwalm(:,io,lm,ias,1),1)
        end do
!$OMP END DO
      end if
    end do
  end do
end do
!$OMP DO
do j=1,nst
  wfmt(1:np,:,j)=0.d0
  do ispn=1,nspinor
    p=0
    do l=0,lmaxo
      do lm=l**2+1,(l+1)**2
        i=npi+lm
        do io=1,apword(l,is)
          p=p+1
          z1=y(p,ispn,j)
          if (ssdph) z1=z1*zq(ispn)
          if (l <= lmaxi) then

            ! SK Begin
            ! print only once the input:
            if (print99) then
              if (not_printed_zfzrf) then
                write(99, '("z1  = (", F10.5, ",", F10.5, ")")') real(z1), aimag(z1)
                write(99, '("lmmaxi = ", I10)') lmmaxi
                write(99, '("lmmaxo = ", I10)') lmmaxo
                !write(92, '("rf array:")')
	      !do   iseb = 1, lrstp
                  !write(92, '(1000F10.5)') apwfr(iseb, :)
                !end do
	      !wri  te(92, '("zf array:")')
	      !do   iseb = 1, ld
                  !write(92, '(1000F10.5)') real(wfmt(iseb, :)), aimag(wfmt(iseb, :))
                !end do
                not_printed_zfzrf = .false.  ! Prevent future prints
              end if
            end if 
	    ! SK End

            call zfzrf(nri,z1,apwfr(1,1,io,l,ias),lmmaxi,wfmt(lm,ispn,j))
          end if
          call zfzrf(nro,z1,apwfr(iro,1,io,l,ias),lmmaxo,wfmt(i,ispn,j))
        end do
      end do
    end do
  end do
end do
!$OMP END DO
if (print99) then
  write(99, '()')
  write(99,'("       [ . . . ]")')
  write(99, '()')
end if 

!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
p=0
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    p=p+1
    i=idxlo(lm,ilo,ias)
    if (tevecsv) then
      do jspn=1,nspnfv
        n=ngp(jspn)
        x(1:nstfv,jspn)=evecfv(n+i,1:nstfv,jspn)
      end do
!$OMP DO
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)
        if (spinpol) then
          jspn=jspnfv(2)
          y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)
        end if
      end do
!$OMP END DO
    else
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=evecfv(ngp(1)+i,k,1)
      end do
    end if
  end do
end do
!$OMP DO
do j=1,nst
  do ispn=1,nspinor
    p=0
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do lm=l**2+1,(l+1)**2
        p=p+1
        i=npi+lm
        z1=y(p,ispn,j)
        if (ssdph) z1=z1*zq(ispn)
        if (l <= lmaxi) then
          call zfzrf(nri,z1,lofr(1,1,ilo,ias),lmmaxi,wfmt(lm,ispn,j))
        end if
        call zfzrf(nro,z1,lofr(iro,1,ilo,ias),lmmaxo,wfmt(i,ispn,j))
      end do
    end do
  end do
end do
!$OMP END DO
! convert to spherical coordinates if required
if (.not.tsh) then
!$OMP DO
  do j=1,nst
    do ispn=1,nspinor
      call zbshtip(nr,nri,wfmt(:,ispn,j))
    end do
  end do
!$OMP END DO
end if
!$OMP END PARALLEL
call freethd(nthd)
return

contains

pure subroutine zfzrf(n, z, rf, ld, zf)
    implicit none
    ! arguments
    integer, intent(in) :: n
    complex(8), intent(in) :: z
    real(8), intent(in) :: rf(lrstp, n)
    integer, intent(in) :: ld
    complex(8), intent(inout) :: zf(ld, n)
      
    ! Perform the computation
    zf(1, :) = zf(1, :) + z * rf(1, :)
        
end subroutine zfzrf

end subroutine
