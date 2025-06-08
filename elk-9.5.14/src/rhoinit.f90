
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,is,ia,ias,nthd
integer nr,nri,nro,nrs,iro,ir
integer nrc,nrci,irco,irc
integer l,lm,i0,i1,ig,ifg
real(8) x,sm,t1,t2
complex(8) z1,z2,z3
! automatic arrays
real(8) ffg(ngvc),wr(nrspmax),jl(0:lmaxi,nrcmtmax)
complex(8) zfmt(npcmtmax)
! allocatable arrays
complex(8), allocatable :: zfft(:)
lmax=min(lmaxi,1)
! compute the superposition of all the atomic density tails
allocate(zfft(ngtot))
zfft(:)=0.d0
call holdthd(nspecies,nthd)

! SK: Added Write-statement.
!write(99, '()')
!write(99, '()')

!write(99,'("-- SK 20.1 -- rhoinit.f90 : We are now initializing the density.")')
!write(99, '()')
!write(99,'("-- SK 20.1 -- rhoinit.f90 : We entering a for-loop over the atomic species now.")')
!write(99, '()')
!flush(99)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ffg,wr,nr,nrs,nro,ig) &
!$OMP PRIVATE(t1,t2,sm,ir,x,ia,ias,ifg) &
!$OMP NUM_THREADS(nthd)

do is=1,nspecies
! SK: Added write statement.

  !$OMP CRITICAL
  !write(99,'(A,I0)') "    -- SK 20.2 -- rhoinit.f90 : Processing species is = ", is
  !flush(99)
  !$OMP END CRITICAL
  
  nr=nrmt(is)
  nrs=nrsp(is)
  nro=nrs-nr+1
  
  ! SK: Added write statements for nrmt, nrsp, and nro.
  !$OMP CRITICAL
  !write(99, '(A,I6,A,I6)') "    -- SK 20.3 -- rhoinit.f90 : nr = nrmt(", is, ") = ", nr
  !write(99, '(A,I6,A,I6)') "    -- SK 20.4 -- rhoinit.f90 : nrs = nrsp(", is, ") = ", nrs
  !write(99, '(A,I0)') "    -- SK 20.5 -- rhoinit.f90 : nro = nrs - nr + 1 = ", nro
  !write(99, '()')
  !flush(99)
  !$OMP END CRITICAL
  

  ! SK: Added write statements.
  !$OMP CRITICAL
  !write(99,'("    -- SK 20.6 -- rhoinit.f90 : We are calling wsplint now with the following arguments:")')
  !write(99, '(A,I6)') "    -- SK 20.6 -- rhoinit.f90 : Value of nro: ", nro
  !write(99, '(A,F12.6)') "    -- SK 20.6 -- rhoinit.f90 : Value of rsp(nr, is): ", rsp(nr, is)
  !write(99, '(A,F12.6)') "    -- SK 20.6 -- rhoinit.f90 : Value of wr(nr): ", wr(nr)
  !write(99, '()')
  !flush(99)
  !$OMP END CRITICAL
  ! determine the weights for the radial integral
  call wsplint(nro,rsp(nr,is),wr(nr))
  !$OMP CRITICAL
  !write(99,'("    -- SK 20.7 -- rhoinit.f90 : We are done calling wsplint now.")')
  !write(99, '(A,F12.6)') "    -- SK 20.6 -- rhoinit.f90 : Value of wr(nr): ", wr(nr)
  !write(99,'("    -- SK 20.8 -- rhoinit.f90 : Question: What did wsplint do?")')

  !write(99, '()')
  !write(99, '()')
  !write(99,'("    -- SK 20.9 -- rhoinit.f90 : We entering another for-loop over ngvc (???).")')
  !write(99, '(A,I6)') "    -- SK 20.10 -- rhoinit.f90 : Value of ngvc: ", ngvc

  !write(99, '()')
  !flush(99)
  !$OMP END CRITICAL
  
  do ig=1,ngvc
    t1=gc(ig)
    !$OMP CRITICAL
    !write(99,'(A,I0)') "        -- SK 20.10.1 -- rhoinit.f90 : Processing ig = ", ig
    !write(99,'(A,F12.6)') "        -- SK 20.10.2 -- rhoinit.f90 : t1 = gc(ig) = ", t1
    !write(99,'(A,F12.6)') "        -- SK 20.10.3 -- rhoinit.f90 : epslat = ", epslat
    !write(99, '()')
    !flush(99)
    !$OMP END CRITICAL
! spherical bessel function j_0(x) times the atomic density tail
    if (t1 > epslat) then
      !$OMP CRITICAL
      !write(99,'("        -- SK 20.10.4 -- rhoinit.f90 : t1 > epslat")')
      !$OMP END CRITICAL
      t2=1.d0/t1
      sm=0.d0
      do ir=nr,nrs
        x=t1*rsp(ir,is)
        sm=sm+t2*sin(x)*rhosp(ir,is)*rsp(ir,is)*wr(ir)
        !$OMP CRITICAL
        !write(99,'(A,I0)') "            -- SK 20.10.4.1 -- rhoinit.f90 : Processing ir = ", ir
        !write(99,'(A,F12.6)') "            -- SK 20.10.4.2 -- rhoinit.f90 : rsp(ir,is) = ", rsp(ir,is)
        !write(99,'(A,F12.6)') "            -- SK 20.10.4.3 -- rhoinit.f90 : x = t1 * rsp(ir,is) = ", x
        !write(99,'(A,F12.6)') "            -- SK 20.10.4.4 -- rhoinit.f90 : rhosp(ir,is) = ", rhosp(ir,is)
        !write(99, '(A, F12.6)') "            -- SK 20.10.4.5 -- rhoinit.f90 : epslat = ", epslat
        !write(99, '()')
        !flush(99)
        !$OMP END CRITICAL
      end do
    else
      sm=sum(rhosp(nr:nrs,is)*(rsp(nr:nrs,is)**2)*wr(nr:nrs))
    end if
! apply low-pass filter
    t1=sm*exp(-4.d0*(gc(ig)/gmaxvr)**2)
    ffg(ig)=(fourpi/omega)*t1
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!$OMP CRITICAL(rhoinit_)
    do ig=1,ngvc
      ifg=igfft(ig)
      zfft(ifg)=zfft(ifg)+ffg(ig)*conjg(sfacg(ig,ias))
    end do
!$OMP END CRITICAL(rhoinit_)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! compute the tails in each muffin-tin
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(jl,zfmt,is,nrc,nrci) &
!$OMP PRIVATE(irco,ig,ifg,irc,x) &
!$OMP PRIVATE(z1,z2,z3,lm,l,i0,i1) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  irco=nrci+1
  zfmt(1:npcmt(is))=0.d0
  do ig=1,ngvc
    ifg=igfft(ig)
    do irc=1,nrc
      x=gc(ig)*rcmt(irc,is)
      call sbessel(lmax,x,jl(:,irc))
    end do
    z1=fourpi*zfft(ifg)*sfacg(ig,ias)
    do l=0,lmax
      z2=z1*zil(l)
      do lm=l**2+1,(l+1)**2
        z3=z2*conjg(ylmg(lm,ig))
        i1=lmmaxi*(nrci-1)+lm
        zfmt(lm:i1:lmmaxi)=zfmt(lm:i1:lmmaxi)+jl(l,1:nrci)*z3
        i0=i1+lmmaxi
        i1=lmmaxo*(nrc-irco)+i0
        zfmt(i0:i1:lmmaxo)=zfmt(i0:i1:lmmaxo)+jl(l,irco:nrc)*z3
      end do
    end do
  end do
  call ztorfmt(nrc,nrci,zfmt,rhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
t1=chgexs/omega
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  i1=lmmaxi*(nri-1)+1
  rhomt(1:i1:lmmaxi,ias)=rhomt(1:i1:lmmaxi,ias)+(t1+rhosp(1:nri,is))/y00
  i0=i1+lmmaxi
  i1=lmmaxo*(nr-iro)+i0
  rhomt(i0:i1:lmmaxo,ias)=rhomt(i0:i1:lmmaxo,ias)+(t1+rhosp(iro:nr,is))/y00
end do
! interstitial density determined from the atomic tails and excess charge
call zfftifc(3,ngridg,1,zfft)
do ir=1,ngtot
  rhoir(ir)=dble(zfft(ir))+t1
! make sure that the density is always positive
  if (rhoir(ir) < 1.d-10) rhoir(ir)=1.d-10
end do
deallocate(zfft)
! deallocate rhosp as it is not used again
deallocate(rhosp)
end subroutine
!EOC

