
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagk
! !INTERFACE:
subroutine rhomagk(ngp,igpig,lock,wppt,occsvp,apwalm,evecfv,evecsv)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   lock   : OpenMP lock for each atom (in,integer(natmtot))
!   wppt   : weight of input p-point (in,real)
!   occsvp : occupation number for each state (in,real(nstsv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density and magnetisation from the
!   eigenvectors at a particular $k$-point. In the muffin-tin region, the
!   wavefunction is obtained in terms of its $(l,m)$-components from both the
!   APW and local-orbital functions. Using a backward spherical harmonic
!   transform (SHT), the wavefunction is converted to real-space and the density
!   obtained from its modulus squared. A similar proccess is used for the
!   intersitial density in which the wavefunction in real-space is obtained from
!   a Fourier transform of the APW functions. See routines {\tt wfmtsv},
!   {\tt genshtmat} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed conversion to spherical harmonics, January 2009 (JKD)
!   Partially de-phased the muffin-tin magnetisation for spin-spirals,
!    February 2009 (FC, FB & LN)
!   Optimisations, July 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
integer(omp_lock_kind), intent(inout) :: lock(natmtot)
real(8), intent(in) :: wppt,occsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,nst,ist
integer is,ias,npc,i,j,k
integer n,igp,nthd
real(8) wo,ts0,ts1
complex(8) z1
! automatic arrays
integer idx(nstsv)
complex(8) wfir(ngtc,nspinor),wfgp(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:)

! SK.
integer :: iseb, jseb, kseb, lseb, mseb
! SK End

call timesec(ts0)


! SK: Added Write-statement.
write(99,'("      Variables, Values and Meanings: ")')
write(99, '()')

write(99, '("      ngkmax = maximal number of G+k-vectors for all k.")')
write(99, '("             ... cf. also findgnkmax.f90.")')
write(99, '("      ngkmax = ", I10)') ngkmax
write(99, '()')

write(99, '("      natmtot = number of atomic species.")')
write(99, '("      natmtot = ", I10)') natmtot
write(99, '()')

write(99, '("      apwordmax = ", I10)') apwordmax
write(99, '()')

write(99, '("      lmmaxapw = (lmaxapw+1)^2.")')
write(99, '("               = ( 8 + 1)^2.")')
write(99, '("      lmmaxapw = ", I10)') lmmaxapw
write(99, '()')

write(99, '("      nmatmax = maximum nmat over all k-points.")')
write(99, '("              ... is determined in init1.f90.")')
write(99, '("      nmatmax = ", I10)') nmatmax
write(99, '()')

write(99, '("      nstfv = number of first-variational states.")')
write(99, '("            = number of computed (not necc. occupied) bands.")')
write(99, '("      nstfv = ", I10)') nstfv
write(99, '()')

write(99, '("      nstsv = number of second-variational states.")')
write(99, '("            = 2 * nstfv")')
write(99, '("      nstsv = ", I10)') nstsv
write(99, '()')


write(99, '("      nspnfv = number of spin-dependent first-variational functions per state (?).")')
write(99, '("             ... seems somehow like a weird variable...")')
write(99, '("      nspnfv = ", I10)') nspnfv
write(99, '()')

write(99, '("      ngp =^ number of G+p-vectors (in,integer(nspnfv)).")')
write(99, '("      ngp array values:")')
write(99, '()')

do iseb = 1, nspnfv
    write(99, '("      ngp(", I3, ") = ", I10)') iseb, ngp(iseb)
end do
write(99, '()')

write(99, '("      wppt = ", F12.6)') wppt
write(99, '()')


!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
write(99,'("      We print igpg now.")')
write(99,'("      Reminder: igpg = : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv)).???")')
write(99, '()')

do jseb = 1, nspnfv
  do iseb = 1, ngkmax
    write(98, '(I10)', advance="no") igpig(iseb, jseb)
  end do
  write(98, *)  ! New line after each row
end do

do jseb = 1, nspnfv
  do iseb = 1, ngp(nspnfv)
    write(96, '(I10)', advance="no") igpig(iseb, jseb)
  end do
  write(96, *)  ! New line after each row
end do

do iseb = 1, nstsv
    write(97, '(F12.6)') occsvp(iseb)  ! Print each value with 6 decimal places
end do

write(99,'("      We entering a for-loop over nstsv now in order to count occupied states")')
flush(99)


!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
! number of and index to occupied states
nst=0
do ist=1,nstsv
  if (abs(occsvp(ist)) < epsocc) cycle
  nst=nst+1
  idx(nst)=ist
end do

write(99,'("      We are done with the for-loop.")')
write(99, '()')


write(99, '("      nst = number of occupied bands / states.")')
write(99, '("      nst = ", I10)') nst
write(99, '()')

write(99, '("      idx = index to the states / bands which are occupied.")')
write(99, '("          = an array which seems extremely useless.")')

do iseb = 1, nst
    write(99, '("idx(", I3, ") = ", I10)') iseb, idx(iseb)
end do
write(99, '()')


write(99, '("      We now allocate the wavefunction-matrix wfmt(npcmtmax,nspinor,nst) with ...")')
write(99, '()')

write(99, '("      npcmtmax = maximum number of points over all packed muffin-tins.")')
write(99, '("               = a variable that is SO WEIRDLY defined, it seems like a joke...")')
write(99, '("                 ... Cf. init0.f90 and genrmesh.f90")')
write(99, '("      npcmtmax = ", I10)') npcmtmax
write(99, '()')

write(99, '("      nspinor = second-variational spinor dimension (1 or 2).")')
write(99, '("              = a variable that seems so weird to introduce...")')
write(99, '("       ... since the second variation takes care of SOC, so how can nspinor be one???")')
write(99, '("      nspinor = ", I10)') nspinor
write(99, '()')

write(99, '("      nst = number of occupied bands / states.")')
write(99, '("      nst = ", I10)') nst
write(99, '()')

allocate(wfmt(npcmtmax,nspinor,nst))


write(99, '("      We entering a for-loop now with ias between 1 and natmtot = ", I10)') natmtot
write(99, '("      Recall for this: idxis(maxatoms*maxspecies = 8*200) = atoms and species indices.")')
write(99, '("      ... Remark: this variable is closely related to another variable.")')
write(99, '()')

write(99, '()')
write(99, '("      We now compute the density- and the magnetization-contributions.")')
write(99, '()')

call holdthd(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,wfgp,ias,is) &
!$OMP PRIVATE(npc,i,j,k,wo) &
!$OMP PRIVATE(ispn,jspn,n,ist) &
!$OMP PRIVATE(z1,igp) &
!$OMP NUM_THREADS(nthd)

! If I put the print-statement from above here, I get it twice in the output-file.
! ---> Can be cured by restricting the print-statement to a single thread.

do ias=1,natmtot
  is=idxis(ias)
  npc=npcmt(is)
!$OMP SINGLE
  
  ! SK Comments
  ! lradstp was one of the convergence-generating-parameters!
  write(99, '("      ias = ", I10)') ias
  write(99, '("      is = idxis(ias) = ", I10)') is
  write(99, '("      npc = npcmt(is) = ", I10)') npc
  write(99, '("      lradstp = ", I10)') lradstp
  write(99, '()')
  !
  write(99,'("      [ Subroutine wfmtsv.f90 started. ]")')
  call wfmtsv(.false.,lradstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,npcmtmax,&
   wfmt)
  write(99,'("      [ Subroutine wfmtsv.f90 finished.]")')
  write(99, '()')
!$OMP END SINGLE
!$OMP DO
  do j=1,nst
    k=idx(j)
    wo=occsvp(k)*wppt
! add to density and magnetisation
    call omp_set_lock(lock(ias))
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        call rmk1(npc,wo,wfmt(:,1,j),wfmt(:,2,j),rhomt(:,ias),magmt(:,ias,1), &
         magmt(:,ias,2),magmt(:,ias,3))
      else
! collinear
        call rmk2(npc,wo,wfmt(:,1,j),wfmt(:,2,j),rhomt(:,ias),magmt(:,ias,1))
      end if
    else
! spin-unpolarised
      call rmk3(npc,wo,wfmt(:,1,j),rhomt(:,ias))
    end if
    call omp_unset_lock(lock(ias))
  end do
!$OMP END DO
end do
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
!$OMP DO
do j=1,nst
  k=idx(j)
  wo=occsvp(k)*wppt/omega
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      n=ngp(jspn)
      wfgp(1:n)=0.d0
      do ist=1,nstfv
        i=(ispn-1)*nstfv+ist
        z1=evecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)) > epsocc) then
          wfgp(1:n)=wfgp(1:n)+z1*evecfv(1:n,ist,jspn)
        end if
      end do
      wfir(:,ispn)=0.d0
      do igp=1,n
        wfir(igfc(igpig(igp,jspn)),ispn)=wfgp(igp)
      end do
! Fourier transform wavefunction to real-space
      call zfftifc(3,ngdgc,1,wfir(:,ispn))
    end do
  else
! spin-unpolarised wavefunction
    wfir(:,1)=0.d0
    do igp=1,ngp(1)
      wfir(igfc(igpig(igp,1)),1)=evecfv(igp,k,1)
    end do
    call zfftifc(3,ngdgc,1,wfir)
  end if
! add to density and magnetisation
!$OMP CRITICAL(rhomagk_)
  if (spinpol) then
! spin-polarised
    if (ncmag) then
! non-collinear
      call rmk1(ngtc,wo,wfir,wfir(:,2),rhoir,magir,magir(:,2),magir(:,3))
    else
! collinear
      call rmk2(ngtc,wo,wfir,wfir(:,2),rhoir,magir)
    end if
  else
! spin-unpolarised
    call rmk3(ngtc,wo,wfir,rhoir)
  end if
!$OMP END CRITICAL(rhomagk_)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
deallocate(wfmt)
call timesec(ts1)
!$OMP ATOMIC
timerho=timerho+ts1-ts0
return

contains

pure subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) wo2,t1,t2
real(8) a1,b1,a2,b2
wo2=2.d0*wo
!$OMP SIMD PRIVATE(a1,b1,a2,b2,t1,t2) SIMDLEN(8)
do i=1,n
  a1=dble(wf1(i)); b1=aimag(wf1(i))
  a2=dble(wf2(i)); b2=aimag(wf2(i))
  t1=a1**2+b1**2; t2=a2**2+b2**2
  mag1(i)=mag1(i)+wo2*(a1*a2+b1*b2)
  mag2(i)=mag2(i)+wo2*(a1*b2-b1*a2)
  mag3(i)=mag3(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
end subroutine

pure subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(8) t1,t2
!$OMP SIMD PRIVATE(t1,t2) SIMDLEN(8)
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  mag(i)=mag(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
end subroutine

pure subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(:)=rho(:)+wo*(dble(wf(:))**2+aimag(wf(:))**2)
end subroutine

end subroutine
!EOC

