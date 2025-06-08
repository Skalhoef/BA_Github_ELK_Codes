
! May 2025, The Half Blood Mathematician 
! Adapted from bandstr.f90 to spintexture1d.f90

!BOP
! !ROUTINE: spintexture
! !INTERFACE:
subroutine spintexture1d
! !USES:
use modmain
use modomp
use modsebbe
! !DESCRIPTION:
!   Produces a spintexture-structure along the path in reciprocal space which connects
!   the vertices in the array {\tt vvlp1d}. The spintexture-structure is obtained from
!   the second-variational eigenstates and is written to the file {\tt C_EVECSV_INFO.OUT}
!   Vertex location lines are written to {\tt BANDLINES.OUT}.
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn,is,ia,ias
integer lmax,lmmax,l,m,lm,iv,nthd
real(8) emin,emax,sm,t1
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:),e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:),apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! initialise universal variables
call init0
call init1
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))


! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr

! generate the spin-orbit coupling radial functions
call gensocfr

emin=1.d5
emax=-1.d5
! begin parallel loop over k-points
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP PRIVATE(dmat,apwalm,ist,ispn) &
!$OMP PRIVATE(ias,l,m,lm,sm) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))

!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL(bandstr_1)
  write(*,'("Info(spintexture): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(bandstr_1)
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
  
  write(91, '("k-point = ", I4)') ik
  
  ! Print evecsv row-by-row
  do iseb = 1, nstsv
    do jseb = 1, nstsv
        write(91, '(2G18.10)', advance='no') real(evecsv(iseb,jseb)), aimag(evecsv(iseb,jseb))
    end do
    write(91,*)  ! Newline after each row
  end do
  

  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ik)=evalsv(ist,ik)-efermi
!$OMP CRITICAL(bandstr_2)
    emin=min(emin,e(ist,ik))
    emax=max(emax,e(ist,ik))
!$OMP END CRITICAL(bandstr_2)
  end do
! end loop over k-points
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)

!$OMP END PARALLEL
call freethd(nthd)





t1=(emax-emin)*0.5d0
emin=emin-t1
emax=emax+t1



! output the band structure for post-processing
if (task == 20 .or. task == 24) then
  open(50,file='BAND.OUT',form='FORMATTED',action='WRITE')
  do ist=1,nstsv
    do ik=1,nkpt
      write(50,'(2G18.10)') dpp1d(ik),e(ist,ik)
    end do
    write(50,*)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(bandstr):")')
  write(*,'(" band structure plot written to BAND.OUT")')
end if



write(*,*)
write(*,'(" Fermi energy is at zero in plot")')
! output the vertex location lines
open(50,file='BANDLINES.OUT',form='FORMATTED',action='WRITE')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,*)
end do
close(50)
write(*,*)
write(*,'(" Vertex location lines written to BANDLINES.OUT")')
deallocate(e)


end subroutine
!EOC
