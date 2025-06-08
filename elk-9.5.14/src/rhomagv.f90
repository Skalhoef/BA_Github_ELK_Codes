
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomagv
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ispn,idm
integer is,ias,n,nthd

! SK 
integer, dimension(7) :: unitseb = [99, 98, 97, 96, 95, 94, 93]
integer :: iseb, jseb, kseb, lseb, mseb
! SK End

! automatic arrays
integer(omp_lock_kind) lock(natmtot)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
! set the charge density and magnetisation to zero
do ias=1,natmtot
  is=idxis(ias)
  rhomt(1:npcmt(is),ias)=0.d0
end do
rhoir(1:ngtc)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    magmt(1:npcmt(is),ias,idm)=0.d0
  end do
  magir(1:ngtc,idm)=0.d0
end do
! initialise the OpenMP locks
do ias=1,natmtot
  call omp_init_lock(lock(ias))
end do
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,ispn) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))

do iseb = 1, size(unitseb)
  write(unitseb(iseb), '("    [ Beginning loop over k-points.]")')
  write(unitseb(iseb), '()')
end do

!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! add to the density and magnetisation
  
  ! SK
  !$OMP CRITICAL
  do iseb = 1, size(unitseb)
    write(unitseb(iseb), '()')
    write(unitseb(iseb), '("     ================")')
    write(unitseb(iseb), '("     k-point = ", I4)') ik
    write(unitseb(iseb), '("     ================")')
    write(unitseb(iseb), '()')
  end do
  !$OMP END CRITICAL

  
  ! Print apwalm
  do jseb = 1, nspnfv
    do mseb = 1, natmtot
      do lseb = 1, lmmaxapw
        do kseb = 1, apwordmax
          do iseb = 1, ngkmax
          write(95, '(ES15.8)', advance="no") apwalm(iseb, kseb, lseb, mseb, jseb)
          end do
          write(95, *)  ! New line after each row
        end do
      end do
    end do
  end do

  ! Print evecfv
  do jseb = 1, nspnfv
    do mseb = 1, nstfv
      do lseb = 1, nmatmax
        write(94, '(ES15.8)', advance="no") evecfv(lseb, mseb, jseb)
      end do
      write(94, *)  ! New line after each row
    end do
  end do

  ! Print evecsv
  do kseb = 1, nstsv
    do lseb = 1, nstsv
      write(93, '(ES15.8)', advance="no") evecsv(kseb, lseb)
    end do
    write(93, *)  ! New line after each row
  end do
  
  ! SK
  write(99,'("     [ Subroutine rhomagk.f90 started. ]")')
  call rhomagk(ngk(:,ik),igkig(:,:,ik),lock,wkpt(ik),occsv(:,ik),apwalm, &
   evecfv,evecsv)
  write(99,'("     [ Subroutine rhomagk.f90 finished.]")')
  write(99, '()')
end do
!$OMP END DO

do iseb = 1, size(unitseb)
  write(unitseb(iseb), '("    [ Finished loop over k-points.]")')
  write(unitseb(iseb), '()')
end do


deallocate(apwalm,evecfv,evecsv)
!$OMP END PARALLEL
call freethd(nthd)
! destroy the OpenMP locks
do ias=1,natmtot
  call omp_destroy_lock(lock(ias))
end do
! convert muffin-tin density/magnetisation to spherical harmonics
call rhomagsh
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(idm) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! symmetrise the density
call symrf(nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax,rhomt,rhoir)
! convert the muffin-tin density from coarse to fine radial mesh
call rfmtctof(rhomt)
! convert the interstitial density from coarse to fine grid
call rfirctof(rhoir,rhoir)
!$OMP SECTION
if (spinpol) then
! symmetrise the magnetisation
  call symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax, &
   magmt,ngtot,magir)
! convert the muffin-tin magnetisation from coarse to fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(:,:,idm))
  end do
! convert the interstitial magnetisation from coarse to fine grid
  do idm=1,ndmag
    call rfirctof(magir(:,idm),magir(:,idm))
  end do
end if
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
! add densities from each process and redistribute
if (np_mpi > 1) then
  n=npmtmax*natmtot
  call mpi_allreduce(mpi_in_place,rhomt,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
  call mpi_allreduce(mpi_in_place,rhoir,ngtot,mpi_double_precision,mpi_sum, &
   mpicom,ierror)
  if (spinpol) then
    n=n*ndmag
    call mpi_allreduce(mpi_in_place,magmt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    n=ngtot*ndmag
    call mpi_allreduce(mpi_in_place,magir,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
  end if
end if
end subroutine

