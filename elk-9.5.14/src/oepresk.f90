
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepresk(ik,vclcv,vclvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: vclcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(in) :: vclvv(nstsv,nstsv,nkpt)
! local variables
integer ist,jst,idm
integer is,ia,ias,ic,m
integer nrc,nrci,npc
real(8) de
complex(8) z1,z2
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),wfcr(:,:)
complex(4), allocatable :: cfmt1(:),cvfmt1(:,:)
complex(4), allocatable :: cfmt2(:,:),cfir2(:)
complex(4), allocatable :: cvfmt2(:,:,:),cvfir2(:,:)
! external functions
complex(8), external :: rcfinp,rcfmtinp
! get the eigenvalues/vectors from file for input k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfir(ngtot,nspinor,nstsv))
call genwfsv_sp(.false.,.false.,nstsv,[0],ngridg,igfft,ngk(1,ik),igkig(:,1,ik),&
 apwalm,evecfv,evecsv,wfmt,ngtot,wfir)
deallocate(apwalm,evecfv,evecsv)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
allocate(wfcr(npcmtmax,2),cfmt1(npcmtmax))
if (spinpol) allocate(cvfmt1(npcmtmax,ndmag))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ic=0
    do ist=1,nstsp(is)
      if (spcore(ist,is)) then
        do m=-ksp(ist,is),ksp(ist,is)-1
          ic=ic+1
! pass in m-1/2 to wavefcr
          call wavefcr(.false.,lradstp,is,ia,ist,m,npcmtmax,wfcr)
          do jst=1,nstsv
            if (evalsv(jst,ik) > efermi) then
              if (spinpol) then
! compute the complex density and magnetisation
                call gencrm(npc,wfcr,wfcr(:,2),wfmt(:,ias,1,jst), &
                 wfmt(:,ias,2,jst),cfmt1,npcmtmax,cvfmt1)
              else
! compute the complex density
                cfmt1(1:npc)=conjg(wfcr(1:npc,1))*wfmt(1:npc,ias,1,jst)
              end if
              z1=conjg(vclcv(ic,ias,jst,ik))
              z2=rcfmtinp(nrc,nrci,wrcmt(:,is),vxmt(:,ias),cfmt1)
              z1=z1-conjg(z2)
              do idm=1,ndmag
                z2=rcfmtinp(nrc,nrci,wrcmt(:,is),bxmt(:,ias,idm),cvfmt1(:,idm))
                z1=z1-conjg(z2)
              end do
              de=evalcr(ist,ias)-evalsv(jst,ik)
              z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! residuals for exchange potential and field
!$OMP CRITICAL(oepresk_)
              call rcadd(npc,z1,cfmt1,dvxmt(:,ias))
              do idm=1,ndmag
                call rcadd(npc,z1,cvfmt1(:,idm),dbxmt(:,ias,idm))
              end do
!$OMP END CRITICAL(oepresk_)
! end loop over jst
            end if
          end do
        end do
! end loop over ist
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(wfcr,cfmt1)
if (spinpol) deallocate(cvfmt1)
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
allocate(cfmt2(npcmtmax,natmtot),cfir2(ngtot))
if (spinpol) then
  allocate(cvfmt2(npcmtmax,natmtot,ndmag),cvfir2(ngtot,ndmag))
end if
do ist=1,nstsv
  if (evalsv(ist,ik) < efermi) then
    do jst=1,nstsv
      if (evalsv(jst,ik) > efermi) then
        if (spinpol) then
! compute the complex density and magnetisation
          call gencfrm(wfmt(:,:,1,ist),wfmt(:,:,2,ist),wfir(:,1,ist), &
           wfir(:,2,ist),wfmt(:,:,1,jst),wfmt(:,:,2,jst),wfir(:,1,jst), &
           wfir(:,2,jst),cfmt2,cfir2,cvfmt2,cvfir2)
        else
! compute the complex density
          call gencrho(.false.,.true.,ngtot,wfmt(:,:,:,ist),wfir(:,:,ist), &
           wfmt(:,:,:,jst),wfir(:,:,jst),cfmt2,cfir2)
        end if
        z1=conjg(vclvv(ist,jst,ik))
        z2=rcfinp(vxmt,vxir,cfmt2,cfir2)
        z1=z1-conjg(z2)
        do idm=1,ndmag
          z2=rcfinp(bxmt(:,:,idm),bxir(:,idm),cvfmt2(:,:,idm),cvfir2(:,idm))
          z1=z1-conjg(z2)
        end do
        de=evalsv(ist,ik)-evalsv(jst,ik)
        z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! add to residuals for exchange potential and field
!$OMP CRITICAL(oepresk_)
        call rcfadd(z1,cfmt2,cfir2,dvxmt,dvxir)
        do idm=1,ndmag
          call rcfadd(z1,cvfmt2(:,:,idm),cvfir2(:,idm),dbxmt(:,:,idm), &
           dbxir(:,idm))
        end do
!$OMP END CRITICAL(oepresk_)
! end loop over jst
      end if
    end do
! end loop over ist
  end if
end do
deallocate(wfmt,wfir,cfmt2,cfir2)
if (spinpol) deallocate(cvfmt2,cvfir2)

contains

pure subroutine rcfadd(za,cfmt,cfir,rfmt,rfir)
implicit none
! arguments
complex(8), intent(in) :: za
complex(4), intent(in) :: cfmt(npcmtmax,natmtot),cfir(ngtot)
real(8), intent(inout) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias
do ias=1,natmtot
  is=idxis(ias)
  call rcadd(npcmt(is),za,cfmt(:,ias),rfmt(:,ias))
end do
call rcadd(ngtot,za,cfir,rfir)
end subroutine

pure subroutine rcadd(n,za,cv,rv)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: za
complex(4), intent(in) :: cv(n)
real(8), intent(out) :: rv(n)
! local variables
real(8) t1
t1=dble(za)
if (abs(t1) > 1.d-12) rv(:)=rv(:)+t1*real(cv(:))
t1=-aimag(za)
if (abs(t1) > 1.d-12) rv(:)=rv(:)+t1*aimag(cv(:))
end subroutine

end subroutine

