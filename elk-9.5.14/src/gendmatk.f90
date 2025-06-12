
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatk(tspndg,tlmdg,lmin,lmax,ias,nst,idx,ngp,apwalm,evecfv, &
 evecsv,ld,dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax,ias
integer, intent(in) :: nst,idx(*)
integer, intent(in) :: ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,nst)
! local variables
integer ispn,jspn,ist,is
integer nrc,nrci,irco,irc
integer l,lma,lmb,lm1,lm2
integer npci,i1,i2
complex(8) zsm
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:)

! Print input arguments
write(99,'("gendmatk: tspndg = ",L1)') tspndg
write(99,'("gendmatk: tlmdg = ",L1)') tlmdg
write(99,'("gendmatk: lmin = ",I10)') lmin
write(99,'("gendmatk: lmax = ",I10)') lmax
write(99,'("gendmatk: ias = ",I10)') ias
write(99,'("gendmatk: nst = ",I10)') nst
write(99,'("gendmatk: ld = ",I10)') ld
write(99,'("gendmatk: idx array:")')
do i1=1,nst
  write(99,'("  idx(",I3,") = ",I10)') i1, idx(i1)
end do
write(99,'("gendmatk: ngp array:")')
do i1=1,nspnfv
  write(99,'("  ngp(",I3,") = ",I10)') i1, ngp(i1)
end do

if (lmin < 0) then
  write(*,*)
  write(*,'("Error(gendmatk): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmax > lmaxo) then
  write(*,*)
  write(*,'("Error(gendmatk): lmax > lmaxo : ",2I8)') lmax,lmaxo
  write(*,*)
  stop
end if

is=idxis(ias)
write(99,'("gendmatk: is = idxis(ias) = ",I10)') is

nrc=nrcmt(is)
write(99,'("gendmatk: nrc = nrcmt(is) = ",I10)') nrc

nrci=nrcmti(is)
write(99,'("gendmatk: nrci = nrcmti(is) = ",I10)') nrci

irco=nrci+1
write(99,'("gendmatk: irco = nrci + 1 = ",I10)') irco

npci=npcmti(is)
write(99,'("gendmatk: npci = npcmti(is) = ",I10)') npci

write(99,'("gendmatk: Allocating wfmt(npcmtmax,nspinor,nst), npcmtmax = ",I10, &
     ", nspinor = ",I10,", nst = ",I10)') npcmtmax, nspinor, nst

! generate the second-variational wavefunctions
allocate(wfmt(npcmtmax,nspinor,nst))
call wfmtsv(.true.,lradstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,npcmtmax, &
 wfmt)
write(99,'("gendmatk: wfmtsv called.")')

! zero the density matrix
dmat(:,:,:,:,:)=0.d0

! Print lmmaxi and lmmaxo before the loop
write(99,'("gendmatk: lmmaxi = ",I10)') lmmaxi
write(99,'("gendmatk: lmmaxo = ",I10)') lmmaxo


! loop over second-variational states
do ist=1,nst
  do ispn=1,nspinor
    do jspn=1,nspinor
      if (tspndg.and.(ispn /= jspn)) cycle
      do l=lmin,lmax
        lma=l**2+1; lmb=lma+2*l
        do lm1=lma,lmb
          do lm2=lma,lmb
            if (tlmdg.and.(lm1 /= lm2)) cycle
            if (l <= lmaxi) then
              zsm=0.d0
              i1=lm1; i2=lm2
              do irc=1,nrci
                write(99,'("gendmatk: l = ",I4,"  i1 = ",I6,"  (irc = ",I6,")")') l, i1, irc
                zsm=zsm+wfmt(i1,ispn,ist)*conjg(wfmt(i2,jspn,ist))*wrcmt(irc,is)
                i1=i1+lmmaxi; i2=i2+lmmaxi
              end do
              do irc=irco,nrc
                write(99,'("gendmatk: l = ",I4,"  i1 = ",I6,"  (irc = ",I6,")")') l, i1, irc
                zsm=zsm+wfmt(i1,ispn,ist)*conjg(wfmt(i2,jspn,ist))*wrcmt(irc,is)
                i1=i1+lmmaxo; i2=i2+lmmaxo
              end do
            else
              zsm=0.d0
              i1=npci+lm1; i2=npci+lm2
              do irc=irco,nrc
                write(99,'("gendmatk: l = ",I4,"  i1 = ",I6,"  (irc = ",I6,")")') l, i1, irc
                zsm=zsm+wfmt(i1,ispn,ist)*conjg(wfmt(i2,jspn,ist))*wrcmt(irc,is)
                i1=i1+lmmaxo; i2=i2+lmmaxo
              end do
            end if
            dmat(lm1,ispn,lm2,jspn,ist)=zsm
          end do
        end do
      end do
    end do
  end do
! end loop over second-variational states
end do
deallocate(wfmt)
end subroutine

