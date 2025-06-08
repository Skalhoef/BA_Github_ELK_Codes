
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwchgk(ik,chgk)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: chgk
! local variables
integer ist,iw,i
real(8) e
complex(8) z1
! allocatable arrays
complex(8), allocatable :: se(:,:,:),gs(:),g(:,:),ge(:,:)
! external functions
complex(8), external :: gwtails
! read the self-energy from file
allocate(se(nstsv,nstsv,0:nwfm))
call getgwsefm(ik,se)
! allocate local arrays
allocate(gs(nstsv),g(nstsv,nstsv),ge(4,nstsv))
chgk=0.d0
do iw=0,nwfm
! compute the diagonal matrix G_s
  do ist=1,nstsv
    e=evalsv(ist,ik)-efermi
    gs(ist)=1.d0/(wfm(iw)-e)
  end do
! compute 1 - G_s Sigma
  do ist=1,nstsv
    z1=-gs(ist)
    g(ist,:)=z1*se(ist,:,iw)
    g(ist,ist)=g(ist,ist)+1.d0
  end do
! invert this matrix
  call zminv(nstsv,g)
! take the trace of G = (1 - G_s Sigma)^(-1) G_s
  do ist=1,nstsv
    g(ist,ist)=g(ist,ist)*gs(ist)
    chgk=chgk+dble(g(ist,ist))
  end do
! store the Green's function at the end point frequencies
  i=0
  if (iw == 0) i=1
  if (iw == 1) i=2
  if (iw == nwfm-1) i=3
  if (iw == nwfm) i=4
  if (i /= 0) then
    do ist=1,nstsv
      ge(i,ist)=g(ist,ist)
    end do
  end if
end do
! add the Matsubara tails analytically
do ist=1,nstsv
  chgk=chgk+dble(gwtails(ge(:,ist)))
end do
chgk=chgk*wkpt(ik)*occmax*kboltz*tempk
deallocate(se,gs,g,ge)
end subroutine

