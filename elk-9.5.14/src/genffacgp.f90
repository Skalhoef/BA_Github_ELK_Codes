
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genffacgp
! !INTERFACE:
pure subroutine genffacgp(is,gpc,ffacgp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is     : species number (in,integer)
!   gpc    : length of G+p-vectors (in,real(ngtot))
!   ffacgp : form factors (out,real(ngtot))
! !DESCRIPTION:
!   Generates the form factors used to determine the smooth characteristic
!   function. See {\tt gencfun} for details.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: gpc(ngtot)
real(8), intent(out) :: ffacgp(ngtot)
! local variables
integer ig
real(8) t1,t2
t1=fourpi/omega
do ig=1,ngtot
  if (gpc(ig) > epslat) then
    t2=gpc(ig)*rmt(is)
    ffacgp(ig)=t1*(sin(t2)-t2*cos(t2))/(gpc(ig)**3)
  else
    ffacgp(ig)=(t1/3.d0)*rmt(is)**3
  end if
end do
end subroutine
!EOC

