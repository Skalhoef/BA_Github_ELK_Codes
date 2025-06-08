
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfirctof(rfirc,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfirc(ngtc)
real(8), intent(out) :: rfir(ngtot)
! local variables
integer ig,ifg
! automatic arrays
complex(8) zfftc(ngtc)
! allocatable arrays
complex(8), allocatable :: zfft(:)
! Fourier transform function on coarse grid to G-space
zfftc(:)=rfirc(:)
call zfftifc(3,ngdgc,-1,zfftc)
! Fourier transform to fine real-space grid
allocate(zfft(nfgrz))
do ifg=1,nfgrz
  ig=igrzf(ifg)
  if (ig <= ngvc) then
    zfft(ifg)=zfftc(igfc(ig))
  else
    zfft(ifg)=0.d0
  end if
end do
call rzfftifc(3,ngridg,1,rfir,zfft)
deallocate(zfft)
end subroutine

