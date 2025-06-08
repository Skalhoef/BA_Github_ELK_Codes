
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine trdmatdu(l,lmmax,dm,trdmat)
use modmain
use moddftu
implicit none
! arguments
integer, intent(in) :: l
integer, intent(in) :: lmmax
complex(8), intent(in) :: dm(lmmax,nspinor,lmmax,nspinor)
complex(8), intent(out) :: trdmat(lmmax,nspinor,lmmax,nspinor,0:1)
! local variables
!integer is,ia,ias
integer ip,inu,p,i,j,itr
integer l1,l2,m1,m2,lm1,lm2,lm1_,lm2_
! allocatable arrays
complex(8), allocatable :: ndmat(:,:,:)
complex(8), allocatable :: ndmatnu(:,:,:,:)
! allocate arrays
allocate(ndmat(lmmax,lmmax,0:3))
allocate(ndmatnu(lmmax,lmmax,0:3,0:1))
! zero matrices
      ndmat(:,:,:)=0.d0
      ndmatnu(:,:,:,:)=0.d0
! calculate the n^p matrix, Sp(\rho \sigma^p) 
      do i=1,lmmax
        do j=1,lmmax
          ndmat(i,j,0)=dm(i,1,j,1)+dm(i,2,j,2)
          ndmat(i,j,1)=dm(i,1,j,2)+dm(i,2,j,1)
          ndmat(i,j,2)=zi*(dm(i,1,j,2)-dm(i,2,j,1))
          ndmat(i,j,3)=dm(i,1,j,1)-dm(i,2,j,2)
        end do
      end do
! calculate  n^{p,\nu} matrix, Eq.59 [Lars notes 16-04-2009]
      p=0
! index ip label Pauli matrices 
      do ip=0,3
         if(ip.gt.0) p=1
         do inu=0,1
            do l1=0,l
               do m1=-l1,l1
                  lm1=l1*(l1+1)+1+m1
                  lm1_=l1*(l1+1)+1-m1
                  do l2=0,l
                     do m2=-l2,l2
                        lm2=l2*(l2+1)+1+m2
                        lm2_=l2*(l2+1)+1-m2
                        ndmatnu(lm1,lm2,ip,inu)=0.5d0*(ndmat(lm1,lm2,ip)+((-1)**(p+inu)) & 
                             *((-1)**(m1+m2))*ndmat(lm2_,lm1_,ip))
                     end do
                  end do
               end do
            end do
         end do
      end do
! construct dmatlunu matrix 
      do inu=0,1
! TR even corresponds to itr=1 and TR odd corresponds to itr=-1
       do i=1,lmmax
          do j=1,lmmax
            trdmat(i,1,j,1,inu)=0.5d0*(ndmatnu(i,j,0,inu)+ndmatnu(i,j,3,inu))
            trdmat(i,1,j,2,inu)=0.5d0*(ndmatnu(i,j,1,inu)-zi*ndmatnu(i,j,2,inu))
            trdmat(i,2,j,1,inu)=0.5d0*(ndmatnu(i,j,1,inu)+zi*ndmatnu(i,j,2,inu))
            trdmat(i,2,j,2,inu)=0.5d0*(ndmatnu(i,j,0,inu)-ndmatnu(i,j,3,inu))
          end do
        end do
      end do
! end loop over atoms and species
deallocate(ndmat,ndmatnu)
return
end subroutine

