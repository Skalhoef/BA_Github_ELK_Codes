
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energytd
use modmain
use modtddft
implicit none
! local variables
integer is,ias
real(8) ca,engya,t1,t2
! external functions
real(8), external :: rfinp
! Coulomb potential energy
engyvcl=rfinp(rhomt,rhoir,vclmt,vclir)
! Madelung term
engymad=0.d0
do ias=1,natmtot
  is=idxis(ias)
  engymad=engymad+0.5d0*spzn(is)*(vclmt(1,ias)-vcln(1,is))*y00
end do
! exchange and correlation energy
engyx=rfinp(rhomt,rhoir,exmt,exir)
engyc=rfinp(rhomt,rhoir,ecmt,ecir)
! vector potential contributions to energy
ca=-1.d0/solsc
engya=ca*dot_product(afieldt(:,itimes),jtot(:))
! constant term (1/2c²)A(t)²Q
ca=0.5d0/solsc**2
t1=sum(afieldt(:,itimes)**2)
t2=sum(chgstot(:))/3.d0
engya=engya+ca*t1*(chgtot-t2)
! electric field contribution (1/2)E(t)²
t1=sum(efieldt(:)**2)
engya=engya+0.5d0*t1
! mass term
t2=ca*(afindpm(0)/afindpm(2))
engya=engya+t1*t2
! total energy
engytot=engykn+0.5d0*engyvcl+engymad+engyx+engyc+engya
end subroutine

