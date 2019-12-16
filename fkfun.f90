subroutine fkfun(x,f,ier2)

use system
use chainsdat
use molecules
use const
use results
use bulk
use kai
use MPI
use fields_fkfun
use kinsol
use conformations
use ematrix
use ellipsoid
use transform
use kaist
implicit none

integer*4 ier2
integer ntot
real*8 x(*),f(*)
real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer jx, jy, jz, jj
real*8 xpot(dimx, dimy, dimz)
! Charge
real*8 psitemp
real*8 MV(3),MU(3),MW(3)
real*8 MVV,MUU,MWW,MVU,MVW,MUW
real*8 psivv,psiuu,psiww, psivu,psivw,psiuw
real*8 psiv(3), epsv(3)

integer, external :: PBCSYMI, PBCREFI

! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err
real*8 avpol_tosend(dimx,dimy,dimz)
real*8 avpol_temp(dimx,dimy,dimz)
real*8 q_tosend, sumgauche_tosend
real*8 gradpsi2
real*8 fv
real*8 histoe2e_tosend(nhist)
real*8 histoe2e_temp(nhist)
!-----------------------------------------------------
! Common variables

shift = 1.0


! Jefe

if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, 2*dimx*dimy*dimz , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

!------------------------------------------------------
! DEBUG
!      if(iter.gt.2000) then
!      do i = 1, n
!      print*,i, x(i)
!      enddo
!      endif


! Recupera xh y psi desde x()

ntot = dimx*dimy*dimz ! numero de celdas
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xh(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
     psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)   
  enddo
 enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Boundary conditions electrostatic potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reflection or PBC, (PBC = 1 or 3)
 
do jx = 0, dimx+1
do jy = 0, dimy+1
do jz = 0, dimz+1

if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)


if ((PBC(1).eq.3).and.(ix.lt.1))ix = PBCREFI(jx,dimx)
if ((PBC(2).eq.3).and.(ix.gt.dimx))ix = PBCREFI(jx,dimx)

if ((PBC(3).eq.3).and.(iy.lt.1))iy = PBCREFI(jy,dimy)
if ((PBC(4).eq.3).and.(iy.gt.dimy))iy = PBCREFI(jy,dimy)

if ((PBC(5).eq.3).and.(iz.lt.1))iz = PBCREFI(jz,dimz)
if ((PBC(6).eq.3).and.(iz.gt.dimz))iz = PBCREFI(jz,dimz)


   psi(jx, jy, jz) = psi(ix, iy, iz)

enddo
enddo
enddo

! Bulk or Wall, PBC = 0 or 2

select case (PBC(1)) ! x = 0
case(0) ! set bulk 
   psi(0,:,:) = 0.0 
case(2)
   psi(0,:,:) = psi(1,:,:) ! zero charge
endselect

select case (PBC(2)) ! x = dimx
case(0) ! set bulk 
   psi(dimx+1,:,:) = 0.0  
case(2)
   psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
endselect

select case (PBC(3)) ! y = 0
case(0) ! set bulk 
   psi(:,0,:) = 0.0  
case(2)
   psi(:,0,:) = psi(:,1,:) ! zero charge
endselect

select case (PBC(4)) ! y = dimy
case(0) ! set bulk 
   psi(:,dimy+1,:) = 0.0
case(2)
   psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
endselect

select case (PBC(5)) ! z = 0
case(0) ! set bulk 
   psi(:,:,0) = 0.0  
case(2)
   psi(:,:,0) = psi(:,:,1) ! zero charge
endselect

select case (PBC(6)) ! z = dimz
case(0) ! set bulk 
   psi(:,:,dimz+1) = 0.0
case(2)
   psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
endselect

! volume fraction and frdir

histoe2e = 0.0

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
    avpol(ix,iy,iz)=0.0
    xpos(ix, iy, iz) = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction 
    xneg(ix, iy, iz) = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
    xOHmin(ix, iy,iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction
    fdis(ix,iy,iz)=1.0 /(1.0 + xHplus(ix,iy,iz)/(K0*xh(ix,iy,iz)) )
   enddo
 enddo  
enddo


! Calculo de xtotal para poor solvent
! en el lattice
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
   xtotal(ix, iy, iz) = 1.0-xpos(ix, iy,iz) - xneg(ix, iy, iz)- xh(ix, iy, iz) - xHplus(ix, iy, iz) - xOHmin(ix, iy, iz) ! xtotal es todo menos solvente e iones
  enddo
 enddo
enddo


! Compute dielectric permitivity
call dielectfcn(xtotal,volprot,epsfcn,Depsfcn)

!------------------------------------------------------------------------
! PDFs polimero
!------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PARALELO: Cada procesador trabaja sobre una cadena...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcula xpot

sttemp = st/(vpol*vsol)

do ix=1,dimx
 do iy=1,dimy
   do iz=1,dimz

     fv = (1.0 - volprot(ix,iy,iz))

     xpot(ix, iy, iz) = xh(ix,iy,iz)**vpol/fdis(ix,iy,iz)*dexp(-psi(ix, iy, iz)*zpol)
     xpot(ix, iy, iz) = xpot(ix,iy,iz)*dexp(voleps(ix,iy,iz))
     
     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 
!     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
!     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

     xpot(ix,iy,iz) = xpot(ix,iy,iz)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol/fv)

     protemp=0.0

     do ax = -Xulimit,Xulimit 
      do ay = -Xulimit,Xulimit
       do az = -Xulimit,Xulimit

            jx = ix+ax
            jy = iy+ay
            jz = iz+az

            if(jx.lt.1) then
            if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jx.gt.dimx) then
            if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jy.lt.1) then
            if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jy.gt.dimy) then
            if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
            endif


            if(jz.lt.1) then
            if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if(jz.gt.dimz) then
            if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
            endif


            if((jx.ge.1).and.(jx.le.dimx)) then
            if((jy.ge.1).and.(jy.le.dimy)) then
            if((jz.ge.1).and.(jz.le.dimz)) then
                fv = (1.0-volprot(jx,jy,jz))
                protemp=protemp + Xu(ax,ay,az)*sttemp*xtotal(jx, jy, jz)*fv
            endif
            endif
            endif

       enddo
      enddo
     enddo

     xpot(ix,iy,iz) = xpot(ix,iy,iz)*dexp(protemp)

   enddo
  enddo
enddo

avpol_tosend = 0.0
q = 0.0
sumgauche = 0.0
histoe2e_tosend = 0.0

do jj = 1, cpp(rank+1)
   ii = cppini(rank+1)+jj

   histoe2e_temp = 0.0

   q_tosend=0.0
   sumgauche_tosend = 0.0
   avpol_temp = 0.0

 do i=1,newcuantas(ii)
   pro(i, jj)=shift
   do j=1,long
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    pro(i, jj) = pro(i, jj) * xpot(ax, ay, az)
   enddo
    pro(i,jj) = pro(i,jj)*exp(-benergy*ngauche(i,ii)) ! energy of gauche bonds
   do j=1,long
   fv = (1.0-volprot(px(i,j, jj),py(i,j, jj),pz(i,j, jj)))
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj))= &
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj))+pro(i, jj)*vpol*vsol/(delta**3)/fv* &
    ngpol(ii)*sc ! ngpol(ii) has the number of chains grafted to the point ii
   enddo

   histoe2e_temp(e2e(i,ii)) = histoe2e_temp(e2e(i,ii)) + pro(i, jj)*ngpol(ii)


   q_tosend=q_tosend+pro(i, jj)
   sumgauche_tosend = sumgauche_tosend+ngauche(i, ii)*pro(i,jj)

 enddo ! i
! norma 
 do ix=1,dimx
  do iy=1,dimy
   do iz=1,dimz
    avpol_tosend(ix,iy,iz)=avpol_tosend(ix, iy, iz) + avpol_temp(ix,iy,iz)/q_tosend
    enddo
   enddo
 enddo
q(ii) = q_tosend ! no la envia ahora
sumgauche(ii) = sumgauche_tosend/q_tosend
histoe2e_tosend = histoe2e_tosend + histoe2e_temp/q_tosend

!print*, rank+1,jj,ii,q(ii)
enddo ! jj

!if(rank.eq.0) then
!do i = 1, 100
!print*, i, pro(i,1), pro(i,100), e2e(i,1), e2e(i,100)
!enddo
!stop
!endif


!------------------ MPI ----------------------------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(histoe2e_tosend, histoe2e, nhist, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err) 
  call MPI_REDUCE(histoe2e_tosend, histoe2e, nhist, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

!!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot

do ix=1,dimx
   do iy=1,dimy
        do iz=1,dimz
  
         fv = (1.0-volprot(ix,iy,iz))

         qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + avpol(ix, iy, iz)*zpol/vpol*fdis(ix,iy,iz) + &
         xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

         qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol    ! OJO

        enddo
   enddo
enddo

! Volume fraction

do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= avpol(ix,iy,iz) + xh(ix,iy,iz) + &
      xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
      xOHmin(ix, iy, iz) -1.000000d0
      enddo
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Poisson equatio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
!

MV(1) = MAT(1,1)
MV(2) = MAT(1,2)  
MV(3) = MAT(1,3)

MU(1) = MAT(2,1)
MU(2) = MAT(2,2)  
MU(3) = MAT(2,3)

MW(1) = MAT(3,1)
MW(2) = MAT(3,2)  
MW(3) = MAT(3,3)

MVV = DOT_PRODUCT(MV,MV)
MUU = DOT_PRODUCT(MU,MU)
MWW = DOT_PRODUCT(MW,MW)

MVU = DOT_PRODUCT(MV,MU)
MVW = DOT_PRODUCT(MV,MW)
MUW = DOT_PRODUCT(MU,MW)

do ix=1,dimx
do iy=1,dimy
do iz=1,dimz

psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)

psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0
psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0
psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0

psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0
psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0
psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0

epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0
epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0
epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0

psitemp = epsfcn(ix,iy,iz)*(MVV*psivv+MUU*psiuu+MWW*psiww+2.0*MVU*psivu+2.0*MVW*psivw+2.0*MUW*psiuw)
psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

! OJO CHECK!!!!

      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)=(psitemp + qtot(ix, iy, iz)*constq)/(-2.0)
      if(electroflag.eq.0)f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)=0.0

enddo
enddo
enddo
 
norma = 0.0

do i = 1, 2*ntot
  norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
if(rank.eq.0)print*,'fkfun:', iter, norma, q(1)
endif

3333 continue
ier2 = 0.0 

return
end
