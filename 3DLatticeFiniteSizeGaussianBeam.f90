!!!!!!!!!!!!!!!!!!!!!1Calculates finite size scaling of 2D lattice system!!!!!!!!!!!!!!!!!!!!!!!!!!1

real*8 function step(x,V0,rc)
 
real*8  :: V0,rc,x

if (abs(x).le.rc) then

step=V0

else
step=0
end if

end function step

real*8 function distance(r1,r2)
!integer :: i1,i2
real*8, dimension(3),intent(in)  :: r1,r2

distance=Sqrt(1.d0*(r1(1)-r2(1))*(r1(1)-r2(1))+1.d0*(r1(2)-r2(2))*(r1(2)-r2(2))+1.d0*(r1(3)-r2(3))*(r1(3)-r2(3)))


end function distance

real*8 function dist2xy(r1,r2)
!integer :: i1,i2
real*8, dimension(3), intent(in)  :: r1,r2

dist2xy=1.d0*r1(1)*r1(1)+1.d0*r1(2)*r1(2)+1.d0*r2(1)*r2(1)+1.d0*r2(2)*r2(2)


end function dist2xy


real*8 function pot(d1,d2,V0,rc) !in units of hbar: V_{ij}/hbar
real*8, intent(in)  :: d1,d2,V0, rc !d1 dist2xy, d2 distance
real*8, parameter   :: Omega=20*2*acos(-1._8)*10**6 !-> t in seconds
real*8, parameter   :: waist=30.d0/0.406 !in micrometer divided by lattice constant in micrometer

!    interface
!     real*8 function dist2xy(r1,r2)
!     real*8, dimension(3), intent(in)   :: r1,r2
!     end function distxy
! 
!     real*8 function distance(r1,r2)
!     real*8, dimension(3),intent(in)  :: r1,r2
! 
!    end interface 

! d1=dist2xy(r1,r2)
! d2=distance(r1,r2)
pot=1963.5*Exp(-d1/(waist*waist))*rc**6/(d2**6+rc**6)


end function pot



function xi(t, Ntot, V0,rc,r)
integer , intent(in)	:: Ntot
real*8, dimension(2)    :: xi
real*8, dimension(Ntot,3), intent(in):: r
real*8 , intent(in) 	:: t, rc
real*8, intent(in)	:: V0
integer	            	::  i,j,k,l
real*8			::JxJy,Jz, Jx2, angle, eps, epstest,pi=acos(-1._8)
real*8	                :: sumJz, sumJxJy, sumJx2, prodJz, prodJxJy, prodJx21, prodJx22
real*8                  ::  theta
real*8                        ::  testik, testij, testjk,testikxy, testjkxy
real*8                        :: temp , Cutoff 

   interface
    real*8 function distance(r1,r2)
    real*8, dimension(3), intent(in)   :: r1,r2
    end function distance

    real*8 function dist2xy(r1,r2)
    real*8, dimension(3), intent(in)  :: r1,r2
    end function dist2xy

    real*8 function pot(d1,d2,V0,rc)
    real*8,intent(in)  :: d1,d2,V0, rc
    end function pot
   end interface 

sumJz=(0.d0,0.d0)
sumJxJy=(0.d0,0.d0)
sumJx2=(0.d0,0.d0)

Cutoff= 15*rc
do i=1,Ntot

  prodJz=(1.d0, 0.d0)    
    do k=1,Ntot
     if(k==i ) cycle
     testik=distance(r(i,:), r(k,:))
     testikxy=dist2xy(r(i,:), r(k,:))

!      if(testik.le.Cutoff) then
     prodJz=prodJz*Cos(t*0.5d0*pot(testik,testikxy, V0, rc) ) 
!      end if

    end do 

    sumJz=sumJz+prodJz
enddo


Jz=-0.5d0*sumJz


!!!!!!!!!!!!!!ab hier Šndern
do i=1,Ntot

 do j=i+1, Ntot
    testij=distance(r(i,:), r(j,:))
   prodJxJy=(1.d0, 0.d0)
   prodJx21=(1.d0, 0.d0)
   prodJx22=(1.d0, 0.d0)
      do k=1,Ntot
!         if (((r(k,1)==r(i,1)) .and. (r(k,2)==r(i,2))) .OR. ((r(k,1)==r(j,1)) .and. (r(k,2)==r(j,2)))) cycle
          if (k==i) cycle
          if (k==j) cycle

            testik=distance(r(i,:), r(k,:))
            testikxy=dist2xy(r(i,:), r(k,:))
            testjk=distance(r(j,:), r(k,:))
            testjkxy=dist2xy(r(j,:), r(k,:))
!         if((testik .le.Cutoff)) then  
! .and. (testij.le.Cutoff)
        prodJxJy=prodJxJy*Cos(t*0.5d0*pot(testikxy,testik, V0, rc)) 
!         end if

!         if((testik.le.Cutoff).and. (testjk.le.Cutoff)) then  
        prodJx21=prodJx21*Cos(t*0.5d0*( pot(testikxy,testik, V0, rc) -pot(testjkxy,testjk, V0, rc)  )  )
        prodJx22=prodJx22*Cos(t*0.5d0*( pot(testikxy, testik,V0, rc) +pot(testjkxy,testjk, V0, rc)  )  )
!         end if


      end do 
!    if(testij.le.Cutoff) then  
   prodJxJy=prodJxJy*Sin(t*0.5d0*pot(dist2xy(r(i,:), r(j,:)),distance(r(i,:), r(j,:)), V0, rc)  )   
!    end if 


sumJxJy=sumJxJy+prodJxJy
sumJx2=sumJx2+prodJx21-prodJx22

 end do
end do


JxJy=-sumJxJy
Jx2=0.25d0*sumJx2


eps=Ntot*(Ntot*0.25d0+0.5d0*Jx2+0.5d0*JxJy)/(Jz*Jz)
xiPi4=Ntot*(Ntot*0.25d0+0.5d0*Jx2+0.5d0*JxJy)/(Jz*Jz)

temp=0
test=2

do i=1,25000
   theta=i*0.5*pi/25000 
!0.5*pi

   
   epstest=Ntot*(Ntot*0.25d0+Cos(theta)*Cos(theta)*Jx2+Cos(theta)*Sin(theta)*JxJy)/(Jz*Jz)
!    temp=4*(1+Cos(theta)*Cos(theta)*3*0.5d0*(1-Cos(V0*2.d0*t)**2)-6*Cos(theta)*Sin(theta)*Sin(V0*t)*Cos(V0*t)**2)/(4*Cos(V0*t)**6)

   if(epstest< eps) then
   eps=epstest
   angle=theta
   end if

! 
!    if(temp<test) then
!    test=temp
!    end if

end do

xi(1)=eps
xi(2)=angle

end function xi

real*8 function xi_prime(t, Ntot, V0,rc, r)
real*8 , intent(in) 	:: t, rc
real*8, dimension(Ntot,3), intent(in):: r
integer , intent(in)	:: Ntot
real*8, intent(in)	:: V0
real*8, parameter        :: dt=0.0000001
real*8, dimension(2)        :: dummy1, dummy2

  interface
    function xi(t, Ntot, V0,rc, r)
    real*8, dimension(2)    :: xi
    real*8, dimension(Ntot,3), intent(in):: r
    real*8 , intent(in) 	:: t, rc
    real*8, intent(in)	:: V0

    end function xi
  end interface

dummy1=xi(t+dt/2, Ntot,V0,rc,r)
dummy2=xi(t-dt/2, Ntot,V0,rc,r)
xi_prime=(dummy1(1)-dummy2(1))/(V0*dt)



end function xi_prime


real*8 FUNCTION rtbis(x1,x2, Ntot, V0,rc,r)
real*8 , intent(in) 	:: rc, x1,x2
real*8, dimension(Ntot,3), intent(in):: r
integer , intent(in)	:: Ntot
real*8, intent(in)	:: V0
INTEGER  :: j
integer, parameter :: jmax=100000
real*8, parameter  ::xacc=0.00000000001
REAL*8   :: dx,f,fmid,xmid

   interface
    real*8 function xi_prime(t, Ntot, V0,rc,r)
    real*8 , intent(in) 	:: t, rc
    real*8, dimension(Ntot,3), intent(in):: r
    integer , intent(in)	:: Ntot
    real*8, intent(in)	:: V0

    end function xi_prime
   end interface


fmid=xi_prime(x2,Ntot,  V0,rc,r)
f=xi_prime(x1,Ntot, V0,rc,r)

! if(f.lt.0.)then
! rtbis=x1
! dx=x2-x1
! else
! rtbis=x2
! dx=x1-x2
! endif
! do  j=1,JMAX
! dx=dx*0.5d0
! xmid=rtbis+dx
! fmid=xi_prime(xmid,Ntotal,Cutoff, V0,rc)
! if(fmid.le.0.)rtbis=xmid
! if(abs(dx).lt.xacc .or. fmid.eq.0.) return
! enddo 
write(10, '(4(a,e15.8,1X))') "fmid", fmid, "f", f, "x1", x1, "x2", x2
write(*, '(4(a,e15.8,1X))') "fmid", fmid, "f", f, "x1", x1, "x2", x2


if(f.lt.0.)then
rtbis=x1
dx=x2-x1
else
rtbis=x2
dx=x1-x2
endif
do  while (abs(fmid) .ge. 0.001)
dx=dx*0.5d0
xmid=rtbis+dx
fmid=xi_prime(xmid,Ntot, V0,rc,r)
write(10, '(4(a,e15.8,1X))') "fmid", fmid, "xmid", xmid
write(*, '(4(a,e15.8,1X))') "fmid", fmid, "xmid", xmid
if(fmid.le.0.)rtbis=xmid
end do

end function rtbis



program expectation_Values

implicit none


!    interface
!       real*8 function step(x,V0,rc)
!       real*8,intent(in)  :: V0,rc,x
!       end function step
! 
!       real*8 function distance(r1,r2)
!       real*8, dimension(2), intent(in)  :: r1,r2
!       end function distance
! 
! 
!       real*8 function pot(x,V0,rc)
!       real*8, intent(in)  :: x,V0, rc
!       end function pot
! 
!    end interface
interface
    function xi(t, Ntot, V0,rc, r)
    real*8, dimension(2)    :: xi
    real*8, dimension(Ntot,3), intent(in):: r
    real*8 , intent(in) 	:: t, rc
    real*8, intent(in)	:: V0

    end function xi


    real*8 function xi_prime(t, Ntot, V0,rc,r)
    real*8 , intent(in) 	:: t, rc
    real*8, dimension(Ntot,3), intent(in):: r
    integer , intent(in)	:: Ntot
    real*8, intent(in)	:: V0

    end function xi_prime


    real*8 FUNCTION rtbis(x1,x2, Ntot, V0,rc,r)
    real*8 , intent(in) 	:: rc, x1,x2
    real*8, dimension(Ntot,3), intent(in):: r
    integer , intent(in)	:: Ntot
    real*8, intent(in)	:: V0
   
    end function rtbis

end interface

integer                       :: i, j, k, l, time, zahl
integer                       :: n,m, num_threads,status, timestep
real*8                        :: rc,x1,x2
integer                       :: Ntot
real*8                        :: sumJz, sumJxJy, sumJx2, prodJz, prodJxJy, prodJx21, prodJx22
real*8, parameter             :: V0=1000, pi=acos(-1._8)
real*8                        :: dt,t, theta, distance, dist2xy  ,testik, testij, testjk,testjkxy, testikxy
real*8                        :: temp, step, temp1, temp2,  pot, Cutoff
! integer, parameter            :: max_time=30
real*8, dimension(:,:), allocatable     :: r
! integer, dimension(100*100,2)     :: r
! real*8, dimension(10*10,2):: r
real*8, dimension(2)          :: results
real*8                        :: JxJy,Jz, Jx2, angle, eps, epstest,xiPi4,test
character( len = 5 )          :: n_str, m_str, rc_str, timestep_str  ! input parameter in string format
character( len = 60 )         :: file_out  ! output filename
integer                       :: error


call getarg( 1, n_str )   ! reads in number of particles as 1st parameter
read( n_str, '(I5)' ) n   ! converts n parameter to integer
call getarg( 2, rc_str )   ! reads in rc as 4th parameter, it must be a real number of type x,yz !!!
read( rc_str, '(F5.2)' ) rc   ! converts rc parameter to float




Ntot=n*n*n

if ((n.ge.1) .and. (n .lt. 10)) then
   write(file_out,'(a24,a2,i2,a1,i1,a9)') '3DFiniteSizeGaussianBeam','rc',int(rc*10),'N',n,'hoch3.dat'
else if((n.ge.10) .and. (n.lt. 100))  then
   write(file_out,'(a24,a2,i2,a1,i2,a9)') '3DFiniteSizeGaussianBeam','rc',int(rc*10),'N',n,'hoch3.dat'
else if((n.ge.100) .and. (n.lt. 1000))  then
   write(file_out,'(a24,a2,i2,a1,i3,a9)') '3DFiniteSizeGaussianBeam','rc',int(rc*10),'N',n,'hoch3.dat'
else
   write(file_out,'(a24,a2,i2,a1,i4,a9)') '3DFiniteSizeGaussianBeam','rc',int(rc*10),'N',n,'hoch3.dat'
endif

write( *, '(3a, F5.2)' ) "Read in parameters: filename = ", file_out, ";  rc = ", rc
write(*,*)  "Ntot= ", n*n*n



! num_threads = 4
! call OMP_SET_NUM_THREADS(num_threads)




!if ( error .ne. 0 ) then
!   write( *, * ) "Fatal error in allocating arrays"
!   stop
!endif

open( unit = 10, FILE = file_out, ACTION = "write", STATUS = "replace" )   ! opens output file for writing

allocate( r(-(Ntot-1)/2:(Ntot-1)/2,3), stat = error )

 Cutoff=15*rc

! Save positions in the lattice

r=0.d0
do i=-(n-1)/2,(n-1)/2
  do j=-(n-1)/2,(n-1)/2
    do k=-(n-1)/2,(n-1)/2

         zahl=i+n*j+n*n*k  

     
         r(zahl,1)=i*1.d0
         r(zahl,2)=j*1.d0
         r(zahl,3)=k*1.d0


! write(10,'(I3,4X,2(F14.8,1X))') k, r(k,:)

   end do
  end do
end do

do i=-(Ntot-1)/2,(Ntot-1)/2

write(*,*) i, r(i,:)

end do

! write(10,'(a5,4X,F14.8,1X)') "erste", r(:,1)
! write(10,'(a6,4X,F14.8,1X)') "zweite", r(:,2)
dt=0.000001


!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(Ntot,r,rc,Jx2,JxJy, Jz,dt,angle,eps,epstest,xiPi4) PRIVATE(sumJx2, sumJz, sumJxJy,prodJx21, prodJx22, prodJz,prodJxJy,t,time,theta)
! do time=1, max_time
!150+70,150+90 for rc=0.5
! t=0.0002152+timestep*dt
!fuer rc=2
! t=0.000068+timestep*dt 
!fuer n=m=50

! 
x1=0.0000001
x2=0.0000259

t=rtbis(x1,x2, Ntot,  V0,rc,r) 
results=xi(t, Ntot, V0,rc,r)

  
write(10,'(3(F14.8,1X),3X,I5,3X,F14.8,3X,F14.8,3X,F14.8,3X,F14.8,3X,F14.8,3X,F14.8 )') t*1963.5, results, Ntot, rc, 1.d0/(3**(1.d0/3.d0)*(0.5d0*Ntot)**(2.d0/3.d0) )






deallocate( r )
close( 10 )

end program expectation_Values
