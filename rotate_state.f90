!Tue 08 Jan 2019 11:06:05 AM MST ! data coming from poppe et.al 2014 
! 13 March 2013
! rmoon_real: radius of the Moon
! dn_sw_real: density in the solar wind (1/m^3)
! mi         : mass of proton (kg)
! e          : electron charge
! gama       : adiabatic index
! mu0        : electromagnetism, permeability
! cs         : sound speed in plasma (Km/s) ! or cs=sqrt((temp_e+((gama)*temp_i))*(e/mi))/1000
!   shared variables
      module constants
        double precision, parameter :: pi=3.14159265358979323846264338328, &
                                       rmoon_real=1737.d0, mi=167.d-29,e=16.d-20, &
                                       gama=5.d0/3.d0,mu0=pi*4.d-7, dn_sw_real=7.d6!, &
                                      ! cs=37.d0
       
        double precision :: rmoon,xmoon,ymoon,dninfty,Bx_sw,By_sw,Bz_sw, &
                            temp_e,temp_i,vx_sw,vy_sw,vz_sw,xmin,xmax,ymin,ymax,z,time

        namelist/charinp/rmoon,xmoon,ymoon,dninfty,Bx_sw,By_sw,Bz_sw, &
        temp_e,temp_i,vx_sw,vy_sw,vz_sw,xmin,xmax,ymin,ymax
      end module constants
!=======================================================================      

! alpha   : angel between B and V_sw
! minfty  : mach number in the solar wind (infinity)
! angel1  :
! angel2  : 

      program rotate
      use constants
      implicit none

!   variables
      integer t,nlines
      !Real(8), Dimension(120) :: x,y,z,cap_radius,norm_distance
      double precision :: alpha,minfty,angle1,angle2,cs,v_sw,B_sw,x,y,zz, &
                          cap_radius,norm_distance,angle3,UT
      double precision :: Ry(3,3),Rz(3,3),R(3,3),pos(3,1),pos_rotate(3,1),B_sw_rotate(3,1), &
                          v_sw_rotate(3,1),B_sw_ART(3,1),v_sw_ART(3,1),Rx(3,3),pos_rotate2(3,1), &
                          pos_new(3,1)

!   Read variables
      open(unit=4,file='char.dat',status='old')      
      read(4,nml=charinp)
     
! z : a constant comes from my model (characteristic model)
      z=(temp_e/(temp_e+temp_i))


! In the paper suggested the temp_i is between (4-6(eV)). I chose 5(eV) for now because I ran FEC with this.
      cs=sqrt((temp_e+((gama)*temp_i))*(e/mi))/1000.d0 ! km/s
      v_sw=sqrt((vx_sw)*(vx_sw) + (vy_sw)*(vy_sw) + (vz_sw)*(vz_sw))
      B_sw=sqrt((bx_sw)*(bx_sw) + (by_sw)*(by_sw) + (bz_sw)*(bz_sw))
      minfty= (v_sw/cs)
      alpha= (180.d0/pi) * acos((vx_sw*bx_sw + vy_sw*by_sw + vz_sw*bz_sw)/(v_sw*B_sw)) ! degree

! these values came from ARTEMIS website, they are in SSE coordinate system which the origin in the Moon.
! My coordinate system's origin is the Sun. So first I have to rotate the coordinate system about y with an angel 180
! Then in the new coordinate x is from Sun toward Moon.
! this rotation means all x componants should be tmes -1.
      B_sw_ART(1,1)=-Bx_sw
      B_sw_ART(2,1)=By_sw
      B_sw_ART(3,1)=Bz_sw
      v_sw_ART(1,1)=-vx_sw
      v_sw_ART(2,1)=vy_sw
      v_sw_ART(3,1)=vz_sw


! ---------------------------------------------------------------------------------------------------------------      
! ***************************************rotation****************************************************************
! ---------------------------------------------------------------------------------------------------------------      

! rotate xyz about y clockwise by angele angle1 (y does not change):  
     ! angle1=atan(vz_sw/vx_sw) ! radian because cos and sin and... accept radian

!   Matrix (Ry):
           !Ry11 Ry12 Ry13 
           !Ry21 Ry22 Ry23 
           !Ry31 Ry32 Ry33 

        !Ry(1,1)= cos(angle1)
        !Ry(2,1)= 0.d0
        !Ry(3,1)= -sin(angle1)

        !Ry(1,2)= 0.d0
        !Ry(2,2)= 1.d0
        !Ry(3,2)= 0.d0
        
        !Ry(1,3)= sin(angle1)
        !Ry(2,3)= 0.d0
        !Ry(3,3)= cos(angle1)


! Now x^ y^ z^ are in a plane which I have to rotate
! rotate x^ y^ z^ about z^ counterclickwise by angele angle2  be in V_sw direction: 

!      angle2=atan(vy_sw/vx_sw) !radian

!   Matrix (Rz):(this is inverse of original Rz)
           !Rz11 Rz12 Rz13 
           !Rz21 Rz22 Rz23 
           !Rz31 Rz32 Rz33 

        Rz(1,1)= v_sw_ART(1,1)/v_sw!vx_sw/v_sw!cos(angle2)
        Rz(2,1)= -v_sw_ART(2,1)/v_sw!-vy_sw/v_sw!-sin(angle2)
        Rz(3,1)= 0.d0

        Rz(1,2)= v_sw_ART(2,1)/v_sw!vy_sw/v_sw!sin(angle2)
        Rz(2,2)= v_sw_ART(1,1)/v_sw!vx_sw/v_sw!cos(angle2)
        Rz(3,2)= 0.d0
        
        Rz(1,3)= 0.d0
        Rz(2,3)= 0.d0
        Rz(3,3)= 1.d0 

        B_sw_rotate=matmul(Rz,B_sw_ART)
        v_sw_rotate=matmul(Rz,v_sw_ART)

        !angle3=atan(B_sw_rotate(2,1)/B_sw_rotate(3,1))
        !angle1=asin(B_sw_rotate(2,1)/sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))

        Rx(1,1)= 1.d0
        Rx(2,1)= 0.d0
        Rx(3,1)= 0.d0

        Rx(1,2)= 0.d0
        Rx(2,2)= B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))!cos(angle3)
        Rx(3,2)= B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))!sin(angle3)
        
        Rx(1,3)= 0.d0
        Rx(2,3)= -B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))!sin(angle3)
        Rx(3,3)= B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))!cos(angle3)     
        

       !r(1,1)=v_sw_ART(1,1)/v_sw!vx_sw/v_sw
       !r(2,1)=(-v_sw_ART(2,1)/v_sw)*(B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))) &
              !(-vy_sw/v_sw)*(B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1))))
       !r(3,1)=(-v_sw_ART(2,1)/v_sw)*(B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1))))
              !(-vy_sw/v_sw)*(B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1))))

      ! r(1,2)=(v_sw_ART(2,1)/v_sw)!vy_sw/v_sw
      ! r(2,2)=v_sw_ART(1,1)/v_sw!vx_sw/v_sw*((B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))))
      ! r(3,2)=(v_sw_ART(1,1)/v_sw)*((B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))))
              !vx_sw/v_sw*((B_sw_rotate(2,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))))

      ! r(1,3)=0.d0
      ! r(2,3)=-((B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))))
      ! r(3,3)=+((B_sw_rotate(3,1)/(sqrt(B_sw_rotate(2,1)*B_sw_rotate(2,1)+B_sw_rotate(3,1)*B_sw_rotate(3,1)))))



! vector a(ax,ay,az) in the new coordinate is a(ax_p,ay_paz_p)=R* a(ax,ay,az)
        
           
!  1.  Read variables


      nlines = 0 
      open(unit=5,file='THB_L1_STATE_27227_finalEdited.txt1')
      DO
      READ(5,*,END=10)
      nlines = nlines + 1
      END DO
10      CLOSE (5)

      !print*, nlines
      
! This code read data (UT,X,Y,Z) of ARTEMIS. Then Calculate the CAP radious and normalized the positions to the cap radious to
! make out 3D data to 2D. after that  roatet them to be in a plane of V_sw and B.

! Rotate:
      open(unit=5,file='THB_L1_STATE_27227_finalEdited.txt1')
      open(unit=55,file='ART_Normalized_Poss_in_V-B_Plane.out')
     
      do t=0,nlines-1
    ! I introduce this "time" to have an idea about the real time of each accepted coordinates(easier for plotting)
       time=t*1.0
         READ(5,*)Ut,x,y,zz

        

        if (y<rmoon_real .and. y>-rmoon_real) then
            cap_radius = sqrt(rmoon_real*rmoon_real - y*y)
	    pos(1,1) = -x/cap_radius 
            pos(2,1) =  zz/cap_radius
            pos(3,1) =  0.d0

            pos_rotate=matmul(Rz,pos)
            pos_new=matmul(Rx,pos_rotate)
            x=pos_new(1,1)
            y=pos_new(2,1)
           ! print*,UT,x,y
            write(55,102)time,UT,x,y


         else
            cap_radius = 0.0
           Print*, t,'no intersection '
         endif
         !print*,t,pos
      enddo


! Calculate cap radius:       

!      open(unit=5,file='THB_L1_STATE_27227_finalEdited.txt1')
!      rmoon = 1737.0
!      do t=1,nlines
!         READ(5,*)x(t),y(t),z(t)
!         if( y(t) < rmoon .and. y(t) > (-rmoon) ) then
!      	   cap_radius(t) = sqrt(rmoon*rmoon - y(t)*y(t))
!     	   norm_distance(t) = (sqrt( x(t)*x(t)+y(t)*y(t)+z(t)*z(t) )) / (cap_radius(t))
!         else
!           cap_radius(t) = 0.0
!          !norm_distance(t) = 100000.0
!           Print*, 'no intersection '
!         endif
!      enddo
!      open(unit=7,file='normalized.out',status='unknown')
!      do t=1,nlines
!         write (7,102)cap_radius(t),norm_distance(t)
!      enddo

 102  format(99es15.5)
 103  format(a,99es15.5)
      stop
      end
