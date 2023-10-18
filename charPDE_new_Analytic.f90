!   By using the book "techniques in partial diffrential equations", 
!   I want to solve the four nonlinear partial diffrential equations.
!   the variables are: U1=ln(n1)=ln1, U2=M||1=mp1, U3=ln(n2)=ln2, U4=M||2=mp2
!   minfty is a constant and is a perpendicular mach number
!   matrix "A" which is the coeeficent of d/dx
!   matrix "B" which is the coeeficent of d/dy 
!   temp_e is electron temperatur and temp_i is ion temperatur
!   in "dgeev"routin:
!       WR and WI   contain the real and imaginary parts,
!                   respectively, of the computed eigenvalues. 
!       VR          (output) DOUBLE PRECISION array, dimension (LDVR,N)
!           If JOBVR = 'V', the right eigenvectors v(j) are stored one
!                   after another in the columns of VR, in the same order
!                   as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!                  the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!                 conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!                 v(j+1) = VR(:,j) - i*VR(:,j+1).
!   var_new(n): updated unknowns, var_new(1)=updated ln1,
!                                 var_new(2)=updated mp1,
!                                 var_new(3)=updated ln2,
!                                 var_new(4)=updated mp2.
!=======================================================================

!   Variable related to comparisons with ARTEMIS
      module ARTEMIS
        double precision, parameter :: pi=3.14159265358979323846264338328, &
                                       rmoon_real=1737.d0, mi=167.d-29,e=16.d-20, &
                                       gama=5.d0/3.d0,mu0=pi*4.d-7, dn_sw_real=7.d6!, &
                                      ! cs=37.d0

        integer:: t,ixb,iyb,nlines
        
        double precision :: time,cap_radius,ww1,ww2,xb,yb,pos_x,pos_y,pos_z, &
                            dn1_inter,dn2_inter,mp1_inter,mp2_inter,xx,yy,zz,ut
                            
        double precision :: Ry(3,3),Rz(3,3),R(3,3),pos(3,1),pos_rotate(3,1),B_sw_rotate(3,1), &
                            v_sw_rotate(3,1),B_sw_ART(3,1),v_sw_ART(3,1),Rx(3,3),pos_rotate2(3,1), &
                            pos_new(3,1)
       
        double precision :: rmoon,xmoon,ymoon,dninfty,Bx_sw,By_sw,Bz_sw, &
                            temp_e,temp_i,vx_sw,vy_sw,vz_sw,xmin,xmax,ymin,ymax,RHS, &
                             alpha,cs,v_sw,B_sw

        namelist/charinp/rmoon,xmoon,ymoon,dninfty,Bx_sw,By_sw,Bz_sw, &
        temp_e,temp_i,vx_sw,vy_sw,vz_sw,xmin,xmax,ymin,ymax,RHS
      end module ARTEMIS


!   Variable related to Characteristic method
      module char_method
      
        integer, parameter :: ny=201,nx=201
      
        double precision :: dx,dy,theta_pos,theta_neg,x_sqrt,y_sqrt,b_pos,b_neg, &
                            yy_pos,yy_neg,x_pos,x_neg,y_test,dn_min,mp1_max, &
                            mp2_min,m_perp,epslon,z
        double precision :: mpinfty,minfty,mp1(0:nx,0:ny),mp2(0:nx,0:ny), &
                            ln1(0:nx,0:ny),ln2(0:nx,0:ny),mp1_p(0:nx,0:ny), &
                            mp2_p(0:nx,0:ny),ln1_p(0:nx,0:ny),ln2_p(0:nx,0:ny), &
                            mp1_anal(0:nx,0:ny),mp2_anal(0:nx,0:ny)
      
      end module char_method
      
      
!   Variable
      module mod1
        integer, parameter :: n=4, LDA=4,itemax=1
        integer:: nr,ir
        double precision :: A(n,n),B(n,n),miu(n),aa(LDA,n),initi_gues, &
                            root_alpha(n),small=1.d-10
      end module
!=======================================================================      

      program char_cell
      use char_method
      use ARTEMIS
      use mod1
      implicit none

!   variables
      CHARACTER :: JOBVL, JOBVR
      integer ix,iy,ix_new,iy_new,i,j,ii,jj,info,j_ite,mini_i,it,iyy
      integer, parameter :: LDVL=1, LDVR=N, Lwork=4*N,nrhs=1
      double precision, allocatable :: tdn(:,:),tdn1(:,:),tdn2(:,:),tdnmp1(:,:), &
                                       tdnmp2(:,:),tdnmp(:,:),tdnm(:,:), &
                                       tdnm1(:,:),tdnm2(:,:)
      double precision :: x,y,x_new,y_new,up0(2),um0(2),up(2),um(2),u0(2),u1(2),u2(2),mpara, &
                          mpara1,mpara2,cm,cp,dn,dn1,dn2,theta1,theta2,c1,c2, & 
                          csalpha,snalpha,cp2,c1m,omega,mp1_min,VL(LDVL,N), VR(LDVR,N), work(Lwork), &
                          WI(n), WR(n),ln1_ite,ln2_ite,mp1_ite,mp2_ite,ty(0:ny),temp2, &
                          beta_coeff(n,n),beta_coeff_t(n,n),B_t(n,n),var_new(n),min_x,max_x,min_dx,ddx, &
                          guess(1:n,0:ny+1)

!   externals
      double precision f3,aazero
      external f3,aazero,ELGS,dgeev,gaussElim

!   Read variables
      open(unit=44,file='char.dat',status='old')      
      read(44,nml=charinp)

!   constant parameters
!     mpinfty: 
!     is parallel velocity which is not constant. for the general case (alpha is not 90)
!     can be defined in the solar wind: mpinfty=vsw*cos(alpha)
!     Minfty: 
!     is the normalized velocity in the solar wind initially NOT just perpendicular velocity.
!     perpendicular velocity is m_perp=vsw*sin(alpha) which is constant in this model


      dx = (xmax-xmin)/(nx) 
      dy = (ymax-ymin)/(ny)     
       
      ! In the paper suggested the temp_i is between (4-6(eV)). I chose 5(eV) for now because I ran FEC with this.
      cs=sqrt((temp_e+((gama)*temp_i))*(e/mi))/1000.d0 ! km/s
      v_sw=sqrt((vx_sw)*(vx_sw) + (vy_sw)*(vy_sw) + (vz_sw)*(vz_sw))
      B_sw=sqrt((bx_sw)*(bx_sw) + (by_sw)*(by_sw) + (bz_sw)*(bz_sw))
      minfty= (v_sw/cs)
      
      alpha= acos((vx_sw*bx_sw + vy_sw*by_sw + vz_sw*bz_sw)/(v_sw*B_sw)) ! rediant
      !(180.d0/pi) * alpha = alpha in degree
      csalpha=cos(alpha)
      snalpha=sin(alpha)
      mpinfty = Minfty*csalpha 
      m_perp  = Minfty*snalpha
      dn_min=1.d-15
      mp1_max=60.d0
      mp2_min=(-1-minfty*tan(alpha)*snalpha) 
      epslon=0.5


      
      !print*,z,temp_e,temp_i,temp_e/(temp_e+temp_i)
      !stop

     ! print*,'_________________________________'
     ! print*,'angle between flow and B:'
     ! print*,'alpha:',alpha*180./pi
     ! print*,'_________________________________'
     ! print*, 'Right Hand Side:'
     ! print*, 'RHS=',RHS

!   allocating the vectors
      allocate(tdn1(0:nx,0:ny))
      allocate(tdn(0:nx,0:ny))
      allocate(tdn2(0:nx,0:ny))

      allocate(tdnmp1(0:nx,0:ny))
      allocate(tdnmp2(0:nx,0:ny))
      allocate(tdnmp(0:nx,0:ny))

      allocate(tdnm1(0:nx,0:ny))
      allocate(tdnm2(0:nx,0:ny))
      allocate(tdnm(0:nx,0:ny)) 

!   open output files
      !open(unit = 15 , file = 'Function_Miu.out',status='unknown')
      !open(unit = 16 , file = 'roots',status='unknown')
      !open(unit = 150 , file = 'test_boundary.out',status='unknown')
      !open(unit = 17 , file = 'compare.out',status='unknown')
      !open(unit = 3 , file = 'testLog_ix0_anal.out',status='unknown')
      !open(unit = 4 , file = 'testLog_ix1_anal.out',status='unknown')

      open(unit = 66 , file = 'new.out',status='unknown')
    !  open(unit = 155 , file = 'analytic.out',status='unknown')

     ! open(unit = 2 , file = 'test_2rmoom_new.out',status='unknown')
     ! open(unit = 3 , file = 'test_3rmoom_new.out',status='unknown')
     ! open(unit = 4 , file = 'test_4rmoom_new.out',status='unknown')
     ! open(unit = 5 , file = 'test_5rmoom_new.out',status='unknown')
     ! open(unit = 16 , file = 'test_6rmoom_new.out',status='unknown')
     ! open(unit = 7 , file = 'test_7rmoom_new.out',status='unknown')
     ! open(unit = 8 , file = 'test_8rmoom_new.out',status='unknown')
     ! open(unit = 9 , file = 'test_9rmoom_new.out',status='unknown')

     ! open(unit = 21 , file = 'test_2rmoom_anal.out',status='unknown')
     ! open(unit = 31 , file = 'test_3rmoom_anal.out',status='unknown')
     ! open(unit = 41 , file = 'test_4rmoom_anal.out',status='unknown')
     ! open(unit = 51 , file = 'test_5rmoom_anal.out',status='unknown')
     ! open(unit = 61 , file = 'test_6rmoom_anal.out',status='unknown')
     ! open(unit = 71 , file = 'test_7rmoom_anal.out',status='unknown')
     ! open(unit = 81 , file = 'test_8rmoom_anal.out',status='unknown')
     ! open(unit = 91 , file = 'test_9rmoom_anal.out',status='unknown')


  


   

!   initialization
      tdn(:,:)=0.d0
      tdn1(:,:)=0.d0
      tdn2(:,:)=0.d0

      tdnmp(:,:)=0.d0
      tdnmp1(:,:)=0.d0  
      tdnmp2(:,:)=0.d0 

      tdnm(:,:)=0.d0
      tdnm1(:,:)=0.d0
      tdnm2(:,:)=0.d0 

      mp1(:,:)=0.d0
      ln1(:,:)=0.d0
      mp2(:,:)=0.d0
      ln2(:,:)=0.d0
      mp1_p(:,:)=0.d0
      ln1_p(:,:)=0.d0
      mp2_p(:,:)=0.d0
      ln2_p(:,:)=0.d0
      mp1_anal(:,:)=0.d0
      mp2_anal(:,:)=0.d0
      
      a(:,:)=0.d0
      b(:,:)=0.d0
      aa(:,:)=0.d0
      miu(:)=0.d0
      beta_coeff(:,:)=0.d0
      beta_coeff_t(:,:)=0.d0
      B_t(:,:)=0.d0
      var_new(:)=0.d0
      ty(:)=0.d0
      VR(:,:)=0.d0
      VL(:,:)=0.d0
      work(:)=0.d0
      WI(:)=0.d0
      WR(:)=0.d0
      guess(:,:)=0.d0

!   Analytic Solution


!   Scann over all mesh points 
        
      do iy=0,ny
        y=ymin+(ymax-ymin)*iy/(ny)
        do ix=0,nx
          x=xmin+(xmax-xmin)*ix/(nx)      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       theta_pos=atan(snalpha/(Minfty+csalpha))
       x_pos=(1.d0/sin(theta_pos))
       b_pos=1.d0/cos(theta_pos)
       yy_pos=tan(theta_pos)*x+b_pos


       theta_neg=atan(-snalpha/(Minfty-csalpha))
       x_neg=(1.d0/sin(theta_neg))
       b_neg=1.d0/cos(theta_neg)
       yy_neg=tan(theta_neg)*x-b_neg

      
       x_sqrt=x*sqrt((x*x)+(y*y)-1.d0)
       y_sqrt=y*sqrt((x*x)+(y*y)-1.d0)
       y_test=tan(alpha)*x-(1.0/csalpha)!  This line is the line that for some angles (check first),
                                 !  above this line 
                                 !  there is no solution for Mp1 (along the magnetic field)


     if (y>=yy_neg) then
      
        mpara1=+1.d0-Minfty*snalpha*((csalpha*(-y+x_sqrt)+snalpha*(+x+y_sqrt)) & 
                                            /(csalpha*(+x+y_sqrt)+snalpha*(y-x_sqrt)))
                                            
!   set a maximum for parallel velocity1???????????????????????????????????????
       ! if(mpara1>=60.d0) then
       !   mpara1=60.d0
       ! endif
!   set a maximum for parallel velocity1???????????????????????????????????????
        mp1(ix,iy)=mpara1
        dn1=dninfty*exp(-abs(mpara1-Minfty*csalpha))
        tdn1(ix,iy)=dn1
 
     else
         mpara1=Minfty*csalpha
         mp1(ix,iy)= mpara1
         dn1=dninfty
         tdn1(ix,iy)=dn1
     endif
     if(y<=yy_pos)then   
         mp2(ix,iy)=-1.d0-Minfty*snalpha*((csalpha*(+y+x_sqrt)+snalpha*(-x+y_sqrt)) & 
                                           /(csalpha*(-x+y_sqrt)-snalpha*(y+x_sqrt)))
         mpara2=mp2(ix,iy)
         dn2=dninfty*exp(-abs(mpara2-Minfty*csalpha))
         tdn2(ix,iy)=dn2
     else
         mpara2=Minfty*csalpha
         mp2(ix,iy)=mpara2
         dn2=dninfty
         tdn2(ix,iy)=dn2
     endif

!******************      ATTENTION ATTENTION   *********************
!******** When there is an angle between B and Vsw this condition is required********
!******************      ATTENTION ATTENTION   *********************
   !  if(y>=y_test )then 
   !     mpara1=mp1(ix,iy-1)
   !     mp1(ix,iy)=mpara1
   !     dn1=dninfty*exp(-abs(mpara1-Minfty*csalpha))! or!dn1=tdn1(ix,iy-1)
   !     tdn1(ix,iy)=dn1
   !  endif
!******************      ATTENTION ATTENTION   *********************
!*******************************************************************
!******************      ATTENTION ATTENTION   *********************
          dn=dn1+dn2
          tdn(ix,iy)=dn

!  set floor for Density
     !   if(tdn1(ix,iy)<=1.d-10)then
     !      tdn1(ix,iy)=1.d-10
     !   endif
     !   if(tdn2(ix,iy)<=1.d-10)then
     !      tdn2(ix,iy)=1.d-10
     !   endif
     !   if(tdn(ix,iy)<=1.d-10)then
     !      tdn(ix,iy)=1.d-10
     !   endif


 
        write(155,102)x,y,tdn1(ix,iy),mp1(ix,iy),tdn2(ix,iy),mp2(ix,iy),tdn1(ix,iy)*mp1(ix,iy), &
                      tdn2(ix,iy)*mp2(ix,iy),tdn1(ix,iy)*mp1(ix,iy)+tdn2(ix,iy)*mp2(ix,iy),tdn1(ix,iy)+ &
                      tdn2(ix,iy)

        enddo
       write(155,*)
      enddo
      
   
      
  !stop
  
!////////////////////////////////////////////////////////////////////////////////
!////////////////// for comparison with measurment//////////////////////////////
!//////////////////////////////////////////////////////////////////////////////



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
        
        



!  5. Read position of satellite from file
      nlines = 0 
      open(unit=555,file='THB_L1_STATE_27227_finalEdited.txt1')
      DO
      READ(555,*,END=10)
      nlines = nlines + 1
      END DO
10      CLOSE (555)

    !  print*, nlines
      
     ! Rotate:
      open(unit=555,file='THB_L1_STATE_27227_finalEdited.txt1')
      open(unit=666,file='interpolated_ARTEMIS_pos.out',status='unknown')
     ! ARTEMIS produces 3D values. My code is a cylinderical 2D. In my code the plane which 
     ! contains B and V_sw are important (the cuts of circle in cylinder) and the direction 
     ! which is perpendicular to this plane is ignorable. 
     ! according to the SSE coordinate system, my y is ARTEMIS z. My x is negative ARTEMIS x.
     ! and ARTEMIS y is ignorale. We are looking for value mostly in the wake. If pos_y is bigger than
     ! ymax or smaller than ymin the position of satellite is outside of the wake. 
     ! if -rmoon<y<rmoon then if -rmoon<z<rmoon it is in the wake otherwise outside of the wake.
     
      do t=1,nlines
   	time = 1.0 * t 
         READ(555,*) ut,xx,yy,zz

          if (yy<rmoon_real .and. yy>-rmoon_real) then
             cap_radius = sqrt(rmoon_real*rmoon_real - yy*yy)
	     pos(1,1) = -xx/cap_radius 
             pos(2,1) =  zz/cap_radius
             pos(3,1) =  0.d0
             
                     ! print*,ut,pos(1,1),pos(2,1)
            
            !print*,'befor rotation',pos(1,1),pos(2,1)
            
             pos_rotate=matmul(Rz,pos)
             pos_new=matmul(Rx,pos_rotate)
             xx=pos_new(1,1)
             yy=pos_new(2,1)

            !print*,'after rotation',xx,yy
           
          
             if( xx>=xmin .and. xx<=xmax ) then
             
           
                ixb=int((xx-xmin)/dx)
                iyb=int((yy-ymin)/dy)
                xb=xmin+((ixb)*dx)
                yb=ymin+((iyb)*dy)

              ! print*, xx,xb,yy,yb
	        ww1=1.d0-((xx-xb)/dx)
                ww2=1.d0-((yy-yb)/dy)
              
          
                dn1_inter= (tdn1(ixb,iyb))*ww1*ww2+ &
                       (tdn1(ixb+1,iyb))*(1.d0-ww1)*ww2+ &
                       (tdn1(ixb+1,iyb+1))*(1.d0-ww1)*(1.d0-ww2)+ &
                       (tdn1(ixb,iyb+1))*ww1*(1.d0-ww2) 

                dn2_inter= (tdn2(ixb,iyb))*ww1*ww2+ &
                       (tdn2(ixb+1,iyb))*(1.d0-ww1)*ww2+ &
                       (tdn2(ixb+1,iyb+1))*(1.d0-ww1)*(1.d0-ww2)+ &
                       (tdn2(ixb,iyb+1))*ww1*(1.d0-ww2) 

                mp1_inter= mp1(ixb,iyb)*ww1*ww2+ &
                       mp1(ixb+1,iyb)*(1.d0-ww1)*ww2+ &
                       mp1(ixb+1,iyb+1)*(1.d0-ww1)*(1.d0-ww2)+ &
                       mp1(ixb,iyb+1)*ww1*(1.d0-ww2) 

                mp2_inter= mp2(ixb,iyb)*ww1*ww2+ &
                       mp2(ixb+1,iyb)*(1.d0-ww1)*ww2+ &
                       mp2(ixb+1,iyb+1)*(1.d0-ww1)*(1.d0-ww2)+ &
                       mp2(ixb,iyb+1)*ww1*(1.d0-ww2) 
                       
               ! print*,time,xb,yb

                !print*,dn1_inter+dn2_inter

                !print*,

                write(666,102) ut,xx,yy,dn1_inter,mp1_inter,dn2_inter,mp2_inter,dn1_inter*mp1_inter, &
                          dn2_inter*mp2_inter,dn1_inter*mp1_inter+dn2_inter*mp2_inter, dn1_inter+dn2_inter
               
            
             else !if( xx>=xmin .and. xx<=xmax ) then
               ! print*,'out of x boundary'
! ?????????????????/I am not sure????????????
              !  dn1_inter=5.d-1
             
                !dn1_inter=0.d0!tdn1(ixb,iyb)!dninfty/2.d0
               !dn2_inter=0.d0!tdn2(ixb,iyb)!0.d0!tdn1(ix,iy)
               print*,'out of x boundary'




              !  mp1_inter=mpinfty
              !  mp2_inter=mpinfty
         !      dn    = dninfty
!   ???????????????????????????????????
             endif !if( xx>=xmin .and. xx<=xmax ) then

          else !(if (yy<rmoon_real .and. yy>-rmoon_real) then)
            cap_radius = 0.0
            Print*, t,'no intersection '
            ! ?????????????????/I am not sure????????????
            dn1_inter=dninfty/2.d0
            dn2_inter=dninfty/2.d0
            mp1_inter=mpinfty
            mp2_inter=mpinfty
          ! dn    = dninfty
!   ???????????????????????????????????
          endif !(if (yy<rmoon_real .and. yy>-rmoon_real) then)
 
        write(666,102) ut,xx,yy,dn1_inter,mp1_inter,dn2_inter,mp2_inter,dn1_inter*mp1_inter, &
                        dn2_inter*mp2_inter,dn1_inter*mp1_inter+dn2_inter*mp2_inter, dn1_inter+dn2_inter
                        

             
      enddo
!////////////////////////////////////////////////////////////////////////////////
!////////////////// for comparison with measurment//////////////////////////////
!//////////////////////////////////////////////////////////////////////////////

102   format(99es15.5)
103   format(a,99es15.5)
      stop
      end
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
!=======================================================================
      subroutine gaussElim(aa,bb,n,nrhs)
      implicit none
!  solve nrhs equations of the form amt*X=bmt. amt is a n*n matix, and
!  bmt is a nrhs column matix with n lines. The solution will be stored
!  in bmt.
!  use simple Gauss elimination with line pivoting.

!  arguments
      integer n,nrhs
      real aa(n,n),bb(n,nrhs)

!  local variables
      integer i,j,k,pivot
      real r,ta(n),tb(nrhs)

!  computation

!  1.  first, scale every row by the largest element in the row
      do i=1,n
        r=abs(aa(i,1))
        do j=2,n
          r=max(r,abs(aa(i,j)))
        enddo
        aa(i,:)=aa(i,:)/r
        bb(i,:)=bb(i,:)/r
      enddo

!  2.  proceed with the Gauss elmination
      do i=1,n
!  2.1 find the line for pivoting
        r=abs(aa(i,i))
        pivot=i
        do j=i+1,n
          if(abs(aa(j,i)) > r) then
            pivot=j
            r=abs(aa(j,i))
          endif
        enddo
        if(pivot > i) then
          ta=aa(i,:)
          tb=bb(i,:)
          aa(i,:)=aa(pivot,:)
          bb(i,:)=bb(pivot,:)
          aa(pivot,:)=ta
          bb(pivot,:)=tb
        endif

!  2.2 eliminate from line i
        r=aa(i,i)
        aa(i,:)=aa(i,:)/r
        bb(i,:)=bb(i,:)/r
        do j=1,i-1
          r=aa(j,i)
          aa(j,:)=aa(j,:)-r*aa(i,:)
          bb(j,:)=bb(j,:)-r*bb(i,:)
        enddo
        do j=i+1,n
          r=aa(j,i)
          aa(j,:)=aa(j,:)-r*aa(i,:)
          bb(j,:)=bb(j,:)-r*bb(i,:)
        enddo

!  2.3 rescale remaining lines
        do k=i+1,n
          r=abs(aa(k,i+1))
          do j=i+2,n
            r=max(r,abs(aa(k,j)))
          enddo
          aa(k,:)=aa(k,:)/r
          bb(k,:)=bb(k,:)/r
        enddo
      enddo

 
      return
      end
!=======================================================================
!   This computer program is written by Tao Pang in conjunction with  
!    his book, "An Introduction to Computational Physics," published   
!    by Cambridge University Press in 1997.
!    Subroutine for evaluating the determinant of a matrix using 
!    the partial-pivoting Gaussian elimination scheme.
!    Copyright (c) Tao Pang 2001.                           

      function f3(x) result(D)
      use mod1
      IMPLICIT NONE
!   local variables
      INTEGER I,J,MSGN,INDX(N)
      Double Precision d,A_MuB(n,n),X
      DO i=1,n
        DO j=1,n
           A_MuB(i,j)=sin(X)*A(i,j)-cos(X)*B(i,j)
        enddo
      enddo
      CALL ELGS(A_MuB,INDX)
      D = 1.0
      DO I = 1, N
         D = D*A_MuB(INDX(I),I)
      END DO
      MSGN = 1
      DO I = 1, N
        DO WHILE (I.NE.INDX(I))
          MSGN = -MSGN
          J = INDX(I)
          INDX(I) = INDX(J)
          INDX(J) = J
        END DO
      END DO
      D = MSGN*D
      do ir=1,nr
      D=D/((x-miu(ir)))
      enddo
      return
      end
!=======================================================================
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.

      SUBROUTINE ELGS (A_MuB,INDX)
      use mod1
      IMPLICIT NONE
!   local variables
      INTEGER I,J,K,ITMP,INDX(N)
      Double Precision C1,PI,PI1,PJ,A_MuB(N,N),C(N)
!   Initialize the index
      DO I = 1, N
         INDX(I) = I
      END DO
!   Find the rescaling factors, one from each row
      DO I = 1, N
         C1= 0.0
        DO J = 1, N
           C1 = AMAX1(C1,ABS(A_MuB(I,J)))
        END DO
        C(I) = C1
      END DO
!   Search the pivoting (largest) element from each column


      DO J = 1, N-1
         PI1 = 0.0
         K=J !!!!!!!!!! I add this line to this routine, Richard suggested
                      !(the code had problem before adding this and INDX(K) was oubound)
                      !but if PI<PI1(the if condition does not satisfy)then K can have any value
                      ! because nothing has assigned to K for when if is not satisfied.
                      ! then we should assign K=J before this condition.
        DO I = J, N
           PI = ABS(A_MuB(INDX(I),J))/C(INDX(I))
           IF (PI.GT.PI1) THEN
              PI1 = PI
              K   = I
           ENDIF
        END DO

!   Interchange the rows via INDX(N) to record pivoting order
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP

        DO I = J+1, N
           PJ  = A_MuB(INDX(I),J)/A_MuB(INDX(J),J)
!   Record pivoting ratios below the diagonal
           A_MuB(INDX(I),J) = PJ
!   Modify other elements accordingly
           DO K = J+1, N
              A_MuB(INDX(I),K) = A_MuB(INDX(I),K)-PJ*A_MuB(INDX(J),K)
           END DO
        END DO
      END DO
      END SUBROUTINE ELGS
!=======================================================================
!  Purpose: find the zero of a REAL function fcn by the secant method.
!  arguments
!  x0: initial guess
!  dx0: initial step
!  dxmin: minimum acceptable step: aazero is accepted if dx.le.dxmin
!  fmin: minmum acceptable function value: aazero is accepted if
!        f.le.fmin*f1
!  xmin: lower bound of the x interval in which to find the root
!  xmax: upper bound of the x interval in which to find the root
!        N.B.: if xmin=xmax, no bounds are assumed.
!  itmax: maximum number of iterations allowable
!  .LT.  meaning <
!  .LE.          <=
!  .GT.          >
!  .GE.          >=
!  .EQ.          =
!  .NE.          /=

      double precision FUNCTION aazero(fcn,x0,dx0,dxmin,fmin,xmin,xmax,itmax)
      IMPLICIT NONE
      INTEGER itmax
      double precision x0,dx0,dxmin,fmin,xmin,xmax

!   Common variables: none

!   Local variables:
      INTEGER it
      double precision dxmax,x1,f1,x2,f2,zfmin,dx,dxinv,one,xl,xr,sgl,sgr,sg
      save one

!   Procedures:
      double precision fcn
      external fcn
      intrinsic abs,sign,min,max

!   Data:
      data one/1.d0/

!   Computation:

      dxmax=100.d0*abs(dx0)
      if(xmin.lt.xmax) then
        x1=max(xmin,min(xmax,x0))
      else
        x1=x0
      endif
      xl=x1
      xr=x1
 
      f1=fcn(x1)
      if(f1.eq.0.d0) then
        x2=x0
        go to 20
      endif
      zfmin=abs(fmin*f1)
      sgl=sign(one,f1)
      sgr=sgl
!     if(x1-xmin.gt.xmax-x1) then
!       dx=abs(dx0)
!     else
!       dx=-abs(dx0)
!     endif
      if(xmin.lt.xmax) then
        if(x1.eq.xmax) then
          dx=abs(dx0)
        elseif(x1.eq.xmin) then
          dx=-abs(dx0)
        else
          dx=-dx0
        endif
      else
        dx=-dx0
      endif
      do 10 it=1,itmax
      x2=x1-dx
      if(xmin.lt.xmax) then
        if(sgl.ne.sgr) then
          x2=max(xl,min(xr,x2))
        else
          x2=max(xmin,min(xmax,x2))
        endif
        if((abs(x2-x1).le.dxmin .and. it.gt.2) .or. x1.eq.x2) go to 30
      endif
      f2=fcn(x2)
      if(abs(f2).le.zfmin) go to 20
      sg=sign(one,f2)
      if(sgl.eq.sgr) then
        xl=min(xl,x2)
        xr=max(xr,x2)
        if(sg.ne.sgl) then
          if(x2.lt.xr) then
            sgl=sg
            xl=x2
          else
            sgr=sg
            xr=x2
          endif
        endif
      endif
      if(sgl.ne.sgr) then
        if(sg.eq.sgl) then
          xl=x2
          sgl=sg
        else
          xr=x2
          sgr=sg
        endif
      endif
 
      dxinv=(f2-f1)/(f2*(x2-x1))
      if(abs(dxinv)*dxmax.le.1.) then
        dx=sign(dxmax,dxinv)
      else
        dx=1./dxinv
        if(abs(dxinv)*dxmin.ge.1. .and. it.gt.1) go to 30
      endif
 
      if(sgl.ne.sgr) then
        if(sg.ne.sgl) then
!   shoot to the left
          if(dx.gt.0.) then
            dx=min(dx,0.8*(x2-xl))
          else
            dx=0.5*(x2-xl)
          endif
        else
!   shoot to the right
          if(dx.lt.0.) then
            dx=max(dx,0.8*(x2-xr))
          else
            dx=0.5*(x2-xr)
          endif
        endif
      endif

      x1=x2
      f1=f2
10    continue

!   no convergence
      WRITE(6,*)' warning: no convergence in aazero'
      aazero=x0
      RETURN

!   good convergence, the function is sufficiently small
20    continue
      aazero=x2
      RETURN

!   good convergence, step size is sufficiently small
30    continue
      aazero=x2-dx
      if(xmin.lt.xmax) then
        aazero=max(xmin,min(xmax,aazero))
      endif
      RETURN
      end


