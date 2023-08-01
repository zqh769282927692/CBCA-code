
C       ######################################
C       MUTILICATION OF TWO VECTORS
C       ###################################
        SUBROUTINE MUTIVT1(C2,C1)
         implicit real*8 (a-h,o-z) 
        real*8 C2(3,3),C1
        C2(3,1)=C2(1,2)*C2(2,3)-C2(1,3)*C2(2,2)
        C2(3,2)=C2(1,3)*C2(2,1)-C2(1,1)*C2(2,3)
        C2(3,3)=C2(1,1)*C2(2,2)-C2(1,2)*C2(2,1)
        C1=C2(3,1)**2+C2(3,2)**2+C2(3,3)**2
!        IF (C1.LT.1.0E-08) GOTO 120
        C1=SQRT(C1)
        C2(3,1)=C2(3,1)/C1
        C2(3,2)=C2(3,2)/C1
        C2(3,3)=C2(3,3)/C1
120     CONTINUE
        RETURN
        END

C       ######################################
C       MUTILICATION OF TWO VECTORS
C       ###################################
        SUBROUTINE MUTIVT(C2,C1)
         implicit real*8 (a-h,o-z) 
        DIMENSION C2(3,3) 
        C2(3,1)=C2(1,2)*C2(2,3)-C2(1,3)*C2(2,2)
        C2(3,2)=C2(1,3)*C2(2,1)-C2(1,1)*C2(2,3)
        C2(3,3)=C2(1,1)*C2(2,2)-C2(1,2)*C2(2,1)
        C1=C2(3,1)**2+C2(3,2)**2+C2(3,3)**2
!        IF (C1.LT.1.0E-08) GOTO 120
        C1=SQRT(C1)
     
        C2(3,1)=C2(3,1)/C1
        C2(3,2)=C2(3,2)/C1
        C2(3,3)=C2(3,3)/C1
          if(c1.lt.1.0E-08) then
!          write(1001,*)'c1=',c1
!          write(1001,*)'==', C2(3,1)**2+C2(3,2)**2+C2(3,3)**2    
          endif
120     CONTINUE
        RETURN
        END

C       #######################
C       FUNCTION DETERMINATION
C       ###########################
        SUBROUTINE FDETEM (B2,T)
        implicit real*8 (a-h,o-z) 
        DIMENSION B2(3,3) 
        T1=B2(1,1)*B2(2,2)*B2(3,3)-B2(1,1)*B2(2,3)*B2(3,2)
        T2=B2(1,2)*B2(2,3)*B2(3,1)-B2(1,2)*B2(2,1)*B2(3,3)
        T3=B2(1,3)*B2(2,1)*B2(3,2)-B2(1,3)*B2(2,2)*B2(3,1)
        T=T1+T2+T3
        RETURN
        END
C       ##########################
C       INTERSECTION OF THREE PLANES
C       #########################
        SUBROUTINE INTERP(P2,E,IW)
       implicit real*8 (a-h,o-z) 
       DIMENSION P2(3,4),E(3),A2(3,3),B2(3,3)

        DO 10 I=1,3
        DO 10 J=1,3
        B2(I,J)=P2(I,J)
10      CONTINUE
        CALL FDETEM(B2,T)
        S=T
        IW=1
        IF (ABS(T).GT.0.000001) GOTO 150
        IW=0
        GOTO 60
150     DO 20 I=1,3
        DO 20 J=1,3
20      A2(I,J)=B2(I,J)
        DO 50 L=1,3
        DO 30 I1=1,3
        DO 30 I2=1,3
30      B2(I1,I2)=A2(I1,I2)
        B2(1,L)=-P2(1,4)
        B2(2,L)=-P2(2,4)
        B2(3,L)=-P2(3,4)
        CALL FDETEM(B2,T)
        E(L)=T/S
50      CONTINUE
60      CONTINUE
        RETURN
        END

C       ##############################
C       A GIVEN PROJECTION SYSTEM
C       #############################
        SUBROUTINE PROJTG(F,A2)
        implicit real*8 (a-h,o-z) 
        DIMENSION F(3),A2(3,3),C2(3,3)
        IF (F(1)**2+F(2)**2.GT.0.0000001) GOTO 200
        A2(3,1)=0.0
        A2(3,2)=0.0
        A2(3,3)=SGN(F(3)+1.2E-07)
        A2(2,1)=0.0
        A2(2,2)=SGN(F(3)+1.2E-07)
        A2(2,3)=0.0
        A2(1,1)=1.0
        A2(1,2)=0.0
        A2(1,3)=0.0
        GOTO 400
C       ---------------------------------
200     FF1=SQRT(F(1)**2+F(2)**2+F(3)**2)
        A2(3,1)=F(1)/FF1
        A2(3,2)=F(2)/FF1
        A2(3,3)=F(3)/FF1
        C2(1,1)=0.0
        C2(1,2)=0.0
        C2(1,3)=1.0
        C2(2,1)=A2(3,1)
        C2(2,2)=A2(3,2)
        C2(2,3)=A2(3,3)
        CALL MUTIVT(C2,C1)
        DO 350 I0=1,3
        A2(1,I0)=C2(3,I0)
        C2(2,I0)=C2(3,I0)
        C2(1,I0)=A2(3,I0)
350     CONTINUE
        CALL MUTIVT(C2,C1)
        A2(2,1)=C2(3,1)
        A2(2,2)=C2(3,2)
        A2(2,3)=C2(3,3)
400     CONTINUE
        RETURN
        END

C       ##################
        FUNCTION SGN(X)
        implicit real*8 (a-h,o-z) 
        IF (X) 10,20,30
10      SGN=-1.0
        GOTO 40
20      SGN=0.0
        GOTO 40 
30      SGN=1.0
40      CONTINUE
        END

C       ##########################
C       COORDINATES TRANSFORM FROM TOTAL TO LOCAL SYSTEM
C       ##########################
        SUBROUTINE TLXYZ(A2,ET,EL)
        implicit real*8 (a-h,o-z) 
        DIMENSION A2(3,3),ET(3),EL(3)
        DO 10 I=1,3
10      EL(I)=0.0
        DO 20 I=1,3
        DO 20 J=1,3
        EL(I)=EL(I)+A2(I,J)*ET(J)
20      CONTINUE
        END

C       #################################
C       COORDINATES TRANSFORM FROM LOCAL TO TOTAL SYSTEM
C       #################################
        SUBROUTINE LTXYZ(A2,ET,EL)
        implicit real*8 (a-h,o-z) 
        DIMENSION A2(3,3),ET(3),EL(3),B2(3,3)
        DO 10 I=1,3
        DO 10 J=1,3
10      B2(I,J)=A2(I,J)
        CALL FDETEM(B2,T)
        S=T
        DO 50 I=1,3
        DO 20 J=1,3
        DO 20 L=1,3
20      B2(J,L)=A2(J,L)
        DO 30 K=1,3
30      B2(K,I)=EL(K)
        CALL FDETEM(B2,T)
        ET(I)=T/S
50      CONTINUE
        END
C       #########################
C       IF POINT (X0,Y0) BE CONTAINED IN A GIVEN BLOCK REGION
C       VERTICES: NVEF; COORDINATES: CVEF(20,2)
C       ########################
        SUBROUTINE POINTC(X0,Y0,NVEF,CVEF,N1)
        implicit real*8 (a-h,o-z) 
        DIMENSION CVEF(800,2),H(800)
        common/Control/Ctl1,Ctl2 
           
        X4=X0
        Y4=Y0
        DO 150 I2=1,NVEF
        X1=CVEF(I2,1)-X4
        Y1=CVEF(I2,2)-Y4
        A1=SQRT(X1**2+Y1**2)
        IF (A1.GT.0.0001) GOTO 140
        N1=0
        GOTO 550
140     H(I2)=ABS(X1)/A1
150     CONTINUE
170     CALL RANDU1(E1)
        DO 210 I2=1,NVEF
        IF (ABS(H(I2)-E1).GT.0.001) GOTO 210
        GOTO 170
210     CONTINUE
        L0=0
        DO 470 I2=1,NVEF
        I3=I2+1
        IF (I3.GT.NVEF) I3=I3-NVEF

        X1=CVEF(I2,1)
        Y1=CVEF(I2,2)
        X2=CVEF(I3,1)
        Y2=CVEF(I3,2)
        A1=Y1-Y2
        B1=X2-X1
        C1=X1*Y2-X2*Y1
        SL1=(A1*X4+B1*Y4+C1)/SQRT((X1-X2)**2+(Y1-Y2)**2) !�����߶εľ���
        X10=E1
        Y10=SQRT(1.-X10**2)
        X3=X4+X10*2000.
        Y3=Y4+Y10*2000.
        SL2=SGN(A1*X3+B1*Y3+C1)
        A1=Y3-Y4
        B1=X4-X3
        C1=X3*Y4-X4*Y3
        SL3=SGN(A1*X1+B1*Y1+C1)
        SL4=SGN(A1*X2+B1*Y2+C1)
        IF (SL3*SL4.GE.0.0) GOTO 470
!                                        <�ݲ��������>
        IF (ABS(SL1).LT.1e-4*Ctl1) GOTO 510  !�ܽ�ʱ��Ϊ��������.JTSEF2���ӳ����У����������ڲ��ĵ㣬ֱ�Ӵ洢���ⲿ�ĵ���Ҫ�󽻡�
!                       �����߶εĵ㵱�������ڵĵ㣬���������ݲ��ı䡣Ȼ����intsb������߶�������߽�Ľ��㡣 2015.5.26
!                        ��prep23dһ�£�1e-8��Ϊ1e-4  2018.1.8
        SL1=SGN(SL1)
        IF (SL1*SL2.GT.0.0) GOTO 470
        L0=L0+1
470     CONTINUE
        N1=L0-INT(L0/2+0.00001)*2
       
        N1=1-N1
        GOTO 550
510     N1=0          !��Ϊ�����ڣ��������ݲ����Ҫƥ�� 2014.8.31
550     CONTINUE
        RETURN
        END
CC####################################################
C      /* mc19: finding a inner point of a loop     */
CC####################################################
        SUBROUTINE MC19(V,J0,X,Y)
       implicit real*8 (a-h,o-z) 
       DIMENSION V(500,2)
        ddx1=9999999.9
        ddy1=9999999.9
        ddx2=-9999999.9
        ddy2=-9999999.9
      
       DO 10 l=1,J0
        if (v(l,1).lt.ddx1) ddx1=v(l,1) 
        if (v(l,2).lt.ddy1) ddy1=v(l,2) 
        if (v(l,1).gt.ddx2) ddx2=v(l,1) 
        if (v(l,2).gt.ddy2) ddy2=v(l,2) 
10     continue
         w0=ddx2-ddx1
	if (ddy2-ddy1.gt.w0) w0=ddy2-ddy1
       w0=w0*10

       x3 =(v(1,1)+v(2,1))/2
       y3 =(v(1,2)+v(2,2))/2
       x00= v(1,2)-v(2,2)         
       y00= v(2,1)-v(1,1)
       a1 = SQRT(x00*x00 + y00*y00)
       x4 = x3 + x00*w0/a1       !!�������߷����֪,�����ظ���Ҳû������
       y4 = y3 + y00*w0/a1

       a2=1 
       DO 100 l=2,J0
      
       x01=v(l,1)           !! v( , )�Ķ�������J0+1��
       y01=v(l,2)
       x2 =v(l+1,1)
       y2 =v(l+1,2)

       x21=x2 - x01
       y21=y2 - y01
       x34=x3 - x4 
       y34=y3 - y4 

       x31=x3 - x01
       y31=y3 - y01

       d3=x21*y34-x34*y21
       d1=x31*y34-x34*y31
       d2=x21*y31-x31*y21
       
       if (abs(d3).LT. .0000001*w0*w0) goto 90   

       t1=d1/d3                !!��01��2 �ϵĴ���ϵ��(����˵�����λ��0��1)
       if (t1.LT.-.00000001.OR.t1.GT.1.00000001) goto 90
       t2=d2/d3                !!��3��4 �ϵĴ���ϵ��(����˵�����λ��0��1)
       if (t2 .LT.0.0000001) goto 90
       if (a2 .LE.  t2) goto 90
       a2=t2
90     CONTINUE
100      CONTINUE 
 
C------/*--------- inner point-----------------*/
       X = x3 + .5 * a2 * x00 * w0/a1
       Y = y3 + .5 * a2 * y00 * w0/a1
      
	RETURN
      END
      
CC####################################################
C   Volume of a tetrahedron  
CC####################################################         
          SUBROUTINE  TET_VOL(B3,T)  
            implicit real*8 (a-h,o-z) 
           DIMENSION B3(4,3),A1(3,3) 
          
          A1(1,1:3) = B3(2,1:3)-B3(1,1:3)
          A1(2,1:3) = B3(3,1:3)-B3(1,1:3)
          A1(3,1:3) = B3(4,1:3)-B3(1,1:3)
          CALL FDETEM(A1,T)
           T=T/6.0     !���������
      ENDSUBROUTINE      
      
C       ########################################
CC--����������̬�ĵ�������������. 2008.6.6  -
C       ########################################
       SUBROUTINE mc12(LJ,SAI)
	 implicit real*8 (a-h,o-z) 
       DIMENSION NDLPZ(1140000,4),NLOC(5500000),CTT(490000,3)
	 DIMENSION LPLN(1140000)
       DIMENSION P1J(90000,9),B2(3,3),V1(800,2),ET(3),EL(3)
       COMMON/NDLPZ/NDLPZ/P1J/P1J/NLOC/NLOC/LPLN/LPLN/CTT/CTT 

	  I1=ABS(NDLPZ(LJ,1))
	  DO 20 J1=1,3
        DO 20 J2=1,3
20       B2(J1,J2)=P1J(I1,(J1-1)*3+J2)

	 DO 80 J=1,NDLPZ(LJ,3)
        J1=NLOC(LPLN(LJ)+J)
	  DO L1=1,3
       ET(L1)=CTT(J1,L1)    	   
        ENDDO  
	  CALL TLXYZ(B2,ET,EL)
        DO L1=1,2
       V1(J,L1)=EL(L1)    	   
        ENDDO
80     CONTINUE

       H0=0
       DO 100 j=1,NDLPZ(LJ,3)
C      /* first node of triangle is (0,0)       */
       x2= V1(J,1)
       y2= V1(J,2)
	 if (j.lt.NDLPZ(LJ,3))then
       x3= V1(J+1,1)
       y3= V1(J+1,2)
	 else
       x3= V1(1,1)
       y3= V1(1,2)
       endif
       f1=(x2*y3-x3*y2)/2
       h0=H0+f1
100    CONTINUE
       SAI=ABS(H0)
      RETURN
	END
 

CC####################################################
C   Triangulation of an arbitrarily-shaped loop          
CC####################################################
        SUBROUTINE Triangulation(V,J0,NNEB,Ntri)
	  IMPLICIT real*8(A-H,O-Z)
       DIMENSION V(800,2),VD(2,2),Ang(800,3),VE(4,2)
       DIMENSION NNEB(3000,3),K(800,3)
       
       IF(J0.EQ.3) THEN
          Ntri=1
          NNEB(Ntri,1)= 1;NNEB(Ntri,2)= 2;NNEB(Ntri,3)= 3	    
           GOTO 1000
	  ENDIF
        
         DO 10 I=1,J0   ! ��·��ÿ����ĺ�㡢�õ㡢ǰ�㣬���δ洢Ϊ1,2,3
 	  K(I,1)=I-1
 	  K(I,2)=I    ! K( ,2)�ǲ�����(���Ϊ-1)��ֱ���õ��γ������κ�����
         K(I,3)=I+1 
 10      continue
 	  K(1,1)=J0
 	  K(J0,3)=1
          
        DO 20 I=1,J0 
        VD(1,1)=V(K(I,2),1)
        VD(1,2)=V(K(I,2),2)
        VD(2,1)=V(K(I,3),1)         !�õ���ǰ����߶�
        VD(2,2)=V(K(I,3),2)
        write(12,*)'check', VD(1,1), VD(1,2)
       CALL dipdirection(VD,d1)
         Ang(K(I,2),3)=d1       !�õ��ǰ��ߵĽǶ�(��X��Ϊ���,��ʱ����ת)
         Ang(K(I,3),1)=d1-180   !�õ��ǰ��ĺ���ߵĽǶ�
	  if(Ang(K(I,3),1).lt.0) Ang(K(I,3),1)=Ang(K(I,3),1)+360
20      continue	  
 
        DO 30 I=1,J0 
         Ang(K(I,2),2)=Ang(K(I,2),1)-Ang(K(I,2),3)     !��·ÿ������ڽ�
       if(Ang(K(I,2),2).lt.0) Ang(K(I,2),2)=Ang(K(I,2),2)+360
30      continue	 

!--------------------------- 2017.2.4  ----
!         goto 70
        IF(J0.EQ.4)THEN          !������С���г�ʱ��
	                           !���ܳ��������λ��Σ���Ӱ�쵽���������廮��
          ang0=-0.1              !�����廮���У��󲿷ֵ�
		                         ! Triangulation���4������Ķ����
        DO 60 I=1,J0
        if(Ang(K(I,2),2).GT.ang0) then    !������ڽ�,�Ӹ��ڽǻ���Ϊ���������� 		
	   ang0=Ang(K(I,2),2)             
	   k0=I             
	  endif
60      continue	
       
        if(k0.gt.1)then
        k0=k0-1          !����ڽǵ�����һ���㣨������ǰ����㣩��Ϊ���� 
        else
         k0=4    
        endif
        
        Ntri=2
        do 65 j=1,Ntri
        if (j.eq.1)goto 66    
         if(k0.le.2)then
          k0=k0+2
        else
          k0=k0-2   
        endif
        
66      j3=K(k0,3)            
        j1=K(k0,1)          
        NNEB(j,1)= J1                  !�γ�������
        NNEB(j,2)= K(k0,2)
        NNEB(j,3)= J3
 
65      continue	      
        RETURN
        ENDIF
!--------------------------        
        
70	  Ntri=0
50        ang0=359.99999 
         i0=0
                 
	 do 51 I=1,J0      !  2015.1.8
         if(K(I,2).eq.-1) goto 51      
51     continue          !  2015.1.8

        DO 100 I=1,J0 
         if(K(I,2).eq.-1) goto 100
 	    i0=i0+1

!---------�ų�����������·����������ཻ������      
  	    ki1=K(i,1)
  	    ki2=K(i,2)
	    ki3=K(i,3)
         VE(1,1)=V(K(I,1),1); VE(1,2)=V(K(I,1),2)
         VE(2,1)=V(K(I,3),1); VE(2,2)=V(K(I,3),2)
        
        DO 80 J=1,J0 
	 if(K(j,2).eq.-1.or.j.eq.i) goto 80
	    kj2=K(j,2)
	    kj3=K(j,3)
          
        if(kj3.eq.ki1.or.kj3.eq.ki2.or.kj3.eq.ki3    !2014.12.5
     # .or.kj2.eq.ki1.or.kj2.eq.ki2.or.kj2.eq.ki3) goto 80 !����
	     !��·ǰ������,kj2-kj3��ߵ�һ������K(i,2)��ǰ��㹲��ʱ
         VE(3,1)=V(Kj2,1); VE(3,2)=V(Kj2,2)
         VE(4,1)=V(kj3,1); VE(4,2)=V(kj3,2)          
         
 	   call intofline(VE,k10)
	  IF(k10.EQ.1) then
!		 write(1001,*)'1-2',VE(1,1), VE(1,2), VE(2,1), VE(2,2)
!           write(1001,*)kj3,kj2,ki1,ki2,ki3
! 	     write(1001,*)'Tringle edge is intersected with loop edge!'
		  GOTO 100
	  endif
80      CONTINUE
                        
!---------�ų�������������ڵĻ�·���С�Ƕ��ཻ���Ա��������λ���  2017.2.8      
!-----  ���ε������ο��ܵ��¸��������붥���γɵ�����������Ϊ��                        
  	    ki1=K(i,1)
  	    ki2=K(i,2)
	    ki3=K(i,3)
         x1=V(K(I,1),1); y1=V(K(I,1),2)
         x2=V(K(I,3),1); y2=V(K(I,3),2)
        
        DO 85 J=1,J0 
	 if(K(j,2).eq.-1.or.j.eq.i) goto 85
	    kj2=K(j,2)
	    kj3=K(j,3)
          
        if(kj3.eq.ki1.and.kj2.ne.ki3.or.
     #     kj2.eq.ki3.and.kj3.ne.ki1) then  
                       !�����ǰ������ڵ����,����,��ֹ���һ�����������˴����ж�
	    
         x3=V(Kj2,1); y3=V(Kj2,2)
         x4=V(kj3,1); y4=V(kj3,2)          
         
        t1=((x1-x2)*(x4-x3)+(y1-y2)*(y4-y3))         ! �߶εļн� 2017.2.8
     #        /sqrt((x1-x2)**2+(y1-y2)**2)/sqrt((x4-x3)**2+(y4-y3)**2)
 	   
	    IF(t1.gt.0.9999985) then     
		            !��С�Ƕ�>3�ȣ�cos(3)=0.9986 cos(1)=0.99848  2017.6.20
                      !   cos(0.1)=0.9999985   2018.2.9
              write(1001,*)'��������ڵĻ�·���С�Ƕ��ཻ',i,kj2,kj3
		  GOTO 100      !����ܻ���ʱ���ڴ˶�  GOTO 100���Ӷ��ò���k0������K(k0,)�������
                          !��˸�Ϊt1.gt.0.9999985   2018.2.9
          endif
         endif 
85      CONTINUE
        
!---�ų�����������ߴ����ص������ݻ�·���ʣ�����ֻ�ܱ���߳�������ߵ�һ����
!--- λ��������ʱ����Ҫ�ų�       
         VE(1,1)=V(K(I,1),1); VE(1,2)=V(K(I,1),2)
         VE(2,1)=V(K(I,3),1); VE(2,2)=V(K(I,3),2)        
        DO 90 J=1,J0 
	 if(K(j,2).eq.-1.or.
     $  K(j,2).eq.K(I,1).or.K(j,2).eq.K(I,2).       !���������I�����ڵ������߶�
     $    or.K(j,2).eq.K(I,3).and.K(j,3).eq.K(I,1)) goto 90 !��Ӧ�㹲��,ʣ��3��
                                                     !�����ù���������ʱ������
	    kj2=K(j,2)
	    kj3=K(j,3)
	   VE(3,1)=V(Kj2,1); VE(3,2)=V(Kj2,2)
         VE(4,1)=V(kj3,1); VE(4,2)=V(kj3,2)          
          
        call overlapline(VE,k10)   
        
	  IF(k10.EQ.1) then
  	     write(1001,*)'Tringle edge overlaps a loop edge!'
		  GOTO 100
	  endif
90      CONTINUE
!---�ų�����ߵ��е�λ�ڻ�·�����������,�������εı߲���λ��������    
	  X0=(V(K(I,1),1)+V(K(I,3),1))/2   
	  Y0=(V(K(I,1),2)+V(K(I,3),2))/2
        
        call POINTC(X0,Y0,J0,V,N1)  !   ,N2,XD,YD)
        if(N1.eq.1) then  !.or.n1.eq.0.and.n2.ne.0  2014.12.7       
! 	     write(1001,*)'Tringle edge is located outside the loop!',X0 ,Y0
	     goto 100   
        endif
 
! ---  ����С�ڽ�    
        if(Ang(K(I,2),2).lt.ang0) then    !����С�ڽ�,���ڽ����ڵ�������������			
	   ang0=Ang(K(I,2),2)               !�γ�������
	   k0=I             
	  endif
100      continue	 
      
!       write(1001,*)'i01=',i0,k0  !,K(k0,2)  
             
        j3=K(k0,3)  
        j1=K(k0,1)  
        if(i0.eq.3) goto 110        !��ֻ��3������δ��ʶΪ-1,ѭ������

 	  K(j1,3)=j3                  ! k0�ĺ����k0ǰ��ֱ�ӽ�����ϵ,�γ�һ�����
	  K(j3,1)=j1

        VD(1,1)=V(j1,1)              !�������γɵ���ߵĽǶ�
        VD(1,2)=V(j1,2)
        VD(2,1)=V(j3,1)      
        VD(2,2)=V(j3,2)	        
        CALL dipdirection(VD,d1)
        Ang(j1,3)=d1         
        Ang(j3,1)=d1-180         
	  if(Ang(j3,1).lt.0) Ang(j3,1)=Ang(j3,1)+360
         Ang(j1,2)=Ang(J1,1)-Ang(J1,3)  ! �µ��ڽ�
       if(Ang(j1,2).lt.0) Ang(j1,2)=Ang(j1,2)+360
         Ang(j3,2)=Ang(J3,1)-Ang(J3,3)  ! �µ��ڽ�
       if(Ang(j3,2).lt.0) Ang(j3,2)=Ang(j3,2)+360

       if(K(k0,2).le.0) goto 1000   

110         Ntri=Ntri+1        
        NNEB(Ntri,1)= J1
        NNEB(Ntri,2)= K(k0,2)
        NNEB(Ntri,3)= J3
!        write(1001,*)'NNEB',  NNEB(Ntri,1:3)
        if(i0.eq.3) goto 1000
	 
	  K(k0,2)=-1  
        goto 50 
     
1000   continue	 
 
	RETURN
        END
        
CC##################### 2011-4-19 #################
       SUBROUTINE dipdirection(VD,d1)
	   IMPLICIT real*8(A-H,O-Z)
	   real*8 VD(2,2)
	  DD=180.0/3.1415926535

	    X1=VD(2,1)-VD(1,1)
          Y1=VD(2,2)-VD(1,2)

          C1=SQRT(X1*X1+Y1*Y1+1.0E-09)
          YC=Y1/C1
       
	   IF (YC.GT.0.9999999) then
            D1=90.0  
            return
         endif
         IF (YC.LT.-0.9999999) then
            D1=270.0
            return
         endif
          
          D1=ATAN(YC/SQRT(1-YC*YC))*DD           !��X��ļн�      

        IF (X1.lt.0.0) D1=180.0-D1    
!--- ���о���180.0-D1�ǶԵġ��ڶ�����ʱ��D1>0�����Ϊ180.0-D1����������ʱ��
!--- ��D1<0����ʱҲΪ180.0-D1��������֮ǰ���и����������ǰ����û���⡣
!--- ��df�У�d1=atan(y12/x12)��������Ϊd1=180.0+d1��2018.5.3
        
        IF (D1.lt.0.0) D1=D1+360.0
     
	RETURN
      END
CC####################################################
C  �ж������߶��Ƿ�ƽ�У�����һ�����Ƿ�λ����һ�߶���2011-4-21  */
CC####################################################
         SUBROUTINE overlapline(V,k10)
	  IMPLICIT real*8(A-H,O-Z)
        real*8 V(4,2)
       
	 k10=0
       x1=v(1,1); y1=v(1,2)
       x2=v(2,1); y2=v(2,2)
 
       x3=v(3,1); y3=v(3,2)	 
	 x4=v(4,1); y4=v(4,2) 
 
       x21=x2 - x1
       y21=y2 - y1
       x34=x3 - x4 
       y34=y3 - y4 

        d3=x21*y34-x34*y21

       if (abs(d3).gt. .000001) goto 90  !2011-10-20,��0.00001Ϊ0.00005
	                               !Խ��Խ�ϸ�2014.11.26

!        SL1=-y21*x3 +x21*y3 +X1*Y2-X2*Y1  ! һ����λ������.����취������,��Ҫ���ϵ�      
!       IF (ABS(SL1).gt.0.0005) GOTO 90  ! λ���߶��ڵ��ж�,��������ʱ,�жϿ��ܲ�׼ȷ
        	 
       d12=sqrt((x2-X1) **2+(Y2-Y1)**2)
       d13=sqrt((x3-X1) **2+(Y3-Y1)**2)
       d23=sqrt((x3-X2) **2+(Y3-Y2)**2)    ! �����������,��ǰ���Ѿ�����
       if(d13+d23-d12.gt.0.000005) GOTO 90    ! �ܶ��ж�,������Ҫ����Ĺ��Ʋ���
	                                       !  ��0.0001Ϊ0.0005
       d14=sqrt((x4-X1) **2+(Y4-Y1)**2)
       d24=sqrt((x4-X2) **2+(Y4-Y2)**2)
       if(d14+d24-d12.gt.0.000005) GOTO 90  !2011-10-20,��0.0001Ϊ0.0005
         k10=1                 !!�����߶��ص�    
90     CONTINUE
 
      END 
CC####################################################
C    /*  �ж������߶��Ƿ���ڽ���  2011-4-21  */
CC####################################################
         SUBROUTINE intofline(V,k10)
	  IMPLICIT real*8(A-H,O-Z)
        real*8 V(4,2)
       
       x1=v(1,1); y1=v(1,2)
       x2=v(2,1); y2=v(2,2)
 
       x3=v(3,1); y3=v(3,2)	 
	 x4=v(4,1); y4=v(4,2) 
 
       x21=x2 - x1
       y21=y2 - y1
       x34=x3 - x4 
       y34=y3 - y4 

       x31=x3 - x1
       y31=y3 - y1

       d3=x21*y34-x34*y21
       d1=x31*y34-x34*y31
       d2=x21*y31-x31*y21
 
	 k10=0
       if (abs(d3).gt.1.0e-8)then  !Խ��Խ�ϸ�  !2015.1.7
	!��С�������߶ιؼ������Ƿ��ཻ���Ƿ�ӽ�ƽ��û����ϵ 2015.1.8
       t1=d1/d3                !!��01��2 �ϵĴ���ϵ��(����˵�����λ��0��1)
       if (t1.LT.0.0000001.OR.t1.GT.0.9999999) goto 90 !�ų��ڶ˵㴦�ཻ2015.1.7
       t2=d2/d3                !!��3��4 �ϵĴ���ϵ��(����˵�����λ��0��1)
       if (t2.LT.0.0000001.OR.t2.GT.0.9999999) goto 90 !�ų��ڶ˵㴦�ཻ2015.1.7
        k10=1;goto 90                 !!�����߶��н���   

	  else
! 	  k10=1                   ! �н���   2014.11. ƽ�ж���Ϊ�н��㲻���� 2017.2.7
        t1=(y1-y2)*x3+(x2-x1)*y3+x1*y2-x2*y1  
	                         ! (x3,y3)����(x1,y1)-(x2,y2)���߶η��̵õ�(x3,y3)��֮�ľ���
        if(abs(t1).lt.1e-6)then           ! 2017.2.7
          k10=1  
        else
           k10=0
        endif   
        endif   

90     CONTINUE
!        if(k10.eq.1)write(1001,*)'intof=',x3,y3,x4,y4,t1,t2,d3
       END
CC####################################################
C   ����������ϸ�֡���������򵥵ġ���ά˫�ؽ��������������߽�����϶��ʱҪ���ӽڵ㣻
C  ��һ����϶�������������У�һ�������ο��ܶ�Ӧ�ж����������Ҫϸ�� 2018.3.28        
C  �����������ӳ���û���⡣����·�Ķ�����10��ʱ�����������߳�Ϊģ�ͳ߶�/10,ϸ�ֳ��ĵ�Ԫ�ܿ쳬��300����
C  �����չNNEB(3000,3)  2018.10.18       
CC####################################################
        SUBROUTINE TriangleFine(V,J0,DT0,VRANGE,NNEB,Ntri)
	  IMPLICIT real*8(A-H,O-Z)
        DIMENSION V(800,2),NNEB(3000,3) 
       
 
209     DT1=0  

       DO 30 I=1,NTRI   
		X1=V(NNEB(I,1),1)
		Y1=V(NNEB(I,1),2)   
		X2=V(NNEB(I,2),1)  
		Y2=V(NNEB(I,2),2)    
		X3=V(NNEB(I,3),1)   
		Y3=V(NNEB(I,3),2)   
      
        DL0=SQRT((X1-X2)**2+(Y1-Y2)**2)   
!           Write(12,*)'X1=',X1,Y1 
!           Write(12,*)'X1=',X2,Y2
!           Write(12,*)'X1=',X3,Y3 
!           Write(12,*)'i=',i 
	   IF(DL0.GT.DT1+0.001) THEN
	    DT1=DL0               !�ҵ����ĵ�Ԫ�߳�
	   
	  X4= X1+(X2-X1)/2;Y4=Y1+(Y2-Y1)/2     !(X1+X2)/2; Y4=(Y1+Y2)/2, 2018.11.15
        ITRI=I; ITRIN=1
         ENDIF
        DL0=SQRT((X2-X3)**2+(Y2-Y3)**2)
	   IF(DL0.GT.DT1+0.001) THEN
	    DT1=DL0
	  
	   X4= X2+(X3-X2)/2;Y4=Y2+(Y3-Y2)/2     !X4=(X2+X3)/2; Y4=(Y2+Y3)/2 
        ITRI=I; ITRIN=2
         ENDIF
        DL0=SQRT((X3-X1)**2+(Y3-Y1)**2)
	   IF(DL0.GT.DT1+0.001) THEN
	    DT1=DL0
	    
	  X4= X3+(X1-X3)/2;Y4=Y3+(Y1-Y3)/2    !X4=(X3+X1)/2.0; Y4=(Y3+Y1)/2.0 
        ITRI=I; ITRIN=3 
         ENDIF
30      CONTINUE
!          Write(12,*)'DT00=',DT1,DT0
	 IF(DT1.LE.DT0) GOTO 309   ! DT0Ϊ�����ĵ�Ԫ���߳���ֵ

        J0=J0+1
	  V(J0,1)=X4; V(J0,2)=Y4    
         
       select case(ITRIN)
       case(1)
 	   ND12=NNEB(ITRI,1)   
	   ND13=NNEB(ITRI,2)
       case(3)
	   ND12=NNEB(ITRI,3)    
	   ND13=NNEB(ITRI,1)
       case(2)
 	   ND12=NNEB(ITRI,2)  
	   ND13=NNEB(ITRI,3)
       end select

        IL2=0
	  DO 31 J=1,NTRI    ! DO 31ѭ�� �������У��ҵ���һ��������
!!            if(j.eq.ITRI)goto 31
	   ND22=NNEB(J,3)
	   ND23=NNEB(J,1)
       IF(ND12.EQ.ND23.AND.ND13.EQ.ND22)THEN       !Ӧ���ǽ����Ӧ���,��˳�ζ�Ӧ��ȣ�
                                                   !�������һ�������ε�3����
	   ITRj=j; IL2=3             !�Ķ����Ӧ. ��ITRj��¼�����α��2013-5-13
       ENDIF 

	   ND22=NNEB(J,1)
	   ND23=NNEB(J,2)
        IF(ND12.EQ.ND23.AND.ND13.EQ.ND22)THEN  !�������һ�������ε�1����
	    ITRj=j; IL2=1                                    !�Ķ����Ӧ
       ENDIF 
	
	   ND22=NNEB(J,2)
	   ND23=NNEB(J,3)
       IF(ND12.EQ.ND23.AND.ND13.EQ.ND22)THEN  !�������һ�������ε�2����
	   ITRj=j;  IL2=2                                  !�Ķ����Ӧ
       ENDIF 
31     CONTINUE      
      
!        write(12,*)'ITRIN=',ITRI,ITRIN
       IF(ITRIN.EQ.1) THEN
        NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=J0;NNEB(NTRI,2)=NNEB(ITRI,2)
	                  NNEB(NTRI,3)=NNEB(ITRI,3)
	  NNEB(ITRI,2)=J0   !ԭ�����ε�һ����Ÿı�       
       ENDIF 

       IF(ITRIN.EQ.3) THEN
	  NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=NNEB(ITRI,1);NNEB(NTRI,2)=NNEB(ITRI,2)
	                  NNEB(NTRI,3)=J0
	  NNEB(ITRI,1)=J0   !ԭ�����ε�һ����Ÿı�     
       ENDIF 

       IF(ITRIN.EQ.2) THEN
	  NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=NNEB(ITRI,1);NNEB(NTRI,2)=NNEB(ITRI,2)
	                  NNEB(NTRI,3)=J0 
	  NNEB(ITRI,2)=J0   !ԭ�����ε�һ����Ÿı� 
       ENDIF 
!--     if(IL2.eq.0) goto 31   !-- IL2.eq.0 ��ʾ�߽��ϣ���Ȼ�Ͳ�ִ������3����� 
       IF(IL2.EQ.3)THEN       
	  NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=NNEB(ITRj,1);NNEB(NTRI,2)=NNEB(ITRj,2)
	                  NNEB(NTRI,3)=J0 
	    NNEB(ITRj,1)=J0   !ԭ�����ε�һ����Ÿı�   
       ENDIF

        IF(IL2.EQ.1)THEN       
        NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=J0;NNEB(NTRI,2)=NNEB(ITRj,2)
	                  NNEB(NTRI,3)=NNEB(ITRj,3) 
	    NNEB(ITRj,2)=J0   !ԭ�����ε�һ����Ÿı�
       ENDIF

        IF(IL2.EQ.2)THEN 
	  NTRI=NTRI+1         !���һ��������
	  NNEB(NTRI,1)=NNEB(ITRj,1);NNEB(NTRI,2)=NNEB(ITRj,2)
	                  NNEB(NTRI,3)=J0 
	    NNEB(ITRj,2)=J0   !ԭ�����ε�һ����Ÿı�
        ENDIF 
        
        goto 209
309     continue  
      
      END 
      
      
!###################################
!###################################
       subroutine DXFHEAD(LV1)

        write(LV1,100)
100     format(1x,'0')
         write(LV1,110)
110      format('SECTION')
        write(LV1,120)
120      format(1x,' 2')
        write(LV1,130)
130      format('ENTITIES')
        end


!###################################
!###################################
       subroutine DXFEND(LV1)
         write(LV1,970)
970      format(1x,' 0')
        write(LV1,980)
980      format('ENDSEC')
        write(LV1,970)
        write(LV1,990)
990      format('EOF')
        END
C       ###############################
C       TO DRAW A LINE  
C       ##############################
        SUBROUTINE DLINE(LV,coor,ict,icc)
         real*8  X1,Y1,Z1,X2,Y2,Z2,coor(2,3)
        WRITE(LV,130)
130     FORMAT(2X,'0')
        WRITE(LV,140)
140     FORMAT('LINE')
        WRITE(LV,150)
150     FORMAT('  8')
        WRITE(LV,160)ict    !ͼ��
160     FORMAT(I1) 
	  WRITE(LV,802)
802      format(' 62')
        WRITE(LV,803) icc    !��ɫ
803      format(i5)

        WRITE(LV,170)
170     FORMAT(' 10')
        WRITE(LV,200) coor(1,1)
        WRITE(LV,180)
180     FORMAT(' 20')
        WRITE(LV,200) coor(1,2)
        WRITE(LV,185)
185     FORMAT(' 30')
        WRITE(LV,200) coor(1,3)
        WRITE(LV,190)
190     FORMAT(' 11')
        WRITE(LV,200) coor(2,1)
        WRITE(LV,195)
195     FORMAT(' 21')
        WRITE(LV,200) coor(2,2)
        WRITE(LV,197)
197     FORMAT(' 31')
        WRITE(LV,200) coor(2,3)
200     FORMAT(F15.6)
        RETURN
        END