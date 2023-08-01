
C       #######################################
c       开挖面切割下的3维块体全空间搜索问题
C       
C        目前数组大小基本接近系统内存限制了
C       #######################################
        program main  !!B3DC1(zmd)  

!!      !MS$ATTRIBUTES DLLEXPORT :: B3DC1

         CHARACTER(10) t1z; CHARACTER(8) d1z
         CHARACTER*60  fl1,fl2
           real*8 t0z,t2z
       COMMON /KRE/KRE/fl1/fl1/fl2/fl2/Ncur0/Ncur
       
!!!注意，主程序和被调用的dll程序都用release 64位编译 2018.10.3
       
        !MS$ ATTRIBUTES DLLIMPORT :: Delaunary   
!         itt=0
!         call Delaunary(itt) 
         
        OPEN(1001,FILE='B3dcDATA.dat')

          OPEN(15,FILE='Inpt_B3DC.DAT',STATUS='OLD')
         READ(15,*) fl1
         READ(15,*) fl2
         close (15)         
         
        OPEN(15,FILE=fl1) 
        READ(15,*) MSIM
        
        CALL INPUTDAT
        DO 200 MS=1,MSIM                 !!随机模拟次数
        CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
	  write(1001,310) "*****The",MS,"time simulation analysis*****"
	  write(1001,*) '---Begining Time:', t1z,'  DATE:', d1z,'---'      
        Read(t1z,*)t0z          !字符型转实型
       write(1001,*)'t0z',t0z
         zRE=1
100      KRE=0
        
        CALL JFZOOM
   
        CALL JTSEF1
	  CALL JTSEF2
	  CALL JTSEF3
	  CALL JTSEF4
        if(Ncur.ne.0) CALL JT_P_Q       !-- 2018.3.28，曲面和平面切割程序合二为一, 2018.9.2
        CALL JTCOTOLO
       
	  CALL INTLIN
        CALL LOOPTOPO
        CALL LOOPDELE
	  CALL LOOPNEXT
        CALL BLOCKCUT3D
        
	  CALL BLOCKVOL
	   
	  CALL OTRST1
        CALL OTRST2
	  
  	  IF(KRE.NE.0) THEN
	   zRE=zRE+1
!	   goto 100  !2018.4.20      2021.11.20
	  ELSE
	WRITE(1001,320)  "Block Search Finished! Search Time is",int(zRE)
	 CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*) '---Finished Time:', t1z,'---' 
        Read(t1z,*)t2z    
        write(1001,*)'t2z',t2z
       t2z=t2z-t0z
       write(1001,'(a31,f9.2)') 'Time-consuming in operation(s):', t2z
        ENDIF
200    CONTINUE
	   zmd=zRE
         close(1001)

310      FORMAT(A8,2X,i2,2X,A29)
320      FORMAT(A38,1X,i2)
	   END  


C       ##################################
C       INPUT PARAMETERS FOR JOINTS
C       ###################################
        SUBROUTINE INPUTDAT        
         implicit real*8 (a-h,o-z) 

        DIMENSION CAZOOM(8,3),SCJFS(6,3),SPJFS(6,3)
        DIMENSION DPJFS(6,3),DDJFS(6,3),SENDC(2,6)
        DIMENSION DPSEND(2,3),DDSEND(2,3),SCSEND(2,3),SPSEND(2,3)
	  DIMENSION DJFDE(100,2),CJFDE(100,3),RJFDE(100)
	  DIMENSION P2E(100,5),P1E(100,9),V3E(100,50,3),V2E(100,50,2)
	  DIMENSION P2D(5000,5),P1D(5000,9),V3D(5000,10,3),V2D(5000,10,2)
        DIMENSION KZ1(100),VV(200,3),ML(1000),MLV(1000,50) 
        DIMENSION F(3),A2(3,3),ET(3),EL(3),M5(12,12),M51(12)
	  DIMENSION  M4(100),etv(50,3,200),itv(50,3),ich(50),NC1(99),NL1(99)
	  DIMENSION DJFDE1(5000,2),CJFDE1(5000,8,3),MJFDE1(5000) 
        DIMENSION D1(50,2),C1(50,8,3),M1(50),ixyz(50),xyz(50,2),bc(50)
        COMMON /P2E/P2E/P1E/P1E/V2E/V2E/V3E/V3E
        COMMON/KZ1/KZ1/ITYPE/ITYPE
        COMMON/ENL/N0,NL0/EML/ML,MLV  
        COMMON/M45/NA9,M4,M5,M51/N8/N8/DK/DK2/N10/N9,N10
        COMMON /DJDE/DJFDE/CJDE/CJFDE/RJDE/RJFDE
        COMMON /NJE0/NJFDE0,NJFDE1/NGF0/NGJF0/NJFDE0V/NJFDE0V
        COMMON /CAOM/CAZOOM/DPFS/DPJFS
	  COMMON /DDFS/DDJFS/SCFS/SCJFS/SPFS/SPJFS
        COMMON/OAXIAL/OAXIAL/MJFDE1/MJFDE1
        COMMON /P2D/P2D/P1D/P1D/V3D/V3D/V2D/V2D/Nfadd/Kfadd,Nfadd
	  COMMON/DDSEND/DPSEND,DDSEND,SCSEND,SPSEND,SENDC/NSEND0/NSEND0
        common/Control/Ctl1,Ctl2/ND0/ND0,NC1,NL1/Ncur0/Ncur
	  DD=3.141592654/180
       
        READ(15,*) OAXIAL
       
        DO 30 I=1,8
        READ(15,*) CAZOOM(I,1),CAZOOM(I,2),CAZOOM(I,3)
30      CONTINUE
        READ(15,*) ITYPE        
!-----        
	 Open (31,FILE='~Bkspe3DV1.dat',STATUS='OLD')        !边坡及洞室模型范围轮廓形态 
       Open (32,FILE='~Bkspe3DP1.dat',STATUS='OLD')        !面的法矢均指向块体内部

        read(31,*) N0    !顶点数
       do 10 i=1,N0
        read(31,*) VV(i,1:3)
10     continue
    
        read(32,*) NL0   !有效面数
       do 25 i=1, NL0
           read(32,*) ML(i)
!      do 20 J=1,ML(i)
           read(32,*) MLV(i, 1:ML(i))
!20      continue
25      continue
          CLOSE (31)
          CLOSE (32)    
       
        N9= NL0         !模型轮廓面数，或复杂边坡的表面数
        N7=N0
         if(ITYPE.eq.2) then            !含有洞室时  
       
       Open (31,FILE='~Bkspe3DV2.dat',STATUS='OLD')  ! 洞室顶点       
         read(31,*) ND0       !洞室和排水孔的总数 2018.7.6 
        do i0=1,ND0
        read(31,*) NC1(i0)    !顶点数        
       do 11 i=1,NC1(i0)
        read(31,*) VV(i+N0,1:3)
11     continue
         N0=N0+NC1(i0)     
        enddo
          CLOSE (31)
          
       Open (32,FILE='~Bkspe3DP2.dat',STATUS='OLD')  !   
        N10=0
        read(32,*) ND0        !洞室和排水孔的总数,每个文件前面都重复写上 
         do i0=1,ND0
        read(32,*) NL1(i0)   !有效面数
       do 26 i=1, NL1(i0)
           read(32,*) ML(i+NL0)        
           read(32,*) MLV(i+NL0,1:ML(i+NL0))
           MLV(i+NL0,1:ML(i+NL0))=MLV(i+NL0,1:ML(i+NL0))+N7 
26     continue   
         NL0=NL0+NL1(i0) 
         N10=N10+NL1(i0)       !N10（洞室的面数之和,N9+N10=NL0），不要用户输入
         N7=N7+NC1(i0)
       enddo 
          CLOSE (32)   
        endif  
!        N0=N0+N1  
!        NL0=NL0+NL1       
!---        
28        READ(15,*) N8  !,N10     !! 读入边界面中临空面的数量N8,以进行块体可动性分析                              
                             
!!--- 读入模型轮廓的产状        
	  DO 40 I=1,N9           
         READ(15,*) KZ1(I),Q1,Q2  ! 根据块体形态分析程序，一些有效面可能是由
	                    ! 定义的切割面一分为二得到.不考虑这种情况，仅考虑一一对应
        Q1=Q1*DD
        Q2=Q2*DD
        P2E(i,1)=SIN(Q1)*SIN(Q2)
        P2E(i,2)=SIN(Q1)*COS(Q2)
        P2E(i,3)=COS(Q1)
!        P2E(i,4)=-P2E(i,1)*VV(MLV(I,1),1)-P2E(i,2)*VV(MLV(I,1),2)  
!     #    -P2E(i,3)*VV(MLV(I,1),3)    !2008.4.5 应该是坐标旋转后大地坐标系下的顶点
	 M4(I) = 0
40      continue            
!---
       if(ITYPE.eq.2) then             !含有洞室时  
        N7=N9
        Open (31,FILE='~BkCavernDir.dat')  ! 洞室形态的产状
         read(31,*) ND0       !洞室和排水孔的总数 2018.7.6 
        do i0=1,ND0
        read(31,*) NL1(i0)     !有效面数,重复读入        
       do 12 i=1,NL1(i0)
        read(31,*) KZ1(I+N7),Q1,Q2
         Q1=Q1*DD
         Q2=Q2*DD
        P2E(i+N7,1)=SIN(Q1)*SIN(Q2)
        P2E(i+N7,2)=SIN(Q1)*COS(Q2)
        P2E(i+N7,3)=COS(Q1)
         M4(I+N7) = 0        
12     continue
         N7=N7+NL1(i0)  
        enddo
          CLOSE (31)
       
      endif  
      
!------!!凹形区域   --------     
        READ(15,*) NA9                 
       If (NA9.EQ.0) GOTO 510 
       
       DO 60 i = 1, NA9
       READ(15,*) M5(i, 12), M51(i)
       DO 50 J = 1, M5(i, 12)
        READ(15,*) M5(i, J)
         j1 = M5(i, J)
         M4(j1) = 1
50       continue
60       continue

!--- 开始读入结构面信息-----

510      READ(15,*) NGJF0
          IF (NGJF0.EQ.0) GOTO 515   
        READ(15,*) ((DPJFS(I,J),J=1,3),I=1,NGJF0)
        READ(15,*) ((DDJFS(I,J),J=1,3),I=1,NGJF0)
        READ(15,*) ((SCJFS(I,J),J=1,3),I=1,NGJF0)
        READ(15,*) ((SPJFS(I,J),J=1,3),I=1,NGJF0)
        
515     READ(15,*) NSEND0            !2009.11.1  沉积岩岩层
        IF (NSEND0.EQ.0) GOTO 520
        READ(15,*)((DPSEND(I,J),J=1,3),I=1,NSEND0)
        READ(15,*)((DDSEND(I,J),J=1,3),I=1,NSEND0)
        READ(15,*)((SCSEND(I,J),J=1,3),I=1,NSEND0)
        READ(15,*)((SPSEND(I,J),J=1,3),I=1,NSEND0)   !真间距
        READ(15,*)((SENDC(I,J),J=1,6),I=1,NSEND0)    !起点坐标,终点坐标 2011-7-29

520     READ(15,*) NJFDE0          !!! 圆形定位结构面
!        write(1001,*)'NJFDE0',NJFDE0     
        IF (NJFDE0.EQ.0) GOTO 530
	  do 70 i=1,NJFDE0
         READ(15,*) DJFDE(I,1),DJFDE(I,2),
     #  CJFDE(I,1),CJFDE(I,2),CJFDE(I,3),RJFDE(I)
70      CONTINUE

530     READ(15,*) NJFDE1        !NJFDE1type !! 凸多边形定位结构面,类型（NJFDE1type=1为实际输入，2为洞室等周围自动加密,
        IF (NJFDE1.EQ.0) GOTO 535 !! 以便渗流分析） 2018.1.31.为了输入定位结构面，也输入自动加密，改为两种方式输入 2018.7.23
!!       if(NJFDE1type.eq.1)then
	  do 80 i=1,NJFDE1
         READ(15,*) DJFDE1(I,1),DJFDE1(I,2),MJFDE1(I) !产状/顶点数 
         
	  DO 80 J=1,MJFDE1(I)              !各顶点坐标
        READ(15,*) CJFDE1(I,J,1),CJFDE1(I,j,2),CJFDE1(I,j,3) 
80      CONTINUE 
        
535     READ(15,*) NJFDE2,dmin      !!dmin最小间距，小于面之间的间距小于该值，则不生成面
         M0=0 
         IF (NJFDE2.EQ.0) GOTO 540
        do 95 i=1,NJFDE2       !! 此时NJFDE2为加密循环次数，沿X，Y，Z方向混合放在一起，
                               !! 也可以重复叠放，可简化编程 2018.1.31
          READ(15,*) D1(I,1),D1(I,2),M1(I)  !产状/顶点数 
    
	   DO J=1,M1(I)              !各顶点坐标
           READ(15,*) C1(I,J,1),C1(I,j,2),C1(I,j,3)         
         enddo      
           READ(15,*)ich(i),ixyz(i),xyz(i,1),xyz(i,2),bc(i) ! ixyz(i)=1，2，3为沿X，Y，Z方向
95      CONTINUE                                     ! 相应方向的起点、终点坐标，步长
        
         etv(1:50,1:3,1:200)=0.0         !已有的xyz取值。第一维孔号，第二维xyz，第三维第几个取值
         itv(1:50,1:3)=0                 !已有的xyz取值的个数。第一维孔号，第二维xyz. 2018.7.4        
      
        do 100 i=1,NJFDE2
         k0=int((xyz(i,2)-xyz(i,1))/bc(i))               !加括号 2018.6.20
        DO 96 J=1,k0+1
         
        if(ixyz(i).eq.1)then                    ! 沿X方向
          t1=C1(I,1,1)+real(j-1)*bc(i)    !x取值
          m3=1
          do it=1,itv(ich(i),1)
            if(abs(t1- etv(ich(i),1,it)).lt.dmin) m3=0             
          enddo
           if(m3.eq.0)goto 96
           itv(ich(i),1)=itv(ich(i),1)+1
           etv(ich(i),1,itv(ich(i),1))=t1
           M0=M0+1
           DJFDE1(M0+NJFDE1,1:2)= D1(I,1:2) 
           MJFDE1(M0+NJFDE1)= M1(I)     
           
          CJFDE1(M0+NJFDE1,1:M1(I),1)= C1(I,1:M1(I),1)+real(j-1)*bc(i)
          CJFDE1(M0+NJFDE1,1:M1(I),2:3)=C1(I,1:M1(I),2:3)     
        endif 
        
         if(ixyz(i).eq.2)then                  ! 沿y方向
          t1=C1(I,1,2)+real(j-1)*bc(i)   !y取值
          m3=1
          do it=1,itv(ich(i),2)
            if(abs(t1- etv(ich(i),2,it)).lt.dmin) m3=0             
          enddo
           if(m3.eq.0)goto 96
           itv(ich(i),2)=itv(ich(i),2)+1
           etv(ich(i),2,itv(ich(i),2))=t1
           M0=M0+1
           DJFDE1(M0+NJFDE1,1:2)= D1(I,1:2) 
           MJFDE1(M0+NJFDE1)= M1(I)     
             
          CJFDE1(M0+NJFDE1,1:M1(I),2)=C1(I,1:M1(I),2)+real(j-1)*bc(i)
          CJFDE1(M0+NJFDE1,1:M1(I),1)=C1(I,1:M1(I),1)  
          CJFDE1(M0+NJFDE1,1:M1(I),3)=C1(I,1:M1(I),3)     
         endif 
         
         if(ixyz(i).eq.3)then                ! 沿z方向
           t1=C1(I,1,3)+real(j-1)*bc(i)   !z取值
           m3=1
          do it=1,itv(ich(i),3)
            if(abs(t1- etv(ich(i),3,it)).lt.dmin) m3=0             
          enddo
           if(m3.eq.0)goto 96
           itv(ich(i),3)=itv(ich(i),3)+1
           etv(ich(i),3,itv(ich(i),3))=t1
           M0=M0+1
           DJFDE1(M0+NJFDE1,1:2)= D1(I,1:2) 
           MJFDE1(M0+NJFDE1)= M1(I)          
               
          CJFDE1(M0+NJFDE1,1:M1(I),3)= C1(I,1:M1(I),3)+real(j-1)*bc(i)
          CJFDE1(M0+NJFDE1,1:M1(I),1:2)=C1(I,1:M1(I),1:2)     
           endif         
 
96       CONTINUE        
100      CONTINUE           

540      NJFDE1=NJFDE1+M0               !凸多边形定位结构面数量,指定输入和加密输入属于一种类型  2018.7.23
!         write(1001,*)'NJFDE1',NJFDE1       
!        endif

        
!!!--- 增加虚拟的圆形定位结构面，用于切割较大块体以提高渗流计算精度 2015.5.5
       READ(15,*) NJFDE0V 
!        write(1001,*)'NJFDE0V',NJFDE0V     
        IF (NJFDE0V.EQ.0) GOTO 550
	  do 71 i=NJFDE0+1,NJFDE0+NJFDE0V    !用同样的数组存储，只需在后面识别一下不同结构面类型
         READ(15,*) DJFDE(I,1),DJFDE(I,2),
     #     CJFDE(I,1),CJFDE(I,2),CJFDE(I,3),RJFDE(I)
71      CONTINUE
        
550      READ(15,*) Ctl1,Ctl2  !容差的系数 
          READ(15,*) DK2    !随机数偏移值,使得初始计算结果不同
          READ(15,*) Ncur   !是否考虑曲面，Ncur=0时不考虑，其它为考虑。
                            !可以暂时输入1，后面在GroundSurf输入正确的Ncur（曲面数）
          READ(15,*) KNFAD    !补充曲面作为边界面，1表示考虑,0则不考虑  2021.12.21
          IF(KNFAD.EQ.1)THEN                             
            READ(15,*)Kfadd,Nfadd   !该面的上方或下方的块体保留（0保留上方，1保留下方），补充的曲面序号
          ENDIF
          CLOSE (15)

C  -------对边坡/洞室形体的顶点坐标进行坐标旋转 

        BS = OAXIAL * DD
        DO 200 I=1,NL0
         DO 110 J=1,3
110      F(J)=P2E(I,J)
        CALL PROJTG(F,A2)
        DO 120 J=1,3
        P1E(I,J)=A2(1,J)
        P1E(I,J+3)=A2(2,J)
        P1E(I,J+6)=A2(3,J)
120      CONTINUE

        DO 130 J=1,ML(i)
	   XX1 = VV(MLV(I,J),1) 
         YY1 = VV(MLV(I,J),2)
         X0 = XX1 * Sin(BS) - YY1 * Cos(BS)
         Y0 = XX1 * Cos(BS) + YY1 * Sin(BS)

	  V3E(I,J,1)=X0
        V3E(I,J,2)=Y0
        V3E(I,J,3)=VV(MLV(I,J),3)

	  ET(1)=V3E(I,J,1)
        ET(2)=V3E(I,J,2)
        ET(3)=V3E(I,J,3)
        CALL TLXYZ(A2,ET,EL)
        V2E(I,J,1)=EL(1)
        V2E(I,J,2)=EL(2)
130      CONTINUE
        P2E(i,4)=-P2E(i,1)*V3E(I,1,1)-P2E(i,2)*V3E(I,1,2)
     #      	  -P2E(i,3)*V3E(I,1,3)
200     CONTINUE

!------对凸多边形定位结构面的数据进行准备08.5.16--- 
        DO 400 I=1,NJFDE1
        Q1=DJFDE1(I,1)*DD
        Q2=DJFDE1(I,2)*DD
        P2D(i,1)=SIN(Q1)*SIN(Q2)
        P2D(i,2)=SIN(Q1)*COS(Q2)
        P2D(i,3)=COS(Q1)
         DO 210 J=1,3
210      F(J)=P2D(I,J)
        CALL PROJTG(F,A2)
        DO 220 J=1,3
        P1D(I,J)=A2(1,J)
        P1D(I,J+3)=A2(2,J)
        P1D(I,J+6)=A2(3,J)
220      CONTINUE

        DO 230 J=1,MJFDE1(I)
	   XX1 = CJFDE1(I,J,1)
         YY1 = CJFDE1(I,J,2)
         X0 = XX1 * Sin(BS) - YY1 * Cos(BS)
         Y0 = XX1 * Cos(BS) + YY1 * Sin(BS)

	  ET(1)=X0
        ET(2)=Y0
        ET(3)=CJFDE1(I,J,3)
        CALL TLXYZ(A2,ET,EL)
	  IF (J.EQ.1) Z10=EL(3)    !!第1点为投影面的基准,这样可以使用户
	                           !!输入时不在同一个面上的点位于同一个面上
        V2D(I,J,1)=EL(1)
        V2D(I,J,2)=EL(2)
        EL(3)= Z10
        CALL LTXYZ(A2,ET,EL)

	  V3D(I,J,1)=ET(1)
        V3D(I,J,2)=ET(2)
        V3D(I,J,3)=ET(3)
230      CONTINUE
        P2D(i,4)=-P2D(i,1)*V3D(I,1,1)-P2D(i,2)*V3D(I,1,2)
     #      	  -P2D(i,3)*V3D(I,1,3)

400      CONTINUE
         RETURN
         END 

C       ############################################
C       PRODUCING 3-D JOINT FACES IN SIMULATION ZOOM
C       #######################################
        SUBROUTINE JFZOOM
         implicit real*8 (a-h,o-z) 	 
        DIMENSION CAZOOM(8,3),SCJFS(6,3),SPJFS(6,3)
        DIMENSION DJFDE(100,2),CJFDE(100,3),RJFDE(100)
        DIMENSION DJFRD(90000,2),CJFRD(90000,3),JTYPE(90000) 
        DIMENSION DPJFS(6,3),DDJFS(6,3)
	  DIMENSION RJFRD(90000),NTJT(8),SENDC(2,6)
	   DIMENSION DPSEND(2,3),DDSEND(2,3),SCSEND(2,3),SPSEND(2,3)
        COMMON /CAOM/CAZOOM/SCFS/SCJFS/SPFS/SPJFS/DJDE/DJFDE
        COMMON/CJDE/CJFDE/RJDE/RJFDE/DJRD/DJFRD/CJRD/CJFRD
        COMMON/RJRD/RJFRD/NGF0/NGJF0/NJE0/NJFDE0,NJFDE1/NJFDE0V/NJFDE0V
        COMMON/NJD0/NJFRD0/DPFS/DPJFS/DDFS/DDJFS
	  COMMON/DDSEND/DPSEND,DDSEND,SCSEND,SPSEND,SENDC/NSEND0/NSEND0
	   COMMON /NTJT/NTJT/JTYPE/JTYPE/VRANGE/VRANGE
	   COMMON /ENL/N0,NL0 
         COMMON/OAXIAL/OAXIAL 

        BS = OAXIAL * 3.141592654/180
         DO I=1,8
	    NTJT(I)=0
         ENDDO

        NJF=0;NJF0=0                     !--   2015.5.5  
        X12=CAZOOM(2,1)-CAZOOM(1,1)
        Y14=CAZOOM(4,2)-CAZOOM(1,2)
        Z15=CAZOOM(5,3)-CAZOOM(1,3)

        VBLOC=ABS(X12*Y14*Z15)
        VRANGE=(X12+Y14+Z15)/3.0
        IF (NGJF0.EQ.0) GOTO 25

        DO 20 I=1,NGJF0

        FN1=SCJFS(I,1)
        FA10=SCJFS(I,2)          !----迹长均值
        FB10=SCJFS(I,3)          !----迹长标准差
  
        FA1=2.0/3.1415926*FA10    !----半径均值
        FB1=2.0/3.1415926*FB10    !----半径标准差 = 根号(方差)

	  EL2=FA1**2+FB1**2        !E(l^2)=D(l)+E(l)^2
        
        FA20=SPJFS(I,2)          !----间距均值
        
	if(FN1.lt.1.5)then       !避免负指数分布时，标准差为0
	                              ! 造成的不准确2011-7-18
        GAMAI=1./(2*3.1415926*FA1*FA1*FA20)   !发现结构面数量少一倍时，
      !块体数量会少很多，这样为了达到相似结果，间距要小一倍。2013-6-28
       else  
	  GAMAI=1./(3.1415926*EL2*FA20)
      endif

        NIG=int(VBLOC*GAMAI)
         NTJT(I+1)=NIG                   !  每一组为一种裂隙类型  
!!	    write(*,*)'dd',i,NIG
        DO 20 J=1,NIG

15	  rq=FV(FA1,FB1,FN1,YEL)
         IF (rq.lt.0.01) GOTO 15   !很小的节理重生成 2005.11.13  
	 
100    CALL RANDU1(YFL)
	
        XC1=CAZOOM(1,1)+X12*YFL
	 CALL RANDU1(YFL)
         YC1=CAZOOM(1,2)+Y14*YFL
        CALL RANDU1(YFL)
         ZC1=CAZOOM(1,3)+Z15*YFL
         IF (NJF.EQ.0) GOTO 19
         TC=0.0
         DO 18 K1=1,NJF
         XC2=CJFRD(K1,1)
         YC2=CJFRD(K1,2)
         ZC2=CJFRD(K1,3)

         A1=abs(XC1-XC2)+abs(YC1-YC2)+abs(ZC1-ZC2)
	
         IF (A1.GT.0.05) GOTO 18
         TC=1.0
      	goto 100
18       CONTINUE

19       NJF=NJF+1
         JTYPE(NJF)=I         !-- 结构面类型，随机节理 2015.5.4
        CJFRD(NJF,1)=XC1
        CJFRD(NJF,2)=YC1
        CJFRD(NJF,3)=ZC1

        FN3=DPJFS(I,1)
        FA3=DPJFS(I,2)
        FB3=DPJFS(I,3)
        FN4=DDJFS(I,1)
        FA4=DDJFS(I,2)
        FB4=DDJFS(I,3)

        DJFRD(NJF,1)=FV(FA3,FB3,FN3,YEL)
        DJFRD(NJF,2)=FV(FA4,FB4,FN4,YEL)
        RJFRD(NJF)=rq

CCC----- 分析倾角、倾向可能的分布和大小情况

        IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
	  END IF
        IF (DJFRD(NJF,2).LT.0) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)+360
	  END IF

        IF (DJFRD(NJF,1).GT.90) THEN
          DJFRD(NJF,1)=180- DJFRD(NJF,1)
          DJFRD(NJF,2)=DJFRD(NJF,2)+180
         IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
         END IF
        END IF
        IF (DJFRD(NJF,1).LT.0) THEN
          DJFRD(NJF,1)= -DJFRD(NJF,1)
          DJFRD(NJF,2)=DJFRD(NJF,2)+180
         IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
         END IF
        END IF
20      CONTINUE
         NJF0=NJF
!!---------2009.11.1 沉积岩层的模拟
25     IF (NSEND0.EQ.0) GOTO 30
         M1=0
        DO 28 I=1,NSEND0
	  Q1=DPSEND(i,2)*3.141592654/180
        Q2=DDSEND(i,2)*3.141592654/180
        P2SE1=SIN(Q1)*SIN(Q2)
        P2SE2=SIN(Q1)*COS(Q2)
        P2SE3=COS(Q1)
        DLIN1=SENDC(i,4)-SENDC(i,1) !X坐标的终点-起点
        DLIN2=SENDC(i,5)-SENDC(i,2)
        DLIN3=SENDC(i,6)-SENDC(i,3)
          DLIN4=SQRT(DLIN1**2+DLIN2**2+DLIN3**2)

       coslp=(P2SE1*DLIN1+P2SE2*DLIN2+P2SE3*DLIN3)/DLIN4
         
        FN1=SPSEND(I,1)     !应该服从负指数分布 
	  FA1=SPSEND(i,2)
	  FB1=SPSEND(i,3) 
      
26 	   cq=FV(FA1,FB1,FN1,YEL)        !真间距的随机模拟
         cq=cq/(abs(coslp)+0.0000001)        !视间距 
         cq0=SPSEND(i,2)/(abs(coslp)+0.0000001)     !视间距均值    
       IF (cq.gt.cq0*5.or.cq.lt.cq0*0.3) GOTO 26 
	                               !间距过大或过小时重新生成 2009.11.1  
	  NJF=NJF+1
          JTYPE(NJF)=NGJF0+I         !-- 结构面类型，层面 2015.5.4
        IF(M1.EQ.0)THEN 
	   CJFRD(NJF,1)=SENDC(i,1)+cq*DLIN1/DLIN4
	   CJFRD(NJF,2)=SENDC(i,2)+cq*DLIN2/DLIN4
	   CJFRD(NJF,3)=SENDC(i,3)+cq*DLIN3/DLIN4
	   M1=1
	   TM1= CJFRD(NJF,1)
	   TM2= CJFRD(NJF,2)
	   TM3= CJFRD(NJF,3)
        ELSE
	    CJFRD(NJF,1)=TM1+cq*DLIN1/DLIN4   ! 2011-7-29
		CJFRD(NJF,2)=TM2+cq*DLIN2/DLIN4  
		CJFRD(NJF,3)=TM3+cq*DLIN3/DLIN4
         TM1=CJFRD(NJF,1)	 
         TM2=CJFRD(NJF,2)
         TM3=CJFRD(NJF,3)
	  ENDIF

!!!        CJFRD(NJF,2)=SENDC(i,2)
!!!        CJFRD(NJF,3)=SENDC(i,3)
!        TJFRD(NJF)=NGJF0+I          !2009.11.1 记录组号
!         write(*,*)  CJFRD(NJF,1),CQ
        FN3=DPSEND(I,1)
        FA3=DPSEND(I,2)
        FB3=DPSEND(I,3)
        FN4=DDSEND(I,1)
        FA4=DDSEND(I,2)
        FB4=DDSEND(I,3)

        DJFRD(NJF,1)=FV(FA3,FB3,FN3,YEL)
        DJFRD(NJF,2)=FV(FA4,FB4,FN4,YEL)
        RJFRD(NJF)=2.0/3.1415926*SCSEND(I,2)  !足够大，不需随机模拟
 
CCC----- 分析倾角、倾向可能的分布和大小情况

        IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
	  END IF
        IF (DJFRD(NJF,2).LT.0) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)+360
	  END IF

        IF (DJFRD(NJF,1).GT.90) THEN
          DJFRD(NJF,1)=180- DJFRD(NJF,1)
          DJFRD(NJF,2)=DJFRD(NJF,2)+180
         IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
         END IF
        END IF
        IF (DJFRD(NJF,1).LT.0) THEN
          DJFRD(NJF,1)= -DJFRD(NJF,1)
          DJFRD(NJF,2)=DJFRD(NJF,2)+180
         IF (DJFRD(NJF,2).GT.360) THEN
          DJFRD(NJF,2)=DJFRD(NJF,2)-360
         END IF
        END IF
	  IF(CJFRD(NJF,1).GT.SENDC(i,4)+0.1.OR.CJFRD(NJF,2).GT.  
     #    SENDC(i,5)+0.1.OR.CJFRD(NJF,3).GT.SENDC(i,6)+0.1) then
             NTJT(NGJF0+1+I)=NJF-NJF0   !前面补上NJF0。不然当没有随机节理而直接进行岩层模拟时，
             NJF0=NJF             ! 这里的NJF0没有赋过值，出现忽大忽小的结果  2015.5.5    每一组为一种裂隙类型  
	   GOTO 28
	  endif
	 GOTO 26
28       CONTINUE

30      IF (NJFDE0+NJFDE0V.EQ.0) GOTO 50     ! 2015.5.5
        DO 40 I=1,NJFDE0+NJFDE0V 
        XI=CJFDE(I,1)
        YI=CJFDE(I,2)
        ZI=CJFDE(I,3)

        NJF=NJF+1
    
        CJFRD(NJF,1)=XI
        CJFRD(NJF,2)=YI
        CJFRD(NJF,3)=ZI
        DJFRD(NJF,1)=DJFDE(I,1)
        DJFRD(NJF,2)=DJFDE(I,2)
        RJFRD(NJF)=RJFDE(I)
        IF(I.LE.NJFDE0)THEN
         JTYPE(NJF)=NGJF0+NSEND0+1         !-- 结构面类型，圆形定位结构面 2015.5.4   
        ELSE
         JTYPE(NJF)=10               !-- 结构面类型，虚拟的圆形定位结构面，编号为10  
        ENDIF
40      CONTINUE
50      NJFRD0=NJF   !!圆形结构面数(定位+随机)

 	 write(1001,*) "Stochastic and fixed joint disc Number:", NJFRD0
 	 write(1001,*) "Fixed joint polygon Number:",NJFDE1
 	 write(1001,*) "Slope/cavern surface Number:",NL0

        NTJT(NGJF0+NSEND0+2)=NJFDE0   !注意,若NGJF0=NSEND0=0,则存在NTJT(2)中
                                     !圆形定位结构面和多边形定位结构面各为一种裂隙类型 2015.6.15
	   if(NJFDE0.eq.0) then    
          NTJT(NGJF0+NSEND0+2)=NJFDE1  
	   else
          NTJT(NGJF0+NSEND0+3)=NJFDE1
         endif
         
         if(NJFDE0.eq.0.and.NJFDE1.EQ.0)THEN    !2015.5.5
             NTJT(NGJF0+NSEND0+2)=NJFDE0V 
          ENDIF   
	   if(NJFDE0.eq.0.and.NJFDE1.NE.0.OR.NJFDE0.NE.0.and.NJFDE1.EQ.0)THEN 
             NTJT(NGJF0+NSEND0+3)=NJFDE0V
         endif 
         if(NJFDE0.NE.0.and.NJFDE1.NE.0)THEN     
             NTJT(NGJF0+NSEND0+4)=NJFDE0V
          ENDIF          
         
          NTJT(1)=NGJF0+NSEND0
	   IF(NJFDE0.NE.0) NTJT(1)=NTJT(1)+1 
	   IF(NJFDE1.NE.0) NTJT(1)=NTJT(1)+1 
         IF(NJFDE0V.NE.0) NTJT(1)=NTJT(1)+1       !2015.5.5
	   !-- NTJT(8)表示最多可以有5组随机结构面(包括层面),圆形定位结构面,多边形定位结构面
          !--- 对节理面圆心进行坐标旋转       
	 DO 500  i=1,NJFRD0
            XX1 = CJFRD(i, 1) 
            YY1 = CJFRD(i, 2)
            CJFRD(i, 1) = XX1 * Sin(BS) - YY1 * Cos(BS)
            CJFRD(i, 2) = XX1 * Cos(BS) + YY1 * Sin(BS)
500    CONTINUE
       RETURN
        END

CC       #####################################
        FUNCTION FV(FA,FB,FN,YFL)
CC       #####################################
	    implicit real*8 (a-h,o-z) 
        IF (FN.LT.1.5) THEN ! =1 负指数分布，标准差不起作用,但在前面
	                      ! 计算结构面数量时起作用2011-7-18       
         CALL RANDU1(YFL)
         FV=-FA*LOG(1-YFL)
        ENDIF
	   
        IF (FN.GT.1.5.AND.FN.LT.2.5) THEN   !=2 正态分布 
         R=0.0
         DO 10 I1=1,36
          CALL RANDU1(YFL)
          R=R+YFL
10       CONTINUE
         FV=(0.577*R-10.392)*FB+FA
        ENDIF

        IF (FN.GT.2.5.AND.FN.LT.3.5) THEN  !=3  均匀分布 
          CALL RANDU1(YFL)
          FV=FA+1.732*FB*(2*YFL-1)
        ENDIF

        IF(FN.GT.3.5.AND.FN.LT.4.5) THEN   !=4  对数正态分布 
	   FA0=LOG(FA/SQRT(FB**2/FA**2+1))	   
	   FB0=SQRT(LOG(FB**2/FA**2+1))

         R=0.0
        DO 20 I1=1,36
         CALL RANDU1(YFL)
20       R=R+YFL
         X=(0.577*R-10.392)*FB0+FA0
         FV=EXP(X)
        END IF

        RETURN
        END

C       ############################
        SUBROUTINE RANDU1(YFL)
         implicit real*8 (a-h,o-z) 
         COMMON/DK/DK2
          CALL RANDOM_NUMBER(YFL)
	    YFL=DK2+YFL
	    IF(YFL.GT.1.0)THEN
           YFL=YFL-1.0
          ENDIF
        END

C       ################################
C          圆形节理面之间的迹线分析
C       ################################
        SUBROUTINE JTSEF1
        implicit real*8 (a-h,o-z) 
        DIMENSION DJFRD(90000,2),CJFRD(90000,3),RJFRD(90000)
        DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3)
        DIMENSION A2(4),CIR(6),CP(2,3),E(3),G(3),B2(3,3)
	  DIMENSION A2A(3,3),B22(3,3),FF(3),CIR1(3)
        DIMENSION X(2,2),Y(4,2),ET(3),EL(3)
        DIMENSION P2J(90000,5),P1J(90000,9) 
	  DIMENSION TA2(4),TCIR(6),TCP(2,3)
        COMMON /DJRD/DJFRD/CJRD/CJFRD/RJRD/RJFRD
        COMMON /NJD0/NJFRD0
        COMMON/P1J/P1J/P2J/P2J 
        COMMON/COLTN/COTO,COLO,NOLO,NJT0
        NOLO(:,:)=0
        
       do 10 i=1,NJFRD0 
	  DD=3.1415926/180.0
        Q1=DJFRD(i,1)*DD
        Q2=DJFRD(i,2)*DD
         X0=CJFRD(I,1)
         Y0=CJFRD(I,2)
         Z0=CJFRD(I,3)
        P2J(I,1)=SIN(Q1)*SIN(Q2)
        P2J(I,2)=SIN(Q1)*COS(Q2)
        P2J(I,3)=COS(Q1)
        P2J(I,4)=-P2J(I,1)*X0-P2J(I,2)*Y0-P2J(I,3)*Z0
C-----------NOTE: 节理圆的法向均指向上  
        DO 5 J=1,3
5      FF(J)=P2J(I,J)

        CALL PROJTG(FF,A2A)
        DO 6 J=1,3
        P1J(I,J)=A2A(1,J)
        P1J(I,J+3)=A2A(2,J)
        P1J(I,J+6)=A2A(3,J)
6      CONTINUE
10      continue

       NJT0=0

       DO 300 I1=1,NJFRD0-1

        DO 22 J1=1,4
22      A2(J1)=P2J(I1,J1)
          D2=SQRT(A2(1)**2+A2(2)**2+A2(3)**2)
        DO 24 J1=1,3
        DO 24 J2=1,3
24      B2(J1,J2)=P1J(I1,(J1-1)*3+J2)
	 
        DO 25 J1=1,3
25         E(J1)=CJFRD(I1,J1)
        CALL TLXYZ(B2,E,G)
         
          CIR1(2)=G(1)
          CIR1(3)=G(2)
         CIR1(1)=RJFRD(I1)
ccc----------很小的节理很难形成封闭回路，不进行交线分析  2005.11.13  
	   if(RJFRD(I1).lt.0.25) goto 300   !由0.05改为0.25  2007.12.22

       DO 100 I2=I1+1,NJFRD0
         DO 28 J1=1,3
        DO 28 J2=1,3
28      B22(J1,J2)=P1J(I2,(J1-1)*3+J2)

        DO 30 J1=1,3
30      E(J1)=CJFRD(I2,J1)
        D1=ABS(A2(1)*E(1)+A2(2)*E(2)+A2(3)*E(3)+A2(4))
        D1=D1/D2
	
	   if(RJFRD(I2).lt.0.25) goto 100  
      
	 IF (D1.Gt.RJFRD(I2)-0.001) GOTO 100   ! 0.05改为0.001，需要形成短迹线，然后被删除，且会留下交点位置， 
        CIR(1)=DJFRD(I2,1)                   ! 以便在LOPADDPOINT找回，若没有形成交点，无法找回  2018.1.16
        CIR(2)=DJFRD(I2,2)                   ! 同时存在或同时不存在短迹线时，正确，否则出错。但所采取的容差
        CIR(3)=RJFRD(I2)                     ! 不能保证正确，可能出现在一个面上形成迹线而另一面上没有
        CIR(4)=CJFRD(I2,1)
        CIR(5)=CJFRD(I2,2)
        CIR(6)=CJFRD(I2,3)

        CALL INTPC(A2,CIR,CP,IW)

CCC     为实现对空间中的2个圆面是否相交的判断，圆面需相互之间互换，并用INTPC
ccc       进行相互分析 

         DO 35 K1=1,4
35          TA2(K1)=P2J(I2,K1)
         TCIR(1)=DJFRD(I1,1)
         TCIR(2)=DJFRD(I1,2)
         TCIR(3)=RJFRD(I1)
         TCIR(4)=CJFRD(I1,1)
         TCIR(5)=CJFRD(I1,2)
         TCIR(6)=CJFRD(I1,3)
        CALL INTPC(TA2,TCIR,TCP,IWT)
       
        IF (IW.EQ.0) GOTO 100
        IF (IWT.EQ.0) GOTO 100
        
        !! 圆形定位结构面相交后，交线有问题（或者是随后分析存在问题）。改了半径后，错误避免，
                                             !!原因查了2个小时没查出来，暂且搁置  2015.6.15 
        
!	 IF(ABS(TCIR(4)-CIR(4)).LT.0.05.AND.ABS(TCIR(5)-CIR(5)).LT.0.05
!     # .AND.ABS(TCIR(6)-CIR(6)).LT.0.05) GOTO 39   !圆心位置太近时(人为给定的
          M11=0                                     !定位结构面容易出现)
	    M12=0
	    M13=1
	  DO 36 I=1,2
	D3=(CP(I,1)-TCIR(4))**2+(CP(I,2)-TCIR(5))**2+(CP(I,3)-TCIR(6))**2
	  D3=SQRT(D3)
 
        IF (D3.Gt.TCIR(3)-0.001) THEN       ! -0.01
	     M11=M11+1
         ENDIF
36       CONTINUE
         IF (M11.LT.2) THEN
           M13=0
        ENDIF
 
	  DO 37 I=1,2
	D3=(TCP(I,1)-CIR(4))**2+(TCP(I,2)-CIR(5))**2+(TCP(I,3)-CIR(6))**2
	  D3=SQRT(D3)
 
         IF (D3.Gt.CIR(3)-0.001) THEN       ! 发现当两个圆的大小相同，交点靠近边缘时，不能形成交线
                                           ! 0.1改为0.001.  2018.1.3
	     M12=M12+1
         ENDIF
37       CONTINUE
         IF (M12.LT.2) THEN
           M13=0
         ENDIF
         
          q1=(CP(1,1)+ CP(2,1))/2.0      !两个圆之间卡在一起，交点都在圆周上
          q2=(CP(1,2)+ CP(2,2))/2.0      !这种情况要保证交线存在 2018.1.16
          q3=(CP(1,3)+ CP(2,3))/2.0          
         D3=(q1-TCIR(4))**2+(q2-TCIR(5))**2+(q3-TCIR(6))**2  !连线的中点与圆心的距离
          D3=SQRT(D3) 
        IF (D3.lt.TCIR(3)-0.001) M13=0 
  
       IF (M13.EQ.1)  GOTO 100
ccc      note: 当2个圆面位于两个面之交棱上的点(分别有2个),其相互之间均位于
ccc       另一圆面交点之外侧时，应排除  

39	  DO 40 J1=1,3
40      E(J1)=CP(1,J1)
        CALL TLXYZ(B2,E,G)
        Z10=G(3)

        X(1,1)=G(1)
        X(1,2)=G(2)
        DO 50 J1=1,3
50      E(J1)=CP(2,J1)
        CALL TLXYZ(B2,E,G)
        X(2,1)=G(1)
        X(2,2)=G(2)
        X11=X(1,1)
        X12=X(1,2)
        X21=X(2,1)
        X22=X(2,2)
         
          DA1=(X11-cir1(2))**2+(X12-cir1(3))**2  
            DA1=sqrt(DA1)
	    DA0=DA1-RJFRD(I1)
      
	  IF (DA0.LT.0.000001) THEN
	     IW1=0
        else
	     IW1=1   
	  ENDIF
         DA1=(X21-cir1(2))**2+(X22-cir1(3))**2  
          DA1=sqrt(DA1)
	   DA0=DA1-RJFRD(I1)
	  IF (DA0.LT.0.000001) THEN
	     IW2=0
        else
	     IW2=1    
        ENDIF
!           if(i1.eq.545.and.i2.eq.556)then
!            write(1001,*)'id',IW1,IW2
!            write(1001,*) CIR1 
!          write(1001,*) X(1,1:2),X(2,1:2) 
!         endif    
        IF (IW1.EQ.1.AND.IW2.EQ.1) THEN      !! 2个点都位于圆外
CC     由于X( , )已是在 I1节理的局部坐标系下的坐标，因此只需输入
cc         半径 CIR1，就可实现迹线与圆周的交点分析
cc      在INTPC中有了IW的判断，因此该线段必定穿过圆，存在2个交点
        CALL INTSB1(CIR1,X,Y,NP)
         IF (NP.Lt.2) write(1001,*)'wrong1'
          
        NJT0=NJT0+1

        DO 52 K1=1,2
        DO 52 K2=1,2

52      COLO(NJT0,(K1-1)*2+K2)=Y(K1,K2)
        
        DO 58 K1=1,2
        DO 53 K2=1,2
        EL(K2)=Y(K1,K2)
53      CONTINUE
        
	    EL(3)=Z10                   !!注意该句不能放错位置
         CALL LTXYZ(B2,ET,EL)
        DO 54 K3=1,3
          COTO(NJT0,(K1-1)*3+K3)=ET(K3)
54        COTO(NJT0+1,(K1-1)*3+K3)=ET(K3)
       
	    CALL TLXYZ(B22,ET,EL)
        DO 55 K2=1,2
55      COLO(NJT0+1,(K1-1)*2+K2)=EL(K2)
         
58      CONTINUE

        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    !! 由于相邻位置要存储2条迹线，在此再增加1
        ELSE IF (IW1.EQ.0.AND.IW2.EQ.0) THEN
   
	 NJT0=NJT0+1

        DO 60 K1=1,2
        DO 60 K2=1,2
60      COLO(NJT0,(K1-1)*2+K2)=X(K1,K2)
        DO 70 K1=1,2
        DO 70 K2=1,3
         COTO(NJT0,(K1-1)*3+K2)=CP(K1,K2)
70      COTO(NJT0+1,(K1-1)*3+K2)=CP(K1,K2)

        DO 78 K1=1,2
        DO 72 K2=1,3
72         ET(K2)=CP(K1,K2)
	    CALL TLXYZ(B22,ET,EL)
        DO 75 K2=1,2
75       COLO(NJT0+1,(K1-1)*2+K2)=EL(K2)
78      CONTINUE

        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    !! 由于相邻位置要存储2条迹线，在此再增加1
        ELSE IF (IW1.EQ.0.AND.IW2.EQ.1) THEN    

	 NJT0=NJT0+1

        COLO(NJT0,1)=X(1,1)
        COLO(NJT0,2)=X(1,2)
 
	 CALL INTSB1(CIR1,X,Y,NP)
        IF (NP.ne.1) then
          write(1001,*)'wrong2' !没有出现这种错误，INTSB1没什么问题
          write(1001,*) CIR1,NP  
          write(1001,*) X(1,1:2),X(2,1:2)
          write(1001,*)Y            
        endif
        COLO(NJT0,3)=Y(1,1)
        COLO(NJT0,4)=Y(1,2)
        DO 80 K1=1,3
        COTO(NJT0,K1)=CP(1,K1)
80      COTO(NJT0+1,K1)=CP(1,K1)
        DO 82 K1=1,2
82      EL(K1)=Y(1,K1)
        EL(3)=Z10
        CALL LTXYZ(B2,ET,EL)
        DO 84 K1=1,3
          COTO(NJT0,K1+3)=ET(K1)
84      COTO(NJT0+1,K1+3)=ET(K1)

        CALL TLXYZ(B22,ET,EL)
        DO 85 K1=1,2
85       COLO(NJT0+1,2+K1)=EL(K1)
        DO 86 K1=1,3
86        ET(K1)=CP(1,K1)
        CALL TLXYZ(B22,ET,EL)
        DO 87 K1=1,2
87       COLO(NJT0+1,K1)=EL(K1)
       
        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    !! 由于相邻位置要存储2条迹线，在此再增加1
        ELSE
	 NJT0=NJT0+1

        COLO(NJT0,3)=X(2,1)
        COLO(NJT0,4)=X(2,2)
       
	 CALL INTSB1(CIR1,X,Y,NP)
         IF (NP.ne.1) write(1001,*)'wrong3'
        COLO(NJT0,1)=Y(1,1)
        COLO(NJT0,2)=Y(1,2)
        DO 90 K1=1,2
90       EL(K1)=Y(1,K1)
        EL(3)=Z10
        CALL LTXYZ(B2,ET,EL)
        DO 92 K1=1,3
        COTO(NJT0,K1+3)=CP(2,K1)
        COTO(NJT0,K1)=ET(K1)
        COTO(NJT0+1,K1+3)=CP(2,K1)
        COTO(NJT0+1,K1)=ET(K1)
92      CONTINUE

        CALL TLXYZ(B22,ET,EL)
        DO 95 K1=1,2
95       COLO(NJT0+1,K1)=EL(K1)
        DO 96 K1=1,3
96        ET(K1)=CP(2,K1)
        CALL TLXYZ(B22,ET,EL)
        DO 97 K1=1,2
97       COLO(NJT0+1,K1+2)=EL(K1)

        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
	   NJT0=NJT0+1 
        END IF
       
100     CONTINUE       
300     CONTINUE
        RETURN
        END

C       ################################
C       多边形定位结构面及临空面与圆形结构面之间形成的迹线分析
C       ################################
        SUBROUTINE JTSEF2
        implicit real*8 (a-h,o-z) 
        DIMENSION DJFRD(90000,2),CJFRD(90000,3),RJFRD(90000)
	  DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3)
        DIMENSION A2(4),CIR(6),CP(2,3),E(3),G(3),B2(3,3)
        DIMENSION X(2,2),EP(50,3),Y(41,2),ET(3),EL(3)
        DIMENSION P2J(90000,5),P1J(90000,9)
	  DIMENSION P2E(100,5),P1E(100,9),V2E(100,50,2),CVEF(800,2)
        DIMENSION B22(3,3),ML(1000),MJFDE1(5000)  
	  DIMENSION P2D(5000,5),P1D(5000,9),V2D(5000,10,2)
        DIMENSION XL(41,2),XL2(41,2),XT(41,3)              
         COMMON /DJRD/DJFRD/CJRD/CJFRD/RJRD/RJFRD
        COMMON /NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1
        COMMON/P1J/P1J/P2J/P2J/MJFDE1/MJFDE1
        COMMON /P2E/P2E/P1E/P1E/V2E/V2E
        COMMON/ENL/N0,NL0/EML/ML,MLV 
        COMMON/COLTN/COTO,COLO,NOLO,NJT0
        COMMON /P2D/P2D/P1D/P1D/V2D/V2D 

        DO 300 I1=NJFRD0+1,NJFRD0+NJFDE1+NL0 
         
        IF(I1.LE.NJFRD0+NJFDE1)THEN   !! I1属于凸多边形定位结构面
        
	  NVEF=MJFDE1(I1-NJFRD0)  
       DO 2 K0=1,NVEF
        CVEF(K0,1)=V2D(I1-NJFRD0,K0,1)
2      CVEF(K0,2)=V2D(I1-NJFRD0,K0,2)

        DO 4 J1=1,4
4      A2(J1)=P2D(I1-NJFRD0,J1)
        D2=SQRT(A2(1)**2+A2(2)**2+A2(3)**2)

        DO 6 J1=1,3
        DO 6 J2=1,3
6      B2(J1,J2)=P1D(I1-NJFRD0,(J1-1)*3+J2)

        DO 8 J1=1,NVEF
        J2=J1+1
        IF (J2.GT.NVEF) J2=J2-NVEF
        EP(J1,1)=CVEF(J1,2)-CVEF(J2,2)
        EP(J1,2)=CVEF(J2,1)-CVEF(J1,1)
        EP(J1,3)=CVEF(J1,1)*CVEF(J2,2)-CVEF(J2,1)*CVEF(J1,2)
8      CONTINUE

	  ENDIF

        IF(I1.GT.NJFRD0+NJFDE1.AND.I1.LE.NJFRD0+NJFDE1+NL0)THEN
	                               !! I1属于多边形开挖面
	  NVEF=ML(I1-NJFRD0-NJFDE1)  
       DO 18 K0=1,NVEF
        CVEF(K0,1)=V2E(I1-NJFRD0-NJFDE1,K0,1)
18      CVEF(K0,2)=V2E(I1-NJFRD0-NJFDE1,K0,2)

        DO 22 J1=1,4
22      A2(J1)=P2E(I1-NJFRD0-NJFDE1,J1)
        D2=SQRT(A2(1)**2+A2(2)**2+A2(3)**2)

        DO 24 J1=1,3
        DO 24 J2=1,3
24      B2(J1,J2)=P1E(I1-NJFRD0-NJFDE1,(J1-1)*3+J2)

        DO 26 J1=1,NVEF
        J2=J1+1
        IF (J2.GT.NVEF) J2=J2-NVEF
        EP(J1,1)=CVEF(J1,2)-CVEF(J2,2)
        EP(J1,2)=CVEF(J2,1)-CVEF(J1,1)
        EP(J1,3)=CVEF(J1,1)*CVEF(J2,2)-CVEF(J2,1)*CVEF(J1,2)
26      CONTINUE

        ENDIF

        DO 200 I2=1,NJFRD0

        DO 28 J1=1,3
        DO 28 J2=1,3
28      B22(J1,J2)=P1J(I2,(J1-1)*3+J2)

        DO 30 J1=1,3
30      E(J1)=CJFRD(I2,J1)
        D1=ABS(A2(1)*E(1)+A2(2)*E(2)+A2(3)*E(3)+A2(4))
        D1=D1/D2
       IF (D1.GE.RJFRD(I2)) GOTO 200
        CIR(1)=DJFRD(I2,1)
        CIR(2)=DJFRD(I2,2)
        CIR(3)=RJFRD(I2)
        CIR(4)=CJFRD(I2,1)
        CIR(5)=CJFRD(I2,2)
        CIR(6)=CJFRD(I2,3)
        CALL INTPC(A2,CIR,CP,IW)   !!

        IF (IW.EQ.0) GOTO 200
         
        DO 40 J1=1,3
40      E(J1)=CP(1,J1)
        CALL TLXYZ(B2,E,G)
        Z10=G(3)
        X(1,1)=G(1)
        X(1,2)=G(2)
        DO 50 J1=1,3
50      E(J1)=CP(2,J1)
        CALL TLXYZ(B2,E,G)
        X(2,1)=G(1)
        X(2,2)=G(2)
        X11=X(1,1)
        X12=X(1,2)
        X21=X(2,1)
        X22=X(2,2)

        CALL POINTC(X11,X12,NVEF,CVEF,IW1)
        CALL POINTC(X21,X22,NVEF,CVEF,IW2)
        if(i2.eq.558.and.I1.eq.NJFRD0+3)then
!              write(1001,*)'1d', IW1,IW2
              endif
        IF (IW1.EQ.1.AND.IW2.EQ.1) THEN      !! 2个点都位于域外
         
	   CALL INTSB(NVEF,CVEF,EP,X,Y,NP)
          IF(NP.GT.2) WRITE(1001,*)'NP=',NP,I1,I2
 	   IF(NP.LT.2) GOTO 200   
            
         M0=0          
301      M0=M0+1   
         M1=M0*2-1; M2=M0*2     !  2015.8.21 
         
         NJT0=NJT0+1

        DO 51 K2=1,2
51      COLO(NJT0,K2)=Y(M1,K2)
        DO 52 K2=1,2
52      COLO(NJT0,2+K2)=Y(M2,K2)
         if(np.gt.2)write(1001,*)Y(M1,:),Y(M2,:)
         DO 53 K2=1,2
        EL(K2)=Y(M1,K2)       !---交点1, 由1改为M1  2015.8.21
53      CONTINUE

	    EL(3)=Z10
         CALL LTXYZ(B2,ET,EL)
        DO 54 K3=1,3
          COTO(NJT0,K3)=ET(K3)
54        COTO(NJT0+1,K3)=ET(K3)
	   CALL TLXYZ(B22,ET,EL)
        DO 55 K2=1,2
55      COLO(NJT0+1,K2)=EL(K2)
         
        DO 56 K2=1,2
        EL(K2)=Y(M2,K2)       !!----  交点2,取编号最大的交点.因为若临空面为凹形面时,存在多个交点。
                       ! 做锦屏的多台阶边坡时发现，“取编号最大的交点”的做法是不对的， 需要严格存为多条线段 2015.8.21
56      CONTINUE
	    EL(3)=Z10
         CALL LTXYZ(B2,ET,EL)
        DO 57 K3=1,3
	    COTO(NJT0,3+K3)=ET(K3)
57        COTO(NJT0+1,3+K3)=ET(K3)
                                  !!可能因定位结构面位置不好，发生POINTC，INTSB不匹配，计算得到交点为空，出现NaN
	   CALL TLXYZ(B22,ET,EL)
        DO 58 K2=1,2
58      COLO(NJT0+1,2+K2)=EL(K2)

        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
	   NJT0=NJT0+1    !! 由于相邻位置要存储2条迹线，在此再增加1
         
         IF(M0.LT.NP/2)GOTO 301     !  2015.8.21 
         
        ELSE IF (IW1.EQ.0.AND.IW2.EQ.0) THEN
!!!!!   将临空面看作无限大，圆与之的交点都在其内部，这样再通过统一的求交点运算，
!!!!!    可以得到节理迹线与边界线的交点，并得到交点的相对位置T1，T2,以便统一进行
!!!!!   回路分析   
!!!!!!   但对于节理圆而言,这样分析夸大了临空面在其上形成的迹线,因此对于节理圆还是分
!!!!!!   4种情况进行处理。且总是存最外侧的交点  

        NJT0=NJT0+1

        DO 60 K1=1,2
        DO 60 K2=1,2
60      COLO(NJT0,(K1-1)*2+K2)=X(K1,K2)
        DO 70 K1=1,2
        DO 70 K2=1,3
         COTO(NJT0,(K1-1)*3+K2)=CP(K1,K2)
70      COTO(NJT0+1,(K1-1)*3+K2)=CP(K1,K2)    !!其实该情况也可能存在交点

        DO 78 K1=1,2
        DO 72 K2=1,3
72       ET(K2)=CP(K1,K2)
	    CALL TLXYZ(B22,ET,EL)
        DO 75 K2=1,2
75       COLO(NJT0+1,(K1-1)*2+K2)=EL(K2)
78      CONTINUE
  
        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1 
	     
        ELSE    !IW1.EQ.0.AND.IW2.EQ.1和IW1.EQ.1.AND.IW2.EQ.0 两种情况合并. 2015.8.27

         CALL INTSB(NVEF,CVEF,EP,X,Y,NP)   
!!       IF (NP.GE.3) WRITE(1001,*)'NP=',NP,I1,I2
         IF (IW1.EQ.0.AND.IW2.EQ.1) THEN   !!这里也要考虑交点为1,3,5,7...个时，按交点排序存为多条线段
           XL(1,1:2)=X(1,1:2)         !  2015.8.27
           DO IT=1,NP
            XL(IT+1,1:2)=Y(IT,1:2) 
           ENDDO    
         ENDIF
         IF (IW1.EQ.1.AND.IW2.EQ.0) THEN 
           XL(1,1:2)=X(2,1:2)         !  2015.8.27
           DO IT=1,NP
            XL(IT+1,1:2)=Y(NP-IT+1,1:2)  !交点的排序从X1算起，这里从X2开始存储，Y()要倒过来

           ENDDO   
         ENDIF
       
         EL(3)=Z10
        DO 401 K1=1,NP+1
        EL(1:2)=XL(K1,1:2)    
           if(i2.eq.558.and.I1.eq.NJFRD0+3)then
              write(1001,*)'1E',  EL(1:2),np
              endif
        CALL LTXYZ(B2,ET,EL)
        XT(K1,1:3)=ET(1:3)  
401     CONTINUE
        DO 402 K1=1,NP+1
        ET(1:3)=XT(K1,1:3)
        CALL TLXYZ(B22,ET,EL)
        XL2(K1,1:2)=EL(1:2)  
402     CONTINUE 
        
          M0=0          
302      M0=M0+1   
         M1=M0*2-1; M2=M0*2      
                 
         NJT0=NJT0+1
        COLO(NJT0,1:2)=XL(M1,1:2)         
        COLO(NJT0,3:4)=XL(M2,1:2)     
        COLO(NJT0+1,1:2)=XL2(M1,1:2)         
        COLO(NJT0+1,3:4)=XL2(M2,1:2)     
        
        COTO(NJT0,1:3)=XT(M1,1:3)  
        COTO(NJT0,4:6)=XT(M2,1:3)
        COTO(NJT0+1,1:3)=XT(M1,1:3)
        COTO(NJT0+1,4:6)=XT(M2,1:3)       
       
        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    
        IF(M0.LT.(NP+1)/2) GOTO 302    
          
        END IF

200     CONTINUE
300     CONTINUE

       RETURN
        END

C       ################################
C       多边形定位结构面之间,多边形定位结构面与临空面形成的迹线分析
C       以定位结构面及临空面为母,定位结构面与之相交
C       ################################
        SUBROUTINE JTSEF3
	  implicit real*8 (a-h,o-z) 
        DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3)
        DIMENSION A2(4),CP(10,3),E(3),G(3),B2(3,3),ML(1000)
        DIMENSION X(2,2),EP(50,3),Y(41,2),ET(3),EL(3),CVEF(800,2)
	  DIMENSION P2E(100,5),P1E(100,9),V2E(100,50,2)
        DIMENSION B22(3,3),MJFDE1(5000),V2D(5000,10,2),MLV(1000,50)  
	  DIMENSION P2D(5000,5),P1D(5000,9),V3D(5000,10,3)
         DIMENSION XL(41,2),XL2(41,2),XT(41,3),V3E(100,50,3)
        COMMON /NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1
        COMMON/MJFDE1/MJFDE1/P2E/P2E/P1E/P1E
        COMMON/ENL/N0,NL0/COLTN/COTO,COLO,NOLO,NJT0/N10/N9,N10
        COMMON /P2D/P2D/P1D/P1D/V3D/V3D/V2E/V2E/V2D/V2D
        COMMON/V3E/V3E/EML/ML,MLV 

        DO 300 I1=NJFRD0+1,NJFRD0+NJFDE1+NL0

        IF(I1.LE.NJFRD0+NJFDE1)THEN   !! I1属于凸多边形定位结构面
        
	  NVEF=MJFDE1(I1-NJFRD0)  
       DO 2 K0=1,NVEF
        CVEF(K0,1)=V2D(I1-NJFRD0,K0,1)
2      CVEF(K0,2)=V2D(I1-NJFRD0,K0,2)

        DO 4 J1=1,4
4       A2(J1)=P2D(I1-NJFRD0,J1)
        D2=SQRT(A2(1)**2+A2(2)**2+A2(3)**2)

        DO 6 J1=1,3
        DO 6 J2=1,3
6      B2(J1,J2)=P1D(I1-NJFRD0,(J1-1)*3+J2)

        DO 8 J1=1,NVEF
        J2=J1+1
        IF (J2.GT.NVEF) J2=J2-NVEF
        EP(J1,1)=CVEF(J1,2)-CVEF(J2,2)
        EP(J1,2)=CVEF(J2,1)-CVEF(J1,1)
        EP(J1,3)=CVEF(J1,1)*CVEF(J2,2)-CVEF(J2,1)*CVEF(J1,2)
8      CONTINUE

	  ENDIF

        IF(I1.GT.NJFRD0+NJFDE1.AND.I1.LE.NJFRD0+NJFDE1+NL0)THEN
	                               !! I1属于多边形开挖面
	  NVEF=ML(I1-NJFRD0-NJFDE1)  
       DO 18 K0=1,NVEF
        CVEF(K0,1)=V2E(I1-NJFRD0-NJFDE1,K0,1)
18      CVEF(K0,2)=V2E(I1-NJFRD0-NJFDE1,K0,2)

        DO 22 J1=1,4
22      A2(J1)=P2E(I1-NJFRD0-NJFDE1,J1)
        D2=SQRT(A2(1)**2+A2(2)**2+A2(3)**2)

        DO 24 J1=1,3
        DO 24 J2=1,3
24      B2(J1,J2)=P1E(I1-NJFRD0-NJFDE1,(J1-1)*3+J2)
        DO 26 J1=1,NVEF
        J2=J1+1
        IF (J2.GT.NVEF) J2=J2-NVEF
        EP(J1,1)=CVEF(J1,2)-CVEF(J2,2)
        EP(J1,2)=CVEF(J2,1)-CVEF(J1,1)
        EP(J1,3)=CVEF(J1,1)*CVEF(J2,2)-CVEF(J2,1)*CVEF(J1,2)
26      CONTINUE

        ENDIF

        IF(I1.LE.NJFRD0+NJFDE1) THEN   !2009.3.31
	   I20=I1+1             ! 避免多边形定位结构面重复求交
         I21=NJFRD0+NJFDE1    ! 多边形结构面与多边形结构面求交   
        ELSE   
	   I20=NJFRD0+1         ! I1为临空面时，I2为多边形结构面，也可以是临空面 
         I21=NJFRD0+NJFDE1+NL0 
        ENDIF

        DO 200 I2=I20,I21       ! I2属于多边形定位结构面,现增加NL0，
             !洞室分析时，洞室形态的临空面和模型范围的轮廓面之间要求得交线 2018.1.3   
        if(i1.GT.NJFRD0+NJFDE1.AND.i1.LE.NJFRD0+NJFDE1+NL0-N10.AND.
     #     i2.GT.NJFRD0+NJFDE1.AND.i2.LE.NJFRD0+NJFDE1+NL0-N10)GOTO 200 
                            !同时属于洞室临空面时，不求交，在JTSEF4直接登记交线 
        if(i1.GT.NJFRD0+NJFDE1+NL0-N10.AND.i1.LE.NJFRD0+NJFDE1+NL0.AND.
     #  i2.GT.NJFRD0+NJFDE1+NL0-N10.AND.i2.LE.NJFRD0+NJFDE1+NL0)GOTO 200
                            !同时属于岩体外轮廓面时，不求交，在JTSEF4直接登记交线 
       
        if(I1.GT.NJFRD0+NJFDE1.and.i2.GT.NJFRD0+NJFDE1)then  !加上if后得到定位结构面与临空面求交 2018.1.30 
        if(i1.ge.i2) goto 200   !同为临空面时，当i1取洞室开挖面时，i2取岩体轮廓面
                                !避免i1取岩体轮廓面时，i2取洞室开挖面（与i1取洞室面，
                                ! i2取轮廓面发生重复） 2018.1.3    
       endif 
        
        if(I2.LE.NJFRD0+NJFDE1)THEN      !--1--!
        DO 28 J1=1,3
        DO 28 J2=1,3
28      B22(J1,J2)=P1D(I2-NJFRD0,(J1-1)*3+J2)

         KJ=0
         DO 30 J1=1,MJFDE1(I2-NJFRD0)   !!线段的两个顶点与面的交点 2008.5.18
         X1=V3D(I2-NJFRD0,J1,1)
         Y1=V3D(I2-NJFRD0,J1,2)
         Z1=V3D(I2-NJFRD0,J1,3)
	   IF(J1.LT.MJFDE1(I2-NJFRD0))THEN
         X2=V3D(I2-NJFRD0,J1+1,1)
         Y2=V3D(I2-NJFRD0,J1+1,2)
         Z2=V3D(I2-NJFRD0,J1+1,3)
         ELSE
         X2=V3D(I2-NJFRD0,1,1)
         Y2=V3D(I2-NJFRD0,1,2)
         Z2=V3D(I2-NJFRD0,1,3)
         ENDIF

         T1=A2(1)*X1+A2(2)*Y1+A2(3)*Z1+A2(4)
	   T2=A2(1)*(X1-X2)+A2(2)*(Y1-Y2)+A2(3)*(Z1-Z2)
	   T1=T1/(T2+1.0e-8)          !!求得交点的待定系数,T2可能为0，需要+1.0e-8, 2018.6.29
 	 IF(T1.GT.0.99999.OR.T1.LT.0.00001) GOTO 30  !! 
	                                          
          KJ=KJ+1
 	   CP(KJ,1)=X1+(X2-X1)*T1
 	   CP(KJ,2)=Y1+(Y2-Y1)*T1
 	   CP(KJ,3)=Z1+(Z2-Z1)*T1        
30       CONTINUE
        endif                           !--1--!
            
       if(I2.GT.NJFRD0+NJFDE1)THEN      !--2--!   2018.1.3
        DO 29 J1=1,3
        DO 29 J2=1,3
29      B22(J1,J2)=P1E(I2-NJFRD0-NJFDE1,(J1-1)*3+J2)
         
         KJ=0
         DO 31 J1=1,ML(i2-NJFRD0-NJFDE1)    ! 2018.1.3
             
         X1=V3E(I2-NJFRD0-NJFDE1,J1,1)      
         Y1=V3E(I2-NJFRD0-NJFDE1,J1,2)
         Z1=V3E(I2-NJFRD0-NJFDE1,J1,3)
	   IF(J1.LT.ML(i2-NJFRD0-NJFDE1))THEN
         X2=V3E(I2-NJFRD0-NJFDE1,J1+1,1)
         Y2=V3E(I2-NJFRD0-NJFDE1,J1+1,2)
         Z2=V3E(I2-NJFRD0-NJFDE1,J1+1,3)
         ELSE
         X2=V3E(I2-NJFRD0-NJFDE1,1,1)
         Y2=V3E(I2-NJFRD0-NJFDE1,1,2)
         Z2=V3E(I2-NJFRD0-NJFDE1,1,3)
         ENDIF
       
         T1=A2(1)*X1+A2(2)*Y1+A2(3)*Z1+A2(4)
	   T2=A2(1)*(X1-X2)+A2(2)*(Y1-Y2)+A2(3)*(Z1-Z2)
	   T1=T1/(T2+1.0e-8)          !!求得交点的待定系数
 	 IF(T1.GT.0.99999.OR.T1.LT.0.00001) GOTO 31  !! 
	                                          
          KJ=KJ+1
 	   CP(KJ,1)=X1+(X2-X1)*T1
 	   CP(KJ,2)=Y1+(Y2-Y1)*T1
 	   CP(KJ,3)=Z1+(Z2-Z1)*T1        
31       CONTINUE
        endif                           !--2--!      
            
            
         IF(KJ.GT.2.or.KJ.eq.1)WRITE(1001,*)"Intersection points between
     # polygons > 2 points or =1 point, ERROR!",KJ           !2018.6.17
	   IF(KJ.EQ.0) GOTO 200

        DO 40 J1=1,3
40      E(J1)=CP(1,J1)
        CALL TLXYZ(B2,E,G)
        Z10=G(3)
        X(1,1)=G(1)
        X(1,2)=G(2)
        DO 50 J1=1,3
50      E(J1)=CP(2,J1)
        CALL TLXYZ(B2,E,G)
        X(2,1)=G(1)
        X(2,2)=G(2)
         X11=X(1,1)
        X12=X(1,2)
        X21=X(2,1)
        X22=X(2,2)

        CALL POINTC(X11,X12,NVEF,CVEF,IW1)
        CALL POINTC(X21,X22,NVEF,CVEF,IW2)
!!            write(1001,*)'COTO',i1,i2
!!!!!! 在此之后与JTSEF2相同 !!!! 2009.3.30
            
        IF (IW1.EQ.1.AND.IW2.EQ.1) THEN      !! 2个点都位于域外
 
	   CALL INTSB(NVEF,CVEF,EP,X,Y,NP)
 	   IF (NP.LT.2) GOTO 200
         
         M0=0          
301      M0=M0+1   
         M1=M0*2-1; M2=M0*2     !  2015.8.27

        NJT0=NJT0+1

        DO 51 K2=1,2
51      COLO(NJT0,K2)=Y(M1,K2)
        DO 52 K2=1,2
52      COLO(NJT0,2+K2)=Y(M2,K2)
        
        DO 53 K2=1,2
        EL(K2)=Y(M1,K2)         
53      CONTINUE

         EL(3)=Z10     
         CALL LTXYZ(B2,ET,EL)
        DO 54 K3=1,3
          COTO(NJT0,K3)=ET(K3)
54        COTO(NJT0+1,K3)=ET(K3)

         CALL TLXYZ(B22,ET,EL)
        DO 55 K2=1,2
55      COLO(NJT0+1,K2)=EL(K2)
        
	  DO 56 K2=1,2
        EL(K2)=Y(M2,K2)       
	                       
56      CONTINUE

	    EL(3)=Z10
         CALL LTXYZ(B2,ET,EL)
        DO 57 K3=1,3
          COTO(NJT0,3+K3)=ET(K3)
57        COTO(NJT0+1,3+K3)=ET(K3)
                                  
	   CALL TLXYZ(B22,ET,EL)
        DO 58 K2=1,2
58      COLO(NJT0+1,2+K2)=EL(K2) 
         
        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    
         IF(M0.LT.NP/2)GOTO 301 
         
        ELSE IF (IW1.EQ.0.AND.IW2.EQ.0) THEN
	  NJT0=NJT0+1
        DO 60 K1=1,2
        DO 60 K2=1,2
60      COLO(NJT0,(K1-1)*2+K2)=X(K1,K2)
        DO 70 K1=1,2
        DO 70 K2=1,3
         COTO(NJT0,(K1-1)*3+K2)=CP(K1,K2)
70      COTO(NJT0+1,(K1-1)*3+K2)=CP(K1,K2)

        DO 78 K1=1,2
        DO 72 K2=1,3
72         ET(K2)=CP(K1,K2)
	    CALL TLXYZ(B22,ET,EL)
        DO 75 K2=1,2
75       COLO(NJT0+1,(K1-1)*2+K2)=EL(K2)
78      CONTINUE

        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    
      
        ELSE !IW1.EQ.0.AND.IW2.EQ.1和IW1.EQ.1.AND.IW2.EQ.0 两种情况合并 2015.8.27
        
          CALL INTSB(NVEF,CVEF,EP,X,Y,NP)
         
          IF (IW1.EQ.0.AND.IW2.EQ.1) THEN                  
           XL(1,1:2)=X(1,1:2)         !  2015.8.27
           DO IT=1,NP
            XL(IT+1,1:2)=Y(IT,1:2) 
           ENDDO    
         ENDIF
         IF (IW1.EQ.1.AND.IW2.EQ.0) THEN 
           XL(1,1:2)=X(2,1:2)         !  2015.8.27
           DO IT=1,NP
            XL(IT+1,1:2)=Y(NP-IT+1,1:2)  !交点的排序从X1算起，这里从X2开始存储，Y()要倒过来
           ENDDO   
         ENDIF        
         
          EL(3)=Z10
        DO 401 K1=1,NP+1
        EL(1:2)=XL(K1,1:2)       
        CALL LTXYZ(B2,ET,EL)
        XT(K1,1:3)=ET(1:3)  
401     CONTINUE
        DO 402 K1=1,NP+1
        ET(1:3)=XT(K1,1:3)
        CALL TLXYZ(B22,ET,EL)
        XL2(K1,1:2)=EL(1:2)  
402     CONTINUE         
        
         M0=0          
302      M0=M0+1   
         M1=M0*2-1; M2=M0*2      
                
	   NJT0=NJT0+1
        COLO(NJT0,1:2)=XL(M1,1:2)         
        COLO(NJT0,3:4)=XL(M2,1:2)     
        COLO(NJT0+1,1:2)=XL2(M1,1:2)         
        COLO(NJT0+1,3:4)=XL2(M2,1:2)     
        
        COTO(NJT0,1:3)=XT(M1,1:3)  
        COTO(NJT0,4:6)=XT(M2,1:3)
        COTO(NJT0+1,1:3)=XT(M1,1:3)
        COTO(NJT0+1,4:6)=XT(M2,1:3)     
    
        NOLO(NJT0,1)=I1
        NOLO(NJT0,2)=I2
        NOLO(NJT0+1,1)=I2
        NOLO(NJT0+1,2)=I1
         NJT0=NJT0+1    
       IF(M0.LT.(NP+1)/2) GOTO 302     
     
        END IF
	     
200     CONTINUE
300     CONTINUE

       RETURN
        END

C       ################################
C       临空面间的交线，看作是迹线
C       ################################
        SUBROUTINE JTSEF4
       implicit real*8 (a-h,o-z) 
       DIMENSION V3E(100,50,3),V2E(100,50,2)
        DIMENSION ML(1000),MLV(1000,50)
	  DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3)
        COMMON/V2E/V2E/V3E/V3E/NJD0/NJFRD0
        COMMON/ENL/N0,NL0/EML/ML,MLV/NJE0/NJFDE0,NJFDE1   
        COMMON/COLTN/COTO,COLO,NOLO,NJT0 

        DO 300 I1=1,NL0-1
            
	   DO 200 J=1,ML(i1)
           K1=MLV(i1,J)
	     j1=j
	      IF(J.LT.ML(i1)) THEN
           K2=MLV(i1,J+1)
	     j2=j+1
	       ELSE
           K2=MLV(i1,1)
	     j2=1
             ENDIF

        DO 100 I2=I1+1,NL0

	   DO 90 k=1,ML(i2)
           K3=MLV(i2,k)
	     j3=k
	      IF(k.LT.ML(i2)) THEN
           K4=MLV(i2,k+1)
	     j4=k+1
	       ELSE
           K4=MLV(i2,1)
	     j4=1
             ENDIF

       IF (K1.EQ.K3.AND.K2.EQ.K4.OR.K1.EQ.K4.AND.K2.EQ.K3) THEN
      
           IF (K1.EQ.K4.AND.K2.EQ.K3) THEN   !曲面与模型边界面交棱相交时，在公共交点分析时发现首尾点不对应，    
               j30=j3                      !此时要交换一下 2018.9.29   
               j3=j4
               j4=j30
           endif
           
        NJT0=NJT0+1

       COLO(NJT0,1)=V2E(I1,j1,1)
       COLO(NJT0,2)=V2E(I1,j1,2)
       COLO(NJT0,3)=V2E(I1,j2,1)
       COLO(NJT0,4)=V2E(I1,j2,2)
	 COLO(NJT0+1,1)=V2E(I2,j3,1)
       COLO(NJT0+1,2)=V2E(I2,j3,2)
       COLO(NJT0+1,3)=V2E(I2,j4,1)
       COLO(NJT0+1,4)=V2E(I2,j4,2)

        DO 70 M2=1,3
        COTO(NJT0,M2)=  V3E(I1,j1,M2)
70      COTO(NJT0,M2+3)=V3E(I1,j2,M2)
        DO 80 M2=1,3
        COTO(NJT0+1,M2)=  V3E(I2,j3,M2)
80      COTO(NJT0+1,M2+3)=V3E(I2,j4,M2)

!           if(i1.eq.8-4.and.i2.eq.9-4)then
!          write(1001,*)'COLO1',COLO(NJT0,1:4)
!         write(1001,*) COLO(NJT0+1,1:4)
!             write(1001,*)COTO(NJT0,1:6)
!         write(1001,*) COTO(NJT0+1,1:6)          
!         endif
        NOLO(NJT0,1)=NJFRD0+NJFDE1+I1
        NOLO(NJT0,2)=NJFRD0+NJFDE1+I2
        NOLO(NJT0+1,1)=NJFRD0+NJFDE1+I2
        NOLO(NJT0+1,2)=NJFRD0+NJFDE1+I1
         NJT0=NJT0+1    !! 由于相邻位置要存储2条迹线，在此再增加1
      
	 ENDIF

90      CONTINUE
100     CONTINUE

200     CONTINUE
300     CONTINUE
        
        write(1001,*)'NJT(I0)=',NJT0
       RETURN
        END

C       #######################################
C   将COTO,COLO,NOLO写成COTOJT,COLOJT,NOLOJT
C    COTO等是杂乱无章地存储在一起，而COTOJT则按照节理次序将同一节理的的极限放在
C    相邻次序上，以便通过NJT的累加(如同“指针”)找到每个节理上的迹线  
c   这样可将原来的三维数组改成二维数组，使得数组空间大为减小2008.3.23
C       #######################################

        SUBROUTINE JTCOTOLO
        implicit real*8 (a-h,o-z) 
        DIMENSION NJT(90000),COLOJT(1300000,4)
        DIMENSION COTOJT(1300000,6),NOLOJT(1300000,3)
        DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3)
         COMMON/NJT/NJT/LOJT/COLOJT/NOLOJT/NOLOJT
        COMMON /TOJT/COTOJT/NJE0/NJFDE0,NJFDE1
	  COMMON/P1J/P1J/P2J/P2J/Ncur0/Ncur
        COMMON/ENL/N0,NL0/NJD0/NJFRD0
        COMMON/COLTN/COTO,COLO,NOLO,NJT0 
        CHARACTER(10) t1z; CHARACTER(8) d1z
          
	   Kij=0
	 DO 500 I0=1,NJFRD0+NJFDE1+NL0+Ncur 
       
        NJT(I0)=0

        DO 100 J=1,NJT0   !! 杂乱地存储在一起的全部迹线数

        IF(NOLO(J,1).EQ.I0) THEN   !!该迹线所属的节理面＝I0
        
	  NJT(I0)=NJT(I0)+1
	  
	  DO 30 M=1,4
         COLOJT(Kij+NJT(I0),M)=COLO(J,M)
30      CONTINUE
        DO 40 M=1,6
         COTOJT(Kij+NJT(I0),M)=COTO(J,M)
40      CONTINUE  
           t1=(COTO(J,1)-COTO(J,4))**2+(COTO(J,2)-COTO(J,5))**2
          t1=sqrt(t1+(COTO(J,3)-COTO(J,6))**2)
          if(t1.lt.0.1)then
              write(1001,*)'t1,i0=',t1,i0
              write(1001,*)COTO(J,1:3)
              write(1001,*)COTO(J,4:6)
          endif
        NOLOJT(Kij+NJT(I0),1:2)=NOLO(J,2:3)   !! 记录是哪个节理面在该节理
	                                      !! 上留下的迹线
            !! 因曲交线之间交点的3个面组成，有两个是一样的。通过  
            !! 记录是其中的哪一段，通过不同段，在LOOPDELE识别顶点整体编号
            !! 并且，在INTLIN中避免求交   2018.4.19
        NOLOJT(Kij+NJT(I0),3)=J    ! 记录总迹线号，以便后面移点时获得不同面上的相同迹线号  2021.9.21
	 ENDIF

100    CONTINUE    
        
        
        Kij=Kij+NJT(I0)
         if(I0.eq.18)write(1001,*)'NJT(I0)',Kij  !查迹线编号，以便在B3DCdraw中画出该面的迹线图
         if(I0.eq.19)write(1001,*)'NJT(I0)',Kij  !<调试入口1>
500    CONTINUE
       CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*) " Total joint traces Number=",Kij,
     $	 '             -*-Time:', t1z

        END 

C       #############################
C       INTERSECTION OF A GIVEN PLANE AND A CIRCLE
C       ##    圆与面（无穷大，位置及产状确定）的2个交点
C       ###########################
        SUBROUTINE INTPC(A2,CIR,CP,IW)
        implicit real*8 (a-h,o-z) 
        DIMENSION A2(4),CIR(6),CP(2,3),B2(3,4),C2(3,3),E(3)
             
        DO 10 I=1,4
10      B2(1,I)=A2(I)
        ALF=CIR(1)*3.1415926535/180.0
        BET=CIR(2)*3.1415926535/180.0
        B2(2,1)=SIN(ALF)*SIN(BET)
        B2(2,2)=SIN(ALF)*COS(BET)
        B2(2,3)=COS(ALF)
        B2(2,4)=-B2(2,1)*CIR(4)-B2(2,2)*CIR(5)-B2(2,3)*CIR(6)
        DO 20 K=1,3
        C2(1,K)=B2(1,K)
        C2(2,K)=B2(2,K)
20      CONTINUE
        CALL MUTIVT(C2,C1)
        IF (C1.GT.0.0001) GOTO 30  !较平行时，当作没有交线
        IW=0
        GOTO 110
30      DO 31 K=1,3
31      B2(3,K)=C2(3,K)

        B2(3,4)=-B2(3,1)*CIR(4)-B2(3,2)*CIR(5)-B2(3,3)*CIR(6)

        CALL INTERP(B2,E,IW1)
        IF (IW1.NE.0) GOTO 51
        IW=0
        GOTO 110
51      D1=(E(1)-CIR(4))**2+(E(2)-CIR(5))**2+(E(3)-CIR(6))**2
        D1=SQRT(D1) 
        IF (D1.LE.CIR(3)) GOTO 50           !2008.1.3  2015.8.22,24,26
!!------保证CIR(6)代表的面跨过2个面之间的交棱
!  采用D1.LE.CIR(3)，不设容差。求得线段，哪怕很短。见<容差分析过程>。2015.8.26

        IW=0
        GOTO 110

50      D2=SQRT(CIR(3)**2-D1**2)
        DO 100 K=1,3
        CP(1,K)=E(K)+D2*C2(3,K)
        CP(2,K)=E(K)-D2*C2(3,K)
100     CONTINUE
	   IW=1 
110     CONTINUE
        END       

        
C       #########################
C       INTERSECTION OF circle with a line
C       #########################
        SUBROUTINE INTSB1(cir1,X,Y,NP)
        implicit real*8 (a-h,o-z) 
        DIMENSION X(2,2),Y(4,2),CIR1(3)
   
          x1=X(1,1)-cir1(2)
          y1=X(1,2)-cir1(3)
          x2=X(2,1)-cir1(2)
          y2=X(2,2)-cir1(3)
	     x0=0                 !可以删除x0，y0
	     y0=0
      	r=cir1(1)   ! +0.0001   2015.8.24
          k0=1
            np=0
        if (abs(y1-y2).LT.0.001)  goto 100   ! 不能太小，否则会出错
	                                  ! 该处是关键的出错位置2008.5.20
         ak= y2-y1

         a=x2-x1
	   b=x1*y2-x2*y1         !-ak*x0
	   c=a**2+ak**2
	   d=2*a*b               !-2*ak**2*y0
	   e=b**2-ak**2*r**2     !+ak**2*y0**2
        
	   f1=(d/2/c)**2-e/c
	 if (f1.lt.0) then
	   write(1001,*)'(d/2/c)**2-e/c=',f1
         write(1001,*)x1,y1,x2,y2
         write(1001,*)x0,y0,r
         write(1001,*)ak,a,b,c,d,e
	   write(1001,*)'d**2-4*c*e=',d**2-4*c*e
        
        if (abs(d**2-4*c*e).lt.0.0000001) then
	   yr1=-d/2/c                    
	   yr2=-d/2/c
	 else
	   yr1=(-d+sqrt(d**2-4*c*e))/2/c   ! 采用补救办法,如果还是
	   yr2=(-d-sqrt(d**2-4*c*e))/2/c   ! sqrt()为负,还会出错   
        endif
        goto 90
	 endif

	   yr1=-d/2/c+sqrt(f1)
	   yr2=-d/2/c-sqrt(f1)

90        xr1=a*yr1/ak+(x1*y2-x2*y1)/ak
         xr2=a*yr2/ak+(x1*y2-x2*y1)/ak
     
         GOTO 110

100     yr1=y1
        yr2=y1
	 
	  xr1=sqrt(r**2-(y1-y0)**2)+X0 
        xr2=-sqrt(r**2-(y1-y0)**2) +X0  

       ! goto 1300        !经检查，下面是对的。可能出现求得的点位于(x1,y1),(x2,y2)的中间或两侧，
                           !只有中间的点是待求的 2015.8.25
110       if (abs(y1-y2).GE.0.001) then
           
          d0=SGN((yr1-y1)*(y2-yr1))           
          if (d0.gT.0.00001) then    !位于y1,y2之间  
	    y(1,1)=xr1+cir1(2)
	    y(1,2)=yr1+cir1(3)
          k0=2
 	     np=1
         endif
          d0=SGN((yr2-y1)*(y2-yr2))              
         if (d0.gT.0.00001) then    !位于y1,y2之间,则有交点,保证交点真实存在
	    y(k0,1)=xr2+cir1(2)
	    y(k0,2)=yr2+cir1(3)
 	   if(k0.eq.1) np=1
         if(k0.eq.2) np=2
	   endif
	 else
          d0=SGN((xr1-x1)*(x2-xr1))              
          if (d0.gT.0.00001) then
	    y(1,1)=xr1+cir1(2)
	    y(1,2)=yr1+cir1(3)
          k0=2
 	     np=1
         endif
          d0=SGN((xr2-x1)*(x2-xr2))              
         if (d0.gT.0.00001) then
	    y(k0,1)=xr2+cir1(2)
	    y(k0,2)=yr2+cir1(3)
         if(k0.eq.1) np=1
         if(k0.eq.2) np=2 
	   endif
       endif


       END

C       #########################
C       INTERSECTION OF A SEGMENT AND THE BOUNDARY OF EF
C       #########################
        SUBROUTINE INTSB(NV,CVEF,EP,X,Y,NP)
        implicit real*8 (a-h,o-z) 
        DIMENSION CVEF(800,2),EP(50,3),X(2,2),A(3),Y(41,2),YT(41,2)
        DIMENSION Yt1(41),km(41)
        common/Control/Ctl1,Ctl2 
!子程序调用时,各数组参量的维数要一致.CVEF(800, ),EP(50, )的维数不同,CVEF是为了
!与POINTC中的一致,而POINTC考虑到在回路中应用时,回路的顶点数可能远大于开挖面的
!顶点数.涉及EP(50,3)的数组,通过NVEF给定了长度,只要前后一致,是没问题的
          
        A(1)=X(1,2)-X(2,2)
        A(2)=X(2,1)-X(1,1)
        A(3)=X(1,1)*X(2,2)-X(2,1)*X(1,2)
        NP=0
       ! ctl1=1.0                     !2015.8.21
            
        DO 100 I=1,NV
        I2=I+1
        IF (I2.GT.NV) I2=I2-NV
        D0=EP(I,1)*A(2)-EP(I,2)*A(1)
        IF (ABS(D0).LT.1.0e-4) GOTO 100     !接近平行
        Y1=EP(I,1)*X(1,1)+EP(I,2)*X(1,2)+EP(I,3)       !点到线段的距离
        Y2=EP(I,1)*X(2,1)+EP(I,2)*X(2,2)+EP(I,3)
        Y3=A(1)*CVEF(I,1)+A(2)*CVEF(I,2)+A(3)
        Y4=A(1)*CVEF(I2,1)+A(2)*CVEF(I2,2)+A(3)
        IF (Y1.LT.-1.0e-4*ctl1.AND.Y2.LT.-1.0e-4*ctl1) GOTO 100  !判断点是否位于线段的同侧  
        IF (Y1.GT. 1.0e-4*ctl1.AND.Y2.GT. 1.0e-4*ctl1) GOTO 100  !容差设得宽一点。线段的顶点靠近区域的边界线时，
        IF (Y3.LT.-1.0e-4*ctl1.AND.Y4.LT.-1.0e-4*ctl1) GOTO 100  !形成交点。这样，相邻面上出现的短线段，也形成交点    
        IF (Y3.GT. 1.0e-4*ctl1.AND.Y4.GT. 1.0e-4*ctl1) GOTO 100  !见笔记本中的图示       
!                 <容差分析过程>：
!    当某个圆盘与面相交，且靠近其棱边时，该圆盘可能切穿到相邻的面， 并在相邻的面上形成短线段,该短线段不能忽略。
! 否则，该棱处两侧面上的线段求交后，一个面上的该棱处有短线段，因此该棱被断开，而另一个面没有短线段，该棱不被断开。
! 这样，回路搜索后，该棱处相关的回路上的棱边不对应，回路将被删除而出错。 2015.8.26    
        
        X0=(EP(I,2)*A(3)-EP(I,3)*A(2))/D0
        Y0=(EP(I,3)*A(1)-EP(I,1)*A(3))/D0
         T1=-y1/D0                       !!交点在X(2,2)表示的线上的相对位置,用于交点排序 2015.8.21
       
!        IF (NP.EQ.0) THEN
        NP=NP+1
        Y(NP,1)=X0
        Y(NP,2)=Y0
         yt1(NP)=T1              !! 2015.8.21
!        ELSE
!        D1=1.0
!        DO 10 J1=1,NP
!        D2=(Y(J1,1)-X0)**2+(Y(J1,2)-Y0)**2
!        IF (D2.GT.0.05) GOTO 10  ! 线段要求有一定长度，2014.7.16调试裂隙搭接时的注释
!        D1=0.0
!10      CONTINUE
!        IF (D1.LT.0.5) GOTO 80
!        NP=NP+1             !!交点数量,当为凹形面时,可以有多个交点

!20      Y(NP,1)=X0
!        Y(NP,2)=Y0
!          yt1(NP)=T1        
!80      CONTINUE
!        END IF
100     CONTINUE
     
!        IF (NP.GE.2) THEN   !! 排序,形成端点1,交点1,交点NP-1,交点NP,
!	                      !! 端点2的排列次序,以便直接取得两侧的交点
!         D1=(X(1,1)-Y(NP,1))**2+(X(1,2)-Y(NP,2))**2
!         D2=(X(1,1)-Y(NP-1,1))**2+(X(1,2)-Y(NP-1,2))**2
!        IF(D1.LT.D2) THEN
!         DO 120 I=1,NP
!	   DO 120 J=1,2
!120	   YT(NP-I+1,J)=Y(I,J)
!         DO 130 I=1,NP
!	   DO 130 J=1,2
!130	   Y(I,J)=YT(I,J)          
!         ENDIF
!	  ENDIF 
     
!!   根据T1的大小对点号进行排序  2015.8.21
         if(np.eq.0) return  !当np=0时下面的Y(m0,1:2)有溢出，造成 write(2111,*)'NP=',np这句不能删除，
                             !但为何有这句，在溢出的情况下也能计算呢   2017.6.2
         m1=0
         do I=1,NP
         km(i)=1   
         enddo
         
!        write(2111,*)'NP=',np   !!为什么 不能删除掉？？？或写成其他语句？nmb,写个2111也行
!        下面写了个 write(1001,*)，又可以删除上面一句了，ri
         
501      yt0=2.001
         DO 120 I=1,NP          
        if(yt1(i).lt.yt0.and.km(i).eq.1)then
          yt0=yt1(i)  
          m0=i
        endif 
120      continue        
        m1=m1+1
        YT(m1,1:2)=Y(m0,1:2) 
!      
        km(m0)=0
        if(m1.lt.NP)goto 501
              
          DO 130 I=1,NP           
130	    Y(I,1:2)= YT(I,1:2)
  
         RETURN
        END

C       ##########################
C       INTERSECTION OF LINES
C       ##########################
        SUBROUTINE INTLIN
        implicit real*8 (a-h,o-z) 
        DIMENSION IA(52000,8),Q(52000,8),P(3800,3)
        DIMENSION P2E(100,5),P1E(100,9),ML(1000),JTYPE(90000)  
        DIMENSION ITOTAL(4),NODC(52000,5),NODCZ(3400000,5),B22(3,3)
        DIMENSION B2(3,3),ET(3),EL(3),NNOD(90000),C2(3,3),NNOD0(90000)
        DIMENSION P2J(90000,5),P1J(90000,9),NODP(52000),ET0(3),EL0(3)
        DIMENSION NDLP(88000),NDLPZ(1140000,4),NDLPQ(5500000),QN(3)
        DIMENSION CNODTZ(3400000,3),CNODL(52000,2),CNODT(52000,3)
	  DIMENSION CVEF(50,2),NOLOJT(1300000,3),P2D(5000,5),P1D(5000,9)
	  DIMENSION NJT(90000),COLOJT(1300000,4),COTOJT(1300000,6)
        DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
        DIMENSION NOTJC(52000),coor(2,3)           
        COMMON /NJT/NJT/LOJT/COLOJT/CNODL/CNODL
        COMMON/IA/IA/Q/Q/NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1 
        COMMON /NODP/NODP/CNODT/CNODT/ITOTAL/ITOTAL
        COMMON/TOJT/COTOJT/NOLOJT/NOLOJT/NNOD/NNOD
	  COMMON/P1J/P1J/P2J/P2J/P2E/P2E/P1E/P1E/P2D/P2D/P1D/P1D
        COMMON /NDLP/NDLP/NDLPZ/NDLPZ/NDLPQ/NDLPQ/NPP/NPP0,NPP
        COMMON/CNODTZ/CNODTZ/NODC/NODC/NODCZ/NODCZ/JTYPE/JTYPE
        COMMON/ENL/N0,NL0/EML/ML,MLV/B2T/B2,Z10/Control/Ctl1,Ctl2
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur 
        CHARACTER(10) t1z; CHARACTER(8) d1z

        ne2=-14    !调试的面编号
        
       OPEN(202,FILE='Loop0.dxf')
       call DXFHEAD(202)
         icc=3;ict=3      
             
       DO 6 J=1,NJFDE1 
       do 2 i=1,9 
2      P1J(NJFRD0+J,I)=P1D(J,I)  
       do 4 i=1,4
4      P2J(NJFRD0+J,I)=P2D(J,I)  
6      CONTINUE
        
       DO 16 J=1,NL0 
       do 12 i=1,9 
12      P1J(NJFRD0+NJFDE1+J,I)=P1E(J,I)  
       do 14 i=1,4
14      P2J(NJFRD0+NJFDE1+J,I)=P2E(J,I)  
16     CONTINUE 
       
      DO 26 J=1,Ncur           !2018.4.22
       t1=0.0
       do k0=1,NDcur(j)
       t1=t1+CoCur(j,K0,3)    
       enddo
       t1=t1/real(NDcur(j))   !平均高程      
        
        P2J(NJFRD0+NJFDE1+NL0+J,1)=0.0
        P2J(NJFRD0+NJFDE1+NL0+J,2)=0.0
        P2J(NJFRD0+NJFDE1+NL0+J,3)=1.0
        P2J(NJFRD0+NJFDE1+NL0+J,4)=-t1
         
        EL(1:3)=P2J(NJFRD0+NJFDE1+NL0+J,1:3)
        CALL PROJTG(EL,B2)
        DO 22 I=1,3
        P1J(NJFRD0+NJFDE1+NL0+J,I)=B2(1,I)
        P1J(NJFRD0+NJFDE1+NL0+J,I+3)=B2(2,I)
        P1J(NJFRD0+NJFDE1+NL0+J,I+6)=B2(3,I)
22      CONTINUE
       
26     CONTINUE 
       
        NPP0=0
	  ncnt=0
	  KN=0 
       
       DO 700 IE=1,NJFRD0+NJFDE1+NL0+Ncur      
 
        IA(:,:)=0
        Q(:,:)=0.0  
        NODP(:)=0      
 	  CNODT(:,:)=0.0	 
 	  CNODL(:,:)=0.0	 
 	  NODC(:,:)=0.0 
        NOTJC(:)=0
!------        
	  NT=0
       
        N1=NJT(IE)

        IF (N1.LT.3) THEN
         ITOTAL(1)=0      !2009.8.19
	    GOTO 510       !迹线数少于3条
        ENDIF
	  KIJ=0
	 IF(IE.NE.1) THEN
	  DO 275 IJ=1,IE-1 
275	 KIJ=KIJ+NJT(IJ)
        ENDIF
 
        DO 270 J1=1,3
        DO 270 J2=1,3
270     B2(J1,J2)=P1j(IE,(J1-1)*3+J2)
!         X1=COTOJT(KIJ+1,1)
!         Y1=COTOJT(KIJ+1,2)                 !采用这几行计算得到的Z10，对于曲面
!         Z1=COTOJT(KIJ+1,3)                 !可能因岩体范围外的KIJ+1算得Z10位于岩体范围外， 
!        Z10=X1*B2(3,1)+Y1*B2(3,2)+Z1*B2(3,3)  造成后面PNT_SUR判断得出回路位于区域外而被删除
!                                块体因此而不被曲面切开而出错   2018.4.24
        Z10=-P2J(ie,4)           !曲面时为平均高程
        EL(3)=Z10 
         
        N11=N1-1
        DO 650 I=1,N11
        NJ1=I
        I2=I+1
        DO 650 J=I2,N1
        NJ2=J
!!             write(1001,*)'NOLOJT',KIJ,i,j
         IZ0=NOLOJT(KIJ+I,1)        !2018.4.19
         JZ0=NOLOJT(KIJ+J,1)         
         KZ0=0       
         MZ0=0  
        if(IZ0.EQ.JZ0) then         !曲交线作为迹线，前后两段不进行交点计算，直接存储
          IZ1=NOLOJT(KIJ+i,2)         !前后两段  
          JZ1=NOLOJT(KIJ+J,2) 
          if(IZ1.EQ.JZ1-1)THEN
          T1=1.0                  !发现采用T1和T2得到的交点坐标不同，并造成回路错误且块体没切开，是T1=0和T2=1写反了 2018.9.20
          T2=0.0
          EL(1:2)=COLOJT(KIJ+J,1:2)
          ET(1:3)=COTOJT(KIJ+J,1:3)
         
           KZ0=NOLOJT(KIJ+i,2)    !曲交线时，取I对应的段号  2018.4.19
           MZ0=NOLOJT(KIJ+J,2)    !增加1维，共5维
        !如1段（顶点1和顶点2构成）和2段（顶点2和顶点3构成）相交，交点顺序号为1
!           write(12,*)'sddd',ie,i,j,IZ0,JZ0,IZ1,JZ1
           CALL PJUDG(NT,NJ1,NJ2,IE,ET,EL,T1,T2,IZ0,JZ0,KZ0,MZ0)   
           NOTJC(NT)=KIJ+J   !对NT点保存的坐标对应的某面上的某节理编号进行记录。因发现后续端点移动后，
                             !移动之前已经保存的坐标没有相应修正，造成交点和端点重合而出错，2021.11.24
                             !需统一将之前已保存的端点坐标再更新一下
             GOTO 650   ! 649
          ENDIF
          
         if(IZ1.eq.NJFRD0+NJFDE1+NL0+Ncur+1)then      !环形交线的首尾相连
            
         t3=sqrt((COTOJT(KIJ+i,1)-COTOJT(KIJ+J,4))**2+    !首段的第1点和尾段的第2点距离
     #           (COTOJT(KIJ+i,2)-COTOJT(KIJ+J,5))**2+
     #           (COTOJT(KIJ+i,3)-COTOJT(KIJ+J,6))**2)          
         if(t3.lt.5e-3)then              !交线首尾相连
             write(12,*)'t3=',t3,i,j 
          T1=0.0          !也改过来 2018.9.20
          T2=1.0
          EL(1:2)=COLOJT(KIJ+j,3:4)
          ET(1:3)=COTOJT(KIJ+j,4:6)
           KZ0=NOLOJT(KIJ+j,2)    !取尾段（j段）对应的段号  2018.9.16
                                  ! 1-2,2-3,3-4,...,9-1相连都取前一段的段号，此时前一段是尾段 
           MZ0=NOLOJT(KIJ+i,2)
           CALL PJUDG(NT,NJ1,NJ2,IE,ET,EL,T1,T2,IZ0,JZ0,KZ0,MZ0)    
              GOTO 650  ! 649  
         endif    
      
         endif 
         GOTO 650
         endif 
       goto 9991   !发现两条远离且接近平行的线段得到交点，因此不用下面if，之间按照以前的算法  2018.9.12 
       
!        if( IE.gt.NJFRD0+NJFDE1+NL0)then    !当为曲面上的两条曲线段时，在三维上是A,B两条曲线段相交，
                                            !在X-Y平面上却是A,C(交点靠近端点），从而造成
                                     !交点位置出错，后续回路连接错，块体没被切开 2018.5.18 
          !发现转为局部坐标系下得到的交点精度很差，按石老师建议改为直接求交
            
           X1=COTOJT(KIJ+i,1);Y1=COTOJT(KIJ+i,2);Z1=COTOJT(KIJ+i,3)
           X2=COTOJT(KIJ+i,4);Y2=COTOJT(KIJ+i,5);Z2=COTOJT(KIJ+i,6)
           X3=COTOJT(KIJ+j,1);Y3=COTOJT(KIJ+j,2);Z3=COTOJT(KIJ+j,3)
           X4=COTOJT(KIJ+j,4);Y4=COTOJT(KIJ+j,5);Z4=COTOJT(KIJ+j,6)
     
            X13=X1-X3;Y13=Y1-Y3;Z13=Z1-Z3
            X21=X2-X1;Y21=Y2-Y1;Z21=Z2-Z1
            X43=X4-X3;Y43=Y4-Y3;Z43=Z4-Z3
            A01=X21*X21+Y21*Y21+Z21*Z21
            B01=-(X43*X21+Y43*Y21+Z43*Z21)
            C01=-(X13*X21+Y13*Y21+Z13*Z21)
            A02=-B01
            B02=-(X43*X43+Y43*Y43+Z43*Z43)
            C02=-(X13*X43+Y13*Y43+Z13*Z43)
            D0=A01*B02-A02*B01
            
            if(abs(D0).lt.1.0e-7) D0=D0+1.0e-7
            T1=(C01*B02-C02*B01)/D0
            T2=(A01*C02-A02*C01)/D0
            if(t1.lt.-1.0e-4.or.t1.gt.1.0001)goto 650
            if(t2.lt.-1.0e-4.or.t2.gt.1.0001)goto 650
!          goto  9999      !下面的else可以不要，只要这里的线段求交点算法
!         else          !当为平面时（其上可有曲线段和直线段）  !-----
                                      
9991      ci1=COLOJT(KIJ+i,1)
          ci2=COLOJT(KIJ+i,2)
          ci3=COLOJT(KIJ+i,3)
          ci4=COLOJT(KIJ+i,4)
           cj1=COLOJT(KIJ+J,1)
          cj2=COLOJT(KIJ+J,2)
          cj3=COLOJT(KIJ+J,3)
          cj4=COLOJT(KIJ+J,4)
         if(i.eq.235.and.ie.eq.10)Then
!            write(1001,*)'ii',COTOJT(KIJ+i,1:3)
!            write(1001,*)COTOJT(KIJ+i,4:6)
         endif 
          if(j.eq.510.and.ie.eq.10.and. KIJ.eq.2408)Then
!            write(1001,*)'jj',COTOJT(KIJ+j,1:3)
!            write(1001,*)COTOJT(KIJ+j,4:6)
         endif 
         X1=CI1
         Y1=CI2
         X2=CI3
         Y2=CI4
         P(I,1)=Y1-Y2
         P(I,2)=X2-X1
         P(I,3)=X1*Y2-X2*Y1  
         X3=CJ1
         Y3=CJ2
         X4=CJ3
         Y4=CJ4
         P(J,1)=Y3-Y4
         P(J,2)=X4-X3
         P(J,3)=X3*Y4-X4*Y3
          
        D0=P(I,1)*P(J,2)-P(I,2)*P(J,1) 
!!	 下句为  <易出错的小量控制> <容差分析过程>
!        IF (ABS(D0).LT.1.0e-5) GOTO 650 !线段接近平行。注意门槛值，
	                            !! 太小会引起交点，远离准确位置 2008.6.3  
        SL1=P(I,1)*CJ1+P(I,2)*CJ2+P(I,3)      
        SL2=P(I,1)*CJ3+P(I,2)*CJ4+P(I,3)
        SL3=P(J,1)*CI1+P(J,2)*CI2+P(J,3)
        SL4=P(J,1)*CI3+P(J,2)*CI4+P(J,3)
            if(i.eq.235.and.j.eq.510.and.ie.eq.10)Then
!                 write(1001,*)'k1',COLOJT(KIJ+i,1:4)
!                 write(1001,*)'k1',COTOJT(KIJ+i,1:6)
               
!                 write(1001,*)'k2',COLOJT(KIJ+j,1:4)
!                 write(1001,*)'k2',COTOJT(KIJ+j,1:6)
                  
!                write(1001,*)'kk',SL1,SL2,SL3,SL4
                endif
!!	 下面为  <易出错的小量控制>从0.05/0.01/0.005都用过       
      IF (SL1.LT.-1.0e-3*Ctl2.AND.SL2.LT.-1.0e-3*Ctl2) GOTO 650   
      IF (SL1.GT. 1.0e-3*Ctl2.AND.SL2.GT. 1.0e-3*Ctl2) GOTO 650   
      IF (SL3.LT.-1.0e-3*Ctl2.AND.SL4.LT.-1.0e-3*Ctl2) GOTO 650  
      IF (SL3.GT. 1.0e-3*Ctl2.AND.SL4.GT. 1.0e-3*Ctl2) GOTO 650   
      !!注意门槛值不能太小，因为其中有临空面的相邻棱线求交点,其交点应该得到
!!        再次将阈值定为0.005，保证临空面在圆上形成的相邻棱线之间的交点能够求得，从而
!!       能够找到正确的回路。保证更大结构面切割得到的可动块体总是包括短小结构面形成的可动块体
!!      另外，本次纠错进一步说明回路删除后的程序没什么问题,2008.9.10.  由0.001改为0.01,2018.10.14
!!!! ------------------
!!  共有3种交点：（1）结构面迹线与结构面迹线、（2）结构面迹线与临空面交线、
!!        （3）临空面交线与临空面交线
!!  其中，（1）过去就遇到过；（2）通过结构面在临空面上的迹线向临空面范围外延伸，这样
!!       就可以很容易获得交点及T1，T2；（3）临空面交线之间的交点，当为相邻的交线时，
!!     计算得到T1，T2等于0或1，交点即为其中的重合顶点（会不会因计算误差而不能完全重合?）
!!!! ------------------
             
!        EL(1)=(P(I,2)*P(J,3)-P(I,3)*P(J,2))/D0
!        EL(2)=(P(I,3)*P(J,1)-P(I,1)*P(J,3))/D0
!        EL(3)=Z10        
!        CALL LTXYZ(B2,ET,EL)     现不这样求	  2018.4.24
      
        T1=SL3/D0                !!交点在NJ1线上的相对位置 0～1 
        T2=-SL1/D0               !!交点在NJ2线上的相对位置 0～1   

!----  !针对曲交线中间的连接点进行坐标移动,理论上可靠，与其他捏点等动作不同，不用担心有问题 2021.9.18
!---  这里有两层问题。1）相邻线段，第一个首先移动了，下一个就不会出现交点紧挨端点；
!-- 2）该迹线（曲交线）在另一个面上，遵循同样的移动规则，能够实现相同移动，且不同面上T1/T2计算结果几乎相等，
!-- 从而保证移动后的坐标几乎相同。通过后续总体编号，同一点只有一个坐标值
            tms=0.025      !10个曲面的计算，修改了该值，0.05改为0.025。0.05时候得到块体没有切开，
                           !采用下面直接goto 9992还有一个为1个面的块体，2021.11.24
              !发现在K(i,j)关系错误，遂改为0.025   2021.11.30
             NLN=NJFRD0+NJFDE1+NL0
            if(ie.le.NLN.AND.NOLOJT(KIJ+i,1).le.NLN.AND.      !位于非曲面，并且形成i线段和j线段的面均为非曲面
     #                       NOLOJT(KIJ+J,1).le.NLN)goto 9992
!------         goto 9992             !不进行坐标移动
            IT=0
            IF(i.GT.1.and.NOLOJT(KIJ+i,1).eq.NOLOJT(KIJ+i-1,1))THEN   !保证调整的是曲交线中间的连接点
             if(t1.lt.tms)then      
             COLOJT(KIJ+I,1:2)=COLOJT(KIJ+I,1:2)+
     #              (COLOJT(KIJ+I,3:4)-COLOJT(KIJ+I,1:2))*(-2*tms)   !当与端点的距离为tms时，坐标移动量是tms的两倍
!-220119!                WRITE(1001,*)'vt101',COTOJT(KIJ+i,1:3)
!-220119!                WRITE(1001,*)'vt101',COTOJT(KIJ+i,4:6)
!-220119!                 WRITE(1001,*)'vt101',COTOJT(KIJ+j,1:3)
!-220119!                WRITE(1001,*)'vt101',COTOJT(KIJ+j,4:6)
             COTOJT(KIJ+I,1:3)=COTOJT(KIJ+I,1:3)+
     #               (COTOJT(KIJ+I,4:6)-COTOJT(KIJ+I,1:3))*(-2*tms)  
!-220119!      WRITE(1001,*)'vt11',COTOJT(KIJ+i,1:3)
 
            COLOJT(KIJ+i-1,3:4)=COLOJT(KIJ+i,1:2)
            COTOJT(KIJ+i-1,4:6)=COTOJT(KIJ+i,1:3)
             IT=1
             ENDIF
            endif  
           IF(i.LT.N1.and.NOLOJT(KIJ+i,1).eq.NOLOJT(KIJ+i+1,1))THEN  
               if(t1.gt.1-tms)then   
             COLOJT(KIJ+I,3:4)=COLOJT(KIJ+I,3:4)+
     #            (COLOJT(KIJ+I,3:4)-COLOJT(KIJ+I,1:2))*(2*tms)
!-220119!       WRITE(1001,*)'vt12',COTOJT(KIJ+i,1:3)
!                WRITE(1001,*)'vt2',COTOJT(KIJ+i,4:6)
             COTOJT(KIJ+I,4:6)=COTOJT(KIJ+I,4:6)+
     #             (COTOJT(KIJ+I,4:6)-COTOJT(KIJ+I,1:3))*(2*tms)   
!               WRITE(1001,*)'vt',COTOJT(KIJ+J,4:6)
            COLOJT(KIJ+i+1,1:2)=COLOJT(KIJ+i,3:4)
            COTOJT(KIJ+i+1,1:3)=COTOJT(KIJ+i,4:6)
!-220119!           WRITE(1001,*)'vt12',COTOJT(KIJ+i+1,1:3)
              IT=2
             endif
             endif  
!-    
          IF(J.GT.1.and.NOLOJT(KIJ+J,1).eq.NOLOJT(KIJ+J-1,1))THEN 
             if(t2.lt.tms)then   
             COLOJT(KIJ+J,1:2)=COLOJT(KIJ+J,1:2)+
     #            (COLOJT(KIJ+J,3:4)-COLOJT(KIJ+J,1:2))*(-2*tms)
!-220119!            WRITE(1001,*)'vt13',COTOJT(KIJ+J,1:3)
!                 WRITE(1001,*)'vt2',COTOJT(KIJ+J,4:6)
             COTOJT(KIJ+J,1:3)=COTOJT(KIJ+J,1:3)+
     #             (COTOJT(KIJ+J,4:6)-COTOJT(KIJ+J,1:3))*(-2*tms)  
!-220119!            WRITE(1001,*)'vt',COTOJT(KIJ+J,1:3)
            COLOJT(KIJ+J-1,3:4)=COLOJT(KIJ+J,1:2)
            COTOJT(KIJ+J-1,4:6)=COTOJT(KIJ+J,1:3)
            IT=3
             ENDIF
          endif  
          
          IF(J.LT.N1.and.NOLOJT(KIJ+J,1).eq.NOLOJT(KIJ+J+1,1))THEN     
            if(t2.gt.1-tms)then   
             COLOJT(KIJ+J,3:4)=COLOJT(KIJ+J,3:4)+
     #             (COLOJT(KIJ+J,3:4)-COLOJT(KIJ+J,1:2))*(2*tms)
             
             COTOJT(KIJ+J,4:6)=COTOJT(KIJ+J,4:6)+
     #             (COTOJT(KIJ+J,4:6)-COTOJT(KIJ+J,1:3))*(2*tms)  
!-220119!             WRITE(1001,*)'vt1',COTOJT(KIJ+J,4:6)
            COLOJT(KIJ+J+1,1:2)=COLOJT(KIJ+J,3:4)
            COTOJT(KIJ+J+1,1:3)=COTOJT(KIJ+J,4:6)
            IT=4
             endif
             endif 
    
          IF(IT.EQ.0) goto 9992
              
          DO 300 IG=1,NJFRD0+NJFDE1+NL0+Ncur    
           IF(IG.EQ.IE) GOTO 300
            DO  J1=1,3
              DO  J2=1,3
             B22(J1,J2)=P1j(iG,(J1-1)*3+J2)
              ENDDO
            ENDDO  
            
             KIJH=0
          IF(IG.NE.1) THEN          
	   DO IJ=1,IG-1 
       	 KIJH=KIJH+NJT(IJ)
          ENDDO  
           ENDIF
          N1H=NJT(IG)
        DO 285 IH=1,N1H      
          IZ0H=NOLOJT(KIJH+IH,1)     
          n02=NOLOJT(KIJH+IH,3)
       IF(it.eq.1.and.NOLOJT(KIJ+i,1).EQ.IG.and.IZ0H.EQ.IE)then 
           n01=NOLOJT(KIJ+i,3)  
          n05= n01-INT(n01/2+0.00001)*2  !奇数或偶数，当为偶数时n05=0
        IF(n05.eq.0.and.n01.eq.n02+1.or.
     #     n05.eq.1.and.n01.eq.n02-1)THEN  !总体交棱存储时，相邻两条分属两个相交的面上,并且总是5和6，9和10这样的相邻两条
                                           !该曲交线在另一个面上也只用考虑为曲交线
            COTOJT(KIJH+IH,1:3)=COTOJT(KIJ+I,1:3)
!-220119!       WRITE(1001,*)'vt22', COTOJT(KIJH+IH,1:3)
            COTOJT(KIJH+IH-1,4:6)=COTOJT(KIJ+i,1:3)  
            
            if(iG.GT.NLN)THEN              !另一个面若为曲面
            COLOJT(KIJH+IH,1:2)=COTOJT(KIJ+I,1:2) 
            COLOJT(KIJH+IH-1,3:4)=COTOJT(KIJ+i,1:2)  
            ELSE                        !另一个面若为平面，局部坐标要通过整体坐标求得 2021.9.24
              ET0(1:3)=COTOJT(KIJ+I,1:3)
              CALL TLXYZ(B22,ET0,EL0)
             COLOJT(KIJH+IH,1:2)=EL0(1:2) 
             COLOJT(KIJH+IH-1,3:4)=EL0(1:2)  
            ENDIF    
!            WRITE(1001,*)'ie',ie,ig,n01,n02,it
!            WRITE(1001,'(a5,2F10.4)')'COT', COLOJT(KIJH+IH,1:2)
!            WRITE(1001,'(a5,3F10.4)')'COT',COTOJT(KIJH+IH,1:3)
!             WRITE(1001,*)'COT',KIJH,IH
            goto 9991     !返回重新计算t1,t2
        endif    
       ENDIF   
       
        IF(it.eq.2.and.NOLOJT(KIJ+i,1).EQ.IG.and.IZ0H.EQ.IE)then 
           n01=NOLOJT(KIJ+i,3)
           n05= n01-INT(n01/2+0.00001)*2 
         IF(n05.eq.0.and.n01.eq.n02+1.or.
     #      n05.eq.1.and.n01.eq.n02-1)THEN 
           
            COTOJT(KIJH+IH,4:6)=COTOJT(KIJ+I,4:6)
            COTOJT(KIJH+IH+1,1:3)=COTOJT(KIJ+i,4:6) 
!-220119!          WRITE(1001,*)'vt23',  COTOJT(KIJH+IH+1,1:3)
            if(iG.GT.NLN)THEN              !另一个面若为曲面
            COLOJT(KIJH+IH,3:4)=COTOJT(KIJ+I,4:5) 
            COLOJT(KIJH+IH+1,1:2)=COTOJT(KIJ+I,4:5)  
            ELSE                        !另一个面若为平面，局部坐标要通过整体坐标求得 2021.9.24
             ET0(1:3)=COTOJT(KIJ+I,4:6)
               CALL TLXYZ(B22,ET0,EL0)
             COLOJT(KIJH+IH,3:4)=EL0(1:2)
             COLOJT(KIJH+IH+1,1:2)=EL0(1:2) 
!             WRITE(1001,*)'ie',ie,ig,n01,n02,it
!            WRITE(1001,'(a5,3F10.4)')'COT',COTOJT(KIJ+I,4:6) 
             endif
            goto 9991
        endif    
        ENDIF   
       
       IF(it.eq.3.and.NOLOJT(KIJ+j,1).EQ.IG.and.IZ0H.EQ.IE)then 
           n01=NOLOJT(KIJ+j,3)
           n05= n01-INT(n01/2+0.00001)*2 
         IF(n05.eq.0.and.n01.eq.n02+1.or.
     #      n05.eq.1.and.n01.eq.n02-1)THEN 
           
            COTOJT(KIJH+IH,1:3)=COTOJT(KIJ+j,1:3)           
            COTOJT(KIJH+IH-1,4:6)=COTOJT(KIJ+j,1:3) 
!-220119!        WRITE(1001,*)'vt43',  COTOJT(KIJH+IH,1:3)
            if(iG.GT.NLN)THEN  
             COLOJT(KIJH+IH,1:2)=COTOJT(KIJ+j,1:2)
             COLOJT(KIJH+IH-1,3:4)=COTOJT(KIJ+j,1:2)  
            ELSE   
               ET0(1:3)=COTOJT(KIJ+j,1:3)
                 CALL TLXYZ(B22,ET0,EL0)
                 COLOJT(KIJH+IH,1:2)=EL0(1:2)
                 COLOJT(KIJH+IH-1,3:4)=EL0(1:2)
             endif
!                WRITE(1001,*)'ie',ie,ig,n01,n02,it
!               WRITE(1001,'(a5,3F10.4)')'COT',COTOJT(KIJ+j,1:3)
            goto 9991
        endif    
       ENDIF   
       
        IF(it.eq.4.and.NOLOJT(KIJ+j,1).EQ.IG.and.IZ0H.EQ.IE)then 
           n01=NOLOJT(KIJ+j,3)
           n05= n01-INT(n01/2+0.00001)*2 
         IF(n05.eq.0.and.n01.eq.n02+1.or.
     #      n05.eq.1.and.n01.eq.n02-1)THEN 
            
            COTOJT(KIJH+IH,4:6)=COTOJT(KIJ+j,4:6)
            COTOJT(KIJH+IH+1,1:3)=COTOJT(KIJ+j,4:6)
!-220119!         WRITE(1001,*)'vt33',COTOJT(KIJH+IH+1,1:3)
            if(iG.GT.NLN)THEN  
             COLOJT(KIJH+IH,3:4)=COTOJT(KIJ+j,4:5)
             COLOJT(KIJH+IH+1,1:2)=COTOJT(KIJ+j,4:5) 
            ELSE
             ET0(1:3)=COTOJT(KIJ+j,4:6)   
               CALL TLXYZ(B22,ET0,EL0)
             COLOJT(KIJH+IH,3:4)=EL0(1:2)
             COLOJT(KIJH+IH+1,1:2)=EL0(1:2)
            endif
!            WRITE(1001,*)'ie',ie,ig,n01,n02,it
!            WRITE(1001,'(a5,3F10.4)')'COT',COTOJT(KIJ+j,4:6)
            goto 9991
        endif    
       ENDIF   
       
285         CONTINUE            
300         CONTINUE              

!----------
           
9992     DO K0=1,2
        EL(K0)=COLOJT(KIJ+J,K0)+(COLOJT(KIJ+J,2+K0)-COLOJT(KIJ+J,K0))*T2
         ENDDO         
        DO K0=1,3
        ET(K0)=COTOJT(KIJ+J,K0)+(COTOJT(KIJ+J,3+K0)-COTOJT(KIJ+J,K0))*T2
        ENDDO    !对于曲面，用平均高程Z10转换为整体坐标系得到的交点位置是不对的，因此这样求 2018.4.24
        
      if(NOLOJT(KIJ+J,2).gt.NJFRD0+NJFDE1+NL0)then !存在曲线-曲线、直线-曲线相交，第二个j总是为曲线.2018.9.17
          KZ0=NOLOJT(KIJ+i,2)      !如果两条线段通过同一段线段，后面判断整体点号时也会出错？这样通过段号是否相同不能完全保证结果正确
        if(KZ0.eq.0)then
          KZ0=int(T1*100+0.5)+100000        ! 直线和曲线相交，直线也得到非零的NOLOJT(KIJ+i,2)。因精度原因，不能*10000  2018.9.29
                                           !造成块体数变少，不能切开的原因，精度问题造成同一点的kz0不同，现在使之模糊，但也可能
                                           !造成挨得近的点认为是相同的点。 <调试入口2> 2018.11.18
        endif   
          MZ0=NOLOJT(KIJ+J,2)  
      endif 
      
         CALL PJUDG(NT,NJ1,NJ2,IE,ET,EL,T1,T2,IZ0,JZ0,KZ0,MZ0)

         
649         IF(IE.EQ.ne2)THEN       !输出某个面的线段交点    <调试入口3>
             write(1001,*)'ii',i,j,NT,NJT(IE)
             WRITE(1001,'(a5,3F10.4)')'TJi', COTOJT(KIJ+I,1:3) 
           WRITE(1001,'(a5,3F10.4)')'TJi',COTOJT(KIJ+I,4:6)
              WRITE(1001,'(a5,3F10.4)')'TJj', COTOJT(KIJ+J,1:3) 
           WRITE(1001,'(a5,3F10.4)')'TJj',COTOJT(KIJ+J,4:6)
           WRITE(1001,'(a5,3F10.4)')'ET0', ET(1:3)
            DO K0=1,3
!        ET(K0)=COTOJT(KIJ+i,K0)+(COTOJT(KIJ+i,3+K0)-COTOJT(KIJ+i,K0))*T1
            ENDDO 
!             WRITE(1001,'(a5,3F10.4)')'ET1', ET(1:3)  !比较交点是否接近,与上面ET0的z坐标有细小差别，因位于不同曲面上
              WRITE(1001,*)T1,t2
           ENDIF
        
650         CONTINUE           
            
         do 750 j1=1,NT        ! 端点移动后，把移动之前直接保存的曲线内交点坐标做相应移动,即后续求交时移动了的端点，
          ij0=NOTJC(j1)        ! 统一将之前已保存的端点坐标再更新一下 2021.11.24
         if(ij0.ne.0)then            
           CNODT(j1,1:3)=COTOJT(ij0,1:3)
           CNODL(j1,1:2)=COLOJT(ij0,1:2)  
         endif    
750      CONTINUE     
           if(ie.eq.ne2)then 
            do j1=1,NT      
          write(1001,'(a3,i6,2x,3f10.4)')'NTt',j1,CNODT(j1,1:3)
           enddo
         endif 
        
         
        ITOTAL(1)=NT
        ITOTAL(2)=N1
        if(ie.eq.ne2)  WRITE(1001,*)'NT', NT
         if(NT.LT.3) GOTO 510 
 
        CALL CPTNET(ie) 
        CALL CDIRET(ie) 
        CALL CBLOCK(ie) 
        
        CALL LOOPLOC(ie) !! 不光对圆形结构面，对开挖面也要分析其回路位置。因为尽管
	                 !! 前面删除了位于面域外的交点，但仍可以形成位于面域外的回路 08.5.4        
 	  CALL LOPADDPOINT

        IO=0
100     IO=IO+1
        K1=NDLP(IO) 
        K2=ABS(K1)
        IF (K1.EQ.0) GOTO 500
C---------------不输出外回路-----------
         IF(K1.LT.0)  goto 480
   
            if (ie.eq.ne2) then               !  <调试入口4> 输出某个面的回路         
             
            write(1001,*)'NDLP(IO)==',ie,NPP0+1,k1                      
            DO  IK=1,K2+1
!          write(1001,*) NDLP(IO+IK),CNODL(NDLP(IO+IK),1:2) 
            write(1001,'(I4,3F10.4)') NDLP(IO+IK),CNODT(NDLP(IO+IK),1:3)
           enddo 
 
              DO  IK=1,K2
              coor(1,1:3)=CNODT(NDLP(IO+IK),1:3)
              coor(2,1:3)=CNODT(NDLP(IO+IK+1),1:3)
               Call DLINE(202,coor,ict,icc) 
              enddo           
            endif
C-----          
        NPP0=NPP0+1
        NDLPZ(NPP0,1)=IE    
        NDLPZ(NPP0,2)=1
        NDLPZ(NPP0,3)=NDLP(IO)
!-------   2015.5.4  ------- !增加结构面类型，以便在渗流分析中识别隙宽，增加NDLPZ(,4)
       if(ie.le.NJFRD0) NDLPZ(NPP0,4)=JTYPE(ie)!随机节理+沉积岩+圆形定位结构面最多为7种,
                                               !其中虚拟的圆形定位结构面为类型10
       if(ie.gt.NJFRD0.and.ie.le.NJFRD0+NJFDE1) NDLPZ(NPP0,4)=8  !多边形定位结构面，第8种
       if(ie.gt.NJFRD0+NJFDE1+NL0) NDLPZ(NPP0,4)=8  !曲面，暂定为第8种，以后再分析 2018.4.20
      if(ie.gt.NJFRD0+NJFDE1.and.ie.le.NJFRD0+NJFDE1+NL0)NDLPZ(NPP0,4)=9
                                                                !临空面，第9种             
C--There are K2+1 nodes in every basic loops,such as n1,n2,n3,n4,n1 
	
	  DO 450 IK=1,K2+1
          NDLPQ(KN+IK)=NDLP(IO+IK)    !!各回路的顶点编号存放在一起
	                                !!存放各节理面上的交点编号
450     CONTINUE
         KN=KN+K2+1 
480      IO=IO+K2+1
        	
        GOTO 100
500     CONTINUE
CC-----------ITOTAL(1) 已在 CPTNET,LOPADDPOINT中发生了变化     
510	  NNOD(IE)=ITOTAL(1)       !各节理面上的交点数     
         IM=NNOD(IE)
	IF (IM.EQ.0) GOTO 700

	 DO 540 I=1,IM
	   
	  DO 520 J=1,3
520      CNODTZ(ncnt+I,J)=CNODT(I,J)   ! 依次存节理面1、节理面2...的交点坐标
                     ! 并用NNOD(IE)记录各节理面上的交点数，起到“指针”的作用

 	  DO 530 J=1,5          
530       NODCZ(ncnt+I,J)=NODC(I,J)  !组成节理面上的交点的3个节理面编号
                             	  !NODCZ(，)数组各维含义与CNODTZ( ,)相同    
540       CONTINUE
          ncnt=ncnt+IM
700       CONTINUE 
       
           NNOD0(:)=NNOD(:)
           NNOD(1)=0
          DO 800 IE=2,NJFRD0+NJFDE1+NL0+Ncur   
           NNOD(ie)=NNOD0(ie-1)+NNOD(ie-1)    !修改NNOD以免后面重复累加2021.10.4 
800       CONTINUE      
         
            call DXFEND(202)
            close(202)
          CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
          write(1001,*) " Intersections Number in CNODTZ=",ncnt,
     #	 '        -*-Time:', t1z

          RETURN
          END

C       ##################################
C       INTERSECTION OF LINES IN 3-D
C       #####################################
        SUBROUTINE PJUDG(NT,NJ1,NJ2,IE,ET,EL,T1,T2,IZ0,JZ0,KZ0,MZ0)
        implicit real*8 (a-h,o-z) 
        DIMENSION IA(52000,8),Q(52000,8),CNODT(52000,3),NODC(52000,5)
        DIMENSION CNODL(52000,2),NODP(52000)
        DIMENSION ET(3),EL(3)
        COMMON /IA/IA/Q/Q/CNODT/CNODT/CNODL/CNODL/NODP/NODP/NODC/NODC

         NT=NT+1
        NL1=NT
 
        DO 120 J1=1,3
120     CNODT(NL1,J1)=ET(J1)
        DO 140  J1=1,2
140     CNODL(NL1,J1)=EL(J1)
 
        NODP(NL1)=IE
         NODC(NL1,1)=IE   !组成交点的3个结构面编号 
	   NODC(NL1,2)=IZ0
	   NODC(NL1,3)=JZ0 
         if(ie.eq.-14)then
          write(1001,'(a2,i6,2x,3f10.4)')'NT',NT,ET(1:3)
         endif    
         if(ie.eq.9.and.IZ0.eq.2.and.JZ0.eq.12)write(1001,*)'d2d'
         if(ie.eq.12.and.IZ0.eq.2.and.JZ0.eq.9)write(1001,*)'d3d'
          NODC(NL1,4)=KZ0
          NODC(NL1,5)=MZ0 
	 
         IF(IE.EQ.7)THEN
!             write(1001,*)'NODC',NL1,NODC(NL1,1:3)
!            write(1001,*)CNODT(NL1,1:3)   
           endif
C       --------------------------
        DO 600 L=1,8
        K1=L
        IF (IA(NL1,L).EQ.NJ1.OR.IA(NL1,L).EQ.0) GOTO 610
600     CONTINUE
610     IA(NL1,K1)=NJ1
        Q(NL1,K1)=T1
C       -------------------
        IF (NJ2.EQ.0) GOTO 750
        DO 620 L=1,8
        K2=L
        IF (IA(NL1,L).EQ.NJ2.OR.IA(NL1,L).EQ.0) GOTO 630
620     CONTINUE
C       ------------------------
630     IA(NL1,K2)=NJ2
        Q(NL1,K2)=T2
750     CONTINUE
        RETURN
        END

C       #######################################
C       COMPUT 1-D NET OF EF--I
C       #####################################
        SUBROUTINE CPTNET(ie)
        implicit real*8 (a-h,o-z) 
        DIMENSION IA(52000,8),Q(52000,8),K(52000,8),CNODL(52000,2)
        DIMENSION IH(800),D(800),NODP(52000),ITOTAL(4)
        DIMENSION CNODT(52000,3),NODC(52000,5)
	  DIMENSION ID0(52000),IA0(52000,8),IHZ(3800,500),CNODT0(52000,3)
	  DIMENSION NODC0(52000,5)
        COMMON/IA/IA/Q/Q/NODK/K/CNODL/CNODL
        COMMON/CNODT/CNODT/NODP/NODP/ITOTAL/ITOTAL/NODC/NODC
        COMMON/COMP0/IHZ,IA0,ID0,CNODT0,NODC0
	 
        K(:,:)=0
  	   
        N1=ITOTAL(2)    !迹线数
        N2=ITOTAL(1)    !交点数

!!      砍树前，对迹线上的交点进行排序，并存储   2007.11.16
 
        DO 40 I=1,N1
        I1=0
        DO 15 J=1,N2
        DO 10 L=1,8
        IF (IA(J,L).NE.I) GOTO 10
        I1=I1+1
        IH(I1)=J
        D(I1)=Q(J,L)
10     CONTINUE
15     CONTINUE
        IF (I1.LT.2) GOTO 40
        IJ1=I1-1
        DO 25 J=1,IJ1
        LJ1=J+1
        DO 20 L=LJ1,I1
        IF (D(J).LE.D(L)) GOTO 20
        NA1=IH(J)
        A2=D(J)
        IH(J)=IH(L)
        D(J)=D(L)
        IH(L)=NA1
        D(L)=A2
20     CONTINUE
25     CONTINUE

       DO 28 J=1,I1
        IHZ(I,J)=IH(J)     !! 记录迹线上交点的排序 2007.11.16
28      CONTINUE
        IHZ(I,500)=I1      !! 记录迹线上交点数  
 
40      CONTINUE

        DO 45 I=1,N2
	   IA0(I,1:8)=IA(I,1:8) 
        CNODT0(I,1:3)=CNODT(I,1:3)
        NODC0(I,1:5)=NODC(I,1:5)
45      CONTINUE

C       ---------------------------
30      N3=0
        DO 230 I=1,N1
        I1=0
        DO 90 J=1,N2
        DO 80 L=1,8
        IF (IA(J,L).NE.I) GOTO 80
        J1=J
        L1=L
        I1=I1+1
80      CONTINUE
90      CONTINUE
        IF (I1.NE.1) GOTO 230
        IA(J1,L1)=0
        N3=1
        I1=0
        DO 200 J=1,8
        IF (IA(J1,J).EQ.0) GOTO 200
        I1=I1+1
        J2=J
200     CONTINUE
        IF (I1.NE.1) GOTO 230
        IA(J1,J2)=0
230     CONTINUE
        IF (N3.EQ.1) GOTO 30
C       -------------

        N4=0
        DO 260 I=1,N2
        I1=0
        DO 240 J=1,8
        IF (IA(I,J).EQ.0) GOTO 240
        I1=I1+1
        J1=J
240     CONTINUE
        IF (I1.NE.1) GOTO 260
        N4=1
        IA(I,J1)=0

260      CONTINUE
        IF (N4.EQ.1) GOTO 30

C       -------------------------------   
        N=0
        DO 400 I=1,N2
        I1=0
        DO 320 J=1,8
        IF (IA(I,J).EQ.0) GOTO 320
        I1=I1+1
320     CONTINUE

        IF (I1.EQ.0) GOTO 400

        N=N+1
	  ID0(N)=I        !砍树后的点号与砍树前的点号对应关系 
        DO 370 J=1,8
        IA(N,J)=IA(I,J)
        Q(N,J)=Q(I,J)
370     CONTINUE
 
       NODP(N)=NODP(I)
       DO 373 J=1,2
373	 CNODL(N,J)=CNODL(I,J)
        DO 374 J=1,3
374     CNODT(N,J)=CNODT(I,J)
        DO 378 J=1,5          !!05.7.20   2018.9.17
378     NODC(N,J)=NODC(I,J)
         
400     CONTINUE
C       ----
	 NB=N+1
	 DO 410 I=NB,N2
	 DO 402 J=1,8
	 IA(I,J)=0
	 Q(I,J)=0.0
402	 CONTINUE
 
	 NODP(I)=0
 
       DO 405 J=1,2
405	 CNODL(I,J)=0
	 DO 406 J=1,3
406	 CNODT(I,J)=0
	 DO 408 J=1,5   
408	 NODC(I,J)=0

410	 CONTINUE

       IF (N.LT.3) GOTO 900
C       ---------------------------
        DO 840 I=1,N1
        I1=0
        DO 550 J=1,N
        DO 540 L=1,8
        IF (IA(J,L).NE.I) GOTO 540
        I1=I1+1
        IH(I1)=J
        D(I1)=Q(J,L)
540     CONTINUE
550     CONTINUE

        IF (I1.EQ.0) GOTO 840
        IJ1=I1-1
        DO 670 J=1,IJ1
        LJ1=J+1
        DO 660   L=LJ1,I1
        IF (D(J).LE.D(L)) GOTO 660
        NA1=IH(J)
        A2=D(J)
        IH(J)=IH(L)
        D(J)=D(L)
        IH(L)=NA1
        D(L)=A2
660     CONTINUE
670     CONTINUE
C       ------
        JI1=I1-1
        DO 830 J=1,JI1
        J1=IH(J)
        J2=IH(J+1)
        DO 750 L=1,8
        L1=L
        IF (K(J1,L).EQ.0.OR.K(J1,L).EQ.J2) GOTO 760
750     CONTINUE
760     K(J1,L1)=J2

        DO 810 L=1,8
        L2=L
        IF (K(J2,L).EQ.0.OR.K(J2,L).EQ.J1) GOTO 820
810     CONTINUE
820     K(J2,L2)=J1
830     CONTINUE
840     CONTINUE

900      CONTINUE
        ITOTAL(1)=N
 
        RETURN
        END

C       ################################
C       COMPUTING DIRECTIONS OF SEGMENTS FROM A NODE (1) IN 2-D
C       ###############################
        SUBROUTINE CDIRET(ie)
        implicit real*8 (a-h,o-z) 
         DIMENSION CNODL(52000,2),IA(52000,8),Q(52000,8)
        DIMENSION K(52000,8),IV(10,2),V(10,2),CNODT(52000,3)
        DIMENSION NODP(52000),ITOTAL(4),VD(2,2)
        COMMON /IA/IA/Q/Q/NODK/K/CNODL/CNODL/NODP/NODP/CNODT/CNODT
        COMMON/ITOTAL/ITOTAL
        N=ITOTAL(1)
	   DD=180.0/3.1415926
       
         Q(:,:)=0
        IA(:,:)=0 
         V(:,:)=0.0
        IV(:,:)=0	   

        DO 510 I=1,N

80      I1=0
        NIE=NODP(I)

        DO 250 J=1,8
        IF (K(I,J).EQ.0) GOTO 250
        J1=K(I,J)
!         if(ie.eq.3.and.i.EQ.32)write(1001,*)'J1=',J1
         VD(1,1:2)=CNODL(I,1:2)
         VD(2,1:2)=CNODL(J1,1:2)         
         
         CALL dipdirection(VD,d1)   !2018.5.4        
     
         I1=I1+1
        IV(I1,1)=J1
        V(I1,1)=D1
250     CONTINUE
C       --------------------
260     JI1=I1-1
        DO 360  J=1,JI1
        LJ=J+1
        DO 360 L=LJ,I1
        IF (V(J,1).LE.V(L,1)) GOTO 360
        JA1=IV(J,1)
        A1=V(J,1)
        IV(J,1)=IV(L,1)
        V(J,1)=V(L,1)
        IV(L,1)=JA1
        V(L,1)=A1
360     CONTINUE
C       ------------------------
        JI1=I1-1
        DO 420 J=1,JI1
        IV(J,2)=IV(J+1,1)
        V(J,2)=V(J+1,1)-V(J,1)
420     CONTINUE
        IV(I1,2)=IV(1,1)
        V(I1,2)=360.0+V(1,1)-V(I1,1)
C       -----------------------
        if (ie.eq.-2) write(1001,'(A9,3f10.4)')'K-CNODL=',CNODT(I,1:3) 
        DO 500 J=1,I1
        IA(I,J)=IV(J,1)
        K(I,J)=IV(J,2)
         if (ie.eq.-2) then
           write(1001,*) i,j, K(I,J)
           write(1001,'(3f10.4)')CNODT(K(I,J),1:3)               
         endif
        Q(I,J)=V(J,2)

500     CONTINUE
510     CONTINUE
        RETURN
        END


C       ############################
C       COMPUTING BLOCKS (1)
C       #############################
        SUBROUTINE CBLOCK(ie)
        implicit real*8 (a-h,o-z) 
        DIMENSION NLP(20,2),CNODL(52000,2),CNODT(52000,3)
        DIMENSION NDLP(188000),IH(8000),IA(52000,8),K(52000,8)
        DIMENSION Q(52000,8),ITOTAL(4)
        COMMON/NLP/NLP/NDLP/NDLP/IA/IA
        COMMON/NODK/K/Q/Q/CNODT/CNODT/CNODL/CNODL 
        COMMON /ITOTAL/ITOTAL
       
        NDLP(:)=0      
	    
        N=ITOTAL(1)

        IO=0
        K0=1
60      DO 130 I=K0,N
        DO 120 J=1,8
        IF (K(I,J).LE.0) GOTO 120
        I1=I
        J1=IA(I,J)
        J0=J
        GOTO 150
120     CONTINUE
130     CONTINUE

        GOTO 530
C       ---------------------
150     K0=I1
        IP=0
        D0=0.0
250     IP=IP+1
        K(I1,J0)=-K(I1,J0)
        IH(IP)=J1
        DO 288 L=1,8
        I0=L
        IF (IA(J1,L).EQ.I1) GOTO 290
288     CONTINUE
290     D0=D0+Q(J1,I0)
        L1=ABS(K(J1,I0))
        I1=J1
        J1=L1
        DO 328 L=1,8
        J0=L
        IF (IA(I1,L).EQ.J1) GOTO 330
328     CONTINUE
330     IF (K(I1,J0).GT.0) GOTO 250
C       ---------------------
350     IO=IO+1

        NDLP(IO)=IP
         if (ie.eq.-8.or.ie.eq.-9) write(1001,*)'IP=', IP,D0
        IF (D0.LT.IP*180.0) GOTO 380
        NDLP(IO)=-IP
380     DO 400 L=1,IP
        NDLP(IO+IP+1-L)=IH(L)   ! 根据角度大小，IH的次序为顺时针，
400     CONTINUE          ! 但存储为回路时，次序相反，因此回路为顺时针
        NDLP(IO+IP+1)=IH(IP)
        IO=IO+IP+1
        NDLP(IO+1)=0

        GOTO 60
530     CONTINUE

C       ---------------------
C         goto 800
        IO=0
100     IO=IO+1
        K1=NDLP(IO)
	  K2=ABS(K1)
 
        IF(K1.EQ.0) GOTO 800
!        IF(K1.LT.0) GOTO 780
         if (ie.eq.-14) then    !这里控制输出回路，以进行纠错检查
           write(1001,*)'NDLP(IO)=', ie,k1
          DO 110 IK=1,K2+1
         write(1001,'(i6,3f12.5)') NDLP(IO+IK),CNODt(NDLP(IO+IK),1:3)
!         write(1001,'(2f12.5)')CNODL(NDLP(IO+IK),1:2)
110     CONTINUE
	  endif
780     IO=IO+K2+1
        GOTO 100
800     CONTINUE 

        RETURN
        END 

C       ############################
C       ANALYSIS LOOP LOCATION 
C       #############################
        SUBROUTINE LOOPLOC(ie) 
        implicit real*8 (a-h,o-z) 
        DIMENSION NDLP(188000),CNODL(52000,2),CNODT(52000,3),VL(500,3)
	  DIMENSION EL(3),ET(3),B2(3,3)
        COMMON/NDLP/NDLP/CNODL/CNODL/B2T/B2,Z10/CNODT/CNODT           

        IO=0
100     IO=IO+1
        K1=NDLP(IO)
	  K2=ABS(K1)
        IF(K1.EQ.0) GOTO 500
        IF(K1.LT.0) GOTO 480
	
      if(ie.le.NJFRD0+NJFDE1+NL0)then
              
	  DO 200 IK=1,K2+1
       VL(IK,1:2)=CNODL(NDLP(IO+IK),1:2) 
200     CONTINUE
       CALL MC19(VL,K1,XL,YL)
        EL(1)=XL
        EL(2)=YL
        EL(3)=Z10       !! 曲面时，是平均高程 2018.4.24 
	          
       CALL LTXYZ(B2,ET,EL) 
       CALL PNT_SUR(ie,ET,M0)          
       
        else               !!针对曲回路不能采用平均高程，要逐点进行分析  2021.12.21
         DO 210 IK=1,K2+1
          ET(1:3)=CNODT(NDLP(IO+IK),1:3)            
           CALL PNT_SUR(ie,ET,M0)    
           IF(M0.EQ.0)GOTO 220               
210     CONTINUE    
        endif
        
220     IF (M0.EQ.0) then
        NDLP(IO)=-NDLP(IO)  !!该回路内的一个点不位于岩体区域内,则该回路
	                      !!不存在.将该回路的点数计为负数(如同外回路).
	  endif	    		  !!后面针对回路点数为负的回路不进行LOPADDPOINT
                            !!分析,也不进行存储  
480     IO=IO+K2+1
        	
        GOTO 100
500     CONTINUE
        RETURN
        END


C       ##################################################################
C       由于回路相邻2个交点间可能有被删除的其它交点，而该交点在其它
C       开挖面（或节理面）上会出现，因此造成不同的面在共棱上的交点不一致而出错.
C       因此需要对这些删除的交点进行补充
C       ###################################################################
        SUBROUTINE LOPADDPOINT
        implicit real*8 (a-h,o-z) 
        DIMENSION NDLP(188000),NDLP0(188000),NODC0(52000,5),IAD(20)
	  DIMENSION ID0(52000),IA0(52000,8),IHZ(3800,500),CNODT0(52000,3)
        DIMENSION CNODT(52000,3),NODC(52000,5),ITOTAL(4),CNODL(52000,2)
        COMMON/NDLP/NDLP/ITOTAL/ITOTAL/CNODT/CNODT/CNODL/CNODL/NODC/NODC
        COMMON/COMP0/IHZ,IA0,ID0,CNODT0,NODC0

          N2=ITOTAL(1)   !交点数

         DO 10 I=1,88000
          NDLP0(I)=NDLP(I)
10        NDLP(I)=0

        JO=0               
        IO=0
100     IO=IO+1

        K1=NDLP0(IO)
        IF (K1.EQ.0) GOTO 500
         K2=ABS(K1)
        IF(K1.LT.0)  goto 480      !不进行外回路分析
         JO=JO+1                     !记录新的NDLP()的编号变化
	   JPO=JO                     ! 记录该回路交点数的存储位置，因为后面
	                              ! 存储编号是逐步增加的而不便查找
	  NDLP(JO)=NDLP0(IO)

	  DO 450 IK=1,K2    
        J1=NDLP0(IO+IK)
        J2=NDLP0(IO+IK+1)
	   JO=JO+1
	   NDLP(JO)=NDLP0(IO+IK)       !! 先存线段的起始点，接着存插入点，
	                               !! 回路最后要存回路的起点
        M1=ID0(J1)                     !! 砍树前的点号 
        M2=ID0(J2)

       DO 250 L1=1,8
	 DO 250 L2=1,8
       IF (IA0(M1,L1).EQ.IA0(M2,L2)) THEN   !!位于同一迹线上
        MJ0=IA0(M1,L1)                      !!迹线号

       GOTO 300
       ENDIF 
250     CONTINUE	  

300     DO 320 L1=1,IHZ(MJ0,500) 
         IF(IHZ(MJ0,L1).EQ.M1) THEN
	     I1=L1 
	   ENDIF
         IF(IHZ(MJ0,L1).EQ.M2) THEN
	     I2=L1 
	   ENDIF
320    CONTINUE  
      
        IF (I1.EQ.I2+1.OR.I1.EQ.I2-1)  GOTO 450 
                               !!相邻的点，中间不存在被删除的点
        IF (I1.GT.I2) THEN
	  DO 330 L1=1,I1-I2-1
!        if(MJ0.eq.247)  write(1001,*)'IAD()=',IHZ(MJ0,I1-L1) 检查插入的点号是否正确
330      IAD(L1)=IHZ(MJ0,I1-L1)    !! 记录中间内插的点号（砍树前的点号）

         I10=I1-I2-1
	  ELSE
	  DO 340 L1=1,I2-I1-1
340      IAD(L1)=IHZ(MJ0,I1+L1)
          I10=I2-I1-1
	  ENDIF

	  DO 400 L1=1,I10
	    N4=0
        DO 360 L2=1,N2
        IF (IAD(L1).EQ.ID0(L2)) THEN   !! 该点已存在
	   N4=1
         L4=L2
	   GOTO 370
	  ENDIF
360     CONTINUE

370     IF(N4.EQ.0) THEN
         N2=N2+1
          ID0(N2)=IAD(L1)       !!插入的新点与砍树前点号对应关系 
	   JO=JO+1
	   NDLP(JO)=N2            !! 插入点为新点   
	    NDLP(JPO)=NDLP(JPO)+1
	    
     	  CNODT(N2,1:3)=CNODT0(IAD(L1),1:3)     !! 插入点的坐标
        NODC(N2,1:5)=NODC0(IAD(L1),1:5)       !! 插入点3个面组成
      
        ELSE
          JO=JO+1
	    NDLP(JO)=L4           !! 插入点为已出现的点
	    NDLP(JPO)=NDLP(JPO)+1
             
	 ENDIF

400    CONTINUE
  
450     CONTINUE
	   JO=JO+1
	   NDLP(JO)=NDLP0(IO+k2+1)     !! 回路最后要存回路的起点

480     IO=IO+K2+1
        	
        GOTO 100
500     CONTINUE

        ITOTAL(1)= N2       !交点数发生变化        

        RETURN
        END

C       ############################################################
C       #######        交点与开挖面的相对位置分析            #######
C       ##  这样，比在形成块体后对块体的位置进行分析更能提高计算速度
C       ############################################################
       
         SUBROUTINE PNT_SUR(ie,ET,M1)      
         implicit real*8 (a-h,o-z) 
         
	   DIMENSION P2E(100,5),P1E(100,9)
	   DIMENSION ET(3),KZ1(100),P(100,4)
 	   DIMENSION M4(100),M5(12,12),M51(12),NC1(99),NL1(99)
         DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100)
         DIMENSION X(5000),XY0(5000,2),NEcur(100)
         COMMON/P2E/P2E/P1E/P1E/KZ1/KZ1/N10/N9,N10/ND0/ND0,NC1,NL1  
         COMMON/ENL/N0,NL0/EML/ML,MLV/ITYPE/ITYPE  
         COMMON/M45/NA9,M4,M5,M51/Control/Ctl1,Ctl2 
         COMMON/NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1/Nfadd/Kfadd,Nfadd
         COMMON/Ncur1/CoCur,XX,NDcur,NEcur
!!       下句为  <易出错的小量控制>  <容差分析过程> 
        Y0=-1.0e-4*Ctl2        

        DO 100 I=1,NL0          ! 简化起见，有效面=定义的面  
	   IF(KZ1(I).EQ.0) THEN
	   DO 10 L1=1,4
10	    P(I,L1)=P2E(I,L1)
         ELSE
	   DO 20 L1=1,4
20	    P(I,L1)=-P2E(I,L1)
      	ENDIF
100    CONTINUE

        M1=1                       ! M1=1时在岩体区域内    
	 DO 200 I=1,N9          ! 首先分析岩体模型范围，洞室时N10 <> 0，边坡时N10=0
        If (M4(i).NE.0) GoTo 200     !在此都是凸体的面，不涉及凹形区 
!        if(i.eq.Nfdel) GOTO 200     !交点位置判断时不用该面作为边界 2021.12.21
                                     !但由于输入的模型边界是通过块体形态分析得到的，面的大小已经限定
          !也就不可能出现在模型边界外的块体了，因此，如果作为补充边界面的曲面部分位于模型范围内
          !必然是模型边界的部分区域会作为实际的边界。因此，该行注销 2021.12.22
        A1=P(I,1)*ET(1)+P(I,2)*ET(2)+P(I,3)*ET(3)+P(I,4)
          If(A1.LT.Y0*100) Then   !Y0绝对值越大，越模糊.针对存在的曲面，因精度而*100.  2021.12.21
	    M1=0  	                !“全部大于才有，一个小于就没有”  
	  RETURN            
          ENDIF
        
200    CONTINUE         
     
        if(Nfadd.eq.0) goto 909                ! 2021.12.21
        if(ie-NJFRD0-NJFDE1-NL0.eq.Nfadd)goto 909    
       
         N1=NDcur(Nfadd)                 
          XY0(1:N1-3,1:2)=CoCur(Nfadd,1:N1-3,1:2)  
          X(1:N1)=XX(Nfadd,1:N1)      
          Xt=ET(1); Yt=ET(2)
        Call CurvPlanePoint(N1,XY0,X,Xt,Yt,Zt)
      
         IF(Kfadd.EQ.0)then
          ZtET=ET(3)-Zt               !保留曲面上方的块体
         else
           ZtET=Zt-ET(3)              !保留曲面下方的块体
         endif
          If(ZtET.LT.Y0*1e3) Then   !考虑到曲面的计算精度，如果曲回路的某个点位于该曲面之上而
                                        !大部分位于之下，需保证能够保留该曲面 2021.12.21          
           M1=0  	                
!           if(ie.eq.8)write(1001,*)'Zt-ET(3)',Zt,ET(1:3)
	    RETURN              
          ENDIF
       
909    if(ITYPE.EQ.2)then
         N7=N9
       do 220 i0=1,ND0               ! 2018.7.6 
          M2=0
       DO 210 I=N7+1,N7+NL1(i0)      ! 针对洞室形态 
        If (M4(i) .NE. 0)  GoTo 210 
        if(ie-NJFRD0-NJFDE1.eq.i)GOTO 910        
	
       A1=P(I,1)*ET(1)+P(I,2)*ET(2)+P(I,3)*ET(3)+P(I,4)
        If(A1.gt.-Y0) GOTO 210  !这里的Y0是否需要*100，还没有研究 2021.12.22
	    M2=1  	           !“全部大于时没有”   2018.1.4          
210    CONTINUE
        if(M2.eq.0) then
          M1=0
          return
        endif
          
         N7=N7+NL1(i0) 
220     CONTINUE
  
         endif
       
      If (NA9.EQ.0) GoTo 910
       
	  DO 500 L = 1,NA9

!!交点与凹区的面无关,即不位于其中的面之上(与块体形态分析不同,没有交点与面有关的情况)
       If (M51(L).NE. 0) Then              !!  2类凹形区
          W1 = 1
		w11 = 1
        DO 300 L1 = 1,M51(L)               !!  对于凹形区的凸面
         i1 = M5(L, L1)
        A1 = p(i1,1)*ET(1)+p(i1,2)*ET(2)+p(i1,3)*ET(3)+p(i1,4)
            If (A1.GT.Y0) GoTo 300     !! 全部满足时,交点存在
               w11 = 0                      !!即交点位于凹形区的全部凸面之前
300      CONTINUE
           If (w11.GT.0.5)   W1 = 0

       DO 310 L2 = M51(L) + 1, M5(L, 12)     !!  ###对于凹形区的凹面
             i2 = M5(L, L2)
         A1 =p(i2, 1)*ET(1)+p(i2,2)*ET(2)+p(i2,3)*ET(3)+p(i2,4)
          If(A1.GT.Y0)   W1 = 0     !!交点位于凹形区凹面之前,则交点存在
310      CONTINUE
!!  '' '########对于凹形区的凹面
        ELse                                !!!  1类凹形区
            W1 = 1
        DO 400 L1 = 1, M5(L, 12)
            i1 = M5(L, L1)
        A1 =p(i1,1)*ET(1)+p(i1, 2)*ET(2)+p(i1,3)*ET(3)+p(i1, 4)
        If (A1.LT.Y0)  GoTo 400    !! 全部满足时,交点不存在,即交点位于
                                         !!  所有的面之后
              W1 = 0       
400      CONTINUE
        End If

       If (W1.LT.0.5)  GoTo 500
       M1 = 0
       GoTo 910
500    CONTINUE
!!--------
910      RETURN
      
         END 
         
C       ##########################
C      因曲交线求解精度问题，2条曲交线在与第3条曲交线相交时，分别在第3条曲交线的相邻线段相交
C     这样就导致后续回路连接不一致的错误。现替换曲回路中的曲交线线段  2021.10.5
C       ##########################
         SUBROUTINE LOOPTOPO
       implicit real*8 (a-h,o-z) 
      
        DIMENSION NDLPZ(1140000,4),NDLPQ(5500000),NDLPZT(1140000,3)
	  DIMENSION NNOD(90000),NODCZ(3400000,5)
	  DIMENSION LPQT(1140000),LPLNT(1140000),coor(2,3)
        DIMENSION NDT(88000),NDLPQT(5500000,2),CNODTZ(3400000,3)
         COMMON /NDLPZ/NDLPZ/NDLPQ/NDLPQ/NDLPZ0/NDLPZ0,NLOC0  
         COMMON/NODCZ/NODCZ/NNOD/NNOD/LPLNT/LPLNT,LPQT 
         COMMON/NPP/NPP0,NPP/NDLPQT/NDLPQT/NDLPZT/NDLPZT 
         COMMON/NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1/ENL/N0,NL0/CNODTZ/CNODTZ
        CHARACTER(10) t1z; CHARACTER(8) d1z	
        
          ne2=-14    !调试的面编号
          
           OPEN(22,FILE='LoopTopo.dxf')
             call DXFHEAD(22) 
           KN=0
           KNT=0
           LPLNT(1)=0
           LPQT(1)=0
           NDLPQT(:,:)=0
           
            do 100 i1=1,NPP0
         
            IE=NDLPZ(I1,1)
            kie=NNOD(ie) 
            J1=NDLPZ(I1,3)
            N1=0
            DO 50 I2=1,J1               
             M0=NDLPQ(KN+I2)
             k01=NODCZ(kie+m0,2) 
             k02=NODCZ(kie+m0,3) 
                      
            if(k01.NE.k02)then    !k01=k02为曲交线的内点（线段连接点）
              N1=N1+1    
              NDLPQT(KNT+N1,1)=M0   !NDLPQT和NDLPQ一样，为了定位并得到CNODTZ和NODCZ数据
               NDT(N1)=i2     
                if(IE==ne2)then
           write(1001,'(a6,2i5,3f10.4)')'nlopt',i1,n1,CNODTZ(kie+m0,1:3)
             endif
            endif           
50          CONTINUE 
            
           IF(NDLPQT(KNT+1,1).NE.NDLPQT(KNT+N1,1))then    !因寻找的初始点原因，可能原回路的首个点为交点，
             N1=N1+1           !则可得到首尾点号相同的拓扑结构回路，否则则得到全部交点，但没有首尾重复
              NDLPQT(KNT+N1,1:2)=NDLPQT(KNT+1,1:2)   !使得拓扑结构回路的首尾相同   
                if(IE==ne2)then
              write(1001,'(a6,2i5,3f10.4)')'nlop0',i1,n1,
     #                                   CNODTZ(kie+NDLPQT(KNT+1,1),1:3)
                endif
            endif 
            
            DO 60 I3=1,N1-1            
              K03=NDT(i3)+1         !对应的原点号的下一点
              M0=NDLPQ(KN+K03)
               k01=NODCZ(kie+m0,2) 
               k02=NODCZ(kie+m0,3)  
             if(k01.eq.k02)then      !对应的原点号的下一点为曲交线的内点
              NDLPQT(KNT+i3,2)=k01  !回路该段是曲交线，形成曲交线的第二个面（第一个面就是IE）       
             endif    
60          CONTINUE     
            
              DO 105 I2=1,n1-1
             if(IE==ne2)then
             icc=3;ict=4
             j3= NDLPQT(KNT+i2,1)               
	       J4=NDLPQT(KNT+i2+1,1)	        
              coor(1,1:3)=CNODTZ(kie+j3,1:3)
              coor(2,1:3)=CNODTZ(kie+J4,1:3)
           Call DLINE(22,coor,ict,icc) 
             endif
105           continue
              
             NDLPZT(i1,3)=N1-1       !此处为N1-1，不是N1
             NDLPZT(i1,1)=IE    
             NDLPZT(i1,2)=1       !没有用到该数据            
             KN=KN+J1+1
             KNT=KNT+N1     !N1已经首尾重复了
             LPLNT(I1+1)=LPLNT(I1)+N1-1
             LPQT(I1+1)=KNT
            
100         CONTINUE
            write(1001,*)" Intersections Number in LPLNT", LPLNT(NPP0+1)
             call DXFEND(22)
              close(22)
         end
         
C       ##########################
C       删除孤立的回路 
C       ##########################
         SUBROUTINE LOOPDELE
       implicit real*8 (a-h,o-z) 
       DIMENSION P2J(90000,5),P1J(90000,9)   
        DIMENSION NDLPZ(1140000,4),NDLPQ(5500000),NDLPZ0(1140000,4)
	  DIMENSION LLOC(5500000),L2NOD(1390000,2),NNOD(90000)
	  DIMENSION NODCZ(3400000,5),NLOC(5500000),NLOC0(5500000)
	  DIMENSION LOP3(5500000),CNODTZ(3400000,3),Cxyz(500,3)   !Cxyz与子程序里面的要长度一致
        dimension Lopt(1140000),PNDLP(1140000,3),KMP0(1140000)
        DIMENSION LPQ(1140000),LPLN(1140000),LPLNT(1140000)
        DIMENSION LPQT(1140000),coor(2,3)
	  DIMENSION NTT(490000,5),CTT(490000,3),LTT(1390000,8)
        DIMENSION NDLPQT(5500000,2),NDLPZT(1140000,3) 
        DIMENSION NCURSEG(4100,420),MCURSEG(11100,5),NLOCT(5500000) 
         COMMON /NDLPZ/NDLPZ/NDLPQ/NDLPQ/NDLPZ0/NDLPZ0,NLOC0  
         COMMON/NODCZ/NODCZ/P1J/P1J/P2J/P2J/LTT/LTT/LPQ/LPQ/LPLN/LPLN 
         COMMON/NLOC/NLOC/NNOD/NNOD/LLOC/LLOC/PNDLP/PNDLP
         COMMON/NPP/NPP0,NPP/CNODTZ/CNODTZ/CTT/CTT/KT/KT
         COMMON/NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1/ENL/N0,NL0
         COMMON/NDLPQT/NDLPQT/NDLPZT/NDLPZT/LPLNT/LPLNT,LPQT/Ncur0/Ncur
        CHARACTER(10) t1z; CHARACTER(8) d1z
	     
          OPEN(22,FILE='LoopExpand.dxf')
            call DXFHEAD(22)  
         write(1001,*) " Loops Number before DELETED:",NPP0 
          
          ne2=-14    !调试的面编号
            
          Lopt(:)=1        
		LTT(:,:)=0  !!总体棱号对应的回路，扩充前有4个.由此在回路删除及相关回路
                  !!分析中只需针对通过相同总体棱号的几个回路进行分析，明显节约时间
	    KT=0
	 
        do 100 i1=1,NPP0
	    IE=NDLPZT(I1,1)
      	kie=NNOD(ie)           
		J1=NDLPZT(I1,3)
        DO 50 I2=1,J1
	  IF(KT.EQ.0) THEN
	     KT=KT+1
	   NLOCT(1)=KT
 
 	    M0=NDLPQT(LPQT(I1)+I2,1)
        NTT(KT,1:5)=NODCZ(kie+m0,1:5)   !NTT(,)：整体编号下各交点的节理面组成
                                        !数组的第1维应明显比CNODTZ( ,)小
        CTT(KT,1:3)=CNODTZ(kie+m0,1:3)   !! 总体交点编号下的坐标 08.5.13     
!          write(1001,*)kt,NODCZ(kie+m0,1:5)
!           write(1001,'(3f10.4)')CNODTZ(kie+m0,1:3)
          GOTO 50
        ENDIF

 	    M0=NDLPQT(LPQT(I1)+I2,1)
	    NZ1=NODCZ(kie+m0,1)  
	    NZ2=NODCZ(kie+m0,2)  
          NZ3=NODCZ(kie+m0,3)
          NZ4=NODCZ(kie+m0,4)
          NZ5=NODCZ(kie+m0,5)    
!        if(NZ1.eq.2.and.NZ2.eq.9.and.NZ3.eq.11.or.
!     #     NZ1.eq.9.and.NZ2.eq.2.and.NZ3.eq.11.or.
!     #      NZ1.eq.11.and.NZ2.eq.2.and.NZ3.eq.9) then
!        write(1001,*)"NODCZ",NODCZ(kie+m0,1:5)         
!        endif
         
	  DO 30 I3=KT,1,-1
 
          NZ6=NTT(I3,1)
          NZ7=NTT(I3,2)
          NZ8=NTT(I3,3)
          NZ9=NTT(I3,4)
          NZ10=NTT(I3,5)   

         if(NZ2.ne.NZ3)then  !NZ2.ne.NZ3包括3种情况:(1)均为平面时,NODC(,4)和NODC(,5)均为0，下述判断等同原判断;(2)两个平面一个曲面，
                             !如曲面与模型边界的交棱相交，此时其实只要用到NODC(,1~3)，由于存储在NODC(,4)的直线段的交点位置相等，
                ! 因此也可以用下述判断；(3)两个曲面和一个平面，此时NODC(,1~3)交错相等而认为是同一点，需要通过NODC(,4~5)是否交错相等
         !进行排除。如两个曲面形成体，又存在平面切一刀的情形，因NODC(,4~5)交错不等而不共点。2018.9.29 
             
	  IF(NZ1.EQ.NZ6.AND.NZ2.EQ.NZ7.AND.NZ3.EQ.NZ8.OR.
     #	 NZ1.EQ.NZ6.AND.NZ2.EQ.NZ8.AND.NZ3.EQ.NZ7.OR.
     #	 NZ1.EQ.NZ7.AND.NZ2.EQ.NZ6.AND.NZ3.EQ.NZ8.OR.
     #	 NZ1.EQ.NZ7.AND.NZ2.EQ.NZ8.AND.NZ3.EQ.NZ6.OR.
     #	 NZ1.EQ.NZ8.AND.NZ2.EQ.NZ6.AND.NZ3.EQ.NZ7.OR.
     #	 NZ1.EQ.NZ8.AND.NZ2.EQ.NZ7.AND.NZ3.EQ.NZ6) THEN  
 
!            if( NZ1.eq.9.and.NZ2.eq.11.and.NZ3.eq.14.or.
!     #           NZ1.eq.11.and.NZ2.eq.9.and.NZ3.eq.14.or.
!     #      NZ1.eq.14.and.NZ2.eq.9.and.NZ3.eq.11) then
!         write(1001,*)"NOc1",i1, NZ1,NZ2,NZ3,NZ4,NZ5
!         write(1001,*)CNODTZ(kie+m0,1:3)
!         write(1001,*)"NOc2", i3, NZ6,NZ7,NZ8,NZ9,NZ10
!         write(1001,*)CTT(I3,1:3)
!            endif 
            
           IF(NZ4.EQ.NZ9.OR.NZ4.EQ.NZ10.OR.              !都为平面时，也能满足，因此判断放在一起
     #	      NZ5.EQ.NZ9.OR.NZ5.EQ.NZ10) THEN    
           if(IE==ne2) then
             write(1001,*)'NODCZ',NODCZ(kie+m0,1:5)
             write(1001,'(3f10.4)')CNODTZ(kie+m0,1:3)
             write(1001,*)'NTT',NTT(I3,1:5)
             write(1001,'(3f10.4)') CTT(I3,1:3)
             endif
!	      NLOCT(LPLNT(I1)+I2)=I3   !  发现回路连接错误，注销。原因是，对于1条曲线穿过一个圆弧的两个交点中，可能第4和5维的数值
!!!            GOTO 50               !  正好有一个相等，则认为是同一点。以前也发现存在这个问题。但当时没有完全研究清楚。2021.11.29
            ENDIF                    !  这样，之前考虑了补充采用距离判断
                        
           t1=(CNODTZ(kie+m0,1)-CTT(I3,1))**2+
     #        (CNODTZ(kie+m0,2)-CTT(I3,2))**2+
     #        (CNODTZ(kie+m0,3)-CTT(I3,3))**2
            t1=sqrt(t1) 
!             write(1001,*)'t1',i1,i3,t1
!             write(1001,*)i1, NZ1,NZ2,NZ3,NZ4,NZ5
!             write(1001,*)i3, NZ6,NZ7,NZ8,NZ9,NZ10   
            if(t1.lt.0.5)then                  !0.3改为0.6，放宽距离的判断。这里要求3个面的构成相同，距离小于给定值
               NLOCT(LPLNT(I1)+I2)=I3          !时认为同一个点  2021.11.22
               GOTO 50
            ENDIF 
          ENDIF    
     
         else         !NZ2=NZ3即分析曲线的前后段，如之前的分析 
             
        IF(NZ1.EQ.NZ6.AND.NZ2.EQ.NZ7.AND.NZ4.EQ.NZ9.OR.     !因曲线段在各个面上保存的排列顺序，只需判断NZ4.EQ.NZ9 
     #	  NZ1.EQ.NZ7.AND.NZ2.EQ.NZ6.AND.NZ4.EQ.NZ9) THEN   
           if(NZ7.EQ.NZ8) THEN  ! NZ2=NZ3时，也要NZ7=NZ8，即比较曲线内的交点  
 
             NLOCT(LPLNT(I1)+I2)=I3
             GOTO 50 
          endif     
         endif  
     
         endif   
         
30       CONTINUE
	   
!!  在平面块体切割中由交点的3个面组成确定是不是同一交点，保证不会出错 
!!  曲面块体切割中，因两条曲线存在多个交点，还需要补充采用距离判断 2021.10.4
         
          KT=KT+1
          NLOCT(LPLNT(I1)+I2)=KT            !!! 总体点号

 	    M0=NDLPQT(LPQT(I1)+I2,1)
          NTT(KT,1:5)=NODCZ(kie+m0,1:5)     
 	    CTT(KT,1:3)=CNODTZ(kie+m0,1:3)    !! 总体交点编号下的坐标
!          write(1001,*)kt,NODCZ(kie+m0,1:5)
!           write(1001,'(3f10.4)')CNODTZ(kie+m0,1:3)
50       CONTINUE
!         KN=KN+J1+1
!         KN1=KN1+J1
100     CONTINUE           
 
        DO i1=1,NPP0
         J1=NDLPZT(I1,3)
        DO I2=1,J1
!         WRITE(1001,*)'NLOCT',I1,I2    
!         WRITE(1001,*)NLOCT(LPLNT(I1)+I2),CTT(NLOCT(LPLNT(I1)+I2),1:3)   
        ENDDO
        ENDDO
        CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
        write(1001,*) " Nodes Quantity in whole
     # Numbering in topologic loops:",KT, '     -*-Time:', t1z
        
C--------- 拓扑结构回路进行扩充 2021.10.7 -----   
      
        
         L0=0                  !所有回路中不重合的曲线段数量
         LPLN(1)=0
       
        DO 200 I1=1,NPP0
            K0=0                !记录扩充后回路i1的顶点数 
          IE=NDLPZT(I1,1)   
         J1=NDLPZT(I1,3) 
      
         DO 180 I2=1,J1           
              K0=K0+1  
         NLOC(LPLN(I1)+K0)= NLOCT(LPLNT(I1)+i2) ! I1回路的i2顶点的整体编号，先保存点i2,紧接着保存可能存在的曲线，形成一个循环  
            if(IE==ne2)   write(1001,*) 'NLOCT',i1,NLOCT(LPLNT(I1)+i2)
        if(NDLPQT(LPQT(I1)+I2,2).eq.0) goto 180  !例如，点1和点2之间有曲线，则NDLPQT(LPQT(I1)+1,2)=1  
        
          N80=0
        IF(L0.NE.0)THEN       !将两点之间的曲线保存在一起，放在每个拓扑结构回路的相应位置，扩充回路
          N1=NLOCT(LPLNT(I1)+i2)     
          IF(I2.LT.J1)THEN
          N2=NLOCT(LPLNT(I1)+i2+1)  
          ELSE
          N2=NLOCT(LPLNT(I1)+1)    
          ENDIF  
          N3=IE                        !面编号
          N4=NDLPQT(LPQT(I1)+I2,2)     !哪个面留下的迹线
!          if(IE==ne2) write(101,*)'n1',N1,N2, N3,N4
              
         DO 115 J11=1,L0               !J11误写为J1
           N5=MCURSEG(J11,4)           
           N6=MCURSEG(J11,5)   
           N7=MCURSEG(J11,2)           !哪个面形成的曲交线留下相同的两个端点
           N8=MCURSEG(J11,3)
!            if(IE==ne2) write(101,*)N5,N6, N7,N8
          IF(N1.EQ.N5.AND.N2.EQ.N6.OR.
     #       N1.EQ.N6.AND.N2.EQ.N5) THEN 
          IF(N3.EQ.N7.AND.N4.EQ.N8.OR.
     #       N3.EQ.N8.AND.N4.EQ.N7) THEN 
           N80=1
           DO J2=1,MCURSEG(J11,1)
              K0=K0+1                       !I1回路的点号+1
              if(N1.EQ.N5)then
             NLOC(LPLN(I1)+K0)=NCURSEG(J11,J2)     
              else
             NLOC(LPLN(I1)+K0)=NCURSEG(J11,MCURSEG(J11,1)+1-J2)    !点号倒过来保存
              endif    
           enddo
           GOTO 180
          ENDIF
          ENDIF
115      continue   
        ENDIF
        
         IF(L0.EQ.0.OR.N80.EQ.0)THEN
        
         IE=NDLPZT(I1,1)
      	kie=NNOD(ie)        
          M0=NDLPQT(LPQT(I1)+I2,1)    !NDLPQT保存各个面上得到的交点,定位CNODTZ和NODCZ
          M1=NDLPQT(LPQT(I1)+I2+1,1)   !因NDLPQT，NDLPQ的每个回路的数据存在首尾重复，就直接+1
          Xn0 = CNODTZ(kie+m0,1)      
          Yn0 = CNODTZ(kie+m0,2)      
          Zn0 = CNODTZ(kie+m0,3)
          XnT = CNODTZ(kie+m1,1)
          YnT = CNODTZ(kie+m1,2)
          ZnT = CNODTZ(kie+m1,3)
          i20=NDLPQT(LPQT(I1)+I2,2)   !留下迹线的i20的面编号更大，但在曲面求交的时候曲面放在前面， 2021.10.11
           if(IE==ne2) then
                write(1001,'(3f10.4)')Xn0,Yn0,Zn0
                write(1001,'(3f10.4)')XnT,YnT,ZnT  
          endif      
          if(ie.le.NJFRD0+NJFDE1+NL0)then
           im0=ie+Ncur          !ie为平面时
          else
           im0=ie-(NJFRD0+NJFDE1+NL0)     !ie为曲面时
          endif
            if(i20.le.NJFRD0+NJFDE1+NL0)then
           im1=i20+Ncur          !ie为平面时
          else
           im1=i20-(NJFRD0+NJFDE1+NL0)     !ie为曲面时
          endif
          
            if(im0.le.im1)then 
              im2=im0;im3=im1 
            else
              im2=im1;im3=im0
            endif    
              
            if(IE==ne2)then
           write(1001,*)'q0',IE,i20,im0,im1,im2,im3 
           endif
          Call CurveLine(Xn0,Yn0,Zn0,XnT,YnT,ZnT,im2,im3,Cxyz,K10)
            if(IE==ne2) then
              do IK=1,K10-1   
               coor(1,1:3)= Cxyz(IK,1:3)  
               coor(2,1:3)= Cxyz(IK+1,1:3)
                 icc=3;ict=3
               write(1001,'(i3,3f10.4)')ik,Cxyz(IK,1:3) 
                Call DLINE(22,coor,ict,icc) 
              enddo
              write(1001,'(i3,3f10.4)')k10,Cxyz(k10,1:3) 
            endif
            
         do IK=1,K10          
            KT=KT+1                        !整体点号接着拓扑结构回路之后，继续增加
!           NTT(KT,1:5)=NODCZ(kie+m0,1:5)   !如同这段线段的首个端点  后面不用NTT了  
              if(IE==ne2)write(1001,*)'KT',KT
           CTT(KT,1:3)= Cxyz(IK,1:3)          
            K0=K0+1                       !I1回路的点号+1
           NLOC(LPLN(I1)+K0)=KT           
         enddo 
         
           L0=L0+1
 
            do IK=1,K10         !L0.EQ.0时，从上面带出来的K10
          NCURSEG(L0,IK)=KT-K10+IK      !存曲线段的总体点号，以便其它回路取用 
           if(IE==ne2)write(1001,*)NCURSEG(L0,IK)
             enddo 
          MCURSEG(L0,1)=K10                     !保存该曲线的顶点数
          MCURSEG(L0,2)=iE
          MCURSEG(L0,3)=i20              
          MCURSEG(L0,4)=NLOCT(LPLNT(I1)+i2)     !保存该曲线两端的整体点号，以便比较点号是否相同
           IF(I2.LT.J1)THEN
          MCURSEG(L0,5)=NLOCT(LPLNT(I1)+i2+1)
           else
          MCURSEG(L0,5)=NLOCT(LPLNT(I1)+1)     
           endif
       
        ENDIF      
         
180      CONTINUE  
         LPLN(i1+1)=LPLN(i1)+k0     
         NDLPZ(I1,3)=k0           !与原回路相比，通过拓扑结构回路扩充形成的新回路的顶点数发生改变
200     CONTINUE     
        
!!!  画出扩充后的回路进行调试分析  
        do 205 i1=1,NPP0  
         J11=NDLPZ(I1,3)
         DO 205 I2=1,J11
             J3=NLOC(LPLN(I1)+I2)   
	   IF(I2.LT.J11) THEN
	       J4=NLOC(LPLN(I1)+I2+1)
	   ELSE
	       J4=NLOC(LPLN(I1)+1)
         ENDIF
        
          if(NDLPZ(I1,1).eq.ne2)then
!            if(i1.eq.3.or.i1.eq.5.or.i1.eq.17)then
           coor(1,1:3)=CTT(j3,1:3)
           coor(2,1:3)=CTT(j4,1:3)
            write(1001,'(a4,3i6,3f10.4)')'ndlp',i1,i2,j3,CTT(j3,1:3)
           icc=4;ict=4
           Call DLINE(22,coor,ict,icc) 
           endif
205      continue         
           call DXFEND(22)
           close(22)         
           
         write(1001,*) " Nodes Quantity in whole Numbering
     # after topologic loop expanded:",KT 
        write(1001,*) " Nodes Quantity in LPLN
     # after topologic loop expanded:",LPLN(NPP0+1) 
        
C---------获得各回路的各棱的总体编号 -----      
         NLT=0
         KN1=0
	  DO 300 i1=1,NPP0
		J10=NDLPZ(I1,3)
        DO 280 I2=1,J10
        IF(NLT.EQ.0) THEN
          NLT=NLT+1                 !! 总体棱号
         LLOC(1)=NLT
	   LTT(NLT,1)=I1
	    L2NOD(NLT,1)=NLOC(1)   !! 该棱2个顶点的总体点号
	    L2NOD(NLT,2)=NLOC(2)   
         GOTO 280
        ENDIF
	  
          L3=NLOC(KN1+I2)
 	  IF(I2.EQ.J10) THEN
	    L4=NLOC(KN1+1)
	    ELSE
          L4=NLOC(KN1+I2+1)
	   ENDIF	

	  DO 250 K1=NLT,1,-1
 
      IF(L2NOD(K1,1).EQ.L3.AND.L2NOD(K1,2).EQ.L4.OR.
     &   L2NOD(K1,1).EQ.L4.AND.L2NOD(K1,2).EQ.L3) THEN  
	    
		  LLOC(KN1+I2)=K1

       DO 260 I4=1,4
	 IF (LTT(K1,I4).EQ.0)THEN
	    LTT(K1,I4)=I1     !存总体棱号对应的回路，即通过该棱的回路数
         GOTO 280           ! 最多有4个. 这里可以存在重复棱. 20080510
       ENDIF
260      CONTINUE

	   GOTO 280
         ENDIF
250      CONTINUE

        NLT=NLT+1              !! 总体棱号
         LLOC(KN1+I2)=NLT
	   LTT(NLT,1)=I1            !!存总体棱号对应的回路 
	    L2NOD(NLT,1)=NLOC(KN1+I2)      !! 该棱2个顶点的总体点号
        IF(I2.EQ.J10) THEN
	    L2NOD(NLT,2)=NLOC(KN1+1)  
	   ELSE
	    L2NOD(NLT,2)=NLOC(KN1+I2+1)  
         ENDIF
280       CONTINUE
          KN1=KN1+J10
300     CONTINUE
       CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*) " Edges Quantity in whole Numbering:",NLT,
     #	 '     -*-Time:', t1z

310 	   M0=0
        DO 500 i1=1,NLT
	   M1=0
        DO 450 I2=1,4
         K1=LTT(i1,i2)
	   IF(LTT(i1,i2).NE.0)THEN
	    IF(Lopt(K1).NE.0)THEN
	     K0=K1
	     M1=M1+1
	    ENDIF
	   ENDIF
450    CONTINUE
        IF(M1.EQ.1) THEN    ! 最后只剩下一个回路通过该棱，删除该回路
	   Lopt(K0)=0         !!!! 标记要删除的回路
         write(1001,*)'The deleted loop is',K0,i1,L2NOD(i1,1:2),
     #     ctt(L2NOD(i1,1),1:3),ctt(L2NOD(i1,2),1:3)
	   M0=1	 
        ENDIF
500     CONTINUE
        IF(M0.EQ.1) GOTO 310

C  ----------重新写数组，将回路编号压缩 ------
        
         LPQ(1)=0     
         LPLN(1)=0
           
        DO 510 I=1,NPP0
	  DO 510 J=1,4
510      NDLPZ0(I,J)=NDLPZ(I,J)

          NLOC0(:)=NLOC(:)
          LOP3(:)=LLOC(:)

	  KM=0
        KN0=0
        KN=0
	  KN11=0
	  KN12=0
       DO 600 I=1,NPP0
	   IF (Lopt(I).Eq.1) THEN
         KM=KM+1
          KMP0(i)=KM             !原回路编号对应新回路编号
        DO 530 J=1,4
            NDLPZ(KM,J)=NDLPZ0(I,J)  
530     CONTINUE
         
        
           IE=NDLPZ(I,1)          
           kie=NNOD(ie)
     
            if (ie.eq.9) then     !--    <调试入口5> 删除回路后，输出某个面的回路
!           write(1001,*)"CN1",i,NDLPZ(I,3)
           DO I2=1,NDLPZ(I,3)
           M0=NDLPQ(KN0+I2)   !不是KN而是KN0，这样才和上面的"CN0"一致 2018.11.15
!           write(1001,*)CNODTZ(kie+m0,1:3)
        enddo
          endif                   !--
          
        PNDLP(KM,1)= P2J(NDLPZ(KM,1),1)     !! 存储法线方向,以便图示中的着色
        PNDLP(KM,2)= P2J(NDLPZ(KM,1),2)
        PNDLP(KM,3)= P2J(NDLPZ(KM,1),3)

	  DO 540 J=1,NDLPZ(I,3)
            NLOC(KN11+J)=NLOC0(KN12+J)
           LLOC(KN11+J)=LOP3(KN12+J)
540      CONTINUE
         KN=KN+NDLPZ0(I,3)+1      !起指针作用
         KN11=KN11+NDLPZ0(I,3) 
         LPQ(KM+1)=KN         !!各回路起始顶点的存储位置,+1后,
	                        ! 第i+1个回路的起点是i结束后位置
         LPLN(KM+1)=KN11       !!各回路起始棱号、点号的存储位置
         ENDIF
          
         
         KN0=KN0+NDLPZ0(I,3)+1     !有些回路被删除，因此KN0>KN
         KN12=KN12+NDLPZ0(I,3)    
600    CONTINUE 

!----写出LTT( ,8)数组，以便在LOOPNEXT中直接对有限的几个有关回路进行分析
!----而不需对全部回路进行循环，明显节约计算时间

       DO 800 I1=1,NLT            !总体棱号在删除回路后不变
        DO 750 I2=1,4
         K1=LTT(i1,i2)
	   IF(K1.NE.0)THEN
	    IF(Lopt(K1).NE.0)THEN
	     LTT(i1,i2)=KMP0(K1)     !原回路编号对应新回路编号
	    ELSE
		 LTT(i1,i2)=0
	    ENDIF
	   ENDIF
750    CONTINUE
800    CONTINUE    
	 CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*) " Loops Number after DELETED (NPP):",KM, 
     #	 '      -*-Time:', t1z

          NPP=KM
C--- 每个回路按2个相反的法线方向扩充成2个编号次序相反的回路 
C--- 输出的NDLPZ(i,j)只到i=NPP，因此后面OTRST1的内容不需更改

       DO 700 I=1,NPP
!            write(1001,*)'npp',i,NDLPZ(I,3)
!            write(1001,*)'NLOC',NLOC(LPLN(I)+1)   !总体点号没有问题 2018.4.24
         NDLPZ(NPP+I,1)=-NDLPZ(I,1)
         NDLPZ(NPP+I,2)=1
         NDLPZ(NPP+I,3)=NDLPZ(I,3)
         NDLPZ(NPP+I,4)=NDLPZ(I,4)             ! 2015.5.4
        PNDLP(NPP+I,1)= -P2J(NDLPZ(i,1),1)
        PNDLP(NPP+I,2)= -P2J(NDLPZ(i,1),2)
        PNDLP(NPP+I,3)= -P2J(NDLPZ(i,1),3)
          
         if (NDLPZ(I,1).eq.2) then     !--
	    DO J=1,NDLPZ(I,3) 
!        write(1001,*)'NLOC',i,j,NLOC(LPLN(I)+j) !,ctt(NLOC(LPLN(I)+j),1:3)
          enddo 
         endif                          !--
CC---------注意NLOC( , )的编号，与NDLPZ( , )不同  
!         NLOC(NPP+I,1)=NLOC(I,1)
         NLOC(KN11+ LPLN(I)+1)=NLOC(LPLN(I)+1) !注意KN11是从上面继承下来的,
	                                         !保证数据向后接续 
	  DO 640 J=2,NDLPZ(I,3)          
         NLOC(KN11+ LPLN(I)+J)=NLOC(LPLN(I)+NDLPZ(I,3)+2-J)
640      CONTINUE
c  --------注意LLOC( , )的编号,回路删除后，总体棱号没有变    
	  DO 650 J=1,NDLPZ(I,3) 
           LLOC(KN11+ LPLN(I)+J)=LLOC(LPLN(I)+NDLPZ(I,3)+1-J)
650      CONTINUE
         LPQ(NPP+I+1)=KN+ LPQ(I)+NDLPZ(I,3)+1
	   LPLN(NPP+I+1)=KN11+LPLN(I)+NDLPZ(I,3)    
700    CONTINUE 
       write(1001,*)" The length of LOOP's vertexes in
     # total:", LPQ(NPP+NPP+1) 

       DO 900 I1=1,NLT
        DO 900 I2=1,4
	  IF(LTT(I1,I2).NE.0)THEN
         LTT(I1,I2+4)=LTT(I1,I2)+NPP !通过总体棱号的最多4个回路扩充成8个
	  ENDIF
900    CONTINUE

        RETURN
        END


C       ########################################
C        ####### 找到回路各边的下一回路    #####
C       ########################################       
           SUBROUTINE LOOPNEXT
        implicit real*8 (a-h,o-z) 
        DIMENSION P2J(90000,5),P1J(90000,9),NLOC(5500000)
        DIMENSION NDLPZ(1140000,4),NLOP(50),LOPNX(5500000)
	  DIMENSION TP(0:50,3),TJ0(3),lopt(1140000),TF(0:50,3)  !,TP0(3),TF0(3)
	  dimension LLOC(5500000),TF1(0:50,3),B2(3,3),ET(3),EL(3),D10(0:50)
	  DIMENSION FP1(50),FP2(50),C2(3,3),LTT(1390000,8),VD(2,2)
        DIMENSION LPQ(1140000),LPLN(1140000),CTT(490000,3)
        DIMENSION CoCur(100,5000,1:3),XX(100,5000),NDcur(100),NEcur(100)
        DIMENSION XY0(5000,3),X(5000),TZ0(3),TZ1(3)
         COMMON /NDLPZ/NDLPZ/NPP/NPP0,NPP/NLOC/NLOC/CTT/CTT 
        COMMON/LOPNX/LOPNX/ENL/N0,NL0/NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1
        COMMON/P1J/P1J/P2J/P2J/LLOC/LLOC/LTT/LTT/LPQ/LPQ/LPLN/LPLN 
        COMMON/Ncur1/CoCur,XX,NDcur,NEcur/DT0/DT0,stp0
	  real*8 T1,T2
        CHARACTER(10) t1z; CHARACTER(8) d1z        
         
   	   lopt(:)=0
        
	   NPP2=NPP*2 

        DO 1000 I1=1,NPP2
	   IE=NDLPZ(I1,1)
      
        DO 990 I2=1,NDLPZ(I1,3)
         N1=NLOC(LPLN(I1)+I2)
         IF(I2.lt.NDLPZ(I1,3)) THEN
	   N2=NLOC(LPLN(I1)+I2+1)
	    ELSE
	   N2=NLOC(LPLN(I1)+1)
	   ENDIF

	  K0=0
	  K2=LLOC(LPLN(I1)+I2)  !! 找到总体棱号，再找到共有8个回路通过该总体棱号
	                  !!总体棱号在删除回路后仍不变
 
	  DO 100 K1=1,8
        J1=LTT(K2,K1)   !!  J1=回路号，最多只需对8个回路进行循环，也包括
	                  !!  该回路本身，以便分析重复棱问题
!      if(i1.eq.3.and.i2.eq.12) write(1001,*)'problems',i1,i2,j1,k1
!        if(i1.eq.126.and.i2.eq.13) write(1001,*)'problems',i1,i2,j1,k1
!        if(i1.eq.55.and.i2.eq.12) write(1001,*)'problems',i1,i2,j1,k1
	  IF(J1.eq.0) goto 100
!      if(i1.eq.138.and.i2.eq.117) write(1001,*)'ps',lopt(j1),NDLPZ(J1,3)
!      if(i1.eq.126.and.i2.eq.13) write(1001,*)'ps',lopt(j1),NDLPZ(J1,3)
!       if(i1.eq.55.and.i2.eq.12)write(1001,*)'ps',lopt(j1),NDLPZ(J1,3)
	   IF(lopt(j1).eq.NDLPZ(J1,3)) goto 100   !如果J1回路作为下一回路已出现了
	                               !与其边数相同的次数，则该回路跳过
	  IF(J1.EQ.(I1+NPP).OR.J1.EQ.(I1-NPP)) GOTO 100  !反回路不作为相关回路

	 DO 80 J2=1,NDLPZ(J1,3)  
 
	   IF(J1.EQ.I1.and.j2.eq.i2) goto 80  !!因为存在多连通回路的重复棱，
	                            !! 重复棱也要找到该棱的下一回路，为该回路本身
         N3=NLOC(LPLN(J1)+J2)
        IF(J2.lt.NDLPZ(J1,3)) THEN
	   N4=NLOC(LPLN(J1)+J2+1)
	    ELSE
	   N4=NLOC(LPLN(J1)+1)
          ENDIF
 
       IF(N1.EQ.N4.AND.N2.EQ.N3) THEN
         K0=K0+1
 	   NLOP(K0)=J1
           if(i1.eq.49.and.i2.eq.101)write(1001,*)'ff',J1
        ENDIF 
80      CONTINUE          
100       CONTINUE         
!         if(i1.eq.226.and.i2.eq.3) write(1001,*)'problems',i1,i2,k0
!         if(i1.eq.230.and.i2.eq.1) write(1001,*)'problems',i1,i2,k0
        if (k0.eq.0.or.k0.gt.3) then
	  write(1001,2100) I1,i2,"'s relative Loops have problems?",k0,ie
        goto 105
	  
        do 101 K1=1,8
	 J1=LTT(K2,K1) 
	  write(1001,*) J1
	IF(J1.eq.0) goto 101
	do J2=1,NDLPZ(J1,3)  

         N3=NLOC(LPLN(J1)+J2)
        IF(J2.lt.NDLPZ(J1,3)) THEN
	   N4=NLOC(LPLN(J1)+J2+1)
	    ELSE
	   N4=NLOC(LPLN(J1)+1)
	   ENDIF
    
       IF(N1.EQ.N4.AND.N2.EQ.N3) THEN
         K0=K0+1
 	   NLOP(K0)=J1
       write(1001,*)'relative loop=', J1 
        ENDIF 
       enddo   
101       CONTINUE         

!!	该情况是由于存在重复棱引起的。通过K00继承上次寻找的结果，也可以得到下一回路，但结果是错的.
!!    在随后块体搜索中一般也不会出现问题.猜测如果以重复棱为起点，得到不在该块体的回路，就
!!		可能会将2个块体同时作为1个块体搜索出来，并进行显示(这样就有问题).
	  endif
105      IF (K0.EQ.1) THEN
	       K00=K0
		 GOTO 900
          ENDIF

 	  TJ0(1)=CTT(N2,1)-CTT(N1,1)  !由回路各顶点对应的总体点号,取得坐标. 080513
	  TJ0(2)=CTT(N2,2)-CTT(N1,2)
	  TJ0(3)=CTT(N2,3)-CTT(N1,3)
        t1=sqrt(TJ0(1)**2+TJ0(2)**2+TJ0(3)**2)
          TJ0(1:3)= TJ0(1:3)/t1
          
       IF(IE.GT.0) THEN       
       DO 110 K1=1,3 
        TP(0,K1)=P2J(IE,K1)             !2018.4.27       
110    CONTINUE
       ELSE
       DO 115 K1=1,3 
        TP(0,K1)=-P2J(-IE,K1)              
115    CONTINUE    
       ENDIF      
        
	  Do 120 K1=1,3
	   C2(1,K1)=TP(0,K1)
120      C2(2,K1)=TJ0(K1) 
        CALL MUTIVT(C2,C1) 
	  DO 130 K1=1,3
        TF(0,K1)=C2(3,K1)
130     CONTINUE     
        
       DO 200 J1=1,K0
	  IL1=NLOP(J1)
        IE1=NDLPZ(IL1,1)
        IF(IE1.GT.0) THEN
       DO 150 K1=1,3 
150	 TP(J1,K1)=P2J(IE1,K1)        
        ELSE
       DO 160 K1=1,3 
160	 TP(J1,K1)=-P2J(-IE1,K1)        
        ENDIF
200    CONTINUE   
       
       do 140 J1=1,K0            !2018.4.27
        Do 121 K1=1,3        
	   C2(1,K1)=TP(j1,K1)
121      C2(2,K1)=TJ0(K1) 
        CALL MUTIVT(C2,C1) 
	  DO 131 K1=1,3
        TF(j1,K1)=-C2(3,K1)     !此处为负
131     CONTINUE
140    CONTINUE          

!曲面时，每段线段旁边的回路法矢不同变化,对tp和tp0进行修正
         do 180 J1=0,K0 
           if (j1.eq.0) then
             IE2=NDLPZ(I1,1)           !此处为IE2而非上面的IE，造成这几天的错误。通过查回路的相关回路，发现有错
           else                        !再查TF, Tp发现矢量有错 2018.4.30
             IE2=NDLPZ(NLOP(J1),1)
           endif
        if (abs(IE2).le.NJFRD0+NJFDE1+NL0)goto 180  
        
         TZ0(1:3)=(CTT(N2,1:3)+CTT(N1,1:3))/2.0   
           
         TZ1(1:3)= TZ0(1:3)+TF(j1,1:3)*Stp0    !沿tf 方向延伸一个步长
!              write(1001,*)'pq',J1,ie2,NLOP(J1)
!        if(i1.eq.55.and.i2.eq.12)write(1001,*)'pf',j1,NLOP(J1),
!     #     TF(j1,1:3), Tp(j1,1:3)
!              write(1001,*)'p0',   TZ0(1:3)
!              write(1001,*)'p1',   TZ1(1:3)
          m1=abs(IE2)-(NJFRD0+NJFDE1+NL0)
         
           N8=NDcur(m1)         !曲面              
           XY0(1:N8-3,1:2)=CoCur(m1,1:N8-3,1:2)  
           X(1:N8)=XX(m1,1:N8)      
            Xt=TZ1(1); Yt=TZ1(2)
        Call CurvPlanePoint(N8,XY0,X,Xt,Yt,Zt)  
             TZ1(3)=Zt 
              
	     C2(1,1:3)=CTT(N2,1:3)-CTT(N1,1:3)
           C2(2,1:3)= TZ1(1:3)-CTT(N1,1:3)
          CALL MUTIVT(C2,C1) 
          IF (j1.eq.0) then
              TP(J1,1:3)=C2(3,1:3)
          ELSE
              TP(J1,1:3)=-C2(3,1:3)
          ENDIF 
         
!            write(1001,*)'pq',TP(J1,1:3)
          if(i1.eq.55.and.i2.eq.12)write(1001,*)'pf',j1,NLOP(J1),
     #     TF(j1,1:3), Tp(j1,1:3)
180      continue       
         
!--   块体没被切开，经查没被切开块体的回路及其相关回路，发现相关回路有错（被错误地用完而
!--   在正确的位置不能用），需要再对TF(0,1:3)进行修正  2018.5.7
         
           C2(1,1:3)=TP(0,1:3)    
           C2(2,1:3)=TJ0(1:3) 
          CALL MUTIVT(C2,C1) 
           TF(0,1:3)=C2(3,1:3)
!!! ---------------
!    经思考，以前这段有些问题，因为为正的不是最大值为夹角最小，为负的也不是最小值
!    为夹角最小。但执行没问题，因事实上最多只有3个相关回路，一个为正，一个接近0，
!    一个为负，下面一段编程也就显得过复杂了
	  M1=0
	  M2=0 
	  DO 300 J1=1,K0
         T1=TF(0,1)*TP(J1,1)+TF(0,2)*TP(J1,2)+TF(0,3)*TP(J1,3)
       IF(T1.GT.1.0e-8) THEN  !改成-0.000001或阈值不同，不影响结果. 正的阈值决定两个基本重合的回路，当为正向夹角时，理应
                   !为相关回路时，却为最后找到。负的阈值时，当为负向夹角时，却作为相关回路。这两种极端情况不会出现 2018.5.7
        M1=M1+1
	  FP1(M1)=J1
	  ELSE
        M2=M2+1
	  FP2(M2)=J1
        ENDIF
300     CONTINUE   

	  IF(M1.EQ.0) GOTO 600
	 if (M1.EQ.1) THEN
	      K00=FP1(M1)
	  GOTO 900
	 ENDIF
         T2=1   
	 DO 420 J1=1,M1
	  MJ1=FP1(J1)
	   T1=TP(MJ1,1)*TP(0,1)+TP(MJ1,2)*TP(0,2)+TP(MJ1,3)*TP(0,3)
	   IF((T1-T2).LT.-0.0000001) THEN
            T2=T1
	      K00=MJ1
		ENDIF
420     CONTINUE   
        GOTO 900

600    if (M2.EQ.1) THEN
	   K00=FP2(M2)
	  GOTO 900
	 ENDIF
         T2=-1  
	 DO 620 J1=1,M2
	   MJ1=FP2(J1)
	   T1=TP(MJ1,1)*TP(0,1)+TP(MJ1,2)*TP(0,2)+TP(MJ1,3)*TP(0,3)
	   IF((T1-T2).GT.0.0000001) THEN
            T2=T1
	      K00=MJ1
		ENDIF
620     CONTINUE   

900      LOPNX(LPLN(I1)+I2)=NLOP(K00)
         lopt(NLOP(K00))=lopt(NLOP(K00))+1        !2007.10.28
!         write(1001,*)'nlop',i1,i2,NLOP(K00)
!     #              ,lopt(NLOP(K00))
         
990     CONTINUE
!          write(*,*)'i1=',i1,NPP2
1000     CONTINUE

          GOTO 2001               !!!!节约时间,不进行检验
 	   do 1200 i=1,2*NPP
         m0=0
        do 1140 i1=1,2*NPP
	   J11=NDLPZ(I1,3)
	  do 1130 j1=1,j11
	   k1=LOPNX(LPLN(I1)+J1)
	  if(k1.eq.i) then    ! 找到与回路i有关系的回路
                            ! 数量为回路i的边数      
	    m0=m0+1 
	if (i.eq.i1) WRITE(1001,*)"Loop i 's next loop is i. ",i
         endif 
1130      continue
1140      continue 

         lopt(i)=m0  
1200     continue
       do 1500 i=1,2*NPP
	if(NDLPZ(i,3).ne.lopt(i))then
       WRITE(1001,2200)'NextLoops =/=edges',
     & 	 i,NDLPZ(i,3),lopt(i),NDLPZ(i,3)-lopt(i)
	endif
1500    continue
        
2001   CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*) " LoopNext Finding finished.          ",
     #	 '               -*-Time:', t1z 

2100    FORMAT(2i7,1X,A32,2x,i2)
2200     FORMAT(a18,I5,1x,i3,1x,i3,1X,i2)   
		RETURN
           END


C       ########################################
C        #######  找到3维块体  ################ 
C       ########################################       
           SUBROUTINE BLOCKCUT3D
        implicit real*8 (a-h,o-z) 
        DIMENSION  LLOC(5500000),NBKC(210000),NDLPZ1(1140000,4)
        DIMENSION NDLPZ(1140000,4),LOPNX(5500000),PNDLP(1140000,3)  
	  DIMENSION TP(20,3),TJ0(3),TP0(3),LPLN(1140000),ET(3),EL(3),B2(3,3)
        DIMENSION LPLN1(1140000),NLOC1(5500000),PNDLP1(1140000,3)
        DIMENSION NDLPK(1140000),LPLP(1140000),NBKC1(210000)
	  integer MLP(4,400000),TMLP(4,800),NB3D(2100000),NB3D1(2100000)
        DIMENSION CTT(490000,3),NLOC(5500000),Knp(1000),NEcur(100)
        DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),P1J(90000,9)
        DIMENSION XY0(5000,3),X(5000),NNEB(3000,3),XYNP(2500,2)
         COMMON /NDLPZ/NDLPZ/NPP/NPP0,NPP/CTT/CTT/LLOC/LLOC
         COMMON /NB3D/NB3D,KB0,NBKC/NLOC/NLOC/PNDLP/PNDLP
         COMMON /LOPNX/LOPNX/KRE/KRE/LPLN/LPLN/KT/KT/P1J/P1J 
         COMMON /NJD0/NJFRD0/NJE0/NJFDE0,NJFDE1/ENL/N0,NL0
         COMMON/Ncur1/CoCur,XX,NDcur,NEcur/IDelu/ANGTOL,IDelu,NMIN
         
        CHARACTER(10) t1z; CHARACTER(8) d1z

	   NPP2=NPP*2
          kd=0
	   KB0=0 
     	   ncnt=0
10        K00=0  

        DO 1000 I1=1,NPP2
         IF(NDLPZ(I1,2).EQ.0) GOTO 1000
        K00=1 
         NDLPZ(I1,2)=0

        M1=NDLPZ(I1,3)  
	  DO 100 J1=1,M1
		MLP(1,J1)=LLOC(LPLN(I1)+J1)
      	MLP(2,J1)=1
       	MLP(3,J1)=I1
	    MLP(4,J1)=J1
100      CONTINUE

105	   K0=0

        DO 200 J2=1,M1
        IF(MLP(2,J2).EQ.1) THEN
	   I20=MLP(3,J2)
         J20=MLP(4,J2)
	  K0=1
       goto 205
	  ENDIF
200      CONTINUE       

        IF(K0.EQ.0) GOTO 700

205      NXL0P=LOPNX(LPLN(i20)+J20)

        IF(NDLPZ(NXL0P,2).EQ.0) then
	     MLP(2,J2)=0       ! J2的值由do 200循环中带出来 
	    goto 105
	  ENDIF

         NDLPZ(NXL0P,2)=0
	   M2=NDLPZ(NXL0P,3)
	  DO 300 J1=1,M2
	   TMLP(1,J1)=LLOC(LPLN(NXL0P)+J1)
         TMLP(2,J1)=1
         TMLP(3,J1)=NXL0P
	   TMLP(4,J1)=J1
300     CONTINUE

      
 	   DO 500 J1=1,M2
		MLP(1,M1+J1)=TMLP(1,J1)
      	MLP(2,M1+J1)=TMLP(2,J1)
       	MLP(3,M1+J1)=TMLP(3,J1)
	    MLP(4,M1+J1)=TMLP(4,J1)
500      CONTINUE
          M1=M1+M2
         GOTO 105
     
700       KB0=KB0+1

           J2=1
         NB3D(J2+ncnt)=MLP(3,1)
          DO 800 J1=2,M1
         IF(NB3D(J2+ncnt).NE.MLP(3,J1)) THEN
		J2=J2+1
		NB3D(J2+ncnt)=MLP(3,J1)     ! 块体1，块体2，..的回路号，
								!   以NBKC(KB0)  为“指针”
         ENDIF
800        CONTINUE
          
           NBKC(KB0)=J2             ! 各块体的回路数         
	       ncnt=ncnt+J2
1000      CONTINUE
        IF(K00.EQ.1) GOTO 10

	 write(1001,*) "Blocks Number (KB0) == :",KB0
        if(ncnt.eq.NPP2) THEN   !所有块体的回路数与总回路数相同
	                          !是保证每个回路只用过1次的基础
           WRITE(1001,2000) "TOTAL Blocks' Face Number:",ncnt," == ",
     # "Loops Number:",NPP2
	  else
           WRITE(1001,2010) "TOTAL Blocks' Face Number:",ncnt," =\= ",
     # "Loops Number:",NPP2,"  ERROR!"	  
        endif  
       
!----------检验所有回路是否用过且只用过1次-----------
        mk=0
	  do 1200 i=1,ncnt
          m0=0
         do 1100 j=1,ncnt
	   if(i.eq.NB3D(j)) then
          m0=m0+1
	   endif
1100    continue
         if (m0.ne.1)  then
	   WRITE(1001,2020)"Loop",i," is NOT used one time! Times are",m0
	     KRE=1
	   mk=1
	endif
1200        continue
       if (mk.eq.0) then
	 WRITE(1001,*)"Block's loops are used only one time, correct!"
       endif   
!     -------检验每个块体是封闭的，即每2个回路间存在公共棱-----
       mk=0
	 DO 1800 I1=1,KB0
	  kib=0
       if (i1.ne.1) then
	  do 1205 j1=1,i1-1
1205        kib=NBKC(j1)+kib
       endif  
    
        IN=NBKC(I1)
       do 1500 K1 = 1,IN  
        LJ1=NB3D(kib+K1)  	!回路号
         J11=NDLPZ(LJ1,3)     !回路的交棱数


	  DO 1420 K2=1,J11
            M1=0
         N1=LLOC(LPLN(LJ1)+K2) 

         DO 1400 L1=1,IN
!!        IF (K1.EQ.L1) GOTO 1400   ! 因为存在多连通回路的重复棱 2007.12.8
	    LJ2=NB3D(kib+L1)	!回路号
		J12=NDLPZ(LJ2,3)

        DO 1380 L2=1,J12	
		N2=LLOC(LPLN(LJ2)+L2)
        if(K1.EQ.L1.and.K2.EQ.L2) goto 1380   !!因为存在多连通回路的重复棱
                                              !!重复棱时，相当于不同回路的不同棱
         IF(N2.EQ.N1) THEN
	    M1=1
	    GOTO 1410
         ENDIF
1380    CONTINUE   
1400     CONTINUE	
1410         IF(M1.EQ.0)  GOTO 1430         	 
1420       CONTINUE
        
1430       IF (M1.EQ.0) THEN
          write(1001,2030)"Block No:",i1," is NOT closed!"
	     mk=1
      	 KRE=1
	  goto 1800
            ENDIF
1500     CONTINUE
1800    CONTINUE        
       CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
        IF(mk.eq.0) write(1001,*)"ALL Blocks are closed!        ",
     #	 '                      -*-Time:', t1z

2000     FORMAT(1X,A27,i7,A4,A13,i7)
2010     FORMAT(1X,A27,i7,A4,A13,i7,A8)
2020     FORMAT(1X,A4,i7,A32,i2)
2030     FORMAT(1X,A10,i7,A15)
!!-------------------------------------------------------         
!!!    针对块体的曲面，进行Delaunary三角划分，然后扩展回路和块体的回路构成 2018.10.7 
!!-----         
          if(IDelu.eq.0) return
          
            N100=NJFRD0+NJFDE1+NL0
           NDLPZ1(:,:)=NDLPZ(:,:)
           LPLN1(:)=LPLN(:)         
           NLOC1(:)=NLOC(:)
           PNDLP1(:,:)=PNDLP(:,:)           
          KM=0    
          LPLN(1)=0          
          LPLP(1)=0
          
       DO 3000 I1=1,NPP  
        IE=NDLPZ1(I1,1)
         
        IF(IE.LE.N100.and.NDLPZ1(I1,3).eq.3)then   !平面上的回路（可能有曲边），当只有3个顶点时，不进行delaunary三角化
        
         KM=KM+1
         NDLPZ(KM,1:4)=NDLPZ1(I1,1:4)     
         PNDLP(KM,1:3)=PNDLP1(I1,1:3)
          DO J1=1,NDLPZ(I1,3)
         NLOC(LPLN(KM)+J1)=NLOC1(LPLN1(I1)+J1)
          enddo 
         LPLN(KM+1)= LPLN(KM)+NDLPZ(I1,3)
         NDLPK(I1)=1              !原回路对应扩展后的回路数
         LPLP(I1+1)=LPLP(I1)+1    !未扩展回路及扩展后的回路依次存放，原回路通过LPLP找到对应的扩展后回路，有1个或多个
              
         ELSE                   !曲面上的回路，以及平面的回路都进行三角化，以采用rhino进行显示   2019.6.1
         open(99,file='Curloop.dat')
          write(99,*)NDLPZ1(I1,3) 
          
        if(IE.LE.N100)then   !!平面上的回路,其顶点需转为局部坐标           
           J02=ABS(NDLPZ1(I1,1))
          DO 1220 K3=1,3
          DO 1220 K4=1,3  
1220      B2(K3,K4)=P1J(J02,(K3-1)*3+K4)  
          
           DO 1280 K3=1,NDLPZ1(I1,3)
            J5=NLOC1(LPLN1(I1)+K3)
             Knp(K3)=J5
          ET(1:3)=CTT(J5,1:3)    	   
          call TLXYZ(B2,ET,EL) 
           Z10=EL(3)
          write(99,'(2f12.6)') EL(1:2)     
1280       continue   
            
        else             
        do 2910 j1=1,NDLPZ1(I1,3)    
         m1=NLOC1(LPLN1(I1)+j1)    
         Knp(j1)=m1              !回路临时编号对应总体点的原编号
        write(99,'(2f12.6)') CTT(m1,1:2)     
2910    continue 
        endif
!         write(99,*) ANGTOL,NMIN
        close(99) 
      
       call Delaunary(ANGTOL,NMIN)    !Delaunary三角化          
        open(98,file='Delau.dat')
        READ(98,*)NP0
        DO K0=1,NP0
         READ(98,*)K1,XYNP(K1,1:2)
        ENDDO
         READ(98,*)NTRI
        DO K0=1,NTRI 
         read(98,*)K1,NNEB(K1,1:3)
        ENDDO
        close(98) 
         write(1001,*)'Delaunary 得到的三角形',i1,ie,NDLPZ1(I1,3),NTRI
         write(*,*)'Delaunary 得到的三角形',i1,ie,NDLPZ1(I1,3),NTRI
      
         if(IE.gt.N100)then  !曲面上的回路
           m1=ABS(IE)-N100        !此处IE>0           
           N8=NDcur(m1)                 
           XY0(1:N8-3,1:2)=CoCur(m1,1:N8-3,1:2)  
           X(1:N8)=XX(m1,1:N8)
         endif   
         
        do 2920 K0=NDLPZ1(I1,3)+1,NP0   !Delaunary剖分时，新增加的点位于NP0的后面 
         KT=KT+1    
         Knp(k0)=KT              !回路内部新增的编号对应总体点的新编号
         
          if(IE.LE.N100)then     !平面上的回路 2019.6.1
            EL(1:2)=XYNP(k0,1:2)    
             EL(3)= Z10
           call LTXYZ(B2,ET,EL) 
            CTT(KT,1:3)=ET(1:3)
          else    
            Xt=XYNP(k0,1); Yt=XYNP(k0,2)
           Call CurvPlanePoint(N8,XY0,X,Xt,Yt,Zt)         
           CTT(KT,1:2)=XYNP(k0,1:2)
           CTT(KT,3)=Zt            !新的总体点号的z坐标
         endif
2920    continue 
        
        DO 2930 K0=1,NTRI       !delaunary后的三角形编入回路（回路只有3个顶点）
          KM=KM+1      
         NDLPZ(KM,1:2)=NDLPZ1(I1,1:2) 
         NDLPZ(KM,3)=3
         NDLPZ(KM,4)=NDLPZ1(I1,4)
         PNDLP(KM,1:3)=PNDLP1(I1,1:3)
         DO J1=1,3              !回路只有3个顶点
         NLOC(LPLN(KM)+J1)=Knp(NNEB(K0,j1))
         enddo
          LPLN(KM+1)= LPLN(KM)+3
2930    continue        
          NDLPK(I1)=NTRI           !原回路对应扩展后的回路数 
          LPLP(I1+1)=LPLP(I1)+NTRI !起到指针作用
         endif 
3000   continue
       
!!!-----每个回路按2个相反的法线方向扩充成2个编号次序相反的回路，与LOOPDELE结尾处相同        
  
         KN11=LPLN(KM+1)     !前面循环带出来的
         KN12=LPLP(NPP+1)     
         do 3050 I1=1,NPP     !对LPLP，NDLPK的扩展与LPLN，NLOC不同，前者针对原回路，后者针对扩展后回路
          NDLPK(NPP+I1)=NDLPK(I1)
          LPLP(NPP+I1+1)=KN12+LPLP(I1)+NDLPK(I1)  
3050     continue   
         
         NPP=KM
      DO 3100 I=1,NPP  
        NDLPZ(NPP+I,1)=-NDLPZ(I,1)
        NDLPZ(NPP+I,2:4)=NDLPZ(I,2:4)
        PNDLP(NPP+I,1:3)= -PNDLP(i,1:3) 
       
        NLOC(KN11+ LPLN(I)+1)=NLOC(LPLN(I)+1)   
        DO J=2,NDLPZ(I,3)  
         NLOC(KN11+ LPLN(I)+J)=NLOC(LPLN(I)+NDLPZ(I,3)+2-J)         
       enddo  
         LPLN(NPP+I+1)=KN11+LPLN(I)+NDLPZ(I,3)  
       
3100  continue         
      
!!------ 
        K10=0
        NBKC1(:)=NBKC(:)  
        NB3D1(:)=NB3D(:)  
        DO 4000 I1=1,KB0 
        J10=NBKC1(I1)      !块体的回路数
  	   kib=0
       if (i1.ne.1) then
	  do ik=1,i1-1
        kib=NBKC1(ik)+kib
        enddo
       endif  
         K20=0              !单个块体的回路数
       DO 3500 K1 = 1,J10 
         LJ=NB3D1(kib+K1)  
        DO ik=1,NDLPK(LJ)     !LJ扩充后的回路数NDLPK(LJ)=1或n
          K10=K10+1
          K20=K20+1
          NB3D(K10)=LPLP(LJ)+ik  !回路LJ前面有LPLP(LJ)个扩展后的回路，然后依次取LJ扩展得到的顺序号即可
        ENDDO  
3500    continue 
        NBKC(I1)=K20         !回路扩充后，块体I1具有K20个回路
!        write(1001,*)'kb0',I1,K20
4000     continue         
 		RETURN
          END

C       ########################################
C        ####### 块体体积计算   ################
C       ########################################       
         SUBROUTINE BLOCKVOL
        implicit real*8 (a-h,o-z) 
        DIMENSION P2J(90000,5),CAZOOM(8,3),P1J(90000,9)
        DIMENSION NDLPZ(1140000,4),NLOC(5500000),CTT(490000,3)
	  DIMENSION NB3D(2100000),NBKC(210000),LPLN(1140000)
        DIMENSION NCONT(210000),NVB(1200),KV5(5)
	  DIMENSION C(3,3),ET(3),EL(3),B2(3,3),NCONT0(1200,2)
	  DIMENSION XC(3),XJ(3),V1(800,2)
        real*8 VOLB(210000),V(4,3)
	  real*8 W7,Z11,Z12,S1,S,vtot,vtot1
         COMMON /NDLPZ/NDLPZ/NB3D/NB3D,KB0,NBKC/VOLB/VOLB
        COMMON/CTT/CTT/P1J/P1J/P2J/P2J/NLOC/NLOC/LPLN/LPLN
        COMMON/CAOM/CAZOOM
        CHARACTER(10) t1z; CHARACTER(8) d1z

	  K5=0 
	 DO 1000 I1=1,KB0
        NCONT(I1)=0
         W7 = 0

	 J10=NBKC(I1)      !块体的回路数
  	   kib=0
       if (i1.ne.1) then
	  do 10 ik=1,i1-1
10       kib=NBKC(ik)+kib
       endif  
!!!------------ 2017.6.2 ---------       
!       Z11 = 999999999.9
!	 DO 100 I2=1,J10
!         J0=NB3D(kib+I2)  !块体的各个回路	  
!        
!         J2=NDLPZ(J0,3)
!	 DO 70 K=1,J2 
!   	  J3=NLOC(LPLN(J0)+K)   !获得总体点号
!        Z12=CTT(J3,3)
!         If (Z12.gt.Z11) GoTo 70
!		Z11=Z12
!70      CONTINUE
!100     CONTINUE
	         
!       do 500 K1 = 1,J10
!         S1 = 0
!         LJ=NB3D(kib+K1)	 
!	   J00=NDLPZ(LJ,1)
!        J0=ABS(J00)
  	   
!        J1 =NLOC(LPLN(LJ)+1)         ! 回路LJ的第1个点
!	  DO 210 L1=1,3
!210     V(1,L1)= CTT(J1,L1)    	   
!       DO 400 K2 = 2,NDLPZ(LJ,3)-1
!	   j2 =NLOC(LPLN(LJ)+K2)
!         J3 =NLOC(LPLN(LJ)+K2+1)
!	  DO 220 L1=1,3
!         V(2,L1)= CTT(J2,L1)
!220      V(3,L1)= CTT(J3,L1)    	   

!        call SUB9400(V,Z11,S)

!        DO 300 k3 = 1,3
!         C(1, k3) = CTT(J2,K3) -CTT(J1,K3)
!         C(2, k3) = CTT(J3,K3) -CTT(J1,K3)
!300     CONTINUE
!          CALL MUTIVT(C,C1)
        
!        A1 = C(3,1)*P2J(J0,1)+C(3,2)*P2J(J0,2)+C(3,3)*P2J(J0,3)
!         if(J00.LT.0) A1=-A1
!        S1 = S1 + SGN(A1) * S
!400     CONTINUE
!!!!        S1 = Abs(S1)
!	  IF(J00.GT.0) THEN
!      	  W7 = W7 - S1
!        ELSE
!	      W7 = W7 + S1
!	  ENDIF
!500    CONTINUE     
!!---------采用第一个回路的第一个点，并在其他回路寻找3点，得到四面体计算体积
!!--------- 计算得到的正负体体积和达到很高的计算精度  2017.6.2 
        
          LJ=NB3D(kib+1)	           !块体的第1个回路	  
!           J00=abs(NDLPZ(LJ,1))      ! 回路所在的面编号
          J4 =NLOC(LPLN(LJ)+1)       !块体的第1个回路的第1个顶点	
          V(4,1:3)= 0  ! CTT(J4,1:3) 第一个面上取一点就忽略了该面形成的体积，当存在曲面时不严格（尤其是存在两个面形成体积时）。
                       ! 因此改为与原点形成四面体进行体积计算   2018.9.12	当曲面有明显起伏时，这样还是不严格，9.13           
           
        DO 700 K1 = 1,J10   ! 2,J10   循环起点改为1，2018.9.12	  
        
         LJ=NB3D(kib+K1)            !块体的各个回路	    
!           J01=abs(NDLPZ(LJ,1))          ! 回路所在的面编号
!         if(J00.eq.J01)goto 700       !同面是跳出，减少工作量.正负体体积和精度减为小数点后6位 
         
         J1 =NLOC(LPLN(LJ)+1)       ! 回路LJ的第1个点
         V(1,1:3)= CTT(J1,1:3)    
         
         do 800 K2 = 2,NDLPZ(LJ,3)-1
           j2 =NLOC(LPLN(LJ)+K2)
           J3 =NLOC(LPLN(LJ)+K2+1)    
           V(2,1:3)= CTT(J2,1:3)
           V(3,1:3)= CTT(J3,1:3)    	    
           call  TET_VOL(V,T)   
           W7 = W7 + T
800      CONTINUE
700     CONTINUE      
        
        VOLB(I1)=W7
       IF( VOLB(I1).LT.0) THEN
        K5=K5+1
       NVB(K5)=I1      !记录负体编号，后面判断只针对负体进行
	 ENDIF
1000    CONTINUE        

!!!---判断腔体，以及未通过复连通回路搜索得到的坑体及环体位于哪个
!!!---正体内部。只需针对负体进行 08.6.6   ------
!        OPEN(50,FILE='BContain1.dat')
           
        X11=(CAZOOM(2,1)-CAZOOM(1,1))*1.2
	 DO 2000 I1=1,K5
       DO 1900 J1=1,KB0
 	 
       IF(VOLB(J1).LT.0) GOTO 1900
       IF(ABS(VOLB(NVB(i1))).GT.VOLB(J1)-0.00001) GOTO 1900

	 J10=NBKC(NVB(i1))       
  	   kib=0
       if (NVB(i1).ne.1) then
	  do ik=1,NVB(i1)-1
         kib=NBKC(ik)+kib
	  enddo
       endif  

       DO 1800 I2=1,J10       !块体NVB(i1)的回路数
!       if(NVB(i1).eq.670) write(50,*) "i2=",i2
        J0=NB3D(kib+I2)   
!	  DO 1700 K21=1,NDLPZ(J0,3)   !需对回路的每个点而不仅是第1个点进行循环
	                        !因为各回路的第1点若恰好与块体J1的回路共面,则不成功 
!还是对第1个点进行循环，并排除可能同面的情况 2009.4.15
   	  J3=NLOC(LPLN(J0)+1)   !获得回路的总体点号
	                     

	 J11=NBKC(J1)          
  	   kib1=0
       if (J1.ne.1) then
	  do ik=1,J1-1
         kib1=NBKC(ik)+kib1
	  enddo
       endif  
        M1=0
       DO 1200 K1=1,J11      !块体J1的回路数
	   LJ=NB3D(kib1+K1)	
	   J01=NDLPZ(LJ,1)
	 DO 1100 K22=1,NDLPZ(LJ,3)
        J4=NLOC(LPLN(LJ)+K22)
        IF(J3.EQ.J4) THEN      !总体点号相同.这样还不能保证完全合理.若点号不同而
	    M1=1                 !同面,则判断可能无交点,因此需要对所有点进行循环而
!	    write(*,*)"the same POINT!"
	    GOTO 1800            !不仅是1个点
        ENDIF

         J02=ABS(J01)
        j00=abs(NDLPZ(J0,1))
!        A1=XC(1)*P2J(J02,1)+XC(2)*P2J(J02,2)+XC(3)*P2J(J02,3)+P2J(J02,4)
!      IF(ABS(A1).LT.0.0001) GOTO 1700 !块体NVB(i1)的点位于块体J1的回路所在的面上
       IF(j00.eq.j02) then  
!	 write(*,*)"the same face!"
	  GOTO 1800          ! 排除同面的情况 2009.4.15
	endif
1100   CONTINUE               
1200   CONTINUE
       
	  DO L1=1,3        
       XC(L1)=CTT(J3,L1)     
        ENDDO          

	 L0=0
       DO 1400 K1=1,J11      !块体J1的回路数
	   LJ=NB3D(kib1+K1)	
	   J01=NDLPZ(LJ,1)
         J02=ABS(J01)
        IF(ABS(P2J(J02,1)).LT.0.00001) GOTO 1400   !避免分母为零，此时无交点
	  A2=XC(1)*P2J(J02,1)+XC(2)*P2J(J02,2)+XC(3)*P2J(J02,3)+P2J(J02,4)
        T1=A2/(P2J(J02,1)*(XC(1)-X11))
         IF(T1.GT.0.99998.OR.T1.LT.0.00002) GOTO 1400 !同面时,通过T1.LT.0.00002
        XJ(1)=XC(1)+(X11-XC(1))*T1                    !避免计入交点数
        XJ(2)=XC(2) 
        XJ(3)=XC(3) 

        DO 1220 K3=1,3
        DO 1220 K4=1,3
1220     B2(K3,K4)=P1J(J02,(K3-1)*3+K4)
        CALL TLXYZ(B2,XJ,EL)
        X0=EL(1)
        Y0=EL(2)
	 DO 1280 K3=1,NDLPZ(LJ,3)
        J5=NLOC(LPLN(LJ)+K3)
	  DO L1=1,3
       ET(L1)=CTT(J5,L1)    	   
        ENDDO
	  CALL TLXYZ(B2,ET,EL)
        DO L1=1,2
       V1(K3,L1)=EL(L1)    	   
        ENDDO
1280     CONTINUE
	  CALL POINTC(X0,Y0,NDLPZ(LJ,3),V1,IW1)
        IF(IW1.EQ.1) THEN
	   GOTO 1400
        ELSE
         L0=L0+1
        ENDIF
1400    CONTINUE
!       if(NVB(i1).eq.670) write(50,*) "i,j=",i1,NVB(i1),j1,L0
        N1=L0-INT(L0/2+0.000001)*2
        N1=1-N1

       IF(N1.EQ.0)THEN        !点位于体内时
        NCONT(NVB(i1))=j1     !负体NVB(i1)位于正体J1内
	  goto 2000
	 ENDIF
!1700   CONTINUE
1800    CONTINUE
1900    CONTINUE 
2000    CONTINUE 
       
	 KCONT=0  
       DO i=1,KB0
	  IF(NCONT(I).NE.0) THEN
	  KCONT=KCONT+1
        NCONT0(KCONT,1)=NCONT(I)     !正体
        NCONT0(KCONT,2)=I            !包含的负体 
        ENDIF
       ENDDO

!!------------------
         vtot=0.0
          vtot1=0.0
        OPEN(30,FILE='B3dcVol.dat')
	   write(30,*) "Block Numer=",KB0
!--  为了在流形法分析中不考虑很小的块体，在此修改了数据格式   2015.4.20 
          write(30,*)"   I=        Vol=       Face Num=     state "
	  DO 2050 i=1,KB0
       if(NBKC(I).gt.3) then
	   IF(VOLB(I).LT.0) THEN
	   write(30,610) I,VOLB(I),NBKC(I) 
     #     ,"Negative Block!"
	   ELSE
           vtot1=vtot1+VOLB(I)
	    write(30,610) I,VOLB(I),NBKC(I),"  true block"
	   ENDIF
	else
	   write(30,610) I,VOLB(I),NBKC(I)
     #     ,"false block!"
       endif
	    vtot=vtot+VOLB(I)
2050    CONTINUE
          write(30,*)" "; write(30,*) "Vol_TOTAL=",vtot
          write(1001,*) "Vol_TOTAL=",vtot
          write(1001,*) "Number of negative blocks contained=",KCONT
          write(30,*) "Positive Block Volume=",vtot1
!!-----------

	 WRITE(30,*)" *--Block_Contained:--"
	 WRITE(30,*)"No.  BlockNumber(+)   Containing(-) "
       DO 2100 i=1,KCONT
	 WRITE(30,670)"I=",i,NCONT0(I,1),NCONT0(I,2)
        vtot1=vtot1+VOLB(NCONT0(I,2))
2100     CONTINUE
      write(30,*) "*--Considering block_containing--*"
       write(30,660) "Positive Block Volume considered
     & block_contained=",vtot1
       write(1001,660) "Positive Block Volume considered
     & block_contained=",vtot1

       DO 2150 I=1,KCONT    !-- 对包含负体的正体体积进行修正
       VOLB(NCONT0(I,1))=VOLB(NCONT0(I,1))+VOLB(NCONT0(I,2))      
2150   CONTINUE

!--------------找出5个最大的正体-----
	 DO I=1,5
	  KV5(i)=0
	 enddo
	  K1=1 
2200      VTT=0
       DO 2300 i=1,KB0
       IF(VOLB(I).LT.0) GOTO 2300
	 IF(VOLB(I).GT.VTT)THEN
	  K3=0
	  DO J=1,K1-1
         IF(I.EQ.KV5(J)) K3=1	
        ENDDO
	 IF(K3.EQ.1) GOTO 2300 
         VTT=VOLB(I)
	  K2=I
	  ENDIF
2300     CONTINUE
	  KV5(K1)=K2
	  K1=K1+1
	 IF(K1.LE.5) GOTO 2200  

       write(30,*)" The biggest 5 block Number,Volume,Percentage(%):"
       write(1001,*)" The biggest 5 block Number,Volume,Percentage(%):"
	 DO I=1,5
        write(30,650) KV5(I),VOLB(KV5(I)),VOLB(KV5(I))/vtot1*100
        write(1001,650) KV5(I),VOLB(KV5(I)),VOLB(KV5(I))/vtot1*100
       ENDDO
       close(30)

       OPEN(40,FILE='BContain.dat')
        write(40,*) KCONT
       DO 2400 i=1,KCONT
	 WRITE(40,*) NCONT0(I,1),NCONT0(I,2)
2400     CONTINUE
        close(40)
       CALL DATE_AND_TIME(TIME = t1z, DATE = d1z)
       write(1001,*)"Blocks volume calculation finished. ",
     #	 '                -*-Time:', t1z

610      FORMAT(I6,2X,F19.6,2X,I6,2X,A16)   !改F16.7为F19.6  2022.7.4
630      FORMAT(1X,A2,I6,2X,A5,F15.5,2X,A9,2X,I6)
640      FORMAT(1X,A23,2X,F12.2)
650      FORMAT(1X,I6,3X,F12.2,3X,F8.4)
660      FORMAT(1X,A50,2X,F12.2)
670      FORMAT(1X,A2,1X,i3,4x,i6,4x,i6)
        Return
        END


C       ###################################
C       OUTPUT CALCULATION RESULTS TO RELEVENT FILES
C       ##################################
        SUBROUTINE OTRST1 
         implicit real*8 (a-h,o-z) 
        CHARACTER*50 FFF1
        DIMENSION NJT(90000),NBKC(210000),PNDLP(1140000,3),ML(1000)
        DIMENSION COTOJT(1300000,6),NDLPZ(1140000,4) 
        DIMENSION CTT(490000,3),NDLPZ0(1140000,4),P2E(100,5)
	  DIMENSION LPLN(1140000),NLOC0(5500000),NLOC(5500000) 
	  DIMENSION NB3D(2100000),V3E(100,50,3),VOLB(210000)
        CHARACTER*60  fl1
        COMMON/NJT/NJT/TOJT/COTOJT/NJE0/NJFDE0,NJFDE1
        COMMON /NJD0/NJFRD0/V3E/V3E/PNDLP/PNDLP/ITYPE/ITYPE
        COMMON/NDLPZ0/NDLPZ0,NLOC0/NLOC/NLOC/CTT/CTT/KT/KT
        COMMON /NDLPZ/NDLPZ/NPP/NPP0,NPP/LPLN/LPLN/Ncur0/Ncur
	  COMMON /NB3D/NB3D,KB0,NBKC/VOLB/VOLB/P2E/P2E
        COMMON/ENL/N0,NL0/EML/ML,MLV/N8/N8/fl1/fl1  
C----------------------
        OPEN(12,FILE='JTB3DC.MAT')
        FFF1='INFORMATION OF JOINT TRACES'
        WRITE(12,100) fl1
        WRITE(12,*) NJFRD0+NJFDE1+NL0+Ncur,N8,NL0-N8      !总面数，临空面数，非临空面数
        WRITE(12,200) (NJT(I),I=1,NJFRD0+NJFDE1+NL0+Ncur)
        NJT0=0
	  do 5 I=1,NJFRD0+NJFDE1+NL0+Ncur
5       NJT0=NJT0+NJT(I)
        
         WRITE(12,*) ((COTOJT(I,K),K=1,6),I=1,NJT0)
	   WRITE(12,300) (((V3E(I,J,K),K=1,3),J=1,ML(i)),I=1,NL0)
        CLOSE (12)
C----------------------
        OPEN(13,FILE='LOOPB3DC.MAT')
        FFF1='INFORMATION OF JOINT LOOPS'
        WRITE(13,100) FFF1
         WRITE(13,600) NPP,KT
         DO 20 I=1,KT                            !前移至此 2015.5.4
          write(13,310) CTT(I,1),CTT(I,2),CTT(I,3)
20       CONTINUE
         
         WRITE(13,200) (LPLN(I),I=1,NPP*2)   

	   DO 10 I=1,NPP
	   WRITE(13,150) I,NDLPZ(I,1),NDLPZ(I,3)      
         WRITE(13,200) (NLOC(LPLN(I)+J),J=1,NDLPZ(I,3))
10        CONTINUE 
        WRITE(13,600) NPP0
          KN=0
	   DO 12 I=1,NPP0
	  WRITE(13,150) I,NDLPZ0(I,1),NDLPZ0(I,3)
        WRITE(13,200) (NLOC0(KN+J),J=1,NDLPZ0(I,3))
         KN=KN+NDLPZ0(I,3)      !不能用LPLN(I)
12       CONTINUE 
          
         DO I=1,NL0
         WRITE(13,610) P2E(i,1:4)   !增加边界面的法向矢量计算结果，以便在渐进失稳分析中计算块体埋深  2017.6.1
         ENDDO 
         CLOSE (13)

C----------------------
        OPEN(20,FILE='BLKCUT3D.MAT')
        FFF1='INFORMATION OF 3D BLOCKS IN SPACE'
        WRITE(20,100) FFF1
        WRITE(20,500) ITYPE            !增加 2018.1.3
        WRITE(20,500) KB0
	  kib=0
	  DO 50 I=1,KB0 
	if (VOLB(I).gt.0) then
	   J10=NBKC(I)
       else
	   j10=-NBKC(I)              !负体时存为负值
	 endif
      	WRITE(20,500) J10
       WRITE(20,200) (NB3D(J),J=kib+1,kib+abs(J10))
	   kib=kib+NBKC(I)
50       continue
       
	  NPP2=NPP*2
	   WRITE(20,500) NPP2
	  DO 60 I=1,NPP2 
	  WRITE(20,155) I,NDLPZ(I,1),NDLPZ(I,3),NDLPZ(I,4)
        WRITE(20,200) (NLOC(LPLN(I)+J),J=1,NDLPZ(I,3))
60        CONTINUE
	  DO 70 I=1,NPP2 
	  WRITE(20,350) PNDLP(i,1),PNDLP(i,2),PNDLP(i,3)  
70        CONTINUE
        CLOSE (20)

100      FORMAT (1X,A50)
150      FORMAT(1X,3(I7,2X))
155      FORMAT(1X,4(I7,2X))
200      FORMAT(1X,10(I8,1X))
300      FORMAT(1X,6(F9.3,1X))
310      FORMAT(1X,3(F16.10,1X))   !增加位数 2016.12.21
350      FORMAT(6(F8.4,1X))
500      FORMAT(1X,I8)
600      FORMAT(1X,4(I6,2X))
610      FORMAT(1X,4(F14.7,1X))
        RETURN
        END

C       ###################################
C       OUTPUT CALCULATION RESULTS TO RELEVENT FILES
C       ##################################
        SUBROUTINE OTRST2
         implicit real*8 (a-h,o-z) 
	   CHARACTER*60 FFF1
        DIMENSION DJFRD(90000,2),CJFRD(90000,3),RJFRD(90000)
        DIMENSION A2(4),E(3),G(3),B2(3,3),NTJT(8),V3D(5000,10,3)
        DIMENSION P2J(90000,5),P1J(90000,9),MJFDE1(5000) 
         COMMON/DJRD/DJFRD/CJRD/CJFRD/RJRD/RJFRD/NJD0/NJFRD0
         COMMON/P1J/P1J/P2J/P2J/NTJT/NTJT/NJE0/NJFDE0,NJFDE1
	   COMMON/MJFDE1/MJFDE1/V3D/V3D
        DD=3.1415926/180

       OPEN(11,FILE='JOINT3D.MAT')
        FFF1='Information of joint circualr disks for 3D display'
        WRITE(11,10) FFF1
        WRITE(11,*) NJFRD0
        WRITE(11,*) NTJT(1)
        WRITE(11,*) (NTJT(I+1),I=1,NTJT(1))
	 
	 DO 100 I1=1,NJFRD0
         DO 22 J1=1,4
22        A2(J1)=P2J(I1,J1)
        DO 24 J1=1,3
        DO 24 J2=1,3
24        B2(J1,J2)=P1J(I1,(J1-1)*3+J2)
        DO 25 J1=1,3
25        E(J1)=CJFRD(I1,J1)
        CALL TLXYZ(B2,E,G)
          X0=G(1)
          Y0=G(2)
         Z10=G(3)
          R=RJFRD(I1)
        DO 60 K=1,36
         X=X0+R*COS(10*K*DD)
	   Y=Y0+R*SIN(10*K*DD)
         G(1)=X
	   G(2)=Y
	   G(3)=Z10
       CALL LTXYZ(B2,E,G)
        WRITE(11,20) (E(J),J=1,3)
60      CONTINUE
100	  CONTINUE
       
	   WRITE(11,*) NJFDE1
        WRITE(11,*) (MJFDE1(I),I=1,NJFDE1)
        DO 200 I=1,NJFDE1
        DO 130 J=1,MJFDE1(I)
         WRITE(11,*) V3D(I,J,1),V3D(I,J,2),V3D(I,J,3) 
130     CONTINUE
200     CONTINUE

        CLOSE (11)
10      FORMAT (1X,A60)
20      FORMAT(1X,3(F9.4,1X))
          RETURN
        END