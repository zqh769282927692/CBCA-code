
          program B3DCdraw
 
cc	   !MS$ATTRIBUTES DLLEXPORT :: B3DCdraw        
c        DIMENSION M0(80,40),M2(40,40),P2(40,5)
c        COMMON /M0EF/M0/M2EF/M2/EQEF/P2/NEF/N5/NODEF/N1
c	  COMMON /KDEF/KDEF
        CALL INRST1
        CALL BLKPLT

        END 


C       ###################################
C       READING CALCULATION RESULTS FROM RELEVENT FILES
C       ##################################
        SUBROUTINE INRST1
        CHARACTER*60 FFF1; CHARACTER*1 y_n

        DIMENSION NJT(90000),COTOJT(1300000,6)
        DIMENSION NDLPZ(1140000,3),NDLPZ0(1140000,3),LPLN(1140000)
	  DIMENSION NDLPZ1(1140000,3),PNDLP(1140000,3),MJFDE1(5000)
	  DIMENSION VV(200,3),ML(100),MLV(100,50),V3D(5000,10,3) 
        DIMENSION NB3D(2100000),V3E(100,50,3),NBKC(210000),NL1(99)
        DIMENSION NLOC0(5500000),NLOC(5500000),NLOC1(5500000)
        DIMENSION vol(210000)
        DIMENSION NTJT(8),V2T(8,3),DXYZ(3240000,3),CTT(490000,3)

        COMMON /NJT/NJT/TOJT/COTOJT/PNDLP/PNDLP
        COMMON /NJD0/NJFRD2,NJFRD0/vol/vol
        COMMON /NTJT/NTJT/COEF/V3E    
        COMMON/V2T/V2T/DXYZ/DXYZ
        COMMON /NDLPZ/NDLPZ,NDLPZ0,NDLPZ1 
	  COMMON/NPP/NPP0,NPP1,NPP2/MJFDE1/MJFDE1
        COMMON /NB3D/NB3D,KB0,NBKC/NJE0/NJFDE1
        COMMON/ENL/N0,NL0/EML/ML,MLV,VV,KZ1   
        COMMON/NJT0/NJT0/LPLN/LPLN/V3D/V3D
        COMMON/NLOC/NLOC,NLOC0,NLOC1/CTT/CTT/KT/KT 
        COMMON/y_n/y_n/Nblock/Nblock,iFaceShow               

C--  输出块体切割后,所有顶点的坐标;每个块体的回路数,回路的顶点数
C--  及点号构成. 见文件BLKCUT3D_reslut.DAT.        -2011.4.24- 
       OPEN(15,FILE='Inpt_B3DC_Draw.DAT')   !改为文件输入 2021.8.4
!       write(*,*)'是否需要形成AutoCad图? Y or N'    !  -2015.8.27- 
       read(15,*) y_n                       !当块体很多时，不形成AutoCad图，
                                           !只用OpengL看，明显节约时间
!       write(*,*)'CAD中画全部块体（选0）或某个块体（选块体编号）? '    !  -2016.11.22- 
       read(15,*) Nblock            !在块体四面体剖分时，调试需要，并发现切割得到的有些块体出错     
       read(15,*) iFaceShow
        close(15)
       OPEN(99,FILE='BLKCUT3D_reslut.DAT')   
        
        OPEN(12,FILE='JTB3DC.MAT')
!        FFF1='INFORMATION OF JOINT TRACES'
        READ(12,110) FFF1
C--------
        OPEN(11,FILE=FFF1)
	   read(11,*) k13
         READ(11,*) Orient
         READ(11,*) ((V2T(I,J),J=1,3),I=1,8)
        CLOSE (11)
cc         添加岩体边界点的坐标旋转
           D1D= 0.0174533
           BS = Orient * D1D
	     DO 500  i=1,8
            XX1 = V2T(i, 1) 
            YY1 = V2T(i, 2)
 	      V2T(i, 1) = XX1 * Sin(BS) - YY1 * Cos(BS)
            V2T(i, 2) = XX1 * Cos(BS) + YY1 * Sin(BS)
500        CONTINUE
c---
       Open (32,FILE='~Bkspe3DP1.dat',STATUS='OLD')    ! 2018.8.31
        read(32,*) NL0   !有效面数
       do 105 i=1, NL0
           read(32,*) ML(i)
!       do 102 J=1,ML(i)
           read(32,*) MLV(i,1:ML(i))  !读入但没有用到
!102      continue
105      continue
          CLOSE (32)     
C---------------------- 
        if(ITYPE.eq.2)then
       Open (32,FILE='~Bkspe3DP2.dat',STATUS='OLD')   
         read(32,*) ND0        !洞室和排水孔的总数 2018.7.6
          do i0=1,ND0
        read(32,*) NL1(i0)   !有效面数
       do 26 i=1, NL1(i0)
           read(32,*) ML(i+NL0)
!       do 21 J=1,ML(i+NL0)
           read(32,*) MLV(i+NL0, 1:ML(i+NL0))
           !     MLV(i+NL0, J)=MLV(i+NL0, J)+N0 与B3dc_sep一致，应该有这句，
           !     但MLV(*)数组没有用到，也没有读入N0，注销这句
!21      continue
26     continue
         NL0=NL0+NL1(i0)
       enddo   
         CLOSE (32)     
        endif
C----------------------
        READ(12,*) NJFRD2,N9,N10     
        READ(12,120) (NJT(I),I=1,NJFRD2)
ccc    (COTOJT(I,J,K) 数组在存储时会发生错误，出现NaN，此问题解决,
ccc  将半径很小的节理在JFZOOM中重生成，并在JTSEF中排除  05.11.13
         NJT0=0
	  do 10 I=1,NJFRD2      
10        NJT0=NJT0+NJT(I)
  
	  READ(12,*) ((COTOJT(I,K),K=1,6),I=1,NJT0)
 	  READ(12,130) (((V3E(I,J,K),K=1,3),J=1,ML(i)),I=1,NL0)  
        CLOSE (12)

C----------------------
        OPEN(13,FILE='LOOPB3DC.MAT')
        FFF1='INFORMATION OF JOINT NETS'
        READ(13,110) FFF1
        READ(13,*)  NPP1,KT
        write(99,*) KT    !! 顶点总数,为回路删除前的顶点数,因此部分可能没有用到
                           !! 但通过CTT获得各点的坐标是正确的
        DO 40 I=1,KT       !! 各点的坐标
	   READ(13,310) CTT(I,1),CTT(I,2),CTT(I,3)
	  write(99,310) CTT(I,1),CTT(I,2),CTT(I,3)
40      CONTINUE
        
        READ(13,120) (LPLN(I),I=1,NPP1*2)   
	   DO 20 I=1,NPP1
	  READ(13,115) I0,NDLPZ1(I,1),NDLPZ1(I,3)
        READ(13,120) (NLOC1(LPLN(I)+J),J=1,NDLPZ1(I,3))
20        CONTINUE 
        READ(13,*)  NPP0                           !! 20071102
	  KN=0
	   DO 30 I=1,NPP0                              
	  READ(13,115) I0,NDLPZ0(I,1),NDLPZ0(I,3)
        READ(13,120) (NLOC0(KN+J),J=1,NDLPZ0(I,3)) 
	   KN=KN+NDLPZ0(I,3)      !不能用LPQ(I)
30        CONTINUE  
         CLOSE (13)

C----------------------
        OPEN(100,FILE='BLKCUT3D.MAT')
        FFF1='INFORMATION OF 3D BLOCKS IN SPACE'
        READ(100,110) FFF1
        READ(100,150) ITYPE    !增加 2018.1.3   
        READ(100,150) KB0
	   write(99,*) KB0         !! 块体总数 
      	kib=0
	  DO 50 I=1,KB0 
	   READ(100,150) J10
	    write(99,150) J10
	    NBKC(I)=J10         !!每个块体的回路数.负体时为负数
	             !!起到“指针”作用,这样NB3D可以连续存储各个块体回路构成编号   
         READ(100,120)(NB3D(J),J=kib+1,kib+abs(J10))
         write(99,120)(NB3D(J),J=kib+1,kib+abs(J10))
            kib=kib+abs(NBKC(I))    !!定位每个块体回路编号在NB3D的起始存储位置
50       continue

	   READ(100,150) NPP2       !!NPP2=NPP1*2
	   write(99,150) NPP2       !!所有回路数
	   write(99,120) (LPLN(I),I=1,NPP1*2)   !!每一个回路第一个顶点的存储位置
                                  !!也如同指针,由NDLPZ(I,3)累加而成. LPLN(1)=0
	  DO 60 I=1,NPP2 
	  READ(100,115) I0,NDLPZ(I,1),NDLPZ(I,3),i1   !!增加i1的读入, 2015.5.5    
        write(99,115) I0,NDLPZ(I,1),NDLPZ(I,3)      !! I0不起作用,NDLPZ(I,1)为
	                        !!回路I所属的结构面编号,NDLPZ(I,3)为回路I顶点数
        READ(100,120) (NLOC(LPLN(I)+J),J=1,NDLPZ(I,3))  
	  write(99,120) (NLOC(LPLN(I)+J),J=1,NDLPZ(I,3))    !! 回路I各顶点编号,
	                                  !! 所有回路编号都存储在一维数组NLOC中
60        CONTINUE 

	  DO 70 I=1,NPP2 
	  READ(100,*) PNDLP(i,1),PNDLP(i,2),PNDLP(i,3) !在此没有用到，也可以不读.
                                                     !现在需要用到2019.6.1
70      CONTINUE

        CLOSE (100)

C----------------------
	  OPEN(20,FILE='JOINT3D.MAT')
        READ(20,110) FFF1
        READ(20,*) NJFRD0
         READ(20,*) NTJT(1)
       READ(20,*) (NTJT(I+1),I=1,NTJT(1))
        NK=NJFRD0*36
         DO 80 I=1,NK
        READ(20,200) (DXYZ(I,J),J=1,3)
80        CONTINUE
	   READ(20,*) NJFDE1
        READ(20,*) (MJFDE1(I),I=1,NJFDE1)
 
        DO 400 I=1,NJFDE1
        DO 350 J=1,MJFDE1(I)
         READ(20,*) V3D(I,J,1),V3D(I,J,2),V3D(I,J,3) 
350     CONTINUE
400     CONTINUE
        CLOSE (20)         
C----------------------
       OPEN(20,FILE='B3dcVol.dat')
         READ(20,110) FFF1
        READ(20,110) FFF1
        do i=1,kb0
        READ(20,*)i0,vol(i0),i1,FFF1
        enddo
           CLOSE (20)

110      FORMAT (1X,A60)
115      FORMAT(1X,4(I7,2X))
120      FORMAT(1X,10(I8,1X))
130      FORMAT(1X,6(F9.3,1X))
150      FORMAT(1X,I8)
300      FORMAT(1X,3(F9.4,1X))
200      FORMAT(1X,3(F9.4,1X))
310      FORMAT(1X,3(F16.10,1X))  !增加位数 2016.12.21

        RETURN
        END

C       ########################
C       TO DRAW EXCAVATION FACES
C       ########################
        SUBROUTINE EFPLT(LV,FFF1,LV1,FOGL,M0,M1)
       CHARACTER*25 FFF1
	 CHARACTER*25 FOGL
        CHARACTER*1 y_n
         DIMENSION VV(200,3),ML(100),MLV(100,50) 
         DIMENSION V2T(8,3),V3E(100,50,3) 
	  COMMON/V2T/V2T/COEF/V3E
	  COMMON/ENL/N0,NL0/EML/ML,MLV,VV,KZ1   
         COMMON/y_n/y_n      
         
        OPEN(LV,FILE=FFF1)
	IF (M0.EQ.1) THEN
        OPEN(LV1,FILE=FOGL)
       ENDIF
  
        WRITE(LV,300)
300     FORMAT(1X,' 0')
        WRITE(LV,304)
304     FORMAT('SECTION')
        WRITE(LV,310)
310     FORMAT(1X,' 2')
        WRITE(LV,320)
320     FORMAT('ENTITIES')

        IF(M1.EQ.0) GOTO 1000
	  KM1=1
           
       DO 400 J=1,4
        X1=V2T(J,1)
        Y1=V2T(J,2)
	  Z1=V2T(J,3)
	  IF (J.LT.4) THEN
        X2=V2T(J+1,1)
        Y2=V2T(J+1,2)
        Z2=V2T(J+1,3)
	 ELSE
        X2=V2T(1,1)
        Y2=V2T(1,2)
        Z2=V2T(1,3)
       ENDIF
       if(y_n.eq.'y'.or.y_n.eq.'Y') CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
        CALL OGLINE(LV1,KM1,X1,Y1,Z1,X2,Y2,Z2)
        
        X2=V2T(J+4,1)
        Y2=V2T(J+4,2)
        Z2=V2T(J+4,3)
       if(y_n.eq.'y'.or.y_n.eq.'Y') CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
        CALL OGLINE(LV1,KM1,X1,Y1,Z1,X2,Y2,Z2)
      
	  X1=V2T(J+4,1)
        Y1=V2T(J+4,2)
	  Z1=V2T(J+4,3)
	  IF (J.LT.4) THEN
        X2=V2T(J+4+1,1)
        Y2=V2T(J+4+1,2)
        Z2=V2T(J+4+1,3)
	 ELSE
        X2=V2T(5,1)
        Y2=V2T(5,2)
        Z2=V2T(5,3)
       ENDIF
       if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
        CALL OGLINE(LV1,KM1,X1,Y1,Z1,X2,Y2,Z2)
         
400     CONTINUE

       KM1=KM1+1

	 DO 600 I=1,NL0
	  DO 500 J=1,ML(i)
        X1=V3E(I,J,1)
        Y1=V3E(I,J,2)
	  Z1=V3E(I,J,3)
	  IF (J.LT.ML(i)) THEN
        X2=V3E(I,J+1,1)
        Y2=V3E(I,J+1,2)
        Z2=V3E(I,J+1,3)
	 ELSE
        X2=V3E(I,1,1)
        Y2=V3E(I,1,2)
        Z2=V3E(I,1,3)
       ENDIF
       if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
        CALL OGLINE(LV1,KM1,X1,Y1,Z1,X2,Y2,Z2)
500     CONTINUE
600     CONTINUE
1000    RETURN
        END

	

C       ###############################
C       TO DRAW A LINE (X1,Y1) TO (X2,Y2)
C       ##############################
        SUBROUTINE DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
        WRITE(LV,130)
130     FORMAT(2X,'0')
        WRITE(LV,140)
140     FORMAT('LINE')
        WRITE(LV,150)
150     FORMAT(2X,'8')
        WRITE(LV,160)
160     FORMAT('0')
        WRITE(LV,170)
170     FORMAT(' 10')
        WRITE(LV,200) X1
        WRITE(LV,180)
180     FORMAT(' 20')
        WRITE(LV,200) Y1
        WRITE(LV,185)
185     FORMAT(' 30')
        WRITE(LV,200) Z1
        WRITE(LV,190)
190     FORMAT(' 11')
        WRITE(LV,200) X2
        WRITE(LV,195)
195     FORMAT(' 21')
        WRITE(LV,200)  Y2
        WRITE(LV,197)
197     FORMAT(' 31')
        WRITE(LV,200)  Z2
200     FORMAT(F15.6)
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   SUBROUTINE OGLINE(LV1,KM,X1,Y1,Z1,X2,Y2,Z2)

        WRITE(LV1,200) KM,X1,Y1,Z1
        WRITE(LV1,201) X2,Y2,Z2

200     FORMAT(I4,3F15.6)
201     FORMAT(3F15.6)
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   SUBROUTINE OGLINE1(LV1,KM2,X1,Y1,Z1)

        WRITE(LV1,200) KM2,X1,Y1,Z1

200     FORMAT(I4,3F15.6)
        RETURN
        END


C       ############################
C       JOINT ROCK SIMULATION PLOT
C       #############################

        SUBROUTINE BLKPLT 
        CHARACTER*25 FFF1,FOGL
         CHARACTER*1 y_n
        DIMENSION COTOJT(1300000,6),NJT(90000) 
        DIMENSION NTJT(8),LJ(6),V3D(5000,10,3) 
	  DIMENSION DXYZ(3240000,3),NBKC(210000),MJFDE1(100)
	  DIMENSION NDLPZ(1140000,3),NDLPZ0(1140000,3),NDLPZ1(1140000,3)
        DIMENSION NLOC0(5500000),NLOC(5500000),NLOC1(5500000)
	  DIMENSION NB3D(2100000),LPLN(1140000),CTT(490000,3)
         DIMENSION B2(3,3),vol(210000)   
        COMMON/NJT/NJT/TOJT/COTOJT 
        COMMON/NJD0/NJFRD2,NJFRD0
        COMMON /DXYZ/DXYZ/NTJT/NTJT
        COMMON /NDLPZ/NDLPZ,NDLPZ0,NDLPZ1
	  COMMON/NPP/NPP0,NPP1,NPP2
        COMMON/V3D/V3D/NJE0/NJFDE1/MJFDE1/MJFDE1
        COMMON /NB3D/NB3D,KB0,NBKC
	  COMMON/NJT0/NJT0/LPLN/LPLN/vol/vol
        COMMON/NLOC/NLOC,NLOC0,NLOC1/CTT/CTT/KT/KT 
        COMMON/y_n/y_n/Nblock/Nblock,iFaceShow  
        
         KM2=3
100     LV=10
        FFF1='JTb3dc.DXF'
	     LV1=11
	  FOGL= '~JTb3dc.DAT'
        CALL EFPLT(LV,FFF1,LV1,FOGL,1,0)
        
	  DO 200 I=1,NJT0
!         if(i.gt.4542.and.i.le.5132)THEN   ! 指定特定的面上显示迹线 2015.8.27 <调试入口1>
!         if(i.gt.2962.and.i.le.3665)THEN   ! 面11
          if(i.gt.8705.and.i.le.9426)THEN   ! 面9
        X1=COTOJT(I,1)         ! 如显示面10，则在B3DC_sep_Cur.for的JTCOTOLO中输入9，10后在B3dcDATA.dat看迹线编号，
                               ! 再将其放在此 2018.1.16 
        Y1=COTOJT(I,2)
        Z1=COTOJT(I,3)
        X2=COTOJT(I,4)
        Y2=COTOJT(I,5)
	  Z2=COTOJT(I,6)
       
	 if(y_n.eq.'y'.or.y_n.eq.'Y') CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
	   CALL OGLINE(LV1,KM2,X1,Y1,Z1,X2,Y2,Z2)
          endif
200     CONTINUE
       
        WRITE(LV,225)
225     FORMAT(1X,' 0')
        WRITE(LV,230)
230     FORMAT('ENDSEC')
        WRITE(LV,225)
        WRITE(LV,240)
240     FORMAT('EOF')
        CLOSE(LV)
		CLOSE(LV1)

CCCCCCC----------------------  !! 回路删除前  20071102

         LV=12
        FFF1='LOOPb3dcB.DXF'  
	 	LV1=13
 	   FOGL= '~LOOPb3dcB.DAT'
        CALL EFPLT(LV,FFF1,LV1,FOGL,1,0)
         
	   KN=0
        DO 450 I=1,NPP0

        J2=NDLPZ0(I,3)
	 DO 410 K=1,J2 
	 J3=NLOC0(KN+K)   
	  IF(K.LT.J2) THEN
	 J4=NLOC0(KN+K+1)
	  ELSE
	 J4=NLOC0(KN+1)
        ENDIF
	  X1=CTT(j3,1)
        Y1=CTT(j3,2)
        Z1=CTT(j3,3)
	  X2=CTT(j4,1)
        Y2=CTT(j4,2)
        Z2=CTT(j4,3)
        if(NDLPZ0(I,1).eq.iFaceShow.or.NDLPZ0(I,1).eq.8)then        !某个面的回路图   2018.4.20 <调试入口2>
                                                             !改为文件输入 2021.8.4
!      if(NDLPZ0(I,1).eq.8.or.NDLPZ0(I,1).eq.12)then  !.or.NDLPZ0(I,1).eq.14
!         if(I==49.or.I==383.or.I==776)then   !.or.I==388-220.or.I==34
            if(I==49)then
                write(1001,*)'NDLPZ0(I,1)',i,NDLPZ0(I,1)  
            endif     
       if(y_n.eq.'y'.or.y_n.eq.'Y') CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
     	 CALL OGLINE(LV1,KM2,X1,Y1,Z1,X2,Y2,Z2)
        endif            
            
!        if(abs(X1-8).lt.0.001.and.abs(y1+14).lt.0.001.
!     #      and.abs(z1+10.6126).lt.0.001)then 
!        write(1001,*)'NDLPZ0(I,1)',NDLPZ0(I,1) 找到通过该点的回路
!         endif
410       CONTINUE
          KN=KN+NDLPZ0(I,3)      !不能用LPLN(I)
450      CONTINUE

        WRITE(LV,225)
        WRITE(LV,230)
        WRITE(LV,225)
        WRITE(LV,240)
        CLOSE (LV)
	  CLOSE(LV1)
CCCCCCC-------------------- !! 回路删除后  20071102

         LV=12
        FFF1='LOOPb3dcA.DXF'   
	 	LV1=13
 	   FOGL= '~LOOPb3dcA.DAT'
        CALL EFPLT(LV,FFF1,LV1,FOGL,1,0)

        DO 451 I=1,NPP1

        J2=NDLPZ1(I,3)
	 DO 411 K=1,J2 

	 J3=NLOC1(LPLN(I)+K)   
	  IF(K.LT.J2) THEN
	 J4=NLOC1(LPLN(I)+K+1)
	  ELSE
	 J4=NLOC1(LPLN(I)+1)
        ENDIF

	  X1=CTT(j3,1)
        Y1=CTT(j3,2)
        Z1=CTT(j3,3)
	  X2=CTT(j4,1)
        Y2=CTT(j4,2)
        Z2=CTT(j4,3)
        if(NDLPZ1(I,1).eq.iFaceShow.or.NDLPZ1(I,1).eq.8)then          ! <调试入口3>
      if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
     	 CALL OGLINE(LV1,KM2,X1,Y1,Z1,X2,Y2,Z2)
        endif
411       CONTINUE

451      CONTINUE

        WRITE(LV,225)
        WRITE(LV,230)
        WRITE(LV,225)
        WRITE(LV,240)
        CLOSE (LV)
	  CLOSE(LV1)
     

CCCCCCCC----------------------
	     LJ(1)=NTJT(2)   
          DO 540 K1=2,NTJT(1)
540       LJ(K1)=LJ(K1-1)+NTJT(K1+1)
       
        KM3=3
           LV=20
        FFF1='JOINT3d.DXF'            
	     LV1=21
	      FOGL= '~JOINT3d.DAT'
         CALL EFPLT(LV,FFF1,LV1,FOGL,1,1)
          K0=1
        DO 600 I=1,NJFRD0
	   IF (I.GT.LJ(K0)) THEN  !不同组,采用不同颜色
          KM3=KM3+1
	     K0=K0+1
	   ENDIF

        DO 580 J1=1,36
	    J2=(I-1)*36+J1
        X1=DXYZ(J2,1)
        Y1=DXYZ(J2,2)
        Z1=DXYZ(J2,3)
       IF (J2.LT.I*36) THEN
	   X2=DXYZ(J2+1,1)
         Y2=DXYZ(J2+1,2)
         Z2=DXYZ(J2+1,3)
       ELSE
	   X2=DXYZ((I-1)*36+1,1)
         Y2=DXYZ((I-1)*36+1,2)
         Z2=DXYZ((I-1)*36+1,3)
	  ENDIF
       if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)   ! 2015.8.22
       CALL OGLINE1(LV1,KM3,X1,Y1,Z1)
580        CONTINUE
600     CONTINUE

       write(LV1,*),"100 -108.0  110.0  -108.0"
	   WRITE(LV1,*) NJFDE1
	    DO I=1, NJFDE1
         WRITE(LV1,*) MJFDE1(I) 
          ENDDO

         KM3=KM3+1
       DO 700 I=1,NJFDE1
	 DO 700 J=1,MJFDE1(I)
        X1=V3D(I,J,1)
        Y1=V3D(I,J,2)
        Z1=V3D(I,J,3)
       IF (J.LT.MJFDE1(I)) THEN
        X2=V3D(I,J+1,1)
        Y2=V3D(I,J+1,2)
        Z2=V3D(I,J+1,3)
       ELSE
        X2=V3D(I,1,1)
        Y2=V3D(I,1,2)
        Z2=V3D(I,1,3)
	  ENDIF
      
        if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
       CALL OGLINE1(LV1,KM3,X1,Y1,Z1)
700     CONTINUE

        WRITE(LV,225)
         WRITE(LV,230)
         WRITE(LV,225)
         WRITE(LV,240)
        CLOSE(LV)
		CLOSE(LV1)
CCCCCCCC----------------------

800     LV=100
        FFF1='BLKCUT3D.DXF'
	 	LV1=101
 	   FOGL= '~BLKCUT3D.DAT'

	  OPEN(LV1,FILE=FOGL)
         write(LV1,*) KB0

	   kib=0
	 DO 850 I=1,KB0 
        WRITE(LV1,*) NBKC(I)          ! 各块体的回路数,以NBKC() 为“指针”  
   	 DO 840 J=1,abs(NBKC(I))
        J0=NB3D(kib+J)                ! 块体1，块体2，..的回路号
        J2=NDLPZ(J0,3)                ! 回路的顶点数
        WRITE(LV1,*) J2
840     continue  
	  kib=kib+abs(NBKC(I))
850     CONTINUE    
    
	  CALL EFPLT(LV,FFF1,LV1,FOGL,0,0)
c---------便于OpenGL图形输出 
       WRITE(LV1,*) KM2,KM2,KM2,KM2    
c---------- 
         OPEN(19,FILE='RhinoBLK.stl')   !-- 2019.6.1
 	   kib=0
         
	  DO 900 I=1,KB0
     
       if(Nblock.eq.0.and.vol(I).gt.0.001.or.I.eq.Nblock)then   !------- 2021.8.10
         write(19,'(A19)')'solid OuterBoundary'  
            endif                                                !--------
         K10=abs(NBKC(I))
!         if(I.eq.5)write(99,*)'K10=',K10
          
        DO 880 J=1,K10 
	  J0=NB3D(kib+J)
         J3=NLOC(LPLN(J0)+1)   
         J4=NLOC(LPLN(J0)+2) 
         J5=NLOC(LPLN(J0)+3)
         B2(1,1:3)=CTT(J5,1:3)- CTT(J3,1:3)
         B2(2,1:3)=CTT(J4,1:3)- CTT(J3,1:3)
           CALL MUTIVT(B2,C1)
             if(Nblock.eq.0.and.vol(I).gt.0.001.or.I.eq.Nblock)then 
         if(vol(I).gt.0.00001)THEN 
           write(19,'(A14,3f10.4)')'  facet normal',B2(3,1:3)  
           write(19,*)'    outer loop'              
           write(19,*)'      vertex',CTT(J3,1:3)  !尽管写出来，但当块体的某个回路顶点不是3个，就不能用rhino显示
           write(19,*)'      vertex',CTT(J5,1:3)
           write(19,*)'      vertex',CTT(J4,1:3)      
         else
           write(19,'(A14,3f10.4)')'  facet normal',-B2(3,1:3)  
           write(19,*)'    outer loop'              
           write(19,*)'      vertex',CTT(J3,1:3)  !尽管写出来，但当块体的某个回路顶点不是3个，就不能用rhino显示
           write(19,*)'      vertex',CTT(J4,1:3)
           write(19,*)'      vertex',CTT(J5,1:3)                
          ENDIF
           endif
      
        J2=NDLPZ(J0,3)
	 DO 870 K=1,J2 
         
	 J3=NLOC(LPLN(J0)+K)   
	  IF(K.LT.J2) THEN
	 J4=NLOC(LPLN(J0)+K+1)
	  ELSE
	 J4=NLOC(LPLN(J0)+1)
        ENDIF
	  X1=CTT(j3,1)
        Y1=CTT(j3,2)
        Z1=CTT(j3,3)
	  X2=CTT(j4,1)
        Y2=CTT(j4,2)
        Z2=CTT(j4,3)
          
      if(Nblock.eq.0)then   !画全部块体 2016.11.22 此时即使前面不形成AutoCad图，在此都输出块体的cad图 2018.7.15
       CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2) !! if(y_n.eq.'y'.or.y_n.eq.'Y')CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)
      ENDIF
       if(I.eq.Nblock)THEN      !画某个块体 2016.11.22
        CALL DLINE(LV,X1,Y1,Z1,X2,Y2,Z2)  !!if(y_n.eq.'y'.or.y_n.eq.'Y')
       ENDIF 
            
     	 CALL OGLINE1(LV1,KM2,X1,Y1,Z1) 
870       CONTINUE
          if(Nblock.eq.0.and.vol(I).gt.0.001.or.I.eq.Nblock)then 
           write(19,*)'    endloop'
           write(19,*)'  endfacet'
          endif
880     CONTINUE 
          if(Nblock.eq.0.and.vol(I).gt.0.001.or.I.eq.Nblock)then 
          write(19,'(A22)')'endsolid OuterBoundary'
         endif
          kib=kib+abs(NBKC(I))
900      CONTINUE
        WRITE(LV,225)
        WRITE(LV,230)
        WRITE(LV,225)
        WRITE(LV,240)
        CLOSE (LV)
	  CLOSE(LV1)
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
