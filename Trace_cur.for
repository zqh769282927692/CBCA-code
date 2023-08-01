 
      
C#####################################
C--   平面与曲面之间的交线分析 2018.3.26
!-- 曲面看作圆盘，因曲面的外轮廓点形成的曲面范围考虑为凸形面
!-- 以曲面为母，多边形平面在其中穿越，分析穿越的交点形成的交线情况      
C#####################################
        SUBROUTINE JT_P_Q

        implicit real*8 (a-h,o-z) 
         
	  DIMENSION COTO(1300000,6),COLO(1300000,4),NOLO(1300000,3) 
        DIMENSION TXY(3000,15,3),CC0(3,3),B2(3,3),QJD(5800,3)
        DIMENSION M3V(3000,15,2),MV(8,3),VM(8,3),MVT(150,500,2)
        DIMENSION MTVC0(5),TVC0(5,1500,3),LS(150),kt0(1500)
	  DIMENSION TVC(150,500,3),MTVC(150),X(5000),XY0(5000,2),CC1(3, 3)  
        DIMENSION CVEF(800,2),V0(800,2),V1(800,3),NNEB(3000,3)
        DIMENSION ET(3),EL(3),V3D(5000,10,3),V2D(5000,10,2),P1D(5000,9) 
        DIMENSION VK(3),Cxyz(500,3),coor(2,3),MJFDE1(5000),P1E(100,9)
        DIMENSION V2E(100,50,2),V3E(100,50,3),MLV(1000,50),ML(1000) 
        DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
        DIMENSION CAZOOM(8,3),a0(3),MV0(2),TXYK(2,3)
        
        COMMON/COLTN/COTO,COLO,NOLO,NJT0/V3D/V3D/V2D/V2D/P1D/P1D        
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur/TXY/TXY
        COMMON/NJE0/NJFDE0,NJFDE1/ENL/N0,NL0/MJFDE1/MJFDE1/P1E/P1E
        COMMON/V2E/V2E/CAOM/CAZOOM/DT0/DT0,stp0/VRANGE/VRANGE 
         COMMON/NJD0/NJFRD0/EML/ML,MLV/V3E/V3E !掉了NJFRD0，造成release和debug的结果不同，release能得到正常结果，估计是
                                     !NJFRD0当成0，debug不执行 DO 200循环，但可以算完，得到结果没有曲线  2019.9.17
        
          X12=CAZOOM(2,1)-CAZOOM(1,1)
          Y14=CAZOOM(4,2)-CAZOOM(1,2)
          Z15=CAZOOM(5,3)-CAZOOM(1,3)       
!          VRANGE=(X12+Y14+Z15)/3.0   
          
          CALL GroundSurf            !2018.3.28          

           K12 = 0        !记录曲面与线段交点的总数量,放在求解之外，以比较所有交点是否接近
          DO 300 I1=1,Ncur              !曲面,作为母面。比较时，该段程序是按照曲面在前编写的, 
                                !在此的顺序为曲面、定位结构面、开挖面.在intlin中的顺序是曲面放在后面
        
	    I11=I1+NJFRD0+NJFDE1+NL0   ! I11为NOLO的存储位置，曲面存在最后，照顾以前切割程序的存储
                               !顺序为定位结构面、开挖面、曲面
           NVEF=NEcur(I1)        
         DO K0=1,NVEF
          CVEF(K0,1:2)= CoCur(i1,K0,1:2)           !三角形与曲面求交，如果两个交点都位于曲面外，则不求交
         enddo
           
         DO 200 I20=I1+1,NJFRD0+NJFDE1+NL0+Ncur    !曲面和平面   没有考虑NJFRD0(=0)，什么时候考虑删掉
             k12=0                                 !二分法求得的交点是否重合，只在求取同一条迹线时进行比较，不同迹线不比较 2018.10.20
           write(12,*)'NcurI20',I20,NJFRD0+NJFDE1+NL0+Ncur
         IF(I20.LE.Ncur)THEN  
          I2=I20                 !曲面
          I22=I2+NJFRD0+NJFDE1+NL0  ! I22为NOLO的存储位置，曲面存在最后
          m2=NEcur(I2)            !定义曲面范围的顶点数
         DO K0=1,m2
          V0(K0,1:2)= CoCur(i2,K0,1:2)
          V1(K0,1:3)= CoCur(i2,K0,1:3)
         enddo
         ELSE
          if(I20.LE.Ncur+NJFDE1)then
          I2=I20-Ncur                   !多边形定位结构面
           I22=I2+NJFRD0         
          m2=MJFDE1(I2)                 !面的顶点数
         DO K0=1,m2
           V0(K0,1:2)=V2D(I2,K0,1:2)    !局部坐标，用作三角形剖分    
           V1(K0,1:3)=V3D(I2,K0,1:3)    !对应的整体坐标
         enddo
          else
           I2=I20-Ncur-NJFDE1           !多边形开挖面
           I22=I2+NJFRD0+NJFDE1   
           m2=ML(I2)                      !面的顶点数     
         DO K0=1,m2
           V0(K0,1:2)=V2E(I2,K0,1:2)  !局部坐标，用作三角形剖分
           V1(K0,1:3)=V3E(I2,K0,1:3)  !对应的整体坐标
         enddo    
         endif
        
         ENDIF
           Write(12,*)'面面相交=',I1,i20
 
       CALL Triangulation(V0,m2,NNEB,Ntri)
         m21=m2
            Write(12,*)'三角形划分后顶点数:',Ntri,m2
            
       CALL TriangleFine(V0,m2,DT0,VRANGE,NNEB,Ntri) !三角形网格细分后，V0和m2增大
         m22=m2
           Write(12,*)'三角形细分后顶点数:',Ntri,m2
          
!!!!         if(m21.eq.m22) goto 20  !网格细分时没有增加三角形，也执行一下下面的语句 
                                     !这样就不需要重复 N1=NDcur(i2) 等
           
        if(I20.LE.Ncur)THEN       !曲面
          N1=NDcur(i2)                !子面
          XY0(1:N1-3,1:2)=CoCur(i2,1:N1-3,1:2)  
          X(1:N1)=XX(i2,1:N1)      
         do K0=m21+1,m22        
           Xt=V0(K0,1) 
           Yt=V0(K0,2)  
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt)
           V1(K0,1:2)= V0(K0,1:2)
           V1(K0,3)=Zt          !网格细分时增加的三角形顶点对应的Z坐标
          enddo    
        else                    !多边形平面
          if(I20.LE.Ncur+NJFDE1)then
        DO 6 J1=1,3
        DO 6 J2=1,3
6       B2(J1,J2)=P1D(I2,(J1-1)*3+J2)
          else
        DO 7 J1=1,3
        DO 7 J2=1,3
7       B2(J1,J2)=P1E(I2,(J1-1)*3+J2)       
         endif 
         
        Z10=V1(1,1)*B2(3,1)+V1(1,2)*B2(3,2)+V1(1,3)*B2(3,3)
        EL(3)=Z10
         
        do K0=m21+1,m22        
        EL(1:2)=V0(K0,1:2)          
        CALL LTXYZ(B2,ET,EL)
        V1(K0,1:3)=ET(1:3)         !增加点的整体坐标         
        enddo
       endif
          
       !!!  每个线段中间内插4个点,即5等分''''''''  
        do 50 i=1,Ntri  
        do 40 J = 1,3                               !ML(I)
         K1 = J
         K2 = K1 + 1
        If (K2.gt.3) K2 = 1 
 
        TXY(I,(J-1)*5+1,1:3)=V1(NNEB(i,j),1:3)     !三角形的3个顶点坐标 
        
        do 30 M1 = 1,4                          
          VK(1:3)=V1(NNEB(i,k1),1:3)+
     #     (V1(NNEB(i,k2),1:3)-V1(NNEB(i,k1),1:3))*real(M1)/5.0
                                     
        if(I20.LE.Ncur)THEN       ! 内插的点也准确求得曲面上的Z坐标
           Xt=VK(1) 
           Yt=VK(2)  
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt) !用前面针对子面定义的N1,XY0,X
           VK(3)=Zt        
        ENDIF  
         TXY(I,(J - 1) * 5 + 1 + M1,1:3)= VK(1:3)  
30      continue 
  
        M3V(I,(J-1)*5+1:(J-1)*5+5,1)=NNEB(i,k1) !10等分后1~5，6~10，11~15线段所在的三角形的2个顶点 
        M3V(I,(J-1)*5+1:(J-1)*5+5,2)=NNEB(i,K2) !用于后面比较线段与曲面的交点是否同一点 
         
40      continue  
         ict=5;icc=5
         do   j1=1,15 
               j2=j1+1
          if(j1.eq.15) j2=1   
          coor(1,1:3)=TXY(i,j1,1:3) 
          coor(2,1:3)=TXY(i,j2,1:3) 
        if(I20.eq.6) write(111,'(2i6,2x,3f10.4)')i1,i,TXY(i,j1,1:3) 
        if(I20.eq.6) call DLINE(21,coor,ict,icc)  !控制哪个面上显示三角形网格        
        enddo   
        
50    continue      
         
!!!!---   
          N1=NDcur(i1)                    !母面              
          XY0(1:N1-3,1:2)=CoCur(i1,1:N1-3,1:2)  
          X(1:N1)=XX(i1,1:N1)      
          
          KQ0=0                    !曲交线数量，面I20所有三角形形成的曲交线总数
         
         do 160 i=1,Ntri 
 
          do 80 j1=1,15               !三角形扩充后为15个顶点
           j2=j1+1
          if(j1.eq.15) j2=1        
          
          Xt1=TXY(i,j1,1); Yt1=TXY(i,j1,2)
        Call CurvPlanePoint(N1,XY0,X,Xt1,Yt1,Zt1)
          Xt2=TXY(i,j2,1); Yt2=TXY(i,j2,2)
        Call CurvPlanePoint(N1,XY0,X,Xt2,Yt2,Zt2) 
        
         If((Zt1-TXY(i,j1,3))*(Zt2-TXY(i,j2,3)).lt.-1.0e-6) Then  ! ''线段的两点在曲面的不同侧        
          GoTo 101
         End If
80       continue
        goto 160
        
101      K03 = 0           ! '''在一个三角形内，累计线段与曲面的交点数
         M6 = 0            ! '''累计顶点串的数量       
   
! 得到K1,K2 /// 第1点位于曲面之下,第2点位于曲面之上,此时的第1点作为起始点.这样得到的顶点串是先出现
                   !位于曲面上的串,再出现位于曲面下的串.即奇次串为曲面上的串,偶次串为曲面下的串
    
          do 110 j=1,15            
         
           If(J.eq.1) Then   !一开始进入循环时，定位K1，K2
          K1 = 1
          K2 = K1+1
          K11 = K1; K22 = K2 
          
201       Xt= TXY(i,K11,1); Yt=TXY(i,K11,2)
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt1)
          Xt= TXY(i,K22,1); Yt=TXY(i,K22,2)
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt2)
       If(Zt1-TXY(i,K11,3).lt.-1.0e-6.or.Zt2-TXY(i,K22,3).gt.1.0e-6)Then
                                                  !     要求第1点位于曲面之下,第2点位于曲面之上
                 !M6.eq.4时，发现第一点位于面下或面上，均不能保证得到正确的两段曲交线,此时将三角形进一步划小 2021.11.19
                 !CurSufInp55-3-1.txt得到重复的曲交线，也是这个原因，也有M6.eq.4，画图发现三角形棱边和母面接近相切。
                 ! 改小三角形边长后ok。   2021.12.1 
           K1 = K1 + 1     
           K2 = K1 + 1
          If (K1.le.15) then
            K11 = K1  
          else
            K11 = K1 - 15   
          endif   
          If (K2.le.15)then
            K22 = K2
         else
            K22 = K2 - 15
          endif    
         GoTo 201
       End If
        
       Else                !K1，K2累加，找下一个不同侧的点位
          K1 = K1 + 1          
          K2 = K1 + 1  
          
           If (K1.le.15) then
            K11 = K1  
          else
            K11 = K1 - 15   
          endif   
          If (K2.le.15)then
            K22 = K2
          else
            if(k2.gt.30)then
                K22 = K2 - 30   
            else
                K22 = K2 - 15
            endif    
          endif    
        End If
        
          Xt= TXY(i,K11,1); Yt=TXY(i,K11,2)
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt1)
           
          Xt= TXY(i,K22,1); Yt=TXY(i,K22,2)         
          Call CurvPlanePoint(N1,XY0,X,Xt, Yt, Zt2) 
         
       If((Zt1- TXY(i,K11,3)) * (Zt2 - TXY(i,K22,3)).lt.-1.0e-8) Then  !  ''线段的两点在曲面的不同侧
            
          CC0(1, 1:3) = TXY(i,K11,1:3)  
          CC0(2, 1:3) = TXY(i,K22,1:3) 
          CC1=CC0                          !ErFen Fa中cc0会随时改变，先把cc0赋值给CC1，用后者进入ErFen Fa子程序
!            write(12,*)'CC1=',CC0(1, 1:3) 
!            write(12,*)'CC2=',CC0(2, 1:3)
         Call ErFenFa(N1,XY0,X,CC1,QJD,k12,0)     ! QJD保存不同面及其不同三角形形成的所有交点坐标 2018.3.30 
                                                ! 以进行距离判断，当距离很近时进行捏点   
         K03 = K03 + 1
         cc0(3,1:3)=cc1(3,1:3)
         VM(K03, 1:3) = CC1(3, 1:3)         ! 暂时认为: K03的最大值不超过8,即一个三角形与曲面最多只有8个交点        
        
!         if(I20.LE.Ncur)THEN
!          a0(1:3)=VM(K03, 1:3)
!          TXYK(1,1:3)= TXY(i,K11,1:3)
!          TXYK(2,1:3)= TXY(i,K22,1:3)
!          write(12,*)'ij=',i,j
!          call JiaoDianXiuZheng(i1,i2,a0,TXYK)
!          VM(K03,1:3)=a0(1:3)
!         endif  

          M6 =K03         ! 每次遇到点K1,K2位于不同侧时, 重新记录1条顶点串.顶点串的数量==曲面与线段的交点数
          
          MV(K03,1:2)=M3V(I,K11,1:2)  !该交点所属线段的所属三角形的两个端点总体编号,不是J，是k11  2018.9.10 
         if(i1.eq.-1.and.i20.eq.6)then
          write(12,*)'MV=',MV(K03,1:2),K11,K22
           write(12,*)'CC1=',CC0(1, 1:3) 
           write(12,*)'CC2=',CC0(2, 1:3)
           write(12,*)'CC3=',CC0(3, 1:3)
           do jj=1,15
            write(12,*)'TXY',TXY(i,jj,1:3)
           enddo 
         endif    
         End If    
110      continue   
 
!----------  获得曲交线   ---------      
      
       if(M6.eq.4)then
         write(12,*)'sandian=',TXY(i,1,1:3), TXY(i,6,1:3), TXY(i,11,1:3)
        do J1 = 1,M6 
          write(12,*) 'MV=',MV(J1,1:2)
        enddo 
        L1=0
        if(MV(1,1).eq.MV(2,1).and.MV(1,2).eq.MV(2,2))L1=1     !不能出现1和2，或3和4在三角形同一边上，否则连线就会
        if(MV(3,1).eq.MV(4,1).and.MV(3,2).eq.MV(4,2))L1=1     !位于三角形外部 2021.12.8
        if(L1.eq.1)then
         MV0(1:2)=MV(1,1:2)             
          a0(1:3)=VM(1,1:3)          !--为了保证都是从“下上”连“上下”，要进行交换
          MV(1,1:2)=MV(3,1:2)        !--第3行传给第1行     2021.12.8
          VM(1,1:3)=VM(3,1:3)  
         MV(3,1:2)=MV0(1:2)          !--第1行传给第3行
         VM(3,1:3)=a0(1:3)
        endif    
         do J1 = 1,M6 
          write(12,*) 'MV0=',MV(J1,1:2)
        enddo 
       endif
       
       do 120 J0 = 1, M6 / 2
         
        KT1 = (J0 - 1) * 2 + 1
        KT2 = (J0 - 1) * 2 + 2      
        
!!       Call DingDianLianJie 该调用直接写在此       
!!    先记录第1个交点,再记录曲线点,再记录第2个交点。通过2个交点,求得中间的连接曲线 
          Xn0 = VM(KT1, 1)        ! ''曲线从第1个交点到第2个交点,第3个交点到第4个交点, ...的连接顺序
          Yn0 = VM(KT1, 2)        !   没有笔记中的‘3-4超长曲线’情况
          Zn0 = VM(KT1, 3)
          XnT = VM(KT2, 1)
          YnT = VM(KT2, 2)
          ZnT = VM(KT2, 3) 
!             write(12,*)'vm1',VM(KT1, 1:3),i  
!             write(12,*)'vm2',VM(KT2, 1:3) 
          CALL POINTC(Xn0,Yn0,NVEF,CVEF,IW1)  
          CALL POINTC(XnT,YnT,NVEF,CVEF,IW2)    
          if(IW1.EQ.1.AND.IW2.EQ.1) goto 120   !2个点都位于域外  2018.9.12
          
          Call CurveLine(Xn0,Yn0,Zn0,XnT,YnT,ZnT,i1,i20,Cxyz,K10)
                    !-- i1,i20传给CurveLine再传给XpieYpie
          M4 = 0 
          KQ0 = KQ0 + 1             ! 曲交线数量+1
         M4 = M4 + 1
         TVC(KQ0,M4,1:3)=VM(KT1,1:3)      !曲交线的第1点的坐标是第一个交点  
         MVT(KQ0,M4,1:2)=MV(KT1,1:2)     !该交点所属线段的所属三角形的两个端点总体编号   
         
          if(I1.eq.6.and.i20.eq.11) write(12,*) 'TVC1=',TVC(KQ0,M4,1:3)
         if(K10.eq.0)goto 149          !2018.11.14
         do IK = 1 , K10
         M4 = M4 + 1
         TVC(KQ0, M4, 1:3) = Cxyz(IK, 1:3)   !  TVC的第一维不是K01，用KQ0表示每段曲交线单独保存     
          if(I1.eq.6.and.i20.eq.11) write(12,*) 'TVC=',TVC(KQ0,M4,1:3)
        enddo 
  
149      M4 = M4 + 1
         TVC(KQ0,M4,1:3)=VM(KT2,1:3)      !曲交线的最后1点的坐标是第二个交点 
         MVT(KQ0,M4,1:2)=MV(KT2,1:2)
         if(I1.eq.6.and.i20.eq.11) write(12,*) 'TVC2=',TVC(KQ0,M4,1:3) !2021.12.1  
         MTVC(KQ0)=M4  
120     continue
                      
160      continue  
           write(12,*) 'KQ0=',KQ0
        if(KQ0.eq.0) goto 200   !跳出，否则在无交线时也会执行下面语句 2018.4.21
!!------ 因子面划分为多个三角形，得到的曲线段需要合并 
!!-- 参照流形法后处理程序DwgTec.for中的回路合并（来源于四面体剖分）
         
         LS(1:KQ0)=1   
         LS(1)=0
         KQ1=1              
         MTVC0(KQ1)=MTVC(1)       !该交线的点数
        TVC0(KQ1,1:MTVC0(KQ1),1:3)= TVC(1,1:MTVC(1),1:3)  
         N1=MVT(1,1,1); N2=MVT(1,1,2)
         N3=MVT(1,MTVC(1),1); N4=MVT(1,MTVC(1),2)  
         
890       I8=0         
         do 900 j0=2,KQ0
         if(LS(j0).eq.0) goto 900
         
         N5=MVT(j0,1,1); N6=MVT(j0,1,2)
         N7=MVT(j0,MTVC(j0),1); N8=MVT(j0,MTVC(j0),2)                 
               
         if(N3.eq.N6.and.N4.eq.N5)then      !线段a--b,c--d 的b，c相连时
          
           do k0=1,MTVC(j0)-1
          TVC0(KQ1,MTVC0(KQ1)+k0,1:3)= TVC(j0,1+k0,1:3)  !j0线段从第2个点开始并入          
           enddo   
            MTVC0(KQ1)= MTVC0(KQ1)+MTVC(j0)-1
            N3=MVT(j0,MTVC(j0),1); N4=MVT(j0,MTVC(j0),2)
           LS(j0)=0           
             I8=1     !存在可以连接的线段，就继续循环       
         endif    
          if(N1.eq.N8.and.N2.eq.N7.and.LS(j0).ne.0)then      !线段a--b,c--d 的b，c相连时
                                              !环形交线，防止接在左端也可以接在右端 2018.9.13
            TVC0(KQ1,1+MTVC(j0)-1:MTVC0(KQ1)+MTVC(j0)-1,1:3)
     #                             = TVC0(KQ1,1:MTVC0(KQ1),1:3)   !数组后移 
            do k0=1,MTVC(j0)-1  
             TVC0(KQ1,k0,1:3)= TVC(j0,k0,1:3)  ! 然后j0线段接入左边      
            enddo   
               MTVC0(KQ1)= MTVC0(KQ1)+MTVC(j0)-1  
              N1=MVT(j0,1,1); N2=MVT(j0,1,2)  
             LS(j0)=0           
             I8=1     !存在可以连接的线段，就继续循环        
          endif     
          
        do k0=1,KQ1
        do k1=1, MTVC0(k0)  
!           write(12,*) TVC0(k0,k1,1:3)         
        enddo
        enddo
900     continue             
             
         IF(I8.EQ.1)GOTO 890     
         
         do 950 j0=2,KQ0
          IF(LS(j0).EQ.1) THEN    !线段合并时，可能形成多段线段（面面相交形成多段线段）
          KQ1=  KQ1+1  
          MTVC0(KQ1)=MTVC(j0)
        TVC0(KQ1,1:MTVC0(KQ1),1:3)= TVC(j0,1:MTVC(j0),1:3)  
          LS(j0)=0
         N1=MVT(j0,1,1); N2=MVT(j0,1,2)
         N3=MVT(j0,MTVC(j0),1); N4=MVT(j0,MTVC(j0),2)             
           GOTO 890
          ENDIF
950      CONTINUE 
           write(12,*) 'KQ11=',KQ1
!!---------前面有捏点，造成KQ0中的线段首尾坐标相同，并带入线段合并中，造成回路分析时因重合点存在
!!-------  而认为是负回路。少了回路造成块体没被切开。需要删除线段中的坐标重合点   2018.5.8    
!!----     将短线段合并 --    2018.11.15  
!          
         do 980 j0=1,KQ1
           kt0(1:1500)=1  
969         M0=0
         do 970 j2=2,MTVC0(j0)-1 
          IF(kt0(j2).EQ.0)GOTO 970     !该点被删除时，前移一点   2018.11.15
           
           j1=j2-1
965         M21=0
           IF(kt0(J1).EQ.0.AND.J1.GT.1)THEN
            j1=J1-1          !继续向后找点    2018.11.15
            M21=1
           ENDIF  
           IF(M21.EQ.1)GOTO 965     
           
           j3=j2+1          !只要一次向前，2021.9.27
!966         M21=0
!           IF(kt0(J3).EQ.0.AND.J3.LT.MTVC0(j0))THEN
!            j3=J3+1          !继续向前找点
!            M21=1
!           ENDIF  
!           IF(M21.EQ.1)GOTO 966       
           
          d1=(TVC0(j0,j1,1)-TVC0(j0,j2,1))**2+
     #       (TVC0(j0,j1,2)-TVC0(j0,j2,2))**2+     
     #       (TVC0(j0,j1,3)-TVC0(j0,j2,3))**2
          d1=sqrt(d1)
          d2=(TVC0(j0,j2,1)-TVC0(j0,j3,1))**2+
     #       (TVC0(j0,j2,2)-TVC0(j0,j3,2))**2+     
     #       (TVC0(j0,j2,3)-TVC0(j0,j3,3))**2
          d2=sqrt(d2)
!           if(j0.eq.1.and.i1.eq.3.and.i20.eq.13)then
!          write(1001,*)'jm', j0,MTVC0(j0),j2,j1,j3,d1
!          write(1001,*)TVC0(j0,j2,1:3)  
!           endif 
           
          if(d1.lt.Stp0*0.1.or.d2.lt.Stp0*0.1.or.d1+d2.lt.Stp0*0.3)then
                                                  !  可改变线段延长的程度 <调试入口6>
              kt0(j2)=0    
              M0=1              
!              GOTO 971      !只判断一个点就跳出，因其前后的点可能被删除，因此不通过循环后移
          endif           
970      CONTINUE 
         
!!!971      IF(M0.EQ.1)GOTO 969
         
          m1=0
          do 975 j1=1,MTVC0(j0)    
         if (kt0(j1).eq.1)then
           m1=m1+1
           TVC0(j0,m1,1:3)=TVC0(j0,j1,1:3)
        endif
975       CONTINUE           
          MTVC0(j0)=m1                
          
980      CONTINUE  
!!-----曲交线的各点是否都位于母面之外 2018.5.15          
         
!!----- 画出曲交线 --      
981        ict=i1;icc=i20
           if(ict.gt.9)ict=9
           if(icc.gt.9)icc=9
        t01=VRANGE
        t02=-1
        do j0=1,KQ1                            ! KQ1为两个面之间的曲交线数量
        do j1=1, MTVC0(j0)-1
       if(I1.eq.-3.and.i20.eq.6) write(1001,*)'j0j1=',j0,j1,MTVC0(j0)
           j2=j1+1 
         coor(1,1:3)=TVC0(j0,j1,1:3) 
         coor(2,1:3)=TVC0(j0,j2,1:3)     
         ict=9;icc=9
           if(I1.eq.3.and.i20.eq.6)Call DLINE(21,coor,ict,icc)   
          
          t1=(coor(1,1)-coor(2,1))**2+(coor(1,2)-coor(2,2))**2
          t1=sqrt(t1+(coor(1,3)-coor(2,3))**2)
          if(t1.lt.0.1)then              
              write(1001,*)coor(1,1:3)
              write(1001,*)coor(2,1:3)
              write(1001,*)'t1=',t1 
          endif
         if(t1.gt.t02) t02=t1
         if(t1.lt.t01) t01=t1
          if(I1.eq.-3.and.i20.eq.6) write(1001,*) TVC0(j0,j1,1:3)      ! 曲交线的坐标
        enddo
         if(I1.eq.-3.and.i20.eq.6) write(1001,*) TVC0(j0,MTVC0(j0),1:3) 
 
        enddo
!!----- 再对曲交线进行存储            
 
        do 1000 j0=1,KQ1
        do 1000 j1=1,MTVC0(j0)-1
         j2=j1+1 
           
        NJT0=NJT0+1
        
        DO 51 K2=1,2
51      COLO(NJT0,K2)=TVC0(j0,j1,k2)   !母面是曲面，三维坐标的前两维就是局部坐标
        DO 52 K2=1,2
52      COLO(NJT0,2+K2)=TVC0(j0,j2,k2)  !Y(M2,K2)
    
        DO 54 K3=1,3
          COTO(NJT0,K3)=TVC0(j0,j1,k3)    !ET(K3)
54        COTO(NJT0+1,K3)=TVC0(j0,j1,k3)  !ET(K3)
     
       if(I20.LE.Ncur)then    !子面是曲面
        COLO(NJT0+1,1:2)=TVC0(j0,j1,1:2)    
       else                   !子面是平面       
         ET(1:3)=TVC0(j0,j1,1:3)
          CALL TLXYZ(B2,ET,EL)
          DO 55 K2=1,2
55        COLO(NJT0+1,K2)=EL(K2)    
       endif 
        DO 57 K3=1,3
          COTO(NJT0,3+K3)=TVC0(j0,j2,k3)   
57      COTO(NJT0+1,3+K3)=TVC0(j0,j2,k3)  
!          if(I20.eq.3)then
!            write(12,*)'COTO',COTO(NJT0,1:3)
!        endif
       if(I20.LE.Ncur)then    !子面是曲面
        COLO(NJT0+1,3:4)=TVC0(j0,j2,1:2)      
       else
        ET(1:3)=TVC0(j0,j2,1:3)
	   CALL TLXYZ(B2,ET,EL)
        DO 58 K2=1,2
58      COLO(NJT0+1,2+K2)=EL(K2)  
        endif
           
        NOLO(NJT0,1)=I11
        NOLO(NJT0,2)=I22
        NOLO(NJT0,3)=NJFRD0+NJFDE1+NL0+Ncur+j1  !该曲交线的第j1段,加上NJFRD0+NJFDE1+NL0.     2018.4.22
                                           !否则可能和I11或I22相等，造成在求总体点号的3个面组成比较时出问题       
        !该问题很隐蔽，删除回路前后，CNODTZ的坐标没问题，但通过NLOC进行的画图有问题，因此定位是求总体点号时出问题
        NOLO(NJT0+1,1)=I22           
        NOLO(NJT0+1,2)=I11
        NOLO(NJT0+1,3)=NJFRD0+NJFDE1+NL0+Ncur+j1    !该曲交线的第j1段 !+Ncur,易错处  2018.5.5
         NJT0=NJT0+1    
         
1000    continue
      
200     CONTINUE               
300     CONTINUE
        call DXFEND(21)
        close(21)
       RETURN
       END       
 
C#####################################   
       SUBROUTINE CurveLine(Xn0,Yn0,Zn0,XnT,YnT,ZnT,i1,i20,Cxyz,K10)
C#####################################      
       implicit real*8 (a-h,o-z)  
       DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
       DIMENSION XY0(5000,2),X(5000), Cxyz(500, 3)
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur
        COMMON/DT0/DT0,stp0
       
         
         If((Xn0-XnT)**2+(Yn0-YnT)**2.lt.(2*Stp0)**2)Then  !起点与终点接近 , 就不求解曲线, 直接取其平分点
            K10 = 0       !1改为0，此时不在其中平分，否则会引起很短的线段，造成intlin中线段求交点错误及随后的回路错误
                          !2018.11.14 
!             Cxyz(K10, 1) = (Xn0 + XnT) / 2.0
!             Cxyz(K10, 2) = (Yn0 + YnT) / 2.0
!             Cxyz(K10, 3) = (Zn0 + ZnT) / 2.0
           GoTo 301
         End If
    
!''''''''''''''''''''''''''先进行公式符号的判断''''''''''''''''
          NK = 0          !  ''' 标记公式中的正负号,为0是为负号,为1是为正号
    
          Call XpieYpie(i1,i20,Xn0,Yn0,Xpie,Ypie)
!!!           write(12,*)'Xn0', Xn0,Yn0,XnT,YnT,Stp0
          Xn1 = Xn0 - 2 * Stp0 * Xpie
          Yn1 = Yn0 - 2 * Stp0 * Ypie
      
      If((Xn1-XnT)**2 +(Yn1-YnT)**2.gt.(Xn0-XnT)**2+(Yn0-YnT)**2)Then ! 一般,距离的变化为单调的
            NK = 1
        End If
 
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
           K10 = 0        ! '''曲线的逼近点数
  
101        K10 = K10 + 1
          Call XpieYpie(i1,i20,Xn0,Yn0,Xpie,Ypie)
!           write(*,*)k10
           If (NK.eq.0) Then
             Xn1 = Xn0 - Stp0 * Xpie
             Yn1 = Yn0 - Stp0 * Ypie
            Else
             Xn1 = Xn0 + Stp0 * Xpie
             Yn1 = Yn0 + Stp0 * Ypie
           End If
!!             write(12,*)'Xn1',Xn1,Yn1
         Xpie1 = Xpie; Ypie1 = Ypie
  
201       Call XpieYpie(i1,i20,Xn1,Yn1,Xpie,Ypie)
  
         Xpie2 = Xpie; Ypie2 = Ypie
         
        If (NK.eq.0) Then
         Xn2 = Xn0 - Stp0 * (Xpie1 + Xpie2) / 2.0          ! ''' 梯形公式进行迭代
         Yn2 = Yn0 - Stp0 * (Ypie1 + Ypie2) / 2.0
        Else
          Xn2 = Xn0 + Stp0 * (Xpie1 + Xpie2) / 2.0
          Yn2 = Yn0 + Stp0 * (Ypie1 + Ypie2) / 2.0
        End If
!          write(*,*)'Xn2',Xn2,Yn2
        If (abs(Xn2-Xn1).gt.0.001.or.abs(Yn2-Yn1).gt.0.001) Then  !0.00001改为0.001
          Xn1 = Xn2
          Yn1 = Yn2
          GoTo 201
        End If
          N1=NDcur(i1)                    !代入母面方程              
          XY0(1:N1-3,1:2)=CoCur(i1,1:N1-3,1:2)  
          X(1:N1)=XX(i1,1:N1)
        
         Call CurvPlanePoint(N1,XY0,X,Xn2, Yn2, Zn2)
!  ''''Zn2 = Cp1 * Xn2 + Cp2 * Yn2 + Cp3  直立面时，出错。不用该式求Zn2 ,
                                   ! '''' 会出现曲线上的点不完全在一个面上,通过图形视觉可以看出误差    
          Xn0 = Xn2
          Yn0 = Yn2    
          Cxyz(K10, 1) = Xn2; Cxyz(K10, 2) = Yn2; Cxyz(K10, 3) = Zn2 
         if(abs(Yn2-80).lt.0.01)write(12,*)'Yn2',stp0,Xn2,Yn2
      If((Xn2-XnT)**2+(Yn2-YnT)**2.gt.(2.0*Stp0)**2.and.K10.lt.499)Then 
                                                      !此处很关键 2021.9.12
                    GoTo 101                 !'与终点较为接近时,停止求解.
          endif                              !''''根据步长,最多进行200次求解?曲线可能很长

301       continue
!          write(12,*) '该交线的交点数量:', K10
           do IK = 1,K10
!          write(12,*) IK, Cxyz(IK, 1:3)
           enddo

        End  
C#####################################  
         SUBROUTINE XpieYpie(i1,i2,Xn,Yn,Xpie,Ypie)
C#####################################  
          implicit real*8 (a-h,o-z)    
     
        DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
        DIMENSION Qx(3),Qy(3),QN(3),Px(3),Py(3)
       DIMENSION CC1(3,3),C2(3,3),P2D(5000,5),P2E(100,5)
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur/P2D/P2D/P2E/P2E 
        COMMON/NJE0/NJFDE0,NJFDE1
!!!---  获得子面的方程参数 --        
         IF(I2.LE.Ncur)THEN          
       
        call PxPy(i2,Xn,Yn,Qx,Qy)       !2018.5.16
        
         C2(1, 1:3) = Qx(1:3)
         C2(2, 1:3) = Qy(1:3)
   
         CALL MUTIVT(C2,C1)     
         QN(1:3) = C2(3, 1:3)            ! 子面（曲面、平面）的单位法矢  
         
        else                             !子面为平面。改变原来通过三角形3个点算法矢的做法 2021.10.5
          if(I2.LE.Ncur+NJFDE1)then
            I21=I2-Ncur                   !多边形定位结构面                
          QN(1:3)=P2D(I21,1:3)            !子面的单位法矢  
          else
           I21=I2-Ncur-NJFDE1           !多边形开挖面
           QN(1:3)=P2E(I21,1:3)          
         endif
      
        endif        
         
!!!----参考书中的公式（3-21）（3-22） 
         call PxPy(i1,Xn,Yn,Px,Py)       !2018.5.16
        
!!!--------
         Xpie = Py(1) * QN(1) + Py(2) * QN(2) + Py(3) * QN(3)
         Xpie = -Xpie
         Ypie = Px(1) * QN(1) + Px(2) * QN(2) + Px(3) * QN(3)
           
        End   
!######################################        
      Subroutine PxPy(i2,Xn,Yn,Qx,Qy)    ! 曲面的x，y向矢量 
!###################################### 
       implicit real*8 (a-h,o-z)    
     
        DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
        DIMENSION Qx(3), Qy(3)   
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur 
        
         N1=NDcur(i2)                !母面
          Qx(1) = 1                  
          Qx(2) = 0
          Qx(3) = XX(i2,2)  
         do J = 1 , N1 - 3
          R0 = (Xn -CoCur(i2,J,1))** 2 + (Yn - CoCur(i2,J,2))** 2
          If (R0.lt.0.0000001) R0 = 0.00001
          Qx(3)=Qx(3)+2*XX(i2,J+3)*(Xn-CoCur(i2,j,1))*(Log(R0)+1)
         enddo 
         
          Qy(1) = 0
          Qy(2) = 1
          Qy(3) =XX(i2,3)
         do J = 1 , N1 - 3
          R0 = (Xn - CoCur(i2,J,1))** 2 + (Yn - CoCur(i2,J,2))** 2
          If (R0.lt.0.0000001) R0 = 0.00001
          Qy(3)=Qy(3)+2*XX(i2,J+3)*(Yn-CoCur(i2,J,2))*(Log(R0)+1)
         enddo
        
         endsubroutine

!############################### 2018.9.8       
       Subroutine JiaoDianXiuZheng(i1,i2,a0,TXYK)       
!###############################     
      implicit real*8 (a-h,o-z) 	 
       DIMENSION X(5000),XY0(5000,2),Qx(3),Qy(3),QN(3),C2(3,3)
       DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
       DIMENSION TXYK(2,3),B2(3,3),ET(3),EL(3),QNL(2) 
       DIMENSION a0(3),a1(3),a2(3),a3(3),CC0(3, 3), QJD(5800, 3)
       COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur 
       
!-- 穿过TXYK(1,)、TXYK(2,)的直立面          
          x21=TXYK(2,1)- TXYK(1,1) 
          y21=TXYK(2,2)- TXYK(1,2)
           write(12,*)'TXYK1=',TXYK(1,1:3)
            write(12,*)'TXYK2=',TXYK(2,1:3)
           write(12,*)'a0=',a0(1:3)
          EL(1)= y21/sqrt(y21**2+x21**2)  !直立面的法矢
          EL(2)=-x21/sqrt(y21**2+x21**2)
          EL(3)=0
!           z10=-(-EL(1)*TXYK(1,1)-EL(2)*TXYK(1,2)-EL(3)*TXYK(1,3))
!          write(12,*)' z10=',z10
          CALL PROJTG(EL,B2)  
          
           n5=1
          Ix1=i1
          Ix2=i2
          mm0=1
!--  求得点a0在母面(i1)（或子面(i2)）上的法矢
999        Xn= a0(1)
           Yn= a0(2)
           call PxPy(Ix1,Xn,Yn,Qx,Qy)  
           C2(1, 1:3) = Qx(1:3)
           C2(2, 1:3) = Qy(1:3)   
           CALL MUTIVT(C2,C1)   
           QN(1:3) = C2(3, 1:3)      ! 点a0在子面的单位法矢  
!          IF(Zt1.GT.a0(3)) QN(1:3)=-QN(1:3)   !置反向
          write(12,*)'QN=',QN(1:3)
          QNL(1)=QN(1)*B2(1,1)+QN(2)*B2(1,2)+QN(3)*B2(1,3)  
          QNL(2)=QN(1)*B2(2,1)+QN(2)*B2(2,2)+QN(3)*B2(2,3)
          t1=sqrt(QNL(1)**2+ QNL(2)**2)
            QNL(1)=QNL(1)/t1          !QN在临时坐标系下的矢量
            QNL(2)=QNL(2)/t1
           write(12,*)'QNL=',QNL(1:2)
           
           CALL TLXYZ(B2,a0,EL)   !点a0转到临时坐标系 
               write(12,*)' EL=', EL(3)
            z10=EL(3)      
            za=QNL(1)
            zb=QNL(2)
            zc=-za*EL(1)-zb*EL(2)
             
99        if(mm0.eq.1)then
               a1(1)=EL(1)+abs(x21)*2.0
           else
               a1(1)=EL(1)-abs(x21)*2.0         !射线方向改为反向  
           endif
           
        
          a1(2)=(-za*a1(1)-zc)/(zb+1.0e-10)  
            a1(3)=z10
            CALL LTXYZ(B2,ET,a1)   !得到的切线延长线上点a1的整体坐标ET 
          
            CC0(1,1:3) = a0(1:3)  
            CC0(2,1:3) = ET(1:3)   
            
            N2=NDcur(Ix2)               ! 点a1与子面(i2)（或母面(i1)）的关系
            XY0(1:N2-3,1:2)=CoCur(Ix2,1:N2-3,1:2)  
            X(1:N2)=XX(Ix2,1:N2)   
           
            if(n5.eq.1)then             !第1次求采用射线求交点时，需要判断射线的方向
            Cx1 = CC0(1, 1)
            Cy1 = CC0(1, 2)
            Cz1 = CC0(1, 3) 
             Call CurvPlanePoint(N2,XY0,X,Cx1, Cy1, Zt1) 
            Cx2 = CC0(2, 1)
            Cy2 = CC0(2, 2)
            Cz2 = CC0(2, 3) 
             Call CurvPlanePoint(N2,XY0,X,Cx2, Cy2, Zt2) 
              write(12,*)'Zt1=',Cx1,Cy1, Cz1,Zt1
              write(12,*)'Zt2=',Cx2,Cy2, Cz2,Zt2
               write(12,*)(Zt1 - Cz1) * (Zt2 -Cz2)
            if ((Zt1 - Cz1) * (Zt2 -Cz2) .gt.1.0e-7) Then   !位于曲面的同侧
!!               a1(1)=EL(1)-abs(x21)*2    !射线方向改为反向  
              mm0=2   
               goto 99
            endif
            endif
            
             Call ErFenFa(N2,XY0,X,CC0,QJD,k12,0)  ! a0--a1射线与子面(i2)（或母面(i1)）相交  k12?????????
              a2(1:3)=CC0(3,1:3)                   ! 求得点a2坐标    
             t1=sqrt((a0(1)-a2(1))**2+(a0(2)-a2(2))**2+
     #       (a0(3)-a2(3))**2)
            write(12,*)'a2=',a2(1:3)      
            write(12,*)'t1=',t1
              if(t1.gt.5.0e-5.and.n5.le.5)then
              
                Ix0=Ix1       !交换         
                Ix1=Ix2
                Ix2=Ix0
                 a0(1:3)=a2(1:3)     !点a2的整体坐标传到点a0，并goto 999    
                 n5=n5+1 
               goto 999  
              else
             write(12,*)'n5=',n5     
              endif 
              
         endsubroutine
           
        
!############################### 2018.5.16       
       Subroutine JiaoDianXiuZheng1(i1,i2,a0,TXYK)       
!###############################     
      implicit real*8 (a-h,o-z) 	 
       DIMENSION X(5000),XY0(5000,2),Qx(3),Qy(3),QN(3),C2(3,3)
       DIMENSION CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100) 
       DIMENSION TXYK(2,3),B2(3,3),ET(3),EL(3),QNL(2) 
       DIMENSION a0(3),a1(3),a2(3),a3(3)
       COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur 
       
!-- 穿过TXYK(1,)、TXYK(2,)的直立面          
          x21=TXYK(2,1)- TXYK(1,1) 
          y21=TXYK(2,2)- TXYK(1,2)
           write(12,*)'TXYK1=',TXYK(1,1:3)
            write(12,*)'TXYK2=',TXYK(2,1:3)
           write(12,*)'a0=',a0(1:3)
          EL(1)= y21/sqrt(y21**2+x21**2)  !直立面的法矢
          EL(2)=-x21/sqrt(y21**2+x21**2)
          EL(3)=0
           z10=-EL(1)*TXYK(1,1)-EL(2)*TXYK(1,2)-EL(3)*TXYK(1,3)
          CALL PROJTG(EL,B2) 
!            write(12,*)'b2',b2(1,1:3)
!            write(12,*)b2(2,1:3)
!            write(12,*)b2(3,1:3)
          n5=1
          Ix=i2                       ! i2为子面
100       N1=NDcur(Ix)                ! 求得点a1坐标    
          XY0(1:N1-3,1:2)=CoCur(Ix,1:N1-3,1:2)  
          X(1:N1)=XX(Ix,1:N1)                       
          Cx1 = a0(1)
          Cy1 = a0(2)           
         Call CurvPlanePoint(N1,XY0,X,Cx1, Cy1, Zt1)
          a1(1)=a0(1)
          a1(2)=a0(2)     
          a1(3)=Zt1
            write(12,*)'a1=',a1(1:3)
!--  求得点a1在子面(i2)或母面(i1)上的法矢
           Xn= a1(1)
           Yn= a1(2)
           call PxPy(Ix,Xn,Yn,Qx,Qy)  
           C2(1, 1:3) = Qx(1:3)
           C2(2, 1:3) = Qy(1:3)   
           CALL MUTIVT(C2,C1)   
           QN(1:3) = C2(3, 1:3)      ! 点a1在子面的单位法矢  
          IF(Zt1.GT.a0(3)) QN(1:3)=-QN(1:3)   !置反向
          write(12,*)'QN=',QN(1:3)
          QNL(1)=QN(1)*B2(1,1)+QN(2)*B2(1,2)+QN(3)*B2(1,3)  
          QNL(2)=QN(1)*B2(2,1)+QN(2)*B2(2,2)+QN(3)*B2(2,3)
          t1=sqrt(QNL(1)**2+ QNL(2)**2)
            QNL(1)=QNL(1)/t1          !QN在临时坐标系下的矢量
            QNL(2)=QNL(2)/t1
           write(12,*)'QNL=',QNL(1:2)
            H0=abs(Zt1-a0(3))
            zL0=abs(H0*QNL(2))
            zLx=zL0*QNL(1)
            zLy=zL0*QNL(2)
            dx=-zLx  
            dy=-H0-zLy  
             write(12,*)H0,zL0,dx,dy
             write(12,*)'s',a1(1:3)
           CALL TLXYZ(B2,a1,EL)   !点a1转到临时坐标系
             
            a2(1)=EL(1)+dx        !临时坐标系下的a2
            a2(2)=EL(2)+dy
            a2(3)=EL(3)
           CALL LTXYZ(B2,ET,a2)   !得到的点a2再转到整体坐标系
            write(12,*)'q',et(1:3)
            
             t1=sqrt((a0(1)-ET(1))**2+(a0(2)-ET(2))**2+
     #       (a0(3)-ET(3))**2)
               a0(1:3)=ET(1:3)     !ET是点a2的整体坐标
            write(12,*)'t1=',i1,i2,t1
              if(t1.gt.5.0e-5.and.n5.le.150)then
              if(ix.eq.i2)then
                Ix0=i1
               else
                 Ix0=i2
               endif
                Ix=Ix0
                 n5=n5+1 
               goto 100  
              else
             write(12,*)'n5=',n5     
              endif 
              
         endsubroutine
           
!################### 求取曲面上的Z坐标 ############        
       SUBROUTINE CurvPlanePoint(N,XY0,X,Xt, Yt, Zt)   
!##################################################
        implicit real*8 (a-h,o-z) 	
         dimension X(5000),XY0(5000,2),Rxy(500),Axy(500)
         
          do K = 1, N - 3
          Rxy(K + 3) = (Xt - XY0(K,1))**2 + (Yt - XY0(K,2))**2
          Axy(K + 3) = Rxy(K + 3) * Log(Rxy(K + 3) + 0.00001)
         enddo
         Axy(1) = 1
         Axy(2) = Xt
         Axy(3) = Yt
    
         Zt = 0.0
 
         do K = 1, N   
          Zt = Zt + Axy(K) * X(K)  
         enddo      
          End  
  
!###############################       
       Subroutine ErFenFa(N1,XY0,X,CC0,QJD,k12,k20)       !'''  求线段与曲面的交点。采用二分法
!###############################     
      implicit real*8 (a-h,o-z) 	 
      DIMENSION X(5000),XY0(5000,2),CC0(3, 3),QJD(5800, 3)  !debug发现800不够，改成5800. 2018.9.17
 
          K11 = 0
          K12 = K12 + 1    !记录曲面与线段交点的总数量
                           !由于块体不同面的同一线段或同一面的不同三角形的起点-终点顺序不同,二分法计算得到结果有所不同,交点有细微的不重合,
                           !通过交点坐标的比较,当充分接近时认为是同一点,以避免二分法计算的误差引起的不共点问题
100      Cx1 = (CC0(1, 1) + CC0(2, 1)) / 2.0
         Cy1 = (CC0(1, 2) + CC0(2, 2)) / 2.0
         Cz1 = (CC0(1, 3) + CC0(2, 3)) / 2.0
         K11 = K11 + 1
          Call CurvPlanePoint(N1,XY0,X,Cx1, Cy1, Zt1)
           Cx2 =CC0(1, 1); Cy2 =CC0(1, 2); Cz2 =CC0(1, 3)
          Call CurvPlanePoint(N1,XY0,X,Cx2, Cy2, Zt2)
!          write(12,*)'C1',Zt1,Cz1,Zt2,Cz2,
!     #        (Zt1 - Cz1)*(Zt2 -Cz2)
!         If ((Zt1 - Cz1) * (Zt2 -Cz2) .lt.1.5e-6) Then  ! 位于曲面的异侧。-0.00001改为-1.0e-6，长期存在的一处曲线跳跃不见了，原因是中分点几乎位于
                                                         ! 曲面上时，造成该循环判断取CC0(2,*)或 CC0(1,*)并最终引起二分法得到的交点偏离了曲面 2018.9.10     
            If ((Zt1 - Cz1).GT.0.0.AND.(Zt2 -Cz2) .lt.0.0.OR.
     #      (Zt1 - Cz1).LT.0.0.AND.(Zt2 -Cz2) .Gt.0.0 ) Then 
          CC0(2, 1) = Cx1
          CC0(2, 2) = Cy1
          CC0(2, 3) = Cz1
         Else
          CC0(1, 1) = Cx1
          CC0(1, 2) = Cy1
          CC0(1, 3) = Cz1
         End If
 
         If(Abs(CC0(1, 1) - CC0(2, 1)) + Abs(CC0(1, 2) - CC0(2, 2)) + 
     #   Abs(CC0(1,3)- CC0(2,3)).gt.0.00001.And.K11.lt.1000) GoTo 100
 
          CC0(3, 1) = (CC0(1, 1) + CC0(2, 1)) / 2.0
          CC0(3, 2) = (CC0(1, 2) + CC0(2, 2)) / 2.0
          CC0(3, 3) = (CC0(1, 3) + CC0(2, 3)) / 2.0
          
          if(k20.eq.0)return  !不进行比较和捏点        2018.9.8
           
!''''''''''' 通过交点坐标的比较,当充分接近时认为是同一点,以避免二分法计算的误差引起的不共点问题''''''''
       If (K12 .eq. 1) Then
        QJD(K12, 1) = CC0(3, 1)                    !''' 将第1个交点进行储存
        QJD(K12, 2) = CC0(3, 2)        ! QJD保存所有的交点坐标以进行距离判断，当距离很近时进行捏点
        QJD(K12, 3) = CC0(3, 3)
       GoTo 101
       End If
  
         do  K = 1 , K12 - 1
          If(Sqrt((QJD(K,1)-CC0(3,1))**2+(QJD(K,2)-CC0(3,2))**2
     #          +(QJD(K,3)-CC0(3,3))**2).lt.1.0e-3) Then  ! <易出错的小量控制>
                                                          !'''限制标准不能小，否则会出错。在此取0.2  0.5
                                        ! 发现取2.0e-4时，两条曲线没有相交，端点坐标距离有些明显，2.0e-4改为1.0e-2后相交。2018.9.20
                                        ! 看来曲线的坐标计算精度不太够. 0.5e-2改为1e-2，计算正确 2018.10.19 
            CC0(3, 1:3) = QJD(K, 1:3)   !'''取原来交点坐标,这样可以将交点合并.  <调试入口>     
            GoTo 101
          End If
        enddo
   
       QJD(K12, 1:3) = CC0(3, 1:3)                !''本次求得的交点与已求交点不重合,存储在QJD(K12, 1),以便下次比较时 
       
101    continue
       End  

    
!######################################        
      Subroutine GroundSurf()    ! 曲面（地面等）求解 
!###################################### 
        implicit real*8 (a-h,o-z) 	
        CHARACTER*60 fl2
       dimension CoCur(100,5000,3),XX(100,5000),NDcur(100),NEcur(100)
       dimension a(5000,5000),b(5000),X(5000),a1(5000,5000),R(5000,5000)
       dimension INDX(5000),CAZOOM(8,3),XY0(5000,2), CVEF(800,2)
       dimension coor(2,3)
        COMMON/Ncur0/Ncur/Ncur1/CoCur,XX,NDcur,NEcur/VRANGE/VRANGE 
        COMMON/CAOM/CAZOOM/DT0/DT0,stp0/fl2/fl2/IDelu/ANGTOL,IDelu,NMIN
         
       OPEN(21,FILE='CurSuf0.dxf')
        call DXFHEAD(21)       
        open(11,FILE=fl2) 
        open(12,file='CurBlock.dat') 
        read(11,*) Ncur
        write(1001,*) "Curved face Number:          ",Ncur
       do 20 I=1,Ncur
          READ(11,*)NDcur(i),NEcur(i)   !每个面有NDcur(i)个坐标点,前NEcur(i)个为面的范围点
       DO 20 j=1,NDcur(i)
        READ(11,*) CoCur(I,j,1:3)
20    CONTINUE
       READ(11,*)DT0,stp0          !多边形平面和曲面区域要进行三角形剖分后再求交线，DT0控制最大的剖分边长
                                   !根据模型大小,估计曲线前进的步长   数学上能否证明总是收敛 
                                   !不同stp造成有的线段和线段交点几乎位于端点上，造成后续回路错，改变stp后可以克服这种错误
                                   !因此通过人为输入，作为纠错的手段 2018.9.24. 改变DT0,stp0时也能使得delaunary容易通过. 2018.10.10
        READ(11,*)IDelu            !是否对块体的曲表面进行delaunary三角化,0否1是
        if(IDelu.eq.1)read(11,*)  ANGTOL,NMIN
        Write(1001,*)'stp0=',stp0
!         READ(11,*), XQ, XZ, XD
!         READ(11,*), YQ, YZ, YD
!         READ(11,*), ZQ
        do 200 i0=1,Ncur
        do I = 1,NDcur(i0)
        do J = 1,NDcur(i0)
         R(I,J) =(CoCur(I0,i,1) - CoCur(I0,J,1))** 2 +
     #      (CoCur(i0,I,2) - CoCur(i0,J,2)) ** 2
       enddo
       enddo
          ! ' 系数矩阵
    
        do I = 1,NDcur(i0)
       a(I, 1) = 1
       a(I, 2) = CoCur(I0,i,1)
       a(I, 3) = CoCur(i0,I,2)
       
      do J = 1,NDcur(i0)
      a(I, J + 3) = R(I, J) * Log(R(I, J) + 0.00001)
      enddo
      enddo 
    
       a(NDcur(i0) + 1, 1:3) =0.0 
       a(NDcur(i0) + 2, 1:3) =0.0 
       a(NDcur(i0) + 3, 1:3) =0.0 
       
      do J = 1,NDcur(i0)
       a(NDcur(i0) + 1, J + 3) = 1.0
       enddo
      do J = 1,NDcur(i0)
       a(NDcur(i0) + 2, J + 3) = CoCur(i0,j,1)
       enddo  
      do J = 1,NDcur(i0)
       a(NDcur(i0) + 3, J + 3) = CoCur(i0,j,2) 
       enddo  
       
      ! '输入已知的方程组的右端向量
      do I = 1, NDcur(i0)
       b(I) = CoCur(i0,i,3)  
       enddo
      b(NDcur(i0) + 1) = 0.0
      b(NDcur(i0) + 2) = 0.0
      b(NDcur(i0) + 3) = 0.0
   
       NDcur(i0) = NDcur(i0) + 3

        do I = 1, NDcur(i0)
        do J = 1 ,NDcur(i0)
            a1(I, J) = a(I, J)
        enddo
        enddo
       Call LUDCMP(a1, NDcur(i0), INDX, P)
       do L = 1 ,NDcur(i0)      
            X(L) = b(L)
        enddo
       
        Call LUBKSB(a1, NDcur(i0), INDX, X)
          do L = 1 ,NDcur(i0)      
            XX(I0,L) = X(L)
        enddo
          Write(12,*)'计算出的方程组的解:'
        do I = 1 , NDcur(i0)
           Write(12,*)'XX(', I,')=', XX(I0,i)
        enddo 
200        continue   
           
!!!       网格表示曲面，插值求得各点坐标并进行图示
          
           XD =VRANGE/ 50.0           ! '''根据模型大小,估计步长     
           yD =VRANGE/ 50.0 
           XQ=CAZOOM(1,1)   
           XZ=CAZOOM(2,1)
           YQ=CAZOOM(1,2)
           YZ=CAZOOM(4,2)
           Kx = Int((XZ - XQ) / XD)
           Ky = Int((YZ - YQ) / YD)
          
           ict=1
          do 300 i1=1,Ncur
           NVEF=NEcur(I1)  
           DO K0=1,NVEF
            CVEF(K0,1:2)= CoCur(i1,K0,1:2)
           enddo
            icc=1
             DO K0=1,NVEF
              K1=K0+1
              if(k1.gt.NVEF) k1=1
             coor(1,1:2)=CVEF(K0,1:2) 
             coor(2,1:2)=CVEF(K1,1:2) 
             coor(1:2,3)=0.0
!              Call DLINE(21,coor,ict,icc) 曲面的xy范围
             enddo
             
            icc=2
            N1=NDcur(i1)                           
            XY0(1:N1-3,1:2)=CoCur(i1,1:N1-3,1:2)  
            X(1:N1)=XX(i1,1:N1)      
          
           do 220 I = 1 , Kx/2 + 1
           do 220 J = 1 , Ky/2
             Xt1 = XQ + XD *2* (I - 1)       !/2和*2，表示面的平行线稀疏一半 2021.12.7
             Yt1 = YQ + YD *2* (J - 1)
             Xt2 = XQ + XD *2* (I - 1)
             Yt2 = YQ + YD *2* J
            CALL POINTC(Xt1,Yt1,NVEF,CVEF,N41)
            Call POINTC(Xt2,Yt2,NVEF,CVEF,N42)  
          
            If (N41.eq.0.and.N42.eq.0) Then                !2个点都在曲面区域内
            Call CurvPlanePoint(N1,XY0,X,Xt1, Yt1, Zt1) 
            Call CurvPlanePoint(N1,XY0,X,Xt2, Yt2, Zt2)
              coor(1,1)=Xt1;coor(1,2)=Yt1;coor(1,3)=Zt1
              coor(2,1)=Xt2;coor(2,2)=Yt2;coor(2,3)=Zt2
        if(i1.eq.3.or.i1.eq.6) Call DLINE(21,coor,ict,icc) 
            End If
     
220        continue 
300       continue      
         
         EndSubroutine
      
!###############################      
      Subroutine LUBKSB(a, N, INDX, b)
!############################### 
       implicit real*8 (a-h,o-z) 	
        dimension INDX(5000),a(5000,5000),b(5000)
        
        II = 0
       do  I = 1, N
        LL = INDX(I)
        Sum = b(LL)
        b(LL) = b(I)
         write(12,*)'b(LL)',b(LL)
        If (II.ne.0) Then       
          do  J = II, I - 1
            Sum = Sum - a(I, J) * b(J)
             write(12,*)'Sum',Sum,a(I, J)
         enddo
        Else
         If (Sum.ne.0.0) II = I
        End If

        b(I) = Sum
        enddo !Next I
        
         do  I = N, 1, -1
            Sum = b(I)
      
        If (I.lt.N) Then
           do  J = I + 1, N
             Sum = Sum - a(I, J) * b(J)
           enddo
        End If
        b(I) = Sum / a(I, I)
        enddo  !Next I
        End 
        
!############################### 
      Subroutine LUDCMP(a, N, INDX, D)
!############################### 
       
       implicit real*8 (a-h,o-z) 	 
        DIMENSION VV(5000),INDX(5000),a(5000,5000),b(5000)
          TINY = 1E-20
         D = 1.0
         do  I = 1, N
          AAMAX = 0.0
         do  J = 1,N
!               write(12,*)'ai0',a(I, J)    debug和release结果一样 2018.9.19
         If (Abs(a(I, J)).GT.AAMAX)  AAMAX = Abs(a(I, J))
        ENDDO
        If (AAMAX.EQ.0.0) WRITE(12,*)'Singular matrix.'
          VV(I) = 1.0 / AAMAX
         ENDDO
         
         do  J = 1, N
         If (J.gt.1) Then
         do  I = 1 , J - 1
           Sum = a(I, J)
         If (I.gt. 1) Then
           do  K = 1, I - 1
!           write(12,*)'ai0',a(I, K) , a(K, J) ,a(I, K) * a(K, J) debug和release结果不一样！ 2018.9.19
              Sum = Sum - a(I, K) * a(K, J)
           enddo
              a(I, J) = Sum
!                 write(12,*)'ai1',a(I, J) 
        End If
         enddo  !-- I
        End If
        AAMAX = 0.0
         do  I = J,N
            Sum = a(I, J)
            If (J.gt. 1) Then
                do  K = 1, J - 1
                  Sum = Sum - a(I, K) * a(K, J)
               enddo
                a(I, J) = Sum
!                   write(12,*)'ai2',a(I, J) 
            End If
            DUM = VV(I) * Abs(Sum)
!           write(12,*)'aij',DUM
            If (DUM.ge.AAMAX) Then
                IMAX = I
                AAMAX = DUM
            End If
        enddo
        If (J.ne.IMAX) Then
             do  K = 1 , N
                DUM1 = a(IMAX, K)
                a(IMAX, K) = a(J, K)
                a(J, K) = DUM1
             enddo  
            D = -D
            VV(IMAX) = VV(J) 
        End If
        INDX(J) = IMAX
        If (J.ne.N) Then
            If (a(J, J).eq.0.0)   a(J, J) = TINY
            DUM = 1.0 / a(J, J)
             do  I = J + 1, N
!                  write(12,*)'ai1',a(I, J),DUM
                a(I, J) = a(I, J) * DUM
!                  write(12,*)'aij',a(I, J),DUM
            enddo
      End If
          enddo  ! Next J
         If (a(N, N).eq.0.0) a(N, N) = TINY
       End  

