 Instructions for code running

All computer codes used in this article can be found online at https://github.com/zqh769282927692/CBCA-code. The codes can be used as follows. 
(1)	Download all files and place them together in a folder.
(2)	There are two vfproj files: B3dc_sep_Cur.vfproj and B3dcDraw_Cur.vfproj. The former is used for analysis and the latter for plotting. Open them using Visual Studio 2010.
(3)	First, run B3dc_sep_Cur.vfproj. This project includes B3DC_sep_Cur.for, Trace_cur.for, and Cal_Geometry.for, as well as DelaunaryDLL.lib, DelaunaryDLL.dll. These files will be loaded automatically when opening this project. Two file names of input data files, 6mian12.dat, and CurSufInp56.txt, are put in Inpt_B3DC.DAT. When running, the program will load Inpt_B3DC.DAT and then load the input data. The meanings of the input data in these files are labeled briefly in the files, as shown in Appendix 1, 2, and 3.
(4)	Then, run B3dcDraw_Cur.vfproj. This project includes B3DCdraw_Cur.for. Its corresponding input data file is Inpt_B3DC_Draw.DAT. 
(5)	Finally, a file named BLKCUT3D.DXF is generated. You can open it using AUTOCAD, as shown in Fig.1. It is exactly Fig. 8 in the article. If the input data are changed, different analysis results will be obtained. Besides, the analysis results are presented in B3dcDATA.dat and B3dcVol.dat.
           
Fig.1 The example of code runing.




Appendix 1:  Input data in CurSufInp56.txt
12                                 'There are 12 surfaces.
27	4                          'There are 27 points defining the first surface. Among them, 
 ' the first 4 points delineate the boundary of this surface.   
-10	-2	70                           'The coordinates of the first point
130	-2	30                       'The coordinates of the second point
130	130	65
-10	130	100
-10	19	75
-10	28	80
-10	36	85
-10	53.3	90
-10	70	95
20	7	65
20	25	75
20	60	85
20	80	91
40	3.4	55
40	18.3	60
40	46.2	70
40	76.8	80
40	100	87
60	9.5	40
60	15.9	50
60	50	60
60	82	70
90	8	33
90	18	40
90	45.4	50
90	80	60
90	98	65

27	4	                                 'Definition of the second surface
-10	-2	40
130	-2	2.0
130	130	35
-10	130	70
-10	19	45
-10	28	50
-10	36	55
-10	53.3	60
-10	70	65
20	7	35
20	25	45
20	60	55
20	80	61
40	3.4	25
40	18.3	30
40	46.2	40
40	76.8	50
40	100	57
60	9.5	10
60	15.9	20
60	50	30
60	82	40
90	8	10
90	18	10
90	45.4	20
90	80	30
90	98	35 
   
6   4
-10 -10  52.5 
130 -10  52
130 130  21
-10 130  30
40  -10 60      
55 45  55.5

6   4
-10  -10  32.5
130  -10  34          
130 130  30
-10 130  35
40  10  44
35 45  43.5	

 6   4
 -10  -10  16.5
130  -10  12
130 130  10
-10 130  15
40  10 32
35 45  25.5

8   4
 -10  -10  66.5
130  -10  62
130 130  40
-10 130  65
40  10    59
35  45   58
60  50  53
55 55  57.5

8   4
 -10  -10  37.5
130  -10  33
130 130  24
-10 130  19
40  10   38
35  45  33
60  50  30
55 55  28

5   4
 -10  -10  76.3
130  -10  71
130 130  41
-10 130  43 
60  50 56
 
5   4
 -10  -10  65.3
130  -10  32
130 130  29.5
-10 130 69.8
60  50   46
 
5   4
 -10  -10  35.3
130  -10  35.9
130 130  79.5
-10 130  89.8
   60  50  66.8

5   4
 -10  -10   -5.3
130  -10  15.9
130 130  89.5
-10 130  39.8  
   60  50  36.8      

 5   4
 -10  -10   63
130  -10  79
130 130  19.5
-10 130  9.8  
 60  50  38      

15 2.4             ! The maximum length of triangles of slave surface, the length of step 
0              !Triangulate the surfaces using Delaunary or not. 0 means no and 1 means yes.
4    8           ！Two control parameters of delaunary.
 
  
Appendix 2:  Input data in 6mian12.dat
1                                         ‘ Frequency of program running
90                                        ‘ Direction of the X axis.
-15.0 -15.0 -15.0                    ‘ Domain of the stochastic fracture network simulation.
135.0 -15.0 -15.0
135.0 135.0 -15.0
-15.0 135.0 -15.0 
-15.0 -15.0 95.0
135.0 -15.0 95.0
135.0 135.0 95.0
-15.0 135.0 95.0 
1                                    ‘ Slope analysis
0                                    ‘The number of excavation faces on the boundary   
     
1  0   0       ‘ The number denoting on or under the boundary surface 1, Dip, Dip direction
0  0   0   
1  90 270         
1  90 90
1  90  0  
1  90  180  
0                         ‘ There is no concave area
0                         ‘ The quantity of stochastic joint groups
0                         ‘ The quantity of parallel formations
0                          ‘The quantity of determinate circular fractures

2                          ‘ The quantity of determinate rectangular fractures
60 180  4                              'Definition of the first fractures
-4  80 95
124 80 95
124  10 -5 
-4    10 -5   

40  0  4                                'Definition of the second fractures
-4  5  83         
124  5  83
124  101 -28 
-4  101 -28 

0   0.2                         
0  
1     1
0.6
  
1          ! Is the surfaces considered？“1” denotes yes.
1          ！The supplementary surface is used as the boundary surface, “1” denotes yes.
1  1       ！The blocks which are located above or below the surface are retained (“1” denotes that the below are retained), the number of the supplementary surface.

 
Appendix 3:  Input data in Inpt_B3DC_Draw.DAT

y              ‘ Whether or not generate the figures in the dxf format. “y” means yes.
0             ‘ The number of the block which will be output. If “0” is input, all blocks 
will be output.

12              ‘The number of face which will be displayed. For debugging the code, 
             ‘ it is more convenient to draw the intersection lines in a specified surface 
 ‘ and then find the problem that needs to be solved. 
  
  
