      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
      
#include <SMAAspUserSubroutines.hdr>


C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
      
      
      REAL :: k_WW, k_rr, k_M0, k_Dist_min, k_C0, k_CC, k_r0, k_rend, 
     1 k_rend_1, k_Volume, k_p, k_Max, k_Min, k_erosion(20000) 
      INTEGER :: K_Number_Elements, k_Number_Neighb, 
     1 K_Number_Elements_0, k_Npt 
      common k_WW, k_rr, k_M0, k_Dist_min, K_Number_Elements, k_p, 
     1 k_C0, k_CC, k_r0, k_Number_Neighb, k_rend, k_rend_1, k_Volume,
     2 k_Max, k_Min, k_Npt, K_Number_Elements_0, k_erosion
      
      
      real k_rd(10)
      pointer(ptrk_rd,k_rd) 
      
      real k_r(10)
      pointer(ptrk_r,k_r) 
         
C------ Surface loading         
       IF (JLTYP.EQ.0) THEN
           
           IF (KINC.LT.50) THEN
               
                 F=110000.0
           
           ELSE
               
               F=110000.0*(1+0.005*(50.0-KINC*1.0))
               
           END IF
           
           
           IF (KINC.GT.250) THEN
               
                 F=0
           
           END IF           
           
      
       ELSE
       
C------ Gravity force              
      F=-22000
      
      IF (KINC.GE.3) THEN
          
      ptrk_rd = SMAFloatArrayCreateSP(21,K_Number_Elements,0.0)
      
      ptrk_r = SMAFloatArrayCreateSP(16,K_Number_Elements,0.0)  
       
       
          IF (k_rd((NOEL-1)*k_Npt+Npt).GT.0.7) THEN
          
              F=-22000
          
C          F=0
          
          ELSE
              
              F=0
              
          END IF
          
            
      END IF
      
      END IF
  
      RETURN
      END
      
      
      
C      -------------------------------------------------------------  
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      
#include <SMAAspUserSubroutines.hdr>
      
C
      DIMENSION TIME(2)
      
      REAL :: k_WW, k_rr, k_M0, k_Dist_min, k_C0, k_CC, k_r0, k_rend, 
     1 k_rend_1, k_Volume, k_p, k_Max, k_Min, k_erosion(20000) 
      INTEGER :: K_Number_Elements, k_Number_Neighb, 
     1 K_Number_Elements_0, k_Npt 
      common k_WW, k_rr, k_M0, k_Dist_min, K_Number_Elements, k_p, 
     1 k_C0, k_CC, k_r0, k_Number_Neighb, k_rend, k_rend_1, k_Volume,
     2 k_Max, k_Min, k_Npt, K_Number_Elements_0, k_erosion
      
       character*256 OUTDIR
       character*256 OUTDIR2
      
       real k_S1(10)
      pointer(ptrS1,k_S1)  
      real k_E1(10)
      pointer(ptrE1,k_E1) 
      
      IF (LOP.EQ.0) THEN

C---- Input    

        CALL GETOUTDIR( OUTDIR, LENOUTDIR )        
        OUTDIR = trim(OUTDIR) // trim('\Input_Data.txt')
        
        open(unit=105, file=OUTDIR)
                
        read (105,*)
        read (105,*) k_p

        read (105,*)
        read (105,*) k_r0    

        read (105,*)
        read (105,*) k_rend

        read (105,*)
        read (105,*) k_C0
        
        read (105,*)
        read (105,*) k_CC
        
C-----  Distance between neighboring points R         
        read (105,*)
        read (105,*) k_dist_min
        
C----- Maximum density          
        read (105,*)
        read (105,*) k_Max
        
C----- Minimum density         
        read (105,*)
        read (105,*) k_Min           
        
        close(105)  
        
          
C---- Input
 
C--        k_p=2  
        
C--         k_r0=0.5
        
C--         k_rend=0.3
        
C--         k_C0=0.001
        
C--         k_CC=1.3
        
C--         k_dist_min=2
        
C--         k_Max=1
        
C--         k_Min=0.0001
        
C-------------------------------        
        k_rend_1=k_rend
      
        K_Number_Elements=1
     
        k_Number_Neighb=1
        
        k_Npt=1
        
      END IF
   
      RETURN
      END       

C      -------------------------------------------------------------  
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
      
#include <SMAAspUserSubroutines.hdr>
      
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      
      REAL :: k_WW, k_rr, k_M0, k_Dist_min, k_C0, k_CC, k_r0, k_rend, 
     1 k_rend_1, k_Volume, k_p, k_Max, k_Min, k_erosion(20000) 
      INTEGER :: K_Number_Elements, k_Number_Neighb, 
     1 K_Number_Elements_0, k_Npt 
      common k_WW, k_rr, k_M0, k_Dist_min, K_Number_Elements, k_p, 
     1 k_C0, k_CC, k_r0, k_Number_Neighb, k_rend, k_rend_1, k_Volume,
     2 k_Max, k_Min, k_Npt, K_Number_Elements_0, k_erosion
     
      real k_w(10)
      pointer(ptrk_w,k_w)
      
      real k_r(10)
      pointer(ptrk_r,k_r) 
      
      real k_Coord_1(10)
      pointer(ptrk_Coord_1,k_Coord_1) 
      
      real k_Coord_2(10)
      pointer(ptrk_Coord_2,k_Coord_2)
      
      real k_Coord_3(10)
      pointer(ptrk_Coord_3,k_Coord_3)  
      
      real k_Neighb(10)
      pointer(ptrk_Neighb,k_Neighb)
      
      real k_rd(10)
      pointer(ptrk_rd,k_rd) 
      
C----- k_C
      real k_C(10)
      pointer(ptrk_C,k_C)
      
C----- E1-E6      
      real k_E1(10)
      pointer(ptrE1,k_E1)   
      real k_E2(10)
      pointer(ptrE2,k_E2) 
      real k_E3(10)
      pointer(ptrE3,k_E3) 
      real k_E4(10)
      pointer(ptrE4,k_E4)             
      real k_E5(10)
      pointer(ptrE5,k_E5)  
      real k_E6(10)
      pointer(ptrE6,k_E6)        
      
      
C     KSTEP=1 
      
      IF (KSTEP.EQ.1) THEN 
      
      FIELD(1)=1
      FIELD(2)=1

      
      END IF
      
C     KSTEP=2 
      
      IF (KSTEP.EQ.2) THEN 
      
C     KINC=1
  
      IF (KINC.EQ.1) THEN
        
        IF (NOEL.GE.K_Number_Elements_0) THEN
           K_Number_Elements_0=NOEL
        END IF
        
        IF (Npt.GE.k_Npt) THEN
           k_Npt=Npt
        END IF
        
        K_Number_Elements=K_Number_Elements_0*k_Npt
        
        FIELD(1)=1
       FIELD(2)=1


      END IF 
      
C     KINC=2
      
      IF (KINC.EQ.2) THEN
    
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,K_Number_Elements,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,K_Number_Elements,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,K_Number_Elements,0.0)
      
      k_Coord_1((NOEL-1)*k_Npt+Npt)=COORD(1)
      k_Coord_2((NOEL-1)*k_Npt+Npt)=COORD(2)
      k_Coord_3((NOEL-1)*k_Npt+Npt)=COORD(3)
         

       FIELD(1)=1

      FIELD(2)=1
      
      END IF 
      
 
C     KINC>2     
      IF (KINC.GT.2) THEN    
    
      ptrk_w = SMAFloatArrayCreateSP(13,K_Number_Elements,0.0)
      
      ptrk_r = SMAFloatArrayCreateSP(16,K_Number_Elements,0.0)
      
      ptrk_Neighb = SMAFloatArrayCreateSP(20,
     1 K_Number_Elements*(k_Number_Neighb+1)*2,0.0)
     
      ptrk_rd = SMAFloatArrayCreateSP(21,K_Number_Elements,0.0)
      
      ptrk_C = SMAFloatArrayCreateSP(14,K_Number_Elements,0.0)      
      
      ptrE1 = SMAFloatArrayCreateSP(7,K_Number_Elements,0.0)
      ptrE2 = SMAFloatArrayCreateSP(8,K_Number_Elements,0.0)
      ptrE3 = SMAFloatArrayCreateSP(9,K_Number_Elements,0.0)
 
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,K_Number_Elements,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,K_Number_Elements,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,K_Number_Elements,0.0)
      
      
       STATEV(1)=k_r((NOEL-1)*k_Npt+Npt)
       STATEV(2)=k_rd((NOEL-1)*k_Npt+Npt)
      STATEV(3)=k_w((NOEL-1)*k_Npt+Npt)
      STATEV(4)=k_C((NOEL-1)*k_Npt+Npt)
 
C      STATEV(5)=k_E1((NOEL-1)*k_Npt+Npt)
      STATEV(6)=k_E2((NOEL-1)*k_Npt+Npt)      
      STATEV(7)=k_E3((NOEL-1)*k_Npt+Npt)

C----- Filtration of density  
    
      IF (STATEV(2).GT.0.677) THEN 
      
           FIELD(1)=1
           STATEV(5)=1
      
      ELSE
          
          FIELD(1)=k_Min
          STATEV(5)=0

      END IF  
      
C------ Location of of the discontinuity       
       IF ((k_Coord_2((NOEL-1)*k_Npt+Npt).LT.0.01).AND.
     1  (k_Coord_2((NOEL-1)*k_Npt+Npt).GT.0.0)) THEN 
      
           FIELD(2)=2

      ELSE
          
           FIELD(2)=1

      END IF      

      END IF
      
      END IF
      
      RETURN
      END
      
C-------------------------------------------------------------
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
C
      INCLUDE 'ABA_PARAM.INC'

#include <SMAAspUserSubroutines.hdr>
 
C 
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
  
      REAL :: k_WW, k_rr, k_M0, k_Dist_min, k_C0, k_CC, k_r0, k_rend, 
     1 k_rend_1, k_Volume, k_p, k_Max, k_Min, k_erosion(20000) 
      INTEGER :: K_Number_Elements, k_Number_Neighb, 
     1 K_Number_Elements_0, k_Npt 
      common k_WW, k_rr, k_M0, k_Dist_min, K_Number_Elements, k_p, 
     1 k_C0, k_CC, k_r0, k_Number_Neighb, k_rend, k_rend_1, k_Volume,
     2 k_Max, k_Min, k_Npt, K_Number_Elements_0, k_erosion
 
      real :: k_ener, k_C1, k_Dist, k_Number_Neighb_0, k_rdd, k_force,
     1 k_smax, k_smin
      INTEGER :: k_i, k_j, k_N, k_F

C----- S1-S6
      real k_S1(10)
      pointer(ptrS1,k_S1)   
      real k_S2(10)
      pointer(ptrS2,k_S2) 
      real k_S3(10)
      pointer(ptrS3,k_S3) 
      real k_S4(10)
      pointer(ptrS4,k_S4)             
      real k_S5(10)
      pointer(ptrS5,k_S5)  
      real k_S6(10)
      pointer(ptrS6,k_S6)  
      
C----- E1-E6      
      real k_E1(10)
      pointer(ptrE1,k_E1)   
      real k_E2(10)
      pointer(ptrE2,k_E2) 
      real k_E3(10)
      pointer(ptrE3,k_E3) 
      real k_E4(10)
      pointer(ptrE4,k_E4)             
      real k_E5(10)
      pointer(ptrE5,k_E5)  
      real k_E6(10)
      pointer(ptrE6,k_E6)  
      
C----- k_w      
      real k_w(10)
      pointer(ptrk_w,k_w)
      
C----- k_C
      real k_C(10)
      pointer(ptrk_C,k_C)
   
C----- k_V
      real k_V(10)
      pointer(ptrk_V,k_V)      
      
C----- k_r
      real k_r(10)
      pointer(ptrk_r,k_r) 

C----- k_Coord      
      real k_Coord_1(10)
      pointer(ptrk_Coord_1,k_Coord_1) 
      
      real k_Coord_2(10)
      pointer(ptrk_Coord_2,k_Coord_2)
      
      real k_Coord_3(10)
      pointer(ptrk_Coord_3,k_Coord_3)
      
C-----   k_Neighb          
      real k_Neighb(10)
      pointer(ptrk_Neighb,k_Neighb)
      
      real k_rd(10)
      pointer(ptrk_rd,k_rd) 
      
       LOVRWRT=1
       
       
      ptrS1 = SMAFloatArrayCreateSP(1,100,0.0)
      ptrS2 = SMAFloatArrayCreateSP(2,100,0.0)
      ptrS3 = SMAFloatArrayCreateSP(3,100,0.0)
      ptrS4 = SMAFloatArrayCreateSP(4,100,0.0)
      ptrS5 = SMAFloatArrayCreateSP(5,100,0.0)
      ptrS6 = SMAFloatArrayCreateSP(6,100,0.0)
   
      ptrE1 = SMAFloatArrayCreateSP(7,100,0.0)
      ptrE2 = SMAFloatArrayCreateSP(8,100,0.0)
      ptrE3 = SMAFloatArrayCreateSP(9,100,0.0)
      ptrE4 = SMAFloatArrayCreateSP(10,100,0.0)
      ptrE5 = SMAFloatArrayCreateSP(11,100,0.0)
      ptrE6 = SMAFloatArrayCreateSP(12,100,0.0)
      
      ptrk_w = SMAFloatArrayCreateSP(13,100,0.0)
      
      ptrk_C = SMAFloatArrayCreateSP(14,100,0.0)
      
      ptrk_V = SMAFloatArrayCreateSP(15,100,0.0)
      
      ptrk_r = SMAFloatArrayCreateSP(16,100,0.0)
      
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,100,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,100,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,100,0.0)
      
      ptrk_Neighb = SMAFloatArrayCreateSP(20,100,0.0)
      
      ptrk_rd = SMAFloatArrayCreateSP(21,100,0.0)
      
      
      IF (KSTEP.EQ.2) THEN
      
C     Finding stress and strain  
      IF (KINC.GT.1) THEN
      
      ptrS1 = SMAFloatArrayCreateSP(1,K_Number_Elements,0.0)
      ptrS2 = SMAFloatArrayCreateSP(2,K_Number_Elements,0.0)
      ptrS3 = SMAFloatArrayCreateSP(3,K_Number_Elements,0.0)
      ptrS4 = SMAFloatArrayCreateSP(4,K_Number_Elements,0.0)
      ptrS5 = SMAFloatArrayCreateSP(5,K_Number_Elements,0.0)
      ptrS6 = SMAFloatArrayCreateSP(6,K_Number_Elements,0.0)
      
      ptrE1 = SMAFloatArrayCreateSP(7,K_Number_Elements,0.0)
      ptrE2 = SMAFloatArrayCreateSP(8,K_Number_Elements,0.0)
      ptrE3 = SMAFloatArrayCreateSP(9,K_Number_Elements,0.0)
      ptrE4 = SMAFloatArrayCreateSP(10,K_Number_Elements,0.0)
      ptrE5 = SMAFloatArrayCreateSP(11,K_Number_Elements,0.0)
      ptrE6 = SMAFloatArrayCreateSP(12,K_Number_Elements,0.0)
      
      ptrk_w = SMAFloatArrayCreateSP(13,K_Number_Elements,0.0)
      
      ptrk_C = SMAFloatArrayCreateSP(14,K_Number_Elements,0.0)
      
      ptrk_V = SMAFloatArrayCreateSP(15,K_Number_Elements,0.0)
      
      ptrk_r = SMAFloatArrayCreateSP(16,K_Number_Elements,0.0)
      
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,K_Number_Elements,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,K_Number_Elements,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,K_Number_Elements,0.0)
      
      END IF

C----- Reading of the results

      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)
C
      N=0
      KK=0
      KKK=0
      

      DO K1=1,999999
         CALL DBFILE(0,JRRAY,JRCD)
        IF (JRCD .NE. 0) GO TO 110
        KEY=JRRAY(1,2)
      
C----- k_V      
      IF (KINC.LE.1) GO TO 110 
      
        IF (KEY.EQ.76) THEN
            KKK=KKK+1
                   
            IF (KKK.LE.K_Number_Elements+1) THEN
                
                    k_V(KKK)=array(3)
C                      k_V(KKK)=0.125
     
            END IF
               
        END IF         
          
C----- k_S      
        IF (KEY.EQ.11) THEN
            N=N+1
C       3 - s11
C       4 - s22 
C       5 - s33
C       6 - s12
C       7 - s13 
C       8 - s23                     
            IF (N.LE.K_Number_Elements+1) THEN
                k_S1(N)=array(3)
                k_S2(N)=array(4)
                k_S3(N)=array(5)
                k_S4(N)=array(6)
                k_S5(N)=array(7)
                k_S6(N)=array(8)
                
   
            END IF
               
        END IF        

C----- k_E        
        IF (KEY.EQ.401) THEN
            KK=KK+1
C       3 - e11
C       4 - e22 
C       5 - e33
C       6 - e12
C       7 - e13 
C       8 - e23                     
            IF (KK.LE.K_Number_Elements+1) THEN
                k_E1(KK)=array(3)
                k_E2(KK)=array(4)
                k_E3(KK)=array(5)

    
            END IF
               
        END IF        
          
      IF (KKK.EQ.K_Number_Elements+2) GO TO 110          
      IF (N.EQ.K_Number_Elements+2) GO TO 110 
      IF (KK.EQ.K_Number_Elements+1) GO TO 110
       
      END DO
      
  110 CONTINUE

C-----KINC=1

      IF (KINC.EQ.1) THEN
      
      ptrk_w = SMAFloatArrayCreateSP(13,K_Number_Elements,0.0)         
      ptrk_C = SMAFloatArrayCreateSP(14,K_Number_Elements,0.0)   
      ptrk_r = SMAFloatArrayCreateSP(16,K_Number_Elements,0.0) 
      ptrk_rd = SMAFloatArrayCreateSP(21,K_Number_Elements,0.0) 
     
       DO j=1, K_Number_Elements
        
           k_C(j)=0
           k_r(j)=k_r0
           k_rd(j)=k_r0
           k_w(j)=0 
           
       END DO

C------ Down       
C       DO j=1, 20000
      
C           k_r(j)=k_erosion(j)
           
C       END DO  
       
       
       
      END IF 
      
C----- KINC=2
C----- Search for neighboring points

      IF (KINC.EQ.2) THEN
      
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,K_Number_Elements,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,K_Number_Elements,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,K_Number_Elements,0.0)
      
C-----  k_Number_Neighb  
       
      k_Number_Neighb=1
       
       DO i=1, K_Number_Elements
        
        k_Number_Neighb_0=1
        k_Dist=0.0
        
          DO j=1, K_Number_Elements
          
            IF (i.NE.j) THEN
       k_Dist=((k_Coord_1(i)-k_Coord_1(j))**2+
     1 (k_Coord_2(i)-k_Coord_2(j))**2+
     2 (k_Coord_3(i)-k_Coord_3(j))**2)**0.5
         
            IF (k_dist_min.GE.k_Dist) then
                k_Number_Neighb_0=k_Number_Neighb_0+1
            END IF
            
           END IF   
         
         END DO
         
         IF (k_Number_Neighb_0.GT.k_Number_Neighb) then
                k_Number_Neighb=k_Number_Neighb_0
         END IF
          
       END DO
       
      ptrk_Neighb = SMAFloatArrayCreateSP(20,
     1 K_Number_Elements*(k_Number_Neighb+1)*2,0.0)

       
C-----  k_Neighb       
       
      DO i=1, K_Number_Elements
      
        k_Neighb((k_Number_Neighb+1)*2*(i-1)+1)=i
        k_Neighb((k_Number_Neighb+1)*2*(i-1)+3)=i
        k_Neighb((k_Number_Neighb+1)*2*(i-1)+4)=k_dist_min*k_V(i)
        
        k_Number_Neighb_0=1
        k_Dist=0.0
        
        DO j=1, K_Number_Elements
        
            IF (i.NE.j) THEN
       k_Dist=((k_Coord_1(i)-k_Coord_1(j))**2+
     1 (k_Coord_2(i)-k_Coord_2(j))**2+
     2 (k_Coord_3(i)-k_Coord_3(j))**2)**0.5
         
            IF (k_dist_min.GE.k_Dist) then
               k_Number_Neighb_0=k_Number_Neighb_0+1
               k_Neighb((k_Number_Neighb+1)*2*(i-1)+
     1          2*k_Number_Neighb_0+1)=j
               k_Neighb((k_Number_Neighb+1)*2*(i-1)+
     1         2*k_Number_Neighb_0+2)=(k_dist_min-k_Dist)*k_V(j)
         
            END IF
            
            END IF
        
      
        END DO
        
        k_Neighb((k_Number_Neighb+1)*2*(i-1)+2)=k_Number_Neighb_0
   
      END DO 
 
 
      DO i=1, K_Number_Elements
        k_Dist=0
        DO j=1,k_Number_Neighb
            k_Dist=k_Dist+
     1       k_Neighb((k_Number_Neighb+1)*2*(i-1)+2*j+2)
        END DO 
        
        DO j=1,k_Number_Neighb
            k_Neighb((k_Number_Neighb+1)*2*(i-1)+2*j+2)=
     1       k_Neighb((k_Number_Neighb+1)*2*(i-1)+2*j+2)/k_Dist
        END DO 
      
      END DO  
               
      k_Volume=0
      
      DO i=1, K_Number_Elements
        k_Volume=k_Volume+k_V(i)
      END DO 
      
C------ Initial density ro
 
      k_i=1
      
       DO j=1, K_Number_Elements      

          k_r(j)=k_Max
          
          
          IF ((k_Coord_2(j).LT.0.01).AND.(k_Coord_2(j).GT.0.0)) THEN 
      
      k_Dist=((k_Coord_1(j)-0.1)**2+
     2 (k_Coord_3(j))**2)**0.5          
              
              IF (k_Dist.LE.0.03) THEN 
      
                 k_r(j)=k_Min 
                  
              END IF    
              
      k_Dist=((k_Coord_1(j))**2+
     2 (k_Coord_3(j))**2)**0.5          
              
              IF (k_Dist.LE.0.03) THEN 
      
                 k_r(j)=k_Min 
                  
              END IF    
              
              
      k_Dist=((k_Coord_1(j)-0.2)**2+
     2 (k_Coord_3(j))**2)**0.5          
              
              IF (k_Dist.LE.0.03) THEN 
      
                 k_r(j)=k_Min 
                  
              END IF                
          
         END IF   
   
C          IF ((k_Coord_3(j).LT.0.01).AND.(k_Coord_2(j).LT.0.01)
C     1    .AND.(k_Coord_2(j).GT.0.0)) THEN 
      IF (k_Coord_3(j).LT.0.005)  THEN 
             
              k_r(j)=k_Min
               
          END IF        
     
     
      END DO
         
      END IF       
      
      
C-----KINC>2

      IF (KINC.GT.2) THEN
   
      ptrS1 = SMAFloatArrayCreateSP(1,K_Number_Elements,0.0)
      ptrS2 = SMAFloatArrayCreateSP(2,K_Number_Elements,0.0)
      ptrS3 = SMAFloatArrayCreateSP(3,K_Number_Elements,0.0)
      ptrS4 = SMAFloatArrayCreateSP(4,K_Number_Elements,0.0)
      ptrS5 = SMAFloatArrayCreateSP(5,K_Number_Elements,0.0)
      ptrS6 = SMAFloatArrayCreateSP(6,K_Number_Elements,0.0)
      
      ptrE1 = SMAFloatArrayCreateSP(7,K_Number_Elements,0.0)
      ptrE2 = SMAFloatArrayCreateSP(8,K_Number_Elements,0.0)
      ptrE3 = SMAFloatArrayCreateSP(9,K_Number_Elements,0.0)

      
      ptrk_w = SMAFloatArrayCreateSP(13,K_Number_Elements,0.0)
      
      ptrk_C = SMAFloatArrayCreateSP(14,K_Number_Elements,0.0)
      
      ptrk_V = SMAFloatArrayCreateSP(15,K_Number_Elements,0.0)
      
      ptrk_r = SMAFloatArrayCreateSP(16,K_Number_Elements,0.0)
      
      ptrk_rd = SMAFloatArrayCreateSP(21,K_Number_Elements,0.0)
      
      ptrk_Coord_1 = SMAFloatArrayCreateSP(17,K_Number_Elements,0.0)
      
      ptrk_Coord_2 = SMAFloatArrayCreateSP(18,K_Number_Elements,0.0)
      
      ptrk_Coord_3 = SMAFloatArrayCreateSP(19,K_Number_Elements,0.0)
      
 
C     Finding M0 - k_M0
  
      DO j=1, K_Number_Elements
          
          k_C(j)=0.0
          
      END DO
      
     
C     k_r
     
      IF (KINC.GT.4) THEN   
         
C------ Critical value      
      DO j=1, K_Number_Elements
          
         IF ((k_Coord_2(j).LT.0.01).AND.(k_Coord_2(j).GT.0.0)) THEN
         
C         k_force=-20000+abs(k_Coord_2(j))/0.03*19700
         
C         k_force=-20000
        k_force=-400000          
         else
             
           IF (k_Coord_2(j).GT.0.01) THEN   
              
             k_force=-20000
             
          ELSE
              
              k_force=-20000
              
          END IF            

             
         END IF
         
           IF (KINC.GT.50) THEN
               
                 k_force=-20000
           
           END IF
         
C----- Finding erosion front         
         k_Dist=0.0
         
         DO i=1, k_Neighb((k_Number_Neighb+1)*2*(j-1)+2)
            
      IF (k_r(k_Neighb((k_Number_Neighb+1)*2*(j-1)+3+2*(i-1))).EQ.k_Min)
     1   THEN  
             
            k_Dist=k_Dist+1.0
      END IF
       
         END DO
 
C------ Finding increment for density    
         
C         k_C1=(k_S1(j)+k_S2(j)+k_S3(j))/3.0
         
C         k_C1=k_E3(j)
        
         k_C1=k_E1(j)
         
C          k_smax=0
C          k_smin=0 
          
C          IF (k_E1(j).LT.0.0) THEN
              
C              k_smax=k_E1(j)
              
C          end if
         
C          IF (k_E2(j).LT.0.0) THEN
              
C              k_smin=k_E2(j)
              
C          end if
          
C          k_C1=k_smax+k_smin
         
         k_w(j)=k_C1
         
         IF (k_Dist.GT.0.0) THEN

             IF (k_C1.GT.k_force) THEN

C             IF (k_C1.GT.0) THEN
              
               k_C(j)=-1.0
 
C              ELSE
                 
C                  k_C(j)=-0.5*(1-k_C1/k_force)
                  
C              END IF
             
              
             END IF
             
               
          
          END IF 
                   
      ENDDO
     
C----- Finding density ro(n+1)       
      
      DO j=1, K_Number_Elements
         
         IF (k_r(j).GE.k_Min) THEN  
              k_r(j)=k_r(j)+k_C(j)  
         END IF
         
          IF (k_r(j).LE.k_Min) THEN
                k_r(j)=k_Min
 
           END IF
   
           IF (k_r(j).GE.k_Max) THEN
                k_r(j)=k_Max

           END IF  
           
            IF ((abs(k_Coord_2(j)).GT.0.09).AND.
     1       (k_Coord_3(j).GT.0.005)) THEN
           
               k_r(j)=k_Max
               
            END IF
            
                
      ENDDO   
             
      END IF
      
C----- filtration of density
      
      DO i=1, K_Number_Elements

            k_rd(i)=0
        
      ENDDO
      
      DO i=1, K_Number_Elements

         DO j=1,k_Number_Neighb

        k_rd(i)=k_rd(i)+k_r(k_Neighb((k_Number_Neighb+1)*2*(i-1)+2*j+1))
     1      *k_Neighb((k_Number_Neighb+1)*2*(i-1)+2*j+2)


        END DO
         
      ENDDO
      
      
      END IF
      
      END IF
      
      RETURN
      END
      
