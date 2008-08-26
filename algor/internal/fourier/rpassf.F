C     SUBROUTINE 'RPASSF' - PERFORMS ONE PASS THROUGH DATA AS PART      
C     OF MULTIPLE REAL FFT (FOURIER SYNTHESIS) ROUTINE                  
C                                                                       
C     A IS FIRST REAL INPUT VECTOR                                      
C         EQUIVALENCE B(1) WITH A (LA*INC1+1)                           
C     C IS FIRST REAL OUTPUT VECTOR                                     
C         EQUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)                       
C     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES                  
C     INC1 IS THE ADDRESSING INCREMENT FOR A                            
C     INC2 IS THE ADDRESSING INCREMENT FOR C                            
C     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A                     
C     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C                    
C     LOT IS THE NUMBER OF VECTORS                                      
C     N IS THE LENGTH OF THE VECTORS                                    
C     IFAC IS THE CURRENT FACTOR OF N                                   
C     LA IS THE PRODUCT OF PREVIOUS FACTORS                             
C     IERR IS AN ERROR INDICATOR:                                       
C              0 - PASS COMPLETED WITHOUT ERROR                         
C              1 - LOT GREATER THAN 64                                  
C              2 - IFAC NOT CATERED FOR                                 
C              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC                   
C     LIPL=.T. => RESULTS ARE RETURNED TO INPUT ARRAY
C              (ONLY VALID IF LA=N/IFAC, I.E. ON LAST PASS)
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE RPASSF(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,   
     *    LA,IERR,LIPL)                                                      
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)                            
      LOGICAL LIPL
C                                                                       
      DATA SIN36/0.587785252292473/,SIN72/0.951056516295154/,           
     *    QRT5/0.559016994374947/,SIN60/0.866025403784437/              
C                                                                       
      M=N/IFAC                                                          
      IINK=LA*INC1                                                      
      JINK=LA*INC2                                                      
      JUMP=(IFAC-1)*JINK                                                
      KSTOP=(N-IFAC)/(2*IFAC)                                           
C                                                                       
      IBASE=0                                                           
      JBASE=0                                                           
      IBAD=0                                                            
C
!     Increase the vector length by fusing the loops if the
!     data layout is appropriate:
      IF (INC1.EQ.LOT.AND.INC2.EQ.LOT.AND.INC3.EQ.1.AND.INC4.EQ.1) THEN
        ILA=1
        ILOT=LA*LOT
        INC21=LA*LOT
      ELSE
        ILA=LA
        ILOT=LOT
        INC21=INC2
      ENDIF
C
      IF (IFAC.EQ.2) THEN
C                                                                       
C     CODING FOR FACTOR 2                                               
C     -------------------                                               
  200 CONTINUE                                                          
      IA=1                                                              
      IB=IA+(2*M-LA)*INC1                                               
      JA=1                                                              
      JB=JA+JINK                                                        
C                                                                       
      IF (LA.NE.M) THEN                                            
C                                                                       
      DO 220 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 210 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      C(JB+J)=A(IA+I)-A(IB+I)                                           
      I=I+INC3                                                          
      J=J+INC4                                                          
  210 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  220 CONTINUE                                                          
      IA=IA+IINK                                                        
      IINK=2*IINK                                                       
      IB=IB-IINK                                                        
      IBASE=0                                                           
      JBASE=JBASE+JUMP                                                  
      JUMP=2*JUMP+JINK                                                  
C
      IF (IA.LT.IB) THEN                                           
      DO 250 K=LA,KSTOP,LA                                              
      KB=K+K                                                            
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      IBASE=0                                                           
      DO 240 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 230 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      D(JA+J)=B(IA+I)-B(IB+I)                                           
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))                 
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))                 
      I=I+INC3                                                          
      J=J+INC4                                                          
  230 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  240 CONTINUE                                                          
      IA=IA+IINK                                                        
      IB=IB-IINK                                                        
      JBASE=JBASE+JUMP                                                  
  250 CONTINUE                                                          
      ENDIF
C
      IF (IA.EQ.IB) THEN                                           
      IBASE=0                                                           
      DO 280 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 270 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)                                                   
      C(JB+J)=-B(IA+I)                                                  
      I=I+INC3                                                          
      J=J+INC4                                                          
  270 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  280 CONTINUE                                                          
      ENDIF                                                         
C                                                                       
      ELSE                !!! Case LA=M                                 
      IF (LIPL) THEN
        DO 294 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 292 IJK=1,ILOT                                                  
        T1=2.0*(A(IA+I)-A(IB+I))                                     
        A(IA+I)=2.0*(A(IA+I)+A(IB+I))                                     
        A(IB+I)=T1
        I=I+INC3                                                          
  292   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  294   CONTINUE                                                          
      ELSE
        DO 298 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 296 IJK=1,ILOT                                                  
        C(JA+J)=2.0*(A(IA+I)+A(IB+I))                                     
        C(JB+J)=2.0*(A(IA+I)-A(IB+I))                                     
        I=I+INC3                                                          
        J=J+INC4                                                          
  296   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  298   CONTINUE                                                          
      ENDIF
      ENDIF
C
      ELSEIF (IFAC.EQ.3) THEN
C                                                                       
C     CODING FOR FACTOR 3                                               
C     -------------------                                               
  300 CONTINUE                                                          
      IA=1                                                              
      IB=IA+(2*M-LA)*INC1                                               
      IC=IB                                                             
      JA=1                                                              
      JB=JA+JINK                                                        
      JC=JB+JINK                                                        
C                                                                       
      IF (LA.NE.M) THEN                                            
C                                                                       
      DO 320 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 310 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      C(JB+J)=(A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))                   
      C(JC+J)=(A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))                   
      I=I+INC3                                                          
      J=J+INC4                                                          
  310 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  320 CONTINUE                                                          
      IA=IA+IINK                                                        
      IINK=2*IINK                                                       
      IB=IB+IINK                                                        
      IC=IC-IINK                                                        
      JBASE=JBASE+JUMP                                                  
      JUMP=2*JUMP+JINK                                                  
C
      IF (IA.LT.IC) THEN                                           
      DO 350 K=LA,KSTOP,LA                                              
      KB=K+K                                                            
      KC=KB+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      IBASE=0                                                           
      DO 340 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 330 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))                                 
      D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))                                 
      C(JB+J)=                                                          
     *    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+I))))
     *   -S1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=                                                          
     *    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+I))))
     *   +C1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=                                                          
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+I))))
     *   -S2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=                                                          
     *    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+I))))
     *   +C2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3                                                          
      J=J+INC4                                                          
  330 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  340 CONTINUE                                                          
      IA=IA+IINK                                                        
      IB=IB+IINK                                                        
      IC=IC-IINK                                                        
      JBASE=JBASE+JUMP                                                  
  350 CONTINUE                                                          
      ENDIF
C
      IF (IA.EQ.IC) THEN                                           
      IBASE=0                                                           
      DO 380 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 370 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      C(JB+J)=(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))                     
      C(JC+J)=-(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))                    
      I=I+INC3                                                          
      J=J+INC4                                                          
  370 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  380 CONTINUE                                                          
      ENDIF                                                         
C                                                                       
      ELSE                !!! Case LA=M                                 
      SSIN60=2.0*SIN60                                                  
      IF (LIPL) THEN
        DO 394 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 392 IJK=1,ILOT                                                  
        T1=(2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))                    
        T2=(2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))                    
        A(IA+I)=2.0*(A(IA+I)+A(IB+I))                                     
        A(IB+I)=T1
        B(IB+I)=T2
        I=I+INC3                                                          
  392   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  394   CONTINUE                                                          
      ELSE
        DO 398 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 396 IJK=1,ILOT                                                  
        C(JA+J)=2.0*(A(IA+I)+A(IB+I))                                     
        C(JB+J)=(2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))                    
        C(JC+J)=(2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))                    
        I=I+INC3                                                          
        J=J+INC4                                                          
  396   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  398   CONTINUE                                                          
      ENDIF
      ENDIF
C
      ELSEIF (IFAC.EQ.4) THEN
C                                                                       
C     CODING FOR FACTOR 4                                               
C     -------------------                                               
  400 CONTINUE                                                          
      IA=1                                                              
      IB=IA+(2*M-LA)*INC1                                               
      IC=IB+2*M*INC1                                                    
      ID=IB                                                             
      JA=1                                                              
      JB=JA+JINK                                                        
      JC=JB+JINK                                                        
      JD=JC+JINK                                                        
C                                                                       
      IF (LA.NE.M) THEN                                            
C                                                                       
      DO 420 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 410 IJK=1,ILOT                                                  
      C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)                                 
      C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)                                 
      C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)                                 
      C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)                                 
      I=I+INC3                                                          
      J=J+INC4                                                          
  410 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  420 CONTINUE                                                          
      IA=IA+IINK                                                        
      IINK=2*IINK                                                       
      IB=IB+IINK                                                        
      IC=IC-IINK                                                        
      ID=ID-IINK                                                        
      JBASE=JBASE+JUMP                                                  
      JUMP=2*JUMP+JINK                                                  
C
      IF (IB.LT.IC) THEN                                           
      DO 450 K=LA,KSTOP,LA                                              
      KB=K+K                                                            
      KC=KB+KB                                                          
      KD=KC+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      C3=TRIGS(KD+1)                                                    
      S3=TRIGS(KD+2)                                                    
      IBASE=0                                                           
      DO 440 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 430 IJK=1,ILOT                                                  
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))                       
      D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))                       
      C(JC+J)=                                                          
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      
     *   -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))                      
      D(JC+J)=                                                          
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      
     *   +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))                      
      C(JB+J)=                                                          
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))                      
     *   -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))                      
      D(JB+J)=                                                          
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))                      
     *   +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))                      
      C(JD+J)=                                                          
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))                      
     *   -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))                      
      D(JD+J)=                                                          
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))                      
     *   +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))                      
      I=I+INC3                                                          
      J=J+INC4                                                          
  430 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  440 CONTINUE                                                          
      IA=IA+IINK                                                        
      IB=IB+IINK                                                        
      IC=IC-IINK                                                        
      ID=ID-IINK                                                        
      JBASE=JBASE+JUMP                                                  
  450 CONTINUE                                                          
      ENDIF
C
      IF (IB.EQ.IC) THEN                                           
      IBASE=0                                                           
      SIN45=SQRT(0.5)                                                   
      DO 480 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 470 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))               
      C(JC+J)=B(IB+I)-B(IA+I)                                           
      C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))              
      I=I+INC3                                                          
      J=J+INC4                                                          
  470 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  480 CONTINUE                                                          
      ENDIF                                                             
C                                                                       
      ELSE                !!! Case LA=M                                 
      IF (LIPL) THEN
        DO 494 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 492 IJK=1,ILOT                                                  
        T1=2.0*((A(IA+I)-A(IC+I))-B(IB+I))                           
        T2=2.0*((A(IA+I)+A(IC+I))-A(IB+I))                           
        T3=2.0*((A(IA+I)-A(IC+I))+B(IB+I))                           
        A(IA+I)=2.0*((A(IA+I)+A(IC+I))+A(IB+I))                           
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        I=I+INC3                                                          
  492   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  494   CONTINUE                                                          
      ELSE
        DO 498 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 496 IJK=1,ILOT                                                  
        C(JA+J)=2.0*((A(IA+I)+A(IC+I))+A(IB+I))                           
        C(JB+J)=2.0*((A(IA+I)-A(IC+I))-B(IB+I))                           
        C(JC+J)=2.0*((A(IA+I)+A(IC+I))-A(IB+I))                           
        C(JD+J)=2.0*((A(IA+I)-A(IC+I))+B(IB+I))                           
        I=I+INC3                                                          
        J=J+INC4                                                          
  496   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  498   CONTINUE                                                          
      ENDIF
      ENDIF
C
      ELSEIF (IFAC.EQ.5) THEN
C                                                                       
C     CODING FOR FACTOR 5                                               
C     -------------------                                               
  500 CONTINUE                                                          
      IA=1                                                              
      IB=IA+(2*M-LA)*INC1                                               
      IC=IB+2*M*INC1                                                    
      ID=IC                                                             
      IE=IB                                                             
      JA=1                                                              
      JB=JA+JINK                                                        
      JC=JB+JINK                                                        
      JD=JC+JINK                                                        
      JE=JD+JINK                                                        
C                                                                       
      IF (LA.NE.M) THEN                                                 
C                                                                       
      DO 520 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 510 IJK=1,ILOT                                                  
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))                                 
      C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I))) 
     *    -(SIN72*B(IB+I)+SIN36*B(IC+I))                                
      C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I))) 
     *    -(SIN36*B(IB+I)-SIN72*B(IC+I))                                
      C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I))) 
     *    +(SIN36*B(IB+I)-SIN72*B(IC+I))                                
      C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I))) 
     *    +(SIN72*B(IB+I)+SIN36*B(IC+I))                                
      I=I+INC3                                                          
      J=J+INC4                                                          
  510 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  520 CONTINUE                                                          
      IA=IA+IINK                                                        
      IINK=2*IINK                                                       
      IB=IB+IINK                                                        
      IC=IC+IINK                                                        
      ID=ID-IINK                                                        
      IE=IE-IINK                                                        
      JBASE=JBASE+JUMP                                                  
      JUMP=2*JUMP+JINK                                                  
C
      IF (IB.LT.ID) THEN                                                
      DO 550 K=LA,KSTOP,LA                                              
      KB=K+K                                                            
      KC=KB+KB                                                          
      KD=KC+KB                                                          
      KE=KD+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      C3=TRIGS(KD+1)                                                    
      S3=TRIGS(KD+2)                                                    
      C4=TRIGS(KE+1)                                                    
      S4=TRIGS(KE+2)                                                    
      IBASE=0                                                           
      DO 540 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 530 IJK=1,ILOT                                                  
C                                                                       
      A10=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))          
     *    +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))                   
      A20=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))          
     *    -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))                   
      B10=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))          
     *    +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))                   
      B20=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))          
     *    -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))                   
      A11=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))               
      A21=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))               
      B11=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))               
      B21=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))               
C                                                                       
      C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))             
      D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))             
      C(JB+J)=C1*(A10-A11)-S1*(B10+B11)                                 
      D(JB+J)=S1*(A10-A11)+C1*(B10+B11)                                 
      C(JE+J)=C4*(A10+A11)-S4*(B10-B11)                                 
      D(JE+J)=S4*(A10+A11)+C4*(B10-B11)                                 
      C(JC+J)=C2*(A20-A21)-S2*(B20+B21)                                 
      D(JC+J)=S2*(A20-A21)+C2*(B20+B21)                                 
      C(JD+J)=C3*(A20+A21)-S3*(B20-B21)                                 
      D(JD+J)=S3*(A20+A21)+C3*(B20-B21)                                 
C                                                                       
      I=I+INC3                                                          
      J=J+INC4                                                          
  530 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  540 CONTINUE                                                          
      IA=IA+IINK                                                        
      IB=IB+IINK                                                        
      IC=IC+IINK                                                        
      ID=ID-IINK                                                        
      IE=IE-IINK                                                        
      JBASE=JBASE+JUMP                                                  
  550 CONTINUE                                                          
      ENDIF
C
      IF (IB.EQ.ID) THEN                                                
      IBASE=0                                                           
      DO 580 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 570 IJK=1,ILOT                                                  
      C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)                                 
      C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) 
     *    -(SIN36*B(IA+I)+SIN72*B(IB+I))                                
      C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     *    -(SIN36*B(IA+I)+SIN72*B(IB+I))                                
      C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) 
     *    -(SIN72*B(IA+I)-SIN36*B(IB+I))                                
      C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     *    -(SIN72*B(IA+I)-SIN36*B(IB+I))                                
      I=I+INC3                                                          
      J=J+INC4                                                          
  570 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  580 CONTINUE                                                          
      ENDIF                                                             
C                                                                       
      ELSE               !!! Case LA=M                                  
      QQRT5=2.0*QRT5                                                    
      SSIN36=2.0*SIN36                                                  
      SSIN72=2.0*SIN72                                                  
      IF (LIPL) THEN
        DO 594 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 592 IJK=1,ILOT                                                  
        T1=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))     
        T2=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))     
        T3=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))     
        T4=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))     
        A(IA+I)=2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))                           
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        I=I+INC3                                                          
  592   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  594   CONTINUE                                                          
      ELSE
        DO 598 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 596 IJK=1,ILOT                                                  
        C(JA+J)=2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))                           
        C(JB+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))     
        C(JC+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))     
        C(JD+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))     
        C(JE+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))                     
     *    +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))     
        I=I+INC3                                                          
        J=J+INC4                                                          
  596   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  598   CONTINUE                                                          
      ENDIF
      ENDIF
C
      ELSEIF (IFAC.EQ.6) THEN
C                                                                       
C     CODING FOR FACTOR 6                                               
C     -------------------                                               
  600 CONTINUE                                                          
      IA=1                                                              
      IB=IA+(2*M-LA)*INC1                                               
      IC=IB+2*M*INC1                                                    
      ID=IC+2*M*INC1                                                    
      IE=IC                                                             
      IF=IB                                                             
      JA=1                                                              
      JB=JA+JINK                                                        
      JC=JB+JINK                                                        
      JD=JC+JINK                                                        
      JE=JD+JINK                                                        
      JF=JE+JINK                                                        
C                                                                       
      IF (LA.NE.M) THEN                                                 
C                                                                       
      DO 620 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 610 IJK=1,ILOT                                                  
      C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))                       
      C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))                       
      C(JB+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))                 
     *    -(SIN60*(B(IB+I)+B(IC+I)))                                    
      C(JF+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))                 
     *    +(SIN60*(B(IB+I)+B(IC+I)))                                    
      C(JC+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))                 
     *    -(SIN60*(B(IB+I)-B(IC+I)))                                    
      C(JE+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))                 
     *    +(SIN60*(B(IB+I)-B(IC+I)))                                    
      I=I+INC3                                                          
      J=J+INC4                                                          
  610 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  620 CONTINUE                                                          
      IA=IA+IINK                                                        
      IINK=2*IINK                                                       
      IB=IB+IINK                                                        
      IC=IC+IINK                                                        
      ID=ID-IINK                                                        
      IE=IE-IINK                                                        
      IF=IF-IINK                                                        
      JBASE=JBASE+JUMP                                                  
      JUMP=2*JUMP+JINK                                                  
C
      IF (IC.LT.ID) THEN                                                
      DO 650 K=LA,KSTOP,LA                                              
      KB=K+K                                                            
      KC=KB+KB                                                          
      KD=KC+KB                                                          
      KE=KD+KB                                                          
      KF=KE+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      C3=TRIGS(KD+1)                                                    
      S3=TRIGS(KD+2)                                                    
      C4=TRIGS(KE+1)                                                    
      S4=TRIGS(KE+2)                                                    
      C5=TRIGS(KF+1)                                                    
      S5=TRIGS(KF+2)                                                    
      IBASE=0                                                           
      DO 640 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 630 IJK=1,ILOT                                                  
C                                                                       
      A11= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))                          
      A20=(A(IA+I)+A(ID+I))-0.5*A11                                     
      A21=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))                   
      B11= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))                          
      B20=(B(IA+I)-B(ID+I))-0.5*B11                                     
      B21=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))                   
C                                                                       
      C(JA+J)=(A(IA+I)+A(ID+I))+A11                                     
      D(JA+J)=(B(IA+I)-B(ID+I))+B11                                     
      C(JC+J)=C2*(A20-B21)-S2*(B20+A21)                                 
      D(JC+J)=S2*(A20-B21)+C2*(B20+A21)                                 
      C(JE+J)=C4*(A20+B21)-S4*(B20-A21)                                 
      D(JE+J)=S4*(A20+B21)+C4*(B20-A21)                                 
C                                                                       
      A11=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))                           
      B11=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))                           
      A20=(A(IA+I)-A(ID+I))-0.5*A11                                     
      A21=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))                   
      B20=(B(IA+I)+B(ID+I))+0.5*B11                                     
      B21=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))                   
C                                                                       
      C(JD+J)=                                                          
     *  C3*((A(IA+I)-A(ID+I))+A11)-S3*((B(IA+I)+B(ID+I))-B11)           
      D(JD+J)=                                                          
     *  S3*((A(IA+I)-A(ID+I))+A11)+C3*((B(IA+I)+B(ID+I))-B11)           
      C(JB+J)=C1*(A20-B21)-S1*(B20-A21)                                 
      D(JB+J)=S1*(A20-B21)+C1*(B20-A21)                                 
      C(JF+J)=C5*(A20+B21)-S5*(B20+A21)                                 
      D(JF+J)=S5*(A20+B21)+C5*(B20+A21)                                 
C                                                                       
      I=I+INC3                                                          
      J=J+INC4                                                          
  630 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  640 CONTINUE                                                          
      IA=IA+IINK                                                        
      IB=IB+IINK                                                        
      IC=IC+IINK                                                        
      ID=ID-IINK                                                        
      IE=IE-IINK                                                        
      IF=IF-IINK                                                        
      JBASE=JBASE+JUMP                                                  
  650 CONTINUE                                                          
      ENDIF
C
      IF (IC.EQ.ID) THEN                                                
      IBASE=0                                                           
      DO 680 L=1,ILA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
!OCL NOVREC                                                   
      DO 670 IJK=1,ILOT                                                  
      C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))                                 
      C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))                                 
      C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I)) 
      C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
      C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+(0.5*(A(IA+I)+A(IC+I))-A(IB+I))   
      C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-(0.5*(A(IA+I)+A(IC+I))-A(IB+I))   
      I=I+INC3                                                          
      J=J+INC4                                                          
  670 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC21                                                  
  680 CONTINUE                                                          
      ENDIF                                                             
C                                                                       
      ELSE                 !!! Case LA=M                                
      SSIN60=2.0*SIN60                                                  
      IF (LIPL) THEN
        DO 694 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 692 IJK=1,ILOT                                                  
        T1=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))                 
     *    -(SSIN60*(B(IB+I)+B(IC+I)))                                   
        T5=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))                 
     *    +(SSIN60*(B(IB+I)+B(IC+I)))                                   
        T2=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))                 
     *    -(SSIN60*(B(IB+I)-B(IC+I)))                                   
        T4=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))                 
     *    +(SSIN60*(B(IB+I)-B(IC+I)))                                   
        T3=(2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))           
        A(IA+I)=(2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))           
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        A(ID+I)=T5
        I=I+INC3                                                          
  692   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  694   CONTINUE                                                          
      ELSE
        DO 698 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 696 IJK=1,ILOT                                                  
        C(JA+J)=(2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))           
        C(JD+J)=(2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))           
        C(JB+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))                 
     *    -(SSIN60*(B(IB+I)+B(IC+I)))                                   
        C(JF+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))                 
     *    +(SSIN60*(B(IB+I)+B(IC+I)))                                   
        C(JC+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))                 
     *    -(SSIN60*(B(IB+I)-B(IC+I)))                                   
        C(JE+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))                 
     *    +(SSIN60*(B(IB+I)-B(IC+I)))                                   
        I=I+INC3                                                          
        J=J+INC4                                                          
  696   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  698   CONTINUE                                                          
      ENDIF
      ENDIF
C
      ELSEIF (IFAC.EQ.8) THEN
C                                                                       
C     CODING FOR FACTOR 8                                               
C     -------------------                                               
  800 CONTINUE                                                          
      IF (LA.NE.M) THEN                                            
        IBAD=3
      ELSE
      IA=1                                                              
      IB=IA+LA*INC1                                                     
      IC=IB+2*LA*INC1                                                   
      ID=IC+2*LA*INC1                                                   
      IE=ID+2*LA*INC1                                                   
      JA=1                                                              
      JB=JA+JINK                                                        
      JC=JB+JINK                                                        
      JD=JC+JINK                                                        
      JE=JD+JINK                                                        
      JF=JE+JINK                                                        
      JG=JF+JINK                                                        
      JH=JG+JINK                                                        
      SSIN45=SQRT(2.0)                                                  
C                                                                       
      IF (LIPL) THEN
        DO 820 L=1,ILA                                                     
        I=IBASE                                                           
!OCL NOVREC
        DO 810 IJK=1,ILOT                                                  
        T2=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))       
        T6=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))       
        T1=2.0*((A(IA+I)-A(IE+I))-B(IC+I))                           
     *    +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))                 
        T5=2.0*((A(IA+I)-A(IE+I))-B(IC+I))                           
     *    -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))                 
        T3=2.0*((A(IA+I)-A(IE+I))+B(IC+I))                           
     *    -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))                 
        T7=2.0*((A(IA+I)-A(IE+I))+B(IC+I))                           
     *    +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))                 
        T4=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))       
        A(IA+I)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))       
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        A(ID+I)=T5
        B(ID+I)=T6
        A(IE+I)=T7  
        I=I+INC3                                                          
  810   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
  820   CONTINUE                                                          
      ELSE
        DO 840 L=1,ILA                                                     
        I=IBASE                                                           
        J=JBASE                                                           
!OCL NOVREC
        DO 830 IJK=1,ILOT                                                  
        C(JA+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))       
        C(JE+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))       
        C(JC+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))       
        C(JG+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))       
        C(JB+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I))                           
     *    +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))                 
        C(JF+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I))                           
     *    -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))                 
        C(JD+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I))                           
     *    -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))                 
        C(JH+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I))                           
     *    +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))                 
        I=I+INC3                                                          
        J=J+INC4                                                          
  830   CONTINUE                                                          
        IBASE=IBASE+INC1                                                  
        JBASE=JBASE+INC21                                                  
  840   CONTINUE                                                          
      ENDIF
C
      ENDIF
C
      ELSE
C
      IBAD=2       !!! Illegal factor
C
      ENDIF
C                                                                       
C     RETURN                                                            
C     ------                                                            
  900 CONTINUE                                                          
      IERR=IBAD                                                         
      RETURN                                                            
      END                                                               
