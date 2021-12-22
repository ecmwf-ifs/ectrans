INTERFACE
 SUBROUTINE SUBLAYER(pz1,pz2,pz3,px1,px2,ldgradps,&
 & pt1,pt2,pw1,pw2,pzs1,pzs2,pzps) 
 USE PARKIND1 ,ONLY : JPIM ,JPRB
 LOGICAL, INTENT(in) :: ldgradps
 REAL(KIND=JPRB), INTENT(in) :: pz1,pz2,pz3,px1,px2,pt1,pt2,pzs1,pzs2
 REAL(KIND=JPRB), INTENT(out) :: pw1,pw2,pzps
 END SUBROUTINE SUBLAYER
END INTERFACE
