C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     (ProjM(2,1)) * C(51,2) * C(52,1)
C     
      SUBROUTINE FFS1C1_0(F2, F1, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP15
      COMPLEX*16 S3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      TMP15 = (F1(3)*F2(3)+F1(4)*F2(4))
      VERTEX = COUP*-CI * TMP15*S3(3)
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     (ProjM(2,1)) * C(51,2) * C(52,1)
C     
      SUBROUTINE FFS1_2C1_0(F2, F1, S3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 S3(*)
      COMPLEX*16 COUP2
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP
      CALL FFS1C1_0(F2,F1,S3,COUP1,VERTEX)
      CALL FFS2C1_0(F2,F1,S3,COUP2,TMP)
      VERTEX = VERTEX + TMP
      END

