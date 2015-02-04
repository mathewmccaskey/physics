ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP()

      IMPLICIT NONE
      DOUBLE PRECISION PI
      LOGICAL READLHA
      PARAMETER  (PI=3.141592653589793D0)

      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      READLHA = .TRUE.
      INCLUDE 'intparam_definition.inc'



      CALL COUP1()
      CALL COUP2()
      CALL COUP3()
C     
C     couplings needed to be evaluated points by points
C     
      CALL COUP4()

      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM()

      IMPLICIT NONE
      DOUBLE PRECISION PI
      LOGICAL READLHA
      PARAMETER  (PI=3.141592653589793D0)

      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      READLHA = .FALSE.
      INCLUDE 'intparam_definition.inc'



C     
C     couplings needed to be evaluated points by points
C     
      CALL COUP4()

      RETURN
      END

