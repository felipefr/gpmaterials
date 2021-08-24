MODULE DETERMINANT
! Used in matrixcoef and in some elements and in writeElemOut, send to elements the specified variables
! Verified if used in jacobian

implicit none
REAL*8  :: DET
INTEGER :: LengthParam,LengthJParam,iElem,IdElemLib
INTEGER :: iWriteElemUnit,Id_Write, OUnit_residualSensitivity , OUnit_residualStress        

!!! good version for SGP
real*8 :: dummy
COMMON dummy,LengthParam,LengthJParam
!~ 		
END MODULE DETERMINANT

MODULE TIME_STEP_ITER
! Used in maingeneric and in some elements, send to elements the specified variables 
       ! REAL*8 :: 
        INTEGER :: Nstep,nstp,Nconverg,Nconvers
END MODULE TIME_STEP_ITER
