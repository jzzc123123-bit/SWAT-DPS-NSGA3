      SUBROUTINE DPS_RBF_EVAL(B,K,N,S,LO,HI,A,W,C,R,
     &                        Z,PHI,RAW,U)
!!    Gaussian policy matching SWAT-DPS-NSGA3 (negative exponent,
!!    one width per state and basis, no factor 1/2 or normalization).
!!    SUBROUTINES/FUNCTIONS CALLED: EXP, MIN, MAX
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: B,K,N
      REAL(KIND=8), INTENT(IN) :: S(B),LO(B),HI(B),A(K)
      REAL(KIND=8), INTENT(IN) :: W(N,K),C(N,B),R(N,B)
      REAL(KIND=8), INTENT(OUT) :: Z(B),PHI(N),RAW(K)
      REAL(KIND=8), INTENT(OUT) :: U(K)
      INTEGER :: IB,JB,KB
      REAL(KIND=8) :: DIST2
      Z=MAX(0.D0,MIN(1.D0,(S-LO)/(HI-LO)))
      DO IB=1,N
        DIST2=0.D0
        DO JB=1,B
          DIST2=DIST2+((Z(JB)-C(IB,JB))/R(IB,JB))**2
        ENDDO
        PHI(IB)=EXP(-DIST2)
      ENDDO
      RAW=A
      DO KB=1,K
        DO IB=1,N
          RAW(KB)=RAW(KB)+W(IB,KB)*PHI(IB)
        ENDDO
      ENDDO
      U=MAX(0.D0,MIN(1.D0,RAW))
      END
