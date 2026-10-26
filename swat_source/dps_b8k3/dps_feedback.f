!!    Subbasin year-start RBF feedback for fertilizer application.
!!    SUBROUTINES/FUNCTIONS CALLED: DPS_FAIL, DPS_RBF_EVAL
      SUBROUTINE DPS_INIT
      USE PARM
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER :: JH,IOP,IOS,NAG,IB,JB,KB,FB,FK,FN,FP,ISB,IPH
      LOGICAL :: PRESENT
      CHARACTER(LEN=128) :: HEADER

      DPS_ENABLED=.FALSE.
      DPS_HAS_STATE=.FALSE.
      DPS_STATE_YEAR=0
      INQUIRE(FILE='dps_policy.in',EXIST=PRESENT)
      IF (.NOT.PRESENT) RETURN
      OPEN(196,FILE='dps_policy.in',STATUS='OLD',IOSTAT=IOS)
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_OPEN')
      READ(196,'(A)',IOSTAT=IOS) HEADER
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_HEADER_READ')
      IF (TRIM(HEADER) /= 'DPS_SUBBASIN_SPLITN_V5') THEN
        CALL DPS_FAIL('POLICY_EXPECTED_DPS_SUBBASIN_SPLITN_V5')
      ENDIF
      READ(196,*,IOSTAT=IOS) DPS_MODE,DPS_START_YEAR,DPS_OUTLET
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_MODE')
      READ(196,*,IOSTAT=IOS) FB,FK,FN,FP
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_DIMENSION_READ')
      IF (FB/=DPS_B.OR.FK/=DPS_K.OR.FN/=DPS_N.OR.FP/=DPS_P)
     &  CALL DPS_FAIL('POLICY_EXPECTED_B8_K3_N12_P231')
      READ(196,*,IOSTAT=IOS) DPS_ANNUAL_LO,DPS_ANNUAL_HI
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_DOSE_BOUNDS')
      READ(196,*,IOSTAT=IOS) DPS_PHASE_LO
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_PHASE_LO')
      READ(196,*,IOSTAT=IOS) DPS_PHASE_HI
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_PHASE_HI')
      READ(196,*,IOSTAT=IOS) DPS_PHASE_REF
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_PHASE_REFERENCE')
      READ(196,*,IOSTAT=IOS) DPS_LO
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_STATE_LO')
      READ(196,*,IOSTAT=IOS) DPS_HI
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_STATE_HI')
      READ(196,*,IOSTAT=IOS) DPS_A
      IF (IOS /= 0) CALL DPS_FAIL('POLICY_BIAS')
      DO IB=1,DPS_N
        READ(196,*,IOSTAT=IOS) (DPS_W(IB,KB),KB=1,DPS_K)
        IF (IOS /= 0) CALL DPS_FAIL('POLICY_WEIGHT')
      ENDDO
      DO IB=1,DPS_N
        READ(196,*,IOSTAT=IOS) (DPS_C(IB,JB),JB=1,DPS_B)
        IF (IOS /= 0) CALL DPS_FAIL('POLICY_CENTER')
      ENDDO
      DO IB=1,DPS_N
        READ(196,*,IOSTAT=IOS) (DPS_R(IB,JB),JB=1,DPS_B)
        IF (IOS /= 0) CALL DPS_FAIL('POLICY_RADIUS')
      ENDDO
      CLOSE(196)
      IF (DPS_MODE<0.OR.DPS_MODE>1) CALL DPS_FAIL('MODE_RANGE')
      IF (DPS_OUTLET<1.OR.DPS_OUTLET>SUBTOT)
     &  CALL DPS_FAIL('OUTLET_RANGE')
      IF (ANY(.NOT.IEEE_IS_FINITE(DPS_LO)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_HI)).OR.
     &    .NOT.IEEE_IS_FINITE(DPS_ANNUAL_LO).OR.
     &    .NOT.IEEE_IS_FINITE(DPS_ANNUAL_HI).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_PHASE_LO)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_PHASE_HI)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_PHASE_REF)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_A)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_W)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_C)).OR.
     &    ANY(.NOT.IEEE_IS_FINITE(DPS_R)))
     &  CALL DPS_FAIL('POLICY_NONFINITE')
      IF (ANY(DPS_HI<=DPS_LO)) CALL DPS_FAIL('STATE_BOUNDS')
      IF (DPS_ANNUAL_LO<0.D0.OR.
     &    DPS_ANNUAL_HI<=DPS_ANNUAL_LO.OR.
     &    DPS_ANNUAL_LO>1.D0.OR.DPS_ANNUAL_HI<1.D0)
     &  CALL DPS_FAIL('DOSE_BOUNDS_MUST_INCLUDE_BASELINE')

      IF (ANY(DPS_PHASE_LO<0.D0).OR.
     &    ANY(DPS_PHASE_HI>1.D0).OR.
     &    ANY(DPS_PHASE_HI<DPS_PHASE_LO).OR.
     &    SUM(DPS_PHASE_LO)>1.D0.OR.SUM(DPS_PHASE_HI)<1.D0)
     &  CALL DPS_FAIL('INFEASIBLE_PHASE_BOUNDS')
      IF (ANY(DPS_PHASE_REF<DPS_PHASE_LO).OR.
     &    ANY(DPS_PHASE_REF>DPS_PHASE_HI).OR.
     &    ABS(SUM(DPS_PHASE_REF)-1.D0)>1.D-12)
     &  CALL DPS_FAIL('INVALID_PHASE_REFERENCE')
      IF (ANY(DPS_A<0.D0).OR.ANY(DPS_A>1.D0))
     &  CALL DPS_FAIL('BIAS_RANGE')
      IF (ANY(ABS(DPS_W)>1.D0)) CALL DPS_FAIL('WEIGHT_RANGE')
      IF (ANY(DPS_C<0.D0).OR.ANY(DPS_C>1.D0))
     &  CALL DPS_FAIL('CENTER_RANGE')
      IF (ANY(DPS_R<0.01D0).OR.ANY(DPS_R>1.D0))
     &  CALL DPS_FAIL('RADIUS_RANGE')
      IF (ALLOCATED(DPS_AGMASK)) DEALLOCATE(DPS_AGMASK)
      IF (ALLOCATED(DPS_STATE)) DEALLOCATE(DPS_STATE)
      IF (ALLOCATED(DPS_AG_AREA)) DEALLOCATE(DPS_AG_AREA)
      IF (ALLOCATED(DPS_AG_COUNT)) DEALLOCATE(DPS_AG_COUNT)
      IF (ALLOCATED(DPS_CHEM_SCALE)) DEALLOCATE(DPS_CHEM_SCALE)
      IF (ALLOCATED(DPS_MAN_SCALE)) DEALLOCATE(DPS_MAN_SCALE)
      ALLOCATE(DPS_AGMASK(NHRU))
      ALLOCATE(DPS_STATE(DPS_B,SUBTOT),DPS_AG_AREA(SUBTOT))
      ALLOCATE(DPS_AG_COUNT(SUBTOT),
     &         DPS_CHEM_SCALE(DPS_K,SUBTOT),
     &         DPS_MAN_SCALE(SUBTOT))
      DPS_AGMASK=.FALSE.
      DPS_STATE=0.D0
      DPS_AG_AREA=0.D0
      DPS_AG_COUNT=0
      DPS_CHEM_SCALE=1.
      DPS_MAN_SCALE=1.
      DO JH=1,NHRU
        DO IOP=1,NOPMX(JH)
          IF (MGTOP(IOP,JH)/=3) CYCLE
          IF (MGT1IOP(IOP,JH)==4.OR.
     &        MGT1IOP(IOP,JH)==44.OR.
     &        MGT1IOP(IOP,JH)==45.OR.
     &        MGT1IOP(IOP,JH)==47.OR.
     &        MGT1IOP(IOP,JH)==51) DPS_AGMASK(JH)=.TRUE.
        ENDDO
        IF (DPS_AGMASK(JH)) THEN
          ISB=HRU_SUB(JH)
          IF (ISB<1.OR.ISB>SUBTOT)
     &      CALL DPS_FAIL('AG_HRU_SUBBASIN_RANGE')
        ENDIF
      ENDDO
      NAG=COUNT(DPS_AGMASK)
      IF (NAG==0) CALL DPS_FAIL('NO_AGRICULTURAL_HRU')
      ALLOCATE(DPS_BASE_CHEM_N(SUBTOT))
      ALLOCATE(DPS_BASE_PHASE_N(DPS_K,SUBTOT))
      ALLOCATE(DPS_CHEM_PHASE(SIZE(MGTOP,1),NHRU))
      DPS_BASE_CHEM_N=0.D0
      DPS_BASE_PHASE_N=0.D0
      DPS_CHEM_PHASE=0
      IF (.NOT.IEEE_IS_FINITE(FMINN(4)+FORGN(4)))
     &  CALL DPS_FAIL('CHEMICAL_N_FRACTION_NONFINITE')
      IF (FMINN(4)+FORGN(4)<=0.)
     &  CALL DPS_FAIL('CHEMICAL_N_FRACTION_INVALID')
      DO JH=1,NHRU
        IF (.NOT.DPS_AGMASK(JH)) CYCLE
        ISB=HRU_SUB(JH)
        DPS_AG_AREA(ISB)=DPS_AG_AREA(ISB)+DBLE(HRU_KM(JH))
        IPH=0
        DO IOP=1,NOPMX(JH)
          IF (MGTOP(IOP,JH)/=3) CYCLE
          IF (MGT1IOP(IOP,JH)/=4) CYCLE
          IPH=IPH+1
          IF (IPH>DPS_K)
     &      CALL DPS_FAIL('MORE_THAN_THREE_CHEMICAL_EVENTS')
          DPS_CHEM_PHASE(IOP,JH)=IPH
          IF (.NOT.IEEE_IS_FINITE(MGT4OP(IOP,JH)))
     &      CALL DPS_FAIL('BASELINE_CHEMICAL_DOSE_NONFINITE')
          IF (MGT4OP(IOP,JH)<=0.)
     &      CALL DPS_FAIL('BASELINE_CHEMICAL_DOSE_INVALID')
          DPS_BASE_PHASE_N(IPH,ISB)=
     &      DPS_BASE_PHASE_N(IPH,ISB)+DBLE(HRU_KM(JH))*
     &      DBLE(MGT4OP(IOP,JH))*DBLE(FMINN(4)+FORGN(4))
        ENDDO
        IF (IPH/=DPS_K)
     &    CALL DPS_FAIL('EXPECTED_THREE_CHEMICAL_EVENTS')
      ENDDO
      DO ISB=1,SUBTOT
        IF (DPS_AG_AREA(ISB)>0.D0) THEN
          DPS_BASE_PHASE_N(:,ISB)=DPS_BASE_PHASE_N(:,ISB)/
     &      DPS_AG_AREA(ISB)
          IF (ANY(DPS_BASE_PHASE_N(:,ISB)<=0.D0))
     &      CALL DPS_FAIL('SUBBASIN_BASE_PHASE_INVALID')
        ENDIF
        DPS_BASE_CHEM_N(ISB)=SUM(DPS_BASE_PHASE_N(:,ISB))
      ENDDO
      DPS_ENABLED=.TRUE.
      OPEN(196,FILE='dps_actions.tsv',STATUS='REPLACE')
      WRITE(196,'(A)') 'year subbasin state_year active'//
     & ' ag_hru_count state_oob_count budget_clip zero_fallback'//
     & ' share_repaired'//
     & ' ag_area_km2 chem_reduction manure_reduction'//
     & ' chem_scale manure_scale'//
     & ' chem_baseline_N_kg_ha_year chem_N_kg_ha_year'//
     & ' baseline_basal_N baseline_topdress1_N baseline_topdress2_N'//
     & ' basal_N topdress1_N topdress2_N'//
     & ' basal_scale topdress1_scale topdress2_scale'//
     & ' trial_basal_share trial_topdress1_share'//
     & ' trial_topdress2_share basal_share topdress1_share'//
     & ' topdress2_share'//
     & ' soil_nitrate_N soil_ammonium_N soil_fresh_org_N'//
     & ' soil_active_org_N soil_stable_org_N aquifer_NO3_N'//
     & ' soil_water shallow_aquifer_water'//
     & ' z1 z2 z3 z4 z5 z6 z7 z8 raw1 raw2 raw3 u1 u2 u3'//
     & ' phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 phi9 phi10'//
     & ' phi11 phi12'
      CLOSE(196)
      OPEN(197,FILE='dps_year_start_state.tsv',STATUS='REPLACE')
      WRITE(197,'(A)') 'year subbasin ag_hru_count ag_area_km2'//
     & ' soil_nitrate_N soil_ammonium_N soil_fresh_org_N'//
     & ' soil_active_org_N soil_stable_org_N aquifer_NO3_N'//
     & ' soil_water shallow_aquifer_water'
      CLOSE(197)
      OPEN(194,FILE='dps_hru_year_start_state.tsv',
     & STATUS='REPLACE')
      WRITE(194,'(A)') 'year hru subbasin area_km2'//
     & ' soil_nitrate_N soil_ammonium_N soil_fresh_org_N'//
     & ' soil_active_org_N soil_stable_org_N aquifer_NO3_N'//
     & ' soil_water shallow_aquifer_water'
      CLOSE(194)
      OPEN(195,FILE='dps_outcome.tsv',STATUS='REPLACE')
      WRITE(195,'(A)') 'year outlet_TN_kg outlet_flow_m3_s'
      CLOSE(195)
      WRITE(*,*) 'DPS shared subbasin RBF B,K,N,P: ',
     & DPS_B,DPS_K,DPS_N,DPS_P
      WRITE(*,*) 'DPS agricultural HRUs: ',NAG
      END

      SUBROUTINE DPS_YEAR_START
      USE PARM
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER :: ACTIVE,OOB,ISB,BCLIP,ZEROFLAG,REPAIR,ITER
      REAL(KIND=8) :: SUMU,TRIAL,TOTAL,AFAC,NTARGET(DPS_K)
      REAL(KIND=8) :: ORIGSHARE(DPS_K),SHARE(DPS_K)
      REAL(KIND=8) :: TAULO,TAUHI,TAU
      IF (.NOT.DPS_ENABLED) RETURN
      IF (.NOT.DPS_HAS_STATE.OR.DPS_STATE_YEAR/=IYR)
     &  CALL DPS_FAIL('YEAR_START_STATE_MISSING')
!!    Each chemical event uses its own fixed annual phase scale.
      DPS_CHEM_SCALE=1.
      DPS_MAN_SCALE=1.
      OPEN(196,FILE='dps_actions.tsv',POSITION='APPEND')
      DO ISB=1,SUBTOT
        DPS_RHO=0.D0
        DPS_FAC=1.D0
        DPS_RAW=0.D0
        DPS_U=0.D0
        DPS_PHI=0.D0
        DPS_Z=0.D0
        NTARGET=DPS_BASE_PHASE_N(:,ISB)
        ORIGSHARE=0.D0
        IF (DPS_BASE_CHEM_N(ISB)>0.D0)
     &    ORIGSHARE=NTARGET/DPS_BASE_CHEM_N(ISB)
        SHARE=ORIGSHARE
        REPAIR=0
        AFAC=1.D0
        OOB=0
        ACTIVE=0
        BCLIP=0
        ZEROFLAG=0
        IF (DPS_AG_COUNT(ISB)>0) THEN
          DPS_Z=MAX(0.D0,MIN(1.D0,
     &      (DPS_STATE(:,ISB)-DPS_LO)/(DPS_HI-DPS_LO)))
          OOB=COUNT(DPS_STATE(:,ISB)<DPS_LO.OR.
     &              DPS_STATE(:,ISB)>DPS_HI)
          IF (DPS_MODE==1.AND.IYR>=DPS_START_YEAR) THEN
            ACTIVE=1
            CALL DPS_RBF_EVAL(DPS_B,DPS_K,DPS_N,
     &        DPS_STATE(:,ISB),DPS_LO,DPS_HI,DPS_A,DPS_W,
     &        DPS_C,DPS_R,DPS_Z,DPS_PHI,DPS_RAW,DPS_U)
            SUMU=SUM(DPS_U)
            TRIAL=DPS_ANNUAL_HI*SUMU
            IF (TRIAL<DPS_ANNUAL_LO) BCLIP=-1
            IF (TRIAL>DPS_ANNUAL_HI) BCLIP=1
            AFAC=MAX(DPS_ANNUAL_LO,MIN(DPS_ANNUAL_HI,TRIAL))
            TOTAL=DPS_BASE_CHEM_N(ISB)*AFAC
            IF (SUMU>0.D0) THEN
              ORIGSHARE=DPS_U/SUMU
            ELSE
!!            The all-zero case uses the regional 40:20:40 reference.
              ZEROFLAG=1
              ORIGSHARE=DPS_PHASE_REF
            ENDIF
            SHARE=ORIGSHARE
!!          Euclidean projection onto a bounded sum-one simplex.
!!          Keep annual total fixed; simultaneous repair of all shares.
            IF (ANY(SHARE<DPS_PHASE_LO).OR.
     &          ANY(SHARE>DPS_PHASE_HI)) THEN
              TAULO=MINVAL(ORIGSHARE-DPS_PHASE_HI)
              TAUHI=MAXVAL(ORIGSHARE-DPS_PHASE_LO)
              DO ITER=1,80
                TAU=0.5D0*(TAULO+TAUHI)
                SHARE=MAX(DPS_PHASE_LO,
     &            MIN(DPS_PHASE_HI,ORIGSHARE-TAU))
                IF (SUM(SHARE)>1.D0) THEN
                  TAULO=TAU
                ELSE
                  TAUHI=TAU
                ENDIF
              ENDDO
              REPAIR=1
            ENDIF
            IF (ABS(SUM(SHARE)-1.D0)>1.D-12.OR.
     &          ANY(SHARE<DPS_PHASE_LO-1.D-12).OR.
     &          ANY(SHARE>DPS_PHASE_HI+1.D-12))
     &        CALL DPS_FAIL('PHASE_PROJECTION_FAILED')
            NTARGET=TOTAL*SHARE
            DPS_FAC=NTARGET/DPS_BASE_PHASE_N(:,ISB)
            DPS_RHO=1.D0-DPS_FAC
          ENDIF
        ENDIF
        IF (ANY(.NOT.IEEE_IS_FINITE(DPS_Z)).OR.
     &      ANY(.NOT.IEEE_IS_FINITE(DPS_PHI)).OR.
     &      ANY(.NOT.IEEE_IS_FINITE(DPS_RAW)).OR.
     &      ANY(.NOT.IEEE_IS_FINITE(DPS_U)).OR.
     &      ANY(.NOT.IEEE_IS_FINITE(DPS_FAC)).OR.
     &      ANY(.NOT.IEEE_IS_FINITE(NTARGET)).OR.
     &      .NOT.IEEE_IS_FINITE(AFAC))
     &    CALL DPS_FAIL('SUBBASIN_ACTION_NONFINITE')
        IF (ANY(NTARGET<0.D0))
     &    CALL DPS_FAIL('NEGATIVE_PHASE_N')
        DPS_CHEM_SCALE(:,ISB)=REAL(DPS_FAC)
        WRITE(196,'(9(I6,1X),56(ES24.16,1X))')
     &    IYR,ISB,DPS_STATE_YEAR,ACTIVE,DPS_AG_COUNT(ISB),OOB,
     &    BCLIP,ZEROFLAG,REPAIR,DPS_AG_AREA(ISB),1.D0-AFAC,0.D0,
     &    AFAC,1.D0,DPS_BASE_CHEM_N(ISB),SUM(NTARGET),
     &    DPS_BASE_PHASE_N(:,ISB),NTARGET,DPS_FAC,
     &    ORIGSHARE,SHARE,
     &    DPS_STATE(:,ISB),DPS_Z,DPS_RAW,DPS_U,DPS_PHI
      ENDDO
      CLOSE(196)
      END

      SUBROUTINE DPS_CAPTURE_OUTLET
      USE PARM
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      IF (.NOT.DPS_ENABLED) RETURN
      DPS_OUTLET_TN=RCHYRO(7,DPS_OUTLET)+
     &  RCHYRO(16,DPS_OUTLET)+RCHYRO(33,DPS_OUTLET)+
     &  RCHYRO(35,DPS_OUTLET)
      DPS_OUTLET_FLOW=RCHYRO(2,DPS_OUTLET)
      IF (.NOT.IEEE_IS_FINITE(DPS_OUTLET_TN).OR.
     &    .NOT.IEEE_IS_FINITE(DPS_OUTLET_FLOW))
     &  CALL DPS_FAIL('OUTLET_NONFINITE')
      OPEN(195,FILE='dps_outcome.tsv',POSITION='APPEND')
      WRITE(195,'(I4,1X,2(ES24.16,1X))')
     & IYR,DBLE(DPS_OUTLET_TN),DBLE(DPS_OUTLET_FLOW)
      CLOSE(195)
      END

      SUBROUTINE DPS_CAPTURE_YEAR_START
      USE PARM
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER :: JH,LYR,ISB
      REAL(KIND=8) :: WT,HH(DPS_B)
      IF (.NOT.DPS_ENABLED) RETURN
      DPS_HAS_STATE=.FALSE.
      DPS_STATE=0.D0
      DPS_AG_AREA=0.D0
      DPS_AG_COUNT=0
      OPEN(194,FILE='dps_hru_year_start_state.tsv',
     & POSITION='APPEND')
      DO JH=1,NHRU
        IF (.NOT.DPS_AGMASK(JH)) CYCLE
        ISB=HRU_SUB(JH)
        IF (ISB<1.OR.ISB>SUBTOT)
     &    CALL DPS_FAIL('STATE_SUBBASIN_RANGE')
        WT=DBLE(HRU_KM(JH))
        IF (.NOT.IEEE_IS_FINITE(WT))
     &    CALL DPS_FAIL('HRU_AREA_NONFINITE')
        IF (WT<=0.D0) CALL DPS_FAIL('HRU_AREA_INVALID')
        IF (SOL_NLY(JH)<1.OR.SOL_NLY(JH)>SIZE(SOL_NO3,1))
     &    CALL DPS_FAIL('SOIL_LAYER_COUNT_INVALID')
        HH=0.D0
        DO LYR=1,SOL_NLY(JH)
          HH(1)=HH(1)+DBLE(SOL_NO3(LYR,JH))
          HH(2)=HH(2)+DBLE(SOL_NH3(LYR,JH))
          HH(3)=HH(3)+DBLE(SOL_FON(LYR,JH))
          HH(4)=HH(4)+DBLE(SOL_AORGN(LYR,JH))
          HH(5)=HH(5)+DBLE(SOL_ORGN(LYR,JH))
        ENDDO
        HH(6)=DBLE(SHALLST_N(JH))
        HH(8)=DBLE(SHALLST(JH))
        IF (.NOT.IEEE_IS_FINITE(SOL_SUMFC(JH)))
     &    CALL DPS_FAIL('SOIL_FIELD_CAPACITY_NONFINITE')
        IF (SOL_SUMFC(JH)>0.) THEN
          HH(7)=DBLE(SOL_SW(JH))/DBLE(SOL_SUMFC(JH))
        ELSE
          CALL DPS_FAIL('SOIL_FIELD_CAPACITY_INVALID')
        ENDIF
        IF (ANY(.NOT.IEEE_IS_FINITE(HH)))
     &    CALL DPS_FAIL('HRU_STATE_NONFINITE')
        DPS_AG_AREA(ISB)=DPS_AG_AREA(ISB)+WT
        DPS_AG_COUNT(ISB)=DPS_AG_COUNT(ISB)+1
        DPS_STATE(:,ISB)=DPS_STATE(:,ISB)+WT*HH
        WRITE(194,'(3(I6,1X),9(ES24.16,1X))')
     &    IYR,JH,ISB,WT,HH
      ENDDO
      CLOSE(194)
      IF (ANY(.NOT.IEEE_IS_FINITE(DPS_AG_AREA)))
     &  CALL DPS_FAIL('SUBBASIN_AREA_NONFINITE')
      DO ISB=1,SUBTOT
        IF (DPS_AG_COUNT(ISB)>0) THEN
          IF (DPS_AG_AREA(ISB)<=0.D0)
     &      CALL DPS_FAIL('SUBBASIN_AREA_INVALID')
          DPS_STATE(:,ISB)=DPS_STATE(:,ISB)/DPS_AG_AREA(ISB)
        ENDIF
      ENDDO
      IF (ANY(.NOT.IEEE_IS_FINITE(DPS_STATE)))
     &  CALL DPS_FAIL('STATE_NONFINITE')
      DPS_STATE_YEAR=IYR
      DPS_HAS_STATE=.TRUE.
      OPEN(197,FILE='dps_year_start_state.tsv',POSITION='APPEND')
      DO ISB=1,SUBTOT
        WRITE(197,'(3(I6,1X),9(ES24.16,1X))')
     &    IYR,ISB,DPS_AG_COUNT(ISB),DPS_AG_AREA(ISB),
     &    DPS_STATE(:,ISB)
      ENDDO
      CLOSE(197)
      END

      SUBROUTINE DPS_FAIL(MSG)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: MSG
      WRITE(*,*) 'DPS_ERROR: ',MSG
      STOP 41
      END
