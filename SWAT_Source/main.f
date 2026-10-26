      program main
!!    this is the main program that reads input, calls the main simulation
!!    model, and writes output.
!!    ~ ~ ~ ~ ~ ~ ~ Calibration parameters  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    Command line: 27 calibration values followed by the output id.
!!    The parameter contract is saved beside this build. Values are
!!    relative or absolute as specified in that contract.
!!    v9 OUTPUTID
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    date        |NA            |date simulation is performed where leftmost
!!                               |eight characters are set to a value of
!!                               |yyyymmdd, where yyyy is the year, mm is the 
!!                               |month and dd is the day
!!    isproj      |none          |special project code:
!!                               |1 test rewind (run simulation twice)
!!    time        |NA            |time simulation is performed where leftmost
!!                               |ten characters are set to a value of
!!                               |hhmmss.sss, where hh is the hour, mm is the 
!!                               |minutes and ss.sss is the seconds and
!!                               |milliseconds
!!    values(1)   |year          |year simulation is performed
!!    values(2)   |month         |month simulation is performed
!!    values(3)   |day           |day in month simulation is performed
!!    values(4)   |minutes       |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    values(5)   |hour          |hour simulation is performed
!!    values(6)   |minutes       |minute simulation is performed
!!    values(7)   |seconds       |second simulation is performed
!!    values(8)   |milliseconds  |millisecond simulation is performed
!!    zone        |NA            |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    prog        |NA            |program name and version
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    i           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: date_and_time
!!    SWAT: getallo, allocate_parms, readfile, readfig
!!    SWAT: readbsn, std1, readwwq, readinpt, std2, storeinitial
!!    SWAT: openwth, headout, simulate, finalbal, writeaa, pestw 
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
      !use license_utils
      implicit none
     
      ! 声明命令行参数变量
      integer :: narg, output_id, iarg
      real :: v1, v2, v3, v4, v5, v6, v7, v8
      character(len=32) :: arg
      !call WELCOME
      !if (.not. check_license()) then
      !    write(*,*) "Error: Invalid license"
      !    write(*,*) "Please contact vendor for a valid license"
      !    stop
      !endif       
      ! 检查命令行参数数量
      narg = command_argument_count()
      if (narg /= 28) then
        write(*,*) "Error: need 27 parameters and one output id"
        write(*,*) "Usage: SWAT664-fixed v1 ... v27 output_id"
        stop 2
      endif

      ! Keep internal slot 5 unused: it was ALPHA_BNK in the
      ! 28-parameter build. External parameters 5-27 map to
      ! internal slots 6-28, preserving all other SWAT injections.
      calp = 0.
      do iarg = 1, 4
        call get_command_argument(iarg, arg)
        read(arg,*) calp(iarg)
      end do
      do iarg = 5, 27
        call get_command_argument(iarg, arg)
        read(arg,*) calp(iarg + 1)
      end do
      v1 = calp(1)
      v2 = calp(2)
      v3 = calp(3)
      v4 = calp(4)
      v5 = 0.
      v6 = calp(6)
      v7 = calp(7)
      v8 = calp(8)
      call get_command_argument(28, arg)
      read(arg,*) output_id

!! BEGIN GENERATED BOUNDS
      if (.not.(calp(1) >= 0.000000000e+00 .and.
     &    calp(1) <= 1.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS CANMX"
        stop 3
      endif
      if (.not.(calp(2) >= -4.500000000e-01 .and.
     &    calp(2) <= 2.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS CN2"
        stop 3
      endif
      if (.not.(calp(3) >= 0.000000000e+00 .and.
     &    calp(3) <= 3.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS CH_N2"
        stop 3
      endif
      if (.not.(calp(4) >= -1.000000000e-02 .and.
     &    calp(4) <= 5.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS CH_K2"
        stop 3
      endif
      if (.not.(calp(6) >= -2.000000000e-01 .and.
     &    calp(6) <= 1.500000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS SOL_AWC"
        stop 3
      endif
      if (.not.(calp(7) >= -2.500000000e-01 .and.
     &    calp(7) <= 3.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS SOL_K"
        stop 3
      endif
      if (.not.(calp(8) >= -1.000000000e-01 .and.
     &    calp(8) <= 1.200000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS SOL_BD"
        stop 3
      endif
      if (.not.(calp(9) >= 0.000000000e+00 .and.
     &    calp(9) <= 1.000000000e+03)) then
        write(*,*) "OUT_OF_BOUNDS CH_W2"
        stop 3
      endif
      if (.not.(calp(10) >= -5.000000000e-02 .and.
     &    calp(10) <= 5.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS CH_L2"
        stop 3
      endif
      if (.not.(calp(11) >= 0.000000000e+00 .and.
     &    calp(11) <= 1.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS ALPHA_BF"
        stop 3
      endif
      if (.not.(calp(12) >= 0.000000000e+00 .and.
     &    calp(12) <= 1.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS ESCO"
        stop 3
      endif
      if (.not.(calp(13) >= 0.000000000e+00 .and.
     &    calp(13) <= 5.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS GW_DELAY"
        stop 3
      endif
      if (.not.(calp(14) >= 0.000000000e+00 .and.
     &    calp(14) <= 3.000000000e+01)) then
        write(*,*) "OUT_OF_BOUNDS CH_D"
        stop 3
      endif
      if (.not.(calp(15) >= -1.000000000e+01 .and.
     &    calp(15) <= 1.000000000e+01)) then
        write(*,*) "OUT_OF_BOUNDS TLAPS"
        stop 3
      endif
      if (.not.(calp(16) >= 5.000000000e-02 .and.
     &    calp(16) <= 9.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS TIMP"
        stop 3
      endif
      if (.not.(calp(17) >= 0.000000000e+00 .and.
     &    calp(17) <= 2.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS GWQMN"
        stop 3
      endif
      if (.not.(calp(18) >= -5.000000000e+00 .and.
     &    calp(18) <= 5.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS SMFMX"
        stop 3
      endif
      if (.not.(calp(19) >= 0.000000000e+00 .and.
     &    calp(19) <= 2.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS GW_REVAP"
        stop 3
      endif
      if (.not.(calp(20) >= 0.000000000e+00 .and.
     &    calp(20) <= 5.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS REVAPMN"
        stop 3
      endif
      if (.not.(calp(21) >= -1.000000000e-03 .and.
     &    calp(21) <= 1.000000000e+01)) then
        write(*,*) "OUT_OF_BOUNDS CH_S2"
        stop 3
      endif
      if (.not.(calp(22) >= -5.000000000e+00 .and.
     &    calp(22) <= 5.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS SMFMN"
        stop 3
      endif
      if (.not.(calp(23) >= 0.000000000e+00 .and.
     &    calp(23) <= 9.000000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS SNO50COV"
        stop 3
      endif
      if (.not.(calp(24) >= 0.000000000e+00 .and.
     &    calp(24) <= 5.000000000e+02)) then
        write(*,*) "OUT_OF_BOUNDS SNOCOVMX"
        stop 3
      endif
      if (.not.(calp(25) >= 0.000000000e+00 .and.
     &    calp(25) <= 2.500000000e-01)) then
        write(*,*) "OUT_OF_BOUNDS SOL_ALB"
        stop 3
      endif
      if (.not.(calp(26) >= 1.000000000e-01 .and.
     &    calp(26) <= 2.000000000e+01)) then
        write(*,*) "OUT_OF_BOUNDS SURLAG"
        stop 3
      endif
      if (.not.(calp(27) >= -1.000000000e+01 .and.
     &    calp(27) <= 5.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS SFTMP"
        stop 3
      endif
      if (.not.(calp(28) >= 0.000000000e+00 .and.
     &    calp(28) <= 1.000000000e+00)) then
        write(*,*) "OUT_OF_BOUNDS EPCO"
        stop 3
      endif
!! END GENERATED BOUNDS















































































      prog = "SWAT Dec 23 2016    VER 2016/Rev 664"
      write (*,1000)
 1000 format(1x,"    SWATBO for Calibration purpose  ",/,             
     &          "               Rev. 2024              ",/,             
     &          "      Soil & Water Assessment Tool    ",/,             
     &          "               WINDOWS Version             ",/,             
     &          " Program reading from file.cio . . . executing",/)
!! process input
		
      call getallo
      call allocate_parms
      call readfile(output_id)
      call readbsn
      call readwwq
      if (fcstyr > 0 .and. fcstday > 0) call readfcst
      call readplant             !! read in the landuse/landcover database
      call readtill              !! read in the tillage database
      call readpest              !! read in the pesticide database
      call readfert              !! read in the fertilizer/nutrient database
      call readurban             !! read in the urban land types database
      call readseptwq            !! read in the septic types database
      call readlup
      call readfig(v1,v2,v3,v4,v5,v6,v7,v8,output_id)
      call readatmodep
      call readinpt
      call calibration_audit(output_id)
      call std1
      call std2
      call openwth
      call headout

      !! convert integer to string for output.mgt file
      subnum = ""
      hruno = ""
      do i = 1, mhru
        write (subnum(i),fmt=' (i5.5)') hru_sub(i)
        write (hruno(i),fmt=' (i4.4)') hru_seq(i)  
      end do

      if (isproj == 2) then 
        hi_targ = 0.0
      end if

!! save initial values
      if (isproj == 1) then
        scenario = 2
        call storeinitial
      else if (fcstcycles > 1) then
        scenario =  fcstcycles
        call storeinitial
      else
        scenario = 1
      endif
        if (iclb /= 4) then
      do iscen = 1, scenario

     
        !! simulate watershed processes
        call simulate

        !! perform summary calculations
        call finalbal
        call writeaa
        call pestw

        !!reinitialize for new scenario
        if (scenario > iscen) call rewind_init
      end do
         end if
      do i = 101, 109       !Claire 12/2/09: change 1, 9  to 101, 109.
        close (i)
      end do
      close(124)
      write (*,1001)
 1001 format (/," Execution successfully completed ")
	
        iscen=1
!! file for Mike White to review to ensure simulation executed normally
      open (9999,file='fin.fin')
      write (9999,*) 'Execution successful'
      close (9999)
      
	stop
      end
