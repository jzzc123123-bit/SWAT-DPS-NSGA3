! Calls the actual published DPS_YEAR_START and DPS_RBF_EVAL routines.
! This is an arithmetic/interface check, not a hydrological simulation.
program policy_oracle
  use parm
  implicit none
  integer :: ncase, ii, ib
  subtot=1
  dps_enabled=.true.
  dps_has_state=.true.
  dps_mode=1
  dps_start_year=1
  dps_annual_lo=.2d0
  dps_annual_hi=1.2d0
  dps_phase_lo=[.3d0,.1d0,.3d0]
  dps_phase_hi=[.5d0,.3d0,.5d0]
  dps_phase_ref=[.4d0,.2d0,.4d0]
  allocate(dps_state(dps_b,1),dps_ag_area(1),dps_ag_count(1))
  allocate(dps_chem_scale(dps_k,1),dps_man_scale(1))
  allocate(dps_base_chem_n(1),dps_base_phase_n(dps_k,1))
  dps_ag_area=1.d0
  dps_ag_count=2
  open(10,file='oracle_input.dat',status='old')
  read(10,*) ncase
  open(196,file='dps_actions.tsv',status='replace')
  close(196)
  do ii=1,ncase
    read(10,*) dps_lo
    read(10,*) dps_hi
    read(10,*) dps_state(:,1)
    read(10,*) dps_base_phase_n(:,1)
    read(10,*) dps_a
    do ib=1,dps_n
      read(10,*) dps_w(ib,:)
    end do
    do ib=1,dps_n
      read(10,*) dps_c(ib,:)
    end do
    do ib=1,dps_n
      read(10,*) dps_r(ib,:)
    end do
    dps_base_chem_n(1)=sum(dps_base_phase_n(:,1))
    iyr=2000+ii
    dps_state_year=iyr
    call dps_year_start
  end do
  close(10)
end program policy_oracle
