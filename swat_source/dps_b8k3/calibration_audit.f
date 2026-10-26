      subroutine calibration_audit(output_id)
      use parm
      implicit none
      integer, intent(in) :: output_id
      integer :: j
      character(len=80) :: aname
      write(aname,'(i0,a)') output_id, 'parameter_audit.txt'
      open(198,file=trim(aname),status='replace')
      write(198,*) 'calp_internal', calp
      write(198,*) 'HRU1_CANMX', canmx(1)
      write(198,*) 'HRU1_CN2', cn2(1)
      write(198,*) 'HRU1_GW', delay(1),alpha_bf(1),gwqmn(1),
     & gw_revap(1),revapmn(1)
      write(198,*) 'HRU1_ESCO', esco(1)
      write(198,*) 'HRU1_EPCO', epco(1)
      write(198,*) 'HRU1_SURLAG', surlag(1)
      write(198,*) 'HRU1_AWC', sol_awc(:,1)
      write(198,*) 'HRU1_BD', sol_bd(:,1)
      write(198,*) 'HRU1_K', sol_k(:,1)
      write(198,*) 'HRU1_ALB', sol_alb(1)
      write(198,*) 'RCH1_W_L_S_N_K', ch_w(2,1), ch_l2(1),
     & ch_s(2,1), ch_n(2,1), ch_k(2,1)
      write(198,*) 'BASIN_SURLAG', surlag_bsn
      write(198,*) 'BASIN_SNOW',sftmp,smfmx,smfmn,timp,
     & snocovmx,sno50cov
      write(198,*) 'SUB1_TLAPS',tlaps(1)
      do j=1,21
        write(198,*) 'REACH',j,ch_w(2,j),ch_d(j),ch_s(2,j),
     &  ch_l2(j),ch_n(2,j),ch_k(2,j),alpha_bnk(j)
      end do
      close(198)
      end
