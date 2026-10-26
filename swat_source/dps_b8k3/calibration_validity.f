      subroutine calibration_finite_check(vals,nvals,label)
      use parm, only: iyr,i,ihru,inum1
      use, intrinsic :: ieee_arithmetic
      implicit none
      integer, intent(in) :: nvals
      real, intent(in) :: vals(nvals)
      character(*), intent(in) :: label
      if (.not.all(ieee_is_finite(vals))) then
        write(*,*) 'NONFINITE_STATE ',label
        write(*,*) 'YEAR,STEP,HRU,REACH ',iyr,i,ihru,inum1
        write(*,*) vals
        stop 91
      endif
      end
