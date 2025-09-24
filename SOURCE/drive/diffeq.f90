!***********************************************************************
!
       subroutine diffeq(dummy,aa,bb,ti,tf,tbar)
!
!***********************************************************************
! * SUBROUTINE:  Analytical integration
! *
! * LANGUAGE:    FORTRAN 77
! *
! * PURPOSE:     Integrate dT/dt = A*T + B over the current time-step
! *              Adapted for HVACSIM+ from the TRNSYS routine DIFFEQ.
! *              See the TRNSYS manual for further details.
! *
! *              The main difference is that there is no need for special
! *              code to deal with the first step time, since initial
! *              values in HVACSIM+ apply to the start of the time-step
! *              and the routine calculates the value of T at the end of
! *              the time-step.
! *
!***********************************************************************
! * inputs
! * ======
! * dummy   : (not need here - retained for compatibility with TRNSYS)
! * aa      : A coefficent, averaged over the time-step
! * bb      : B coefficent, averaged over the time-step
! * ti      : value of T at the beginning of the time-step
! *
! * outputs
! * =======
! * tf      : value of T at the end of the time-step
! * tbar    : value of T averaged over the time-step
! *
!***********************************************************************
!
   use modsim_head
   implicit none
   integer,parameter    :: pd=selected_real_kind(15) ! double precision
   real(kind=pd)        :: tbard,aad,bbd,tid,tfd,tstepd
   real                 :: tbar,tf,ti,aa,bb,dummy

! *** Solution depends on whether the A coefficient is zero

   if(abs(aa) == 0.0) then

! *** A coefficient is zero - linear solution

     tf = bb*tstep + ti
     tbar = (tf + ti)/2.0
   else

! *** A coefficient is non-zero - exponential solution

     aad = aa
     bbd = bb
     tstepd = tstep
     tid = ti
     tfd = (tid + bbd/aad)*dexp(aad*tstepd) - bbd/aad
     tbard = (tid + bbd/aad)/aad/tstepd*(dexp(aad*tstepd)- 1.)&
              - bbd/aad
     tf = real(tfd)
     tbar = real(tbard)
   endif

   return
   end subroutine diffeq

