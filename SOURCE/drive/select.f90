!***********************************************************************
!
      subroutine select(iu,xin,out)
!
! ----------------------------------------------------------------------
!
!     SELECT : Call type subroutines specified by ITYPE.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated to Fortran 90.  June 6, 2007 Cheol Park, NIST
!
!***********************************************************************

    use modsim_head
    implicit none

    integer                   :: iu,i,ivar,itype,nn,nnn
    integer,dimension(minoiu) :: iostat
    real,dimension(minoiu)    :: xin,out

    do i=1,nin(iu)
      ivar=in(iu,i)
      iostat(i)=icheck(ivar)
    enddo

    itype=iunits(iu)
    nn=jpar(iu)
    nnn=isaved(iu)

    select case(itype)

! -----------------------------------------------------------------------------

!  Fan or pump

      case(1)
        call type1(xin,out,par(nn),saved(nnn),iostat)

!  Conduit (duct or pipe)

      case(2)
        call type2(xin,out,par(nn),saved(nnn),iostat)

!  Inlet conduit (duct or pipe)

      case(3)
        call type3(xin,out,par(nn),saved(nnn),iostat)

!  Flow merge

      case(4)
        call type4(xin,out,par(nn),saved(nnn),iostat)

!  Damper or valve

      case(5)   
        call type5(xin,out,par(nn),saved(nnn),iostat)

!  Flow split

      case(6)   
        call type6(xin,out,par(nn),saved(nnn),iostat)

!  Temperature sensor

      case(7)   
        call type7(xin,out,par(nn),saved(nnn),iostat)

!  Proportional-integral ontroller

      case(8)   
        call type8(xin,out,par(nn),saved(nnn),iostat)

!  Linear valve with pneumatic actuator

      case(9)   
        call type9(xin,out,par(nn),saved(nnn),iostat)

!  Hot water to air heating coil -- simple model

      case(10)   
        call type10(xin,out,par(nn),saved(nnn),iostat)

!  Hot water to air heating coil -- detailed model

      case(11)   
        call type11(xin,out,par(nn),saved(nnn),iostat)

!  Cooling or dehumidifying coil (nr)

      case(12)   
        call type12(xin,out,par(nn),saved(nnn),iostat)

!  Three-way valve with actuator

      case(13)   
        call type13(xin,out,par(nn),saved(nnn),iostat)

!  Evaporative humidifier

      case(14)   
        call type14(xin,out,par(nn),saved(nnn),iostat)

!  Room with constant loads

      case(15)   
        call type15(xin,out,par(nn),saved(nnn),iostat)

!  'Sticky' proportional controller

      case(16)   
        call type16(xin,out,par(nn),saved(nnn),iostat)

!  Mixing dampers and merge

      case(17)   
        call type17(xin,out,par(nn),saved(nnn),iostat)

!  Plenum

      case(18)   
        call type18(xin,out,par(nn),saved(nnn),iostat)

!  Flow balance control

      case(19)   
        call type19(xin,out,par(nn),saved(nnn),iostat)

!  High or low limit controller

      case(20)   
        call type20(xin,out,par(nn),saved(nnn),iostat)

!  'Grounded' flow split

      case(21)   
        call type21(xin,out,par(nn),saved(nnn),iostat)

!  Steam spray humidifier

      case(22)   
        call type22(xin,out,par(nn),saved(nnn),iostat)

!  Steam nozzle

      case(23)   
        call type23(xin,out,par(nn),saved(nnn),iostat)

!  Ideal gas nozzle

      case(24)   
        call type24(xin,out,par(nn),saved(nnn),iostat)

!  Steam to air heat exchanger

      case(25)   
        call type25(xin,out,par(nn),saved(nnn),iostat)

!  Control signal inverter

      case(26)   
        call type26(xin,out,par(nn),saved(nnn),iostat)

!  Moist air flow merge

      case(27)   
        call type27(xin,out,par(nn),saved(nnn),iostat)

!  Constant flow resistance

      case(28)   
        call type28(xin,out,par(nn),saved(nnn),iostat)

!  Inlet onstant flow resistance

      case(29)   
        call type29(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Cooling or dehumidifying oil (#2)

      case(30)   
        call type30(xin,out,par(nn),saved(nnn),iostat)

!  Moist air mixing dampers and merge

      case(33)   
        call type33(xin,out,par(nn),saved(nnn),iostat)

!  Multiplier

      case(34)   
        call type34(xin,out,par(nn),saved(nnn),iostat)

!  Mean values of temperatures and humidity ratios

      case(35)   
        call type35(xin,out,par(nn),saved(nnn),iostat)

!  Summer mode of supply air temperature reset schedule

      case(36)   
        call type36(xin,out,par(nn),saved(nnn),iostat)

!  Time of day controller with zone demand reset

      case(39)   
        call type39(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Zone envelope

      case(50)   
        call type50(xin,out,par(nn),saved(nnn),iostat)

!  Building surface temperatures

      case(51)   
        call type51(xin,out,par(nn),saved(nnn),iostat)

!  Zone air temperature and humidity ratio

      case(52)   
        call type52(xin,out,par(nn),saved(nnn),iostat)

!  Weather input

      case(53)   
        call type53(xin,out,par(nn),saved(nnn),iostat)

!  Zone envelope (repeat)

      case(54)   
        call type54(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Hot water boiler with a domestic hot water heating coil

      case(62)
        call type62(xin,out,par(nn),saved(nnn),iostat)

!  Hot water coil with constant wall temperature

      case(63)
        call type63(xin,out,par(nn),saved(nnn),iostat)

!  Boiler burner and circulating pump control

      case(64)
        call type64(xin,out,par(nn),saved(nnn),iostat)

!  data storage allocation for input for HVAC-Cx AHU
        
      case(109)
        call type109(xin,out,par(nn),saved(nnn),iostat)
! -----------------------------------------------------------------------------

!  Static boiler (IEA Annex 10)

      case(122)
        call type122(xin,out,par(nn),saved(nnn),iostat)

!  Chiller (IEA Annex 17)

      case(124)
        call type124(xin,out,par(nn),saved(nnn),iostat)

!  Cooling tower (IEA Annex 17)

      case(143)
        call type143(xin,out,par(nn),saved(nnn),iostat)

!  Chiller sump

      case(144)
        call type144(xin,out,par(nn),saved(nnn),iostat)

!  Cooling tower (IEA Annex 17) - diffeq

      case(145)
        call type145(xin,out,par(nn),saved(nnn),iostat)

!  Chiller sump - diffeq

      case(146)
        call type146(xin,out,par(nn),saved(nnn),iostat)

!  Cooling tower controller (IEA Annex 17)

      case(179)
        call type179(xin,out,par(nn),saved(nnn),iostat)

!  Mixing of water flows

      case(200)
        call type200(xin,out,par(nn),saved(nnn),iostat)

!  Adding electric powers

      case(201)
        call type201(xin,out,par(nn),saved(nnn),iostat)

!  Integrator of thermal power for energy

      case(202)
        call type202(xin,out,par(nn),saved(nnn),iostat)

!  Relative humidity and wet-bulb temprature

      case(203)
        call type203(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Temperature sensor

      case(301)   
        call type301(xin,out,par(nn),saved(nnn),iostat)

!  Humidity sensor

      case(302)   
        call type302(xin,out,par(nn),saved(nnn),iostat)

!  Flow rate sensor

      case(303)   
        call type303(xin,out,par(nn),saved(nnn),iostat)

!  Total pressure sensor

      case(304)   
        call type304(xin,out,par(nn),saved(nnn),iostat)

!  Static pressure sensor

      case(305)   
        call type305(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Motor-driven actuator

      case(321)   
        call type321(xin,out,par(nn),saved(nnn),iostat)

!  Damper

      case(322)   
        call type322(xin,out,par(nn),saved(nnn),iostat)

!  Damper - calculates flow rate

      case(323)   
        call type323(xin,out,par(nn),saved(nnn),iostat)

!  Mixing box, mixed flow as input

      case(324)   
        call type324(xin,out,par(nn),saved(nnn),iostat)

!  Mixing box, calculates mixed flow

      case(325)   
        call type325(xin,out,par(nn),saved(nnn),iostat)

!  Mixing box, min oa damper, mixed flow as input

      case(326)   
        call type326(xin,out,par(nn),saved(nnn),iostat)

!  Mixing box with minimum outside air damper

      case(327)   
        call type327(xin,out,par(nn),saved(nnn),iostat)

!  Two port control valve

      case(328)   
        call type328(xin,out,par(nn),saved(nnn),iostat)

!  Two port valve - calculates flow

      case(329)   
        call type329(xin,out,par(nn),saved(nnn),iostat)

!  Three port mixing valve

      case(330)   
        call type330(xin,out,par(nn),saved(nnn),iostat)

!  Variable speed drive

      case(333)   
        call type333(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Fluid resistance

      case(341)   
        call type341(xin,out,par(nn),saved(nnn),iostat)

!  Fluid resistance - calculates flow rate

      case(342)   
        call type342(xin,out,par(nn),saved(nnn),iostat)

!  Flow split

      case(345)   
        call type345(xin,out,par(nn),saved(nnn),iostat)

!  Asymmetric flow split

      case(346)   
        call type346(xin,out,par(nn),saved(nnn),iostat)

!  Flow merge

      case(348)   
        call type348(xin,out,par(nn),saved(nnn),iostat)

!  Room air mass balance

      case(349)   
        call type349(xin,out,par(nn),saved(nnn),iostat)

!  Fan or pump

      case(350)   
        call type350(xin,out,par(nn),saved(nnn),iostat)

!  Fan or pump (implicit flow)

      case(351)   
        call type351(xin,out,par(nn),saved(nnn),iostat)

!  Fan or pump - temperature rise

      case(352)   
        call type352(xin,out,par(nn),saved(nnn),iostat)

!  Fan or pump (implicit flow) - temperature rise

      case(353)   
        call type353(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Heating/cooling coil

      case(362)   
        call type362(xin,out,par(nn),saved(nnn),iostat)

!  Ideal heating of fluid stream

      case(366)   
        call type366(xin,out,par(nn),saved(nnn),iostat)

!  Mixing of two moist air streams

      case(367)   
        call type367(xin,out,par(nn),saved(nnn),iostat)

!  Mixing of six moist air streams

      case(368)   
        call type368(xin,out,par(nn),saved(nnn),iostat)

!  Supply and return flow rates

      case(369)   
        call type369(xin,out,par(nn),saved(nnn),iostat)

!  Modify pressure value MAG 240915

      case(391)   
        call type391(xin,out,par(nn),saved(nnn),iostat)

!  Modify flow value MAG 240915

      case(392)   
        call type392(xin,out,par(nn),saved(nnn),iostat)

!  Modify temperature value MAG 240915

      case(393)   
        call type393(xin,out,par(nn),saved(nnn),iostat)

!  Modify control value MAG 240915

      case(394)   
        call type394(xin,out,par(nn),saved(nnn),iostat)

!  Modify other value MAG 240915

      case(395)   
        call type395(xin,out,par(nn),saved(nnn),iostat)

!  Modify energy value MAG 240915

      case(396)   
        call type396(xin,out,par(nn),saved(nnn),iostat)

!  Modify power value MAG 240915

      case(397)   
        call type397(xin,out,par(nn),saved(nnn),iostat)

!  Modify humidity value MAG 240915

      case(398)   
        call type398(xin,out,par(nn),saved(nnn),iostat)

!  Modify variable value by multiplier MAG 240915

      case(399)   
        call type399(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Room with interzone flows (no plenum)

      case(401)   
        call type401(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum return and interzone flows

      case(402)   
        call type402(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum and ducted return and interzone flows

      case(403)   
        call type403(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum and ducted return and interzone flows (explicit integ)

      case(404)   
        call type404(xin,out,par(nn),saved(nnn),iostat)

!  Room (no plenum, no interzone flows)

      case(411)   
        call type411(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum return (no interzone flows)

      case(412)   
        call type412(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum and ducted return (no interzone flows)

      case(413)   
        call type413(xin,out,par(nn),saved(nnn),iostat)

!  Room with plenum and ducted return (no interzone flows) explicit integ

      case(414)   
        call type414(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  PID controller

      case(441)   
        call type441(xin,out,par(nn),saved(nnn),iostat)

!  Fanger pmv and ppd

      case(461)   
        call type461(xin,out,par(nn),saved(nnn),iostat)

!  Heat meter

      case(462)   
        call type462(xin,out,par(nn),saved(nnn),iostat)

!  Energy meter

      case(463)   
        call type463(xin,out,par(nn),saved(nnn),iostat)

!  E-51 supply fan static pressure control

      case(481)   
        call type481(xin,out,par(nn),saved(nnn),iostat)

!  E-51 flow difference control of return fan

      case(482)   
        call type482(xin,out,par(nn),saved(nnn),iostat)

!  E-51 return fan reset control

      case(483)   
        call type483(xin,out,par(nn),saved(nnn),iostat)

!  E-51 minimun outside air damper control

      case(484)   
        call type484(xin,out,par(nn),saved(nnn),iostat)

!  E-51 modulated mixed air damper control

      case(485)   
        call type485(xin,out,par(nn),saved(nnn),iostat)

!  E-51 supply air temperature control

      case(486)   
        call type486(xin,out,par(nn),saved(nnn),iostat)

!  E-51 economizer control

      case(487)   
        call type487(xin,out,par(nn),saved(nnn),iostat)

!  E-51 low temperature control

      case(488)   
        call type488(xin,out,par(nn),saved(nnn),iostat)

!  E-51 supply air temperature reset

      case(489)   
        call type489(xin,out,par(nn),saved(nnn),iostat)

!  VAV room temperature control with reheat

      case(490)   
        call type490(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!  Heating/cooling coil + two port valve

      case(521)   
        call type521(xin,out,par(nn),saved(nnn),iostat)

!  Coil + two port valve, calculates flow

      case(522)   
        call type522(xin,out,par(nn),saved(nnn),iostat)

!  Heating/cooling coil + three port valve

      case(523)   
        call type523(xin,out,par(nn),saved(nnn),iostat)

!  Coil + three port valve, calculates flow

      case(524)   
        call type524(xin,out,par(nn),saved(nnn),iostat)

!  Motorized pressure-independent vav box

      case(525)   
        call type525(xin,out,par(nn),saved(nnn),iostat)

!  Pressure-independent vav box (calculates flow)

      case(526)   
        call type526(xin,out,par(nn),saved(nnn),iostat)

!  Pressure-independent vav box (implicit flow)

      case(527)   
        call type527(xin,out,par(nn),saved(nnn),iostat)

!  Flow split and pressure-independent vav box

      case(528)   
        call type528(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

!
!  Holder for data storage allocation

      case(107)
        call type107(xin,out,par(nn),saved(nnn),iostat)

!  Holder for data storage allocation

      case(108)
        call type108(xin,out,par(nn),saved(nnn),iostat)

!  Holder for data storage allocation

      case(109)
        call type109(xin,out,par(nn),saved(nnn),iostat)

!  Mixing of three moist air streams

      case(365)
        call type365(xin,out,par(nn),saved(nnn),iostat) 

!  Mixing box

      case(371)
        call type371(xin,out,par(nn),saved(nnn),iostat)

!  Room

      case(428)
        call type428(xin,out,par(nn),saved(nnn),iostat)

!  Supply fan

      case(471)
        call type471(xin,out,par(nn),saved(nnn),iostat)

!  Return fan

      case(480)
        call type480(xin,out,par(nn),saved(nnn),iostat)
!  Supply air temperature control

      case(492)
        call type492(xin,out,par(nn),saved(nnn),iostat)

!  Supply air temperature reset control

      case(493)
        call type493(xin,out,par(nn),saved(nnn),iostat)

!  Economizer control

      case(495)
        call type495(xin,out,par(nn),saved(nnn),iostat)

!  VAV room temperature control with reheat

      case(496)
        call type496(xin,out,par(nn),saved(nnn),iostat)

!  Modulated mixed air damper control w/occupancy

      case(497)   
        call type497(xin,out,par(nn),saved(nnn),iostat)

!  Air damper control

      case(499)
        call type499(xin,out,par(nn),saved(nnn),iostat)

!  Read inputs from a file

      case(504)   
        call type504(xin,out,par(nn),saved(nnn),iostat)

!  Dynamic or steady state coil and three port valve

      case(530)
        call type530(xin,out,par(nn),saved(nnn),iostat)

! -----------------------------------------------------------------------------

      case default

!  No such type (probably commented out)

        print *,' Error in subroutine select: call to type',itype,&
                ' not found.'
        stop
    end select

!  Set variable status (icheck) from unit output status (iostat).

    do i=1,nout(iu)
      ivar=iout(iu,i)
      if(ivar>0) then
        if(i<=nde(iu)) then
          idechk(inde(iu,i))=iostat(i)
        endif
        if((icheck(ivar)>=-2).and.(icheck(ivar)<=1)) then
          if(icheck(ivar)>=0) then
            icheck(ivar)=iostat(i)
          elseif(iostat(i)==0) then
            icheck(ivar)=-1
          endif
        endif
      endif
    enddo
   
    !if((iu >= 40) .and. (iu <= 45)) then
    !   	  write(44, fmt="(f7.0)", advance="no") time
	!      write(44, fmt="(a)", advance="no") ", "
    !   	  write(44, fmt="(i3)", advance="no") iu
	!      write(44, fmt="(a)", advance="no") ", "
	!      write(44, '(9999(e16.8, a))') (state(i), ",", i=1,nstate)
    !endif

    return
    end subroutine select

