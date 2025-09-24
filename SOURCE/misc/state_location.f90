! ****************************************************************************
!
!        State_location.f90
!
!        To create a reference table of state variables for variable indices
!
!        Created:  Nov. 27, 1998
!        Changed:  Jan. 19, 2001
!        Modified: March 28, 2007
!        Cheol Park, NIST
!
! ***************************************************************************

      program state_location
      implicit none
      integer, parameter           :: max_var = 600
      integer, dimension(max_var)  :: state, pressure, flow, temp,&
                                      control, other, energy, power,&
                                      humidity
      integer  :: nstate,npressure,nflow,ntemp,ncontrol,nother,&
                  nenergy,npower,nhumidity,i

      namelist /nam_input/ nstate,npressure,nflow,ntemp,ncontrol,nother,&
                           nenergy,npower,nhumidity

      open(unit=3,file='state_location.inp')
      open(unit=7,file='state_location.out')

      read(3,nam_input)
      do i=1,npressure
         write(7,*)'state(',i,')  ', '  p(',i,')  '
         print *,  'state(',i,')  ', '  p(',i,')  '
      end do

      do i=1,nflow
         write(7,*) 'state(',i+npressure,')  ', '  m(',i,')'
         print *,   'state(',i+npressure,')  ', '  m(',i,')'
      end do
      do i=1,ntemp
          write(7,*) 'state(',i+npressure+nflow,')  ', '  t(',i,')'
          print *,   'state(',i+npressure+nflow,')  ', '  t(',i,')'
      end do
      do i=1,ncontrol
          write(7,*) 'state(',i+npressure+nflow+ntemp,')  ','  c(',i,')'
          print *,   'state(',i+npressure+nflow+ntemp,')  ','  c(',i,')'
      end do
      do i=1,nother
          write(7,*) 'state(',i+npressure+nflow+ntemp+ncontrol,')  ',&
     		      '  o(',i,')'
          print *,   'state(',i+npressure+nflow+ntemp+ncontrol,')  ',&
     		      '  o(',i,')'
      end do
      do i=1,nenergy
          write(7,*)&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother,&
           ')  ',     '  e(',i,')'
          print *,&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother,&
           ')  ',     '  e(',i,')'
      end do
      do i=1,npower
          write(7,*)&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother+nenergy&
           ,')  ','  q(',i,')'
          print *,&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother+nenergy&
           ,')  ','  q(',i,')'
      end do
      do i=1,nhumidity
          write(7,*)&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother+nenergy+&
            npower,')  ','  h(',i,')'
          print *,&
           'state(',i+npressure+nflow+ntemp+ncontrol+nother+nenergy+&
            npower,')  ','  h(',i,')'
      end do
      stop
      end program state_location

! ***************************************************************************
! namelist input file: state_location.inp
!
!&nam_input nstate=32, npressure=3, nflow=1, ntemp=12, ncontrol=5,
!           nother=0,  nenergy=0,  npower=8, nhumidity=3/
! ***************************************************************************


















