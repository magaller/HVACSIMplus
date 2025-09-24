! ******************************************************************************
!
!     UPD_INFO:   Update the model information file that is an output file
!                 of HVACGEN with option VIEW ALL
!
!     Note that the UNIT numbers in the supplemental information file
!     must be located at columns 1 - 3 in right-justified format.
!
!     Created:    Feb. 6, 1995   Cheol Park
!     Updated:    May 25, 2007
!                 converted to Fortran 90
!
! ******************************************************************************

      program upd_info
      implicit none

      character(len=80)      :: data_char, unit_info
      character(len=50)      :: unit_desc
      character(len=30)      :: in_file_1, in_file_2, out_file
      integer                :: ios, i

      print *,' Enter the name of view output file of HVACGEN (viewsave.txt)'
      write(*,fmt='(a4)',advance='no') ' => '
      read (*, fmt='(a30)') in_file_1

      print *,' Enter the supplemental unit information file (hvacsim.inf)'
      write(*,fmt='(a4)',advance='no') ' => '
      read (*, fmt='(a30)') in_file_2

      print *,' Enter the name of model information file (hvacsim.model)'
      write(*,fmt='(a4)',advance='no') ' => '
      read (*, fmt='(a30)') out_file

      open (unit=10, file=in_file_1, status='old')
      open (unit=11, file=in_file_2, status='old')
      open (unit=12, file=out_file, status='unknown')

      i = 0
      do
         read (10, fmt='(a80)', iostat=ios) data_char
         if(ios < 0) then
            print *,'End of file'
            stop
         endif

         i = i +1
         if (data_char(2:5) == 'UNIT') then                                   
             print *, 'Line No. ', i
             rewind 11
             20  read (11, fmt='(a50)') unit_desc
             if (unit_desc(1:3) == data_char(7:9)) then                       
                print *, '----- matched data ----'                            
                unit_info = data_char(1:23)//' ----- '//unit_desc(5:50)       
             else                                                             
                goto 20                                                       
             endif                                                            
             write (12, fmt='(a80)') unit_info                                
         else                                                                 
             write (12, fmt='(a80)') data_char                                
         endif                                                                
      enddo                                                                
      end program upd_info

