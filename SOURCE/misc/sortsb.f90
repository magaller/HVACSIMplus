! **********************************************************************
!
!   sortsb: Sort the output data for each superblock in order to
!           use with a plotting routine
!
!           The first line lists the row and column numbers of
!           the selected superblock data set.  Output can be made
!           with skipping a specified number of rows, nskip.
!
! -------------------------------------------------------------------
!
!   created :    June 25,1985   Cheol Park
!   updated to Fortran 90:    Jan. 3, 2007
!
! **********************************************************************

    program sortsb
    implicit none

    integer,parameter                  :: maxcha=60*15       !  4/6/2008
    character(len=1)                   :: ans
    character(len=2)                   :: number
    character(len=4)                   :: label
    character(len=46)                  :: infile,outfile
    character(len=80)                  :: indata
    character(len=1),dimension(maxcha) :: charx
    character(len=80),dimension(8)     :: showf
    character(len=2),dimension(10)     :: nsblk=(/' 1',' 2',' 3',' 4',' 5',&
                                                  ' 6',' 7',' 8',' 9','10'/)
    integer                            :: isblk,nskip,icount,i,ndata,ncol,io,&
                                          irow,nrow,nrow_n,nrow_new,kk,iskip
    integer                            :: ifile1=7,ifile2=8,ifile3=9
    real                               :: scale,time,times
    real,dimension(50)                 :: data

    print *,' Enter input file name '
    write(*,fmt='(a4)',advance='no') ' => '
    read (*,fmt='(a46)') infile
10  print *,' Enter output file name '
    write(*,fmt='(a4)',advance='no') ' => '
    read (*,fmt='(a46)') outfile
    open(ifile1,file=infile, status='old')
    open(ifile2,file=outfile,status='replace')
    open(ifile3, status='scratch')

    print *,' Superblock # ?'
    write(*,fmt='(a4)',advance='no') ' => '
    read *, isblk
    print *, ' Number of output lines to be skipped =  '
    write(*,fmt='(a4)',advance='no') ' => '
    read *, nskip

20  read(ifile1,500) indata
    500 format(a80)
    if(indata(1:4)=='SUPE'.and.indata(11:12)==nsblk(isblk)) then
      icount=0
      read (ifile1,fmt='(80a1)') charx
      do i=1,maxcha
        if(charx(i)=='.') then
          icount=icount+1
        elseif(charx(i)=='S') then
          ndata=icount
          exit
        endif
      enddo
      rewind ifile1
    else
      goto 20
    endif

    read(ifile1,500) showf
    write(*,fmt='(1x,a80)') showf

    do i=1,8
      backspace ifile1
    enddo

    print *,' Number of seconds per unit time?'
    write(*,fmt='(a4)',advance='no') ' => '
    read *, scale
    ncol = ndata + 1

    irow = 0
40  read(ifile1,500,iostat=io) indata
    if(io == 0) then
      if(indata(1:4)=='SUPE'.and.indata(11:12)==nsblk(isblk)) then
        backspace ifile1                                              
        read(ifile1,fmt='(a4,6x,a2,f12.2)') label,number,time         
        read(ifile1,fmt='(5g15.6)') (data(i),i=1,ndata)               
        times=time/scale                                              
        write(ifile3,800) times,(data(i),i=1,ndata)                   
        800 format(50g14.6)                                           
        irow = irow + 1                                               
      endif                                                           
      goto 40
    endif
    nrow = irow

    rewind ifile3
    nrow_new = nrow / (nskip + 1)
    write(ifile2,fmt='(2i10)') nrow_new, ncol ! rows and columns of data set
    kk = 0
    iskip = 0
    do while ( kk < nrow )
       read(ifile3,800)  times,(data(i),i=1,ndata)
       if ( iskip == nskip ) then
         write(ifile2,800) times,(data(i),i=1,ndata)
         iskip = 0
       else
         iskip = iskip + 1
       endif
       kk = kk + 1
    end do
    close (ifile2)

    print *,' Extract another superblock? <n>  '
    write(*,fmt='(a4)',advance='no') ' => '
    read(*,fmt='(a1)') ans
    if((ans=='y').or.(ans=='Y')) then
      rewind ifile1
      rewind ifile3
      go to 10
    endif

    stop  '---- End of sortsb -----'
    end program sortsb

