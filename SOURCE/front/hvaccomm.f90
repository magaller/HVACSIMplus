!=======================================================================

module hvacsim_par

!=======================================================================


!------- The following integer,parameters are used to define the array sizes
!        in hvacgen, slimcon and modsim.

!           maxsbk = maximum superblocks in the simulation plus 1
!           maxblk = maximum blocks in the simulation
!           maxunt = maximum units in the simulation
!           maxstv = maximum state variables in the simulation
!           maxdeq = maximum differential equations in the simulation
!           maxpar = maximum unit integer,parameters in the simulation
!           maxsav = maximum saved variables in the simulation
!           maxbnd = maximum boundary variables in the simulation
!           mdmbnd = dimension of a boundary variable calculation array
!           mblkis = maximum blocks in a superblock
!           muntib = maximum units in a block
!           mdeqis = maximum differential equations in a superblock
!           mdeqiu = maximum differential equations in a unit
!           mseqis = maximum simultaneous equations in a superblock
!           mseqib = maximum simultaneous equations in a block
!           mbndis = maximum boundary conditions in a single superblock
!           minoib = maximum inputs or outputs in a block
!           minoiu = maximum inputs or outputs in a unit
!           mrptis = maximum reported variables in a superblock
!           mpariu = maximum integer,parameters in a unit
!           mdmrep = dimension of a reported variable calculation array
!           maxtyp = maximum number of types
!           nlscre = number of lines on screen - used for scrolling

!           maxlbl = maximum number of labels
!           minois = maximum inputs or outputs in a superblock
!           munitb = maximum units in a block
!           msumary = max(mseqis,mseqib,minoiu,mrptis,maxlbl)

integer,parameter  :: maxsbk=40,maxblk=50,maxunt=400,maxstv=3000,maxdeq=90
integer,parameter  :: maxpar=5000,maxsav=9000,maxbnd=50,mdmbnd=maxbnd*8
integer,parameter  :: mblkis=20,muntib=40
integer,parameter  :: mdeqis=50,mdeqiu=10,mseqis=20,mseqib=75
integer,parameter  :: mbndis=50,minoib=200,minoiu=50,mrptis=60,maxlbl=8
integer,parameter  :: mpariu=30,mdmrep=mrptis*maxsbk
integer,parameter  :: maxtyp=600
integer,parameter  :: nlscre=19
integer,parameter  :: minois=mblkis*minoib
integer,parameter  :: msumary=75

!     For Building Shell Model

!           maxzn  = maximum number of zones in the simulation
!           maxns  = maximum number of surfaces in a zone
!           maxstr = maximum number of constructs
!           maxnrf = maximum number of terms in a conduction transfer function
!           maxord = maximum order of a conduction transfer function

integer,parameter  :: maxzn=10,maxns=12
integer,parameter  :: maxstr=10,maxnrf=20,maxord=5

!     the maximum time (seconds) for writing data on the file, hvacsim.sum,
!     is set by wrtmax.  by changing this value, the size of output can be
!     adjusted.  refer to subroutine report in modino.f90.

real,parameter     :: wrtmax=3600.0

!     typar.dat must be present in current directory

character (len=50) :: typarff = 'typar.dat '

end module hvacsim_par

!===============================================================================

module hvaccomm

!=======================================================================

use hvacsim_par

!---------type declarations---------------------------------------------

character(len=3)  :: extnsn = '  '
character(len=47) :: filnme = &
                     '                                              '
character(len=7),dimension(8) :: iounit = &
           (/'(kPa)  ','(kg/s) ','(C)    ','(-)    ', '(-)    ', &
             '(kJ)   ','(kW)   ','(kg/kg)'/)
character(len=80)                   :: titlsm,dscrib
character(len=80),dimension(minoiu) :: inunit,otunit,imssge,omssge
character(len=80),dimension(mpariu) :: pmssge

logical :: svmode,change
real    :: rtolx,atolx,xtol,ttime

integer :: nnin,nnout,nnpar,unit,icheck
integer :: irdflg,blkflg,block,icnt
integer :: supflg,blkcnt,supunt,supcnt,untcnt
integer :: simflg,simsup,simunt,simblk
integer :: pres,flow,temp,cntr,othr,enrg,powr,ahum
integer :: nunit,nunct
integer :: stcnt,bcnt,decblk,decunt,decsup
integer :: istat
!------- Logical unit declaration
integer :: lutypr=10,luf=11,lufb=12,luview=13

real,dimension(maxsbk)           :: repinv
real,dimension(maxstv)           :: state
real,dimension(maxunt,mpariu)    :: pvalue
real,dimension(muntib,mpariu)    :: valtem

integer,dimension(maxsbk)        :: scan,freeze,simcnt,supnum
integer,dimension(maxblk)        :: blkunt
integer,dimension(maxunt)        :: untnum,blkin,blkout,blkpar,untyp
integer,dimension(mblkis)        :: blknum
integer,dimension(mbndis)        :: scnt,scnt2,scnt3
integer,dimension(muntib)        :: numtem,typtem,tin,tout,tpar
integer,dimension(minoiu)        :: iident,iodent
integer,dimension(mpariu)        :: parnum
integer,dimension(mdmbnd)        :: bndry

integer,dimension(maxsbk,mdmrep) :: ntyp,nindx,repvar
integer,dimension(maxunt,minoiu) :: indxin,indxot,intype,ottype
integer,dimension(muntib,minoiu) :: xintem,xottem

end module hvaccomm

