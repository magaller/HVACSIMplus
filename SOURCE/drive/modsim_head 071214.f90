! ***********************************************************************

 module precision1

   implicit none
   save
   integer,parameter :: pp = selected_real_kind(6)! single precision
!   integer,parameter :: pp = selected_real_kind(15)! double precision

 end module precision1

! ***********************************************************************

 module modsim_head

!------- The following parameters are used to define the array sizes
!        in hvacgen, slimcon and modsim.

!    maxsbk = maximum superblocks in the simulation plus 1
!    maxblk = maximum blocks in the simulation
!    maxunt = maximum units in the simulation
!    maxstv = maximum state variables in the simulation
!    maxdeq = maximum differential equations in the simulation
!    maxpar = maximum unit parameters in the simulation
!    maxsav = maximum saved variables in the simulation
!    maxbnd = maximum boundary variables in the simulation
!    mdmbnd = dimension of a boundary variable calculation array
!    mblkis = maximum blocks in a superblock
!    muntib = maximum units in a block
!    mdeqis = maximum differential equations in a superblock
!    mdeqiu = maximum differential equations in a unit
!    mseqis = maximum simultaneous equations in a superblock
!    mseqib = maximum simultaneous equations in a block
!    mbndis = maximum boundary conditions in a single superblock
!    minoib = maximum inputs or outputs in a block
!    minoiu = maximum inputs or outputs in a unit
!    mrptis = maximum reported variables in a superblock
!    mpariu = maximum parameters in a unit
!    mdmrep = dimension of a reported variable calculation array
!    maxtyp = maximum number of TYPES
!    nlscre = number of lines on screen - used for scrolling

!    maxlbl = maximum number of labels
!    minois = maximum inputs or outputs in a superblock
!    munitb = maximum units in a block
!    msumary = MAX(MSEQIS,MSEQIB,MINOIU,MRPTIS,MAXLBL)

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

! The maximum time (seconds) for writing data on the file, hvacsim.sum,
! is set by wrtmax.  By changing this value, the size of output can be
! adjusted.  Refer to subroutine report in modino.f90.

 real                              :: wrtmax=3600.0
                                                                                     
 logical                           :: frz                                            
 character(len=80)                 :: title                                          
 character(len=8)                  :: typarff='typar.dat'                            
                                                                                     
 integer                           :: nblock,nbound                                  
 integer                           :: npar,iprint                                    
 integer                           :: nd,mtime,mtstep,itime                                
 integer                           :: nsblok,nu,nstate,init,nsaved                   
 integer,dimension(maxsbk)         :: korder,ndsb,msbt,mstep,ideflg,&                
                                      ifzopt,insopt,ninsb,jsslv,nreprt,&             
                                      njsslv,nsuper,nssolv,kstep                     
 integer,dimension(maxsbk,7)       :: ip                                             
 integer,dimension(maxsbk,8)       :: ipa                                            
 integer,dimension(maxsbk,mblkis)  :: isuper                                         
 integer,dimension(maxsbk,minois)  :: insb                                           
 integer,dimension(maxsbk,mseqis)  :: issolv,jssolv                                  
 integer,dimension(maxsbk,mdeqis)  :: indesb                                         
 integer,dimension(maxsbk,mrptis)  :: ireprt                                         
 integer,dimension(maxsbk,mrptis,2):: ident                                          
 integer,dimension(maxblk)         :: noutb,ninb,njsolv,nsolve,&                     
                                      nunits,ibstat                                  
 integer,dimension(maxblk,mseqib)  :: isolve,jsolve                                  
 integer,dimension(maxblk,minoib)  :: ioutb,inb                                      
 integer,dimension(maxblk,muntib)  :: iblock                                         
 integer,dimension(maxunt)         :: nin,nout,nde,jpar,isaved,iunits                
 integer,dimension(maxunt,minoiu)  :: in,iout                                        
 integer,dimension(maxunt,mdeqiu)  :: inde                                           
 integer,dimension(maxdeq)         :: idevar                                         
 integer,dimension(maxstv)         :: isign                                          
 integer,dimension(maxbnd)         :: ibound                                         
 integer,dimension(maxlbl)         :: ndent,ndnt
 integer,dimension(maxstv)         :: icheck                                         
 integer,dimension(maxdeq)         :: idechk                                         
 integer,dimension(mblkis)         :: ufzflg                                         
 real                              :: time,tstep,ttime,tmin,&                  
                                      rtolx,atolx,tprton,tprtof,xtola
 real,dimension(maxsbk)            :: dnk,fnk,tnext,sbtime,step,treprt
 real,dimension(maxsbk,7)          :: ah,gh,a,g                                      
 real,dimension(maxsbk,8,8)        :: am                                             
 real,dimension(maxdeq)            :: spred                                          
 real,dimension(maxdeq,7)          :: spast                                          
 real,dimension(maxpar)            :: par                                            
 real,dimension(maxsav)            :: saved                                          
 real,dimension(maxstv)            :: stold,state,tstate                             
! ------------------ shell model ----------------------------------------
 integer                           :: ndn,nsky,nhor,ntoa,nwoa,npoa,nvwind            
 integer,dimension(maxstr)         :: norda,nctfa                                    
 real                              :: tshell,tzero,hdn,hsky,hhor,&                   
                                      ssazm,ssalt,toa,woa,poa,vwind                  
 real,dimension(maxstr)            :: u                                              
 real,dimension(0:maxnrf,maxstr)   :: xc,yc,zc
 real,dimension(maxord,maxstr)     :: rj                                             
 real,dimension(maxzn)             :: fqrad                                          
 real,dimension(maxns,maxzn)       :: hisc,hisr,as,fview,solint,transm,&             
                                      sc,absis                                       
! ----------------- BLOCK DATA items ------------------------------------
 character(len=4),dimension(maxlbl):: lbls=(/'pres','flow','temp','ctrl',&           
                                             'othr','enrg','powr','humt'/)           
 integer                           :: ifile1=7,ifile2=8,ifile3=9,ifile4=10,&         
                                      ifile5=11,ifile6=12,inp=5                      
 integer                           :: ifile7=13,ifile8=14                            
 integer                           :: lastit=0                                       
 data   isign/maxstv*1/                                                              
                                                                                     
 end module modsim_head                                                              

! ***********************************************************************

 module cool_plant

   implicit none
   logical     :: cooling
   real        :: start_chiller =  3600.0,&
                  stop_chiller  = 1209600.0   ! 3600 s/h * 24 h/day * 7 days

 end module cool_plant


