!===============================================================================

 module slim

!=======================================================================

 use hvacsim_par

 character (len=80)                 :: title
 integer                            :: luf=7,lud=8,luty=10,ioerr
 logical                            :: algbraic

 real,dimension(maxstv)             :: state
 real,dimension(maxpar)             :: par
 real,dimension(maxsbk)             :: treprt
 real                               :: atolx,rtolx,ttime,xtol

 integer,dimension(maxsbk)          :: nssolv,nreprt,ifzopt,insopt,nsuper, &
                                       ndeqsb,ntdbvs
 integer,dimension(maxblk)          :: nunits,nsolve,maxibl,maxobl
 integer,dimension(maxunt)          :: nin,isaved,jpar,nout,nde,maxpun
 integer,dimension(minoiu)          :: iodent,iident
 integer,dimension(maxlbl)          :: ndent,ndnt
 integer,dimension(maxdeq)          :: idevar
 integer,dimension(maxbnd)          :: ibound
 integer,dimension(maxsbk,mrptis)   :: ireprt
 integer,dimension(maxsbk,mblkis)   :: isuper
 integer,dimension(maxsbk,mseqis)   :: issolv
 integer,dimension(maxblk,mseqib)   :: isolve
 integer,dimension(maxblk,muntib)   :: iblock
 integer,dimension(maxunt,minoiu)   :: in,iout
 integer,dimension(maxunt,mdeqiu)   :: inde
 integer,dimension(maxsbk,mrptis,2) :: ident

 end module slim

