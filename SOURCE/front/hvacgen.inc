!===============================================================================

! hvacgen.inc : include file

!===============================================================================

character(len=80) :: item
character(len=80) :: line,messag
logical           :: same,abort,reask
real              :: anumbr
integer           :: report,intnum,status
integer           :: crtflg,edflg
integer           :: rtn_crunit,rtn_crblk,rtn_crsup,rtn_crsim,rtn_recrt
integer           :: rtn_fsave,rtn_readin,rtn_rdunt,rtn_rdblk,rtn_rdsim, &
                     rtn_savunt,rtn_savblk,rtn_savsup,rtn_savsim
integer           :: rtn_vewunt,rtn_vewblk,rtn_vewsim,rtn_vewstruct, &
                     rtn_varval,rtn_bound,rtn_rptvar,rtn_errors,rtn_revew, &
                     rtn_struct
integer           :: rtn_edtitl,rtn_str,rtn_edval,rtn_edbnd,rtn_insert, &
                     rtn_delete,rtn_edrep,rtn_rpchng,rtn_rpinrt,rtn_rpdelt, &
                     rtn_ederr,rtn_reedt,rtn_edstr,rtn_extblk
integer           :: rtn_inssim,rtn_insunt,rtn_insblk,rtn_recalc,rtn_typinf, &
                     rtn_delsim,rtn_delunt,rtn_delblk,rtn_repsim,rtn_repunt
integer           :: i,j,k,l,m,n,ii,jj,kk,ll,mm,nn
integer           :: lincnt,iok,inswer,nsuper,iscnt                      

