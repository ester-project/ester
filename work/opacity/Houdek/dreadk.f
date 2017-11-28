      subroutine readk(lun,kurtab,iz,md,mt,rlk,tlk,opk)
      implicit double precision (a-h,o-z)
c
c    read the Kurucz-tables (new) from
c    the file kurtab, which was produced
c    by the program 'kur2liv.f'
c
c    History:
c
c       12.3.1993 creation
c
c       13/11/95: modified for OPAL95 tables
c
c    Last modification: 13/11/95
c       
c    input variables:
c    ----------------
c    lun  ..... logival unit number for file
c               kurtab
c    kurtab ... character variable, which defines 
c               the filename to be read
c    iz     ... used for first index-dimension 
c               of array opk
c               nr. of different Z values (13)
c    md     ... used for  array dimension of rlk
c               nr. of different density values (27) 
c               pro temperature value
c    mt    .... used for array dimension of tlk
c               nr. of different temperature values (30)
c    
c    output parameters:
c    ------------------
c    rlk ...... array(md) containing the density values
c               actual rlk = log10(rho/T6^3)
c    tlk ...... array(mt) containing the temperature values
c    opk ...... array(iz,md*mt) containing the
c               opacity values
c
      parameter(nz=13)
      character*(*) kurtab
      integer lun,iz,md,mt
      dimension rlk(iz,md*mt),tlk(mt),opk(iz,md*mt)
      dimension zval(nz)
c
      open(lun,file=kurtab,status='old',form='unformatted',err=9030)
c
      read(lun)izr,mdr,mtr
c
c     check array boundary
c
      if(izr.gt.nz)goto 9700
      if(izr.ne.iz)goto 9710
      if(mtr.ne.mt)goto 9720
      if(mdr.ne.md)goto 9730
c
c     read Z-values
c
      read(lun)(zval(i),i=1,izr)
c
c     read temperature values
c
      read(lun)(tlk(i),i=1,mtr)
c
c     read density and opacity values
c
      do 2001 i=1,izr
        read(lun) (rlk(i,j),opk(i,j),j=1,mdr*mtr)
2001  continue
c
      close(lun)
c
      return
c
9030  print *,'readk: open error of ',kurtab
      stop
9700  print *,'readk: array-dimension error of zval,nz= ',nz 
      stop
9710  print *,'readk: inconsistency-error of opk,iz= ',iz 
      stop
9720  print *,'readk: inconsistency-error of tlk,mt= ',mt 
      stop
9730  print *,'readk: inconsistency-error of rlk:md,mdr= ',md,mdr 
      stop
c
      end
