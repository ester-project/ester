      subroutine readkx(lun,kurtab,ix,iz,md,mt,rlk,tlk,opkx)
      implicit double precision (a-h,o-z)
c
c    read the X-extrapolated Kurucz-tables from
c    the file kurtab, which was produced
c    by the program 'dopxext.f'
c    derived from s/r 'dreadk.f' (introduced new argument ix)
c
c    History:
c
c       17.5.1993 creation
c
c       13/11/95: modified fort OPAL95 tables
c
c    Last modification: 13/11/95
c       
c    input variables:
c    ----------------
c    lun  ..... logival unit number for file
c               kurtab
c    kurtab ... character variable, which defines 
c               the filename to be read
c    ix     ... index-dimension of array opkx
c               (= ntab ... nr. of different X-values (7))
c    iz     ... used for index-dimension 
c               of array opkx
c               nr. of different Z values (13)
c    md     ... used for  array dimension of rlk
c               nr. of different density values (nxir=19) 
c               pro temperature value
c    mt    .... used for array dimension of tlk
c               nr. of different temperature values (30)
c    
c    output parameters:
c    ------------------
c    rlk ...... array(md) containing the density values
c               actual rlk = log10(rho/T6^3)
c    tlk ...... array(mt) containing the temperature values
c    opkx...... array(md,mt,ix,iz) containing the
c               opacity values
c
      parameter(nz=13)
      character*(*) kurtab
      integer lun,ix,iz,md,mt
      dimension rlk(md),tlk(mt),opkx(md,mt,ix,iz)
c
c     table dimension for s/r opintc{f} and opints
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyifi,mdi,
     +                nti,iali
c
c
      open(lun,file=kurtab,status='old',form='unformatted',err=9030)
c
      read(lun)ixr,izr,mdr,mtr
c
c     check array boundary
c
      if(ixr.ne.ntabi)goto 9705
      if(izr.gt.nzvai)goto 9710
      if(mtr.ne.mt)goto 9720
      if(mdr.ne.mdi)goto 9730
c
c     read temperature values
c
      read(lun)(tlk(j),j=1,mtr)
c
c     read density values
c
      read(lun) (rlk(i),i=1,mdr)
c
c     read density and opacity values
c
      do 2001 l=1,izr
        do 2010 itab=1,ixr
         read(lun) ((opkx(i,j,itab,l),i=1,mdr),j=1,mtr)
2010    continue
2001  continue
c
      close(lun)
c
      return
c
9030  print *,'readkx: open error of ',kurtab
      stop
9705  print *,'readkx: inconsistency-error of opkx,ix= ',ix 
      stop
9710  print *,'readkx: inconsistency-error of opkx,iz= ',iz 
      stop
9720  print *,'readkx: inconsistency-error of tlk,mt= ',mt 
      stop
9730  print *,'readkx: inconsistency-error of rlk,md= ',md 
      stop
c
      end
