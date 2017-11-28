      subroutine readal(lun,alxtab,ix,iz,ir,it,
     +                  rlgalx,tlgalx,opaalx)
      implicit double precision (a-h,o-z)
c
c    read Alexander-tables from
c    the file alxtab, which was produced
c    by the program 'alex2b'
c
c    History:
c
c       06/04/95  creation
c
c    Last modification: 19/08/97
c       
c    input variables:
c    ----------------
c    lun  ..... logical unit number for input-file
c    kurtab ... input-filename
c    ix     ... desired nr. of different X-values (max=7))
c    iz     ... desired nr. of different Z values (max=13)
c    ir     ... desired nr. of different density values rlalx (max=17) 
c    it    .... desired nr. of different temperature values 
c               tlalx (max=23)
c    
c    output parameters:
c    ------------------
c    rlgalx ... array(ir,ix,iz) containing the density values
c    tlgalx ... array(it,ix,iz) containing the temperature values
c    opaalx.... array(ir,it,ix,iz) containing the opacity values
c
      character*(*) alxtab
      integer lun,ix,iz,ir,it
      dimension rlgalx(ir,ix,iz),tlgalx(it,ix,iz),
     +          opaalx(ir,it,ix,iz)
c
      dimension xtab(7),ztab(13)
c
      open(lun,file=alxtab,status='old',form='unformatted',err=9030)
c
c-----read xtab & ztab
      read(lun)izr,(ztab(l),l=1,izr)
      read(lun)ixr,(xtab(k),k=1,ixr)
      if(ixr.gt.ix)goto 9705
      if(izr.gt.iz)goto 9700
      do 4001 l=1,izr
         do 5001 k=1,ixr
c-----------read rlgagx
            read(lun)irr,(rlgalx(i,k,l),i=1,irr)
            if(irr.ne.ir)goto 9730
c-----------read tlgalx
            read(lun)itr,(tlgalx(j,k,l),j=1,itr)
            if(itr.ne.it)goto 9720
c-----------read opaalx
            read(lun)((opaalx(i,j,k,l),i=1,irr),j=1,itr)
 5001    continue
 4001 continue
c
      close(lun)
c
      return
c
9030  print *,'readal: open error of ',alxtab
      stop
9700  print *,'readal: inconsistency-error of iz= ',iz 
      stop
9705  print *,'readal: inconsistency-error of ix= ',ix 
      stop
9720  print *,'readal: inconsistency-error of it= ',it 
      stop
9730  print *,'readal: inconsistency-error of ir= ',ir
      stop
c
      end
