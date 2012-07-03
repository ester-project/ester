      program opdos
      implicit double precision (a-h,o-z)
c
c     History:
c
c     8.2.1993: creation, test-driver for interpolation over a 
c                         given domain 
c
c     9.2.1993: use s/r opintd to get partial derivative
c               values oph & opz
c
c    10.2.1993: included iexp flag, iexp > 0 if interpolated
c               point lies in the extrapolation domain
c
c    3.3.1993:  changed from dlg->rlg
c
c    5.3.1993:  changed filename-treatment
c               added new argument to opinit, tabnam
c               which defines filename to be used
c
c   12.3.1993:  new argument ier in s/r opintd
c               if ier.gt.0 then result is not valid !
c
c   22.5.1993   use new s/r opintf (fast)
c               new argument 'iorder' (= nr. of used tables for Akima-int.)
c
c   8.11.1993:  inserted ieee-floating handling
c               new argument imode to select algorithm:
c               if imode = 0 -> 677 only (minimum norm)  (s/r opint{c,f})
c                  imode = 1 -> 677 + birational splines (s/r opint{c,f} + opints)
c
c   6.1.1994:   binary outputfiles (unformatted output)
c
c   Last modification:
c   20/08/97
c
c     external myhandleri
c     external myhandleru
c     external myhandlero
c     external myhandlerd
c
      character*79 tabnam
c      
      parameter(nn=200)
      dimension tabnam(2)
      dimension xi(nn),yi(nn),
     +          opalg(nn,nn),opr(nn,nn),opt(nn,nn),
     +          opx(nn,nn),opz(nn,nn)
c
      data iexp/0/
      data tabnam /'OPINTPATH_95','OPINTPATH_AX'/
      data iorder /4/
      data imode  /2/     ! rat. splies + electron conduction
c
c     setting up floating-point exception handler
c
c     ieeer = ieee_handler ('set', 'invalid', myhandleri)
c     if (ieeer .ne. 0) print *, 'ieee_handler i cant be set !!'
c     ieeer = ieee_handler ('set', 'division', myhandlerd)
c     if (ieeer .ne. 0) print *, 'ieee_handler d cant be set !!'
c     ieeer = ieee_handler ('set', 'overflow', myhandlero)
c     if (ieeer .ne. 0) print *, 'ieee_handler o cant be set !!'
c     ieeer = ieee_handler ('set', 'underflow', myhandleru)
c     if (ieeer .ne. 0) print *, 'ieee_handler u cant be set !!'
c
c     print *,'get machine precision'
      call maceps(drelpr)
c
c
      itab=1
      write(6,'(A,$)') 
     .        ' Enter low-T tables, 1=Kuru91, 2=Alex94 : '
      read(5,'(i2)',err=9991)itab
      write(6,'(A,$)')' Enter rlg start: '
      read(5,'(d12.5)',err=9991)xa
      write(6,'(A,$)')' Enter rlg end  : '
      read(5,'(d12.5)',err=9991)xe
      write(6,'(A,$)')' Enter # of rlg points: '
      read(5,'(i4)',err=9991)nxi
      write(6,'(A,$)')' Enter tlg start : '
      read(5,'(d12.5)',err=9991)ya
      write(6,'(A,$)')' Enter tlg end   : '
      read(5,'(d12.5)',err=9991)ye
      write(6,'(A,$)')' Enter # of tlg  points: '
      read(5,'(i4)',err=9991)nyi
 321  write(6,'(A,$)')' Enter X-Value [0.-.8 ]: '
      read(5,'(d12.5)',err=9991)xval
      if(xval.lt.0.d0.or.xval.gt.8d-1)goto 321
 322  write(6,'(A,$)')' Enter Z-Value [0.-.1]: '
      read(5,'(d12.5)',err=9991)zval
      if(xval.lt.0.d0.or.xval.gt.8d-1)goto 322

c     initialize optables (read them)
c     print *,'initialize tables'
      call opinit(drelpr,iorder,tabnam(itab),imode)
c
c     set grid coordinates
c
      dx = ( xe  - xa)/(nxi-1)
      dy = ( ye -  ya)/(nyi-1)
      open(8,form='unformatted')
      open(9,form='unformatted')
      do 2005 k=1,nxi
         xi(k)=xa+(k-1)*dx
c        write(8,'(f9.4,$)')xi(k)
 2005 continue
      do 2200 k=1,nyi
         yi(k)=ya+(k-1)*dy
c        write(9,'(f9.4,$)')yi(k)
 2200 continue
      write(8) nxi,(xi(k),k=1,nxi)
      write(9) nyi,(yi(k),k=1,nyi)
      close(8)
      close(9)
c
c now do the interpolation
c
      do 5001 j=1,nyi
         do 5010 i=1,nxi
c
            rlg=xi(i)
c           rlg=xi(i)-3.d0*yi(j)+18.d0
            tlg=yi(j)
c
c           now call the actual interpolation-s/r 
c           (no extrapolation checking: ier=0 always)
            call opints(xval,zval,tlg,rlg,opalg(i,j),opr(i,j),opt(i,j),
     +                  opx(i,j),opz(i,j),iexp,ier)
c
            if(ier.gt.0)write(6,'(A,i3)')
     +                  'opdof: error in interpolation s/r,ier= ',ier
c           check if point lies in extrapolation domain
c           if(iexp.gt.0)print 6000,'extrapolated points at ( ',
c    +                           'rlg=',xi(i),' tlg=',tlg,
c    +                           ') : ',iexp
c
c           output the results
c
c
 5010    continue
 5001 continue
      open(80,form='unformatted')
      open(81,form='unformatted')
      open(82,form='unformatted')
      open(83,form='unformatted')
      open(84,form='unformatted')
c
      write(80)((opalg(i,j),i=1,nxi),j=1,nyi)
      write(81)((opr(i,j),i=1,nxi),j=1,nyi)
      write(82)((opt(i,j),i=1,nxi),j=1,nyi)
      write(83)((opx(i,j),i=1,nxi),j=1,nyi)
      write(84)((opz(i,j),i=1,nxi),j=1,nyi)
c
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
c
      stop
 6000 format(2A,F7.3,A,F7.3,A,i2)
c
 9991 print *,'opdos: input error'
      stop
      end
c
d     integer function myhandlerd (sig, code, context)
d     integer sig, code, context(5)
d     print 1
d     myhandlerd=0
d     return
d 1   format(' division by zero !!')
d     end
c
d     integer function myhandleru (sig, code, context)
d     integer sig, code, context(5)
d     print 1
d     myhandleru=0
d     return
d 1   format(' underflow !!')
d     end
c
d     integer function myhandlero (sig, code, context)
d     integer sig, code, context(5)
d     print 1
d     myhandlero=0
d     return
d 1   format(' overflow !!')
d     end
c
d     integer function myhandleri (sig, code, context)
d     integer sig, code, context(5)
d     print 1
d     myhandleri=0
d     return
d 1   format(' invalid !!')
d     end
