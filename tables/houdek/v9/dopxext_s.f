      program opxexs
      implicit double precision (a-h,o-z)
c
c     History:
c
c     13.5.1993: creation
c                evaluate the X-derivatives [d(opa)/dX]
c                at the canonical opal-tablepoints AND
c                extrapolate these X-derivatives for the 
c                Kurucz-temperature-tablepoints using
c                Akima-Interpolation (linear extrapolation);
c                following by an integration process
c                to evaluate Kurucz-opacity-values
c                for other X-values than the already
c                given value X=0.7.
c                Integration is done with a Runge-Kutta-Fehlberg
c                4-5th order, where for the evaluation of the
c                RHS Akima-Interpolation is used.
c
c     19.5.1993: use Shepard's Interpolation-formula
c                for the extrapolation-process
c 
c     20.5.1993: iwk-array (triangle indices) will only be
c                written one (triangulation is the same for
c                all tables; independend from X and Z)
c
c     25.5.1993: changed s/r intrpl[d] (orig. akima) -->
c                        s/r uvip3p[d] (impr. akima)
c
c      8.6.1993: mt=22 -> mt=31 (index of new Kurucz fittingpoint tk_fit=4.05)
c                (must be equal to ylo - output of kur2l2g)
c
c                introduced new parameter nyis, defining
c                the starting # of OPAL-temperature, for 
c                fit with Kurucz-table at tk_fit=4.05
c                (nyis=7 for t6=0.012 ->  tl_fit=4.079)
c
c     28.1.1994: mt=31 -> mt=30 (index of new Kurucz fittingpoint tk_fit=4.0)
c                (must be equal to ylo - output of kur2l2g)
c
c                fit with Kurucz-table at tk_fit=4.0 -> mt=30
c                (nyis=6 for t6=0.011 ->  tl_fit=4.041)
c
c    Last modification:
c    29.1.1994
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c              parameters ntab,nzva,nxir,mt                !
c     MUST have same values than in the s/r f (RHS)        !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      parameter(ntab=7,nzva=10)
      parameter(nxir=17,nyir=50,nyis=6)
      parameter(niwk=58000,nwk=32500)
      parameter(ndat=1686)
      parameter(lun=31,nout=6)
      parameter(md=25,mt=30)
      parameter(ncp=7)
c
c     definitions for the rksuite stuff
c     ---------------------------------
c
      integer           neq, lenwrk, method
      double precision  zero, one, two, four
      parameter         (neq=1,lenwrk=32*neq,method=2)
      parameter         (zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
d     double precision  waste, hnext, stpsok, totf, stpcst
      double precision  hstart, t, tend, const,
     &                  tlast, tol, tstart, twant
      integer           l, uflag
      logical           errass, mesage
      double precision  thres(neq), work(lenwrk), y(neq), ymax(neq),
     &                  yp(neq), ystart(neq),xtab(ntab)
c     .. external subroutines ..
      external          f, setup, stat, ut
c     .. intrinsic functions ..
      intrinsic         atan, cos, sin
c
c     end of definitions for the rksuite stuff
c     ----------------------------------------
c
      character*79 optabe,paderi,ivadat,kurtab,kurtabx
c
      dimension xd(ndat),yd(ndat),zd(ndat),
     +          iwk(niwk),wk(nwk)
c
c     dimension for the X-derivatives
c
      dimension xderiv(ntab,nyir,nxir,nzva),dummy(ntab),tlgtab(nyir),
     +          zz(ntab),zitab(nxir),opasol(ntab),
     +          rlgtab(nxir),otab(2),zotab(2),xv(7),yv(7),zv(7)
c
c     dimension for the Livermore-tables
c
      dimension rlg(nxir,ntab,nzva)
      dimension tlg(nyir,ntab,nzva)
      dimension opa(nxir,nyir,ntab,nzva)
      dimension ival(nyir,ntab,nzva)
c
c     dimension for the Kurucz-tables
c
      dimension opk(nzva,md*mt)
      dimension opkx(nxir,mt,ntab,nzva)
      dimension tlk(mt)
      dimension rlk(nzva,md*mt)
c
      data xtab
     +    /0.0d0,0.1d0,0.2d0,0.35d0,
     +     0.5d0,0.7d0,0.8d0/
c
      common /xderix/ xderix(ntab,mt,nxir,nzva)
      common /tables/ xtab
      common /tabsel/ ls,is,js
c
c     table dimension for s/r readi & readkx
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyifi,mdi,
     +                nti,iali
c
c     pointers for the interpolation routine masube
c     nt = # of triangles
      common /jpoint/ nt,jipt,jipl,jiwp,jind,jngp,jigp
      common /ipoint/ ipd,ial,ibe,iga,iei,ials,ibes,igas,izx,izy

c
      nzvai=nzva
      ntabi=ntab
      nxiri=nxir
      nyiri=nyir
      mdi  =md
c
      call maceps(drelpr)
c
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat,kurtab.dat)
c     first read the pathname of the file, which declares
c     the the absolute pathnames of the above mentioned input files
c     you can assign a logical name to OPINTPATH; eg. in a csh:
c     setenv OPINTPATH <absolute_path_name>
c
      open(29,file='OPINTPATH_92X',
     +form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      read(29,'(A)')paderi
      read(29,'(A)')ivadat
      read(29,'(A)')kurtab
      read(29,'(A)')kurtabx
      close(29)
      print *,'Using tables:'
      print *,optabe
      print *,paderi
      print *,ivadat
      print *,kurtab
      print *,kurtabx
c
c     write(nout,'(A,$)')' Enter tension parameter: '
c     read(5,'(d12.5)')tp
c
      tp=0.d0
c
c     read the opacity-tables (Livermore + Kurucz)
      call readi(lun,optabe,ntab,nxir,nyir,nzva,rlg,tlg,opa,ival)
      call readk(lun+2,kurtab,nzva,md,mt,rlk,tlk,opk)
c
c     setup the opal-temperature values
c
      do 3510 i=1,nyir
       tlgtab(i)=tlg(i,ntab,nzva)
3510  continue
c
c     setup the opal-density values
c
      do 3515 i=1,nxir
       rlgtab(i)=rlg(i,ntab,nzva)
3515  continue
c
c
c
c--------------------------------------------------------
c----             perform the derivation of X         ---
c--------------------------------------------------------
c
      print '(/,a)','evaluate X-derivatives of opal-tables'
      print *      ,'and extrapolate them for Kurucz-tables'
      print '(a,/)','--------------------------------------'
c     repeat for nzva tables (normaly: nzva=10)
      do 4010 l=1,nzva
       print *,'processing table Z= ',l
c      repeat for nxir density points (normaly: nxir=17)
       do 4020 i=1,nxir
c       repeat for (nyir-nyis) temp. points (normaly: nyir=50)
        do 4030 j=nyis,nyir
c        fillup the array for the numerical derivation of X
         do 4040 itab=1,ntab
           zz(itab)=opa(i,j,itab,l)
4040     continue
c
c----    now evaluate the X-derivatives (xderiv) using the ntab different
c----    X-values (xtab) as the abszisse-table-points using a modified 
c----    Akima-Interpolation algorithm
c
c        call intrpd(6,ntab,xtab,zz,ntab,xtab,dummy,
c    +                 xderiv(1,j,i,l)                 )
c
         call uvip3d(3,ntab,xtab,zz,ntab,xtab,dummy,
     +                 xderiv(1,j,i,l)                 )
c
4030    continue
4020   continue
c
c----------------------------------------------------------------
c----   extrapolation of the X-derivatives using Shepard's method
c----------------------------------------------------------------
c
c          ---> rlg (i)
c    (jk)
c    tlg
c     ^    * * * * * * * * * * * *
c     |    * 1 * * * * * * * * * *
c     |    2 3 4 * * * * * * * * *
c     |    5 6 7 * * * * * * * * *
c            x
c
c * ..... given table-points
c x ..... new extrapolated point
c 1-7 ... ncp (=7) used table-points for extrapolation of x
c
c 
c
c       do the extrapolation for all ntab X-values
c
        do 5030 itab=1,ntab
c        mt temperature values
         do 5020 jk=1,mt
c         nxir-2 density values (1-st and last with Akima)
          do 5010 i=2,nxir-1          
            xv(1)=rlgtab(i)
            yv(1)=tlgtab(nyis+2)
            zv(1)=xderiv(itab,nyis+2,i,l)
            do 5040 icp=2,4
              xv(icp)=rlgtab(i-1+icp-2)
              yv(icp)=tlgtab(nyis+1)
              zv(icp)=xderiv(itab,nyis+1,i-1+icp-2,l)
5040        continue
            do 5045 icp=5,ncp
              xv(icp)=rlgtab(i-1+icp-5)
              yv(icp)=tlgtab(nyis)
              zv(icp)=xderiv(itab,nyis,i-1+icp-5,l)
5045        continue
c
            x0=rlgtab(i)
            y0=tlk(mt-jk+1)
            xderix(itab,jk,i,l) = shep(ncp,xv,yv,zv,x0,y0)
5010      continue
c
c         use Akima-Interpolation (linear extrapolation)
c         for the first and last point
c
          n=2
          otab(1)=rlgtab(1)
          otab(2)=rlgtab(nxir)
          do 5050 ii=1,nxir-2
             zitab(ii)=xderix(itab,jk,ii+1,l)
5050      continue
c         using 5-th degree
          nd =  5
          call uvip3p(nd,nxir-2,rlgtab(2),zitab,
     +                   n,     otab,     zotab)
          xderix(itab,jk,1,l)=zotab(1)
          xderix(itab,jk,nxir,l)=zotab(2)
5020     continue
5030    continue
4010  continue
c
c--------------------------------------------------------
c---- end of the derivation of X -> results in xderiv ---
c----                      and                        ---
c---- end of the extrapolation of d(opa)/dX -> xderix ---
c--------------------------------------------------------
c
c
c======================================================================
c
c     just for testing with idl
c     output one extrapolated array d(opa)/dX [xderix]
c
d     l=1
d     itab=1
d     print *,'output test-table l=',l
d     open(8)
d     open(9)
d     open(90)
d     open(91)
d     write(8,'(17f9.5)')(rlg(i,itab,l),i=1,nxir)
d     write(9,'(50f9.5,$)')(tlk(j),j=1,mt)
d     write(9,'(50f9.5)')(tlgtab(j),j=nyis,nyir)
c
d     do 120 itab=1,ntab
c
d        do 112 j=1,mt
d         write(90,'(17f10.5),$')(xderix(itab,j,i,l),i=1,nxir)
d         write(91,'(17f10.5)')(91.d0,i=1,nxir)
d 112    continue
d        do 111 j=nyis,nyir
d         write(90,'(17f10.5)')(xderiv(itab,j,i,l),i=1,nxir)
d         write(91,'(17f10.5)')(96.d0,i=1,nxir)
d 111    continue
c
d 120 continue
c
d     close(8)
d     close(9)
d     close(90)
d     close(91)
c
c======================================================================
c
c     ---------------------------------------------------
c     now do the integration process using the rksuite-lib
c     ---------------------------------------------------
c
c
c  set the initial conditions.  note that tend is taken well past
c  the last output point, tlast.  when this is possible, and it
c  usually is, it is good practice.
c
c  set error control parameters.
c
      tol = 5.0d-5
      do 20 l = 1, neq
         thres(l) = 1.0d-10
   20 continue
c
      print '(/,a)','integrate X-derivatives of Kurucz-tables'
      print '(a,/)','----------------------------------------'
c
c----- start of triple-loop ---------
c
c     rlk-loop is running faster than tlk-loop !
c     (look in s/r kur2l2g or s/r readk)
c
      do 6001 ls=1,nzva
        print *,'processing table: Z= ',ls
        do 6010 js=1,mt
          do 6020 is=1,nxir
c
c     set index-pointer for opk-array
      ip = (js-1)*md+is
c
d     write(*,'(/a,i4,2f6.3)') 'integrate at triple: ',
d    +                         ls,rlk(ls,is),tlk(js)
c
      tstart    = xtab(1)
      ystart(1) = zero
      tlast     = xtab(ntab)
      tend      = xtab(ntab)+0.1d0
      opasol(1) = zero
c
c
c  call the setup routine. because messages are requested, mesage = .true.,
c  there is no need later to test values of flags and print out explanations.
c  in this variant no error assessment is done, so errass is set .false..
c  by setting hstart to zero, the code is told to find a starting (initial)
c  step size automatically .
c
      mesage = .true.
      errass = .false.
      hstart = zero
      call setup(neq,tstart,ystart,tend,tol,thres,method,'usual task',
     &           errass,hstart,work,lenwrk,mesage)
c
c
c  compute answers at ntab output points for constant
c  Z (1<=ls<=nzva), rlg (1<=is<=nxir), tlg (1<=js<=mt).
c
c
      do 40 itab = 2, ntab
         twant = xtab(itab)
         call ut(f,twant,t,y,yp,ymax,work,uflag)
         opasol(itab)=y(1)
c
         if (uflag.gt.2) go to 60
c
c        success. t = twant. output computed solution components.
c        write (*,'(1x,f6.3,3x,f9.4)') t, y(1)
   40 continue
c
c  now add the integration constant.
c  the solution must fit the point given by the table for constant
c  X=0.7 (itab=6) for each triple (rlg,tlg,Z).
c
      itab  = 6
      const = opk(ls,ip) - opasol(itab)
d     write(*,'(/,a,f6.3,/)') 'integration const.: ',const
c
      do 2700 itab=1,ntab
         opkx(is,js,itab,ls) = opasol(itab) + const
c        write (*,'(1x,f6.3,3x,f9.4)') xtab(itab), opkx(is,js,itab,ls)
 2700 continue
c
c  the integration is complete or has failed in a way reported in a
c  message to the standard output channel.
c
   60 continue
c
c
c  ymax(l) is the largest magnitude computed for the solution component
c  y(l) in the course of the integration from tstart to the last t.  it
c  is used to decide whether thres(l) is reasonable and to select a new
c  thres(l) if it is not.
c
d     write (*,'(a/)') '             ymax(l) '
d     do 80 l = 1, neq
d        write (*,'(13x,1pe8.2)')    ymax(l)
d  80 continue
c  the subroutine stat is used to obtain some information about the progress
c  of the integration. totf is the total number of calls to f made so far
c  in the integration; it is a machine-independent measure of work.  at present
c  the integration is finished, so the value printed out refers to the overall
c  cost of the integration.
c
d     call stat(totf,stpcst,waste,stpsok,hnext)
d     write (*,'(/a,i10)')
d    &  ' the cost of the integration in evaluations of f is', totf
c
c
c-----------end of triple-loop ------
c
 6020     continue
 6010   continue
 6001 continue
c
c  ------------------ idl test-output for a certain table ----
c
c     select test-array
c
d     ls=8
c
d     open(92)
d     do 7030 itab=1,ntab
d      do 7020 js=1,mt
d       write(92,'(17f10.5)')(opkx(is,js,itab,ls),is=1,nxir)
d7020  continue
d7030 continue
d     close(92)
c
c  -------- end of idl test-output ------------------------
c
c  --------  Partial derivatives calculation --------------
c
c  Now the usual evaluation of the partial derivatives
c  for the actual opacity-interpolation will be done
c  (this is the same scheme as in the s/r pderiv)
c
c     open file for partial derivatives
      open(lun,file=paderi,status='UNKNOWN',form='unformatted')
c
      print '(/,a)','evaluate partial derivatives opal+Kurucz-tables'
      print '(a,/)','------------------------------------------------'
      write(nout,8610)
c     repeat for nzva tables (normaly: nzva=10)
      do 8010 l=1,nzva
c
c       define xd- ,yd- and zd- values (each of size ndat)
c       and define the xx- and yy- arrays
c       ( nxi x nyi coordinates at the desired z-values)
        do 8001 itab=1,ntab
          nd = 0
          do 1002 j=nyis,nyir
           do 1001 i=1,nxir
               nd = nd + 1
               zd(nd) = opa(i,j,itab,l)
               xd(nd) = rlg(i,itab,l)
               yd(nd) = tlg(j,itab,l)
 1001      continue
 1002     continue
c
c         now fillup the Kurucz-tables opkx (extrapolated above)
c
          do 1225 js=1,mt
           do 1221 is=1,nxir
            nd=nd+1
            zd(nd) = opkx(is,js,itab,l)
            xd(nd) = rlk(l,is)
            yd(nd) = tlk(js)
 1221      continue
 1225     continue
c
          write(nout,8600)l,itab,nd
c
c
c         interpolate at xi=-1.,yi=1. (dummy values)
c         and evaluate partial derivatives
c         stored in the wk array
c
          call masube(drelpr,1,0,nd,xd,yd,zd,tp,1,1,-1.d0,1.d0,z,1,
     +                nout,iwk,wk)
c
c         output of pointers, indices of triangles and derivative-values
          if((itab.eq.1).and.(l.eq.1)) then
             write(lun)nzva
             write(lun)ntab
             write(lun)nt
             write(lun)ial
             write(lun)(iwk(m),m=8,3*nt+7)
          endif
          write(nout,'(2i5)')nt,ial
          write(lun)(wk(m),m=5,ial-1)
c
8001    continue
8010  continue
c
c
      close(lun)
c
c-------------------------------------------------------------
c
c     now output the new interpolated values
c
      print*,'****************************************'
      print*,'output Kurucz-tables to file ',kurtabx
      print*,'****************************************'
c
      open(lun+1,file=kurtabx,status='unknown',
     +           form='unformatted')
c
      write(lun+1) ntab,nzva,nxir,mt
c
c     output temperature values
c
      write(lun+1) (tlk(j),j=1,mt)
c
c     output density values
c
      write(lun+1) (rlk(10,i),i=1,nxir)
c
c     output density and opacity values
c
      do 4501 l=1,nzva
       do 4510 itab=1,ntab
        write(lun+1)((opkx(is,js,itab,l),is=1,nxir),js=1,mt)
4510   continue
4501  continue
      close(lun+1)
c
c-------------------------------------------------------------
c
      stop
c
8600  format('evaluate partial derivative values of table:',
     +        2i3,i5,$)
8610  format(46x,'z  x  nd   nt   ial')
9030  print *,'dopinit: error in opening of OPINTPATH'
      print *,'dopinit: try to assign a logical name to OPINTPATH'
      print *,'dopinit: with setenv OPINTPATH <absolute_path_name>'
      print *,'dopinit: in a csh where <absolute_path_name> is the'
      print *,'dopinit: file containing the absolute pathnames of'
      print *,'dopinit: the input files optabe.bin pderivs.dat'
      print *,'dopinit: and ival.dat'
      stop
c
      end
c
c ------------------END of Main program ---------------------
c
      subroutine f(t,y,yp)
      implicit double precision (a-h,o-z)
c
      double precision  t
      double precision  y(*), yp(*)
c
c     this is the right hand side (RHS)
c     of the Runge Kutta integration-sheme.
c     continuously values of the depend variable yp
c     are achieved trough Akima Interpolation between
c     the table points
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                     parameter values                     !
c     MUST have same values than in the calling s/r opxext !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      parameter(ntab=7,nzva=10)
      parameter(nxir=17,nyir=50,nyis=6)
      parameter(md=25,mt=30)

      common /xderix/ xderix(ntab,mt,nxir,nzva)
      common /tables/ xtab(ntab)
      common /tabsel/ ls,is,js
c
c     using 3-th degree
      nd =  3
c
      call uvip3p(nd,ntab,xtab,xderix(1,js,is,ls),
     +               1   ,t   ,yp(1)              )

      return
      end
