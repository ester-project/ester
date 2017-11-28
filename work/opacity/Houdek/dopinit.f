c      algorithm 677, collected algorithms from acm.
c      this work published in transactions on mathematical software,
c      vol. 15, no. 4, pp. 365-374.
c
      subroutine opinit(eps,iord,tabnam,imode)
      implicit double precision (a-h,o-z)
c
c  Initialize (read) Opazity-tables ,partial derivative values
c  and index-coordinates
c
c     History:
c     10.6.1992: creation
c     25.8.1992: reduce of iwk-      & wk-arrays to
c                          iwk(3*nt) & wk(2*n) respectively
c                where nt = # of triangles (1568) and n = # of datapoints
c                (nxi*nyi = 17*50 = 850)
c     30.9.1992: included 'status='old'' in open statement for VAX/VMS
c
c     20.12.1992 modified for z-domain interpolation
c                nzva different z-values
c                (z = mass fraction of heavy elements)
c
c                modified opacity values input from
c                'formatted' to 'unformatted'
c                first you have to run the program 'a2b'
c                (ascii-to-binary conversion of opacity-tables)
c
c      25.1.1993 changed variable t6 -> tlg = log10(t)
c                table-values changed to tlg in da2b
c
c       5.2.1993 modified pathname-treatment of inputfiles.
c                A single file will be used to define the
c                absolute pathnames of the 3 inputfiles
c                (optabe.bin, pderivs.dat, ival.dat) line by line
c                in the above order !
c                This file has to be placed in the
c                current working directory and must have
c                the filename OPINTPATH, or may be
c                assigned to a logical name: eg. 
c                under UNIX in a csh:
c              
c                setenv OPINTPATH /dir1/dir2/opacity.paths
c
c       5.3.1993 redefined pathname-treatment
c                defined new argument tabnam, which
c                defines pathname of file, which table
c                has to be used
c
c      11.3.1993 introduced parameter statement
c                for array-dimensions
c                removed ntab from common /tablex/
c                removed nzva from common /tablez/
c                add new common /arrdim/ for
c                consistency-check of the
c                common - array defs, because parameters
c                can not be passed through common-statement
c
c     12.3.1993  modified pderiv.dat read- handling
c                for use for the Kurucz-tables
c                modified common /jpoint/
c                         common /ipoint/
c
c     13.3.1993  use 'D' in 1-st column for
c                debugging commands
c
c     17.5.1993  use d(opa)/dX-extrapolated Kurucz-tables for all X-values
c                (s/r readkx, array opkx)
c                changed parameter md = 17
c
c     20.5.1993: iwk-array (triangle indices) will only be
c                read once (triangulation is the same for
c                all tables; independend from X and Z)
c
c     25.5.1993: included variable 'iorder' in argument list
c                and added iorder to common /tablex/
c
c      6.6.1993: change read-sequence of input-files (tabnam)
c
c      8.6.1993: mt=22 -> mt=31 (index of new Kurucz fittingpoint tk_fit=4.05)
c                (must be equal to ylo - output of kur2l2g)
c
c                introduced new parameter nyis, defining
c                the starting # of OPAL-temperature, for 
c                fit with Kurucz-table at tk_fit=4.05
c                (nyis=7 for t6=0.012 ->  tl_fit=4.079)
c
c                parameter nt   -> 2368 (from dopxext_s.f output)
c                          ndat -> 1275
c
c      7.11.1993: added H. Spath's birational spline algorithm
c                 on a regulary grid s/r rat2d & rbivpd (s/r dopints)
c                 !!!!!!! B U T  !!!!!!!!!
c                 the stored coefficients array ra will allocate
c                 > 10MBte ! (17*75*4*4*7*10)*8 byte
c
c                 new argument parameter imode
c
c                 if imode = 0 -> only minimum norm (677) algorithm (s/r opint{c,f})
c                    imode > 0 -> birational splines + 677 (s/r opints + opint{c,f})
c
c      8.11.1993: change drelpr->10.d0*drelpr for imode >0  prevent 'division by zero'
c                 for the analytically evaluated partial derivatives wrt x
c
c      28.1.1994:  fit with Kurucz-table at tk_fit=4.0 (mt=30)
c                 (nyis=6 for t6=0.011 ->  tl_fit=4.041)
c
c      13/11/95: modified for OPAL95 tables
c
c      23/05/97: new array ztabl = log10(ztab), used for 
c                logarithmic interpolation in Z
c
c     17/08/97: include electron conduction (new variable lec [logical])
c
c     19/08/97: include Alexander95 low-temperature tables
c
c     Last modification:
c     21/08/97
c
c  input parameters
c  
c  eps ............... relative machine precision
c  iord  ............. nr. of used tables for Akima-int. [2<=iorder<=4]
c  tabnam ............ filename, which tables to be used
c  imode ............. selection: if|imode|= 1 -> minimum norm (677) only
c                                   |imode|= 2 -> birational splines only
c                                   |imode|> 2 -> birational splines + 677
c                                    imode < 0 -> disable electron conduction
c
c     setup array-dimensions (max values for OPAL95 tables)
c
      parameter(ntab=7,nzva=13)
      parameter(nxir=19,nyir=70,nyis=5)     ! OPAL95, nyis=7 for KURU91
      parameter(ntal=19)                    ! Alexander95
      parameter(md=19,mt=30)                ! Kurucz91
c     parameter(ndat=1786)
      parameter(ndat=(nyir-nyis+1+mt)*nxir)
c     nyii=nyir-nyis+1+mt
c     nt=2*(nxir+nyii-2)+2*(nxir*nyii-2*(nxir+nyii-2)-1)
c     parameter(nt=3348,ial=2*ndat+5)       !ial=3235 for OPAL95+ALEX94
      parameter(nt=3348)                    ! nt=3024 for OPAL95+ALEX94
c     niwkl=3*nt,nwkl=2*nd=ial-5
      parameter(niwkl=3*nt,nwkl=2*ndat)
c
c---  dimension for s/2 rat2d
c 
      parameter(nyif=nyir+mt-nyis+1)
c
      real ra
      logical*1 lec                  ! logical for electron conduction
      character*80 tabnam
      character*79 optabe,pderiv,ivadat,lowtab
      dimension ival(nyir,ntab,nzva) ! for OPAL tables
c
c     dimension for Alexander-tables
c
      dimension opaal(nxir,ntal,ntab,nzva)
      dimension rlgal(nxir,ntab,nzva)
      dimension tlgal(ntal,ntab,nzva)
c
c     dimension for Kurucz-tables
c
      dimension opkx(md,mt,ntab,nzva)
      dimension tlk(mt)
      dimension rlk(md)
c
c --- definitions for s/r rat2d (rational splines)
c
      dimension p(nxir,nyif),
     +          q(nxir,nyif),
     +          r(nxir,nyif),
c    +          ra(nxir,nyif,4,4),
     +          rsdx(nxir),rsdy(nyif),
     +          rsax(nxir),rsbx(nxir),rscx(nxir),rsrx(nxir),
     +          rsay(nyif),rsby(nyif),rscy(nyif),rsry(nyif),
     +          ops(nxir,nyif,ntab,nzva)
c
c     partial derivatives and triangulation indices
      common /pderiv/ iwk(niwkl),wk(nwkl,ntab,nzva)
c
      common /opadat/ opa(nxir,nyir,ntab,nzva),
     +                rlg(nxir,ntab,nzva),
     +                tlg(nyir,ntab,nzva)
      common /valdat/ ivalo(nyir,ntab,nzva)
      common /tablex/ xtab(ntab),iorder,lec
      common /tablez/ ztab(nzva),ztabl(nzva)
      common /xyzdat/ xd(ndat,nzva),
     +                yd(ndat,nzva),
     +                zd(ndat,ntab,nzva)
      common /machin/ drelpr,toll,eps10
c
c
c---- birational spline coefficients for s/r rbivpd
c---- evaluated in s/r rat2d (s/r opinit)
c
      common /birasp/ pp,tls(nyif),
     +                ra(nxir,nyif,4,4,ntab,nzva)
c
c     common for array-consistency-check of common-arrays
c     between s/r opinit and s/r opintd
d     common /arrdim/ nxiv,nyiv,ndatv,ntabv,nzvav,niwkv,nwkv
c
c     table dimension for s/r opintc{f} and opints
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,iyifi,mdi,
     +                nti,iali
c
      data xtab(1),xtab(2),xtab(3),xtab(4),
     +     xtab(5),xtab(6),xtab(7)
     +    /0.0d0,0.1d0,0.2d0,0.35d0,
     +     0.5d0,0.7d0,0.8d0/
c
      data ztab(1),ztab(2),ztab(3),ztab(4),ztab(5),
     +     ztab(6),ztab(7),ztab(8),ztab(9),ztab(10),
     +     ztab(11),ztab(12),ztab(13)
     +    /0.0d0,  1.0d-4 , 3.0d-4 , 1.0d-3 , 2.0d-3,
     +     4.0d-3, 1.0d-2 , 2.0d-2 , 3.0d-2 , 4.0d-2,
     +     6.0d-2, 8.0d-2 , 1.0d-1/
c
c     compute log10(ztab)
      ztabl(1)=-7.0d0
      do 5 i=2,nzva
        ztabl(i)=dlog10(ztab(i))
    5 continue
c
c     tension parameter for rational splines
      pp = 0.d0
c
c     setup variables for common /arrdim/
d     nxiv=nxir
d     nyiv=nyir
d     ndatv=ndat
d     ntabv=ntab
d     nzvav=nzva
d     niwkv=niwkl
d     nwkv=nwkl
c 
c     setup iorder for common /tablex/
      iorder = iord
c---  relativemachine precision
c---  in case using linear interpolation for part. deriv.
c     drelpr = eps
c     toll   = 1.d0 - eps
c     eps10  = 10.d0 * eps
c---  in case of using analytical part.deriv.
      drelpr = 10.d0*eps
      toll   = 1.d0 - drelpr
      eps10  = drelpr
c
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat,lowtab.dat)
c
      open(29,file=tabnam,form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      read(29,'(A)')lowtab
      read(29,'(A)')pderiv
      read(29,'(A)')ivadat
      close(29)
c      print '(/,a)',' using tables:'
c      print *,optabe
c      print *,lowtab
c      print *,pderiv
c      print *,ivadat
c
c     check which low-T tables we have to deal with
c
      lun=31
      open(lun+2,file=lowtab,status='old',form='unformatted')
      read(lun+2),iii
      close(lun+2)
c
c     read pointer and partial derivatives for subroutine masubi(e)
      open(33,file=pderiv,status='old',
     +        form='unformatted',err=9020)
c
c     # of different z-values (mass fraction of heavy elements)
      read(33)nzvai
      if(nzvai.gt.nzva)then
        print *,'opinit:error in Z-numbers'
        stop
      endif
c
c     # of opacity-tables
      read(33)ntabi
      if(ntabi.gt.ntab)then
        print *,'opinit:error in X-numbers'
        stop
      endif
c
      if (abs(imode).eq.2) goto 4567    ! only rational splines
      print *,'read partial derivative values and triangle indices...'
c
c     read # of triangles (nt), ial=2*nd+5 and triangualtion-indices (iwk)
      read(33)nti
      if(nti.gt.nt)then
        print *,'error in triangle-numbers nt= ',nt,' nti= ',nti
        stop
      endif
c
      read(33)iali
      read(33)(iwk(n),n=1,3*nti)
c
c     read partial derivative values (wk ; n=2*nd=ial-5)
      do 2500 l=1,nzvai
       do 2301 m=1,ntabi
        read(33)(wk(n,m,l),n=1,iali-1-4)
 2301  continue
 2500 continue
c
 4567 close(33)
c
c     read ival-array from file 
c      print *, 'read ival.dat for extrapolation domain checking'
      open(32,file=ivadat,status='old',err=9010)
      do 3001 l=1,nzvai
       if(nzvai.eq.13)then
c        OPAL95
         nxiri=nxir
         nyiri=nyir
         nyisi=nyis                  !=5:AlEX94<->OPAL95 fit at log(T)=3.90<->3.95
         if(iii.ne.nzva)nyisi=nyisi+2!=7:Kuru91<->OPAL95 fit at log(T)=4.00<->4.05
         mdi=md
         nyifi=nyiri+mt-nyisi+1
         read(32,'(70i3)')((ivalo(m,k,l),m=1,nyiri),k=1,ntabi)
       else
c        OPAL92
         nxiri=17
         nyiri=50
         nyisi=6
         mdi=17
         nyifi=nyiri+mt-nyisi+1
         read(32,'(50i3)')((ivalo(m,k,l),m=1,nyiri),k=1,ntabi)
       endif
 3001 continue
      close(32)
c
c     read the livermore opacity-tables
      lun =  31
      if(nzvai.eq.13)then
c        print *, 'read OPAL95 tables ....'
        call rdi95(lun,optabe,ntab,nxir,nyir,nzva,rlg,tlg,opa,ival)
      else
        print *, 'read OPAL92 tables ....'
        call readi(lun,optabe,ntab,nxir,nyir,nzva,rlg,tlg,opa,ival)
      endif
c
      if(iii.eq.nzva)then
c        print *, 'read Alexander low-temperature tables ....'
        call readal(lun+2,lowtab,ntab,nzva,nxir,ntal,rlgal,tlgal,opaal)
        iyifi= nyifi-mt+ntal                ! total nr. of log(T) values
      else
        print *, 'read Kurucz low-temperature tables ....'
        call readkx(lun+2,lowtab,ntab,nzva,md,mt,rlk,tlk,opkx)
        iyifi= nyifi                        ! total nr. of log(T) values
      endif
c
c     do we use electron conduction ?
      if(imode.gt.0)then
         lec=.true.
c         print '(a)',' electron conduction enabled'
      else
         lec=.false.
         print '(a)',' electron conduction disabled'
      endif
c
c     load arrays with opacity-values suitable for the s/r masub & rbival
c
      do 4010 l=1,nzvai
c       define xd- and yd- values (each of size nd=(nyir-nyis)*nxir)
        do 4001 itab=1,ntabi
          nd = 0
          do 1002 j=nyisi,nyiri
           jj = iyifi-nyiri+j
           do 1001 i=1,ival(j,itab,l)
c                      ^^^^^^^^^^^^^^
c                      ival is determined in 'readi' from
c                      the extrapolated tab-values (should be '17' or '19')
               nd = nd + 1
               zd(nd,itab,l) = opa(i,j,itab,l)
               xd(nd,l) = rlg(i,itab,l)
               yd(nd,l) = tlg(j,itab,l)
               ops(i,jj,itab,l) = opa(i,j,itab,l)
 1001      continue
           tls(jj) = tlg(j,itab,l)
 1002     continue
c
c         load low-temperature opacity-values
c
          if(iii.eq.nzva)then
            do js=1,ntal                        ! load Alexander tables
             do is=1,nxiri                      ! of size nd=19x19=361 ...OPAL95
              nd=nd+1
              zd(nd,itab,l)     = opaal(is,js,itab,l)
              xd(nd,l)          = rlgal(is,itab,l)
              yd(nd,l)          = tlgal(js,itab,l)
              ops(is,js,itab,l) = opaal(is,js,itab,l)
             enddo
             tls(js)=tlgal(js,itab,l)
            enddo 
          else
            do 1225 js=1,mt                     ! load Kurucz tables
             do 1221 is=1,nxiri                 ! of size nd=19x30=570 ...OPAL95
              nd=nd+1                           !         nd=17x30=510 ...OPAL92
              zd(nd,itab,l) = opkx(is,js,itab,l)
              xd(nd,l) = rlk(is)
              yd(nd,l) = tlk(js)
              ops(is,js,itab,l) = opkx(is,js,itab,l)
 1221        continue
             tls(js)=tlk(js)
 1225       continue
          endif
c
c      calculate birational spline coefficients ra if |imode|>1
c
        if (abs(imode).gt.1) then
c          if((l.eq.1).and.(itab.eq.1))
c     .              print *,'compute birational spline coefficients...'
          ir=1             ! compute derivatives at boundaries
          call rat2d(nxiri,iyifi,nxir,nyif,
     .               rlg(1,itab,l),tls,ops(1,1,itab,l),
     .               ir,pp,drelpr,
     .               p,q,r,ra(1,1,1,1,itab,l),iflag,
     .               rsdx,rsax,rsbx,rscx,rsrx,
     .               rsdy,rsay,rsby,rscy,rsry)
          if(iflag.ne.0)then
           print *,'opinit: ERROR in s/r rat2d: iflag= ',iflag
          endif
        endif
4001    continue
4010  continue      
c
c
      return
c
9010  print *,'opinit: ERROR in opening ',ivadat
      stop
9020  print *,'opinit: ERROR in opening ',pderiv
      stop
9030  print *,'opinit: error in opening ',tabnam
      stop
c
      end
