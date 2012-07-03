      program kur2l2g
      implicit real*8(a-h,o-z)
c
c     read new kurucz-tables and
c     convert them to the canonical
c     Livermore Z-values using univariate
c     interpolation in Z (improved Akima interpolation)
c     and interpolate density values to given 
c     rlt-Livermore-gridpoints (at constant t)
c
c     History
c
c     14.3.1993: creation
c
c      8.6.1993: changed tlg_max (fitting temperature with opal
c                tables) to 4.05
c
c     28.1.1994: changed tlg_max (fitting temperature with opal
c                tables) to 4.00
c                new array ltl4, to determine, from which table
c                the opacity value at tlg=4.0 should be used
c                in the final table
c
c     31.1.1994: for the last 4 table-points at constant density
c                spline interpolation is used to achieve smooth
c                transition between the Kurucz and Livermore-tables
c
c                     T     I    I    I    T       (T...tables point
c              tlg=  3.91 3.94 3.97 4.0 4.0416      I... interpolated point)
c
c     parameter for Livermore-tables
c
      parameter(ntab=7,nzva=10)
      parameter(nxir=17,nyir=50)
      parameter(lunl=31,ntlgc=100)
c
c     parameter for Kurucz-tables
c
      parameter(it=19,iz=10,irlt=25)
      parameter(mt=56,mp=38,iv=3,nkur=11,lkur=mp*mt)
c
c     dimension for Kurucz-tables
c
      logical ltl4(nxir)
      double precision opk(it,lkur),opl(iz,lkur)
      double precision dek(it,lkur),del(iz,lkur)
      double precision zku(it)     ,zli(iz)
      double precision rlg(mp)  ,ops(mp)
      double precision rlt(irlt),opa(irlt,mt,iz)
      double precision opi(mt+nyir),tli(mt+nyir),y2a(mt+nyir)
      double precision tlf(3),opf(3)
c
      character*12 tabnam,outtab
      character*20 infile(it)
c
c     table dimension for s/r opintc{f} and opints
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyifi,mdi,
     +                nti,iali
c
      common/kurtab/tl(mt),var(mt,mp,iv),zv,nt,np(mt)
c
c ... var(i,j,1)=log p, var(i,j,2)=log op, var(i,j,3)=log rho
c     i ... 1,...,mt
c     j ....1,...,mp
c
c
c     dimension for Livermore-tables
c
      double precision rlgl(nxir,ntab,nzva)
      double precision tlgl(nyir,ntab,nzva)
      double precision opal(nxir,nyir,ntab,nzva)
      double precision tlgc(ntlgc)
      integer ival(nyir,ntab,nzva)
      character*79 optabe
c
c     define fitting point (temperature) to opal tables
      data tlm /4.0d0/
c
      data zsun /1.886d-2/
      data tabnam /'KURTABS'/
      data outtab /'kuru91.bin'/
      data lun /21/
      data zli /0.0d0 , 1.0d-4 , 3.0d-4 , 1.0d-3 , 2.0d-3,
     +          4.0d-3, 1.0d-2 , 2.0d-2 , 3.0d-2 , 4.0d-2/
      data rlt /-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0,
     +          -4.5d0,-4.0d0,-3.5d0,-3.0d0,-2.5d0,
     +          -2.0d0,-1.5d0,-1.0d0,-0.5d0, 0.0d0,
     +           0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0,
     +           3.0d0, 3.5d0, 4.0d0, 4.5d0, 5.0d0/
      data ltl4 /.true. ,.false.,.true. ,.true. ,.true. ,
     +           .true. ,.false.,.false.,.false.,.false.,
     +           .false.,.false.,.false.,.false.,.false.,
     +           .false.,.false./
c
      nzvai=nzva
      ntabi=ntab
      nxiri=nxir
      nyiri=nyir
c
c======================================================================
c
c     read all tables listed in tabnam
c
      open(lun+2,file=tabnam,form='formatted',status='old',err=9030)
      read(lun+2,'(i3)') itab
      if(itab.gt.it)goto 9040 
      read(lun+2,'(A)') (infile(i),i=1,itab)
      close(lun+2)
      do 2010 n=1,itab
        call kurset(lun,infile(n))
        ip=0
        do 2020 i=1,mt
           do 2030 j=1,mp
              ip=ip+1
              opk(n,ip)=var(i,j,2)
              dek(n,ip)=var(i,j,3)
2030       continue
2020    continue
        if(ip.ne.lkur)goto 9050
        zku(n)=(10.d0**zv)*zsun
2010  continue
c
c======================================================================
c
c     read all Livermore opacity-tables
c
      open(29,file='OPINTPATH_92X',
     +form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      close(29)
      print *,'Using table:'
      print *,optabe
      call readi(lunl,optabe,ntab,nxir,nyir,nzva,rlgl,tlgl,opal,ival)
c
c======================================================================
c
c     do univariate interpolation in Z
c     using improved Akima Interpolation
c     ACM/Toms Alg. 697
c
c     using 3-rd degree
      nd =  3
      do 3001 i=1,lkur
c
c       interpolate the opacity values
c
        call uvip3p(nd,itab,zku,opk(1,i),iz,zli,
     +                          opl(1,i)         )
c
c       interpolate the density values
c
        call uvip3p(nd,itab,zku,dek(1,i),iz,zli,
     +                          del(1,i)         )
3001  continue
c
c======================================================================
c
c     search index for the maximum temperature
c     tlm=tlg_max=4.0 (fitting point to Livermore-tabs)
c     evaluate Z-domain (mass fraction of heavy element):
c     search in z-table ztab (bisection-method)
c
      jlo=1
      call hunt(tl,mt,tlm,jlo)
      print *,'index at fitting point tlm= ',tlm
      print *,'upper index of Kurucz-tables ',
     +        'extrapoint jlo= ',jlo
c
c     search index for minimum temperature of
c     Livermore-tables, to fit at tlm
c
      llo=1
      call hunt(tlgl(1,1,1),nyir,tlm,llo)
      llo=llo+1
      print *,'lower index of Livermore-tables llo= ',llo
      print *,'tlgl_min= ',tlgl(llo,1,1)
c
c======================================================================
c
c
c     fill up a continuous temperature array tlgc (Kurucz+Livermore)
c
      do 4001 i=1,jlo
         tlgc(i)=tl(i)
4001  continue
      do 4002 i=llo,nyir
         tlgc(jlo+i-llo+1)=tlgl(i,1,1)
4002  continue
      itlgc=jlo+nyir-llo+1
      if(itlgc.ge.ntlgc)then
        print *,'array dimension error of tlgc: ',itlgc
      endif
c
c======================================================================
c
c     interpolate the density values at constant t 
c     (tlg_min --> tlg_max) at the
c     irlt given gridpoints rlt and set last point
c     at tlg=4.0 to either Kurucz- or Livermore-value
c     according to ltl4 (.true. = Livermore-value)
c
c     rlg = log10(density/t6**3)
c     rlg = del -3*rl +18
c
c     setup  5 intermediate temperature interpolation points
c     opacities are defined later and 'if' fitting points, to be interpolated
      if=3
      do 5003 i=1,jlo-if
         tli(i)=tl(i)
5003  continue
      do 5004 i=1,nyir-llo+1
         tli(jlo-if+i)=tlgl(llo+i-1,1,1)
5004  continue
c     setup fitting temperature table points
      tlf(1)=tl(jlo-2)
      tlf(2)=tl(jlo-1)
      tlf(3)=tl(jlo)
c
      do 5001 l=1,iz
c
c        for use with idl-precedure kurucz.pro
c        do 5010 i=1,mt
c
         do 5010 i=1,jlo
           do 5020 j=1,mp
              rlg(j)=del(l,(i-1)*mp+j)-3.d0*tl(i)+18.d0
              ops(j)=opl(l,(i-1)*mp+j)
5020       continue
           call uvip3p(nd,mp,rlg,ops,irlt,rlt,
     +                           opa(1,i,l)   )
5010     continue
c        set table point at tlg(jlo)=4.0 to Livermore-value
         do 5005 m=1,nxir
            if(ltl4(m).eqv..true.) opa(m,jlo,l)=opal(m,llo-1,6,l) 
c                                                           ^
c                                                           X=0.7
c           interpolate intermediate values for 
c           last 4 table points at constant density
            if(m.ge.(nxir-3)) then
               do 5013 n=1,jlo-if
                  opi(n)=opa(m,n,l)
5013           continue
               do 5015 n=1,nyir-llo+1
                  opi(jlo-if+n)=opal(m,llo+n-1,6,l)
5015           continue
               ni=jlo-if+nyir-llo+1
               call spline(tli,opi,ni,opi(1),opi(ni),y2a)
               do 5017 n=1,if
                  call splint(tli,opi,y2a,ni,tlf(n),opf(n),ydx)
                  opa(m,jlo-if+n,l)=opf(n)
5017           continue
            endif
c
5005     continue
5001  continue
c
c======================================================================
c
c     now output the new interpolated values
c     only from tl_min -> tl_max + inserted extrapoint
c
      print*,'*********************************'
      print*,'output tables to file ',outtab
      print*,'*********************************'
c
      open(lun+1,file=outtab,status='unknown',
     +           form='unformatted')
c
      write(lun+1) iz,irlt,jlo
      write(lun+1) (zli(i),i=1,iz)
      write(lun+1) (tl(i),i=1,jlo+1)
      do 6001 l=1,iz
         write(lun+1) ((rlt(j),opa(j,i,l),j=1,irlt),i=1,jlo)
6001  continue
      close(lun+1)
c
c======================================================================
c
c     just for testing with idl (kurucz.pro)
c     output one interpolated table
c     
      l=8
      print *,'output test-table Z=',zli(l)
      open(8,form='unformatted')
      open(9,form='unformatted')
      open(10,form='unformatted')
      open(80,form='unformatted')
      open(81,form='unformatted')
      write(8)irlt,(rlt(j),j=1,irlt)
      write(9)mt,(tl(i),i=1,mt)
      write(10)nyir,(tlgl(i,6,l),i=1,nyir)
c
      write(80)((opa(j,i,l),j=1,irlt),i=1,mt)
      write(81)((opal(j,i,6,l),j=1,nxir),i=1,nyir)
c
      close(8)
      close(9)
      close(10)
      close(80)
      close(81)
c
c======================================================================
c
      stop
c
9030  print *,'kur2l2g: error in opening ',tabnam
      stop
9040  print *,'kur2l2g: array-dim error of it'
      stop
9050  print *,'kue2l2g: wrong nr. of tablepoints',ip,lkur
      stop
c
      end
