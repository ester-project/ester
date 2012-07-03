      program exakop95
      implicit double precision (a-h,o-z)
c
c  Extrapolate OPAL95 tables using Akima's method
c
c  History:
c
c     16/08/97: created from exop95.f
c
c  last modification: 16/08/97
c
      parameter(ntab=7,nzva=13)
      parameter(nxi=19,nyi=70)
      parameter(ndat=100)
      parameter(lun=31,nout=6)
c
      character*80 optabe,paderi,ivadat
c
c
      dimension xd(ndat),zd(ndat)
c
c
      dimension rlg(nxi,ntab,nzva)
      dimension tlg(nyi,ntab,nzva)
      dimension opa(nxi,nyi,ntab,nzva)
      dimension ival(nyi,ntab,nzva)
c
      dimension ireg1(3),ireg2(3),ireg3(5),ireg4(5),ireg5(5),ireg6(5)
      dimension itar1(4),itar2(1),itar3(3),itar4(7),itar5(9),itar6(10)
      dimension xi(10),di(10)
c
c     define indices of table-points, used for the input-data
c     of the extrapolation s/r uvip3d
c
      data nreg1,ireg1
     +     /3,
     +      13,14,15/
c
      data nreg2,ireg2
     +     /3,
     +      67,68,69/
c
      data nreg3,ireg3
     +     /5,
     +      63,64,67,68,69/
c
      data nreg4,ireg4
     +     /5,
     +      59,60,67,68,69/
c
      data nreg5,ireg5
     +     /5,
     +      57,58,67,68,69/
c
      data nreg6,ireg6
     +     /5,
     +      56,57,67,68,69/
c
c     target indices of extrapolated data
c
      data ntar1,itar1
     +     /4,
     +      16,17,18,19/
c
      data ntar2,itar2
     +     /1,
     +      70/
c
      data ntar3,itar3
     +     /3,
     +      65,66,70/
c
      data ntar4,itar4
     +     /7,
     +      61,62,63,64,65,66,70/
c
      data ntar5,itar5
     +     /9,
     +      59,60,61,62,63,64,65,66,70/
c
      data ntar6,itar6
     +     /10,
     +      58,59,60,61,62,63,64,65,66,70/
c
c 
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat)
c
      open(29,file='OPINTPATH_95X',
     +     form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      read(29,'(A)')paderi
      read(29,'(A)')ivadat
      close(29)
c
c     read the opacity-tables and open outputfile on lun+1
c
      call rde95(lun,optabe,ntab,nxi,nyi,nzva,rlg,tlg,opa,ival)
c
c     output of ival-array (index of last valid opacity value/line)
c
      open(lun,file=ivadat,status='UNKNOWN',err=9010)
      do 101 l=1,nzva
        write(lun,'(70i3)')((ival(m,k,l),m=1,nyi),k=1,ntab)
  101 continue
      close(lun)
c
      np = 3                ! degree of polynomial
      ni = 1                ! used for interp. of one point only
c
c     repeat ntab tables (ntab=7)
      do 6020 itab=1,ntab
c       repeat for nzva tables (nzva=13 for OPAL95)
        do 6010 l=1,nzva
c
c       inter/extrapolation is performed with Akima's univariate method
c
        write(nout,6000)l,itab
c
c       region 1 (constant T)
c       ---------------------
        j=67                                     ! log(T)=8.1
        do i=1,nreg1
           xd(i) = rlg(ireg1(i),  itab,l)
           zd(i) = opa(ireg1(i),j,itab,l)
        enddo
        do i=1,ntar1
           xi(i) = rlg(itar1(i),  itab,l)
        enddo
        call uvip3d(np,nreg1,xd,zd,
     .                 ntar1,xi,opa(itar1(1),j,itab,l),di)
        do j=68,69                               ! log(T)=[8.3,8.5,8.7]
           do i=2,nreg1
              xd(i-1) = rlg(ireg1(i),  itab,l)
              zd(i-1) = opa(ireg1(i),j,itab,l)
           enddo
           call uvip3d(np,nreg1-1,xd,zd,
     .                    ntar1,xi,opa(itar1(1),j,itab,l),di)
        enddo
c
c       region 2 (constant R)
c       ---------------------
        j=15                                     ! log(R) = -1.
        do i=1,nreg2
           xd(i) = tlg(  ireg2(i),itab,l)
           zd(i) = opa(j,ireg2(i),itab,l)
        enddo
        xi(ntar2) = tlg(  itar2(ntar2),itab,l)   ! only one,
        call uvip3d(np,nreg2,xd,zd,              ! last (max T) point
     .  ni,xi(ntar2),opa(j,itar2(ntar2),itab,l),di)
c
c       region 3 (constant R)
c       ---------------------
        j=16                                     ! log(R) = -0.5
        do i=1,nreg3
           xd(i) = tlg(   ireg3(i),itab,l)
           zd(i) = opa(j, ireg3(i),itab,l)
        enddo
        do i=1,ntar3
           xi(i) = tlg(   itar3(i),itab,l)
           call uvip3d(np,nreg3,xd,zd,
     .     ni,xi(i),opa(j,itar3(i),itab,l),di)
        enddo
c
c       region 4 (constant R)
c       ---------------------
        j=17                                     ! log(R) = 0.0
        do i=1,nreg4
           xd(i) = tlg(   ireg4(i),itab,l)
           zd(i) = opa(j, ireg4(i),itab,l)
        enddo
        do i=1,ntar4
           xi(i) = tlg(   itar4(i),itab,l)
           call uvip3d(np,nreg4,xd,zd,
     .     ni,xi(i),opa(j,itar4(i),itab,l),di)
        enddo
c
c       region 5 (constant R)
c       ---------------------
        j=18                                     ! log(R) = 0.5
        do i=1,nreg5
           xd(i) = tlg(   ireg5(i),itab,l)
           zd(i) = opa(j, ireg5(i),itab,l)
        enddo
        do i=1,ntar5
           xi(i) = tlg(   itar5(i),itab,l)
           call uvip3d(np,nreg5,xd,zd,
     .     ni,xi(i),opa(j,itar5(i),itab,l),di)
        enddo
c
c       region 6 (contant R)
c       --------------------
        j=19                                     ! log(R) = 1.0
        do i=1,nreg6
           xd(i) = tlg(   ireg6(i),itab,l)
           zd(i) = opa(j, ireg6(i),itab,l)
        enddo
        do i=1,ntar6
           xi(i) = tlg(   itar6(i),itab,l)
           call uvip3d(np,nreg6,xd,zd,
     .     ni,xi(i),opa(j,itar6(i),itab,l),di)
        enddo
c
c       output of rlg-values
c       --------------------
c
        write(lun+1) (rlg(i,itab,l),i=1,nxi)
c       write(77,'(/,7X,19f7.3,/)') (rlg(i,7,l),i=1,nxi)
c
c       output tlg + extrapolated opacity-values
c
        do j=1,nyi
           write(lun+1)tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
c          write(77,'(20f7.3)')tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
        enddo
c
 6010   continue
 6020 continue
c
c     closing of outputfile, which has been opened in 'reade'
      close(lun+1)
c
      stop
c
6000  format('extrapolation of table:',2i3)
9010  print *,'exop95: error in opening of file ',ivadat
      stop
9030  print *,'exop95: error in opening of OPINTPATH_95'
      stop
c
      end
