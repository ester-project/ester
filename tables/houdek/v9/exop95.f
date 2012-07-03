      program exop95
      implicit double precision (a-h,o-z)
c
c  Extrapolate of OPAL95 tables using Shepard's method
c
c  History:
c
c     13/11/95: created from dexak2.f             
c
c  last modification: 13/11/95
c
      parameter(ntab=7,nzva=13)
      parameter(nxi=19,nyi=70)
      parameter(ndat=100)
      parameter(lun=31,nout=6)
c
      character*80 optabe,paderi,ivadat
c
c
      dimension xd(ndat),yd(ndat),zd(ndat)
c
c
      dimension rlg(nxi,ntab,nzva)
      dimension tlg(nyi,ntab,nzva)
      dimension opa(nxi,nyi,ntab,nzva)
      dimension ival(nyi,ntab,nzva)
c
      dimension ireg1( 5,2),ireg2( 8,2),ireg3(11,2)
      dimension itar1( 5,2),itar2(12,2),itar3(25,2)
c
c     define inideces of table-points, used for the input-data
c     of the extrapolation s/r shepard
c
      data nreg1,ireg1
     +     /5,
     +      19,18,17,17,17,
     +      57,58,58,59,60/
c
      data nreg2,ireg2
     +     /8,
     +      19,18,17,16,16,16,16,16,
     +      60,60,60,60,61,62,63,64/
c
      data nreg3,ireg3
     +     /11,
     +      19,18,17,16,15,15,15,15,15,15,14,
     +      64,64,64,64,64,65,66,67,68,69,70/
c
c     target indeces of extrapolated data
c
      data ntar1,itar1
     +     /5,
     +      19,19,18,19,18,
     +      58,59,59,60,60/
c
      data ntar2,itar2
     +     /12,
     +      19,18,17,19,18,17,19,18,17,19,18,17,
     +      61,61,61,62,62,62,63,63,63,64,64,64/
c
      data ntar3,itar3
     +     /25,
     +      19,18,17,16,19,18,17,16,19,18,17,16,19,18,17,16,
     +      19,18,17,16,19,18,17,16,15,
     +      65,65,65,65,66,66,66,66,67,67,67,67,68,68,68,68,
     +      69,69,69,69,70,70,70,70,70/
c
c 
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat)
c     first read the pathname of the file, which declares
c     the the absolute pathnames of the above mentioned input files
c     you can assign a logical name to OPINTPATH; eg. in a csh:
c     setenv OPINTPATH <absolute_path_name>
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
c     repeat ntab tables (ntab=7)
      do 6020 itab=1,ntab
c       repeat for nzva tables (nzva=13 for OPAL95)
        do 6010 l=1,nzva
c
c       Extrapolation is performed with Shepard's method
c
        write(nout,6000)l,itab
c
c       extrapolate region 1
c
        do 1005 i=1,nreg1
            xd(i) = rlg(ireg1(i,1),itab,l)
            yd(i) = tlg(ireg1(i,2),itab,l)
            zd(i) = opa(ireg1(i,1),ireg1(i,2),itab,l)
 1005   continue
c
        do 1010 i=1,ntar1
           opa(itar1(i,1),itar1(i,2),itab,l)=
     +     shep(nreg1,xd,yd,zd,
     +          rlg(itar1(i,1),itab,l),tlg(itar1(i,2),itab,l))
 1010   continue    
c
c       extrapolate region 2
c
        do 1015 i=1,nreg2
            xd(i) = rlg(ireg2(i,1),itab,l)
            yd(i) = tlg(ireg2(i,2),itab,l)
            zd(i) = opa(ireg2(i,1),ireg2(i,2),itab,l)
 1015   continue
c
        do 1020 i=1,ntar2
           opa(itar2(i,1),itar2(i,2),itab,l)=
     +     shep(nreg2,xd,yd,zd,
     +          rlg(itar2(i,1),itab,l),tlg(itar2(i,2),itab,l))
 1020   continue    
c
c       extrapolate region 3
c
        do 1025 i=1,nreg3
            xd(i) = rlg(ireg3(i,1),itab,l)
            yd(i) = tlg(ireg3(i,2),itab,l)
            zd(i) = opa(ireg3(i,1),ireg3(i,2),itab,l)
 1025   continue
c
        do 1030 i=1,ntar3
           opa(itar3(i,1),itar3(i,2),itab,l)=
     +     shep(nreg1,xd,yd,zd,
     +          rlg(itar3(i,1),itab,l),tlg(itar3(i,2),itab,l))
 1030   continue    
c
c       output of rlg-values
c
        write(lun+1) (rlg(i,itab,l),i=1,nxi)
c       write(77,'(/,7X,19f7.3,/)') (rlg(i,7,l),i=1,nxi)
c
c       output tlg + extrapolated opacity-values
c
        do 1100,j=1,nyi
          write(lun+1)tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
c         write(77,'(20f7.3)')tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
 1100   continue
c
 6010   continue
 6020 continue
c     closing of outputfile, which has been opened in 'reade'
      close(lun+1)
c
      stop
c
6000  format('extrapolation of table:',2i3)
9010  print *,'exop95: error in opening of file ',ivadat
      stop
9030  print *,'exop95: error in opening of OPINTPATH_95'
c
      stop
c
      end
