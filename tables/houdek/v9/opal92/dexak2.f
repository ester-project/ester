      program exak2
      implicit double precision (a-h,o-z)
c
c  Extrapolation of the opacity tables using akima method
c        extrapolate 1-dimensional for const. temperature 
c        (tlg>6., ....)
c
c  History:
c
c     13.2.1993: use improved Akima-Interpolation
c                (ACM TOMS 17.3 Alg. 697)
c
c      11.3.1993 introduced parameter statement
c                for array-dimensions
c
c
c  last modification: 11.3.1993
c
      parameter(ntab=7,nzva=10)
      parameter(nxi=17,nyi=50)
      parameter(ndat=1000,ndai=100)
      parameter(lun=31,nout=6)
c
      character*80 optabe,paderi,ivadat
c
c
      dimension yd(ndat),zd(ndat)
c
c
      dimension dlg(nxi,ntab,nzva)
      dimension tlg(nyi,ntab,nzva)
      dimension opa(nxi,nyi,ntab,nzva)
      dimension ival(nyi,ntab,nzva)
      dimension yval(nxi),zval(nxi)
c
      call maceps(drelpr)
c
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat)
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
      close(29)
c
      print *,'using tables:'
      print *,optabe
      print *,paderi
      print *,ivadat
c
c     read the opacity-tables and open outputfile on lun+1
      call reade(lun,optabe,ntab,nxi,nyi,nzva,dlg,tlg,opa,ival)
c
c     output of ival-array (index of last valid opacity value/line)
c     VAX/VMS needs recl=1024 (?)
c     open(lun,file=ivadat,status='UNKNOWN',recl=1024,err=9010)
      open(lun,file=ivadat,status='UNKNOWN',err=9010)
      do 101 l=1,nzva
        write(lun,'(50i3)')((ival(m,k,l),m=1,nyi),k=1,ntab)
  101 continue
      close(lun)
c
c     repeat for nzva tables (normaly: nzva=10)
      do 6010 l=1,nzva
c
c       using Akima - Interpolation (Extrapolation)
c
c       special treatment of table itab=1 (X=0.0)
c
        itab = 1
c
        write(nout,6000)l,itab
c
c       output of dlg-values
c       write(lun+1) (dlg(i,itab,l),i=1,nxi)
        write(lun+1) (dlg(i,7,l),i=1,nxi)
c        write(77,'(7X,17f7.4)') (dlg(i,7,l),i=1,nxi)
c
        do 1001 j=1,nyi
          if(ival(j,itab,l).lt.nxi) then
            nd=ival(j,itab,l)
            do 1005 i=1,nd
               yd(i) = dlg(i,itab,l)
               zd(i) = opa(i,j,itab,l)
 1005       continue
c
            call uvip3p(5,nd,yd,zd,nxi-nd,
     +                   dlg(nd+1,7,l),
     +                   opa(nd+1,j,itab,l))
          endif
c
c         output tlg + extrapolated opacity-values
c
          write(lun+1)tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
c          write(77,'(18f7.4)')tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
 1001   continue
c
c       end of special treatment of tables with X=0.0
c
c       repeat for ntab X-values (ntab=7)
        do 6001 itab=2,ntab
c
          write(nout,6000)l,itab
c
c         output of dlg-values
c
          write(lun+1) (dlg(i,itab,l),i=1,nxi)
c          write(77,'(7X,17f7.4)') (dlg(i,7,l),i=1,nxi)
c
c         first extrapolate for 6.0<=tlg<=8.0 (33<=j<=50)
c                              -0.5<=rlg<=1.0 (15<=i<=17)
c
          do 2001 j=33,nyi 
            if(ival(j,itab,l).lt.nxi) then
              nd=ival(j,itab,l)
c             nd should be 14
              do 2005 i=1,nd
                 yd(i) = dlg(i,itab,l)
                 zd(i) = opa(i,j,itab,l)
 2005         continue
c
              call uvip3p(5,nd,yd,zd,nxi-nd,
     +                    dlg(nd+1,itab,l),
     +                    opa(nd+1,j,itab,l))
            endif
 2001     continue
c
c         now interpolate for rlg =1.0 (i=17), 
c                      4.5xx<=tlg<=6.0 (j=16..32)
c
          ni=0
          nd=0
          i=17
          do 3001 j=1,15
               nd=nd+1
               zd(nd) = opa(i,j,itab,l)
               yd(nd) = tlg(j,itab,l)
 3001     continue
          do 3010 j=16,32
               ni=ni+1
               yval(j-15) = tlg(j,itab,l)
 3010     continue
          do 3020 j=33,nyi
               nd=nd+1
               zd(nd) = opa(i,j,itab,l)
               yd(nd) = tlg(j,itab,l)
 3020     continue
c
c         do interpolation
c
          call uvip3p(5,nd,yd,zd,ni,yval,zval)
c
c         fill opa-array with interpolated values zval
c
          do 3030 j=1,17
            opa(i,j+15,itab,l)=zval(j)
 3030     continue
c
c         end of interpolation for rlg =1.0 (i=17),
c                           4.5xx<=tlg<=6.0 (j=16...32)
c
c         now do interpolation for rlg=-0.5...0.5 (i=14...16)
c                           4.5xx<=tlg<=6.0       (j=16...32)
c                     
c
          do 4001 j=16,32
            if(ival(j,itab,l).lt.nxi) then
c             nd should be 13
              nd=ival(j,itab,l)
              do 4005 i=1,nd
                 yd(i) = dlg(i,itab,l)
                 zd(i) = opa(i,j,itab,l)
 4005         continue
c
c             add last point (rlg=1.0,[i=17]), to the yd/zd-Array
c
              i=17
              nd=nd+1
              yd(nd) = dlg(i,itab,l)
              zd(nd) = opa(i,j,itab,l)
c
c             do final interpolation
c
              call uvip3p(5,nd,yd,zd,nxi-nd,
     +                    dlg(nd,itab,l),
     +                    opa(nd,j,itab,l))
 
            endif
 4001     continue
c
c         output tlg + extrapolated opacity-values
c
          do 5001 j=1,nyi
            write(lun+1)tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
c            write(77,'(18f7.4)')tlg(j,itab,l),(opa(i,j,itab,l),i=1,nxi)
 5001     continue
 6001   continue
 6010 continue
c     closing of outputfile, which has been opened in 'reade'
      close(lun+1)
c
      stop
c
6000  format('extrapolation of table:',2i3)
9010  print *,'dexak2: error in opening of file ',ivadat
      stop
9030  print *,'dexak2: error in opening of OPINTPATH'
      print *,'dexak2: try to assign a logical name to OPINTPATH'
      print *,'dexak2: with setenv OPINTPATH <absolute_path_name> in '
      print *,'dexak2: a csh where <absolute_path_name> is the file '
      print *,'dexak2: containing the absolute pathnames of the '
      print *,'dexak2: input files optabe.bin, pderivs.dat, ival.dat'
      stop
c
      end
