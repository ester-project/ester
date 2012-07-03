	program alex94bext
        implicit double precision (a-h,o-z)
c
c	purpose:
c		read alexander95 opacity tables extrapolate for 
c               log(R)=[-8,-7.5)  and convert values from ascii 
c               to binarie format .
c
c
c       History: 
c
c       18/08/97:  creation from s/r alex2b.f
c
c       last modification: 21/08/97
c
c	Variables:
c	lun ....... logical unit number of opacity-table-file
c	ntab ...... nr. of different x-values
c	nval ...... max. nr. of op.values pro line (= # of rlg-values)
c	nlin ...... max. nr. of op.table lines (= # of log(T)=values)
c       nzva ...... nr. of different z-values
c                   (z = mass fraction of heavy elements)
c
c	rlg ....... array[1..nval,1..ntab], 
c		    decade log of r=density(gm/cm**3)/t6**3
c	tlg ....... array[1..nlin,1..ntab].
c                   decade log of temperature
c	opa ....... array[1..nvar,1..nlin,1..ntab],
c                   opacity values
c
c       setup array-dimensions
c
        parameter(ntab=7,nzva=13)
        parameter(nval=19,nlin=23)
        parameter(lun=21)
c
        dimension       ir (     ntab,nzva)
	dimension	rlg(nval,ntab,nzva)
	dimension	tlg(nlin,ntab,nzva)
	dimension	opa(nval,nlin,ntab,nzva)
        dimension       xtab(ntab),ztab(nzva)
        dimension       xi(2),di(2)

	character*132 	line
	character*17	ifname(ntab),ofname
c
c      initialize input filenames
       data ifname 
     .     /'alex_x00_g93.tab', 'alex_x01_g93.tab','alex_x02_g93.tab',
     .      'alex_x035_g93.tab','alex_x05_g93.tab','alex_x07_g93.tab',
     .      'alex_x08_g93.tab'/
c
c-----initialize output filename
      data ofname
     .     /'alex94.bin'/
c
      np = 3               ! degree of polynomial for s/r uvip3d
      ni = 2
      xi(1) = -8.0d0       ! x-values for extrapolated points
      xi(2) = -7.5d0
c
c     open outputfile
      lun1=lun+1
      open(lun1,file=ofname,status='unknown',
     .          form='unformatted',err=9011)
c
c
c
      do 1015 k=1,ntab
c         open inputfiles
          open(lun,file=ifname(k),status='old',err=9001)
c
c---------read tables
          do 3001 l=1,nzva
c-----------read 1st line 
            read(lun,'(a)',err=9002) line
            print *,'read ',ifname(k),line(1:37)
            print *,'and extrapolate table for log(R)=[-8.,-7.5]'
            read(line,'(19x,f4.2,4x,f6.4)')xtab(k),ztab(l)
c-----------read data
            do 2001, j=1,nlin
c--------------read tlg,rlg and opacity-values
               read(lun,9999,err=9002)
     +                it,ir(k,l),tlg(it,k,l),rlg(3,k,l),
     +                (opa(i+2,it,k,l),i=1,ir(k,l))
 9999          format(i2,i3,f6.3,f5.1,
     +                8f8.3,/,10f8.3,/,3f8.3)
               if(j.ne.it)then
                     print *,' error in temperature index....'
                     stop
               endif
               ir(k,l) = ir(k,l)+2                   ! should be 19
               if(ir(k,l).ne.nval) then
                    print *,' error in log(R) index....'
                    stop
               endif
               do i=1,nval
                  rlg(i,k,l) = -8.0 + 0.5d0*(i-1)    ! fill up new log(R)-array
               enddo
               ns=5
               if(j .gt. 20) ns=2 
               call uvip3d(np,ns,rlg(3,k,l),opa(3,j,k,l),
     .                        ni,xi        ,opa(1,j,k,l),di)
 2001       continue
 3001     continue
          close(lun)
 1015   continue
c
c
c       output data in binary format
c       --------------------------------
c
        write(lun1) nzva,(ztab(l),l=nzva,1,-1)
        write(lun1) ntab,(xtab(k),k=1,ntab)
        do 4001 l=nzva,1,-1
           do 5001 k=1,ntab
              print '(a,f4.2,a,f6.4)', 
     +              ' write X= ',xtab(k),' Z= ',ztab(l)
              write(lun1)ir(k,l),(rlg(i,k,l),i=1,ir(k,l))
              write(lun1)nlin-4 ,(tlg(j,k,l),j=nlin,5,-1)!'-4' ignore 4 largest tlg
              write(lun1)((opa(i,j,k,l),i=1,ir(k,l)),j=nlin,5,-1)
 5001      continue
 4001   continue
c
        close(lun1)
	stop
c
 9001   print *,'alex94bext: error in opening ',ifname(k),' lun= ',lun
	stop
 9011   print *,'alex94bext: error in opening ',ofname,' lun= ',lun
	stop
 9002	print *,'alex94bext: error in reading: line: ',j
	stop
c
	end
