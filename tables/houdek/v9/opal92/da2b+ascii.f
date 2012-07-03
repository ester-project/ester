	program da2bpa
        implicit double precision (a-h,o-z)
c
c	purpose:
c		read the opacity values of the
c		livermoor opacity tables and convert 
c               the values from ascii to binarie format .
c 	reference: 
c		rogers & iglesias, ap.j.supplement 79, april(1992)
c
c	written: g.w. houdek
c
c       History: 
c
c      13.12.1992 creation
c
c      25.1.1993: changed binary-output of 
c                           t6 -> log10(t) = 6.+log10(t6)
c
c      11.3.1993 introduced parameter statement
c                for array-dimensions
c
c      18.8.1993: added asci-output of tlg values
c                
c	last modification: 18.8.1992
c
c	Variables:
c	lun ....... logical unit number of opacity-table-file
c	ntab ...... nr. of opacity tables/file
c	nval ...... max. nr. of op.values pro line (const. t6)
c	nlin ...... max. nr. of op.table lines
c       nzva ...... Nr. of input files (different z-values)
c                   (z = mass fraction of heavy elements)
c
c	rlg ....... array[1..nval,1..ntab], 
c		    decade log of r=density(gm/cm**3)/t6**3
c	t6 ........ array[1..nlin,1..ntab].
c                   temperature/10**6
c	opa ....... array[1..nvar,1..nlin,1..ntab],
c                   opacity values
c
c       setup array-dimensions
c
        parameter(ntab=7,nzva=10)
        parameter(nval=17,nlin=50)
        parameter(lun=21)
c
	dimension	rlg(nval,nzva)
	dimension	t6 (nlin,nzva)
	dimension	opa(nval,nlin,nzva)

	character*132 	line
	character*11	ifname(10),ofname
	character*7	ct6
c
c	initialize input filenames
	data   ifname(1)    , ifname(2)    , ifname(3),
     +         ifname(4)    , ifname(5)    , ifname(6),
     +         ifname(7)    , ifname(8)    , ifname(9),
     +         ifname(10)
     +       /'z0000.tab','z0001.tab','z0003.tab',
     +        'z0010.tab','z0020.tab','z0040.tab',
     +        'z0100.tab','z0200.tab','z0300.tab',
     +        'z0400.tab'/
c
c	initialize output filename
	data	ofname
     +         /'optab.bin'/
c
c	  open outputfile
	  open(lun+1,file=ofname,status='unknown',
     +               form='unformatted',err=9011)
c
c
c
	do 1015 l=1,nzva
c	  open inputfiles
	  open(lun,file=ifname(l),status='old',err=9001)
c
c         read ntab tables/input-file
	  do 3001 k=1,ntab
c	  read 1st 6 lines 
	  do 1001 i=1,6
	   read(lun,'(a)') line
           if (i.eq.3) print *,ifname(l),line(44:62)
 1001     continue	
c
c	  read/write string 'T6' + rlg-values
	  read(lun,   '(a7,17f7.3)') ct6,(rlg(i,k),i=1,nval)
	  write(lun+1              )     (rlg(i,k),i=1,nval)
c
c	  read blank line
	  read(lun,'(a)') line
c
c	  read/write t6/log10(t) & opacity-values
          do 1010 j=1,nlin
      	    read(lun,'(f7.3,17f7.3)',err=9002)
     +           t6 (j,k)           ,(opa(i,j,k),i=1,nval)
     	    write(lun+1,             err=9020)
     +           6.d0+log10(t6(j,k)),(opa(i,j,k),i=1,nval)
            if(k.eq.1.and.l.eq.1)write(99,'(f12.8)',err=9020) 
     +                           6.d0+log10(t6(j,k))
 1010     continue
c     	  print *,'da2b: ',j-1,' lines read in file ',ifname(k)
 3001     continue
	  close(lun)
 1015   continue
	close(lun+1)
c
	stop
c
 9001   print *,'da2b: error in opening ',ifname(k),'lun= ',lun
	stop
 9011   print *,'da2b: error in opening ',ofname,'lun= ',lun+1
	stop
 9002	print *,'da2b: error in reading: line: ',j
	stop
 9020   print *,'da2b: error in writing ',ofname,'lun= ',lun+1 
        stop
c
	end
