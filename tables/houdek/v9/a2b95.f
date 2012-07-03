	program a2b95
        implicit double precision (a-h,o-z)
c
c	purpose:
c		read the OPAL95 of the
c		livermoor opacity tables and convert 
c               the values from ascii to binarie format .
c 	reference: 
c		rogers & iglesias, ap.j.supplement 79, april(1992)
c
c	written: g.w. houdek
c
c       History: 
c
c      13.11.1995 creation, derived from da2b.f
c
c	last modification: 13.11.1995
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
        parameter(ntab=7,nzva=13)
        parameter(nval=19,nlin=70)
        parameter(lun=21)
c
	dimension	rlg(nval,nzva)
	dimension	tlg(nlin,nzva)
	dimension	opa(nval,nlin,nzva)

	character*132 	line
	character*10	ifname,ofname
	character*6	ct6
c
c	initialize input filename
	data   ifname
     +         /'opal95.tab'/
c
c	initialize output filename
	data	ofname
     +          /'opal95.bin'/
c
c	  open outputfile
	  open(lun+1,file=ofname,status='unknown',
     +               form='unformatted',err=9011)
c
c
c
c	open inputfile
	open(lun,file=ifname,status='old',err=9001)
c
c       read ntab tables/Z-value
        do 3001 k=1,ntab
c
	do 1015 l=1,nzva
c	  read 1st 5 lines 
	  do 1001 i=1,5
	   read(lun,'(a)') line
           if (i.eq.2) print *,line(1:60)
 1001     continue	
c
c	  read/write string 'logT' + rlg-values
	  read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k),i=1,nval)
	  write(lun+1              )     (rlg(i,k),i=1,nval)
c
c	  read blank line
	  read(lun,'(a)') line
c
c	  read/write logT & opacity-values
          do 1010 j=1,nlin
      	    read(lun,'(f5.3,19f7.3)',err=9002)
     +           tlg(j,k)           ,(opa(i,j,k),i=1,nval)
     	    write(lun+1,             err=9020)
     +           tlg(j,k)           ,(opa(i,j,k),i=1,nval)
 1010     continue
c     	  print *,'a2b95: ',j-1,' lines read in file ',ifname
 1015     continue
 3001   continue
	close(lun+1)
c
	stop
c
 9001   print *,'a2b95: error in opening ',ifname,'lun= ',lun
	stop
 9011   print *,'a2b95: error in opening ',ofname,'lun= ',lun+1
	stop
 9002	print *,'a2b95: error in reading: line: ',j
	stop
 9020   print *,'a2b95: error in writing ',ofname,'lun= ',lun+1 
        stop
c
	end
