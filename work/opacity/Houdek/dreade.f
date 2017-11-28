	subroutine reade (lun,optabe,ntab,nval,nlin,nzva,rlg,
     +                    tlg,opa,ival)
        implicit double precision (a-h,o-z)
        character*(*) optabe
c
c	purpose:
c		read the opacity values of the
c		livermoor opacity tables.
c 	reference: 
c		rogers & iglesias, ap.j.supplement 76, april(1992)
c
c	written: g.w. houdek
c
c       History: 13.6.1992 creation
c
c		 30.9.1992 inserted 'status='unknown'' for VAX/VMS
c
c     19.12.1992 modified for z-domain interpolation
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
c       7.2.1993 modified pathname-treatment of inputfiles.
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
c                the character variable opinte is defined new
c                in the list of the subroutine arguments
c
c      7.2.1993  modified array-dimension-definitions to
c                fixed numbers, because s/r-arguments have to
c                be declared identical as in the calling routine.
c
c     11.3.1993  modified array-dimension-definitions again to
c                argument-values. In the calling routine
c                the parameter statement has benn introduced 
c    
c	last modification: 8.4.1993
c
c	input values:
c	lun ....... logical unit number of opacity-table-file
c       optabe .... absolute filename path of the (extrap.) opacity table
c	ntab ...... nr. of opacity tables
c	nval ...... max. nr. of op.values pro line (const. tlg)
c	nlin ...... max. nr. of op.table lines
c       nzva ...... max. nr. of different z-values [normaly z=10]
c                   (z = mass fraction of heavy elements)
c
c       output values:
c	rlg ....... array[1..nval,1..ntab,1..nzva], 
c		    decade log of r=density(gm/cm**3)/t6**3
c	tlg ....... array[1..nlin,1..ntab,1..nzva].
c                   log10(temperature)
c	opa ....... array[1..nval,1..nlin,1..ntab,1..nzva],
c                   log10(opacity) values
c       ival ...... array[1..nlin,1..ntab,1..nzva],
c                   contains for each read line, the nr. of the last
c                   valid (given) opacity value for extrapolation domain
c                   cheking (if enabled)
c
	integer 	lun,ntab,nval,nlin,nzva
	dimension	rlg(nval,ntab,nzva)
	dimension	tlg(nlin,ntab,nzva)
	dimension	opa(nval,nlin,ntab,nzva)
        dimension       ival(nlin,ntab,nzva)

	character*10	ifname,ofname
c
c	initialize input filename
	data ifname /'opal92.bin'/
c
c
c	  open inputfile/outputfile
	  open(lun,  file=ifname,status='old',
     +               form='unformatted',    err=9001)
	  open(lun+1,file=optabe,status='unknown',
     +               form='unformatted',    err=9011)
c
c
c       read nzva-times ntab tables
	do 1015 l=1,nzva
c
c	  read ntab tables
	  do 3001 k=1,ntab
c
c	   read line of rlg-values
	   read(lun) (rlg(i,k,l),i=1,nval)
c
c	   read 'tlg' + opacity-values
           do 1010 j=1,nlin
      	    read(lun,end=1010,err=9002) 
     +           tlg(j,k,l),(opa(i,j,k,l),i=1,nval)
c
            do 1005 m=nval,1,-1
               if(opa(m,j,k,l).ne.0.0d0)goto 1007
 1005       continue
 1007       ival(j,k,l) = m
 1010      continue
 3001	  continue
 1015   continue
        print *,'reade: last table ',j-1,' lines read'
c       close input file lun
        close(lun)
c       close of output file lun+1 in calling routine (dextra.f)
c
	return
c
 9001   print *,'reade: error in opening ',ifname,'lun= ',lun
	stop
 9011   print *,'reade: error in opening ',ofname,'lun= ',lun+1
	stop
 9002	print *,'reade: error in read: line: ',j
	stop
c
	end
