	subroutine rdi95(lun,optabe,ntab,nval,nlin,nzva,rlg,
     +                    tlg,opa,ival)
        implicit double precision (a-h,o-z)
        character*(*) optabe
c
c	purpose:
c		read extrapolated, binary OPAL95 tables
c 	reference: 
c		rogers & iglesias, ap.j.supplement 76, april(1992)
c
c	written: g.w. houdek
c
c       History: 13.6.1992 creation
c
c                19.12.1992 modified for z-domain interpolation
c                           nzva different z-values
c                           (z = mass fraction of heavy elements)
c
c                           modified opacity values input from
c                           'formatted' to 'unformatted'
c                           first you have to run the program 'a2b'
c                           (ascii-to-binary conversion of opacity-tables)
c
c       25.1.1993 changed variable t6 -> tlg = log10(t)
c                 table-values changed to tlg in da2b
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
c                setenv OPINTPATH /dir1/dir2/opacity.pathes
c
c                the character variable opinte is defined new
c                in the list of the subroutine arguments
c
c      13/11/95: modified for OPAL95 tables
c
c      last modification: 13/11/95
c
c	input values:
c	lun ....... logical unit number of opacity-table-file
c       optabe .... absolute filename path of the opacity table
c	ntab ...... nr. of different x-values [ntab=7]
c	nval ...... max. nr. of op.values per line (const. tlg)
c	nlin ...... max. nr. of op.table lines
c       nzva ...... max. nr. of different z-values [normaly z=13]
c                   (z = mass fraction of heavy elements)
c
c       output values:
c	rlg ....... array[1..nval,1..ntab,1..nzva], 
c		    decade log of r=density(gm/cm**3)/t6**3
c	tlg ....... array[1..nlin,1..ntab,1..nzva].
c                   temperature/10**6
c	opa ....... array[1..nvar,1..nlin,1..ntab,1..nzva],
c                   opacity values
c       ival ...... array[1..nlin,1..ntab,1..nzva],
c                   contains for each read line, the nr. of the last
c                   valid (given) opacity value
c
	integer 	lun,ntab,nval,nlin,nzva
	dimension	rlg(nval,ntab,nzva)
	dimension	tlg(nlin,ntab,nzva)
	dimension	opa(nval,nlin,ntab,nzva)
        dimension       ival(nlin,ntab,nzva)
c
c     table dimension for s/r opintc{f} and opints
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyifi,mdi,
     +                nti,iali
c
c	character*10	ifname
c
c	initialize input filenames
c	data  ifname /'opal95e.bin'/
c
c	  open inputfile
	  open(lun,file=optabe,form='unformatted',
     +             status='old',err=9001)
c
c       read ntabi tables
	do 3001 k=1,ntabi
c
c         read nzva-times ntab tables
          do 1015 l=1,nzvai
c
c          read line of rlg-values
	   read(lun) (rlg(i,k,l),i=1,nxiri)
c
c	   read 'tlg' + opacity-values
           do 1010 j=1,nyiri
      	    read(lun,end=1010,err=9002)
     +           tlg(j,k,l),(opa(i,j,k,l),i=1,nxiri)
c
            do 1005 m=nxiri,1,-1
               if(opa(m,j,k,l).ne.0.0d0)goto 1007
 1005       continue
 1007       ival(j,k,l) = m
 1010      continue
 1015     continue
 3001   continue
c    	print *,'rdi95: last table ',j-1,' lines read'
c       close input file lun
        close(lun)
c
	return
c
 9001   print *,'readi: error in open statement'
	stop
 9002	print *,'readi: error in read: line: ',j
	stop
c
	end
