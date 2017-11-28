	subroutine rde95(lun,optabe,ntab,nval,nlin,nzva,rlg,
     +                    tlg,opa,ival)
        implicit double precision (a-h,o-z)
        character*(*) optabe
c
c	purpose:
c		read the binary OPAL95 values
c 	reference: 
c		rogers & iglesias, ap.j.supplement 76, april(1992)
c
c	written: g.w. houdek
c
c       History: 
c
c       13/11/95: creation, derived from dreade.f
c
c	last modification: 13/11/95
c
c	input values:
c	lun ....... logical unit number of opacity-table-file
c       optabe .... absolute filename path of the (extrap.) opacity table
c	ntab ...... nr. of opacity tables
c	nval ...... max. nr. of op.values pro line (const. tlg)
c	nlin ...... max. nr. of op.table lines
c       nzva ...... max. nr. of different z-values [z=13 for OPAL95]
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

	character*10	ifname
        character*11    ofname
c
c	initialize input filename
	data ifname /'opal95.bin'/
c
c
c	  open inputfile/outputfile
	  open(lun,  file=ifname,status='old',
     +               form='unformatted',    err=9001)
	  open(lun+1,file=optabe,status='unknown',
     +               form='unformatted',    err=9011)
c
c       read ntab tables
	do 3001 k=1,ntab
c
c         read nzva-times ntab tables
	  do 1015 l=1,nzva
c
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
 1015	  continue
 3001   continue
        print *,'rde95: last table ',j-1,' lines read'
c       close input file lun
        close(lun)
c       close of output file lun+1 in calling routine (exop95.f)
c
	return
c
 9001   print *,'rde95: error in opening ',ifname,'lun= ',lun
	stop
 9011   print *,'rde95: error in opening ',ofname,'lun= ',lun+1
	stop
 9002	print *,'rde95: error in read: line: ',j
	stop
c
	end
