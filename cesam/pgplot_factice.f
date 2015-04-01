
c***********************************************************

	SUBROUTINE pgbegin(n,device,n1,n2)
	 INTEGER, INTENT(in) :: n, n1, n2
	 CHARACTER (len=10), INTENT(in) :: device
	 RETURN
	END SUBROUTINE pgbegin	

c***********************************************************

	SUBROUTINE pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 REAL, INTENT(in) :: xtick,ytick		
	 INTEGER, INTENT(in) :: nxsub,nysub
	 CHARACTER (len=*), INTENT(in) :: xopt,yopt	
	 RETURN
	END SUBROUTINE pgbox 	

c***********************************************************

	SUBROUTINE pgdraw(x,y)
	 REAL, INTENT(in) :: x,y		
	 RETURN
	END SUBROUTINE pgdraw

c***********************************************************
	
	SUBROUTINE pgend
	 RETURN
	END SUBROUTINE pgend

c***********************************************************

	SUBROUTINE pgenv(xmin,xmax,ymin,ymax,n1,n2)
	 REAL, INTENT(in) :: xmax, xmin, ymax, ymin
	 INTEGER :: n1, n2
	 RETURN
	END SUBROUTINE pgenv

c***********************************************************

	SUBROUTINE pglabel(text,t,titre)
	 CHARACTER (len=*), INTENT(in) :: text,t,titre	
	 RETURN
	END SUBROUTINE pglabel		
	
c***********************************************************

	SUBROUTINE  pgline(n,x,y)
	 REAL, DIMENSION(:), INTENT(in) :: x, y 
	 INTEGER, INTENT(in) :: n
	 RETURN
	END SUBROUTINE pgline

c***********************************************************

	SUBROUTINE pgmove(x,y)
	 REAL, INTENT(in) :: x,y		
	 RETURN
	END SUBROUTINE pgmove
	
c***********************************************************

	SUBROUTINE pgmtext(side,disp,coord,fjust,text) 
	 REAL, INTENT(in) :: disp,coord,fjust		
	 CHARACTER (len=*), INTENT(in) :: side,text
	 RETURN
	END SUBROUTINE pgmtext	
	
c***********************************************************

	SUBROUTINE pgpoint(n,x,y,m)
	 REAL, INTENT(in), DIMENSION(:) :: x,y		
	 INTEGER, INTENT(in) :: n,m
	 RETURN
	END SUBROUTINE pgpoint	
	
c***********************************************************

	SUBROUTINE pgpoly(n,x,y)
	 REAL, INTENT(in), DIMENSION(:) :: x,y		
	 INTEGER, INTENT(in) :: n
	 RETURN
	END SUBROUTINE pgpoly	

c***********************************************************

	SUBROUTINE pgscf(n)
	 INTEGER, INTENT(in) :: n
	 RETURN
	END SUBROUTINE pgscf

c***********************************************************

	SUBROUTINE pgsch(s)
	 REAL, INTENT(in) :: s		
	 RETURN
	END SUBROUTINE pgsch	
	
c***********************************************************

	SUBROUTINE pgsci(i)
	 INTEGER, INTENT(in) ::i		
	 RETURN
	END SUBROUTINE pgsci
		
c***********************************************************

	SUBROUTINE pgscr(i,x,y,z)
	 REAL, INTENT(in) ::x,y,z 		
	 INTEGER, INTENT(in) ::i 
	 RETURN
	END SUBROUTINE pgscr	

c***********************************************************

	SUBROUTINE pgsfs(i) 
	 INTEGER, INTENT(in) :: i
	 RETURN
	END SUBROUTINE pgsfs
	
c***********************************************************

	SUBROUTINE pgsls(i)		
	 INTEGER, INTENT(in) ::i 
	 RETURN
	END SUBROUTINE pgsls
	
c***********************************************************

	SUBROUTINE pgslw(i)
	 INTEGER, INTENT(in) :: i
	 RETURN
	END SUBROUTINE pgslw
	
c***********************************************************

	SUBROUTINE pgtext(x,y,t)
	 REAL, INTENT(in) :: x,y 		
	 CHARACTER (len=*), INTENT(in) :: t
	 RETURN
	END SUBROUTINE 	pgtext
	
c***********************************************************

	SUBROUTINE pgvsize(x,x1,y,y1)
	 REAL, INTENT(in) :: x,x1,y,y1		
	 RETURN
	END SUBROUTINE 	pgvsize
	
c***********************************************************

	SUBROUTINE pgwindow(xmin,xmax,ymin,ymax) 
	 REAL, INTENT(in) :: xmin,xmax,ymin,ymax		
	 RETURN
	END SUBROUTINE pgwindow
	
c***********************************************************








