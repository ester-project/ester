      program opdalex94
      implicit double precision (a-h,o-z)
c
c    computes partial derivatives for OPAL95+ALEX94
c    as used for the minimum norm interpolation scheme
c
c    History:
c    19/08/97: creation from dopxext_s_95.f
c
c    Last modification:
c    21/08/97
c
      parameter(ntab=7,nzva=13)
      parameter(nxir=19,nyir=70,nyis=5)
      parameter(ntal=19)
      parameter(ndat=nxir*(nyir-nyis+1+ntal))
      parameter(niwk=70000,nwk=60000)
      parameter(lun=31,nout=6)
c
c
      character*79 optabe,paderi,ivadat,lowtab
c
c     for s/r masube
      dimension xd(ndat),yd(ndat),zd(ndat),
     +          iwk(niwk),wk(nwk)
c
c     dimension for OPAL95-tables
c
      dimension rlg(nxir,ntab,nzva)
      dimension tlg(nyir,ntab,nzva)
      dimension opa(nxir,nyir,ntab,nzva)
      dimension ival(nyir,ntab,nzva)
c
c     dimension for Alexander-tables
c
      dimension opaal(nxir,ntal,ntab,nzva)
      dimension rlgal(nxir,ntab,nzva)
      dimension tlgal(ntal,ntab,nzva)
c
c     pointers for the interpolation routine masube
c     nt = # of triangles
      common /jpoint/ nt,jipt,jipl,jiwp,jind,jngp,jigp
      common /ipoint/ ipd,ial,ibe,iga,iei,ials,ibes,igas,izx,izy
c
c     table dimension for s/r rdi95
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyif,mdi,
     .                nti,iali
c
      nzvai=nzva      ! define values for s/r rdi95
      ntabi=ntab
      nxiri=nxir
      nyiri=nyir
      nyisi=nyis
c
      call maceps(drelpr)
c
c     read the absolute pathnames of the tables
c
      open(29,file='OPINTPATH_AX',
     +form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      read(29,'(A)')lowtab
      read(29,'(A)')paderi
      read(29,'(A)')ivadat
      close(29)
      print *,'Using tables:'
      print *,optabe
      print *,lowtab
      print *,paderi
      print *,ivadat
c
      tp=0.d0                 ! tension parameter for s/r masube
c
c     read the opacity-tables (OPAL95 + ALEX94)
      print *,'read OPAL95 tables...'
      call rdi95(lun,optabe,ntab,nxir,nyir,nzva,rlg,tlg,opa,ival)
      print *,'read ALEX94 tables...'
      call readal(lun+2,lowtab,ntab,nzva,nxir,ntal,rlgal,tlgal,opaal)
c
c  --------  Partial derivatives calculation --------------
c
c     open file for partial derivatives
      open(lun,file=paderi,status='UNKNOWN',form='unformatted')
c
      print '(/,a)','evaluate partial derivatives for OPAL95+ALEX94'
      print '(a,/)','----------------------------------------------'
      write(nout,8610)
c     repeat for nzva tables (nzva=13 for OPAL95)
      do 8010 l=1,nzva
c
c       define xd- ,yd- and zd- values (each of size ndat)
c       and define the xx- and yy- arrays
c       ( nxi x nyi coordinates at the desired z-values)
        do 8001 itab=1,ntab
          nd = 0
          do 1002 j=nyis,nyir
           do 1001 i=1,nxir
               nd = nd + 1
               zd(nd) = opa(i,j,itab,l)
               xd(nd) = rlg(i,itab,l)
               yd(nd) = tlg(j,itab,l)
 1001      continue
 1002     continue
c
c         load Alex94-table
c
          do js=1,ntal
           do is=1,nxir
            nd=nd+1
            zd(nd) = opaal(is,js,itab,l)
            xd(nd) = rlgal(is,itab,l)
            yd(nd) = tlgal(js,itab,l)
           enddo
          enddo
c
          write(nout,8600)l,itab,nd
c
c         interpolate at xi=-1.,yi=5. (dummy values)
c         and evaluate partial derivatives
c         stored in the wk array
c
          call masube(drelpr,1,0,nd,xd,yd,zd,tp,1,1,-1.d0,5.d0,z,1,
     +                nout,iwk,wk)
c
c         output of pointers, indices of triangles and derivative-values
          if((itab.eq.1).and.(l.eq.1)) then
             write(lun)nzva
             write(lun)ntab
             write(lun)nt
             write(lun)ial
             write(lun)(iwk(m),m=8,3*nt+7)
          endif
          write(nout,'(2i5)')nt,ial
          write(lun)(wk(m),m=5,ial-1)
c
8001    continue
8010  continue
c
c
      close(lun)
c
      stop
c
c-------------------------------------------------------------
c
8600  format('evaluate partial derivative values for table:',
     +        2i3,i5,$)
8610  format(47x,'z  x  nd   nt   ial')
9030  print *,'dopdalex94: ERROR in opening of OPINTPATH'
      stop
c
      end
c
c ------------------END of Main program ---------------------
