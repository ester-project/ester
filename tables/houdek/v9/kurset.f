      subroutine kurset(ir,infile)
c======================================================================
c
      implicit real*8 (a-h,o-z)
      character*(*) infile
      character*80 headr1,headr2
c
      parameter(mt=56,mp=38,iv=3,nkur=11,lkur=2128)
c
c ... var(i,j,1)=log p, var(i,j,2)=log op, var(i,j,3)=log rho
c
c ... vark(1...nkur) (of original Kurucz data file):
c ... logT logP 0km/s 1km/s 2km/s 4km/s 8km/s log Ne log Na log rho IPloCM-1
c
      common/kurtab/tl(mt),var(mt,mp,iv),zv,nt,np(mt)
      dimension vark(nkur)
c
      open(unit=ir,status='old',file=infile,err=9030)
c
      read(ir,'(a80)') headr1
      read(ir,'(a80)') headr2
      write(6,*) headr1
      read(headr1(42:45),'(f4.1)') zv
c     write(6,*) headr2
c
      tlgold = 0.
      it = 0
      ip = 0
      do 100 ikur=1,3000
      read(ir,9000,end=1000) (vark(j),j=1,nkur)
      lines = ikur
      tlgnew = vark(1)
      if(tlgnew.ne.tlgold) then
         it = it + 1
         ip = 0
         tlgold = tlgnew
      end if
      nt = it
      ip = ip + 1
      np(it) = ip
      tl(it) = tlgnew
      var(it,ip,1) = vark(2)
      var(it,ip,2) = vark(3)
      var(it,ip,3) = vark(10)
 100  continue
c
1000  if(lines.ne.lkur) then
         write(6,*) 'kurucz opacity: wrong line number in table.'
         write(6,*) 'lines expected, read: ',lkur,lines
      else
         write(6,'(2a,2i3,i5)') 'table successfully read. nt,np(nt),',
     +                         'lines = ',nt,np(nt),lines
      end if
c
      return
 9000 format(2f5.2,5f7.3,3f9.5,f8.3)
 9030 print *,'error in opening ',infile
      stop
c
      end
