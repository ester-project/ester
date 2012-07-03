      subroutine inttwo(x,n,y,m,v,w,i,j,iflag)
      implicit double precision(a-h,o-z)
      dimension x(n),y(m)
      iflag=0
      if(v.lt.x(1).or.v.gt.x(n).or.w.lt.y(1).or.w.gt.y(m))then
        iflag=3
        return
      endif
      if(i.lt.1.or.i.ge.n.or.j.lt.1.or.j.ge.m)then
        i=1
        j=1
      endif
      if(v.lt.x(i))  goto 10
      if(v.le.x(i+1))goto 40
      l=n
      goto 30
  10  l=i
      i=1
  20  k=(i+l)/2
      if(v.lt.x(k))then
        l=k
      else
        i=k
      endif
  30  if(l.gt.i+1) goto 20
  40  if(w.lt.y(j))goto 50
      if(w.le.y(j+1))return
      l=m
      goto 70
  50  l=j
      j=1
  60  k=(j+l)/2
      if(w.lt.y(k))then
        l=k
      else
        j=k
      endif
  70  if(l.gt.j+1) goto 60
      return
      end
