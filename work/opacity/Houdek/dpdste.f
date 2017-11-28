      subroutine pdste(n,x,y,z,nt,ipt,pd,ipd)
      implicit double precision (a-h,o-z)
c it estimates the first order partial derivative values
c at the data points following the klucewicz method.
c
c the input parameters are
c     n   = number of data points,
c     x,y,z = arrays of dimension n containing the x,y and z
c             coordinates of the data points,
c     nt  = number of triangles,
c     ipt = integer array of dimension 3*nt containing the indices
c           of the vertexes of the triangles.
c
c the output parameter is
c     pd  = array of dimension 2*n containing the partial derivative
c           values at the data points.
c
c the other parameter is
c     ipd = integer array of dimension n used internally as work area.
c
c declaration statement.
      dimension x(*),y(*),z(*),ipt(*),pd(*),ipd(*)
c preliminary processing.
      n2=n+n
        do 10 i=1,n2
          pd(i)=0.d0
          ipd(i)=0
10      continue
c estimates for each triangle the slopes of the plane through the
c function's value at the vertexes.
      do 20 i=1,nt
        l=3*i
        i1=ipt(l-2)
        i2=ipt(l-1)
        i3=ipt(l)
        x21=x(i2)-x(i1)
        x31=x(i3)-x(i1)
        y21=y(i2)-y(i1)
        y31=y(i3)-y(i1)
        z21=z(i2)-z(i1)
        z31=z(i3)-z(i1)
        c=y21*x31-x21*y31
        dx=(y21*z31-z21*y31)/c
        dy=(z21*x31-x21*z31)/c
c  updates the ipd and pd arrays.
        ipd(i1)=ipd(i1)+1
        i1=i1+i1
        pd(i1-1)=pd(i1-1)+dx
        pd(i1)=pd(i1)+dy
        ipd(i2)=ipd(i2)+1
        i2=i2+i2
        pd(i2-1)=pd(i2-1)+dx
        pd(i2)=pd(i2)+dy
        ipd(i3)=ipd(i3)+1
        i3=i3+i3
        pd(i3-1)=pd(i3-1)+dx
        pd(i3)=pd(i3)+dy
20    continue
c  averages the derivative values stored in the pd array.
      j=0
        do 30 i=2,n2,2
          j=j+1
          den=ipd(j)
          pd(i-1)=pd(i-1)/den
          pd(i)=pd(i)/den
30      continue
      return
      end
