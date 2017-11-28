      subroutine intrp(drelpr,tol,x,y,z,ipt,pd,iti,itpv,xii,yii,zii)
      implicit double precision (a-h,o-z)
c it carries out punctual interpolation, i.e., it determines
c the z value at a given point in a triangle by means of the
c 9-parameter discretized version of nielson's scheme.
c
c the input parameters are
c    drelpr,tol = double relative precision and tolerance,
c    x,y,z= arrays of dimension n containing the x,y and z coordinates
c         of the data points, where n is the number of the data points,
c    ipt = integer array of dimension 3*nt, where nt is the number of
c          triangles, containing the indices of the vertexes of the
c          triangles-themselves,
c    pd  = array of dimension 2*n containing the partial derivative
c          values at the data points,
c    iti = index of the triangle where the point for which
c          interpolation has to be performed, lies,
c    itpv= index of the triangle considered in the previous call,
c    xii,yii = x and y coordinates of the point for which inter-
c              polation has to be performed.
c
c the output parameter is
c    zii = interpolated z value.
c
c declaration statements.
      dimension x(*),y(*),z(*),ipt(*),pd(*)
      save
      ito=iti
        if(ito.eq.itpv)go to 10
c  selects the triangle containing the point (xii,yii).
          itpv=ito
          ito3=3*ito
          ind1=ipt(ito3-2)
          ind2=ipt(ito3-1)
          ind3=ipt(ito3)
c  calculates the basic quantities relatives to the selected triangle.
          x12=x(ind1)-x(ind2)
          x13=x(ind1)-x(ind3)
          x23=x(ind2)-x(ind3)
          y12=y(ind1)-y(ind2)
          y13=y(ind1)-y(ind3)
          y23=y(ind2)-y(ind3)
          d=x13*y23-x23*y13
          index=2*ind1
          d3v1=-x13*pd(index-1)-y13*pd(index)
          d2v1=-x12*pd(index-1)-y12*pd(index)
          index=2*ind2
          d1v2=x12*pd(index-1)+y12*pd(index)
          d3v2=-x23*pd(index-1)-y23*pd(index)
          index=2*ind3
          d2v3=x23*pd(index-1)+y23*pd(index)
          d1v3=x13*pd(index-1)+y13*pd(index)
          e1=x23*x23+y23*y23
          e2=x13*x13+y13*y13
          e3=x12*x12+y12*y12
          e=2.d0*e1
          alfa21=(e2+e1-e3)/e
          alfa31=(e3+e1-e2)/e
          e=2.d0*e2
          alfa12=(e1+e2-e3)/e
          alfa32=(e3+e2-e1)/e
          e=2.d0*e3
          alfa13=(e1+e3-e2)/e
          alfa23=(e2+e3-e1)/e
c  calculates the remaining quantities necessary for the zi evaluation
c  depending from the interpolation point.
   10   b1=((xii-x(ind3))*y23-(yii-y(ind3))*x23)/d
        b2=((yii-y(ind1))*x13-(xii-x(ind1))*y13)/d
        b3=1.d0-b1-b2
          if(b1.ge.tol)go to 30
            if(b2.ge.tol)go to 40
              if(b1+b2.le.drelpr)go to 50
        w=b1*b2*b3/(b1*b2+b1*b3+b2*b3)
        qb=b1*b1
        wb=w*b1
        zii=z(ind1)*(qb*(3.d0-2.d0*b1)+6.d0*wb*(b3*alfa12+b2*alfa13))+
     *      d3v1*(qb*b3+wb*(3.d0*b3*alfa12+b2-b3))+
     *      d2v1*(qb*b2+wb*(3.d0*b2*alfa13+b3-b2))
        qb=b2*b2
        wb=w*b2
        zii=zii+z(ind2)*(qb*(3.d0-2.d0*b2)+6.d0*wb*(b1*alfa23+b3*
     *      alfa21))+
     *      d1v2*(qb*b1+wb*(3.d0*b1*alfa23+b3-b1))+
     *      d3v2*(qb*b3+wb*(3.d0*b3*alfa21+b1-b3))
        qb=b3*b3
        wb=w*b3
        zii=zii+z(ind3)*(qb*(3.d0-2.d0*b3)+6.d0*wb*(b2*alfa31+b1*
     *      alfa32))+
     *      d2v3*(qb*b2+wb*(3.d0*b2*alfa31+b1-b2))+
     *      d1v3*(qb*b1+wb*(3.d0*b1*alfa32+b2-b1))
      return
   30 zii=z(ind1)
      return
   40 zii=z(ind2)
      return
   50 zii=z(ind3)
      return
      end
