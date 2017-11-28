      subroutine pdmin(drelpr,n,x,y,z,nout,iadve,nadve,index,tp,pd,
     *                 alfa,beta,dgamma,eijq,alfas,betas,gammas,zx,zy)
      implicit double precision (a-h,o-z)
c it estimates the first order partial derivative values at
c the data points by means of a global method based on a minimum
c norm under tension network .
c
c the input parameters are
c     drelpr = double relative precision
c     n   = number of data points,
c     x,y,z = array of dimension n containing the x,y and z
c             coordinates of the data points,
c     nout = logical unit number for the standard output unit of the
c            system,
c     iadve = integer array of dimension 6*n-12 containing the indices
c             of the vertexes adjacent to each vertex in the
c             triangulation,
c     nadve = integer array of dimension n+1 containing the number of
c             the vertexes adjacent to each vertex in the triangulation
c     tp  = tension parameter,
c     pd  = array of dimension 2*n containing an initial evaluation
c           of the partial derivative values at the data points.
c
c the output parameter is
c     pd  = array of dimension 2*n containing the partial derivative
c           values at the data points.
c
c the other parameters are
c     index = integer array of dimension 6*n-15 used internally as
c             work area,
c     alfa,beta,gamma,eijq = arrays of dimension 3*n-6 used internally
c                            as work areas,
c     alfas,betas,gammas,zx,zy = arrays of dimension n used internally
c                                as work areas.
c
c  the  reler  constant in the data  initialization statement  is a
c  relative error tolerance to stop the iterative method.
c  therefore it is machine dependent;
c  the absolute error tolerance  tau  is then obtained by
c         tau=reler*amax1(abs(pd(i)),i=1,2*n)+2*n*drelpr.
c
c declaration statements.
      dimension x(*),y(*),z(*),iadve(*),nadve(*),pd(*)
      dimension index(*),alfa(*),beta(*),dgamma(*),eijq(*),
     *          alfas(*),betas(*),gammas(*),zx(*),zy(*)
      data reler/1.d-15/
c  calculates the part of matrix coefficients independent
c  from the tension parameter tp.
      k=0
      pdm=0.0d0
        do 60 i=1,n
          j=i+i
          pdm=dmax1(pdm,dabs(pd(j)),dabs(pd(j-1)))
          zx(i)=0.d0
          zy(i)=0.d0
          alfas(i)=0.d0
          betas(i)=0.d0
          gammas(i)=0.d0
          jin=nadve(i)+1
          jfin=nadve(i+1)
            do 50 j=jin,jfin
              ind=iadve(j)
              dx=x(i)-x(ind)
              dy=y(i)-y(ind)
              dz=z(i)-z(ind)
              dxq=dx*dx
              dyq=dy*dy
              alunq=dxq+dyq
              al3=alunq*dsqrt(alunq)
              zx(i)=zx(i)+dz*dx/al3
              zy(i)=zy(i)+dz*dy/al3
                if(ind.gt.i)go to 30
                  lin=nadve(ind)+1
                  lfin=nadve(ind+1)
                    do 10 l=lin,lfin
                      if(i.eq.iadve(l))go to 20
   10               continue
   20             index(j)=index(l)
                  go to 40
   30         k=k+1
              index(j)=k
              al3p2=al3+al3
              eijq(k)=alunq
              alfa(k)=dxq/al3p2
              beta(k)=dx*dy/al3p2
              dgamma(k)=dyq/al3p2
   40       alfas(i)=alfas(i)+alfa(index(j))
            betas(i)=betas(i)+beta(index(j))
            gammas(i)=gammas(i)+dgamma(index(j))
   50       continue
          zx(i)=3.d0*zx(i)/2.d0
          zy(i)=3.d0*zy(i)/2.d0
          alfas(i)=2.d0*alfas(i)
          betas(i)=2.d0*betas(i)
          gammas(i)=2.d0*gammas(i)
   60   continue
      if(tp.eq.0)go to 100
c  calculates the part of matrix coefficients depending from
c  the tension parameter tp.
      tpq=tp*tp
      tpq60=tpq/60.d0
      tpq40=tpq/40.d0
      tpq15=tpq/15.d0
        do 80 i=1,n
          al=0.d0
          be=0.d0
          ga=0.d0
          z1=0.d0
          z2=0.d0
          jin=nadve(i)+1
          jfin=nadve(i+1)
            do 70 j=jin,jfin
              k=index(j)
              ind=iadve(j)
              dz=z(i)-z(ind)
              al=al+alfa(k)*eijq(k)
              be=be+beta(k)*eijq(k)
              ga=ga+dgamma(k)*eijq(k)
              sqeijq=dsqrt(eijq(k))
              z1=z1+(x(i)-x(ind))*dz/sqeijq
              z2=z2+(y(i)-y(ind))*dz/sqeijq
   70       continue
          alfas(i)=alfas(i)+tpq15*al
          betas(i)=betas(i)+tpq15*be
          gammas(i)=gammas(i)+tpq15*ga
          zx(i)=zx(i)+tpq40*z1
          zy(i)=zy(i)+tpq40*z2
   80   continue
        m=nadve(n+1)/2
          do 90 i=1,m
            alfa(i)=alfa(i)-tpq60*alfa(i)*eijq(i)
            beta(i)=beta(i)-tpq60*beta(i)*eijq(i)
            dgamma(i)=dgamma(i)-tpq60*dgamma(i)*eijq(i)
   90     continue
c  calculates the solutions of the system following the gauss-seidel
c  method.
  100 iter=1
      tau=reler*pdm+2*n*drelpr
  110 erq=0.d0
        do 130 i=1,n
        adx=0.d0
        bdx=0.d0
        bdy=0.d0
        gdy=0.d0
        jin=nadve(i)+1
        jfin=nadve(i+1)
          do 120 j=jin,jfin
            ind=2*iadve(j)
            k=index(j)
            adx=adx+alfa(k)*pd(ind-1)
            bdy=bdy+beta(k)*pd(ind)
            bdx=bdx+beta(k)*pd(ind-1)
            gdy=gdy+dgamma(k)*pd(ind)
  120     continue
        det=alfas(i)*gammas(i)-betas(i)*betas(i)
        s1=adx+bdy-zx(i)
        s2=bdx+gdy-zy(i)
        pdxn=(-gammas(i)*s1+betas(i)*s2)/det
        pdyn=(betas(i)*s1-alfas(i)*s2)/det
        ipi=i+i
        erx=pdxn-pd(ipi-1)
        ery=pdyn-pd(ipi)
        erq=erq+erx*erx+ery*ery
        pd(ipi-1)=pdxn
        pd(ipi)=pdyn
  130 continue
c  checks whether convergence is reached with the prescribed tolerance.
      if(erq.lt.tau)return
        if(iter.eq.20)go to 150
          iter=iter+1
          go to 110
c
c  error exit.
  150 write(nout,1)
      return
    1 format(1x,26hminimization not completed)
      end
