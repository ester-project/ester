      subroutine adjac(nt,ipt,n,iadve,nadve)
      implicit double precision (a-h,o-z)
c it estimates the adjacencies monodimensional array containing
c for each vertex the indices of the vertexes adjacent in the
c triangulation.
c
c the input parameters are
c     nt  = number of triangles,
c     ipt = integer array of dimension 3*nt containing the indices
c           of the vertexes of the triangles,
c     n   = number of data points.
c
c the output parameters are
c     iadve = integer array of dimension  6*n-12  containing for each
c             vertex the indices of the vertexes adjacent in the
c             triangulation,
c     nadve = integer array of dimension  n+1  containing for each
c             vertex the number of the vertexes adjacent in the
c             triangulation.
c
c  declaration statement.
      dimension ipt(*),iadve(*),nadve(*),item(1000)
      nadve(1)=0
      kin=0
      nt3=3*nt
        do 70 i=1,n
          i2=0
c  stores the indices of the adjacent vertexes.
          do 30 j1=1,nt3,3
            j2=j1+2
            do 10 j=j1,j2
              if(i.eq.ipt(j))go to 20
   10       continue
          go to 30
   20       i1=i2+1
            i2=i2+2
            item(i1)=ipt(j1+mod(j,3))
            item(i2)=ipt(j1+mod(j+1,3))
   30     continue
c  discards the indices that have been stored twice.
          jin=kin+1
          kin=kin+2
          jfin=kin
          iadve(jin)=item(1)
          iadve(jfin)=item(2)
            if(i2.eq.2)go to 60
              do 50 j=3,i2
                do 40 l=jin,jfin
                  if(item(j).eq.iadve(l))go to 50
   40           continue
                kin=kin+1
                iadve(kin)=item(j)
                jfin=kin
   50         continue
   60       nadve(i+1)=kin
   70   continue
      return
      end
