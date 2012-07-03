      subroutine trdisb(n,a,b,c,d)
      implicit double precision(a-h,o-z)
      dimension a(n),b(n),c(n),d(n)
      d(1)=d(1)/b(1)
      do 10 k=2,n
         km1=k-1
         d(k)=(d(k)-a(km1)*d(km1))/b(k)
  10  continue
      do 20 k=n-1,1,-1
         d(k)=d(k)-c(k)*d(k+1)
  20  continue
      return
      end
