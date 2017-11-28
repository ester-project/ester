      subroutine trdisa(n,a,b,c,eps,iflag)
      implicit double precision(a-h,o-z)
      dimension a(n),b(n),c(n)
      iflag=0
      n1=n-1
      do 10 k=2,n
         km1=k-1
         h=b(km1)
         if(abs(h).lt.eps)then
           iflag=2
           return
         endif
         h=a(km1)/h
         c(km1)=h
         b(k)=b(k)-h*a(km1)
  10  continue
      return
      end
