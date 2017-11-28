      subroutine ratmat(h,pp,b)
      implicit double precision (a-h,o-z)
      dimension b(4,4)
      fa=h/(pp+3.d0)
      b(1,2)= b(1,1)*fa
      b(2,2)= b(2,1)*fa
      b(3,2)=-b(1,2)
      b(4,2)=-b(2,2)
      b(1,4)= b(4,2)
      b(2,4)= b(2,2)*(pp+2.d0)
      b(3,4)= b(2,2)
      b(4,4)=-b(2,4)
      return
      end
