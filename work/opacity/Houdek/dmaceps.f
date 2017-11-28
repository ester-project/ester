      subroutine maceps(drelpr)
      implicit double precision (a-h,o-z)
c it carries out an approximation of the double relative precision
c constant
      drelpr=1.0d0
10    drelpr=0.5d0*drelpr
c the do loop statement is necessary when the arithmetic unit has more
c bits than in storage (intel 8087 family of arithmetic units is of
c this type), infact a do loop involves a store from register to memory
c of the value for drelpr.
      do 20 i=1,2
20    continue
      if (drelpr+1.0d0.gt.1.0d0) goto 10
      drelpr=2.0d0*drelpr
      return
      end
