      subroutine solve_rs(f,b)
c
c      solve 4x4 system b*x = f
c      replaces x in place of f
c
      implicit none
      integer              :: i,j
      real, dimension(4,4)      :: b
      real, dimension(4)        :: f
c
      real     :: d1,d2,d3,d4

      b(1,1)     = 1./b(1,1)
      b(1,2)     = b(1,2)*b(1,1)
      b(1,3)     = b(1,3)*b(1,1)
      b(1,4)     = b(1,4)*b(1,1)
c
      b(2,2)     = 1./(b(2,2) -b(2,1)*b(1,2))
      b(2,3)     = (b(2,3)-b(2,1)*b(1,3))*b(2,2)
      b(2,4)     = (b(2,4)-b(2,1)*b(1,4))*b(2,2)
      b(3,2)     = b(3,2)-b(3,1)*b(1,2)
      b(3,3)     = 1./(b(3,3) -b(3,1)*b(1,3) - b(3,2)*b(2,3))
      b(3,4)     = (b(3,4) -b(3,1)*b(1,4) -b(3,2)*b(2,4))*b(3,3)
      b(4,2)     = b(4,2) -b(4,1)*b(1,2)
      b(4,3)     = b(4,3) -b(4,1)*b(1,3) -b(4,2)*b(2,3)
      b(4,4)     = 1./(b(4,4)-b(4,1)*b(1,4)-b(4,2)*b(2,4)-b(4,3)*b(3,4))
c
c            fwd substitution
c
      d1         = f(1)*b(1,1)
      d2         = (f(2) -b(2,1)*d1)*b(2,2)
      d3         = (f(3) -b(3,1)*d1-b(3,2)*d2)*b(3,3)
      d4         = (f(4) -b(4,1)*d1-b(4,2)*d2-b(4,3)*d3)*b(4,4)
c
c            bwd substitution
c
      f(4)       = d4
      f(3)       = d3 -b(3,4) *f(4)
      f(2)       = d2 -b(2,3)*f(3)-b(2,4)*f(4)
      f(1)       = d1-b(1,2)*f(2)-b(1,3)*f(3)- b(1,4)*f(4)
c
      return
c
      end
