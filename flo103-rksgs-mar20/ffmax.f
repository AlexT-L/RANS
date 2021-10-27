      subroutine ffmax(ls,le,fs,ys,fx,yx,lx)
c
c     ******************************************************************
c
      dimension fs(1),ys(1)
c
c     ******************************************************************
c
      do l=ls,le
         if (fs(l).ge.fx) then
            fx = fs(l)
            yx = ys(l)
            lx = l
         end if
      end do

      if (lx.gt.ls.and.lx.lt.le) then
         y1 = ys(lx-1)
         y2 = ys(lx)
         y3 = ys(lx+1)
         a1 = fs(lx-1)/((y1-y2)*(y1-y3))
         a2 = fs(lx  )/((y2-y1)*(y2-y3))
         a3 = fs(lx+1)/((y3-y1)*(y3-y2))
         yx = (a1*(y2+y3) + a2*(y1+y3) +a3*(y1+y2))/(2.*(a1+a2+a3))
         fx = a1*(yx-y2)*(yx-y3)+a2*(yx-y1)*(yx-y3)+a3*(yx-y1)*(yx-y2)
      end if

      return

      end
