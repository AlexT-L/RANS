      integer function idmin(n,sx,incx)
c
c     ******************************************************************
c     *                                                                *
c     *   find the index of element having min value                   *
c     *                                                                *
c     ******************************************************************
c
      integer   i,incx,ix,n
c
c     *****************************************************************
c
      real      sx(n),smin
c
c     *****************************************************************
c
      idmin     = 0
      if (n.lt.1) return

      idmin     = 1
      if (n.eq.1) return

      if (incx.eq.1) go to 11
c
c     code for increment not equal to 1
c
      ix        = 1
      smin      = sx(1)
      ix        = ix + incx
      do i=2,n
         if (sx(ix).lt.smin) then
            idmin     = i
            smin      = sx(ix)
         end if
         ix        = ix + incx
      end do

      return
c
c      code for increment equal to 1
c
   11 smin      = sx(1)
      do i=2,n
         if (sx(i).lt.smin) then
            idmin     = i
            smin      = sx(i)
         end if
      end do

      return

      end
