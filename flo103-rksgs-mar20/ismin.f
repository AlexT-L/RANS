      function ismin(n,sx,incx)
c
c     find the largest value of the
c     elements of sx..if n<=0, return 0
c
      dimension sx(n)
c
      ismin = 0
      if(n .le. 0) return
      sxmin = 1000000000
      do 100 i=1,n,incx
      if (sx(i) .lt. sxmin) then
      sxmin = sx(i)
      ismin = i
      end if
  100 continue
c
      return
      end
