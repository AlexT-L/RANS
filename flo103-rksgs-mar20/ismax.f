      function ismax(n,sx,incx)
c
c     find the largest value of the
c     elements of sx..if n<=0, return 0
c
      dimension sx(n)
c
      ismax = 0
      if(n .le. 0) return
      sxmax = 0
      do 100 i=1,n,incx
      if (sx(i) .gt. sxmax) then
      sxmax = sx(i)
      ismax = i
      end if
  100 continue
c
      return
      end
