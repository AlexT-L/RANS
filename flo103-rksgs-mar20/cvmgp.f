      function cvmgp(x1,x2,x3)
c
c     tests for positive or zero.  x1 is returned if x3 >= 0.
c     x2 is returned if x3 < 0.
c
      if(x3.lt.0.0) then
        cvmgp=x2
      else
        cvmgp=x1
      endif
c
      return
      end
