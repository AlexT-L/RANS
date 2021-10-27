      function cutime(arg)
c
c     ******************************************************************
c     *                                                                *
c     *   cumulative computing time                                    *
c     *                                                                *
c     ******************************************************************
c
      real*4 tt(2)
c     nsec      = mclock()
c     cutime    = .01*float(nsec)
      cutime    = etime(tt)
      return
      end
