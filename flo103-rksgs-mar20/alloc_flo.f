      subroutine alloc_flo
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the flow solution                        *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use flo_var
c
c     ******************************************************************
c
      allocate (w(0:ib,0:jb,4))
      allocate (p(0:ib,0:jb))

      allocate (rlv(0:ib,0:jb))
      allocate (rev(0:ib,0:jb))
      allocate (aev(0:ib,0:jb))
      allocate (bev(0:ib,0:jb))

      allocate (tw(ie))
      allocate (cfric(ie))

      allocate (dsti(ie))
      allocate (ynot(ie))
      allocate (ssmax(ie))

      return

      end
