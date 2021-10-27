      subroutine alloc_out
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the output                               *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use out_var
c
c     ******************************************************************
c
      allocate (xp(il))
      allocate (yp(il))
      allocate (cp(il))
      allocate (cf(il))

      allocate (xq(ie,je))
      allocate (yq(ie,je))
      allocate (qc(ie,je))

      allocate (icf(ie,je))

      return

      end
