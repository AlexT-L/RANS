      subroutine alloc_psm
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the psmw solution                        *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use psm_var
c
c     ******************************************************************
c
      allocate (rs(ie,je,4))

      allocate (r00(ie,je))
      allocate (p00(ie,je))
      allocate (u00(ie,je))
      allocate (v00(ie,je))

      return

      end
