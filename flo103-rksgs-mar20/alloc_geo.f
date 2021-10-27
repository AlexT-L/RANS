      subroutine alloc_geo
c
c     ******************************************************************
c     *                                                                *
c     *   allocate all arrays for the mesh generator                   *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use geo_var
c
c     ******************************************************************
c
      allocate (a0(il))
      allocate (a1(il))
      allocate (b0(jl))
      allocate (s0(il))

      return

      end
