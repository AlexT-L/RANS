      subroutine alloc_mesh
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the mesh                                 *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use mesh_var
c
c     ******************************************************************
c
      allocate (x(il,jl,2))
      allocate (xc(0:ib,0:jb,2))

      allocate (vol(0:ib,0:jb))

      allocate (pori(il,jl))
      allocate (porj(il,jl))

      allocate (fint(il,jl))

      return

      end
