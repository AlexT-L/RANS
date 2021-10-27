      subroutine alloc_solv
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the flow and adjoint solutions           *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use solv_var
c
c     ******************************************************************
c
      allocate (wn(ie,je,4))
      allocate (w1(ie,je,4))
      allocate (wr(ie,je,4))

      allocate (dw(0:ib,0:jb,4))
      allocate (fw(0:ib,0:jb,4))
      allocate (vw(0:ib,0:jb,3))

      allocate (fs(il,jl,4))

      allocate (dtl(ie,je))
      allocate (rfl(ie,je))

      allocate (rfli(ie,je))
      allocate (rflj(ie,je))

      allocate (radi(ie,je))
      allocate (radj(ie,je))

      return

      end
