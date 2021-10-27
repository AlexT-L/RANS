      subroutine alloc_mg
c
c     ******************************************************************
c     *                                                                *
c     *   allocate arrays for the multigrid scheme                     *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use mg_var
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: ix,jx,idx,idx2,idx4
c
c     ******************************************************************
c
      ix        = nx
      jx        = ny
      idx       = (ix  +4)*(jx  +4)

   11 if (mod(ix,2).ne.0.or.mod(jx,2).ne.0) go to 21

      ix        = ix/2
      jx        = jx/2
      idx       = idx  +(ix  +4)*(jx  +4)
      go to 11

   21 idx2      = 2*idx
      idx4      = 4*idx

      allocate (ww(idx4))
      allocate (ww1(idx4))
      allocate (wwr(idx4))

      allocate (dtlc(idx))
      allocate (rflc(idx))
      allocate (rlvc(idx))
      allocate (revc(idx))
      allocate (radic(idx))
      allocate (radjc(idx))

      allocate (xx(idx2))
      allocate (xxc(idx2))
      allocate (volc(idx))

      return

      end
