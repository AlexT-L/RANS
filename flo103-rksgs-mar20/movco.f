      subroutine movco
c
c     ******************************************************************
c     *                                                                *
c     *   move to a coarser mesh                                       *
c     *                                                                *
c     ******************************************************************
c
      use dims
      use dimsc
c
c     ******************************************************************
c
      use flo_var
      use mesh_var
      use solv_var
      use mg_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
      use geo_param
      use mg_param
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
      integer  :: i,j,k
      integer  :: nc,ile,ite
c
c     ******************************************************************
c
      real     :: rqq
c
c     ******************************************************************
c
c     set the coarse mesh dimensions
c
      iil       = nx/2  +1
      jjl       = ny/2  +1
      iie       = nx/2  +2
      jje       = ny/2  +2
      iib       = nx/2  +3
      jjb       = ny/2  +3
c
c     transfer the solution and residuals to a coarser mesh
c
      if (mode.ge.0)
     .call collc (ww(kw),wwr(kw),rlvc(k1),revc(k1))
c
c     expand the mesh
c
      if (mode.lt.0.or.ncyc.eq.1) call xpand (xx(kx))
c
c     exchange variables between the coarse and fine meshes
c
      call shift (ww(kw),ww1(kw),wwr(kw),xx(kx),xxc(kx),volc(k1),
     .            rlvc(k1),revc(k1))
c
c     reset the pointers to those of the coarse mesh
c
      nc        = (iib  +1)*(jjb  +1)
      k1        = k1  +nc
      kx        = kx  +2*nc
      kw        = kw  +4*nc
c
c     set the mesh dimensions to those of the coarse mesh
c
      nx        = nx/2
      ny        = ny/2
      il        = nx  +1
      jl        = ny  +1
      ie        = nx  +2
      je        = ny  +2
      ib        = nx  +3
      jb        = ny  +3
c
c     set the limits of the profile
c
      ite       = .50000005*xte*(il -1)
      ile       = il/2  +1
      itl       = ile  -ite
      itu       = ile  +ite

      if (mode.lt.0) return
c
c     calculate the cell areas
c
      if (ncyc.eq.1) call metric
c
c     flag the wall and far field boundary points
c
      call bound
c
c     calculate the pressure
c
      do j=2,je
      do i=1,ie
         rqq       = .5*(w(i,j,2)**2  +w(i,j,3)**2)/w(i,j,1)
         p(i,j)    = (gamma  -1.)*dim(w(i,j,4),rqq)
      end do
      end do
c
c     update the surface pressure
c
      call bcwall
c
c     option to update the far field
c
      if (fbc.gt.0.) call bcfar
c
c     set extended values of the flow variables
c     in the halo surrounding the mesh
c
      call halo

      return

      end
