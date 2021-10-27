      subroutine setup
c
c     ******************************************************************
c     *                                                                *
c     *   sets up the metric and initial data                          *
c     *                                                                *
c     ******************************************************************
c
      use mg_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
c     flag the wall and far field boundary points
c
      call bound
c
c     calculate the cell areas
c
      call metric
c
c     define the initial data
c
      if (mode.lt.0) call init
c
c     update the surface pressure
c
      call bcwall
c
c     set extended values of the flow variables
c     in the halo surrounding the mesh
c
      call halo

      return

      end
