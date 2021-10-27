      subroutine gmesh
c
c     ******************************************************************
c     *                                                                *
c     *   generates the mesh                                           *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use geo_var
      use mesh_var
c
c     ******************************************************************
c
      use flo_param
      use geo_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      common/tit/ title
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
      common/ims/ imesh,iin
c
c     ******************************************************************
c
      integer  :: imesh,iin
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: iwrit
      integer  :: i,j
      integer  :: ile,ite
c
c     ******************************************************************
c
      real     :: xmax,xmin
c
c     ******************************************************************
c
      iwrit     = 6
c
c     set the mesh dimensions
c
      il        = nx  +1
      jl        = ny  +1
      ie        = nx  +2
      je        = ny  +2
      ib        = nx  +3
      jb        = ny  +3
c
c     set the limits of the profile
c
      ite       = .5000005*xte*(il  -1)
      ile       = il/2  +1
      itl       = ile  -ite
      itu       = ile  +ite
c
c     if imesh.ge.1 read the mesh
c
      if (imesh.ge.1) go to 21
c
c     set the limits of the outer mesh for a viscous simulation
c
      if (kvis.gt.1) then
         nbl       = jl/2
         ny        = ny  -nbl
         jl        = jl  -nbl
      end if
c
c     define point distributions in each coordinate direction
c
      call coord
c
c     define the transformation for generating the mesh
c
      call geom
c
c     generate the mesh
c
      call mesh
c
c     print the mesh parameters
c
      write (iwrit,600)
      write (iwrit,12)
   12 format(1x,'mapped coordinates and far field function'/
     .       3x,'     a0     ','     s0     ')

      do i=1,il
         write (iwrit,610) a0(i),s0(i)
      end do

      write (iwrit,600)
      write (iwrit,14)
   14 format(3x,'     b0     ')

      do j=1,jl
         write (iwrit,610) b0(j)
      end do

      write (iwrit,18)
   18 format(3x,'   boundx   ','   boundy   ','    bunch   ',
     .          '     xte    ','    ylim1   ','    ylim2   ')
      write (iwrit,610) boundx,boundy,bunch,xte,ylim1,ylim2

      write (iwrit,20)
   20 format(3x,'     ax     ','     ay     ','     sy     ')
      write (iwrit,610) ax,ay,sy
c
c     option to plot the euler mesh
c
      if (imesh.lt.-1) call grid(title)
c
c     insert an inner sublayer in the mesh for viscous simulations
c
      if (kvis.gt.0) then
         jlinv     = jl
         nyinv     = ny
         ny        = ny  +nbl
         jl        = jl  +nbl
         call vmesh
      end if

      write (6,1000) kvis,il,jl,nbl
 1000 format(1x,'gmesh,kvis,il,jl,nbl',4i6)
c
c     sangho's modification for non-dimentionalization
c
      xmax       = x(itl,1,1)
      xmin       = x(itl,1,1)

      do i=itl,itu
         xmin       = min(xmin,x(i,1,1))
      end do

      scal       = xmax  -xmin

      do i=1,il
      do j=1,jl
         x(i,j,1)   = x(i,j,1)/scal
         x(i,j,2)   = x(i,j,2)/scal
      end do
      end do

      xmax       = xmax/scal
      xmin       = xmin/scal
      scal       = 1.

      chord      = scal*(xmax  -xmin)
      xm         = scal*xmin  +.25*chord
      ym         = scal*x(itl,1,2)
      write (6,1006) itl,itu
      write (6,1008) xte,scal,chord
c
c     option to plot the viscous mesh
c
      if (imesh.le.-1) call grid(title)

      return
c
c     read the mesh from mesh.d
c
   21 open  (iin,file='mesh.d',access='sequential',form='unformatted')

      read  (iin) il,jl,itl,itu
      read  (iin) ((x(i,j,1),i=1,il),j=1,jl),
     .            ((x(i,j,2),i=1,il),j=1,jl)

      itu        = il  -itl  +1
      xte        = float(itu  -itl)/float(il  -1)
      scal       = 1.
      xmax       = x(itl,1,1)
      xmin       = x(itl,1,1)

      do i=itl,itu
         xmin       = min(xmin,x(i,1,1))
      end do

      chord      = scal*(xmax  -xmin)
      xm         = scal*xmin  +.25*chord
      ym         = scal*x(itl,1,2)

      write (6,1002) iin
 1002 format(1x,'mesh read from unit ',i6)
      write (6,1004) il,jl
 1004 format(1x,'il,jl',2i6)
      write (6,1006) itl,itu
 1006 format(1x,'itl,itu',2i6)
      write (6,1008) xte,scal,chord
 1008 format(1x,'xte,scal,chord',3f15.8)

c     call wplot

      return
c
c     ******************************************************************
c
  600 format(1x)
  610 format(1x,6f12.4)

      end
