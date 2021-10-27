      subroutine  zoomgrid(title,i1,i2,j1,j2)
c
c     ******************************************************************
c     *                                                                *
c     *   plots the mesh                                               *
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
      implicit none
c
c     ******************************************************************
c
      common/psp/ npgrid,npcp,npcnv
c
c     ******************************************************************
c
      integer  :: npgrid,npcp,npcnv
      integer  :: i1,i2,j1,j2
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j
      integer  :: kp
c
c     ******************************************************************
c
      real     :: xp,yp
      real     :: xmax,xmin,ymax,ymin,scale
c
c     ******************************************************************
c
      character(80) :: b
c
c     ******************************************************************
c
c     open a new page
c
c     open  (18,file='GRID.PLOT',access='append')
      open  (18,file='GRID.PLOT',position='append')

      npgrid = npgrid +1
      call initpl(2.,1.25,1.,npgrid)
c
c     write the subtitles
c
      write (b,12)  title
   12 format(a80)
      call symbol(0.,-.25,.20,b,0.,80)
      write(b,14) nx,ny
   14 format('Grid ',i4,' x ',i4)
      call symbol(0.,-.5,.18,b,0.,16)
      write(b,16) i1,i2,j1,j2
   16 format('I1 ',i4,'  I2 ',i4,'  J1 ',i4,'  J2 ',i4)
      call symbol(0.,-.75,.18,b,0.,34)
c
c     calculate the scale
c
      xmax      = x(i1,j1,1)
      xmin      = x(i1,j1,1)
      ymax      = x(i1,j1,2)
      ymin      = x(i1,j1,2)

      do j=j1,j2
      do i=i1,i2
         xmax      = max(x(i,j,1),xmax)
         xmin      = min(x(i,j,1),xmin)
         ymax      = max(x(i,j,2),ymax)
         ymin      = min(x(i,j,2),ymin)
      end do
      end do

      scale     = 5.0/max((xmax  -xmin),(ymax  -ymin))
c
c     draw the grid lines in the i direction
c
      do j=j1,j2
         kp        = 3
      do i=i1,i2
         xp        = scale*(x(i,j,1)  -xmin)
         yp        = scale*(x(i,j,2)  -ymin)  +1.5
         call plot(xp,yp,kp)
         kp        = 2
      end do
      end do
c
c     draw the grid lines in the j direction
c
      do i=i1,i2
         kp        = 3
      do j=j1,j2
         xp        = scale*(x(i,j,1)  -xmin)
         yp        = scale*(x(i,j,2)  -ymin)  +1.5
         call plot(xp,yp,kp)
         kp        = 2
      end do
      end do
c
c     close the plot
c
      call endplt
      endfile 18
      close (18)

      return

      end
