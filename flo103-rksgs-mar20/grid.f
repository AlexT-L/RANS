      subroutine  grid(title)
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
      real     :: xmax,xmin,scale
      real     :: xfmax,xfmin,yfmax,yfmin
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
      call symbol(0.,-.5,.20,b,0.,80)
      write(b,14) nx,ny
   14 format('Grid ',i4,' x ',i4)
      call symbol(0.,-.75,.14,b,0.,16)
c
c     calculate the scale
c
      xmax      = x(itl,1,1)
      xmin      = x(itl,1,1)

      do i=itl+1,itu
         xmax      = max(x(i,1,1),xmax)
         xmin      = min(x(i,1,1),xmin)
      end do

      scale     = 2.5/(xmax  -xmin)
c
c     draw the grid lines in the i direction
c
      do j=1,jl
         kp        = 3
      do i=1,il
         xp        = scale*(x(i,j,1)  -xmin)  +1.
         yp        = scale*x(i,j,2)  +4.5
         call plot(xp,yp,kp)
         kp        = 2
      end do
      end do
c
c     draw the grid lines in the j direction
c
      do i=1,il
         kp        = 3
      do j=1,jl
         xp        = scale*(x(i,j,1)  -xmin)  +1.
         yp        = scale*x(i,j,2)  +4.5
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
c
c     find the mesh limits
c
      xfmax     = x(itl,1,1)
      xfmin     = x(itl,1,1)
      yfmax     = x(itl,1,2)
      yfmin     = x(itl,1,2)

      do j=1,jl
      do i=1,il
         xfmax     = max(x(i,j,1),xfmax)
         xfmin     = min(x(i,j,1),xfmin)
         yfmax     = max(x(i,j,2),yfmax)
         yfmin     = min(x(i,j,2),yfmin)
      end do
      end do

      write (6,1002) nx,ny,xmax,xmin,xfmax,xfmin,yfmax,yfmin
 1002 format(1x,'nx,ny',2i4,
     .       2x,'xmax,xmin,xfmax,xfmin,yfmax,yfmin',6f8.2)

      return

      end
