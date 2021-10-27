      subroutine  vmesh
c
c     ******************************************************************
c     *                                                                *
c     *   mofifies the mesh for viscous simulations                    *
c     *   by inserting additional mesh intervals in a sublayer         *
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
      use flo_param
      use geo_param
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
      integer  :: i,j
      integer  :: nt,ind
      integer  :: ile
c
c     ******************************************************************
c
      real     :: xb0,yb0
      real     :: dx,dy,ds
      real     :: scut,dt,dta,dtb
      real     :: achord,rex,cfx,a,dbl1
c
c     ******************************************************************
c
      real, dimension(jl)          :: s,t
      real, dimension(jl)          :: xn,x1,x2,x3,xv
      real, dimension(jl)          :: yn,y1,y2,y3,yv
c
c     ******************************************************************
c
      jbl       = nbl  +ncut  +1
      ile       = il/2  +1

      write (6,1101) jlinv,nbl,ncut,jbl,jl,ile
 1101 format(1x,'jlinv,nbl,ncut,jbl,jl,ile',6i6)

c
c     estimate the mesh interval for yplus = aplus at 10 percent chord
c
      achord    = x(itu,1,1)  -x(ile,1,1)
      rex       = .1*re*achord
      a         = log(.06*rex)
      cfx       = .455/(a*a)
      dbl1      = aplus*achord/(re*sqrt(.5*cfx))

      write (6,1102) dbl1
 1102 format(1x,'dbl1',e12.4)

      do 10 i=1,il
c
c     save the surface mesh points
c
      xb0       = x(i,1,1)
      yb0       = x(i,1,2)
c
c     copy the coordinates along the mesh line
c
      do j=1,jlinv
         xn(j)     = x(i,j,1)
         yn(j)     = x(i,j,2)
      end do
c
c     calculate the arc length along an outgoing mesh line
c
      s(1)      = 0.

      do j=1,nyinv
         dx        = xn(j+1)  -xn(j)
         dy        = yn(j+1)  -yn(j)
         ds        = sqrt(dx*dx  +dy*dy)
         s(j+1)    = s(j)  +ds
      end do

      scut      = s(ncut+2)
c
c     fit splines to the x and y coordinates
c
      call splif (1,jlinv,s,xn,x1,x2,x3,2,0.,2,0.,0,0.,ind)
      call splif (1,jlinv,s,yn,y1,y2,y3,2,0.,2,0.,0,0.,ind)
c
c     shift the mesh indices to the outer inviscid range
c
      do j=jlinv,1,-1
         x(i,j+nbl,1) = x(i,j,1)
         x(i,j+nbl,2) = x(i,j,2)
         t(j+nbl)  = s(j)
      end do
c
c     redistribute the inner boundary layer mesh spacing
c
      nt        = jbl  +1
      dta       = dbl1/scut
      dtb       = (t(jbl+1)  -t(jbl))/scut

      call tanhds(nt,dta,dtb,t,ind)
c
c     scale the new distribution to match the original distribution
c
      do j=1,nt
         t(j)      = scut*t(j)
      end do
c
c     interpolate the mesh coordinates to the new spacing
c     for both the inner and the outer mesh
c
      call intpl (2,ny,t,xv,1,jlinv,s,xn,x1,x2,x3,0)
      call intpl (2,ny,t,yv,1,jlinv,s,yn,y1,y2,y3,0)

      dta       = scut*dta
      dtb       = scut*dtb
      write (6,1103) i,ind,scut,dta,dtb
 1103 format(1x,'i,iter',2i6,2x,'scut,dta,dtb',3f15.8)

      if (i.eq.ile) then

         do j=1,jlinv
            write (6,1104) i,j,s(j)
         end do

         do j=1,ny
            dt       = t(j+1)  -t(j)
            write (6,1105) i,j,t(j),dt
         end do

      end if

 1104 format(1x,'i,j',2i6,2x,'s',f15.8)
 1105 format(1x,'i,j',2i6,2x,'t,dt',2f15.8)

c
c     copy the surface coordinates
c
      x(i,1,1) = xb0
      x(i,1,2) = yb0
c
c     copy the coordinates along the mesh line
c
      do j=2,ny
         x(i,j,1) = xv(j)
         x(i,j,2) = yv(j)
      end do

   10 continue

      return

      end
