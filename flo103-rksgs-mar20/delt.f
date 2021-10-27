      subroutine delt
c
c     *****************************************************************
c     *                                                               *
c     *   calculates the boundary layer thickness                     *
c     *                                                               *
c     *****************************************************************
c
      use dims
c
c     ******************************************************************
c
      use flo_var
      use mesh_var
      use solv_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
c
c     ******************************************************************
c
      real, dimension(je)               :: qs,dn,ut
c
c     ******************************************************************
c
      logical locke
c
c     ******************************************************************
c
c     js        = 2*jl/3  -2
      js        = .75*(ny  -4)

      do 20 i=2,il

      qs(1)     = 0.
      ut(1)     = 0.
      j         = js
      xy        = .5*(x(i,j,1)  -x(i,j-1,1)
     .               +x(i-1,j,1)  -x(i-1,j-1,1))
      yy        = .5*(x(i,j,2)  -x(i,j-1,2)
     .               +x(i-1,j,2)  -x(i-1,j-1,2))
      qs(j)     = (yy*w(i,j,2)  -xy*w(i,j,3))/(w(i,j,1))
      si        = sign(1.,qs(js))

      do j=2,js
         xy        = .5*(x(i,j,1)  -x(i,j-1,1)
     .                  +x(i-1,j,1)  -x(i-1,j-1,1))
         yy        = .5*(x(i,j,2)  -x(i,j-1,2)
     .                  +x(i-1,j,2)  -x(i-1,j-1,2))
         dsi       = 1./sqrt(xy**2  +yy**2)
         qs(j)     = si*(yy*w(i,j,2)  -xy*w(i,j,3))
         dn(j)     = 1./dsi
         ut(j)     = qs(j)*dsi
      end do

      dsti(i)   = 0.
      ynot(i)   = 0.
      ssmax(i)  = 0.
      cdu       = .98
      jmax      = idmax(js,ut,1)
      fx        = .6*ut(jmax)
      lend      = 2
      lbig      = 2
      locke     = .false.

      do  j=3,js
         if ( ut(j-1) .lt. 0. .and. ut(j) .ge. .0) then
            lbig      = j
         end if
         if (.not. locke) then
            if ( ut(j-1) .ge. cdu*ut(j) .and. ut(j) .gt. fx) then
               locke     = .true.
               lend      = j
            end if
         end if
      end do

      uinf      = 1./ut(lend)
      do j=lbig,lend
         dsti (i)  = dsti(i) + (ut(lend)*dn(j) - qs(j))
         ra        = w(i,lend,1)/w(i,j,1)
         ssmax(i)  = ssmax(i) + ra*ut(j)*uinf*(dn(j)-qs(j)*uinf)
      end do

      dsti(i)   = dsti(i)*uinf
      dsti (i)  = max(dsti(i),1.e-6)
      ra        = w(i,lend,1)/w(i,2,1)
      ssmax(i)  = max(ssmax(i),ra*qs(2)*uinf)
      fc        = .95*ut(lend)

      jse       = lend
      do j=2,jse
         lend      = j  -1
         if (ut(j).gt.fc) go to 11
      end do

   11 xbi       = .5*(x(i,1,1)+x(i-1,1,1))
      ybi       = .5*(x(i,1,2)+x(i-1,1,2))
      ycorr     = sqrt((xc(i,lend,1) - xbi)**2+(xc(i,lend,2)-ybi)**2)
      ynot(i)   = 1.5*(ycorr  +dn(lend)*
     .                (fc  -ut(lend))/(ut(lend+1)  -ut(lend)))

   20 continue

      return

      end
