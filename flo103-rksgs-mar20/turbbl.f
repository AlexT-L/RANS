      subroutine turbbl
c
c     *****************************************************************
c     *                                                               *
c     *   calculates the eddy viscosity                               *
c     *   using baldwin and lomax model                               *
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
      integer  :: i,j,i1,ii,js,jmax,jmin
      integer  :: idmin
c
c     ******************************************************************
c
      real     :: gm1,scf,sscf
      real     :: aplus,ccp,cwk,clk,cck,cmutm
      real     :: rho1,ru1,rv1,p1,u1,v1,t1
      real     :: rho2,ru2,rv2,p2,u2,v2,t2
      real     :: dxi,dxj,dyi,dyj,dui,duj,dvi,dvj,dti,dtj,dsj
      real     :: dux,duy,dvx,dvy,dtx,dty
      real     :: rlva,reva,rmu,rlam,rk,term
      real     :: sigx,sigy,tauxy
      real     :: dx13,dy13,dx24,dy24
      real     :: du13,dv13,du24,dv24
      real     :: dsij,dvdx,dudy
      real     :: avora,yplus,ex
      real     :: rhow,rmuw,wl,xbi,ybi
      real     :: udiff,fin,emumax,pex
c
c     ******************************************************************
c
      real, dimension(ie,je)            :: u,v,t,avor
      real, dimension(ie)               :: fmax,ymax,fwake
      real, dimension(je)               :: f,dnor,vel,emuin,emuout
      real, dimension(je)               :: fkleb
c
c     *****************************************************************
c
      logical     eul,lam,tur,rng,ftran,outer
c
c     *****************************************************************
c
      data aplus,ccp,cwk,clk,cck,cmutm/26.,1.6,1.,.4,.0168,14./
c
c     *****************************************************************
c
      ftran     = .false.
      gm1       = gamma - 1.
c     scf       = re/(sqrt(gamma)*rm)
      scf       = (scal*re/chord)/(sqrt(gamma)*rm)
      sscf      = sqrt(scf)
      js        = 2*jl/3.

      do j=1,je
      do i=1,ie
         u(i,j)    = w(i,j,2)/w(i,j,1)
         v(i,j)    = w(i,j,3)/w(i,j,1)
         t(i,j)    = p(i,j)/(gm1*w(i,j,1))
      end do
      end do
c
c     calculate the wall stress tauw
c
      do 10 i=2,il

      rho1    = .25*(w(i,1,1)  +w(i+1,1,1)  +w(i,2,1)  +w(i+1,2,1))
      ru1     = .25*(w(i,1,2)  +w(i+1,1,2)  +w(i,2,2)  +w(i+1,2,2))
      rv1     = .25*(w(i,1,3)  +w(i+1,1,3)  +w(i,2,3)  +w(i+1,2,3))
      p1      = .25*(p(i,1)    +p(i+1,1)    +p(i,2)    +p(i+1,2))
      u1      = ru1/rho1
      v1      = rv1/rho1
      t1      = p1/(gm1*rho1)

      rho2    = .25*(w(i-1,1,1)  +w(i,1,1)  +w(i-1,2,1)  +w(i,2,1))
      ru2     = .25*(w(i-1,1,2)  +w(i,1,2)  +w(i-1,2,2)  +w(i,2,2))
      rv2     = .25*(w(i-1,1,3)  +w(i,1,3)  +w(i-1,2,3)  +w(i,2,3))
      p2      = .25*(p(i-1,1)    +p(i,1)    +p(i-1,2)    +p(i,2))
      u2      = ru2/rho2
      v2      = rv2/rho2
      t2      = p2/(gm1*rho2)

      dxi       = x(i,1,1)  -x(i-1,1,1)
      dyi       = x(i,1,2)  -x(i-1,1,2)
      dxj       = xc(i,2,1)  -xc(i,1,1)
      dyj       = xc(i,2,2)  -xc(i,1,2)
      dsj       = 1./(dxi*dyj  -dyi*dxj)

      dui       = u1  -u2
      dvi       = v1  -v2
      dti       = t1  -t2
      duj       = u(i,2)  -u(i,1)
      dvj       = v(i,2)  -v(i,1)
      dtj       = t(i,2)  -t(i,1)

      dux       = (dui*dyj  -duj*dyi)*dsj
      dvx       = (dvi*dyj  -dvj*dyi)*dsj
      dtx       = (dti*dyj  -dtj*dyi)*dsj
      duy       = (duj*dxi  -dui*dxj)*dsj
      dvy       = (dvj*dxi  -dvi*dxj)*dsj
      dty       = (dtj*dxi  -dti*dxj)*dsj

      rlva      = .5*(rlv(i,2)  +rlv(i,1))
      reva      = .5*(rev(i,2)  +rev(i,1))
      rmu       = rlva  +reva
      rk        = gamma*(rlva/prn  +reva/prt)
      rlam      = -2.*rmu/3.

      term      = rlam*(dux  +dvy)
      sigx      = term  +2.*rmu*dux
      sigy      = term  +2.*rmu*dvy
      tauxy     = rmu*(duy  +dvx)
      tw(i)     = (dxi*dyi*(sigx  -sigy)  -(dxi**2  -dyi**2)*tauxy)/
     .            (dxi**2  +dyi**2)

   10 continue
c
c     calculate the vorticity for each cell
c
      do j=1,jl
      do i=1,il
         term      = rlam*(dux  +dvy)
         dx13      = xc(i,j,1)   - xc(i+1,j+1,1)
         dy13      = xc(i,j,2)   - xc(i+1,j+1,2)
         dx24      = xc(i+1,j,1) - xc(i,j+1,1)
         dy24      = xc(i+1,j,2) - xc(i,j+1,2)
         du13      = u(i,j) - u(i+1,j+1)
         dv13      = v(i,j) - v(i+1,j+1)
         du24      = u(i+1,j) - u(i,j+1)
         dv24      = v(i+1,j) - v(i,j+1)
         dsij      = 1./(dx13*dy24 - dx24*dy13)
         dvdx      =  dsij * (dv13*dy24 - dv24*dy13)
         dudy      = -dsij * (du13*dx24 - du24*dx13)
         avor(i,j) = abs(dudy-dvdx)
         rev(i,j)  = 0.
      end do
      end do
c
c     calculate the turbulent viscosity at interior points
c
      do 20 i=2,il

      fmax(i)   = 0.
      ymax(i)   = 0.
      fwake(i)  = 0.

      do j=1,js
         f(j)      = 0.
         vel(j)    = 0.
         emuin(j)  = 0.
         emuout(j) = 0.
         dnor(j)   = 0.
         fkleb(j)  = 0.
      end do

      rhow      = .5*(w(i,1,1)  +w(i,2,1))
      rmuw      = .5*(rlv(i,1)  +rlv(i,2))
      wl        = sqrt(abs(tw(i))*rhow)/rmuw

      xbi       = .5*(x(i-1,1,1)  +x(i,1,1))
      ybi       = .5*(x(i-1,1,2)  +x(i,1,2))

      if (i.gt.itl.and.i.le.itu) then

      vel(1)    = 0.

      do j=2,js
         avora     = .25*(avor(i-1,j-1)  +avor(i-1,j)
     .                   +avor(i,j-1)    +avor(i,j))
         dnor(j)   = sqrt((xc(i,j,1)  -xbi)**2  +(xc(i,j,2)  -ybi)**2)
         yplus     = wl*dnor(j)*sscf
         ex        = -yplus/aplus
         term      = exp(ex)
         emuin(j)  = scf*w(i,j,1)*avora*(clk*dnor(j)*(1.  -term))**2
         f(j)      = dnor(j)*avora*(1.  -term)
         vel(j)    = sqrt(u(i,j)**2  +v(i,j)**2)
      end do

      else

      vel(1)    = sqrt((.5*(u(i,1)+u(i,2)))**2
     .                +(.5*(v(i,1)+v(i,2)))**2)

      do j=2,js
         avora     = .25*(avor(i-1,j-1)  +avor(i-1,j)
     .                   +avor(i,j-1)    +avor(i,j))
         dnor(j)   = sqrt((xc(i,j,1)  -xbi)**2  +(xc(i,j,2)  -ybi)**2)
         emuin(j)  = scf*w(i,j,1)*avora*(clk*dnor(j))**2
         f(j)      = dnor(j)*avora
         vel(j)    = sqrt(u(i,j)**2  +v(i,j)**2)
      end do

      jmin      = idmin(js,vel,1)
      vel(1)    = vel(jmin)

      end if

      f(1)      = 0.

      call ffmax(2,js,f,dnor,fmax(i),ymax(i),jmax)

      udiff     = vel(jmax) - vel(1)
      if (fmax(i) .gt. 0.) then
         fwake(i)  = cwk*ymax(i)*udiff**2/fmax(i)
         fwake(i)  = min(fwake(i),ymax(i)*fmax(i))
      end if

      ymax(i)   = max(ymax(i),dnor(2))

      do j=2,js
         fin       = .3*dnor(j)/ymax(i)
         fkleb(j)  = 1./(1. + 5.5*fin**6)
         emuout(j) = scf*cck*ccp*w(i,j,1)*fwake(i)*fkleb(j)
      end do

      outer     = .false.
      emumax    = 0.
      do j=2,js
         if (emuin(j) .ge. emuout(j) .or. outer ) then
            rev(i,j)  = emuout(j)
            emumax    = max(emumax,rev(i,j))
            outer     = .true.
         else
            rev(i,j)  = emuin(j)
            emumax    = max(emumax,rev(i,j))
         end if
      end do
c
c     simulate transition  or fix transition
c
      if (.not. ftran .and. emumax .le. cmutm ) then
         do j=2,jl
            rev(i,j)  = 0.
         end do
      end if

      if (ftran .and. xc(i,2,1) .le. xtran) then
         do j=2,jl
            rev(i,j)  = 0.
         end do
      end if

   20 continue
c
c     adjust the near wake
c
      ii        = ie
      do i=2,itl+1
      ii        = ii  -1
      do j=2,je
         pex       = -(xc(i,2,1)-xc(itl+1,2,1))/(20.*ymax(itl+1))
         rev(i,j)  = rev(i,j)+(rev(itl+1,j)-rev(i,j))*exp(pex)
         pex       = -(xc(ii,2,1)-xc(itu,2,1))/(20.*ymax(itu))
         rev(ii,j) = rev(ii,j)+(rev(itu,j)-rev(ii,j))*exp(pex)
      end do
      end do
c
c     set boundary values along the cut in the c mesh
c
      ii        = il  +1
      do i=2,il
         ii        = ii  -1
         rev(i,1)  = rev(ii,2)
      end do
c
c     set boundary values along the profile
c
      i1        = itl  +1
      do i=i1,itu
         rev(i,1)  = -rev(i,2)
      end do
c
c     set boundary values at j=je
c
      do i=2,il
         rev(i,je) =  rev(i,jl)
      end do
c
c     set boundary values at i=1 and i=ie
c
      do j=1,je
         rev(1,j)  = rev(2,j)
         rev(ie,j) = rev(il,j)
      end do

      return

      end
