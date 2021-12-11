      subroutine bcfar(il, jl, ie, je, itl, itu,
     & w, p, rlv, rev, 
     & x, xc,
     & cp, cf,
     & gamma,rm,rho0,p0,h0,c0,u0,v0,ca,sa,re,prn,prt,scal,chord,xm,
     & ym,kvis,
     & bc,
     & mode)
c
c     ******************************************************************
c     *                                                                *
c     *   far field boundary condition                                 *
c     *                                                                *
c     ******************************************************************
c
c     w(i,j,1)  = density
c     w(i,j,2)  = momentum in x direction
c     w(i,j,3)  = momentum in y direction
c     w(i,j,4)  = total energy
c
c     ******************************************************************
c
C      use dims
c
c     ******************************************************************
c
c      use flo_var
C      use mesh_var
C      use out_var
c
c     ******************************************************************
c
c      use flo_param
c      use solv_param
c
c     ******************************************************************
c     dims
      integer, intent(in) :: il, jl, ie, je, itl, itu

c     flo_var
      real(8), intent(inout), dimension(:,:,:) :: w
      real(8), intent(inout), dimension(:,:)   :: p
      real(8), intent(inout), dimension(:,:)   :: rlv, rev

c     mesh_var
      real(8), intent(in), dimension(:,:,:) :: x,xc

c     out_var
      real(8), intent(inout), dimension(:)     :: cp
      real(8), intent(inout), dimension(:)     :: cf

c     flo_param
      real(8), intent(in)      :: gamma,rm,rho0,p0,h0,c0,u0,v0,ca,sa
      real(8), intent(in)      :: re,prn,prt
      real(8), intent(in)      :: scal,chord,xm,ym

      integer, intent(in)   :: kvis

c     solv_param
      real(8), intent(in)      :: bc

c     mg_param
      integer, intent(in)   :: mode

      implicit none
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j,l
c
c     ******************************************************************
c
      real(8)     :: pi,gm,gmg,gmm,s0,x0,beta,circ
      real(8)     :: xx,yx,xy,yy,d,xa,ya,r,angl,c,s
      real(8)     :: qn,qv,u,v,ufr,vfr,cfr
      real(8)     :: er,fr,qt,qq,cc,a,b,t
      real(8)     :: rmax2,omega,omega1,rq2

c
c     ******************************************************************
c
      real(8)     :: dx,dy,dcl,dcd
      real(8)     :: gm1,gpm,scf,clvis,cdvis
      real(8)     :: rho1,ru1,rv1,p1,u1,v1,t1
      real(8)     :: rho2,ru2,rv2,p2,u2,v2,t2
      real(8)     :: dxi,dxj,dyi,dyj,dui,duj,dvi,dvj,dsj
      real(8)     :: dux,duy,dvx,dvy
      real(8)     :: rlva,reva,rmu,rlam,rk,term
      real(8)     :: sigx,sigy,tauxy
c
c     ******************************************************************
c
      real(8), dimension(il)               :: gs2,gs3,gs4

      real(8) :: cl,cd,cm,cdv,clv
c
c     ******************************************************************
c
      pi        = 4.*atan(1.)
      gmg       = (gamma  -1.)/gamma
      gm        = 1./(gamma  -1.)
      gmm       = gamma*rm**2
      a         = (2.*gm  +1.)*gm
      s0        = rho0**gamma/p0
c
c     calculate the pressure and lift coefficients
c
      do i=2,il
         cp(i)     = ((p(i,2)  +p(i,1))/p0  -2.)/gmm
      end do

      cl        = 0.
      cd        = 0.
      cm        = 0.
      gmm       = gamma*rm**2

      do i=2,il
         cp(i)     = ((p(i,2)  +p(i,1))/p0  -2.)/gmm
      end do

      do i=itl+1,itu
         dx        = scal*(x(i,1,1)  -x(i-1,1,1))/chord
         dy        = scal*(x(i,1,2)  -x(i-1,1,2))/chord
c
c     this seems to be a bug.
c     when generating mesh xm is based on input xp.
c     when reading mesh xm is based on x & y.
c     i modified forcf.f gmesh.f st. xm is based on x & y.
c
c        xa        = (.5*scal*(x(i,1,1)  +x(i-1,1,1))  -xm)/chord
c        ya        = (.5*scal*(x(i,1,2)  +x(i-1,1,2))  -ym)/chord
c
         xa        = scal*(.5*(x(i,1,1)  +x(i-1,1,1))  -xm)/chord
         ya        = scal*(.5*(x(i,1,2)  +x(i-1,1,2))  -ym)/chord
         dcl       = -cp(i)*dx
         dcd       = cp(i)*dy
         cl        = cl  +dcl
         cd        = cd  +dcd
         cm        = cm  +dcd*ya  -dcl*xa
      end do

      if (kvis.eq.0) then
         clvis     = 0.
         cdvis     = 0.
      else
c
c     evaluation of viscous corrections to cl and cd.
c
      clvis   = 0.
      cdvis   = 0.

      gm1       = gamma - 1.
      gpm       = 2./(gamma*rm**2)
      scf       = sqrt(gamma)*rm/(scal*re/chord)

      do j=1,2
      do i=1,ie
         u(i,j)    = w(i,j,2)/w(i,j,1)
         v(i,j)    = w(i,j,3)/w(i,j,1)
c        t(i,j)     = p(i,j)/(gm1*w(i,j,1))
      end do
      end do

      if (mode.ne.0) then
         do i=itl+1,itu
            rev(i,1)  = -rev(i,2)
            u(i,1)    = -u(i,2)
            v(i,1)    = -v(i,2)
         end do
      end if
c
c     evaluation of viscous flux terms in i and j directions
c
      do 10 i=itl+1,itu

      rho1    = .25*(w(i,1,1)  +w(i+1,1,1)  +w(i,2,1)  +w(i+1,2,1))
      ru1     = .25*(w(i,1,2)  +w(i+1,1,2)  +w(i,2,2)  +w(i+1,2,2))
      rv1     = .25*(w(i,1,3)  +w(i+1,1,3)  +w(i,2,3)  +w(i+1,2,3))
      p1      = .25*(p(i,1)    +p(i+1,1)    +p(i,2)    +p(i+1,2))
      u1      = ru1/rho1
      v1      = rv1/rho1
c     t1      = p1/(gm1*rho1)

      rho2    = .25*(w(i-1,1,1)  +w(i,1,1)  +w(i-1,2,1)  +w(i,2,1))
      ru2     = .25*(w(i-1,1,2)  +w(i,1,2)  +w(i-1,2,2)  +w(i,2,2))
      rv2     = .25*(w(i-1,1,3)  +w(i,1,3)  +w(i-1,2,3)  +w(i,2,3))
      p2      = .25*(p(i-1,1)    +p(i,1)    +p(i-1,2)    +p(i,2))
      u2      = ru2/rho2
      v2      = rv2/rho2
c     t2      = p2/(gm1*rho2)

      dxi       = x(i,1,1)  -x(i-1,1,1)
      dyi       = x(i,1,2)  -x(i-1,1,2)
      dxj       = xc(i,2,1)  -xc(i,1,1)
      dyj       = xc(i,2,2)  -xc(i,1,2)
      dsj       = 1./(dxi*dyj  -dyi*dxj)

      dui       = u1  -u2
      dvi       = v1  -v2
c     dti       = t1  -t2
      duj       = u(i,2)  -u(i,1)
      dvj       = v(i,2)  -v(i,1)
c     dtj       = t(i,2)  -t(i,1)

      dux       = (dui*dyj  -duj*dyi)*dsj
      dvx       = (dvi*dyj  -dvj*dyi)*dsj
c     dtx       = (dti*dyj  -dtj*dyi)*dsj

      duy       = (duj*dxi  -dui*dxj)*dsj
      dvy       = (dvj*dxi  -dvi*dxj)*dsj
c     dty       = (dtj*dxi  -dti*dxj)*dsj

      rlva      = .5*(rlv(i,2)  +rlv(i,1))
      reva      = .5*(rev(i,2)  +rev(i,1))
      rmu       = rlva  +reva
      rk        = gamma*(rlva/prn  +reva/prt)
      rlam      = -2.*rmu/3.

      term      = rlam*(dux  +dvy)
      sigx      = term  +2.*rmu*dux
      sigy      = term  +2.*rmu*dvy
      tauxy     = rmu*(duy  +dvx)
c     qx        = -rk*dtx
c     qy        = -rk*dty

c
c     residuals in "nsflux.f" have extra 2, compensated later in "euler.f"
c     so, for viscous stresses in lines below:                   -2. ==> -1.
c     also, gs2 and gs3 -are forces on the body (not fluid), so: -1. ==>  1.
c
      gs2(i)    = (-sigx*dyi  +tauxy*dxi)*scf*gpm
      gs3(i)    = ( sigy*dxi  -tauxy*dyi)*scf*gpm

c     gs2(i)    = (-sigx*dyi  +tauxy*dxi)*scf
c     gs3(i)    = ( sigy*dxi  -tauxy*dyi)*scf

      cdvis     = cdvis  +gs2(i)
      clvis     = clvis  +gs3(i)

      cf(i)     = tauxy*scf*gpm

   10 continue

      end if

      dcl       = cl*ca  -cd*sa
      cd        = cl*sa  +cd*ca
      cl        = dcl
      cdv       = clvis*sa  +cdvis*ca
      clv       = clvis*ca  -cdvis*sa

c     cd        = cd +clvis*sa  +cdvis*ca
c     cl        = cl +clvis*ca  -cdvis*sa

      return

      end


      if (rm.ge.1.) go to 31
c
c     calculate the vortex contribution to the far field
c
      x0        = x(itl,1,1)  -.5*chord/scal
      beta      = sqrt(1.  -rm**2)
      circ      = .25*chord*cl*rm*c0*beta/(pi*scal)

      if (bc.ge.0.) go to 11
c
c     use simple extrapolation if bc.lt.0.
c
      do i=2,il
         xx        = x(i,jl,1)  -x(i-1,jl,1)
         yx        = x(i,jl,2)  -x(i-1,jl,2)
         qn        = xx*w(i,jl,3)  -yx*w(i,jl,2)
         if (qn.le.0.) then
            xa        = .5*(x(i,jl,1)  +x(i-1,jl,1))  -x0
            ya        = .5*(x(i,jl,2)  +x(i-1,jl,2))
            r         = sqrt(xa**2  +ya**2)
            angl      = atan2(ya,xa)
            c         = cos(angl)
            s         = sin(angl)
            qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
            ufr       = u0  +qv*s
            vfr       = v0  -qv*c
            p(i,je)   = p(i,jl)
            w(i,je,1) = p(i,je)/(gmg*(h0  -.5*(ufr**2  +vfr**2)))
            w(i,je,2) = w(i,je,1)*ufr
            w(i,je,3) = w(i,je,1)*vfr
            w(i,je,4) = w(i,je,1)*h0  -p(i,je)
         else
            u         = w(i,jl,2)/w(i,jl,1)
            v         = w(i,jl,3)/w(i,jl,1)
            w(i,je,1) = p0/(gmg*(h0  -.5*(u**2  +v**2)))
            w(i,je,2) = w(i,je,1)*u
            w(i,je,3) = w(i,je,1)*v
            w(i,je,4) = w(i,je,1)*h0  -p0
         end if
      end do

      do j=2,jl
         xx        = x(1,j,1)  -x(1,j-1,1)
         yx        = x(1,j,2)  -x(1,j-1,2)
         qn        = xx*w(2,j,3)  -yx*w(2,j,2)
         if (qn.le.0.) then
            xa        = .5*(x(1,j,1)  +x(1,j-1,1))  -x0
            ya        = .5*(x(1,j,2)  +x(1,j-1,2))
            r         = sqrt(xa**2  +ya**2)
            angl      = atan2(ya,xa)
            c         = cos(angl)
            s         = sin(angl)
            qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
            ufr       = u0  +qv*s
            vfr       = v0  -qv*c
            p(1,j)   = p(2,j)
            w(1,j,1) = p(2,j)/(gmg*(h0  -.5*(ufr**2  +vfr**2)))
            w(1,j,2) = w(1,j,1)*ufr
            w(1,j,3) = w(1,j,1)*vfr
            w(1,j,4) = w(1,j,1)*h0  -p(1,j)
         else
            u        = w(2,j,2)/w(2,j,1)
            v        = w(2,j,3)/w(2,j,1)
            w(1,j,1) = p0/(gmg*(h0  -.5*(u**2  +v**2)))
            w(1,j,2) = w(1,j,1)*u
            w(1,j,3) = w(1,j,1)*v
            w(1,j,4) = w(1,j,1)*h0  -p0
         end if
      end do

      do j=jl,2,-1
         xx        = x(il,j-1,1)  -x(il,j,1)
         yx        = x(il,j-1,2)  -x(il,j,2)
         qn        = xx*w(2,j,3)  -yx*w(2,j,2)
         if (qn.le.0.) then
            xa        = .5*(x(il,j,1)  +x(il,j-1,1))  -x0
            ya        = .5*(x(il,j,2)  +x(il,j-1,2))
            r         = sqrt(xa**2  +ya**2)
            angl      = atan2(ya,xa)
            c         = cos(angl)
            s         = sin(angl)
            qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
            ufr       = u0  +qv*s
            vfr       = v0  -qv*c
            p(ie,j)   = p(il,j)
            w(ie,j,1) = p(il,j)/(gmg*(h0  -.5*(ufr**2  +vfr**2)))
            w(ie,j,2) = w(ie,j,1)*ufr
            w(ie,j,3) = w(ie,j,1)*vfr
            w(ie,j,4) = w(ie,j,1)*h0  -p(ie,j)
         else
            u        = w(il,j,2)/w(il,j,1)
            v        = w(il,j,3)/w(il,j,1)
            w(ie,j,1) = p0/(gmg*(h0  -.5*(u**2  +v**2)))
            w(ie,j,2) = w(ie,j,1)*u
            w(ie,j,3) = w(ie,j,1)*v
            w(ie,j,4) = w(ie,j,1)*h0  -p0
         end if
      end do

      go to 41
c
c     extrapolate the outgoing riemann invariant normal to the boundary
c     and fix the incoming riemann invariant if bc.eq.0.
c
   11 if (bc.gt.0.) go to 21

      do i=2,il
         xa        = .5*(x(i,jl,1)  +x(i-1,jl,1))  -x0
         ya        = .5*(x(i,jl,2)  +x(i-1,jl,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(i,jl,1)  -x(i-1,jl,1)
         yx        = x(i,jl,2)  -x(i-1,jl,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(i,jl,2)/w(i,jl,1)
         v         = w(i,jl,3)/w(i,jl,1)
         c         = sqrt(gamma*p(i,jl)/w(i,jl,1))
         er        = xx*v   -yx*u   +2.*gm*c
         fr        = xx*vfr  -yx*ufr  -2.*gm*cfr
         c         = .25*(er  -fr)/gm
         qn        = .5*(er  +fr)
         qt        = xx*ufr  +yx*vfr
         s         = s0
         if (qn.gt.0.) then
            qt        = xx*u   +yx*v
            s         = w(i,jl,1)**gamma/p(i,jl)
         end if
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(i,je,1) = (s*cc)**gm
         w(i,je,2) = w(i,je,1)*u
         w(i,je,3) = w(i,je,1)*v
         p(i,je)   = w(i,je,1)*cc
         w(i,je,4) = w(i,je,1)*h0  -p(i,je)
      end do

      do j=2,jl
         xa        = .5*(x(1,j,1)  +x(1,j-1,1))  -x0
         ya        = .5*(x(1,j,2)  +x(1,j-1,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(1,j,1)  -x(1,j-1,1)
         yx        = x(1,j,2)  -x(1,j-1,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(2,j,2)/w(2,j,1)
         v         = w(2,j,3)/w(2,j,1)
         c         = sqrt(gamma*p(2,j)/w(2,j,1))
         er        = xx*v   -yx*u   +2.*gm*c
         fr        = xx*vfr  -yx*ufr  -2.*gm*cfr
         c         = .25*(er  -fr)/gm
         qn        = .5*(er  +fr)
         qt        = xx*ufr  +yx*vfr
         s         = s0
         if (qn.gt.0.) then
            qt        = xx*u   +yx*v
            s         = w(2,j,1)**gamma/p(2,j)
         end if
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(1,j,1) = (s*cc)**gm
         w(1,j,2) = w(1,j,1)*u
         w(1,j,3) = w(1,j,1)*v
         p(1,j)   = w(1,j,1)*cc
         w(1,j,4) = w(1,j,1)*h0  -p(1,j)
      end do

      do j=jl,2,-1
         xa        = .5*(x(il,j,1)  +x(il,j-1,1))  -x0
         ya        = .5*(x(il,j,2)  +x(il,j-1,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(il,j-1,1)  -x(il,j,1)
         yx        = x(il,j-1,2)  -x(il,j,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(il,j,2)/w(il,j,1)
         v         = w(il,j,3)/w(il,j,1)
         c         = sqrt(gamma*p(il,j)/w(il,j,1))
         er        = xx*v   -yx*u   +2.*gm*c
         fr        = xx*vfr  -yx*ufr  -2.*gm*cfr
         c         = .25*(er  -fr)/gm
         qn        = .5*(er  +fr)
         qt        = xx*ufr  +yx*vfr
         s         = s0
         if (qn.gt.0.) then
            qt        = xx*u   +yx*v
            s         = w(il,j,1)**gamma/p(il,j)
         end if
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(ie,j,1) = (s*cc)**gm
         w(ie,j,2) = w(ie,j,1)*u
         w(ie,j,3) = w(ie,j,1)*v
         p(ie,j)   = w(ie,j,1)*cc
         w(ie,j,4) = w(ie,j,1)*h0  -p(ie,j)
      end do

      go to 41
c
c     extrapolate the upstream riemann invariant normal to the boundary
c     and fix the total enthalpy if bc.gt.0.
c
   21 do i=2,il
         xa        = .5*(x(i,jl,1)  +x(i-1,jl,1))  -x0
         ya        = .5*(x(i,jl,2)  +x(i-1,jl,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(i,jl,1)  -x(i-1,jl,1)
         yx        = x(i,jl,2)  -x(i-1,jl,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(i,jl,2)/w(i,jl,1)
         v         = w(i,jl,3)/w(i,jl,1)
         c         = sqrt(gamma*p(i,jl)/w(i,jl,1))
         qn        = xx*v   -yx*u
         s         = s0
         qt        = xx*ufr  +yx*vfr
         r         = xx*v   -yx*u  +2.*gm*c
         t         = -1.
         if (qn.gt.0.) then
            s         = w(i,jl,1)**gamma/p(i,jl)
            qt        = xx*u   +yx*v
            r         = xx*vfr  -yx*ufr  -2.*gm*cfr
            t         = 1.
         end if
         b         = -t*gm*r
         c         = .5*(r**2  +qt**2)  -h0
         c         = (sqrt(b**2  -a*c)  +b)/a
         qn        = r  +2.*t*gm*c
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(i,je,1) = (s*cc)**gm
         w(i,je,2) = w(i,je,1)*u
         w(i,je,3) = w(i,je,1)*v
         p(i,je)   = w(i,je,1)*cc
         w(i,je,4) = w(i,je,1)*h0  -p(i,je)
      end do

      do j=2,jl
         xa        = .5*(x(1,j,1)  +x(1,j-1,1))  -x0
         ya        = .5*(x(1,j,2)  +x(1,j-1,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(1,j,1)  -x(1,j-1,1)
         yx        = x(1,j,2)  -x(1,j-1,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(2,j,2)/w(2,j,1)
         v         = w(2,j,3)/w(2,j,1)
         c         = sqrt(gamma*p(2,j)/w(2,j,1))
         qn        = xx*v   -yx*u
         s         = s0
         qt        = xx*ufr  +yx*vfr
         r         = xx*v   -yx*u  +2.*gm*c
         t         = -1.
         if (qn.gt.0.) then
            s         = w(2,j,1)**gamma/p(2,j)
            qt        = xx*u   +yx*v
            r         = xx*vfr  -yx*ufr  -2.*gm*cfr
            t         = 1.
         end if
         b         = -t*gm*r
         c         = .5*(r**2  +qt**2)  -h0
         c         = (sqrt(b**2  -a*c)  +b)/a
         qn        = r  +2.*t*gm*c
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(1,j,1) = (s*cc)**gm
         w(1,j,2) = w(1,j,1)*u
         w(1,j,3) = w(1,j,1)*v
         p(1,j)   = w(1,j,1)*cc
         w(1,j,4) = w(1,j,1)*h0  -p(1,j)
      end do

      do j=jl,2,-1
         xa        = .5*(x(il,j,1)  +x(il,j-1,1))  -x0
         ya        = .5*(x(il,j,2)  +x(il,j-1,2))
         r         = sqrt(xa**2  +ya**2)
         angl      = atan2(ya,xa)
         c         = cos(angl)
         s         = sin(angl)
         qv        = circ/(r*(1.  -(rm*(s*ca  -c*sa))**2))
         ufr       = u0  +qv*s
         vfr       = v0  -qv*c
         cfr       = sqrt((gamma  -1.)*(h0  -.5*(ufr**2  +vfr**2)))
         xx        = x(il,j-1,1)  -x(il,j,1)
         yx        = x(il,j-1,2)  -x(il,j,2)
         d         = 1./sqrt(xx**2  +yx**2)
         xx        = xx*d
         yx        = yx*d
         u         = w(il,j,2)/w(il,j,1)
         v         = w(il,j,3)/w(il,j,1)
         c         = sqrt(gamma*p(il,j)/w(il,j,1))
         qn        = xx*v   -yx*u
         s         = s0
         qt        = xx*ufr  +yx*vfr
         r         = xx*v   -yx*u  +2.*gm*c
         t         = -1.
         if (qn.gt.0.) then
            s         = w(il,j,1)**gamma/p(il,j)
            qt        = xx*u   +yx*v
            r         = xx*vfr  -yx*ufr  -2.*gm*cfr
            t         = 1.
         end if
         b         = -t*gm*r
         c         = .5*(r**2  +qt**2)  -h0
         c         = (sqrt(b**2  -a*c)  +b)/a
         qn        = r  +2.*t*gm*c
         u         = xx*qt  -yx*qn
         v         = xx*qn  +yx*qt
         cc        = c**2/gamma
         w(ie,j,1) = (s*cc)**gm
         w(ie,j,2) = w(ie,j,1)*u
         w(ie,j,3) = w(ie,j,1)*v
         p(ie,j)   = w(ie,j,1)*cc
         w(ie,j,4) = w(ie,j,1)*h0  -p(ie,j)
      end do

      go to 41
c
c     extrapolate all quantities at the outflow boundary
c     and fix all quantities at the inflow boundary
c     if the flow is supersonic
c
   31 do i=2,il
         xx        = x(i,jl,1)  -x(i-1,jl,1)
         yx        = x(i,jl,2)  -x(i-1,jl,2)
         qn        = xx*v0  -yx*u0
         if (qn.gt.0.) then
            w(i,je,1) = w(i,jl,1)
            w(i,je,2) = w(i,jl,2)
            w(i,je,3) = w(i,jl,3)
            w(i,je,4) = w(i,jl,4)
            p(i,je)   = p(i,jl)
         else
            w(i,je,1) = rho0
            w(i,je,2) = rho0*u0
            w(i,je,3) = rho0*v0
            w(i,je,4) = rho0*h0  -p0
            p(i,je)   = p0
         end if
      end do

      do j=2,jl

         w(1,j,1) = w(2,j,1)
         w(1,j,2) = w(2,j,2)
         w(1,j,3) = w(2,j,3)
         w(1,j,4) = w(2,j,4)
         p(1,j)   = p(2,j)

         w(ie,j,1) = w(il,j,1)
         w(ie,j,2) = w(il,j,2)
         w(ie,j,3) = w(il,j,3)
         w(ie,j,4) = w(il,j,4)
         p(ie,j)   = p(il,j)

      end do
c
c     outflow bc for viscous calculations
c
   41 if (kvis.gt.0) then

      do j=2,jl

         w(1,j,1)  = 2.*w(2,j,1) - w(3,j,1)
         w(1,j,2)  = 2.*w(2,j,2) - w(3,j,2)
         w(1,j,3)  = 2.*w(2,j,3) - w(3,j,3)
         rmax2     = w(2,j,2)**2/(gamma*p(2,j)*w(2,j,1))
         omega1    = gamma*rmax2/(1.+(gamma-1)*rmax2)
         omega     = min(1.,omega1)
         p(1,j)    = omega*(2.*p(2,j)-p(3,j))+(1.-omega)*p(2,j)
         rq2       = w(1,j,2)**2+w(1,j,3)**2
         w(1,j,4)  = p(1,j)*gm + .5*rq2/w(1,j,1)

         w(ie,j,1) = 2.*w(il,j,1) - w(il-1,j,1)
         w(ie,j,2) = 2.*w(il,j,2) - w(il-1,j,2)
         w(ie,j,3) = 2.*w(il,j,3) - w(il-1,j,3)
         rmax2     = w(il,j,2)**2/(gamma*p(il,j)*w(il,j,1))
         omega1    = gamma*rmax2/(1.+(gamma-1)*rmax2)
         omega     = min(1.,omega1)
         p(ie,j)   = omega*(2.*p(il,j)-p(il-1,j))+(1.-omega)*p(il,j)
         rq2       = w(ie,j,2)**2+w(ie,j,3)**2
         w(ie,j,4) = p(ie,j)*gm + .5*rq2/w(ie,j,1)

      end do

      end if
c
c     enforce bc at corner cells
c
      w(1,je,1)   = w(2,je,1)
      w(1,je,2)   = w(2,je,2)
      w(1,je,3)   = w(2,je,3)
      w(1,je,4)   = w(2,je,4)
      p(1,je)     = p(2,je)

      w(ie,je,1)  = w(il,je,1)
      w(ie,je,2)  = w(il,je,2)
      w(ie,je,3)  = w(il,je,3)
      w(ie,je,4)  = w(il,je,4)
      p(ie,je)    = p(il,je)

      w(1,1,1)    = w(ie,2,1)
      w(1,1,2)    = w(ie,2,2)
      w(1,1,3)    = w(ie,2,3)
      w(1,1,4)    = w(ie,2,4)
      p(1,1)      = p(ie,2)

      w(ie,1,1)   = w(1,2,1)
      w(ie,1,2)   = w(1,2,2)
      w(ie,1,3)   = w(1,2,3)
      w(ie,1,4)   = w(1,2,4)
      p(ie,1)     = p(1,2)

      end
