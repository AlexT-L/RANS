      subroutine bcfar
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
      use dims
c
c     ******************************************************************
c
      use flo_var
      use mesh_var
      use out_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
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
c
c     ******************************************************************
c
      real     :: pi,gm,gmg,gmm,s0,x0,beta,circ
      real     :: xx,yx,xy,yy,d,xa,ya,r,angl,c,s
      real     :: qn,qv,u,v,ufr,vfr,cfr
      real     :: er,fr,qt,qq,cc,a,b,t
      real     :: rmax2,omega,omega1,rq2
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

      call forcf

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
