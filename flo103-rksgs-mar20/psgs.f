      subroutine  psgs
c
c     Symmetric Gauss-Seidel preconditioner
c
c     Programmed by Antony Jameson January 2013
c     using the transformation to the symmetric Euler Jacobian matrix
c     alternative formulation to the scheme proposed by Cord Rossow,
c     Charlie Swanson and Eli Turkel
c
      use dims
c
c     ******************************************************************
c
      use flo_var
      use solv_var
      use mesh_var
      use psm_var
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
      real, dimension(4,4)      :: xx
      real, dimension(4)        :: rhs
c
c     ******************************************************************
c
      integer  :: i,j,ia,ibf,ja,jbf,idir,jdir,m,n
      integer  :: iter,niter
c
c     ******************************************************************
c
      real     :: rhoa,sar,sar2,svol,dtv2
      real     :: rmuel,rmuet,sgrmrei,scale
      real     :: sxa,sya,vnx,vny
      real     :: rl,pl,ul,vl,rr,pr,ur,vr
      real     :: rl1,rr1,rb1,rt1,rl2,rr2,rb2,rt2
      real     :: rl3,rr3,rb3,rt3,rl4,rr4,rb4,rt4
      real     :: a1,a2,a3,a4
c
c     ******************************************************************
c
      sgrmrei   = 0.
      if (kvis.gt.0) sgrmrei = sqrt(gamma)*rm/re
c
      niter = 2
c
      if (nstage .eq. 1) then
c
c     save the velocities,density and pressure
c
      do j=1,je
      do i=1,ie
        r00(i,j) = w(i,j,1)
        u00(i,j) = w(i,j,2)/w(i,j,1)
        v00(i,j) = w(i,j,3)/w(i,j,1)
        p00(i,j) = p(i,j)
      end do
      end do
c
      end if
c
c     ******************************************************************
c     *                                                                *
c     *   setup and solve linear implicit system                       *
c     *                                                                *
c     ******************************************************************
c
c     transform the residuals to the symmetrizing variables
c
      call res_cons_to_sym
c
c     set the residuals to zero at the boundaries
c
      call bcdw
c
c     set the modified residuals to zero
c
      do j=1,je
      do i=1,ie
        rs(i,j,1) = 0.
        rs(i,j,2) = 0.
        rs(i,j,3) = 0.
        rs(i,j,4) = 0.
      end do
      end do
c
c     ******************************************************************
c     *                                                                *
c     *   solution by symmetric gauss-seidel iteration                 *
c     *                                                                *
c     ******************************************************************
c
      idir =  -1
c
      do 500 iter = 1,niter
c
      idir = -idir
      jdir = idir
c
      call bcrs
c
      if (idir.eq.1) then
        ia  = 2
        ja  = 2
        ibf = il
        jbf = jl
      else
        ia  = il
        ja  = jl
        ibf = 2
        jbf = 2
      end if
c
      do 400 j=ja,jbf,jdir
      do 300 i=ia,ibf,idir
c
        do m=1,4
        do n=1,4
          xx(m,n) = 0.
        end do
        end do
c
        dtv2   = eps*cfl*dtl(i,j)
c
c       left variables
c
        rl     = r00(i,j)
        pl     = p00(i,j)
        ul     = u00(i,j)
        vl     = v00(i,j)
c
c       left i interface
c
        a1     = rs(i-1,j,1)
        a2     = rs(i-1,j,2)
        a3     = rs(i-1,j,3)
        a4     = rs(i-1,j,4)
c
        sxa    = x(i-1,j-1,2)   -x(i-1,j,2)
        sya    = x(i-1,j,1)     -x(i-1,j-1,1)
        sar2   = sxa*sxa  +sya*sya
        sar    = sqrt(sar2)
        vnx    = sxa/sar
        vny    = sya/sar
c
c       right variables
c
        rr     = r00(i-1,j)
        pr     = p00(i-1,j)
        ur     = u00(i-1,j)
        vr     = v00(i-1,j)
c
c       viscous coefficients
c
        rhoa   = .5*(rl  +rr)
        svol   = .5*(vol(i-1,j)  +vol(i,j))
        scale  = sgrmrei*sar2/(svol*rhoa)
        if (kvis.eq.0) scale = 0.
        rmuel  = .5*scale*(rlv(i-1,j)  +rlv(i,j))
        rmuet  = .5*scale*(rev(i-1,j)  +rev(i,j))
c
        call rhs_face(vnx, vny, sar,
     #                dtv2, rmuel, rmuet,
     #                rl, pl, ul, vl,
     #                rr, pr, ur, vr,
     #                a1, a2, a3, a4,
     #                rl1, rl2, rl3, rl4,
     #                xx)
c
c       right i interface
c
        a1     = rs(i+1,j,1)
        a2     = rs(i+1,j,2)
        a3     = rs(i+1,j,3)
        a4     = rs(i+1,j,4)
c
        sxa    = x(i,j,2)    -x(i,j-1,2)
        sya    = x(i,j-1,1)  -x(i,j,1)
        sar2   = sxa*sxa  +sya*sya
        sar    = sqrt(sar2)
        vnx    = sxa/sar
        vny    = sya/sar
c
c       right variables
c
        rr     = r00(i+1,j)
        pr     = p00(i+1,j)
        ur     = u00(i+1,j)
        vr     = v00(i+1,j)
c
c       viscous coefficients
c
        rhoa   = .5*(rl  +rr)
        svol   = .5*(vol(i,j)  +vol(i+1,j))
        scale  = sgrmrei*sar2/(svol*rhoa)
        if (kvis.eq.0) scale = 0.
        rmuel  = .5*scale*(rlv(i,j)  +rlv(i+1,j))
        rmuet  = .5*scale*(rev(i,j)  +rev(i+1,j))
c
        call rhs_face(vnx, vny, sar,
     #                dtv2, rmuel, rmuet,
     #                rl, pl, ul, vl,
     #                rr, pr, ur, vr,
     #                a1, a2, a3, a4,
     #                rr1, rr2, rr3, rr4,
     #                xx)
c
c       bottom j interface
c
        a1     = rs(i,j-1,1)
        a2     = rs(i,j-1,2)
        a3     = rs(i,j-1,3)
        a4     = rs(i,j-1,4)
c
        sxa    = x(i,j-1,2)    -x(i-1,j-1,2)
        sya    = x(i-1,j-1,1)  -x(i,j-1,1)
        sar2   = sxa*sxa  +sya*sya
        sar    = sqrt(sar2)
        vnx    = sxa/sar
        vny    = sya/sar
c
c       right variables
c
        rr     = r00(i,j-1)
        pr     = p00(i,j-1)
        ur     = u00(i,j-1)
        vr     = v00(i,j-1)
c
c       viscous coefficients
c
        rhoa   = .5*(rl  +rr)
        svol   = .5*(vol(i,j-1)  +vol(i,j))
        scale  = sgrmrei*sar2/(svol*rhoa)
        if (kvis.eq.0) scale = 0.
        rmuel  = .5*scale*(rlv(i,j-1)  +rlv(i,j))
        rmuet  = .5*scale*(rev(i,j-1)  +rev(i,j))
c
        call rhs_face(vnx, vny, sar,
     #                dtv2, rmuel, rmuet,
     #                rl, pl, ul, vl,
     #                rr, pr, ur, vr,
     #                a1, a2, a3, a4,
     #                rb1, rb2, rb3, rb4,
     #                xx)
c
c       top j interface
c
        a1     = rs(i,j+1,1)
        a2     = rs(i,j+1,2)
        a3     = rs(i,j+1,3)
        a4     = rs(i,j+1,4)
c
        sxa    = x(i-1,j,2)  -x(i,j,2)
        sya    = x(i,j,1)    -x(i-1,j,1)
        sar2   = sxa*sxa  +sya*sya
        sar    = sqrt(sar2)
        vnx    = sxa/sar
        vny    = sya/sar
c
c       right variables
c
        rr     = r00(i,j+1)
        pr     = p00(i,j+1)
        ur     = u00(i,j+1)
        vr     = v00(i,j+1)
c
c       viscous coefficients
c
        rhoa   = .5*(rl  +rr)
        svol   = .5*(vol(i,j)  +vol(i,j+1))
        scale  = sgrmrei*sar2/(svol*rhoa)
        if (kvis.eq.0) scale = 0.
        rmuel  = .5*scale*(rlv(i,j)  +rlv(i,j+1))
        rmuet  = .5*scale*(rev(i,j)  +rev(i,j+1))
c
        call rhs_face(vnx, vny, sar,
     #                dtv2, rmuel, rmuet,
     #                rl, pl, ul, vl,
     #                rr, pr, ur, vr,
     #                a1, a2, a3, a4,
     #                rt1, rt2, rt3, rt4,
     #                xx)
c
c************assemble  complete right-hand-side
c
        rhs(1) = dw(i,j,1) - rl1 - rr1 - rb1 - rt1
        rhs(2) = dw(i,j,2) - rl2 - rr2 - rb2 - rt2
        rhs(3) = dw(i,j,3) - rl3 - rr3 - rb3 - rt3
        rhs(4) = dw(i,j,4) - rl4 - rr4 - rb4 - rt4

c
c--------*----------*---------*---------*---------*---------*---------*---------*
c
        xx(1,1) = 1. + xx(1,1)
        xx(2,2) = 1. + xx(2,2)
        xx(3,3) = 1. + xx(3,3)
        xx(4,4) = 1. + xx(4,4)
c
        call solve_rs(rhs,xx)
c
        rs(i,j,1) = rhs(1)
        rs(i,j,2) = rhs(2)
        rs(i,j,3) = rhs(3)
        rs(i,j,4) = rhs(4)
c
  300 continue
  400 continue
c
  500 continue                ! end gs
c
c     transform the residuals back to the conservative variables
c
      call res_sym_to_cons

      return

      end
c
c     ******************************************************************
c
      subroutine  rhs_face(vnx, vny, sar,
     #                     dtv2, rmuel, rmuet,
     #                     rl, pl, ul, vl,
     #                     rr, pr, ur, vr,
     #                     a1, a2, a3, a4,
     #                     rhs_f1, rhs_f2, rhs_f3, rhs_f4,
     #                     xx)
c
c     calculate the face contributions to the right hand side
c     and preconditioning matrix
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
      real, dimension(4,4)      :: xx
c
c     ******************************************************************
c
      real     :: gm1,dlim
      real     :: vnx,vny,sar,dtv2,rmuel,rmuet
      real     :: rl,pl,ul,vl,rr,pr,ur,vr,a1,a2,a3,a4
      real     :: rhoa,ua,va,pa,qn,cc,c
      real     :: e1,e2,e3,a,q1,q2,q3,r1,r2,aqsa
      real     :: epsv1,epsv2,amu,amupr,vq1,vr1
      real     :: a11,a12,a13,a14,a22,a23,a33,a41,a44
      real     :: rhs_f1,rhs_f2,rhs_f3,rhs_f4
c
c     ******************************************************************
c
      gm1     = gamma - 1.
      dlim    = .10
      if (kvis.eq.0) dlim = .25
      if (ncyc.lt.5) dlim = .5
c
c     arithmetic averaging
c
      rhoa  = .5*(rl  +rr)
      ua    = .5*(ul  +ur)
      va    = .5*(vl  +vr)
      pa    = .5*(pl  +pr)
      qn    = vnx*ua  +vny*va
      cc    = gamma*pa/rhoa
      c     = sqrt(cc)
c
c     viscous coefficients
c
      epsv1 = 1.
      epsv2 = .25
      amu   = epsv1*dtv2*(rmuel + rmuet)
      amupr = epsv2*dtv2*(rmuel/prn + rmuet/prt)
c
c     eigenvalues
c
      e1    = abs(qn)
      e2    = abs(qn  +c)
      e3    = abs(qn  -c)
      a     = dlim*c
      if (e1.lt.a) e1 = .5*(a  +e1*e1/a)
      if (e2.lt.a) e2 = .5*(a  +e2*e2/a)
      if (e3.lt.a) e3 = .5*(a  +e3*e3/a)
      a     = .5*sar*dtv2
      aqsa  = a*diag*abs(qn)
c
c     form the right hand side with the negative jacobian
c
      q1    = a*(qn  -e1)
      q2    = a*(qn  +c  -e2)
      q3    = a*(qn  -c  -e3)
      r1    = .5*(q2  +q3)  -q1
      r2    = .5*(q2  -q3)

      vq1   = q1  -amu
      vr1   = r1  -amu/3.

      a11   = r1  +q1  -gm1*amupr  -aqsa
      a12   = vnx*r2
      a13   = vny*r2
      a14   = -amupr
      a22   = vnx*vnx*vr1  +vq1  -aqsa
      a23   = vnx*vny*vr1
      a33   = vny*vny*vr1  +vq1  -aqsa
      a41   = -gm1*amupr
      a44   = q1  -amupr  -aqsa

      rhs_f1   = a11*a1  +a12*a2  +a13*a3  +a14*a4
      rhs_f2   = a12*a1  +a22*a2  +a23*a3
      rhs_f3   = a13*a1  +a23*a2  +a33*a3
      rhs_f4   = a41*a1  +a44*a4
c
c     form the preconditioning matrix
c
      xx(1,1)  = xx(1,1)  -a11
      xx(1,2)  = xx(1,2)  -a12
      xx(1,3)  = xx(1,3)  -a13
      xx(1,4)  = xx(1,4)  -a14

      xx(2,1)  = xx(2,1)  -a12
      xx(2,2)  = xx(2,2)  -a22
      xx(2,3)  = xx(2,3)  -a23

      xx(3,1)  = xx(3,1)  -a13
      xx(3,2)  = xx(3,2)  -a23
      xx(3,3)  = xx(3,3)  -a33

      xx(4,1)  = xx(4,1)  -a14
      xx(4,4)  = xx(4,4)  -a44
 
      return

      end
c
c     ******************************************************************
c
      subroutine  bcdw
c
c     set the residuals to zero at the boundaries
c     and match the residuals across the cut in the c-mesh
c
      use dims
c
c     ******************************************************************
c
      use solv_var
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      integer  :: i,j,ii
c
c     ******************************************************************
c
      do 10 i=1,ie
        dw(i,1,1) = 0.0
        dw(i,1,2) = 0.0
        dw(i,1,3) = 0.0
        dw(i,1,4) = 0.0
c
        dw(i,je,1) = 0.0
        dw(i,je,2) = 0.0
        dw(i,je,3) = 0.0
        dw(i,je,4) = 0.0
   10 continue
c
      do 20 j=1,je
        dw(1,j,1) = 0.0
        dw(1,j,2) = 0.0
        dw(1,j,3) = 0.0
        dw(1,j,4) = 0.0
c
        dw(ie,j,1) = 0.0
        dw(ie,j,2) = 0.0
        dw(ie,j,3) = 0.0
        dw(ie,j,4) = 0.0
   20 continue
c
      do 30 i=1,ie
         if (i.le.itl .or. i.gt.itu) then
            ii    = ie-i+1
            dw(i,1,1) = dw(ii,2,1)
            dw(i,1,2) = dw(ii,2,2)
            dw(i,1,3) = dw(ii,2,3)
            dw(i,1,4) = dw(ii,2,4)
         end if
         if (i.gt.itl .or. i.lt.itl) then
            dw(i,1,1) = dw(i,2,1)
            dw(i,1,2) =-dw(i,2,2)
            dw(i,1,3) =-dw(i,2,3)
c           dw(i,1,2) = dw(i,2,2)
c           dw(i,1,3) = dw(i,2,3)
            dw(i,1,4) = dw(i,2,4)
         end if
   30 continue
c
      return
      end
c
c     ******************************************************************
c
      subroutine  bcrs
c
c     set the modified residuals to zero at the boundaries
c     and match the modified residuals across the cut in the c-mesh
c
      use dims
c
c     ******************************************************************
c
      use psm_var
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      integer  :: i,j,ii
c
c     ******************************************************************
c
      do 10 i=1,ie
        rs(i, 1,1) = 0.0
        rs(i, 1,2) = 0.0
        rs(i, 1,3) = 0.0
        rs(i, 1,4) = 0.0
c
        rs(i,je,1) = 0.0
        rs(i,je,2) = 0.0
        rs(i,je,3) = 0.0
        rs(i,je,4) = 0.0
   10 continue
c
      do 20 j=1,je
        rs( 1,j,1) = 0.0
        rs( 1,j,2) = 0.0
        rs( 1,j,3) = 0.0
        rs( 1,j,4) = 0.0
c
        rs(ie,j,1) = 0.0
        rs(ie,j,2) = 0.0
        rs(ie,j,3) = 0.0
        rs(ie,j,4) = 0.0
   20 continue
c
      do 30 i=1,ie
         if (i.le.itl .or. i.gt.itu) then
            ii    = ie-i+1
            rs(i,1,1) = rs(ii,2,1)
            rs(i,1,2) = rs(ii,2,2)
            rs(i,1,3) = rs(ii,2,3)
            rs(i,1,4) = rs(ii,2,4)
         end if
   30 continue
c
      return
      end
c
c     ******************************************************************
c
      subroutine  res_cons_to_sym
c
c     conservative residuals dw changed to nonconservative residuals
c
      use dims
c
c     ******************************************************************
c
      use flo_var
      use solv_var
c
c     ******************************************************************
c
      use flo_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      integer  :: i,j
c
c     ******************************************************************
c
      real     :: gm1
      real     :: r1,r2,r3,r4
      real     :: ua,va,qq,cc,c
c
c     ******************************************************************
c
      gm1   = gamma-1.

      do j=2,jl
      do i=2,il

        r1     = dw(i,j,1)
        r2     = dw(i,j,2)
        r3     = dw(i,j,3)
        r4     = dw(i,j,4)

        ua     = w(i,j,2)/w(i,j,1)
        va     = w(i,j,3)/w(i,j,1)
        qq     =.5*(ua*ua  +va*va)
        cc     = gamma*p(i,j)/w(i,j,1)
        c      = sqrt(cc)

        dw(i,j,1) = gm1*(r4  +qq*r1  -ua*r2  -va*r3)/cc
        dw(i,j,2) = (r2  -ua*r1)/c
        dw(i,j,3) = (r3  -va*r1)/c
        dw(i,j,4) = dw(i,j,1)  -r1

      end do
      end do

      return

      end
c
c     ******************************************************************
c
      subroutine  res_sym_to_cons
c
c     nonconservative residuals dw changed to conservative residuals
c
      use dims
c
c     ******************************************************************
c
      use solv_var
      use flo_var
      use psm_var
c
c     ******************************************************************
c
      use flo_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      integer  :: i,j
c
c     ******************************************************************
c
      real     :: gm1
      real     :: r1,r2,r3,r4
      real     :: ua,va,ha,qq,cc,c
c
c     ******************************************************************
c
      gm1   = gamma-1.

      do j=2,jl
      do i=2,il

        r1     = rs(i,j,1)
        r2     = rs(i,j,2)
        r3     = rs(i,j,3)
        r4     = rs(i,j,4)

        ua     = w(i,j,2)/w(i,j,1)
        va     = w(i,j,3)/w(i,j,1)
        ha     = (w(i,j,4)  +p(i,j))/w(i,j,1)
        qq     =.5*(ua*ua  +va*va)
        cc     = gamma*p(i,j)/w(i,j,1)
        c      = sqrt(cc)

        dw(i,j,1) = r1  -r4
        dw(i,j,2) = ua*(r1  -r4)  +c*r2
        dw(i,j,3) = va*(r1  -r4)  +c*r3
        dw(i,j,4) = ha*r1  +c*(ua*r2  +va*r3) -qq*r4

      end do
      end do

      return

      end
c
c     ******************************************************************
c
      subroutine solve_rs(f,b)
c
c      solve 4x4 system b*x = f
c      replaces x in place of f
c
      implicit none
c
c     ******************************************************************
c
      integer              :: i,j
c
c     ******************************************************************
c
      real, dimension(4,4)      :: b
      real, dimension(4)        :: f
c
c     ******************************************************************
c
      real     :: d1,d2,d3,d4
c
c     ******************************************************************
c
      b(1,1)     = 1./b(1,1)
      b(1,2)     = b(1,2)*b(1,1)
      b(1,3)     = b(1,3)*b(1,1)
      b(1,4)     = b(1,4)*b(1,1)
c
      b(2,2)     = 1./(b(2,2) -b(2,1)*b(1,2))
      b(2,3)     = (b(2,3)-b(2,1)*b(1,3))*b(2,2)
      b(2,4)     = (b(2,4)-b(2,1)*b(1,4))*b(2,2)
      b(3,2)     = b(3,2)-b(3,1)*b(1,2)
      b(3,3)     = 1./(b(3,3) -b(3,1)*b(1,3) - b(3,2)*b(2,3))
      b(3,4)     = (b(3,4) -b(3,1)*b(1,4) -b(3,2)*b(2,4))*b(3,3)
      b(4,2)     = b(4,2) -b(4,1)*b(1,2)
      b(4,3)     = b(4,3) -b(4,1)*b(1,3) -b(4,2)*b(2,3)
      b(4,4)     = 1./(b(4,4)-b(4,1)*b(1,4)-b(4,2)*b(2,4)-b(4,3)*b(3,4))
c
c            fwd substitution
c
      d1         = f(1)*b(1,1)
      d2         = (f(2) -b(2,1)*d1)*b(2,2)
      d3         = (f(3) -b(3,1)*d1-b(3,2)*d2)*b(3,3)
      d4         = (f(4) -b(4,1)*d1-b(4,2)*d2-b(4,3)*d3)*b(4,4)
c
c            bwd substitution
c
      f(4)       = d4
      f(3)       = d3 -b(3,4) *f(4)
      f(2)       = d2 -b(2,3)*f(3)-b(2,4)*f(4)
      f(1)       = d1-b(1,2)*f(2)-b(1,3)*f(3)- b(1,4)*f(4)
c
      return
c
      end
