      subroutine  rhs_face(vnx, vny, sar, dtv2, rmuem,
     #                     rl, pl, ul, vl,
     #                     rr, pr, ur, vr,
     #                     a1, a2, a3, a4,
     #                     rhs_f1, rhs_f2, rhs_f3, rhs_f4,
     #                     xx)
c
c******  computation of right-hand-side face-coefficients for sgs********
c
c
c     ******************************************************************
c
      use flo_param
c
c     ******************************************************************
c
      use solv_param
c
c     ******************************************************************
c
      implicit none
      real, dimension(4,4)      :: xx

      real     :: gm1,dlim
      real     :: vnx,vny,sar,dtv2,rmuem
      real     :: rl,pl,ul,vl,rr,pr,ur,vr,a1,a2,a3,a4
      real     :: rhoa,us,vs,ps,qs,ccs,cs
      real     :: e1,e2,e3,a,q1,q2,q3,r1,r2,aqsa
      real     :: amu,amupr,vq1,vr1
      real     :: a11,a12,a13,a14,a22,a23,a33,a41,a44
      real     :: rhs_f1,rhs_f2,rhs_f3,rhs_f4
c
      gm1     = gamma - 1.
      dlim    = .10
      if (kvis.eq.0) dlim = .25
      if (ncyc.lt.5) dlim = .5
c
c---- arithmetic averaging
c
      rhoa  = 0.5*(rl  +rr)
      us    = 0.5*(ul  +ur)
      vs    = 0.5*(vl  +vr)
      ps    = 0.5*(pl  +pr)
      qs    = vnx*us  +vny*vs
      ccs   = gamma*ps/rhoa
      cs    = sqrt(ccs)
c
c---- viscous coefficients
c
      amu   = dtv2*rmuem
      amupr = amu/(prn  +prt)
c
c---- eigenvalues
c
      e1    = abs(qs)
      e2    = abs(qs  +cs)
      e3    = abs(qs  -cs)
      a     = dlim*cs
      if (e1.lt.a) e1 = .5*(a  +e1*e1/a)
      if (e2.lt.a) e2 = .5*(a  +e2*e2/a)
      if (e3.lt.a) e3 = .5*(a  +e3*e3/a)
      a     = sar*dtv2
      aqsa  = a*diag*abs(qs)
c
c---- form the right hand side with the negative jacobian
c
      q1    = a*(qs  -e1)
      q2    = a*(qs  +cs  -e2)
      q3    = a*(qs  -cs  -e3)
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
c---- form the preconditioning matrix wtth the positive jacobian
c
      q1    = a*(qs  +e1)
      q2    = a*(qs  +cs  +e2)
      q3    = a*(qs  -cs  +e3)
      r1    = .5*(q2  +q3)  -q1
      r2    = .5*(q2  -q3)

      vq1   = q1  +amu
      vr1   = r1  +amu/3.

      xx(1,1)   =  r1  +q1  +gm1*amupr  +aqsa  + xx(1,1)
      xx(1,2)   =  vnx*r2                      + xx(1,2)
      xx(1,3)   =  vny*r2                      + xx(1,3)
      xx(1,4)   =  amupr                       + xx(1,4)

      xx(2,1)   =  vnx*r2                      + xx(2,1)
      xx(2,2)   =  vnx*vnx*vr1  +vq1  +aqsa    + xx(2,2)
      xx(2,3)   =  vnx*vny*vr1                 + xx(2,3)

      xx(3,1)   =  vny*r2                      + xx(3,1)
      xx(3,2)   =  vnx*vny*vr1                 + xx(3,2)
      xx(3,3)   =  vny*vny*vr1  +vq1  +aqsa    + xx(3,3)

      xx(4,1)   =  gm1*amupr                   + xx(4,1)
      xx(4,4)   =  q1  +amupr  +aqsa           + xx(4,4)

      return

      end
