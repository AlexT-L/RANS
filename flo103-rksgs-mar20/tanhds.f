c
c***********************************************************************
c
      subroutine tanhds (ns, sp1, sp2, s, status)
c
c***************
c
c    Program by:  jc vassberg    copyright (c) 2000
c    Written on:  19 apr 2000
c    Updated on:  19 apr 2000
c
c***************
c
c    TANHDS computes a normalized arclength distribution based on the
c    hyperbolic tangent function.  Spacings can be specified at either
c    or both ends.  See Joe Thompson's grid generation book.
c
c    Nomenclature
c    ------------
c        ns:  number of points in distribution     (input)
c       sp1:  normalized initial spacing           (input)
c       sp2:  normalized final spacing             (input)
c    status:  error status, ge.0 = no errors       (output)
c         s:  normalized arclength distribution    (output)
c
c    Notes
c    -----
c      * A two-sided distribution results for sp1 > 0.0 and sp2 > 0.0
c        with initial spacing (sp1) and final spacing (sp2)
c
c      * A one-sided distribution results for sp1 > 0.0 and sp2 = 0.0
c        with initial spacing (sp1)
c
c      * A one-sided distribution results for sp1 = 0.0 and sp2 > 0.0
c        with final spacing (sp2)
c
c      * Equal spacing results if sp1 = sp2 = 0.0
c
c      * Iterates until (delta sp < tol*sp) or (iter = itmx)
c
c      * TANHDS calls:  functions ZSINHV & ATANHV
c
c***************
c
      implicit none
c
      integer
     .   ns, status,
     .   i, j, iter, itmx, n, nm1, nsmid
c
      real*8
     .   sp1, sp2, s(ns), si, sj,
     .   a, atanhv, aterm, b, bterm, c, cterm, 
     .   del, del1, del2, delcmp, delta, delta1, 
     .   denom, divfac, zsinhv, f, fp, fpa, fpb, 
     .   fpc, pct, tol,
     .   pct1, pct2, scmp, scmp1, scmp2, 
     .   sp, u, u1, u2, xden, xintnh, 
     .   xnm1, xnum, xnum1, xnum2, xnumer
c
c***************
c
c    Initialize
      status    = 0
      itmx      = 100
      tol       = 0.005d0
c     tol       = 0.0005d0
c     tol       = 0.0001d0
      s(1)      = 0.d0
      s(ns)     = 1.d0
      nm1       = ns - 1
      xnm1      = dble(nm1)
      nsmid     = (ns + 1) / 2
c
c    Done
      if (ns .eq. 2) then
       goto 999
      endif
c
c    Bad input
      if (ns .lt. 2) then
         status    = -1
         goto 999
      endif
c
c    Bad input
      if ((sp1 .lt. 0.d0) .or. (sp2 .lt. 0.d0)) then
         status    = -1
         goto 999
      endif
c
c    Hyperbolic tangent distribution (two-sided)
      if ((sp1 .gt. 0.d0) .and. (sp2 .gt. 0.d0)) then
c
c       Bad input
         if ((sp1+sp2) .gt. 1.d0) then
            status    = -1
            goto 999
         endif
c
c       Initial guess for delta
         b         = 1.d0 / (xnm1*sqrt(sp1*sp2))
         delta     = zsinhv(b)
         delta1    = delta
c
c       Bad input
         if (b .eq. 1.d0) then
            status    = -1
            goto 999
         endif
c
         iter      = 1
  200    continue
            if (b .gt. 1.d0) then
               c         = (tanh(delta*((1.d0/xnm1)-0.5d0))
     .                     /tanh(0.5d0*delta))
               a         = (1.d0+c-sp1-(c*sp1)) / (sp1-c*sp1)
               denom     = sp2*(a-1.d0) - a
               xnumer    = (tanh(0.5d0*delta))*(a-sp2-a*sp2)
               xintnh    = atanhv(-xnumer/denom)
               divfac    = ((xnm1-1.d0)/xnm1) - 0.5d0
               delcmp    = xintnh / divfac
               xnum1     = tanh(delcmp*(1.d0/xnm1-0.5d0))
               xnum2     = tanh(delcmp*((xnm1-1.d0)/xnm1-0.5d0))
               xden      = tanh(0.5d0*delcmp)
            else
     .      if (b .lt. 1.d0) then
               c         = (tan(delta*((1.d0/xnm1)-0.5d0))
     .                     /tan(0.5d0*delta))
               a         = (1.d0+c-sp1-(c*sp1)) / (sp1-c*sp1)
               denom     = sp2*(a-1.d0) - a
               xnumer    = (tan(0.5d0*delta))*(a-sp2-a*sp2)
               xintnh    = atan(-xnumer/denom)
               divfac    = ((xnm1-1.d0)/xnm1) - 0.5d0
               delcmp    = xintnh / divfac
               xnum1     = tan(delcmp*(1.d0/xnm1-0.5d0))
               xnum2     = tan(delcmp*((xnm1-1.d0)/xnm1-0.5d0))
               xden      = tan(0.5d0*delcmp)
            endif
c
c          Check for convergence
            u1        = 0.5d0*(1.d0+xnum1/xden)
            scmp1     = u1 / (a+(1.d0-a)*u1)
            u2        = 0.5d0*(1.d0+xnum2/xden)
            scmp2     = 1.d0 - (u2/(a+(1.d0-a)*u2))
            del1      = abs(sp1-scmp1)
            del2      = abs(sp2-scmp2)
            pct1      = del1 / sp1
            pct2      = del2 / sp2
c
c          Not converged
            if ((pct1 .gt. tol) .or. (pct2 .gt. tol)) then
               delta     = delcmp
               iter      = iter + 1
               if (iter .lt. itmx) then
                  goto 200
               else
                  delta     = delta1
                  if (b .gt. 1.d0) then
                     c         = (tanh(delta*((1.d0/xnm1)-0.5d0))
     .                           /tanh(0.5d0*delta))
                  else
     .            if (b .lt. 1.d0) then
                     c         = (tan(delta*((1.d0/xnm1)-0.5d0))
     .                           /tan(0.5d0*delta))
                  endif
                  a         = (1.d0+c-sp1-(c*sp1))/(sp1-c*sp1)
                  delcmp    = delta
               endif
            endif
c        end-loop-200
c
c       Converged
         status    = iter
         delta     = delcmp
         do 210  i = 2,nm1
            if (b .gt. 1.d0) then
               xnum      = tanh(delta*(dble(i-1)/xnm1-0.5d0))
               xden      = tanh(0.5d0*delta)
            else
     .      if (b .lt. 1.d0) then
               xnum      = tan(delta*(dble(i-1)/xnm1-0.5d0))
               xden      = tan(0.5d0*delta)
            endif
            u         = 0.5d0*(1.d0+xnum/xden)
            s(i)      = u / (a+(1.d0-a)*u)
  210    continue
c
c    Hyperbolic tangent distribution (one-sided)
      else
     .if ((sp1 .gt. 0.d0) .or. (sp2 .gt. 0.d0)) then
c
c       Bad input
         if ((sp1 .gt. 0.5d0) .or. (sp2 .gt. 0.5d0)) then
            status    = -1
            goto 999
         endif
c
         if (sp1 .gt. 0.d0) sp = sp1
         if (sp2 .gt. 0.d0) sp = sp2
c
c       Initial guess for delta
         b         = 1.d0 / (xnm1*sp)
         delta     = zsinhv(b)
         delta1    = delta
         iter      = 1.d0
         aterm     = 0.5d0*((1.d0/xnm1)-1.d0)
         bterm     = 0.5d0
         cterm     = sp - 1.d0
  100    continue
            if (b .gt. 1.d0) then
               f         = tanh(aterm*delta)/tanh(bterm*delta) - cterm
               fpa       = aterm*tanh(bterm*delta)/cosh(aterm*delta)**2
               fpb       = bterm*tanh(aterm*delta)/cosh(bterm*delta)**2
               fpc       = tanh(bterm*delta)**2
            else
     .      if (b .lt. 1.d0) then
               f         = tan(aterm*delta)/tan(bterm*delta) - cterm
               fpa       = aterm*tan(bterm*delta)/cos(aterm*delta)**2
               fpb       = bterm*tan(aterm*delta)/cos(bterm*delta)**2
               fpc       = tan(bterm*delta)**2
            endif
            fp        = (fpa - fpb) / fpc
            delcmp    = delta - (f/fp)
            if (b .gt. 1.d0) then
               c         = tanh((0.5d0*delcmp)*((1.d0/xnm1)-1.d0))
               a         = tanh(0.5d0*delcmp)
            else
     .      if (b .lt. 1.d0) then
               c         = tan((0.5d0*delcmp)*((1.d0/xnm1)-1.d0))
               a         = tan(0.5d0*delcmp)
            endif
            scmp      = 1.d0 + c/a
            del       = abs(sp - scmp)
            pct       = del / sp
c
c          Not converged
            if (pct .gt. tol) then
               delta     = delcmp
               iter      = iter + 1
               if (iter .lt. itmx) then
                  goto 100
               else
                  delta     = delta1
                  if( b .gt. 1.d0 )then
                     a      = tanh(0.5d0*delta)
                  else
     .            if (b .lt. 1.d0) then
                     a      = tan(0.5d0*delta)
                  endif
                  delcmp = delta
               endif
            endif
c        end-loop-100
c
c       Converged
         status    = iter
         delta     = delcmp
         do 110  i = 2,nm1
           if (b .gt. 1.d0) then
              c         = tanh((0.5d0*delta)*((dble(i-1)/xnm1)-1.d0))
              s(i)      = 1.d0 + c/a
           else
     .     if (b .lt. 1.d0) then
              c         = tan((0.5d0*delta)*((dble(i-1)/xnm1)-1.d0))
              s(i)      = 1.d0 + c/a
           endif
  110    continue
c
c       Swap direction if required
         if (sp2 .gt. 0.d0) then
            do 120  i = 1,nsmid
               j         = ns + 1 - i
               si        = 1.d0 - s(i)
               sj        = 1.d0 - s(j)
               s(i)      = sj
               s(j)      = si
  120       continue
         endif
c
c    Equal spacing distribution (no-sided)
      else
         do 10   n = 2,nm1
            s(n)       = dble(n-1) / xnm1
   10    continue
      endif
c
  999 continue
c
      return
      end
c
c**********************************************************************
c
      real*8 function zsinhv(b)
c
c***************
c
c    program by:  jc vassberg    copyright (c) 2000
c    written on:  19 apr 2000
c    updated on:  19 apr 2000
c
c***************
c
c    ZSINHV finds x for an approximate zero of:
c       sin(x)/x - b  for  0.0 < b < 1.0
c      sinh(x)/x - b  for  1.0 < b
c
c***************
c
      implicit none
c
c***************
c
      real*8
     .   b, b1, pi, v, w, c1, c2, c3, c4, c5
c
c***************
c
c    Approximate root of sin(x)/x - b
      if ((b .gt. 0.d0) .and. (b .lt. 1.d0)) then
         if (b .le. 0.26938972d0) then
            pi        = 4*atan(1.d0)
            c1        = 11.726095d0
            c2        = 13.205501d0
            c3        =  6.794732d0
            zsinhv    = pi*((((((c1*b-c2)*b+c3)*b
     .                  -(1.d0+pi*pi/6.d0))*b+1.d0)*b-1.d0)*b+1.d0)
         else
            b1        = 1.d0 - b
            c1        = 0.075845134d0
            c2        = 0.053337753d0
            c3        = 0.048774238d0
            c4        = 0.057321429d0
            c5        = 0.150000000d0
            zsinhv    = sqrt(6.d0*b1)
     .                  *(((((c1*b1-c2)*b1+c3)*b1+c4)*b1+c5)*b1+1.d0)
         endif
c
c    Approximate root of sinh(x)/x - b
      else
     .if (b .gt. 1.d0) then
         if (b .le. 2.7829681d0) then
            b1        = b - 1.0
            c1        = 0.0010794123d0
            c2        = 0.0077424461d0
            c3        = 0.0249072950d0
            c4        = 0.0573214290d0
            c5        = 0.1500000000d0
            zsinhv    = -sqrt(6.d0*b1)
     .                  *(((((c1*b1-c2)*b1+c3)*b1-c4)*b1+c5)*b1-1.d0)
         else
            v         = log(b)
            w         = 1.d0/b - 0.028527431d0
            c1        = 8.56795911d0
            c2        = 2.62945470d0
            c3        = 1.94964430d0
            c4        = 0.24902722d0
            c5        = 0.02041793d0
            zsinhv    = v + (1.d0+1.d0/v)*log(2.d0*v)
     .                    + (((c1*w-c2)*w+c3)*w+c4)*w - c5
         endif
      else
         zsinhv    = 0.d0
      endif
c
      return
      end
c
c**********************************************************************
c
      real*8 function atanhv(x)
c
c***************
c
c    program by:  jc vassberg    copyright (c) 2000
c    written on:  19 apr 2000
c    updated on:  19 apr 2000
c
c***************
c
c    ATANHV finds the arc-hyperbolic-tangent of x.
c
c***************
c
      implicit none
c
c***************
c
      real*8
     .   x
c
c***************
c
      if (abs(x) .ge. 1.d0) then
         atanhv    = 0.d0
      else
         atanhv    = 0.5d0*log((1.d0+x)/(1.d0-x))
      endif
c
      return
      end
