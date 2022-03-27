      subroutine  turb2(il,jl,ie,je,ib,jb,itl,itu,
     & w,p,rev,x,vol,
     & gamma,rm,re,xtran,ncyc)
c
c     **********************************************************************
c     *                                                                    *
c     *   baldwin-lomax turbulence model:  modtur = 2                      *
c     *                                                                    *
c     *   calculates turbulent viscosity at the cell faces and then        *
c     *   averages to obtain cell center values                            *
c     *   fully vectorized routine                                         *
c     *                                                                    *
c     **********************************************************************
c
c     use dims
c
c     ******************************************************************
c
c     use flo_var
c
c     ******************************************************************
c
c     use solv_var
c
c     ******************************************************************
c
c     use mesh_var
c
c     ******************************************************************
c
c     use psm_var
c
c     ******************************************************************
c
c     use flo_param
c
c     ******************************************************************
c
c     use solv_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
c     input variables
c
c     ******************************************************************
c     from dims
      integer, intent(inout) :: il,jl,ie,je,ib,jb,itl,itu

c     from flo_var
      real(8), intent(in), dimension(0:ib,0:jb,4) :: w
      real(8), intent(in), dimension(0:ib,0:jb) :: p
      real(8), intent(inout), dimension(0:ib,0:jb) :: rev

c     from mesh_var
      real(8), intent(in), dimension(1:il,1:jl,2) :: x
      real(8), intent(in), dimension(0:ib,0:jb) :: vol

c     from flo_param
      real(8), intent(in)     :: gamma,rm,re,xtran

c     from solv_param
      integer, intent(in)     :: ncyc
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer :: i,j,i2,j2,jlm,jstop,itlp,iwrit,modbl,restarr,nturbw
      integer :: itr1,itr2,itrp1,jmyv,jwrit,itup,kwrite1,iloc
      integer :: jmaxv,jmaxut,jminut,jmaxyv,jcros,jcrosm,icross,ivect
      real(8) :: aplusi,cwk1,ckleb,ccp,rey,rein,sgam,sgrm,sgrmi
      real(8) :: xxa,yxa,uy,vy,uavg,vavg,vor1,vor2,vor3,vort,utotal
      real(8) :: xyw,xye,yyw,yye,volawi,volaei,uxe,vxe,uxw,vxw
      real(8) :: tur1,tur2,tur3,volai,amub,avor1,avora,avorb,avorc,avor2
      real(8) :: xc2,yc2,xyc,yyc,ylen1,ylenm1,dyvm,dyvp,y1,damp
      real(8) :: udiff,udiff1,coeff,fwake1,coeff2,fwake2,fwake,fkleb0
      real(8) :: fkleb1,tscali,amuti1,cmuti,cmuto
      real(8) :: ystop,amudif,amutc,scale,delta1,dstar1
      real(8) :: ue1,tauw1,aylen,yplus,cf1,atauw,atauwr,ustar,ustari
      real(8) :: yval,ayval,uval,uplus
      integer, dimension(1:il) :: jedge
      real(8), dimension(1:il)      :: uedge,tauw,yscal,scalf,avor
      real(8), dimension(1:il)      :: avorm,ravg,amuto,amuti,yvor
      real(8), dimension(1:il)      :: yvorm,utotm,ylenm,fkleb
      real(8), dimension(1:il)      :: utmin,utmax,utot1,fcros
      real(8), dimension(0:ib,0:jb) :: amu,u,v,vor,rinv,utot,vola,t,ylen
      real(8), dimension(1:ie,1:je) :: amut

      integer, external :: ismin,ismax
      real(8), external :: cvmgp

c     dimension       uedge(il),tauw(il),yscal(il),scalf(il),avor(il),
c    .                avorm(il),ravg(il),amuto(il),amuti(il),yvor(il),
c    .                yvorm(il),utotm(il),ylenm(il),fkleb(il),jedge(il)
c     dimension       utmin(il),utmax(il),utot1(il),fcros(il)
c     dimension       amu(ib,jb),u(ib,jb),v(ib,jb),vor(ib,jb),
c    .                rinv(ib,jb),utot(ib,jb),vola(ib,jb),
c    .                t(ib,jb),ylen(ib,jb)
c     dimension       amut(ie,je)
cs
      i2        = ie
      j2        = je
co    il        = i2- 1
co    jl        = j2- 1
      jlm       = jl- 1
c
      jstop     = 3* (j2- 2)/5
co    if (cmesh .lt. 0.0) jstop = jl- 1
      itlp      = itl+ 1
      iwrit     = 6
c
      aplusi    = 1./26.
      cwk1      = 1.0
      ckleb     = 0.3
      ccp       = 1.6
      modbl     = 3
      restarr   = 0
c
      rey       = re
      rein      = 1.0/rey
      sgam      = sqrt(gamma)
      sgrm      = sgam*rm
      sgrmi     = 1.0/sgrm
c
c     **********************************************************************
c     *                                                                    *
c     *    next 2 lines of code is activated to "test" turbulence model    *
c     *      effect on convergence                                         *
c     *                                                                    *
c     *    turbulent viscosity is frozen after ncyct cycles                *
c     *                                                                    *
c     **********************************************************************
c
c     ncyct     = 10
c     if (ncyc .gt. ncyct) return
c
      do 5 j=1,j2
      do 5 i=1,i2
        rinv(i,j) = 1.0/w(i,j,1)
        t(i,j)    = p(i,j)* rinv(i,j)
        u(i,j)    = w(i,j,2)* rinv(i,j)
        v(i,j)    = w(i,j,3)* rinv(i,j)
        amu(i,j)  = t(i,j)
        amut(i,j) = 0.0
    5 continue
c
c     **********************************************************************
c     *                                                                    *
c     *   determination of eddy viscosity                                  *
c     *                                                                    *
c     *   turbulence model:                                                *
c     *     wall boundary layer --- baldwin-lomax model                    *
c     *     wake region         --- baldwin-lomax model (cwake= 1.0)       *
c     *                                                                    *
c     *   calculate vorticity and total velocity                           *
c     *                                                                    *
c     **********************************************************************
c
      do 10 j=1,jl
      do 10 i=1,i2
        vola(i,j) = 0.5* (vol(i,j)+ vol(i,j+1))
   10 continue
c
      do 15 i=2,il
        xxa       = x(i,1,1)-x(i-1,1,1)
        yxa       = x(i,1,2)-x(i-1,1,2)
        uy        = u(i,2)- u(i,1)
        vy        = v(i,2)- v(i,1)
        uavg      = .5* (u(i,1)+u(i,2))
        vavg      = .5* (v(i,1)+v(i,2))
        vor1      = (xxa*uy + yxa*vy)/vola(i,1)
        vor2      = 0.0
        vor3      = 0.0
        vort      = vor1-vor2-vor3
        vor(i,1)  = abs(vort)
        utotal    = uavg*uavg + vavg*vavg
        utot(i,1) = sqrt(utotal)
   15 continue
c
      do 20 j=2,jlm
      do 20 i=2,il
        xxa       = x(i,j,1)-x(i-1,j,1)
        yxa       = x(i,j,2)-x(i-1,j,2)
        uy        = u(i,j+1)- u(i,j)
        vy        = v(i,j+1)- v(i,j)
        uavg      = 0.5* (u(i,j)+ u(i,j+1))
        vavg      = 0.5* (v(i,j)+ v(i,j+1))
c
c     thin-layer navier-stokes contribution to vorticity
c
        vor1      = (xxa*uy + yxa*vy)/vola(i,j)
c
c     additional contributions to vorticity
c
        xyw       = 0.5* (x(i-1,j+1,1)- x(i-1,j-1,1))
        xye       = 0.5* (x(i,j+1,1)- x(i,j-1,1))
        yyw       = 0.5* (x(i-1,j+1,2)- x(i-1,j-1,2))
        yye       = 0.5* (x(i,j+1,2)- x(i,j-1,2))
        volawi    = 2.0/(vola(i,j)+ vola(i-1,j))
        volaei    = 2.0/(vola(i,j)+ vola(i+1,j))
        uxe       = 0.5* (u(i+1,j)+u(i+1,j+1)) - uavg
        vxe       = 0.5* (v(i+1,j)+v(i+1,j+1)) - vavg
        uxw       = uavg - 0.5* (u(i-1,j)+u(i-1,j+1))
        vxw       = vavg - 0.5* (v(i-1,j)+v(i-1,j+1))
        vor2      = 0.5* (xye* volaei* uxe+ xyw* volawi* uxw)
        vor3      = 0.5* (yye* volaei* vxe+ yyw* volawi* vxw)
        vort      = vor1- vor2- vor3
        vor(i,j)  = abs(vort)
        utotal    = uavg* uavg+ vavg* vavg
        utot(i,j) = sqrt(utotal)
   20 continue
c
c     determine transition index
c
      itr1      = 0
      itr2      = 0
      j         = 1
      do i=1,il
        if (x(i,j,1) .le. xtran) then
          itr1      = i - 1
          go to 21
        end if
      end do
   21 continue
c
      itrp1     = itr1 + 1
      do i=itrp1,il
        if (x(i,j,1) .ge. xtran) then
          itr2      = i
          go to 22
        end if
      end do
   22 continue
c
      do 30 i=2,il
c
        do 25 j=1,jlm
          avor(j)   = vor(i,j)
          utot1(j)  = utot(i,j)
   25   continue
c
c
c     effect of using jlm or jstop needs to be checked;
c     it does not seem to make a difference
c
c       jmaxv     = ismax(jstop,avor,1)
        jmaxv     = ismax(jlm,avor,1)
        if (jmaxv .eq. 0) jmaxv = 1
cc
        jminut    = ismin(jlm,utot1,1)
cc
        jmaxut    = ismax(jstop,utot1,1)
cc
        avorm(i)  = avor(jmaxv)
        utmin(i)  = utot1(jminut)
        utmax(i)  = max(utot1(jmaxut),1.e-3)
        utotm(i)  = utmax(i)
        yscal(i)  = 1000000.
   30 continue
c
      if (modbl .eq. 1) then
        tur1    = 1.0
        tur2    = 0.
        tur3    = 0.
      else if (modbl .eq. 2) then
        tur1    = 0.
        tur2    = 1.0
        tur3    = 0.
      else
        tur1    = 0.
        tur2    = 0.
        tur3    = 1.0
      end if
c
      do 35 i=itlp,itu
        xxa       = x(i,1,1)- x(i-1,1,1)
        yxa       = x(i,1,2)- x(i-1,1,2)
        volai     = 1.0/vola(i,1)
        uy        = 2.0* u(i,2)
        amub      = .5* (amu(i,1)+ amu(i,2))
        tauw(i)   = amub* (xxa* uy)* volai
        avor1     = vor(i,1)
        avora     = avor1
        avorb     = 0.5* (avor1+ avorm(i))
        avorc     = avorm(i)
        avor2     = tur1*avora + tur2*avorb + tur3*avorc
        yscal(i)  = sqrt(rey* sgrmi* amub* avor2* w(i,2,1))/(26.*amub)
   35 continue
c
c     **********************************************************************
c     *                                                                    *
c     *   compute normal distance ylen(i,j) and function  yvor             *
c     *   (yvor = y* vorticity)                                            *
c     *                                                                    *
c     **********************************************************************
c
co    if (ncyc .eq. ncyci1) then
        do 39 i=2,il
          ylen(i,1) = 0.0
          do 37 j=2,jl
            xc2       = .50* (x(i,j,1)+ x(i-1,j,1)
     1                       -x(i,j-1,1)- x(i-1,j-1,1))
            yc2       = .50* (x(i,j,2)+ x(i-1,j,2)
     1                       -x(i,j-1,2)- x(i-1,j-1,2))
            xyc       = xc2
            yyc       = yc2
            scalf(j)  = sqrt(xyc*xyc + yyc*yyc)
            ylen(i,j) = ylen(i,j-1)+ scalf(j)
   37     continue
   39   continue
co    end if
c
      do 50 i=2,il
        ylen1     = 0.5* ylen(i,2)
        do 40 j=1,jstop
          y1        = yscal(i)* ylen(i,j)
          damp      = 1.0- exp(-y1)
          yvor(j)   = ylen(i,j)* vor(i,j)* damp
   40   continue
c
        jmaxyv    = ismax(jstop,yvor,1)
        jmaxyv    = max(jmaxyv,2)
        jedge(i)  = jmaxyv
c
c     next line of code replaced because it caused convergence
c     stall when m = 0.001 - 12/10/05 (check this further !!!!)
c
        yvorm(i)  = max(yvor(jmaxyv),1.e-6)
        ylenm(i)  = max(ylen(i,jmaxyv),ylen1)
c
        if (jedge(i) .lt. jstop) then
          ylenm1  = ylenm(i)
c
          if (ncyc.ge.10 .or. restarr.eq.1) then
            jmyv    = jedge(i)
            dyvm    = yvor(jmyv)-yvor(jmyv-1)
            dyvp    = yvor(jmyv)-yvor(jmyv+1)
c
            if (yvor(jmyv-1) .lt. yvor(jmyv+1)) then
              ylenm(i) = ylen(i,jmyv)+ .5*(ylen(i,jmyv+1)- ylen(i,jmyv))
     .                   *(1.- dyvp/dyvm)
            else
              ylenm(i) = ylen(i,jmyv)- .5*(ylen(i,jmyv)- ylen(i,jmyv-1))
     .                   *(1.- dyvm/dyvp)
            end if
c
          else
            ylenm(i)  = ylenm1
          end if
c
        end if
   50 continue
c
c     **********************************************************************
c     *                                                                    *
c     *   compute outer eddy viscosity                                     *
c     *                                                                    *
c     *   outer do loop                                                    *
c     *                                                                    *
c     **********************************************************************
c
      do 200 i=2,il
c
        udiff     = abs(utmax(i)- utmin(i))
        udiff1    = cwk1* udiff
        do 60 j=2,jstop
          ravg(j)   = 0.5* (w(i,j,1)+ w(i,j+1,1))
          coeff     = 0.0168* ccp
          fwake1    = coeff* yvorm(i)* ylenm(i)
          coeff2    = coeff* cwk1* cwk1
          fwake2    = coeff2* ylenm(i)* udiff* udiff/yvorm(i)
          fwake     = min(fwake1,fwake2)
          fkleb0    = ckleb* ylen(i,j)/ylenm(i)
          fkleb1    = min(fkleb0,1.e5)
          fkleb(j)  = 1.0/(1.0+ 5.5* fkleb1**6)
          amuto(j)  = rey* sgrmi* ravg(j)* fwake* fkleb(j)
          amuto(j)  = abs(amuto(j))
   60   continue
        amuto(1)  = amuto(2)
c
c     **********************************************************************
c     *                                                                    *
c     *   compute inner eddy viscosity                                     *
c     *                                                                    *
c     **********************************************************************
c
        do 70 j=2,jstop
          y1        = yscal(i)* ylen(i,j)
          damp      = 1.0- exp(-y1)
          tscali    = 0.4* ylen(i,j)* damp
          amuti1    = tscali* tscali* vor(i,j)
          amuti(j)   = rey* sgrmi* ravg(j)* amuti1
          amuti(j)   = abs(amuti(j))
   70   continue
        amuti(1)  = 0.0
        if (i.le.itl .or. i.gt.itu) amuti(1) = amuti(2)
c
c        load viscosity coeffs. into array, use inner value until
c        match point is reached
c        scalar coding
c
        ivect     = 1
        if (ivect .eq. 0) then
          icross    = 0
          amut(i,1) = amuti(1)
          do 75 j=2,jstop
            if (amuti(j).le.amuto(j) .and. icross.eq.0) then
              amut(i,j) = amuti(j)
            else
              icross    = 1
              amut(i,j) = amuto(j)
            end if
   75     continue
        else
c
          amut(i,1) = amuti(1)
          ystop     = jstop
          do 80 j=1,jstop
            amudif    = amuti(j)- amuto(j)
            fcros(j)  = cvmgp(real(j,8),real(1000.,8),amudif)
   80     continue
          jcros    = ismin(jstop,fcros,1)
          if (jcros .eq. 1) jcros = 2
          jcrosm   = jcros- 1
c
          do 90 j=1,jcrosm
            amut(i,j) = amuti(j)
   90     continue
c
          do 100 j=jcros,jstop
            amut(i,j) = amuto(j)
  100     continue
        end if
c
c     **********************************************************************
c     *                                                                    *
c     *   compute turbulent viscosity at cell center                       *
c     *                                                                    *
c     **********************************************************************
c
        do 110 j=2,jstop
          amutc     = 0.5* (amut(i,j)+ amut(i,j-1))
          amu(i,j)  = amutc
  110   continue
c
        do 120 j=2,jstop
          amut(i,j) = amu(i,j)
  120   continue
c
        if (i.gt.itr1 .and. i.le.itr2) then
          do 130 j=2,jstop
            amut(i,j) = 0.
  130     continue
        end if
c
c     **********************************************************************
c     *                                                                    *
c     *   debugging check (activate with jwrit = 1)                        *
c     *                                                                    *
c     **********************************************************************
c
        jwrit     = 0
        if (jwrit .eq. 1) then
          write(6,8000) i
 8000     format(5x,'i = ',i5)
          write(6,8050)
 8050     format(5x,'j,y1,fkleb,ravg,muti,muto,mut')
 8100     format(i5,6e15.5)
          do 195 j=2,jstop
            y1        = yscal(i)*ylen(i,j)
            cmuti     = amuti(j)
            cmuto     = amuto(j)
            write(6,8100) j,y1,fkleb(j),ravg(j),cmuti,cmuto,amut(i,j)
  195     continue
        end if
c
c       finish of outer loop
  200 continue
c
c       copy amut to rev
c
      scale     = 1.
      do j=1,je
      do i=1,ie
         rev(i,j)  = scale*amut(i,j)
      end do
      end do
      return
c
c     **********************************************************************
c     *                                                                    *
c     *   printout for boundary layer details                              *
c     *   nturbw  -  input for frequency of printing b.l. data             *
c     *              (currently set at large value)                        *
c     *                                                                    *
c     **********************************************************************
c
      nturbw = 5000
c
      if (mod(ncyc,nturbw) .eq. 0) then
        do 250 i=2,il
          if(i .eq. 2) write(iwrit,840)
          write (iwrit,850) i,jedge(i),ylenm(i),utotm(i),yvorm(i),
     .                      yscal(i)
  250   continue
c
        do 300 i=itlp,itu
          if(i.eq.itlp) write(6,860)
          je      = jedge(i)
          delta1  = ylenm(i)
          dstar1  = 1.0
          uedge(i)= w(i,je,2)/w(i,je,1)
          ue1     = uedge(i)* sgrmi
          tauw1   = abs(tauw(i))
          aylen   = abs(ylen(i,2))
          yplus   = yscal(i)* aylen* 26.
          yplus   = 0.5* yplus
          cf1     = 2.0* rein* sgrmi* tauw1
          write(6,870) i,je,ue1,delta1,dstar1,tauw1,cf1,yplus
  300   continue
c
        iloc    = itu
        atauw   = abs(tauw(iloc))
        atauwr  = atauw* rinv(iloc,2)
        ustar   = sqrt(rein* rm* sgam* atauwr)
        if(ustar.eq.0.0) ustar = 1.0
        ustari  = 1.0/ustar
c
        do 310 j=2,j2
          yval    = 0.5* (ylen(iloc,j-1)+ ylen(iloc,j))
          ayval   = abs(yval)
          yplus   = yscal(iloc)* ayval* 26.
          uval    = w(iloc,j,2)* rinv(iloc,j)
          uplus   = uval* ustari
          write(6,880) j,yplus,uplus,ylen(iloc,j),ylen(iloc,j-1),yval,
     .                 yscal(iloc),yplus
  310   continue
      end if
c
c
c     **********************************************************************
c     *                                                                    *
c     *   additional checks of routine (activate with kwrite = 1)          *
c     *                                                                    *
c     **********************************************************************
c
      kwrite1   = 0
      if (kwrite1.eq.1) then
        do i=itl,itl-4,-4
        do j=1,10
          write(970,971) i,j,amut(i,j)
  971     format (5x,'i,j,amut = ',2i5,e15.6)
        end do
        end do
c
        itup      = itu + 1
        do i=itup,itup+4,4
        do j=1,10
          write(970,971) i,j,amut(i,j)
        end do
        end do
c
c---->    finish of kwrite if
        end if
c
      return
c
  800 format(1h ,'ylen,ravg,amuto1,amuto2,amuto,udiff,fkleb')
  810 format(1h ,5x,i5,7e15.6)
  820 format(1h ,'ylen,utot,vor,yvor,amuti,amuto,fkleb')
  830 format(1h ,5x,i5,7e15.6)
  840 format(1h ,'ylenm,utotm,yvorm,yscal')
  850 format(1h ,5x,2i5,4e15.6)
  860 format(1h ,1x,4h  i ,4h je ,1x,
     .  10h    ue    ,10h   delta  ,10h   dstar  ,
     .  10h   tauw   ,12h     cf     ,10h   yplus  )
  870 format(1x,2i4,3f10.6,2e15.6,f10.6)
  880 format(1h ,5x,i5,7f12.4)

      end subroutine

      integer function ismin(n,sx,incx)
c
c     find the largest value of the
c     elements of sx..if n<=0, return 0
c
      integer, intent(in) :: n,incx
      real(8), dimension(n) :: sx
      real(8) :: sxmin
c
      ismin = 0
      if(n .le. 0) return
      sxmin = 1000000000
      do 100 i=1,n,incx
      if (sx(i) .lt. sxmin) then
      sxmin = sx(i)
      ismin = i
      end if
  100 continue
c
      return
      end function
c
c     ******************************************************************
c
      integer function ismax(n,sx,incx)
c
c     find the largest value of the
c     elements of sx..if n<=0, return 0
c
      integer, intent(in) :: n,incx
      real(8), dimension(n) :: sx
      real(8) :: sxmax
c
      ismax = 0
      if(n .le. 0) return
      sxmax = 0
      do 100 i=1,n,incx
      if (sx(i) .gt. sxmax) then
      sxmax = sx(i)
      ismax = i
      end if
  100 continue
c
      return
      end function
c
c     ******************************************************************
c
      real(8) function cvmgp(x1,x2,x3)
c
c     tests for positive or zero.  x1 is returned if x3 >= 0.
c     x2 is returned if x3 < 0.
c
      real(8), intent(in) :: x1,x2,x3

      if(x3.lt.0.0) then
        cvmgp=x2
      else
        cvmgp=x1
      endif
c
      return
      end function
c
c     ******************************************************************
c
      
