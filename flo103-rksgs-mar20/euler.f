      subroutine euler
c
c     ******************************************************************
c     *                                                                *
c     *   multistage time stepping scheme with residual averaging      *
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
      integer  :: i,j,n
      integer  :: ismoop,lmax
c
c     ******************************************************************
c
      real     :: fn,dt,rqq,rh,f,qq,h,rt
c
c     ******************************************************************
c
      real, dimension(il,jl)   :: res
c
c     ******************************************************************
c
c     if ksmoop.gt.0 smooth the residuals at all stages
c     if ksmoop.lt.0 smooth the residuals at alternate stages
c
      ismoop    = iabs(ksmoop)
      if (mod(mstage,2).eq.0) ismoop = ksmoop
      rfl0      = .5*cfl/cflim
      csmoop    = 0.
c
c     ******************************************************************
c     *                                                                *
c     *   save the initial values on first entry to a coarser mesh     *
c     *   during the multigrid cycle                                   *
c     *                                                                *
c     ******************************************************************
c
      if (kode.eq.2) then

      do n=1,4
         do j=1,je
         do i=1,ie
            w1(i,j,n) = w(i,j,n)
         end do
         end do
      end do

      end if
c
c     set initial values for the multistage scheme
c
      do n=1,4
         do j=2,jl
         do i=2,il
            wn(i,j,n) = w(i,j,n)
            fw(i,j,n) = 0.
         end do
         end do
      end do
c
c     ******************************************************************
c     *                                                                *
c     *   loop over the stages of the multistage scheme                *
c     *                                                                *
c     ******************************************************************
c
      do 10 nstage=1,mstage

      fn        = .5*cstp(nstage)*cfl
      rfil      = cdis(nstage)
      if (ksmoop.lt.0) ismoop = -ismoop
c
c     calculate the artificial dissipation if rfil.ne.0.
c     use simple dissipation on coarse grids
c
      if (mode.eq.0.and.rfil.ne.0.) call dflux
      if (mode.ne.0.and.rfil.ne.0.) call dfluxc
c
c     calculate the viscous terms if kvis.gt.0 and rfil.ne.0.
c
      if (kvis.gt.0.and.rfil.ne.0.) call nsflux
c
c     calculate the euler terms
c
      call eflux
c
c     add the dissipative terms
c
      do n=1,4
         do j=2,jl
         do i=2,il
            dw(i,j,n) = dw(i,j,n)  +fw(i,j,n)
         end do
         end do
      end do

      if (kvis.gt.0) then
c
c     add the viscous terms
c
      do n=1,3
         do j=2,jl
         do i=2,il
            dw(i,j,n+1) = dw(i,j,n+1)  +vw(i,j,n)
         end do
         end do
      end do

      end if
c
c     ******************************************************************
c     *                                                                *
c     *   modification of the residuals on any coarser mesh            *
c     *   in the multigrid cycle                                       *
c     *                                                                *
c     ******************************************************************
c
      if (mode.gt.0) then

      if (kode.eq.2.and.nstage.eq.1) then

         do n=1,4
            do j=2,jl
            do i=2,il
               wr(i,j,n) = fcoll*wr(i,j,n)  -dw(i,j,n)
            end do
            end do
         end do

      end if

         do n=1,4
            do j=2,jl
            do i=2,il
               dw(i,j,n) = dw(i,j,n)  +wr(i,j,n)
            end do
            end do
         end do

      end if
c
c     return if the residuals are being calculated for transfer
c     to a coarser mesh
c
      if (kode.lt.0) return
c
c     ******************************************************************
c
c
c     multiply by the time step
c
      do j=2,jl
      do i=2,il
      if (iprec.eq.0) then
         dt        = fn*rfl(i,j)*dtl(i,j)
      else
         dt        = fn*dtl(i,j)
      end if
         res(i,j)  = dw(i,j,1)/vol(i,j)
         dw(i,j,1) = dt*dw(i,j,1)
         dw(i,j,2) = dt*dw(i,j,2)
         dw(i,j,3) = dt*dw(i,j,3)
         dw(i,j,4) = dt*dw(i,j,4)
      end do
      end do
c
c     option to smooth the residuals
c
      if (iprec.eq.0) then
        if (ismoop.gt.0) call psmoo
      else
        call psgs
      end if
c
c     update the flow variables
c
      do j=2,jl
      do i=2,il
         w(i,j,1)  = wn(i,j,1)  -dw(i,j,1)
         w(i,j,2)  = wn(i,j,2)  -dw(i,j,2)
         w(i,j,3)  = wn(i,j,3)  -dw(i,j,3)
         w(i,j,4)  = wn(i,j,4)  -dw(i,j,4)
         rqq       = .5*(w(i,j,2)**2  +w(i,j,3)**2)/w(i,j,1)
         p(i,j)    = (gamma  -1.)*dim(w(i,j,4),rqq)
      end do
      end do

      if (nstage.lt.mstage) then
c
c     update the surface pressure
c
      call bcwall
c
c     set extended values in the halo
c
      call halo

      end if

   10 continue
c
c     option for enthalpy damping on the fine mesh
c     enthalpy damping should not be used on the coarser meshes
c     of the multigrid cycle with a cell centered scheme
c
      if (kvis.eq.0.and.mode.eq.0.and.hm.gt.0.) then

      if (rm.gt.1.) hm = -abs(hm)

      do j=2,jl
      do i=2,il
         rh        = w(i,j,4)  +p(i,j)  -w(i,j,1)*h0
         f         = 1./(1.  +max(hm,0.)*rh/w(i,j,1))
         w(i,j,1)  = f*w(i,j,1)
         w(i,j,2)  = f*w(i,j,2)
         w(i,j,3)  = f*w(i,j,3)
         w(i,j,4)  = rh/(1.  +abs(hm))  -p(i,j)  +w(i,j,1)*h0
         rqq       = .5*(w(i,j,2)**2  +w(i,j,3)**2)/w(i,j,1)
         p(i,j)    = (gamma  -1.)*dim(w(i,j,4),rqq)
      end do
      end do

      end if
c
c     update the surface pressure
c
      call bcwall
c
c     update the far field boundary condition
c
      if (mode.eq.0.or.fbc.gt.0.) call bcfar
c
c     set extended values in the halo
c
      call halo

      if (mode.gt.0) return
c
c     option to check convergence
c
      hrms      = 0.
      hmax      = 0.
      ih        = 0
      jh        = 0
      rtrms     = 0.
      rtmax     = 0.
      irt       = 0
      jrt       = 0
      lmax      = nsup
      nsup      = 0

      do j=2,jl
      do i=2,il
         qq        = w(i,j,2)**2  +w(i,j,3)**2
         nsup      = nsup  +min(1.,qq/max(qq,(gamma*p(i,j)*w(i,j,1))))
         h         = (w(i,j,4)  +p(i,j))/w(i,j,1)  -h0
         hrms      = hrms  +h**2
c        dt        = cfl*rfl(i,j)*dtl(i,j)*vol(i,j)
c        dt        = dtl(i,j)*vol(i,j)
c        rt        = dw(i,j,1)/dt
         rt        = res(i,j)
         rtrms     = rtrms  +rt**2
      end do
      end do

      rtrms     = sqrt(rtrms/float(nx*ny))
      hrms      = sqrt(hrms/float(nx*ny))

      if (lmax.lt.0) return

      do j=2,jl
      do i=2,il
         h         = (w(i,j,4)  +p(i,j))/w(i,j,1)  -h0
         if (abs(h).gt.abs(hmax)) then
            hmax      = h
            ih        = i
            jh        = j
         end if
c        dt        = cfl*rfl(i,j)*dtl(i,j)*vol(i,j)
c        dt        = dtl(i,j)*vol(i,j)
c        rt        = dw(i,j,1)/dt
         rt        = res(i,j)
         if (abs(rt).gt.abs(rtmax)) then
            rtmax     = rt
            irt       = i
            jrt       = j
         end if
      end do
      end do

      return

      end
