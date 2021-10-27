      module mg_param
c
c     ******************************************************************
c     *                                                                *
c     *   parameters controlling the multigrid process                 *
c     *                                                                *
c     ******************************************************************
c
      integer   :: kode,mode,lcyc
      integer   :: k1,kx,kw

      real      :: fcoll,fadd,fbc
      real      :: cflf,cflc,cfl0,hmf,hmc
      real      :: epsf,epsc

      end module mg_param
