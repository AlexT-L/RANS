      module solv_param
c
c     ******************************************************************
c     *                                                                *
c     *   parameters controlling the solution process                  *
c     *                                                                *
c     ******************************************************************
c
      integer   :: mcyc,ncyc,ntim,nout,nprnt,iprnt,lprnt,nmesh
      integer   :: irt,jrt,ih,jh,nsup
      integer   :: mstage,nstage,ksmoop
      integer   :: nres,imon
      integer   :: lmesh,mmesh
      integer   :: isym

      integer   :: iprec

      real      :: cstp(6),cdis(6)
      real      :: cfl,cflim,vt,hm,vis0,vis2,vis4,adis,qdis,bc,rfil,gcl
      real      :: smoopi,smoopj,rfl0,csmoop

      real      :: eps,diag

      real      :: rtrms,rtmax,hrms,hmax
      real      :: tot,rtrms0

      real      :: fcyc(4),fprnt(4),fout(4),ftim(4)
      real      :: gprnt(4),hprnt(4),hmesh(4)

      end module solv_param
