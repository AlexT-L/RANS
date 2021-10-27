      subroutine syslock
c
c     ******************************************************************
c     *                                                                *
c     *   check license expiration date set in gettd.c at compile time *
c     *                                                                *
c     *   check hardware address of the machine this program is        *
c     *    being run on.  allowed machines are set in haddr.c          *
c     *                                                                *
c     *   terminate program execution if necessary                     *
c     *                                                                *
c     ******************************************************************
c
      implicit none
      integer  lock,iproc
c
c     ******************************************************************
c
c     check license termination date
c
      iproc = 1
      call gettd(lock,iproc)
c
      if (lock.ne.0) stop
c
c     check machine hardware address
c
      call haddr(lock,iproc)
c
      if (lock.ne.0) stop
c
      return
      end
