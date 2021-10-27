      module solv_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the solution             *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:,:,:)    , allocatable :: wn,dw,fw,vw,fs
      real, dimension(:,:,:)    , allocatable :: w1,wr
      real, dimension(:,:)      , allocatable :: dtl,rfl
      real, dimension(:,:)      , allocatable :: rfli,rflj
      real, dimension(:,:)      , allocatable :: radi,radj

      end module solv_var
