      subroutine  res_cons_to_sym
c
c     conservative residuals dw changed to nonconservative residuals
c
      use dims
c
c     ******************************************************************
c
      use flo_var
c
c     ******************************************************************
c
      use solv_var
c
c     ******************************************************************
c
      use flo_param
c
c     ******************************************************************
c
      implicit none
      integer  :: i,j
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
        dw(i,j,4) = r1  -dw(i,j,1)

      end do
      end do

      return

      end
