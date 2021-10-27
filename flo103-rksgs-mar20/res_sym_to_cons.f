      subroutine  res_sym_to_cons
c
c     nonconservative residuals dw changed to conservative residuals
c
      use dims
c
c     ******************************************************************
c
      use solv_var
c
c     ******************************************************************
c
      use flo_var
c
c     ******************************************************************
c
      use psm_var
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

        dw(i,j,1) = r1  +r4
        dw(i,j,2) = ua*(r1  +r4)  +c*r2
        dw(i,j,3) = va*(r1  +r4)  +c*r3
        dw(i,j,4) = ha*r1  +c*(ua*r2  +va*r3) +qq*r4

      end do
      end do

      return

      end
