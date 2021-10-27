      subroutine  bcresj
c
c*****   set boundary values to zero and implement cut bc's
c*****  for face values of residuals
c
      use dims
      use psm_var

      implicit none
      integer  :: i,j,ii
c
      do 10 i= 1,ie
        rs(i, 1,1) = 0.0
        rs(i, 1,2) = 0.0
        rs(i, 1,3) = 0.0
        rs(i, 1,4) = 0.0
c
        rs(i,je,1) = 0.0
        rs(i,je,2) = 0.0
        rs(i,je,3) = 0.0
        rs(i,je,4) = 0.0
   10 continue
c
      do 20 j= 1,je
        rs( 1,j,1) = 0.0
        rs( 1,j,2) = 0.0
        rs( 1,j,3) = 0.0
        rs( 1,j,4) = 0.0
c
        rs(ie,j,1) = 0.0
        rs(ie,j,2) = 0.0
        rs(ie,j,3) = 0.0
        rs(ie,j,4) = 0.0
   20 continue
c
      do 30 i= 1,ie
         if (i.le.itl .or. i.gt.itu) then
            ii    = ie-i+1
            rs(i,1,1) = rs(ii,2,1)
            rs(i,1,2) = rs(ii,2,2)
            rs(i,1,3) = rs(ii,2,3)
            rs(i,1,4) = rs(ii,2,4)
         end if
   30 continue
c
      return
      end
