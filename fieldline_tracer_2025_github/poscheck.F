      subroutine poscheck (pos0,pos1,bbox)
c check if position is outside volume for periodic data
c if so produce new position back inside volume
c assumes in adjacent periodic domain 
c Peter Wyper June 202
c updated to pos0(9) and pos1(9) for Tassev Q method - June 2023
c
c inputs
      real      pos0(10),bbox(2,3)
c outputs
      real      pos1(10)

c internal
c      integer   xlo,xhi,ylo,yhi,zlo,zhi
c      real      dvec(3),tolbfd2
c            
c * * *
c
      pos1(:) = pos0(:)

        if (bbox(1,1) .gt. pos0(1)) then
           pos1(1) = pos0(1) + (bbox(2,1)-bbox(1,1))
        endif
        if (bbox(2,1) .lt. pos0(1)) then
           pos1(1) = pos0(1) - (bbox(2,1)-bbox(1,1))
        endif
        if (bbox(1,2) .gt. pos0(2)) then
           pos1(2) = pos0(2) + (bbox(2,2)-bbox(1,2))
        endif
        if (bbox(2,2) .lt. pos0(2)) then
           pos1(2) = pos0(2) - (bbox(2,2)-bbox(1,2))
        endif
        if (bbox(1,3) .gt. pos0(3)) then
           pos1(3) = pos0(3) + (bbox(2,3)-bbox(1,3))
        endif
        if (bbox(2,3) .lt. pos0(3)) then
           pos1(3) = pos0(3) - (bbox(2,3)-bbox(1,3))
        endif

        return


      end subroutine poscheck
