      subroutine boundary_per (out,pos0,pos1,bbox,xp,yp,zp)
c check if position is outside volume
c if so place on boundary
c Peter Wyper Feb 2015
c update for periodic boundaries - June 2020
c update to pos(9) for qsl - June 2023
c      use bfield_common
c inputs
      integer   out
      real      pos0(9),pos1(9),bbox(2,3)
      logical   xp,yp,zp

c internal
      integer   xlo,xhi,ylo,yhi,zlo,zhi
      real      dvec(9),tolbfd2
c            
c * * *
c
        out = 0
c test whether the new point is outside volume or on boundaries 
        xlo = 0
        xhi = 0
        ylo = 0
        yhi = 0
        zlo = 0
        zhi = 0

        if (bbox(1,1) .gt. pos1(1)) then
           if (xp) then
c periodic in x-direction - do nothing 
           else
           xlo = 1
           endif
        endif
        if (bbox(2,1) .lt. pos1(1)) then
           if (xp) then
           else
           xhi = 1
           endif
        endif
        if (bbox(1,2) .gt. pos1(2)) then
           if (yp) then
           else
           ylo = 1
           endif
        endif
        if (bbox(2,2) .lt. pos1(2)) then
           if (yp) then
           else
           yhi = 1
           endif
        endif
        if (bbox(1,3) .gt. pos1(3)) then
           if (zp) then
           else
           zlo = 1
           endif
        endif
        if (bbox(2,3) .lt. pos1(3)) then
           if (zp) then
           else
           zhi = 1
           endif
        endif

        
        if (xlo+xhi+ylo+yhi+zlo+zhi .eq. 0) then
           return
        endif



        if (xlo+xhi+ylo+yhi+zlo+zhi .gt. 0) then
           out = 1
           dvec(1:3) = pos1(1:3)-pos0(1:3)

           if (xlo  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(1,1)-pos1(1))/dvec(1)*dvec(1:3) !outside - interpolate
           endif
           if (xhi  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(2,1)-pos1(1))/dvec(1)*dvec(1:3)
           endif

           if (ylo  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(1,2)-pos1(2))/dvec(2)*dvec(1:3)
           endif
           if (yhi  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(2,2)-pos1(2))/dvec(2)*dvec(1:3)
           endif
             
           if (zlo  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(1,3)-pos1(3))/dvec(3)*dvec(1:3)
           endif

           if (zhi  .eq. 1) then
           pos1(1:3) = pos1(1:3) + (bbox(2,3)-pos1(3))/dvec(3)*dvec(1:3)
           endif

        endif



      end subroutine boundary_per
