      subroutine next_block (pos1,bsearch,rsearch,nb)
c
c     decide where to place footpoint for fbfield
c     Peter Wyper Jan 2015
c
c inputs
      real      pos1(3), bsearch(2,3)
c outputs
      real      rsearch(3)
      integer   nb
c internal 
      real      xb,yb,zb
      integer   xlo,xhi,ylo,yhi,zlo,zhi
c
c * * *
c
c test whether the new point is in a block on a face, vertex or corner
        xlo = 0
        xhi = 0
        ylo = 0
        yhi = 0
        zlo = 0
        zhi = 0
        
        if (bsearch(1,1) .gt. pos1(1)) then
           xlo = 1
        endif
        if (bsearch(2,1) .lt. pos1(1)) then
           xhi = 1
        endif
        if (bsearch(1,2) .gt. pos1(2)) then
           ylo = 1
        endif
        if (bsearch(2,2) .lt. pos1(2)) then
           yhi = 1
        endif
        if (bsearch(1,3) .gt. pos1(3)) then
           zlo = 1
        endif
        if (bsearch(2,3) .lt. pos1(3)) then
           zhi = 1
        endif


        if (xlo+xhi+ylo+yhi+zlo+zhi .eq. 0) then
c still in the old block
           nb = 0
           return
        endif


        if (xlo+xhi+ylo+yhi+zlo+zhi .eq. 1) then
c    put footpoint to feed into fbfield on box face (opposite new point)
           if (xlo+xhi .eq. 1) then
              yb = pos1(2)
              zb = pos1(3)
              if (xlo .eq. 1) then
                 xb = bsearch(1,1)
              else
                 xb = bsearch(2,1)
              endif
           endif
           if (ylo+yhi .eq. 1) then
              xb = pos1(1)
              zb = pos1(3)
              if (ylo .eq. 1) then
                 yb = bsearch(1,2)
              else
                 yb = bsearch(2,2)
              endif
           endif
           if (zlo+zhi .eq. 1) then
              xb = pos1(1)
              yb = pos1(2)
              if (zlo .eq. 1) then
                 zb = bsearch(1,3)
              else
                 zb = bsearch(2,3)
              endif
           endif

           rsearch(1) = xb
           rsearch(2) = yb
           rsearch(3) = zb
           nb = 1
           return
        endif
c
        if (xlo+xhi+ylo+yhi+zlo+zhi .eq. 2) then
c     put footpoint to feed to fbfield on box vertex (next to new point)           

           if (xlo+xhi .eq. 0) then
              xb = pos1(1)
              
              if (ylo .eq. 1) then
                 yb = bsearch(1,2)
              else
                 yb = bsearch(2,2)
              endif
              
              if (zlo .eq. 1) then
                 zb = bsearch(1,3)
              else
                 zb = bsearch(2,3)
              endif
           endif
           
           if (ylo+yhi .eq. 0) then
              yb = pos1(2)
              
              if (xlo .eq. 1) then
                 xb = bsearch(1,1)
              else
                 xb = bsearch(2,1)
              endif
              if (zlo .eq. 1) then
                 zb = bsearch(1,3)
              else
                 zb = bsearch(2,3)
              endif
           endif

           if (zlo+zhi .eq. 0) then 
              zb = pos1(3)
              
              if (xlo .eq. 1) then
                 xb = bsearch(1,1)
              else
                 xb = bsearch(2,1)
              endif
              if (ylo .eq. 1) then
                 yb = bsearch(1,2)
              else
                 yb = bsearch(2,2)
              endif
           endif

           rsearch(1) = xb
           rsearch(2) = yb
           rsearch(3) = zb
           nb = 1
           return
        endif

        if (xlo+xhi+ylo+yhi+zlo+zhi .eq. 3) then
c put footpoint to feed to fbfield in box corner
           
           if (xlo .eq. 1) then
              xb = bsearch(1,1)
           else
              xb = bsearch(2,1)
           endif

           if (ylo .eq. 1) then
              yb = bsearch(1,2)
           else
              yb = bsearch(2,2)
           endif

           if (zlo .eq. 1) then
              zb = bsearch(1,3)
           else
              zb = bsearch(2,3)
           endif

           rsearch(1) = xb
           rsearch(2) = yb
           rsearch(3) = zb
           nb = 1
           return
        endif
 
 

      end subroutine next_block
