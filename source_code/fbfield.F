      subroutine fbfield (jsearch, lsearch, rsearch, bsearch, fail)
c
c     find HelioSpace bfield data.
c
      use bfield_common
c
      integer   jsearch, lsearch
      integer   fail
c
      real      rsearch(3), bsearch(2,3)
c
      integer   lbt,    lbn
      integer   ibound, jbound, kbound, lbound
      integer   nrlb,   nrpe
      integer   chlb,   chpe
      integer   ptlb,   ptpe
c
 4001 format (' fbfield finds neighboring leaf & grandparent blocks',
     .        ' ... terminating')
 4002 format (' fbfield finds no neighbor to parent block',
     .        ' ... terminating')
c
c * * *
c
c     start a search from scratch.
c
      fail = 0 

        if (jsearch .eq. 0) then
          lsearch = -1
c
c     knucklehead search of the tree data.
c
          do lbt = 1, ntotl
            if (type(lbt) .eq. 1) then
              if (bndbox(1,1,lbt) .gt. rsearch(1)) go to 1000
              if (bndbox(2,1,lbt) .lt. rsearch(1)) go to 1000
              if (bndbox(1,2,lbt) .gt. rsearch(2)) go to 1000
              if (bndbox(2,2,lbt) .lt. rsearch(2)) go to 1000
              if (bndbox(1,3,lbt) .gt. rsearch(3)) go to 1000
              if (bndbox(2,3,lbt) .lt. rsearch(3)) go to 1000

c
c     success!
c
              lsearch      = lbt
              bsearch(1,1) = bndbox(1,1,lbt)
              bsearch(2,1) = bndbox(2,1,lbt)
              bsearch(1,2) = bndbox(1,2,lbt)
              bsearch(2,2) = bndbox(2,2,lbt)
              bsearch(1,3) = bndbox(1,3,lbt)
              bsearch(2,3) = bndbox(2,3,lbt)
c
              go to 1001
c
            end if
 1000       continue
          end do
 1001     continue
c

          return
c
        else if (jsearch .eq. 1) then
c
c     smart search for neighbor block.
c
          ibound = 0
          jbound = 0
          kbound = 0
c
          if (abs(rsearch(1)-bsearch(1,1)) .lt. tolbfd) then
            ibound = -1
c     found a domain boundary.  success!
            if (neigh(1,1,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          else if (abs(rsearch(1)-bsearch(2,1)) .lt. tolbfd) then
            ibound = +1
c     found a domain boundary.  success!
            if (neigh(1,2,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          end if
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
c     found a domain boundary.  success!
            if (neigh(1,3,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
c     found a domain boundary.  success!
            if (neigh(1,4,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
c     found a domain boundary.  success!
            if (neigh(1,5,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
c     found a domain boundary.  success!
            if (neigh(1,6,lsearch) .le. -20) then
              lsearch = 0
              return
            end if
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
c     1st coord face search.
c
          if (abs(ibound) .eq. 1) then
            lbt     = lsearch
            lsearch = -1
c
c     left.
c
            if (ibound .eq. -1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,1,lbt) .gt. 0) then
                nrlb = neigh(1,1,lbt)
                nrpe = neigh(2,1,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage2 = 0.5 * (bndbox(1,2,lbn) + bndbox(2,2,lbn))
                  baverage3 = 0.5 * (bndbox(1,3,lbn) + bndbox(2,3,lbn))
                  if (rsearch(3) .le. baverage3) then
                    if (rsearch(2) .le. baverage2) then
                      chlb = chld(1,2,lbn)
                      chpe = chld(2,2,lbn)
                    else
                      chlb = chld(1,4,lbn)
                      chpe = chld(2,4,lbn)
                    end if
                  else
                    if (rsearch(2) .le. baverage2) then
                      chlb = chld(1,6,lbn)
                      chpe = chld(2,6,lbn)
                    else
                      chlb = chld(1,8,lbn)
                      chpe = chld(2,8,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,1,lbp)
                nrpe = neigh(2,1,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                  fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
c     right.
c
            else if (ibound .eq. +1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,2,lbt) .gt. 0) then
                nrlb = neigh(1,2,lbt)
                nrpe = neigh(2,2,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage2 = 0.5 * (bndbox(1,2,lbn) + bndbox(2,2,lbn))
                  baverage3 = 0.5 * (bndbox(1,3,lbn) + bndbox(2,3,lbn))
                  if (rsearch(3) .le. baverage3) then
                    if (rsearch(2) .le. baverage2) then
                      chlb = chld(1,1,lbn)
                      chpe = chld(2,1,lbn)
                    else
                      chlb = chld(1,3,lbn)
                      chpe = chld(2,3,lbn)
                    end if
                  else
                    if (rsearch(2) .le. baverage2) then
                      chlb = chld(1,5,lbn)
                      chpe = chld(2,5,lbn)
                    else
                      chlb = chld(1,7,lbn)
                      chpe = chld(2,7,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,2,lbp)
                nrpe = neigh(2,2,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          ibound = 0
c
          if (abs(rsearch(2)-bsearch(1,2)) .lt. tolbfd) then
            jbound = -1
          else if (abs(rsearch(2)-bsearch(2,2)) .lt. tolbfd) then
            jbound = +1
          else
            jbound = 0
          end if
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                   fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
            end if
c
          end if
c
c     2nd coord face search.
c
          if (abs(jbound) .eq. 1) then
            lbt     = lsearch
            lsearch = -1
c
c     left.
c
            if (jbound .eq. -1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,3,lbt) .gt. 0) then
                nrlb = neigh(1,3,lbt)
                nrpe = neigh(2,3,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage1 = 0.5 * (bndbox(1,1,lbn) + bndbox(2,1,lbn))
                  baverage3 = 0.5 * (bndbox(1,3,lbn) + bndbox(2,3,lbn))
                  if (rsearch(3) .le. baverage3) then
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,3,lbn)
                      chpe = chld(2,3,lbn)
                    else
                      chlb = chld(1,4,lbn)
                      chpe = chld(2,4,lbn)
                    end if
                  else
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,7,lbn)
                      chpe = chld(2,7,lbn)
                    else
                      chlb = chld(1,8,lbn)
                      chpe = chld(2,8,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,3,lbp)
                nrpe = neigh(2,3,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                   fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
c     right.
c
            else if (jbound .eq. +1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,4,lbt) .gt. 0) then
                nrlb = neigh(1,4,lbt)
                nrpe = neigh(2,4,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage1 = 0.5 * (bndbox(1,1,lbn) + bndbox(2,1,lbn))
                  baverage3 = 0.5 * (bndbox(1,3,lbn) + bndbox(2,3,lbn))
                  if (rsearch(3) .le. baverage3) then
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,1,lbn)
                      chpe = chld(2,1,lbn)
                    else
                      chlb = chld(1,2,lbn)
                      chpe = chld(2,2,lbn)
                    end if
                  else
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,5,lbn)
                      chpe = chld(2,5,lbn)
                    else
                      chlb = chld(1,6,lbn)
                      chpe = chld(2,6,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,4,lbp)
                nrpe = neigh(2,4,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          jbound = 0
c
          if (abs(rsearch(3)-bsearch(1,3)) .lt. tolbfd) then
            kbound = -1
          else if (abs(rsearch(3)-bsearch(2,3)) .lt. tolbfd) then
            kbound = +1
          else
            kbound = 0
          end if
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                   fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
            end if
c
          end if
c
c     3rd coord face search.
c
          if (abs(kbound) .eq. 1) then
            lbt     = lsearch
            lsearch = -1
c
c     left.
c
            if (kbound .eq. -1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,5,lbt) .gt. 0) then
                nrlb = neigh(1,5,lbt)
                nrpe = neigh(2,5,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage1 = 0.5 * (bndbox(1,1,lbn) + bndbox(2,1,lbn))
                  baverage2 = 0.5 * (bndbox(1,2,lbn) + bndbox(2,2,lbn))
                  if (rsearch(2) .le. baverage2) then
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,5,lbn)
                      chpe = chld(2,5,lbn)
                    else
                      chlb = chld(1,6,lbn)
                      chpe = chld(2,6,lbn)
                    end if
                  else
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,7,lbn)
                      chpe = chld(2,7,lbn)
                    else
                      chlb = chld(1,8,lbn)
                      chpe = chld(2,8,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,5,lbp)
                nrpe = neigh(2,5,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                   fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
c     right.
c
            else if (kbound .eq. +1) then
c
c     a neighbor block exists at the same level.
c
              if (neigh(1,6,lbt) .gt. 0) then
                nrlb = neigh(1,6,lbt)
                nrpe = neigh(2,6,lbt)
                lbn  = lbno(nrlb,nrpe+1)
c
c     the neighbor block is a leaf block.  success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a parent block.  what child is this?
c
                else if (type(lbn) .eq. 2) then
                  baverage1 = 0.5 * (bndbox(1,1,lbn) + bndbox(2,1,lbn))
                  baverage2 = 0.5 * (bndbox(1,2,lbn) + bndbox(2,2,lbn))
                  if (rsearch(2) .le. baverage2) then
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,1,lbn)
                      chpe = chld(2,1,lbn)
                    else
                      chlb = chld(1,2,lbn)
                      chpe = chld(2,2,lbn)
                    end if
                  else
                    if (rsearch(1) .le. baverage1) then
                      chlb = chld(1,3,lbn)
                      chpe = chld(2,3,lbn)
                    else
                      chlb = chld(1,4,lbn)
                      chpe = chld(2,4,lbn)
                    end if
                  end if
c
c     success!
c
                  lbc          = lbno(chlb,chpe+1)
                  lsearch      = lbc
                  bsearch(1,1) = bndbox(1,1,lbc)
                  bsearch(2,1) = bndbox(2,1,lbc)
                  bsearch(1,2) = bndbox(1,2,lbc)
                  bsearch(2,2) = bndbox(2,2,lbc)
                  bsearch(1,3) = bndbox(1,3,lbc)
                  bsearch(2,3) = bndbox(2,3,lbc)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     the neighbor block is a grandparent block (not allowed; stop).
c
                else if (type(lbn) .eq. 2) then
                  write (*,4001)
                  stop 4001
                end if
c
c     no neighbor block exists at the same level.  find the parent's neighbor.
c
              else
                ptlb = prnt(1,lbt)
                ptpe = prnt(2,lbt)
                lbp  = lbno(ptlb,ptpe+1)
c
                nrlb = neigh(1,6,lbp)
                nrpe = neigh(2,6,lbp)
                lbn  = lbno(nrlb,nrpe+1)
c
c     success!
c
                if (type(lbn) .eq. 1) then
                  lsearch      = lbn
                  bsearch(1,1) = bndbox(1,1,lbn)
                  bsearch(2,1) = bndbox(2,1,lbn)
                  bsearch(1,2) = bndbox(1,2,lbn)
                  bsearch(2,2) = bndbox(2,2,lbn)
                  bsearch(1,3) = bndbox(1,3,lbn)
                  bsearch(2,3) = bndbox(2,3,lbn)
c
          kbound = 0
c
          lbound = abs(ibound) + abs(jbound) + abs(kbound)
c
                  if (lbound .eq. 0) return
c
c     oops.
c
                else
                   fail = 1
c                  write (*,4002)
c                  stop 4002
                end if
              end if
c
            end if
c
          end if
c
        end if
c
      return
      end subroutine fbfield
