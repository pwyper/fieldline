      module bfield_common

      save
      parameter (mbpe=300,mpes=2048)
      parameter (mblks=mbpe*mpes)
      parameter (lvars=3)
      parameter (mx=9, my=9, mz=9)

      logical   lcartsn, lsphrcl, lsphexp

      real      bndbox(2,3,mblks)
      real      vars(lvars)
!      real      varb(mx,my,mz,lvars,mblks)
      real      coord1(mx), coord2(my), coord3(mz)

      real      time
      real      pi, pio2
      real      tolbfd

      integer   lbpe(2,mblks)
      integer   type(mblks)
      integer   prnt(2,mblks)
      integer   chld(2,8,mblks)
      integer   lrefine(mblks)
      integer   which_child(mblks)
      integer   neigh(2,6,mblks)
      integer   lbno(mbpe,mpes)
      integer   data(mblks)
      integer   nbpe(mpes)

      integer   nleaf, ntotl
      integer   nvars, ivars
      integer   n1p, n2p, n3p
      integer   ndim

      end module bfield_common
