      subroutine rbfield (jstep,varb)
c
c     read HelioSpace bfield data.
c
      use bfield_common
c
      integer   jstep
c
      character*1  cdum
      character*14 hsname
c
      integer   lbt
      integer   i, j, k, l, m
      real varb(mx,my,mz,lvars,mblks)
c
 1001 format ('bfield.000000',i1)
 1002 format ('bfield.00000',i2)
 1003 format ('bfield.0000',i3)
 1004 format ('bfield.000',i4)
 1005 format ('bfield.00',i5)
 1006 format ('bfield.0',i6)
 1007 format ('bfield.',i7)
 2001 format (a1)
 2002 format (10x,a1)
 2003 format (i2)
 3001 format (' nx = ',i2,' > mx = ',i2,' ... rbfield terminates',/)
 3002 format (' ny = ',i2,' > my = ',i2,' ... rbfield terminates',/)
 3003 format (' nz = ',i2,' > mz = ',i2,' ... rbfield terminates',/)
 3004 format (' npes = ',i4,' > mpes = ',i4,' ... rbfield terminates',/)
 3005 format (' nbpe = ',i3,' > mbpe = ',i3,' ... rbfield terminates',/)
c
c * * *
c
c     allocate bfield array
c         allocate (varb(mx,my,mz,lvars,mblks))
c
#ifdef AXISYM
        ndim = 2
#else
        ndim = 3
#endif
c
        open (unit=2, file='bfield.hdr', status='old')
c
c     read the header file.
c
        read (2,2001) cdum
        read (2,2001) cdum
c
        if (cdum .eq. 'C') then
          lcartsn = .true.
          lsphrcl = .false.
        else if (cdum .eq. 'S') then
          lsphrcl = .true.
          lcartsn = .false.
        end if
c
        lsphexp = .false.
        if (lsphrcl) then
          backspace (2)
          read (2,2002,err=100) cdum
          if (cdum .eq. 'E') then
            lsphexp = .true.
          end if
  100     continue
        end if
c
        if (lsphrcl) then
          pio2 = asin (1.0)
          pi   = 2.0 * pio2
        end if
c
        read (2,2003) n1p
        read (2,2003) n2p
        read (2,2003) n3p
c
        nstop = 0
        if (n1p .gt. mx) then
          write (*,3001) n1p, mx
          nstop = 1
        end if
        if (n2p .gt. my) then
          write (*,3002) n2p, my
          nstop = 1
        end if
        if (n3p .gt. mz) then
          write (*,3003) n3p, mz
          nstop = 1
        end if
c
        if (nstop .eq. 1) stop 3000
c
        nvars = 3
c
        if (jstep .le. 9) then
          write (hsname,1001) jstep
        else if (jstep .le. 99) then
          write (hsname,1002) jstep
        else if (jstep .le. 999) then
          write (hsname,1003) jstep
        else if (jstep .le. 9999) then
          write (hsname,1004) jstep
        else if (jstep .le. 99999) then
          write (hsname,1005) jstep
        else if (jstep .le. 999999) then
          write (hsname,1006) jstep
        else if (jstep .le. 9999999) then
          write (hsname,1007) jstep
        end if
c
        open (unit=3, file=hsname, form='unformatted',
     .        status='old',CONVERT='BIG_ENDIAN')
c
        npes           =  0
        lbpe(:,:)      = -1
        type(:)        = -1
        prnt(:,:)      = -1
        chld(:,:,:)    = -1
        lrefine(:)     = -1
        which_child(:) = -1
        neigh(:,:,:)   = -1
        lbno(:,:)      = -1
        nbpe(:)        =  0
c
c     read the data file.
c
        read (3) time
        read (3) ntotl, nleaf
c
c     read the tree data.
c
        do lbt = 1, ntotl
          read (3) (lbpe(l,lbt),l=1,2),
     .             type(lbt),
     .             (prnt(l,lbt),l=1,2),
     .             ((chld(l,m,lbt),l=1,2),m=1,8),
     .             lrefine(lbt),
     .             which_child(lbt),
     .             ((neigh(l,m,lbt),l=1,2),m=1,6)
          read (3) ((bndbox(l,m,lbt),l=1,2),m=1,3)
c
          npes                = max (lbpe(2,lbt)+1, npes)
          nbpe(lbpe(2,lbt)+1) = max (lbpe(1,lbt), nbpe(lbpe(2,lbt)+1))
c
c     set the reverse lookup array to find blocks in this list.
c
          lbno(lbpe(1,lbt),lbpe(2,lbt)+1) = lbt
c
c     read and process the leaf-block data.
c
          if (type(lbt) .eq. 1) then
c
            do k = 1, n3p
              do j = 1, n2p
                do i = 1, n1p
                  read (3) (vars(l), l=1,nvars)
                  varb(i,j,k,1,lbt) = vars(1)
                  varb(i,j,k,2,lbt) = vars(2)
                  varb(i,j,k,3,lbt) = vars(3)
                end do
              end do
            end do
c
          end if
        end do
c
        nbpemx = 0
        do lpes = 1, npes
          nbpemx = max (nbpemx, nbpe(lpes))
        end do
c
        nstop = 0
        if (npes .gt. mpes) then
          write (*,3004) npes, mpes
          nstop = 1
        end if
c
        if (nbpemx .gt. mbpe) then
          write (*,3005) nbpemx, mbpe
          nstop = 1
        end if
c
        if (nstop .eq. 1) stop 3000
c
        close (3)
        close (2)
c
      return
      end subroutine rbfield
