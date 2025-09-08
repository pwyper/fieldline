      subroutine rk4 (pos,n1,len,lmax,ep,bbox,varb,xp,yp,zp,bv1,sph,apt)
c
c     4th order Runge Kutta field line tracer - periodic boundaries
c                      
      use bfield_common
c  

      real      rsearch(3), bsearch(2,3)
      real      xmin, ymin, zmin
      real      xmax, ymax, zmax
      real      xfpt, yfpt, zfpt
      real      xfpt2, yfpt2, zfpt2
      real      pos1(10),pos(n1,10)
      real      xb,yb,zb
      logical   sph, apt
c                     
      integer   jsearch, lsearch, ls0
      integer   jstep
      integer   lbt,nb
      integer   n,n1
      real      bbox(2,3),ep
      integer   out
      integer   fail
c
      integer   xlo,xhi,ylo,yhi,zlo,zhi
      real      eighth,eight,one,two,half,six,sixth
      real      dv(3), corner(3), bv(3)
      real      blen,len1,len,lenb,lmax
      real      k1(9),k2(9),k3(9),k4(9)
      real      delta(3),err1
      real      uvec(3), vvec(3),bhat(3)
c
      logical   xp,yp,zp
      real      posp(10),bv1(3)
      real      tmp, r
c
      real varb(mx,my,mz,lvars,mblks)
c
c * * *
c
        eight = real (8.0)
        six = real (6.0)
        one = real (1.0)
        two = real (2.0)
        four = real(4.0)
        zero = real(0.0)
        eighth = one/eight
        half = one/two
        sixth = one/six
c
        jsearch = 0
        lsearch = -1
        bsearch(:,:) = 0.0
        fail = 0
c 
c should be in volume from initialisation - don't need poscheck if non-periodic
c
        call poscheck (pos(1,:),posp,bbox)
c posp - insde the volume for purposes of interpolation
c pos(1,1:3) - still outside to give correct value for Q etc
        call fbfield (0, lsearch, posp(1:3), bsearch,fail)
        fail = 0
        jsearch = 1
c block position for 1st position
        ls0 = lsearch 
        dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
        corner(:) = bsearch(1,:)
c------------------------------------------------------------------
c                   initialise U and V vectors
c------------------------------------------------------------------
        call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp(1:3),bv)
        blen = sqrt(bv(1)**2.0+bv(2)**2.0+bv(3)**2.0)
        bhat(:) = bv(:)/blen
c set so vvec.bhat = 0 while avoiding zero length for vvec
        if (abs(bhat(3)) .eq. maxval(abs(bhat))) then
           vvec(1) = bhat(3)
           vvec(2) = zero
           vvec(3) = -bhat(1)
        else
           vvec(1) = bhat(2)
           vvec(2) = -bhat(1)
           vvec(3) = zero
        endif
c normalise
        vvec(:) = vvec(:)/sqrt(vvec(1)**2.0+vvec(2)**2.0+vvec(3)**2.0)
c set uvec = cross(vvec,bhat)
        uvec(1) = vvec(2)*bhat(3)-vvec(3)*bhat(2)
        uvec(2) = vvec(3)*bhat(1)-vvec(1)*bhat(3)
        uvec(3) = vvec(1)*bhat(2)-vvec(2)*bhat(1)
c normalise to be sure (should be unit length anyway)
        uvec(:) = uvec(:)/sqrt(uvec(1)**2.0+uvec(2)**2.0+uvec(3)**2.0)
c place into pos vector to be advanced
c pos(:,1:3) - position
c pos(:,4:6) - uvec
c pos(:,7:9) - vvec
        pos(1,4:6) = uvec(:)
        pos(1,7:9) = vvec(:)
        posp(4:6) = uvec(:)
        posp(7:9) = vvec(:)
c----------------------------------------------------------------
c               initialise twist calc
c----------------------------------------------------------------
c pos(:,10) = curl(B).B/(4pi*B^2)
        if (sph) then
           call tw_sph (varb(:,:,:,:,lsearch),dv,corner,posp)
           pos(1,10) = posp(10)           
        else
           call tw_cart (varb(:,:,:,:,lsearch),dv,corner,posp)
           pos(1,10) = posp(10)
        endif
c----------------------------------------------------------------
c
c step size set by: 
c lmax - externally set max
c lenb - 4 cells/half block wide
c len  - starting length
c len1 - above combined +  4rth/3rd order check 
c stops jumping blocks and loosing high res information

        lenb = minval(dv*four)
        len1 = len !+1000.0

        out = 0
      do n=1,n1-1
c         print *,'n = ',n 
         if (out .eq. 0) then
c ---------------------------------------------------------------
c set step length
            if (apt) then
               len1 = minval([abs(len1),lenb,lmax])*sign(1.0,len)
            endif
c     len1 = minval([lmax])*sign(1.0,len)
c --------------------------------------------------------------- 
      call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp(1:3),bv)

      blen = sqrt(bv(1)**2.0+bv(2)**2.0+bv(3)**2.0)
      if (blen .ne. 0.0) then
c     cartesian: dx/dl = B/|B|, sphercial: dr/dl = Br/|B|, dth/dl = Bth/r|B|, dphi/dl=Bphi/(rsin(th)|B|)
c     But tracing in ln(r) space: dr/dl = dr/du*du/dl = r*du/dl = Br/|B| where u = ln(r) so 1/r factor 
c         print *,'posp(1), sin(posp(2)) = ',posp(1), sin(posp(2))
         if (sph) then
            r = exp(posp(1))
            k1(1:3) = len1*bv(:)/blen/r
            k1(3) = k1(3)/sin(posp(2))
         else
            k1(1:3) = len1*bv(:)/blen
         endif
      else
         k1(1:3) = 0.0
      endif

c find k1 for uvec and vvec
      if (sph) then
         call qsl_sph (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k1)
      else
         call qsl_cart (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k1)
      endif

      pos1(1:9) = pos(n,1:9)+half*k1(1:9)
      call boundary_per (out,pos(n,:),pos1(:),bbox,xp,yp,zp)   
      if (out .eq. 1) then
         go to 4000
      endif

      call poscheck (pos1,posp,bbox)
      call next_block(posp(1:3),bsearch,rsearch,nb)
      if (nb .eq. 1) then
         call fbfield (1, lsearch, rsearch, bsearch,fail)

         if (posp(1) .gt. bsearch(2,1)) fail = 1 !rsearch
         if (posp(1) .lt. bsearch(1,1)) fail = 1
         if (posp(2) .gt. bsearch(2,2)) fail = 1
         if (posp(2) .lt. bsearch(1,2)) fail = 1
         if (posp(3) .gt. bsearch(2,3)) fail = 1
         if (posp(3) .lt. bsearch(1,3)) fail = 1

         if (fail .eq. 1) then
            call fbfield (0, lsearch, posp(1:3), bsearch,fail)
            fail = 0
         endif
         if (lsearch .eq. 0) then
            out = 1
            go to 4000
         endif
         dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
         corner(:) = bsearch(1,:)
      endif
c ---------------------------------------------------------------
      call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp(1:3),bv)
      blen = sqrt(bv(1)**2.0+bv(2)**2.0+bv(3)**2.0)
      if (blen .ne. 0.0) then
c         k2(1:3) = bv(:)*len1/blen
         if (sph) then
            r = exp(posp(1))
            k2(1:3) = len1*bv(:)/blen/r
            k2(3) = k2(3)/sin(posp(2))
c     k2(1) = len1*bv(1)/blen/r
c            k2(2) = len1*bv(2)/blen/r
c            k2(3) = len1*bv(3)/blen/(r*sin(posp(2)))
         else
            k2(1:3) = len1*bv(:)/blen
         endif
      else
         k2(1:3) = 0.0
      endif

c find k2 for uvec and vvec                                       
      if (sph) then
         call qsl_sph (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k2)
      else
         call qsl_cart (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k2)
      endif

      pos1(1:9) = pos(n,1:9)+half*k2(1:9)

      call boundary_per (out,pos(n,:),pos1(:),bbox,xp,yp,zp)
      if (out .eq. 1) then
         go to 4000
      endif

      call poscheck (pos1,posp,bbox)
      call next_block(posp(1:3),bsearch,rsearch,nb)
      if (nb .eq. 1) then
         call fbfield (1, lsearch, rsearch, bsearch,fail)
         if (posp(1) .gt. bsearch(2,1)) fail = 1 !rsearch
         if (posp(1) .lt. bsearch(1,1)) fail = 1
         if (posp(2) .gt. bsearch(2,2)) fail = 1
         if (posp(2) .lt. bsearch(1,2)) fail = 1
         if (posp(3) .gt. bsearch(2,3)) fail = 1
         if (posp(3) .lt. bsearch(1,3)) fail = 1
         if (fail .eq. 1) then
            call fbfield (0, lsearch, posp(1:3), bsearch,fail)
            fail = 0
         endif
         if (lsearch .eq. 0) then
            out = 1
            go to 4000
         endif
         dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
         corner(:) = bsearch(1,:)
      endif
c ---------------------------------------------------------------
      call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp,bv)
      blen = sqrt(bv(1)**2.0+bv(2)**2.0+bv(3)**2.0)

      if (blen .ne. 0.0) then
c         k3(1:3) = bv(:)*len1/blen
         if (sph) then
            r = exp(posp(1))
            k3(1:3) = len1*bv(:)/blen/r
            k3(3) = k3(3)/sin(posp(2))
c     k3(1) = len1*bv(1)/blen/r
c            k3(2) = len1*bv(2)/blen/r
c            k3(3) = len1*bv(3)/blen/(r*sin(posp(2)))
         else
            k3(1:3) = len1*bv(:)/blen
         endif
      else
         k3(1:3) = 0.0
      endif

c find k3 for uvec and vvec                                       
      if (sph) then
         call qsl_sph (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k3)
      else
         call qsl_cart (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k3)
      endif

      pos1(1:9) = pos(n,1:9) + k3(1:9)
      call boundary_per (out,pos(n,:),pos1(:),bbox,xp,yp,zp)
      if (out .eq. 1) then
         go to 4000
      endif

      call poscheck (pos1,posp,bbox)
      call next_block(posp(1:3),bsearch,rsearch,nb)

      if (nb .eq. 1) then
         call fbfield (1, lsearch, rsearch, bsearch,fail)

      if (posp(1) .gt. bsearch(2,1)) fail = 1 !rsearch(1) etc
      if (posp(1) .lt. bsearch(1,1)) fail = 1
      if (posp(2) .gt. bsearch(2,2)) fail = 1
      if (posp(2) .lt. bsearch(1,2)) fail = 1
      if (posp(3) .gt. bsearch(2,3)) fail = 1
      if (posp(3) .lt. bsearch(1,3)) fail = 1

      if (fail .eq. 1) then
         call fbfield (0, lsearch, posp(1:3), bsearch,fail)
         fail = 0 
      endif

         if (lsearch .eq. 0) then
            out = 1
            go to 4000
         endif
         dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
         corner(:) = bsearch(1,:)
      endif
c ---------------------------------------------------------------
      call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp(1:3),bv)

      blen = sqrt(bv(1)**2.0+bv(2)**2.0+bv(3)**2.0)
      if (blen .ne. 0.0) then
c         k4(1:3) = bv(:)*len1/blen
         if (sph) then
            r = exp(posp(1))
            k4(1:3) = len1*bv(:)/blen/r
            k4(3) = k4(3)/sin(posp(2))
c            k4(1) = len1*bv(1)/blen/r
c            k4(2) = len1*bv(2)/blen/r
c            k4(3) = len1*bv(3)/blen/(r*sin(posp(2)))
         else
            k4(1:3) = len1*bv(:)/blen
         endif
      else
         k4(1:3) = 0.0
      endif

c find k4 for uvec and vvec                                       
      if (sph) then
         call qsl_sph (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k4)
      else
         call qsl_cart (varb(:,:,:,:,lsearch),dv,corner,posp,len1,k4)
      endif

      pos1(1:9) = pos(n,1:9) + (k1(1:9)+two*k2(1:9))*sixth
      pos1(1:9) = pos1(1:9) + (two*k3(1:9)+k4(1:9))*sixth 


      call boundary_per (out,pos(n,:),pos1(:),bbox,xp,yp,zp)
      if (out .eq. 1) then
         go to 4000
      endif

c set as next position unless on boundary 
      pos(n+1,:) = pos1(:)
      call poscheck (pos(n+1,:),posp,bbox)
      call next_block(posp(1:3),bsearch,rsearch,nb)
      if (nb .eq. 1) then
         call fbfield (1, lsearch, rsearch, bsearch,fail)
         if (posp(1) .gt. bsearch(2,1)) fail = 1 !rsearch
         if (posp(1) .lt. bsearch(1,1)) fail = 1
         if (posp(2) .gt. bsearch(2,2)) fail = 1
         if (posp(2) .lt. bsearch(1,2)) fail = 1
         if (posp(3) .gt. bsearch(2,3)) fail = 1
         if (posp(3) .lt. bsearch(1,3)) fail = 1
         if (fail .eq. 1) then
            call fbfield (0, lsearch, posp(1:3), bsearch,fail)
            fail = 0
         endif
         if (lsearch .eq. 0) then
            out = 1
            go to 4000
         endif
         dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
         corner(:) = bsearch(1,:)
      endif

c ---------------------------------------------------------------
c check 4th order against 3rd order for error estimation
c ---------------------------------------------------------------
      if (apt) then
        call interp_block (varb(:,:,:,:,lsearch),dv,corner,posp(1:3),bv)
        blen = sqrt(bv(1)**two+bv(2)**two+bv(3)**two)

        if (blen .ne. 0.0) then
           delta(:) = (k4(1:3)-len1*bv(:)/blen)/six
           err1=sqrt(delta(1)**2.0+delta(2)**2.0+delta(3)**2.0)
         
           if (err1 .gt. ep) then
              len1 = half*len1
           else
              len1 = two*len1
           endif
        else
           len1 = half*len1
        endif
      endif         
c-----------------------------------------
c          twist calculation
c-----------------------------------------
c pos(:,10) = curl(B).B/(4pi*B^2)
        if (sph) then
           call tw_sph (varb(:,:,:,:,lsearch),dv,corner,posp)
           pos(n+1,10) = posp(10)
        else
           call tw_cart (varb(:,:,:,:,lsearch),dv,corner,posp)
           pos(n+1,10) = posp(10)
        endif
c---------------------------------------------------------------
      else
c if outside the volume put the rest of the points on the boundary
 4000    continue
      pos(n+1,:) = pos1(:)
      endif


      end do
c find B at end point (for bends in fieldline_per.F)
c make sure in the volume for non-periodic directions
      if (xp) then
      else
         if (pos(n1,1) .lt. bbox(1,1)) then
            pos(n1,1) = bbox(1,1)
         endif
         if (pos(n1,1) .gt. bbox(2,1)) then
            pos(n1,1) = bbox(2,1)
         endif
      endif
      if (yp) then
      else
         if (pos(n1,2) .lt. bbox(1,2)) then
            pos(n1,2) = bbox(1,2)
         endif
         if (pos(n1,2) .gt. bbox(2,2)) then
            pos(n1,2) = bbox(2,2)
         endif
      endif
      if (zp) then
      else
         if (pos(n1,3) .lt. bbox(1,3)) then
            pos(n1,3) = bbox(1,3)
         endif
         if (pos(n1,3) .gt. bbox(2,3)) then
            pos(n1,3) = bbox(2,3)
         endif
      endif

c put in volume for periodic directions
      call poscheck (pos(n1,:),posp,bbox)
      ls = -1
      call fbfield (0, ls, posp(1:3), bsearch,fail)

      dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
      corner(:) = bsearch(1,:)
      call interp_block(varb(:,:,:,:,ls),dv,corner,posp(1:3),bv1)

c find the projected values of uvec and vvec at the end points
      blen = sqrt(bv1(1)**2.0+bv1(2)**2.0+bv1(3)**2.0)
      if (blen .eq. 0.0) then
      else
      bhat(:) = bv1(:)/blen
      endif

      tmp = pos(n1,4)*bhat(1)+pos(n1,5)*bhat(2)+pos(n1,6)*bhat(3)
      pos(n1,4:6) = pos(n1,4:6)-bhat(:)*tmp

      tmp = pos(n1,7)*bhat(1)+pos(n1,8)*bhat(2)+pos(n1,9)*bhat(3)
      pos(n1,7:9) = pos(n1,7:9)-bhat(:)*tmp
      
      return

      end subroutine rk4
