      program fieldline_per
c
c     trace a grid of field lines
c     
c     outputs: 
c     the fieldlines
c     the start and end grid positions (the mapping) - forward and backward traced
c     the b field values at the start and end points
c     the Q map
c     the twist map
c
c     Peter Wyper Feb 2015
c     update for periodic boundaries - June 2020
c     major update - June 2023
c     removed tmap
c     added forward/back tracing outputs
c     added volumetric Q using Tassev et al. 17 method
c     added Tw calculation
c
      use bfield_common
c
      character*23  fbfname
      character*23  bbfname
      character*21  fmfname
      character*21  bmfname
      character*22  bfname
      character*17  qfname
      character*19  fbendnm
      character*19  bbendnm
      character*17  tname
      character*17  flname
      character*18  twname
c
      real, dimension(:), allocatable :: pos(:,:)
      real, dimension(:), allocatable :: gridfl(:,:,:)
      real, dimension(:), allocatable :: fmapping(:,:,:,:)
      real, dimension(:), allocatable :: bmapping(:,:,:,:)
      real, dimension(:), allocatable :: bbends(:,:,:,:)
      real, dimension(:), allocatable :: fbends(:,:,:,:)
      real, dimension(:), allocatable :: qmap(:,:)
      real, dimension(:), allocatable :: df1(:),df2(:),df3(:),df4(:)
      real, dimension(:), allocatable :: flines(:,:,:,:)
      real, dimension(:), allocatable :: flen(:,:)
      real, dimension(:), allocatable :: twist(:,:)
c
      real, dimension(:), allocatable :: varb(:,:,:,:,:)
c
      integer   jstep
      integer   lmin, lmax
      integer   n1,n2,n3,nn,nsnap
      integer   nn1,nn2,ns
      logical   sp
      real      bbox(2,3)
      real      rlo(3),rhi(3)
      integer   i,j,fail,kk
c
      real      eighth,eight,one,two,half,six,sixth
      real      len,len2,lmx
      real      bsearch(2,3),corn(3),dv(3),bv(3)
      integer   ls, switch
      real      tolb
      real      dx,tmp,pi1,tmp1
c
      logical   trnmap,savefls,saveflen,savetwist
      logical   oldqmap, apt
      real      y0,z0
      real      top,bot
      logical   xp,yp,zp
      real      posp(3),bv1(3)   
c
      real      uvf(3),uvb(3),vvf(3),vvb(3)
      real      bb,fb,b0
c
 1001  format (' rbfield reads bfield data at jstep = ',i7,/)
 1002  format (' rbfield reads nbtt = ',i12,/,
     .        '               nblf = ',i12,/,
     .        '               lmin = ',i12,/,
     .        '               lmax = ',i12,//,
     .        '               time = ',1p,1e12.4,/,
     .        '               xmin = ',1e12.4,/,
     .        '               xmax = ',1e12.4,/,
     .        '               ymin = ',1e12.4,/,
     .        '               ymax = ',1e12.4,/,
     .        '               zmin = ',1e12.4,/,
     .        '               zmax = ',1e12.4,/)

c
 1003    format (i4,' = nsteps',/,i4, ' = n1',/,i4,' = n2')
 1004    format ( 1e12.4,' = rlox',/
     .        ,1e12.4,' = rhix',/
     .        ,1e12.4,' = rloy',/
     .        ,1e12.4,' = rhiy',/
     .        ,1e12.4,' = rloz',/
     .        ,1e12.4,' = rhiz',/)
 1005    format (1e12.4,' = lstart',/,1e12.4,' = lmax',/,1e12.4,' = ep')
c
 2001     format ('PFLS/fieldline.000000',i1)
 2002     format ('PFLS/fieldline.00000',i2)
 2003     format ('PFLS/fieldline.0000',i3)
 2004     format ('PFLS/fieldline.000',i4)
 2005     format ('PFLS/fieldline.00',i5)
 2006     format ('PFLS/fieldline.0',i6)
 2007     format ('PFLS/fieldline.',i7)

 3001     format ('PFLS/fmapping.000000',i1)
 3002     format ('PFLS/fmapping.00000',i2)
 3003     format ('PFLS/fmapping.0000',i3)
 3004     format ('PFLS/fmapping.000',i4)
 3005     format ('PFLS/fmapping.00',i5)
 3006     format ('PFLS/fmapping.0',i6)
 3007     format ('PFLS/fmapping.',i7)
 3501     format ('PFLS/bmapping.000000',i1)
 3502     format ('PFLS/bmapping.00000',i2)
 3503     format ('PFLS/bmapping.0000',i3)
 3504     format ('PFLS/bmapping.000',i4)
 3505     format ('PFLS/bmapping.00',i5)
 3506     format ('PFLS/bmapping.0',i6)
 3507     format ('PFLS/bmapping.',i7)
 4001     format ('PFLS/qmap.000000',i1)
 4002     format ('PFLS/qmap.00000',i2)
 4003     format ('PFLS/qmap.0000',i3)
 4004     format ('PFLS/qmap.000',i4)
 4005     format ('PFLS/qmap.00',i5)
 4006     format ('PFLS/qmap.0',i6)
 4007     format ('PFLS/qmap.',i7)
 5001     format ('PFLS/fbends.000000',i1)
 5002     format ('PFLS/fbends.00000',i2)
 5003     format ('PFLS/fbends.0000',i3)
 5004     format ('PFLS/fbends.000',i4)
 5005     format ('PFLS/fbends.00',i5)
 5006     format ('PFLS/fbends.0',i6)
 5007     format ('PFLS/fbends.',i7)
 5501     format ('PFLS/bbends.000000',i1)
 5502     format ('PFLS/bbends.00000',i2)
 5503     format ('PFLS/bbends.0000',i3)
 5504     format ('PFLS/bbends.000',i4)
 5505     format ('PFLS/bbends.00',i5)
 5506     format ('PFLS/bbends.0',i6)
 5507     format ('PFLS/bbends.',i7)
 7001     format ('PFLS/flen.000000',i1)
 7002     format ('PFLS/flen.00000',i2)
 7003     format ('PFLS/flen.0000',i3)
 7004     format ('PFLS/flen.000',i4)
 7005     format ('PFLS/flen.00',i5)
 7006     format ('PFLS/flen.0',i6)
 7007     format ('PFLS/flen.',i7)
 8001     format ('PFLS/twist.000000',i1)
 8002     format ('PFLS/twist.00000',i2)
 8003     format ('PFLS/twist.0000',i3)
 8004     format ('PFLS/twist.000',i4)
 8005     format ('PFLS/twist.00',i5)
 8006     format ('PFLS/twist.0',i6)
 8007     format ('PFLS/twist.',i7)
c     
c * * *
c
c     allocate the main array for bfield data
         allocate (varb(mx,my,mz,lvars,mblks))
c         
c global variable to be used in fbfield.F 
        tolbfd = EPSILON(real(1.0)) 


c local variables
        eight = real (8.0)
        six = real (6.0)
        one = real (1.0)
        two = real (2.0)
        eighth = one/eight
        half = one/two
        sixth = one/six
        fail = 0
        pi1 = 3.14159265359
c
c     open files. 
c
        open (unit=15, file='fieldline.cnt', status='old')
c
c    read in data from fieldline.cnt
        read  (15,*) sp
        read  (15,*) xp,yp,zp
        read  (15,*) savefls,saveflen,savetwist
        read  (15,*) n1,n2,n3
        read  (15,*) rlo(1),rlo(2),rlo(3)
        read  (15,*) rhi(1),rhi(2),rhi(3)
        read  (15,*) apt
        read  (15,*) ns,len,lmx,ep
c        read  (15,*) top, bot

        read  (15,*) nsnap

        if (n1 .eq. 0) then
           if (savefls)then
           allocate (flines(n2+1,n3+1,2*ns,3))
           endif

           if (saveflen)then
           allocate (flen(n2+1,n3+1))
           endif
           if (savetwist)then
           allocate (twist(n2+1,n3+1))
           endif

           allocate (gridfl(n2+1,n3+1,3))
           nn1 = n2
           nn2 = n3
        endif

        if (n2 .eq. 0) then
           if (savefls) then
           allocate (flines(n1+1,n3+1,2*ns,3))
           endif
           if (saveflen) then
           allocate (flen(n1+1,n3+1))
           endif
           if (savetwist)then
           allocate (twist(n1+1,n2+1))
           endif
           allocate (gridfl(n1+1,n3+1,3))
           nn1 = n1
           nn2 = n3
        endif

        if (n3 .eq. 0) then
           if (savefls) then
           allocate (flines(n1+1,n2+1,2*ns,3))          
           endif
           if (saveflen) then
           allocate (flen(n1+1,n2+1))
           endif
           if (savetwist)then
           allocate (twist(n1+1,n2+1))
           endif
           allocate (gridfl(n1+1,n2+1,3))
           nn1 = n1
           nn2 = n2
        endif

c --------------------------------------------------------------------
c  set starting points on plane 

        if (n1 .eq. 0) then
           gridfl(:,:,1) = gridfl(:,:,1)*0.0 + rlo(1) !+ tol
           
           do i=1,n2+1
              gridfl(i,1,2) = (i-1)*(rhi(2)-rlo(2))/n2 +rlo(2)
           end do
           
           do j=1,n3+1
              gridfl(1,j,3) = (j-1)*(rhi(3)-rlo(3))/n3 +rlo(3)
              gridfl(:,j,2) = gridfl(:,1,2)
           end do

           do i=1,n2+1
              gridfl(i,:,3) = gridfl(1,:,3)
           end do
        endif
         

        if (n2 .eq. 0) then
           gridfl(:,:,2) = gridfl(:,:,2)*0.0 + rlo(2) !+ tol

           do i=1,n1+1
              gridfl(i,1,1) = (i-1)*(rhi(1)-rlo(1))/n1 +rlo(1)
           end do

           do j=1,n3+1
              gridfl(1,j,3) = (j-1)*(rhi(3)-rlo(3))/n3 +rlo(3)
              gridfl(:,j,1) = gridfl(:,1,1)
           end do

           do i=1,n1+1
              gridfl(i,:,3) = gridfl(1,:,3)
           end do
        endif


        if (n3 .eq. 0) then
           gridfl(:,:,3) = gridfl(:,:,3)*0.0 + rlo(3) ! + tol

           do i=1,n1+1
              gridfl(i,1,1) = (i-1)*(rhi(1)-rlo(1))/n1 +rlo(1)
           end do

           do j=1,n2+1
              gridfl(1,j,2) = (j-1)*(rhi(2)-rlo(2))/n2 +rlo(2)
              gridfl(:,j,1) = gridfl(:,1,1)
           end do

           do i=1,n1+1
              gridfl(i,:,2) = gridfl(1,:,2)
           end do
        endif
c --------------------------------------------------------------------
c        open (unit=19,file='PFLS/fieldline.prt', status='unknown')
c        write (19,1003) nstep, nn1+1, nn2+1
c        write (19,1005) len, lmx, ep
c        write (19,1004) rlo(1),rhi(1),rlo(2),
c     .                  rhi(2),rlo(3),rhi(3) 

c        open (unit=22,file='PFLS/jsteps.prt',status='new')

        allocate (pos(ns,10))
        allocate (fmapping(nn1+1,nn2+1,2,3))
        allocate (bmapping(nn1+1,nn2+1,2,3))
        allocate (bbends(nn1+1,nn2+1,2,3))
        allocate (fbends(nn1+1,nn2+1,2,3))
        allocate (qmap(nn1+1,nn2+1))
        allocate (df1(nn1+1))
        allocate (df2(nn2+1))

        do nn=1,nsnap

        read  (15,*) jstep
c        write (22,*) jstep

        print *,'jstep = ',jstep
c  read bfield data at this time step
        call rbfield (jstep,varb)

        if (sp) then
c -------------------------------------------------------------------------------
c spericals
c -------------------------------------------------------------------------------
c r=log(r) in restart files - rescale the start positions to log(r)
           gridfl(:,:,1) = alog(gridfl(:,:,1)) 
c theta in restart files saved in lattitude - convert to co-lattitude 
           bndbox(:,2,:) = pi1*0.5 - bndbox(:,2,:)
c start theta positions in lattitude and degrees - convert to co-lattitude & radians
           gridfl(:,:,2) = pi1*0.5 - gridfl(:,:,2)*pi1/180.0
c start phi positions in degrees - convert to radians
           gridfl(:,:,3) = gridfl(:,:,3)*pi1/180.0
c 
           fmapping = 0.0*fmapping
           fmapping(:,:,1,:) = gridfl(:,:,:)
           bmapping = 0.0*bmapping
           bmapping(:,:,1,:) = gridfl(:,:,:)
c
        lmin = +1000
        lmax = -1000
        bbox(1,1) = +1000.
        bbox(2,1) = -1000.
        bbox(1,2) = +1000.
        bbox(2,2) = -1000.
        bbox(1,3) = +1000.
        bbox(2,3) = -1000.
c               
        do lbt = 1, ntotl
          lmin = min (lmin, lrefine(lbt))
          lmax = max (lmax, lrefine(lbt))
          bbox(1,1) = min (bbox(1,1), bndbox(1,1,lbt))
          bbox(2,1) = max (bbox(2,1), bndbox(2,1,lbt))
          bbox(1,2) = min (bbox(1,2), bndbox(1,2,lbt))
          bbox(2,2) = max (bbox(2,2), bndbox(2,2,lbt))
          bbox(1,3) = min (bbox(1,3), bndbox(1,3,lbt))
          bbox(2,3) = max (bbox(2,3), bndbox(2,3,lbt))
       end do


c loop over each starting point and find field line
      do i=1,nn1+1
         print *,'i = ',i
        do j=1,nn2+1
c if the footpoint is outside the volume -  put on the edge
           if (gridfl(i,j,1) .lt. bbox(1,1)) then
              gridfl(i,j,1) = bbox(1,1)
           endif
           if (gridfl(i,j,1) .gt. bbox(2,1)) then
              gridfl(i,j,1) = bbox(2,1)
           endif
c
           if (gridfl(i,j,2) .lt. bbox(1,2)) then
              gridfl(i,j,2) = bbox(1,2)
           endif
           if (gridfl(i,j,2) .gt. bbox(2,2)) then
              gridfl(i,j,2) = bbox(2,2)
           endif
c
           if (gridfl(i,j,3) .lt. bbox(1,3)) then
              gridfl(i,j,3) = bbox(1,3)
           endif
           if (gridfl(i,j,3) .gt. bbox(2,3)) then
              gridfl(i,j,3) = bbox(2,3)
           endif
c find b at starting point
           ls = -1
           call fbfield (0, ls, gridfl(i,j,:), bsearch,fail)
           dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
           corn(:) = bsearch(1,:)
         call interp_block (varb(:,:,:,:,ls),dv,corn,gridfl(i,j,:),bv)
         fbends(i,j,1,:)=bv(:)
         bbends(i,j,1,:)=bv(:)

c trace in the direction of the bfield  
         len2 = len

c--------------------------------------------
c trace field line forwards (in log(r) space)
c--------------------------------------------
        pos(1,1:3) = gridfl(i,j,:)
        call rk4(pos(:,:),ns,len2,lmx,ep,bbox,varb,xp,yp,zp,bv1,sp,apt)
        fmapping(i,j,2,:)=pos(ns,1:3)
c find b at end point
        fbends(i,j,2,:)=bv1(:)

c store field line positions if asked for
        if (savefls) then
           flines(i,j,ns+1:2*ns,:) = pos(:,1:3)
        endif

        if (saveflen) then
c convert r component back
           pos(:,1) = exp(pos(:,1))
c theta stays in co-lattitude and radians
c phi stays in radians
c dl = sqrt(dr^2+r^2*dtheta^2+(rsin(theta))^2*dphi^2)
           tmp1 = 0.0
           do kk=1,ns-1
          tmp = (pos(kk,1)*sin(pos(kk,2)))**2*(pos(kk+1,3)-pos(kk,3))**2
            tmp = tmp + (pos(kk+1,1)-pos(kk,1))**2
            tmp = tmp + pos(kk,1)**2*(pos(kk+1,2)-pos(kk,2))**2
            tmp1 = tmp1 + sqrt(tmp)
           end do
           flen(i,j) = tmp1
           pos(:,1) = alog(pos(:,1))
        endif

        if (savetwist) then
c convert r component back                            
           pos(:,1) = exp(pos(:,1))
c theta stays in co-lattitude and radians             
c phi stays in radians
c dl = sqrt(dr^2+r^2*dtheta^2+(rsin(theta))^2*dphi^2) 
           tmp1 = 0.0
           do kk=1,ns-1
          tmp = (pos(kk,1)*sin(pos(kk,2)))**2*(pos(kk+1,3)-pos(kk,3))**2
            tmp = tmp + (pos(kk+1,1)-pos(kk,1))**2
            tmp = tmp + pos(kk,1)**2*(pos(kk+1,2)-pos(kk,2))**2
            tmp1 = tmp1 + sqrt(tmp)*(pos(kk,10)+pos(kk+1,10))*half
           end do
           twist(i,j) = tmp1
           pos(:,1) = alog(pos(:,1))
        endif

c store uvec and vvec at end point
        uvf(:) = pos(ns,4:6)
        vvf(:) = pos(ns,7:9)

c--------------------------------------------   
c trace field line backwards (in log(r) space)   
c--------------------------------------------   
        pos(1,1:3) = gridfl(i,j,:)
        call rk4(pos(:,:),ns,-len2,lmx,ep,bbox,varb,xp,yp,zp,bv1,sp,apt)
        bmapping(i,j,2,:)=pos(ns,1:3)
c find b at end point 
        bbends(i,j,2,:)=bv1(:)

c store field line positions if asked for 
        if (savefls) then
           do kk=1,ns
              flines(i,j,ns-kk+1,:) = pos(kk,1:3)
           end do
        endif

        if (saveflen) then
c convert r component back  
           pos(:,1) = exp(pos(:,1))
c           print *,'(A1) = ',pos(ns,3)
c theta stays in co-lattitude and radians  
c phi stays in radians    
c dl = sqrt(dr^2+r^2*dtheta^2+(rsin(theta))^2*dphi^2) 
           tmp1 = 0.0
           do kk=1,ns-1
          tmp = (pos(kk,1)*sin(pos(kk,2)))**2*(pos(kk+1,3)-pos(kk,3))**2
            tmp = tmp + (pos(kk+1,1)-pos(kk,1))**2
            tmp = tmp + pos(kk,1)**2*(pos(kk+1,2)-pos(kk,2))**2
            tmp1 = tmp1 + sqrt(tmp)
           end do
           flen(i,j) = flen(i,j) + tmp1
           pos(:,1) = alog(pos(:,1))
        endif

        if (savetwist) then
c convert r component back         
           pos(:,1) = exp(pos(:,1))
c theta stays in co-lattitude and radians      
c phi stays in radians  
c dl = sqrt(dr^2+r^2*dtheta^2+(rsin(theta))^2*dphi^2)
           tmp1 = 0.0
           do kk=1,ns-1
          tmp = (pos(kk,1)*sin(pos(kk,2)))**2*(pos(kk+1,3)-pos(kk,3))**2
            tmp = tmp + (pos(kk+1,1)-pos(kk,1))**2
            tmp = tmp + pos(kk,1)**2*(pos(kk+1,2)-pos(kk,2))**2
            tmp1 = tmp1 + sqrt(tmp)*(pos(kk,10)+pos(kk+1,10))*half
           end do
           twist(i,j) = twist(i,j) + tmp1
           pos(:,1) = alog(pos(:,1))
        endif

c store uvec and vvec at end point
        uvb(:) = pos(ns,4:6)
        vvb(:) = pos(ns,7:9)

c equation 11 from Tassev et al. 17 
c Q = Bf*Bb/B0^2*(af^2*bb^2+ab^2*bf^2-2*(ab dot ab)(af dot bf)                                                                                               
        tmp = uvf(1)**2+uvf(2)**2+uvf(3)**2
        tmp = tmp*(vvb(1)**2+vvb(2)**2+vvb(3)**2)

        tmp1 = uvb(1)**2+uvb(2)**2+uvb(3)**2
        tmp1 = tmp1*(vvf(1)**2+vvf(2)**2+vvf(3)**2)

        tmp = tmp + tmp1

        tmp1 = uvb(1)*vvb(1)+uvb(2)*vvb(2)+uvb(3)*vvb(3)
        tmp1 = tmp1*(uvf(1)*vvf(1)+uvf(2)*vvf(2)+uvf(3)*vvf(3))

        tmp = tmp - 2.0*tmp1

       bb=sqrt(bbends(i,j,2,1)**2+bbends(i,j,2,2)**2+bbends(i,j,2,3)**2)
       fb=sqrt(fbends(i,j,2,1)**2+fbends(i,j,2,2)**2+fbends(i,j,2,3)**2)

c B0^2
       b0=fbends(i,j,1,1)**2+fbends(i,j,1,2)**2+fbends(i,j,1,3)**2                                                                                                                                                                                                                                      
       if (b0 .eq. 0.0) then
       qmap(i,j) = 2.0
       else
       qmap(i,j) = bb*fb*tmp/b0
       endif

         end do
       end do

c---------------------------
c
c rescale the radial direction back again 
           fmapping(:,:,:,1) = exp(fmapping(:,:,:,1)) !all the points along the line now
           bmapping(:,:,:,1) = exp(bmapping(:,:,:,1)) 
           gridfl(:,:,1) = exp(gridfl(:,:,1))
c convert theta direction back to lattitude
           bndbox(:,2,:) = pi1*0.5 - bndbox(:,2,:)
c convert theta fieldline positions back to lattitude and degrees
           fmapping(:,:,:,2) = (pi1*0.5-fmapping(:,:,:,2))*180./pi1
           bmapping(:,:,:,2) = (pi1*0.5-bmapping(:,:,:,2))*180./pi1
           gridfl(:,:,2) = (pi1*0.5-gridfl(:,:,2))*180./pi1
c convert phi fieldline positions to degrees
           fmapping(:,:,:,3) = fmapping(:,:,:,3)*180./pi1
           bmapping(:,:,:,3) = bmapping(:,:,:,3)*180./pi1
           gridfl(:,:,3) = gridfl(:,:,3)*180./pi1
c
c rescale and convert file line positions if asked for
           if (savefls) then
              flines(:,:,:,1)=exp(flines(:,:,:,1))
              flines(:,:,:,2)=(pi1*0.5-flines(:,:,:,2))*180./pi1
              flines(:,:,:,3)=flines(:,:,:,3)*180./pi1
           endif
c
      else
c -------------------------------------------------------------------------------
c cartesian
c -------------------------------------------------------------------------------
c
        fmapping(:,:,1,:) = gridfl(:,:,:)
        bmapping(:,:,1,:) = gridfl(:,:,:)
c
        lmin = +1000
        lmax = -1000
        bbox(1,1) = +1000.
        bbox(2,1) = -1000.
        bbox(1,2) = +1000.
        bbox(2,2) = -1000.
        bbox(1,3) = +1000.
        bbox(2,3) = -1000.
c            
        do lbt = 1, ntotl
          lmin = min (lmin, lrefine(lbt))
          lmax = max (lmax, lrefine(lbt))
          bbox(1,1) = min (bbox(1,1), bndbox(1,1,lbt))
          bbox(2,1) = max (bbox(2,1), bndbox(2,1,lbt))
          bbox(1,2) = min (bbox(1,2), bndbox(1,2,lbt))
          bbox(2,2) = max (bbox(2,2), bndbox(2,2,lbt))
          bbox(1,3) = min (bbox(1,3), bndbox(1,3,lbt))
          bbox(2,3) = max (bbox(2,3), bndbox(2,3,lbt))
       end do

c loop over each starting point and find field line                                                                                                                                     
      do i=1,nn1+1
         print *,'i = ',i
        do j=1,nn2+1
c find b at starting point                                                                                                                                                              
           ls = -1
           call fbfield (0, ls, gridfl(i,j,:), bsearch,fail)
           dv(:) = (bsearch(2,:)-bsearch(1,:))*eighth
           corn(:) = bsearch(1,:)
        call interp_block (varb(:,:,:,:,ls),dv,corn,gridfl(i,j,:),bv)
        fbends(i,j,1,:)=bv(:)
        bbends(i,j,1,:)=bv(:)

c trace in bfield direction 
        len2  = len
c-------------------------
c trace field line forward
c-------------------------
        pos(1,1:3) = gridfl(i,j,:)
        call rk4(pos(:,:),ns,len2,lmx,ep,bbox,varb,xp,yp,zp,bv1,sp,apt)
        fmapping(i,j,2,:)=pos(ns,1:3)
c save b at end point        
         fbends(i,j,2,:)=bv1(:)

        if (savefls) then
           flines(i,j,ns+1:2*ns,:) = pos(:,1:3)
        endif

        if (saveflen) then
           tmp1 = 0.0
           do kk=1,ns-1
              tmp = (pos(kk+1,1)-pos(kk,1))**2
              tmp = tmp + (pos(kk+1,2)-pos(kk,2))**2
              tmp = tmp + (pos(kk+1,3)-pos(kk,3))**2
              tmp1 = tmp1 + sqrt(tmp)
           end do
           flen(i,j) = tmp1
        endif

        if (savetwist) then
           tmp1 = 0.0
           do kk=1,ns-1
              tmp = (pos(kk+1,1)-pos(kk,1))**2
              tmp = tmp + (pos(kk+1,2)-pos(kk,2))**2
              tmp = tmp + (pos(kk+1,3)-pos(kk,3))**2
              tmp1 = tmp1 + sqrt(tmp)*(pos(kk,10)+pos(kk+1,10))*half
           end do
           twist(i,j) = tmp1
        endif

c store uvec and vvec at end point
        uvf(:) = pos(ns,4:6)
        vvf(:) = pos(ns,7:9)

c--------------------------- 
c trace field line backwards 
c--------------------------- 
        pos(1,1:3) = gridfl(i,j,:)
        call rk4(pos(:,:),ns,-len2,lmx,ep,bbox,varb,xp,yp,zp,bv1,sp,apt)
        bmapping(i,j,2,:)=pos(ns,1:3)
c save b at end point              
         bbends(i,j,2,:)=bv1(:)

        if (savefls) then
           do kk=1,ns
              flines(i,j,ns-kk+1,:) = pos(kk,1:3)
           end do
        endif

        if (saveflen) then
           tmp1 = 0.0
           do kk=1,ns-1
              tmp = (pos(kk+1,1)-pos(kk,1))**2
              tmp = tmp + (pos(kk+1,2)-pos(kk,2))**2
              tmp = tmp + (pos(kk+1,3)-pos(kk,3))**2
              tmp1 = tmp1 + sqrt(tmp)
           end do
           flen(i,j) = flen(i,j) + tmp1
        endif

        if (savetwist) then
           tmp1 = 0.0
           do kk=1,ns-1
              tmp = (pos(kk+1,1)-pos(kk,1))**2
              tmp = tmp + (pos(kk+1,2)-pos(kk,2))**2
              tmp = tmp + (pos(kk+1,3)-pos(kk,3))**2
              tmp1 = tmp1 + sqrt(tmp)*(pos(kk,10)+pos(kk+1,10))*half
           end do
           twist(i,j) = twist(i,j) + tmp1
        endif

c store uvec and vvec at end point
        uvb(:) = pos(ns,4:6)
        vvb(:) = pos(ns,7:9)

c equation 11 from Tassev et al. 17
c Q = Bf*Bb/B0^2*(af^2*bb^2+ab^2*bf^2-2*(Ab dot Bb)(af dot bf)
        tmp = uvf(1)**2+uvf(2)**2+uvf(3)**2
        tmp = tmp*(vvb(1)**2+vvb(2)**2+vvb(3)**2)

        tmp1 = uvb(1)**2+uvb(2)**2+uvb(3)**2
        tmp1 = tmp1*(vvf(1)**2+vvf(2)**2+vvf(3)**2)

        tmp = tmp + tmp1
        
        tmp1 = uvb(1)*vvb(1)+uvb(2)*vvb(2)+uvb(3)*vvb(3)
        tmp1 = tmp1*(uvf(1)*vvf(1)+uvf(2)*vvf(2)+uvf(3)*vvf(3))

        tmp = tmp - 2.0*tmp1
 
       bb=sqrt(bbends(i,j,2,1)**2+bbends(i,j,2,2)**2+bbends(i,j,2,3)**2)
       fb=sqrt(fbends(i,j,2,1)**2+fbends(i,j,2,2)**2+fbends(i,j,2,3)**2) 
       b0=sqrt(fbends(i,j,1,1)**2+fbends(i,j,1,2)**2+fbends(i,j,1,3)**2) 
       if (b0 .eq. 0.0) then
       qmap(i,j) = 2.0
       else
       qmap(i,j) = bb*fb*tmp/(b0**2)
       endif

        end do
      end do

      endif
c ----------------------------------------------------------------------- 
c      close (unit=19)
c      close (unit=22)
c -----------------------------------------------------------------------

c -----------------------------------------------------------------------
      oldqmap = .false.
c old version of qmap calculation
      if (oldqmap) then
c calculate Q map - may need to trim if field lines hit side walls
      if (sp) then
c -----------------------------------------------------------------------
c sphericals - eqns 22 & 23 Titov 2007
c note - requires sin(theta) non-zero (i.e. no north or south pole)
c -----------------------------------------------------------------------
      qmap = 0.0*qmap
      if (n1 .eq. 0) then
        allocate (df3(nn1+1))
        allocate (df4(nn2+1))

c convert theta to co-lattitiude & radians
         fmapping(:,:,:,2) = pi1*0.5 - fmapping(:,:,:,2)*pi1/180.0
c convert phi to radians
         fmapping(:,:,:,3) = fmapping(:,:,:,3)*pi1/180.0

         dx = fmapping(2,1,1,2)-fmapping(1,1,1,2) !dtheta

         do j=1,nn2+1
c dTHETA/dtheta  
            call diff4 (fmapping(:,j,2,2),dx,nn1+1,df1(:))
c dPHI/dtheta
            call diff4 (fmapping(:,j,2,3),dx,nn1+1,df3(:))
            
            do i=1,nn1+1
               qmap(i,j) = (sin(fmapping(i,j,2,2))*df3(i) )**two !(sin(THETA)*dPHI/dtheta)^2
               qmap(i,j) = qmap(i,j) + df1(i)**two ! + (dTHETA/dtheta)^2
            end do
         end do

         dx = fmapping(1,2,1,3)-fmapping(1,1,1,3) !dphi
         do i=1,nn1+1
c dTHETA/dphi 
            call diff4 (fmapping(i,:,2,2),dx,nn2+1,df2(:))
c dPHI/dphi 
            call diff4 (fmapping(i,:,2,3),dx,nn2+1,df4(:))

            do j=1,nn2+1 
              qmap(i,j)=qmap(i,j)+ (df2(j)/sin(fmapping(i,j,1,2)) )**two ! + [(1/sin(theta)*dTHETA/dphi)]^2

              tmp = sin(fmapping(i,j,2,2))*df4(j)/sin(fmapping(i,j,1,2)) 
               qmap(i,j)=qmap(i,j)+ tmp**two ! + [(sin(THETA)/sin(theta))*dPHI/dphi]^2
            end do
         end do

c         qmap(:,:) = qmap(:,:)*(mapping(:,:,2,1)/mapping(:,:,1,1))**two ! *(R*/R0)^2
c         qmap(:,:) = qmap(:,:)*abs(bends(:,:,2,1)/bends(:,:,1,1)) ! *(Br*/Br)

         do i=1,nn1+1
            do j=1,nn2+1

          if (abs(fmapping(i,j,1,1)) .lt. 1e-6) then
          qmap(i,j) = qmap(i,j)
          else
        qmap(i,j) = qmap(i,j)*(fmapping(i,j,2,1)/fmapping(i,j,1,1))**two ! *(R*/R0)^2                                                                                 
          end if

          if (abs(fbends(i,j,1,1)) .lt. 1e-6) then
          qmap(i,j) = 2.0
          else
          qmap(i,j) = qmap(i,j)*abs(fbends(i,j,2,1)/fbends(i,j,1,1)) ! *(Br*/Br)                                                                                        
          end if

            end do
         end do


c convert theta fieldline positions back to lattitude and degrees
           fmapping(:,:,:,2) = (pi1*0.5-fmapping(:,:,:,2))*180./pi1
c convert phi fieldline positions to degrees
           fmapping(:,:,:,3) = fmapping(:,:,:,3)*180./pi1

         deallocate (df3,df4)

      endif

      else
c -----------------------------------------------------------------------
c cartesian
c -----------------------------------------------------------------------
      qmap = 0.0*qmap
      if (n1 .eq. 0) then

         dx = fmapping(2,1,1,2)-fmapping(1,1,1,2)

c dY/dy         
c dZ/dy
         do j=1,nn2+1
         call diff4 (fmapping(:,j,2,2),dx,nn1+1,df1(:))
         qmap(:,j) =  df1(:)**two
         call diff4 (fmapping(:,j,2,3),dx,nn1+1,df1(:))
         qmap(:,j) =  qmap(:,j) + df1(:)**two
         end do

         dx = fmapping(1,2,1,3)-fmapping(1,1,1,3)
c dY/dz
c dZ/dz
         do i=1,nn1+1
         call diff4 (fmapping(i,:,2,2),dx,nn2+1,df2(:))
         qmap(i,:)=qmap(i,:)+ df2(:)**two
         call diff4 (fmapping(i,:,2,3),dx,nn2+1,df2(:))
         qmap(i,:)=qmap(i,:)+ df2(:)**two
         end do

         qmap(:,:) = qmap(:,:)*abs(fbends(:,:,2,1)/fbends(:,:,1,1))
      endif

      endif

      endif
c ----------------------------------------------------------------------
c   save data
c ----------------------------------------------------------------------    

c bfield data at end points
      if (jstep .le. 9) then
         write (fbendnm,5001) jstep
      else if (jstep .le. 99) then
         write (fbendnm,5002) jstep
      else if (jstep .le. 999) then
         write (fbendnm,5003) jstep
      else if (jstep .le. 9999) then
         write (fbendnm,5004) jstep
      else if (jstep .le. 99999) then
         write (fbendnm,5005) jstep
      else if (jstep .le. 999999) then
         write (fbendnm,5006) jstep
      else if (jstep .le. 9999999) then
         write (fbendnm,5007) jstep
      end if

      open (unit=21, file=fbendnm, status='new'
     .     ,form='unformatted')
      write (21) fbends                                                                                                                                                                                                                                                                                       
      close (unit=21)

      if (jstep .le. 9) then
         write (bbendnm,5501) jstep
      else if (jstep .le. 99) then
         write (bbendnm,5502) jstep
      else if (jstep .le. 999) then
         write (bbendnm,5503) jstep
      else if (jstep .le. 9999) then
         write (bbendnm,5504) jstep
      else if (jstep .le. 99999) then
         write (bbendnm,5505) jstep
      else if (jstep .le. 999999) then
         write (bbendnm,5506) jstep
      else if (jstep .le. 9999999) then
         write (bbendnm,5507) jstep
      end if

      open (unit=21, file=bbendnm, status='new'                                                                                                                                                                                                                                                               
     .     ,form='unformatted')                                                                                                                                                                                                                                                                               
      write (21) bbends                                                                                                                                                                                                                                                                                       
      close (unit=21)  

c field lines positions
      if (savefls) then
        if (jstep .le. 9) then
          write (bfname,2001) jstep
        else if (jstep .le. 99) then
          write (bfname,2002) jstep
        else if (jstep .le. 999) then
          write (bfname,2003) jstep
        else if (jstep .le. 9999) then
          write (bfname,2004) jstep
        else if (jstep .le. 99999) then
          write (bfname,2005) jstep
        else if (jstep .le. 999999) then
          write (bfname,2006) jstep
        else if (jstep .le. 9999999) then
          write (bfname,2007) jstep
        end if

        open (unit=23, file=bfname, status='new' 
     .                  ,form='unformatted')
        write (23) flines        
        close (unit=23)
      endif

c field line lengths                                                                                                                                                                                  
      if (saveflen) then
        if (jstep .le. 9) then
          write (flname,7001) jstep
        else if (jstep .le. 99) then
          write (flname,7002) jstep
        else if (jstep .le. 999) then
          write (flname,7003) jstep
        else if (jstep .le. 9999) then
          write (flname,7004) jstep
        else if (jstep .le. 99999) then
          write (flname,7005) jstep
        else if (jstep .le. 999999) then
          write (flname,7006) jstep
        else if (jstep .le. 9999999) then
          write (flname,7007) jstep
        end if

        open (unit=25, file=flname, status='new'
     .                  ,form='unformatted')
        write (25) flen
        close (unit=25)
      endif

c twist
      if (savetwist) then
        if (jstep .le. 9) then
          write (twname,8001) jstep
        else if (jstep .le. 99) then
          write (twname,8002) jstep
        else if (jstep .le. 999) then
          write (twname,8003) jstep
        else if (jstep .le. 9999) then
          write (twname,8004) jstep
        else if (jstep .le. 99999) then
          write (twname,8005) jstep
        else if (jstep .le. 999999) then
          write (twname,8006) jstep
        else if (jstep .le. 9999999) then
          write (twname,8007) jstep
        end if

        open (unit=26, file=twname, status='new'
     .                  ,form='unformatted')
        write (26) twist
        close (unit=26)
      endif

c mapping
        if (jstep .le. 9) then
           write (fmfname,3001) jstep
        else if (jstep .le. 99) then
           write (fmfname,3002) jstep
        else if (jstep .le. 999) then
           write (fmfname,3003) jstep
        else if (jstep .le. 9999) then
           write (fmfname,3004) jstep
        else if (jstep .le. 99999) then
           write (fmfname,3005) jstep
        else if (jstep .le. 999999) then
           write (fmfname,3006) jstep
        else if (jstep .le. 9999999) then
           write (fmfname,3007) jstep
        end if

        open (unit=19,file=fmfname,status='new'                                                                                                                                                                                                                                                               
     .                  ,form='unformatted')                                                                                                                                                                                                                                                                  
        write (19) fmapping                                                                                                                                                                                                                                                                                   
        close (19) 


        if (jstep .le. 9) then
           write (bmfname,3501) jstep
        else if (jstep .le. 99) then
           write (bmfname,3502) jstep
        else if (jstep .le. 999) then
           write (bmfname,3503) jstep
        else if (jstep .le. 9999) then
           write (bmfname,3504) jstep
        else if (jstep .le. 99999) then
           write (bmfname,3505) jstep
        else if (jstep .le. 999999) then
           write (bmfname,3506) jstep
        else if (jstep .le. 9999999) then
           write (bmfname,3507) jstep
        end if

        open (unit=19,file=bmfname,status='new'                                                                                                                                                                                                                                                               
     .                  ,form='unformatted')                                                                                                                                                                                                                                                                  
        write (19) bmapping                                                                                                                                                                                                                                                                                   
        close (19)    


        if (oldqmap) then
        if (n1 .eq. 0) then
c Qmap - only for x or r = photosphere  
           if (jstep .le. 9) then
              write (qfname,4001) jstep
           else if (jstep .le. 99) then
              write (qfname,4002) jstep
           else if (jstep .le. 999) then
              write (qfname,4003) jstep
           else if (jstep .le. 9999) then
              write (qfname,4004) jstep
           else if (jstep .le. 99999) then
              write (qfname,4005) jstep
           else if (jstep .le. 999999) then
              write (qfname,4006) jstep
           else if (jstep .le. 9999999) then
              write (qfname,4007) jstep
           end if
           open (unit=20,file=qfname,status='new'
     .             ,form='unformatted')
           write (20) qmap
           close (20)
        endif
        else
c Qmap - in the volume also
           if (jstep .le. 9) then
              write (qfname,4001) jstep
           else if (jstep .le. 99) then
              write (qfname,4002) jstep
           else if (jstep .le. 999) then
              write (qfname,4003) jstep
           else if (jstep .le. 9999) then
              write (qfname,4004) jstep
           else if (jstep .le. 99999) then
              write (qfname,4005) jstep
           else if (jstep .le. 999999) then
              write (qfname,4006) jstep
           else if (jstep .le. 9999999) then
              write (qfname,4007) jstep
           end if
           open (unit=20,file=qfname,status='new'
     .             ,form='unformatted')
           write (20) qmap
           close (20)
        endif

      end do
      
      close (unit=15)
c
        deallocate (pos,gridfl,fmapping,bmapping,fbends,bbends)
        deallocate (qmap,df1,df2) 
c       deallocate (pos,gridfl,mapping,bends,qmap,df1,df2)
c
c *************************** 
c deallocate the bfield array 
      deallocate (varb)
c *************************** 
        if (savefls) then
           deallocate (flines)
        endif
        if (saveflen) then
           deallocate (flen)
        endif
        if (savetwist) then
           deallocate (twist)
        endif

c
        stop
c
      end program fieldline_per

      subroutine diff4 (f,dx,n,dfdx)
c central difference derivative for Q map
      integer n
      real f(n),dfdx(n),dx
      real diff, eght
      integer i

      diff = real(1.0)/real(12.0)
      eght = real(8.0)

      do i=3,n-1
         dfdx(i)=(-f(i-2)+f(i-1)*eght-f(i+1)*eght+f(i+2))*diff
      enddo

      dfdx(1:2)   = real(0.0)
      dfdx(n:n+1) = real(0.0)

      dfdx = dfdx/dx

      end subroutine diff4











