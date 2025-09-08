      subroutine interp_block (bvec,dv,corner,pos,bv)
c  
c     linear interpolation of B field within a block
c     Peter Wyper Jan 2015
c
c inputs
      real pos(3)
      real dv(3),corner(3)
      real bvec(9,9,9,3)
c internal
      real one
      real dx,dy,dz
      integer i,j,k
      real tmp1(3),tmp2(3),tmp3(3),tmp4(3)
c output
      real bv(3)

      one = real(1.0)

      i = int((pos(1)-corner(1))/dv(1))+1
      j = int((pos(2)-corner(2))/dv(2))+1
      k = int((pos(3)-corner(3))/dv(3))+1

      dx = mod((pos(1)-corner(1)),dv(1))/dv(1)
      dy = mod((pos(2)-corner(2)),dv(2))/dv(2)
      dz = mod((pos(3)-corner(3)),dv(3))/dv(3)

c      print *,'corner = ',corner
c      print *,'i,j,k = ',i,j,k
c      print *,'dx,dy,dz = ',dx,dy,dz


      if (i .eq. 9) then
c top corner
         if (j .eq. 9) then
            if (k .eq. 9) then
               bv(:) = bvec(i,j,k,:)
               return
            endif
         endif

c top xy vertex
         if (j .eq. 9) then
            bv(:) = dz*bvec(i,j,k+1,:)+(one-dz)*bvec(i,j,k,:)
            return
         endif
c top xz vertex
         if (k .eq. 9) then
            bv(:) = dy*bvec(i,j+1,k,:)+(one-dy)*bvec(i,j,k,:) 
            return
         endif
         
      endif

c top yz vertex
      if (j .eq. 9) then
         if (k .eq. 9) then
            bv(:) = dx*bvec(i+1,j,k,:)+(one-dx)*bvec(i,j,k,:)
            return
         endif
      endif

c top x face
      if (i .eq. 9) then
         tmp1(:) = dz*bvec(i,j,k+1,:)+(one-dz)*bvec(i,j,k,:)
         tmp2(:) = dz*bvec(i,j+1,k+1,:)+(one-dz)*bvec(i,j+1,k,:)

         bv(:) = dy*tmp2(:)+(one-dy)*tmp1(:)
         return
      endif

c top y face
      if (j .eq. 9) then
         tmp1(:) = dz*bvec(i,j,k+1,:)+(one-dz)*bvec(i,j,k,:)
         tmp2(:) = dz*bvec(i+1,j,k+1,:)+(one-dz)*bvec(i+1,j,k,:)

         bv(:) = dx*tmp2(:)+(one-dx)*tmp1(:)
         return
      endif

c top z face
      if (k .eq. 9) then
         tmp1(:) = dy*bvec(i,j+1,k,:)+(one-dy)*bvec(i,j,k,:)
         tmp2(:) = dy*bvec(i+1,j+1,k,:)+(one-dy)*bvec(i+1,j,k,:)

         bv(:) = dx*tmp2(:)+(one-dx)*tmp1(:)
         return
      endif

c main volume

      tmp1(:) = dz*bvec(i,j,k+1,:)+(one-dz)*bvec(i,j,k,:)
      tmp2(:) = dz*bvec(i+1,j,k+1,:)+(one-dz)*bvec(i+1,j,k,:)
      tmp3(:) = dx*tmp2(:)+(one-dx)*tmp1(:)

      tmp1(:) = dz*bvec(i,j+1,k+1,:)+(one-dz)*bvec(i,j+1,k,:)
      tmp2(:) = dz*bvec(i+1,j+1,k+1,:)+(one-dz)*bvec(i+1,j+1,k,:)
      tmp4(:) = dx*tmp2(:)+(one-dx)*tmp1(:)

      bv(:) = dy*tmp4(:)+(one-dy)*tmp3(:)
      return


      end subroutine interp_block

