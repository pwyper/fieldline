      subroutine tw_cart (bvec,dv,corner,pos)
c saves curl(B).B/(4Pi*B^2) into pos(10)
c integrated along field line later to give Tw (eqn. 7, Liu et al. 2016)
c Peter Wyper June 2023                  
c inputs 
      real pos(10)
      real dv(3),corner(3)
      real bvec(9,9,9,3)

c internal 
      real bv(2,2,2,3)
      real a(3),b(3),c(3),d(3),e(3),f(3),g(3),h(3)
      real dx,dy,dz
      integer i,j,k
      real dbdx(3),dbdy(3),dbdz(3)
      real bpos(3),curlb(3),blen
      real pi1
c
      pi1 = 3.14159265359
c
c position in the block
      i = int((pos(1)-corner(1))/dv(1))+1
      j = int((pos(2)-corner(2))/dv(2))+1
      k = int((pos(3)-corner(3))/dv(3))+1
c position as a fraction of the cell
      dx = mod((pos(1)-corner(1)),dv(1))/dv(1)
      dy = mod((pos(2)-corner(2)),dv(2))/dv(2)
      dz = mod((pos(3)-corner(3)),dv(3))/dv(3)
c if on side boundary (numerical error in pos takes it over) adjust the index
      if (i .eq. 9) then
         i = 8
         dx = real(1.0)
      endif
      if (j .eq. 9) then
         j = 8
         dy = real(1.0)
      endif
      if (k .eq. 9) then
         k = 8
         dz = real(1.0)
      endif

c cell itself
      bv(:,:,:,:) = bvec(i:i+1,j:j+1,k:k+1,:)

c construct the interpolation coefficients 
c eqn. 1 in Haynes and Parnell 07 but in vector form (for 3 field components)
      a(:) = bv(2,1,1,:)
      b(:) = bv(2,1,1,:)-bv(1,1,1,:) 
      c(:) = bv(1,2,1,:)-bv(1,1,1,:)
      d(:) = bv(2,2,1,:)-bv(2,1,1,:)-bv(1,2,1,:)+bv(1,1,1,:)
      e(:) = bv(1,1,2,:)-bv(1,1,1,:)
      f(:) = bv(2,1,2,:)-bv(2,1,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      g(:) = bv(1,2,2,:)-bv(1,2,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      h(:) = bv(2,2,2,:)-bv(2,2,1,:)-bv(2,1,2,:)-bv(1,2,2,:)
      h(:) = h(:) + bv(2,1,1,:)+bv(1,2,1,:)+bv(1,1,2,:)-bv(1,1,1,:)

c bfield at the position
      bpos(:) = a(:)+b(:)*dx+c(:)*dy+d(:)*dx*dy+e(:)*dz+f(:)*dx*dz
      bpos(:) = bpos(:)+g(:)*dy*dz+h(:)*dx*dy*dz
c Blength
      blen = sqrt(bpos(1)**2+bpos(2)**2+bpos(3)**2)

c analytical expressions for dBx/dx, dBx/dy etc found by differentiating the interpolation
c dbdx = (dbxdx,dbydx,dbzdx) and so on
      dbdx(:) = (b(:) + d(:)*dy + f(:)*dz + h(:)*dy*dz)/dv(1)
      dbdy(:) = (c(:) + d(:)*dx + g(:)*dz + h(:)*dx*dz)/dv(2)
      dbdz(:) = (e(:) + f(:)*dx + g(:)*dy + h(:)*dx*dy)/dv(3)

c dv factor since dx etc are scaled from 0 to 1 in each direction
c so e.g. x' = 0 -> 1, with x' = (x-x(x'=0))/(x(1)-x(0))
c then db/dx = db/dx'*dx'/dx = db/dx'/(x(1)-x(0))

c find curlB
      curlb(1) = dbdy(3)-dbdz(2)
      curlb(2) = dbdz(1)-dbdx(3)
      curlb(3) = dbdx(2)-dbdy(1)

      pos(10) = curlb(1)*bpos(1)+curlb(2)*bpos(2)+curlb(3)*bpos(3)
      pos(10) = pos(10)/(real(4.0)*pi1*blen**2)

      end subroutine tw_cart


      subroutine tw_sph (bvec,dv,corner,pos)
c saves curl(B).B/(4Pi*B^2) into pos(10) 
c integrated along field line later to give Tw (eqn. 7, Liu et al. 2016) 
c Peter Wyper June 2023                                                        
c inputs                                                                       
      real pos(10)
      real dv(3),corner(3)
      real bvec(9,9,9,3)
c internal                                                                     
      real bv(2,2,2,3)
      real a(3),b(3),c(3),d(3),e(3),f(3),g(3),h(3)
      real dr,dt,dp
      integer i,j,k
      real dbdr(3),dbdt(3),dbdp(3)
      real r, th, bpos(3), curlb(3)
      real pi1
c                        
      pi1 = 3.14159265359

c position in the block                                                        
      i = int((pos(1)-corner(1))/dv(1))+1
      j = int((pos(2)-corner(2))/dv(2))+1
      k = int((pos(3)-corner(3))/dv(3))+1
c position as a fraction of the cell                                           
      dr = mod((pos(1)-corner(1)),dv(1))/dv(1)
      dt = mod((pos(2)-corner(2)),dv(2))/dv(2)
      dp = mod((pos(3)-corner(3)),dv(3))/dv(3)
c if on side boundary (numerical error in pos takes it over) adjust the index  
      if (i .eq. 9) then
         i = 8
         dr = real(1.0)
      endif
      if (j .eq. 9) then
         j = 8
         dt = real(1.0)
      endif
      if (k .eq. 9) then
         k = 8
         dp = real(1.0)
      endif

c cell itself                                                                              
      bv(:,:,:,:) = bvec(i:i+1,j:j+1,k:k+1,:)

c construct the interpolation coefficients                                                 
c eqn. 1 in Haynes and Parnell 07 but in vector form (for 3 field components)              
      a(:) = bv(1,1,1,:)
      b(:) = bv(2,1,1,:)-bv(1,1,1,:)
      c(:) = bv(1,2,1,:)-bv(1,1,1,:)
      d(:) = bv(2,2,1,:)-bv(2,1,1,:)-bv(1,2,1,:)+bv(1,1,1,:)
      e(:) = bv(1,1,2,:)-bv(1,1,1,:)
      f(:) = bv(2,1,2,:)-bv(2,1,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      g(:) = bv(1,2,2,:)-bv(1,2,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      h(:) = bv(2,2,2,:)-bv(2,2,1,:)-bv(2,1,2,:)-bv(1,2,2,:)
      h(:) = h(:) + bv(2,1,1,:)+bv(1,2,1,:)+bv(1,1,2,:)-bv(1,1,1,:)

c interpolation being done in u = ln(r) space  
c find r from pos
      r = exp(pos(1))
      th = pos(2)

c find bv using interpolation cooefficients
      bpos(:) = a(:)+b(:)*dr+c(:)*dt+d(:)*dr*dt+e(:)*dp+f(:)*dr*dp
      bpos(:) = bpos(:)+g(:)*dt*dp+h(:)*dr*dt*dp

c Blength                                          
      blen = sqrt(bpos(1)**2+bpos(2)**2+bpos(3)**2)

c analytical expressions for dBr/dr, dBr/dt etc found by differentiating the interpolation 
c dbdr = (dbrdr,dbtdr,dbpdr) and so on 
      dbdr(:) = ((b(:) + d(:)*dt + f(:)*dp + h(:)*dt*dp)/dv(1))/r
      dbdt(:) = (c(:) + d(:)*dr + g(:)*dp + h(:)*dr*dp)/dv(2)
      dbdp(:) = (e(:) + f(:)*dr + g(:)*dt + h(:)*dr*dt)/dv(3)

c dv factor since dr etc are scaled from 0 to 1 in each direction          
c so e.g. t' = 0 -> 1, with t' = (t-t(t'=0))/(t(1)-t(0))                   
c then db/dt = db/dt'*dt'/dt = db/dt'/(t(1)-t(0))                          
c also db/dr = db/du'*du'/du*du/dr = db/du'/dv/r, where u = ln(r)

c find curlB 
      curlb(1) = (dbdt(3)*sin(th)+bpos(3)*cos(th)-dbdp(2))/(r*sin(th))
      curlb(2) = (dbdp(1)/sin(th)-r*dbdr(3)-bpos(3))/r
      curlb(3) = (dbdr(2)*r+bpos(2)-dbdt(1))/r

      curlb(:) = curlb(:)/blen

      pos(10) = curlb(1)*bpos(1)+curlb(2)*bpos(2)+curlb(3)*bpos(3)
      pos(10) = pos(10)/(real(4.0)*pi1*blen)


      end subroutine tw_sph

