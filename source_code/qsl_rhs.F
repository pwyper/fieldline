      subroutine qsl_cart (bvec,dv,corner,pos,len,kvec)
c returns the kvector for uvec and vvec advancement in rk4
c Peter Wyper June 2023
c updated to correct form for Q (Aslanyan et al. 24) April 2025
c inputs 
      real pos(10)
      real dv(3),corner(3)
      real bvec(9,9,9,3)
      real kvec(9), len

c internal 
      real bv(2,2,2,3)
      real b(3),c(3),d(3),e(3),f(3),g(3),h(3)
      real dx,dy,dz
      integer i,j,k
      real dbdx(3),dbdy(3),dbdz(3)
      real bv1(3), blen

c find |B| at pos
      call interp_block (bvec,dv,corner,pos,bv1)
      blen = sqrt(bv1(1)**2+bv1(2)**2+bv1(3)**3)
      
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
      b(:) = bv(2,1,1,:)-bv(1,1,1,:) 
      c(:) = bv(1,2,1,:)-bv(1,1,1,:)
      d(:) = bv(2,2,1,:)-bv(2,1,1,:)-bv(1,2,1,:)+bv(1,1,1,:)
      e(:) = bv(1,1,2,:)-bv(1,1,1,:)
      f(:) = bv(2,1,2,:)-bv(2,1,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      g(:) = bv(1,2,2,:)-bv(1,2,1,:)-bv(1,1,2,:)+bv(1,1,1,:)
      h(:) = bv(2,2,2,:)-bv(2,2,1,:)-bv(2,1,2,:)-bv(1,2,2,:)
      h(:) = h(:) + bv(2,1,1,:)+bv(1,2,1,:)+bv(1,1,2,:)-bv(1,1,1,:)

c analytical expressions for dBx/dx, dBx/dy etc found by differentiating the interpolation
c dbdx = (dbxdx,dbydx,dbzdx) and so on
      dbdx(:) = (b(:) + d(:)*dy + f(:)*dz + h(:)*dy*dz)/dv(1)
      dbdy(:) = (c(:) + d(:)*dx + g(:)*dz + h(:)*dx*dz)/dv(2)
      dbdz(:) = (e(:) + f(:)*dx + g(:)*dy + h(:)*dx*dy)/dv(3)

c dv factor since dx etc are scaled from 0 to 1 in each direction
c so e.g. x' = 0 -> 1, with x' = (x-x(x'=0))/(x(1)-x(0))
c then db/dx = db/dx'*dx'/dx = db/dx'/(x(1)-x(0))

c create kvec for uvec and vvec
c kvec = 1/|B|*(U dot del)B * len - eqns 15 to 17 Aslanyan et al. 24 
c pos(4:6) = kvec for vvec 
c pos(7:9) = kvec for uvec
      len = len/blen
      kvec(4:6) = (pos(4)*dbdx(:)+pos(5)*dbdy(:)+pos(6)*dbdz(:))*len
      kvec(7:9) = (pos(7)*dbdx(:)+pos(8)*dbdy(:)+pos(9)*dbdz(:))*len

      end subroutine qsl_cart


      subroutine qsl_sph (bvec,dv,corner,pos,len,kvec)
c returns the kvector for uvec and vvec advancement in rk4 in spherical coordinates                    
c Peter Wyper June 2023
c updated to correct form for Q (Aslanyan et al. 24) April 2025
c inputs                                                                       
      real pos(10)
      real dv(3),corner(3)
      real bvec(9,9,9,3)
      real kvec(9), len
c internal                                                                     
      real bv(2,2,2,3)
      real a(3),b(3),c(3),d(3),e(3),f(3),g(3),h(3)
      real dr,dt,dp
      integer i,j,k
      real dbdr(3),dbdt(3),dbdp(3)
      real r, th, bpos(3)
      real blen
      
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

c spherical (u.grad)bhat requires bv at pos
c find bv using interpolation cooefficients
      bpos(:) = a(:)+b(:)*dr+c(:)*dt+d(:)*dr*dt+e(:)*dp+f(:)*dr*dp
      bpos(:) = bpos(:)+g(:)*dt*dp+h(:)*dr*dt*dp

c find blen
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

c create kvec for uvec and vvec 
c kvec = 1/|B|*(U dot del)B * len - eqns 18 to 20 Aslanyan et al. 24    
c pos(4:6) = delta xvec for vvec
c pos(7:9) = delta xvec for uvec

c (u.grad)b -> each component
      kvec(4:6)=pos(4)*dbdr(:)+pos(5)*dbdt(:)/r
      kvec(4:6)=kvec(4:6)+pos(6)*dbdp(:)/(r*sin(th))
c geometrical factors
c dVth/dl --> + (Br*Vth-Vr*Bth)/r       
      kvec(5) = kvec(5) + (bpos(1)*pos(5)-pos(4)*bpos(2))/r
c dVphi/dl --> + (Br*Vphi-Vr*Bphi)/r + (Bth*Vphi-Vth*Bphi)/(r*tan(th))
      kvec(6) = kvec(6) + (bpos(1)*pos(6)-pos(4)*bpos(3))/r
      kvec(6) = kvec(6) + (bpos(2)*pos(6)-pos(5)*bpos(3))/(r*tan(th))
      
c (v.grad)b -> each component                  
      kvec(7:9)=pos(7)*dbdr(:)+pos(8)*dbdt(:)/r
      kvec(7:9)=kvec(7:9)+pos(9)*dbdp(:)/(r*sin(th))
c geometrical factors
c dUth/dl --> + (Br*Uth-Ur*Bth)/r       
      kvec(8) = kvec(8) + (bpos(1)*pos(8)-pos(7)*bpos(2))/r
c dUphi/dl --> + (Br*Uphi-Ur*Bphi)/r + (Bth*Uphi-Uth*Bphi)/(r*tan(th))
      kvec(9) = kvec(9) + (bpos(1)*pos(9)-pos(7)*bpos(3))/r
      kvec(9) = kvec(9) + (bpos(2)*pos(9)-pos(8)*bpos(3))/(r*tan(th))


      kvec(4:9) = kvec(4:9)*len/blen

      end subroutine qsl_sph

