MODULE global
  INTEGER, PARAMETER :: dp = 8 ! kind of real variables (double precision)
  INTEGER, PARAMETER :: n = 10 ! Lattice size
  !INTEGER, DIMENSION(n,n) :: s
  real(kind=dp) :: pi
  
END MODULE global

program onsager
 use global
 implicit none
 
 Real(kind=dp) :: E,k,T
 integer :: lattice_type
 real(kind=dp) :: alpha
 real(kind=dp) :: beta
 real(kind=dp) :: x, tc
 real(kind=dp), parameter ::j = 1.0_dp 
 pi = 2*acos(0.0)

 read(*,*) T
 T = T/100.0_dp
 k = J/T
 
 !write (*,*) "Enter lattice type 1 for square 2 for triangular"
 read(*,*) lattice_type
 
 
 x = exp(2*k)
 if ( lattice_type == 1) then
   alpha = 2*(x**4 - 6*x*x + 1)
   alpha = alpha / (pi*(x**4 -1))
   
   beta = 16*(x*x - 1)**2*x*x/((x*x+1)** 4)
 else if (lattice_type == 2) then
   !Got a triangular lattice need to find out if T < Tc 
   tc = 4.0_dp/log(3.0)
 !  write(*,*) "TC=",tc
   if ((x**2 .gt. 3) )then !.and. (t < tc)) then
   !  Write(*,*) "Opting for T<Tc, X^2 > 3"
     alpha = 4*x**2*(x**2-3)*(x**2+1)
     alpha = alpha /(pi * sqrt((x**2 - 1)**5*(x**2+3)))
     
   !  beta = 16*(x**2 - 1)**2*x**2
     !beta = beta /((x**2-1)**3*(x**2+3))
     beta = 16*x**2/((x**2-1)**3*(x**2+3))
   else if ((x**2 .lt. 3) )then !.and. (t > tc)) then
     alpha = x*(x**2-3)*(x**2+1)
     alpha = alpha /(pi * (x**2 - 1))
     
     beta = (x**2-1)**3*(x**2+3)
     beta = beta / (16*(x**2))

   else
     Write (*,*) "Solution not known for x^2=",x**2; Stop
   end if
 else 
   Write (*,*) "Invalid lattice type"; Stop
 end if
 !write (*,*) "A=",alpha, " B =", beta
 
 if(abs(beta-1._dp) > 1e-9)then
    E = -(J/tanh(2.0*K) + alpha*cel(sqrt(1.0-beta),1._dp,  1._dp, 1._dp))
    else 
    E = -(J/tanh(2.0*K))
    end if
	
	
 write (*,*) T, E
 ! Write (*,*) 1._dp/tanh(2._dp*k), alpha, beta, cel(sqrt(1.0_dp-beta), 1._dp, 1._dp, 1._dp)
 contains 
   !Function to evaluate the Elliptic Integral  1/(1-Bsin^2(x)) dx 
   !between 0 and 2pi
   real(kind=dp) Function I(b)
    real(kind=dp) , intent(in) :: b
    real(kind=dp) :: errtol, third, c1, c2,c3,c4
    real(kind=dp) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
    
    third = 1.0_dp/3.0_dp
    c1 = 1.0_dp/2.4_dp
    c2 = 0.1_dp
    c3 = 3.0_dp / 44.0_dp
    c4 = 1.0_dp /14.0_dp
    errtol = 0.0025
    
    
    
    xt = 0.0_dp
    yt = 1.0_dp - b
    zt = 1.0_dp
    do
      sqrtx = sqrt(xt)
      sqrty = sqrt(yt)
      sqrtz = sqrt(zt)
      alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt = 0.25*(xt+alamb)
      yt = 0.25*(yt+alamb)
      zt = 0.25*(zt+alamb)
      ave = third*(xt+yt+zt)
      delx = (ave-xt)/ave
      dely = (ave-yt)/ave
      delz = (ave-zt)/ave
      if (max(max(abs(delx), abs(dely)), abs(delz)) < errtol) then
        exit
      end if
	  
    end do
    e2 = delx*dely - delz*delz
    e3 = delx*dely*delz
    I = (1.0+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave)
   end function I
   
   real(kind=dp) Function el2(x, kc,aa,bb)
     real(kind=dp), intent(in) :: x, kc, aa,bb
     real(kind=dp) ::  qc, b,a, d, p,z,eye,f,l,em,g,y,c,e
     
     real(kind=dp), parameter :: CA=0.0003, cb=1.e-9
     if(x.eq. 0 ) then
       el2 = 0.0_dp
     else if (kc .ne.0) then
       qc = kc
       a=aa
       b=bb
       c = x**2
       d=1._dp+c
       p = sqrt((1.+qc**2*c)/D)
       D = x/D
       c = d/(2._dp*p)
       z = a-b
       eye = a
       a = 0.5_dp*(b+a)
       y = abs(1._dp/x)
       f = 0.0_dp
       l = 0.0_dp
       em = 1.0_dp
       qc = abs(qc)
       do
         b = eye*qc+b
         e = em*qc
         g = e/p
         d=f*g+d
         f=c
         eye=a
         p=g+p
         c = 0.5_dp*(d/p+c)
         g = em
         em = qc+em
         a = 0.5_dp*(b/em+a)
         y = -e/y+y
         if (y.eq.0._dp) y=sqrt(e)*cb
           if(abs(g-qc).gt.ca*g)then
             qc=sqrt(e)*2._dp
             
             l=l+l
           if(y.lt.0._dp) l=l+1
         else 
           exit
         end if
       end do
       if(y.lt.0._dp)l=l+1
       E = (atan(em/y)+pi*l)*a/em
       if(x.lt.0._dp)e=-E
       el2=e+c*z
     end if
   end function el2
   
   real(kind=dp) function cel(qqc, pp,aa,bb)
     real(kind=dp), intent(in) :: qqc, aa, pp, bb
     real(kind=dp) :: qc,a,b,p,em,e,f,g,q
     real(kind = dp), parameter :: ca=.0003
     real(kind=dp) :: pio2 
     pio2 = pi/2
     if (qqc.eq.0) stop 'failuer in cel'
     qc = abs(qqc)
     a = aa
     b = bb
     p=pp
     e = qc
     em = 1
     if (p.gt.0._dp) then
     p = sqrt(p)
     b=b/p
     else 
     f = qc*qc
     q = 1._dp - f
     g = 1._dp -p
     f = f-p
     q=q*(b-a*p)
     p = sqrt(f/g)
     a = (a-b)/g
     b = -q/(g*g*p)+a*p
     end if
     do
       f = a
       a=a+b/p
       g = e/p
       b = b+f*g
       b = b+b
       p=g+p
       g=em
       em=qc+em
       if(abs(g-qc).gt. g*ca) then
         qc=sqrt(e)
         qc=qc+qc
         e = qc*em
       else 
         exit
       end if
     end do
     cel = pio2*(b+a*em)/(em*(em+p))
     end function cel
       
       
         
   
   
   
end program onsager
