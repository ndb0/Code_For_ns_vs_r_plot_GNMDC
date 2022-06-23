


module mod1

!use precision
!use ieee_arithmetic
implicit none
!integer, parameter :: dl =SELECTED_REAL_KIND(P=15)
integer, parameter :: dl =SELECTED_REAL_KIND(R=200)
integer, parameter :: dls =SELECTED_REAL_KIND(R=200)

!double precision, parameter:: fg=2.22072,b= 1.0998,c=0.266252
!double precision, parameter::d= 2.67098933990!139463318799079794,nu=.55293312,cof=4.5E-10
!double precision, parameter::nu=.55293312,cof=4.5E-10
double precision, parameter:: nef=75.0d0,phi_ini=9.10d0,a00=4.8e-5,PI=4.D0*DATAN(1.D0)
double precision, parameter:: M_pl=1.d0,k_end=4e16,dhv00=1e-10

double precision :: error(50)
double precision :: css(800000),cssn(800000),ctt(800000),zns(800000),znt(800000)
integer :: error_index
!double precision m=7.147E-6 !mass of the inflation
contains
		
!((a*phi**2)/nu**2+(b*phi**4)/nu**4+(c*phi**6)/nu**6)/(1+(d*phi**2)/nu**2)**2


!p2=(nu*(1+(d*phi**2)/nu**2)**3)
!(2*phi*(a+(2*b*phi**2)/nu**2-(a*d*phi**2)/nu**2+(*c*phi**4*(3+(d*phi**2)/nu**2))/nu**4))/p2


!p4=(3*a*d**2*phi**4)/nu**4+(4*c*d*phi**6)/nu**6+(*c*d**2*phi**8)/nu**8
!p3=a+(6*b*phi**2)/nu**2-(8*a*d*phi**2)/nu**2+(15*c*phi**4)/nu**4
!p2=(p3-(6*b*d*phi**4)/nu**4+p4)
!p1=(1+(d*phi**2)/nu**2)**4
!(2*p2)/p1








	real(dl)   function  V(coff,phi)  
!potential function
		implicit none
		real(dl)   y,m,phi0,z,phi
		real(dl) coff(25)
		real(dl) fg,b,c,d,nu,cof
		cof=coff(6)
		!!print **,fg,b,c,d,cof,nu
		V=cof*0.5*phi**2
		
	return
	end function
	
	
	
		!################################option to chose#########################################
! G inflation term : arXiv:2001.05909
	real(dl)   function  theta(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25)
		
		theta=theta_sin(coff,phi)+theta_mono(coff,phi) 
		

        return
	

end function

!1st derivative
	real(dl)   function  thetaphi(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25)
		
		thetaphi=thetaphi_sin(coff,phi)+thetaphi_mono(coff,phi)
		

        return
	

end function

!2nd derivative of G_inf wrt
	real(dl)   function  thetaphiphi(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25)
		
		thetaphiphi=thetaphiphi_sin(coff,phi)+thetaphiphi_mono(coff,phi)
		

        return
	end function
	
	
		real(dl)   function  thetaphiphiphi(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25)
		
		thetaphiphiphi=0.0!thetaphiphiphi_bessl(phi,theta)
		

        return
	end function
	

	
!________________________________________________________________________________________________!
	

	!################################sin/()#########################################
! G inflation term : arXiv:2001.05909
	real(dl)   function  theta_sin(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,t,f1,factor
		
		A0=coff(18)/coff(6)*coff(23)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
                t=coff(25)
                
               
                
                theta_sin=(A0*Sin(c1+b1*(1.d0+s2*Log(phic/phi))*(phi-  &
		1.d0*phic)**m))/(1.d0+(phi-1.d0*phic)**p/sigma)**(1.d0*q)
		
        return
	

end function

real(dl)   function  theta_mono(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,t,f1,factor,A1		
		
		A1=coff(16)/coff(6)
	!	phic=coff(1)*coff(15)
                t=coff(25)
                
         
         
                
                theta_mono=A1*(phi)**t
		
        return
	

end function

	real(dl)   function  theta_only_sin(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,t,f1
		
		A0=coff(18)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
                t=coff(25)
               !f1=coff(18)
                !f1=0.4
                
                theta_only_sin=Sin(c1+b1*(phi-1.d0*phic)**m*(1.d0+  &
		s2*Log(phic/phi)))/(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(1.d0*q)
                !theta_sin=A0*((phi/phic)**t+(f1*Sin(c1+b1*(1.d0+s2*Log(phic/phi))*(phi  &
                !-1.d0*phic)**m))/(1.d0+(phi-1.d0*phic)**p/sigma)**(1.d0*q))

		
!            theta_sin=A0*(((Sin(c1+b1*(phi-1.d0*phic)**m*(1.d0+  &
 !   s2*Log(phic/phi))))/(1.d0+(phi-  &
!    1.d0*phic)**p/sigma)**(1.d0*q)))+0.4*A0*phi**4/phic**4
        return
	

	end function


		real(dl)   function  thetaphi_sin(coff,phi)  
		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,f1,t,factor
		A0=coff(18)/coff(6)*coff(23)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
              !  f1=coff(18)
               !f1=0.4
                t=coff(25)
               thetaphi_sin=(-1.d0*A0*p*q*Sin(c1+b1*(1.d0+s2*Log(phic/phi))*(phi-  &
		1.d0*phic)**m)*(phi-1.d0*phic)**(-1.d0+p)*(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(-1.d0-1.d0*q))/sigma+(A0*Cos(c1+  &
		b1*(1.d0+s2*Log(phic/phi))*(phi-1.d0*phic)**m)*(b1*m*(1.d0+  &
		s2*Log(phic/phi))*(phi-1.d0*phic)**(-1.d0+m)-  &
		(1.d0*b1*s2*(phi-1.d0*phic)**m)/phi))/(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(1.d0*q)
                
                
               
	return
	end function
	
	real(dl)   function  thetaphi_mono(coff,phi)  
		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,f1,t,factor,A1
		A1=coff(16)/coff(6)
	!	phic=coff(1)*coff(15)
                t=coff(25)
               thetaphi_mono=(A1*t*(phi)**(-1.d0+t))       
               
	return
	end function
	
	
	!2nd derivative of G_inf wrt
		real(dl)   function  thetaphiphi_sin(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,t,f1,factor
		A0=coff(18)/coff(6)*coff(23)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4)
		sigma=coff(3)
		p=coff(20)
	        s2=coff(19)
                q=coff(21)        
                m=coff(22)
               !f1=coff(18)
      		!f1=0.4
                t=coff(25)

		thetaphiphi_sin=(-1.d0*A0*p**2*(-1.d0-1.d0*q)*q*Sin(c1+b1*(1.d0+  &
		s2*Log(phic/phi))*(phi-1.d0*phic)**m)*(phi-  &
		1.d0*phic)**(-2.d0+2.d0*p)*(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(-2.d0-1.d0*q))/sigma**2-  &
		(2.d0*A0*p*q*Cos(c1+b1*(1.d0+s2*Log(phic/phi))*(phi-  &
		1.d0*phic)**m)*(phi-1.d0*phic)**(-1.d0+p)*(b1*m*(1.d0+  &
		s2*Log(phic/phi))*(phi-1.d0*phic)**(-1.d0+m)-  &
		(1.d0*b1*s2*(phi-1.d0*phic)**m)/phi)*(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(-1.d0-1.d0*q))/sigma-(1.d0*A0*(-1.d0+  &
		p)*p*q*Sin(c1+b1*(1.d0+s2*Log(phic/phi))*(phi-  &
		1.d0*phic)**m)*(phi-1.d0*phic)**(-2.d0+p)*(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(-1.d0-1.d0*q))/sigma+(A0*Cos(c1+  &
		b1*(1.d0+s2*Log(phic/phi))*(phi-1.d0*phic)**m)*(b1*(-1.d0+  &
		m)*m*(1.d0+s2*Log(phic/phi))*(phi-1.d0*phic)**(-2.d0+m)-  &
		(2.d0*b1*m*s2*(phi-1.d0*phic)**(-1.d0+m))/phi+(b1*s2*(phi-  &
		1.d0*phic)**m)/phi**2))/(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(1.d0*q)-(1.d0*A0*Sin(c1+b1*(1.d0+  &
		s2*Log(phic/phi))*(phi-1.d0*phic)**m)*(b1*m*(1.d0+  &
		s2*Log(phic/phi))*(phi-1.d0*phic)**(-1.d0+m)-  &
		(1.d0*b1*s2*(phi-1.d0*phic)**m)/phi)**2)/(1.d0+(phi-  &
		1.d0*phic)**p/sigma)**(1.d0*q)
            
       ! if ( thetaphiphi != thetaphiphi ) then
                     !  write(*,*) 'thetaphiphi is NaN'
                      ! !print **, omega,phic,sigma,dentst
                      ! !print **, "omega,phic,sigma"
                      ! !print **, omega,phic,sigma
                      ! read(*,*) rff
        !              thetaphiphi=0.d0
        !endif
		
	return
	end function
	
	real(dl)   function  thetaphiphi_mono(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m,t,f1,factor,A1
		A1=coff(16)/coff(6)
		phic=coff(1)*coff(15)
                t=coff(25)

		thetaphiphi_mono=(A1*(-1.d0+t)*t*(phi)**(-2.d0+t))
       ! if ( thetaphiphi != thetaphiphi ) then
                     !  write(*,*) 'thetaphiphi is NaN'
                      ! !print **, omega,phic,sigma,dentst
                      ! !print **, "omega,phic,sigma"
                      ! !print **, omega,phic,sigma
                      ! read(*,*) rff
        !              thetaphiphi=0.d0
        !endif
		
	return
	end function

	
		!3rd derivative of G_inf wrt
		real(dl)   function  thetaphiphiphi_sin(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma
		integer p
		A0=coff(16)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4)
		sigma=coff(3)
		p=4
		
                thetaphiphiphi_sin=0.0
            
       ! if ( ieee_is_nan(thetaphiphiphi) ) then
                     !  write(*,*) 'thetaphiphiphi is NaN'
                      ! !print **, omega,phic,sigma,dentst
                      ! !print **, "omega,phic,sigma,dentist"
                      ! !print **, omega,phic,sigma,dentist
                      ! read(*,*) rff
        !              thetaphiphiphi=0.d0
        !endif            
		
	return
	end function
	
	
!________________________________________________________________________________________________!
	
!##########################################Bessel-fn-part################################

	real(dl)   function  theta_bessl(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m
		
		A0=coff(16)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)

		theta_bessl=(-1.d0*A0*bessel_JN(1,b1+c1*exp((c1*(-1.d0*phi+  &
            phic))/phic)))/(b1+c1*exp((c1*(-1.d0*phi+phic))/phic))
        return
	
 
end function

	real(dl)   function  thetaphi_bessl(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m
		
		A0=coff(16)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
		thetaphi_bessl=(-1.d0*A0*c1**2*exp(c1)*bessel_JN(2,b1+c1*exp(c1-  &
            (1.d0*c1*phi)/phic)))/(c1*exp(c1)*phic+  &
            b1*exp((c1*phi)/phic)*phic)
        return
	

end function


	real(dl)   function  thetaphiphi_bessl(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m
		
		A0=coff(16)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
		thetaphiphi_bessl=(A0*c1**3*exp(c1)*(c1*(b1*exp(c1)+c1*exp(c1*(2.d0-  &
            (1.d0*phi)/phic)))*bessel_JN(1,b1+c1*exp(c1-  &
            (1.d0*c1*phi)/phic))+(-2.d0*c1*exp(c1)+  &
            b1*exp((c1*phi)/phic))*bessel_JN(2,b1+c1*exp(c1-  &
            (1.d0*c1*phi)/phic))))/((c1*exp(c1)+  &
            b1*exp((c1*phi)/phic))**2*phic**2)
                    return
	

end function


	real(dl)   function  thetaphiphiphi_bessl(coff,phi)  

		implicit none
		real(dl)   phi
		real(dl) coff(25),A0,phic,b1,c1,cof,sigma,s2,p,q,m
		
		A0=coff(16)/coff(6)
		b1=coff(2)
		phic=coff(1)*coff(15)
		c1=coff(4) !coff(3)
		sigma=coff(3)
		p=coff(20)
		s2=coff(19)
                q=coff(21)
                m=coff(22)
		thetaphiphiphi_bessl=0.0
        return
	

end function

!________________________________________________________________________________________________!



!1st partial derivative of potential function wrt efolds
	  real(dl) function V_phi(coff,phi)  

		implicit none
		real(dl) y1,m,phi,H,phi_n,p1,p2
		real(dl) coff(25)
		real(dl) cof
		cof=coff(6)
		V_phi=cof*phi
                return
	end function

	!1st derivative of G_inf wrt phi

	
	!2nd partial derivative of potential function wrt phi
	!2nd partial derivative of potential function wrt phi
	real(dl)  function  V_phi_phi(coff,phi)  

	implicit none
		real(dl)   y2,m,phi,p1,p2,p3,p4
		real(dl) coff(25)
		real(dl) cof
		cof=coff(6)
		V_phi_phi=cof
	return
	end function


	
		real(dl)  function  V_phi_phi_phi(coff,phi)  

	implicit none
		real(dl)   phi
		real(dl) coff(25)
		real(dl) cof
		cof=coff(6)


	        V_phi_phi_phi=0.0
	return
	end function
	

	
	real(dl) function wind(k,R)

		implicit none

		real(dl) :: k,R
		!R=1e15
		!!print **,"kr",k,r
		wind=exp(-((k*R)**2)/2.d0)
	return
	end function

	    subroutine write_info(x0,coff,info,filename)
!use mod1
implicit none
		real(dl) coff(25),info(40),eff_efold,x0
		integer i
		character(90) filename
                		OPEN(unit=1555,file=filename ,action='write',position='append')
		write(1555,*)"##############################New set #############################################"
		write(1555,*)"phi_ini=coff(15)=",coff(15),"efold=coff(9)=",coff(9),"x0=",x0 
		write(1555,*) "nonminimal coupling term=coff(16)=",coff(16)
		eff_efold=coff(9)-info(33)
        write(1555,*) "efold at k=.05=info(33):",info(33), "     " , "effective no of efols :", eff_efold
		write(1555,*)"a,b,c,d,nu,V0=(coff(1-6)=",(coff(i),i=1,6)
		write(1555,*)"a00=coff(8)",coff(8),"constr=info(17)",info(17)
		write(1555,*)".............................................."
		write(1555,*)".............................................."
		write(1555,*) "coff(12)=f_PSn(Press-Sc-numerical)=",coff(12),"coff(14)=f_Peakn(peaks-th-numerical)=",coff(14)
		write(1555,*) "info(10)=max_m1(max-pbh-PS-mass)=",info(10),"info(12)=max_m1(max-pbh-peaks-mass)=",info(12)
		write(1555,*) "info(1)=ns05(from slow-roll-para)",info(1),"info(19)=ns002(from slow-roll-para)=",info(19)
		write(1555,*) "info(14)=ns05(from slope)=",info(14),"info(22)=ns002(from-slope)",info(22)
		write(1555,*)".............................................."
		write(1555,*)".............................................."
		
		write(1555,*) "coff(11)=f_PSn1(slowroll)=",coff(11),"coff(13)=f_Peakn1(slowroll)=",coff(13)
		write(1555,*) "info(11)=max_betaN(max_beta-ps)",info(11),"info(13)=max_peakN(max_beta-peaks)=",info(13)
		write(1555,*) "info(6)=max_ps(max-power)",info(6),"info(8)=max_i(max-k)",info(8)
		write(1555,*) "info(7)=min_ps(min-power)",info(7),"info(9)=min_i(min-k)",info(9)
		
		write(1555,*) "info(2)=ep05",info(2),"info(20)=ep002",info(20)
		write(1555,*) "info(3)=eta05",info(3),"info(21)=eta002",info(21)
		write(1555,*) "info(4)=ns(linear-fit)=",info(4),"info(5)=ns(linear-fit2)",info(5)

		
		write(1555,*) "info(15,16,18)=nss05,nt05,tp_to_sp_ratio05(slope_of_ns,temsor_index,tp05/sp05)=",info(15),info(16),info(18)
		write(1555,*) "info(23-15)=nss002,nt002,tp_to_sp_ratio002(slope_of_ns,temsor_index,tp002/sp002)=",info(23),info(24),info(25)
		write(1555,*) "minimum value of phi_N=info(26)=",info(26),"at info(28)=k_max=",info(28)
		write(1555,*) "maximum power=info(27)",info(27)
		write(1555,*) "pbh-mass=info(27)",info(29)
		write(1555,*) "at min phi_n, phi, h, n",info(30),info(31),info(32)
		close(1555)
                
              return
      end subroutine        
		
end module mod1





