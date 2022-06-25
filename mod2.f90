module mod2
!use precision
use mod1
implicit none
contains


!2nd derivative of scaler field phi wrt efolds(n)
		real(dl)  function  phi_nn(coff,phi201,H201,phi_n201)
					use mod1
			
		implicit none
		real(dl)   phi201,H201,phi_n201
		real(dl) coff(25)
		!gterm=-((H*phi_n)**2)*Gphi_inf(coff,phi)
		!gterm1=1.d0-2.d0*G_inf(coff,phi)
		
		
        real(dl) theta1201,v_phi1201,thetaphi1201

		!H_n1=H_n(coff,phi,H,phi_n)
		theta1201=theta(coff,phi201)
		v_phi1201=V_phi(coff,phi201)
		thetaphi1201=thetaphi(coff,phi201)
		!phi_nn1=phi_nn(coff,phi,H,phi_n)
		!thetaphiphi1=thetaphiphi(coff,phi)
		!v_phi_phi1=V_phi_phi(coff,phi)
		!h_nn1=h_nn(coff,phi,phi_n,H)
		!phi_nnn1=phi_nnn(coff,phi,phi_n,H)
		!V_phi_phi_phi1=V_phi_phi_phi(coff,phi)
		!thetaphiphiphi1=thetaphiphiphi(coff,phi)
		
		
		
		phi_nn=(0.5d0*(V_phi1201*(-4.d0+ 6.d0*H201**2*theta1201*phi_N201**2)-  &
        1.d0*H201**2*phi_N201*(12.d0-2.d0*phi_N201**2-108.d0*H201**4*theta1201**2*phi_N201**2+  &
        2.d0*H201**2*thetaphi1201*phi_N201*(3.d0+phi_N201**2)+  &
                3.d0*H201**2*theta1201*(12.d0- 14.00d0*phi_N201**2+  &
                3.d0*H201**2*thetaphi1201*phi_N201**3))))/(H201**2*(2.d0+  &
                9.d0*H201**4*theta1201**2*phi_N201**2-1.d0*H201**2*theta1201*(-6.d0+phi_N201**2)))
		return
		end function

!1st derivative of Hubble parameter H  wrt efolds(n)
		real(dl)  function  H_n(coff,phi202,H202,phi_n202)
						use mod1
			
		implicit none
		real(dl)   phi202,H202,phi_n202
		real(dl) coff(25)
        real(dl) theta1202,v_phi1202,thetaphi1202

		!H_n1=H_n(coff,phi,H,phi_n)
		theta1202=theta(coff,phi202)
		v_phi1202=V_phi(coff,phi202)
		thetaphi1202=thetaphi(coff,phi202)
		
		!gterm5=(-G_inf(coff,phi) +1.d00)
		H_n=(-1.d0*H202*phi_N202*(phi_N202+  27.d0*H202**4*theta1202**2*phi_N202-  &
                1.d0*H202**2*thetaphi1202*phi_N202**2+2.d0*theta1202*(6.d0*H202**2*phi_N202+  &
                V_phi1202)))/(2.d0+ 9.d0*H202**4*theta1202**2*phi_N202**2-  &
                1.d0*H202**2*theta1202*(-6.d0+phi_N202**2))
		return 
		end function






! the function z=R/q=a.phi_t/H

		real(dl) function zzzz(coff,phi_n203,a203,H203,phi203)
					use mod1
			
		implicit none
		real(dl) phi_n203,a203,H203,phi203
		real(dl) coff(25)
		
        real(dl) theta1203,v_phi1203,thetaphi1203

		!H_n1=H_n(coff,phi,H,phi_n)
		theta1203=theta(coff,phi203)
		v_phi1203=V_phi(coff,phi203)
		thetaphi1203=thetaphi(coff,phi203)
		
		
		!gterm2=sqrt(1.d0-2.d0*G_inf(coff,phi))
		
		
		zzzz=a203*sqrt((-1.0*phi_n203**2*(-2.0+  &
		H203**2*theta1203*phi_n203**2)*(2.0+  &
		9.0*H203**4*theta1203**2*phi_n203**2-  &
		1.0*H203**2*theta1203*(-6.0+  &
		phi_n203**2)))/(2.0-  &
		3.0*H203**2*theta1203*phi_n203**2)**2)
		! !print **,zzzz,zzzz,phi_n,a,H
		!
		return
		end function


!1st derivative of (aH)  wrt efolds(n)
!		real(dl) function aH_n(coff,a,H,phi_n,phi)
!					use mod1
!			
!		implicit none
!		real(dl) phi_n,a,H,H_n2,phi
!		real(dl) coff(25)
	!	H_n2=(theta(coff,phi) -1.d00)*(phi_n**2)*H/2.d0
!		aH_n=a*H+a*H_n2
!		return
!		end

 

end module mod2
