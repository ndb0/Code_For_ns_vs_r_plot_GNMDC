module mod3
!use precision
use mod1
use mod2
!use ieee_arithmetic
implicit none
contains

!!check github
!this subroutine provides the derivative values for ode solver in array format in background part 
	subroutine  f(coff,ff251,r251,t251,k251,hn)
			use mod1
			use mod2
			
		    implicit none
		    real(dl)  ,dimension(3)::r251
		    real(dl)  ,dimension(3)::ff251
		    real(dl)  :: t251,k251,hn,phi2251,phi_n2251,fphi251,fphi_n251,fH251,H2251
                    real(dl) coff(25)
		    phi2251 = r251(1)
		    phi_n2251 = r251(2)
		    H2251=r251(3)
		    fphi251 = phi_n2251
		    fphi_n251=  phi_nn(coff,phi2251,H2251,phi_n2251)
                   
		    fH251=H_n(coff,phi2251,H2251,phi_n2251)
		    
		    ff251=[fphi251,fphi_n251,fH251]
		   ! !print **, fff, ff
		   ! !print **, r1, r
!		return
		end subroutine

		!i 2nd derivative of hubble parameter wrt efolds
		real(dl) function h_nn(coff,phi252,phi_n252,H252)
					use mod1
			use mod2
		implicit none
		real(dl) phi_n252,H252,phi252
		real(dl) coff(25)
		
        real(dl) h_n1252,theta1252,v_phi1252,thetaphi1252,phi_nn1252,thetaphiphi1252
		real(dl) v_phi_phi1252

		H_n1252=H_n(coff,phi252,H252,phi_n252)
		theta1252=theta(coff,phi252)
		v_phi1252=V_phi(coff,phi252)
		thetaphi1252=thetaphi(coff,phi252)
		!phi_nn1=phi_nn(coff,phi,H,phi_n)
		thetaphiphi1252=thetaphiphi(coff,phi252)
		v_phi_phi1252=V_phi_phi(coff,phi252)
		!h_nn1=h_nn(coff,phi,phi_n,H)
		!phi_nnn1=phi_nnn(coff,phi,phi_n,H)
		!V_phi_phi_phi1=V_phi_phi_phi(coff,phi)
		!thetaphiphiphi1=thetaphiphiphi(coff,phi)
		
				h_nn= (0.5d0*(13122.d0*H252**14*theta1252**6*phi_N252**6-  &
            162.d0*H252**10*theta1252**5*phi_N252**4*(27.d0*H252**2*(2.d0  &
            +phi_N252**2)+2.d0*phi_N252**2*V_phi_phi1252)+  &
            18.d0*H252**6*theta1252**4*phi_N252**2*(-6.d0*V_phi1252**2*phi_N252**2  &
            +9.d0*H252**4*(72.d0-66.d0*phi_N252**2+(-38.d0+  &
            thetaphi1252*V_phi1252)*phi_N252**4)-  &
            2.d0*H252**2*phi_N252*(V_phi1252*(72.d0+  &
            39.d0*phi_N252**2)-2.d0*phi_N252*(-6.d0+  &
            phi_N252**2)*V_phi_phi1252)+  &
            9.d0*H252**6*phi_N252**5*(21.d0*thetaphi1252+  &
            phi_N252*thetaphiphi1252))-  &
            1.d0*H252**4*theta1252**3*phi_N252*(405.d0*H252**8*thetaphi1252**2*phi_N252**7  &
            +4.d0*V_phi1252**2*phi_N252*(48.d0+phi_N252**2)-  &
            36.d0*H252**4*phi_N252*(360.d0+6.d0*(-23.d0+  &
            thetaphi1252*V_phi1252)*phi_N252**2+(-35.d0+  &
            4.d0*thetaphi1252*V_phi1252)*phi_N252**4)+  &
            4.d0*H252**2*(6.d0*V_phi1252*(-72.d0+54.d0*phi_N252**2+  &
            11.d0*phi_N252**4)+phi_N252*(36.d0+24.d0*phi_N252**2+  &
            phi_N252**4)*V_phi_phi1252)+  &
            36.d0*H252**6*phi_N252**4*(-9.d0*thetaphi1252*(10.d0+  &
            3.d0*phi_N252**2)+phi_N252*(-6.d0+  &
            phi_N252**2)*thetaphiphi1252))+  &
            theta1252**2*(27.d0*H252**10*thetaphi1252**2*phi_N252**6*(-20.d0  &
            +phi_N252**2)-  &
            2.d0*H252**6*phi_N252**2*(thetaphi1252*V_phi1252*(-36.d0  &
            -72.d0*phi_N252**2+phi_N252**4)+27.d0*(-96.d0+20.d0*phi_N252**2  &
            +phi_N252**4))+48.d0*H252**2*V_phi1252**2+  &
            4.d0*H252**4*phi_N252*(V_phi1252*(324.d0-  &
            48.d0*phi_N252**2+phi_N252**4)+4.d0*phi_N252*(-6.d0+  &
            phi_N252**2)*V_phi_phi1252)+  &
            2.d0*H252**8*phi_N252**3*(3.d0*thetaphi1252*(-180.d0+  &
            312.d0*phi_N252**2+7.d0*phi_N252**4)+phi_N252*(36.d0+  &
            24.d0*phi_N252**2+phi_N252**4)*thetaphiphi1252))+  &
            4.d0*phi_N252*(V_phi1252*(4.d0-  &
            10.d0*H252**2*thetaphi1252*phi_N252)+H252**2*phi_N252*(12.d0-  &
            1.d0*phi_N252**2+  &
            H252**4*thetaphi1252**2*phi_N252**2*(-15.d0+phi_N252**2)+  &
            2.d0*H252**2*phi_N252*(-15.d0*thetaphi1252+  &
            phi_N252*thetaphiphi1252)))+  &
            2.d0*theta1252*(4.d0*H252**2*V_phi1252*phi_N252*(36.d0  &
            -2.d0*phi_N252**2+3.d0*H252**2*thetaphi1252*phi_N252*(-4.d0  &
            +phi_N252**2))+8.d0*V_phi1252**2-  &
            1.d0*H252**2*phi_N252**2*(-1.d0*H252**2*(432.d0-54.d0*phi_N252**2+  &
            phi_N252**4)+H252**6*thetaphi1252**2*phi_N252**2*(90.d0+  &
            60.d0*phi_N252**2+phi_N252**4)+8.d0*V_phi_phi1252+  &
            4.d0*H252**4*phi_N252*(thetaphi1252*(90.d0-  &
            33.d0*phi_N252**2)+phi_N252*(-6.d0+  &
            phi_N252**2)*thetaphiphi1252)))))/(H252*(2.d0+  &
            9.d0*H252**4*theta1252**2*phi_N252**2-  &
            1.d0*H252**2*theta1252*(-6.d0+phi_N252**2))**3)
		return
		end function
		
!3rd derivative of scalar field wrt efolds
		real(dl) function phi_nnn(coff,phi253,phi_n253,H253)
					use mod1
			use mod2
		implicit none
		real(dl) phi_n253,H253,phi253
		real(dl) coff(25)
		
        real(dl) h_n1253,theta1253,v_phi1253,thetaphi1253,phi_nn1253
		real(dl) v_phi_phi1253,thetaphiphi1253

		H_n1253=H_n(coff,phi253,H253,phi_n253)
		theta1253=theta(coff,phi253)
		v_phi1253=V_phi(coff,phi253)
		thetaphi1253=thetaphi(coff,phi253)
		!phi_nn1=phi_nn(coff,phi,H,phi_n)
		thetaphiphi1253=thetaphiphi(coff,phi253)
		v_phi_phi1253=V_phi_phi(coff,phi253)
		!h_nn1=h_nn(coff,phi,phi_n,H)
		!phi_nnn1=phi_nnn(coff,phi,phi_n,H)
		!V_phi_phi_phi1=V_phi_phi_phi(coff,phi)
		!thetaphiphiphi1=thetaphiphiphi(coff,phi)
		
phi_nnn=(0.5d0*(8.d0*theta1253*V_phi1253**2*phi_N253*(-8.d0 &
                    +54.d0*H253**6*theta1253**3*phi_N253**4-  &
                    3.d0*H253**4*theta1253**2*phi_N253**2*(-6.d0+phi_N253**2)+  &
                    10.d0*H253**2*theta1253*(-6.d0+phi_N253**2))+  &
                    2.d0*V_phi1253*(24.d0-20.d0*phi_N253**2+  &
                    81.d0*H253**8*theta1253**3*phi_N253**4*(-2.d0*thetaphi1253*phi_N253*(2.d0  &
                    +phi_N253**2)+theta1253*(30.d0+31.d0*phi_N253**2))+  &
                    4.d0*H253**2*(4.d0*thetaphi1253*phi_N253*(3.d0+  &
                    2.d0*phi_N253**2)+9.d0*theta1253*(4.d0-  &
                    14.000000000000002d0*phi_N253**2+phi_N253**4))+  &
                    12.d0*H253**6*theta1253**2*phi_N253**2*(thetaphi1253*phi_N253*(18.d0  &
                    -15.d0*phi_N253**2+phi_N253**4)+6.d0*theta1253*(-72.d0+  &
                    27.d0*phi_N253**2+phi_N253**4))-  &
                    1.d0*H253**4*theta1253*(-6.d0+  &
                    phi_N253**2)*(8.d0*thetaphi1253*phi_N253*(3.d0+  &
                    5.d0*phi_N253**2)+theta1253*(36.d0-504.d0*phi_N253**2+  &
                    13.d0*phi_N253**4))+  &
                    486.d0*H253**10*theta1253**4*phi_N253**6*(-1.d0*thetaphi1253*phi_N253  &
                    +9.d0*theta1253))+  &
                    phi_N253*(-16.d0*V_phi_phi1253+4.d0*H253**2*(36.d0+  &
                    3.d0*phi_N253**4-24.d0*theta1253*V_phi_phi1253+  &
                    2.d0*phi_N253**2*(-12.d0+  &
                    5.d0*theta1253*V_phi_phi1253))+  &
                    729.d0*H253**14*theta1253**4*phi_N253**6*(2.d0*thetaphi1253**2*phi_N253**2  &
                    +72.d0*theta1253**2-  &
                    1.d0*theta1253*phi_N253*(18.d0*thetaphi1253+  &
                    phi_N253*thetaphiphi1253))+  &
                    243.d0*H253**12*theta1253**3*phi_N253**4*(2.d0*thetaphi1253**2*phi_N253**2*(6.d0  &
                    +phi_N253**2)+36.d0*theta1253**2*(6.d0+7.d0*phi_N253**2)  &
                    -3.d0*theta1253*phi_N253*(thetaphi1253*(18.d0+  &
                    17.d0*phi_N253**2)+2.d0*phi_N253*thetaphiphi1253))-  &
                    2.d0*H253**4*(3.d0*theta1253*(-216.d0+420.d0*phi_N253**2-  &
                    70.d0*phi_N253**4+phi_N253**6)+2.d0*theta1253**2*(36.d0-  &
                    12.d0*phi_N253**2+7.d0*phi_N253**4)*V_phi_phi1253+  &
                    2.d0*phi_N253*(thetaphi1253*(-54.d0-51.d0*phi_N253**2+  &
                    4.d0*phi_N253**4)+2.d0*phi_N253*(3.d0+  &
                    phi_N253**2)*thetaphiphi1253))+  &
                    27.d0*H253**10*theta1253**2*phi_N253**2*(36.d0*thetaphi1253**2*phi_N253**2*(2.d0  &
                    +phi_N253**2)+54.d0*theta1253**2*(-36.d0+  &
                    40.d0*phi_N253**2+13.d0*phi_N253**4)+  &
                    18.d0*theta1253**3*phi_N253**4*V_phi_phi1253+  &
                    theta1253*phi_N253**2*(-72.d0*thetaphi1253*phi_N253*(7.d0  &
                    +phi_N253**2)+(-36.d0-12.d0*phi_N253**2+  &
                    phi_N253**4)*thetaphiphi1253))+  &
                    2.d0*H253**6*(2.d0*thetaphi1253**2*phi_N253**2*(36.d0+  &
                    6.d0*phi_N253**2+phi_N253**4)-3.d0*theta1253**2*(-648.d0  &
                    +3060.d0*phi_N253**2-894.d0*phi_N253**4+7.d0*phi_N253**6)+  &
                    3.d0*theta1253**3*phi_N253**2*(-36.d0+36.d0*phi_N253**2+  &
                    phi_N253**4)*V_phi_phi1253+  &
                    2.d0*theta1253*phi_N253*(thetaphi1253*(324.d0+  &
                    450.d0*phi_N253**2-159.d0*phi_N253**4+2.d0*phi_N253**6)+  &
                    phi_N253*(-36.d0-15.d0*phi_N253**2+  &
                    2.d0*phi_N253**4)*thetaphiphi1253))+  &
                    H253**8*theta1253*(2.d0*thetaphi1253**2*phi_N253**2*(216.d0 &
                    +324.d0*phi_N253**2+90.d0*phi_N253**4-1.d0*phi_N253**6)+  &
                    18.d0*theta1253**2*(216.d0-  &
                    2916.0000000000005d0*phi_N253**2+  &
                    1458.0000000000002d0*phi_N253**4+77.d0*phi_N253**6)-  &
                    108.d0*theta1253**3*phi_N253**4*(-3.d0+  &
                    phi_N253**2)*V_phi_phi1253+  &
                    theta1253*phi_N253*(27.d0*thetaphi1253*(72.d0+  &
                    132.d0*phi_N253**2-190.d0*phi_N253**4+phi_N253**6)-  &
                    2.d0*phi_N253*(108.d0+216.d0*phi_N253**2+9.d0*phi_N253**4+  &
                    phi_N253**6)*thetaphiphi1253)))))/(H253**2*(2.d0+  &
                    9.d0*H253**4*theta1253**2*phi_N253**2-  &
                    1.d0*H253**2*theta1253*(-6.d0+phi_N253**2))**3)
		return
		end function

!i 1st derivative of z wrt efolds
		real(dl) function Z_n(coff,phi254,phi_n254,H254,a254)
					use mod1
			use mod2
		implicit none
		real(dl) phi_n254,a254,H254,phi254,z_n_a254
		real(dl) coff(25)
        real(dl) h_n1254,theta1254,v_phi1254,thetaphi1254,thetaphiphi1254
		real(dl) v_phi_phi1254,V_phi_phi_phi1254,thetaphiphiphi1254,phi_nn1254

		!H_n1254=H_n(coff,phi254,H254,phi_n254)
		theta1254=theta(coff,phi254)
		v_phi1254=V_phi(coff,phi254)
		thetaphi1254=thetaphi(coff,phi254)
		!phi_nn1254=phi_nn(coff,phi254,H254,phi_n254)
		!thetaphiphi1254=thetaphiphi(coff,phi254)
		!v_phi_phi1254=V_phi_phi(coff,phi254)
		!h_nn1=h_nn(coff,phi,phi_n,H)
		!phi_nnn1=phi_nnn(coff,phi,phi_n,H)
		!V_phi_phi_phi1254=V_phi_phi_phi(coff,phi254)
		!thetaphiphiphi1254=thetaphiphiphi(coff,phi254)
		!!!!!!!!!!!!!!!!!!!!!!here Z_n =Z'/Z to be used in mukhanov eqn !!!!!!!!!!!!
		
		Z_n= (v_phi1254*(-16.0+27.0*H254**10*theta1254**5*phi_n254**8-  &
		9.0*H254**8*theta1254**4*phi_n254**6*(12.0+phi_n254**2)+  &
		16.0*H254**2*theta1254*(-3.0+2.0*phi_n254**2)+  &
		8.0*H254**6*theta1254**3*phi_n254**4*(24.0+5.0*phi_n254**2)-  &
		8.0*H254**4*theta1254**2*phi_n254**2*(18.0+7.0*phi_n254**2))+  &
		H254**2*phi_n254*(972.0*H254**12*theta1254**6*phi_n254**8+  &
		27.0*H254**10*theta1254**5*phi_n254**6*(-84.0+phi_n254**2)+8.0*(-4.0+  &
		phi_n254**2)+3.0*H254**8*theta1254**3*phi_n254**4*(thetaphi1254*phi_n254**5+  &
		theta1254*(936.0+54.0*phi_n254**2-20.0*phi_n254**4))-  &
		4.0*H254**2*(-2.0*thetaphi1254*phi_n254**3+theta1254*(48.0-36.0*phi_n254**2  &
		+5.0*phi_n254**4))-  &
		1.0*H254**6*theta1254**2*phi_n254**2*(2.0*thetaphi1254*phi_n254**3*(24.0+  &
		5.0*phi_n254**2)+theta1254*(720.0-132.0*phi_n254**2-220.0*phi_n254**4+  &
		3.0*phi_n254**6))+2.0*H254**4*theta1254*(2.0*thetaphi1254*phi_n254**3*(24.0  &
		+phi_n254**2)+theta1254*(-144.0+60.0*phi_n254**2-132.0*phi_n254**4+  &
		7.0*phi_n254**6))))/(H254**2*phi_n254*(-2.0+  &
		H254**2*theta1254*phi_n254**2)*(-2.0+3.0*H254**2*theta1254*phi_n254**2)*(2.0+  &
		9.0*H254**4*theta1254**2*phi_n254**2-1.0*H254**2*theta1254*(-6.0+  &
		phi_n254**2))**2)
		
		
		
		
		
		
     !   if ( ieee_is_nan(z_N) ) then
                    !write(*,*) 'z_N is NaN'
     !               open (unit = 1586, file = "error.log" ,action='write',position='append')
     !               write(1586,*)  "error kind : z_N_NaN : mod3",z_N
                    
                   ! write(1586,*)x0,(coff(is),is=1,16)
     !               close(1586)
     !   endif
		return 
    end	function	
		
		!first slow roll parameter
		real(dl) function epsln2(coff,phi_n_ep,H_ep,phi_ep)
					use mod1
                    use mod2
		implicit none
		real(dl) phi_n_ep,H_ep,phi_ep
		real(dl) coff(25)
		!gterm4=( -G_inf(coff,phi) +1.d00)
		!!print **, in epsilon,gterm4,G_inf(coff,phi),coff(16),phi
		epsln2=-H_n(coff,phi_ep,H_ep,phi_n_ep)/H_ep
    !    if ( epsln2.lt.0.0d0 ) then
            !write(*,*) 'epsilon is Negetive'
    !        open (unit = 1586, file = "error.log" ,action='write',position='append')
	!		write(1586,*)  "error kind : epsilon_negetive : mod3",epsln2
			
			!write(1586,*)x0,(coff(is),is=1,16)
	!		close(1586)
    !    endif
		
		return
		end function
		
end module mod3
