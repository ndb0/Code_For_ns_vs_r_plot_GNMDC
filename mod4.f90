

module mod4
!use precision
use mod1
use mod2
use mod3
use mod3_5
implicit none
contains


!2nd derivative of q-k wrt efolds 
!!!!!!!!!!!!!!!!!!!!!!!!! here q is representing R.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real(dl) function qnn(coff,phi0_q2,H0_q2,phi_n0_q2,a0_q2,k_q2,q_n0_q2,q0_q2)
		use mod1
		use mod2
		use mod3
		use mod3_5
		implicit none
		real(dl) phi0_q2,H0_q2,phi_n0_q2,a0_q2,k_q2,q_n0_q2,q0_q2,H_n0_q2,cs1_q2
		real(dl) coff(25),z1_q2, term1_q2, term2_q2,term0_q2
		
        real(dl) h_n1_q2,z_n1_q2,z_nn1_q2,cs_n1_q2
		!real(dl) phi_nnn1,h_nn1,v_phi_phi1,V_phi_phi_phi1,thetaphiphiphi1
		!real(dl) aHbyCs_N,d2z_deta2,d2Q_deta2
       ! external aHbyCs_N,d2z_deta2,d2Q_deta2
        
		H_n1_q2=H_n(coff,phi0_q2,H0_q2,phi_n0_q2)
		!theta1=theta(coff,phi0)
		!v_phi1=V_phi(coff,phi0)
	!	thetaphi1=thetaphi(coff,phi)
		!phi_nn1=phi_nn(coff,phi0,phi_n0,H0)
		cs1_q2=cs(coff,phi0_q2,phi_n0_q2,H0_q2)
		cs_n1_q2=cs_n(coff,phi0_q2,phi_n0_q2,H0_q2)
		z_n1_q2=Z_n(coff,phi0_q2,phi_n0_q2,H0_q2,a0_q2)
		! z_nn1_q2=z_NN(coff,phi0_q2,phi_n0_q2,H0_q2,a0_q2)  !!  not required in eqn for R
		!aHbyCs_N1=aHbyCs_N(coff,phi0,phi_n0,H0,a0)
		!d2Q_deta21=d2Q_deta2(coff,phi0,phi_n0,H0,a0,k,q0)
		!d2z_deta21=d2z_deta2(coff,phi0,phi_n0,H0,a0)
		!z1_q2=zzzz(coff,phi_n0_q2,a0_q2,H0_q2,phi0_q2)
		
		
		!term0_q2=1.d0+H_n1_q2/H0_q2+cs1_q2*cs_n1_q2
		!term1_q2=((k_q2*cs1_q2/(a0_q2*H0_q2))**2-z_nn1_q2/z1_q2-z_n1_q2/z1_q2*term0_q2)*q0_q2
		!term2_q2=q_n0_q2*term0_q2
		!Print *, "H_n,H0,a0,cs1,cs_n1",H_n1,H0,a0,cs1,cs_n1
		!Print *, "***********",term0,term1,term2
		qnn=-(1.0+H_n1_q2/H0_q2+2.0*z_n1_q2)*q_n0_q2-((k_q2*cs1_q2/(a0_q2*H0_q2))**2)*q0_q2
		!Print *, "***********",term0,term1,term2,q0
		!qnn=-(cs1/(a0*H0))*q_n0*aHbyCs_N1+((cs1/(a0*H0))**2)*d2Q_deta21
		! !print **,qnn,qnn,bbbb,eeee,Z_nn(phi0,phi_n0,H0,a0),a0	
	return
	end function
	

!2nd derivative of v-k wrt efolds (for the tensor modes :: v''+{1+2(zt'/zt)+H'/H}v'+((ct*k)/(aH))^2*v=0 )
	real(dl) function vnn(coff,phi0_v2,H0_v2,phi_n0_v2,a0_v2,k_v2,v_n0_v2,v0_v2)
		use mod1
		use mod2
		use mod3
		implicit none
		real(dl) phi0_v2,H0_v2,phi_n0_v2,a0_v2,k_v2,v_n0_v2,v0_v2,H_n0_v2,ct_v2,z1tbyz_v2
		real(dl) coff(25)	
		ct_v2=cst(coff,phi0_v2,phi_n0_v2,H0_v2,a0_v2)
        H_n0_v2=H_n(coff,phi0_v2,H0_v2,phi_n0_v2)
        z1tbyz_v2=z1tbyzt(coff,phi0_v2,phi_n0_v2,H0_v2,a0_v2)
        vnn=-(1.0+H_n0_v2/H0_v2+2.0*z1tbyz_v2)*v_n0_v2-((k_v2*ct_v2/(a0_v2*H0_v2))**2)*v0_v2
    !!vnn=-(1.d0+H_n0_v2/H0_v2)*v_n0_v2-v0_v2*((k_v2/(a0_v2*H0_v2))**2)+(2.d0+H_n0_v2/H0_v2)*v0_v2
		!!print **,vnn,vnn,phi0,H0,phi_n0,a0,k,v_n0,v0	
	return
	end function

	subroutine qvnn(coff,rs,k,pain,a,qnnr,qnnc,vnnr,vnnc)
	use mod1
	use mod2
	use mod3
	implicit none
	
	real(dl) rs(11),coff(25)
	real(dl) qnnr,qnnc,vnnr,vnnc
	
	real(dl) phi,H,phi_n
	real(dl) qr,qc,q_nr,q_nc,vr,vc,v_nr,v_nc
	real(dl) cs1,cs_n1,z_n1,h_n1, z1tbyzt1,ct,k,a
	
	integer pain
	
    phi = rs(1)
    phi_n= rs(2)
    H=rs(3)
	
    qr=rs(4)
    qc=rs(5)
    q_nr=rs(6)
    q_nc=rs(7)
    vr=rs(8)
    vc=rs(9)
    v_nr=rs(10)
    v_nc=rs(11)
    
    !a=coff(8)
    
    H_n1=H_n(coff,phi,H,phi_n)

    cs1=css(pain) !cs(coff,phi,phi_n,H)
    cs_n1=cssn(pain)!cs_n(coff,phi,phi_n,H)
    z_n1=zns(pain)!Z_n(coff,phi,phi_n,H,a)
    
    ct=ctt(pain)!cst(coff,phi,phi_n,H,a)
 
    z1tbyzt1=znt(pain)!z1tbyzt(coff,phi,phi_n,H,a)
    
    
    qnnr=-(1.0+H_n1/H+2.0*z_n1)*q_nr-((k*cs1/(a*H))**2)*qr
    qnnc=-(1.0+H_n1/H+2.0*z_n1)*q_nc-((k*cs1/(a*H))**2)*qc
    vnnr=-(1.0+H_n1/H+2.0*z1tbyzt1)*v_nr-((k*ct/(a*H))**2)*vr
    vnnc=-(1.0+H_n1/H+2.0*z1tbyzt1)*v_nc-((k*ct/(a*H))**2)*vc
    
	
	end subroutine qvnn
	
	subroutine assignvals(coff,r,a,css1,cssn1,ctt1,zns1,znt1)
    use mod1
	use mod2
	use mod3
	implicit none
	
	real(dl) r(3),coff(25),h_n1
	real(dl) css1,cssn1,ctt1,zns1,znt1
		real(dl) phi,H,phi_n,a
		
    phi = r(1)
    phi_n= r(2)
    H=r(3)
	

    

    
    H_n1=H_n(coff,phi,H,phi_n)

    css1=cs(coff,phi,phi_n,H)
    cssn1=cs_n(coff,phi,phi_n,H)
    zns1=Z_n(coff,phi,phi_n,H,a)
    
    ctt1=cst(coff,phi,phi_n,H,a)
 
    znt1=z1tbyzt(coff,phi,phi_n,H,a)
	
	
	
	end subroutine assignvals
	

	real(dl) function eta(coff,phi,phi_n,H)
		use mod1
		use mod2
		use mod3
		implicit none
		real(dl)   y4,phi,H,phi_n,phi_nn1,H_n1
		real(dl) coff(25)
		H_n1=H_n(coff,phi,H,phi_n)
		phi_nn1=phi_nn(coff,phi,H,phi_n)
		eta=-(phi_nn1/phi_n+H_n1/H)
	return
	end function

subroutine set_initial(x0,coff,phi0,phi_n0,H0)
use mod1
use mod3
implicit none


		real(dl) coff(25),info(40),eff_efold,x0,phi0,phi_n0,H0, mini,phi1,phimin,thet,factor
		integer i,k
		character(90) filename
		real(dl) ep,correct_dhv,phic
		
      correct_dhv=coff(17)
       
       phi0=coff(15)
       phic=coff(1)*coff(15)
   587  H0=sqrt(V(coff,phi0)/3.0*(1.0+correct_dhv))
   
 !  go to 245
   	! setting gnmdc function
		phi1=phi0
		mini=0.0
		do while(phi1.lt.phi0*coff(1)*1.50)
			thet=theta_only_sin(coff,phi1) 
			if (mini.gt.thet)then
				mini=thet
				phimin=phi1
 !               print *, phi1,phimin
			end if
			phi1=phi1*(1.001)
		end do	
		
       factor= (-coff(16)/mini)*(phimin/coff(15))**coff(25)
       !factor=(coff(16))*((phimin/phic)**coff(25))*(-1.0/mini)
		coff(18)=factor
!		print*, "phimin,min,coff(18),phi1=",phimin,mini, coff(18),phi1
		
   	
   	
   !	 setting gnmdc function
	
		
        
       
        
    245    phi_n0=-Sqrt(2.0*(3.0*H0**2-V(coff,phi0))/(H0**2  &
		             +9.0*(H0**4)*theta(coff,phi0)))
		ep=epsln2(coff,phi_n0,H0,phi0) 
		
        if((1.0+9.0*(H0**2)*theta(coff,phi0)).lt.0.0)then
              
               print *, "phi_N imaginary: necessary to change theta amplitude"
                !print *, "Changing from coff(16)=",coff(16)
                 print *, "Changing from coff(18)=",coff(18)
                 !coff(16)=coff(16)*0.99
                  coff(18)=coff(18)*0.99
                  OPEN(unit=1522,file="../NidraB",action='write',position='append')
                  !write(*,*) "phi_N imaginary: set_ini: coff(16)=",coff(16)
                  write(*,*) "phi_N imaginary: set_ini: coff(16)=",coff(18)
                  close(1522)
              print *, "to coff(16)=",coff(16)
              print *, "to coff(18)=",coff(18)
                go to 587
         endif
        if(ep.ge.0.99)then
               
                print *, "ep> 1: necessary to change slow roll proximity"
              print *, "Changing from correct_dhv =",correct_dhv
                 correct_dhv=correct_dhv*0.9
                
                go to 587
         endif 
		 
		 coff(17)=correct_dhv
		! coff(18)=ep
		
	!	print *, "to correct_dhv=",correct_dhv,ep,phi_n0,H0,phi0,V(coff,phi0),theta(coff,phi0)
		!stop
		!print*, "gnmdc term fixed" 
		! Print *, "+++++So we are initiating with epsilon,phi_n=",ep,phi_n0,"++++++++"
      return
      end subroutine 
      
      
      
      subroutine touch_initial(x0,coff,phi0,phi_n0,H0)
use mod1
use mod3
implicit none


		real(dl) coff(25),info(40),eff_efold,x0,phi0,phi_n0,H0
		integer i
		character(90) filename
		real(dl) ep,correct_dhv
		
		correct_dhv=coff(17)
        phi0=coff(15)
        H0=sqrt(V(coff,phi0)/3.0*(1.0+correct_dhv))
        
     !   if((1.0+9.0*(H0**2)*theta(coff,phi0)).lt.0.0)then
               
     !       print *, "Error: phi_N imaginary: necessary to change theta amplitude"
     !    endif       
        
        phi_n0=-Sqrt(2.0*(3.0*H0**2-V(coff,phi0))/(H0**2  &
		             +9.0*(H0**4)*theta(coff,phi0)))
		ep=epsln2(coff,phi_n0,H0,phi0)             
      !  if(ep.ge.0.99)then
               
     !           print *, "Error: ep> 1: necessary to change slow roll proximity"
   
     !    endif 
		 
		 
      return
      end subroutine 
      
      


end module mod4

	

	
real(dl) function aHbyCs_N(coff,phi,phi_n,H,a)
		use mod2
		use mod1
		use mod3
		use mod3_5
		implicit none
		real(dl) phi_n,a,H,phi,cs1,a1,n,cs_n1
		real(dl) coff(25)

        real(dl) h_n1,theta1,v_phi1,thetaphi1,phi_nn1,thetaphiphi1
		real(dl) phi_nnn1,h_nn1,v_phi_phi1,V_phi_phi_phi1,thetaphiphiphi1

		H_n1=H_n(coff,phi,H,phi_n)
		theta1=theta(coff,phi)
		v_phi1=V_phi(coff,phi)
	!	thetaphi1=thetaphi(coff,phi)
		phi_nn1=phi_nn(coff,phi,H,phi_n)
		cs1=cs(coff,phi,phi_n,H)
		!a1=coff(8)*exp(n)
		cs_n1=cs_n(coff,phi,phi_n,H)
		aHbyCs_N=a*(H_n1/cs1+H*cs_n1+H/cs1)
		
		return
		end function
		
    real(dl) function d2z_deta2(coff,phi,phi_n,H,a)
		use mod2
		use mod1
		use mod3
		use mod3_5
		implicit none
		real(dl) phi_n,a,H,phi,cs1,a1,n,z_n1,z_nn1,aHbyCs_N
		real(dl) coff(25)

        real(dl) h_n1,theta1,v_phi1,thetaphi1,phi_nn1,thetaphiphi1
		real(dl) phi_nnn1,h_nn1,v_phi_phi1,V_phi_phi_phi1,thetaphiphiphi1

		H_n1=H_n(coff,phi,H,phi_n)
		theta1=theta(coff,phi)
		v_phi1=V_phi(coff,phi)
	!	thetaphi1=thetaphi(coff,phi)
		phi_nn1=phi_nn(coff,phi,H,phi_n)
		cs1=cs(coff,phi,phi_n,H)
		z_n1=Z_n(coff,phi,phi_n,H,a)
		z_nn1=z_NN(coff,phi,phi_n,H,a)
		!a1=coff(8)*exp(n)
		d2z_deta2=((a*H/cs1)**2)*z_nn1+(a*H/cs1)*Z_n1*aHbyCs_N(coff,phi,phi_n,H,a)
		
		return
		end function
		
        real(dl) function d2Q_deta2(coff,phi,phi_n,H,a,k,q)
		use mod2
		use mod1
		use mod3
		use mod3_5
		implicit none
		real(dl) phi_n,a,H,phi,cs1,d2z_deta2,a1,n,zrat,q,k
		real(dl) coff(25)

        real(dl) h_n1,theta1,v_phi1,thetaphi1,phi_nn1,thetaphiphi1
		real(dl) phi_nnn1,h_nn1,v_phi_phi1,V_phi_phi_phi1,thetaphiphiphi1
        
        zrat=d2z_deta2(coff,phi,phi_n,H,a)/zzzz(coff,phi_n,a,H,phi)
		
		d2Q_deta2=(k**2-zrat)*q
		
        return
        end function

