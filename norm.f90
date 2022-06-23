module spectral
use mod1
use mod2
use mod3
use mod4
use mod5
implicit none
contains

!########################################################################################################
!In this subroutine for a particular potential function and for given potential parameters we
!find a suitable starting point value for scalar field φ which satisfies the fixed number of efolds
!requirement and also takes care that the system always remains in inflationary phase throughout
!the whole evolution.
!Starting from a given initial φ using background equations system is evolved.Here we assume
!that initially it starts from slow roll.In each step the first Hubble slow role parameter (epsilon ) is
!evaluated and checked if it reaches the end of inflation condition, whenever the condition is
!achieved,we exit the program and obtain number of efolds if the number of efolds is less than
!70, we run the process again taking a largely increased value of φ. Alternatively if  epsilon stays
!within the inflationary range, we evolve the system upto 70 efolds and check the value of  H .The
!process is continued with adaptively decreasing initial φ values , until at 70 efolds we get the
!value of epsilon inside the error bar, fixed previously. 
!##########################################################################################################

   subroutine norm(phi000,nef,coff)

		implicit none

		real(dl) :: a0,H0,phi0,phi_n0,temp,odef,ndef
		real(dl)  :: phi00,phi_n00,H00,N,hn,t,ee,a,ddf,zz,k,nef,phi000,dk
		real(dl) , dimension(3)::r
		real(dl) coff(25)
		character(3) rff
		k=.1
		a0=a00
		phi00=phi_ini
		phi_n00=phi_n000
		H00=dsqrt(V(coff,phi00)/(3.0-(phi_n00**2)/2.d0))
		r=[phi00,phi_n00,H00]
		hn=1E-3
		n=400
		ee=epsln2(phi_n00,H00)
		dk=.01
		temp=100
		!print *,'n,r,ee', N,r,ee
		do while (abs(n-nef).ge. 15)
			!print *,'n,r,ee', N,phi00,ee

		      	


67			if((n.le.1) .or. (dk.le.1e-10))then
				Print *,"we have hit a wall:to proceed say 'yes'"
				read (*,*)rff
				if (rff=="yes")then
					print *,"fine then"
				else
					STOP
				endif;endif	
			odef=nef-temp
			ndef=nef-n
			temp=n
			if(abs(odef+ndef).eq.(abs(odef)+abs(ndef))) then
					print *, "ok"
			else
				dk=dk/2.d0
			endif
			if (n.le. nef) then
				
				if(coff(5).ge.0.0)then
					
					coff(5)=coff(5)+dk
				else
					coff(5)=coff(5)-dk
				end if
				 print *,"le1",n,coff(5),dk

			else 
				if(coff(5).ge.0.0)then
					
					coff(5)=coff(5)-dk
				else
					
					coff(5)=coff(5)+dk
				end if	
				print *,"ge1",n,coff(5),dk	
			endif
			go to 56

		
			
			
			
56			phi_n00=1E-2 
 
			H00=dsqrt(V(coff,phi00)/(3.0-(phi_n00**2)/2.d0))
			r=[phi00,phi_n00,H00]
			N=0.0
			ee=epsln2(phi_n00,H00)
			do while (ee.le.1.00)
			!print *,'n1111,r,ee', N,r,ee
				call odesolve(r,N,hn,coff)
				phi0=r(1)
				phi_n0=r(2)
				H0=r(3)
				ee=epsln2(phi_n0,H0)
				n=n+hn
				if (n.ge.200) then
				     
				     Go to 67
				end if     
				!print *,'n,phi0,ee', N,phi0,ee,coff(5)
71			enddo
	
		enddo
		print *, "final E-fold no",n
		phi000=phi00
      return
      end

!######################################################################################################
!This subroutine is used to evolve the system in sub-horizon limit, for a given initial value of φ
!obtained from find phi ini subroutine.
!Here also we start with slow roll condition,evolve the system with efolds using back ground
!equations and in each step the horizon limit is checked with ddf parameter. Whenever the
!sub-horizon approximation wears off, we exit the program. Passing all the variables
!#######################################################################################################


subroutine subhorizon_evolution(k, phi0,phi_n0,H0,N,hn,a0,coff)
		
		implicit none

		real(dl) :: a0,H0,phi0,phi_n0,pi,k,hn,ee,N,a,ddf,lddf
		real(dl) , dimension(3)::r
		real(dl):: coff(25)
		r=[phi0,phi_n0,H0]
		ee=epsln2(phi_n0,H0)
		a=a0*dexp(N)
		ddf=k/(a*H0)
		! print *,"ddfs",ddf,k,a,H0
		lddf=10
		do while (ddf.ge.lddf)

			! write (11,*) N,r,ee
		

			call odesolve(r,N,hn,coff)
			N=N+hn
			if(abs(ddf-lddf).le.0.05)then
			print *,"am in stage one"
			endif
			phi0=r(1)
			phi_n0=r(2)
			H0=r(3)
			ee=epsln2(phi_n0,H0)
			a=a0*dexp(N-hn)
			ddf=k/(a*H0)
			!print *,"ddf",ddf,k,a,H0
		enddo
return 
end



!###################################################################################################
!This subroutine is used to implement Bunch Davies vacuum condition , as a starting point of
!perturbative evolution.
!Whenever the sub-horizon limit is crossed, the curvature perturbation is initialized both for
!scalar and tensor part , with 4 additional complex variables (q,v and their first derivatives).
!####################################################################################################
subroutine perturbation_initiation(qr0,qc0,q_nr0,q_nc0,vr0,vc0,v_nr0,v_nc0,k,a,H0)
		
                        real(dl) qr0,qc0,q_nr0,q_nc0,vr0,vc0,v_nr0,v_nc0,k,a,H0
			!PI=4.D0*DATAN(1.D0)
			qr0=1/dsqrt(2*k)
			qc0=0.0
			q_nc0=-1*k/dsqrt(2*k)/(a*H0)
			q_nr0=0.0
			vr0=1/dsqrt(2*k)  
			v_nc0=-1*k/dsqrt(2*k)/(a*H0)
			vc0=0.0
	                v_nr0=0.0
			print *,"am in stage two"
return
end



!#####################################################################################################
!This is the computational bottleneck of the program, for the required accuracy we need to take
!very fine efolds binning but, This program essentially solve 11 coupled differential equations and
!takes long time for small efold-steps. After mode initiation, each of these variables being complex
!, their real and imaginary part is evolved separately.On the top of it the back ground is evolved.
!The error is compared with the function values and whenever error increases significantly warning
!messeges are print out. Here also in each step the perturbation length scale is compared with
!Horizon scale, whenever the super horizon limit is achieved, we stop the program and pass the
!q,v values as a function of efolds along with z also as a function of efolds.
!#######################################################################################################



subroutine perturbation_evolution(phi0,phi_n0,H0,qr0,qc0,q_nr0,q_nc0,vr0,vc0,v_nr0,v_nc0,N,a0,k,hs,coff)

	implicit none
	real(dl) :: a0,H0,phi0,phi_n0,pi,hs,uddf
	real(dl)  :: q_nr0,q_nc0,vr0,vc0,v_nr0,v_nc0,k,qr0,qc0,N,hn,t,ee,a,ddf,zz
	real(dl)  , dimension(11)::rs
	real(dl):: coff(25)
		uddf=1e-6
                rs=[phi0,phi_n0,H0,qr0,qc0,q_nr0,q_nc0,vr0,vc0,v_nr0,v_nc0]
		a=a0*dexp(N)
		ddf=k/(a*H0)
		!hn=1e-4
		do while (ddf.ge.uddf)	
			if(abs(ddf-uddf).le.1e-9)then
			print *,"am in stage three"
			endif
			call odesolve2(rs,N,hs,k,coff)
			N=N+hs
		
			phi0=rs(1)
			phi_n0=rs(2)
			H0=rs(3)
			ee=epsln2(phi_n0,H0)
			if (ee.ge.1.0) then
				     STOP
			endif
			a=a0*dexp(N)
			ddf=k/(a*H0)
			qr0=rs(4)
			qc0=rs(5)
			q_nr0=rs(6)
			q_nc0=rs(7)
			vr0=rs(8)
			vc0=rs(9)
			v_nr0=rs(10)
			v_nc0=rs(11)
			zz=zzzz(coff,phi_n0,a,H0)



		enddo
return
end

end module spectral
