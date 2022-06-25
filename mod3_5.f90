module mod3_5
!use precision
use mod1
use mod2
use mod3
!use ieee_arithmetic
implicit none
contains
!d(aH/Cs)/dN

        
        
	real(dl) function cs(coff,phi_cs,phi_n_cs,H_cs)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n_cs,H_cs,phi_cs,thetaphi1_cs,cs2_cs
		real(dl) coff(25)

        real(dl) h_n1_cs,theta1_cs
		real(dl) phi_nn1_cs
        integer is
		H_n1_cs=H_n(coff,phi_cs,H_cs,phi_n_cs)
		theta1_cs=theta(coff,phi_cs)
		!v_phi1=V_phi(coff,phi)
		thetaphi1_cs=thetaphi(coff,phi_cs)
		phi_nn1_cs=phi_nn(coff,phi_cs,H_cs,phi_n_cs)

		
		
		cs2_cs=((h_n1_cs*(0.8888888888888888d0-  &
            3.111111111111111d0*h_cs**2*theta1_cs*phi_N_cs**2+  &
            0.6666666666666666d0*h_cs**4*theta1_cs**2*phi_N_cs**4+  &
            0.3333333333333333d0*h_cs**6*theta1_cs**3*phi_N_cs**6)+  &
            h_cs**3*phi_N_cs*(thetaphi1_cs*(-0.4444444444444444d0*phi_N_cs**2-  &
            0.4444444444444444d0*h_cs**2*theta1_cs*phi_N_cs**4+  &
            0.3333333333333333d0*h_cs**4*theta1_cs**2*phi_N_cs**6)+  &
            theta1_cs*(h_cs**2*theta1_cs*phi_N_cs**2*(-0.8888888888888888d0*phi_N_cs  &
            -0.8888888888888888d0*phi_NN1_cs)+  &
            h_cs**4*theta1_cs**2*phi_N_cs**4*(1.3333333333333333d0*phi_N_cs+  &
            0.6666666666666666d0*phi_NN1_cs)-  &
            0.8888888888888888d0*phi_NN1_cs)))/(h_cs*phi_N_cs**2*(-2.d0+  &
            h_cs**2*theta1_cs*phi_N_cs**2)*(0.2222222222222222d0+  &
            1.d0*h_cs**4*theta1_cs**2*phi_N_cs**2+  &
            h_cs**2*theta1_cs*(0.6666666666666666d0-  &
            0.1111111111111111d0*phi_N_cs**2))))
            if(cs2_cs.ge.0.d0)then
            cs=sqrt(cs2_cs)
            else
              !  !print **, "error kind 18 : main.f rates line 591"
                error_index=18
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
             !   write(1586,*)  "error kind 18: mod3_5 line no 55"
                write(1586,*)"cs^2 negetive:",cs2_cs
            !    write(1586,*)(coff(is),is=1,16)
             write(1586,*) "!######################################################!"
                close(1586)
            !cs=0.0
           ! !print **, "*****cs2 negetive error ******",cs
            endif
      return
      end function


        real(dl) function cs_new_sqr(coff,phi11,phiN11,H11,a11)
                use mod2
                use mod1
                use mod3
                implicit none
                real(dl) phiN11,H11,phi11,thetaphi1_cs,cs2_cs,thetN11
                real(dl) coff(25),a11

        real(dl) hN11,thet11
        real(dl) phiNN11
        integer is
                HN11=H_n(coff,phi11,H11,phiN11)
                thet11=theta(coff,phi11)
                thetN11=thetaphi(coff,phi11)*phiN11
                !v_phi1=V_phi(coff,phi)
                thetaphi1_cs=thetaphi(coff,phi11)
                phiNN11=phi_nn(coff,phi11,H11,phiN11)

    cs_new_sqr =(HN11*(8.d0+H11**2*phiN11**2*thet11*(-28.000000000000004d0+  &
        3.d0*H11**2*phiN11**2*thet11*(2.d0+  &
       H11**2*phiN11**2*thet11)))+  &
       H11**3*phiN11*(-8.d0*phiNN11*thet11-  &
       8.d0*H11**2*phiN11**2*(phiN11+phiNN11)*thet11**2+  &
       6.d0*H11**4*phiN11**4*(2.d0*phiN11+phiNN11)*thet11**3+  &
       phiN11*(-2.d0+H11**2*phiN11**2*thet11)*(2.d0+  &
       3.d0*H11**2*phiN11**2*thet11)*thetN11))/(H11*phiN11**2*(-2.d0  &
       +H11**2*phiN11**2*thet11)*(2.d0+H11**2*thet11*(6.d0+  &
       phiN11**2*(-1.d0+9.d0*H11**2*thet11))))
     return
      end function

      real(dl) function cst_new(coff,phi,phiN,H,a)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phiN,a,H,phi
		real(dl) coff(25)

        real(dl) theta0    !h_n720,theta720,v_phi720,thetaphi720,phi_nn720,thetaphiphi720
        
        theta0=theta(coff,phi)
		!v_phi720=V_phi(coff,phi720)
		!thetaphi720=thetaphi(coff,phi720)
		
		cst_new =  Sqrt(-1.0 + 4.0/(2.0 - H**2*phiN**2*theta0))
			
			
	return
	end function

	real(dl) function cs_n(coff,phi352,phi_n352,H352)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n352,H352,phi352,termden352
		real(dl) coff(25)
        integer is
        real(dl) h_n1352,theta1352,v_phi1352,thetaphi1352,phi_nn1352,thetaphiphi1352,v_phi_phi1352
		!real(dl) phi_nnn1,h_nn1,v_phi_phi1!,V_phi_phi_phi1,thetaphiphiphi1

		H_n1352=H_n(coff,phi352,H352,phi_n352)
		theta1352=theta(coff,phi352)
		v_phi1352=V_phi(coff,phi352)
		thetaphi1352=thetaphi(coff,phi352)
		phi_nn1352=phi_nn(coff,phi352,H352,phi_n352)
		thetaphiphi1352=thetaphiphi(coff,phi352)
		v_phi_phi1352=V_phi_phi(coff,phi352)
		!h_nn1=h_nn(coff,phi,phi_n,H)
		!phi_nnn1=phi_nnn(coff,phi,phi_n,H)
		!V_phi_phi_phi1=V_phi_phi_phi(coff,phi)
		!thetaphiphiphi1=thetaphiphiphi(coff,phi)
		termden352=(-0.09876543209876543d0  &
            +4.333333333333333d0*H352**10*theta1352**5*phi_n352**6+  &
            H352**2*theta1352*(-0.5925925925925924d0+  &
            0.5925925925925926d0*theta1352*V_phi1352*phi_N352+  &
            0.24691358024691357d0*phi_n352**2)+  &
            H352**8*theta1352**3*phi_n352**4*(-0.07407407407407407d0*thetaphi1352*phi_n352**3  &
            +theta1352*(-8.666666666666666d0+  &
            0.9629629629629631d0*phi_n352**2))+  &
            H352**4*theta1352*(-0.5925925925925924d0*theta1352**2*V_phi1352*phi_n352**3  &
            -0.2962962962962963d0*thetaphi1352*phi_n352**3+  &
            theta1352*(-0.8888888888888888d0+  &
            2.4691358024691357d0*phi_n352**2-  &
            0.17283950617283947d0*phi_n352**4))+  &
            H352**6*theta1352**2*phi_n352**2*(0.14814814814814814d0*theta1352**2*V_phi1352*phi_n352**3  &
            +0.2962962962962962d0*thetaphi1352*phi_n352**3+  &
            theta1352*(5.185185185185186d0-3.0123456790123453d0*phi_n352**2+  &
            0.03703703703703703d0*phi_n352**4)))/((-2.d0+  &
            H352**2*theta1352*phi_n352**2)*(0.2222222222222222d0+  &
            1.d0*H352**4*theta1352**2*phi_n352**2+  &
            H352**2*theta1352*(0.6666666666666666d0-  &
            0.1111111111111111d0*phi_n352**2))**2)
            
            if(termden352.le.0.d0)then
            
                print *, "error kind 19 : main.f rates line 591"
                error_index=19
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
               ! write(1586,*)  "error kind 19: mod3_5 line no 117"
                write(1586,*)"(derivative of cs)^2 negetive:",termden352
               ! write(1586,*)(coff(is),is=1,16)
                write(1586,*) "!######################################################!"
                close(1586)
            !termden352=0.0!termden352
            endif

					cs_n= (1.4210854715202004e-14*H352**20*theta1352**10*phi_N352**12+  &
            H352**2*thetaphi1352*phi_N352*(-6.1679056923619804e-18-  &
            0.004877305288827921d0*phi_N352**2+  &
            0.014631915866483767d0*H352**2*thetaphi1352*phi_N352**3)+  &
            H352**16*theta1352**9*phi_N352**10*(-1.5789838572446672e-15*H352**4*thetaphi1352*phi_N352**3  &
            +H352**2*(17.3333333333333d0+5.777777777777776d0*phi_N352**2)+  &
            phi_N352*(0.44444444444444503d0*V_phi1352+  &
            0.07407407407407407d0*phi_N352*V_phi_phi1352))+  &
            H352**12*theta1352**8*phi_N352**8*(0.041152263374485604d0*V_phi1352**2*phi_N352**2  &
            +H352**4*(-112.d0-11.555555555555555d0*phi_N352**2+  &
            (3.1111111111111116-  &
            0.03703703703703706d0*thetaphi1352*V_phi1352)*phi_N352**4)+  &
            H352**2*phi_N352*(V_phi1352*(-0.09876543209876587d0+  &
            0.6995884773662548d0*phi_N352**2)+  &
            phi_N352*(-0.3456790123456789d0-  &
            0.016460905349794237d0*phi_N352**2)*V_phi_phi1352)+  &
            H352**6*(thetaphi1352*(3.1579677144893344e-15*phi_N352**3-  &
            1.4074074074074074d0*phi_N352**5)-  &
            0.037037037037037035d0*phi_N352**6*thetaphiphi1352))+  &
            H352**10*theta1352**7*phi_N352**6*(0.09259259259259259d0*H352**8*thetaphi1352**2*phi_N352**8  &
            +V_phi1352**2*phi_N352**2*(-0.1481481481481481d0+  &
            0.002743484224965707d0*phi_N352**2)+  &
            H352**4*(175.40740740740742d0-44.246913580246925d0*phi_N352**2+  &
            (-11.390946502057613d0+  &
            0.17283950617283958d0*thetaphi1352*V_phi1352)*phi_N352**4+  &
            (0.4554183813443072d0-  &
            0.05761316872427983d0*thetaphi1352*V_phi1352)*phi_N352**6)+  &
            H352**2*phi_N352*(V_phi1352*(-9.744855967078195d0-  &
            2.716049382716049d0*phi_N352**2+  &
            0.08413351623228171d0*phi_N352**4)+  &
            phi_N352*(0.32921810699588466d0+  &
            0.12071330589849107d0*phi_N352**2+  &
            0.0009144947416552353d0*phi_N352**4)*V_phi_phi1352)+  &
            H352**6*(thetaphi1352*(2.565848768022584e-15*phi_N352**3+  &
            5.777777777777775d0*phi_N352**5-  &
            0.38271604938271603d0*phi_N352**7)+  &
            phi_N352**6*(0.17283950617283944d0+  &
            0.008230452674897118d0*phi_N352**2)*thetaphiphi1352))+  &
            H352**4*theta1352**4*phi_N352*(V_phi1352**2*phi_N352*(-0.9218106995884773d0  &
            -0.04389574759945132d0*phi_N352**2)+  &
            H352**8*thetaphi1352**2*phi_N352**7*(0.49382716049382686d0-  &
            0.18106995884773663d0*phi_N352**2-  &
            0.0027434842249657045d0*phi_N352**4)+  &
            H352**4*phi_N352*(20.543209876543212d0+(-24.669410150891625d0+  &
            0.19753086419753085d0*thetaphi1352*V_phi1352)*phi_N352**2+  &
            (-0.3511659807956104d0-  &
            0.6584362139917695d0*thetaphi1352*V_phi1352)*phi_N352**4+  &
            (0.2853223593964334d0+  &
            0.040237768632830406d0*thetaphi1352*V_phi1352)*phi_N352**6+  &
            6.853228547068867e-19*phi_N352**8)+  &
            H352**2*(V_phi1352*(2.984910836762688d0-  &
            8.881572930955645d0*phi_N352**2-  &
            0.029263831732967795d0*phi_N352**4+  &
            0.017070568510897725d0*phi_N352**6)+  &
            phi_N352*(-0.2633744855967078d0+  &
            0.08779149519890256d0*phi_N352**2-  &
            0.0731595793324188d0*phi_N352**4)*V_phi_phi1352)+  &
            H352**6*(thetaphi1352*phi_N352**2*(3.947459643111668e-16+  &
            13.563786008230453d0*phi_N352**2-  &
            6.811156835848193d0*phi_N352**4+  &
            0.1438805060204237d0*phi_N352**6+  &
            0.0033531473860691934d0*phi_N352**8)+  &
            phi_N352**5*(0.19753086419753088d0-  &
            0.1316872427983539d0*phi_N352**2-  &
            0.018289894833104704d0*phi_N352**4)*thetaphiphi1352))+  &
            H352**8*theta1352**6*phi_N352**4*(V_phi1352**2*phi_N352**2*(-0.1097393689986283d0  &
            -0.02011888431641518d0*phi_N352**2)+  &
            H352**8*thetaphi1352**2*phi_N352**8*(-0.4320987654320987d0+  &
            0.0020576131687242813d0*phi_N352**2)+  &
            H352**4*(-101.92592592592594d0+119.04526748971193d0*phi_N352**2+  &
            (4.960219478737992d0-  &
            0.16460905349794233d0*thetaphi1352*V_phi1352)*phi_N352**4+  &
            (-2.2569730224051203d0+  &
            0.26886145404663914d0*thetaphi1352*V_phi1352)*phi_N352**6+  &
            (0.02011888431641518d0-  &
            0.0004572473708276187d0*thetaphi1352*V_phi1352)*phi_N352**8)+  &
            H352**2*phi_N352*(V_phi1352*(25.283950617283946d0+  &
            0.7681755829903978d0*phi_N352**2-  &
            0.4755372656607224d0*phi_N352**4+  &
            0.0009144947416552344d0*phi_N352**6)+  &
            phi_N352*(0.3950617283950615-  &
            0.30727023319615904d0*phi_N352**2-  &
            0.009144947416552352d0*phi_N352**4)*V_phi_phi1352)+  &
            H352**6*phi_N352**5*(thetaphi1352*(-3.2921810699588416+  &
            1.3580246913580238d0*phi_N352**2-  &
            0.023776863283036135d0*phi_N352**4)+  &
            phi_N352*(-0.16460905349794233d0-  &
            0.06035665294924554d0*phi_N352**2-  &
            0.00045724737082761767d0*phi_N352**4)*thetaphiphi1352))+  &
            H352**6*theta1352**5*phi_N352**2*(V_phi1352**2*phi_N352**2*(0.9218106995884771d0  &
            +0.051211705532693164d0*phi_N352**2)+  &
            H352**8*thetaphi1352**2*phi_N352**8*(0.41152263374485587d0+  &
            0.023319615912208505d0*phi_N352**2+  &
            0.0004572473708276178d0*phi_N352**4)+  &
            H352**4*(19.753086419753092d0-88.23045267489711d0*phi_N352**2+  &
            (21.24554183813442d0-  &
            0.19753086419753085d0*thetaphi1352*V_phi1352)*phi_N352**4+  &
            (3.2775491540923625d0-  &
            0.19753086419753105d0*thetaphi1352*V_phi1352)*phi_N352**6+  &
            (-0.1286389269928364d0-  &
            0.0027434842249657006d0*thetaphi1352*V_phi1352)*phi_N352**8-  &
            1.7133071367672168e-19*phi_N352**10)+  &
            H352**2*phi_N352*(V_phi1352*(-19.40192043895747d0+  &
            7.65249199817101d0*phi_N352**2+  &
            0.7852461515012954d0*phi_N352**4-  &
            0.0067062947721383895d0*phi_N352**6)+  &
            phi_N352*(-0.39506172839506176d0+  &
            0.2633744855967078d0*phi_N352**2+  &
            0.03657978966620941d0*phi_N352**4)*V_phi_phi1352)+  &
            H352**6*(thetaphi1352*phi_N352**3*(1.5789838572446672e-15-  &
            11.237311385459531d0*phi_N352**2+  &
            0.5523548239597625d0*phi_N352**4+  &
            0.07925621094345381d0*phi_N352**6-  &
            0.0004572473708276172d0*phi_N352**8)+  &
            phi_N352**6*(-0.19753086419753074d0+  &
            0.15363511659807952d0*phi_N352**2+  &
            0.004572473708276176d0*phi_N352**4)*thetaphiphi1352))+  &
            theta1352**3*(H352**2*V_phi1352**2*(0.0877914951989026d0-  &
            0.014631915866483726d0*phi_N352**2)+  &
            H352**10*thetaphi1352**2*phi_N352**6*(-0.49382716049382724d0+  &
            0.36213991769547327d0*phi_N352**2+  &
            0.0036579789666209375d0*phi_N352**4)+  &
            H352**6*phi_N352**2*(7.374485596707821d0-  &
            2.13625971650663d0*phi_N352**2-  &
            0.2341106538637403d0*phi_N352**4-  &
            6.853228547068867e-19*phi_N352**6+  &
            thetaphi1352*V_phi1352*(0.13168724279835403d0+  &
            1.0096021947873801d0*phi_N352**2-  &
            0.13900320073159575d0*phi_N352**4))+  &
            H352**4*phi_N352*(V_phi1352*(2.077732053040695d0-  &
            0.8193872885230906d0*phi_N352**2-  &
            0.014631915866483769d0*phi_N352**4)+  &
            phi_N352*(-0.17558299039780517+  &
            0.07315957933241882d0*phi_N352**2)*V_phi_phi1352)+  &
            H352**8*(thetaphi1352*(-9.868649107779167e-17*phi_N352-  &
            3.248285322359396d0*phi_N352**3+  &
            7.776863283036121d0*phi_N352**5-  &
            0.8632830361225419*phi_N352**7-  &
            0.008535284255448855d0*phi_N352**9)+  &
            phi_N352**4*(0.1316872427983539d0-  &
            0.04389574759945128d0*phi_N352**2+  &
            0.0365797896662094d0*phi_N352**4)*thetaphiphi1352))+  &
            theta1352*phi_N352*(V_phi1352*(0.009754610577655842d0-  &
            0.10242341106538635d0*H352**2*thetaphi1352*phi_N352)+  &
            H352**2*phi_N352*(0.029263831732967545d0+  &
            H352**4*thetaphi1352**2*phi_N352**2*(-0.06584362139917695d0-  &
            0.02194787379972564d0*phi_N352**2)+  &
            H352**2*phi_N352*(thetaphi1352*(-0.39018442310623375d0+  &
            0.0024386526444139553d0*phi_N352**2)+  &
            0.014631915866483767d0*phi_N352*thetaphiphi1352)))+  &
            theta1352**2*(H352**2*V_phi1352*phi_N352*(0.3901844231062337d0-  &
            0.004877305288827898d0*phi_N352**2+  &
            H352**2*thetaphi1352*phi_N352*(-0.2633744855967078d0+  &
            0.1975308641975308d0*phi_N352**2))+  &
            0.029263831732967534d0*V_phi1352**2+  &
            H352**2*phi_N352*(H352**2*(0.994970278920896d0*phi_N352+  &
            0.019509221155311715*phi_N352**3)+  &
            H352**6*thetaphi1352**2*phi_N352**3*(-0.3292181069958848d0-  &
            0.18655692729766812d0*phi_N352**2+  &
            0.007315957933241881d0*phi_N352**4)-  &
            0.029263831732967534d0*phi_N352*V_phi_phi1352+  &
            H352**4*(thetaphi1352*(-9.868649107779167e-17-  &
            2.2094192958390493d0*phi_N352**2+  &
            1.0925163846974546d0*phi_N352**4+  &
            0.0073159579332418715*phi_N352**6)+  &
            phi_N352**3*(0.0877914951989026d0-  &
            0.03657978966620941d0*phi_N352**2)*thetaphiphi1352))))/((-2.d0+  &
            H352**2*theta1352*phi_N352**2)**2*(0.2222222222222222d0+  &
            1.d0*H352**4*theta1352**2*phi_N352**2+  &
            H352**2*theta1352*(0.6666666666666666d0-  &
            0.1111111111111111d0*phi_N352**2))**4*Sqrt(termden352))
            
          !  if ( ieee_is_nan(cs_N) ) then
          !             write(*,*) 'cs_N is NaN'
          !             !print **, coff,phi352,phi_n352,H352
          !             !print **, " h_n1,theta1,v_phi1,thetaphi1,phi_nn1,thetaphiphi1,v_phi_phi1"
          !             !print **, h_n1352,theta1352,v_phi1352,thetaphi1352,phi_nn1352,thetaphiphi1352,v_phi_phi1352
                      ! read(*,*) rff
          !      endif
		
		return
		end function
		
		! speed of tensor modes.
		real(dl) function cst(coff,phi720,phi_n720,H720,a720)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n720,a720,H720,phi720
		real(dl) coff(25)

        real(dl) h_n720,theta720,v_phi720,thetaphi720,phi_nn720,thetaphiphi720
        
        theta720=theta(coff,phi720)
		v_phi720=V_phi(coff,phi720)
		thetaphi720=thetaphi(coff,phi720)
		
		
		
		!!!!! corrected definition of cst ::: 17-02-2022.
		cst = Sqrt(-1.0 + 4.0/(2.0 - H720**2*phi_N720**2*theta720))
   ! sqrt(-1.0+1.0/(0.50+  &
			!.250*H720**2*theta720*phi_n720**2))
			
			
	return
	end function
	
	! derivative of tensor speed  :::::::::::::::::::::::: corrected 22-06-22 
	real(dl) function cs1t(coff,phi721,phi_n721,H721,a721)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n721,a721,H721,phi721
		real(dl) coff(25)

        real(dl) h_n721,theta721,v_phi721,thetaphi721,phi_nn721,thetaphiphi721
        
        theta721=theta(coff,phi721)
		v_phi721=V_phi(coff,phi721)
		thetaphi721=thetaphi(coff,phi721)
		
		cs1t =(phi_n721*(H721**2*theta721*(-2.6666666666666665d0-  &
        8.d0*H721**2*theta721)*phi_n721+H721**2*(0.4444444444444444d0*thetaphi721+  &
        0.4444444444444444d0*theta721**2*v_phi721)*phi_n721**2  &
        +H721**4*theta721**2*(4.d0+12.d0*H721**2*theta721)*phi_n721**3-  &
        0.2222222222222222d0*H721**4*theta721*thetaphi721*phi_n721**4-  &
        0.8888888888888888d0*theta721*v_phi721))/((2.d0-  &
        1.d0*H721**2*theta721*phi_n721**2)**2*Sqrt(-1.d0+1/(0.5d0-  &
        0.25d0*H721**2*theta721*phi_n721**2))*(0.2222222222222222d0+  &
        H721**2*theta721*(0.6666666666666666d0+(-0.1111111111111111d0+  &
        1.d0*H721**2*theta721)*phi_n721**2)))
                
                
                
                
                
                
                
               
        
	return
	end function
	! zt for tensors. :::::::::::::::::::::::; corrected 22-06-22
	real(dl) function zt(coff,phi723,phi_n723,H723,a723) 
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n723,a723,H723,phi723
		real(dl) coff(25)

        real(dl) h_n723,theta723,v_phi723,thetaphi723,phi_nn723,thetaphiphi723
        
        	theta723=theta(coff,phi723)
		v_phi723=V_phi(coff,phi723)
		thetaphi723=thetaphi(coff,phi723)
		
		
		zt=0.353553390593273730*a723*sqrt(1.0-  &
		0.50*H723**2*theta723*phi_n723**2)
	return
	end function
	

	
	!z1t/zt for tensors i.e zt'/zt for tensors in mukhanov eqn ::::::::::::::::::::;corrected 22-06-22
	real(dl) function z1tbyzt(coff,phi722,phi_n722,H722,a722)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n722,a722,H722,phi722
		real(dl) coff(25)
		real(dl) h_n722,theta722,v_phi722,thetaphi722,phi_nn722,thetaphiphi722
		
		
		theta722=theta(coff,phi722)
		v_phi722=V_phi(coff,phi722)
		thetaphi722=thetaphi(coff,phi722)

        
		z1tbyzt=(-0.4444444444444445d0-  &
        0.22222222222222224d0*theta722*v_phi722*phi_n722+  &
        4.000000000000001d0*H722**6*theta722**3*phi_n722**4+  &
        H722**4*theta722*phi_n722**2*(-0.05555555555555556d0*thetaphi722*phi_n722**3+  &
        theta722*(-3.333333333333334d0+0.888888888888889d0*phi_n722**2))+  &
        H722**2*(-0.2222222222222223d0*theta722*phi_n722**2+  &
        (0.11111111111111112d0*thetaphi722+  &
        0.11111111111111112d0*theta722**2*v_phi722)*phi_n722**3  &
        -1.3333333333333335d0*theta722))/((-2.d0+  &
        1.d0*H722**2*theta722*phi_n722**2)*(0.2222222222222222d0+  &
        H722**2*theta722*(0.6666666666666666d0+(-0.1111111111111111d0+  &
        1.d0*H722**2*theta722)*phi_n722**2)))
			
	return
	end function
	
			
!2nd derivative of z wrt efolds
	real(dl) function Z_nn(coff,phi353,phi_n353,H353,a353)
		use mod2
		use mod1
		use mod3
		implicit none
		real(dl) phi_n353,a353,H353,phi353
		real(dl) coff(25)

        real(dl) h_n1353,theta1353,v_phi1353,thetaphi1353,phi_nn1353,thetaphiphi1353
		real(dl) phi_nnn1353,h_nn1353,v_phi_phi1353,V_phi_phi_phi1353,thetaphiphiphi1353

		h_n1353=H_n(coff,phi353,H353,phi_n353)
		theta1353=theta(coff,phi353)
		v_phi1353=V_phi(coff,phi353)
		thetaphi1353=thetaphi(coff,phi353)
		phi_nn1353=phi_nn(coff,phi353,H353,phi_n353)
		thetaphiphi1353=thetaphiphi(coff,phi353)
		V_phi_phi1353=V_phi_phi(coff,phi353)
		h_nn1353=h_nn(coff,phi353,phi_n353,H353)
		phi_nnn1353=phi_nnn(coff,phi353,phi_n353,H353)
		V_phi_phi_phi1353=V_phi_phi_phi(coff,phi353)
		thetaphiphiphi1353=thetaphiphiphi(coff,phi353)
		
	!!!!!!!!!!!!!!!!! Z_nn is not corrected for now for the new def of Z=Sqrt(2*Qs)*a	
       Z_nn=-(0.13042371500610728d0*a353*phi_N353**6*(-0.75d0*(H353**11*theta1353**5*(2.d0*h_n1353  &
        -9.d0*H353-117.d0*H353**3*theta1353)*thetaphi1353*phi_N353**12+  &
        (7.111111111111111d0+42.666666666666664d0*H353**2*theta1353+  &
        64.d0*H353**4*theta1353**2)*phi_nn1353+phi_N353*(7.111111111111111d0  &
        +21.333333333333332d0*H353*theta1353*h_n1353+  &
        64.d0*H353**3*theta1353**2*h_n1353+  &
        H353**2*theta1353*(42.666666666666664d0-  &
        53.333333333333336d0*theta1353*v_phi1353*phi_nn1353)+  &
        64.d0*H353**4*theta1353**2)+  &
        H353*theta1353*phi_N353**3*(10.666666666666666d0*h_n1353-  &
        103.11111111111111d0*H353**2*theta1353*h_n1353-  &
        512.d0*H353**4*theta1353**2*h_n1353+H353*(-32.d0-  &
        21.333333333333332d0*thetaphi1353*v_phi1353-  &
        10.666666666666666d0*theta1353*V_phi_phi1353)+  &
        H353**3*(64.d0*theta1353**2*v_phi1353*phi_nn1353+  &
        37.33333333333333d0*thetaphi1353*phi_nn1353-  &
        263.1111111111111d0*theta1353)-  &
        501.3333333333333d0*H353**5*theta1353**2)+  &
        H353**6*theta1353**2*phi_N353**8*(H353*theta1353**2*(10.666666666666668d0*H353*phi_nn1353  &
        +237.33333333333326d0*H353**3*theta1353*phi_nn1353+  &
        936.d0*H353**5*theta1353**2*phi_nn1353-  &
        13.333333333333327d0*theta1353*h_n1353*v_phi1353+  &
        53.333333333333336d0*H353*theta1353*v_phi1353)+(9.333333333333336d0  &
        +8.d0*H353*theta1353*h_n1353+0.8888888888887974d0*H353**2*theta1353-  &
        780.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**2*phi_N353**4*(H353*theta1353**2*(26.666666666666654d0*H353*phi_nn1353  &
        +462.22222222222206d0*H353**3*theta1353*phi_nn1353+  &
        1317.333333333333d0*H353**5*theta1353**2*phi_nn1353-  &
        32.d0*theta1353*h_n1353*v_phi1353+128.d0*H353*theta1353*v_phi1353)+  &
        (5.333333333333333d0+21.333333333333332d0*H353*theta1353*h_n1353-  &
        30.22222222222222d0*H353**2*theta1353-  &
        256.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**8*theta1353**3*phi_N353**10*(H353**2*theta1353**2*(-2.d0*phi_nn1353-  &
        234.d0*H353**4*theta1353**2*phi_nn1353+  &
        theta1353*(-52.d0*H353**2*phi_nn1353-8.d0*v_phi1353))+(-2.d0-  &
        6.666666666666664d0*H353*theta1353*h_n1353+  &
        16.666666666666703d0*H353**2*theta1353+  &
        468.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**4*theta1353*phi_N353**6*(H353*theta1353**2*(-22.22222222222223d0*H353*phi_nn1353  &
        -444.44444444444434d0*H353**3*theta1353*phi_nn1353-  &
        1560.d0*H353**5*theta1353**2*phi_nn1353+  &
        48.d0*theta1353*h_n1353*v_phi1353-128.d0*H353*theta1353*v_phi1353)+  &
        (-13.333333333333336d0-16.d0*H353*theta1353*h_n1353+  &
        15.999999999999945d0*H353**2*theta1353+  &
        658.6666666666665d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353*phi_N353**2*(-512.d0*H353**5*theta1353**3*phi_nn1353-  &
        21.333333333333332d0*theta1353**2*h_n1353*v_phi1353+  &
        H353*(-21.333333333333332d0*theta1353*phi_nn1353-  &
        42.666666666666664d0*theta1353**2*v_phi1353+  &
        10.666666666666666d0*thetaphi1353)+  &
        H353**3*theta1353*(-234.66666666666663d0*theta1353*phi_nn1353+  &
        32.d0*thetaphi1353))+  &
        1.d0*H353**12*theta1353**5*phi_N353**13*thetaphiphi1353+  &
        H353**10*theta1353**3*phi_N353**11*(theta1353**2*(-2.d0+  &
        thetaphi1353*(5.d0*H353**2*phi_nn1353-2.d0*v_phi1353))+  &
        H353**3*theta1353**4*(-234.d0*h_n1353-234.d0*H353)+  &
        theta1353**3*(-26.d0*H353*h_n1353-2.d0*V_phi_phi1353-52.d0*H353**2)+  &
        3.333333333333332d0*thetaphi1353**2-  &
        6.666666666666668d0*theta1353*thetaphiphi1353)+  &
        H353**5*theta1353*phi_N353**7*(18.66666666666667d0*theta1353**2*h_n1353  &
        -126.2222222222224d0*H353**2*theta1353**3*h_n1353-  &
        1560.d0*H353**4*theta1353**4*h_n1353+  &
        H353**3*theta1353**2*(26.66666666666667d0*theta1353**2*v_phi1353*phi_nn1353  &
        +56.00000000000001d0*thetaphi1353*phi_nn1353-  &
        636.4444444444443d0*theta1353)-1840.d0*H353**5*theta1353**4+  &
        H353*(theta1353**2*(-40.888888888888886d0-8.d0*thetaphi1353*v_phi1353)  &
        -32.d0*theta1353**3*V_phi_phi1353+8.d0*thetaphi1353**2-  &
        16.d0*theta1353*thetaphiphi1353))+  &
        H353**3*phi_N353**5*(-26.66666666666667d0*theta1353**2*h_n1353+  &
        159.9999999999999d0*H353**2*theta1353**3*h_n1353+  &
        1317.333333333333d0*H353**4*theta1353**4*h_n1353+  &
        1418.6666666666665d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-48.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        64.d0*thetaphi1353*phi_nn1353+604.4444444444443d0*theta1353)+  &
        H353*(theta1353**2*(53.333333333333336d0+16.d0*thetaphi1353*v_phi1353)  &
        +32.d0*theta1353**3*V_phi_phi1353+  &
        5.333333333333333d0*thetaphi1353**2+  &
        5.333333333333333d0*theta1353*thetaphiphi1353))+  &
        H353**7*theta1353**2*phi_N353**9*(-4.d0*theta1353**2*h_n1353+  &
        86.66666666666674d0*H353**2*theta1353**3*h_n1353+  &
        936.d0*H353**4*theta1353**4*h_n1353+1092.d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-6.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        26.66666666666667d0*thetaphi1353*phi_nn1353+  &
        301.3333333333333d0*theta1353)+  &
        H353*(theta1353**2*(14.666666666666666d0+  &
        6.666666666666664d0*thetaphi1353*v_phi1353)+  &
        13.333333333333334d0*theta1353**3*V_phi_phi1353-  &
        12.d0*thetaphi1353**2+  &
        16.d0*theta1353*thetaphiphi1353)))*(-234.d0*H353**13*theta1353**7*h_n1353*phi_N353**11  &
        +H353*theta1353*h_n1353*phi_N353*(21.333333333333332d0-  &
        21.333333333333332d0*theta1353*v_phi1353*phi_N353+  &
        10.666666666666664d0*phi_N353**2)+  &
        H353**11*theta1353**5*h_n1353*phi_N353**9*(2.d0*thetaphi1353*phi_N353**3  &
        +theta1353*(936.d0-26.d0*phi_N353**2))+  &
        H353**9*theta1353**4*h_n1353*phi_N353**7*(-6.6666666666666625d0*thetaphi1353*phi_N353**3  &
        +theta1353*(-1560.d0+86.66666666666666d0*phi_N353**2))+  &
        H353**3*theta1353*h_n1353*phi_N353*(-32.00000000000001d0*theta1353**2*v_phi1353*phi_N353**3  &
        +21.333333333333332d0*thetaphi1353*phi_N353**3+theta1353*(64.d0-  &
        103.1111111111111d0*phi_N353**2-  &
        26.666666666666668d0*phi_N353**4))+  &
        H353**7*theta1353**3*h_n1353*phi_N353**5*(-13.333333333333334d0*theta1353**2*v_phi1353*phi_N353**3  &
        +8.000000000000004d0*thetaphi1353*phi_N353**3+  &
        theta1353*(1317.333333333333d0-126.22222222222226d0*phi_N353**2-  &
        3.9999999999999982d0*phi_N353**4))+  &
        H353**5*theta1353**2*h_n1353*phi_N353**3*(48.00000000000001d0*theta1353**2*v_phi1353*phi_N353**3  &
        -15.99999999999999d0*thetaphi1353*phi_N353**3+theta1353*(-512.d0+  &
        160.d0*phi_N353**2+18.666666666666664d0*phi_N353**4))+  &
        7.111111111111111d0*phi_nn1353+  &
        H353**14*theta1353**6*phi_N353**10*(-117.d0*thetaphi1353*phi_N353**2-  &
        234.d0*theta1353*phi_nn1353)+  &
        H353**2*(thetaphi1353*phi_N353**2*(10.666666666666666d0-  &
        21.333333333333332d0*theta1353*v_phi1353*phi_N353+  &
        5.333333333333332d0*phi_N353**2)+  &
        theta1353*((42.666666666666664d0-  &
        21.333333333333336d0*phi_N353**2)*phi_nn1353+  &
        theta1353*(-53.333333333333336d0*v_phi1353*phi_N353*phi_nn1353-  &
        10.666666666666666d0*phi_N353**3*V_phi_phi1353)))+  &
        H353**12*theta1353**5*phi_N353**8*(-52.d0*theta1353*phi_N353**2*phi_nn1353  &
        +936.d0*theta1353*phi_nn1353+thetaphi1353*phi_N353**2*(468.d0-  &
        13.d0*phi_N353**2+5.d0*phi_N353*phi_nn1353)+  &
        1.d0*phi_N353**5*thetaphiphi1353)+  &
        H353**6*theta1353*phi_N353**2*(8.000000000000002d0*thetaphi1353**2*phi_N353**5  &
        +  &
        theta1353**2*(-8.000000000000004d0*thetaphi1353*v_phi1353*phi_N353**5  &
        +(-511.9999999999999d0+462.22222222222206d0*phi_N353**2-  &
        22.22222222222223d0*phi_N353**4)*phi_nn1353)+  &
        theta1353**3*(-47.999999999999986d0*v_phi1353*phi_N353**3*phi_nn1353  &
        -32.d0*phi_N353**5*V_phi_phi1353)+  &
        theta1353*phi_N353**2*(thetaphi1353*(-256.d0+80.d0*phi_N353**2+  &
        9.333333333333332d0*phi_N353**4-  &
        63.99999999999999d0*phi_N353*phi_nn1353)-  &
        16.d0*phi_N353**3*thetaphiphi1353))+  &
        H353**10*theta1353**3*phi_N353**6*(3.333333333333333d0*thetaphi1353**2*phi_N353**5  &
        +theta1353**2*(-2.d0*thetaphi1353*v_phi1353*phi_N353**5+(-1560.d0+  &
        237.33333333333337d0*phi_N353**2-2.d0*phi_N353**4)*phi_nn1353)+  &
        theta1353**3*(-6.d0*v_phi1353*phi_N353**3*phi_nn1353-  &
        2.d0*phi_N353**5*V_phi_phi1353)+  &
        theta1353*phi_N353**2*(thetaphi1353*(-780.d0+  &
        43.33333333333333d0*phi_N353**2-  &
        26.666666666666668d0*phi_N353*phi_nn1353)-  &
        6.666666666666666d0*phi_N353**3*thetaphiphi1353))+  &
        H353**4*(5.333333333333333d0*thetaphi1353**2*phi_N353**5+  &
        theta1353**2*(15.99999999999999d0*thetaphi1353*v_phi1353*phi_N353**5  &
        +(64.d0-234.66666666666666d0*phi_N353**2+  &
        26.666666666666668d0*phi_N353**4)*phi_nn1353)+  &
        theta1353**3*(63.99999999999999d0*v_phi1353*phi_N353**3*phi_nn1353  &
        +32.d0*phi_N353**5*V_phi_phi1353)+  &
        theta1353*phi_N353**2*(thetaphi1353*(32.d0-  &
        51.55555555555556d0*phi_N353**2-  &
        13.333333333333334d0*phi_N353**4+  &
        37.33333333333333d0*phi_N353*phi_nn1353)+  &
        5.333333333333333d0*phi_N353**3*thetaphiphi1353))+  &
        H353**8*theta1353**2*phi_N353**4*(-12.000000000000002d0*thetaphi1353**2*phi_N353**5  &
        +  &
        theta1353**2*(6.6666666666666625d0*thetaphi1353*v_phi1353*phi_N353**5  &
        +(1317.333333333333d0-444.44444444444434d0*phi_N353**2+  &
        10.666666666666668d0*phi_N353**4)*phi_nn1353)+  &
        theta1353**3*(26.666666666666668d0*v_phi1353*phi_N353**3*phi_nn1353  &
        +13.333333333333334d0*phi_N353**5*V_phi_phi1353)+  &
        theta1353*phi_N353**2*(thetaphi1353*(658.6666666666666d0-  &
        63.11111111111113d0*phi_N353**2-  &
        1.9999999999999991d0*phi_N353**4+  &
        55.999999999999986d0*phi_N353*phi_nn1353)+  &
        16.d0*phi_N353**3*thetaphiphi1353)))+1.d0*phi_N353*(2.d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)*(0.6666666666666666d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)*(1.3333333333333333d0-  &
        58.5d0*H353**10*theta1353**5*phi_N353**6+H353**2*theta1353*(8.d0-  &
        8.d0*theta1353*v_phi1353*phi_N353-3.333333333333334d0*phi_N353**2)  &
        +H353**8*theta1353**3*phi_N353**4*(1.d0*thetaphi1353*phi_N353**3+  &
        theta1353*(117.d0-13.d0*phi_N353**2))+  &
        H353**6*theta1353**2*phi_N353**2*(-2.d0*theta1353**2*v_phi1353*phi_N353**3  &
        -4.d0*thetaphi1353*phi_N353**3+theta1353*(-70.d0+  &
        40.666666666666664d0*phi_N353**2-0.5d0*phi_N353**4))+  &
        H353**4*theta1353*(8.d0*theta1353**2*v_phi1353*phi_N353**3+  &
        4.d0*thetaphi1353*phi_N353**3+theta1353*(12.d0-  &
        33.333333333333336d0*phi_N353**2+  &
        2.3333333333333335d0*phi_N353**4)))*(H353**11*theta1353**5*(2.d0*h_n1353  &
        -9.d0*H353-117.d0*H353**3*theta1353)*thetaphi1353*phi_N353**12+  &
        (7.111111111111111d0+42.666666666666664d0*H353**2*theta1353+  &
        64.d0*H353**4*theta1353**2)*phi_nn1353+phi_N353*(7.111111111111111d0  &
        +21.333333333333332d0*H353*theta1353*h_n1353+  &
        64.d0*H353**3*theta1353**2*h_n1353+  &
        H353**2*theta1353*(42.666666666666664d0-  &
        53.333333333333336d0*theta1353*v_phi1353*phi_nn1353)+  &
        64.d0*H353**4*theta1353**2)+  &
        H353*theta1353*phi_N353**3*(10.666666666666666d0*h_n1353-  &
        103.11111111111111d0*H353**2*theta1353*h_n1353-  &
        512.d0*H353**4*theta1353**2*h_n1353+H353*(-32.d0-  &
        21.333333333333332d0*thetaphi1353*v_phi1353-  &
        10.666666666666666d0*theta1353*V_phi_phi1353)+  &
        H353**3*(64.d0*theta1353**2*v_phi1353*phi_nn1353+  &
        37.33333333333333d0*thetaphi1353*phi_nn1353-  &
        263.1111111111111d0*theta1353)-  &
        501.3333333333333d0*H353**5*theta1353**2)+  &
        H353**6*theta1353**2*phi_N353**8*(H353*theta1353**2*(10.666666666666668d0*H353*phi_nn1353  &
        +237.33333333333326d0*H353**3*theta1353*phi_nn1353+  &
        936.d0*H353**5*theta1353**2*phi_nn1353-  &
        13.333333333333327d0*theta1353*h_n1353*v_phi1353+  &
        53.333333333333336d0*H353*theta1353*v_phi1353)+(9.333333333333336d0  &
        +8.d0*H353*theta1353*h_n1353+0.8888888888887974d0*H353**2*theta1353-  &
        780.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**2*phi_N353**4*(H353*theta1353**2*(26.666666666666654d0*H353*phi_nn1353  &
        +462.22222222222206d0*H353**3*theta1353*phi_nn1353+  &
        1317.333333333333d0*H353**5*theta1353**2*phi_nn1353-  &
        32.d0*theta1353*h_n1353*v_phi1353+128.d0*H353*theta1353*v_phi1353)+  &
        (5.333333333333333d0+21.333333333333332d0*H353*theta1353*h_n1353-  &
        30.22222222222222d0*H353**2*theta1353-  &
        256.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**8*theta1353**3*phi_N353**10*(H353**2*theta1353**2*(-2.d0*phi_nn1353-  &
        234.d0*H353**4*theta1353**2*phi_nn1353+  &
        theta1353*(-52.d0*H353**2*phi_nn1353-8.d0*v_phi1353))+(-2.d0-  &
        6.666666666666664d0*H353*theta1353*h_n1353+  &
        16.666666666666703d0*H353**2*theta1353+  &
        468.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**4*theta1353*phi_N353**6*(H353*theta1353**2*(-22.22222222222223d0*H353*phi_nn1353  &
        -444.44444444444434d0*H353**3*theta1353*phi_nn1353-  &
        1560.d0*H353**5*theta1353**2*phi_nn1353+  &
        48.d0*theta1353*h_n1353*v_phi1353-128.d0*H353*theta1353*v_phi1353)+  &
        (-13.333333333333336d0-16.d0*H353*theta1353*h_n1353+  &
        15.999999999999945d0*H353**2*theta1353+  &
        658.6666666666665d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353*phi_N353**2*(-512.d0*H353**5*theta1353**3*phi_nn1353-  &
        21.333333333333332d0*theta1353**2*h_n1353*v_phi1353+  &
        H353*(-21.333333333333332d0*theta1353*phi_nn1353-  &
        42.666666666666664d0*theta1353**2*v_phi1353+  &
        10.666666666666666d0*thetaphi1353)+  &
        H353**3*theta1353*(-234.66666666666663d0*theta1353*phi_nn1353+  &
        32.d0*thetaphi1353))+  &
        1.d0*H353**12*theta1353**5*phi_N353**13*thetaphiphi1353+  &
        H353**10*theta1353**3*phi_N353**11*(theta1353**2*(-2.d0+  &
        thetaphi1353*(5.d0*H353**2*phi_nn1353-2.d0*v_phi1353))+  &
        H353**3*theta1353**4*(-234.d0*h_n1353-234.d0*H353)+  &
        theta1353**3*(-26.d0*H353*h_n1353-2.d0*V_phi_phi1353-52.d0*H353**2)+  &
        3.333333333333332d0*thetaphi1353**2-  &
        6.666666666666668d0*theta1353*thetaphiphi1353)+  &
        H353**5*theta1353*phi_N353**7*(18.66666666666667d0*theta1353**2*h_n1353  &
        -126.2222222222224d0*H353**2*theta1353**3*h_n1353-  &
        1560.d0*H353**4*theta1353**4*h_n1353+  &
        H353**3*theta1353**2*(26.66666666666667d0*theta1353**2*v_phi1353*phi_nn1353  &
        +56.00000000000001d0*thetaphi1353*phi_nn1353-  &
        636.4444444444443d0*theta1353)-1840.d0*H353**5*theta1353**4+  &
        H353*(theta1353**2*(-40.888888888888886d0-8.d0*thetaphi1353*v_phi1353)  &
        -32.d0*theta1353**3*V_phi_phi1353+8.d0*thetaphi1353**2-  &
        16.d0*theta1353*thetaphiphi1353))+  &
        H353**3*phi_N353**5*(-26.66666666666667d0*theta1353**2*h_n1353+  &
        159.9999999999999d0*H353**2*theta1353**3*h_n1353+  &
        1317.333333333333d0*H353**4*theta1353**4*h_n1353+  &
        1418.6666666666665d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-48.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        64.d0*thetaphi1353*phi_nn1353+604.4444444444443d0*theta1353)+  &
        H353*(theta1353**2*(53.333333333333336d0+16.d0*thetaphi1353*v_phi1353)  &
        +32.d0*theta1353**3*V_phi_phi1353+  &
        5.333333333333333d0*thetaphi1353**2+  &
        5.333333333333333d0*theta1353*thetaphiphi1353))+  &
        H353**7*theta1353**2*phi_N353**9*(-4.d0*theta1353**2*h_n1353+  &
        86.66666666666674d0*H353**2*theta1353**3*h_n1353+  &
        936.d0*H353**4*theta1353**4*h_n1353+1092.d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-6.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        26.66666666666667d0*thetaphi1353*phi_nn1353+  &
        301.3333333333333d0*theta1353)+  &
        H353*(theta1353**2*(14.666666666666666d0+  &
        6.666666666666664d0*thetaphi1353*v_phi1353)+  &
        13.333333333333334d0*theta1353**3*V_phi_phi1353-  &
        12.d0*thetaphi1353**2+16.d0*theta1353*thetaphiphi1353)))+3.d0*(2.d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)*(0.6666666666666666d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)*(1.3333333333333333d0-  &
        58.5d0*H353**10*theta1353**5*phi_N353**6+H353**2*theta1353*(8.d0-  &
        8.d0*theta1353*v_phi1353*phi_N353-3.333333333333334d0*phi_N353**2)  &
        +H353**8*theta1353**3*phi_N353**4*(1.d0*thetaphi1353*phi_N353**3+  &
        theta1353*(117.d0-13.d0*phi_N353**2))+  &
        H353**6*theta1353**2*phi_N353**2*(-2.d0*theta1353**2*v_phi1353*phi_N353**3  &
        -4.d0*thetaphi1353*phi_N353**3+theta1353*(-70.d0+  &
        40.666666666666664d0*phi_N353**2-0.5d0*phi_N353**4))+  &
        H353**4*theta1353*(8.d0*theta1353**2*v_phi1353*phi_N353**3+  &
        4.d0*thetaphi1353*phi_N353**3+theta1353*(12.d0-  &
        33.333333333333336d0*phi_N353**2+  &
        2.3333333333333335d0*phi_N353**4)))*phi_nn1353*(H353**11*theta1353**5*(2.d0*h_n1353  &
        -9.d0*H353-117.d0*H353**3*theta1353)*thetaphi1353*phi_N353**12+  &
        (7.111111111111111d0+42.666666666666664d0*H353**2*theta1353+  &
        64.d0*H353**4*theta1353**2)*phi_nn1353+phi_N353*(7.111111111111111d0  &
        +21.333333333333332d0*H353*theta1353*h_n1353+  &
        64.d0*H353**3*theta1353**2*h_n1353+  &
        H353**2*theta1353*(42.666666666666664d0-  &
        53.333333333333336d0*theta1353*v_phi1353*phi_nn1353)+  &
        64.d0*H353**4*theta1353**2)+  &
        H353*theta1353*phi_N353**3*(10.666666666666666d0*h_n1353-  &
        103.11111111111111d0*H353**2*theta1353*h_n1353-  &
        512.d0*H353**4*theta1353**2*h_n1353+H353*(-32.d0-  &
        21.333333333333332d0*thetaphi1353*v_phi1353-  &
        10.666666666666666d0*theta1353*V_phi_phi1353)+  &
        H353**3*(64.d0*theta1353**2*v_phi1353*phi_nn1353+  &
        37.33333333333333d0*thetaphi1353*phi_nn1353-  &
        263.1111111111111d0*theta1353)-  &
        501.3333333333333d0*H353**5*theta1353**2)+  &
        H353**6*theta1353**2*phi_N353**8*(H353*theta1353**2*(10.666666666666668d0*H353*phi_nn1353  &
        +237.33333333333326d0*H353**3*theta1353*phi_nn1353+  &
        936.d0*H353**5*theta1353**2*phi_nn1353-  &
        13.333333333333327d0*theta1353*h_n1353*v_phi1353+  &
        53.333333333333336d0*H353*theta1353*v_phi1353)+(9.333333333333336d0  &
        +8.d0*H353*theta1353*h_n1353+0.8888888888887974d0*H353**2*theta1353-  &
        780.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**2*phi_N353**4*(H353*theta1353**2*(26.666666666666654d0*H353*phi_nn1353  &
        +462.22222222222206d0*H353**3*theta1353*phi_nn1353+  &
        1317.333333333333d0*H353**5*theta1353**2*phi_nn1353-  &
        32.d0*theta1353*h_n1353*v_phi1353+128.d0*H353*theta1353*v_phi1353)+  &
        (5.333333333333333d0+21.333333333333332d0*H353*theta1353*h_n1353-  &
        30.22222222222222d0*H353**2*theta1353-  &
        256.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**8*theta1353**3*phi_N353**10*(H353**2*theta1353**2*(-2.d0*phi_nn1353-  &
        234.d0*H353**4*theta1353**2*phi_nn1353+  &
        theta1353*(-52.d0*H353**2*phi_nn1353-8.d0*v_phi1353))+(-2.d0-  &
        6.666666666666664d0*H353*theta1353*h_n1353+  &
        16.666666666666703d0*H353**2*theta1353+  &
        468.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**4*theta1353*phi_N353**6*(H353*theta1353**2*(-22.22222222222223d0*H353*phi_nn1353  &
        -444.44444444444434d0*H353**3*theta1353*phi_nn1353-  &
        1560.d0*H353**5*theta1353**2*phi_nn1353+  &
        48.d0*theta1353*h_n1353*v_phi1353-128.d0*H353*theta1353*v_phi1353)+  &
        (-13.333333333333336d0-16.d0*H353*theta1353*h_n1353+  &
        15.999999999999945d0*H353**2*theta1353+  &
        658.6666666666665d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353*phi_N353**2*(-512.d0*H353**5*theta1353**3*phi_nn1353-  &
        21.333333333333332d0*theta1353**2*h_n1353*v_phi1353+  &
        H353*(-21.333333333333332d0*theta1353*phi_nn1353-  &
        42.666666666666664d0*theta1353**2*v_phi1353+  &
        10.666666666666666d0*thetaphi1353)+  &
        H353**3*theta1353*(-234.66666666666663d0*theta1353*phi_nn1353+  &
        32.d0*thetaphi1353))+  &
        1.d0*H353**12*theta1353**5*phi_N353**13*thetaphiphi1353+  &
        H353**10*theta1353**3*phi_N353**11*(theta1353**2*(-2.d0+  &
        thetaphi1353*(5.d0*H353**2*phi_nn1353-2.d0*v_phi1353))+  &
        H353**3*theta1353**4*(-234.d0*h_n1353-234.d0*H353)+  &
        theta1353**3*(-26.d0*H353*h_n1353-2.d0*V_phi_phi1353-52.d0*H353**2)+  &
        3.333333333333332d0*thetaphi1353**2-  &
        6.666666666666668d0*theta1353*thetaphiphi1353)+  &
        H353**5*theta1353*phi_N353**7*(18.66666666666667d0*theta1353**2*h_n1353  &
        -126.2222222222224d0*H353**2*theta1353**3*h_n1353-  &
        1560.d0*H353**4*theta1353**4*h_n1353+  &
        H353**3*theta1353**2*(26.66666666666667d0*theta1353**2*v_phi1353*phi_nn1353  &
        +56.00000000000001d0*thetaphi1353*phi_nn1353-  &
        636.4444444444443d0*theta1353)-1840.d0*H353**5*theta1353**4+  &
        H353*(theta1353**2*(-40.888888888888886d0-8.d0*thetaphi1353*v_phi1353)  &
        -32.d0*theta1353**3*V_phi_phi1353+8.d0*thetaphi1353**2-  &
        16.d0*theta1353*thetaphiphi1353))+  &
        H353**3*phi_N353**5*(-26.66666666666667d0*theta1353**2*h_n1353+  &
        159.9999999999999d0*H353**2*theta1353**3*h_n1353+  &
        1317.333333333333d0*H353**4*theta1353**4*h_n1353+  &
        1418.6666666666665d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-48.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        64.d0*thetaphi1353*phi_nn1353+604.4444444444443d0*theta1353)+  &
        H353*(theta1353**2*(53.333333333333336d0+16.d0*thetaphi1353*v_phi1353)  &
        +32.d0*theta1353**3*V_phi_phi1353+  &
        5.333333333333333d0*thetaphi1353**2+  &
        5.333333333333333d0*theta1353*thetaphiphi1353))+  &
        H353**7*theta1353**2*phi_N353**9*(-4.d0*theta1353**2*h_n1353+  &
        86.66666666666674d0*H353**2*theta1353**3*h_n1353+  &
        936.d0*H353**4*theta1353**4*h_n1353+1092.d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-6.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        26.66666666666667d0*thetaphi1353*phi_nn1353+  &
        301.3333333333333d0*theta1353)+  &
        H353*(theta1353**2*(14.666666666666666d0+  &
        6.666666666666664d0*thetaphi1353*v_phi1353)+  &
        13.333333333333334d0*theta1353**3*V_phi_phi1353-  &
        12.d0*thetaphi1353**2+16.d0*theta1353*thetaphiphi1353)))+  &
        5.d0*H353*phi_N353**2*(2.d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)*(1.3333333333333333d0-  &
        58.5d0*H353**10*theta1353**5*phi_N353**6+H353**2*theta1353*(8.d0-  &
        8.d0*theta1353*v_phi1353*phi_N353-3.333333333333334d0*phi_N353**2)  &
        +H353**8*theta1353**3*phi_N353**4*(1.d0*thetaphi1353*phi_N353**3+  &
        theta1353*(117.d0-13.d0*phi_N353**2))+  &
        H353**6*theta1353**2*phi_N353**2*(-2.d0*theta1353**2*v_phi1353*phi_N353**3  &
        -4.d0*thetaphi1353*phi_N353**3+theta1353*(-70.d0+  &
        40.666666666666664d0*phi_N353**2-0.5d0*phi_N353**4))+  &
        H353**4*theta1353*(8.d0*theta1353**2*v_phi1353*phi_N353**3+  &
        4.d0*thetaphi1353*phi_N353**3+theta1353*(12.d0-  &
        33.333333333333336d0*phi_N353**2+  &
        2.3333333333333335d0*phi_N353**4)))*(1.d0*H353*thetaphi1353*phi_N353**2  &
        +theta1353*(2.d0*h_n1353*phi_N353+  &
        2.d0*H353*phi_nn1353))*(H353**11*theta1353**5*(2.d0*h_n1353-9.d0*H353-  &
        117.d0*H353**3*theta1353)*thetaphi1353*phi_N353**12+  &
        (7.111111111111111d0+42.666666666666664d0*H353**2*theta1353+  &
        64.d0*H353**4*theta1353**2)*phi_nn1353+phi_N353*(7.111111111111111d0  &
        +21.333333333333332d0*H353*theta1353*h_n1353+  &
        64.d0*H353**3*theta1353**2*h_n1353+  &
        H353**2*theta1353*(42.666666666666664d0-  &
        53.333333333333336d0*theta1353*v_phi1353*phi_nn1353)+  &
        64.d0*H353**4*theta1353**2)+  &
        H353*theta1353*phi_N353**3*(10.666666666666666d0*h_n1353-  &
        103.11111111111111d0*H353**2*theta1353*h_n1353-  &
        512.d0*H353**4*theta1353**2*h_n1353+H353*(-32.d0-  &
        21.333333333333332d0*thetaphi1353*v_phi1353-  &
        10.666666666666666d0*theta1353*V_phi_phi1353)+  &
        H353**3*(64.d0*theta1353**2*v_phi1353*phi_nn1353+  &
        37.33333333333333d0*thetaphi1353*phi_nn1353-  &
        263.1111111111111d0*theta1353)-  &
        501.3333333333333d0*H353**5*theta1353**2)+  &
        H353**6*theta1353**2*phi_N353**8*(H353*theta1353**2*(10.666666666666668d0*H353*phi_nn1353  &
        +237.33333333333326d0*H353**3*theta1353*phi_nn1353+  &
        936.d0*H353**5*theta1353**2*phi_nn1353-  &
        13.333333333333327d0*theta1353*h_n1353*v_phi1353+  &
        53.333333333333336d0*H353*theta1353*v_phi1353)+(9.333333333333336d0  &
        +8.d0*H353*theta1353*h_n1353+0.8888888888887974d0*H353**2*theta1353-  &
        780.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**2*phi_N353**4*(H353*theta1353**2*(26.666666666666654d0*H353*phi_nn1353  &
        +462.22222222222206d0*H353**3*theta1353*phi_nn1353+  &
        1317.333333333333d0*H353**5*theta1353**2*phi_nn1353-  &
        32.d0*theta1353*h_n1353*v_phi1353+128.d0*H353*theta1353*v_phi1353)+  &
        (5.333333333333333d0+21.333333333333332d0*H353*theta1353*h_n1353-  &
        30.22222222222222d0*H353**2*theta1353-  &
        256.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**8*theta1353**3*phi_N353**10*(H353**2*theta1353**2*(-2.d0*phi_nn1353-  &
        234.d0*H353**4*theta1353**2*phi_nn1353+  &
        theta1353*(-52.d0*H353**2*phi_nn1353-8.d0*v_phi1353))+(-2.d0-  &
        6.666666666666664d0*H353*theta1353*h_n1353+  &
        16.666666666666703d0*H353**2*theta1353+  &
        468.d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353**4*theta1353*phi_N353**6*(H353*theta1353**2*(-22.22222222222223d0*H353*phi_nn1353  &
        -444.44444444444434d0*H353**3*theta1353*phi_nn1353-  &
        1560.d0*H353**5*theta1353**2*phi_nn1353+  &
        48.d0*theta1353*h_n1353*v_phi1353-128.d0*H353*theta1353*v_phi1353)+  &
        (-13.333333333333336d0-16.d0*H353*theta1353*h_n1353+  &
        15.999999999999945d0*H353**2*theta1353+  &
        658.6666666666665d0*H353**4*theta1353**2)*thetaphi1353)+  &
        H353*phi_N353**2*(-512.d0*H353**5*theta1353**3*phi_nn1353-  &
        21.333333333333332d0*theta1353**2*h_n1353*v_phi1353+  &
        H353*(-21.333333333333332d0*theta1353*phi_nn1353-  &
        42.666666666666664d0*theta1353**2*v_phi1353+  &
        10.666666666666666d0*thetaphi1353)+  &
        H353**3*theta1353*(-234.66666666666663d0*theta1353*phi_nn1353+  &
        32.d0*thetaphi1353))+  &
        1.d0*H353**12*theta1353**5*phi_N353**13*thetaphiphi1353+  &
        H353**10*theta1353**3*phi_N353**11*(theta1353**2*(-2.d0+  &
        thetaphi1353*(5.d0*H353**2*phi_nn1353-2.d0*v_phi1353))+  &
        H353**3*theta1353**4*(-234.d0*h_n1353-234.d0*H353)+  &
        theta1353**3*(-26.d0*H353*h_n1353-2.d0*V_phi_phi1353-52.d0*H353**2)+  &
        3.333333333333332d0*thetaphi1353**2-  &
        6.666666666666668d0*theta1353*thetaphiphi1353)+  &
        H353**5*theta1353*phi_N353**7*(18.66666666666667d0*theta1353**2*h_n1353  &
        -126.2222222222224d0*H353**2*theta1353**3*h_n1353-  &
        1560.d0*H353**4*theta1353**4*h_n1353+  &
        H353**3*theta1353**2*(26.66666666666667d0*theta1353**2*v_phi1353*phi_nn1353  &
        +56.00000000000001d0*thetaphi1353*phi_nn1353-  &
        636.4444444444443d0*theta1353)-1840.d0*H353**5*theta1353**4+  &
        H353*(theta1353**2*(-40.888888888888886d0-8.d0*thetaphi1353*v_phi1353)  &
        -32.d0*theta1353**3*V_phi_phi1353+8.d0*thetaphi1353**2-  &
        16.d0*theta1353*thetaphiphi1353))+  &
        H353**3*phi_N353**5*(-26.66666666666667d0*theta1353**2*h_n1353+  &
        159.9999999999999d0*H353**2*theta1353**3*h_n1353+  &
        1317.333333333333d0*H353**4*theta1353**4*h_n1353+  &
        1418.6666666666665d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-48.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        64.d0*thetaphi1353*phi_nn1353+604.4444444444443d0*theta1353)+  &
        H353*(theta1353**2*(53.333333333333336d0+16.d0*thetaphi1353*v_phi1353)  &
        +32.d0*theta1353**3*V_phi_phi1353+  &
        5.333333333333333d0*thetaphi1353**2+  &
        5.333333333333333d0*theta1353*thetaphiphi1353))+  &
        H353**7*theta1353**2*phi_N353**9*(-4.d0*theta1353**2*h_n1353+  &
        86.66666666666674d0*H353**2*theta1353**3*h_n1353+  &
        936.d0*H353**4*theta1353**4*h_n1353+1092.d0*H353**5*theta1353**4+  &
        H353**3*theta1353**2*(-6.d0*theta1353**2*v_phi1353*phi_nn1353-  &
        26.66666666666667d0*thetaphi1353*phi_nn1353+  &
        301.3333333333333d0*theta1353)+  &
        H353*(theta1353**2*(14.666666666666666d0+  &
        6.666666666666664d0*thetaphi1353*v_phi1353)+  &
        13.333333333333334d0*theta1353**3*V_phi_phi1353-  &
        12.d0*thetaphi1353**2+16.d0*theta1353*thetaphiphi1353)))+  &
        1.d0*phi_N353*(-2.d0+  &
        H353**2*theta1353*phi_N353**2)*(-0.6666666666666666d0+  &
        1.d0*H353**2*theta1353*phi_N353**2)*(1.3333333333333333d0-  &
        58.5d0*H353**10*theta1353**5*phi_N353**6+H353**2*theta1353*(8.d0-  &
        8.d0*theta1353*v_phi1353*phi_N353-3.333333333333334d0*phi_N353**2)  &
        +H353**8*theta1353**3*phi_N353**4*(1.d0*thetaphi1353*phi_N353**3+  &
        theta1353*(117.d0-13.d0*phi_N353**2))+  &
        H353**6*theta1353**2*phi_N353**2*(-2.d0*theta1353**2*v_phi1353*phi_N353**3  &
        -4.d0*thetaphi1353*phi_N353**3+theta1353*(-70.d0+  &
        40.666666666666664d0*phi_N353**2-0.5d0*phi_N353**4))+  &
        H353**4*theta1353*(8.d0*theta1353**2*v_phi1353*phi_N353**3+  &
        4.d0*thetaphi1353*phi_N353**3+theta1353*(12.d0-  &
        33.333333333333336d0*phi_N353**2+  &
        2.3333333333333335d0*phi_N353**4)))*(-21.333333333333332d0*theta1353**2*h_n1353**2*v_phi1353*phi_N353**2  &
        +theta1353*h_n1353**2*phi_N353*(21.333333333333332d0+  &
        10.666666666666666d0*phi_N353**2)+  &
        H353**13*theta1353**6*phi_N353**10*(-3276.d0*h_n1353*thetaphi1353*phi_N353**2  &
        +theta1353*(-234.d0*phi_N353*h_nn1353+h_n1353*(-3276.d0*phi_N353-  &
        5850.d0*phi_nn1353)))+7.111111111111111d0*phi_nn1353+  &
        7.111111111111111d0*phi_nnn1353+  &
        H353*(h_n1353*thetaphi1353*phi_N353**2*(42.666666666666664d0+  &
        21.333333333333332d0*phi_N353**2)+  &
        theta1353*(phi_N353*(21.333333333333332d0+  &
        10.666666666666666d0*phi_N353**2)*h_nn1353+  &
        h_n1353*(85.33333333333333d0*phi_N353+(-64.d0-  &
        85.33333333333333d0*thetaphi1353*v_phi1353)*phi_N353**3+  &
        106.66666666666666d0*phi_nn1353-  &
        10.666666666666664d0*phi_N353**2*phi_nn1353))+  &
        theta1353**2*phi_N353*(-21.333333333333332d0*v_phi1353*phi_N353*h_nn1353  &
        +h_n1353*(v_phi1353*(-85.33333333333333d0*phi_N353-  &
        149.33333333333331d0*phi_nn1353)-  &
        42.666666666666664d0*phi_N353**2*V_phi_phi1353)))+  &
        H353**14*theta1353**5*phi_N353**9*(-702.d0*thetaphi1353**2*phi_N353**4+  &
        theta1353**2*(-2340.d0*phi_nn1353**2+phi_N353*(-2574.d0*phi_nn1353  &
        -234.d0*phi_nnn1353))+  &
        theta1353*(thetaphi1353*phi_N353**2*(-1638.d0*phi_N353-  &
        3042.0000000000005d0*phi_nn1353)-  &
        117.d0*phi_N353**4*thetaphiphi1353))+  &
        H353**5*theta1353*phi_N353**2*(16.d0*h_n1353*thetaphi1353**2*phi_N353**5  &
        +theta1353**2*(phi_N353*(-512.d0+159.9999999999999d0*phi_N353**2+  &
        18.66666666666667d0*phi_N353**4)*h_nn1353+  &
        h_n1353*(-3008.d0*phi_N353+3626.6666666666656d0*phi_N353**3+  &
        (-245.33333333333331d0+144.d0*thetaphi1353*v_phi1353)*phi_N353**5  &
        -4608.d0*phi_nn1353+3573.333333333332d0*phi_N353**2*phi_nn1353-  &
        2.6666666666666856d0*phi_N353**4*phi_nn1353))+  &
        theta1353**3*phi_N353**4*(48.d0*v_phi1353*h_nn1353+  &
        h_n1353*(-768.d0*v_phi1353-144.d0*phi_N353*V_phi_phi1353))+  &
        theta1353*(-16.d0*thetaphi1353*phi_N353**4*h_nn1353+  &
        h_n1353*phi_N353**2*(thetaphi1353*(-3071.9999999999995d0+  &
        575.9999999999993d0*phi_N353**2+  &
        112.00000000000003d0*phi_N353**4-480.d0*phi_N353*phi_nn1353)-  &
        112.d0*phi_N353**3*thetaphiphi1353)))+  &
        H353**9*theta1353**3*phi_N353**6*(6.666666666666664d0*h_n1353*thetaphi1353**2*phi_N353**5  &
        +theta1353**2*(phi_N353*(-1560.d0+  &
        86.66666666666674d0*phi_N353**2)*h_nn1353+  &
        h_n1353*(-18400.d0*phi_N353+3013.333333333333d0*phi_N353**3+  &
        (-20.d0-20.d0*thetaphi1353*v_phi1353)*phi_N353**5-  &
        26520.d0*phi_nn1353+3153.333333333333d0*phi_N353**2*phi_nn1353-  &
        20.d0*phi_N353**4*phi_nn1353))+  &
        theta1353**3*h_n1353*phi_N353**3*(v_phi1353*(-80.d0*phi_N353-  &
        60.d0*phi_nn1353)-20.d0*phi_N353**2*V_phi_phi1353)+  &
        theta1353*(-6.666666666666664d0*thetaphi1353*phi_N353**4*h_nn1353+  &
        h_n1353*phi_N353**2*(thetaphi1353*(-15600.d0+  &
        600.0000000000008d0*phi_N353**2-  &
        333.3333333333333d0*phi_N353*phi_nn1353)-  &
        73.33333333333334d0*phi_N353**3*thetaphiphi1353)))+  &
        H353**11*theta1353**4*phi_N353**8*(10.d0*h_n1353*thetaphi1353**2*phi_N353**5  &
        +theta1353**2*(phi_N353*(936.d0-26.d0*phi_N353**2)*h_nn1353+  &
        h_n1353*(13103.999999999998d0*phi_N353-624.d0*phi_N353**3+  &
        19656.d0*phi_nn1353-910.d0*phi_N353**2*phi_nn1353))+  &
        theta1353*(2.d0*thetaphi1353*phi_N353**4*h_nn1353+  &
        h_n1353*phi_N353**2*(thetaphi1353*(11232.d0-264.d0*phi_N353**2+  &
        84.d0*phi_N353*phi_nn1353)+  &
        14.000000000000002d0*phi_N353**3*thetaphiphi1353)))+  &
        H353**3*(42.666666666666664d0*h_n1353*thetaphi1353**2*phi_N353**5+  &
        theta1353**2*(phi_N353*(64.d0-103.11111111111111d0*phi_N353**2-  &
        26.66666666666667d0*phi_N353**4)*h_nn1353+  &
        h_n1353*(256.d0*phi_N353-1052.4444444444443d0*phi_N353**3+  &
        (213.33333333333331d0-32.d0*thetaphi1353*v_phi1353)*phi_N353**5+  &
        320.d0*phi_nn1353-1248.d0*phi_N353**2*phi_nn1353-  &
        26.666666666666742d0*phi_N353**4*phi_nn1353))+  &
        theta1353**3*phi_N353**3*(-32.d0*v_phi1353*phi_N353*h_nn1353+  &
        h_n1353*(v_phi1353*(512.d0*phi_N353+128.d0*phi_nn1353)+  &
        96.d0*phi_N353**2*V_phi_phi1353))+  &
        theta1353*(21.333333333333332d0*thetaphi1353*phi_N353**4*h_nn1353+  &
        h_n1353*phi_N353**2*(thetaphi1353*(256.d0-  &
        327.11111111111114d0*phi_N353**2-  &
        106.66666666666669d0*phi_N353**4+  &
        234.66666666666663d0*phi_N353*phi_nn1353)+  &
        42.666666666666664d0*phi_N353**3*thetaphiphi1353)))+  &
        H353**7*theta1353**2*phi_N353**4*(-72.d0*h_n1353*thetaphi1353**2*phi_N353**5  &
        +theta1353**2*(phi_N353*(1317.333333333333d0-  &
        126.2222222222224d0*phi_N353**2-4.d0*phi_N353**4)*h_nn1353+  &
        h_n1353*(11349.333333333332d0*phi_N353-  &
        5091.555555555555d0*phi_N353**3+(117.33333333333334d0-  &
        13.333333333333327d0*thetaphi1353*v_phi1353)*phi_N353**5+  &
        17125.33333333333d0*phi_nn1353-  &
        4439.111111111111d0*phi_N353**2*phi_nn1353+  &
        49.333333333333336d0*phi_N353**4*phi_nn1353))+  &
        theta1353**3*phi_N353**3*(-13.333333333333327d0*v_phi1353*phi_N353*h_nn1353  &
        +h_n1353*(v_phi1353*(426.6666666666667d0*phi_N353+  &
        106.66666666666674d0*phi_nn1353)+  &
        93.33333333333334d0*phi_N353**2*V_phi_phi1353))+  &
        theta1353*(8.d0*thetaphi1353*phi_N353**4*h_nn1353+  &
        h_n1353*phi_N353**2*(thetaphi1353*(10538.666666666664d0-  &
        497.7777777777792d0*phi_N353**2-32.d0*phi_N353**4+  &
        512.d0*phi_N353*phi_nn1353)+136.d0*phi_N353**3*thetaphiphi1353)))+  &
        H353**12*theta1353**4*phi_N353**7*(-3042.0000000000005d0*theta1353**3*h_n1353**2*phi_N353**4  &
        +theta1353**2*(7488.d0*phi_nn1353**2-  &
        520.d0*phi_N353**2*phi_nn1353**2+phi_N353**3*(-572.d0*phi_nn1353  &
        -52.d0*phi_nnn1353)+phi_N353*(9828.d0*phi_nn1353+  &
        936.d0*phi_nnn1353))+thetaphi1353*phi_N353**4*(thetaphi1353*(2340.d0  &
        -45.d0*phi_N353**2+25.d0*phi_N353*phi_nn1353)+  &
        5.d0*phi_N353**3*thetaphiphi1353)+  &
        theta1353*(phi_N353**4*thetaphiphi1353*(468.d0-9.d0*phi_N353**2+  &
        18.d0*phi_N353*phi_nn1353)+  &
        thetaphi1353*phi_N353**2*(-312.d0*phi_N353**3+10296.d0*phi_nn1353+  &
        phi_N353*(6552.d0+55.00000000000001d0*phi_nn1353**2)+  &
        phi_N353**2*(-420.d0*phi_nn1353+5.d0*phi_nnn1353))+  &
        1.d0*phi_N353**7*thetaphiphiphi1353))+  &
        H353**6*(-93.3333333333333d0*theta1353**5*h_n1353**2*v_phi1353*phi_N353**8  &
        +8.d0*thetaphi1353**3*phi_N353**8+  &
        theta1353*thetaphi1353*phi_N353**5*(thetaphi1353*(-512.d0+  &
        31.99999999999989d0*phi_N353**2+  &
        18.66666666666667d0*phi_N353**4-72.d0*phi_N353*phi_nn1353)-  &
        16.d0*phi_N353**3*thetaphiphi1353)+  &
        theta1353**3*phi_N353*(-1504.d0*phi_N353*phi_nn1353+  &
        3022.2222222222217d0*phi_N353**3*phi_nn1353-  &
        286.2222222222222d0*phi_N353**5*phi_nn1353-  &
        1024.d0*phi_nn1353**2+  &
        1848.8888888888882d0*phi_N353**2*phi_nn1353**2-  &
        133.33333333333337d0*phi_N353**4*phi_nn1353**2-  &
        512.d0*phi_N353*phi_nnn1353+  &
        462.22222222222206d0*phi_N353**3*phi_nnn1353-  &
        22.22222222222223d0*phi_N353**5*phi_nnn1353+  &
        thetaphi1353*(v_phi1353*phi_N353**5*(-512.d0*phi_N353-  &
        248.d0*phi_nn1353)+phi_N353**7*(56.00000000000001d0*h_n1353**2-  &
        136.d0*V_phi_phi1353))-8.d0*v_phi1353*phi_N353**7*thetaphiphi1353)  &
        +  &
        theta1353**2*phi_N353**3*(-24.d0*thetaphi1353**2*v_phi1353*phi_N353**5  &
        +phi_N353**2*thetaphiphi1353*(-256.d0+  &
        15.999999999999945d0*phi_N353**2+  &
        9.333333333333336d0*phi_N353**4-176.d0*phi_N353*phi_nn1353)+  &
        thetaphi1353*(1813.3333333333328d0*phi_N353**3-  &
        122.66666666666666d0*phi_N353**5-2560.d0*phi_nn1353+  &
        8.d0*phi_N353**4*phi_nn1353+phi_N353*(-1504.d0-  &
        320.d0*phi_nn1353**2)+  &
        phi_N353**2*(1482.6666666666658d0*phi_nn1353-  &
        64.d0*phi_nnn1353))-  &
        16.d0*phi_N353**5*thetaphiphiphi1353)+  &
        theta1353**4*phi_N353**4*(h_n1353**2*(9221.333333333332d0*phi_N353  &
        -883.5555555555568d0*phi_N353**3-  &
        28.000000000000004d0*phi_N353**5)-  &
        272.d0*phi_N353**2*V_phi_phi1353*phi_nn1353-  &
        240.d0*v_phi1353*phi_nn1353**2+  &
        v_phi1353*phi_N353*(-768.d0*phi_nn1353-48.d0*phi_nnn1353)-  &
        128.d0*phi_N353**3*V_phi_phi1353-  &
        32.d0*phi_N353**4*V_phi_phi_phi1353))+  &
        H353**2*(-96.d0*theta1353**3*h_n1353**2*v_phi1353*phi_N353**4-  &
        21.333333333333332d0*thetaphi1353**2*v_phi1353*phi_N353**4+  &
        thetaphi1353*phi_N353*(42.666666666666664d0*phi_N353-  &
        32.d0*phi_N353**3+64.d0*phi_nn1353)+  &
        phi_N353**3*(10.666666666666666d0+  &
        5.333333333333333d0*phi_N353**2)*thetaphiphi1353+  &
        theta1353*(42.666666666666664d0*phi_nn1353-  &
        96.d0*phi_N353**2*phi_nn1353-  &
        42.666666666666664d0*phi_N353*phi_nn1353**2+  &
        42.666666666666664d0*phi_nnn1353-  &
        21.333333333333332d0*phi_N353**2*phi_nnn1353+  &
        thetaphi1353*(v_phi1353*phi_N353**2*(-85.33333333333333d0*phi_N353  &
        -170.66666666666666d0*phi_nn1353)+phi_N353**4*(64.d0*h_n1353**2  &
        -42.666666666666664d0*V_phi_phi1353))-  &
        21.333333333333332d0*v_phi1353*phi_N353**4*thetaphiphi1353)+  &
        theta1353**2*(h_n1353**2*(192.d0*phi_N353-  &
        309.33333333333337d0*phi_N353**3-  &
        80.00000000000001d0*phi_N353**5)-  &
        85.33333333333333d0*phi_N353**2*V_phi_phi1353*phi_nn1353-  &
        53.333333333333336d0*v_phi1353*phi_nn1353**2+  &
        v_phi1353*phi_N353*(-85.33333333333333d0*phi_nn1353-  &
        53.333333333333336d0*phi_nnn1353)-  &
        42.666666666666664d0*phi_N353**3*V_phi_phi1353-  &
        10.666666666666666d0*phi_N353**4*V_phi_phi_phi1353))+  &
        H353**10*theta1353**2*phi_N353**5*(9.999999999999996d0*thetaphi1353**3*phi_N353**7  &
        +theta1353*thetaphi1353*phi_N353**4*(thetaphi1353*(-3120.d0+  &
        66.66666666666681d0*phi_N353**2-  &
        70.00000000000003d0*phi_N353*phi_nn1353)-  &
        20.000000000000004d0*phi_N353**3*thetaphiphi1353)+  &
        theta1353**3*(-12880.d0*phi_N353*phi_nn1353+  &
        2712.d0*phi_N353**3*phi_nn1353-22.d0*phi_N353**5*phi_nn1353-  &
        9360.d0*phi_nn1353**2+  &
        1898.666666666666d0*phi_N353**2*phi_nn1353**2-  &
        20.d0*phi_N353**4*phi_nn1353**2-1560.d0*phi_N353*phi_nnn1353+  &
        237.33333333333326d0*phi_N353**3*phi_nnn1353-  &
        2.d0*phi_N353**5*phi_nnn1353+  &
        thetaphi1353*(v_phi1353*phi_N353**5*(-48.d0*phi_N353-  &
        57.99999999999999d0*phi_nn1353)+phi_N353**7*(22.d0*h_n1353**2-  &
        14.000000000000002d0*V_phi_phi1353))-  &
        2.d0*v_phi1353*phi_N353**7*thetaphiphi1353)+  &
        theta1353**2*phi_N353**2*(-10.d0*thetaphi1353**2*v_phi1353*phi_N353**5  &
        +phi_N353**2*thetaphiphi1353*(-780.d0+  &
        16.666666666666703d0*phi_N353**2-  &
        100.00000000000001d0*phi_N353*phi_nn1353)+  &
        thetaphi1353*(1506.6666666666665d0*phi_N353**3-10.d0*phi_N353**5-  &
        14040.d0*phi_nn1353-10.d0*phi_N353**4*phi_nn1353+  &
        phi_N353*(-9200.d0-240.00000000000006d0*phi_nn1353**2)+  &
        phi_N353**2*(1353.3333333333333d0*phi_nn1353-  &
        26.66666666666667d0*phi_nnn1353))-  &
        6.666666666666668d0*phi_N353**5*thetaphiphiphi1353)+  &
        theta1353**4*phi_N353**3*(h_n1353**2*(10296.d0*phi_N353-  &
        286.d0*phi_N353**3)-  &
        28.000000000000004d0*phi_N353**2*V_phi_phi1353*phi_nn1353-  &
        54.d0*v_phi1353*phi_nn1353**2+  &
        v_phi1353*phi_N353*(-80.d0*phi_nn1353-6.d0*phi_nnn1353)-  &
        8.d0*phi_N353**3*V_phi_phi1353-  &
        2.d0*phi_N353**4*V_phi_phi_phi1353))+  &
        H353**8*theta1353*phi_N353**3*(-24.d0*thetaphi1353**3*phi_N353**7+  &
        theta1353*thetaphi1353*phi_N353**4*(thetaphi1353*(1975.9999999999993d0  &
        +2.6666666666663916d0*phi_N353**2-6.d0*phi_N353**4+  &
        60.d0*phi_N353*phi_nn1353)+24.d0*phi_N353**3*thetaphiphi1353)+  &
        theta1353**3*(7093.333333333333d0*phi_N353*phi_nn1353-  &
        4455.11111111111d0*phi_N353**3*phi_nn1353+  &
        132.d0*phi_N353**5*phi_nn1353+5269.333333333332d0*phi_nn1353**2  &
        -2666.666666666666d0*phi_N353**2*phi_nn1353**2+  &
        85.33333333333334d0*phi_N353**4*phi_nn1353**2+  &
        1317.333333333333d0*phi_N353*phi_nnn1353-  &
        444.44444444444434d0*phi_N353**3*phi_nnn1353+  &
        10.666666666666668d0*phi_N353**5*phi_nnn1353+  &
        thetaphi1353*(v_phi1353*phi_N353**5*(266.6666666666667d0*phi_N353+  &
        193.33333333333334d0*phi_nn1353)+  &
        phi_N353**7*(-59.99999999999998d0*h_n1353**2+  &
        73.33333333333333d0*V_phi_phi1353))+  &
        6.666666666666664d0*v_phi1353*phi_N353**7*thetaphiphi1353)+  &
        theta1353**2*phi_N353**2*(26.666666666666654d0*thetaphi1353**2*v_phi1353*phi_N353**5  &
        +phi_N353**2*thetaphiphi1353*(658.6666666666665d0+  &
        0.8888888888887974d0*phi_N353**2-2.d0*phi_N353**4+  &
        200.d0*phi_N353*phi_nn1353)+  &
        thetaphi1353*(-2545.7777777777774d0*phi_N353**3+  &
        58.666666666666664d0*phi_N353**5+  &
        9221.333333333332d0*phi_nn1353+  &
        22.66666666666667d0*phi_N353**4*phi_nn1353+  &
        phi_N353*(5674.666666666666d0+392.d0*phi_nn1353**2)+  &
        phi_N353**2*(-1770.666666666667d0*phi_nn1353+  &
        56.00000000000001d0*phi_nnn1353))+  &
        16.d0*phi_N353**5*thetaphiphiphi1353)+  &
        theta1353**4*phi_N353**3*(h_n1353**2*(-14040.d0*phi_N353+  &
        780.0000000000007d0*phi_N353**3)+  &
        146.66666666666669d0*phi_N353**2*V_phi_phi1353*phi_nn1353+  &
        186.66666666666669d0*v_phi1353*phi_nn1353**2+  &
        v_phi1353*phi_N353*(426.6666666666667d0*phi_nn1353+  &
        26.66666666666667d0*phi_nnn1353)+  &
        53.333333333333336d0*phi_N353**3*V_phi_phi1353+  &
        13.333333333333334d0*phi_N353**4*V_phi_phi_phi1353))+  &
        H353**4*(240.d0*theta1353**4*h_n1353**2*v_phi1353*phi_N353**6+  &
        thetaphi1353*phi_N353**3*(thetaphi1353*(32.d0-  &
        30.22222222222222d0*phi_N353**2-  &
        13.333333333333336d0*phi_N353**4+  &
        63.99999999999999d0*phi_N353*phi_nn1353)+  &
        16.d0*phi_N353**3*thetaphiphi1353)+theta1353**2*(64.d0*phi_nn1353-  &
        789.3333333333333d0*phi_N353**2*phi_nn1353+  &
        266.66666666666663d0*phi_N353**4*phi_nn1353-  &
        469.33333333333326d0*phi_N353*phi_nn1353**2+  &
        106.66666666666663d0*phi_N353**3*phi_nn1353**2+  &
        64.d0*phi_nnn1353-234.66666666666663d0*phi_N353**2*phi_nnn1353+  &
        26.666666666666654d0*phi_N353**4*phi_nnn1353+  &
        thetaphi1353*(v_phi1353*phi_N353**4*(384.d0*phi_N353+  &
        272.d0*phi_nn1353)+phi_N353**6*(-80.d0*h_n1353**2+  &
        112.d0*V_phi_phi1353))+  &
        16.d0*v_phi1353*phi_N353**6*thetaphiphi1353)+  &
        theta1353*phi_N353*(32.d0*thetaphi1353**2*v_phi1353*phi_N353**5+  &
        phi_N353**2*thetaphiphi1353*(32.d0-  &
        30.22222222222222d0*phi_N353**2-  &
        13.333333333333336d0*phi_N353**4+  &
        63.99999999999999d0*phi_N353*phi_nn1353)+  &
        thetaphi1353*(-526.2222222222222d0*phi_N353**3+  &
        106.66666666666666d0*phi_N353**5+192.d0*phi_nn1353-  &
        26.6666666666667d0*phi_N353**4*phi_nn1353+phi_N353*(128.d0+  &
        111.99999999999999d0*phi_nn1353**2)+  &
        phi_N353**2*(-590.2222222222222d0*phi_nn1353+  &
        37.33333333333333d0*phi_nnn1353))+  &
        5.333333333333333d0*phi_N353**5*thetaphiphiphi1353)+  &
        theta1353**3*phi_N353**2*(h_n1353**2*(-2560.d0*phi_N353+  &
        799.9999999999994d0*phi_N353**3+  &
        93.33333333333336d0*phi_N353**5)+  &
        224.d0*phi_N353**2*V_phi_phi1353*phi_nn1353+  &
        192.d0*v_phi1353*phi_nn1353**2+  &
        v_phi1353*phi_N353*(512.d0*phi_nn1353+64.d0*phi_nnn1353)+  &
        128.d0*phi_N353**3*V_phi_phi1353+  &
        32.d0*phi_N353**4*V_phi_phi_phi1353)))))/((0.6666666666666666d0  &
        -1.d0*H353**2*theta1353*phi_N353**2)**10*((-1.d0*phi_N353**4*(-2.d0+  &
        H353**2*theta1353*phi_N353**2)*(1.3333333333333333d0-  &
        58.5d0*H353**10*theta1353**5*phi_N353**6+H353**2*theta1353*(8.d0-  &
        8.d0*theta1353*v_phi1353*phi_N353-3.333333333333334d0*phi_N353**2)  &
        +H353**8*theta1353**3*phi_N353**4*(1.d0*thetaphi1353*phi_N353**3+  &
        theta1353*(117.d0-13.d0*phi_N353**2))+  &
        H353**6*theta1353**2*phi_N353**2*(-2.d0*theta1353**2*v_phi1353*phi_N353**3  &
        -4.d0*thetaphi1353*phi_N353**3+theta1353*(-70.d0+  &
        40.666666666666664d0*phi_N353**2-0.5d0*phi_N353**4))+  &
        H353**4*theta1353*(8.d0*theta1353**2*v_phi1353*phi_N353**3+  &
        4.d0*thetaphi1353*phi_N353**3+theta1353*(12.d0-  &
        33.333333333333336d0*phi_N353**2+  &
        2.3333333333333335d0*phi_N353**4))))/(0.6666666666666666d0-  &
        1.d0*H353**2*theta1353*phi_N353**2)**4)**1.75d0)
        
    !    if ( ieee_is_nan(z_NN) ) then
    !                   write(*,*) 'z_NN is NaN'
    !                   !print **, coff,phi353,phi_n353,H353,a353
                      ! read(*,*) rff
    !    endif
          
	return
	end function
	
	
	end module mod3_5
