module mod5
!use precision
use mod1
use mod2
use mod3
use mod4
implicit none
contains


!ode solver for the perturbation part

		subroutine  odesolve2(rs501,t501,hn,k501,coff,info,x0,pain)
		use mod1
		use mod2
		use mod3
		use mod4
		implicit none
		
		real(dl)   ,dimension(11)::k1501
		real(dl)  ,dimension(11)::k2501
		real(dl)  ,dimension(11)::k3501
		real(dl)  ,dimension(11)::k4501
		real(dl)  ,dimension(11)::rs501
		real(dl)  ,dimension(11)::rserr501
		real(dl)  ,dimension(11)::gg501
		real(dl)   k501,t501,hn,rees501,check501
		real(dl)  ,dimension(11)::rstemp501
		real(dl):: coff(25),info(40),x0
		integer pain


    integer i501,is501
    REAL(dl), PARAMETER :: A2=0.2_dl,A3=0.3_dl,A4=0.6_dl,A5=1.0_dl,&
    A6=0.875_dl,B21=0.2_dl,B31=3.0_dl/40.0_dl,B32=9.0_dl/40.0_dl,&
    B41=0.3_dl,B42=-0.9_dl,B43=1.2_dl,B51=-11.0_dl/54.0_dl,&
    B52=2.5_dl,B53=-70.0_dl/27.0_dl,B54=35.0_dl/27.0_dl,&
    B61=1631.0_dl/55296.0_dl,B62=175.0_dl/512.0_dl,&
    B63=575.0_dl/13824.0_dl,B64=44275.0_dl/110592.0_dl,&
    B65=253.0_dl/4096.0_dl,C1=37.0_dl/378.0_dl


    REAL(dl), PARAMETER ::C3=250.0_dl/621.0_dl,C4=125.0_dl/594.0_dl

    REAL(dl), PARAMETER ::C6=512.0_dl/1771.0_dl,DC1=C1-2825.0_dl/27648.0_dl,&
    DC3=C3-18575.0_dl/48384.0_dl,DC4=C4-13525.0_dl/55296.0_dl,&
    DC5=-277.0_dl/14336.0_dl,DC6=C6-0.25_dl
    real(dl) ,dimension(11):: ak1501,ak2501,ak3501,ak4501,ak5501,ak6501,ree501


		call g(coff,gg501,rs501,t501,k501,hn,pain)
		ak1501=gg501
		rstemp501=rs501+B21*hn*ak1501
		
		call g(coff,gg501,rstemp501,t501+A2*hn,k501,hn,pain)
		 ak2501 = gg501
		
		rstemp501=rs501+hn*(B31*ak1501+B32*ak2501)

		
		call g(coff,gg501,rstemp501,t501+A3*hn,k501,hn,pain)
		ak3501 = gg501
		rstemp501=rs501+hn*(B41*ak1501+B42*ak2501+B43*ak3501)
	          
				
		call g(coff,gg501,rstemp501,t501+A4*hn,k501,hn,pain)
	          ak4501 = gg501
		rstemp501=rs501+hn*(B51*ak1501+B52*ak2501+B53*ak3501+B54*ak4501)		
		call g(coff,gg501,rstemp501,t501+A5*hn,k501,hn,pain)
		   ak5501 = gg501
		rstemp501=rs501+hn*(B61*ak1501+B62*ak2501+B63*ak3501+B64*ak4501+B65*ak5501)
		call g(coff,gg501,rstemp501,t501+A6*hn,k501,hn,pain)
		   ak6501=gg501
		rs501=rs501+hn*(C1*ak1501+C3*ak3501+C4*ak4501+C6*ak6501)
		rserr501=hn*(DC1*ak1501+DC3*ak3501+DC4*ak4501+DC5*ak5501+DC6*ak6501)
		ree501=rserr501/rs501
		rees501=0.d0
		do 10 i501 = 1, 11
                rees501=rees501+ree501(i501)**2
10		continue
		check501=sqrt(rees501)/11
		!!print **, rees
		if (check501.ge.1e-2) then
			!print **, "ode_solve2::Truncation error too large, reduce stepsize."

                        !print **, "error kind 9 : mod5.f odesolve2 line 82"
			error_index=9
			error(error_index)=error(error_index)+1
			open (unit = 1586, file = "error.log" ,action='write',position='append')
			write(1586,*) "error kind and count :" ,error_index,error(error_index)
			write(1586,*)  "error kind 9 : mod5.f odesolve2 line 82"
			write(1586,*) "Truncation error too large, reduce stepsize. Toterror%=",check501
			write(1586,*)x0,(coff(is501),is501=1,16),(info(is501),is501=1,25)
			close(1586)



		endif
		return
		end subroutine


!ode solver for the background part
		subroutine  odesolve(r502,t502,hn,coff,info,x0)
		use mod1
		use mod2
		use mod3
		use mod4
		implicit none
		
		real(dl)   ,dimension(3)::k1502
		real(dl)  ,dimension(3)::k2502
		real(dl)  ,dimension(3)::k3502
		real(dl)  ,dimension(3)::k4502
		real(dl)  ,dimension(3)::r502
		real(dl)  ,dimension(3)::rserr502
		real(dl)  ,dimension(3)::ff502
		real(dl)   k502,t502,hn,rees502,check502
		real(dl)  ,dimension(3)::rstemp502
		real(dl):: coff(25),info(40),x0
		
integer i502,is502
REAL(dl), PARAMETER :: A2=0.2_dl,A3=0.3_dl,A4=0.6_dl,A5=1.0_dl,&
A6=0.875_dl,B21=0.2_dl,B31=3.0_dl/40.0_dl,B32=9.0_dl/40.0_dl,&
B41=0.3_dl,B42=-0.9_dl,B43=1.2_dl,B51=-11.0_dl/54.0_dl,&
B52=2.5_dl,B53=-70.0_dl/27.0_dl,B54=35.0_dl/27.0_dl,&
B61=1631.0_dl/55296.0_dl,B62=175.0_dl/512.0_dl,&
B63=575.0_dl/13824.0_dl,B64=44275.0_dl/110592.0_dl,&
B65=253.0_dl/4096.0_dl,C1=37.0_dl/378.0_dl


REAL(dl), PARAMETER ::C3=250.0_dl/621.0_dl,C4=125.0_dl/594.0_dl

REAL(dl), PARAMETER ::C6=512.0_dl/1771.0_dl,DC1=C1-2825.0_dl/27648.0_dl,&
DC3=C3-18575.0_dl/48384.0_dl,DC4=C4-13525.0_dl/55296.0_dl,&
DC5=-277.0_dl/14336.0_dl,DC6=C6-0.25_dl
real(dl) ,dimension(3):: ak1502,ak2502,ak3502,ak4502,ak5502,ak6502,ree502

		call f(coff,ff502,r502,t502,k502,hn)
		ak1502=ff502
		rstemp502=r502+B21*hn*ak1502
		
		call f(coff,ff502,rstemp502,t502+A2*hn,k502,hn)
		 ak2502 = ff502
		
		rstemp502=r502+hn*(B31*ak1502+B32*ak2502)

		
		call f(coff,ff502,rstemp502,t502+A3*hn,k502,hn)
		ak3502 = ff502
		rstemp502=r502+hn*(B41*ak1502+B42*ak2502+B43*ak3502)
	          
				
		call f(coff,ff502,rstemp502,t502+A4*hn,k502,hn)
	          ak4502 = ff502
		rstemp502=r502+hn*(B51*ak1502+B52*ak2502+B53*ak3502+B54*ak4502)		
		call f(coff,ff502,rstemp502,t502+A5*hn,k502,hn)
		   ak5502 = ff502
		rstemp502=r502+hn*(B61*ak1502+B62*ak2502+B63*ak3502+B64*ak4502+B65*ak5502)
		call f(coff,ff502,rstemp502,t502+A6*hn,k502,hn)
		   ak6502=ff502
		r502=r502+hn*(C1*ak1502+C3*ak3502+C4*ak4502+C6*ak6502)
		rserr502=hn*(DC1*ak1502+DC3*ak3502+DC4*ak4502+DC5*ak5502+DC6*ak6502)
		ree502=rserr502/r502
		rees502=0.d0
		do 10 i502 = 1, 3
                rees502=rees502+ree502(i502)**2
10		continue
		check502=sqrt(rees502)/3
		!!print **, rees
		if (check502.ge.1e-2) then
			!print **, "ode_solve1::Truncation error too large, reduce stepsize."
                        !print **, "error kind 10 : mod5.f odesolve1 line 169"
			error_index=10
			error(error_index)=error(error_index)+1
			open (unit = 1586, file = "error.log" ,action='write',position='append')
			write(1586,*) "error kind and count :" ,error_index,error(error_index)
			write(1586,*)  "error kind 10 : mod5.f odesolve1 line 169"
			write(1586,*) "Truncation error too large, reduce stepsize. Toterror%=",check502
			write(1586,*)x0,(coff(is502),is502=1,16),(info(is502),is502=1,25)
			close(1586)


		endif
		return
		end subroutine



end module mod5


!this provides the derivatives to ode solver in array format for perturbation part
             
	subroutine g(coff,gg503,rs503,t503,k503,hn,pain)
	    !use precision
		use mod1
		use mod2
		use mod3
		use mod4
		implicit none
		real(dl) phi_n503,a503,H503,phi503,fphi503,fphi_n503,fH503,fqr503
		real(dl) fqc503,fq_nr503,fq_nc503,fvr503,fvc503,fv_nr503,fv_nc503
		real(dl) k503,qr503,qc503,vr503,vc503,a0503,t503
		real(dl) hn,q_nr503,q_nc503,v_nr503,v_nc503
		real(dl)  ,dimension(11)::rs503
		!real(dl)  ,dimension(3)::f
		!real(dl)  ,dimension(11)::ree
		real(dl)  ,dimension(11)::gg503
		real(dl):: coff(25),qnnr,qnnc,vnnr,vnnc
		integer pain 
		
		a0503=coff(8)
		    a503=a0503*exp(t503)
		    phi503 = rs503(1)
		    phi_n503= rs503(2)
		    H503=rs503(3)
		    fphi503 = phi_n503
		    fphi_n503=phi_nn(coff,phi503,H503,phi_n503)
		    fH503=H_n(coff,phi503,H503,phi_n503)
		    qr503=rs503(4)
		    qc503=rs503(5)
		    q_nr503=rs503(6)
		    q_nc503=rs503(7)
		    vr503=rs503(8)
		    vc503=rs503(9)
		    v_nr503=rs503(10)
		    v_nc503=rs503(11)
		    
		    call qvnn(coff,rs503,k503,pain,a503,qnnr,qnnc,vnnr,vnnc)
		    fqr503=q_nr503
		    fqc503=q_nc503
		    fq_nr503=qnnr!qnn(coff,phi503,H503,phi_n503,a503,k503,q_nr503,qr503)
		    fq_nc503=qnnc!(coff,phi503,H503,phi_n503,a503,k503,q_nc503,qc503)
		    fvr503=v_nr503
		    fvc503=v_nc503
		    fv_nr503=vnnr!0.0!vnn(coff,phi503,H503,phi_n503,a503,k503,v_nr503,vr503)
		    fv_nc503=vnnc!0.0!vnn(coff,phi503,H503,phi_n503,a503,k503,v_nc503,vc503)
		    
		    
		    
    gg503=[fphi503,fphi_n503,fH503,fqr503,fqc503,fq_nr503,fq_nc503,fvr503,fvc503,fv_nr503,fv_nc503]

		! !print **,"gg",a0,t
		!return
		end subroutine




