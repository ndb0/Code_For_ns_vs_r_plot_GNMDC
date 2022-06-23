module spectral
!use precision
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

subroutine efold_optimized_old(coff,n,info,x0)
use mod1
use mod2
use mod3
use mod4
use mod5
implicit none

		real(dl) :: a0,H0,phi0,phi_n0,temp,odef,ndef,nu, start_nu
		real(dl)  :: phi00,phi_n00,H00,N,hn,t,ee,a,ddf,zz,k,phi000,dk,fac1,dk1,keep,nu2
		real(dl) , dimension(3)::r
		real(dl) coff(25),info(40),x0,min_phi_n
		integer i
		character(3) rff
                integer ad1,ad_in,is
		!print **, "spectral begin",nef
		k=.1d0
		a0=coff(8)
		
		call set_initial(x0,coff,phi00,phi_n00,H00)
		
		r=[phi00,phi_n00,H00]
		hn=1E-3
		min_phi_n=.10d0
		!n=200.d0
		ee=epsln2(coff,phi_n00,H00,phi00)
		dk=.01d0
		temp=100.d0
		start_nu=coff(5)
526		coff(5)=start_nu
!		coff(5)=.5
	       !  coff(5)=coff(5)/2.0d0
36		         nu=coff(5)
		        dk=.1d0
 			keep=coff(5)
			a0=coff(8)
			
			call set_initial(x0,coff,phi00,phi_n00,H00)

			r=[phi00,phi_n00,H00]
		!	!print **,"okkkkk",r,coff(5),n
			N=0.d0
			ee=epsln2(coff,phi_n00,H00,phi00)
			do while (ee.le.1.00)
			!!print **,'n1111,r,ee', N,r,ee
				call odesolve(r,N,hn,coff,info,x0)
				phi00=r(1)
				phi_n00=r(2)
				H00=r(3)
				ee=epsln2(coff,phi_n00,H00,phi00)
				n=n+hn
				


25				if (n.ge.(nef*1.3d0)) then
				      keep=coff(5)
				      coff(5)=coff(5)*.9d0
				      nu=coff(5)
				      coff(6)=coff(6)*(.9d0)**2
				      dk=(keep-nu)/(2.d0*nu)
				      go to 36
				  end if     
				!!print **,'n,phi0,ee', N,phi0,ee,coff(5)
		enddo
!	        !print **,"need now" ,coff(5),n,nef
		temp=n
		odef=nef-temp
		dk1=dk
!		!print **, "spectral loop begins",n,nef,dk
		!!print **,'n,r,ee', N,r,ee
		do while ((abs(n-nef)/nef).ge. 0.010d0)
			!!print **,'n,r,ee', N,phi00,ee

			keep=coff(5)
!                        !print **, "value is saved",keep
                        coff(5)=nu*dk+nu
!##########################################coff(5)=(nu*dk+nu)######################################################
 !                       !print **, "am in add",keep,dk
!			ad1=int(1.d0/dk)
 !                       ad_in=1
!		!	 !print **, "am in add",ad1
  !                      coff(5)=dk*keep
!			Do while (ad_in.le.ad1)
!			  coff(5)=coff(5)+keep/ad1
!			  ad_in=ad_in+1
                          !!print **, ad_in
!                      end do

!##################################################################################################################
			
		!dk1=dk
56 		Print *,"n,coff5",n,coff(5)
        call set_initial(x0,coff,phi00,phi_n00,H00)
			r=[phi00,phi_n00,H00]
			N=0.d0
			ee=epsln2(coff,phi_n00,H00,phi00)
			min_phi_n=.10d0
			do while (ee.le.1.00)
				call odesolve(r,N,hn,coff,info,x0)
				phi00=r(1)
				phi_n00=r(2)
				H00=r(3)
				ee=epsln2(coff,phi_n00,H00,phi00)
				n=n+hn

				if (n.ge.(nef*3.d0)) then
				      go to 67
				  end if  
71			enddo
67 if((abs(dk1).le.(1.0e-17)).or.(abs(dk).ge.(1e4)).or. (abs(coff(5)).le.(1.0e-7)).or.(abs(coff(5)).ge.(1.e9)))then
				Print *,"error : in phi_ini line 56: we have hit a wall:to proceed say 'yes'",n,dk
					!print **, "error kind 2 : spectral.f find_phi_ini line 67 "
					!print **, n,dk,dk1,coff(5)
					error_index=2
					error(error_index)=error(error_index)+1
					open (unit = 1586, file = "error.log" ,action='write',position='append')
					write(1586,*) "error kind and count :" ,error_index,error(error_index)
					write(1586,*) "error kind 2 : spectral.f find_phi_ini line 67, n,dk= ",n,dk
					write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
					close(1586)
			
				if (start_nu.ge.1.0)then
					start_nu=.115
					go to 526
				else
					start_nu=10.d0
					go to 526
				endif


                        !print **, "want to continue?"
			read (*,*)rff
			if (rff=="yes")then
			!print **,"fine then"
			else
			STOP
			endif
				endif	
!			!print **,"need now" ,coff(5),n,nef	
			!odef=nef-temp
			ndef=nef-n
			temp=n
			
			if(((odef)/abs(odef)).eq.((ndef)/abs(ndef)))then
					!if (abs(odef).le.abs(ndef))then
					!	dk=dk*1.1
					!else
					!!print **, "c1"
                                       ! !print **, "loop1"
					dk=dk+dk1
			else		
				!!print **, "c2"
				nu=keep
				nu2=coff(5)
				coff(5)=nu
			!	!print **, "loop2",nu,nu2
				dk=((nu2-nu)/(5.d0*nu))
				dk1=dk
				!fac1=fac1*.9
				go to 56
			endif
		enddo
               coff(9)=n
		!print **, "final E-fold no",n
		coff(10)=phi00
		
      return
      end subroutine

      
      subroutine set_efolds(coff,nstar,info,x0,nch)  
    use mod1
	    	
	implicit none
    
        
		real(dl) , dimension(3)::r
		real(dl) coff(25),info(40),x0,min_phi_n,temp
		
		character(3) rff
		integer nch
        real(dl) detnorm,n,phi_end,co5p,co5n,co5star,nstar,hn,neta,init
        
        
        
        !nch=15
        !########################################################################################### 
        temp=coff(nch)
        hn=1d-3
        init=coff(nch)
       ! coff(nch)=init/1.1
        call run_back(coff,n,x0,phi_end,hn)
        do while(n.le.nef)
      !  print *,"n<nef: n,coff(10),x0,coff(15)",n,coff(5),x0,coff(15)
          coff(nch)=coff(nch)/1.01d0
          call run_back(coff,n,x0,phi_end,hn)
        enddo

       ! coff(nch)=init/1.1
       ! print *, "in set_efold temp,nef=",temp,nef
        
       ! call run_back(coff,n,x0,phi_end,hn)
        
        
        
        do while(n.ge.nef*3.d0)
     !   print *,"n>3nef: n,coff(10),x0,coff(15)",n,coff(5),x0,coff(15)
          co5p=coff(nch)
          coff(nch)=coff(nch)*1.02d0
          call run_back(coff,n,x0,phi_end,hn)
        enddo  
        
        detnorm=n-nef
        
        if(detnorm.le.0.d0)then
          co5n=coff(nch)
         ! goto 908
          !call run_back(coff,n,x0,phi_end,hn)
          do while((phi_end.le.coff(1)*coff(nch)).and.(n.lt.3.d0*nef))
         !  print *,"phi_end>phi0: ,n,coff(15),phi_end,coff(15)*coff(1)",n,coff(15),phi_end,coff(15)*coff(1)
            coff(nch)=coff(nch)/1.01d0
            call run_back(coff,n,x0,phi_end,hn)
            detnorm=n-nef
          enddo  
      !    print *,"phi_end>phi0: ,n,coff(15),phi_end,coff(15)*coff(1)",n,coff(15),phi_end,coff(15)*coff(1)
          if(detnorm.le.0.d0)then
             co5n=coff(nch)
          else
             co5n=coff(nch)*1.01d0
            co5p=coff(nch)
          endif   
      !    908 co5n=coff(5)
        else
        co5p=coff(nch)
    588 coff(nch)=coff(nch)*1.02d0
        call run_back(coff,n,x0,phi_end,hn)
        detnorm=n-nef
        if(detnorm.ge.0.d0)then
            co5p=coff(nch)
            go to 588
            else
            co5n=coff(nch)
            endif
        endif
    

    
        co5star=coff(nch)
        nstar=n
    !    neta=nch
            do while (abs(nstar-nef)/nef.ge.1e-3)
             !   print *, "before starting get_efolds : ", co5star,nstar,co5n,co5p
                call get_efolds(co5star,nstar,co5n,co5p,nch,info,x0,coff,hn) 
             !   print *, "after getefolds : ",co5star,nstar,co5n,co5p
               ! read(*,*)  rff
            enddo
            coff(15)=co5star
          !  print *, "final E-fold no",nstar
    return
    end subroutine
    
subroutine get_efolds(co5star,nstar,co5p,co5n,nch,info,x0,coff,hn)  
      	use mod1
	    use OMP_LIB	
	implicit none



	
	
	
	
	
	
   
	character(90) rff,filename
	real(dl):: coff(25),info(40),x0,neta
	


	real(dl) co5star,ntar,rati,nstar,hn
	
	integer ipara,iparacount,nch

	real(dl) copara(50),diff(50),efolds(50),phi_end(50)
	real(dl) diffmin,diffp,diffn,co5p,co5n
	
	real(dl) nw,phi_endw
	
	real(dl) pricof(25)
    real(dl) co5i,co5
    real(dl) dry,co5temp
    integer ercount1,ercount2,ICOR
	!pragma omp parallel
		


	ICOR=OMP_GET_NUM_PROCS()

	neta=float(ICOR)-2.0
!	PRINT *, icor,neta
	dry=(co5p-co5n)/neta
	co5=co5n   
	ipara=1
	!print *, "initiating get_efolds ",co5p,co5n
    do while (co5.lt.(co5p))
                    copara(ipara)=co5
                    !print **,ipara,copara(ipara)
                    ipara=ipara+1
                    
                    co5=co5+dry
    end do
    copara(ipara)=co5p
	iparacount=ipara
	!!CALL OMP_SET_NUM_THREADS(16)
	!OMP parallel private(nw,phi_endw,pricof)
	
    !$OMP PARALLEL DO private(nw,phi_endw,pricof)
    !Ipara is private by default
        DO Ipara=1,iparacount
                
               ! !print **, ipara,copara(ipara)
                pricof=coff
                pricof(nch)=copara(ipara)
                call run_back(pricof,nw,x0,phi_endw,hn)
                !kkw(ipara)=kpara2
                efolds(ipara)=nw
                phi_end(ipara)=phi_endw
                diff(ipara)=nw-nef
                !!print **, kk(i),power(i)
                !i=i+1
                

        ENDDO
        !$OMP END PARALLEL DO
        !OMP end parallel
        ercount1=0
        ercount2=0
        DO Ipara=1,iparacount
            if(efolds(ipara).lt.(nef*3.d0))then
              ercount1=ercount1+1
            if(phi_end(ipara).gt.(copara(ipara)*coff(1)))then
                ercount2=ercount2+1
            endif

             else
              ercount2=ercount2+1       
            endif  

        !    print *,"cof15,n,phi_,cof15*cof1 :", copara(ipara),efolds(ipara),phi_end(ipara),copara(ipara)*coff(1)
        ENDDO
        
        if((ercount1.ge.1).and.(ercount2.ge.1))then
            diffmin=3.d0*nef
            diffp=3.d0*nef
            diffn=-3.d0*nef
            DO Ipara=1,iparacount
                if(diffmin.ge.abs(diff(ipara)))then
                    diffmin=abs(diff(ipara))
                    co5star=copara(ipara)
                    nstar=efolds(ipara)
                endif
                if((diff(ipara).ge.0.d0).and.(diffp.ge.diff(ipara)))then
                    diffp=diff(ipara)
                    co5n=copara(ipara)
                    if(efolds(ipara).ge.(nef*3.d0))then
                    
                        exit
                    endif   
                endif
                if((diff(ipara).le.0.d0).and.(diffn.le.diff(ipara)))then
                    diffn=diff(ipara)
                   ! if(phi_end(ipara).lt.(copara(ipara)*x0))then
                            co5p=copara(ipara)
                   ! endif
                endif
             ! print *, copara(ipara),efolds(ipara)
            ENDDO
        else if((ercount1.ge.1).and.(ercount2.lt.1))then
            co5n=copara(iparacount)/1.01d0
            co5p=copara(iparacount)
         else if((ercount1.lt.1).and.(ercount2.ge.1))then   
            co5n=copara(1)
            co5p=copara(1)*1.01d0
         else
            open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :in get_efold unrecoverable: existing"

                 write(1586,*) "!######################################################!"
               close(1586)
            read(*,*) rff
         endif   
	!print **, "concluding get_efolds"
    
    end subroutine
    
subroutine run_back(coff,n,x0,phi_end,hn)
        use mod1
        use mod2
        use mod3
        use mod4
        use mod5
        implicit none
        real(dl) , dimension(3)::r
		real(dl) coff(25),info(40),x0
		real(dl) phi_end
		real(dl) phi00,phi_n00,h00,n,hn,ee
		real(dl) need(25)
		
		!print *, "in run back"
        call set_initial(x0,coff,phi00,phi_n00,H00)
        r=[phi00,phi_n00,H00]
        N=0.d0
        ee=epsln2(coff,phi_n00,H00,phi00)
      
        do while (ee.lt.1.00)
				call odesolve(r,N,hn,coff,info,x0)
				phi00=r(1)
				phi_n00=r(2)
				H00=r(3)
				ee=epsln2(coff,phi_n00,H00,phi00)
				n=n+hn

				if (n.ge.(nef*3.d0)) then
				      exit
				  end if  
71		enddo
        phi_end=phi00
        coff(10)=phi_end
    return
    end subroutine
    
    subroutine run_back_assignvalues(coff,hn,nstart,nend)
        use mod1
        use mod2
        use mod3
        use mod4
        use mod5
        implicit none
        real(dl) , dimension(3)::r
		real(dl) coff(25),info(40),x0
		real(dl) phi_end
		real(dl) phi00,phi_n00,h00,n,hn,ee
		real(dl) need(25),hsn,aa
		integer i
		real(dl) css1,cssn1,ctt1,zns1,znt1,nstart,nend
		character(150) rff,filename
		
		
		i=1
		!print *, "in run back"
        call touch_initial(x0,coff,phi00,phi_n00,H00)
        r=[phi00,phi_n00,H00]
        N=0.d0
        aa=coff(8)*Exp(n)
        ee=epsln2(coff,phi_n00,H00,phi00)
        call assignvals(coff,r,aa,css1,cssn1,ctt1,zns1,znt1)
      css(i)=css1
      cssn(i)=cssn1
      ctt(i)=ctt1
      zns(i)=zns1
      znt(i)=znt1
          WRITE(filename,'(a,i4.4,a)') "../data/asign_back.dat"     
            OPEN(unit=785,file=filename ,action='write',position='append')
        do while (n.lt.40.0)
				call odesolve(r,N,hn,coff,info,x0)
				phi00=r(1)
				phi_n00=r(2)
				H00=r(3)
			!	ee=epsln2(coff,phi_n00,H00,phi00)
				n=n+hn
				aa=coff(8)*Exp(n)
				i=i+1
				if((n.ge.nstart-2.0).and.(n.le.nend+2.0))then
				call assignvals(coff,r,aa,css1,cssn1,ctt1,zns1,znt1)
                
                css(i)=css1
                cssn(i)=cssn1
                ctt(i)=ctt1
                zns(i)=zns1
                znt(i)=znt1
                
                write(785,*) n,css1,ctt1
                endif
                
       
             
71		enddo
      close(785)
       ! phi_end=phi00
       ! coff(10)=phi00
    return
    end subroutine
    
    
      
subroutine get_N(coff,n,info,x0,jin)
        use mod1
        use mod2
        use mod3
        use mod4
        use mod5
        use mod3_5
        implicit none

		real(dl) :: a0,H0,phi0,phi_n0,temp,odef,ndef,nu, start_nu
		real(dl)  :: phi00,phi_n00,H00,N,hn,t,ee,a,ddf,zz,k,phi000,dk,fac1,dk1,keep,nu2
		real(dl) , dimension(3)::r
		real(dl) coff(25),info(40),x0,min_phi_n,phi_nn00,H_n00,z_00,z_n00,cs_00,cs_n00,thet,new,new1
		integer i,jin
		character(150) rff,filename
        integer ad1,ad_in,is,iss
                
        integer erflag
                
        real(dl) a_gn ,ct_00,ct_n00,diff_ct_c,cs_sqr1 ,ct_new,thetphi,thetad
       real(dl)  v1, v2, v3, v4, v5, pot, vprime, factor
        real(dl) v6,v7,v8,PR,PT,ratio,Qs,Qt,ff,gg,epD,ct,c_scalar
        
        
        
       iss=0 
     !  print *, "iss=",iss
        error(25)=0
        
		!print **, "N____:",nef
		k=.1d0
		a0=coff(8)
        call set_initial(x0,coff,phi00,phi_n00,H00)
		!H00=(-0.05555555555555556d0*(-6.d0+phi_N00**2+sqrt(36.d0-  &
         !       12.d0*(1.d0+6.d0*theta(coff,phi00)*V(coff,phi00)*phi_N00**4  &
         !       +phi_N00**4)))/(theta(coff,phi00)*phi_N00**2))
		
		!H00=sqrt(V(coff,phi00)/(3.d0-(phi_n00**2)/2.d0))
		

		r=[phi00,phi_n00,H00]
		hn=1E-3
		min_phi_n=.10d0
		!print **,jin,phi00,phi_n00,H00,V(coff,phi00)
		!WRITE(filename,'(a,i4.4,a)') "../data/ini_back",jin,".dat"
        phi_nn00=phi_nn(coff,phi00,H00,phi_n00)
        h_n00=h_n(coff,phi00,H00,phi_n00)
        WRITE(filename,'(a,i4.4,a)') "../data/ini_back",jin,".TXT"
        OPEN(unit=785,file=filename )
        close(785,STATUS='DELETE' )
		OPEN(unit=785,file=filename ,action='write',position='append')
		
		
               
		N=0.d0
        !print **,"okkkkk",r,coff(5),n
        !print **, n,h00,phi00,phi_n00,phi_nn00,h_n00
        ee=epsln2(coff,phi_n00,H00,phi00)
      !  print *, "***************starting value of epsilon******************** epsilon=",ee
        do while (ee.lt.1.d0)
			
			!!print **,'n1111,r,ee', N,r,ee
			    phi_nn00=phi_nn(coff,phi00,H00,phi_n00)
			    h_n00=h_n(coff,phi00,H00,phi_n00)
				call odesolve(r,N,hn,coff,info,x0)
				phi00=r(1)
				phi_n00=r(2)
				H00=r(3)
				ee=epsln2(coff,phi_n00,H00,phi00)
				n=n+hn
				a_gn=a0*Exp(n)
				z_00=zzzz(coff,phi_n00,a_gn,H00,phi00)
                z_n00=Z_n(coff,phi00,phi_n00,H00,a_gn)
				cs_00=cs(coff,phi00,phi_n00,H00)
				cs_n00=cs_n(coff,phi00,phi_n00,H00)
				ct_00=cst(coff,phi00,phi_n00,H00,a_gn)
				ct_n00=cs1t(coff,phi00,phi_n00,H00,a_gn)
				diff_ct_c=ct_00-1.0
				 thet=theta(coff,phi00)
				 thetphi=thetaphi(coff,phi00)
				 new=z_n00/(a_gn*H00)
                cs_sqr1=cs_new_sqr(coff,phi00,phi_n00,H00,a_gn)
                ct_new=cst_new(coff,phi00,phi_N00,H00,a_gn)
                thetad=theta_diped(coff,phi00)
                factor=(0.05*cs_00)/(a_gn*H00)
                pot=V(coff,phi00)
                vprime=V_phi(coff,phi00)
                epD=1.5*thet*(H00*phi_n00)**2     !eqn 2.19 of 1910.00622
                ff=(1.0-epD/3.0)/(1.0-epD)
                gg=(epD/3.0)*(1+3.0*(H00**2)*thet*(1.0+epD)/(1.0-epD/3.0))     !eqn 3.2 of 1910.00622
                Qs=(ff**2)*gg/(thet*H00**2)
                Qt=(1.0-epD/3.0)/4.0
                ct=sqrt(1.0+2.0*epD/3.0)
                c_scalar=(epD/3.0*gg)*(1+epD+3.0*(H00**2)*thet*(1.0+epD+4.0*epD/(9.0*ff))+ &
                6.0*H00*h_n00*thet*(1.0-epD/3.0))/(1.0-epD/3.0)        ! eq 3.3  
                PR=(H00**2)/(8.0*(Pi**2)*Qs*c_scalar**3)        ! 3.4
                PT=(H00**2)/(2.0*(Pi**2)*Qt*ct**3)
                ratio=PT/PR
                if (factor.ge.0.001) then
                    v1=0.5*(vprime/pot)**2
                    v2=1+3.0*(H00**2)*thet+thetphi*phi_n00*H00**2
                    v3=16.0*v1/v2
                    v4=phi00
                    v5=coff(16)
                    v6=PR
                    v7=PT
                    v8=ratio
                endif
				! print *, "******************** epsilon=",ee
				!!print **, n,h00,phi00,phi_n00,phi_nn00,h_n00
           
            write(785,*) n,a_gn,phi00,phi_n00,phi_nn00,H00,h_n00,ee,z_00,z_n00,cs_00,cs_n00,ct_00,ct_n00 &
                                        ,diff_ct_c,thet,new,cs_sqr1,ct_new,thetphi,thetad
				
    if((cs_00 /=cs_00).or.(cs_n00 /= cs_n00).or.(z_00 /=z_00).or.(z_n00 /= z_n00))then
                error_index=22
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
             write(1586,*)  "error kind 22: spectral line no 504"
             write(1586,*)"Some of the values are NaN"
            write(1586,*) "N,z_00,z_n00,cs_00,cs_n00", N,z_00,z_n00,cs_00,cs_n00
            write(1586,*)(coff(is),is=1,8)
            write(1586,*) "!######################################################!"
               close(1586)
     endif           
	     if(	ee.lt.-1e2)then
            error_index=23
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
             write(1586,*)  "error kind 23: spectral line no 504"
             write(1586,*)"1st slow roll param is negetive"
            write(1586,*) "N,ee", N,ee
            write(1586,*)(coff(is),is=1,8)
            write(1586,*) "!######################################################!"
               close(1586)
               
             !  print *, "***************error******************** epsilon=",ee
               exit
	       endif		
        if(n.ge.nef*2.0)then
                exit
        endif        



				
	if((cs_00 /=cs_00).or.(cs_n00 /= cs_n00).or.(z_00 /=z_00).or.(z_n00 /= z_n00))then !.or.(	ee.le.0.0)
          if((1.0+9.0*(H00**2)*theta(coff,phi00)).lt.0.0)then
                if(iss.lt.1)then
                coff(24)=9.0*(H00**2)*theta(coff,phi00)
               ! print *, "NAN predicted",iss,coff(24)
                endif
                iss=iss+1
            !    print *, "get_N: error predicted ==",iss,coff(24)
            else     
             !   print *, "get_N: No error predicted ##",iss
       endif                
				error(25)=1
			!	print *, "N,ee,thet,cs_00,cs_n00,z_00,z_n00,phi,V:::",N,ee,thet,cs_00,cs_n00 &
             !                           ,z_00,z_n00,phi00,V(coff,phi00)
				!exit
	endif	
	
	   
		
				
        enddo   
        close(785)
        WRITE(filename,'(a,i4.4,a)') "../data/analytical",jin,".TXT"
        OPEN(unit=785,file=filename )
        close(785,STATUS='DELETE' )

        OPEN(unit=999,file=filename ,action='write',position='append')
        write(999,*) v5,v1,v2,v3,v4,v6,v7,v8
        close(999)
        info(36)=h00
      
 	!	print *, "final E-fold no",n
		coff(10)=phi00
		!stop
      return
      end  subroutine   
      
      
      

!######################################################################################################
!This subroutine is used to evolve the system in sub-horizon limit, for a given initial value of φ
!obtained from find phi ini subroutine.
!Here also we start with slow roll condition,evolve the system with efolds using back ground
!equations and in each step the horizon limit is checked with ddf parameter. Whenever the
!sub-horizon approximation wears off, we exit the program. Passing all the variables
!#######################################################################################################


subroutine subhorizon_evolution(k_sh,phi0_sh,phi_n0_sh,H0_sh,N_sh,hn,coff,info,x0,li)
        use mod1
        use mod2
        use mod3
        use mod4
        use mod5
        implicit none

		real(dl) :: a0_sh,H0_sh,phi0_sh,phi_n0_sh,k_sh,hn,ee_sh,N_sh
		real(dl) :: a_sh,ddf_sh,lddf_sh,h_n0_sh,phi_nn0_sh,cs1sub
		real(dl) , dimension(3)::r_sh
		real(dl):: coff(25),info(40),x0
		
		!real(dl)  , dimension(5000000,20):: psdata
		integer subnum_sh,is_sh,rff
		!integer, dimension(5) :: ii
		integer li
		


        a0_sh=coff(8)
        call touch_initial(x0,coff,phi0_sh,phi_n0_sh,H0_sh)

		!H00=(-0.05555555555555556d0*(-6.d0+phi_N00**2+sqrt(36.d0-  &
         !       12.d0*(1.d0+6.d0*theta(coff,phi00)*V(coff,phi00)*phi_N00**4  &
         !       +phi_N00**4)))/(theta(coff,phi00)*phi_N00**2))
		
		!H0_sh=sqrt(V(coff,phi0_sh)/3.0*(1.0+dhv00))
		!!print **, phi0_sh,phi_n0_sh,h0_sh,V(coff,phi0_sh)
		
		
		r_sh=[phi0_sh,phi_n0_sh,H0_sh]
		ee_sh=epsln2(coff,phi_n0_sh,H0_sh,phi0_sh)
		a_sh=a0_sh*exp(N_sh)
		cs1sub=cs(coff,phi0_sh,phi_n0_sh,H0_sh)
		ddf_sh=k_sh*cs1sub/(a_sh*H0_sh)
		!print **, "Enetering stage 1 ::","at N=",N_sh,"k/(aH)=",ddf_sh
		! !print **,"ddfs",ddf,k,a,H0
		lddf_sh=100.d0
		subnum_sh=0
		!!print **, "before sub_starts",r
		!!print **, "before sub_starts",n,a
		!read(*,*)rff
		!OPEN(unit=784,file="../data/subhor_back.dat" ,action='write',position='append')
		do while (ddf_sh.ge.lddf_sh)
			
			! write (11,*) N,r,ee
			subnum_sh=subnum_sh+1
			
			call odesolve(r_sh,N_sh,hn,coff,info,x0)
			N_sh=N_sh+hn
			li=li+1
			!pain=pain+1
			phi0_sh=r_sh(1)
			phi_n0_sh=r_sh(2)
			H0_sh=r_sh(3)
			
			!!print **, "check",j
           ! !print **,"subH:", k,n,h0,phi0,phi_n0
			ee_sh=epsln2(coff,phi_n0_sh,H0_sh,phi0_sh)
			a_sh=a0_sh*exp(N_sh-hn)
			cs1sub=cs(coff,phi0_sh,phi_n0_sh,H0_sh)
			ddf_sh=k_sh/(a_sh*H0_sh/cs1sub)
			!h_n0_sh=H_n(coff,phi0_sh,H0_sh,phi_n0_sh)
			!phi_nn0_sh=phi_nn(coff,phi0_sh,H0_sh,phi_n0_sh)
			!!print **, n,phi0,phi_n0,h0,ddf
			!write(784,*) n,phi0,phi_n0,h0,ee,phi_nn0,h_n0,ddf
			!j=j+1
                    !do 10 i = 1, 3
                    !psdata(j,i)=r(i)
                    !10  continue
                    !psdata(j,13)=N
                    !psdata(j,17)=ddf
                    !psdata(j,18)=k
			!!print **,"ddf",ddf,k,a,H0
		enddo
		! close(784)
		!print **, "number of loops in subhorizon limit",subnum_sh
		if (subnum_sh.le.2) then
 print *, "*******Not satisfied subhorizon condition*********",subnum_sh

				!print **, "error kind 3 : spectral.f subhorizon_evolution line 229 "
				error_index=3
				error(error_index)=error(error_index)+1
				open (unit = 1586, file = "error.log" ,action='write',position='append')
				write(1586,*) "error kind and count :" ,error_index,error(error_index)
				write(1586,*)  "error kind 3 : spectral.f subhorizon_evolution line 229 "
				write(1586,*) "not satisfied subhorizon condition change a0",subnum_sh
				write(1586,*)x0,(coff(is_sh),is_sh=1,16),(info(is_sh),is_sh=1,25)
				close(1586)



				!   STOP
			endif
		
		
return 
end subroutine



!###################################################################################################
!This subroutine is used to implement Bunch Davies vacuum condition , as a starting point of
!perturbative evolution.
!Whenever the sub-horizon limit is crossed, the curvature perturbation is initialized both for
!scalar and tensor part , with 4 additional complex variables (q,v and their first derivatives).
!####################################################################################################
subroutine perturbation_initiation(coff,qr0_pi,qc0_pi,q_nr0_pi,q_nc0_pi,vr0_pi,vc0_pi, &
                    v_nr0_pi,v_nc0_pi,k_pi,a_pi,H0_pi,phi0_pi,phi_n0_pi)
        !!!!!!!!!!!!!!!!q is representing R everywhere!!!!!!!!!!!!!!!
	use mod1
        use mod2
        use mod3
        use mod4
        use mod5
        use mod3_5
        implicit none		
        real(dl) qr0_pi,qc0_pi,q_nr0_pi,q_nc0_pi,vr0_pi,vc0_pi,v_nr0_pi,v_nc0_pi,k_pi,a_pi,H0_pi
        real(dl) cs1ini,cs_n1ini,phi0_pi,phi_n0_pi,coff(25)
        real(dl) z1_ini,z_n1_ini,zt_ini,zt1byzt_ini,ct_ini,ct_n1ini
        
			!PI=4.D0*DATAN(1.D0)
			cs1ini=cs(coff,phi0_pi,phi_n0_pi,H0_pi)
			cs_n1ini=cs_n(coff,phi0_pi,phi_n0_pi,H0_pi)
			z1_ini=zzzz(coff,phi_n0_pi,a_pi,H0_pi,phi0_pi)   ! z
			z_n1_ini=Z_n(coff,phi0_pi,phi_n0_pi,H0_pi,a_pi)! this is actually  z'/z
			zt_ini=zt(coff,phi0_pi,phi_n0_pi,H0_pi,a_pi)
			zt1byzt_ini =z1tbyzt(coff,phi0_pi,phi_n0_pi,H0_pi,a_pi)
			ct_ini=cst(coff,phi0_pi,phi_n0_pi,H0_pi,a_pi)
			ct_n1ini=cs1t(coff,phi0_pi,phi_n0_pi,H0_pi,a_pi)
			
			qr0_pi=1.0/(z1_ini*sqrt(2.0*k_Pi*cs1ini))
			qc0_pi=0.0
			
			q_nc0_pi=-(k_pi*cs1ini)/(sqrt(2.0*k_pi*cs1ini)*z1_ini*a_pi*H0_pi)
			q_nr0_pi=-(k_pi*cs_n1ini)/(z1_ini*(k_pi*cs1ini*2.0)**1.5)-z_n1_ini/(sqrt(2.0*k_pi*cs1ini)*z1_ini)
			!q_nc0=-1*k/sqrt(2*k)/(a*H0)
			!q_nr0=0.0d0
			!vr0_pi=1/sqrt(2*k_pi*cs1ini)  
			vr0_pi=1.0/(zt_ini*sqrt(2.0*k_Pi*ct_ini))
			!v_nc0_pi=-1*k_pi/sqrt(2*k_pi)/(a_pi*H0_pi)
			vc0_pi=0.0
            !v_nr0_pi=0.0
			v_nc0_pi=-(k_pi*ct_ini)/(sqrt(2.0*k_pi*ct_ini)*zt_ini*a_pi*H0_pi)
			v_nr0_pi=-(k_pi*ct_n1ini)/(zt_ini*(k_pi*ct_ini*2.0)**1.5)-zt1byzt_ini/(sqrt(2.0*k_pi*ct_ini)*zt_ini)
			!print **,"Completed stage 2"
return
end subroutine



!#####################################################################################################
!This is the computational bottleneck of the program, for the required accuracy we need to take
!very fine efolds binning but, This program essentially solve 11 coupled differential equations and
!takes long time for small efold-steps. After mode initiation, each of these variables being complex
!, their real and imaginary part is evolved separately.On the top of it the back ground is evolved.
!The error is compared with the function values and whenever error increases significantly warning
!messeges are !print *out. Here also in each step the perturbation length scale is compared with
!Horizon scale, whenever the super horizon limit is achieved, we stop the program and pass the
!q,v values as a function of efolds along with z also as a function of efolds.
!#######################################################################################################



subroutine perturbation_evolution(phi0_pe,phi_n0_pe,H0_pe,qr0_pe,qc0_pe,q_nr0_pe,q_nc0_pe,vr0_pe, &
        vc0_pe,v_nr0_pe,v_nc0_pe, N_pe,k_pe,hs,coff,srsp_pe, info,x0,li)
    use mod1
    use mod2
    use mod3
    use mod4
    use mod5
    use mod3_5
    implicit none
        
	real(dl) :: a0_pe,H0_pe,phi0_pe,phi_n0_pe,hs,uddf_pe,r_pe(3),sp1_pe,tp1_pe
	real(dl)  :: q_nr0_pe,q_nc0_pe,vr0_pe,vc0_pe,v_nr0_pe,v_nc0_pe,qr0_pe,qc0_pe,N_pe,hn
	real(dl) ee_pe,a_pe,ddf_pe,zz_pe,srsp_pe,k_pe
	real(dl)  , dimension(11)::rs_pe
	real(dl):: coff(25),x0,info(40),z_nn1_pe,z_n1_pe,cs1_pe,cs_n1_pe,h_nn1_pe,phi_nnn1_pe
	
	!real(dl)  , dimension(5000000,20):: psdata
	integer supnum_pe,count2,is_pe
	integer, dimension(5) :: ii
               character(3) rff 
               
    integer li,pain
    
    pain=int(N_pe/hs)
		!phi0=r(1)
		!phi_n0=r(2)
		!h0=r(3)
		a0_pe=coff(8)
		uddf_pe=1e-3
    rs_pe=[phi0_pe,phi_n0_pe,H0_pe,qr0_pe,qc0_pe,q_nr0_pe,q_nc0_pe,vr0_pe,vc0_pe,v_nr0_pe,v_nc0_pe]
		a_pe=a0_pe*exp(N_pe)
		cs1_pe=cs(coff,phi0_pe,phi_n0_pe,H0_pe)
		ddf_pe=k_pe/(a_pe*H0_pe/cs1_pe)
		!hn=1e-4
		supnum_pe=0
		ee_pe=epsln2(coff,phi_n0_pe,H0_pe,phi0_pe)
		!count2_pe=0
		!print **,"Enetering 3.1 ::","at N=",N_pe,"k/(aH)=",ddf_pe
		!OPEN(unit=785,file="../data/ini_pert05.dat" ,action='write',position='append')
		do while ((ddf_pe.ge.1.0).and.(ee_pe.lt.1.0))	
			supnum_pe=supnum_pe+1
            if (ee_pe.ge.1.0) then
					!print **, "spectral : 284 reach ee 1 "
					!print **, "why ?",ee_pe,phi0_pe,phi_n0_pe,H0_pe,coff(16)
				    ! STOP
			endif
			call odesolve2(rs_pe,N_pe,hs,k_pe,coff,info,x0,pain)
			N_pe=N_pe+hs
		    li=li+1    
		    pain=pain+1
			phi0_pe=rs_pe(1)
			phi_n0_pe=rs_pe(2)
			H0_pe=rs_pe(3)
			ee_pe=epsln2(coff,phi_n0_pe,H0_pe,phi0_pe)
			!!print **,"supH:",k, n,h0,phi0,phi_n0

			a_pe=a0_pe*exp(N_pe)
			cs1_pe=cs(coff,phi0_pe,phi_n0_pe,H0_pe)
			ddf_pe=k_pe/(a_pe*H0_pe/cs1_pe)
			qr0_pe=rs_pe(4)
			qc0_pe=rs_pe(5)
			q_nr0_pe=rs_pe(6)
			q_nc0_pe=rs_pe(7)
			vr0_pe=rs_pe(8)
			vc0_pe=rs_pe(9)
			v_nr0_pe=rs_pe(10)
			v_nc0_pe=rs_pe(11)
			!h_nn1_pe=h_nn(coff,phi0,phi_n0,H0)
			!phi_nnn1_pe=phi_nnn(coff,phi0,phi_n0,H0)
			!cs_n1=cs_n(coff,phi0,phi_n0,H0,a0)
			!z_nn1=Z_nn(coff,phi0,phi_n0,H0,a0)
			
			!z_n1=Z_n(coff,phi0,phi_n0,H0,a0)
			zz_pe=zzzz(coff,phi_n0_pe,a_pe,H0_pe,phi0_pe)
            !sP1=((qr0**2/zz**2+qc0**2/zz**2))*k**3.0/2.0/pi**2
            !tP1=((vr0**2/a**2+vc0**2/a**2))*k**3.0/2.0/pi**2
       !    !print **,"3.1: k,N,h,phi,phi_n=",k,N,H0,phi0,phi_n0
       !   !print **,"3.1: z_nn1,z_n1,cs1,cs_n1=",z_nn1,z_n1,cs1,cs_n1
		!	!print **,"sp1,tp1",sp1,tp1
		!	!print **,"qr,qc,z,sp1", qr0,qc0,zz,sp1
			!!print **, "check",j
			!j=j+1
			!ee=epsln2(coff,phi_n0,H0,phi0)
            !read(*,*) rff
         !   write(785,*) n,phi0,phi_n0,h0,ee,zz,z_n1,z_nn1,phi_nnn1,h_nn1,cs1,cs_n1,qr0,qc0,q_nr0,q_nc0
		enddo
		!close(785)
		!Stop
		!print **,"Enetering 3.2 ::","at N=",N_pe,"k/(aH)=",ddf_pe
		srsp_pe=(h0_pe**2/phi_n0_pe**2)/(4.d0*pi**2)
		
		do while ((ddf_pe.ge.uddf_pe).and.(ee_pe.lt.1.0))	
			supnum_pe=supnum_pe+1
			if (ee_pe.ge.1.0) then
					!print **, "spectral : 316 reach ee 1 "
				    ! STOP
			endif
			call odesolve2(rs_pe,N_pe,hs,k_pe,coff,info,x0,pain)
			N_pe=N_pe+hs
		    pain=pain+1
		    li=li+1
			phi0_pe=rs_pe(1)
			phi_n0_pe=rs_pe(2)
			H0_pe=rs_pe(3)
			ee_pe=epsln2(coff,phi_n0_pe,H0_pe,phi0_pe)

		!	!print **,"supH:", k,n,h0,phi0,phi_n0
			a_pe=a0_pe*exp(N_pe)
			cs1_pe=cs(coff,phi0_pe,phi_n0_pe,H0_pe)
			ddf_pe=k_pe/(a_pe*H0_pe/cs1_pe)
			qr0_pe=rs_pe(4)
			qc0_pe=rs_pe(5)
			q_nr0_pe=rs_pe(6)
			q_nc0_pe=rs_pe(7)
			vr0_pe=rs_pe(8)
			vc0_pe=rs_pe(9)
			v_nr0_pe=rs_pe(10)
			v_nc0_pe=rs_pe(11)
			!cs_n1=cs_n(coff,phi0,phi_n0,H0,a0)
			!z_nn1=Z_nn(coff,phi0,phi_n0,H0,a0)
			!cs1=cs(coff,phi0,phi_n0,H0,a0)
			!z_n1=Z_n(coff,phi0,phi_n0,H0,a0)
			zz_pe=zzzz(coff,phi_n0_pe,a_pe,H0_pe,phi0_pe)
           ! sP1=((qr0**2/zz**2+qc0**2/zz**2))*k**3.0/2.0/pi**2
           ! tP1=((vr0**2/a**2+vc0**2/a**2))*k**3.0/2.0/pi**2
        !  !print **,"3.2: k,N,h,phi,phi_n=",k,N,H0,phi0,phi_n0
        !    !print **,"3.2: z_nn1,z_n1,cs1,cs_n1=",z_nn1,z_n1,cs1,cs_n1
		!	!print **,"sp1,tp1",sp1,tp1
		!	!print **,"qr,qc,z,sp1", qr0,qc0,zz,sp1
			!!print **, "check",j
			!j=j+1
			!ee=epsln2(coff,phi_n0,H0,phi0)
           ! read(*,*) rff
		enddo
		!print **,"Completed 3.2 ::","at N=",N_pe,"k/(aH)=",ddf_pe
!		do while (ee.lt.1.0)	
!			supnum=supnum+1
!			
!			if (ee.gt.1.0) then
!					
!				     !print **, "spectral : 342 reach ee 1 "
!				     if (count2.ge.1) then
!				    	STOP
!				    endif
!				    count2=count2+1
!			endif
!			a=a0*exp(N)
!			ddf=k/(a*H0)
!			qr0=rs(4)
!			qc0=rs(5)
!			q_nr0=rs(6)
!			q_nc0=rs(7)
!			vr0=rs(8)
!			vc0=rs(9)
!			v_nr0=rs(10)
!			v_nc0=rs(11)
!			zz=zzzz(coff,phi_n0,a,H0)
!			!ee=epsln2(coff,phi_n0,H0)
!			call odesolve2(rs,N,hs,k,coff)
!			N=N+hs
!			phi0=rs(1)
!			phi_n0=rs(2)
!			H0=rs(3)
!			ee=epsln2(coff,phi_n0,H0)
!		enddo
		!print **, "number of loops in superhorizon limit",supnum_pe
		if (supnum_pe.le.2) then
	print *, "************Superhorizon condition not statisfied*************",supnum_pe

				!print **, "error kind 4:spectral.f perturbation_evolution line 412"
				error_index=4
				error(error_index)=error(error_index)+1
				open (unit = 1586, file = "error.log" ,action='write',position='append')
				write(1586,*) "error kind and count :" ,error_index,error(error_index)
				write(1586,*)  "error kind 4:spectral.f perturbation_evolution line 412"
				write(1586,*) "not satisfied superhorizon condition change a0",supnum_pe
				write(1586,*)x0,(coff(is_pe),is_pe=1,16)
				close(1586)

				!   STOP
			endif
		
		
return
end subroutine

end module spectral
