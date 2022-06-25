program main

!implicit none
!contains
!subroutine get_scalar_pow(V00,phi00,sigma0,coff23,p_in,kk,power)

	use mod1
	use mod2
	use mod3
	use mod4
	use mod5
	!use solve
	use spectral
        use OMP_LIB
	implicit none
	
        real(dl) R,R_end,vdc,svdc,rat,Ms,peak,vc,beta,zai,zai_end,intgn1,nn,vdc1,vdc2,intgn2,intgn,xxx(8),ddd(8),flags(50)
        real(dl) beta2,beta1,srsp,W,dk,k,svdc1,peak1,vc1,rat1,zaic,mass,cont,mass1,dm,eta1,constr,max_m, max_m_final
	real(dl) :: power(500)
	
	real(dl) :: power1(500),xx(500),yy(500),fsrp(500),phi_exit(500),tpp(500)
	real(dl):: kk(500)
	real(dl):: coff(25), info(40)
	real(dl):: x0,b,d,tp,n1,ee,zz,qr,qc,vr,hs,xbar,ybar,top,bot,ampl,slope,sp00,temp,odef,ndef,dot_def0,dot_def1,c
	real(dl):: s1,s2,s3,s4,s5,slope1,scalar_in,k1,phi000,oslp,diff22,max_ps,min_ps
	real(dl) :: ep_pivot,eta_pivot,ns_pivot,hn,aa,phi_n00,h00,n11,arch(3)
	real(dl) :: ep_pivot1,eta_pivot1,ns_pivot1,ns_pivot_alt,dd,dd2,keep11
	real(dl) :: pbhpeak,max_peak,max_mass,max_r
	integer i,j,n,max_i,min_i,jin,deltac,r_in,p_in,loop_in,iii,mx_iii,ij,is
	real(dl):: fsp(50000),ssp(50000),si(50000),ksp(50000)
	real(dl) dgap,oldef,newdef,temp1,oldef0,newdef0,er_count
	character(90) rff,filename
	real(dl) tg,ae,wreh,nreh,greh,fac1,temp2,dgap1,efold,dgapf,dgap0,eff_efold,min_frac,loop_in_sp
	real(dl):: f_PS,f_PSN,f_Peak,f_PeakN,f_PS1,f_PSN1,f_Peak1,f_PeakN1
	real, dimension(2,169) :: array1
	integer ljin,ICOR,iss0
	real(dl) bbb,aaa,ccc,brat
    CHARACTER(100) :: num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,timed
    integer date_time(8),date_time1(8)
    character*10 btime(3)
   
    real(dl) V00,phi00,sigma0,coff23

    
    !To run :   ./primospectra .98 1e-11 20e-11 1.0
    
	iss0=0

	jin=1
	coff(16)=1e-8
	do while (coff(16).le.100.0)
	coff(16)=coff(16)*1.4
	jin=jin+1
   !CALL GET_COMMAND_ARGUMENT(1,num1)  
    !CALL GET_COMMAND_ARGUMENT(2,num2)  
    !   CALL GET_COMMAND_ARGUMENT(3,num3)  
 !      CALL GET_COMMAND_ARGUMENT(4,num4) 
        !CALL GET_COMMAND_ARGUMENT(4,num4) 
       ! CALL GET_COMMAND_ARGUMENT(6,num6) 
      !  CALL GET_COMMAND_ARGUMENT(7,num7) 
      !  CALL GET_COMMAND_ARGUMENT(8,num8) 
      !  CALL GET_COMMAND_ARGUMENT(9,num9) 
      !  CALL GET_COMMAND_ARGUMENT(10,num10)
    !   print *, num1,num3,num5 
     ! coff(1)=phi00                   !then, convert them to REALs
     ! coff(3)=sigma0       !then, convert them to REALs
     !  coff(6)=V00       !then, convert them to REALs
    !coff(23)=coff23
    
   ! coff(16)=1000.0
     coff(6)=4.3252818563344161E-010         !1.9975201487454334E-009
    
    !   READ(num1,*)jin
    !    READ(num2,*)coff(16)
      !  READ(num3,*)coff(6)
!        READ(num4,*)coff(4)
        !READ(num4,*)coff(4)
      ! READ(num8,*)coff(21)
      !  READ(num9,*)coff(22)
      !  READ(num10,*)jin
        
        coff(1)=0.92   ! the value of phi corresponding to wiggles 
         coff(2)=15.0    ! b1 in theta function
	   !  coff(4)=-0.1       ! phase term 'c1' in theta function
       	coff(3)=5.0e-5
         coff(20)=4.0
         coff(21)=0.5
         coff(22)=1.0
         coff(19)=0.0
         coff(25)=4.0     
         coff(18)=350 ! the factor in front of monomial term of gnmdc :: fixed by hand 

         coff(23)=0.0

       call date_and_time(btime(1), btime(2), btime(3), date_time)
    ICOR=OMP_GET_NUM_PROCS()
   ! OPEN(unit=1522,file="../Nidra",action='write',position='append')
   !  write(1522,*) "--------------Primo initiated----------------------------" 
   ! write(1522,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"Date:",&
   !     date_time(3),".",date_time(2),".",date_time(1) &
   !     ,"Time:",date_time(5),".",date_time(6),".",date_time(7),".",date_time(8) &
   !     ,"(hh.mm.ss.milisec)"
   !     write(1522,*) "------------------------------------------------"
    write(*,*) " initiated primo-sin ncore+inputs:",icor,coff(1),coff(3),coff(6)
   ! close(1522)
    

			!print **, "nl1"
	!		open(120, file="constr44")
			! read in values
			!print **, "nl2"
	!		read(120,*) array1
	!		close(120)

        
	call system('mkdir -p ../plot/' )	
	call system('mkdir -p ../data/' )	
	do 274 error_index = 1, 50, 1
         error(error_index)=0
  274  continue		
			
	max_m_final=1e-3
        
	zaic=.414
     !   coff(19)=10.0

  !	jin=50
	deltac=20
	dgap0=.01d0
!       CALL GET_COMMAND_ARGUMENT(1,num1char)  
!   	CALL GET_COMMAND_ARGUMENT(2,num2char)
 ! 	CALL GET_COMMAND_ARGUMENT(3,num3char)  
  !	CALL GET_COMMAND_ARGUMENT(4,num4char)
  !     READ(num1char,*)coff(1)                   !then, convert them to REALs
  !   READ(num2char,*)coff(3)
  !     READ(num3char,*)coff(4)                  !then, convert them to REALs
     !  READ(num4char,*)coff(3)
       


!	coff(6)=3.67e-11      ! the overall coefficient in potential i.e. m^2 here.
	!coff(16)=550.870890562614612e-1    ! (A0 of the theta function *) 
	


	!coff(15)=4.1436729025735670
        coff(15)= 2.6911415955835922          !3.0         !1.1112404375013905!1.05296979904174805                                        !17.2907278468742379026076582001494    ! initial value of phi
	
        coff(17)=dhv00
	
	
               

		temp2=a00
		dgap1=0.4d0
		coff(8)=a00

		

	!        call solveqn(x0,d,c,coff)
		temp1=12.d0
		
		dgap=dgap0
		nreh=0.0
		wreh=1.0/3.0	
		cont=100.0
		constr=.3d0
		min_frac=1.0e-15
		ij=0

		ljin=1
  
                       
344		error(25)=0
                call get_pbh(x0,coff,power,power1,kk,zaic,temp2,dgap1,nreh,wreh,j,jin,deltac,tpp,efold,info,ljin)
     
        if(error(25).gt.0)then
               
 
      !  OPEN(unit=1522,file="../Nidra",action='write',position='append')
        if((iss0.lt.1).and.(abs(coff(24)).gt.1.2))then
               brat=abs(abs(coff(24))-0.15)
               !coff(16)=coff(16)/(abs(brat))
                coff(18)=coff(18)/(abs(brat))
       !        write(1522,*) "Nan-Negetive Err: in brat:",iss0,coff(24),brat,coff(16)
               
        else
                !coff(16)=coff(16)*0.9
                coff(18)=coff(18)*0.9
           !write(*,*) "Nan-Negetive Err: in else:",iss0,coff(24),brat,coff(16)
            write(*,*) "Nan-Negetive Err: in else:",iss0,coff(24),brat,coff(18)
        endif 
             iss0=iss0+1   
        !write(1522,*) "executed: NaN negetive error is there.: ratio:",coff(24),brat
             !write(1522,*) "new/nextvalue of coff(16) is " ,coff(16)    
            write(*,*) "new/nextvalue of coff(16) is " ,coff(16)
         !   close(1522)
                go to 344
       endif



        90		p_in=j
		WRITE(filename,'(a,i4.4,a)') "../data/final_power",jin,".TXT"
	 		OPEN(unit=152,file=filename ,action='write',position='append')
                        write(152,*)"#",(coff(i),i=1,23)
	 				i=1
       					do while(i.le.p_in-1)
					write(152,*) x0,kk(i),power(i),tpp(i),power1(i),tpp(i)/power(i)
					i=i+1
					end do
			close(152)	

		open (unit = 100, file = "../coff.dat",action='write',position='append')
	        Write(100,*) (coff(i),i=1,23)
                close(100)

             !  open (unit = 111, file = "../data/ifo.dat",action='write',position='append')
		

		
		
		filename="../Nidra-human-readable"
         Call write_info(x0,coff,info,filename)

			
	        er_count=0
		do 170 i = 1, 50
                er_count=er_count+error(i)
170		continue

			open (unit = 1586, file = "error.log" ,action='write',position='append')
			write(1586,*) "Total error count=",er_count
			write(1586,*) "Description of errors::::\n"
			write(1586,*) ((i),i=1,15)
			close(1586)
         call date_and_time(btime(1), btime(2), btime(3), date_time1)
              !  OPEN(unit=1522,file="../Nidra",action='write',position='append')
              !   write(*,*) "------------------------------------------------"
              write(*,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"start-Date:",&
        date_time(3),".",date_time(2),".",date_time(1) &
        ,"Time:",date_time(5),".",date_time(6),".",date_time(7),".",date_time(8) &
        ,"(hh.mm.ss.milisec)"
write(*,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"end-Date:",&
        date_time1(3),".",date_time1(2),".",date_time1(1) &
        ,"Time:",date_time1(5),".",date_time1(6),".",date_time1(7),".",date_time1(8) &
        ,"(hh.mm.ss.milisec)"
          !    write(1522,*) "coff(1),coff(3),coff(15),coff(16),code(17),coff(6)=="
         !       write(1522,*) coff(1),coff(3),coff(15),coff(16),coff(17),coff(6)
           !     write(1522,*) "=======================Primordial code completed==========================="
          !      close(1522)

enddo
end program main






subroutine get_pbh(x0,coff,power,power1,kk,zaic,temp2,dgap1,nreh,wreh,j,jin,deltac,tpp,efold,info,ljin)

	use mod1
	use mod2
	use mod3
	use mod4
	use mod5
	!use solve
	use spectral
	implicit none


      real(dl) R,R_end,vdc,svdc,rat,Ms,peak,vc,beta,zai,zai_end,intgn1,nn,vdc1,vdc2,intgn2,intgn,k_infl,k_ini
        real(dl) beta2,beta1,srsp,W,dk,k,svdc1,peak1,vc1,rat1,zaic,mass,cont,mass1,dm,eta1
	real(dl) :: power(500)
	real(dl) :: power1(500),xx(500),yy(500),fsrp(500),phi_exit(500),tpp(500),err_16_mem
	real(dl):: kk(500)
	real(dl):: coff(25),info(40)
	 integer is
	real(dl):: x0,b,d,tp,n1,ee,zz,qr,qc,vr,hs,xbar,ybar,top,bot,ampl,slope,sp00,temp,odef,ndef,keep25,dk98,ov
	real(dl):: s1,s2,s3,s4,s5,slope1,scalar_in,k1,phi000,oslp,diff22,max_ps,min_ps
	real(dl) :: ep_pivot,eta_pivot,ns_pivot,hn,aa,phi_n00,phi00,h00,n11,arch(3)
	real(dl) :: ep_pivot1,eta_pivot1,ns_pivot1,ns_pivot_alt
	real(dl) :: pbhpeak,max_peak,max_mass,max_r
	integer i,j,n,max_i,min_i,jin,deltac,r_in,p_in,ljin
	real(dl):: fsp(50000),ssp(50000),si(50000),ksp(50000)
	real(dl) dgap,oldef,newdef,temp1,oldef0,newdef0
	character(90) rff,filename
	real(dl) tg,ae,wreh,nreh,greh,fac1,temp2,dgap1,efold,efold_pivot
	real(dl):: f_PS,f_PSN,f_Peak,f_PeakN,f_PS1,f_PSN1,f_Peak1,f_PeakN1
	real(dl) co6star,spstar
	real(dl) detnorm,co6p,co6n,sp
	integer aloopi
         integer date_time(8),date_time1(8)
         character*10 btime(3)
    real(dl) sp1,sp2,sp3,sp4, tp1,tp2,tp3,tp4,tp_to_sp,scalar_index

  call date_and_time(btime(1), btime(2), btime(3), date_time)
              write(*,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"start-getpbh:",&
        date_time(3),".",date_time(2),".",date_time(1) &
        ,"Time:",date_time(5),".",date_time(6),".",date_time(7),".",date_time(8) &
        ,"(hh.mm.ss.milisec)"



       ! error(25)=0
	!ljin=ljin+1
	!coff(4)=2.6836816975470263
	hs=1e-4
 ! 	write(*,*)"After solve the loop",x0,(coff(i),i=1,6)
	!coff(6)=2.9d-16
	!print **,"coff6",coff(6)


85	i=1
       ! coff(6)=2.8e-10
	phi000=coff(15)
	!go to 1235
242	coff(5)=10.11542896334306045247
243     call get_N(coff,efold,info,x0,jin) !	call set_efolds(coff,efold,info,x0)
        print *," efold=",efold, "coff(15)=",coff(15),"coff(10)=",coff(10),"coff(16)=",coff(16),nef,"coff(18)=",coff(18)
        
        
        if(error(25).gt.0)then
        
                print *, "NaN negetive error is there"
                go to 570
        endif
        
        if((efold.gt.(nef+1.0)).or.(efold.le.(nef-1.0)))then
        print *, "in set efolds"
         call set_efolds(coff,efold,info,x0,15)
         print *, coff(15),efold
                go to 243
        endif


       ! coff(16)=coff(16)*coff(23)
	
	!###########################################################################################
  
    
   ! go to 859 ! going directly to pwr
    
    
    !########################################################################################### 
	
	
	
1248	aloopi=0
	i=1
	dgap1=.4

	
258	sp00=2.09052*1e-9
	coff(6)=20.0e-11                     !1.8951315855878927E-009
	dk=.1
	call  run_back_assignvalues(coff,1e-4_dl)
	call get_spectra(.05_dl,tp,sp,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
	print *,"sp at .o5",sp
	if(abs(sp-sp00)/sp00.le.1e-2)then
	Print *,"skipping tp going to Normalization and initial scale factor optimization"
	  go to 410
	endif  
	coff(6)=((5.5*1e-5)**2)/sp*1.8951315855878927E-009

222	call  run_back_assignvalues(coff,1e-4_dl)
   call get_spectra(.05_dl,tp,sp,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
!###########################################################################################
  ! k_infl=1e15
  !  go to 859 ! going directly to pwr
    
    
    !########################################################################################### 
        detnorm=sp-sp00
        if(detnorm.le.0.d0)then
        co6n=coff(6)
    589  coff(6)=coff(6)*1.5d0
        call  run_back_assignvalues(coff,1e-4_dl)
        call get_spectra(.05_dl,tp,sp,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
        detnorm=sp-sp00
        if(detnorm.le.0.d0)then
            co6n=coff(6)
            go to 589
            else
            co6p=coff(6)
            endif
        else
        co6p=coff(6)
    588  coff(6)=coff(6)/1.5d0
        call  run_back_assignvalues(coff,1e-4_dl)
        call get_spectra(.05_dl,tp,sp,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
        detnorm=sp-sp00
        if(detnorm.ge.0.d0)then
            co6p=coff(6)
            go to 588
            else
            co6n=coff(6)
            endif
        endif
        
        

        co6star=coff(6)
        spstar=sp
            do while (abs(spstar-sp00)/sp00.ge.1e-2)
            print *, "before starting loop_normal : ", co6star,spstar,co6star,spstar,co6n,co6p
          !  call loop_normal(co6star,spstar,co6p,co6n,info,x0,ljin,jin,coff)
            call loop_normal(co6star,spstar,co6p,co6n,info,x0,ljin,jin,coff,sp00)
            print *, "after loop normal : ",co6star,spstar,co6n,co6p
            !read(*,*)  rff
            enddo
     !################################################################################################	       
            coff(6)=co6star
            sp=spstar
   ! print *, ")))))))))))))))))))))))))))))))))))))))))))))loopcheck2((((((((((((((((((((((((((((((((((((((((((((((((((("
        !call get_spectra(.05d0,tp,sp,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff)
   410     print *,"sp at .o5",sp, "for cof6=", coff(6)	
	
       !##################################################################
        ! Note that there are mo normalization loop here !!!!

       !##################################################################
	WRITE(filename,'(a,i4.4,a,i4.4,a)') '../data/reheating_d_check',jin,"ln",ljin,".TXT"
	open (unit = 1144, file = filename ,action='write',position='append')
        write(1144,*) "--------------New Case----------------------"
	write(1144,*) x0,(coff(i),i=1,16)
        write(1144,*) "------------------------------------"
	write(1144,*) "||",coff(15),efold,n11,h00,ae,coff(8)
	close(1144)	
	!###########################################################################################################################
	! Reheating Part
	call touch_initial(x0,coff,phi00,phi_n00,H00)
	k_ini=coff(8)*H00*1e3
        	call rates_ini(coff,info,p_in,jin,x0,n11,h00,ljin,efold_pivot)
     !   n11=efold
      info(36)=h00
  !  print *, "after normalization efold number=", n11, "H_end=",h00
	k_infl=coff(8)*H00*Exp(n11)*1e-4
	
	info(33)=efold_pivot
!###################################################################################
        if((efold-n11)/efold.ge.0.10d0)then

			Write(*,*) "efold,n11 values dont match",efold,n11
                        !print **, "error kind 5 : main.f get_pbh line 394"
			error_index=5
			error(error_index)=error(error_index)+1
			open (unit = 1586, file = "error.log" ,action='write',position='append')
			write(1586,*) "error kind and count :" ,error_index,error(error_index)
			write(1586,*) "error kind 5 : main.f get_pbh line 394"
			write(1586,*)"efold,n11 values dont match",efold,n11
			write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
			close(1586)

 			read(*,*) info(25)
	else
         !   Write(*,*) "efold,n11 values do match",efold,n11
 		!	read(*,*) info(25)
 		n11=efold
        endif
    !    k_infl=coff(8)*H00*Exp(n11)*1e-4

		!close(25)
		!print **, "H0",H00,"N",n11
		!H00=H00!*(1.22*(1.0e19))**2
289		Tg = 2.7255*8.61733*(1e-14)*8.18*(1e-20)
		greh = 106.75
		fac1=(3.d0*(H00**2)*Exp(-3.d0*nreh*(1.d0 +wreh)))**(-1.0/4.0)*1.90776*(10**(57.0-30.004))/Exp(n11-70.0)
		ae=((43.0/11.d0/greh)**(1.d0/3.d0))*(((Pi**2)*greh/30.d0)**(1.0/4.0))*Tg*fac1*Exp(-nreh)
	!	print *,"ae=",ae,"where a0=",coff(8)



                !################################## initial scaling of a0 ####################
	        if(((abs(ae/coff(8)).ge.10.d0).or.(abs(ae/coff(8)).lt.0.1d0)).and.(aloopi.le.0))then
                   coff(8)=ae*1.2
                   aloopi=aloopi+1
                   go to 222
                endif
	        !###############################################################################
		
		oldef0=coff(8)-temp2
		newdef0=coff(8)-ae
		temp2=ae
			if(abs(oldef0+newdef0).eq.(abs(oldef0)+abs(newdef0))) then
			!		print *, "ok"
					dgap1=dgap1
			else
			!	print *, "before taking 1/2 th dgap",dgap 
				dgap1=dgap1/2.d0
			!	print *, "after taking 1/2 th dgap",dgap 
			endif
	if(ae.ge.(coff(8)*1.1))then
		   coff(8)=coff(8)*(1.d0+dgap1)
		   go to 222
	else if(ae.le.(coff(8)*0.9))then
		   coff(8)=coff(8)*(1.d0-dgap1)
		   go to 222
	endif  

	!############################################################################################################################
338	i=1
    j=1
	!call infl_evolve(coff,phi00,H00,Phi_n00,n11,efold,info,x0)
	
go to 859
!#############################################################################################
	if (info(26).ge.(4e-5))then
				        !print **, "error kind 13 (Warning) : main.f get_pbh line 549 "
					error_index=13
					error(error_index)=error(error_index)+1
					open (unit = 1586, file = "error.log" ,action='write',position='append')
					write(1586,*) "error kind and count :" ,error_index,error(error_index)
					write(1586,*) "no inflection feature: info(26)=min_phi_n=",info(26)
					write(1586,*) "info(27)=max_pow=",info(27)
					write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
					close(1586)
	coff(12)=0.d0
	!STOP
	go to 570
	end if

        if (info(26).le.(1.9e-5))then
                   !print **, "error kind 13_alt (Warning) :main.f get_pbh line 733"
                   open (unit = 1586, file = "error.log",action='write',position='append')
                                       
                    write(1586,*) "over inflection feature:phi_n_min",info(26)
                    write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
                    close(1586)
                    coff(12)=10.d0+info(27)
        !STOP
                    go to 570
        end if

!#############################################################################################

!#############################################################################################
	if ((info(26).le.(4e-5)).and.(info(29).le.(1e-22)))then
				        !print **, "error kind 14 (Warning) : main.f get_pbh line 565 "
					error_index=14
					error(error_index)=error(error_index)+1
					open (unit = 1586, file = "error.log" ,action='write',position='append')
					write(1586,*) "error kind and count :" ,error_index,error(error_index)
					write(1586,*) "inflection feature is far: info(26)=min_phi_n=",info(26)
					write(1586,*) "info(27)=max_pow=",info(27)
					write(1586,*) "info(28)=k_maxx=",info(28)
					write(1586,*) "info(29)=massn_max=",info(29)
					write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
					close(1586)
	coff(12)=0.d0
	!coff(15)=coff(15)-0.01d0
	!STOP
	go to 570
	end if
!#############################################################################################
	err_16_mem=error(16)
  !  Write(*,*) "coffs=",(coff(i),i=1,16)
 859 call date_and_time(btime(1), btime(2), btime(3), date_time)
              write(*,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"start-pwr:",&
        date_time(3),".",date_time(2),".",date_time(1) &
        ,"Time:",date_time(5),".",date_time(6),".",date_time(7),".",date_time(8) &
        ,"(hh.mm.ss.milisec)"
   call  run_back_assignvalues(coff,1e-4_dl)
   
   call date_and_time(btime(1), btime(2), btime(3), date_time)
              write(*,'(A10,i02.2,1(A1,i02.2),A1,i4,A10,i02.2,2(A1,i02.2),A1,i3,A20)')"startpwr1:",&
        date_time(3),".",date_time(2),".",date_time(1) &
        ,"Time:",date_time(5),".",date_time(6),".",date_time(7),".",date_time(8) &
        ,"(hh.mm.ss.milisec)"
    call pwr(kk,power,power1,phi_exit,fsrp,j,coff,tpp,k_infl,info,x0,ljin,jin,k_ini)

    call get_spectra(.05_dl,tp1,sp1,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
    call get_spectra(.05_dl+1e-4,tp2,sp2,N1,ee,eta1,zz,qr,qc,vr,vc,hs,srsp,coff,info,x0,ljin,jin)
    tp_to_sp=tp1/sp1
    scalar_index=1.0+(log(sp2)-log(sp1))/(log(.05_dl+1e-4)-log(.05_dl))
    print *, "scalar index=",scalar_index, "tp_to_sp=", tp_to_sp
	
	 OPEN(unit=1522,file="../ns_vs_r.dat",action='write',position='append')
	  write(1522,*) "#coff(16),coff(25),coff(6),coff(15),tp_to_sp,scalar_index, P_R_.05, P_h_.05"
	 write(1522,*) coff(16),coff(25),coff(6),coff(15),scalar_index,tp_to_sp,sp1,tp1
	close(1522)
	
	  570 return

end subroutine






subroutine rates_ini(coff,info,p_in,jin,x0,n11s,h00s,ljin,efold_pivot)
	use mod1
	use mod2
	use mod3
	use mod4
	use mod5
	!use solve
	use spectral
	implicit none
	
	
	


	
	real(dl) a_eq,gdof,mp,gev,tg
    integer  i,is,jin,p_in,ljin
     real(dl):: coff(25), info(40),x0
     
    real(dl) phi_n00s,phi00s,h00s,archs(3),n11s,aas,ee_s
    character(90) rff,filename
    real(dl) :: ep_pivot,eta_pivot,ns_pivot,ns_pivot_alt,massn_max,a_form1,hn,ns_pivot1
	real(dl) ep002,eta002,eta05,ep05,k_maxx,n_s002,n_ss002,n_t002,tp_to_sp_ratio002,efold_pivot
	
	real(dl) max_pow,min_phi_n,n01,phi01,h01,min_ee_s
	real(dl) epV,etaV,Anmdc,ns_nmdc,ns_nmdc05,ns_nmdc002
	
	
	a_eq=1.d0/3300.d0
	gdof=106.75d0
	mp=1.22*1e19
	gev=2.488767d0*1e37*2.d0*Pi
	tg=2.7255*8.61733*1e-14/mp
	



!###############################################################################################################	
!   scalar index part	
	i=1
	hn=1d-3
	
	
	
	
			call set_initial(x0,coff,phi00s,phi_n00s,H00s)

			archs=[phi00s,phi_n00s,H00s]
			N11s=0.d0
			aas=coff(8)
			ee_s=epsln2(coff,phi_n00s,H00s,phi00s)
			!WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/indexes01_slow_",jin,"ln",ljin,".dat"
			!open (unit = 25, file = filename ,action='write',position='append')	
			!WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/back",jin,"ln",ljin,".dat"
			!open (unit = 106, file =filename  ,action='write',position='append')
			min_phi_n=2.1d0
			min_ee_s=1.d0
			do while (ee_s.le.1.d0)


!#######################################################################################################
		if(n11s.ge.(nef*1.5d0))then
			Print *,"a problem : rates-ini  :to procee_sd say 'yes'"

                        !print **, "error kind 12: main.f rates line 673"
			error_index=12
			error(error_index)=error(error_index)+1
			open (unit = 1586, file = "error.log" ,action='write',position='append')
			write(1586,*) "error kind and count :" ,error_index,error(error_index)
			write(1586,*)  "error kind 7 : main.f rates line 591"
			write(1586,*)"efold no n11s> (nef*1.5d0) : n11s=",n11s
			write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
			close(1586)


			read (*,*)rff
						if (rff=="yes")then
							!print **,"fine then"
						else
							STOP
						endif
		 endif	
			!!print **,'n11s11,r,ee_s', N,r,ee_s
!######################################################################################################	

			
				call odesolve(archs,N11s,hn,coff,info,x0)
				phi00s=archs(1)
				phi_n00s=archs(2)
				H00s=archs(3)
				n11s=n11s+hn
				aas=coff(8)*exp(n11s)
				ee_s=epsln2(coff,phi_n00s,H00s,phi00s)
				!fsp(i)=((phi_n00s/2.d0)**2)*2.d0
				!ssp(i)=eta(coff,phi00s,phi_n00s,H00s)
				!ksp(i)=aas*h00s
				!si(i)=1.d0-4.d0*(fsp(i))+2*(ssp(i))
				ep_pivot=((phi_n00s/2.d0)**2)*2.d0
				eta_pivot=eta(coff,phi00s,phi_n00s,H00s)
				ns_pivot=1.d0-4.d0*ep_pivot+2.d0*eta_pivot
				ns_pivot1=1.d0-4.d0*abs(ep_pivot)+2.0*abs(eta_pivot)
				ns_pivot_alt=1.d0-2.d0*ep_pivot+eta_pivot
                epV=0.5d0*(V_phi(coff,phi00s)/V(coff,phi00s))**2 
                etaV=V_phi_phi(coff,phi00s)/V(coff,phi00s)
                Anmdc=1.d0+3.d0*(H00s**2)*theta(coff,phi00s)  
				ns_nmdc=1.d0-(2.d0*epV*(4.d0-1.d0/Anmdc)-2.d0*etaV)/Anmdc
				
				!if(aas*H00s.ge.1e4)then
					if(abs(ee_s).le.min_ee_s)then		
						min_phi_n=abs(phi_n00s)
						max_pow=(h00s**2/phi_n00s**2)/(4.d0*pi**2)
						min_ee_s=abs(ee_s)
						k_maxx=aas*h00s
						phi01=phi00s
						h01=h00s
						n01=n11s
					endif
				!endif
					if((aas*H00s).le.0.002)then
						eta002=eta_pivot
						ep002=ep_pivot
						ns_nmdc002=ns_nmdc
					endif
					if((aas*H00s).le.0.05)then
						eta05=eta_pivot
						ep05=ep_pivot
						efold_pivot=n11s
						ns_nmdc05=ns_nmdc
					!	write(25,*) aas*h00s,ep_pivot,eta_pivot,ns_pivot
					endif
		
			        ! write(106,*)aas*h00s,n11s,phi00s,h00s,phi_n00s,V(coff,phi00s)
			enddo
			!close(106)
			!close(25)
		

	!print **,"pivot slow roll parametrs",ep_pivot,eta_pivot,aas*h00s,n11s
	a_form1=(.07*.14/k_maxx)*a_eq
	massn_max=.2*4.0*Pi*(a_form1/k_maxx)*gev*mp**2*1.602*1e-10/(3.d0*1e8)**2/2.d0/1e30
	info(1)=ns_nmdc05
	info(2)=ep05
	info(3)=eta05
	info(19)=ns_nmdc002
	info(20)=ep002
	info(21)=eta002
	info(26)=min_phi_n
	info(27)=max_pow
	info(28)=k_maxx
	info(29)=massn_max
	info(30)=phi01
	info(31)=h01
	info(32)=n01
end subroutine


	



subroutine pwr(kpara,powerw,power1w,phi_exitw,fsrpw,j,coff,tppw,k_infl,info,x0,ljin,jin,k_ini)


	use mod1
	use spectral
        use OMP_LIB        
	implicit none


!##########################################################################################################
!k_ini=the lowest value of k,k_end=the highest value of k, sp=scaler power spectra 
!,tp=tensor power spectra,coff(15)=arbitary initial value ,nef=no of efolds requirement for inflation
!ns=scalar spectral index,nt=tensor spectral index,phi=scalar field, phi_n=1st dervative of scalar wrt efolds
!H-hubble parameter , N=efolds,
!##########################################################################################################
!module 1 specifies the potential function changing it we can shift to diffrent models
!in first part taking any arbitary value of coff(15) will lead to find a particular value for required efolds
!then we evaluate the scalar and tensor power specta along with spectral indexes and save the data to files
!###########################################################################################################





   ! real(dl) sp,tp,sd,ns,nt,k,dk
	!real(dl) N1,ee,zz,qr,qc,vr,vc,hs,srsp,phi00,phi000,eta1,spmem
	!,i
	
	
	
	
	
	
	
	
	
	integer ::i,j,ljin,jin
	character(90) rff,filename
	real(dl):: coff(25),k_infl,info(40),x0,k_ini
	real(dl):: tp11,sp11,srsp11,k11,kstar,r11
	real(dl) powerw(500),fsrpw(500),ssrpw(500),phi_exitw(500),tppw(500)
	real(dl) :: power1w(500),nlast(500)
	
	
	integer ipara,iparacount,is
	real(dl) kpara(500),kpara1
	
	real(dl) kpara2,tpw,spw,N1w,eew,eta1w,zzw,qrw,qcw,vrw,vcw,hs,srspw
	
	real(dl) pricof(25)
	!pragma omp parallel
	kstar=info(28)*1e-2
	

	!phi000=coff(15)
	!call find_phi_ini(phi000,nef,coff)
	!!print **,"starting value of phi:",phi000
	kpara1=1e-6

        k_infl=10.0
	!phi00=phi000
	!i=1
	
	!j=1
     !      k_infl=1.0e1
	!open (unit = 12, file = "logforplot.dat",action='write',position='append')
	ipara=1
	hs=1e-4
	
	do while (kpara1.le.k_infl)
			
		if ((kpara1.ge.5e-3).and.(kpara1.le.0.1)) then	
	!		!dk=k*0.1
			kpara(ipara)=kpara1
	!		!print **,ipara, kpara(ipara)
			ipara=ipara+1
			
			kpara1=kpara1*1.2d0
			
	        
			!k=k*1.05
		
            else if ((kpara1.ge.5e-4).and.(kpara1.le.5e-3)) then
			!dk=k*0.4
            kpara(ipara)=kpara1
            !print **,ipara, kpara(ipara)
			ipara=ipara+1
			

			kpara1=kpara1*1.07d0
			
		!	!print **,"u"
		else if((kpara1.le.5e-4).and. (kpara1.ge.5e-5)) then
			!dk=k*0.4
            kpara(ipara)=kpara1
            !print **,ipara, kpara(ipara)
			ipara=ipara+1
			

			kpara1=kpara1*1.07d0
			
		!	!print **,"u"
		
        else
            			!dk=k*0.4
            kpara(ipara)=kpara1
            !print **,ipara, kpara(ipara)
			ipara=ipara+1
			

			kpara1=kpara1*1.4d0
        
			
		endif
		
	end do
	iparacount=ipara
	j=iparacount-1
!!	CALL OMP_SET_NUM_THREADS(48)
	!OMP parallel private(kpara2,tpw,spw,N1w,eta1w,zzw,qrw,qcw,vrw,vcw,srspw)
	
    !$OMP PARALLEL DO private(kpara2,tpw,spw,N1w,eta1w,zzw,qrw,qcw,vrw,vcw,srspw,pricof)
    !Ipara is private by default
    DO Ipara=1,iparacount-1
            kpara2=kpara(ipara)
            !print **, ipara,kpara(ipara)
            pricof=coff
            !pricof(4)=1e-5
            call get_spectra(kpara(ipara),tpw,spw,N1w,eew,eta1w,zzw,qrw,qcw,vrw,vcw,hs,srspw,pricof,info,x0,ljin,jin)
            !kkw(ipara)=kpara2
            powerw(ipara)=spw
            power1w(ipara)=srspw
            fsrpw(ipara)=eew
            phi_exitw(ipara)=coff(7)
            tppw(ipara)=tpw
            nlast(ipara)=N1w
          print *, kpara(ipara),spw, tpw, tpw/spw, OMP_GET_THREAD_NUM()
            !i=i+1
            
            !spmem=sp
            if(spw.ge.0.1d0)then
                
                Print *,"a problem in pwr scalar power spec > 0.10  :"

                            !print **, "error kind 16 : main.f rates line 591"
                error_index=16
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
                write(1586,*)  "error kind 7 : pwr, main.f90  line 1280"
                write(1586,*)"sp> 0.1 : sp=",spw
                write(1586,*)x0,(coff(is),is=1,16),(info(is),is=1,25)
                close(1586)
                !exit
            endif
        ENDDO
        !$OMP END PARALLEL DO
        !OMP end parallel
      !  WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/power",jin,".",ljin,".dat"
      !  open (unit = 11, file = filename,action='write',position='append')
        DO Ipara=1,iparacount-1
        k11=kpara(ipara)
        sp11=powerw(ipara)
        srsp11=power1w(ipara)
        tp11=tppw(ipara)
        r11=tp11/sp11
       write (*,*) k11,sp11,tp11,srsp11,r11,nlast(ipara)
      ! write (11,*) k11,sp11,tp11,srsp11,r11
        if(sp11 /=sp11)then
                !sp11=0.0
                Print *,"a problem in pwr scalar power spec =NaN  :"

                            !print **, "error kind 16 : main.f rates line 591"
                error_index=21
                error(error_index)=error(error_index)+1
                open (unit = 1586, file = "error.log" ,action='write',position='append')
                write(1586,*) "error kind and count :" ,error_index,error(error_index)
                write(1586,*)  "error kind 21 : pwr, main.f90  line 1421"
                write(1586,*)"sp=NaN : sp=",sp11
                write(1586,*) k11,sp11,srsp11,(coff(is),is=1,8)
                close(1586)
                
                
        endif
      !  write (*,*) k11,sp11
        ENDDO
      !  close (11)
        !close(12)
        !close(13)
        !close(14)
        
	
	
	end subroutine



subroutine loop_normal(co6star,spstar,co6p,co6n,info,x0,ljin,jin,coff,sp00)


	use mod1
	use spectral	
        use OMP_LIB
	implicit none



	
	
	
	
	
	
	integer ::i,j,ljin,jin
	character(90) rff,filename
	real(dl):: coff(25),k_infl,info(40),x0
	real(dl):: tp11,sp11,srsp11,k11,kstar
	real(dl) powerw(500),fsrpw(500),ssrpw(500),phi_exitw(500),tppw(500)
	real(dl) :: power1w(500)
	real(dl) co6star,spstar,rati,sp00
	
	integer ipara,iparacount

	real(dl) copara(50),kpara2,diff(50),dry,diffmin,diffp,diffn,co6p,co6n
	
	real(dl) tpw,spw,N1w,eew,eta1w,zzw,qrw,qcw,vrw,vcw,hs,srspw
	
	real(dl) pricof(25)
    real(dl) co6i,co6
	real(dl) sp05,icore
	!pragma omp parallel
		


	icore=OMP_GET_NUM_PROCS()
      !  print *, "number of available threads",icore
	sp05=sp00
	hs=1e-4
	dry=(co6p-co6n)/(icore-1.0)
	co6=co6n   
	ipara=1
	!print **, "initiating loop normal"
    do while (co6.le.co6p)
                    copara(ipara)=co6
                    !print **,ipara,copara(ipara)
                    ipara=ipara+1
                    
                    co6=co6+dry
    end do
    copara(ipara)=co6p
	!open (unit = 12, file = "logforplot.dat",action='write',position='append')
	!ipara=1	

	
	

	iparacount=ipara
	!!print **, "lets try="
	!pricof=coff
	!!print **, pricof,x0,jin,ljin
	!call get_spectra(0.05d0,tpw,spw,N1w,eew,eta1w,zzw,qrw,qcw,vrw,vcw,hs,srspw,pricof,info,x0,ljin,jin)
	!!print **,"sp=",spw
	!read(*,*) rff
	!!CALL OMP_SET_NUM_THREADS(48)
	!OMP parallel private(kpara2,tpw,spw,N1w,eta1w,zzw,qrw,qcw,vrw,vcw,srspw,pricof)
	
    !$OMP PARALLEL DO private(kpara2,tpw,spw,N1w,eta1w,zzw,qrw,qcw,vrw,vcw,srspw,pricof)
    !Ipara is private by default
    DO Ipara=1,iparacount
            
            !print **, ipara,copara(ipara)
            pricof=coff
            pricof(6)=copara(ipara)
            call  run_back_assignvalues(coff,1e-4_dl)
            call get_spectra(0.05_dl,tpw,spw,N1w,eew,eta1w,zzw,qrw,qcw,vrw,vcw,hs,srspw,pricof,info,x0,ljin,jin)
            !kkw(ipara)=kpara2
            powerw(ipara)=spw
            power1w(ipara)=srspw
            fsrpw(ipara)=eew
            phi_exitw(ipara)=coff(7)
            tppw(ipara)=tpw
            diff(ipara)=spw-sp05
            !!print **, kk(i),power(i)
            !i=i+1
            

        ENDDO
        !$OMP END PARALLEL DO
        !OMP end parallel
        diffmin=1.d0
        diffp=1.d0
        diffn=-1.d0
        DO Ipara=1,iparacount
          if(diffmin.ge.abs(diff(ipara)))then
              diffmin=abs(diff(ipara))
              co6star=copara(ipara)
              spstar=powerw(ipara)
          endif
         if((diff(ipara).ge.0.d0).and.(diffp.ge.diff(ipara)))then
              diffp=diff(ipara)
              co6p=copara(ipara)
          endif
        if((diff(ipara).le.0.d0).and.(diffn.le.diff(ipara)))then
              diffn=diff(ipara)
              co6n=copara(ipara)
          endif
          !print **, copara(ipara),powerw(ipara)
        ENDDO
        
	!print **, "concluding loop normal"
	
	end subroutine














subroutine get_spectra(k_gs,tp_gs,sp_gs,N_gs,ee,eta1,zz_gs,qr0_gs,qc0_gs,vr0_gs,vc0_gs,hs,srsp_gs, &
                coff,info,x0,ljin,jin)


!######################################################################################################
!M_pl=plank mass,a0=initial scale factor ,coff(15)=starting value of scalar field, 
!H=hubble parameter ,phi_n=1st dervative of scalar wrt efolds, N=efolds, phi= scalar field ,
!H0=initial H,pi=pi
!######################################################################################################
!this subroutine evaluates the scalar and tensor power spectra calling different subroutines 
!for diffterent conditions and also provides various other parameters at super-horizon limit 
!######################################################################################################




	use mod1
	use mod2
	use mod3
	use mod4
	use mod5
	use spectral
	implicit none

	real(dl) :: a0,hs,hn
	real(dl)  :: q_nr0_gs,q_nc0_gs,vr0_gs,vc0_gs,v_nr0_gs,v_nc0_gs,k_gs,qr0_gs,qc0_gs,N_gs,zz_gs
	
	real(dl) ee,eta1
	real(dl) :: tp_gs,sp_gs,a_gs,srsp_gs,H0_gs,phi0_gs,phi_n0_gs
	real(dl)  , dimension(11)::rs_gs
	real(dl) , dimension(3)::r_gs
	real(dl):: coff(25),info(40),x0
    integer is,i,j,l,ljin,jin
	!integer, dimension(5) :: ii
		integer jinn_gs
	character(90) rff1,filename1,filename
	
	integer li
	
	!!print **, "in get_spectra"
	
	jinn_gs=int(abs(log(k_gs)*1000))
	!jinn=jinn+1
	!print **,"log10(k)=", log10(k_gs)
	!WRITE(filename1,'(a,i8.4,a)') "scalardata",jinn,".dat"
 	!OPEN(unit=jinn,file=filename1 ,action='write',position='append')
    
    li=1
	i=1
	j=1
	!k=1
	l=1
	!PI=4.D0*DATAN(1.D0)
	!M_pl=1.d0
	a0=coff(8)
	
	call touch_initial(x0,coff,phi0_gs,phi_n0_gs,H0_gs)

	
	
	
	!print **, "line 1591: get_spectra :",a0,phi0_gs,phi_n0_gs,H0_gs
	r_gs=[phi0_gs,phi_n0_gs,H0_gs]
!	do 10 i = 1, 3
	!!print **, "check",i
  	!psdata(j,i)=r(i)
 ! 10  continue
  	
	hn=1E-3
	N_gs=0.d0
	!psdata(j,13)=N
	!psdata(j,18)=k
	call subhorizon_evolution(k_gs, phi0_gs,phi_n0_gs,H0_gs,N_gs,hn,coff,info,x0,li)
	!!print **, "check",j
	r_gs=[phi0_gs,phi_n0_gs,H0_gs]
	a_gs=a0*exp(N_gs)
	call perturbation_initiation(coff,qr0_gs,qc0_gs,q_nr0_gs,q_nc0_gs,vr0_gs,vc0_gs,v_nr0_gs,v_nc0_gs,k_gs,a_gs,H0_gs, &
	                           phi0_gs,phi_n0_gs)
	r_gs=[phi0_gs,phi_n0_gs,H0_gs]
	!!print **,"======#######=====at N H",H0,N
  call perturbation_evolution(phi0_gs,phi_n0_gs,H0_gs,qr0_gs,qc0_gs,q_nr0_gs,q_nc0_gs,vr0_gs,vc0_gs,v_nr0_gs,v_nc0_gs, &
             N_gs,k_gs,hs,coff,srsp_gs,info,x0,li)
	r_gs=[phi0_gs,phi_n0_gs,H0_gs]
	a_gs=a0*exp(N_gs)
	!!print **,"=====###############======at N H",H0,N
	rs_gs=[phi0_gs,phi_n0_gs,H0_gs,qr0_gs,qc0_gs,q_nr0_gs,q_nc0_gs,vr0_gs,vc0_gs,v_nr0_gs,v_nc0_gs]
	!zz_gs=zzzz(coff,phi_n0_gs,a_gs,H0_gs,phi0_gs)
	
      !  print "efold number",n_gs	
	qr0_gs=rs_gs(4)
	qc0_gs=rs_gs(5)
	vr0_gs=rs_gs(8)
	vc0_gs=rs_gs(9)
	!ee=((phi_n0_gs)**2)/2.d0
	!eta1=eta(coff,phi0_gs,phi_n0_gs,H0_gs)
	!srsp=H0**2/(ee*8*PI**2)
	sP_gs=((qr0_gs**2+qc0_gs**2))*k_gs**3.0/2.0/pi**2
	tP_gs=((vr0_gs**2+vc0_gs**2))*k_gs**3.0/2.0/pi**2
	!tP_gs=((vr0_gs**2/a_gs**2+vc0_gs**2/a_gs**2))*k_gs**3.0/2.0/pi**2
	!coff(7)=phi0_gs
	
	
	!do 11 l = 1,j
	
  	!write(jinn,*) (psdata(l,i),i=1,18)
  
  !11  continue
	!close (jinn)
	
	
!	Print *, "k=" ,k_gs,"srsp=",srsp_gs,"sp=",sp_gs
	!!print **,"nnnnnn",N
	end subroutine


subroutine get_f(kk,power,power1,zaic,deltac,jin,j,nreh,wreh,coff,info,x0,ljin)


	use mod1
	use spectral	
	implicit none


!##########################################################################################################
!k_ini=the lowest value of k,k_end=the highest value of k, sp=scaler power spectra 
!,tp=tensor power spectra,coff(15)=arbitary initial value ,nef=no of efolds requirement for inflation
!ns=scalar spectral index,nt=tensor spectral index,phi=scalar field, phi_n=1st dervative of scalar wrt efolds
!H-hubble parameter , N=efolds,
!##########################################################################################################
!module 1 specifies the potential function changing it we can shift to diffrent models
!in first part taking any arbitary value of coff(15) will lead to find a particular value for required efolds
!then we evaluate the scalar and tensor power specta along with spectral indexes and save the data to files
!###########################################################################################################




        real(dl) R,R_end,vdc,svdc,rat,Ms,peak,vc,beta,zai,zai_end,intgn1,nn,vdc1,vdc2,intgn2,intgn
        real(dl) beta2,beta1,srsp,W,dk,k,svdc1,peak1,vc1,rat1,zaic,mass,cont,mass1,dm,eta1,x0
	integer i,j,deltac,jin,is,ljin
	real(dl) :: power(500),nreh,wreh,betaN,peakN,betaN1,peakN1,cont_peak,cont_ps,cont_peakn,cont_psn,betaform
	real(dl) :: power1(500),massn,massn1,dmn,a_form,a_eq,a_end,a0_rad,a0_infl,intgnn,intgnn1,r_in,max_peakN,max_m,max_betaN
	real(dl):: kk(500),mp,gev,Tg,gdof,Noef,info(25),max_betaform
	real(dl):: coff(25),f_PS,f_PSN,f_Peak,f_PeakN,f_PS1,f_PSN1,f_Peak1,f_PeakN1,hs,cont_peakn1,cont_psn1,max_m1,sp00
	character(90) rff,filename
	
	
			cont_peak=0.d0
			cont_ps=0.d0
			cont_peakn=0.d0
			cont_psn=0.d0
			cont_peakn1=0.d0
			cont_psn1=0.d0
		!	max_betaN=0.d0
		!	max_peakn=0.d0
	hs=1E-4
	a_eq=1.d0/3300.d0
	gdof=106.75
	mp=1.22*1e19
	gev=2.488767*1e37*2.d0*Pi
	tg=2.7255*8.61733*1e-14/mp
!	a0_infl=coff(8)
!	noef=coff(9)
!	a_end=a0_infl*Exp(noef)
!	a0_rad=a_end*Exp(nreh)
	
	
	! do while (zaic.le..414)
		 
		!!print **,"=====VVVVV1===#################==========================="
		WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/mass",jin,"ln",ljin,".dat"
		!WRITE(filename,'(a,i4.4,a)') "mass",jin,".TXT"
	 	OPEN(unit=deltac,file=filename ,action='write',position='append')
        WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/betamax",jin,"ln",ljin,".dat"
		!WRITE(filename,'(a,i4.4,a)') "mass",jin,".dat"
	 	OPEN(unit=15,file=filename ,action='write',position='append')
		WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/cont11",jin,"ln",ljin,".dat"
	 	!WRITE(filename,'(a,i4.4,a)') "cont11",jin,".TXT"
		open (unit = 19, file = filename ,action='write',position='append')
		WRITE(filename,'(a,i4.4,a,i4.4,a)') "../data/mass25_",jin,"ln",ljin,".dat"
		!WRITE(filename,'(a,i4.4,a)') "mass25_",jin,".TXT"
		open (unit = 18, file = filename ,action='write',position='append')
		!WRITE(filename,'(a,i4.4,a)') "mass_frac",jin,".dat"
		WRITE(filename,'(a,i4.4,a)') "../data/mass_frac",jin,".TXT"
		open(unit=1234, file=filename)
		close(1234, status='delete')
		open (unit = 123, file = filename ,action='write',position='append')
		i=1
	
	
	


	
	
	
	
	
		R=1e-22
		R_end=1e-4
		cont=0.d0
		max_peakN=0.d0
		max_m=0.d0
		max_betaN=-0.d0
                max_betaform=-0.d0
		max_m1=0.d0
                r_in=0
		do while (R.le.R_end)
		r_in=r_in+1
			!a_form=a0_rad*Sqrt(R)
			a_form=(R*.07*.14)*a_eq
			vdc=0.d0
			vdc1=0.d0
			vdc2=0.d0
			svdc=0.d0
			svdc1=0.d0
			nn=0.0
			i=1
			!!print **,"checkpoint1"
			do while (i.le.j-2)
				k=kk(i)
				dk=kk(i+1)-k
				!k=kk(i)
				sp00=power(i)
				srsp=power1(i)
				W=wind(k,R)
				vdc=vdc+16.d0/81.d0*(k**3)*(R**4)*sp00*(W**2)*dk
				vdc1=vdc1+16.d0/81.d0*(k**3)*(R**4)*srsp*(W**2)*dk
				svdc=svdc+16.d0/81.d0*(k**5)*(R**4)*sp00*(W**2)*dk
				svdc1=svdc1+16.d0/81.d0*(k**5)*(R**4)*srsp*(W**2)*dk
				i=i+1
				!!print **,kk(i),power(i),power1(i),vdc,vdc1,svdc,svdc1
				!!print **,"=====VVVVV1===#################==========================="
				
			end do
			!write (16,*) R,vdc,vdc1,svdc,svdc1
			!!print **, "vdc=====",vdc
			rat=svdc/vdc
			rat1=svdc1/vdc1
			zai=zaic
			zai_end=100000
			intgn1=0.d0
			intgn=0.d0
			do while (zai.le.zai_end)
					intgn=intgn+2.0*exp(-.5*zai**2/vdc+56.8)*R_end/R/sqrt(2*Pi)/sqrt(vdc)*zai*.1
					intgnN=intgnN+2.0*exp(-.5*zai**2/vdc)/sqrt(2*Pi)/sqrt(vdc)*zai*.1
					intgn1=intgn1+2.0*exp(-.5*zai**2/vdc1+56.8)*R_end/R/sqrt(2*Pi)/sqrt(vdc1)*zai*.1
					intgnN1=intgnN1+2.0*exp(-.5*zai**2/vdc1)/sqrt(2*Pi)/sqrt(vdc1)*zai*.1
					zai=zai*1.1
					!write (*,*)"intgn1,2,3,nn",intgn,intgn1
				
			end do	

			
			
			vc=zaic/sqrt(vdc)
			vc1=zaic/sqrt(vdc)
			beta=intgn
			beta1=intgn1
			betaform=Erfc(vc/sqrt(2.0))
			betaN=Erfc(vc/sqrt(2.0))/a_form*a_eq

			betaN1=Erfc(vc1/sqrt(2.0))/a_form*a_eq
			Mass=.2*3.d0/2.d0*((3.d0/100.d0)**(1/3))*((.07*.14*R)**2)*(7.0e17/1.98)
			massn=.2*4.0*Pi*(R*a_form)*gev*mp**2*1.602*1e-10/(3.d0*1e8)**2/2.d0/1e30
            if(max_betaform.le.betaform)then
				max_betaform=betaform
				!max_m1=massn
			endif
			if(max_betaN.le.betaN)then
				max_betaN=betaN
				max_m1=massn
			endif	
			Ms=(2.d0*Pi)**(3.0/2.0)*(R)**3
			peak=(((rat/3.d0)**(3.0/2.0))*(exp(-.5*vc**2))*(vc**2 - 1.d0)/(2.d0*Pi)**2)*Ms*exp(56.8d0)*R_end/R
			peakn=(((rat/3.d0)**(3.0/2.0))*(exp(-.5*vc**2))*(vc**2 - 1.d0)/(2.d0*Pi)**2)*ms/a_form*a_eq
			peak1=(((rat1/3.d0)**(3.0/2.0))*(exp(-.5*vc1**2))*(vc1**2 - 1.d0)/(2.d0*Pi)**2)*Ms*exp(56.8d0)*R_end/R
			peakn1=(((rat1/3.d0)**(3.0/2.0))*(exp(-.5*vc1**2))*(vc1**2 - 1.d0)/(2.d0*Pi)**2)*ms/a_form*a_eq
			if(max_peakN.le.peakN)then
				max_peakN=peakN
				max_m=massn
			endif	
			!call contraints(Mass)
			!beta=exp(-.5*vc**2+0.5d0*150**2)/sqrt(2.d0*Pi)/vc
			!write (*,*) R,Ms,beta,rat,vc,peak
			!write (14,*) R,Mass,beta/.42,peak/.42,beta1/.42,peak1/.42
			
		

			write(deltac,*) massn,betan,betan1,peakn,peakn1
			write(123,*) mass,Massn,betan,betan1,betaform,peakn		
			write(15,*) R,massn, betaform,max_betaform,betan,vc	    !write (*,*) R,beta,vdc
			!write (*,*) "===============================///////======================================================"
			!write (*,*) R, beta,peak
			!write (*,*) "===============================///////======================================================"
			!close(14)	
			!close(16)
			
			R=R*1.01
			Mass1=0.2*3.d0/2.d0*((3.d0/100.d0)**(1/3))*((.07*.14*R)**2)*(7.0e17/1.98)
			massn1=.2*4.0*Pi*(R*a_form)*gev*mp**2*1.602*1e-10/(3.d0*1e8)**2/2.d0/1e30
			dm=mass1-mass
			dmn=massn1-massn
			cont_peak=cont_peak+dm*peak/mass
			cont_ps=cont_ps+dm*beta/mass
			cont_peakn=cont_peakn+dmn*peakn/massn
			cont_psn=cont_psn+dmn*betan/massn
			cont_peakn1=cont_peakn1+dmn*peakn1/massn
			cont_psn1=cont_psn1+dmn*betan1/massn
			write (18,*) zaic,cont_peakn,cont_psn,cont_peakn1,cont_psn1			
			

		enddo
		f_Peakn1=cont_peakn1/.42
		f_psn1=cont_psn1/.42
		f_peakn=cont_peakn/.42
		f_psn=cont_psn/.42
		coff(11)=f_PSn1
		coff(12)=f_PSN
		!coff(12)=max_betaform
		coff(13)=f_Peakn1
		coff(14)=f_PeakN
		info(10)=max_m1
		info(11)=max_betaN
		info(12)=max_m
		info(13)=max_peakN
		
		!Print *, f_ps,f_psn,f_peak,f_peakn

		!zaic=zaic+.01
		close(123)
		close(18)
		close(15)
		close(deltac)
		write (*,*) "jao pakhi balo hawa chalochalo abchaya janlar kanch ami ki amake hariyechi bake rupkathar anach kanach:/"
		write (*,*) f_PeakN,f_Peakn1,f_PSN,f_PSn1
		write (19,*)  f_PeakN,f_Peakn1,f_PSN,f_PSn1
		close(19)
!	enddo
	
	
	end subroutine

        subroutine print_Matrix(array, n, m)
	use mod1
	implicit none
	real, intent(in) :: array(n,m)
	integer, intent(in) :: n,m
	integer :: i,j
	do i = 1,m
	!do j = 1,n
	!print **, (array(j,i),j = 1,n)
	!end do
	!print **, "\n"
	end do
	end subroutine

!end module scalar_pow
