module solve
use mod1
implicit none

contains
	






subroutine solveqn(x,d,b,coff)
        use mod1
	implicit none
	integer,parameter::n=2
	!real(dl),dimension (n)::col,y
	real(dl),dimension(n,n)::a,a1
	integer::i,j,k,l,m
	real(dl)::z,k1,k2,k3,k4,k5
	real(dl) b,d,x
	real(dl)::coff(25)
	CHARACTER:: rff
!(2 x (a + 2 b x^2 - a d x^2 + c x^4 (3 + d x^2)))/(1 + d x^2)^3
!(2 (a + 6 b x^2 - 8 a d x^2 + 15 c x^4 - 6 b d x^4 + 3 a d^2 x^4 + 4 c d x^6 + c d^2 x^8))/(1 + d x^2)^4
	!!print **, "we r in solve::type anything followed by enter to proceed"
	!!print **, "x,d,b,y(1),y(2)",x,d,b,y(1),y(2)
	!read(*,*)rff
	!x=1.000d0
	!d=2.67033d0
	!b=1.0998d0

!	y(1)=-2.0
!	y(2)=-2.0

!	DO while ( (y(1).le.0.0) .OR. (y(1).ge.1.0e34) .OR.  (y(2).le.0.0) .OR. (y(2).ge.1.0e34) )
	!!print **, "d,b,y(1),y(2)",x,d,b,y(1),y(2)
	!!print **, "we r in while loop::type anything followed by enter to proceed"
	!!print **, "x,d,b,y(1),y(2)",x,d,b,y(1),y(2)
	!read(*,*)rff
!64			if((y(1).ge.1.0e34) .or. (y(2).ge.1.0e34))then
!				Print *,"we have hit a wall in solve:to proceed say 'yes'"
!				read (*,*)rff
!				if (rff=="yes")then
!					!print **,"fine then"
!				else
!					STOP
!				endif;endif	
	
	
	
	
	
!	a(1,1)=(1.d0-d*x**2)
!	a(1,2)=2.d0*x**2
!	!a(1,2)=(x**4)*(3.d0 + d*x**2)
!	a(2,1)=(1.d0-8.d0*d*x**2+ 3.d0*d**2 *x**4)
!	a(2,2)=(6.d0*x**2 - 6.d0* d* x**4)
!	!a(2,2)=(15.d0*x**4+4.d0*d*x**6+d**2*x**8)
!	col(1)=-b*x**4 *(3.d0 + d*x**2)
!	!col(1)=-2.d0*b*x**2
!	col(2)=-(15.d0*x**4+4.d0*d*x**6+d**2*x**8)*b
!	!col(2)=-(6.d0*b*x**2- 6.d0*b*d*x**4)
!	!data a/1,-1,1,1,2,-1,-3,0,1/
!	a1(1,1)=1.d0
!	a1(1,2)=0.0d0
!	a1(2,1)=0.0d0
!	a1(2,2)=1.d0



	!data a1/1.d0,0.d0,0.d0,1.d0/
	!data b/9,6,-5/
!	call !print *m(a,a1,n)
	!read(*,*)rff
!	do i=1,n
!			if (a(i,i)==0)then
!			!print **, i
!				l=i+1
!58				do m=1,n
!						k1=a(l,m)
!						k2=a1(l,m)
!						a(l,m)=a(i,m)
!						a1(l,m)=a1(i,m)
!						a(i,m)=k1
!						a1(i,m)=k2
!				end do
!				if (a(i,i)==0)then
!				l=l-2
!				go to 58
!				end if
!			endif
!		
!		do j=i+1,n
!			if (a(i,i)==0)then
!				!print **,"in trouble"
!			endif
!			k3=a(j,i)/a(i,i)
!			do k=1,n
!			
!			!!print **,"i,j,k",i,j,k
!				a(j,k)=a(j,k)-k3*a(i,k)
!				a1(j,k)=a1(j,k)-k3*a1(i,k)
!			!call !print *m(a,a1,n)
!			
!			end do
!			call !print *m(a,a1,n)
!			!print **,"i,j,k",i,j,k
!
!		end do
!	end do
!
!	call !print *m(a,a1,n)
!	!!print **,"********************" 
!	do i=1,n
!	k4=a(i,i)
!	do k=1,n
!		a(i,k)=a(i,k)/k4
!		a1(i,k)=a1(i,k)/k4
!	enddo;enddo


!	call !print *m(a,a1,n)
!	do i=1,n-1
!	do j=i+1,n
!	k5=a(i,j)
!	do k=1,n
!		a(i,k)=a(i,k)-k5*a(j,k)
!		a1(i,k)=a1(i,k)-k5*a1(j,k)
!		
!	enddo;enddo;enddo
!	call !print *m(a,a1,n)
!	!print **,"####################"



!	CALL MULTI(a1,col,y)
!	do i=1,n
!		!print **,"x",i,"=",y(i)
!	enddo
!	10 format(2x,a,i1,a,f10.5)




!	coff(1)=1.d0                                    ! a: free parameter
!	coff(2)=(-4.d0*fg*nu)/(3.d0*phi)                 !b in terms of a
!	coff(3)=(fg*nu**2)/(2.d0*phi**2)                  !c in terms of a
!	b=b*1.1
	!d=d*1.02
!	enddo
	!b=b/1.1


end subroutine


subroutine printm(a,a1,n)
        use mod1
        implicit none
	integer::p,q,n
	real(dl),dimension(n,n)::a,a1
	do p=1,n
	!	write(*,*)(a(p,q),q=1,n) , (a1(p,q),q=1,n)
		60 format(2x,3(f10.6),10x,3(f10.6))
	enddo
	!!print **,"********************" 
end subroutine


subroutine multi(a1,col,y)
        use mod1
	implicit none
	integer,parameter::n=2
	real(dl),dimension (n)::col,y
	real(dl),dimension(n,n)::a1
	integer::i,j,k
	do i=1,n
		do j=1,n
			y(i)=0
			do k=1,n
			y(i)=y(i)+a1(i,k)*col(k)
			enddo
		enddo
	enddo
end subroutine


end module solve
