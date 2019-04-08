program testfem

use inout2 
use matrix_solve
use hybrd_module
use minpackfunctions
	
implicit none

    integer, parameter :: pr = selected_real_kind(15) 
    character(30) :: path !variable for path for fetching data
    integer dims(2) !variable for fetching dimensions 
    real(pr), dimension(:,:), allocatable :: Points,  Triangles, Edges

    real(pr), dimension(:,:), allocatable :: Ku,Kv,Kmatrix,Kmatrix_copy,Ku_copy,Lower_Left_K, K_initial
    real(pr), dimension(:), allocatable :: Fu,Fv,F,c,F_copy,Fv_copy,f_initial_copy,c0,f_initial
	real(pr) :: x1,x2,x3,y1,y2,y3,Ak,fixed_sum,dphi1_dr,dphi2_dr,dphi3_dr,dphi1_dz,dphi2_dz,dphi3_dz
	real(pr) :: sigma_ur,sigma_uz,sigma_vr,sigma_vz,V_mu_ref,E_a_vmu_ref,V_mfv_ref
	real(pr) :: E_a_vmfv_ref,K_mu,K_mv,K_mfu,r_q,rho_u,rho_v,p_atm
	real(pr)::R_g,V_mu,V_mfv,Cu_amb,Cv_amb
	real(pr) :: Temper,T_stan,eta_u,eta_v,T_ref
	real(pr) :: k_e1_1, k_e2_2, k_e3_3, k_e1_2, k_e1_3, k_e2_3
	real(pr) :: ke11,ke12,ke22,fe1,fe2, length_edge,finish,start

	integer triangle_index,triangle1,triangle2,triangle3, edge_index, edge1, edge2
    	integer T,E,M

	!hybrd related variables
   	integer maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
   	real(pr) xtol,epsfcn,factor
    	real(pr), dimension(:), allocatable:: fvec,diag, r,qtf, wa1, wa2, wa3, wa4
    	real(pr), dimension(:,:), allocatable:: fjac
	
	external fcn !needed for nonlinear function calculation by hybrd

	!integer(pr), dimension(:,:), allocatable :: Triangles, Edges
    !real(pr), dimension(:), allocatable :: Fu,Cu
	
    call cpu_time(start)

    !Read the required mesh information

    !POINTS/NODES INLADEN
    path='Pear2/Points'
    call matrix_dims(dims,path)	
    !print *, 'Dimensies van nodesmatrix: ', dims(1), ',', dims(2)
    allocate(Points(dims(1),dims(2)))
    call read_matrix(Points,path)
    !print  *, 'De 4e waarde van de tweede rij', Points(2,4)

    M = dims(2)

    !TRIANGLES INLADEN 
    path='Pear2/Triangles'
    call matrix_dims(dims,path)	
    !print *, 'Dimensies van trianglesmatrix: ', dims(1), ',', dims(2)
    allocate(Triangles(dims(1),dims(2)))
    call read_matrix(Triangles,path)
    !print  *, 'De 4e waarde van de derde rij', Triangles(3,4)

    T = dims(2)	

    ! EDGES INLADEN  
    path='Pear2/Edges' 
    call matrix_dims(dims,path)	
    !print *, 'Dimensies van edgesmatrix: ', dims(1), ',', dims(2)
    allocate(Edges(dims(1),dims(2)))
    call read_matrix(Edges,path)
    !print  *, 'De 4e waarde van de tweede rij', Edges(2,4)

    E = dims(2)

!------------- ALLES OPGEHAALD
!-------CONSTANTEN DEFINIEREN 
 	T_ref=293.15
    T_stan = 273.15
 	Temper=T_stan+(-1.0)
    eta_u = 0.02
    eta_v =0.007
	
	sigma_ur= 2.8*10.0**(-10.0)
    sigma_uz= 1.10*10.0**(-9.0)
    sigma_vr= 2.32*10.0**(-9.0)
    sigma_vz= 6.97*10.0**(-9.0)
  

   	V_mu_ref=2.39*10.0**(-4.0)
    E_a_vmu_ref = 80200.0
    V_mfv_ref= 1.61*10.0**(-4.0)
    E_a_vmfv_ref = 56700


 	K_mu=0.4103
    K_mv=27.2438
    K_mfu=0.1149

 	r_q=0.97;

    rho_u=7.0*10.0**(-7.0)
    rho_v=7.5*10.0**(-7.0)

    p_atm=101300.0
    R_g= 8.314

    V_mu=V_mu_ref*exp((E_a_vmu_ref/R_g)*((1.0/T_ref)-(1.0/Temper)))
    V_mfv = V_mfv_ref*exp((E_a_vmfv_ref/R_g)*((1.0/T_ref)-(1.0/Temper)))

    Cu_amb = p_atm*eta_u/(R_g*Temper)
    Cv_amb = p_atm*eta_v/(R_g*Temper)


!--------------------------

   allocate(Ku(M,M))
   allocate(Kv(M,M))
   allocate(Ku_copy(M,M))
   allocate(Fu(M))
   allocate(Fv(M))
   allocate(Kmatrix(2*M,2*M))
   allocate(Kmatrix_copy(2*M,2*M))
   allocate(F(2*M))
   allocate(F_copy(2*M))
   allocate(c(2*M))
   allocate(Fv_copy(M))
   allocate(Lower_Left_K(M,M))
   allocate(K_initial(2*M,2*M))	
   allocate(f_initial_copy(2*M))
   allocate(c0(2*M))
   allocate(f_initial(2*M))

   !Forloop over triangles
	do triangle_index=1,T
		triangle1=int(Triangles(1,triangle_index)) !slechte benaming dit is trianglenode1
		triangle2=int(Triangles(2,triangle_index))
		triangle3=int(Triangles(3,triangle_index))
		
		x1= Points(1,triangle1)		
		y1= Points(2,triangle1)		
		x2= Points(1,triangle2)		
		y2= Points(2,triangle2)		
		x3= Points(1,triangle3)		
		y3= Points(2,triangle3)

        	Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2.0 !Oppervlakte 

		!Aanpassingen beginnen

		fixed_sum = (1.0/6.0)*(x1+x2+x3)
		dphi1_dr = -(1.0/(2.0*Ak))*(y3-y2)
    		dphi2_dr = -(1.0/(2.0*Ak))*(y1-y3)
    		dphi3_dr = -(1.0/(2.0*Ak))*(y2-y1)

    		dphi1_dz = -(1.0/(2.0*Ak))*(x2-x3)
    		dphi2_dz = -(1.0/(2.0*Ak))*(x3-x1)
    		dphi3_dz = -(1.0/(2.0*Ak))*(x1-x2)

		!Ku adjustments:
		k_e1_1 = ((sigma_ur*dphi1_dr*dphi1_dr)+(sigma_uz*dphi1_dz*dphi1_dz))*fixed_sum*2.0*Ak 
		k_e2_2 = ((sigma_ur*dphi2_dr*dphi2_dr)+(sigma_uz*dphi2_dz*dphi2_dz))*fixed_sum*2.0*Ak
		k_e3_3 = ((sigma_ur*dphi3_dr*dphi3_dr)+(sigma_uz*dphi3_dz*dphi3_dz))*fixed_sum*2.0*Ak


		k_e1_2 = ((sigma_ur*dphi1_dr*dphi2_dr)+(sigma_uz*dphi1_dz*dphi2_dz))*fixed_sum*2.0*Ak

		k_e1_3 = ((sigma_ur*dphi1_dr*dphi3_dr)+(sigma_uz*dphi1_dz*dphi3_dz))*fixed_sum*2.0*Ak

		k_e2_3 = ((sigma_ur*dphi2_dr*dphi3_dr)+(sigma_uz*dphi2_dz*dphi3_dz))*fixed_sum*2.0*Ak


		Ku(triangle1,triangle1) = Ku(triangle1,triangle1) + k_e1_1
		Ku(triangle1,triangle2) = Ku(triangle1,triangle2) + k_e1_2
		Ku(triangle1,triangle3) = Ku(triangle1,triangle3) + k_e1_3
		
		Ku(triangle2,triangle1) = Ku(triangle2,triangle1) + k_e1_2 !symm
		Ku(triangle2,triangle2) = Ku(triangle2,triangle2) + k_e2_2
		Ku(triangle2,triangle3) = Ku(triangle2,triangle3) + k_e2_3
		
		Ku(triangle3,triangle1) = Ku(triangle3,triangle1) + k_e1_3 !symm
		Ku(triangle3,triangle2) = Ku(triangle3,triangle2) + k_e2_3 !symm
		Ku(triangle3,triangle3) = Ku(triangle3,triangle3) + k_e3_3

		!Kv adjustments:
		k_e1_1 = ((sigma_vr*dphi1_dr*dphi1_dr)+(sigma_vz*dphi1_dz*dphi1_dz))*fixed_sum*2.0*Ak
		k_e2_2 = ((sigma_vr*dphi2_dr*dphi2_dr)+(sigma_vz*dphi2_dz*dphi2_dz))*fixed_sum*2.0*Ak
		k_e3_3 = ((sigma_vr*dphi3_dr*dphi3_dr)+(sigma_vz*dphi3_dz*dphi3_dz))*fixed_sum*2.0*Ak

		k_e1_2 = ((sigma_vr*dphi1_dr*dphi2_dr)+(sigma_vz*dphi1_dz*dphi2_dz))*fixed_sum*2.0*Ak
		k_e1_3 = ((sigma_vr*dphi1_dr*dphi3_dr)+(sigma_vz*dphi1_dz*dphi3_dz))*fixed_sum*2.0*Ak
		k_e2_3 = ((sigma_vr*dphi2_dr*dphi3_dr)+(sigma_vz*dphi2_dz*dphi3_dz))*fixed_sum*2.0*Ak

		Kv(triangle1,triangle1) = Kv(triangle1,triangle1) + k_e1_1
		Kv(triangle1,triangle2) = Kv(triangle1,triangle2) + k_e1_2
		Kv(triangle1,triangle3) = Kv(triangle1,triangle3) + k_e1_3
		
		Kv(triangle2,triangle1) = Kv(triangle2,triangle1) + k_e1_2 !symm
		Kv(triangle2,triangle2) = Kv(triangle2,triangle2) + k_e2_2
		Kv(triangle2,triangle3) = Kv(triangle2,triangle3) + k_e2_3
		
		Kv(triangle3,triangle1) = Kv(triangle3,triangle1) + k_e1_3 !symm
		Kv(triangle3,triangle2) = Kv(triangle3,triangle2) + k_e2_3 !symm
		Kv(triangle3,triangle3) = Kv(triangle3,triangle3) + k_e3_3

	enddo !All triangles have been handled

	!Forloop over edges
	do edge_index=1,E
		edge1 = int(Edges(1,edge_index))
		edge2 = int(Edges(2,edge_index))

		x1= Points(1,edge1)		
		y1= Points(2,edge1)		
		x2= Points(1,edge2)		
		y2= Points(2,edge2)


		if ((x1 .GT. 0.0017) .OR. (x2 .GT. 0.0017)) then !check if on gamma 2 edge then
			
		    	length_edge = sqrt((x2-x1)**2.0 + (y2-y1)**2.0);
				
				!K_adjustments:
				!Ku
   				ke11 = (1.0/12.0)*(3.0*x1 + x2)*rho_u*length_edge
    			ke12 = (1.0/12.0)*(x1 + x2)*rho_u*length_edge
   				ke22 = (1.0/12.0)*(x1 + 3.0*x2)*rho_u*length_edge
			
				Ku(edge1,edge1) = Ku(edge1,edge1)  + ke11
				Ku(edge1,edge2) = Ku(edge1,edge2)  + ke12
				Ku(edge2,edge1) = Ku(edge2,edge1)  + ke12
				Ku(edge2,edge2) = Ku(edge2,edge2)  + ke22
				!Kv
				ke11 = (1.0/12.0)*(3.0*x1 + x2)*rho_v*length_edge
    			ke12 = (1.0/12.0)*(x1 + x2)*rho_v*length_edge
   				ke22 = (1.0/12.0)*(x1 + 3.0*x2)*rho_v*length_edge
				
				Kv(edge1,edge1) = Kv(edge1,edge1)  + ke11
				Kv(edge1,edge2) = Kv(edge1,edge2)  + ke12
				Kv(edge2,edge1) = Kv(edge2,edge1)  + ke12
				Kv(edge2,edge2) = Kv(edge2,edge2)  + ke22

				!F_adjustments
				!Fu
				fe1 = -length_edge*(Cu_amb*rho_u)*(1.0/6.0)*(2.0*x1 + x2) 
       			fe2 = -length_edge*(Cu_amb*rho_u)*(1.0/6.0)*(x1 + 2.0*x2) 

				Fu(edge1)= Fu(edge1) - fe1
				Fu(edge2)= Fu(edge2) - fe2
				
				!Fv
				fe1 = -length_edge*(Cv_amb*rho_v)*(1.0/6.0)*(2.0*x1 + x2)
       			fe2 = -length_edge*(Cv_amb*rho_v)*(1.0/6.0)*(x1 + 2.0*x2)

				Fv(edge1)= Fv(edge1) - fe1
				Fv(edge2)= Fv(edge2) - fe2

		end if
	enddo !loop over edges is done


	Kmatrix = 0
	Kmatrix(1:M,1:M) = Ku
	Kmatrix(M+1:2*M,M+1:2*M) = Kv
	Kmatrix_copy=Kmatrix
	F(1:M) = Fu
	F(M+1:2*M) = Fv

	!F IS OK; gecheckt ZELFDE ALS MATLAB
	!K is OK; gecheckt zelfde als matlab

	F_copy = F
	call gauss_pivot(Kmatrix,F,2*M) !in F komt dan c te zitten
	
	Kmatrix = Kmatrix_copy !restore Kmatrix
	c = F
	!print *, c
	!c pure diffusie is OK; gecheckt zelfde als matlab
	F = F_copy
	
	Ku_copy = Ku
	Fv_copy = Fv

	do triangle_index=1,T
	
		triangle1=int(Triangles(1,triangle_index))
		triangle2=int(Triangles(2,triangle_index))
		triangle3=int(Triangles(3,triangle_index))
		
		x1= Points(1,triangle1)		
		y1= Points(2,triangle1)		
		x2= Points(1,triangle2)		
		y2= Points(2,triangle2)		
		x3= Points(1,triangle3)		
		y3= Points(2,triangle3)

		Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2.0 !Oppervlakte 

		k_e1_1 =  (1.0/60.0)*(3.0*x1 + x2 + x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak
	    
	    	k_e1_2 = (1.0/60.0)*(x1 + x2 + (1.0/2.0)*x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak
	   
	    	k_e1_3 = (1.0/60.0)*(x1 + (1.0/2.0)*x2 + x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak
	   
	    	k_e2_2 =  (1.0/60.0)*(x1 + 3.0*x2 + x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak
	    
	    	k_e2_3 = (1.0/60.0)*((1.0/2.0)*x1 + x2 + x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak
	    
	    	k_e3_3 =  (1.0/60.0)*(x1 + x2 + 3.0*x3) * ((V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak

		Ku_copy(triangle1,triangle1) = Ku_copy(triangle1,triangle1) + k_e1_1
		Ku_copy(triangle1,triangle2) = Ku_copy(triangle1,triangle2) + k_e1_2
		Ku_copy(triangle1,triangle3) = Ku_copy(triangle1,triangle3) + k_e1_3

		Ku_copy(triangle2,triangle1) = Ku_copy(triangle2,triangle1) + k_e1_2
		Ku_copy(triangle2,triangle2) = Ku_copy(triangle2,triangle2) + k_e2_2
		Ku_copy(triangle2,triangle3) = Ku_copy(triangle2,triangle3) + k_e2_3

		Ku_copy(triangle3,triangle1) = Ku_copy(triangle3,triangle1) + k_e1_3
		Ku_copy(triangle3,triangle2) = Ku_copy(triangle3,triangle2) + k_e2_3
		Ku_copy(triangle3,triangle3) = Ku_copy(triangle3,triangle3) + k_e3_3


		Fv_copy(triangle1) = Fv_copy(triangle1) - V_mfv*(- x1/12.0 - x2/24.0- x3/24.0)*2.0*Ak
		Fv_copy(triangle2) = Fv_copy(triangle2) - V_mfv*(- x1/24.0 - x2/12.0- x3/24.0)*2.0*Ak
		Fv_copy(triangle3) = Fv_copy(triangle3) - V_mfv*(- x1/24.0 - x2/24.0- x3/12.0)*2.0*Ak

		k_e1_1 =  -(1.0/60.0)*(3.0*x1 + x2 + x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;
    
    		k_e1_2 = -(1.0/60.0)*(x1 + x2 + (1.0/2.0)*x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;
    		k_e1_3 = -(1.0/60.0)*(x1 + (1.0/2.0)*x2 + x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;
    
    		k_e2_2 = - (1.0/60.0)*(x1 + 3.0*x2 + x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;
    
    		k_e2_3 = -(1.0/60.0)*((1.0/2.0)*x1 + x2 + x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;
    		k_e3_3 =  -(1.0/60.0)*(x1 + x2 + 3.0*x3) * ((r_q*V_mu)/(((K_mu)*(1.0+Cv_amb/K_mv))))*2.0*Ak;

		Lower_Left_K(triangle1,triangle1) = Lower_Left_K(triangle1,triangle1) + k_e1_1
		Lower_Left_K(triangle1,triangle2) = Lower_Left_K(triangle1,triangle2) + k_e1_2
		Lower_Left_K(triangle1,triangle3) = Lower_Left_K(triangle1,triangle3) + k_e1_3

		Lower_Left_K(triangle2,triangle1) = Lower_Left_K(triangle2,triangle1) + k_e1_2
		Lower_Left_K(triangle2,triangle2) = Lower_Left_K(triangle2,triangle2) + k_e2_2
		Lower_Left_K(triangle2,triangle3) = Lower_Left_K(triangle2,triangle3) + k_e2_3

		Lower_Left_K(triangle3,triangle1) = Lower_Left_K(triangle3,triangle1) + k_e1_3
		Lower_Left_K(triangle3,triangle2) = Lower_Left_K(triangle3,triangle2) + k_e2_3
		Lower_Left_K(triangle3,triangle3) = Lower_Left_K(triangle3,triangle3) + k_e3_3
		
	enddo

   K_initial = 0
   K_initial(1:M,1:M) = Ku_copy
   K_initial(M+1:2*M,1:M) = Lower_Left_K
   K_initial(M+1:2*M,M+1:2*M) = Kv
   f_initial(1:M) = Fu
   f_initial(M+1:2*M) = Fv_copy

   f_initial_copy = f_initial
   call gauss_pivot(K_initial,f_initial,2*M) !in f_initial komt dan c0 te zitten
   c0 = f_initial

   !print *,c0	
   !c0_initial is OK; gecheckt zelfde als matlab 

   !nu moeten we het niet-lineair stelsel gaan bouwen en oplossen. Daartoe gebruiken we de library minpack van fortran en specifiek, de functie om een systeem van niet-lineaire vergelijkingen op te lossen aan de hand van een benaderde jacobiaan: hybrd.

	c=c0 ! initialise c as c0

	!on input c contains the initial estimate
	!on output c contains the final solution after running through this



	ml=2*M !assume jacobian not to be banded
	mu=2*M !assume jacobian not to be banded
	xtol=10.0**(-13)
        maxfev=100
        epsfcn=10.0**(-12)
        mode=1 !variables are scaled internally
        factor=.2 !recommended value
        nprint=1 !print toeter
	ldfjac=2*M
	lr=2*M*(2*M+1)/2

     	allocate(fvec(2*M))
        allocate(diag(2*M))
        allocate(fjac(2*M,2*M))
        allocate(qtf(2*M))
        allocate(wa1(2*M))
        allocate(wa2(2*M))
        allocate(wa3(2*M))
        allocate(wa4(2*M))
        allocate(r(lr))

	!c0(1:M)=10000.0
	!c0(M+1:2*M) = 20000.0
	!print *,'zever c'
	!print *, c0
	!Test om te kijken wat oproep van fcn zoal bekomt voor bepaalde c
        !call fcn(2*M,c0,fvec,0,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q) 

	!print *,'f evaluation with initial c'
	!print *, fvec
        
	call hybrd(fcn,2*M,c,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,factor, &
& nprint,info, nfev,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)

	!PROBLEEM: fcn evaluatie levert sowieso kleine f op. Daardoor doet hybrd ook amper aanpassingen. Is onze evaluatie F(x) verkeerd?

	!c bevat de oplossing na uitwerken van het niet-lineaire systeem
	!print *,c	
	!print *,info

	deallocate(fvec)
        deallocate(diag)
        deallocate(fjac)
        deallocate(qtf)
        deallocate(wa1)
        deallocate(wa2)
        deallocate(wa3)
        deallocate(wa4)
        deallocate(r)

   deallocate(Kmatrix)
   deallocate(Kmatrix_copy)
   deallocate(c)
   deallocate(F_copy)
   deallocate(F)
   deallocate(Ku_copy)
   deallocate(Fv_copy)
   deallocate(Lower_Left_K)
   deallocate(K_initial)
   deallocate(f_initial)
   deallocate(c0)
   deallocate(f_initial_copy)

   deallocate(Ku)
   deallocate(Kv)
   deallocate(Fu)
   deallocate(Fv)

	call cpu_time(finish)
	!print*,'time ', finish-start
		
end program testfem

subroutine fcn(n,x,fvec,iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q) 
	implicit none
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd.
!         in this case set iflag to a negative integer.
!	  n -> Amount of variables (nodes*2)
!	  
!	THIS FUNCTION EVALUATES F(CU,CV)

        integer, parameter :: fcn_pr = selected_real_kind(15) 
	integer n,iflag, M, triangle_index, T
	integer triangle1, triangle2, triangle3
	real(fcn_pr) :: x(n), fvec(n), Kmatrix(n,n), F(n), Points(2,n/2),Triangles(3,T)
	real(fcn_pr) :: c_u(n/2),c_v(n/2)
	real(fcn_pr) :: x1,x2,x3,y1,y2,y3,c_u_coeff,c_v_coeff,V_mu,K_mv,K_mu,fixed_sum_one_third,r_u
	real(fcn_pr) :: fixed_sum_one_fifth, V_mfv,r_q,K_mfu, r_v_const, Ak
	

	!if (iflag .e. 0) go to 100

      	M = n/2
	fvec = matmul(Kmatrix,x)-F !makkelijke stuk
	c_u =x(1:M)
	c_v =x(M+1:n)
	
	!Nu Hu en Hv benaderen

	do triangle_index=1,T

		triangle1=int(Triangles(1,triangle_index))
		triangle2=int(Triangles(2,triangle_index))
		triangle3=int(Triangles(3,triangle_index))
		
		x1= Points(1,triangle1)		
		y1= Points(2,triangle1)		
		x2= Points(1,triangle2)		
		y2= Points(2,triangle2)		
		x3= Points(1,triangle3)		
		y3= Points(2,triangle3)
		
		Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2.0 !Oppervlakte 

		!------------------Deel vanuit H_u

		!part1
		c_u_coeff = c_u(triangle1)*(1.0/3.0) + c_u(triangle2)*(1.0/3.0) + c_u(triangle3)*(1.0/3.0) 
		c_v_coeff = c_v(triangle1)*(1.0/3.0) + c_v(triangle2)*(1.0/3.0) + c_v(triangle3)*(1.0/3.0) 

		r_u = (V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		fixed_sum_one_third = (1.0/3.0)*(x1 + x2 + x3) 

		fvec(triangle1) =fvec(triangle1) + (-27.0/96.0)*fixed_sum_one_third*r_u*(1.0/3.0)*2.0*Ak

		fvec(triangle2) = fvec(triangle2) + (-27.0/96.0)*fixed_sum_one_third*r_u*(1.0/3.0)*2.0*Ak

		fvec(triangle3) = fvec(triangle3) + (-27.0/96.0)*fixed_sum_one_third*r_u*(1.0/3.0)*2.0*Ak

		!part2

		c_u_coeff = c_u(triangle1)*(3.0/5.0) + c_u(triangle2)*(1.0/5.0) + c_u(triangle3)*(1.0/5.0) 
		c_v_coeff = c_v(triangle1)*(3.0/5.0) + c_v(triangle2)*(1.0/5.0) + c_v(triangle3)*(1.0/5.0) 

		r_u = (V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		fixed_sum_one_fifth = (3.0/5.0)*x1 + (1.0/5.0)*x2 + (1.0/5.0)*x3 
		

		fvec(triangle1) = fvec(triangle1) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(3.0/5.0)*2.0*Ak

		fvec(triangle2) = fvec(triangle2) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		fvec(triangle3) = fvec(triangle3) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		!part3

		c_u_coeff = c_u(triangle1)*(1.0/5.0) + c_u(triangle2)*(1.0/5.0) + c_u(triangle3)*(3.0/5.0) 
		c_v_coeff = c_v(triangle1)*(1.0/5.0) + c_v(triangle2)*(1.0/5.0) + c_v(triangle3)*(3.0/5.0) 

		r_u = (V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		fixed_sum_one_fifth = (1.0/5.0)*x1 + (1.0/5.0)*x2 + (3.0/5.0)*x3 
		
		fvec(triangle1) = fvec(triangle1) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		fvec(triangle2) = fvec(triangle2) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		fvec(triangle3) = fvec(triangle3) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(3.0/5.0)*2.0*Ak

		!part4

		c_u_coeff = c_u(triangle1)*(1.0/5.0) + c_u(triangle2)*(3.0/5.0) + c_u(triangle3)*(1.0/5.0) 
		c_v_coeff = c_v(triangle1)*(1.0/5.0) + c_v(triangle2)*(3.0/5.0) + c_v(triangle3)*(1.0/5.0) 

		r_u = (V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		fixed_sum_one_fifth = (1.0/5.0)*x1 + (3.0/5.0)*x2 + (1.0/5.0)*x3 
		
		fvec(triangle1) = fvec(triangle1) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		fvec(triangle2) = fvec(triangle2) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(3.0/5.0)*2.0*Ak

		fvec(triangle3) = fvec(triangle3) + (25.0/96.0)*fixed_sum_one_fifth*r_u*(1.0/5.0)*2.0*Ak

		
		!---------------------------Deel vanuit H_v

		!part1
		c_u_coeff = c_u(triangle1)*(1.0/3.0) + c_u(triangle2)*(1.0/3.0) + c_u(triangle3)*(1.0/3.0) 
		c_v_coeff = c_v(triangle1)*(1.0/3.0) + c_v(triangle2)*(1.0/3.0) + c_v(triangle3)*(1.0/3.0) 

		r_u=(V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		r_v_const=V_mfv/(1.0+c_u_coeff/K_mfu)

		fixed_sum_one_third = (1.0/3.0)*(x1 + x2 + x3) 		

		fvec(triangle1+M) =fvec(triangle1+M) + (-27.0/96.0)*&
&((-fixed_sum_one_third)*r_q*r_u*(1.0/3.0)*2.0*Ak -&
&(fixed_sum_one_third*r_v_const*(1.0/3.0)*2.0*Ak))

		fvec(triangle2+M) =fvec(triangle2+M) + (-27.0/96.0)* &
& ((fixed_sum_one_third)*r_q*r_u*(1.0/3.0)*2.0*Ak-&
& (fixed_sum_one_third*r_v_const*(1.0/3.0)*2.0*Ak))
		
		fvec(triangle3+M) =fvec(triangle3+M) +(-27.0/96.0)* &
& ((-fixed_sum_one_third)*r_q*r_u*(1.0/3.0)*2.0*Ak -&
& (fixed_sum_one_third*r_v_const*(1.0/3.0)*2.0*Ak))

		!part2
		c_u_coeff = c_u(triangle1)*(3.0/5.0) + c_u(triangle2)*(1.0/5.0) + c_u(triangle3)*(1.0/5.0) 
		c_v_coeff = c_v(triangle1)*(3.0/5.0) + c_v(triangle2)*(1.0/5.0) + c_v(triangle3)*(1.0/5.0) 

		r_u=(V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		r_v_const=V_mfv/(1.0+c_u_coeff/K_mfu)

		fixed_sum_one_fifth = (3.0/5.0)*x1 + (1.0/5.0)*x2 + (1.0/5.0)*x3 

		fvec(triangle1+M) =fvec(triangle1+M) + (25.0/96.0)* &
&((-fixed_sum_one_fifth)*r_q*r_u*(3.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(3.0/5.0)*2.0*Ak))

		fvec(triangle2+M) =fvec(triangle2+M) + (25.0/96.0)* &
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))
		
		fvec(triangle3+M) =fvec(triangle3+M) + (25.0/96.0)* &
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))


		!part3

		c_u_coeff = c_u(triangle1)*(1.0/5.0) + c_u(triangle2)*(1.0/5.0) + c_u(triangle3)*(3.0/5.0) 
		c_v_coeff = c_v(triangle1)*(1.0/5.0) + c_v(triangle2)*(1.0/5.0) + c_v(triangle3)*(3.0/5.0) 

		r_u=(V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		r_v_const=V_mfv/(1.0+c_u_coeff/K_mfu)

		fixed_sum_one_fifth = (1.0/5.0)*x1 + (1.0/5.0)*x2 + (3.0/5.0)*x3 

		fvec(triangle1+M) =fvec(triangle1+M) + (25.0/96.0)*& 
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))

		fvec(triangle2+M) =fvec(triangle2+M) + (25.0/96.0)*&
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
& (fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))
		
		fvec(triangle3+M) =fvec(triangle3+M) + (25.0/96.0)*&
&((-fixed_sum_one_fifth)*r_q*r_u*(3.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(3.0/5.0)*2.0*Ak))


		!part4

		c_u_coeff = c_u(triangle1)*(1.0/5.0) + c_u(triangle2)*(3.0/5.0) + c_u(triangle3)*(1.0/5.0) 
		c_v_coeff = c_v(triangle1)*(1.0/5.0) + c_v(triangle2)*(3.0/5.0) + c_v(triangle3)*(1.0/5.0) 

		r_u=(V_mu*c_u_coeff)/((1.0+c_v_coeff/K_mv)*(K_mu+c_u_coeff))
		r_v_const=V_mfv/(1.0+c_u_coeff/K_mfu)

		fixed_sum_one_fifth = (1.0/5.0)*x1 + (3.0/5.0)*x2 + (1.0/5.0)*x3 

		fvec(triangle1+M) =fvec(triangle1+M) + (25.0/96.0)*&
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))

		fvec(triangle2+M) =fvec(triangle2+M) + (25.0/96.0)* &
&((-fixed_sum_one_fifth)*r_q*r_u*(3.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(3.0/5.0)*2.0*Ak))
		
		fvec(triangle3+M) =fvec(triangle3+M) + (25.0/96.0)* &
&((-fixed_sum_one_fifth)*r_q*r_u*(1.0/5.0)*2.0*Ak- &
&(fixed_sum_one_fifth*r_v_const*(1.0/5.0)*2.0*Ak))


	enddo

	!print *,"laten we fvec printen na deze fcn oproep"
	!print *,fvec
	
	!na afloop zit in fvec de volledige functie evaluatie voor een zekere x=(cu;cv)

end subroutine fcn








