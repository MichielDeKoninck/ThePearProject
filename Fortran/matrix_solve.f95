module matrix_solve
	implicit none
	save
	!private
	!public eig
	integer, parameter :: wp = selected_real_kind(15) !define working precision which can be adapted. (6 for single, 15 for double)


contains
	subroutine gaussexplicitUandL(A,U,L) !This subroutine uses seperate memory to calculate U and L
		integer m
		real(wp), dimension (:,:), intent(in) :: A
		integer k
		integer j
		integer f
		real(wp), dimension(:,:) :: U
		real(wp), dimension(:,:) :: L

		m=size(A(:,1)) ! m is set to be the amount of rows

		!make L to be identity matrix
		do f=1,m 
			L(f,f)=1
		end do  
		
		U=A

		do  k=1,m-1 !All columns
			do  j=k+1,m !All rows within each column
 				L(j,k)=U(j,k)/U(k,k)
				U(j,k:m)=U(j,k:m)-L(j,k)*U(k,k:m)
			end do
		end do
		!The calculation of Ax=v is not executed within this subroutine, it was implemented just to compare to results of the gauss(A,v,m) subroutine. 
	end subroutine gaussexplicitUandL

	subroutine gauss(A,v) !all L and U calculations happen within memory of A
		integer m
		real(wp), dimension (:,:) :: A
		real(wp), dimension(:) :: v 
		integer k
		integer j
		integer i

		m=size(A(:,1))

		do  k=1,m-1 !All columns
				do  j=k+1,m !All rows within each column
	 				A(j,k)=A(j,k)/A(k,k)
					A(j,k+1:m)=A(j,k+1:m)-A(j,k)*A(k,k+1:m)
				end do
		end do

		!Now let's calculate x for the equation Ax=v
		!First we calculate y using Forward sub:
		do i=1,m 
			do j=1,i-1
				v(i)= v(i) - A(i,j)*v(j) !forward substitution
			end do
		end do

		!Forward substitution is completed, let's now do backward to find x
		!The result will once more be saved in v
		
		do i=m,1,-1
			do j=i+1,m
				v(i)= (v(i) - A(i,j)*v(j)) !backward substitution
			end do
			v(i)= v(i)/(A(i,i))
		end do
		
	end subroutine gauss

	subroutine gauss_pivot(A,v) !After execution v contains the solution of Ax=v where x is the solution. A will contain LU.
		integer :: m,n
		real(wp), dimension (:,:) :: A
		real(wp), allocatable :: P(:,:)
		real(wp), dimension(:) :: v 

		m=size(A(:,1))
		n=size(A(1,:))
		
		allocate(P(m,n))
		
		!implementation is split into 2 different parts
		call LUdecomp(A,P) !calculates LU of A and delivers permutation matrix
		call triangularsystemsolver(A,v,P) !solves system LUx=v

		deallocate(P)

	end subroutine gauss_pivot

	subroutine LUdecomp(A,P) !calculates the LU-decomposition of A
		integer m,n
		real(wp), dimension (:,:) :: A
		integer k,j,i
		real(wp) bigval
		integer swaprow
		real(wp), dimension (:,:) :: P
		real(wp), allocatable :: tempholder(:)
		
		m=size(A(:,1))
		n=size(A(1,:))
		
		allocate(tempholder(n))
		
		
		P=0
		do i=1,m
			P(i,i)=1
		end do
		
		do k=1,m-1
			bigval = 0 !initialise as 0 for each column 
			do i=k,m 
				if (abs(A(i,k)) > bigval) then
					bigval=abs(A(i,k))
					swaprow=i	!swaprow will contain the number of the row with the biggest element for this colomn (will always be at least k) 
				endif 
			end do
			
			!Swaps within A are executed
			tempholder = A(k,1:m)
			A(k,1:m) = A(swaprow,1:m)
			A(swaprow, 1:m) = tempholder
			!Swaps within permutation matrix are executed
			tempholder = P(k,:)
			P(k,:) = P(swaprow,:)
			P(swaprow,:) = tempholder

			!now let's do gaussian elimination with the permuted matrix:
			do  j=k+1,m 
		 		A(j,k)=A(j,k)/A(k,k)
				A(j,k+1:m)=A(j,k+1:m)-A(j,k)*A(k,k+1:m)
			end do
		end do

		deallocate(tempholder)

	end subroutine LUdecomp

	subroutine triangularsystemsolver (A,v,P) !solves LUx=v where A contains LU for x and v is already permutated
		real(wp), dimension (:,:) :: A,P
		real(wp), dimension(:) :: v
		integer m,j,i
		
		v=matmul(P,v)
		
		m=size(A(:,1))
		
		!now let's calculate x for the equation Ax=v
		!first we calculate y using Forward sub:
		do i=1,m 
			do j=1,i-1
			v(i)= v(i) - A(i,j)*v(j) !forward substitution
			end do
		end do

		!Forward substitution is completed, let's now do backward to find x
		!The result will once more be saved in v
			
		do i=m,1,-1
			do j=i+1,m
				v(i)= (v(i) - A(i,j)*v(j)) !backward substitution
			end do
			v(i)= v(i)/(A(i,i))
		end do

	end subroutine triangularsystemsolver
	
	
		
	subroutine schur(G,A,gg,hh) ! this should solve the special matrix system described in the task gg= g en hh= h
		integer m,n,i
		real(wp), dimension (:,:) :: A,G
		real(wp), dimension(:) :: gg,hh
		real(wp), allocatable :: P(:,:), A_copy_T(:,:), G_copy(:,:), temp(:,:), X(:,:)
		real(wp), allocatable :: gg_copy(:), hh_copy(:), column(:)

		!we need the original versions of the vectors and matrices so copies need to be made
		! A is mxn, G is nxn 
		m=size(A(:,1)) ! = m=size(hh)
		n=size(A(1,:)) ! = n=size(gg) 
		
		allocate(P(n,n)) !the permutation matrix is square and same size as G
		allocate(A_copy_T(n,m))
		allocate(X(n,m))
		allocate(G_copy(n,n))
		allocate(temp(m,m))
		allocate(gg_copy(n))
		allocate(column(n))

		A_copy_T=transpose(A)
		G_copy=G
		gg_copy=gg
		!----------- lambda* -----------------
		
		!1. calculate G-1AT
		
		!It would be inefficient to let gauss_pivot calculate the LU factorisation each time, therefore we make two different functions, one which calculates LU and one which uses LU to calculate x for Ax=v.
		call LUdecomp(G_copy,P) !calculate the LU decomp of G, P will be the permutation matrix
		do i=1,m !handle each column of AT
			call triangularsystemsolver(G_copy,A_copy_T(:,i),P)
		end do

		!2. calculate G-1g using the LU of G
		call triangularsystemsolver(G_copy,gg_copy,P) !!gg_copy contains G-1g

		!3. solve the eventualsystem

		temp= matmul(A,A_copy_T) ! AG-1AT

		hh= matmul(A,gg_copy)-hh !AG-1g-h 
			
		call gauss_pivot(temp,hh) 

		!the result lambda is now in hh

		!----------- p -----------------

		gg= (matmul(transpose(A),hh)-gg )
		call triangularsystemsolver(G_copy,gg,P) ! gg now contains p
		
		deallocate(P)
		deallocate(A_copy_T)
		deallocate(G_copy)
		deallocate(gg_copy)
		deallocate(temp)
		deallocate(X)
		deallocate(column)
		
	end subroutine schur


	subroutine gausspivot_kkt(G,A,gg,hh,kkt_matrix)
		real(wp), dimension (:,:), intent(out) :: kkt_matrix
		real(wp), dimension (:,:), intent(in) :: A,G
		real(wp), dimension (:) :: gg,hh
		real(wp), allocatable :: v(:)
		real(wp), allocatable :: copy_of_kkt_matrix(:,:)
		integer :: n,m

		m = size(A(:,1))
		n = size(A(1,:))

		allocate(v(m+n))
		allocate(copy_of_kkt_matrix(m+n,m+n))

		kkt_matrix=0
		kkt_matrix(1:n,1:n) = G
		kkt_matrix(n+1:n+m,1:n) = A
		kkt_matrix(1:n,n+1:n+m) = transpose(A)
		copy_of_kkt_matrix=kkt_matrix

		v(1:n) = gg
		v(n+1:n+m)= hh
		
		call gauss_pivot(copy_of_kkt_matrix,v) !solve using gauss

		gg=-v(1:n) !we want to find -p
		hh=v(n+1:n+m)
		
		deallocate(v)
		deallocate(copy_of_kkt_matrix)

	end subroutine 

	



end module matrix_solve

