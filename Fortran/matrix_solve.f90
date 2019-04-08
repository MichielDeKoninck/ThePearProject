module matrix_solve
	implicit none
	save
	!private
	!public eig
	integer, parameter :: prec = selected_real_kind(15) !define working precision which can be adapted. (6 for single, 15 for double)
	integer k,j,i

contains
	subroutine gauss(A, v, m)
		
		
		integer, intent(in) :: m
		real(prec), dimension(m,m) :: A
		real(prec), dimension(m) :: v
		real(prec), dimension(m) :: y
		real(prec), dimension(m,m) :: U
		real(prec), dimension(m,m) :: L
		real(prec), dimension(m,m) :: copy_A
	
		copy_A = A

		! Ik heb U en L ook berekend om zo de efficiënt berekende A te kunnen vergelijken met L en U.
		U = A
		do i = 1,m
			do j = 1,m
				L(i,j) = 0
			enddo
		enddo
		! Ik maak van L op het begin een eenheidsmatrix.
		do i = 1,m
			L(i,i) = 1
		enddo

		do k = 1,m-1
			do j = k+1,m
				L(j,k) = U(j,k) / U(k,k)
				U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
			enddo
		enddo

		! k in de laatste lijn (voor efficiënte A) aangepast naar k+1 (zie discussieforum).
		do k = 1,m-1
			do j = k+1,m
				A(j,k) = A(j,k) / A(k,k)
				A(j,k+1:m) = A(j,k+1:m) - A(j,k)*A(k,k+1:m)
			enddo
		enddo

		! y stel ik op het begin gelijk aan v voor de voorwaartse substitutie toe te passen.
		y = v
		! Voorwaartse substitutie waarbij we i = 1 niet moeten doen, want y_1 = v_1.
		do i = 2,m
			do j = 1,i-1
				y(i) = y(i) - A(i,j)*y(j)
			enddo
		enddo
		
		! De berekende x, wat dus v is hier, stel ik op het begin gelijk aan y/a_ii voor de 	
		! achterwaartse substitutie toe te passen.
		do i = 1,m
			v(i) = y(i)/A(i,i)
		enddo
		
		! Achterwaartse substitutie toepassen om de uiteindelijk oplossing van Ax = v te berekenen.
		do i = m-1,1,-1
			do j = i+1,m
				v(i) = v(i) - (A(i,j)*v(j))/A(i,i)
			enddo
		enddo
		
	end subroutine gauss

	subroutine gauss_pivot(A, v, m)
	
		integer(prec) i
		integer(prec) j
		integer(prec) k
		integer m
		integer(prec) g
		real(prec), dimension(m,m) :: A
		real(prec), dimension(m) :: v
		real(prec), dimension(m) :: y
		real(prec), dimension(m,m) :: copy_A
		real(prec), dimension(m) :: v_copy
		real(prec), dimension(m) :: TEMP_SWAP_A
		integer(prec), dimension(m) :: p
		real(prec) maxwaarde
		real(prec) temp
		integer(prec) ix
		integer(prec) p_temp

		copy_A = A

		! Ik stel de vector p van lengte m die de nuttige informatie van de permutatiematrix bevat.
		p  = (/ (1.0 * g, g = 1,m) /)

		! Ik zoek de rij waarin de gezochte pivot zit.
		do k = 1,m-1
			! Ik moet wel voor elke nieuwe kolom de oude maxwaarde van de vorige kolom resetten.
			maxwaarde = 0
			do i = k,m
				temp = abs(A(i,k))
				if (temp > maxwaarde) then 
					maxwaarde = temp
					! ix is de rij waar grootste pivot zit in bepaalde kolom van  						! submatrix. Omdat ik de maxwaarde bij elke nieuwe kolom al reset,  						! moet ik ix niet elke keer opnieuw resetten.		
					ix = i  	
				endif
			enddo
			! Hier swap ik rij ix met rij k van de efficiënt berekende matrix A. 
			TEMP_SWAP_A = A(ix,1:m)
			A(ix,1:m) = A(k,1:m)
			A(k,1:m) = TEMP_SWAP_A
			! Hier swap ik rij ix met rij k van de permutatiematrix.
			p_temp = p(k)
			p(k) = p(ix)
			p(ix) = p_temp

			! k in de laatste lijn (voor efficiënte A) aangepast naar k+1 (zie discussieforum).
			do j = k+1,m
				A(j,k) = A(j,k)/A(k,k)
				A(j,k+1:m) = A(j,k+1:m) - A(j,k)*A(k,k+1:m)
			enddo
		enddo

		v_copy = v

		! Vanaf hier bereken ik de nieuwe v door de oude v te permuteren met de bekomen 		! permutatiematrix, namelijk P*v = v'. Nu kan ik exact op dezelfde manier werken als bij het 			! Gauss-algoritme hierboven.
		do i = 1,m
			v(i) = v_copy(p(i))
		enddo

		! y stel ik op het begin gelijk aan v voor de voorwaartse substitutie toe te passen.
		y = v

		! Voorwaartse substitutie waarbij we i = 1 niet moeten doen, want y_1 = v_1.
		do i = 2,m
			do j = 1,i-1
				y(i) = y(i) - A(i,j)*y(j)
			enddo
		enddo

		! De berekende x, wat dus v is hier, stel ik op het begin gelijk aan y/a_ii voor de 	
		! achterwaartse substitutie toe te passen.
		do i = 1,m
			v(i) = y(i)/A(i,i)
		enddo
		
		! Achterwaartse substitutie toepassen om de uiteindelijk oplossing van Ax = v te berekenen.
		do i = m-1,1,-1
			do j = i+1,m
				v(i) = v(i) - (A(i,j)*v(j))/A(i,i)
			enddo
		enddo

	end subroutine gauss_pivot

	
end module matrix_solve

