module hybrd_module

use minpackfunctions

contains

!*==HYBRD.spg  processed by SPAG 6.72Dc at 21:34 on  6 Apr 2019
      SUBROUTINE HYBRD(FCN,N,X,Fvec,Xtol,Maxfev,Ml,Mu,Epsfcn,Diag,Mode, &
                     & Factor,Nprint,Info,Nfev,Fjac,Ldfjac,R,Lr,Qtf,Wa1,&
                     & Wa2,Wa3,Wa4,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
      IMPLICIT NONE
!*--HYBRD6
      INTEGER N , Maxfev , Ml , Mu , Mode , Nprint , Info , Nfev ,      &
            & Ldfjac , Lr, T
      DOUBLE PRECISION Xtol , Epsfcn , Factor, V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q
      DOUBLE PRECISION X(N) , Fvec(N) , Diag(N) , Fjac(Ldfjac,N) ,      &
                     & R(Lr) , Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) ,      &
                     & Wa4(N), Kmatrix(n,n), F(n), Points(2,n/2),Triangles(3,T)

      EXTERNAL FCN
!     **********
!
!     subroutine hybrd
!
!     the purpose of hybrd is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
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
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least maxfev
!         by the end of an iteration.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , jm1 , l , msum , ncfail , ncsuc ,  &
            & nslow1 , nslow2
      INTEGER iwa(1)
      LOGICAL jeval , sing
      DOUBLE PRECISION actred , delta , epsmch , fnorm , fnorm1 , one , &
                     & pnorm , prered , p1 , p5 , p001 , p0001 , ratio ,&
                     & sum , temp , xnorm , zero
      !DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p001 , p0001 , zero/1.0D0 , 1.0D-1 , 5.0D-1 ,&
         & 1.0D-3 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
!
!     check the input parameters for errors.
!
      IF ( N.LE.0 .OR. Xtol.LT.zero .OR. Maxfev.LE.0 .OR. Ml.LT.0 .OR.  &
         & Mu.LT.0 .OR. Factor.LE.zero .OR. Ldfjac.LT.N .OR.            &
         & Lr.LT.(N*(N+1))/2 ) GOTO 300
      IF ( Mode.EQ.2 ) THEN
         DO j = 1 , N
            IF ( Diag(j).LE.zero ) GOTO 300
         ENDDO
      ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      !print *,"Ik ben voor de 1e call van FCN"
      !print *,Fvec
      CALL FCN(N,X,Fvec,iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q) !added variables required for specific fcn calculation
     ! print *,"Ik ben na de 1e call van FCN"
     ! print *,Fvec
      Nfev = 1
      IF ( iflag.LT.0 ) GOTO 300
      fnorm = ENORM(N,Fvec)
!
!     determine the number of calls to fcn needed to compute
!     the jacobian matrix.
!
      msum = MIN0(Ml+Mu+1,N)
!
!     initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     beginning of the outer loop.
!
 100  jeval = .TRUE.
!
!        calculate the jacobian matrix.
!
      iflag = 2
     ! print *,"Ik ben voor de 1e call van FDJAC1"
     ! print *,Fvec
      CALL FDJAC1(FCN,N,X,Fvec,Fjac,Ldfjac,iflag,Ml,Mu,Epsfcn,Wa1,Wa2,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
      Nfev = Nfev + msum
     ! print *,"Ik ben na de 1e call van FDJAC1"
      !print *,Fvec
      IF ( iflag.LT.0 ) GOTO 300
!
!        compute the qr factorization of the jacobian.
!
      CALL QRFAC(N,N,Fjac,Ldfjac,.FALSE.,iwa,1,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      IF ( iter.EQ.1 ) THEN
         IF ( Mode.NE.2 ) THEN
            DO j = 1 , N
               Diag(j) = Wa2(j)
               IF ( Wa2(j).EQ.zero ) Diag(j) = one
            ENDDO
         ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         DO j = 1 , N
            Wa3(j) = Diag(j)*X(j)
         ENDDO
         xnorm = ENORM(N,Wa3)
         delta = Factor*xnorm
         IF ( delta.EQ.zero ) delta = Factor
      ENDIF
!
!        form (q transpose)*fvec and store in qtf.
!
      DO i = 1 , N
         Qtf(i) = Fvec(i)
      ENDDO
      DO j = 1 , N
         IF ( Fjac(j,j).NE.zero ) THEN
            sum = zero
            DO i = j , N
               sum = sum + Fjac(i,j)*Qtf(i)
            ENDDO
            temp = -sum/Fjac(j,j)
            DO i = j , N
               Qtf(i) = Qtf(i) + Fjac(i,j)*temp
            ENDDO
         ENDIF
      ENDDO
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .FALSE.
      DO j = 1 , N
         l = j
         jm1 = j - 1
         IF ( jm1.GE.1 ) THEN
            DO i = 1 , jm1
               R(l) = Fjac(i,j)
               l = l + N - i
            ENDDO
         ENDIF
         R(l) = Wa1(j)
         IF ( Wa1(j).EQ.zero ) sing = .TRUE.
      ENDDO
!
!        accumulate the orthogonal factor in fjac.
!
      CALL QFORM(N,N,Fjac,Ldfjac,Wa1)
!
!        rescale if necessary.
!
      IF ( Mode.NE.2 ) THEN
         DO j = 1 , N
            Diag(j) = DMAX1(Diag(j),Wa2(j))
         ENDDO
      ENDIF
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  IF ( Nprint.GT.0 ) THEN
         iflag = 0
	 !print *,"Ik ben voor de 2e call van FCN"
        ! print *,Fvec
         IF ( MOD(iter-1,Nprint).EQ.0 ) CALL FCN(N,X,Fvec,iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
	 !print *,"Ik ben na de 2e call van FCN"
         !print *,Fvec
         IF ( iflag.LT.0 ) GOTO 300
      ENDIF
!
!           determine the direction p.
!
      CALL DOGLEG(N,R,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      DO j = 1 , N
         Wa1(j) = -Wa1(j)
         Wa2(j) = X(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      ENDDO
      pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      IF ( iter.EQ.1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      !print *,"Ik ben voor de 3e call van FCN"
      !print *,Fvec
      CALL FCN(N,Wa2,Wa4,iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
      !print *,"Ik ben na de 3e call van FCN"
      !print *,Fvec
      Nfev = Nfev + 1
      IF ( iflag.GE.0 ) THEN
         fnorm1 = ENORM(N,Wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         IF ( fnorm1.LT.fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         DO i = 1 , N
            sum = zero
            DO j = i , N
               sum = sum + R(l)*Wa1(j)
               l = l + 1
            ENDDO
            Wa3(i) = Qtf(i) + sum
         ENDDO
         temp = ENORM(N,Wa3)
         prered = zero
         IF ( temp.LT.fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         IF ( prered.GT.zero ) ratio = actred/prered
!
!           update the step bound.
!
         IF ( ratio.GE.p1 ) THEN
            ncfail = 0
            ncsuc = ncsuc + 1
            IF ( ratio.GE.p5 .OR. ncsuc.GT.1 )                          &
               & delta = DMAX1(delta,pnorm/p5)
            IF ( DABS(ratio-one).LE.p1 ) delta = pnorm/p5
         ELSE
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         ENDIF
!
!           test for successful iteration.
!
         IF ( ratio.GE.p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
            DO j = 1 , N
               X(j) = Wa2(j)
               Wa2(j) = Diag(j)*X(j)
               Fvec(j) = Wa4(j)
            ENDDO
            xnorm = ENORM(N,Wa2)
            fnorm = fnorm1
            iter = iter + 1
         ENDIF
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         IF ( actred.GE.p001 ) nslow1 = 0
         IF ( jeval ) nslow2 = nslow2 + 1
         IF ( actred.GE.p1 ) nslow2 = 0
!
!           test for convergence.
!
         IF ( delta.LE.Xtol*xnorm .OR. fnorm.EQ.zero ) Info = 1
         IF ( Info.EQ.0 ) THEN
!
!           tests for termination and stringent tolerances.
!
            IF ( Nfev.GE.Maxfev ) Info = 2
            IF ( p1*DMAX1(p1*delta,pnorm).LE.epsmch*xnorm ) Info = 3
            IF ( nslow2.EQ.5 ) Info = 4
            IF ( nslow1.EQ.10 ) Info = 5
            IF ( Info.EQ.0 ) THEN
!
!           criterion for recalculating jacobian approximation
!           by forward differences.
!
               IF ( ncfail.EQ.2 ) GOTO 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               DO j = 1 , N
                  sum = zero
                  DO i = 1 , N
                     sum = sum + Fjac(i,j)*Wa4(i)
                  ENDDO
                  Wa2(j) = (sum-Wa3(j))/pnorm
                  Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
                  IF ( ratio.GE.p0001 ) Qtf(j) = sum
               ENDDO
!
!           compute the qr factorization of the updated jacobian.
!
               CALL R1UPDT(N,N,R,Lr,Wa1,Wa2,Wa3,sing)
               CALL R1MPYQ(N,N,Fjac,Ldfjac,Wa2,Wa3)
               CALL R1MPYQ(1,N,Qtf,1,Wa2,Wa3)
!
!           end of the inner loop.
!
               jeval = .FALSE.
!
!        end of the outer loop.
!
               GOTO 200
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 300  IF ( iflag.LT.0 ) Info = iflag
      iflag = 0
      !print *,"Ik ben voor de 4e call van FCN"
      !print *,Fvec
      IF ( Nprint.GT.0 ) CALL FCN(N,X,Fvec,iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
      !print *,"Ik ben na de 4e call van FCN"
     ! print *,Fvec
!
!     last card of subroutine hybrd.
!
      END SUBROUTINE HYBRD

end module hybrd_module

