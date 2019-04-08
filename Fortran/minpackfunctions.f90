module minpackfunctions

!module containing all subroutines needed to use hybrd

contains

!____________________________________________________________________________________________________________________________
!*==DOGLEG.spg  processed by SPAG 6.72Dc at 21:44 on  6 Apr 2019
      SUBROUTINE DOGLEG(N,R,Lr,Diag,Qtb,Delta,X,Wa1,Wa2)
      IMPLICIT NONE
!*--DOGLEG4
      INTEGER N , Lr
      DOUBLE PRECISION Delta
      DOUBLE PRECISION R(Lr) , Diag(N) , Qtb(N) , X(N) , Wa1(N) , Wa2(N)
!     **********
!
!     subroutine dogleg
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta, the
!     problem is to determine the convex combination x of the
!     gauss-newton and scaled gradient directions that minimizes
!     (a*x - b) in the least squares sense, subject to the
!     restriction that the euclidean norm of d*x be at most delta.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization of a. that is, if a = q*r, where q has
!     orthogonal columns and r is an upper triangular matrix,
!     then dogleg expects the full upper triangle of r and
!     the first n components of (q transpose)*b.
!
!     the subroutine statement is
!
!       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an input array of length lr which must contain the upper
!         triangular matrix r stored by rows.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       x is an output array of length n which contains the desired
!         convex combination of the gauss-newton direction and the
!         scaled gradient direction.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jj , jp1 , k , l
      DOUBLE PRECISION alpha , bnorm , epsmch , gnorm , one , qnorm ,   &
                     & sgnorm , sum , temp , zero
      !DOUBLE PRECISION DPMPAR , ENORM
      DATA one , zero/1.0D0 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
!     first, calculate the gauss-newton direction.
!
      jj = (N*(N+1))/2 + 1
      DO k = 1 , N
         j = N - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         IF ( N.GE.jp1 ) THEN
            DO i = jp1 , N
               sum = sum + R(l)*X(i)
               l = l + 1
            ENDDO
         ENDIF
         temp = R(jj)
         IF ( temp.EQ.zero ) THEN
            l = j
            DO i = 1 , j
               temp = DMAX1(temp,DABS(R(l)))
               l = l + N - i
            ENDDO
            temp = epsmch*temp
            IF ( temp.EQ.zero ) temp = epsmch
         ENDIF
         X(j) = (Qtb(j)-sum)/temp
      ENDDO
!
!     test whether the gauss-newton direction is acceptable.
!
      DO j = 1 , N
         Wa1(j) = zero
         Wa2(j) = Diag(j)*X(j)
      ENDDO
      qnorm = ENORM(N,Wa2)
      IF ( qnorm.GT.Delta ) THEN
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
         l = 1
         DO j = 1 , N
            temp = Qtb(j)
            DO i = j , N
               Wa1(i) = Wa1(i) + R(l)*temp
               l = l + 1
            ENDDO
            Wa1(j) = Wa1(j)/Diag(j)
         ENDDO
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
         gnorm = ENORM(N,Wa1)
         sgnorm = zero
         alpha = Delta/qnorm
         IF ( gnorm.NE.zero ) THEN
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
            DO j = 1 , N
               Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
            ENDDO
            l = 1
            DO j = 1 , N
               sum = zero
               DO i = j , N
                  sum = sum + R(l)*Wa1(i)
                  l = l + 1
               ENDDO
               Wa2(j) = sum
            ENDDO
            temp = ENORM(N,Wa2)
            sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
            alpha = zero
            IF ( sgnorm.LT.Delta ) THEN
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
               bnorm = ENORM(N,Qtb)
               temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
               temp = temp - (Delta/qnorm)*(sgnorm/Delta)               &
                    & **2 + DSQRT((temp-(Delta/qnorm))                  &
                    & **2+(one-(Delta/qnorm)**2)*(one-(sgnorm/Delta)**2)&
                    & )
               alpha = ((Delta/qnorm)*(one-(sgnorm/Delta)**2))/temp
            ENDIF
         ENDIF
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
         temp = (one-alpha)*DMIN1(sgnorm,Delta)
         DO j = 1 , N
            X(j) = temp*Wa1(j) + alpha*X(j)
         ENDDO
      ENDIF
!
!     last card of subroutine dogleg.
!
      END SUBROUTINE DOGLEG

! _____________________________________________________________________________________________________________________________

!*==R1UPDT.spg  processed by SPAG 6.72Dc at 21:52 on  6 Apr 2019
      SUBROUTINE R1UPDT(M,N,S,Ls,U,V,W,Sing)
      IMPLICIT NONE
!*--R1UPDT4
      INTEGER M , N , Ls
      LOGICAL Sing
      DOUBLE PRECISION S(Ls) , U(M) , V(N) , W(M)
!     **********
!
!     subroutine r1updt
!
!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that
!
!                   t
!           (s + u*v )*q
!
!     is again lower trapezoidal.
!
!     this subroutine determines q as the product of 2*(n - 1)
!     transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     where gv(i), gw(i) are givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.
!
!     the subroutine statement is
!
!       subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of s.
!
!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.
!
!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.
!
!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.
!
!       u is an input array of length m which must contain the
!         vector u.
!
!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more,
!     john l. nazareth
!
!     **********
      INTEGER i , j , jj , l , nmj , nm1
      DOUBLE PRECISION cos , cotan , giant , one , p5 , p25 , sin ,     &
                     & tan , tau , temp , zero
      !DOUBLE PRECISION DPMPAR
      DATA one , p5 , p25 , zero/1.0D0 , 5.0D-1 , 2.5D-1 , 0.0D0/
!
!     giant is the largest magnitude.
!
      giant = DPMPAR(3)
!
!     initialize the diagonal element pointer.
!
      jj = (N*(2*M-N+1))/2 - (M-N)
!
!     move the nontrivial part of the last column of s into w.
!
      l = jj
      DO i = N , M
         W(i) = S(l)
         l = l + 1
      ENDDO
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
      nm1 = N - 1
      IF ( nm1.GE.1 ) THEN
         DO nmj = 1 , nm1
            j = N - nmj
            jj = jj - (M-j+1)
            W(j) = zero
            IF ( V(j).NE.zero ) THEN
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
               IF ( DABS(V(N)).GE.DABS(V(j)) ) THEN
                  tan = V(j)/V(N)
                  cos = p5/DSQRT(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               ELSE
                  cotan = V(N)/V(j)
                  sin = p5/DSQRT(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  IF ( DABS(cos)*giant.GT.one ) tau = one/cos
               ENDIF
!
!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.
!
               V(N) = sin*V(j) + cos*V(N)
               V(j) = tau
!
!        apply the transformation to s and extend the spike in w.
!
               l = jj
               DO i = j , M
                  temp = cos*S(l) - sin*W(i)
                  W(i) = sin*S(l) + cos*W(i)
                  S(l) = temp
                  l = l + 1
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     add the spike from the rank 1 update to w.
!
      DO i = 1 , M
         W(i) = W(i) + V(N)*U(i)
      ENDDO
!
!     eliminate the spike.
!
      Sing = .FALSE.
      IF ( nm1.GE.1 ) THEN
         DO j = 1 , nm1
            IF ( W(j).NE.zero ) THEN
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
               IF ( DABS(S(jj)).GE.DABS(W(j)) ) THEN
                  tan = W(j)/S(jj)
                  cos = p5/DSQRT(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               ELSE
                  cotan = S(jj)/W(j)
                  sin = p5/DSQRT(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  IF ( DABS(cos)*giant.GT.one ) tau = one/cos
               ENDIF
!
!        apply the transformation to s and reduce the spike in w.
!
               l = jj
               DO i = j , M
                  temp = cos*S(l) + sin*W(i)
                  W(i) = -sin*S(l) + cos*W(i)
                  S(l) = temp
                  l = l + 1
               ENDDO
!
!        store the information necessary to recover the
!        givens rotation.
!
               W(j) = tau
            ENDIF
!
!        test for zero diagonal elements in the output s.
!
            IF ( S(jj).EQ.zero ) Sing = .TRUE.
            jj = jj + (M-j+1)
         ENDDO
      ENDIF
!
!     move w back into the last column of the output s.
!
      l = jj
      DO i = N , M
         S(l) = W(i)
         l = l + 1
      ENDDO
      IF ( S(jj).EQ.zero ) Sing = .TRUE.
!
!     last card of subroutine r1updt.
!
      END SUBROUTINE R1UPDT





! _____________________________________________________________________________________________________________________________





!*==R1MPYQ.spg  processed by SPAG 6.72Dc at 21:51 on  6 Apr 2019
      SUBROUTINE R1MPYQ(M,N,A,Lda,V,W)
      IMPLICIT NONE
!*--R1MPYQ4
      INTEGER M , N , Lda
      DOUBLE PRECISION A(Lda,N) , V(N) , W(N)
!     **********
!
!     subroutine r1mpyq
!
!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine r1mpyq(m,n,a,lda,v,w)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.
!
!     subroutines called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , nmj , nm1
      DOUBLE PRECISION cos , one , sin , temp
      DATA one/1.0D0/
!
!     apply the first set of givens rotations to a.
!
      nm1 = N - 1
      IF ( nm1.GE.1 ) THEN
         DO nmj = 1 , nm1
            j = N - nmj
            IF ( DABS(V(j)).GT.one ) cos = one/V(j)
            IF ( DABS(V(j)).GT.one ) sin = DSQRT(one-cos**2)
            IF ( DABS(V(j)).LE.one ) sin = V(j)
            IF ( DABS(V(j)).LE.one ) cos = DSQRT(one-sin**2)
            DO i = 1 , M
               temp = cos*A(i,j) - sin*A(i,N)
               A(i,N) = sin*A(i,j) + cos*A(i,N)
               A(i,j) = temp
            ENDDO
         ENDDO
!
!     apply the second set of givens rotations to a.
!
         DO j = 1 , nm1
            IF ( DABS(W(j)).GT.one ) cos = one/W(j)
            IF ( DABS(W(j)).GT.one ) sin = DSQRT(one-cos**2)
            IF ( DABS(W(j)).LE.one ) sin = W(j)
            IF ( DABS(W(j)).LE.one ) cos = DSQRT(one-sin**2)
            DO i = 1 , M
               temp = cos*A(i,j) + sin*A(i,N)
               A(i,N) = -sin*A(i,j) + cos*A(i,N)
               A(i,j) = temp
            ENDDO
         ENDDO
      ENDIF
!
!     last card of subroutine r1mpyq.
!
      END SUBROUTINE R1MPYQ


! _____________________________________________________________________________________________________________________________


!*==QRFAC.spg  processed by SPAG 6.72Dc at 21:50 on  6 Apr 2019
      SUBROUTINE QRFAC(M,N,A,Lda,Pivot,Ipvt,Lipvt,Rdiag,Acnorm,Wa)
      IMPLICIT NONE
!*--QRFAC4
      INTEGER M , N , Lda , Lipvt
      INTEGER Ipvt(Lipvt)
      LOGICAL Pivot
      DOUBLE PRECISION A(Lda,N) , Rdiag(N) , Acnorm(N) , Wa(N)
!     **********
!
!     subroutine qrfac
!
!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!     the subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dmax1,dsqrt,min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jp1 , k , kmax , minmn
      DOUBLE PRECISION ajnorm , epsmch , one , p05 , sum , temp , zero
      !DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p05 , zero/1.0D0 , 5.0D-2 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
!     compute the initial column norms and initialize several arrays.
!
      DO j = 1 , N
         Acnorm(j) = ENORM(M,A(1,j))
         Rdiag(j) = Acnorm(j)
         Wa(j) = Rdiag(j)
         IF ( Pivot ) Ipvt(j) = j
      ENDDO
!
!     reduce a to r with householder transformations.
!
      minmn = MIN0(M,N)
      DO j = 1 , minmn
         IF ( Pivot ) THEN
!
!        bring the column of largest norm into the pivot position.
!
            kmax = j
            DO k = j , N
               IF ( Rdiag(k).GT.Rdiag(kmax) ) kmax = k
            ENDDO
            IF ( kmax.NE.j ) THEN
               DO i = 1 , M
                  temp = A(i,j)
                  A(i,j) = A(i,kmax)
                  A(i,kmax) = temp
               ENDDO
               Rdiag(kmax) = Rdiag(j)
               Wa(kmax) = Wa(j)
               k = Ipvt(j)
               Ipvt(j) = Ipvt(kmax)
               Ipvt(kmax) = k
            ENDIF
         ENDIF
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = ENORM(M-j+1,A(j,j))
         IF ( ajnorm.NE.zero ) THEN
            IF ( A(j,j).LT.zero ) ajnorm = -ajnorm
            DO i = j , M
               A(i,j) = A(i,j)/ajnorm
            ENDDO
            A(j,j) = A(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
            jp1 = j + 1
            IF ( N.GE.jp1 ) THEN
               DO k = jp1 , N
                  sum = zero
                  DO i = j , M
                     sum = sum + A(i,j)*A(i,k)
                  ENDDO
                  temp = sum/A(j,j)
                  DO i = j , M
                     A(i,k) = A(i,k) - temp*A(i,j)
                  ENDDO
                  IF ( .NOT.(.NOT.Pivot .OR. Rdiag(k).EQ.zero) ) THEN
                     temp = A(j,k)/Rdiag(k)
                     Rdiag(k) = Rdiag(k)*DSQRT(DMAX1(zero,one-temp**2))
                     IF ( p05*(Rdiag(k)/Wa(k))**2.LE.epsmch ) THEN
                        Rdiag(k) = ENORM(M-j,A(jp1,k))
                        Wa(k) = Rdiag(k)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         Rdiag(j) = -ajnorm
      ENDDO
!
!     last card of subroutine qrfac.
!
      END SUBROUTINE QRFAC








! _____________________________________________________________________________________________________________________________







!*==QFORM.spg  processed by SPAG 6.72Dc at 21:50 on  6 Apr 2019
      SUBROUTINE QFORM(M,N,Q,Ldq,Wa)
      IMPLICIT NONE
!*--QFORM4
      INTEGER M , N , Ldq
      DOUBLE PRECISION Q(Ldq,M) , Wa(M)
!     **********
!
!     subroutine qform
!
!     this subroutine proceeds from the computed qr factorization of
!     an m by n matrix a to accumulate the m by m orthogonal matrix
!     q from its factored form.
!
!     the subroutine statement is
!
!       subroutine qform(m,n,q,ldq,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a and the order of q.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       q is an m by m array. on input the full lower trapezoid in
!         the first min(m,n) columns of q contains the factored form.
!         on output q has been accumulated into a square matrix.
!
!       ldq is a positive integer input variable not less than m
!         which specifies the leading dimension of the array q.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       fortran-supplied ... min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jm1 , k , l , minmn , np1
      DOUBLE PRECISION one , sum , temp , zero
      DATA one , zero/1.0D0 , 0.0D0/
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
      minmn = MIN0(M,N)
      IF ( minmn.GE.2 ) THEN
         DO j = 2 , minmn
            jm1 = j - 1
            DO i = 1 , jm1
               Q(i,j) = zero
            ENDDO
         ENDDO
      ENDIF
!
!     initialize remaining columns to those of the identity matrix.
!
      np1 = N + 1
      IF ( M.GE.np1 ) THEN
         DO j = np1 , M
            DO i = 1 , M
               Q(i,j) = zero
            ENDDO
            Q(j,j) = one
         ENDDO
      ENDIF
!
!     accumulate q from its factored form.
!
      DO l = 1 , minmn
         k = minmn - l + 1
         DO i = k , M
            Wa(i) = Q(i,k)
            Q(i,k) = zero
         ENDDO
         Q(k,k) = one
         IF ( Wa(k).NE.zero ) THEN
            DO j = k , M
               sum = zero
               DO i = k , M
                  sum = sum + Q(i,j)*Wa(i)
               ENDDO
               temp = sum/Wa(k)
               DO i = k , M
                  Q(i,j) = Q(i,j) - temp*Wa(i)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
!     last card of subroutine qform.
!
      END SUBROUTINE QFORM






! _____________________________________________________________________________________________________________________________





!*==FDJAC1.spg  processed by SPAG 6.72Dc at 21:49 on  6 Apr 2019
      SUBROUTINE FDJAC1(FCN,N,X,Fvec,Fjac,Ldfjac,Iflag,Ml,Mu,Epsfcn,Wa1,&
                      & Wa2,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
      IMPLICIT NONE
!*--FDJAC15
      INTEGER N , Ldfjac , Iflag , Ml , Mu, T
      DOUBLE PRECISION Epsfcn,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q
      DOUBLE PRECISION X(N) , Fvec(N) , Fjac(Ldfjac,N) , Wa1(N) , Wa2(N), Kmatrix(n,n), F(n), Points(2,n/2),Triangles(3,T)
!     **********
!
!     subroutine fdjac1
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
!                         wa1,wa2)
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
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
!         least n, then the jacobian is considered dense, and wa2 is
!         not referenced.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , k , msum
      DOUBLE PRECISION eps , epsmch , h , temp , zero
      !DOUBLE PRECISION DPMPAR
      DATA zero/0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      eps = DSQRT(DMAX1(Epsfcn,epsmch))
      msum = Ml + Mu + 1
      IF ( msum.LT.N ) THEN
!
!        computation of banded approximate jacobian.
!
         DO k = 1 , msum
            DO j = k , N , msum
               Wa2(j) = X(j)
               h = eps*DABS(Wa2(j))
               IF ( h.EQ.zero ) h = eps
               X(j) = Wa2(j) + h
            ENDDO
            CALL FCN(N,X,Wa1,Iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
            IF ( Iflag.LT.0 ) GOTO 99999
            DO j = k , N , msum
               X(j) = Wa2(j)
               h = eps*DABS(Wa2(j))
               IF ( h.EQ.zero ) h = eps
               DO i = 1 , N
                  Fjac(i,j) = zero
                  IF ( i.GE.j-Mu .AND. i.LE.j+Ml ) Fjac(i,j)            &
                     & = (Wa1(i)-Fvec(i))/h
               ENDDO
            ENDDO
         ENDDO
      ELSE
!
!        computation of dense approximate jacobian.
!
         DO j = 1 , N
            temp = X(j)
            h = eps*DABS(temp)
            IF ( h.EQ.zero ) h = eps
            X(j) = temp + h
            CALL FCN(N,X,Wa1,Iflag,Kmatrix,F,Points,Triangles,T,V_mu,K_mv,K_mu,K_mfu,V_mfv,r_q)
            IF ( Iflag.LT.0 ) GOTO 99999
            X(j) = temp
            DO i = 1 , N
               Fjac(i,j) = (Wa1(i)-Fvec(i))/h
            ENDDO
         ENDDO
      ENDIF
!
!     last card of subroutine fdjac1.
!
99999 END SUBROUTINE FDJAC1


! _____________________________________________________________________________________________________________________________





!*==ENORM.spg  processed by SPAG 6.72Dc at 21:48 on  6 Apr 2019
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      IMPLICIT NONE
!*--ENORM4
      INTEGER N
      DOUBLE PRECISION X(N)
!     **********
!
!     function enorm
!
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i
      DOUBLE PRECISION agiant , floatn , one , rdwarf , rgiant , s1 ,   &
                     & s2 , s3 , xabs , x1max , x3max , zero
      DATA one , zero , rdwarf , rgiant/1.0D0 , 0.0D0 , 3.834D-20 ,     &
         & 1.304D19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = N
      agiant = rgiant/floatn
      DO i = 1 , N
         xabs = DABS(X(i))
         IF ( xabs.GT.rdwarf .AND. xabs.LT.agiant ) THEN
!
!           sum for intermediate components.
!
            s2 = s2 + xabs**2
         ELSEIF ( xabs.LE.rdwarf ) THEN
!
!              sum for small components.
!
            IF ( xabs.LE.x3max ) THEN
               IF ( xabs.NE.zero ) s3 = s3 + (xabs/x3max)**2
            ELSE
               s3 = one + s3*(x3max/xabs)**2
               x3max = xabs
            ENDIF
!
!              sum for large components.
!
         ELSEIF ( xabs.LE.x1max ) THEN
            s1 = s1 + (xabs/x1max)**2
         ELSE
            s1 = one + s1*(x1max/xabs)**2
            x1max = xabs
         ENDIF
      ENDDO
!
!     calculation of norm.
!
      IF ( s1.NE.zero ) THEN
         ENORM = x1max*DSQRT(s1+(s2/x1max)/x1max)
      ELSEIF ( s2.EQ.zero ) THEN
         ENORM = x3max*DSQRT(s3)
      ELSE
         IF ( s2.GE.x3max ) ENORM = DSQRT(s2*(one+(x3max/s2)*(x3max*s3))&
                                  & )
         IF ( s2.LT.x3max ) ENORM = DSQRT(x3max*((s2/x3max)+(x3max*s3)))
      ENDIF
!
!     last card of function enorm.
!
      END FUNCTION ENORM






! _____________________________________________________________________________________________________________________________








	!*==DPMPAR.spg  processed by SPAG 6.72Dc at 21:47 on  6 Apr 2019
      DOUBLE PRECISION FUNCTION DPMPAR(I)
      IMPLICIT NONE
!*--DPMPAR4
      INTEGER I
!     **********
!
!     Function dpmpar
!
!     This function provides double precision machine parameters
!     when the appropriate set of data statements is activated (by
!     removing the c from column 1) and all other data statements are
!     rendered inactive. Most of the parameter values were obtained
!     from the corresponding Bell Laboratories Port Library function.
!
!     The function statement is
!
!       double precision function dpmpar(i)
!
!     where
!
!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. If the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are
!
!         dpmpar(1) = b**(1 - t), the machine precision,
!
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,
!
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
!
!     Argonne National Laboratory. MINPACK Project. November 1996.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
!
!     **********
      INTEGER mcheps(4)
      INTEGER minmag(4)
      INTEGER maxmag(4)
      DOUBLE PRECISION dmach(3)
      EQUIVALENCE (dmach(1),mcheps(1))
      EQUIVALENCE (dmach(2),minmag(1))
      EQUIVALENCE (dmach(3),maxmag(1))
!
!     Machine constants for the IBM 360/370 series,
!     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
!     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
!
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /
!     data minmag(1),minmag(2) / z00100000, z00000000 /
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
!
!     Machine constants for the Honeywell 600/6000 series.
!
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
!
!     Machine constants for the CDC 6000/7000 series.
!
!     data mcheps(1) / 15614000000000000000b /
!     data mcheps(2) / 15010000000000000000b /
!
!     data minmag(1) / 00604000000000000000b /
!     data minmag(2) / 00000000000000000000b /
!
!     data maxmag(1) / 37767777777777777777b /
!     data maxmag(2) / 37167777777777777777b /
!
!     Machine constants for the PDP-10 (KA processor).
!
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
!
!     Machine constants for the PDP-10 (KI processor).
!
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
!
!     Machine constants for the PDP-11.
!
!     data mcheps(1),mcheps(2) /   9472,      0 /
!     data mcheps(3),mcheps(4) /      0,      0 /
!
!     data minmag(1),minmag(2) /    128,      0 /
!     data minmag(3),minmag(4) /      0,      0 /
!
!     data maxmag(1),maxmag(2) /  32767,     -1 /
!     data maxmag(3),maxmag(4) /     -1,     -1 /
!
!     Machine constants for the Burroughs 6700/7700 systems.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o7770000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o7777777777777777 /
!
!     Machine constants for the Burroughs 5700 system.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o0000000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o0007777777777777 /
!
!     Machine constants for the Burroughs 1700 system.
!
!     data mcheps(1) / zcc6800000 /
!     data mcheps(2) / z000000000 /
!
!     data minmag(1) / zc00800000 /
!     data minmag(2) / z000000000 /
!
!     data maxmag(1) / zdffffffff /
!     data maxmag(2) / zfffffffff /
!
!     Machine constants for the Univac 1100 series.
!
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
!
!     Machine constants for the Data General Eclipse S/200.
!
!     Note - it may be appropriate to include the following card -
!     static dmach(3)
!
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
!     data mcheps/32020k,3*0/
!
!     Machine constants for the Harris 220.
!
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /
!     data minmag(1),minmag(2) / '20000000, '00000201 /
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /
!
!     Machine constants for the Cray-1.
!
!     data mcheps(1) / 0376424000000000000000b /
!     data mcheps(2) / 0000000000000000000000b /
!
!     data minmag(1) / 0200034000000000000000b /
!     data minmag(2) / 0000000000000000000000b /
!
!     data maxmag(1) / 0577777777777777777777b /
!     data maxmag(2) / 0000007777777777777776b /
!
!     Machine constants for the Prime 400.
!
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
!
!     Machine constants for the VAX-11.
!
!     data mcheps(1),mcheps(2) /   9472,  0 /
!     data minmag(1),minmag(2) /    128,  0 /
!     data maxmag(1),maxmag(2) / -32769, -1 /
!
!     Machine constants for IEEE machines.
!
      DATA dmach(1)/2.22044604926D-16/
      DATA dmach(2)/2.22507385852D-308/
      DATA dmach(3)/1.79769313485D+308/
!
      DPMPAR = dmach(I)
!
!     Last card of function dpmpar.
!
      END FUNCTION DPMPAR

end module minpackfunctions

