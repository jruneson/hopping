      module maths
         use types
         implicit real(dp) (a-h,o-z)

         complex(dpc), parameter :: iu = (0.0_dp,1.0_dp)
      contains

         subroutine symevp (a,lda,n,d,ierr)
            implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     This subroutine uses LAPACK DSYEV to
c     diagonalise a real symmetric matrix.
c     -----------------------------------------------------------------
c
            dimension a(lda,n),d(n)
            allocatable :: work(:)
c
            lwork = 34*n
            allocate (work(lwork))
            call dsyev ('V','U',n,a,lda,d,work,lwork,ierr)
            deallocate (work)
         end subroutine

         subroutine hermevp (a,lda,n,d,ierr)
            implicit double precision (a-h,o-z)
            complex(dpc) a,work
c
c     -----------------------------------------------------------------
c     This subroutine uses LAPACK ZHEEV to
c     diagonalise a complex hermitian matrix.
c     First does a workspace query.
c     -----------------------------------------------------------------
c
            dimension a(lda,n),d(n)
            allocatable :: work(:), rwork(:)
c
            allocate(rwork(3*n-2))
            ! Workspace query
            lwork = -1
            allocate (work(1))
            call zheev ('V','U',n,a,lda,d,work,lwork,rwork,ierr)
            lwork = int(work(1))
            deallocate (work)
            ! Set work space
            allocate (work(lwork))
            call zheev ('V','U',n,a,lda,d,work,lwork,rwork,ierr)
            deallocate (work,rwork)
         end subroutine

         subroutine symevpd (a,lda,n,d,ierr)
c
c     ------------------------------------------------------------------
c     This subroutine uses LAPACK DSYEVD to
c     diagonalise a real symmetric matrix.
c     ------------------------------------------------------------------
c
            dimension a(lda,n),d(n)
            allocatable :: work(:),iwork(:)
c
            lwork = 1+6*n+2*n*n
            liwork = 3+5*n
            allocate (work(lwork),iwork(liwork))
            call dsyevd ('V','U',n,a,lda,d,work,lwork,iwork,liwork,ierr)
            deallocate (work,iwork)
         end subroutine


         subroutine gasdev (g,n)
c
c     ------------------------------------------------------------------
c     Generates an array of n normal deviates.
c     ------------------------------------------------------------------
c
            dimension g(n)
            allocatable :: r(:)
            allocate (r(n+1))
            call random_number(r)
c
            twopi = 2*dacos(-1.d0)
            do j = 1,n,2
               x = 1.d0-r(j)! r(j) is in [0,1) and can be exactly 0...
               y = r(j+1)
               rho = dsqrt(-2*dlog(x))! ...in which case rho = 0 (OK).
               phi = twopi*y
               g(j) = rho*dcos(phi)
               if (j .eq. n) exit
               g(j+1) = rho*dsin(phi)
            enddo
            deallocate (r)
         end subroutine

         subroutine sampsun (n,c)
            complex(dpc) c
c
c     ------------------------------------------------------------------
c     Samples a random SU(N) electronic state
c     ------------------------------------------------------------------
c
            dimension c(n)
            allocatable :: p(:),q(:)
            allocate (p(n),q(n))
            call gasdev (p,n)
            call gasdev (q,n)
            s = 0.0
            do j = 1,n
               s = s+p(j)**2+q(j)**2
            enddo
            s = 1.0/sqrt(s)
            do j = 1,n
               c(j) = s*dcmplx(p(j),q(j))
            enddo
            deallocate (p,q)
         end subroutine

         subroutine rnseed(mode)
            implicit integer (a-z)
c
c     ------------------------------------------------------------------
c     Random initial seed for the random_number generator.
c     mode >= 0: use system clock
c     mode < 0: use a fixed seed
c     ------------------------------------------------------------------
c
            allocatable :: seed(:)
c
            call random_seed (size=n)
            allocate (seed(n))
            if (mode.le.0) then
               call system_clock (count)
            else if (mode.gt.0) then
               count = 123456789 * mode
            end if
            do i = 1,n
               seed(i) = i*count
            enddo
            call random_seed (put=seed)
            deallocate (seed)
         end subroutine

         function trace (A)
            dimension A(:,:)
            trace = 0.d0
            n = size(A,1)
            do j=1,n
               trace = trace + A(j,j)
            end do
         end function

         function trnorm (A)
            complex(dpc) :: A, B
            dimension A(:,:)
            allocatable :: B(:,:)
            n = size(A,1)
            allocate(B(n,n))
            B = matmul(conjg(transpose(A)),A)
            trnorm = sqrt(trace(real(B)))
         end function

      end module
