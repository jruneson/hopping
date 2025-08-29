module tul
   use types
   implicit double precision (a-h,o-z)
   public :: pot,force
   real(dp) :: A, B, C, D, E
   integer :: model_number ! 1, 2 or 3
   logical :: branching

contains

   subroutine inittul(model_number_,branching_)
      use system, only : rmass, gam, sys_pot=>pot, sys_grad=>grad
      use system, only : nf, ns, sampnuc_opt, sampnuc_sys=>sampnuc
      logical :: branching_
      if (nf.ne.1) stop 'Tully models require nf=1'
      if (ns.ne.2) stop 'Tully models require ns=2'

      model_number = model_number_
      branching = branching_
      rmass = 2000.d0
      gam = 0.5d0
      if (model_number.eq.1) then
         A = 0.01d0
         B = 1.6d0
         C = 0.005d0
         D = 1.0d0
         sys_pot => pot_tully1
         sys_grad=>grad_tully1   
      else if (model_number.eq.2) then
         A = 0.1d0
         B = 0.28d0
         C = 0.015d0
         D = 0.06d0
         E = 0.05d0
         sys_pot => pot_tully2
         sys_grad=>grad_tully2
      else if (model_number.eq.3) then
         A = 6.d-4
         B = 0.10
         C = 0.90
         sys_pot => pot_tully3
         sys_grad=>grad_tully3   
      else
         stop 'Invalid model number'
      end if
      sampnuc_opt = 'wig'
      sampnuc_sys => sampnuc

      ! call print_potential()

   end subroutine

   subroutine print_potential()
      use system, only : potmat,nac
      real(dp) :: q(1), v(2)
      complex(dpc) :: u(2,2), dnac(1)
      qmin = -10.d0
      qmax = 10.d0
      open(10,file='tully_pot.dat')
      do i=1,100
         q(1) = qmin + (qmax-qmin)*(i-1)/99
         call potmat(q,u,v)
         call nac(q,u,v,1,2,dnac)
         write(10,*) q(1), v, real(dnac)
      end do
      close(10)
   end subroutine

   subroutine pot_tully1(q,V)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: V(:,:)
      V = 0.d0
      q1 = q(1)
      do n=1,2
         if (n == 1) then
            s = sign(1.d0,q1)
         else
            s = -sign(1.d0,q1)
         end if
         V(n,n) = s*A*(1-exp(-B*abs(q1)))
      end do
      V(1,2) = C*exp(-D*q1**2)
      V(2,1) = V(1,2)
   end subroutine 

   subroutine grad_tully1(q,dVdq)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: dVdq(:,:,:)
      dVdq = 0.d0
      dVdq(1,1,1) = A*B*exp(-B*abs(q(1)))
      dVdq(1,2,2) = -dVdq(1,1,1)
      dVdq(1,1,2) = -2.d0*D*q(1)*C*exp(-D*q(1)**2)
      dVdq(1,2,1) = dVdq(1,1,2)
   end subroutine 


   subroutine pot_tully2(q,V)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: V(:,:)
      V = 0.d0
      q1 = q(1)
      V(1,1) = 0.0
      V(2,2) = -A*exp(-B*q1**2) + E
      V(1,2) = C*exp(-D*q1**2)
      V(2,1) = V(1,2)
   end subroutine

   subroutine grad_tully2(q,dVdq)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: dVdq(:,:,:)
      dVdq = 0.d0
      q1 = q(1)
      dVdq(1,2,2) = 2.d0*A*B*q1*exp(-B*q1**2)
      dVdq(1,1,2) = -2.d0*C*D*q1*exp(-D*q1**2)
      dVdq(1,2,1) = dVdq(1,1,2)
   end subroutine
   

   subroutine pot_tully3(q,V)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: V(:,:)
      V = 0.d0
      V(1,1) = A
      V(2,2) = -A
      if (q(1) < 0.d0) then
         V(1,2) = B*exp(C*q(1))
      else
         V(1,2) = B*(2.d0 - exp(-C*q(1)))
      end if
   end subroutine

   subroutine grad_tully3(q,dVdq)
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: dVdq(:,:,:)
      dVdq = 0.d0
      if (q(1) < 0.d0) then
         dVdq(1,2,1) = B*C*exp(C*q(1))
      else
         dVdq(1,2,1) = B*C*exp(-C*q(1))
      end if
      dVdq(1,1,2) = dVdq(1,2,1)
   end subroutine

   subroutine sampnuc (p,q)
      use maths, only : gasdev
      use system, only : nf, sampnuc_opt, gam
!
!     ------------------------------------------------------------------
!     Samples the nuclear momenta (p) and coordinates (q)
!     ------------------------------------------------------------------
!
      real(dp), intent(out) :: p(:),q(:)
!
!     Sample p and q from the Boltzmann distribution
!
      call gasdev (p,nf)
      call gasdev (q,nf)
      pmean = 10.d0
      qmean = -15.0
      ! gamma = 0.5d0
      gamma = gam(1)
      do i=1,nf
         if (sampnuc_opt.eq.'wig') then
            ! Wigner
            deltap = sqrt(gamma/2.0)
            deltaq = sqrt(1.0/(2.0*gamma))
         else
            stop 'Unknown sampnuc_opt'
         end if
         p(i) = deltap*p(i) + pmean
         q(i) = deltaq*q(i) + qmean
      end do
   end subroutine

end module tul
