module fmo
   use types
   implicit real(dp) (a-h,o-z)

   ! FMO system
   integer :: nf_bath
   allocatable :: Vconst(:,:) ! Constant part of potential matrix
   allocatable :: omega(:)    ! Phonon frequencies
   allocatable :: rkappa(:)   ! On-diagonal coupling strengths

   ! Constants 
   real(dp), parameter :: cmm1 = 4.556335d-6


contains

   subroutine initfmo()
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc, sys_pot=>pot
      use system, only : sys_nac=>nac, sys_h0=>h0, sys_grad=>grad
      use system, only : ns, nf, rmass, gam
      ! dimension :: Vconst_(ns_,ns_), omega_(nf_), kappa_(nf_)
!
!     ------------------------------------------------------------------ 
!     Initialise FMO module
!     ------------------------------------------------------------------ 
!
      if (ns .ne. 7) stop 'FMO requires 7 states'
      if (mod(nf,ns) .ne. 0) stop 'nf must be divisible by ns'
      nf_bath = nf/ns
      allocate(Vconst(ns,ns),omega(nf_bath),rkappa(nf_bath))
      Vconst(1,2) = -87.7
      Vconst(1,3) =   5.5
      Vconst(1,4) = - 5.9
      Vconst(1,5) =   6.7
      Vconst(1,6) = -13.7
      Vconst(1,7) = - 9.9
      Vconst(2,3) =  30.8
      Vconst(2,4) =   8.2
      Vconst(2,5) =   0.7
      Vconst(2,6) =  11.8
      Vconst(2,7) =   4.3
      Vconst(3,4) = -53.5
      Vconst(3,5) = - 2.2
      Vconst(3,6) = - 9.6
      Vconst(3,7) =   6.0
      Vconst(4,5) = -70.7
      Vconst(4,6) = -17.0
      Vconst(4,7) = -63.3
      Vconst(5,6) =  81.1
      Vconst(5,7) = - 1.3
      Vconst(6,7) =  39.7
      Vconst = Vconst + transpose(Vconst)
      Vconst(1,1) = 12410.
      Vconst(2,2) = 12530.
      Vconst(3,3) = 12210.
      Vconst(4,4) = 12320.
      Vconst(5,5) = 12480.
      Vconst(6,6) = 12630.
      Vconst(7,7) = 12440.
      shift = 12210.
      do n=1,ns
         Vconst(n,n) = Vconst(n,n) - shift
      end do
      Vconst = Vconst*cmm1
!
!     Discretize bath
!
      rlamda = 35.0*cmm1
      omegac = 106.14*cmm1
      fac = sqrt(2*rlamda/nf_bath)
      rkappa = 0.0
      do i = 1,nf_bath
         omega(i) = omegac*tan(0.5*pi*(i-0.5)/nf_bath)
         rkappa(i) = fac*omega(i)
      end do
!
!     Mass and width parameter
!
      rmass = 1.d0
      do i=1,nf
         j = mod(i-1,nf_bath)+1
        gam(i) = rmass(i) * omega(j)
         ! TEST: classical density
         ! gam(i) = beta * 0.5 * rmass(i) * omega(j)**2
         ! ! Wigner 
         ! x = beta*omega(j)/2.0
         ! y = tanh(x)/x
         ! gam(i) = gam(i) * y
         ! END TEST
      end do

!
      sys_force0=>force0
      sys_force=>force
      sys_sampnuc=>sampnuc
      sys_pot=>pot
      sys_grad=>gradient_dia
      sys_nac=>nac
      sys_h0=>h0
   end subroutine

   subroutine force0 (q,f0)
!
!     ------------------------------------------------------------------ 
!     Calculates the state-independent (CPA) force f0 = - omega**2 * q
!     ------------------------------------------------------------------ 
!  
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: f0(:)
      call force0_reshaped(q,f0)
   end subroutine

   subroutine force0_reshaped(q,f0)
      use system, only : ns
      real(dp), intent(in) :: q(nf_bath,ns)
      real(dp), intent(out) :: f0(nf_bath,ns)
!
!     Compute state-independent (CPA) force for special case of a Frenkel-exciton model
!
      do n=1,ns
         f0(:,n) = - omega**2 * q(:,n)
      end do
   end subroutine

   subroutine force (q,u,f)
!
!     ------------------------------------------------------------------ 
!     Calculates the MASH forces on the nuclei in adiabatic state u.
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(in) :: u(:)
      real(dp), intent(out) :: f(:)
!
      unused = q(1)
      call force_reshaped(u,f)
   end subroutine

   subroutine force_reshaped(u,f)
      use system, only : ns
      complex(dpc), intent(in) :: u(:)
      real(dp), intent(out) :: f(nf_bath,ns)
!
!     Compute gradient of adiabatic potential a for special case of a Frenkel-exciton model
!
      do n=1,ns
         f(:,n) = -abs(u(n))**2 * rkappa(:)
      end do
   end subroutine

   subroutine sampnuc (p,q)
      use maths, only : gasdev
      use system, only : nf, beta, sampnuc_opt
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
      do i=1,nf
         j = mod(i-1,nf_bath)+1
         if (sampnuc_opt.eq.'cla') then
            ! Classical
            deltap = sqrt(1.0/beta) 
         else if (sampnuc_opt.eq.'wig') then
            ! (Thermal) Wigner 
            deltap = sqrt(omega(j)/(2.0*tanh(beta*omega(j)/2.0))) 
         else 
            stop 'Unknown sampnuc_opt'
         end if
         p(i) = deltap*p(i)
         deltaq = deltap/omega(j)
         q(i) = deltaq*q(i)
      end do
   end subroutine

   subroutine pot (q,v)
      use maths, only : symevp
!
!     ------------------------------------------------------------------ 
!     Diabatic potential energy matrix V(q)
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: v(:,:)
!
      v = Vconst
      call addbath(q,v)
   end subroutine

   subroutine addbath(q, V)
      use system, only : ns
      real(dp), intent(in) :: q(nf_bath,ns)
      complex(dpc), intent(inout) :: V(ns,ns)
!
!     Potential for Frenkel-exciton model. 
!     Reshaping of q is done automatically when called with shape q(nf)
!
      do n=1,ns
         V(n,n) = V(n,n) + sum(rkappa(:)*q(:,n))
      end do
   end subroutine

   subroutine gradient_dia(q,g)
      use system, only : ns, nf
!
!     ------------------------------------------------------------------
!     Diabatic potential gradient (excluding harmonic part)
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: g(:,:,:) ! (nf,ns,ns)
      unused = q(1)
      g = 0.d0
      do n=1,ns
         do i=1,nf
            j = mod(i-1,nf_bath)+1
            g(j,n,n) = rkappa(j)
         end do
      end do
      end subroutine

   subroutine nac (q,u,v,ia,ib,d)
      use system, only : ns,nf
!
!     ------------------------------------------------------------------ 
!     Calculate element ia,ib of the nonadiabatic coupling vector
!     from the transformation matrix (u) and adiabatic energies (v).
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:), v(:)
      complex(dpc), intent(in) :: u(:,:)
      integer, intent(in) :: ia,ib
      complex(dpc), intent(out) :: d(:)
      complex(dpc), allocatable :: d_reshaped(:,:)
      unused = q(1)
      allocate(d_reshaped(nf_bath,ns))
      call nac_reshaped(u,v,ia,ib,d_reshaped)
      d = reshape(d_reshaped,(/nf/))
      deallocate(d_reshaped)
   end subroutine

   subroutine nac_reshaped(u,v,ia,ib,d)
      use system, only : ns
      complex(dpc), intent(in) :: u(:,:)
      real(dp), intent(in) :: v(:)
      integer, intent(in) :: ia,ib
      complex(dpc), intent(out) :: d(nf_bath,ns)
!
!     Compute element ia,ib of the nonadiabatic coupling vector for special case of a Frenkel-exciton model
!     NOTE: Written for a real u matrix
!
      do n=1,ns
         d(:,n) = u(n,ia)*u(n,ib)*rkappa(:)/(v(ib)-v(ia))
      end do
   end subroutine

   function h0 (p,q)
      use system, only : nf
!
!     ------------------------------------------------------------------ 
!     Calculates state-independent (bath) energy.
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: p(:),q(:)
      h0 = 0.d0
      do i=1,nf
         j = mod(i-1,nf_bath)+1
         h0 = h0 + (p(i)**2+(omega(j)*q(i))**2)/2
      end do
   end function


end module
