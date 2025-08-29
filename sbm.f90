module sbm
   use types
   implicit real(dp) (a-h,o-z)

   allocatable :: Vconst(:,:) ! Constant part of potential matrix
   allocatable :: omega(:)    ! Phonon frequencies
   allocatable :: rkappa(:)    ! On-diagonal coupling strengths


contains

   subroutine initsbm(eps,Delta,rlamda,omegac)
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc, sys_pot=>pot
      use system, only : sys_nac=>nac, sys_h0=>h0, sys_grad=>grad
      use system, only : ns, nf, gam, rmass
!
!     ------------------------------------------------------------------ 
!     Initialise Spin-boson module
!     ------------------------------------------------------------------ 
!
      if (ns .ne. 2) stop 'Spin-Boson requires 2 states'
      allocate(Vconst(ns,ns),omega(nf),rkappa(nf))
      Vconst = 0.0
      Vconst(1,1) = eps
      Vconst(2,2) = -eps
      Vconst(1,2) = Delta
      Vconst(2,1) = Delta
      ! Vconst = Vconst*cmm1
!
!     Discretize bath
!
      ! rlamda = rlamda*cmm1
      ! omegac = omegac*cmm1
      fac = sqrt(0.5*rlamda/nf)
      rkappa = 0.0
      do i = 1,nf
         omega(i) = omegac*tan(0.5*pi*(i-0.5)/nf)
         rkappa(i) = fac*omega(i)
      end do
!
!     Mass and width parameter
!
      rmass = 1.d0
      ! print*, 'beta', beta
      gam = rmass * omega !* tanh(beta*omega/2.d0)
!
      sys_force0=>force0
      sys_force=>force
      sys_grad=>gradient_dia
      sys_sampnuc=>sampnuc
      sys_pot=>pot
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
      f0 = - omega**2 * q
   end subroutine

   subroutine force (q,u,f)
!
!     ------------------------------------------------------------------ 
!     Calculates the forces on the nuclei in complex state u.
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(in) :: u(:)
      real(dp), intent(out) :: f(:)
!
      unused = q(1)
      f = -abs(u(1))**2 * rkappa(:)
      f = f + abs(u(2))**2 * rkappa(:)
   end subroutine

   subroutine gradient_dia(q,g)
!
!     ------------------------------------------------------------------
!     Diabatic potential gradient (excluding harmonic part)
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: g(:,:,:) ! (nf,ns,ns)
      unused = q(1)
      g = 0.d0
      g(:,1,1) = rkappa
      g(:,2,2) = -rkappa
   end subroutine

   subroutine sampnuc (p,q)
      use maths, only : gasdev
      use system, only : nf, beta, sampnuc_opt, rmass
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
         if (sampnuc_opt.eq.'cla') then
            ! Classical
            deltap = sqrt(rmass(i)/beta)
            deltaq = deltap/(rmass(i)*omega(i)) 
         else if (sampnuc_opt.eq.'wig') then
            ! (Thermal) Wigner 
            deltap = sqrt(rmass(i)*omega(i)/(2.0*tanh(beta*omega(i)/2.0))) 
            deltaq = deltap/(rmass(i)*omega(i))
         else 
            print*, sampnuc_opt
            stop 'Unknown sampnuc_opt'
         end if
         p(i) = deltap*p(i)
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
      v(1,1) = v(1,1) + sum(rkappa(:)*q)
      v(2,2) = v(2,2) - sum(rkappa(:)*q)
   end subroutine

   subroutine nac (q,u,v,ia,ib,d)
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
!
!     Returns components d(i,b,a) for a given adiabatic state a
!
      unused = q(1)
      d(:) = (u(1,ia)*u(1,ib) - u(2,ia)*u(2,ib)) * rkappa(:) / (v(ib)-v(ia))
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
         h0 = h0 + (p(i)**2+(omega(i)*q(i))**2)/2
      end do
   end function


end module
