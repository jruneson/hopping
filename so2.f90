module so2
   use types
   implicit real(dp) (a-h,o-z)

   allocatable :: w(:) ! Frequencies
   allocatable :: Vconst(:,:)
   allocatable :: Vlin(:,:,:)
   allocatable :: Vsoc(:,:) ! (Real in this model)

   allocatable :: Vconst_SF(:,:), Vlin_SF(:,:,:)


contains
   subroutine initso2(model_number)
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc
      use system, only : sys_h0=>h0, sys_pot=>pot, sys_grad=>grad
      use system, only : ns, nf, nmult, rmass, gam, sampnuc_opt
!
!     ------------------------------------------------------------------
!     Initialise SO2 module.
!     Input parameters from Plasser et al, PCCP 2019
!     at CIS level (model_number=1) or at CISD level (model_number=2)
!     ------------------------------------------------------------------
!
      allocatable :: rkappa(:,:), rlamda(:,:,:), eps(:)
      if (nf.ne.3) stop 'SO2 requires nf=3'
      if (ns.ne.13) stop 'SO2 requires ns=13'
      nsing = 4
      ntrip = 3
      nmult = 7
      if (nsing+ntrip.ne.nmult) stop 'nsing and ntrip must add to nmult'
      allocate(w(nf),Vconst(ns,ns),Vlin(nf,ns,ns),Vsoc(ns,ns))
      allocate(Vconst_SF(nmult,nmult),Vlin_SF(nf,nmult,nmult))
      allocate(rkappa(nf,ns),rlamda(nf,ns,ns),eps(nmult))

      ! Constants
      ev_to_au = 1./27.2113961
      cmm1 = 4.556335d-6

      ! Initialize
      rkappa = 0.0
      rlamda = 0.0
      Vconst = 0.0
      Vlin = 0.0
      Vconst_SF = 0.0
      Vlin_SF = 0.0
      sampnuc_opt = 'wig'

      ! Frequencies
      w = (/ 0.0023635993,  0.0053089059,  0.0064027455 /)
   
      ! Masses
      rmass = 1.d0/w

      ! Width parameter
      gam = 1.d0

      select case (model_number)
      case(1)
         ! CIS
         ! Vertical transition energies (SOC free states)
                               ! Term = SHARC notation
         eps(1) = 0.0000000000 ! 1A1  = 1  1
         eps(2) = 0.1640045037 ! 1B1  = 1  2
         eps(3) = 0.1781507393 ! 1A2  = 1  3
         eps(4) = 0.2500917150 ! 1B2  = 1  4
         eps(5) = 0.1341247878 ! 3B1  = 3  1
         eps(6) = 0.1645714941 ! 3B2  = 3  2
         eps(7) = 0.1700512159 ! 3A2  = 3  3
   
         ! On-diagonal coupling constants (mode,state)
         rkappa(1,2) =  1.19647e-03
         rkappa(2,2) =  1.21848e-02
         rkappa(1,3) = -7.98767e-03
         rkappa(2,3) =  1.85494e-02
         rkappa(1,4) = -4.45712e-03
         rkappa(2,4) =  1.77413e-02
         rkappa(1,5) =  2.87702e-03
         rkappa(2,5) =  8.19136e-03
         rkappa(1,6) = -5.03345e-03
         rkappa(2,6) =  1.83566e-02
         rkappa(1,7) = -8.05428e-03
         rkappa(2,7) =  1.86062e-02

         ! Off-diagonal coupling constantds (mode, state, state)
         rlamda(3,1,4) = -1.83185e-02
         rlamda(3,2,3) =  7.32022e-03
         rlamda(3,5,7) =  5.71746e-03

      case(2)
         ! CISD
         ! Vertical transition energies (SOC free states)
                               ! Term = SHARC notation
         eps(1) = 0.0000000000 ! 1A1  = 1  1
         eps(2) = 0.1553277839 ! 1B1  = 1  2
         eps(3) = 0.1688652931 ! 1A2  = 1  3
         eps(4) = 0.3088375445 ! 1B2  = 1  4
         eps(5) = 0.1231069775 ! 3B1  = 3  1
         eps(6) = 0.1547471701 ! 3B2  = 3  2
         eps(7) = 0.1600894729 ! 3A2  = 3  3

         ! On-diagonal coupling constants (mode,state)
         rkappa(1,1) = -6.48289e-04
         rkappa(2,1) = -9.15365e-04
         rkappa(1,2) =  1.06541e-03
         rkappa(2,2) =  1.08185e-02
         rkappa(1,3) = -7.92795e-03
         rkappa(2,3) =  1.63364e-02
         rkappa(1,4) = -4.80026e-03
         rkappa(2,4) =  1.70393e-02
         rkappa(1,5) =  2.60893e-03
         rkappa(2,5) =  6.95032e-03
         rkappa(1,6) = -5.24793e-03
         rkappa(2,6) =  1.71110e-02
         rkappa(1,7) = -7.98845e-03
         rkappa(2,7) =  1.62033e-02
   
         ! Off-diagonal coupling constantds (mode, state, state)
         rlamda(3,1,4) =  1.65920e-02
         rlamda(3,2,3) = -5.55909e-03
         rlamda(3,5,7) = -3.66103e-03
      case default
         stop 'so2.f90: Unknown model_number'
      end select

      ! Fill Vconst_SF and Vlin_SF
      nmult = nsing + ntrip
      do i=1,nmult
         Vconst_SF(i,i) = eps(i)
         Vlin_SF(:,i,i) = rkappa(:,i)
         do i2=i+1,nmult
            Vlin_SF(:,i,i2) = rlamda(:,i,i2)
            Vlin_SF(:,i2,i) = rlamda(:,i,i2)
         end do
      end do

      ! Fill Vconst and Vlin
      do i=1,nsing ! Singlets
         Vconst(i,i) = eps(i)
         Vlin(:,i,i) = rkappa(:,i)
         do i2=i+1,nsing ! S-S coupling
            Vlin(:,i,i2) = rlamda(:,i,i2)
         end do
      end do
      do i=1,ntrip ! Triplets
         do j=1,3 ! Multiplicities
            k = nsing + (j-1)*ntrip + i
            Vconst(k,k) = eps(nsing+i)
            Vlin(:,k,k) = rkappa(:,nsing+i)
            do i2=1,nsing
               Vlin(:,i2,k) = rlamda(:,i2,nsing+i)
               Vlin(:,k,i2) = rlamda(:,i2,nsing+i)
            end do
            do i2=i+1,ntrip
               do j2=1,3
                  k2 = nsing + (j2-1)*ntrip + i2
                  Vlin(:,k,k2) = rlamda(:,nsing+i,nsing+i2)
                  Vlin(:,k2,k) = rlamda(:,nsing+i,nsing+i2)
               end do
            end do
         end do
      end do

      open(10,file='soc_re.in')
      do i=1,ns
         read(10,*) Vsoc(i,:) ! Read SOC
         ! read(10,*) Vsoc(:,i) ! Read SOC
      end do
      close(10)

      sys_force0=>force0
      sys_force=>force
      sys_sampnuc=>sampnuc
      sys_pot=>potenl_dia
      sys_grad=>gradient_dia
      sys_h0=>h0

      deallocate(rkappa,rlamda,eps)
   end subroutine

   subroutine force0 (q,f0)
!
!     ------------------------------------------------------------------
!     Calculates the state-independent (CPA) force f0 = - w * q
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: f0(:)
      f0 = - w * q
   end subroutine

   subroutine force (q,u,f)
      use system, only : nf, ns
!
!     ------------------------------------------------------------------
!     Calculates the state-dependent force on the nuclei in adiabatic state u.
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(in) :: u(:)
      real(dp), intent(out) :: f(:)
      complex(dpc), allocatable :: g(:,:,:)
      allocate(g(nf,ns,ns))
      call gradient_dia(q,g)
      do i=1,nf
         f(i) = -real(sum(conjg(u)*matmul(g(i,:,:),u)))
      end do
      deallocate(g)
   end subroutine

   subroutine potenl_dia(q,u)
      use system, only : ns
!
!     ------------------------------------------------------------------
!     Diabatic potential
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: u(:,:)
      u = Vconst + Vsoc
      do n=1,ns
         do m=1,ns
            u(n,m) = u(n,m) + sum(Vlin(:,n,m)*q)
         end do
      end do
   end subroutine

   subroutine potenl_sf(q,u)
      use system, only : nmult
!
!     ------------------------------------------------------------------
!     Spin-free potential
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: u(:,:)
      u = Vconst_SF
      do n=1,nmult
         do m=1,nmult
            u(n,m) = u(n,m) + sum(Vlin_SF(:,n,m)*q)
         end do
      end do
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
      g = Vlin
   end subroutine

   subroutine gradient_SF(q,g)
!
!     ------------------------------------------------------------------
!     Spin-free potential gradient (excluding harmonic part)
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: g(:,:,:)
      unused = q(1)
      g = Vlin_SF
   end subroutine

   subroutine sampnuc (p,q)
      use maths, only : gasdev
      use system, only : nf, sampnuc_opt
!
!     ------------------------------------------------------------------
!     Samples the nuclear momenta (p) and coordinates (q)
!     ------------------------------------------------------------------
!
      real(dp), intent(out) :: p(:),q(:)
      call gasdev (p,nf)
      call gasdev (q,nf)
      do i=1,nf
         if (sampnuc_opt.eq.'cla') then
            ! Wigner
            deltap = sqrt(0.0)
            deltaq = sqrt(0.0)
            p(i) = deltap*p(i)
            q(i) = deltaq*q(i)
         else if (sampnuc_opt.eq.'wig') then
            ! Wigner
            deltap = sqrt(0.5)
            deltaq = sqrt(0.5)
            p(i) = deltap*p(i)
            q(i) = deltaq*q(i)
         else
            stop 'Unknown sampnuc_opt'
         end if
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
         h0 = h0 + w(i)*(p(i)**2+q(i)**2)/2.d0
      end do
   end function


end module

