module ful
   use types
   implicit real(dp) (a-h,o-z)

   allocatable :: w(:) ! Frequencies
   allocatable :: Vconst(:,:)
   allocatable :: Vlin(:,:,:)

contains
   subroutine initful()
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc
      use system, only : sys_h0=>h0, sys_pot=>pot, sys_grad=>grad
      use system, only : ns, nf, rmass, gam, sampnuc_opt
      !
      !     ------------------------------------------------------------------
      !     Initialise fulvene module.
      !     Input parameters from Gomez et al PCCP 26, 1829 (2024)
      !     ------------------------------------------------------------------
      !
      allocatable :: rkappa(:,:), rlamda(:,:,:), eps(:)
      if (nf.ne.30) stop 'fulvene requires nf=30'
      if (ns.ne.2) stop 'fulvene requires ns=2'
      allocate(w(nf),Vconst(ns,ns),Vlin(nf,ns,ns))
      allocate(rkappa(ns,nf),rlamda(ns,ns,nf),eps(ns))

      ! Constants
      ev = 1./27.2113961
      cmm1 = 4.556335d-6

      ! Initialize
      rkappa = 0.0
      rlamda = 0.0
      Vconst = 0.0
      Vlin = 0.0
      sampnuc_opt = 'wig'

      ! Frequencies
      w(1 ) =  0.02610
      w(2 ) =  0.04530
      w(3 ) =  0.06174
      w(4 ) =  0.07760
      w(5 ) =  0.08560
      w(6 ) =  0.08750
      w(7 ) =  0.09489
      w(8 ) =  0.09766
      w(9 ) =  0.10603
      w(10) =  0.11074
      w(11) =  0.11154
      w(12) =  0.11242
      w(13) =  0.11739
      w(14) =  0.12811
      w(15) =  0.12909
      w(16) =  0.14578
      w(17) =  0.14629
      w(18) =  0.16851
      w(19) =  0.18011
      w(20) =  0.18285
      w(21) =  0.19361
      w(22) =  0.20310
      w(23) =  0.20964
      w(24) =  0.22047
      w(25) =  0.41409
      w(26) =  0.42437
      w(27) =  0.42759
      w(28) =  0.42790
      w(29) =  0.43081
      w(30) =  0.43248

      w = w * ev

      ! Masses
      rmass = 1.d0/w

      ! Wavepacket widths
      gam = 1.d0

      ! Vertical transition energies 
      eps(1) = 0.00000 ! X
      eps(2) = 4.16142 ! A

      ! On-diagonal coupling constants (state,mode)
      rkappa(1,6 ) = -0.00799
      rkappa(1,13) =  0.00608
      rkappa(1,15) =  0.00487
      rkappa(1,17) = -0.01463
      rkappa(1,20) = -0.01265
      rkappa(1,21) = -0.00096
      rkappa(1,22) = -0.00421
      rkappa(1,24) = -0.00090
      rkappa(1,25) =  0.00464
      rkappa(1,28) =  0.01287
      rkappa(1,30) =  0.05083
      rkappa(2,6 ) = -0.16463
      rkappa(2,13) =  0.16902
      rkappa(2,15) =  0.20870
      rkappa(2,17) = -0.19129
      rkappa(2,20) =  0.12386
      rkappa(2,21) = -0.06968
      rkappa(2,22) = -0.38302
      rkappa(2,24) =  0.49293
      rkappa(2,25) = -0.05766
      rkappa(2,28) =  0.01640
      rkappa(2,30) =  0.03264
      

      ! Off-diagonal coupling constantds (state, state, mode)
      rlamda(1,2,2 ) = -0.02235
      rlamda(1,2,9 ) = -0.14659
      rlamda(1,2,14) =  0.04999
      rlamda(1,2,16) =  0.02678
      rlamda(1,2,18) = -0.07546
      rlamda(1,2,19) = -0.14017
      rlamda(1,2,23) =  0.04228
      rlamda(1,2,26) =  0.00122
      rlamda(1,2,27) =  0.00416
      rlamda(1,2,29) = -0.00788

      ! Fill Vconst and Vlin
      do i=1,ns 
         Vconst(i,i) = eps(i)
         Vlin(:,i,i) = rkappa(i,:)
         do j=i+1,ns
            Vlin(:,i,j) = rlamda(i,j,:)
            Vlin(:,j,i) = rlamda(i,j,:) ! Symmetrize Vlin
         end do
      end do

      Vconst = Vconst * ev
      Vlin = Vlin * ev

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
      u = Vconst
      do n=1,ns
         do m=1,ns
            u(n,m) = u(n,m) + sum(Vlin(:,n,m)*q)
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


