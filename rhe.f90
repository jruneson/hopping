module rhe
   use types
   implicit real(dp) (a-h,o-z)

   ! Rhenium(I) carbonyl alpha-diimine complex

   allocatable :: w(:) ! Frequencies
   allocatable :: Vconst(:,:)
   allocatable :: Vlin(:,:,:)
   complex(dpc), allocatable :: Vsoc(:,:)


contains

   subroutine initrhe()
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc
      use system, only : sys_h0=>h0, sys_pot=>pot, sys_grad=>grad
      use system, only : ns, nf, rmass, gam, sampnuc_opt, potmat
!
!  ------------------------------------------------------------------
!  Initialise Rhenium(I) carbonyl alpha-diimine complex
!  Input data from Fumanal et al. JCTC 13, 1293 (2017)
!  ------------------------------------------------------------------
!
      complex(dpc) :: u
      allocatable :: eps(:), rkappa(:,:), rlamda(:,:,:), soc_re(:,:), soc_im(:,:)
      allocatable :: q(:), u(:,:), v(:)
      if (nf.ne.15) stop 'Rhenium(I) carbonyl alpha-diimine complex requires nf=15'
      if (ns.ne.14) stop 'Rhenium(I) carbonyl alpha-diimine complex requires ns=14'
      allocate(w(nf),Vconst(ns,ns),Vlin(nf,ns,ns),Vsoc(ns,ns),eps(ns),rkappa(6,nf),rlamda(6,6,nf),soc_re(ns,ns),soc_im(ns,ns))
      sampnuc_opt = 'wig'

      Vconst = 0.d0
      Vlin = 0.d0
      Vsoc = 0.d0

      ! Frequencies (Table S4)
      ! a' modes
      w(1)  =   93.0
      w(2)  =  235.0
      w(3)  =  439.0
      w(4)  =  498.0
      w(5)  =  552.0
      w(6)  =  637.0
      w(7)  = 1174.0
      w(8)  = 1336.0
      w(9)  = 1444.0
      w(10) = 1554.0
      w(11) = 1623.0
      w(12) = 1660.0
      ! a'' modes
      w(13) =   90.0
      w(14) =  475.0
      w(15) =  626.0

      ! Convert from cmm1 to au
      cmm1 = 4.556335d-6
      w = w * cmm1

      ! Masses
      rmass = 1.d0/w

      ! Width parameters
      gam = 1.d0

      ! Vertical transition energy (without SOC) (Table S1)
      ! Order S1, S2, T1, T2, T3, T4
      eps(1) = 0.1145475957
      eps(2) = 0.1248374022
      eps(3) = 0.1094394418
      eps(4) = 0.1129673755
      eps(5) = 0.119214758
      eps(6) = 0.1256091377

      ! Fill Vconst with vertical transition energies.
      ! Order of states is S1, S2, [T1,T2,T3,T4]-, [T1,T2,T3,T4]_0, [T1,T2,T3,T4]+
      Vconst(1,1) = eps(1)
      Vconst(2,2) = eps(2)
      do i=1,3
         ii = 3 + (i-1)*4
         Vconst(ii,ii) = eps(3)
         Vconst(ii+1,ii+1) = eps(4)
         Vconst(ii+2,ii+2) = eps(5)
         Vconst(ii+3,ii+3) = eps(6)
      end do

      ! Linear spin-independent couplings (state,mode) (Table S4)
      rkappa = 0.d0
      rkappa(1,1) =  0.00083788
      rkappa(2,1) = -0.00069824
      rkappa(3,1) =  0.00063944
      rkappa(4,1) = -0.00061004
      rkappa(5,1) = -0.00015067
      rkappa(6,1) = -0.00010657
      rkappa(1,2) = -0.00210574
      rkappa(2,2) = -0.00142587
      rkappa(3,2) = -0.00164269
      rkappa(4,2) = -0.00106206
      rkappa(5,2) = -0.00084156
      rkappa(6,2) = -0.00065781
      rkappa(1,3) =  0.00034912
      rkappa(2,3) = -0.00008452
      rkappa(3,3) =  0.00045202
      rkappa(4,3) = -0.00022785
      rkappa(5,3) = -0.00223436
      rkappa(6,3) =  0.00043732
      rkappa(1,4) = -0.00327804
      rkappa(2,4) = -0.00288850
      rkappa(3,4) = -0.00285542
      rkappa(4,4) = -0.00206531
      rkappa(5,4) = -0.00116495
      rkappa(6,4) = -0.00324496
      rkappa(1,5) =  0.00054756
      rkappa(2,5) = -0.00006982
      rkappa(3,5) =  0.00040792
      rkappa(4,5) = -0.00089668
      rkappa(5,5) = -0.00041159
      rkappa(6,5) =  0.00018007
      rkappa(1,6) = -0.00109145
      rkappa(2,6) =  0.00137810
      rkappa(3,6) = -0.00087831
      rkappa(4,6) =  0.00115760
      rkappa(5,6) =  0.00032707
      rkappa(6,6) = -0.00006982
      rkappa(1,7) =  0.00265698
      rkappa(2,7) =  0.00139280
      rkappa(3,7) =  0.00217923
      rkappa(4,7) =  0.00172354
      rkappa(5,7) = -0.00239605
      rkappa(6,7) =  0.00198079
      rkappa(1,8) =  0.00391380
      rkappa(2,8) =  0.00383663
      rkappa(3,8) =  0.00411225
      rkappa(4,8) =  0.00557487
      rkappa(5,8) =  0.00554915
      rkappa(6,8) =  0.00410122
      rkappa(1,9) =  0.00303917
      rkappa(2,9) =  0.00338461
      rkappa(3,9) =  0.00337726
      rkappa(4,9) =  0.00426659
      rkappa(5,9) =  0.00663693
      rkappa(6,9) =  0.00297302
      rkappa(1,10) = -0.00675452
      rkappa(2,10) = -0.00460836
      rkappa(3,10) = -0.00634293
      rkappa(4,10) = -0.00619961
      rkappa(5,10) =  0.00198079
      rkappa(6,10) = -0.00554915
      rkappa(1,11) =  0.00475536
      rkappa(2,11) =  0.00247323
      rkappa(3,11) =  0.00418575
      rkappa(4,11) =  0.00333316
      rkappa(5,11) = -0.00296934
      rkappa(6,11) =  0.00336991
      rkappa(1,12) = -0.00183379
      rkappa(2,12) =  0.00005880
      rkappa(3,12) = -0.00040424
      rkappa(4,12) =  0.00210206
      rkappa(5,12) =  0.00869856
      rkappa(6,12) = -0.00082318

      ! ev_to_au = 1./27.2113961
      ! print*, 'test'
      ! do i = 1, 12
      !    print '(6F10.4)', (rkappa(j,i)/ev_to_au, j = 1, 6)
      ! end do

      ! Fill Vlin with rkappa
      Vlin(:,1,1) = rkappa(1,:)
      Vlin(:,2,2) = rkappa(2,:)
      do i=1,3
         ii = 3 + (i-1)*4
         Vlin(:,ii,ii) = rkappa(3,:) ! T1
         Vlin(:,ii+1,ii+1) = rkappa(4,:) ! T2
         Vlin(:,ii+2,ii+2) = rkappa(5,:) ! T3
         Vlin(:,ii+3,ii+3) = rkappa(6,:) ! T4
      end do

      ! Off-diagonal linear spin-independent couplings (state,state',mode) (Table S4)
      rlamda = 0.d0
      rlamda(4,6,1) = 0.00027929
      rlamda(3,5,2) = 0.00034912
      rlamda(3,5,3) = 0.00069456
      rlamda(3,5,4) = 0.00050714
      rlamda(4,6,4) = 0.00083053
      rlamda(3,5,5) = 0.00067986
      rlamda(4,6,5) = 0.00075336
      rlamda(3,5,6) = 0.00069456
      rlamda(3,5,7) = 0.00152142
      rlamda(3,5,9) = 0.00137075
      rlamda(4,6,9) = 0.00069456
      rlamda(3,5,10) = 0.00169782
      rlamda(3,5,11) = 0.00176764
      rlamda(3,5,12) = 0.00301712
      rlamda(1,2,13) = 0.00058064
      rlamda(3,4,13) = 0.00022417
      rlamda(5,6,14) = 0.00063576
      rlamda(1,2,14) = 0.00071661
      rlamda(3,4,14) = 0.00138912
      rlamda(5,6,14) = 0.00027929
      rlamda(1,2,15) = 0.00144425
      rlamda(3,4,15) = 0.00083053
      rlamda(5,6,15) = 0.00128255

      ! Fill Vlin with rlamda
      Vlin(:,1,2) = rlamda(1,2,:)
      do i=1,3
         ii = 3 + (i-1)*4
         Vlin(:,ii,ii+1) = rlamda(3,4,:) ! T1-T2
         Vlin(:,ii,ii+2) = rlamda(3,5,:) ! T1-T3
         Vlin(:,ii,ii+3) = rlamda(3,6,:) ! T1-T4
         Vlin(:,ii+1,ii+2) = rlamda(4,5,:) ! T2-T3
         Vlin(:,ii+1,ii+3) = rlamda(4,6,:) ! T2-T4
         Vlin(:,ii+2,ii+3) = rlamda(5,6,:) ! T3-T4
      end do

      ! Symmetrize Vlin
      do i=1,ns
         do j=i+1,ns
            Vlin(:,j,i) = Vlin(:,i,j)
         end do
      end do

      ! Spin-orbit coupling (state,state'), from file
      open(10,file='soc_re.in')
      open(11,file='soc_im.in')
      do i=1,ns
         read(10,*) soc_re(i,:) ! Read real part
         read(11,*) soc_im(i,:) ! Read imaginary part
      end do
      close(10)
      close(11)
      Vsoc = dcmplx(soc_re,soc_im)
      
      sys_force0=>force0
      sys_force=>force
      sys_sampnuc=>sampnuc
      sys_pot=>potenl_dia
      sys_grad=>gradient_dia
      sys_h0=>h0

      allocate(q(nf),u(ns,ns),v(ns))
      q = 0.0
      call potmat(q,u,v)
      ev_to_au = 1./27.2113961
      ! print '(F10.4)', v/ev_to_au
      ! do i=1,ns
      !    print '(14F10.4)', (abs(u(j,i))**2, j=1,ns)
      ! end do
      deallocate(q,u,v)

      deallocate(eps,rkappa,rlamda,soc_re,soc_im)

      ! stop 'check energies'
   end subroutine

   subroutine force0 (q,f0)
!
!     ------------------------------------------------------------------
!     Calculates the state-independent (CPA) force f0 = - m*w**2 * q
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
!
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
!
!     Sample p and q from the Boltzmann distribution
!
      call gasdev (p,nf)
      call gasdev (q,nf)
      do i=1,nf
         if (sampnuc_opt.eq.'cla') then
            ! Classical
            deltap = 0.d0
            deltaq = 0.d0
         else if (sampnuc_opt.eq.'wig') then
            ! Wigner
            deltap = sqrt(0.5)
            deltaq = sqrt(0.5)
         else
            stop 'Unknown sampnuc_opt'
         end if
         p(i) = deltap*p(i)
         q(i) = deltaq*q(i)
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
