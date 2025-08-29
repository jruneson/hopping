module dma
   use types
   implicit real(dp) (a-h,o-z)

   allocatable :: w(:) ! Frequencies
   allocatable :: Vconst(:,:)
   allocatable :: Vlin(:,:,:)

contains
   subroutine initdma()
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc
      use system, only : sys_h0=>h0, sys_pot=>pot, sys_grad=>grad
      use system, only : ns, nf, rmass, gam
      !
      !     ------------------------------------------------------------------
      !     Initialise DMABN module.
      !     Input parameters from Gomez et al PCCP 26, 1829 (2024)
      !     ------------------------------------------------------------------
      !
      allocatable :: rkappa(:,:), rlamda(:,:,:), eps(:)
      if (nf.ne.57) stop 'DMABN requires nf=57'
      if (ns.ne.3) stop 'DMABN requires ns=3'
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

      ! Frequencies
      w = (/ 0.00622, &
             0.00909, &
             0.01033, &
             0.01673, &
             0.02188, &
             0.02317, &
             0.03186, &
             0.03575, &
             0.04213, &
             0.05324, &
             0.05919, &
             0.06109, &
             0.06234, &
             0.07081, &
             0.07103, &
             0.08237, &
             0.08386, &
             0.09314, &
             0.10163, &
             0.10377, &
             0.10575, &
             0.12272, &
             0.12294, &
             0.12408, &
             0.12707, &
             0.13479, &
             0.14102, &
             0.14110, &
             0.14298, &
             0.14786, &
             0.14935, &
             0.15535, &
             0.16070, &
             0.16529, &
             0.17258, &
             0.17568, &
             0.17782, &
             0.18109, &
             0.18219, &
             0.18324, &
             0.18354, &
             0.18661, &
             0.18749, &
             0.19708, &
             0.20350, &
             0.21190, &
             0.29621, &
             0.37412, &
             0.37509, &
             0.38359, &
             0.38373, &
             0.39291, &
             0.39407, &
             0.39931, &
             0.39943, &
             0.40299, &
             0.40309 /)
      w = w * ev

      ! Masses
      rmass = 1.d0/w

      ! Widths
      gam = 1.d0

      ! Vertical transition energies 
      eps(1) = 0.00000 ! X
      eps(2) = 4.93397 ! A
      eps(3) = 5.31142 ! B

      ! On-diagonal coupling constants (state,mode)
      rkappa(1,3 ) = -0.00036
      rkappa(1,4 ) = -0.00082
      rkappa(1,5 ) = -0.00045
      rkappa(1,7 ) =  0.00032
      rkappa(1,8 ) =  0.00037
      rkappa(1,9 ) = -0.00043
      rkappa(1,11) = -0.00049
      rkappa(1,12) = -0.00071
      rkappa(1,13) = -0.00038
      rkappa(1,15) =  0.00018
      rkappa(1,17) = -0.00019
      rkappa(1,19) = -0.00044
      rkappa(1,22) =  0.00033
      rkappa(1,23) = -0.00038
      rkappa(1,25) = -0.00026
      rkappa(1,26) = -0.00015
      rkappa(1,29) =  0.00014
      rkappa(1,31) = -0.00033
      rkappa(1,32) =  0.00046
      rkappa(1,36) = -0.00066
      rkappa(1,37) = -0.00025
      rkappa(1,41) =  0.00036
      rkappa(1,43) =  0.00014
      rkappa(1,46) =  0.00012
      rkappa(1,48) = -0.00014
      rkappa(1,49) =  0.00011
      rkappa(1,53) = -0.00038
      rkappa(1,54) = -0.00018
      rkappa(1,56) = -0.00028
      rkappa(1,57) = -0.00027
      rkappa(2,1 ) =  0.01981
      rkappa(2,2 ) =  0.00114
      rkappa(2,3 ) =  0.00679
      rkappa(2,4 ) = -0.00129
      rkappa(2,5 ) = -0.01080
      rkappa(2,6 ) = -0.00051
      rkappa(2,7 ) =  0.00074
      rkappa(2,8 ) =  0.01364
      rkappa(2,9 ) =  0.01326
      rkappa(2,10) =  0.00025
      rkappa(2,11) = -0.00067
      rkappa(2,12) = -0.02516
      rkappa(2,13) = -0.01572
      rkappa(2,14) =  0.00029
      rkappa(2,15) =  0.00120
      rkappa(2,16) = -0.00017
      rkappa(2,17) =  0.00630
      rkappa(2,18) =  0.00095
      rkappa(2,19) = -0.09430
      rkappa(2,21) = -0.00292
      rkappa(2,22) = -0.03331
      rkappa(2,23) =  0.04050
      rkappa(2,25) = -0.00984
      rkappa(2,26) = -0.00036
      rkappa(2,27) = -0.01563
      rkappa(2,28) =  0.00299
      rkappa(2,29) = -0.00022
      rkappa(2,30) = -0.01587
      rkappa(2,31) =  0.00569
      rkappa(2,32) =  0.07393
      rkappa(2,33) =  0.00135
      rkappa(2,34) = -0.00099
      rkappa(2,35) =  0.00021
      rkappa(2,36) =  0.07246
      rkappa(2,37) =  0.00079
      rkappa(2,38) = -0.00031
      rkappa(2,39) =  0.00674
      rkappa(2,40) =  0.01043
      rkappa(2,41) = -0.00021
      rkappa(2,42) = -0.00014
      rkappa(2,43) =  0.02299
      rkappa(2,44) =  0.00436
      rkappa(2,45) = -0.00038
      rkappa(2,46) =  0.10678
      rkappa(2,47) = -0.02123
      rkappa(2,48) = -0.00070
      rkappa(2,49) = -0.01697
      rkappa(2,50) =  0.00086
      rkappa(2,51) =  0.00237
      rkappa(2,52) = -0.00024
      rkappa(2,53) = -0.00153
      rkappa(2,54) = -0.02192
      rkappa(2,55) =  0.00756
      rkappa(2,56) = -0.00137
      rkappa(2,57) =  0.03034
      rkappa(3,1 ) =  0.00797
      rkappa(3,2 ) = -0.00096
      rkappa(3,3 ) = -0.00024
      rkappa(3,4 ) = -0.00094
      rkappa(3,5 ) = -0.00628
      rkappa(3,6 ) = -0.00062
      rkappa(3,7 ) =  0.00110
      rkappa(3,8 ) =  0.00920
      rkappa(3,9 ) = -0.00695
      rkappa(3,10) = -0.00019
      rkappa(3,11) = -0.00010
      rkappa(3,12) = -0.01937
      rkappa(3,13) = -0.01093
      rkappa(3,14) = -0.00014
      rkappa(3,15) = -0.00166
      rkappa(3,17) = -0.00475
      rkappa(3,18) = -0.00549
      rkappa(3,19) = -0.07216
      rkappa(3,20) =  0.00024
      rkappa(3,21) =  0.00189
      rkappa(3,22) = -0.03457
      rkappa(3,23) =  0.04590
      rkappa(3,25) = -0.01089
      rkappa(3,26) = -0.00013
      rkappa(3,27) = -0.01695
      rkappa(3,28) =  0.00394
      rkappa(3,29) =  0.00141
      rkappa(3,30) = -0.02115
      rkappa(3,31) = -0.07638
      rkappa(3,32) =  0.04662
      rkappa(3,33) =  0.00065
      rkappa(3,34) =  0.00011
      rkappa(3,36) = -0.00036
      rkappa(3,37) = -0.00018
      rkappa(3,38) = -0.00029
      rkappa(3,39) =  0.00487
      rkappa(3,40) = -0.00584
      rkappa(3,41) =  0.00064
      rkappa(3,43) = -0.01030
      rkappa(3,44) =  0.06696
      rkappa(3,45) =  0.00062
      rkappa(3,46) = -0.12549
      rkappa(3,47) = -0.10953
      rkappa(3,48) =  0.00058
      rkappa(3,49) =  0.02129
      rkappa(3,50) =  0.00104
      rkappa(3,51) =  0.00253
      rkappa(3,52) = -0.00031
      rkappa(3,53) =  0.00905
      rkappa(3,54) = -0.01492
      rkappa(3,55) =  0.00526
      rkappa(3,56) = -0.00034
      rkappa(3,57) =  0.00257
      

      ! Off-diagonal coupling constantds (state, state, mode)
      rlamda(1,2,1 ) =  0.00023
      rlamda(1,2,2 ) =  0.00716
      rlamda(1,2,3 ) = -0.00237
      rlamda(1,2,4 ) =  0.00655
      rlamda(1,2,5 ) = -0.00235
      rlamda(1,2,6 ) =  0.02356
      rlamda(1,2,7 ) = -0.01311
      rlamda(1,2,8 ) = -0.00098
      rlamda(1,2,9 ) =  0.00015
      rlamda(1,2,10) = -0.02099
      rlamda(1,2,11) =  0.04703
      rlamda(1,2,12) = -0.00047
      rlamda(1,2,13) = -0.00013
      rlamda(1,2,14) = -0.05253
      rlamda(1,2,15) = -0.00131
      rlamda(1,2,16) =  0.01382
      rlamda(1,2,17) =  0.00073
      rlamda(1,2,18) = -0.00052
      rlamda(1,2,20) = -0.00319
      rlamda(1,2,23) = -0.00066
      rlamda(1,2,24) =  0.00046
      rlamda(1,2,25) = -0.00050
      rlamda(1,2,26) = -0.03919
      rlamda(1,2,27) =  0.00331
      rlamda(1,2,28) =  0.01165
      rlamda(1,2,29) = -0.07418
      rlamda(1,2,30) = -0.00155
      rlamda(1,2,31) =  0.00043
      rlamda(1,2,32) =  0.00221
      rlamda(1,2,33) = -0.23453
      rlamda(1,2,34) = -0.26477
      rlamda(1,2,35) =  0.33483
      rlamda(1,2,36) = -0.00165
      rlamda(1,2,37) =  0.05187
      rlamda(1,2,38) =  0.02173
      rlamda(1,2,39) = -0.00409
      rlamda(1,2,40) =  0.01359
      rlamda(1,2,41) =  0.13492
      rlamda(1,2,42) =  0.06252
      rlamda(1,2,43) = -0.00181
      rlamda(1,2,44) =  0.00021
      rlamda(1,2,45) =  0.03665
      rlamda(1,2,46) = -0.00047
      rlamda(1,2,47) =  0.00018
      rlamda(1,2,48) = -0.01867
      rlamda(1,2,49) =  0.00063
      rlamda(1,2,50) =  0.00470
      rlamda(1,2,51) = -0.00202
      rlamda(1,2,52) = -0.00113
      rlamda(1,2,54) = -0.00052
      rlamda(1,2,55) = -0.00103
      rlamda(1,2,56) =  0.00272
      rlamda(1,2,57) =  0.00041
      rlamda(1,3,1 ) = -0.06569
      rlamda(1,3,3 ) = -0.02458
      rlamda(1,3,5 ) =  0.02188
      rlamda(1,3,6 ) =  0.00162
      rlamda(1,3,7 ) = -0.00029
      rlamda(1,3,8 ) = -0.02782
      rlamda(1,3,9 ) = -0.02332
      rlamda(1,3,10) =  0.00032
      rlamda(1,3,11) = -0.00013
      rlamda(1,3,12) = -0.00192
      rlamda(1,3,13) =  0.01311
      rlamda(1,3,14) =  0.00024
      rlamda(1,3,15) = -0.00651
      rlamda(1,3,16) =  0.00185
      rlamda(1,3,17) = -0.05622
      rlamda(1,3,18) = -0.00342
      rlamda(1,3,19) =  0.01769
      rlamda(1,3,20) = -0.00020
      rlamda(1,3,21) =  0.00219
      rlamda(1,3,22) = -0.02000
      rlamda(1,3,23) =  0.02829
      rlamda(1,3,25) = -0.08859
      rlamda(1,3,26) = -0.00015
      rlamda(1,3,27) =  0.03954
      rlamda(1,3,28) = -0.00879
      rlamda(1,3,29) = -0.00048
      rlamda(1,3,30) =  0.02818
      rlamda(1,3,31) =  0.08258
      rlamda(1,3,32) = -0.04090
      rlamda(1,3,33) =  0.00041
      rlamda(1,3,34) =  0.00041
      rlamda(1,3,35) = -0.00069
      rlamda(1,3,36) = -0.13269
      rlamda(1,3,37) = -0.00241
      rlamda(1,3,38) =  0.00037
      rlamda(1,3,39) = -0.01557
      rlamda(1,3,40) =  0.04191
      rlamda(1,3,41) = -0.00465
      rlamda(1,3,42) =  0.00085
      rlamda(1,3,43) = -0.04983
      rlamda(1,3,44) =  0.09946
      rlamda(1,3,46) =  0.16533
      rlamda(1,3,47) =  0.05411
      rlamda(1,3,48) =  0.00277
      rlamda(1,3,49) =  0.07859
      rlamda(1,3,50) = -0.00343
      rlamda(1,3,51) = -0.00810
      rlamda(1,3,52) = -0.00015
      rlamda(1,3,53) =  0.01106
      rlamda(1,3,54) = -0.00814
      rlamda(1,3,55) =  0.00304
      rlamda(1,3,57) = -0.00324
      rlamda(2,3,2 ) =  0.00195
      rlamda(2,3,4 ) =  0.00617
      rlamda(2,3,5 ) = -0.00123
      rlamda(2,3,6 ) =  0.01004
      rlamda(2,3,7 ) = -0.00077
      rlamda(2,3,9 ) =  0.00014
      rlamda(2,3,10) = -0.00338
      rlamda(2,3,11) =  0.01277
      rlamda(2,3,12) =  0.00013
      rlamda(2,3,14) =  0.02641
      rlamda(2,3,15) =  0.00099
      rlamda(2,3,16) =  0.00530
      rlamda(2,3,17) =  0.00031
      rlamda(2,3,20) = -0.00208
      rlamda(2,3,24) = -0.00278
      rlamda(2,3,25) = -0.00030
      rlamda(2,3,26) = -0.01234
      rlamda(2,3,27) =  0.00099
      rlamda(2,3,28) =  0.00252
      rlamda(2,3,29) = -0.01659
      rlamda(2,3,30) = -0.00026
      rlamda(2,3,31) = -0.00020
      rlamda(2,3,32) =  0.00057
      rlamda(2,3,33) = -0.05087
      rlamda(2,3,34) = -0.00962
      rlamda(2,3,35) =  0.03118
      rlamda(2,3,36) = -0.00021
      rlamda(2,3,37) =  0.01499
      rlamda(2,3,38) =  0.00439
      rlamda(2,3,39) = -0.00069
      rlamda(2,3,40) =  0.00299
      rlamda(2,3,41) =  0.02879
      rlamda(2,3,42) =  0.02474
      rlamda(2,3,43) = -0.00032
      rlamda(2,3,44) =  0.00058
      rlamda(2,3,45) = -0.09530
      rlamda(2,3,46) = -0.00035
      rlamda(2,3,47) = -0.00014
      rlamda(2,3,48) = -0.00565
      rlamda(2,3,49) =  0.00019
      rlamda(2,3,50) =  0.00185
      rlamda(2,3,51) = -0.00073
      rlamda(2,3,52) = -0.00199
      rlamda(2,3,54) =  0.00635
      rlamda(2,3,55) =  0.01709
      rlamda(2,3,56) = -0.00350

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
!     Calculates the state-independent force f0 = - w * q
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
         else if (sampnuc_opt.eq.'foc') then
            ! Focused (circle)
            q(i) = sin(p(i))
            p(i) = cos(p(i))
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


