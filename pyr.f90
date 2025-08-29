module pyr
   use types
   implicit real(dp) (a-h,o-z)

   allocatable :: w(:) ! Frequencies
   allocatable :: Vconst(:,:)
   allocatable :: Vlin(:,:,:)
   allocatable :: Vbilin(:,:,:,:)

   logical :: mass_rescaled = .false.


contains

   subroutine initpyr()
      use system, only : sys_force0=>force0, sys_force=>force
      use system, only : sys_sampnuc=>sampnuc
      use system, only : sys_h0=>h0, sys_pot=>pot, sys_grad=>grad
      use system, only : ns, nf, rmass, gam
!
!     ------------------------------------------------------------------
!     Initialise pyrazine module.
!     24D Input data from Raab, Worth, Meyer, Cederbaum
!     J. Chem. Phys. 110, 936 (1999)
!     ------------------------------------------------------------------
!
      allocatable :: ai(:) ! Linear coefficients (surface 1)
      allocatable :: bi(:) ! Linear coefficients (surface 2)
      allocatable :: ci(:) ! Linear coupling coefficients
      allocatable :: aij(:,:) ! Bilinear coefficients (surface 1)
      allocatable :: bij(:,:) ! Bilinear coefficients (surface 2)
      allocatable :: cij(:,:) ! Bilinear coupling coefficients

      if (nf.ne.24 .and. nf.ne.3) stop 'Pyrazine requires nf=24 or nf=3'
      if (ns.ne.2) stop 'Pyrazine requires ns=2'
      allocate(Vconst(ns,ns),Vlin(nf,ns,ns),Vbilin(nf,nf,ns,ns))
      allocate(w(nf), ai(nf), bi(nf), ci(nf), aij(nf,nf), bij(nf,nf), cij(nf,nf))

      w(:) = 0.d0
      ai(:) = 0.d0
      bi(:) = 0.d0
      ci(:) = 0.d0
      aij(:,:) = 0.d0
      bij(:,:) = 0.d0
      cij(:,:) = 0.d0
      
      if (nf.eq.3) then
         ! Reduced 3D model (Schneider and Domcke CPL 1988)
         w(1) = 0.126   ! 1
         w(2) = 0.074   ! 61
         w(3) = 0.118   ! 10a
         ai = (/ -0.037, -0.105, 0.0   /)
         bi = (/ -0.254,  0.149, 0.0   /)
         ci = (/  0.0,    0.0,   0.262 /)
         delta = (4.84-3.94)/2.0
      else
         ! Frequencies (Table I, experiment)
         w(1) = 0.0739   ! 6a
         w(2) = 0.1258   ! 1
         w(3) = 0.1525   ! 9a
         w(4) = 0.1961   ! 8a
         w(5) = 0.3788   ! 2
         w(6) = 0.1139   ! 10a
         w(7) = 0.0937   ! 4
         w(8) = 0.1219   ! 5
         w(9) = 0.0873   ! 6b
         w(10) = 0.1669  ! 3
         w(11) = 0.1891  ! 8b
         w(12) = 0.3769  ! 7b
         w(13) = 0.0423  ! 16a
         w(14) = 0.1190  ! 17a
         w(15) = 0.1266  ! 12
         w(16) = 0.1408  ! 18a
         w(17) = 0.1840  ! 19a
         w(18) = 0.3734  ! 13
         w(19) = 0.1318  ! 18b
         w(20) = 0.1425  ! 14
         w(21) = 0.1756  ! 19b
         w(22) = 0.3798  ! 20b
         w(23) = 0.0521  ! 16b
         w(24) = 0.0973  ! 11


         ! Linear tuning coefficients (Table II, CIS)
         ai(1) = -0.0981
         ai(2) = -0.0503
         ai(3) =  0.1452
         ai(4) = -0.0445
         ai(5) =  0.0247

         bi(1) =  0.1355
         bi(2) = -0.171 ! Table V
         bi(3) =  0.0375
         bi(4) =  0.0168
         bi(5) =  0.0162


         ! Diagonal bilinear tuning coefficients - surface 1.
         ! (Table III, CIS - note that Au modes breaks order in the table)
         aij(1,1) = 0.d0 ! Comment (2) for modes 1-5
         aij(2,2) = 0.d0
         aij(3,3) = 0.d0
         aij(4,4) = 0.d0
         aij(5,5) = 0.d0
         aij(6,6) = -0.01159
         aij(7,7) = -0.02252 ! Table V
         aij(8,8) = -0.01825

         aij(9,9) = -0.00741
         aij(10,10) =  0.05183
         aij(11,11) = -0.05733
         aij(12,12) = -0.00333

         aij(13,13) =  0.01145 ! Table V
         aij(14,14) = -0.02040

         aij(15,15) = -0.04819
         aij(16,16) = -0.00792
         aij(17,17) = -0.02429
         aij(18,18) = -0.00492
         aij(19,19) = -0.00277 ! Table V
         aij(20,20) =  0.03924 ! Table V
         aij(21,21) =  0.00992
         aij(22,22) = -0.00110
         aij(23,23) = -0.02176
         aij(24,24) =  0.00315


         ! Off-diagonal bilinear tuning coefficients - surface 1.
         aij(1,2) =  0.00108
         aij(1,3) = -0.00204
         aij(1,4) = -0.00135
         aij(1,5) = -0.00285

         aij(2,3) =  0.00474
         aij(2,4) =  0.00154
         aij(2,5) = -0.00163

         aij(3,4) =  0.00872
         aij(3,5) = -0.00474

         aij(4,5) = -0.00143

         aij(7,8) = -0.00049

         aij(9,10) =  0.01321
         aij(9,11) = -0.00717
         aij(9,12) =  0.00515

         aij(10,11) = -0.03942
         aij(10,12) =  0.00170

         aij(11,12) = -0.00204

         aij(13,14) =  0.00100

         aij(15,16) =  0.00525
         aij(15,17) = -0.00485
         aij(15,18) = -0.00326
         aij(16,17) =  0.00852
         aij(16,18) =  0.00888
         aij(17,18) = -0.00443

         aij(19,20) =  0.00016 ! Table V
         aij(19,21) = -0.00250
         aij(19,22) =  0.00357

         aij(20,21) = -0.00197
         aij(20,22) = -0.00355
         aij(21,22) =  0.00623

         aij(23,24) = -0.00624


         ! Diagonal bilinear tuning coefficients - surface 2.
         bij(1,1) = 0.d0 ! Comment (2) for modes 1-5
         bij(2,2) = 0.d0
         bij(3,3) = 0.d0
         bij(4,4) = 0.d0
         bij(5,5) = 0.d0
         bij(6,6) = -0.01159
         bij(7,7) = -0.03445 ! Table V
         bij(8,8) = -0.00265
         bij(9,9)   = -0.00385
         bij(10,10) =  0.04842
         bij(11,11) = -0.06332
         bij(12,12) = -0.00040
         bij(13,13) = -0.01459 ! Table V
         bij(14,14) = -0.00618
         bij(15,15) = -0.00840
         bij(16,16) =  0.00429
         bij(17,17) = -0.00734
         bij(18,18) =  0.00062
         bij(19,19) = -0.01179 ! Table V
         bij(20,20) =  0.04000 ! Table V
         bij(21,21) =  0.01246
         bij(22,22) =  0.00069
         bij(23,23) = -0.02214
         bij(24,24) = -0.00496


         ! Off-diagonal bilinear tuning coefficients - surface 2.
         bij(1,2) = -0.00298
         bij(1,3) = -0.00189
         bij(1,4) = -0.00203
         bij(1,5) = -0.00128

         bij(2,3) =  0.00155
         bij(2,4) =  0.00311
         bij(2,5) = -0.00600

         bij(3,4) =  0.01194
         bij(3,5) = -0.00334

         bij(4,5) = -0.00713

         bij(7,8) = 0.00911

         bij(9,10)  = -0.00661
         bij(9,11)  =  0.00429
         bij(9,12)  = -0.00246

         bij(10,11) = -0.03034
         bij(10,12) = -0.00185

         bij(11,12) = -0.00388

         bij(13,14) = -0.00091

         bij(15,16) =  0.00536
         bij(15,17) = -0.00097
         bij(15,18) =  0.00034

         bij(16,17) =  0.00209
         bij(16,18) = -0.00049

         bij(17,18) =  0.00346

         bij(19,20) = -0.00884 ! Table V
         bij(19,21) =  0.07000 ! Table V
         bij(19,22) = -0.01249

         bij(20,21) = -0.05000 ! Table V
         bij(20,22) =  0.00265

         bij(21,22) = -0.00422

         bij(23,24) = -0.00261
      

         ! Linear coupling coefficient 
         ci(6) = 0.2080

         ! Bilinear coupling coefficient
         cij(1,6) = -0.01000 ! Table V
         cij(2,6) = -0.00551
         cij(3,6) =  0.00127
         cij(4,6) =  0.00799
         cij(5,6) = -0.00512

         cij(7,9)  = -0.01372
         cij(7,10) = -0.00466
         cij(7,11) =  0.00329
         cij(7,12) = -0.00031

         cij(8,9)  =  0.00598
         cij(8,10) = -0.00914
         cij(8,11) =  0.00961
         cij(8,12) =  0.00500

         cij(13,15) = -0.01056
         cij(13,16) =  0.00559
         cij(13,17) =  0.00401
         cij(13,18) = -0.00226

         cij(14,15) = -0.01200
         cij(14,16) = -0.00213
         cij(14,17) =  0.00328
         cij(14,18) = -0.00396

         cij(19,23) =  0.00118
         cij(20,23) = -0.00009
         cij(21,23) = -0.00285
         cij(22,23) = -0.00095

         cij(19,24) =  0.01281
         cij(20,24) = -0.01780
         cij(21,24) =  0.00134
         cij(22,24) = -0.00481


         ! Half energy splitting between S1 and S2 (experiment)
         delta = 0.423d0


         ! Symmetrize
         do i = 1, 23
            do j = i+1,24
               aij(j,i) = aij(i,j)
               bij(j,i) = bij(i,j)
               cij(j,i) = cij(i,j)
            enddo
         enddo
      end if

      ! Convert units from eV to Hartree.
      ev_to_au = 1.0/27.2113961
      w(:) = w(:) * ev_to_au
      ai(:) = ai(:) * ev_to_au
      bi(:) = bi(:) * ev_to_au
      ci(:) = ci(:) * ev_to_au
      aij(:,:) = aij(:,:) * ev_to_au
      bij(:,:) = bij(:,:) * ev_to_au
      cij(:,:) = cij(:,:) * ev_to_au
      delta = delta * ev_to_au
      

      ! Define mass as the inverse of the frequencies
      rmass = 1.d0  / w

      ! Width parameter is one with the above mass definition
      gam = 1.d0

      if (mass_rescaled) then
         ! Definition in Mannouch & Richardson JCP 158, 104111 (2023)
         ! (gives equivalent dynamics)
         ai = ai * sqrt(w)
         bi = bi * sqrt(w)
         ci = ci * sqrt(w)
         do i=1,nf
            aij(:,i) = aij(:,i) * sqrt(w) * sqrt(w(i))
            bij(:,i) = bij(:,i) * sqrt(w) * sqrt(w(i))
            cij(:,i) = cij(:,i) * sqrt(w) * sqrt(w(i))
         end do
         rmass = 1.d0
         gam = rmass * w
      end if


      ! Put pieces together in Vconst, Vlin and Vbilin
      Vconst = 0.d0
      Vconst(1,1) = -delta
      Vconst(2,2) = delta
      Vlin(:,1,1) = ai(:)
      Vlin(:,2,2) = bi(:)
      Vlin(:,1,2) = ci(:)
      Vlin(:,2,1) = ci(:)
      Vbilin(:,:,1,1) = aij
      Vbilin(:,:,2,2) = bij
      Vbilin(:,:,1,2) = cij
      Vbilin(:,:,2,1) = cij


      sys_force0=>force0
      sys_force=>force
      sys_sampnuc=>sampnuc
      sys_pot=>potenl_dia
      sys_grad=>gradient_dia
      sys_h0=>h0

      deallocate(ai,bi,ci,aij,bij,cij)
   end subroutine

   subroutine force0 (q,f0)
      use system, only : rmass
!
!     ------------------------------------------------------------------
!     Calculates the state-independent force f0 = - w * q
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: f0(:)
      f0 = - rmass*w**2 * q
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
      use system, only : ns, nf
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
            do j=1,nf
               u(n,m) = u(n,m) + sum(Vbilin(:,j,n,m)*q(j)*q(:))
            end do
         end do
      end do
   end subroutine

   subroutine gradient_dia(q,g)
      use system, only : nf
!
!     ------------------------------------------------------------------
!     Diabatic potential gradient (excluding harmonic part)
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: g(:,:,:) ! (nf,ns,ns)
      do i=1,nf
         g(i,:,:) = Vlin(i,:,:)
         do j=1,nf
            g(i,:,:) = g(i,:,:) + 2.d0*Vbilin(i,j,:,:)*q(j)
         end do
      end do
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
            if (mass_rescaled) then
               deltap = sqrt(0.5*w(i))
               deltaq = sqrt(0.5/w(i))
            end if
         else
            stop 'Unknown sampnuc_opt'
         end if
         p(i) = deltap*p(i)
         q(i) = deltaq*q(i)
      end do
   end subroutine

   function h0 (p,q)
      use system, only : nf,rmass
!
!     ------------------------------------------------------------------ 
!     Calculates state-independent (bath) energy.
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: p(:),q(:)
      h0 = 0.d0
      do i=1,nf
         ! h0 = h0 + w(i)*(p(i)**2+q(i)**2)/2.d0
         h0 = h0 + 0.5*(p(i)**2/rmass(i) + rmass(i)*w(i)**2*q(i)**2)
      end do
   end function
end module
