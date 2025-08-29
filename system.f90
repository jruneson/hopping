module system
   use types
   implicit real(dp) (a-h,o-z)
!
!     ------------------------------------------------------------------
!     Module for a generic system.
!     ------------------------------------------------------------------
!
   integer, protected :: ns ! Number of states
   integer, protected :: nf ! Number of nuclear degrees of freedom
   integer :: nh ! Number of elements in h
   integer :: nv ! Number of dimensions of velocity
   integer :: nmult ! Number of spin-free states (each multiplet counts as one)
   real(dp) :: beta = 0.0  ! Inverse temperature
   real(dp) :: gamma = 0.0 ! Friction coefficient
   allocatable :: rmass(:) ! Masses
   allocatable :: gam(:)   ! Width parameters
   allocatable :: gammaT(:) ! Thermalization rate of Langevin thermostat
   character(3) :: model
   character(3) :: sampnuc_opt = 'wig'
   procedure (formh_interface), pointer :: formh
   procedure (image_interface), pointer :: image
   procedure (force0_interface), pointer :: force0
   procedure (force_interface), pointer :: force
   procedure (hmult_interface), pointer :: hmult
   procedure (vmult_interface), pointer :: vmult
   procedure (sampnuc_interface), pointer :: sampnuc
   procedure (potmat_interface), pointer :: potmat
   procedure (potii_interface), pointer :: pot_ii
   procedure (nac_interface), pointer :: nac
   procedure (pot_interface), pointer :: pot
   procedure (grad_interface), pointer :: grad
   procedure (h0_interface), pointer :: h0
   procedure (get_omega_interface), pointer :: get_omega
   procedure (sampdis_interface), pointer :: sampdis

   interface
      subroutine formh_interface(q,h)
         ! Diabatic potential matrix
         import dp
         real(dp), intent(in) :: q(:)
         real(dp), intent(out) :: h(:)
      end subroutine

      subroutine image_interface(h,c,d)
         ! Image of the electronic Hamiltonian
         import dp
         real(dp), intent(in) :: h(:)
         real(dp), intent(in) :: c(:)
         real(dp), intent(out) :: d(:)
      end subroutine

      subroutine force0_interface(q,f0)
         ! State-independent force
         import dp
         real(dp), intent(in) :: q(:)
         real(dp), intent(out) :: f0(:)
      end subroutine

      ! subroutine eforce_interface(c,f0)
      !    ! Ehrenfest force
      !    import dp,dpc
      !    complex(dpc), intent(in) :: c(:)
      !    real(dp), intent(out) :: f0(:)
      ! end subroutine

      subroutine force_interface(q,u,f0)
         ! Force in complex state vector u
         import dp, dpc
         real(dp), intent(in) :: q(:)
         complex(dpc), intent(in) :: u(:)
         real(dp), intent(out) :: f0(:)
      end subroutine

      subroutine hmult_interface(h,c,d)
         ! Forms |d> = H*|c>
         import dp,dpc
         real(dp), intent(in) :: h(:)
         complex(dpc), intent(in) :: c(:)
         complex(dpc), intent(out) :: d(:)
      end subroutine

      subroutine vmult_interface(h,c,d,idim)
         ! Forms |d> = V*|c>
         import dp,dpc
         optional :: idim
         real(dp), intent(in) :: h(:)
         complex(dpc), intent(in) :: c(:)
         complex(dpc), intent(out) :: d(:)
      end subroutine

      subroutine sampnuc_interface(p,q)
         ! Samples a random nuclear state
         import dp,dpc
         real(dp), intent(out) :: p(:), q(:)
      end subroutine

      subroutine potmat_interface(q,u,v)
         import dp,dpc
         real(dp), intent(in) :: q(:)
         complex(dpc), intent(out) :: u(:,:)
         real(dp), intent(out) :: v(:)
      end subroutine

      subroutine potii_interface(q,ua,va,u,v,l)
         import ns,dp
         real(dp), intent(in) :: q(:), ua(:), va
         real(dp), intent(out) :: u(ns,l), v(l)
      end subroutine

      subroutine nac_interface(q,u,v,ia,ib,d)
         ! Column ia of the non-adiabatic coupling vector
         import dp,dpc
         real(dp), intent(in) :: q(:), v(:)
         complex(dpc), intent(in) :: u(:,:)
         integer, intent(in) :: ia, ib
         complex(dpc), intent(out) :: d(:)
      end subroutine

      subroutine pot_interface(q,u)
         ! Diabatic potential matrix, shape (ns,ns)
         import dp,dpc
         real(dp), intent(in) :: q(:)
         complex(dpc), intent(out) :: u(:,:)
      end subroutine

      subroutine grad_interface(q,g)
         ! Diabatic gradient tensor, shape (nf,ns,ns)
         import dp,dpc
         real(dp), intent(in) :: q(:)
         complex(dpc), intent(out) :: g(:,:,:)
      end subroutine

      real(dp) function h0_interface(p,q)
         import dp
         real(dp), intent(in) :: p(:),q(:)
      end function

      real(dp) function get_omega_interface()
         ! Output omega (for SSH and PPP)
         import dp
      end function

      subroutine sampdis_interface()
         ! import n,dp
      end subroutine

   end interface

contains
   subroutine init_system(nf_,ns_,model_,nuc_opt_,beta_,gamma_)
      integer, intent(in) :: nf_,ns_
      character(3) :: model_, nuc_opt_
      optional :: beta_, gamma_
!
!     ------------------------------------------------------------------
!     Initialises the system
!     ------------------------------------------------------------------
!
      nf = nf_
      ns = ns_
      nh = ns ! Default: number of states = number of sparse Hamiltonian elements
      nv = 1 ! Default: number of velocity dimensions = 1
      model = model_
      sampnuc_opt = nuc_opt_
      if (present(beta_)) beta = beta_
      if (present(gamma_)) gamma = gamma_
      allocate(rmass(nf))
      rmass = 1.d0 ! Unit masses as default
      allocate(gam(nf))
      gam = 0.d0 ! No width as default
      allocate(gammaT(nf))
      gammaT = 0.d0 ! No thermalization rate as default
      potmat => potmat_default
      force0 => force0_default
      force => mforce_default
      nac => nac_default
      h0 => h0_default

   end subroutine

   subroutine getcvc (q,c,cvc)
      complex (dpc) c,d
!
!     ------------------------------------------------------------------
!     Calculates cvc = <c|V(q)|c>.
!     ------------------------------------------------------------------
!
      dimension q(nf),c(ns)
      allocatable :: h(:),d(:)
      allocate (h(nh),d(ns))
!
      print*, 'getcvc'
      call formh (q,h)
      print*, 'formh'
      call hmult (h,c,d)
      print*, 'hmult'
      cvc = real(dot_product(c,d))
      deallocate (h,d)
   end subroutine

   subroutine force0_default(q,f0)
!
!     ------------------------------------------------------------------
!     Calculates the state-independent force on the nuclei.
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: f0(:)
      f0 = 0.d0 * q
   end subroutine

   subroutine mforce_default (q,u,f)
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
      call grad(q,g)
      do i=1,nf
         f(i) = -real(sum(conjg(u)*matmul(g(i,:,:),u)))
      end do
      deallocate(g)
   end subroutine


   subroutine potmat_default (q,u,v)
      use maths, only : hermevp
!
!     ------------------------------------------------------------------ 
!     Diagonalises the electronic potential energy matrix V(q)
!     to obtain its eigenvectors U(q) and eigenvalues v(q).
!     Excludes harmonic part of potential.
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:)
      complex(dpc), intent(out) :: u(:,:)
      real(dp), intent(out) :: v(:)
!
      call pot(q,u)
      call hermevp (u,ns,ns,v,ierr)
      if (ierr .ne. 0) stop 'symevp'
   end subroutine


   subroutine naccol (q,u,v,ib,d)
!
!     ------------------------------------------------------------------ 
!     Calculate column ia of the nonadiabatic coupling vector
!     from the transformation matrix (u) and adiabatic energies (v).
!     Outputs d(ia,j) = <ia|dV/dq_j|ib> 
!     ------------------------------------------------------------------ 
!
      real(dp), intent(in) :: q(:), v(:)
      complex(dpc), intent(in) :: u(:,:)
      integer, intent(in) :: ib
      complex(dpc), intent(out) :: d(:,:) ! NOTE: (ns,nf)
      do ia=1,ns
         if (ia.eq.ib) then
            d(ia,:) = 0.d0
         else
            call nac(q,u,v,ia,ib,d(ia,:))
         end if
      end do
   end subroutine

   subroutine nacrow (q,u,v,ia,d)
!
!     ------------------------------------------------------------------ 
!     Calculate row ia of the nonadiabatic coupling vector
!     from the transformation matrix (u) and adiabatic energies (v).
!     Outputs d(j,ib) = <ia|dV/dq_j|ib>
!     ------------------------------------------------------------------
      real(dp), intent(in) :: q(:), v(:)
      complex(dpc), intent(in) :: u(:,:)
      integer, intent(in) :: ia
      complex(dpc), intent(out) :: d(:,:) ! NOTE: (nf,ns)
      do ib=1,ns
         if (ib.eq.ia) then
            d(:,ib) = 0.d0
         else
            call nac(q,u,v,ia,ib,d(:,ib))
         end if
      end do
   end subroutine

   subroutine grad_ab(q,u,ia,ib,g)
!
!     ------------------------------------------------------------------
!     Offdiagonal adiabatic gradient element
!     ------------------------------------------------------------------
!
      real(dp) :: q(:)
      complex(dpc) :: u(:,:)
      integer :: ia, ib
      complex(dpc) :: g(:)
      complex(dpc), allocatable :: gdia(:,:,:)
      allocate(gdia(nf,ns,ns))
      call grad(q,gdia)
      do i=1,nf
         g(i) = dot_product(u(:,ia),matmul(gdia(i,:,:),u(:,ib)))
      end do
      deallocate(gdia)
   end subroutine

   subroutine nac_default(q,u,v,ia,ib,d)
!
!     ------------------------------------------------------------------
!     Non-adiabatic coupling vector
!     ------------------------------------------------------------------
!
      real(dp), intent(in) :: q(:), v(:)
      complex(dpc), intent(in) :: u(:,:)
      integer, intent(in) :: ia, ib
      complex(dpc), intent(out) :: d(:)
!
      complex(dpc), allocatable :: g(:,:,:)
      allocate(g(nf,ns,ns))
      call grad(q,g)
      do i=1,nf
         d(i) = dot_product(u(:,ia),matmul(g(i,:,:),u(:,ib))) / (v(ib)-v(ia))
      end do
      deallocate(g)
   end subroutine

   real(dp) function h0_default(p,q)
      real(dp), intent(in) :: p(:),q(:)
      unused = q(1)
      h0_default = 0.5*sum(p**2/rmass)
   end function


end module
