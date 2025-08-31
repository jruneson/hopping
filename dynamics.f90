module dynamics
   use types
   implicit real(dp) (a-h,o-z)
!
!     ------------------------------------------------------------------
!     Module for running mixed quantum-classical dynamics.
!     Currently supports FSSH, MASH, Ehrenfest, CPA (=no back-action).
!     ------------------------------------------------------------------
!
   character :: method 
   character :: hop_opt ! 'd' for NACV probability, 't' for TDC probability with finite difference,
                        ! 'e' for Ehrenfest
   character :: rescale_opt
   character :: dec_opt 
   real(dp) :: dec_param 
   real(dp) :: alphan ! MASH constant

contains

   subroutine initdyn (method_,hop_opt_,rescale_opt_,dec_opt_,dec_param_)
      use system, only : ns
!
!     ------------------------------------------------------------------
!     Initialize method ('m' for MASH, 'f' for FSSH).
!     ------------------------------------------------------------------
!
      character :: method_
      character :: hop_opt_
      character :: rescale_opt_
      character :: dec_opt_
!
      method = method_
      hop_opt = hop_opt_
      rescale_opt = rescale_opt_
      dec_opt = dec_opt_
      dec_param = dec_param_

      ! Set MASH constant
      alphan = get_alpha(ns)
   end subroutine

   real(dp) function get_alpha (n)
      if (n.eq.1) then 
         get_alpha = 1
      else
         hn = 0.0
         do i = 1,n
            hn = hn + 1.0/i
         end do
         get_alpha = (n-1.0)/(hn-1.0)
      end if
   end function

!
!  =================== SINGLE TRAJECTORY SUBROUTINES ===================
!

   subroutine ftraj (pt,qt,ct,cat,iat,ut,vt,et,nt,dt)
      use system, only : nf, ns, potmat, force0, force, h0, model
      use tul, only : branching
      complex(dpc) c,ca,ct,cat,u,ut,T
!
!     ------------------------------------------------------------------
!     Single FSSH trajectory starting from pt(0),qt(0),ct(0),iat(0)
!     ------------------------------------------------------------------
!
      dimension pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt)
      dimension iat(0:nt),ut(ns,ns,0:nt),vt(ns,0:nt),et(0:nt)
      allocatable :: p(:),q(:),c(:),ca(:),u(:,:),v(:),f0(:),f(:),T(:,:)
      if (method.ne.'f') stop 'ftraj is only to be used with method="f"'
      allocate (p(nf),q(nf),c(ns),ca(ns),u(ns,ns),v(ns),f0(nf),f(nf),T(ns,ns))
!
      p = pt(:,0)
      q = qt(:,0)
      c = ct(:,0)
      ia = iat(0)
      call potmat (q,u,v)
      call force0 (q,f0)
      if (hop_opt.eq.'e') then
         call force (q,c,f)
      else
         call force (q,u(:,ia),f)
      end if
!
      dt2 = dt/2
      ca = matmul(c,conjg(u))
      do it = 0,nt
         if (it.gt.0) then
            call fstep (p,q,c,ca,u,v,f0,f,T,ia,dt)
         end if
         pt(:,it) = p
         qt(:,it) = q
         ct(:,it) = c
         cat(:,it) = ca
         iat(it) = ia
         ut(:,:,it) = u
         vt(:,it) = v
         if (hop_opt.eq.'e' .or. dec_opt.eq.'e') then
            et(it) = sum(abs(ca)**2 * v) + h0(p,q)
         else
            et(it) = v(ia) + h0(p,q)
         end if
         ! HACK FOR EARLY TERMINATION OF TULLY BRANCHING
         if (model.eq.'tul' .and. branching) then
            if (abs(q(1)).gt.20) then
               qt(:,nt) = q
               cat(:,nt) = ca
               iat(nt) = ia
               return
            end if
         end if
         ! END TULLY BRANCHING
      end do
      deallocate (p,q,c,ca,u,v,f0,f,T)
   end subroutine

   subroutine mtraj (pt,qt,ct,cat,iat,ut,vt,et,nt,dt)
      use system, only : nf, ns, potmat, force0, force, h0
      complex(dpc) c,ct,cat,u,u0,ut
!
!     ------------------------------------------------------------------
!     Single MASH trajectory starting from pt(0),qt(0),ct(0),iat(0)
!     ------------------------------------------------------------------
!
      dimension pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt)
      dimension iat(0:nt),ut(ns,ns,0:nt),vt(ns,0:nt),et(0:nt)
      allocatable :: p(:),q(:),c(:),u0(:,:),u(:,:),v(:),f0(:),f(:)
      if (method.ne.'m') stop 'mtraj is only to be used with method="m"'
      allocate (p(nf),q(nf),c(ns),u0(ns,ns),u(ns,ns),v(ns),f0(nf),f(nf))
!
      p = pt(:,0)
      q = qt(:,0)
      c = ct(:,0)
      ia = iat(0)
      call potmat (q,u,v)
      call force0 (q,f0)
      call force (q,u(:,ia),f)
!
      dt2 = dt/2
      do it = 0,nt
         if (it.gt.0) then
            call mstep (p,q,c,u,v,f0,f,ia,dt)
         end if
         pt(:,it) = p
         qt(:,it) = q
         ct(:,it) = c
         cat(:,it) = matmul(c,u)
         iat(it) = ia
         ut(:,:,it) = u
         vt(:,it) = v
         et(it) = v(ia) + h0(p,q)
      end do
      deallocate (p,q,c,u0,u,v,f0,f)
   end subroutine

!
!  ========================= COMMON SUBROUTINES ========================
!

   subroutine sampad (v,ia)
      use system, only : ns, beta
!
!     ------------------------------------------------------------------
!     Samples the adiabatic state (ia) from a Boltzmann distribution.
!     Uses the eigenvalues of the Hamiltonian (v).
!     ------------------------------------------------------------------
!
      dimension v(ns)
      allocatable :: f(:)
      allocate(f(ns))
!
!     Sample an electronic state ia from the Boltzmann distribution
!
      call random_number (r)
      z = 0.0
      do j = 1,ns
         f(j) = exp(-beta*(v(j)-v(1)))
         z = z+f(j)
      enddo
      f = f/z
      s = 0.0
      ia = 0
      do j = 1,ns
         ia = j
         s = s+f(j)
         if (s .ge. r) exit
      enddo
      deallocate(f)
   end subroutine

!
!  ====================== MASH and FSSH SUBROUTINES ====================
!

   subroutine evolve (p,q,c,u,v,f0,f,ia,dt)
      use system, only : potmat, force0, force, nf, ns, rmass
      complex (dpc) c,u
!
!     ------------------------------------------------------------------
!     Evolves p, q, and c for a time interval dt.
!     Central q version for MASH and FSSH 
!     (involves diagonalization once per step)
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),c(ns),u(ns,ns),v(ns),f0(nf),f(nf)
!
!     Evolve c for dt/2
!
      dt2 = dt/2
      call cstep_exp (c,u,v,dt2)
!
!     Evolve p for dt/2
!
      p = p + dt2*(f0+f)
!
!     Evolve q for dt
!
      q = q+dt*p/rmass
!
!     Calculate the new adiabatic states and the new nuclear forces
!
      call potmat (q,u,v)
      call force0 (q,f0)
      if (hop_opt.eq.'e') then
         ! Ehrenfest force
         call force (q,c,f)
      else
         ! Force on the active state
         call force (q,u(:,ia),f)
      end if
!
!     Evolve p for dt/2
!
      p = p + dt2*(f0+f)
!
!     Evolve c for dt/2
!
      call cstep_exp (c,u,v,dt2)
   end subroutine

   subroutine cstep_exp(c,u,v,dt)
      use system, only : ns
      use maths, only : iu
      complex (dpc) c,u
!
!     ------------------------------------------------------------------
!     Uses the exponential propagator to evolve c through dt.
!     ------------------------------------------------------------------
!
      dimension c(ns),u(ns,ns),v(ns)
!
      ! Convert to eigenbasis
      c = matmul(conjg(transpose(u)),c)
      ! Propagate
      c = exp(-iu*v*dt)*c
      ! Convert back to diabatic basis
      c = matmul(u,c)
   end subroutine


   subroutine fstep (p,q,c,ca,u,v,f0,f,T,ia,dt)
      use system, only : nf, ns, nacrow, rmass
      complex (dpc) c,ca,u,u0,d,T,phase
!
!     ------------------------------------------------------------------
!     FSSH step
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),c(ns),ca(ns),u(ns,ns),v(ns),f0(nf),f(nf),T(ns,ns)
      allocatable :: u0(:,:),d(:,:)
      allocate(u0(ns,ns),d(nf,ns))
      u0 = u
      ! Velocity Verlet step
      call evolve(p,q,c,u,v,f0,f,ia,dt)

      ! Fix phase of eigenvectors before calculating NACV and TDC
      do i=1,ns
         ! First version based on real Hamiltonians:
         ! if (real(dot_product(u0(:,i),u(:,i))).lt.0.d0) u(:,i) = -u(:,i)
         ! General fix by Akimov, JPCL 2018:
         phase = dot_product(u0(:,i),u(:,i))
         phase = phase / abs(phase)
         u(:,i) = u(:,i) * conjg(phase)
      end do

      ! Calculate nonadiabatic coupling vector d(:,ia,ib) = <\psi_ia | d/dq | \psi_ib>
      ! (we only need the ia-th row)
      call nacrow (q,u,v,ia,d)

      ! Calculate time-derivative couplings T(i,j) = <\psi_i | d/dt | \psi_j>
      if (hop_opt.eq.'t') then
         ! Calculate overlap matrix 
         do i=1,ns
            do j=1,ns
               T(i,j) = dot_product(u0(:,i),u(:,j)) - dot_product(u(:,i),u0(:,j))
            end do
         end do
         T = T/(2*dt)
      else
         ! Calculate T from the NACV
         do ib=1,ns
            T(ia,ib) = sum(d(:,ib)*p/rmass)
            T(ib,ia) = -T(ia,ib)
         end do
      end if

      ! Hop and decohere
      ca = matmul(transpose(conjg(u)),c)
      call fs_hop (p,q,ca,u,v,f,ia,d,T,dt)
      call decohere(ca,p,q,u,v,f,ia,d,T,dt)
      c = matmul(u,ca)

      deallocate(u0,d)
   end subroutine

   subroutine mstep (p,q,c,u,v,f0,f,ia,dt)
      use system, only : nf, ns
      complex (dpc) c,u,u0
!
!     ------------------------------------------------------------------
!     MASH step
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),c(ns),u(ns,ns),v(ns),f0(nf),f(nf)
      allocatable :: u0(:,:)
      allocate (u0(ns,ns))

      u0 = u
      call evolve (p,q,c,u,v,f0,f,ia,dt)
      call smash (p,q,c,u,v,f,ia)
      deallocate (u0)
   end subroutine
   
   subroutine fs_hop (p,q,ca,u,v,f,ia,d,T,dt)
      use system, only : nf,ns,force,nacrow,rmass
      complex (dpc) ca,u,d,T
!
!     ------------------------------------------------------------------
!     FSSH hopping algorithm.
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),ca(ns),u(ns,ns),v(ns),f(nf),d(nf,ns),T(ns,ns)
      allocatable :: prob(:)

      allocate(prob(ns))
      ! Store initial state
      ia0 = ia
      ! Calculate probabilities
      prob = 0.0
      do ib=1,ns
         if (ib.eq.ia) cycle
         select case (hop_opt)
         case ('d','t')
            prob(ib) = max(0.d0,2*dt*real(ca(ib)/ca(ia)*T(ia,ib)))
         case ('e')
            ! Ehrenfest: don't hop
            exit
         case default 
            stop 'Unknown hopping option'
         end select
      end do
      ! Try hops
      call random_number(r)
      cumprob = 0.0
      do ib=1,ns
         cumprob = cumprob + prob(ib)
         if (r .lt. cumprob) then
            ! Hop to ib if energetically allowed
            call rescale (p,ca,v,d(:,ib),ia,ib)
            exit
         end if
      end do
      ! Update force and nac if hop occurred
      if (ia.ne.ia0) then
         call force(q,u(:,ia),f)
         select case (dec_opt)
         case ('o','i','v')
            ! Update NAC vector when needed for rescaling
            call nacrow(q,u,v,ia,d)
            if (hop_opt.eq.'d') then
               ! Update T
               do ib=1,ns
                  T(ia,ib) = sum(d(:,ib)*p/rmass)
                  T(ib,ia) = -T(ia,ib)
               end do
            end if
         end select
      end if
      deallocate(prob)
   end subroutine

   subroutine rescale (p,ca,v,d,ia,ib)
      use system, only : nf,ns,nac,rmass,force,grad_ab
      complex (dpc) ca,d
!
!     ------------------------------------------------------------------
!     Velocity rescaling in FSSH.
!     ------------------------------------------------------------------
!
      dimension :: p(nf),ca(ns),v(ns),d(nf)
      allocatable :: dir(:), pnac(:)
      if (ia.eq.ib) return
      dv = v(ib) - v(ia)
      select case(rescale_opt)
      case('d')
         ! Transform momentum to mass-scaled coordinates
         p = p/sqrt(rmass)
         ! Find direction of NAC vector in mass-scaled coordinates
         allocate(dir(nf),pnac(nf))
         dir = real(conjg(ca(ia))*d*ca(ib))
         dir = dir / sqrt(rmass)
         dir = dir / sqrt(sum(abs(dir)**2))
         ! Find momentum projection along NAC vector
         pnac = sum(p*dir)*dir
         ! Check that there is enough kinetic energy to hop
         Eknac = 0.5*sum(pnac**2)
         if (Eknac.gt.dv) then
            ! Succesful hop: Rescale along NAC vector
            p = p - pnac
            p = p + sqrt(2.d0*(Eknac-dv)) * pnac / sqrt(sum(pnac**2))
            ! Update active state
            ia = ib
         else
            ! Frustrated hop: Reverse along NAC vector
            p = p - 2*pnac
         end if
         ! Undo mass scaling
         p = p*sqrt(rmass)
         deallocate(dir,pnac)
      case('h') 
         !
         ! Momentum rescaling along d as in Hammes-Schiffer & Tully 1994
         ! (standard choice in many FSSH implementations)
         ! NOTE: only for real Hamiltonians
         !
         a = 0.5*sum(abs(d)**2/rmass)
         b = sum(p/rmass*real(d))
         c = dv
         det = b**2 - 4*a*c
         if (det.lt.0.0) then
            ! Frustrated hop, reversal
            alph = b/a
            p = p - alph*real(d)
            return 
         end if
         if (b.lt.0) then 
            alph = (b+sqrt(det))/(2*a)
         else 
            alph = (b-sqrt(det))/(2*a)
         end if
         p = p - alph*real(d)
         ia = ib
      case ('v','p')
         !
         !  Rescale along entire velocity
         !
         ekin = 0.5*sum(p**2/rmass)
         ekin_new = ekin - dv
         if (ekin_new.ge.0) then
            p = p * sqrt(ekin_new/ekin)
            ia = ib
         else
            p = -p
         end if
      case ('n')
         ! No rescaling, just accept all hops
         ia = ib
      case default 
         stop 'Unknown rescaling option'
      end select
   end subroutine

   subroutine decohere (c,p,q,u,v,f,ia,dvec,Tmat,dt)
      use system, only : nf, ns, nac, force, rmass, grad, gam
      use maths, only : pi
      complex (dpc) c,u,d,dvec,Tmat
!
!     ------------------------------------------------------------------
!     Apply decoherence correction to the active state.
!     c -- input: adiabatic wavefunction
!     ------------------------------------------------------------------
!
      intent(in) :: ia
      dimension c(ns),p(nf),q(nf),u(ns,ns),v(ns),f(nf),dvec(nf,ns),Tmat(ns,ns)
      allocatable :: fb(:),d(:),dir(:),pnac(:),dvel(:),dmom(:)
      logical :: allowed
      if (dec_opt.eq.'n') return ! No decoherence
      ! Store initial density matrix: R for off-diagonal and S for diagonal elements
!
      allocate(fb(nf),d(nf),dir(nf),pnac(nf),dvel(nf),dmom(nf))   
!
!     Perform decoherence correction
!
      do ib=1,ns
         if (ib.eq.ia) cycle
         dv = v(ib)-v(ia)
         r = abs(Tmat(ia,ib))/abs(dv)
         select case (dec_opt)
         case ('e','c')
            ! Energy-based decoherence correction
            Ekin = 0.5*sum(p**2/rmass)
            tau = (1.0+dec_param/Ekin)/abs(dv)
            c(ib) = exp(-dt/tau)*c(ib)
         case ('i')
            ! Instantaneous decoherence
            if (r.gt.dec_param .and. dec_param.gt.0.0) cycle
            c(ib) = 0.0
         case ('o')
            ! Check that states are uncoupled before applying decoherence
            if (r.gt.dec_param .and. dec_param.gt.0) cycle

            if (abs(Tmat(ia,ib))*dt.lt.1.d-15) then
               ! If states are entirely uncoupled, then it is unsafe to do 
               ! rescaling. In this situation, we instead want to decohere instantaneously
               c(ib) = 0.0
               cycle
            end if

            ! Estimate decoherence rate based on overlap decay
            ! exp(-k1*t - (k2*t)^2)
            !
            ! First calculate momentum difference that would arise from a hop
            ! Use mass-scaled variables
            dir = dreal(dvec(:,ib)) / sqrt(rmass)
            dir = dir / sqrt(sum(abs(dir)**2))
            ! Calculate momentum along NAC vector
            pnac = sum(p*dir/sqrt(rmass))*dir 
            ! Check that there is enough energy to hop
            Eknac = 0.5*sum(pnac**2)
            pnac = pnac * sqrt(rmass)
            dv = v(ib)-v(ia)
            if (Eknac.gt.dv) then  
               ! Rescale along NAC vector
               dmom = pnac*(sqrt(1-dv/Eknac)-1) 
               allowed = .true.
            else
               dmom = -2*pnac 
               allowed = .false.
            end if
            if (allowed) then
               ! Calculate force at state b
               call force(q,u(:,ib),fb)
            else 
               ! Use force at state a
               fb = f
            end if
            ! k1 decoherence rate: (1/gam) * dmom * dF
            b = 0.0
            do i=1,nf
               b = b + (1.d0/gam(i))*dmom(i)*(fb(i)-f(i)) 
            end do
            ! k2 decoherence rate: a = 0.5*(gam * dmom^2/m^2 + 1/gam * dF^2)
            a = 0.0
            do i=1,nf
               a = a + gam(i)*(dmom(i)/rmass(i))**2 + 1.d0/gam(i)*(fb(i)-f(i))**2 
            end do
            a = 0.5*a

            ! Fit exp(-b*t-a*t^2) to exp(-2kt)(2-exp(-2kt)) w.r.t. k
            y = b/(2*sqrt(a))
            tmp = -b/(2*a) + exp(-y**2)/(sqrt(pi*a)*erfc(y))
            rate = 7.0/(12.0*tmp)

            ! Apply decoherence
            c(ib) = exp(-rate*dt)*c(ib)
         case ('v')
            if (r.gt.dec_param .and. dec_param.gt.0) cycle
            ! Estimate decoherence rate based on just velocity difference over wavepacket width
            ! Calculate direction of NAC vector
            call nac (q,u,v,ia,ib,d)
            dir = dreal(d)
            ! Use mass-scaled variables
            p = p / sqrt(rmass)
            dir = dir / sqrt(rmass)
            dir = dir / sqrt(sum(abs(dir)**2))
            ! Calculate momentum along NAC vector
            pnac = sum(p*dir)*dir
            ! Undo mass-scaling
            p = p * sqrt(rmass)
            ! Check that there is enough energy to hop
            Eknac = 0.5*sum(pnac**2)
            dv = v(ib)-v(ia)
            if (Eknac.gt.dv) then  
               ! Rescale along NAC vector
               dvel = pnac*(sqrt(1-dv/Eknac)-1) / sqrt(rmass)
            else
               ! Reverse along NAC vector
               dvel = -2*pnac / sqrt(rmass)
            end if
            ! Calculate decoherence rate
            rate = 7.0/12.0 * sqrt(pi/2.0*sum(gam*dvel**2)) 
            
            c(ib) = exp(-rate*dt)*c(ib)
         case default
            stop 'Unknown decoherence option'
         end select
      end do
      c = c/sqrt(sum(abs(c)**2))
!
      deallocate(fb,d,dir,pnac,dvel,dmom)
   end subroutine


   subroutine smash (p,q,c,u,v,f,ia)
      use system, only : nf, ns, force
      complex (dpc) c,cad,u
!
!     ------------------------------------------------------------------
!     Multi-state MASH hopping algorithm.
!     Find new most populated state and hop to it if possible.
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),c(ns)
      dimension u(ns,ns),v(ns),f(nf)
      allocatable :: cad(:)
      call maxpop (u,c,ib) ! Find most populated adiabatic state
      if (ib .eq. ia) return
      allocate(cad(ns))
      cad = matmul(conjg(transpose(u)),c) ! Transform to adiabatic basis
      call mrescale (p,q,cad,u,v,ia,ib) ! Attempt to hop to ib
      if (ib .eq. ia) call force (q,u(:,ib),f) ! Update force if hop occurred
      deallocate (cad)
   end subroutine


   subroutine mrescale (p,q,cad,u,v,ia,ib)
      use system, only : nf, ns, naccol, rmass
      complex (dpc) cad,u,a,b
!
!     ------------------------------------------------------------------
!     MASH rescaling algorithm.
!     Check if hop from ia to ib is energetically allowed.
!     If yes, rescale momentum along d. 
!     If not, reverse momentum along d.
!     ------------------------------------------------------------------
!
      dimension p(nf),q(nf),cad(ns)
      dimension u(ns,ns),v(ns)
      allocatable :: a(:,:),b(:,:),d(:),pp(:),po(:)
      if (ia.eq.ib) return ! No hop necessary
      allocate (a(ns,nf),b(ns,nf),d(nf),pp(nf),po(nf))
!
!     Construct columns ia and ib of the non-adiabatic coupling matrix
!
      call naccol (q,u,v,ia,a)
      call naccol (q,u,v,ib,b)
!
!     Construct the direction d of momentum rescaling/reversal
!
      do j = 1,nf
         d(j) =  sum(real(conjg(cad(:))*a(:,j)*cad(ia))) &
               - sum(real(conjg(cad(:))*b(:,j)*cad(ib)))
      end do
      d = d/sqrt(sum(d**2))
!
!     Decompose p into components parallel and orthogonal to d
!
      p = p/sqrt(rmass)
      pd = sum(p*d)
      pp = pd * d
      po = p - pp
!
!     Check if there is enough kinetic energy parallel to d to hop
!
      told = sum(pp**2)
      if (told.gt.2*(v(ib)-v(ia))) then
!
!        If so, then hop
!
         tnew = told - 2*(v(ib)-v(ia))
         scale = sqrt(tnew/told)
         pp = scale * pp
         p = po + pp
         ia = ib
      else
!
!        If not, then reverse the direction of pp instead
!
         p = po - pp
      end if
      p = p*sqrt(rmass)
      deallocate (a,b,d,pp,po)
   end subroutine

   subroutine maxpop(u,c,ia)
      use system, only : ns
      complex(dpc) c(ns), u(ns,ns)
!
!     ------------------------------------------------------------------
!     Find the most populated adiabatic state.
!     ------------------------------------------------------------------
!
      popmax = 0.0
      ia = 0
      do ib = 1,ns
         pop = abs(sum(u(:,ib)*c))**2
         if (pop>popmax) then
            popmax = pop
            ia = ib
         end if
      end do
      if (ia.eq.0) stop 'No valid adiabatic states found in maxpop'
   end subroutine

   subroutine sampcap (ia,c)
      use maths, only : sampsun
      use system, only : ns
      complex (dpc) c, cswap
!
!     ------------------------------------------------------------------
!     Sample an electronic state (c) from the cap of an
!     adiabatic state ia of the Hamiltonian.
!     ------------------------------------------------------------------
!
      dimension c(ns)
!
!     Generate a random SU(N) electronic state
!
      call sampsun (ns,c)
!
!     Check which state (k) is currently most populated
!
      k = 0
      r = 0.0
      do j = 1,ns
         t = abs(c(j))**2
         if (t .gt. r) then
            k = j
            r = t
         end if
      end do
!
!     Permute state indices such that ia becomes the most populated state.
!
      if (ia .ne. k) then
         cswap = c(k)
         c(k) = c(ia)
         c(ia) = cswap
      endif
   end subroutine



end module
