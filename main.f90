program main
!
!     ------------------------------------------------------------------
!     Main script for nonadiabatic dynamics.
!     Usage: polaron.x [input.in]
!     ------------------------------------------------------------------
!
   use maths, only : rnseed
   use types
   use system, only : model
   use tul, only : branching
   implicit real(dp) (a-h,o-z)
!
!  Initialise the random_number generator for reproducibility. 
!  To use wall clock, call rnseed(-1)
!
   call rnseed(42)
!
!  Read input and initialize module
!
   call init(ntraj,nt,dt)
!
!  Run what you want
!
   select case (model)
   case ('tul')
      ! Tully model
      if (branching) then
         call tully_branch (ntraj,dt)
      else
         call tully (ntraj,nt,dt)
      end if
   case default
      ! Population dynamics (default)
      call populations (ntraj,nt,dt)
   end select
!
!  Done!
!
end program

subroutine init(ntraj,nt,dt)
   use types
   use dynamics, only : initdyn
   use fmo, only : initfmo
   use sbm, only : initsbm
   use pyr, only : initpyr
   use dma, only : initdma
   use ful, only : initful
   use so2, only : initso2
   use rhe, only : initrhe
   use tul, only : inittul
   use system, only : init_system
   implicit real(dp) (a-h,o-z)
   character(1000) :: filename, buffer, label
   integer, parameter :: iu = 10
!
   character :: method
   character(3) :: model, nuc_opt
   character(1) :: hop_opt
   character(1) :: rescale_opt
   character(1) :: dec_opt
   logical :: branching
!
!  Constants
!
   fs = 41.34136d0
   rkB = 3.166829d-6
!
!  Set default inputs
!
   hop_opt = 'd'
   rescale_opt = 'd'
   dec_opt = 'n'
   dec_param = 0.0
   nuc_opt = 'wig'
   Delta = 1.0
   temp = 0.0
   n = 0
   branching = .false.
!
!  Open input file
!
   call getarg(1,filename)
   open(iu,file=filename,status='old',action='read',iostat=ios)
   if (ios.ne.0) stop 'Error opening input file'
!
!  Read input file line by line
!
   do while (.true.)
      read(iu, '(A)', iostat=ios) buffer
      if (ios.ne.0) exit
!
!     Find the first space in the line. Split label and values.
!
      ipos = scan(buffer, ' ')
      label = buffer(1:ipos)
      buffer = buffer(ipos+1:)
      buffer = TRIM(ADJUSTL(buffer))
      select case(label)
      case ("n")
         read (buffer,*,iostat=ios) n
      case ("ntraj")
         read (buffer,*,iostat=ios) ntraj
      case ("dt")
         read (buffer,*,iostat=ios) dt
      case ("nt")
         read (buffer,*,iostat=ios) nt
      case ("model")
         read (buffer,'(A)',iostat=ios) model
      case ("model_number")
         read (buffer,*,iostat=ios) model_number
      case ("branching")
         branching = .true.
      case ("method")
         read (buffer,*,iostat=ios) method
      case ("temp")
         read (buffer,*,iostat=ios) temp
      case ("frac")
         read (buffer,*,iostat=ios) frac
      case ("hop_opt","hopopt","hopping_opt","hopping_option")
         read (buffer,'(A)',iostat=ios) hop_opt
      case ("res_opt","resopt","rescale_opt","rescale_option")
         read (buffer,'(A)',iostat=ios) rescale_opt
      case ("dec_opt","decopt","decoherence_opt","decoherence_option")
         read (buffer,'(A)',iostat=ios) dec_opt
      case ("dec_param","decparam","decoherence_param","decoherence_parameter")
         read (buffer,*,iostat=ios) dec_param
      case ("nuc_opt","nucopt","sampnuc_opt")
         read (buffer,'(A)',iostat=ios) nuc_opt
      case ("eps","epsilon")
         read (buffer,*,iostat=ios) eps
      case ("Delta","delta")
         read (buffer,*,iostat=ios) Delta
      case ("lambda","lamda","rlamda","rlambda")
         read (buffer,*,iostat=ios) rlamda
      case ("omegac")
         read (buffer,*,iostat=ios) omegac
      case ("#", "!")
         continue
      case default
         print*, 'Unknown keyword: ', label
      end select
   end do

   if (temp.gt.0.0) then
      beta = 1.0/temp
   else 
      beta = 0.0
   end if
   close(iu)
!
!  Setup model
!
   select case (model)
   case('fmo')
   !
   !  FMO parameters
   ! 
      ns = 7
      nf = ns * n
      dt = dt * fs
      temp = temp * rkB
      beta = 1.0/temp
      call init_system(nf,ns,model,nuc_opt,beta)
      call initfmo()
   case ('sbm')
   !
   !  Spin-boson model
   !
      ns = 2
      nf = n
      call init_system(nf,ns,model,nuc_opt,beta)
      call initsbm(eps,Delta,rlamda,omegac)
   case ('pyr')
   !
   !  Pyrazine (3 or 24 modes)
   !
      ns = 2
      if (n.eq.0) n = 24
      nf = n
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt)
      call initpyr()
   case ('dma')
   !
   !  DMABN (57 modes)
   !
      ns = 3
      nf = 57
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt)
      call initdma()
   case ('ful')
   !
   !  Fulvene (30 modes)
   !
      ns = 2
      nf = 30
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt)
      call initful()
   case ('so2')
   !
   !  SO2 
   !
      model_number = n
      ns = 13
      nf = 3
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt)
      call initso2(model_number)   
   case ('rhe')
   !
   !  Rhenium(I) carbonyl alpha-diimine complex
   !
      ns = 14
      nf = 15
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt)
      call initrhe()
   case ('tul')
   !
   !  Tully model
   !
      ns = 2
      nf = 1
      dt = dt * fs
      call init_system(nf,ns,model,nuc_opt,beta)
      call inittul(model_number,branching)
   case default
      print*, 'Unknown model: ', model
      stop
   end select
!
!  Initialise dynamics module
!
   call initdyn(method,hop_opt,rescale_opt,dec_opt,dec_param)
!
!  Init done!
!
end subroutine


subroutine populations (ntraj,nt,dt)
!
!     ------------------------------------------------------------------
!     Calculate populations of diabatic states.
!     ------------------------------------------------------------------
!
   use types
   use dynamics, only : ftraj, mtraj, method, alphan
   use system, only : ns, nf, sampnuc, potmat, model
   implicit real(dp) (a-h,o-z)
   complex(dpc) ct,cat,ut
   allocatable :: pt(:,:),qt(:,:),ct(:,:),cat(:,:),iat(:)
   allocatable :: ut(:,:,:),vt(:,:),et(:),det(:)
   allocatable :: popd1(:,:),dpd1(:,:),popd2(:,:),dpd2(:,:),popa1(:,:),dpa1(:,:),popa2(:,:),dpa2(:,:)
   ! character(len=4) :: basis = 'site'
   allocate(pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt),iat(0:nt))
   allocate(ut(ns,ns,0:nt),vt(ns,0:nt),et(0:nt),det(0:nt))
   allocate(popd1(ns,0:nt),dpd1(ns,0:nt),popd2(ns,0:nt),dpd2(ns,0:nt),popa1(ns,0:nt),dpa1(ns,0:nt),popa2(ns,0:nt),dpa2(ns,0:nt))
   ! Setup initial density matrix
   init = 1
   if (model.eq.'pyr') init = 2
   if (model.eq.'so2') init = 2
   if (model.eq.'rhe') init = 2
   if (model.eq.'dma') init = 3
   if (model.eq.'ful') init = 2
   ! Exciton transformation
   ! qt(:,0) = 0.0
   ! call potmat (qt(:,0),ux,vx)
   ! if (basis.eq.'exci') then
   !    rho = real(matmul(conjg(ux),matmul(rho,transpose(ux))))
   ! end if
   ! Loop over trajectories
   popd1 = 0.0
   popd2 = 0.0
   popa1 = 0.0
   popa2 = 0.0
   Rt = 0.0
   St = 0.0
   et = 0.0
   m = 0
!$omp parallel do default(shared) private(pt,qt,ct,cat,iat,ut,vt,det, &
!$omp& dpd1,dpd2,dpa1,dpa2)
   do k = 1,ntraj
      ! Sample initial conditions
      call sampnuc(pt(:,0),qt(:,0)) ! Sample nuclear positions and momenta
      if (method.eq.'f') then
         !
         !   FSSH
         !
         call sample_adiabat(qt(:,0),ct(:,0),iat(0),init) 
         call ftraj(pt,qt,ct,cat,iat,ut,vt,det,nt,dt)
         call getpops(ut,cat,iat,dpd1,dpd2,dpa1,dpa2,nt)
      else if (method.eq.'m') then
         !
         !   Multi-MASH 
         !
         call sample_adiabat(qt(:,0),ct(:,0),iat(0),init) 
         call sample_mash(qt(:,0),ct(:,0),iat(0),init)
         call mtraj (pt,qt,ct,cat,iat,ut,vt,det,nt,dt)
         !
         !   Phi observable
         !
         ! if (basis.eq.'exci') then
         !    do it = 0,nt
         !       ct(:,it) = matmul(transpose(ux),ct(:,it))
         !    end do
         ! end if
         dpd1 = alphan*(abs(ct)**2  - 1.0/ns) + 1.0/ns
         dpa1 = alphan*(abs(cat)**2 - 1.0/ns) + 1.0/ns
      end if
!$omp critical
      m = m+1
      et = et + det
      popd1 = popd1 + dpd1
      popd2 = popd2 + dpd2
      popa1 = popa1 + dpa1
      popa2 = popa2 + dpa2
      if (ntraj.ge.20 .and. mod(m,ntraj/10).eq.0) then
         write (6,601) (100*m)/ntraj
601      format(1x,'Calculation',i4,'% done')
         call flush (6)
         call printdata(nt,dt,popd1/m,'popd1.out ',ns)
         call printdata(nt,dt,popa1/m,'popa1.out',ns)
         call printdata(nt,dt,et/m,'et.out   ',1)
         if (method.eq.'f') then
            call printdata(nt,dt,popd2/m,'popd2.out ',ns)
            call printdata(nt,dt,popa2/m,'popa2.out',ns)
         end if
      endif
      if (ntraj.eq.1) then
         ! Print iat, ut, vt as well
         open (unit=40,file='iat.out')
         open (unit=50,file='vt.out')
         open (unit=60,file='ut.out')
         do it = 0,nt
            write (40,*) it*dt,iat(it)
            write (50,*) it*dt,vt(:,it)
            write (60,*) it*dt,ut(:,:,it)
         end do
      end if   
!$omp end critical
   enddo
!$omp end parallel do
   et = et/ntraj
   popd1 = popd1/ntraj
   popd2 = popd2/ntraj
   popa1 = popa1/ntraj
   popa2 = popa2/ntraj
   call printdata(nt,dt,popd1,'popd1.out',ns)
   call printdata(nt,dt,popa1,'popa1.out',ns)
   call printdata(nt,dt,et,'et.out   ',1)
   if (method.eq.'f') then
      call printdata(nt,dt,popd2,'popd2.out',ns)
      call printdata(nt,dt,popa2,'popa2.out',ns)
   end if
   deallocate (pt,qt,ct,cat,iat,ut,vt,det)
   deallocate (et,popd1,dpd1,popd2,dpd2,popa1,dpa1,popa2,dpa2)
end subroutine

subroutine sample_adiabat(q,c,ia,init)
!
!        Adiabatic state sampling
!
   use types
   use system, only : ns, nf, potmat, model
   implicit real(dp) (a-h,o-z)
   complex(dpc) c,u
   dimension q(nf),c(ns)
   allocatable :: u(:,:),v(:),rho(:,:),Pa(:)
   allocate(u(ns,ns),v(ns),rho(ns,ns),Pa(ns))
   rho = 0.0
   rho(init,init) = 1.0
   call potmat (q,u,v) ! Adiabatic transformation
   do ib = 1,ns
      Pa(ib) = real(dot_product(u(:,ib),matmul(rho,u(:,ib))))
   end do
   ! Sample ia with probability according to Pa
   call random_number(r)
   Psum = 0.0
   do ib = 1,ns
      Psum = Psum + Pa(ib)
      if (r.lt.Psum) then
         ia = ib
         exit
      end if
   end do
   c = 0.0
   c(init) = 1.0
   if (model.eq.'ura') c = u(:,ia)
   deallocate(u,v,rho,Pa)
end subroutine

subroutine sample_mash(q,c,ia,init)
!
!        Diabatic cap sampling
!
   use types 
   use system, only : ns, nf, potmat
   use dynamics, only : sampcap, dec_opt, maxpop
   use maths, only : sampsun, pi
   implicit real(dp) (a-h,o-z)
   complex(dpc) c,u
   dimension :: q(nf),c(ns)
   allocatable :: u(:,:),v(:),list(:)
   allocate(u(ns,ns),v(ns),list(ns))
   call potmat(q,u,v) ! Adiabatic transformation
!
!     Diabatic sampling
!
   select case (dec_opt)
   case ('c','n')
      ! Cap
      call sampcap (init,c) ! Sample diabatic c vector
   case ('s')
      ! Sphere
      call sampsun (ns,c) ! Sample diabatic c vector
   end select
   call maxpop(u,c,ia) ! Find most populated adiabatic state
!
!     Adiabatic cap sampling
!
   ! call potmat (q,u,v) ! Adiabatic transformation
   ! rho = 0.0
   ! rho(init,init) = 1.0 ! Setup initial density matrix in diabatic basis
   ! do ia = 1,ns
   !    Pa(ia) = dreal(dot_product(u(:,ia),matmul(rho(:,:),u(:,ia)))) ! Initial adiabatic populations
   ! end do
   ! ! Sample ia with probability according to Pa
   ! call random_number(r)
   ! Psum = 0.0
   ! do ia = 1,ns
   !    Psum = Psum + Pa(ia)
   !    if (r.lt.Psum) then
   !       exit
   !    end if
   ! end do
   ! call sampcap (ia,c) ! Sample adiabatic c vector
   ! c = matmul(u,c) ! Transform to diabatic basis
!
   deallocate(u,v,list)
end subroutine

subroutine getpops (ut,cat,iat,dpd1,dpd2,dpa1,dpa2,nt)
!
!  Calculate population observables
!
   use types
   use system, only : ns
   implicit real(dp) (a-h,o-z)
   complex(dpc) :: ut(ns,ns,0:nt),cat(ns,0:nt)
   dimension :: iat(0:nt), dpd1(ns,0:nt), dpd2(ns,0:nt), dpa1(ns,0:nt), dpa2(ns,0:nt)
   complex(dpc), allocatable :: c(:),u(:,:)
   allocate(c(ns),u(ns,ns))
   do it = 0,nt
      ia = iat(it)
      u = ut(:,:,it)
      c = cat(:,it)
      ! if (basis.eq.'exci') then
         ! u = matmul(transpose(ux),u) ! Measure in exciton basis
      ! end if
!
!           Pi-style observables
!
      dpd1(:,it) = abs(u(:,ia))**2
      do ib=1,ns
         do ic=1,ns
            if (ib.eq.ic) cycle
            dpd1(:,it) = dpd1(:,it) + real(conjg(u(:,ic))*u(:,ib)  &
               * conjg(c(ic))*c(ib))
         end do
      end do
!
!           WF-style observables
!
      dpd2(:,it) = abs(matmul(u,c))**2
!
!           Adiabatic observables
!
      dpa1(:,it) = 0.0
      dpa1(ia,it) = 1.d0
      dpa2(:,it) = abs(c)**2
   end do
   deallocate(c,u)
end subroutine

subroutine tully_branch (ntraj,dt)
!
!  Branching ratios in Tully models starting from set of initial momenta
!
   use types
   use dynamics, only : method, ftraj
   use system, only : ns, nf, sampnuc, potmat, rmass, gam
   use tul, only : model_number
   implicit real(dp) (a-h,o-z)
   character(len=100) :: fmt
   complex(dpc) ct,cat,ut
   allocatable :: pt(:,:),qt(:,:),ct(:,:),cat(:,:),iat(:)
   allocatable :: ut(:,:,:),vt(:,:),et(:,:),det(:)
   allocatable :: pinit(:),refl1(:,:),refl2(:,:),dr1(:),dr2(:)
   allocatable :: tran1(:,:),tran2(:,:),dt1(:),dt2(:)
   if (model_number.eq.1) then
      pmin = 5.0
      pmax = 30.0
      np = 26
      init = 1
   else if (model_number.eq.2) then
      pmin = 10.0
      pmax = 50.0
      np = 40 ! Number of initial momenta
      init = 1
   else if (model_number.eq.3) then
      pmin = 3.0
      pmax = 30.0
      np = 28 ! Number of initial momenta
      init = 2 ! Lower adiabat = second diabat in this case
   end if
   allocate(pinit(np),refl1(np,ns),refl2(np,ns),dr1(ns),dr2(ns))
   allocate(tran1(np,ns),tran2(np,ns),dt1(ns),dt2(ns))
   ! Loop over trajectories
   refl1 = 0.0
   refl2 = 0.0
   tran1 = 0.0
   tran2 = 0.0
   ! Set up pinit grid
   do ip=1,np
      pinit(ip) = pmin + (ip-1.0)/(np-1.0) * (pmax-pmin)
   end do
   do ip=1,np
      nt = int(100.0/(dt*pinit(ip)/rmass(1)))
      gam = 2*(pinit(ip)/20.)**2
      print*, 'pinit', pinit(ip), nt
      allocate(pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt),iat(0:nt))
      allocate(ut(ns,ns,0:nt),vt(ns,0:nt),et(np,0:nt),det(0:nt))   
!omp parallel do default(shared) private(pt,qt,iat,ut,vt,ca,cat,det,dr1,dr2,dt1,dt2)
      do k=1,ntraj
         ! call sampnuc(pt(:,0),qt(:,0))
         ! qt(:,0) = qt(:,0) - 15.0
         ! pt(:,0) = pt(:,0) + pinit(ip)
         pt(:,0) = pinit(ip)
         qt(:,0) = -15.0
         ct(:,0) = 0.0
         ct(init,0) = 1.0
         iat(0) = 1 ! Start on lower adiabatic state
         if (method.eq.'f') then
            call ftraj (pt,qt,ct,cat,iat,ut,vt,det,nt,dt)
         end if
         dr1 = 0.0
         dt1 = 0.0
         dr2 = 0.0
         dt2 = 0.0
         if (qt(1,nt).lt.0.0) then
            dr1(iat(nt)) = 1.0
            dr2(:) = abs(cat(:,nt))**2
         else
            dt1(iat(nt)) = 1.0
            dt2(:) = abs(cat(:,nt))**2
         end if
!omp critical
         refl1(ip,:) = refl1(ip,:) + dr1(:)
         refl2(ip,:) = refl2(ip,:) + dr2(:)
         tran1(ip,:) = tran1(ip,:) + dt1(:)
         tran2(ip,:) = tran2(ip,:) + dt2(:)
         et(ip,:) = et(ip,:) + det
!omp end critical
      end do
!omp end parallel do
      deallocate(pt,qt,ct,cat,iat,ut,vt,et,det)
   end do
   refl1 = refl1/ntraj
   refl2 = refl2/ntraj
   tran1 = tran1/ntraj
   tran2 = tran2/ntraj
   ncol = 2*ns
   write(fmt, '(A, I0, A)') '(f10.3, 1x, ', ncol, 'f10.6)'
   open (unit=30,file='branch1.out')
   open (unit=31,file='branch2.out')
   do ip = 1,np
      write (30,fmt) pinit(ip),refl1(ip,:),tran1(ip,:)
      write (31,fmt) pinit(ip),refl2(ip,:),tran2(ip,:)
   end do
   close (unit=30)
   close (unit=31)
   deallocate(pinit,refl1,refl2,tran1,tran2,dr1,dr2,dt1,dt2)
end subroutine

subroutine tully (ntraj,nt,dt)
!
!  Run population and coherence dynamics for Tully models, starting from a wavepacket
!
   use types
   use dynamics, only : method, ftraj, mtraj, alphan
   use system, only : ns, nf, sampnuc, potmat
   use tul, only : model_number
   implicit real(dp) (a-h,o-z)
   complex(dpc) ct,cat,c,u,ut
   allocatable :: pt(:,:),qt(:,:),ct(:,:),cat(:,:),iat(:)
   allocatable :: ut(:,:,:),vt(:,:),wt(:),et(:),det(:)
   allocatable :: u(:,:),c(:)
   allocatable :: pop1(:,:),pop2(:,:),dp1(:,:),dp2(:,:),xy(:,:),dxy(:,:)
   allocate(pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt),iat(0:nt))
   allocate(ut(ns,ns,0:nt),vt(ns,0:nt),wt(0:nt),et(0:nt),det(0:nt))
   allocate(u(ns,ns),c(ns))
   allocate(pop1(ns,0:nt),pop2(ns,0:nt),dp1(ns,0:nt),dp2(ns,0:nt),xy(ns,0:nt),dxy(ns,0:nt))
   ! Loop over trajectories
   pop1 = 0.0
   pop2 = 0.0
   xy = 0.0
   et = 0.0
   pinit = 25.0
   init = 1
   if (model_number.eq.3) then
      init = 2
   end if
!$omp parallel do default(shared) private(pt,qt,iat,ut,vt,wt,ct,cat,ia,ib,ic,u,v,c,& 
!$omp& det,dp1,dp2,dxy)
   do k = 1,ntraj
      call sampnuc(pt(:,0),qt(:,0))
      ct(:,0) = 0.0
      ct(init,0) = 1.0
      iat(0) = 1
      if (method.eq.'f') then
!
!        FSSH trajectory
!
         call ftraj (pt,qt,ct,cat,iat,ut,vt,det,nt,dt)
!
!        Measure populations
!
         do it = 0,nt
            ia = iat(it)
            u = ut(:,:,it)
            c = cat(:,it)
!
!           Standard FSSH style observables
!
            dp1(:,it) = 0.0
            dp1(ia,it) = 1.0
            dp2(:,it) = 0.0
            dp2(:,it) = abs(c)**2
            dxy(1,it) = dreal(conjg(c(1))*c(2))
            dxy(2,it) = aimag(conjg(c(1))*c(2))
         end do
      else if (method.eq.'m') then
         call sample_mash(qt(:,0),ct(:,0),iat(0),init)
         call mtraj (pt,qt,ct,cat,iat,ut,vt,det,nt,dt)
         do it = 0,nt
            ia = iat(it)
            u = ut(:,:,it)
            c = cat(:,it)
!
!           Standard MASH style observables
!
            dp1(:,it) = 0.0
            dp1(ia,it) = 1.0
            dp2(:,it) = alphan*(abs(c)**2 - 1.0/ns) + 1.0/ns
            dxy(1,it) = alphan*dreal(conjg(c(1))*c(2))
            dxy(2,it) = alphan*aimag(conjg(c(1))*c(2))
         end do
      end if
!$omp critical
      et = et + det
      pop1 = pop1 + dp1
      pop2 = pop2 + dp2
      xy = xy + dxy
!$omp end critical
   end do
!$omp end parallel do

   et = et/ntraj
   pop1 = pop1/ntraj
   pop2 = pop2/ntraj
   xy = xy/ntraj
   call printdata(nt,dt,pop1(:,:),'popa1.out',ns)
   call printdata(nt,dt,et,'et.out   ',1)
   call printdata(nt,dt,pop2(:,:),'popa2.out',ns)
   call printdata(nt,dt,xy(:,:),'xy.out   ',2)
!
   deallocate (pt,qt,ct,cat,iat,ut,vt,et,det,u,c,pop1,pop2,dp1,dp2,xy,dxy)
end subroutine


subroutine maxpop(u,c,ia)
   use types
   use system, only : ns
   implicit real(dp) (a-h,o-z)
   complex(dpc) c(ns), u(ns,ns)
   integer :: ia
   real(dp) :: mpop
   mpop = 0.0
   do ib = 1,ns
      pop = abs(sum(u(:,ib)*c))**2
      if (pop>mpop) then
         mpop = pop
         ia = ib
      end if
   end do
end subroutine

subroutine printdata(nt,dt,data,fname,ncol)
   use types
   use system, only : model
   implicit real(dp) (a-h,o-z)
   dimension data(ncol,0:nt)
   character(len=100) :: fmt
   character(len=9) :: fname
   write(fmt, '(A, I0, A)') '(f10.3, 1x, ', ncol, 'f16.8)'
   fs = 41.34136d0
   if (model.eq.'sbm') then
      fs = 1.0
   end if
   open (unit=30,file=fname)
   do it = 0,nt
      t = it*dt/fs
      write (30,fmt) t,data(:,it)
   end do
   close (unit=30)
end subroutine

subroutine output(nt,dtin,xmsd,xmsd2,ct,et,ekt,pnt,method)
!   ------------------------------------------------------------------
!   Output MSD(t), C(t), D(t) and E(t) to file.
!   ------------------------------------------------------------------
   use types
   use system, only : model
   implicit real(dp) (a-h,o-z)
   dimension xmsd(0:nt),xmsd2(0:nt),ct(0:nt),et(0:nt),ekt(0:nt),pnt(0:nt)
   character :: method
   open (unit=10,file='Ct.out')
   open (unit=20,file='Dt.out')
   open (unit=30,file='Et.out')
   open (unit=40,file='pnt.out')
   dci = 0.0
   dt = dtin
   fs = 41.34136d0
   if (model.eq.'shi') dt = dtin / fs
   do it = 1,nt
      if (method.eq.'e'.or.method.eq.'c'.or.method.eq.'m') then
!
!        Numerically integrate C(t) to get D(t)
!
         dci = dci+(dt/2)*(ct(it-1)+ct(it))
      end if
!
!     Numerically differentiate MSD(t) to get 2*D(t)
!
      if (it.eq.nt) then
         dc = (xmsd(it)-xmsd(it-1))/dt
      else
         dc = (xmsd(it+1)-xmsd(it-1))/(2*dt)
      end if
      dc = 0.5*dc
      t = it*dt
      ! k = max(1,it/100)
      ! if (mod(it,k).eq.0 .and. it.lt.nt) then
         write (10,*) it*dt,xmsd(it),xmsd2(it),ct(it)
         write (20,*) it*dt,dc,xmsd2(it)/(2*t),dci
         write (30,*) it*dt,et(it),ekt(it)
         write (40,*) it*dt,pnt(it)
      ! end if
   enddo
   close (unit=10)
   close (unit=20)
   close (unit=30)
   close (unit=40)
end subroutine

subroutine runtrj(nt,dt)
!
!     ------------------------------------------------------------------
!     Calculate single trajectory and output to file.
!     ------------------------------------------------------------------
!
   use dynamics, only : mtraj, sampad, method
   use system, only : ns, nf, sampnuc, sampdis, potmat
   use types
   implicit real(dp) (a-h,o-z)
   complex(dpc) ct,cat,u,ut,uut
   allocatable :: pt(:,:),qt(:,:),ct(:,:),iat(:),ut(:,:),vt(:,:),et(:)
   allocatable :: xt(:), u(:,:), v(:), cat(:,:), v0t(:), pnt(:), uut(:,:,:)
   character(len=100) :: fmt
   allocate(pt(nf,0:nt),qt(nf,0:nt),ct(ns,0:nt),cat(ns,0:nt),pnt(0:nt))
   allocate(iat(0:nt),xt(0:nt),ut(ns,0:nt),vt(ns,0:nt),et(0:nt),v0t(0:nt))
   allocate(u(ns,ns),v(ns),uut(ns,ns,0:nt))
   call sampnuc (pt(:,0),qt(:,0))
   call potmat (qt(:,0),u,v)
   call sampad (v,iat(0))
   ct(:,0) = u(:,iat(0))
   if (method.eq.'m') then
      call mtraj (pt,qt,ct,cat,iat,uut,vt,et,nt,dt)
   ! TODO: fill in ftraj
   else
      stop 'Error: method not recognised'
   end if
!
!     Set output format for big arrays
!
   write(fmt, '(A, I0, A)') '(f10.3, 1x, ', ns, 'f8.4)'
!
!     Output to file
!
   !  General output
   open (unit=40,file='et.out')
   open (unit=50,file='vt.out')
   open (unit=80,file='xt.out')
   do it = 0,nt
      write (40,*) it*dt,et(it)
      ! v0 = 0.0
      ! nf_bath = nf/ns
      ! do i=1,nf
         ! j = mod(i-1,nf_bath)+1
         ! v0 = v0 + 0.5*omega(j)**2*qt(i,it)**2
      ! end do
      write (50,fmt) it*dt,vt(:,it) !+ v0 !+ 0.5*sum(omega**2*qt(:,it)**2)
      write (80,*) it*dt,xt(it)
   enddo
   close (unit=40)
   close (unit=50)
   close (unit=80)
   if (method.eq.'l') then
      ! Output for LZSH
      open (unit=30,file='pt.out')
      open (unit=40,file='qt.out')
      open (unit=50,file='v0t.out')
      open (unit=60,file='ut.out')
      open (unit=70,file='iat.out')
      open (unit=90,file='pnt.out')
      do it = 0,nt
         write (30,fmt) it*dt,pt(:,it)
         write (40,fmt) it*dt,qt(:,it)
         write (50,*) it*dt,v0t(it)
         write (60,fmt) it*dt,ut(:,it)
         write (70,*) it*dt,iat(it)
         write (90,*) it*dt,pnt(it)
      end do
      close (unit=30)
      close (unit=40)
      close (unit=50)
      close (unit=60)
      close (unit=70)
      close (unit=90)
   else
      ! Output for CPA, Ehrenfest and MASH
      open (unit=31,file='Pdt.out')
      do it = 0,nt
         write (31,*) it*dt,abs(ct(:,it))**2
      end do
      close (unit=31)
   end if
   if (method.eq.'m') then
      ! Output for MASH
      open (unit=32,file='Pat.out')
      open (unit=60,file='ut.out')
      open (unit=70,file='iat.out')
      do it = 0,nt
         write (32,*) it*dt,abs(cat(:,it))**2
         write (60,*) it*dt,real(uut(:,iat(it),it))
         write (70,*) it*dt,iat(it)
      end do
      close (unit=32)
      close (unit=60)
      close (unit=70)
   end if
   deallocate (pt,qt,ct,cat,iat,ut,vt,et,u,v,xt,uut)
end subroutine