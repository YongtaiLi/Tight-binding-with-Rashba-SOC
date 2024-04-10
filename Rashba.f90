program Rashba
  implicit none
  integer, parameter :: kind = 8 ! For double-precision

  real(kind), parameter :: pi = 3.1415926535897932384626433_kind
  real(kind), parameter :: e  = 2.7182818285904523536028747_kind
  real(kind), parameter :: zeror = 0.000000000000000000000000_kind
  real(kind), parameter :: oner = 1.000000000000000000000000_kind
  complex(kind), parameter :: ii=(0.0_kind,1.0_kind)
  complex(kind), parameter :: zeroc=(0.0_kind,0.0_kind)
  complex(kind), parameter :: onec=(1.0_kind,0.0_kind)
  
  real(kind) :: delta, eta, sp_BR
  integer :: nwn, kcx, kcy, kcz, kc  ! The grids
  integer :: testkx, testky ! The n-th grid(s) to print the Green's functions.
  integer :: ik, jk, kk, i, j, k, n ! loops and iterations
  integer :: ndim ! spatial dimensions

  real(kind), allocatable :: kx(:), ky(:), kz(:), omega(:) ! Momentum scales and energy scales.
  real(kind), allocatable :: Eps_1D(:), Eps(:,:), Eps_3D(:,:,:), alpha_so_real(:), alpha_so_imag(:), & ! Dispersion and Rashba
                             A_1D(:,:), A_up(:,:,:), A_dn(:,:,:), A_up_3D(:,:,:,:), A_dn_3D(:,:,:,:), & ! Spectral, 2D and 3D
                             DOS_up(:), DOS_dn(:), DOS_tot(:) ! DOS
  complex(kind), allocatable :: alpha_so(:,:), alpha_so_conjg(:,:), deltemp(:,:,:), deltemp_3D(:,:,:,:), & ! Energetics
                                Ginv_1D(:,:), G_1D(:,:), & ! Green's functions, 1D
                                Ginv(:,:,:,:,:), G(:,:,:,:,:), G_up(:,:,:), G_dn(:,:,:), & ! Green's functions, 2D
                                Ginv_3D(:,:,:,:,:,:), G_3D(:,:,:,:,:,:), G_up_3D(:,:,:,:), G_dn_3D(:,:,:,:) ! Green's functions, 3D
  
  complex(kind) :: pauli(2,2,3) ! Pauli matrix.
  real(kind) :: identity(2,2) !identity matrix
  
  ! For zgetri and zgetrf from LAPACK
  integer, parameter :: lwork = 100 ** 2
  integer :: info
  complex(kind), allocatable :: work(:)
  integer, allocatable :: ipvt(:)

  ! For files writing
  character(len=20) :: Inverse_Greens, Greens_up, Greens_dn, TotalDOS_2D, TotalDOS_3D, TotalDOS_1D
  integer :: io, io1, io2, io3, io4
  !****************************************************************************************************
  ! This practice program simulates the Green's functions, spectral functions, and the density of states of
  ! tight-binding electrons on a 2D square lattice with Rashba SO-coupling.
  !
  ! Here are the explanations of some variables:
  !
  ! delta        : Energy windows, so the energy varies from [-delta, +delta].
  ! kx(:), ky(:) : x- and y-components of momentum
  ! omega(:)     : Energy, omega \in [-delta, +delta].
  ! sp_BR        : Rashba SO-coupling strength.
  ! eta          : The broadening parameter that keeps the Green's functions analytic in the upper-half plane.
  !
  ! nwn             : Number of frequency points in the range (0, +delta], same for the negative half.
  ! kcx, kcy, kcz   : Number of momentum points (from [-pi, pi]).
  ! ik, jk, kk      : iteration indices for momentum-oriented loops.
  ! n               : iteration index for energy-oriented loops.
  !
  ! vareps(:,:)                  : dispersion of tight-binding, non-interacting electrons.
  ! alpha_so(:,:)                : k-dependent prefactor for Rashba (to be distinguished from "alp")
  ! G_up(:,:,:), G_dn(:,:,:)     : Up-spin Green's function (single-particle)
  ! G(:,:,:,:,:)                 : Green's function
  ! deltemp(:,:,:)               : The (omega + i eta - vareps) thing (only those three)
  ! sqmod_alpha_so(:,:)          : The square-modulus of alpha_so(:,:)
  ! Rash(:,:,:)                  : The term to add to a Green's function due to Rashba effect.
  !
  ! A_up(:,:,:), A_dn(:,:,:)     : The up-spin spectral function.
  ! DOS_up(:), DOS_dn(:,:,:)     : The density of states projected to up-spin electrons.
  ! 
  ! Here, we set t = 1, and we suppress t in this code.
  !****************************************************************************************************
  ! Specify some parameters. One can either make it user-friendly by asking them to specify the value as this below.
  print *, "Enter the spatial/reciprocal dimension, ndim."
  read(*,*) ndim
  print *, "Enter the energy window, delta."
  read(*,*) delta
  print *, "Enter the broadening parameter, eta."
  read(*,*) eta
  print *, "Enter the number of frequency points in positive half range, (0,delta], nwn"
  read(*,*) nwn
  print *, "So that the actual number of frequency points is (2*nwn+1), that's", (2*nwn+1)
  print *, "Enter the KGrid, kc."
  print *, "Note: For 3D, kc~20, for 2D, kc~100"
  read(*,*) kc
  print *, "Enter the strength of Rashba So-coupling, sp_BR."
  print *, "Note that sp_BR ranges from 0 to 1, most commonly around 0.2."
  read(*,*) sp_BR

  ! Or to set those values ourselves.
  kcx = kc; kcy = kc; kcz = kc

  print *, "In this program, the K-grids are, kc =", kcx
  print *, "The energy windows and grids are, delta = ", delta, ",and nwn =", nwn
  print *, "The strength of Rashba SOC is, sp_BR =", sp_BR
  print *, "And the broadening parameters are, eta =", eta
  print *, "The dimension is, ndim =", ndim
  if(ndim==2.or.ndim==3)then
    print *, "Valid dimensions, we continue."
  else if(ndim==1)then
    print *, "You specified the one-dimensional space simulation, where no Rashba SOC applies. Are you sure?"
    print *, "If you are, we can continue."
  else
    print *, "Invalid dimensions. STOP!"
    stop
  end if

  ! Making the momentum and energy axes.
  allocate(kx(kcx),ky(kcy),omega(-nwn:nwn))
  kx = zeror; ky = zeror; omega = zeror
  do ik = 1, kcx ! making kx(:)
    kx(ik) = -pi + 2*pi*(ik-1)/(kcx-1)
  end do
  do jk = 1, kcy ! making ky(:)
    ky(jk) = -pi + 2*pi*(jk-1)/(kcy-1)
  end do
  if(ndim==3)then
    allocate(kz(kcz))
    kz = zeror
    do kk = 1, kcz ! making kz(:)
      kz(kk) = -pi + 2*pi*(kk-1)/(kcz-1)
    end do
  end if
  do n = -nwn, nwn ! making omega(:)
      omega(n) = delta*(n)/nwn
  end do
  
  ! Making the identity and Pauli matrices, pauli(2,2,3)
  ! The 1st coord. is the row number, the 2nd coord. is the column number, the 3rd coord is x, y, z
  identity(1,:)=[oner, zeror]; identity(2,:)=[zeror, oner]! identity matrix
  pauli(1,:,1)=[zeror, oner];  pauli(2,:,1)=[oner, zeror] ! sigma_x
  pauli(1,:,2)=[zeroc, -ii];    pauli(2,:,2)=[ii, zeroc]  ! sigma_y
  pauli(1,:,3)=[oner, zeror];  pauli(2,:,3)=[zeror, -oner]! sigma_z

  ! The kinetic energy of tight-binding electrons on a 2D square lattice, Eps(:,:)
  if(ndim==2)then
    allocate(Eps(kcx,kcy))
    Eps = zeror
    do ik = 1, kcx
      do jk = 1, kcy
        Eps(ik,jk) = -2*(cos(kx(ik))+cos(ky(jk)))
      end do
    end do
  else if(ndim==3)then
    allocate(Eps_3D(kcx,kcy,kcz))
    Eps_3D = zeror
    do ik = 1, kcx
      do jk = 1, kcy
        do kk = 1, kcz
          Eps_3D(ik,jk,kk) = -2*(cos(kx(ik))+cos(ky(jk))+cos(kz(kk)))
        end do
      end do
    end do
  else if(ndim==1)then
    allocate(Eps_1D(kcx))
    do ik = 1, kcx
      Eps_1D(ik) = -2*cos(kx(ik))
    end do
  end if

  ! The (omega + i eta - vareps) stuff, deltemp(:,:,:) (Work for 2D and 3D. If 1D, then it's trivially the inverse G)
  if(ndim==2)then
    allocate(deltemp(kcx,kcy,-nwn:nwn))
    deltemp = zeroc
    do n = -nwn,nwn
      do ik = 1, kcx
        do jk = 1, kcy
          deltemp(ik,jk,n) = dcmplx(omega(n)-Eps(ik,jk), eta)
        end do
      end do
    end do
  else if(ndim==3)then
    allocate(deltemp_3D(kcx,kcy,kcz,-nwn:nwn))
    deltemp_3D = zeroc
    do n = -nwn,nwn
      do ik = 1, kcx
        do jk = 1, kcy
          do kk = 1, kcz
            deltemp_3D(ik,jk,kk,n) = dcmplx(omega(n)-Eps_3D(ik,jk,kk), eta)
          end do
        end do
      end do
    end do
  else
    continue
  end if

  ! The real and imaginary part of the Rashba coef., alpha_so_real(:,:), alpha_so_imag(:,:)
  if(ndim==2.or.ndim==3)then
    allocate(alpha_so_real(kcy),alpha_so_imag(kcx))
    alpha_so_real=zeror; alpha_so_imag=zeror

    do jk = 1, kcy
      alpha_so_real(jk) = -sp_BR*sin(ky(jk))
    end do
    do ik = i, kcx
      alpha_so_imag(ik) = -sp_BR*sin(kx(ik))
    end do
  else
    continue
  end if

  ! Assemble the inverse Green's function using Pauli matrices (Note, if 1D, then no Pauli matrices are needed)
  if(ndim==2)then
    allocate(Ginv(2,2,kcx,kcy,-nwn:nwn))
    Ginv=zeroc
    do n = -nwn,nwn 
      do ik = 1,kcx
        do jk = 1, kcy
              Ginv(:,:,ik,jk,n) = deltemp(ik,jk,n)*identity(:,:) - alpha_so_real(jk)*pauli(:,:,1) & 
                                  + alpha_so_imag(ik)*pauli(:,:,2)
        end do
      end do
    end do

  else if(ndim==3)then
    allocate(Ginv_3D(2,2,kcx,kcy,kcz,-nwn:nwn))
    Ginv_3D=zeroc
    do n = -nwn,nwn
      do ik = 1,kcx
        do jk = 1, kcy
          do kk = 1, kcz
            Ginv_3D(:,:,ik,jk,kk,n) = deltemp_3D(ik,jk,kk,n)*identity(:,:) - alpha_so_real(jk)*pauli(:,:,1) &
                                      + alpha_so_imag(ik)*pauli(:,:,2)
          end do
        end do
      end do
    end do

  else if(ndim==1)then
    allocate(Ginv_1D(kcx,-nwn:nwn))
    Ginv_1D = zeroc
    do n = -nwn,nwn
      do ik = 1, kcx
        Ginv_1D(ik,n) = dcmplx(omega(n)-Eps_1D(ik), eta)
      end do
    end do
  end if
   
  ! Matrix inversion to obtain G (For 1D, trivial inversion)
  !#########################################################################
  ! When compiling with the library LAPACK, type
  ! (...)$ [compiler] [source_file_name].f90 -llapack -lblas
  !#########################################################################
  if(ndim==2)then
    allocate(G(2,2,kcx,kcy,-nwn:nwn))
    allocate(ipvt(2), work(lwork))
    G = Ginv
    do n = -nwn,nwn
      do jk = 1,kcy
        do ik = 1,kcx
          call zgetrf(2,2,G(1,1,ik,jk,n),2,ipvt,info)
          if (info==0.d0) then
            call zgetri(2,G(1,1,ik,jk,n),2,ipvt,work,lwork,info)
          else
            print *, "info \=0. Matrix inversion failed for Ginv -> G: Info=", info
            stop
          end if
        end do
      end do
    end do

  else if(ndim==3)then
    allocate(G_3D(2,2,kcx,kcy,kcz,-nwn:nwn))
    allocate(ipvt(2), work(lwork))
    G_3D = Ginv_3D
    do n = -nwn,nwn
      do kk = 1, kcz
        do jk = 1,kcy
          do ik = 1,kcx
            call zgetrf(2,2,G_3D(1,1,ik,jk,kk,n),2,ipvt,info)
            if (info==0.d0) then
              call zgetri(2,G_3D(1,1,ik,jk,kk,n),2,ipvt,work,lwork,info)
            else
              print *, "info \=0. Matrix inversion failed for Ginv_3D -> G_3D: Info=", info
              stop
            end if
          end do
        end do
      end do
    end do

  else ! For 1D, just trivial one-over-Ginv flipping. 
    allocate(G_1D(kcx,-nwn:nwn))
    do n = -nwn,nwn
      do ik = 1, kcx
        G_1D(ik,n) = 1/Ginv_1D(ik, n)
      end do
    end do
  end if
  
  ! Distribute the G to G_up and G_dn. Note, this only happens in 3D and 2D.
  if(ndim==2)then
    allocate(G_up(kcx,kcy,-nwn:nwn),G_dn(kcx,kcy,-nwn:nwn))
    do n = -nwn,nwn
      do jk = 1,kcy
        do ik = 1,kcx
          G_up(ik,jk,n) = G(1,1,ik,jk,n)
          G_dn(ik,jk,n) = G(2,2,ik,jk,n)
        end do
      end do
    end do
  
  else if(ndim==3)then
    allocate(G_up_3D(kcx,kcy,kcz,-nwn:nwn),G_dn_3D(kcx,kcy,kcz,-nwn:nwn))
    do n = -nwn,nwn
      do kk = 1, kcz
        do jk = 1,kcy
          do ik = 1,kcx
            G_up_3D(ik,jk,kk,n) = G_3D(1,1,ik,jk,kk,n)
            G_dn_3D(ik,jk,kk,n) = G_3D(2,2,ik,jk,kk,n)
          end do
        end do
      end do
    end do
  else
    continue
  end if
  
  ! Form the spectral functions, A_up and A_dn for 2D and 3D, and only A_1D for 1D.
  if(ndim==2)then
    allocate(A_up(kcx,kcy,-nwn:nwn),A_dn(kcx,kcy,-nwn:nwn))
    A_up=zeror; A_dn=zeror
    do n = -nwn,nwn
      do jk = 1,kcy
        do ik = 1,kcx
          A_up(ik,jk,n) = - dimag(G_up(ik,jk,n))/pi
          A_dn(ik,jk,n) = - dimag(G_dn(ik,jk,n))/pi
        end do
      end do
    end do
  else if(ndim==3)then
    allocate(A_up_3D(kcx,kcy,kcz,-nwn:nwn),A_dn_3D(kcx,kcy,kcz,-nwn:nwn))
    A_up_3D=zeror; A_dn_3D=zeror
    do n = -nwn,nwn
      do kk = 1, kcz
        do jk = 1, kcy
          do ik = 1, kcx
            A_up_3D(ik,jk,kk,n) = -dimag(G_up_3D(ik,jk,kk,n))/pi
            A_dn_3D(ik,jk,kk,n) = -dimag(G_dn_3D(ik,jk,kk,n))/pi
          end do
        end do
      end do
    end do
  else
    allocate(A_1D(kcx,-nwn:nwn))
    A_1D=zeror
    do n= -nwn,nwn
      do ik = 1, kcx
        A_1D(ik,n) = -dimag(G_1D(ik,n))/pi
      end do
    end do
  end if
  
  ! Write the density of states (DOS). 
  if(ndim==2)then
    allocate(DOS_up(-nwn:nwn),DOS_dn(-nwn:nwn),DOS_tot(-nwn:nwn))
    DOS_up = zeror; DOS_dn = zeror; DOS_tot = zeror
    do n = -nwn,nwn
      do ik = 1,kcx
        do jk = 1,kcy
          DOS_up(n) = DOS_up(n) + A_up(ik,jk,n)
          DOS_dn(n) = DOS_dn(n) + A_dn(ik,jk,n)
        end do
      end do
    end do
    
    ! We scale the density of states by 2*kcx*kcy.
    DOS_up(:) = DOS_up(:)/dble(2*kcx*kcy); DOS_dn(:) = DOS_dn(:)/dble(2*kcx*kcy)
    DOS_tot(:) = DOS_up(:) + DOS_dn(:)

    ! Write the (total) density of states.
    print *, "Start writing the total density of states for 2D system..."
    TotalDOS_2D = "DOS_tot.dat"
    open(file=trim(TotalDOS_2D), newunit=io)
    do n = -nwn,nwn
      write(io,*) omega(n), DOS_tot(n)
    end do
    close(io)
  else if(ndim==3)then
    allocate(DOS_up(-nwn:nwn),DOS_dn(-nwn:nwn),DOS_tot(-nwn:nwn))
    DOS_up = zeror; DOS_dn = zeror; DOS_tot = zeror
    do n = -nwn,nwn
      do ik = 1,kcx
        do jk = 1,kcy
          do kk = 1,kcz
             DOS_up(n) = DOS_up(n) + A_up_3D(ik,jk,kk,n)
             DOS_dn(n) = DOS_dn(n) + A_dn_3D(ik,jk,kk,n)
          end do
        end do
      end do
    end do
    
    ! We scale the density of states by 2*kcx*kcy*kcz.
    DOS_up(:) = DOS_up(:)/dble(2*kcx*kcy*kcz); DOS_dn(:) = DOS_dn(:)/dble(2*kcx*kcy*kcz)
    DOS_tot(:) = DOS_up(:) + DOS_dn(:)

    ! Write the (total) density of states.
    print *, "Start writing the total density of states for 3D system..."
    TotalDOS_3D = "DOS_tot.dat"
    open(file=trim(TotalDOS_3D), newunit=io)
    do n = -nwn,nwn
      write(io,*) omega(n), DOS_tot(n)
    end do
    close(io)
  else
    allocate(DOS_tot(-nwn:nwn))
    DOS_tot=zeror
    do n = -nwn,nwn
      do ik = 1, kcx
        DOS_tot(n) = DOS_tot(n) + A_1D(ik,n)
      end do
    end do

    ! We scale the density of states by kcx.
    DOS_tot(:) = DOS_tot(:)/dble(kcx)
    
    ! Write the DOS
    print *, "Start writing the total density of states for 3D system..."
    TotalDOS_1D = "DOS_tot.dat"
    open(file=trim(TotalDOS_1D), newunit=io)
    do n = -nwn,nwn
      write(io,*) omega(n), DOS_tot(n)
    end do
    close(io)
  end if

  print *, "The total density of states has been written in DOS_tot.dat"
  print *, "The first column is the frquency, and the second column is the DOS_tot."
  if(ndim==2) print *, "Note that the total DOS is normalized by 2*kcx*kcy =", 2*kcx*kcy
  if(ndim==3) print *, "note that the total DOS is normalized by 2*kcx*kcy*kcz =", 2*kcx*kcy*kcz
  if(ndim==1) print *, "Note that the total DOS is normalized by kcx =", kcx
end program Rashba

