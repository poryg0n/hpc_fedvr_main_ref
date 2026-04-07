      module observables
      use constants, only : ci
      use constants, only : ppi
      use constants, only : kapp
      use exploit_parameters
      use math_util
      use space_time_ops
      implicit none
      integer :: mode_k = 0                ! 1 for symmetric
      contains

!     subroutine compute_observables(...)
!        call spectral_obs(...)
!        call spatial_obs(...)
!     end subroutine
      

      subroutine compute_dyn_observables(nmax, psi, phi,           &
                                        xx, wx, jac,               &
                                        eigval, eigvec,            &
                                        norm_1, norm_2,            &
                                        p0, pexc, pion,            &
                                        dipole, momentum, energy)
      
        implicit none
        integer, intent(in) :: nmax
        complex(8), intent(in) :: psi(nmax), phi(nmax)
        real(8), intent(in) :: eigval(nmax), eigvec(nmax,nmax)
        real(8), intent(in) :: xx(nmax), wx(nmax), jac
      
        real(8), intent(out) :: norm_1, norm_2 
        real(8), intent(out) :: p0, pexc, pion
        complex(8), intent(out) :: energy, dipole, momentum
      
        complex(8) :: psi_x(nmax), p_psi(nmax)
        complex(8) :: phi_x(nmax), p_phi(nmax)
        integer :: i
      

        ! --- transform to DVR ---
        psi_x = matmul(eigvec, psi)
        phi_x = matmul(eigvec, phi)
        psi_x = psi_x/wx/dsqrt(jac)
        phi_x = phi_x/wx/dsqrt(jac)
        
        ! --- norm ---
        norm_1 = sum(abs(psi_x)**2 * (wx**2) * jac)
        norm_2 = sum(abs(phi_x)**2 * (wx**2) * jac)
        
        ! --- populations ---
        p0   = abs(psi(1))**2
        pexc = sum(abs(psi(2:nmax))**2)
        pion = 1.d0 - norm_1
        
        ! --- dipole ---
        dipole = sum(conjg(psi_x) * xx * psi_x * (wx**2) * jac)
        
        ! --- energy (field-free) ---
        energy = sum(abs(psi)**2 * eigval)
        
        ! --- momentum (CLEAN) ---
        call apply_momentum_operator(nmax, eigvec, xx, wx, jac,        &
                                                     psi, p_psi, 0)
        momentum = sum(conjg(psi) * p_psi)
      
      end subroutine


      subroutine compute_pemd_zrp(nmax, krange, t_end,           &
                           xx, wx, jacc,                               &
                           eigvec, eigval, wf0, wf,                 &
                           k_max, kk, p_ion, p0, a0, ak)
      
        implicit none
        integer, intent(in) :: nmax, krange
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        real(8), intent(out) :: p_ion, p0
        real(8), intent(out) :: kk(krange)
        complex(8), intent(in) :: wf0(nmax)
        complex(8), intent(in) :: wf(nmax)
        complex(8), intent(out) :: ak(krange)
        complex(8), intent(out) :: a0
      
        ! locals
        integer :: j, k, ij
        integer :: p, n_cont
        real(8) :: aux1, aux2, dk
        real(8) :: auxr1(krange/2), auxr2(krange/2)
        complex(8) :: auxc_
        complex(8) :: auxc(nmax), wfc_k(nmax)
        complex(8) :: wfc0(nmax)
        complex(8) :: wfc(nmax)
        complex(8), parameter :: ci = (0.d0,1.d0)
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

      
        p=0
        do p=1,nmax
           if (eigval(p).gt.0.0d0) then
              exit
           end if
        enddo
!       write(*,*) p

        ak = (0.d0, 0.d0)
        dk = k_max/(krange/2-1)

        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf, wfc)
      
        do j=1,krange
      
           if (j.le.(krange/2-p+1)) then
              ij = krange/2 + 1 - j
              kk(j)= -k_max + (j-1)*dk
      
           else if (j.ge.(krange/2+p)) then
              ij = j - krange/2
              kk(j) = (ij-1)*dk
      
           else
              kk(j)=0.d0
           end if

!          call build_wfc_k(xx, kk(j), kapp, mode_k, wfc_k)
      
           wfc_k = exp(ci*kk(j)*xx) +                                   &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                  &
                   exp(-ci*abs(kk(j)*xx))
      
           auxc = conjg(wfc_k) * wfc * wx*wx*jacc
           ak(j) = sum(auxc)
           ak(j) = exp(ci*(0.5d0*kk(j)**2)*t_end) * ak(j)
      
        end do

      
        ! --- probabilities ---
        p_ion = 0.d0
        call integr_over_range(krange, kk, ak, auxc_)
        p_ion = real(auxc_) / (2.d0*ppi)
        write(*,*) p_ion


        a0 = (0.d0, 0.d0)
        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0, wfc0)
        auxc = conjg(wfc0) * wfc * wx*wx*jacc
  
        a0 = sum(auxc)
        a0 = exp(ci*eigval(1)*t_end)*a0
      
        p0 = abs(a0)**2
      
      end subroutine



      subroutine compute_phi_elems(workdir,                  &
                           nmax, krange, t_end,              &
                           xx, wx, jacc,                             &
                           eigvec, eigval,                   &
                           wf0_1, wf_1, wf0_2, wf_2,                 &
                           omega, k_max, kk, a0, ak, b0w, bkw)
      
        implicit none
        integer, intent(in) :: nmax, krange
        real(8), intent(in) :: jacc, omega
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        complex(8), intent(in) :: a0

        character(255) :: workdir

        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: wf0_1(nmax)
        complex(8), intent(in) :: wf0_2(nmax)
        complex(8), intent(in) :: wf_1(nmax)
        complex(8), intent(in) :: wf_2(nmax)
        complex(8), intent(in) :: ak(krange)

        complex(8), intent(out) :: b0w
        complex(8), intent(out) :: bkw(krange)
      
        ! locals
        integer :: j, k, l, ij
        integer :: p, n_cont
        real(8) :: E0, Ek, Ekp
        real(8) :: aux1, aux2, dk, delta_kk
        real(8) :: auxr1(krange/2), auxr2(krange/2)
        complex(8) :: factor, denom
        complex(8) :: vec_0
        complex(8) :: b0wT, aux
        complex(8) :: wfc_k(nmax), wfc_k_(nmax)
        complex(8) :: dwfc_k(nmax), pwfc_k(nmax)
        complex(8) :: dwfc_0(nmax), pwfc_0(nmax)
        complex(8) :: wfc0_1(nmax)
        complex(8) :: wfc0_2(nmax)
        complex(8) :: wfc_1(nmax)
        complex(8) :: wfc_2(nmax)
        complex(8) :: auxc(nmax)
        complex(8) :: vec_1(krange)
        complex(8) :: vec_2(krange)
        complex(8) :: vec_k(krange)
        complex(8) :: bkwT(krange)
        complex(8), parameter :: ci = (0.d0,1.d0)
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

        real(8) :: eta = 1.d-6
        complex(8) :: dk0_, dkk_
        complex(8) :: pk0_, pkk_
        complex(8) :: pk0(krange), p0k(krange), pkk(krange)

        integer :: unit_pk0, unit_pkk
        open(newunit=unit_pk0, file=trim(workdir)//"/pk0.dat", status="replace")
        open(newunit=unit_pkk, file=trim(workdir)//"/pkk.dat", status="replace")
      
        p=0
        do p=1,nmax
           if (eigval(p).gt.0.0d0) then
              exit
           end if
        enddo
!       write(*,*) p

        dk = k_max/(krange/2-p+1)

        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0_1, wfc0_1)
        call differentiate(xx, wfc0_1, dwfc_0)
        pwfc_0 = -ci * dwfc_0

        E0  = 0.5d0 * kapp**2

        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf_1, wfc_1)
        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf_2, wfc_2)
     
!       kk = 0.d0 
!       do j=1,krange
!     
!          if (j.le.(krange/2-p+1)) then
!             ij = krange/2 + 1 - j
!             kk(j)= -k_max + (j-1)*dk
!     
!          else if (j.ge.(krange/2+p)) then
!             ij = j - krange/2
!             kk(j) = (ij-1)*dk
!     
!          end if
!       enddo

        do j=1,krange
!          call build_wfc_k(xx, kk(j), kapp, mode_k, wfc_k)
      
           wfc_k = exp(ci*kk(j)*xx) +                                   &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                  &
                   exp(-ci*abs(kk(j)*xx))
       
           Ek  = 0.5d0 * kk(j)**2 
 
           do l=1,krange
!             call build_wfc_k(xx, kk(l), kapp, mode_k, wfc_k_)
      
              wfc_k_ = exp(ci*kk(l)*xx) +                              &
                      (ci*kapp/(-abs(kk(l)) - ci*kapp)) *              &
                      exp(-ci*abs(kk(l)*xx))
     
              Ekp = 0.5d0 * kk(l)**2 

              call differentiate(xx, wfc_k_, dwfc_k)
              pwfc_k = -ci * dwfc_k

              auxc = conjg(wfc_k) * pwfc_k * wx*wx*jacc
              pkk(l) = sum(auxc)

              auxc = conjg(wfc_k) * xx * wfc_k_ * wx*wx*jacc
              dkk_ = sum(auxc)
              dkk_ = -ci* ( Ek - Ekp ) * dkk_

              ! *** analytical formula
              if(j.eq.l) then
                 delta_kk = 1.d0/dk
              else
                 delta_kk = 0.d0
              end if

              denom = (kk(j)**2 - kk(l)**2 + ci*eta)
              factor = ( kk(l)*abs(kk(j)) / (abs(kk(j)) -ci*kapp)      &
                      -  kk(j)*abs(kk(l)) / (abs(kk(l)) +ci*kapp)   )

              pkk_ = 2.d0*ppi*kk(j) * delta_kk                         &
                                   - 2.d0*kapp* factor/denom
     
              write(unit_pkk, *) kk(j), kk(l), pkk(l), dkk_, pkk_ 
              vec_2(l) = pkk(l) * ak(j) / ( Ek + omega - Ekp + ci*eta )
              vec_2(l) = exp(ci * ( Ek + omega - Ekp ) * t_end ) * vec_2(j)

           enddo

           call integr_over_range(krange, kk, vec_2, vec_k(j))

           auxc = conjg(wfc_k) * pwfc_0 * wx*wx*jacc
           pk0(j) = sum(auxc)
           p0k(j) = conjg(pk0(j))


           auxc = conjg(wfc_k) * xx * wfc0_1 * wx*wx*jacc
           dk0_ = sum(auxc)
           dk0_ = -ci* ( Ek - E0 ) * dk0_

           denom = ( kapp**2 + kk(j)**2 )

           pk0_ = 2.d0 * kk(j) * kapp**(3.d0/2) /denom
     
           write(unit_pk0, *) kk(j), pk0(j), dk0_, pk0_

           vec_1(j) = p0k(j) * ak(j) / ( E0 + omega - Ek + ci*eta )
           vec_1(j) = exp(ci * ( E0 + omega - Ek ) * t_end ) * vec_1(j)
        enddo

        call integr_over_range(krange, kk, vec_1, vec_0)

        auxc = conjg(wfc0_1) * wfc_2 * wx*wx*jacc
        b0wT = sum(auxc)
        b0wT = exp(ci*eigval(1)*t_end) * b0wT

        auxc = conjg(wfc_k) * wfc_2 * wx*wx*jacc
        bkwT = sum(auxc)
        bkwT = exp(ci*eigval(1)*t_end) * bkwT
       
        vec_1 = pk0 * a0 / ( Ek + omega - E0 + ci*eta )
        vec_1 = exp(ci * ( Ek + omega - E0 ) * t_end ) * aux

        b0w = b0wT + vec_0
        bkw = bkwT + vec_1 +  vec_k

        close(unit_pk0)
        close(unit_pkk)
      
      end subroutine


      subroutine compute_hhg_ip(nt, time, x_t, wsteps, omg,          &
                                               hhg_1, hhg_2)
      
        implicit none
        integer, intent(in) :: nt
        integer, intent(in) :: wsteps
        real(8), intent(in) :: time(nt)
        complex(8), intent(in) :: x_t(nt)
        real(8), allocatable, intent(out) :: omg(:)
        real(8), allocatable, intent(out) :: hhg_1(:)
        real(8), allocatable, intent(out) :: hhg_2(:)
      
        complex(8) :: auxc(nt), integral, amp
        complex(8), parameter :: ci = (0.d0, 1.d0)
        integer :: nw, i, it
        real :: w

        nw = wsteps+1
        allocate(omg(nw), hhg_1(nw), hhg_2(nw))

        do i=1,nw
           omg(i) = wmin + (i-1)*dw0
!          write(*,*) i, omg(i), wmin
        enddo

      
        do i = 1, nw
      
           ! --- integral term ---
           auxc = exp(ci*omg(i)*time) * x_t
           call composite_simpson_18c(nt, time, auxc, amp)
      
!          hhg_2(i) = omg(i)**2 * abs(amp)**2
           ! --- integration by parts ---
           amp = exp(ci*omg(i)*time(nt)) * x_t(nt)   &
               - exp(ci*omg(i)*time(1))  * x_t(1)    &
               - ci*omg(i) * amp
      
           ! --- spectrum ---
           hhg_1(i) = abs(amp)**2



           do it=1,nt
              w = sin(ppi*(time(it)-time(1))/(time(nt)-time(1)))**2
              auxc(it) = exp(ci*omg(i)*time(it))* x_t(it) * w
           enddo
           call composite_simpson_18c(nt, time, auxc, amp)
      
!          hhg_2(i) = abs(amp)**2
           hhg_2(i) = omg(i)**2 * abs(amp)**2
        end do
      
      end subroutine
      
      end module
