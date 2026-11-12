program test
! e-(k1) e+(k2) -> phi+(k3) phi-(k4) in scalar QED (model sqed_ct).
!
! Checks at two phase-space points (sqrt(s) = 500 GeV, all masses kept):
!  1. LO against the analytic Born of the s-channel photon;
!  2. the double pole vanishes (all charged particles are massive);
!  3. the single pole of the renormalized amplitude equals the soft-photon
!     pole of massive charged legs, -sum_{j,k} q_j q_k A_jk (units of
!     alpha/2pi, relative to LO): the UV poles cancel against the
!     counterterms;
!  4. the finite part depends on mu^2 as the single pole does;
!  5. ir_subtraction (QED insertion operator, with the masses of the
!     initial-state electrons) gives the same poles.
! The finite part is printed but not compared: there is no independent
! reference for it.
use eepp_config, only: ki, nlo_prefactors
use eepp_matrix, only: initgolem, exitgolem, samplitude, ir_subtraction
use eepp_model, only: me, mW, e
implicit none

integer, parameter :: logf = 27
real(ki), parameter :: pi = 3.14159265358979323846264338327950288_ki
real(ki), parameter :: sqrts = 500.0_ki
! tolerances (observed agreement: LO ~1e-15, single pole ~1e-13)
real(ki), parameter :: eps_lo = 1.0e-10_ki, eps_pole = 1.0e-9_ki
real(ki), parameter :: eps_dpole = 1.0e-10_ki, eps_mu = 1.0e-8_ki

real(ki), dimension(3, 2) :: dirs
real(ki), dimension(4, 4) :: vecs
real(ki), dimension(0:3) :: amp, amp4
real(ki), dimension(2:3) :: irp
real(ki) :: scale2, born, pole, dmu, d
integer :: ipt, prec, ch
integer, dimension(2) :: channels
logical :: success

! directions of phi+ (normalized below)
dirs(:,1) = (/ -35.903073612782187_ki, 188.40993425041188_ki, 138.75289043803963_ki /)
dirs(:,2) = (/ 145.99739685097936_ki, 74.345616944054157_ki, 170.87184665106653_ki /)

channels = (/ logf, 6 /)
open(file="test.log", unit=logf)
success = .true.

call initgolem()
nlo_prefactors = 0

do ipt = 1, 2
   call kinematics(dirs(:,ipt), vecs)
   scale2 = sqrts**2
   call samplitude(vecs, scale2, amp, prec)
   call samplitude(vecs, 4.0_ki*scale2, amp4, prec)
   call ir_subtraction(vecs, scale2, irp)

   born = born_analytic(vecs)
   pole = soft_pole(vecs)
   dmu = (amp4(1) - amp(1)) / amp(0) / log(4.0_ki)

   do ch = 1, 2
      write(channels(ch),'(A,I1)') "phase-space point ", ipt
      write(channels(ch),'(A,4(1x,ES23.15))') "  k3 =", vecs(3,:)
      write(channels(ch),'(A,2(1x,ES23.15))') "  LO           GoSam, analytic:", amp(0), born
      write(channels(ch),'(A,1x,ES23.15)')    "  finite part  GoSam:          ", amp(1)/amp(0)
      write(channels(ch),'(A,2(1x,ES23.15))') "  single pole  GoSam, soft:    ", amp(2)/amp(0), pole
      write(channels(ch),'(A,1x,ES23.15)')    "  double pole  GoSam:          ", amp(3)/amp(0)
      write(channels(ch),'(A,2(1x,ES23.15))') "  ir_subtraction single, double:", irp(2)/amp(0), irp(3)/amp(0)
      write(channels(ch),'(A,1x,ES23.15)')    "  d(finite)/d(ln mu^2):        ", dmu
   end do

   d = abs(amp(0) - born) / abs(born)
   call check(d .le. eps_lo, "LO", d)
   d = abs(amp(3) / amp(0))
   call check(d .le. eps_dpole, "double pole", d)
   d = abs(amp(2)/amp(0) - pole) / abs(pole)
   call check(d .le. eps_pole, "single pole", d)
   d = abs(dmu - amp(2)/amp(0)) / abs(amp(2)/amp(0))
   call check(d .le. eps_mu, "scale dependence", d)
   d = abs(irp(2) - amp(2)) / abs(amp(2))
   call check(d .le. eps_pole, "ir_subtraction, single pole", d)
   d = abs(irp(3) / amp(0))
   call check(d .le. eps_dpole, "ir_subtraction, double pole", d)
end do

call exitgolem()

if (success) then
   write(unit=logf,fmt="(A15)") "@@@ SUCCESS @@@"
   write(unit=6,fmt="(A15)") "@@@ SUCCESS @@@"
else
   write(unit=logf,fmt="(A15)") "@@@ FAILURE @@@"
   write(unit=6,fmt="(A15)") "@@@ FAILURE @@@"
end if
close(unit=logf)

contains

subroutine check(ok, what, d)
   logical, intent(in) :: ok
   character(len=*), intent(in) :: what
   real(ki), intent(in) :: d
   if (.not. ok) then
      success = .false.
      do ch = 1, 2
         write(channels(ch),'(A,A,A,ES10.3)') "==> comparison of ", what, &
            & " failed, relative difference ", d
      end do
   end if
end subroutine check

! e- along +z, e+ along -z, phi+ along dir, phi- opposite, in the cms
subroutine kinematics(dir, vecs)
   real(ki), dimension(3), intent(in) :: dir
   real(ki), dimension(4, 4), intent(out) :: vecs
   real(ki) :: E, pe, pf
   E = sqrts / 2.0_ki
   pe = sqrt(E**2 - me**2)
   pf = sqrt(E**2 - mW**2)
   vecs(1,:) = (/ E, 0.0_ki, 0.0_ki,  pe /)
   vecs(2,:) = (/ E, 0.0_ki, 0.0_ki, -pe /)
   vecs(3,1) = E
   vecs(3,2:4) = pf * dir / sqrt(sum(dir**2))
   vecs(4,1) = E
   vecs(4,2:4) = -vecs(3,2:4)
end subroutine kinematics

pure function dot(p, q)
   real(ki), dimension(4), intent(in) :: p, q
   real(ki) :: dot
   dot = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
end function dot

! spin-averaged |M|^2 = e^4/s^2 [ 2 (p1.k)(p2.k) - k^2 (p1.p2 + me^2) ],
! k = p3 - p4
function born_analytic(vecs) result(m2)
   real(ki), dimension(4, 4), intent(in) :: vecs
   real(ki) :: m2, s
   real(ki), dimension(4) :: k
   k = vecs(3,:) - vecs(4,:)
   s = dot(vecs(1,:) + vecs(2,:), vecs(1,:) + vecs(2,:))
   m2 = e**4 / s**2 * (2.0_ki*dot(vecs(1,:), k)*dot(vecs(2,:), k) &
      & - dot(k, k)*(dot(vecs(1,:), vecs(2,:)) + me**2))
end function born_analytic

! 1/eps pole of 2 Re(M0^* M1)/|M0|^2 in units of alpha/(2 pi) from soft
! photons between massive charged legs: -sum_{j,k} q_j q_k A_jk with
! q = +Q (outgoing), -Q (incoming), A_jj = 1,
! A_jk = ln((1+v)/(1-v))/(2v), v = sqrt(1 - m_j^2 m_k^2/(p_j.p_k)^2)
function soft_pole(vecs) result(pole)
   real(ki), dimension(4, 4), intent(in) :: vecs
   real(ki) :: pole, r, v, a
   real(ki), dimension(4), parameter :: q = (/ 1.0_ki, -1.0_ki, 1.0_ki, -1.0_ki /)
   real(ki), dimension(4) :: m2
   integer :: j, k
   m2 = (/ me**2, me**2, mW**2, mW**2 /)
   pole = 0.0_ki
   do j = 1, 4
      do k = 1, 4
         if (j .eq. k) then
            a = 1.0_ki
         else
            ! masses from the model, not p.p: for the electron
            ! E^2 - |p|^2 would lose all digits of m_e^2
            r = m2(j)*m2(k) / dot(vecs(j,:), vecs(k,:))**2
            v = sqrt(1.0_ki - r)
            ! 1 - v = r/(1 + v), stable for v -> 1
            a = log((1.0_ki + v)**2 / r) / (2.0_ki*v)
         end if
         pole = pole - q(j)*q(k)*a
      end do
   end do
end function soft_pole

end program test
