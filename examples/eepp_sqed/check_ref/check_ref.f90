program check_ref
! Compare GoSam (model sqed_ct) with reference values for
!   e-(p1) e+(p2) -> pi-(p3) pi+(p4) [ + A(k) ]
! computed independently with a photon mass lambda^2 = 1e-30 GeV^2.
!
! GoSam ordering: e-, e+, phi+, phi-, [A]  -> pi+ = reference p4, pi- = p3.
! IR mapping (all charged legs massive, no double pole):
!   photon-mass result = GoSam finite part at mu^2 = lambda^2,
! i.e. 1/eps_IR -> ln(lambda^2/mu^2). The on-shell UV renormalisation
! (alpha(0)) makes the result independent of the UV scale.
! nlo_prefactors = 2: amp(1) includes alpha/(2 pi), as the reference.
! Conventions of the reference found by this comparison:
!  - spin SUM over e-, e+ (GoSam averages): factor 4;
!  - no vacuum polarisation in the s-channel photon. GoSam includes the
!    loops of e, mu, tau, phi (alpha(0) scheme); 2 Re Pi_ren(s) * Born
!    is subtracted analytically (closed formulas below).
! For 2->3 the subtraction is NOT applied: the photon virtuality differs
! between diagrams there (see the printed lines).
use PROC_config, only: ki, nlo_prefactors
use PROC_matrix, only: initgolem, exitgolem, samplitude
use PROC_model, only: alpha, me, mmu, mtau, mW
#ifdef NOVP
! process generated with PROC_novp.rc (vacpol=false, in virtual_novp/)
use PROC_ew_ct, only: vacpol
#endif
implicit none
integer, parameter :: N = NLEG
real(ki), parameter :: lambda2 = 1.0e-30_ki
real(ki), dimension(N, 4) :: ref, vecs
real(ki), dimension(0:3) :: amp, ampS
real(ki) :: mpi2, s, born_ref, virt_ref, fin_shift, vp
integer :: prec

! ---- reference momenta (E, px, py, pz) ----
#include "momenta.inc"

! pion mass from the momenta (E^2 - |p|^2 is safe for the pion)
mpi2 = ref(3,1)**2 - sum(ref(3,2:4)**2)

alpha = 0.72973525694968E-02_ki
me    = 0.510998950e-3_ki
mmu   = 105.658375e-3_ki
mW    = sqrt(mpi2)          ! phi mass

vecs = ref
vecs(3,:) = ref(4,:)   ! phi+ = pi+
vecs(4,:) = ref(3,:)   ! phi- = pi-  (5th leg, if any: the photon)

#ifdef NOVP
if (vacpol) stop "virtual_novp was not generated with vacpol=false"
write(*,'(A)') " variant: NO vacuum polarisation (virtual_novp, vacpol=false)"
#endif
call initgolem()        ! init_functions + ew_ren with the values above
nlo_prefactors = 2

s = 2.0_ki*(ref(1,1)*ref(2,1) - sum(ref(1,2:4)*ref(2,2:4))) + 2.0_ki*me**2
call samplitude(vecs, lambda2, amp, prec)
call samplitude(vecs, s, ampS, prec)
! same thing via mu^2 = s plus the pole times ln(lambda^2/s)
fin_shift = ampS(1) + ampS(2)*log(lambda2/s)

write(*,'(A,ES23.15)') " m_pi    = ", sqrt(mpi2)
write(*,'(A,ES23.15)') " m_tau   = ", mtau
write(*,'(A,I3)')      " prec    = ", prec
write(*,'(A,2ES23.15,ES11.2)') " BORN     GoSam, ref, rel.diff: ", amp(0), born_ref, amp(0)/born_ref-1
write(*,'(A,2ES23.15,ES11.2)') " VIRTUAL  GoSam, ref, rel.diff: ", amp(1), virt_ref, amp(1)/virt_ref-1
#ifdef NOVP
vp = 0.0_ki
if (.false.) then
#else
if (N == 4) then
#endif
vp = 2.0_ki*(pif(me) + pif(mmu) + pif(mtau) + pis(mW))
write(*,'(A,ES23.15)')         " 2 Re Pi(s), e+mu+tau+pi:      ", vp
write(*,'(A,2ES23.15,ES11.2)') " 4*BORN   GoSam, ref, rel.diff: ", 4*amp(0), born_ref, 4*amp(0)/born_ref-1
write(*,'(A,2ES23.15,ES11.2)') " 4*VIRT   GoSam, ref, rel.diff: ", 4*amp(1), virt_ref, 4*amp(1)/virt_ref-1
write(*,'(A,2ES23.15,ES11.2)') " 4*(VIRT-VP*BORN)  vs ref:      ", 4*(amp(1)-vp*amp(0)), virt_ref, &
   & 4*(amp(1)-vp*amp(0))/virt_ref-1
else
write(*,'(A,2ES23.15,ES11.2)') " 4*BORN   GoSam, ref, rel.diff: ", 4*amp(0), born_ref, 4*amp(0)/born_ref-1
write(*,'(A,2ES23.15,ES11.2)') " 4*VIRT   GoSam, ref, rel.diff: ", 4*amp(1), virt_ref, 4*amp(1)/virt_ref-1
end if
write(*,'(A,ES23.15)')         " VIRTUAL  via mu^2=s + pole*ln(lambda^2/s): ", fin_shift
write(*,'(A,2ES23.15)')        " poles (1/eps, 1/eps^2) at mu^2=lambda^2: ", amp(2), amp(3)
call exitgolem()

contains

! Re of the on-shell renormalised photon self-energy, Pi(s) - Pi(0),
! for s > 4 m^2, sign such that alpha_eff = alpha/(1 - Pi):
! Dirac fermion of charge 1 and complex scalar of charge 1
function pif(m) result(p)
   real(ki), intent(in) :: m
   real(ki) :: p, r, b
   real(ki), parameter :: pi = 3.14159265358979323846264338327950288_ki
   r = m**2/s
   b = sqrt(1.0_ki - 4.0_ki*r)
   p = alpha/(3.0_ki*pi)*(-5.0_ki/3.0_ki - 4.0_ki*r &
     & + b*(1.0_ki + 2.0_ki*r)*log((1.0_ki + b)**2/(4.0_ki*r)))
end function pif

function pis(m) result(p)
   real(ki), intent(in) :: m
   real(ki) :: p, r, b
   real(ki), parameter :: pi = 3.14159265358979323846264338327950288_ki
   r = m**2/s
   b = sqrt(1.0_ki - 4.0_ki*r)
   p = alpha/(12.0_ki*pi)*(-8.0_ki/3.0_ki + 8.0_ki*r &
     & + b**3*log((1.0_ki + b)**2/(4.0_ki*r)))
end function pis

end program check_ref
