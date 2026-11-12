[%' vim: syntax=golem
'%]module     [% process_name asprefix=\_ %]_ew_ct
!
! On-shell counterterms of scalar QED: photon, charged leptons (e, mu, tau)
! and a stable charged scalar phi+- (mass mW, no width), for the model sqed_ct.
! Installed as common/ew_ct.f90 instead of the Standard-Model ew_ct.f90
! when the model is sqed_ct (template.xml, option internal_ct_model=sqed).
!
! Same interface as the Standard-Model module:
!  - ew_ren() fills the common block /gctvalues/ of dzdecl.h, which the
!    counterterm diagrams (diagramsct.f90) read;
!  - ddrr is the charge-renormalization shift of the G_mu scheme, zero here:
!    the charge is renormalized in the Thomson limit, alpha(0);
!  - each counterterm is an array (1:2): (1) = finite part, (2) = coefficient
!    of 1/eps, both at mu = 1 GeV. matrix.f90 restores the scale with
!    log(scale2) times the pole part. UV and IR poles are not separated.
!
! Conventions: Denner, Fortsch. Phys. 41 (1993) 307, as in ew_ct.f90:
!  - self-energies with alpha/(4 pi) and the loop functions of OneLOop;
!  - lepton self-energies in dimensional reduction (as ew_ct.f90 with
!    dred=.true.): Sigma^{L,R} = -a 2 B1, Sigma^S = -a 4 B0 (a = alpha/4pi);
!    the photon and the phi self-energies are the same in CDR and DRED;
!  - photon self-energy from leptons and phi, derivative at k^2 = 0:
!    dZ_AA = -a ( 4/3 sum_l B0(0,ml,ml) + 1/3 B0(0,M,M) ), dZ_e = -dZ_AA/2;
!  - phi self-energy (photon-phi bubble; the seagull tadpole is scaleless):
!    Sigma_phi(p2) = -a ( 4 p2 (B0 + B1)(p2; 0, M) + A0(M) ),
!    dM^2 = Sigma_phi(M^2), dZ_phi = -Sigma_phi'(M^2) (IR divergent);
!  - counterterm vertices of the model (sqed_ct.hh):
!      AA        C1 = gctA   = dZ_AA,                  C2 = 0
!      phi phi   C1 = gctphi1 = dZ_phi,                C2 = gctphi = M^2 dZ_phi + dM^2
!      l l       CL = CR = dZ_l,                       CP = CM = m dZ_l + dm_l
!      A phi phi gctAPP  = -( dZ_e + dZ_AA/2 + dZ_phi )
!      AA phi phi gctAAPP = 2 ( 2 dZ_e + dZ_AA + dZ_phi )
!      A l l     gctAl = gctAr = -Q ( dZ_e + dZ_AA/2 + dZ_l ),  Q = -1
!    No A-Z mixing, no tadpole, no phi^4 vertex.
!
      implicit none
      private
      complex*16, dimension(1:2) :: ddrr = (0.d0, 0.d0)
      ! phi field renormalization, coefficient C1 of the phi phi counterterm;
      ! not in dzdecl.h because the Standard Model has no such counterterm
      complex*16, dimension(1:2) :: gctphi1 = (0.d0, 0.d0)
      public :: ew_ren, ddrr, gctphi1
contains

      subroutine ew_ren()
      use avh_olo
      use [% process_name asprefix=\_ %]_model, only: me, mmu, mtau, mW, alpha
      implicit none
      include 'dzdecl.h'

      real*8, parameter :: pi = 3.1415926535897932384626433832795029d0
      complex*16, parameter :: cone = (1.d0, 0.d0), zero = (0.d0, 0.d0)
      ! charge of the leptons in units of e (field em = lepton)
      real*8, parameter :: ql = -1.d0

      integer, save :: ini = 0
      complex*16 :: a, mphi2
      complex*16, dimension(0:1) :: dza, dze, dzphi, dmphi2, db0, db1
      complex*16, dimension(0:1,3) :: dzl, dml
      complex*16, dimension(0:2) :: b0, b1, b11, b00, a0
      real*8 :: ml(3), m
      integer :: k, i1

      if (ini .ne. 0) return
      ini = 1

      call olo_onshell(1.d-10)
      call olo_scale(1.d0)

      a = cone * dble(alpha) / 4.d0 / pi
      ml(1) = dble(me)
      ml(2) = dble(mmu)
      ml(3) = dble(mtau)
      mphi2 = cone * dble(mW)**2

      ! photon wave function, charge
      dza(:) = zero
      do k = 1, 3
         if (ml(k) .gt. 0.d0) then
            call olo(b0, zero, cone*ml(k)**2, cone*ml(k)**2)
            dza(:) = dza(:) - a * 4.d0/3.d0 * ql**2 * b0(0:1)
         end if
      end do
      if (abs(mphi2) .gt. 0.d0) then
         call olo(b0, zero, mphi2, mphi2)
         dza(:) = dza(:) - a * 1.d0/3.d0 * b0(0:1)
      end if
      dze(:) = -0.5d0 * dza(:)

      ! charged leptons: mass and wave function (photon loop only)
      do k = 1, 3
         m = ml(k)
         if (m .gt. 0.d0) then
            ! B(m^2; m, 0)
            call olo_b11(b11, b00, b1, b0, cone*m*m, cone*m*m, zero)
            ! IR-singular derivatives at p^2 = m^2 (eqs. 5.51, 5.54 of
            ! Bardin-Passarino, as b0pir and b1pir of ew_ct.f90)
            db0(0) = -(2.d0 - log(m*m)) / (2.d0*m*m) * cone
            db0(1) = -1.d0 / (2.d0*m*m) * cone
            db1(0) =  (3.d0 - log(m*m)) / (2.d0*m*m) * cone
            db1(1) =  1.d0 / (2.d0*m*m) * cone
            dml(:,k) = -a * ql**2 * m * (2.d0*b1(0:1) + 4.d0*b0(0:1))
            dzl(:,k) =  a * ql**2 * (2.d0*b1(0:1) &
                     &   + m*m*(4.d0*db1(:) + 8.d0*db0(:)))
         else
            dml(:,k) = zero
            dzl(:,k) = zero
         end if
      end do

      ! charged scalar: mass and wave function
      dmphi2(:) = zero
      dzphi(:) = zero
      if (abs(mphi2) .gt. 0.d0) then
         ! B(M^2; 0, M): photon first, as in the derivation above
         call olo_b11(b11, b00, b1, b0, mphi2, zero, mphi2)
         call olo(a0, mphi2)
         dmphi2(:) = -a * (4.d0*mphi2*(b0(0:1) + b1(0:1)) + a0(0:1))
         ! IR-singular B0'(M^2; 0, M), as b0pir of ew_ct.f90
         db0(0) = -(2.d0*cone - log(mphi2)) / (2.d0*mphi2)
         db0(1) = -1.d0 / (2.d0*mphi2)
         ! Sigma_phi(p2) = -a ( 2 (p2 + M^2) B0 - A0 ), derivative at M^2
         dzphi(:) = a * (2.d0*b0(0:1) + 4.d0*mphi2*db0(:))
      end if

      ! counterterm couplings: (1) finite, (2) pole
      do i1 = 1, 2
         gctA(i1)    = dza(i1-1)
         gctphi1(i1) = dzphi(i1-1)
         gctphi(i1)  = mphi2*dzphi(i1-1) + dmphi2(i1-1)
         gctAPP(i1)  = -(dze(i1-1) + 0.5d0*dza(i1-1) + dzphi(i1-1))
         gctAAPP(i1) = 2.d0*(2.d0*dze(i1-1) + dza(i1-1) + dzphi(i1-1))

         gctCLe(i1) = dzl(i1-1,1)
         gctCRe(i1) = dzl(i1-1,1)
         gctCPe(i1) = ml(1)*dzl(i1-1,1) + dml(i1-1,1)
         gctCMe(i1) = ml(1)*dzl(i1-1,1) + dml(i1-1,1)
         gctAle(i1) = -ql*(dze(i1-1) + 0.5d0*dza(i1-1) + dzl(i1-1,1))
         gctAre(i1) = gctAle(i1)

         gctCLmu(i1) = dzl(i1-1,2)
         gctCRmu(i1) = dzl(i1-1,2)
         gctCPmu(i1) = ml(2)*dzl(i1-1,2) + dml(i1-1,2)
         gctCMmu(i1) = ml(2)*dzl(i1-1,2) + dml(i1-1,2)
         gctAlmu(i1) = -ql*(dze(i1-1) + 0.5d0*dza(i1-1) + dzl(i1-1,2))
         gctArmu(i1) = gctAlmu(i1)

         gctCLtau(i1) = dzl(i1-1,3)
         gctCRtau(i1) = dzl(i1-1,3)
         gctCPtau(i1) = ml(3)*dzl(i1-1,3) + dml(i1-1,3)
         gctCMtau(i1) = ml(3)*dzl(i1-1,3) + dml(i1-1,3)
         gctAltau(i1) = -ql*(dze(i1-1) + 0.5d0*dza(i1-1) + dzl(i1-1,3))
         gctArtau(i1) = gctAltau(i1)
      end do

      end subroutine ew_ren

end module [% process_name asprefix=\_ %]_ew_ct
