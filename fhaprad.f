      SUBROUTINE fhaprad(Ebeam, xc, q2c, zc, ptc, phic, m2th, sib, sig,
     &                   delta, tail)

      IMPLICIT NONE

      DOUBLE PRECISION Ebeam, epsphir, epstau, epsrr, x, y, z, pt
      DOUBLE PRECISION phi_had, xc, q2c, zc, ptc, phic, rcfac
      DOUBLE PRECISION sib, sig, delinf, delta, tai(3), tail
      DOUBLE PRECISION mp, mpi, mn, Sx, nu, Mx2, raddeg, m2th
      INTEGER ilep, iphi_rad

      DATA mp/0.93827d0/, mpi/0.13957d0/, mn/0.93956536d0/
      DATA raddeg/57.2957795131d0/

C Input for HAPRAD
      ilep = 1       ! ilep - registered lepton: 1 - electron,  2 - muon
      iphi_rad = 0   ! iphi_rad - integration over phi_{rad} (1) or approximation (0)
      epsphir = 0.01 ! relative accuracies of integration over phi_r, tau, R
      epstau = 0.001
      epsrr = 0.01
      
      x = xc
      y = -q2c       ! y or -Q2 if negative
      z = zc
      pt = ptc       ! pt or t if negative
      phi_had = phic / raddeg
      
      nu = q2c / (2.d0 * mp * xc)
      Sx = 2.d0 *mp *nu
      Mx2 = mp ** 2 +Sx * (1.d0 - zc) + ptc ! Mass square of unobservable hadrons
      
      IF (Mx2 .GT. m2th) THEN
         CALL ihaprad(Ebeam, ilep, iphi_rad, epsphir, epstau, epsrr, 
     &                x, y, z, pt, phi_had, sib, sig, delinf, delta,
     &                tai)
         rcfac = sig / sib
         sib = sib * 1.d-3
         sig = sig * 1.d-3
         tail = tai(2) * 1.d-3
      ELSE
         sib = 0.d0
         sig = 0.d0
         delta = 0.d0
         tail = 0.d0
      ENDIF
      
      RETURN
      END
