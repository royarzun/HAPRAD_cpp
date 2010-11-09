      PROGRAM rcdat
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'

      DOUBLE PRECISION Ebeam, x, q2, z, t, phi, nu
      DOUBLE PRECISION sib, sig, delta, tail
      DOUBLE PRECISION mn, DSIampS, rc, s
      DOUBLE PRECISION phi_min, phi_max, t_min, t_max, z_min, z_max
      DOUBLE PRECISION m2th, x_min, x_max

      DATA mn/0.93956536d0/

      amhh = 0.1395675d0      ! mass of detected hadron
      amhu = amp!0.93956536d0 ! mass of undetected hadron for exclusive radiative tail 
      m2th = (mn + ampi) ** 2 ! minimal invariant mass square of undetected hadron
C     m2th = 1.4d0 ** 2       ! another limit on  minimal invariant mass square of undetected hadron 
      Ebeam = 6d0             ! energy of lepton beam
      s = 2d0 * amp * Ebeam
      q2 = 2.5d0
      x_min = q2 / (2d0 * amp * Ebeam)
      x_max = q2 / (q2 + m2th - amp ** 2)
      phi_min = 0d0
      phi_max = 360d0
      x = 0.32d0
      nu = q2 / (2d0 * amp * x)
      z_min = amhh / nu
      z_max = 1d0
      phi = 140d0
      z = 0.1d0

      t_min = max(m2th - amp ** 2 - q2 / x* (1d0 - z) + 1d-4,
     &            -q2 + amhh ** 2 -2d0 * (nu ** 2 * z +
     &            sqrt((nu ** 2 + q2) * ((z * nu) ** 2 - amhh ** 2))))

      t_max = -q2 + amhh ** 2 - 2d0 * (nu ** 2 * z -
     &        sqrt((nu ** 2 + q2) * ((z*nu)**2-amhh**2)))

      t = 0.5d0 * (t_max + t_min)

      OPEN (8,FILE='res.dat') 
      WRITE (8,'(9a12)') 'ebeam     ', 'q2      ', 'x      ',
     &                   'phi      ', 'z     ', 't     ', 'sib    ',
     &                   'sig  ', 'tail'

      CALL fhaprad(Ebeam, x, q2, z, t, phi, m2th, sib, sig, delta, tail)

      WRITE (8,'(9g12.4)') ebeam, q2, x, phi, z, t, sib, sig, tail
      WRITE (*,'(9g12.4)') ebeam, q2, x, phi, z, t, sib, sig, tail
      CLOSE(8)
      END
