      SUBROUTINE ihaprad(bmom, ilepm, iphi_radm, epsphirm, epstaum,
     &                   epsrrm, xmas, ymas, zmas, tmas, phimas, sib,
     &                   sig, delinf, delta, tai)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'epsil.inc'
      INCLUDE 'epsmarch.inc'

      DOUBLE PRECISION bmom, epsphirm, epstaum, epsrrm, xmas, ymas
      DOUBLE PRECISION zmas, tmas, phimas, sib, sig, delinf, delta
      DOUBLE PRECISION tai(3), tmom, snuc, yma, ymi, sqnuq, p22max
      DOUBLE PRECISION tdmax, epspl

      INTEGER ilepm, iphi_radm
      EXTERNAL sphih
      
      ipol = 0       ! ipol - 1 -- long; 2 -- tran; 0 -- unpol <-- type of polarization
      iphi_had = 0   ! iphi_had - integration over phi_{had} (1) or not (0)
      iphi_rad = iphi_radm
      ilep = ilepm
      epsphir = epsphirm
      epstau = epstaum
      epsrr = epsrrm
      
      
      isf1 = 1
      isf2 = isf20
      isf3 = 1
      un = 1.d0
      pl = 0.d0
      pn = 0.d0
      qn = 0.d0

      tmom = 0.d0     ! tmom - momentum per nucleon
C     snuc = 2.d0*(sqrt(tmom**2+amh**2)*sqrt(bmom**2+aml2)+bmom*tmom)
C     snuc = 2.d0*(sqrt(tmom**2+amh**2)*sqrt(bmom**2)+bmom*tmom)
      snuc = 2.d0*amp*bmom
C     PRINT *, 'tmas', tmas
      
      xs = xmas
      zdif = zmas
      tdif = tmas
      IF (ymas .GE. 0.0) THEN
         ys = ymas
         y = snuc * xs * ys      !  q2
      ELSE
         y = -ymas      !  q2
         ys = y / (snuc * xs)
      ENDIF
      yma = 1.d0 / (1.d0 + amp * *2 * xs / snuc)
      amc2 = (amp + amhh) ** 2
      ymi = (amc2 - amp ** 2) / (snuc * (1.d0 - xs))
      IF (ys .GT. yma .OR.
     &         ys .LT. ymi .OR.
     &         xs .GT. 1.d0 .OR.
     &         xs .LT. 0.d0) THEN
         print *, ' Warning! Wrong kinematics!!!! skeep the point!'
         print *, ' ys =  ', ys
         print *, ' xs =  ', xs
         RETURN
      ENDIF
      CALL conkin(snuc)
      ehad = anu * zdif
      sqnuq = SQRT(anu ** 2 + y)
      
      IF (ehad .LT. amhh) THEN
         PRINT *, ' Warning! Wrong kinematics!!!! skeep the point!'
         PRINT *, ' ehad  = ', ehad
         RETURN
      ENDIF
      
      pph = SQRT(ehad ** 2 - amhh ** 2)
      
      IF (tdif .GE. 0.d0) THEN
         pth = tdif
         IF (pph .LT. pth) THEN
            PRINT *, ' Warning! Wrong kinematics!!!! skeep the point!'
            PRINT *, ' pph  = ', pph
            PRINT *, ' pth  = ', pth
            RETURN
         ENDIF

         plh = SQRT(pph ** 2 - pth ** 2)
         IF (pph .GT. pth) THEN
            an = an * sqly / 2.d0 / amp / plh
         ELSE
            an = 0.d0
         ENDIF

         tdif = amhh ** 2 - y + 2.d0 * (sqnuq * plh - anu * ehad)
         PRINT *, 'pl: ', plh, tdif, (plh - (tdif + y - amhh ** 2 +
     &                                2.d0 * anu * ehad) / 2.d0 / sqnuq)
      ELSE
         plh = (tdif + y - amhh ** 2 + 2.d0 * anu * ehad) / 2.d0 / sqnuq
         IF (pph .LT. ABS(plh)) THEN
            epspl = SQRT((tdif * epsmarch / sqnuq) ** 2 +
     &              (2.d0 * amhh ** 2 * epsmarch / sqnuq) ** 2 +
     &              2.d0 * (2.d0 * anu * ehad * epsmarch / sqnuq) ** 2 +
     &              ((tdif + y - amhh ** 2 + 2.d0 * anu * ehad) /
     &              sqnuq * epsmarch) ** 2) / 2.d0
            IF (ABS(pph - ABS(plh)) .GT. epspl) THEN
               PRINT *, ' Warning! Wrong kinematics! Skeep the point!'
               PRINT *, ' pph  = ', pph
               PRINT *, ' plh  = ', plh, (pph - ABS(plh))
C              IF (ABS(pph - plh) .LT. 1.d-9) STOP
               RETURN
            ELSE
               PRINT *, 'Zero pt!', plh, (pph - ABS(plh)), epspl
               plh = SIGN(1.d0, plh) * pph
            ENDIF
         ENDIF
         pth = SQRT(pph ** 2 - plh ** 2)
C        WRITE (*, '(4f22.19, f7.4)') ehad, pph, pth, plh, anu
      ENDIF

C     p22max = (sqrt(w2) - amp) ** 2
      p22max = w2 - amhh ** 2
      p22 = amp2 + sx * (1.d0 - zdif) + tdif
      IF (p22 .LT. amc2 .OR. p22 .GT. p22max) THEN
         PRINT *, ' Warning! Wrong kinematics! Skeep the point!'
         PRINT *, ' p22  = ', p22
         PRINT *, ' amc2 = ', amc2
         PRINT *, ' p22m = ', p22max
         RETURN
      ENDIF

      tdmin = amhh ** 2 - y + 2.d0 * ( sqnuq * pph - anu * ehad)
      tdmax = amhh ** 2 - y + 2.d0 * (-sqnuq * pph - anu * ehad)
      
C     PRINT *, tdif, tdmin, tdmax, (tdif - tdmin)
      
C     STOP
C     tdif = tdif + tdmin

      IF ((tdif - tdmin) .GT. epsmarch .OR. tdif .LT. tdmax) THEN
         PRINT *, ' Warning! Wrong kinematics! Skeep the point!'
         PRINT *, ' tdif  = ', tdif
         PRINT *, ' tdmax = ', tdmax
         PRINT *, ' tdmin = ', tdmin, (tdif - tdmin)
         RETURN
      ENDIF
      
      
C     plh = (tdif + y - amhh ** 2 + 2.d0 * anu * ehad) / 2.d0 / sqnuq
C     plh = (amhh ** 2 - y + 2.d0 * (sqnuq * SQRT(pph ** 2 -
C    &         pth ** 2) - anu * ehad) + y - amhh ** 2 + 2.d0 * anu *
C    &         ehad) / (2.d0 * sqnuq)
C     dum = SQRT(pph ** 2 - pth ** 2)
C     plh = dum - (amp * dum + ehad - ehad) / amp
C     PRINT *, 'diff: ', plh, (plh - dum), sqnuq, tdif, anu, ehad
      
      
     
C     PRINT *, y, xs, zdif, pth, pph, plh, tdif, anu, ehad, sqnuq
      
      CALL sphih(phimas, sib, sig, delinf, delta, tai)
      RETURN
      END



      SUBROUTINE sphih(phih, sib, sig, delinf, delta, tai)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      DOUBLE PRECISION phih, sib, sibt, sig, delinf, delta, tai(3)
      DOUBLE PRECISION costs, sints, phk1, costx, sintx, phk2, tr
      DOUBLE PRECISION extai1
      INTEGER il, infin, in

      phidif = phih
      costs = (s * (s - x) + 2.d0 * amp2 * y) / sqls / sqly

      IF ((s * x * y - amp2 * y ** 2 - aml2 * aly) .GT. 0.d0) THEN
         sints = 2.d0 * amp * SQRT(s * x * y - amp2 * y ** 2 -
     &                        aml2 * aly) / sqls / sqly
      ELSE
         PRINT *, 'sphih: sints = NaN ', (s * x * y - amp2 * y ** 2 -
     &                                    aml2 * aly)
         sints = 0.d0
      ENDIF
      vv10 = (s * ehad - sqls *
     &        (costs * plh + sints * pth * COS(phidif))) / amp
      phk1 = 0.5d0 * (s * ehad - sqls * costs * plh) / amp
      phk12 = -0.5d0 * sqls * sints * pth / amp
      costx = (x * (s - x) - 2.d0 * amp2 * y) / sqlx / sqly
      IF ((s * x * y - amp2 * y ** 2 - aml2 * aly) .GT. 0.d0) THEN
         sintx = 2.d0 * amp * SQRT(s * x * y - amp2 * y ** 2 -
     &           aml2 * aly) / sqlx / sqly
      ELSE
         PRINT *, 'sphih: sintx = NaN ', (s * x * y - amp2 * y ** 2 -
     &                                    aml2 * aly)
         sintx = 0.d0
      ENDIF
      vv20 = (x * ehad - sqlx * (costx * plh + sintx * pth *
     &        COS(phidif))) / amp
      phk2 = 0.5d0 * (x * ehad - sqlx * costx * plh) / amp
      phkp = phk1 + phk2
      phkm = phk2 - phk1
C     PRINT  *, sints ** 2 + costs ** 2
C     PRINT  *, sintx ** 2 + costx ** 2

C Delta is factorizing part of virtual AND real leptonic bremsstrahlung
      CALL deltas(delta, delinf, tr)
      DO 10 il = 1, 1
      infin = 3
      IF (ipol .EQ. 0) infin = 1
      DO 10 in = 1, infin, 2
         DO 30 ita = 1, 2
C           WRITE (9, '(10(1h*), '' ita  =  '', i2, 10(1h*))') ita
            WRITE (*, '(10(1h*), '' ita  =  '', i2, 10(1h*))') ita
            IF (ita .EQ. 1) THEN
               CALL bornin(sib)
               CALL bornintest(sibt)
               PRINT *, 'sib1', sib
               PRINT *, 'sibt', sibt
C              STOP
               IF (sib .EQ. 0.0)  THEN
                  tai(1) = 0.d0
                  GOTO 30
               ENDIF
C              STOP
            ENDIF
            CALL qqt(tai(ita))
            WRITE (*, '(a5, i1, a2, g12.4)') 'tai(', ita, ') = ',
     &                                       tai(ita)
   30    CONTINUE
         delinf = 0d0
         extai1 = EXP(alfa / pi * delinf)
         sig = sib * extai1 * (1.d0 + alfa / pi * (delta - delinf)) +
     &         tai(1) + tai(2)
C        PRINT  *, 'dd ', ((1.d0 - zdif) * sx + tdif), sib, tai(2),
C    &          (tai(2) / sib)
   10 CONTINUE
      RETURN
      END



CDECK  ID>,  CONKIN.
      SUBROUTINE conkin(snuc)
C Set of kinematical constants
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'pol.inc'
      DOUBLE PRECISION snuc, axy, sqn
      
      amt = amp
      aml2 = amlep(ilep)
      aml2 = amlep(ilep) ** 2
      al2 = 2.d0 * aml2
      pi2 = pi ** 2
      amh = amp
      ap = 2.d0 * amp
      amp2 = amp ** 2
      ap2 = 2.d0 * amp2
      
      s = snuc * amp / amh
      x = s * (1.d0 - ys)
      sx = s - x
      sxp = s + x
      ym = y + al2
      tpl = s ** 2 + x ** 2
      tmi = s ** 2 - x ** 2
      w2 = amp2 + s - y - x
      als = s * s - al2 * ap2
      alx = x * x - al2 * ap2
      alm = y * y + 4.d0 * aml2 * y
      aly = sx ** 2 + 4.d0 * amp2 * y
      IF (als .LT. 0.d0)  PRINT *, 'conkin: als<0 ', als
      sqls = SQRT(MAX(0.d0, als))
      IF (alx .LT. 0.d0)  PRINT *, 'conkin: alx<0 ', alx
      sqlx = SQRT(MAX(0.d0, alx))
      IF (aly .LT. 0.d0)  PRINT *, 'conkin: aly<0 ', aly
      sqly = SQRT(MAX(0.d0, aly))
      IF (alm .LT. 0.d0)  PRINT *, 'conkin: alm<0 ', alm
      sqlm = SQRT(MAX(0.d0, alm))
      allm = LOG((sqlm + y) / (sqlm - y)) / sqlm
      anu = sx / ap
      axy = pi * (s - x)
C     an = 2. * alfa ** 2 / sqls * axy * barn * amh / amp
      an = pi * alfa ** 2 * ys * sx * amp / 2.d0 / sqly * barn
C     print *, 'an1', an
C     tamin = (sx - sqly) / ap2
      tamax = (sx + sqly) / ap2
      tamin =  -y / amp2 / tamax
      as = s / 2.d0 / aml / sqls
      bs = 0.d0
      cs =  -aml / sqls
      IF (ipol .NE. 2) THEN
         ae = amp / sqls
         be = 0.d0
         ce =  - s / ap / sqls
      ELSE
         IF ((s * x * y - aly * aml2 - amp2 * y * y)  .GT. 0.d0) THEN
            sqn = SQRT(s * x * y - aly * aml2 - amp2 * y * y)
         ELSE
            PRINT *, 'conkin: sqn = NaN ', (s *x * y - aly * aml2 -
     &      amp2 * y * y)
            sqn = 0.d0
         ENDIF
         ae = (-s * x + ap2 * ym) / sqls / sqn / 2.d0
         be = sqls / sqn / 2.d0
         ce =  -(s * y + al2 * sx) / sqls / sqn / 2.d0
      ENDIF
      apq =  - y * (ae - be) + ce * sx
      apn = (y + 4.d0 * aml2) * (ae + be) + ce * sxp
      dk2ks = as * ym + al2 * bs + cs * x
      dksp1 = as * s + bs * x + cs * ap2
      dapks = 2.d0 * (al2 * (as * ae + bs * be) + ap2 * cs * ce + ym *
     &        (as * be + bs * ae)+ s * (as * ce + cs * ae) + x *
     &        (bs * ce + cs * be))
      RETURN
      END



CDECK  ID>,  BORNIN.
      SUBROUTINE bornintest(sibort)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'phi.inc'
      DOUBLE PRECISION sibort, sfm0(8), q2, Eb, H1, H2, H3, H4, Eh, rt
      DOUBLE PRECISION rtz, gev2mb, semi
C     CALL strf(0.d0, 0.d0, 0.d0, sfm0)
      q2 = y
      Eb = s / 2.d0 / amp
      CALL semi_inclusive_model(q2, xs, ys, zdif, pth**2, p22, plh, H1,
     &                          H2, H3, H4)
      Eh = zdif * sx / 2.d0 / amp
      rt = 1.d0 - ys - amp * xs * ys / (2.d0 * Eb)
      rtz = SQRT(rt / (1.d0 + 2.d0 * amp * xs / (ys * Eb)))
      gev2mb = 389.37929d0
C     semi = (1.0d + 3) * gev2mb * (4.d0 * pi * alfa ** 2 * amp * Eb) /
C    &       (q2 ** 2) * Eh / plh
      semi = (16.d0 * an * amp * Eb) / (q2 *  * 2) * Eh / ys / sx * 
     &       (xs * ys ** 2 * H1 + rt * H2 + pth / SQRT(q2) *
     &       (2.d0 - ys) * rtz * COS(phidif) * H3 + pth ** 2 / q2 *
     &       rtz ** 2 * COS(2.d0 * phidif) * H4)
C     PRINT *, 'born ', semi * 1.d-3, (xs * ys ** 2 * H1 + rt * H2),
C    &      (pth / sqrt(q2) * (2.d0 - ys) * rtz * cos(phidif) * H3),
C    &      (pth ** 2 / q2 * rtz ** 2 * COS(2.d0 * phidif) * H4), phidif
C     PRINT *, 'btest ', Eb, Eh, plh, q2, xs, zdif, pth, phidif,
C    &      H1, H2, H3, H4
C     STOP
      sibort = semi
      RETURN
      END



      SUBROUTINE bornin(sibor)
C     sibor is born cross section with polarized initial
C     lepton AND polarized target
C     siamm is contribution of anomalous magnetic moment.
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'pol.inc'
      INCLUDE 'print.inc'
      DOUBLE PRECISION sibor, sfm0(8), tm(8), ssum, aa, hi2
      INTEGER isf
      ipri1 = 1
      CALL strf(0.d0, 0.d0, 0.d0, sfm0)
      ipri1 = 0

      hi2 = 0.25d0

      tm(1) = y
      tm(2) = (s * x - amp2 * y) / 2.d0
      tm(3) = (vv10 * vv20 - amhh ** 2 * y) / 2.d0
      tm(4) = (s * vv20 + x * vv10 - zdif * sx * y) / 2.d0
      aa = sx * (zdif - 2.d0 * amp * plh / sqly) / 2.d0 / amp2
      ssum = 0.d0
      DO 1 isf = isf1, isf2, isf3
         ssum = ssum + tm(isf) * sfm0(isf)
C        PRINT *, isf
C        PRINT *, ssum, tm(isf), sfm0(isf)
C        PRINT *, ssum
    1 CONTINUE
C     PRINT *, 'an2', an, ssum, y, tm(2), sfm0(2)
      sibor = ssum * an / y / y * 2.d0
      RETURN
      END



CDECK  ID>,  DELTAS.
      SUBROUTINE deltas(delta, delinf, tr)
C Delta is factorizing part of virtual AND real leptonic bremsstrahlung
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'phi.inc'
      DOUBLE PRECISION delta, delinf, tr, del1, del2, sum, aj0, deltai
      DOUBLE PRECISION ssh, xxh, alss, alxx, sqlss, sqlxx, allss, allxx
      DOUBLE PRECISION dlm, delta0, delta_old, fspen, vacpol, sfpr

      del1 = -ym * (alm * allm *  * 2 / 2.d0 + 2.d0 *
     &       fspen(2.d0 * sqlm / (y + sqlm)) - pi2 / 2.d0) / sqlm
      del2 = (3.d0 * y / 2.d0 + 4.d0 * aml2) * allm - 2.d0

      sum = vacpol(y)

      aj0 = 2.d0 * (ym * allm - 1.d0)
      deltai = aj0 * LOG((p22 - amc2) / aml / SQRT(p22))

      ssh = x + y - vv20
      xxh = s - y - vv10
      alss = ssh ** 2 - 2.d0 * p22 * al2
      alxx = xxh ** 2 - 2.d0 * p22 * al2
      IF (alss .LT. 0.d0) PRINT *, 'deltas: alss<0 ', alss
      sqlss = SQRT(MAX(0.d0, alss))
      IF (alxx .LT. 0.d0) PRINT *, 'deltas: alxx<0 ', alxx
      sqlxx = SQRT(MAX(0.d0, alxx))
      allss = LOG((sqlss + ssh) / ( - sqlss + ssh)) / sqlss
      allxx = LOG((sqlxx + xxh) / ( - sqlxx + xxh)) / sqlxx
      dlm = LOG(y / aml2)
      sfpr = dlm ** 2 / 2.d0 - dlm * LOG(ssh * xxh / aml2 / p22) -
     &       (LOG(ssh / xxh)) ** 2 / 2.d0 +
     &       fspen(1.d0 - p22 * y / ssh / xxh) - pi2 / 3.d0
      delta0 = (ssh * allss + xxh * allxx) / 2.d0 + sfpr
      delta_old = deltai + delta0 + del1 + del2 + sum
      delinf = (dlm - 1.d0) * LOG((p22 - amc2) ** 2 / ssh / xxh)
      tr = alfa / pi * (dlm - 1.d0)
      delta = delinf + sum + (1.5d0 * dlm - 2.d0 - 0.5d0 *
     &        LOG(xxh / ssh) ** 2 +
     &        fspen(1.d0 - p22 * y / ssh / xxh) - pi2 / 6.d0)
C     WRITE (*, '(6hdelta:, 2g12.4)') delta_old, delta
      RETURN
      END



CDECK  ID>,  VACPOL.
      DOUBLE PRECISION FUNCTION vacpol(t)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      DOUBLE PRECISION t, am2(3), suml, a2, sqlmi
      DOUBLE PRECISION allmi, aaa, bbb, ccc, sumh
      INTEGER i
C am2 : squared masses of charge leptons
      DATA am2/0.26110d-6, 0.111637d-1, 3.18301d0/

      suml = 0.d0
      DO 10 i = 1, 3
         a2 = 2.d0 * am2(i)
         sqlmi = SQRT(t * t + 2.d0 * a2 * t)
         allmi = LOG((sqlmi + t) / (sqlmi - t)) / sqlmi
   10 suml = suml + 2.d0 * (t + a2) * allmi / 3.d0 - 10.d0 / 9.d0 +
     &       4.d0 * a2 * (1.d0 - a2 * allmi) / 3.d0 / t
      IF (t .LT. 1.d0) THEN
         aaa  =  -1.345d-9
         bbb  =  -2.302d-3
         ccc  =   4.091d0
      ELSEIF (t .LT. 64.d0) THEN
         aaa  =  -1.512d-3
         bbb  =  -2.822d-3
         ccc  =   1.218d0
      ELSE
         aaa  =  -1.1344d-3
         bbb  =  -3.0680d-3
         ccc  =   9.9992d-1
      ENDIF
      sumh  =  -(aaa + bbb * LOG(1.d0 + ccc * t)) * 2.d0 *pi / alfa
      vacpol = suml + sumh
      RETURN
      END



CDECK  ID>,  QQT.
      SUBROUTINE qqt(tai)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'epsil.inc'
      DOUBLE PRECISION tai, qqtphi, rv2tr, phiar(4), tar(6), ta1, ta2
      DOUBLE PRECISION ot, otr, am(2), bm(2), wrk(500), rere, re, re2
      INTEGER iph, ico, id, mir, ma, i
      EXTERNAL qqtphi
      EXTERNAL rv2tr
     
      IF (ita .EQ. 1) THEN
         IF (iphi_rad .EQ. 1) THEN
            CALL simpsx(0.d0, (2.d0*pi), 150, epsphir, qqtphi, tai)
            tai = an * alfa * tai / pi ** 2 / 4.d0 / sqly
         ELSEIF (iphi_rad .EQ. 0) THEN
            tai = an * alfa / pi * qqtphi(0.d0) / 2.d0 / sqly
         ENDIF
C        STOP
      ELSE
         ta1 =  - y / s     
         ta2 = y / x     
         phiar(1) = 0.d0
         phiar(2) = 0.01d0 * pi
         phiar(3) = 2.d0 * pi - 0.01d0 * pi
         phiar(4) = 2.d0 * pi
         tar(1) = tamin
         tar(2) = ta1 - 0.15d0 * (ta1 - tamin)
         tar(3) = ta1 + 0.15d0 * (ta2 - ta1)
         tar(4) = ta2 - 0.15d0 * (ta2 - ta1)
         tar(5) = ta2 + 0.15d0 * (tamax - ta2)
         tar(6) = tamax
         ot = 1.d-3

         rere = 0.d0
         DO iph = 1, 3
            DO ico = 1, 5
               am(1) = tar(ico)
               bm(1) = tar(ico+1)
               am(2) = phiar(iph)
               bm(2) = phiar(iph+1)
               IF (am(2) .GT. bm(2)) WRITE (*, *)' am(2)<bm(2)'
               IF (am(1) .GT. bm(1)) WRITE (*, *)' am(1)<bm(2)'
               id = 1
               mir = 10000
               ma = 10 * mir
C              PRINT *, 'start d01fce integration!'
               CALL d01fce(2, am, bm, mir, ma, rv2tr, ot, otr, 500,
     &                     wrk, re, id)
C              WRITE (9, '(1x, ''tai:'', 2i3, g13.4, f8.4, i9, i3)')
C    &               ico, iph, re, otr, mir, id
               WRITE (*, '(1x, ''tai:'', 2i3, g15.6, f8.4, i9, i3)')
     &               ico, iph, re, otr, mir, id
               rere = rere + re
            ENDDO
         ENDDO
         tai =  - alfa / (64.d0 * pi ** 5 * sqly * amp) * an * rere
      ENDIF
      RETURN
      END



****************** rv2tr **************************************
      DOUBLE PRECISION FUNCTION rv2tr(ndim, arg)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'ex.inc'
      DOUBLE PRECISION arg(15), tm(8, 6), sfm(4), vv, fwiw, ta, dmu, r
      DOUBLE PRECISION tldq2, tldw2, tldtm, podinlz, pres
      INTEGER ndim, isf, irr
      IF (ndim .LT. 2 .OR. ndim .GT. 15) STOP
      phirad = arg(2)
      ta = arg(1)
      amh2 = amhh ** 2
      amu2 = amhu ** 2
      vv = (1.d0 - zdif) * sx + tdif + amp2 - amu2
      dmu = (ehad - plh * (sx - ta * ap2) / sqly - ap2 * pth *
     &      SQRT((ta - tamin) * (tamax - ta)) *
     &      COS(phidif - phirad) / sqly) / amp     
      fwiw = 1.d0 + ta - dmu
      r = vv / fwiw
C     tt =  - y + amhad(ivec) ** 2 - vv1 + vv2
C     print *, 't', tt
C     print *, 't', tdif
      tldq2 = y + r * ta
      tldw2 = w2 - r * (1.d0 + ta)
      tldtm = tdif - r * (ta - dmu)
      CALL sffun(sfm, tldq2, tldw2, tldtm)
      CALL tails(ta, tm, dmu)
      podinlz = 0.d0
      DO  isf = 1, 4
         DO  irr = 1, 3
            pres = sfm(isf) * tm(isf, irr) / tldq2 ** 2 * r ** (irr - 2)
            podinlz = podinlz + pres
         ENDDO
      ENDDO
      rv2tr = podinlz / fwiw
      RETURN
      END



      SUBROUTINE sffun(sfm, q2, w2, t)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'ex.inc'
      DOUBLE PRECISION sfm(4), q2, w2, t, sqw2, sx, aly, sxt, tq, sqlw
      DOUBLE PRECISION sqll, cspion, st, sl, stt, slt, sltp, sfm10
      DOUBLE PRECISION sfm20, sfm2tl, sfm4tl, sfm4tt, sfm3tt, sfm2tt
      DOUBLE PRECISION sfm5tl, coetr
      sqw2 = SQRT(w2)
      sx = w2 + q2 - amp2
      aly = sx ** 2 + 4.d0 * amp2 * q2
      sxt = sx + t + amp2 - amu2
      tq = t + q2 - amh2

      IF (((W2 - amu2 - amh2) ** 2 - 4.d0 * amu2 * amh2) .LT. 0.d0)
     &   PRINT *, 'sffun: sqlw = NaN ', ((W2 - amu2 - amh2) ** 2 -
     &         4.d0 * amu2 * amh2)
      sqlw = SQRT(MAX(0.d0,
     &                ((W2 - amu2 - amh2) ** 2 - 4.d0 * amu2 * amh2)))
      IF ((q2*sxt**2-sxt*sx*tq-amp2*tq**2-amh2*aly) .LT. 0.d0)
     &   PRINT *, 'ssffun: qll = NaN ', (q2 * sxt ** 2 -
     &         sxt * sx * tq - amp2 * tq ** 2 - amh2 * aly)
      sqll = SQRT(MAX(0.d0, (q2*sxt**2-sxt*sx*tq-amp2*tq**2-amh2*aly)))
      cspion = (2.d0*tq*w2+(sx-2.d0*q2)*(w2+amh2-amu2))/sqlw/SQRT(aly)

C     Eh = sxt / 2. / amp
C     anu = sx / 2. / amp
C     qm = sqrt(anu ** 2 + q2)
C     Eh_cm = (w2 + amh2 - amp2) / (2.d0 * sqw2)
C     ph_cm = sqrt(Eh_cm ** 2 - amh2)
C     bet = qm / (anu + amp)
C     gam = (anu + amp) / sqw2
C     c_th_cm = (Eh - gam * Eh_cm) / (bet * gam * ph_cm)
C     PRINT *, 'costhcm =  ', cspion, c_th_cm
C     STOP

C Exclusive peak model (cross sections sigma_L, T, LT... from MAID2003)
      CALL exclusive_model(q2, sqw2, cspion, st, sl, stt, slt, sltp)
      
C Structure functions
      IF (aly .GT. 0.0 .AND. sqll .GT. 0.0 .AND. sqlw .GT. 0.0) THEN
         sfm10 = (st - stt)
         sfm20 = 4.d0 * (st + sl) * q2 / aly
         sfm2tl = 2.d0 * slt * sqrt(q2) * (-sx * tq + 2 * q2 * sxt) /
     &            (aly * sqll)
         sfm4tl = -slt * sqrt(q2) / sqll
         sfm4tt = -2.d0 * stt *(-sx * tq + 2 * q2 * sxt) / sqll ** 2
         sfm3tt = 2.d0 * stt * aly / sqll ** 2
         sfm2tt = 2.d0 *stt *((-sx * tq + 2 * q2 * sxt) ** 2 -
     &            2 * q2 * sqll **2 ) / (aly * sqll ** 2)
         sfm5tl = -sltp * sqrt(q2) / sqll
         coetr = 16.d0 *pi * (w2 - amp2) * w2 / (alfa * sqlw) / barn *
     &           1.d3
         sfm(1) = coetr * sfm10
         sfm(2) = coetr * (sfm20 + sfm2tl + sfm2tt)
         sfm(3) = coetr * sfm3tt
         sfm(4) = coetr * (sfm4tl + sfm4tt)
      ELSE
         sfm(1) = 0.d0
         sfm(2) = 0.d0
         sfm(3) = 0.d0
         sfm(4) = 0.d0
      ENDIF
       
      IF (sfm(2) .NE. sfm(2)) PRINT *, 'sffun: ', coetr, st, sl, stt,
     &         slt, sltp, q2, aly, sx, sqw2, cspion, sxt, sqll
      RETURN
      END



      DOUBLE PRECISION FUNCTION qqtphi(phi)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'epsil.inc'
      DOUBLE PRECISION phi, tlm(4), ta1, ta2, tar(6), res, re, ep
      INTEGER i
      external rv2ln
      phirad = phi
      ep = 1.d-12
      ta1 =  -y / s  
      ta2 = y / x  
      tar(1) = tamin
      tar(2) = ta1 - 0.15d0 * (ta1 - tamin)
      tar(3) = ta1 + 0.15d0 * (ta2 - ta1)
      tar(4) = ta2 - 0.15d0 * (ta2 - ta1)
      tar(5) = ta2 + 0.15d0 * (tamax - ta2)
      tar(6) = tamax
      res = 0.d0
      DO i = 1, 5
         CALL simptx(LOG(xs + tar(i)) + ep, LOG(xs + tar(i+1)) - ep,
     &               100, epstau, rv2ln, re)
         res = res + re
      ENDDO 
      qqtphi = res
      RETURN
      END



CDECK  ID>,  TAILS.
      SUBROUTINE tails(ta, tm, amu)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'pol.inc'
      DOUBLE PRECISION ta, tm(8, 6), phka, phkb, b2, b1, c1, c2, bb
      DOUBLE PRECISION sc1, sc2, bi12, bi1pi2, bis, bir, b1i, b11i
      DOUBLE PRECISION sqrtmb, z1, z2, hi2, amh2, zh, vvp, vvm, amu

      IF (iphi_rad .EQ. 0 .AND. ita .EQ. 1) THEN
         b2 = ( - aly * ta + sxp * sx * ta + 2.d0 * sxp * y) / 2.d0
         b1 = ( - aly * ta - sxp * sx * ta - 2.d0 * sxp * y) / 2.d0
         c1 = -(4.d0 * (amp2 * ta ** 2-sx * ta-y) * aml2 -
     &        (s * ta + y) ** 2)
         c2 = -(4.d0 * (amp2 * ta ** 2 - sx * ta - y) * aml2 -
     &        (ta * x - y) ** 2)
         bb = 1.d0!/sqly
         IF (c1 .LT. 0.d0) PRINT *, 'tails: sc1 = NaN ', c1
         sc1 = SQRT(MAX(0.d0, c1))
         IF (c2 .LT. 0.d0) PRINT *, 'tails: sc2 = NaN ', c2
         sc2 = SQRT(MAX(0.d0, c2))
         bi12 = sqly * (sxp * (sx * ta + 2.d0 * y)) /
     &          (sc1 * sc2 * (sc1 + sc2))
         bi1pi2 = sqly / sc2 + sqly / sc1
         bis = sqly * ( - b1 / sc1 / c1 + b2 / sc2 / c2) * aml2
         bir = sqly * (b2 / sc2 / c2 + b1 / sc1 / c1) * aml2
         b1i =  - sqly * b1 / aly / sqly
         b11i = sqly * (3.d0 * b1 ** 2 - aly * c1) / 2.d0 /
     &          aly ** 2 / sqly
      ELSE
         IF ((ta - tamin) * (tamax - ta) *
     &         (s * x * y - y ** 2 * amp2 - aml2 * aly) .GT. 0.d0) THEN
            sqrtmb = SQRT((ta - tamin) * (tamax - ta) *
     &               (s * x * y - y ** 2 * amp2 - aml2 * aly))
         ELSE
            PRINT *, 'tails: sqrtmb = NaN ',
     &            ((ta - tamin) * (tamax - ta) *
     &            (s * x * y - y ** 2 * amp2 - aml2 * aly))
            sqrtmb = 0.d0
         ENDIF

         z1 = (y * sxp + ta * (s * sx + ap2 * y) - ap *
     &        COS(phirad) * sqrtmb) / aly
         z2 = (y * sxp + ta * (x * sx - ap2 * y) - ap *
     &        COS(phirad) * sqrtmb) / aly
         bb = 1.d0!/sqly/pi
         bi12 = bb / (z1 * z2)
         bi1pi2 = bb / z2 + bb / z1
         bis = (bb / z2 ** 2 + bb / z1 ** 2) * aml2
         bir = (bb / z2 ** 2 - bb / z1 ** 2) * aml2
         b1i = bb * z1
         b11i = bb * z1 ** 2
      ENDIF
      hi2 = bis - ym * bi12
      amh2 = amhh ** 2
      zh = zdif
      vvp = (vv10 + vv20) / 2.d0
      vvm = (vv10 - vv20) / 2.d0
      tm(1,1) = 4.d0 * y * hi2
      tm(1,2) = 4.d0 * ta * hi2
      tm(1,3) =  - 2.d0 * (bi12 * ta ** 2 + 2.d0 * bb)
      tm(2,1) = 2.d0 * hi2 * (s * x - amp2 * y)
      tm(2,2) = (bi1pi2 * sx * sxp - bi12 * ta * sxp ** 2 + 2.d0 *
     &          bir * sxp + 2.d0 * hi2 * (sx - 2 * amp2 * ta)) / 2.d0
      tm(2,3) = (bi12 * ta * (2.d0 * amp2 * ta - sx) - bi1pi2 * sxp +
     &          4.d0 * amp2 * bb) / 2.d0
      tm(3,1) = 2.d0 * hi2 * (vv10 * vv20 - amh2 * y)
      tm(3,2) = -2.d0 * ((amh2 * ta - amu * vvm) * hi2 - bir * amu *
     &          vvp - bi1pi2 * vvm * vvp + bi12 * ta * vvp ** 2)
      tm(3,3) = bi12 * ta * (amh2 * ta - amu * vvm) - bi1pi2 * amu *
     &          vvp + 2.d0 * amh2 * bb
      tm(4,1) = -2.d0 * (vv10 * sx - 2.d0 * s * vvp + y * sx * zh) *
     &          hi2
      tm(4,2) = -2.d0 * bi12 * ta * vvp * sxp + bi1pi2 *
     &          (sxp * vvm + sx * vvp) + bir *
     &          (amu * sxp + 2.d0 * vvp) + hi2 *
     &          (sx * (amu - 2.d0 * ta * zh) + 2.d0 * vvm)
      tm(4,3) = bi12 * ta * (sx * (ta * zh - amu / 2d0) - vvm) -
     &          bi1pi2 * (amu / 2d0 * sxp + vvp) + 2.d0 * sx * zh * bb
     
      RETURN
      END



CDECK  ID>,  RV2LN.
      DOUBLE PRECISION FUNCTION rv2ln(taln)
C Integrand (over ta)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'epsil.inc'
      INCLUDE 'amf2.inc'
      DOUBLE PRECISION taln, podinl, ta, costk, sintk, d2kvir, phka
      DOUBLE PRECISION phkb,factor, rmin, rmax, res, rv, vyv(0:1000)
      DOUBLE PRECISION DSIMPS, aval, bval
      INTEGER i
      LOGICAL nonzero
      EXTERNAL podinl
      
      ta = EXP(taln) - y / sx
      taa = ta
      costk = (sx - ap2 * ta) / sqly
      IF (ABS(costk) .LE. 1.d0) THEN
         sintk = SQRT(1.d0 - costk ** 2)
      ELSE
         PRINT *, 'rv2ln: costk>1 ', costk
      sintk = 0.d0
      ENDIF
      d2kvir = (ehad - plh * costk - pth * sintk *
     &         COS(phirad - phidif)) / amp
      phka = 0.5d0 * (ehad - plh * costk) / amp
      phkb = 0.5d0 * ( - pth * sintk) / amp
      daa = d2kvir
      factor = 1.d0 + ta - d2kvir
      CALL tails(ta, tm, d2kvir)

      IF (ita .EQ. 1) THEN
         CALL strf(0.d0, 0.d0, 0.d0, sfm0)
         rmin = 1.d-8
         rmax = (p22 - amc2) / factor
      
C Check minimum r
C        print*, 'rmin1: ', rmin, rmax
C        nonzero = .false.
C        DO i = 1, 4
C        IF(sfm0(i).ne.0.d0) nonzero = .true.
C        enddo
C        IF(nonzero) THEN
C           IF((ap*d2kvir-2.d0*ehad).gt.0.d0) THEN
C              aval = (2.d0*ehad*sx-ap*(amhad(ivec)**2-y-tdif))/
C    &         (ap*d2kvir-2.d0*ehad)
C              bval = pph**2*2.d0*sqrt(sx**2+(ap*y)**2)/
C    &         (ap*d2kvir-2.d0*ehad)
C              rmin = aval-bval
C              rmin = (abs(aval)+abs(bval))*1.d-7
C           ELSE
C              aval = (2.d0*ehad*sx-ap*(amhad(ivec)**2-y-tdif))
C              bval = pph**2*2.d0*sqrt(sx**2+(ap*y)**2)
C              rmin = (abs(aval)+abs(bval))*1.d-7
C           ENDIF
C        ENDIF
C        PRINT *, 'rmin2: ', rmin, rmax
      
      
C        CALL qvnc8(podinl, rmin, rmax, 1.d-4, 1.d-9, res, er, nn, fl, 3500)
C        CALL dqn32(rmin, rmax, podinl, res)

         CALL simpux(rmin, rmax, 100, epsrr, podinl, res)

C        PRINT *, 'start integration '
C        DO i = 0, 100
C           rv = rmin+dble(i)*(rmax-rmin)/1.d2
C           print*, 'start integration: ', i
C           vyv(i) = podinl(rv)
C           print*, rv, vyv(i)
C        ENDDO
C        res = DSIMPS(vyv, rmin, rmax, 100)
C        PRINT *, '-----------------------------', res
      
C        PRINT *, '###################INTEGRATE !'
      
      ELSEIF (ita .EQ. 2) THEN
         res = podinl((p22 - amp2) / factor) / factor
C     ELSE
C        aa = amt/amp
C        IF (ita .EQ. 3) aa = 1.
C        res = podinl((sx-y/aa)/(1.d0+ta/aa))/(1.+ta/aa) /aa**2
      ENDIF
      rv2ln = res * (y / sx + ta)
C     PRINT *, rv2ln, res, rmin, rmax, epsrr
      RETURN
      END



CDECK  ID>,  PODINL.
      DOUBLE PRECISION FUNCTION podinl(r)
C Integrand (over r )
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'amf2.inc'
      DOUBLE PRECISION r, pp, pres, sfm(8), ta, d2kvir, rm
      INTEGER isf, irr, i
      logical reget

      ta = taa
      d2kvir = daa
      CALL strf(ta, d2kvir, r, sfm)
      podinl = 0.d0
      
C Check lower bound
      IF (ita .EQ. 1)  THEN
         rm = r
 1010    reget = .false.
         DO i = 1, 4
            IF (sfm(i) .EQ. 0.d0 .AND.
     &            sfm0(i) .NE. 0.d0 .AND.
     &            r .LT. 1.d-6) reget = .true.
            IF (sfm(i) .EQ. 0.d0 .AND.
     &            sfm0(i) .NE. 0.d0 .AND.
     &            r .LT. 1.d-6) PRINT *, i, sfm(i), sfm0(i)
         ENDDO
         IF (reget) THEN
            rm = rm + 1.d-8
            PRINT *, 'regetting ', rm, r
            CALL strf(ta, d2kvir, rm, sfm)
            STOP
            GOTO 1010
         ENDIF
      ENDIF
      
C     PRINT *, 'podinl start ', isf1, isf2, isf3, i1(isf)+i2(isf)-1
      
      DO 11 isf = isf1, isf2, isf3
         DO 1 irr = 1, (i1(isf) + i2(isf) - 1)
            pp = sfm(isf)
            IF (irr .EQ. 1 .AND. ita .EQ. 1) pp = pp - sfm0(isf) *
     &            (1.d0+ r * ta / y) ** 2
            pres = pp * r ** (irr - 2) / (y + r * ta) ** 2
            podinl = podinl - tm(isf, irr) * pres
C           IF (irr .EQ. 1) PRINT *, 'pod: ', isf, irr, sfm0(isf), sfm(isf), pp, pres, r**(irr-2)
    1    CONTINUE
   11 CONTINUE
C     PRINT *, 'podinl: ', r, podinl
      RETURN
      END



CDECK  ID>,  STRF.
      SUBROUTINE strf(ta, d2kvir, rr, sfm)
C The programm calculates deep inelastic (ita = 1), 
C elastic (ita = 2),  quasielastic (ita = 3) structure FUNCTIONs
C in kinematical point (ta, rr).
C    rr = sx-tt, 
C    ta = (t-y)/rr, 
C where tt = t+amf2-amp2,  amf2 is invarint mass of final hadrons
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'print.inc'
      INCLUDE 'epsmarch.inc'
      DOUBLE PRECISION ta, d2kvir, rr, sfm(8), tldq2, tldtd, tldaly
      DOUBLE PRECISION tldsqly, phq, tldnu, aks, zh, tldsq, tldplh, b1
      DOUBLE PRECISION b2, b3, b4, epspt2, epsnu, epst, epsphq, epsq
      DOUBLE PRECISION epspl, tldpt2, H1z, H2z, H3z, H4z, epsi, aa, h1
      DOUBLE PRECISION h2, h3, h4, dum, tldp22, a, tlde, tldy
      INTEGER i
     
      DO i = 1, 8
         sfm(i) = 0.d0
      ENDDO

      tldq2 = y + rr * ta
      tldtd = tdif - rr * (ta - d2kvir)
      tldaly = (sx - rr) ** 2 + 4.d0 * amp2 * tldq2
      tldsqly = SQRT(tldaly)
      tldp22 = p22 - rr * (1.d0 + ta - d2kvir)
      phq = (amhh ** 2 - tldq2 - tldtd) / 2.d0

C     PRINT *, '---------------Begin-----------------'
C     WRITE (*, '(3f22.19)') sx, rr, ap
      
      tldnu = (sx - rr) / ap
      aks = tldq2 / ap / tldnu
      zh = ehad / tldnu
      tldsq = SQRT(tldq2 + tldnu ** 2)
C     tldplh = (tldtd+tldq2-amhad(ivec)**2+2.d0*tldnu*ehad)/2.d0/tldsq
      tldplh = (ehad * tldnu - phq) / tldsq
      tldpt2 = pph ** 2 - tldplh ** 2
      
      epsnu = SQRT((sx * epsmarch / ap) ** 2 +
     &        (rr * epsmarch / ap) ** 2 +
     &        ((sx - rr) / ap * epsmarch) ** 2)

      epst = SQRT((tdif * epsmarch) ** 2 +
     &       ((ta - d2kvir) * rr * epsmarch) ** 2 +
     &       (rr * ta * epsmarch) ** 2 + (rr * d2kvir * epsmarch) ** 2)

      epsphq = SQRT((2.d0 * amhh ** 2 * epsmarch) ** 2 +
     &         (tldq2 * epsmarch) ** 2 + epst ** 2) / 2.d0

      epsq = 1.d0 / 2.d0 / SQRT(tldq2 + tldnu ** 2) *
     &       SQRT((tldq2 * epsmarch) ** 2 +
     &       (2.d0 * tldnu ** 2 * epsnu) ** 2)

      epspl = SQRT((ehad * epsmarch * tldnu / tldsq) ** 2 +
     &        (ehad * epsnu / tldsq) ** 2 + (epsphq / tldsq) ** 2 +
     &        ((ehad * tldnu - phq) / (tldsq ** 2) * epsq) ** 2)

      epspt2 = 2.d0 * SQRT(pph ** 4 * epsmarch ** 2 + tldplh ** 4 *
     &                     epspl ** 2)
C     PRINT *, 'pt2 =  ', tldpt2, epspt2, (tldpt2**2-epspt2**2), epspl, tldplh, pph
C     IF (tldpt2 .LT. 0.d0) RETURN
      IF (tldpt2 .LT. 0.d0 .AND.
     &         (tldpt2 ** 2 - epspt2 ** 2) .GT. 0.d0) RETURN
C     b1 = 0.d0
C     b2 = 0.d0
C     b3 = 0.d0
C     b4 = 0.d0
C     PRINT *, 'q2', tldq2
C Call semi-inclusive model (H_i defined in Mulders)
      
C     dum = tldq2/ap/aks
C     IF(tldnu.ne.dum) THEN
C        print*, 'redefine'
C        aks = tldq2/ap/tldnu
C        zh = ehad/tldnu
C        tldsq = sqrt(tldq2+tldnu**2)
C        tldplh = (ehad*tldnu-phq)/tldsq
C        tldpt2 = pph**2-tldplh**2
C        IF(tldpt2.lt.0.d0) RETURN
C        WRITE (*, '(10f22.19)') tldq2, aks, zh, tldpt, tldnu, Ehad, pph, 
C    &         dum, (tldq2/ap/tldnu), (zh*tldnu)
C     ENDIF
      
C     WRITE (*, '(10f22.19)') tldq2, aks, zh, tldpt, tldnu, Ehad, pph, 
C    &      dum, (tldq2/ap/tldnu), (zh*tldnu)
      
C Recover y
      a = s / ap * (s / ap - anu) * tldq2 / y
      tlde = (tldnu + SQRT(tldnu ** 2 + 4.d0 * a)) / 2.d0
      tldy = tldnu / tlde
      CALL semi_inclusive_model(tldq2, aks, tldy, zh, tldpt2, tldp22,
     &                          tldplh, H1z, H2z, H3z, H4z)
C     PRINT *, '----------------------------------------------'
C     PRINT *, 'Hs ', tldq2, aks, zh, tldpt, H1z, H2z, H3z, H4z
C     PRINT *, 'semi-inclusive'
C     PRINT *, Eb, tldq2, aks, zh, sqrt(tldpt2)
C     PRINT *, H1z, H2z, H3z, H4z
     
C No photon emission (k-->0)
C     aa = sx * (zdif - 2. * amp * plh / sqly) / 2.d0 / amp2
C     h(1) = zdif*(2.d0*y*xs*sfm0(1)-pth**2*sfm0(4))/amp/y/xs
C     h(2) = 2.*zdif*(2.*xs*(xs*sfm0(2)+aa*sfm0(3))*aly+
C    &sfm0(4)*(aa**2*aly-2.*pth**2*y))/amp/xs/y/aly
C     h(3) = 2.*sfm0(4)*zdif/amp/y/xs
C     h4 = -2.*zdif*(xs*sfm0(3)+aa*sfm0(4))/amp/y/xs
      
C Including kinematic shift due to photon emission
      epsi = ap2 / (sx - rr)
      aa = (sx - rr) * (zh - 2.d0 * amp * tldplh / tldsqly) /
     &2.d0 / amp2
      h1 = zh * (2.d0 * tldq2 * aks * h1z - tldpt2 * h4z) /
     &amp/tldq2/aks
      h2 = 2.d0 * zh * (2.d0 * aks * (aks * h2z + aa * h3z) *
     &tldaly + h4z * (aa ** 2 * tldaly - 2.d0 * tldpt2 * tldq2)) /
     &amp / aks / tldq2 / tldaly
      h3 = 2.d0 * h4z * zh / amp / tldq2 / aks
      h4 = -2.d0 * zh * (aks * h3z + aa * h4z) / amp / tldq2 / aks
C     PRINT *, epsi, aa, h1, h2, h3, h4
      sfm(1) = un * h1
      sfm(2) = un * h2
      sfm(3) = un * h3
      sfm(4) = un * h4
C     sfm(5) = epsi**2*b1
C     sfm(6) = epsi**3*(b2/3.d0+b3+b4)
C     sfm(7) = epsi*(b2/3.d0-b3)
C     sfm(8) = epsi**2*(b2/3.d0-b4)
      
      
C     PRINT *, 'strf: ', tldq2, aks, zh, tldpt2, tldpt, H1z, H2z, H3z, H4z, 
C    &tldnu, sx, ap, ehAD
      
      RETURN
      END
