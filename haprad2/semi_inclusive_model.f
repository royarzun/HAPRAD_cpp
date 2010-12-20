      SUBROUTINE semi_inclusive_model(q2, x, y, z, pt2, mx2, pl,
     &                                H1z, H2z, H3z, H4z)
c     Semi - inclusive cross section from the product of parton
c     distribution and fragmentation functions

      IMPLICIT NONE

      INCLUDE 'partons.inc'
      INCLUDE 'constants8.inc'

      DOUBLE PRECISION q2, x, y, z, pt2, mx2, pl
      DOUBLE PRECISION XD, ZD, Q2D, xv, zv, q2v
      DOUBLE PRECISION aa, H1z, H2z, H3z, H4z
      DOUBLE PRECISION SCALE, UPV, DNV, USEA, DSEA, STR
      DOUBLE PRECISION CHM, BOT, TOP, GL
      DOUBLE PRECISION uff(2), dff(2), sff(2), cff(2), bff(2), gff
      DOUBLE PRECISION bs, GTMD, H1, H2, H3m, H4m, nu, Eh, ph, ph_long
      DOUBLE PRECISION rt, rtz, Rb, igev2mb, sgmpt, rlt
      DOUBLE PRECISION uq, dq, sq, cq, bq, tq, gg
      DOUBLE PRECISION ac, bc, cc, dc, ec, fc
      DOUBLE PRECISION h3, h4, qm, t, pi_thresh, r, xi
      DOUBLE PRECISION Ebeam, cterm, m_cos_phi, m_cos_2phi

      INTEGER ih, nc, GPDF, SPDF, IFINI, ISET, ICHARGE

      DATA ac/1.2025D - 10/, bc/ - 5.2703D - 02/, cc/3.7467D - 01/
      DATA dc/6.5397D - 02/, ec/ - 2.2136D - 01/, fc/ - 1.0621D - 01/
      DATA rlt/0.14d0/

      COMMON /FRAGINI/ IFINI

c NLO PDFs
c      DATA GPDF/5/, SPDF/7/  ! GRV 94 HO NLL DIS
c      DATA GPDF/3/, SPDF/36/ ! MRS H NLL DIS
c      DATA GPDF/4/, SPDF/28/ ! CTEQ 2pD NLL DIS
c      DATA GPDF/2/, SPDF/8/  ! DFLM 260 NLL DIS

c LO PDFs
      DATA GPDF/5/, SPDF/5/   ! GRV 94 LO
c     DATA GPDF/3/, SPDF/72/  ! MRST c - g LO
c     DATA GPDF/4/, SPDF/32/  ! CTEQ 4L LO
c     DATA GPDF/2/, SPDF/5/   ! DFLM ca LO

c Pi+ FF
      DATA ISET/1/, ICHARGE/1/
      DATA nc/0/
      DATA Rb/1.4142135624d0/
      DATA igev2mb/389.379292d0/

c Init
      H1  = 0.d + 0
      H2  = 0.d + 0
      H3m = 0.d + 0
      H4m = 0.d + 0
      H1z = 0.d + 0
      H2z = 0.d + 0
      H3z = 0.d + 0
      H4z = 0.d + 0

c Check kinematics
      if (x .le. 0.d0 .or. x .ge. 1.d0) return
      if (z .le. 0.d0 .or. z .ge. 1.d0) return
      nc = nc + 1
c     print *, nc

c Obtain Parton Distribution Functions
c     XD = x
      r = sqrt(1.d0 + (2.d0 * mp * x) ** 2 / q2)
      xi = 2.d0 * x / (1.d0 + r)
      XD = xi

c QCD scale
      if(q2 .gt. 1.d0) then
         SCALE = sqrt(q2)
      else
         SCALE = 1.0000001
      endif

c Init PDFs
      if (nc .eq. 1) CALL init_pdf(GPDF,SPDF)
      CALL STRUCTM(XD,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)

c Obtain Parton Fragmentation Functions
      if (z .gt. 0.01) then
         ZD = z
      else
         ZD = 0.01000001
      endif

      if(q2 .gt. 1.d0) then
         Q2D = q2
      else
         Q2D = 1.0000001
      endif

c Init FFs
      IFINI = 1
      if (nc .eq. 1) IFINI = 0
      CALL PKHFF(ISET,ICHARGE,ZD,Q2D,uff,dff,sff,cff,bff,gff)

c Parton transverse momentum distribution (Gaussian fit)
      sgmpt = ac + bc * x + cc * z + dc * x ** 2 + ec * z ** 2 +
     &        fc * x * z

      if(sgmpt .lt. 0.02) sgmpt = 0.02
      if(sgmpt .gt. 0.15) sgmpt = 0.15

      if(pl .gt. 0.0) then
         GTMD = exp( - pt2 / (2.d0 * sgmpt)) / (2.d0 * pi * sgmpt)
      else
         GTMD = exp( - (pt2 + 2.d0 * pl ** 2) / (2.d0 * sgmpt)) /
     &          (2.d0 * pi * sgmpt)
      endif

c Kinematics
c     nu = q2 / (2.d0 * mp * x)
c     qm = sqrt(q2 + nu ** 2)
c     Eh = z * nu
c     if(Eh .lt. mpi) return
c     ph = sqrt(Eh ** 2 - mpi ** 2)
c     if(ph .lt. pt) return
c     ph_long = sqrt(ph ** 2 - pt ** 2)
      
      
c     write(*, '(10f22.19)') q2, x, z, pt, nu, Eh, ph,
c    &                       (q2 / (2.d0 * mp * x)),
c    &                       (q2 / (2.d0 * mp * nu)),
c    &                       (z * nu)
      
      
c     t = mpi ** 2 - q2 - 2.d0 * nu ** 2 * z + 2.d0 * ph_long * qm
c     mx2 = mp ** 2 + 2.d0 * mp * nu * (1.d0 - z) + t
      if(mx2 .lt. ((mp + mpi) ** 2)) return
      
c      print *, mx2, ((mp + mpi) ** 2)
      
      pi_thresh = sqrt(1.d0 - (mp + mpi) ** 2 / mx2)

c Semi - inclusive structure functions at LO and LT
      uq = eu2 * ((UPV + USEA) * uff(1) + USEA * uff(2))
      dq = ed2 * ((DNV + DSEA) * dff(1) + DSEA * dff(2))

      sq = es2 * (STR * sff(1) + STR * sff(2))
      cq = ec2 * (CHM * cff(1) + CHM * cff(2))
      bq = eb2 * (BOT * bff(1) + BOT * bff(2))

      tq = 0.0d + 0
      gg = 0.0d + 0

      H2 = (uq + dq + sq + cq + bq + tq + gg) * GTMD * pi_thresh
      H1 = H2 / (2.d0 * x * (1.d0 + rlt)) *
     &          (1.d0 + 4.d0 * mp ** 2 * x ** 2 / q2)
c     H1 = H2 / (2.0 * x)

      xv = x
      if(x .lt. 0.1) xv = 0.1

      q2v = q2
      if(q2 .lt. 1.0) q2v = 1.0

      zv = z
      if(z .lt. 0.1) zv = 0.1

      H3m = h3(x,q2,z) * GTMD * pi_thresh
      H4m = h4(x,q2,z) * GTMD * pi_thresh
      
c Check that <cos(phi)> and <cos(2phi)> are < 1
      Ebeam = q2 / (2.d0 * mp * x * y)

      rt = 1.d0 - y-mp * x * y / (2.d0 * Ebeam)
      rtz = sqrt(rt / (1.d0 + 2.d0 * mp * x / (y * Ebeam)))
      cterm = x * y ** 2 * H1 + rt * H2
      m_cos_phi = sqrt(pt2 / q2) * (2.d0 - y) *
     &            rtz * H3m / (2.d0 * cterm)

c     if(Ebeam .eq. 27.d0) print * ,'H3,H4 ',H3m,H4m
c     if(Ebeam .eq. 27.d0) print * ,'<cos(phi)> ',m_cos_phi

      if (abs(4.d0 * m_cos_phi) .gt. 0.9d0) stop
      if (abs(4.d0 * m_cos_phi) .gt. 0.9d0) then
         H3m = 0.9d0 * sign(1.d0,m_cos_phi) * (2.d0 * cterm) /
     &         (sqrt(pt2 / q2) * (2.d0 - y) * rtz) / 4.d0
c        if(Ebeam .eq. 27.d0) print *, 'Correct H3= ', H3m
      endif
      m_cos_2phi = pt2 / q2 * rtz ** 2 * H4m / (2.d0 * cterm)
c     if(Ebeam .eq. 27.d0) print *, '<cos(phi)> ', m_cos_2phi

      if (abs(4.d0 * m_cos_2phi) .gt. 0.9d0) then
         H4m = 0.9d0 * sign(1.d0,m_cos_2phi) * (2.d0 * cterm) /
     &         (pt2 / q2 * rtz ** 2) / 4.d0
c        if (Ebeam .eq. 27.d0) print *, 'Correct H4= ', H4m
      endif

c     if(Ebeam .eq. 27.d0) then
c        print *, 'final ', (x * y ** 2 * H1 + rt * H2),
c    &                      (sqrt(pt2 / q2) * (2.d0 - y) * rtz * H3m),
c    &                      (pt2 / q2 * rtz ** 2 * H4m)
c     endif
      
c     H1 = 0.
c     H2 = 0.
c     H3m = 0.
c     H4m = 0.

c Keep in hadronic tensor of Levlet & Mulders, PRD49, 96 (1994)
      H1z = H1
      H2z = H2
      H3z = H3m
      H4z = H4m

c Convert into new hadronic tensor of Akushevich, Eur.Phys.J. C10, 681 (1999)
c     aa = (Eh * qm - ph_long * nu) / (mp * qm)
c     H1z = 2.d0 * z / mp * (H1 - pt ** 2 / (q2 * 2.d0 * x) * H4m)
c     H2z = 2.d0 * z / (nu * mp ** 2) *
c    &            (H2 + aa / x * H3m + nu ** 2 / q2 *
c    &               (2.d0 * aa ** 2 * mp ** 2 / q2 - (pt / qm) ** 2) *
c    &               H4m)
c     H3z = 4.d0 * z * nu / (q2 ** 2) * H4m
c     H4z = - 2.0 * z / (mp * q2) * (H3m + aa / x * H4m)

c     print *, 'new kinem: ', q2, x, z, pt
c     print *, H1z, H2z, H3z, H4z, GTMD, pi_thresh, UPV, uff(1)
                                                          
      RETURN
      END
