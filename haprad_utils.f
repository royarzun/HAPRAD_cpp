      SUBROUTINE simpu(a1, b1, h1, reps1, aeps1, funct, x, ai, aih,
     &                 aiabs)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
      DIMENSION f(7), p(5)
      h = SIGN(h1, b1 - a1)
      s = SIGN(1.d0, h)
      a = a1
      b = b1
      ai = 0.d0
      aih = 0.d0
      aiabs = 0.d0
      p(2) = 4.d0
      p(4) = 4.d0
      p(3) = 2.d0
      p(5) = 1.d0
      IF (b - a) 1, 2, 1
    1 reps = ABS(reps1)
      aeps = ABS(aeps1)
      DO 3 k = 1, 7
    3 f(k) = 10.d16
      x = a
      c = 0.d0
      f(1) = funct(x) / 3.d0
C     PRINT *, 'c ', x, f(1)
    4 x0 = x
      IF ((x0 + 4.d0 * h - b) * s) 5, 5, 6
    6 h = (b - x0) / 4.d0
      IF(h) 7, 2, 7
    7 DO 8 k = 2, 7
    8 f(k) = 10.d16
      c = 1.d0
    5 di2 = f(1)
      di3 = ABS(f(1))
      DO 9 k = 2, 5
         x = x + h
         IF ((x - b) * s) 23, 24, 24
   24    x = b
   23    IF (f(k) - 10.d16) 10, 11, 10
   11    f(k) = funct(x) / 3.d0
C        PRINT *, 'd'
C        WRITE (*, '(2f19.15)') x, f(k)
   10    di2 = di2 + p(k) * f(k)
    9 di3 = di3 + p(k) * ABS(f(k))
      di1 = (f(1) + 4.d0 * f(3) + f(5)) * 2.d0 * h
      di2 = di2 * h
      di3 = di3 * h
      IF (reps) 12, 13, 12
   13 IF (aeps) 12, 14, 12
   12 eps = ABS((aiabs + di3) * reps)
      IF (eps - aeps) 15, 16, 16
   15 eps = aeps
   16 delta = ABS(di2 - di1)
      IF (delta - eps) 20, 21, 21
   20 IF (delta - eps / 8.d0) 17, 14, 14
   17 h = 2.d0 * h
      f(1) = f(5)
      f(2) = f(6)
      f(3) = f(7)
      DO 19 k = 4, 7
   19 f(k) = 10.d16
      GO TO 18
   14 f(1) = f(5)
      f(3) = f(6)
      f(5) = f(7)
      f(2) = 10.d16
      f(4) = 10.d16
      f(6) = 10.d16
      f(7) = 10.d16
   18 di1 = di2 + (di2 - di1) / 15.d0
      ai = ai + di1
      aih = aih + di2
      aiabs = aiabs + di3
      GO TO 22
   21 h = h / 2.d0
      f(7) = f(5)
      f(6) = f(4)
      f(5) = f(3)
      f(3) = f(2)
      f(2) = 10.d16
      f(4) = 10.d16
      x = x0
      c = 0.d0
      GO TO 5
   22 IF (c) 2, 4, 2
    2 RETURN
      END



CDECK  ID>, SIMPS.
      SUBROUTINE simpsz(a1, b1, h1, reps1, aeps1, funct, x, ai, aih,
     &                  aiabs)
C simps
C
C a1, b1       - the limits of integration
C h1            - an initial step of integration
C reps1, aeps1 - relative AND absolute PRECISION of integration
C funct         - a name of FUNCTION subprogram for calculation of integrand +
C x             - an argument of the integrand
C ai            - the value of integral
C aih           - the value of integral with the step of integration
C aiabs         - the value of integral for module of the integrand
C
C This subrogram calculates the definite integral with the relative or
C absolute PRECISION by simpson's method with the automatical choice of
C the step of integration.
C If aeps1 is very small (like 1 * E-17), THEN calculation of integral
C with reps1, AND IF reps1 is very small (like 1 * E-10), THEN calculation
C of integral with aeps1. When aeps1 = reps1 = 0, THEN calculation with
C the constant step h1.
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      DIMENSION f(7), p(5)
      h = SIGN(h1, b1 - a1)
      s = SIGN(1.d0, h)
      a = a1
      b = b1
      ai = 0.d0
      aih = 0.d0
      aiabs = 0.d0
      p(2) = 4.d0
      p(4) = 4.d0
      p(3) = 2.d0
      p(5) = 1.d0
      IF (b - a) 1, 2, 1
    1 reps = ABS(reps1)
      aeps = ABS(aeps1)
      DO 3 k = 1, 7
    3 f(k) = 10.d16
      x = a
      c = 0.d0
      f(1) = funct(x) / 3.d0
    4 x0 = x
      IF ((x0 + 4.d0 * h - b) * s) 5, 5, 6
    6 h = (b - x0) / 4.d0
      IF(h) 7, 2, 7
    7 DO 8 k = 2, 7
    8 f(k) = 10.d16
      c = 1.d0
    5 di2 = f(1)
      di3 = ABS(f(1))
      DO 9 k = 2, 5
         x = x+h
         IF ((x - b) * s) 23, 24, 24
   24    x = b
   23    IF (f(k) - 10.d16) 10, 11, 10
   11    f(k) = funct(x) / 3.d0
   10    di2 = di2 + p(k) * f(k)
    9 di3 = di3 + p(k) * ABS(f(k))
      di1 = (f(1) + 4.d0 * f(3) + f(5)) * 2.d0 * h
      di2 = di2 * h
      di3 = di3 * h
      IF (reps) 12, 13, 12
   13 IF (aeps) 12, 14, 12
   12 eps = ABS((aiabs + di3) * reps)
      IF (eps - aeps) 15, 16, 16
   15 eps = aeps
   16 delta = ABS(di2-di1)
      IF (delta - eps) 20, 21, 21
   20 IF (delta - eps / 8.d0) 17, 14, 14
   17 h = 2.d0 * h
      f(1) = f(5)
      f(2) = f(6)
      f(3) = f(7)
      DO 19 k = 4, 7
   19 f(k) = 10.d16
      GO TO 18
   14 f(1) = f(5)
      f(3) = f(6)
      f(5) = f(7)
      f(2) = 10.d16
      f(4) = 10.d16
      f(6) = 10.d16
      f(7) = 10.d16
   18 di1 = di2 + (di2 - di1) / 15.d0
      ai = ai + di1
      aih = aih + di2
      aiabs = aiabs + di3
      GO TO 22
   21 h = h/2.d0
      f(7) = f(5)
      f(6) = f(4)
      f(5) = f(3)
      f(3) = f(2)
      f(2) = 10.d16
      f(4) = 10.d16
      x = x0
      c = 0.d0
      GO TO 5
   22 IF(c) 2, 4, 2
    2 RETURN
      END



      SUBROUTINE simpsx(a, b, np, ep, func, res)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      EXTERNAL func
      step = (b - a) / DBLE(np)
      CALL simpsz(a, b, step, ep, 1d-18, func, ra, res, r2, r3)
      RETURN
      END



      SUBROUTINE simptx(a, b, np, ep, func, res)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      EXTERNAL func
      step = (b - a) / dble(np)
      CALL simpt(a, b, step, ep, 1d-18, func, ra, res, r2, r3)
      RETURN
      END



      SUBROUTINE simpux(a, b, np, ep, func, res)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      EXTERNAL func
      step = (b - a) / dble(np)
      CALL simpu(a, b, step, ep, 1d-18, func, ra, res, r2, r3)
      RETURN
      END



      SUBROUTINE simpt(a1, b1, h1, reps1, aeps1, funct, x, ai, aih,
     &                 aiabs)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
      dimension f(7), p(5)
      h = sign(h1, b1 - a1)
      s = sign(1.d0, h)
      a = a1
      b = b1
      ai = 0.d0
      aih = 0.d0
      aiabs = 0.d0
      p(2) = 4.d0
      p(4) = 4.d0
      p(3) = 2.d0
      p(5) = 1.d0
      IF (b - a) 1, 2, 1
    1 reps = ABS(reps1)
      aeps = ABS(aeps1)
      DO 3 k = 1, 7
    3 f(k) = 10.d16
      x = a
      c = 0.d0
      f(1) = funct(x) / 3.d0
C     PRINT *, 'a ', x, f(1)
      
    4 x0 = x
      IF ((x0 + 4.d0 * h - b) * s) 5, 5, 6
    6 h = (b - x0) / 4.d0
      IF (h) 7, 2, 7
    7 DO 8 k = 2, 7
    8 f(k) = 10.d16
      c = 1.d0
    5 di2 = f(1)
      di3 = ABS(f(1))
      DO 9 k = 2, 5
         x = x + h
         IF ((x - b) * s) 23, 24, 24
   24    x = b
   23    IF (f(k) - 10.d16) 10, 11, 10
   11    f(k) = funct(x) / 3.d0
C        PRINT *, 'b ', x, f(k)
   10    di2 = di2 + p(k) * f(k)
    9 di3 = di3 + p(k) * ABS(f(k))
      di1 = (f(1) + 4.d0 * f(3) + f(5)) * 2.d0 * h
      di2 = di2 * h
      di3 = di3 * h
      IF (reps) 12, 13, 12
   13 IF (aeps) 12, 14, 12
   12 eps = ABS((aiabs + di3) * reps)
      IF (eps - aeps) 15, 16, 16
   15 eps = aeps
   16 delta = ABS(di2 - di1)
      IF (delta - eps) 20, 21, 21
   20 IF (delta - eps / 8.d0) 17, 14, 14
   17 h = 2.d0 * h
      f(1) = f(5)
      f(2) = f(6)
      f(3) = f(7)
      DO 19 k = 4, 7
   19 f(k) = 10.d16
      GO TO 18
   14 f(1) = f(5)
      f(3) = f(6)
      f(5) = f(7)
      f(2) = 10.d16
      f(4) = 10.d16
      f(6) = 10.d16
      f(7) = 10.d16
   18 di1 = di2 + (di2 - di1) / 15.d0
      ai = ai + di1
      aih = aih + di2
      aiabs = aiabs + di3
      GO TO 22
   21 h = h / 2.d0
      f(7) = f(5)
      f(6) = f(4)
      f(5) = f(3)
      f(3) = f(2)
      f(2) = 10.d16
      f(4) = 10.d16
      x = x0
      c = 0.d0
      GO TO 5
   22 IF (c) 2, 4, 2
    2 RETURN
      END



CDECK  ID>, D01FCE.
      SUBROUTINE d01fce(ndim, a, b, minpts, maxpts, functn, eps, 
     * acc, lenwrk, wrkstr, finval, ifail)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
C Mark 8 Release. NAG Copyright 1979.
C
C Adaptive multidimensional integration SUBROUTINE
C
C Input parameters:
C
C ndim   - INTEGER number of variables, must exceed 1 but not exceed 15.
C a      - real array of lower limits, with dimension ndim
C b      - real array of upper limits, with dimension ndim
C minpts - INTEGER minimum number of integrand values to be allowed, 
C          which must not exceed maxpts.
C maxpts - INTEGER maximum number of integrand values to be
C          allowed, which must be at least 2 ** ndim + 2 * ndim ** 2 + 2 * ndim + 1.
C functn - externally declared user defined real FUNCTION integrand.
C          It must have parameters (ndim, z), where z is a real array
C          of dimension ndim.
C eps    - real required relative accuracy, must be greater than zero
C lenwrk - INTEGER length of array wrkstr, must be at least 2 * ndim + 4.
C ifail  - INTEGER nag failure parameter
C             ifail = 0 for hard fail
C             ifail = 1 for soft fail
C
C Output parameters:
C
C minpts - INTEGER number of integrand values used by the routine
C wrkstr - real array of working storage of dimension (lenwrk).
C acc    - real estimated relative accuracy of finval
C finval - real estimated value of integral
C ifail  - ifail = 0 for normal exit, when estimated relative less
C              integaccuracy rand values used.
C          ifail = 1 IF ndim .lt. 2, ndim .gt. 15, minpts .gt. maxpts, 
C              maxpts .lt. 2**ndim+2*ndim*(ndim+1)+1, eps.LE.0
C              or lenwrk .lt. 2*ndim+4.
C          ifail = 2 IF maxpts was too small for d01fce to obtain
C              the required relative accuracy eps. In this case
C              d01fce RETURNs a value of finval with estimated
C              relative accuracy acc.
C          ifail = 3 IF lenwrk too small for maxpts integrand
C              values. In this case d01fce RETURNs a value of
C              finval with estimated accuracy acc using the working
C              storage available, but acc will be greater than eps.
C
C ********************************************************************
C
C Scalar arguments ..
      DOUBLE PRECISION eps, finval, acc
      INTEGER ifail, lenwrk, maxpts, minpts, ndim
C Aray arguments ..
      DOUBLE PRECISION a, b, wrkstr
      DIMENSION a(ndim), b(ndim), wrkstr(lenwrk)
C Function arguments ..
      DOUBLE PRECISION functn
C Local scalars ..
      CHARACTER*8 srname
      DOUBLE PRECISION 
     & abserr, df1, df2, difmax, f1, f2, f3, f4, half, lamda2, 
     & lamda4, lamda5, one, ratio, rgncmp, rgnerr, rgnert, rgnval, 
     & rgnvlt, rgnvol, rlndim, sum1, sum2, sum3, sum4, sum5, two, 
     & twondm, weit1, weit2, weit3, weit4, weit5, weitp1, weitp2, 
     & weitp3, weitp4, zero
      INTEGER dvaxes, dvaxis, dvflag, funcls, ierror, j, k, maxaxs, 
     & mxrgns, pointr, rgncls, rulcls, sbrgns, subrgn, subtmp, 
     & tpontp, tpontr
C Local arrays ..
      DIMENSION center(15), dif(15), oldcnt(15), width(15), z(15)
      INTEGER dvcntl(15), dvcntr(15)
C Function references ..
      DOUBLE PRECISION x02aae
      INTEGER p01aae, x02bbe

      DATA srname /'  d01fce'/
      DATA zero, one, two, half /0.d0, 1.d0, 2.d0, 0.5d0/

C Subroutine initialisation and parameter checking
      IF (ndim .LT. 2 .OR. ndim .GT. 15) GO TO 560
      IF (minpts .GT. maxpts) GO TO 560
      IF (eps .LE. zero) GO TO 560
      IF (lenwrk .LT. 2 * ndim + 4) GO TO 560
      funcls = 0
      finval = zero
      abserr = zero
      twondm = two ** ndim
      rgnvol = twondm
      dvflag = 1
      fffff1 = DBLE(x02bbe(one))
      fffff2 = 1.d0 / x02aae(0.0d0)
      maxaxs = INT(MIN(fffff1, fffff2))
C     maxaxs = int(amin1(float(x02bbe(one)), 1.0 / x02aae(0.0d0)))
      maxaxs = (maxaxs - ndim) / (ndim + 1)
      mxrgns = lenwrk / (2 * ndim + 4)
      sbrgns = 0
      rgnvlt = zero
      rgnert = zero
      DO 20 j = 1, ndim
         center(j) = (a(j) + b(j)) * half
         dif(j) = zero
         width(j) = (b(j) - a(j)) * half
         dvcntl(j) = 1
         dvcntr(j) = 1
         oldcnt(j) = center(j)
         rgnvol = rgnvol * width(j)
   20 CONTINUE
C End subroutine initialisation

C Basic rule initialisation
      rulcls = 2 ** ndim + 2 * ndim * ndim + 2 * ndim + 1
      funcls = rulcls
      IF (maxpts .LT. rulcls) GO TO 560
      rlndim = ndim
      lamda2 = SQRT(9.d0 / 70.d0)
      lamda4 = SQRT(9.d0 / 10.d0)
      lamda5 = SQRT(9.d0 / 19.d0)
      weit1 = (12824.d0 - 9120.d0 * rlndim + 400.d0 * rlndim * rlndim) /
     &        19683.d0
      weit2 = 980.d0/6561.d0
      weit3 = (1820.d0 - 400.d0 * rlndim) / 19683.d0
      weit4 = 200.d0/19683.d0
      weit5 = 6859.d0 / 19683.d0 / twondm
      weitp1 = (729.d0 - 950.d0 * rlndim + 50.d0 * rlndim ** 2) / 729.d0
      weitp2 = 245.d0 / 486.d0
      weitp3 = (265.d0 - 100.d0 * rlndim) / 1458.d0
      weitp4 = 25.d0/729.d0
      ratio = (lamda2 / lamda4) ** 2
C End basic rule initialisation

      GO TO 100

C Divide subregion with largest error AND prepare to use basic rule on
C each portion
   40 subrgn = 1
      pointr = wrkstr(1)
      rgncls = rulcls
      rgnvol = twondm
      tpontr = pointr + 2
      DO 60 j = 1, ndim
         tpontr = tpontr + 2
         center(j) = wrkstr(tpontr - 1)
         width(j) = wrkstr(tpontr)
         dvcntr(j) = 1
         dvcntl(j) = 1
         oldcnt(j) = center(j)
         rgnvol = rgnvol * width(j)
   60 CONTINUE
      dvaxes = wrkstr(pointr + 2)
      IF (dvaxes .LT. 0) GO TO 600
   80 dvaxis = dvaxes
      dvaxes = dvaxis / (ndim + 1)
      dvaxis = dvaxis  -  (ndim + 1) * dvaxes
      dvcntl(dvaxis) = 2 * dvcntl(dvaxis)
      rgncls = rgncls * 2
      IF (dvaxes .GT. 0) GO TO 80
      IF (funcls + rgncls .GT. maxpts) GO TO 580
      IF (rgncls / rulcls + sbrgns - 1 .GT. mxrgns) dvflag = 2
      funcls = funcls + rgncls
C     PRINT *, funcls
      abserr = abserr - wrkstr(pointr)
      finval = finval  -  wrkstr(pointr + 1)

C Begin basic rule
  100 DO 120 j = 1, ndim
         z(j) = center(j)
  120 CONTINUE
      sum1 = functn(ndim, z)
      sum2 = zero
      sum3 = zero
      DO 140 j = 1, ndim
         z(j) = center(j) - lamda2 * width(j)
         f1 = functn(ndim, z)
         z(j) = center(j) + lamda2 * width(j)
         f2 = functn(ndim, z)
         z(j) = center(j) -  lamda4 * width(j)
         f3 = functn(ndim, z)
         z(j) = center(j) + lamda4*width(j)
         f4 = functn(ndim, z)
         sum2 = sum2 + f1 + f2
         sum3 = sum3 + f3 + f4
         df1 = f1 + f2 - two * sum1
         df2 = f3 + f4 - two * sum1
         dif(j) = dif(j) + ABS(df1 - ratio * df2)
         z(j) = center(j)
  140 CONTINUE
      sum4 = zero
      DO 200 j = 2, ndim
         z(j-1) = center(j-1) - lamda4 * width(j-1)
         DO 160 k = j, ndim
            z(k) = center(k) - lamda4 * width(k)
            sum4 = sum4 + functn(ndim, z)
            z(k) = center(k) + lamda4 * width(k)
            sum4 = sum4 + functn(ndim, z)
            z(k) = center(k)
  160    CONTINUE
         z(j-1) = center(j-1) + lamda4 * width(j-1)
         DO 180 k = j, ndim
            z(k) = center(k) - lamda4 * width(k)
            sum4 = sum4 + functn(ndim, z)
            z(k) = center(k) + lamda4 * width(k)
            sum4 = sum4 + functn(ndim, z)
            z(k) = center(k)
  180    CONTINUE
         z(j-1) = center(j-1)
  200 CONTINUE
      sum5 = zero
      DO 220 j = 1, ndim
         z(j) = center(j) - lamda5 * width(j)
  220 CONTINUE
  240 DO 260 j = 2, ndim
         IF (z(j-1) .LT. center(j-1) + width(j-1)) GO TO 280
         z(j-1) = center(j-1) - lamda5 * width(j-1)
         z(j) = z(j) + two * lamda5 * width(j)
  260 CONTINUE
      IF (z(ndim) .GT. center(ndim) + width(ndim)) GO TO 300
  280 sum5 = sum5 + functn(ndim, z)
      z(1) = z(1) + two * lamda5 * width(1)
      GO TO 240
  300 rgnval = rgnvol * (weit1 * sum1 + weit2 * sum2 + weit3 * sum3 +
     &                   weit4 * sum4 + weit5 * sum5)
      rgncmp = rgnvol * (weitp1 * sum1 + weitp2 * sum2 +
     &                   weitp3 * sum3 + weitp4 * sum4)
      rgnerr = ABS(rgnval - rgncmp)
C End basic rule

C Store results of basic rule application
      rgnvlt = rgnvlt + rgnval
      rgnert = rgnert + rgnerr
      finval = finval + rgnval
      abserr = abserr + rgnerr
      IF (dvflag .EQ. 0) GO TO 340
      IF (dvflag .EQ. 2) GO TO 500
      pointr = mxrgns + sbrgns * (2 * ndim + 3) + 1
      sbrgns = sbrgns + 1
      wrkstr(sbrgns) = pointr
      subrgn = sbrgns
      tpontr = pointr + 2
      DO 320 j = 1, ndim
         tpontr = tpontr + 2
         wrkstr(tpontr-1) = center(j)
         wrkstr(tpontr) = width(j)
  320 CONTINUE
  340 wrkstr(pointr) = rgnert
      wrkstr(pointr+1) = rgnvlt
C Determine axis along which fourth difference is largest
      difmax = zero
      DO 380 j = 1, ndim
         IF (difmax .GT. dif(j)) GO TO 360
         difmax = dif(j)
         dvaxis = j
  360    dif(j) = zero
  380 CONTINUE
      tpontr = pointr + 2 * (dvaxis + 1)
      wrkstr(tpontr) = width(dvaxis) * half
      wrkstr(tpontr-1) = center(dvaxis) - wrkstr(tpontr)
      IF (dvflag .NE. 2) GO TO 400
      dvaxes = wrkstr(pointr + 2)
      IF (dvaxes .GT. maxaxs) dvaxes = -1
      dvaxis = dvaxis + (ndim + 1) * dvaxes
  400 wrkstr(pointr+2) = dvaxis
      IF (dvflag .EQ. 1) GO TO 460
C Determine the position in the parially ordered list of the subregion
C which replaces most recently divided subregion
  420 subtmp = 2 * subrgn
      IF (subtmp .GT. sbrgns) GO TO 480
      tpontr = wrkstr(subtmp)
      IF (subtmp .EQ. sbrgns) GO TO 440
      tpontp = wrkstr(subtmp + 1)
      IF (wrkstr(tpontr) .GE. wrkstr(tpontp)) GO TO 440
      subtmp = subtmp + 1
      tpontr = tpontp
  440 IF (rgnert .GE. wrkstr(tpontr)) GO TO 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      GO TO 420
C When working storage is not used up, determine the position in the
C partially ordered list for the description of other portion(s) of most
C recently divided subregion
  460 subtmp = subrgn/2
      IF (subtmp .lt. 1) GO TO 480
      tpontr = wrkstr(subtmp)
      IF (rgnert.LE.wrkstr(tpontr)) GO TO 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      GO TO 460
  480 rgnvlt = zero
      rgnert = zero
      IF (dvflag .EQ. 2) GO TO 540
      dvflag = 1 - dvflag
C Count to determine the next part of the recently divided subregion for
C application of the basic rule
  500 center(1) = center(1) + two * width(1)
      dvcntr(1) = dvcntr(1) + 1
      DO 520 j = 2, ndim
         IF (dvcntr(j-1) .LE. dvcntl(j-1)) GO TO 100
         dvcntr(j-1) = 1
         center(j-1) = oldcnt(j-1)
         dvcntr(j) = dvcntr(j) + 1
         center(j) = center(j) + two * width(j)
  520 CONTINUE
      IF (dvcntr(ndim) .LE. dvcntl(ndim)) GO TO 100
      center(ndim) = oldcnt(ndim)
      IF (dvflag .EQ. 2) GO TO 340
C End ordering of basic rule results

C Make checks for possible termination of routine
  540 acc = abserr / ABS(finval)
      IF (acc .GT. eps .OR. funcls .LT. minpts) GO TO 40
C Loop back to apply basic rule
C Termination point, set ifail AND RETURN
      ierror = 0
      GO TO 620
  560 ierror = 1
      GO TO 620
  580 ierror = 2
      GO TO 620
  600 ierror = 3
  620 minpts = funcls
      ifail = p01aae(ifail, ierror, srname)
      RETURN
      END



C Cern library routine e104 (interpolation)
      FUNCTION dfint(narg, arg, nent, ent, table)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      DIMENSION arg(narg), nent(narg), ent(1000), table(100500)
      DIMENSION d(narg), ncomb(narg), ient(narg)
      kd = 1
      m = 1
      ja = 1
      DO 5 i = 1, narg
         ncomb(i) = 1
         jb = ja - 1+nent(i)
         DO 2 j = ja, jb
            IF (arg(i) .LE. ent(j)) GO TO 3
    2    CONTINUE
         j = jb
    3    IF (j .NE. ja) GO TO 4
         j = j + 1
    4    jr = j - 1
         d(i) = (ent(j) - arg(i)) / (ent(j) - ent(jr))
         ient(i) = j - ja
         kd = kd + ient(i) * m
         m = m * nent(i)
    5 ja = jb + 1
      dfint = 0.d0
   10 fac = 1.d0
      iadr = kd
      ifadr = 1
      DO 15 i = 1, narg
         IF (ncomb(i) .EQ. 0) GO TO 12
         fac = fac * (1.d0 - d(i))
         GO TO 15
   12    fac = fac * d(i)
         iadr = iadr - ifadr
   15 ifadr = ifadr * nent(i)
      dfint = dfint + fac * table(iadr)
      il = narg
   40 IF (ncomb(il) .EQ. 0) GO TO 80
      ncomb(il) = 0
      IF (il .EQ. narg) GO TO 10
      il = il + 1
      DO 50 k = il, narg
   50 ncomb(k) = 1
      GO TO 10
   80 il = il - 1
      IF(il .NE. 0) GO TO 40
      RETURN
      END



CDECK  ID>, FSPEN.
      DOUBLE PRECISION FUNCTION fspen(x)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
      DATA f1/1.644934d0/

      IF (x) 8, 1, 1
    1 IF (x - .5d0) 2, 2, 3
    2 fspen = fspens(x)
      RETURN
    3 IF (x - 1.d0) 4, 4, 5
    4 fspen = f1 - log(x) * log(1.d0 - x + 1d-10) - fspens(1.d0 - x)
      RETURN
    5 IF (x - 2.d0) 6, 6, 7
    6 fspen = f1 - 0.5d0 * log(x) * log((x - 1.d0) ** 2 / x) +
     &        fspens(1.d0 - 1.d0 / x)
      RETURN
    7 fspen = 2.d0 * f1 - 0.5d0 * log(x) ** 2 - fspens(1.d0 / x)
      RETURN
    8 IF (x + 1.d0) 10, 9, 9
   9  fspen = -0.5d0 * log(1.d0 - x) ** 2 - fspens(x / (x - 1.d0))
      RETURN
  10  fspen = -0.5d0 * log(1. - x) * log(x ** 2 / (1.d0 - x)) - f1 +
     &        fspens(1.d0 / (1.d0 - x))
      RETURN
      END



CDECK  ID>, FSPENS.
      DOUBLE PRECISION FUNCTION fspens(x)
C Spence function
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
      f = 0.d0
      a = 1.d0
      an = 0.d0
      tch = 1.d-16
    1 an = an + 1.d0
      a = a * x
      b = a / an ** 2
      f = f + b
      IF (b - tch) 2, 2, 1
    2 fspens = f
      RETURN
      END



      DOUBLE PRECISION FUNCTION x02aae(x)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
C NAG Copyright 1975
C Mark 4.5 Release
C if = IBM.
C for ibm/360/370/3090
C     DATA z/z3380000000000000/
C     x02aae = z
C for sun
      DATA z/1.1d-16/
      x02aae = z
C * eps *
C Returns the value eps where eps is the smallest positive number such
C that 1.0 + eps > 1.0
C The x parameter is not used
C for icl 1900
C     x02aae = 2.0 ** (-37.0)
C if = pc.
C for pdp11
C     x02aae = 2.d0 ** (-23.d0)
      RETURN
      END



      INTEGER  FUNCTION x02bbe(x)
      IMPLICIT DOUBLE PRECISION(a-h, o-z)
C NAG Copyright 1975
C Mark 4.5 Release
C Real x
C * maxint *
C Returns the largest INTEGER representable on the computer
C the x parameter is not used
C for icl 1900
C     x02bbe = 8388607
C for ibm, sun, vax, ibm pc/386/486
      x02bbe = 2147483647
C for pdp11
C     x02bbe = 32767
      RETURN
      END



      INTEGER FUNCTION p01aae(ifail, error, srname)
C Mark 1 Release. NAG Copyright 1971
C Mark 3 revised
C Mark 4a revised, ier-45
C Mark 4.5 revised
C Mark 7 revised (dec 1978)
C Returns the value of error or terminates the program.
      INTEGER error, ifail, nout
      CHARACTER*8 srname
C Test if no error detected
      IF (error .EQ. 0) GO TO 20
C Determine output unit for message
      CALL x04aae (0, nout)
C Test for soft failure
      IF (MOD(ifail, 10) .EQ. 1) GO TO 10
C     hard failure
      WRITE (nout, 99999) srname, error
C Stopping mechanism may also differ
      STOP
C Soft fail
C Test IF error messages suppressed
   10 IF (MOD(ifail / 10, 10) .EQ. 0) GO TO 20
      WRITE (nout, 99999) srname, error
   20 p01aae = error
      RETURN
99999 FORMAT (1h0, 38herror detected by nag library routine , a8, 
     &        11h - ifail = , i5//)
      END



      SUBROUTINE x04aae(i, nerr)
C Mark 7 Release. NAG Copyright 1978
C Mark 7c Revised ier-190 (may 1979)
C If i = 0, sets nerr to current error message unit number (stored in
C nerr1).
C If i = 1, changes current error message unit number to value specified
C by nerr.
C
C *** note ***
C this routine assumes that the value of nerr1 is saved between calls.
C in some implementations it may be necessary to store nerr1 in a
C labelled common block /ax04aa/ to achieve this.
C
C Scalar arguments
      INTEGER i, nerr
C Local scalars
      INTEGER nerr1
      DATA nerr1 /5/
      IF (i .EQ. 0) nerr = nerr1
      IF (i .EQ. 1) nerr1 = nerr
      RETURN
      END
