      DOUBLE PRECISION FUNCTION h3(x,  q2,  z)
      IMPLICIT NONE

      DOUBLE PRECISION x, q2, z
      DOUBLE PRECISION q0, lambda
      DOUBLE PRECISION a, a1, a2, b1, b2, bb

      DATA q0, lambda/1.d0, 0.25d0/
      DATA a/-0.36544D-03/, a1/-2.1855d0/, a2/3.4176d0/
      DATA b1/-1.7567d0/, b2/1.1272d0/, bb/8.9985d0/

      IF (q2 .GT. q0) THEN
         h3 = a * x ** a1 * (1.d0 - x) ** a2 * z ** b1 *
     &        (1.d0-z) ** b2 * (LOG(q2 / (lambda ** 2)) /
     &                          LOG(q0 / (lambda ** 2))) ** bb
      ELSE
         h3 = a * x ** a1 * (1.d0 - x) ** a2 * z ** b1 *
     &        (1.d0 - z) ** b2
      ENDIF
      RETURN
      END
