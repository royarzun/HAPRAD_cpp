      DOUBLE PRECISION FUNCTION h4(x, q2, z)
      IMPLICIT NONE

      DOUBLE PRECISION x, q2, z
      DOUBLE PRECISION q0, lambda
      DOUBLE PRECISION a, a1, a2, b1, b2, bb

      DATA q0, lambda/1.d0, 0.25d0/
      DATA a/0.10908D-02/, a1/-0.35265D-06/, a2/0.30276D-07/
      DATA b1/-0.66787d0/, b2/3.5868d0/, bb/6.8777d0/

      IF (q2 .GT. q0) THEN
         h4 = a * x ** a1 * (1.d0 - x) ** a2 * z ** b1 *
     &        (1.d0 - z) ** b2 * (LOG(q2 / (lambda ** 2)) /
     &                            LOG(q0 /( lambda ** 2))) ** (bb / x)
      ELSE
         h4 = a * x ** a1 * (1.d0 - x) ** a2 * z ** b1 *
     &        (1.d0 - z) ** b2
      ENDIF
      RETURN
      END
