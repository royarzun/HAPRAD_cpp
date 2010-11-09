      SUBROUTINE init_pdf(g, s)
      IMPLICIT NONE
      DOUBLE PRECISION val(20)
      INTEGER g, s
      CHARACTER*20 parm(20)

C Init PDFSET
      parm(1) = 'Init0'
      val(1)  = 0.d + 0
      CALL PDFSET(parm,val)

C Choise of the distribution
      parm(1) = 'Nptype'
      val(1)  = 1.d + 0
      parm(2) = 'Ngroup'
      val(2)  = FLOAT(g)
      parm(3) = 'Nset'
      val(3)  = FLOAT(s)
      CALL PDFSET(parm,val)
      RETURN
      END
