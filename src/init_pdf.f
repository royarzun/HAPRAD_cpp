      subroutine init_pdf(g,s)
      implicit none
      double precision VAL(20)
      integer g,s
      character*20 PARM(20)
c-----Init PDFSET
      PARM(1)='Init0'
      VAL(1)=0.d+0
      call PDFSET(PARM,VAL)
c-----Choise of the distribution
      PARM(1)='Nptype'
      VAL(1)=1.d+0
      PARM(2)='Ngroup'
      VAL(2)=float(g)
      PARM(3)='Nset'
      VAL(3)=float(s)
      call PDFSET(PARM,VAL)
      return
      end

      subroutine exec_structm(XD,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      implicit none
      double precision XD,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
      call STRUCTM(XD,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      return
      end
      
      subroutine exec_pkhff(ISET,ICHARGE,ZD,Q2D,uff,dff,sff,cff,bff,gff)
      implicit none 
      double precision ISET,ICHARGE,ZD,Q2D,uff(3),dff(3),sff(3),cff(3),bff(3),gff(3) 
      call PKHFF(ISET,ICHARGE,ZD,Q2D,uff,dff,sff,cff,bff,gff)
      return
      end
