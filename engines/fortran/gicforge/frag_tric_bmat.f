C ORACLE GICForge fragment-coordinate analytic Wilson B rows.
C
C This file is strict Fortran77-style code.  It is the shared
C Fortran implementation for fragment-center and TRIC/geomeTRIC-style
C inter-fragment rotation coordinates.  It is meant to be called by the
C imported legacy GICForge MkBNew path instead of duplicating fragment
C derivatives locally.
C
C Public entry points use one-based atom indices:
C   ORCFCDI : center-center distance B row
C   ORCFCAD : fragment-center / atom distance B row
C   ORCFTRN : center-center Cartesian translation component B row
C   ORCFROT : exponential-map fragment rotation component B row
C   ORCGSPC : classify ORACLE special protected primitive families
C   ORCGSEL : protected-first modified Gram-Schmidt rank selection
C   ORCGLCB : frozen GIC linear-combination B row
C
      Subroutine ORCGSPC(FAMILY,ISPEC)
      Character*(*) FAMILY
      Integer ISPEC
C
      ISPEC=0
      If(FAMILY.eq.'FRAG_DISTANCE') ISPEC=1
      If(FAMILY.eq.'FRAG_CENTER_ATOM_DISTANCE') ISPEC=1
      If(FAMILY.eq.'FRAG_TRANSLATION') ISPEC=1
      If(FAMILY.eq.'FRAG_ORIENTATION') ISPEC=1
      If(FAMILY.eq.'CENTER_ATOM_DISTANCE') ISPEC=1
      Return
      End

      Subroutine ORCGSEL(NCAND,NCOL,BMAT,IPROT,TARGET,TOL,ISEL,
     $  NSEL,RANK,Q,WORK,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NCAND,NCOL,IPROT(*),TARGET,ISEL(*),NSEL,RANK,IFAIL
      Dimension BMAT(NCOL,*),Q(NCOL,*),WORK(*)
C
C     Protected-first non-redundant reduction.  BMAT(:,I) is the
C     analytic B row for candidate I.  IPROT(I).ne.0 marks
C     CLASS=SPECIAL_PROTECTED.  IFAIL=2 means that independent
C     protected rows alone exceed TARGET.
C
      NSEL=0
      RANK=0
      IFAIL=0
      If(TARGET.lt.0) then
       IFAIL=1
       Return
      EndIf
C
      Do 10 IC=1,NCAND
       If(IPROT(IC).eq.0) Goto 10
       If(RANK.ge.TARGET) then
        Call ORCGIND(NCOL,BMAT(1,IC),RANK,Q,TOL,WORK,IADD)
        If(IADD.ne.0) then
         IFAIL=2
         Return
        EndIf
        Goto 10
       EndIf
       Call ORCGADD(NCOL,BMAT(1,IC),RANK,Q,TOL,WORK,IADD)
       If(IADD.ne.0) then
        NSEL=NSEL+1
        ISEL(NSEL)=IC
       EndIf
   10 Continue
C
      If(RANK.ge.TARGET) Return
C
      Do 20 IC=1,NCAND
       If(IPROT(IC).ne.0) Goto 20
       If(RANK.ge.TARGET) Goto 30
       Call ORCGADD(NCOL,BMAT(1,IC),RANK,Q,TOL,WORK,IADD)
       If(IADD.ne.0) then
        NSEL=NSEL+1
        ISEL(NSEL)=IC
       EndIf
   20 Continue
   30 Continue
      If(RANK.ne.TARGET) IFAIL=3
      Return
      End

      Subroutine ORCGADD(NCOL,ROW,RANK,Q,TOL,WORK,IADD)
      Implicit Real*8(A-H,O-Z)
      Integer NCOL,RANK,IADD
      Dimension ROW(*),Q(NCOL,*),WORK(*)
C
      Call ORCGIND(NCOL,ROW,RANK,Q,TOL,WORK,IADD)
      If(IADD.eq.0) Return
      RANK=RANK+1
      RN=DSQRT(ORCGDOT(NCOL,WORK,WORK))
      Do 10 J=1,NCOL
       Q(J,RANK)=WORK(J)/RN
   10 Continue
      Return
      End

      Subroutine ORCGIND(NCOL,ROW,RANK,Q,TOL,WORK,IADD)
      Implicit Real*8(A-H,O-Z)
      Integer NCOL,RANK,IADD
      Dimension ROW(*),Q(NCOL,*),WORK(*)
C
      IADD=0
      Do 10 J=1,NCOL
       WORK(J)=ROW(J)
   10 Continue
      RN=DSQRT(ORCGDOT(NCOL,WORK,WORK))
      If(RN.le.TOL) Return
      Do 20 J=1,NCOL
       WORK(J)=WORK(J)/RN
   20 Continue
      Do 30 K=1,RANK
       PROJ=0.0D0
       Do 40 J=1,NCOL
        PROJ=PROJ+WORK(J)*Q(J,K)
   40  Continue
       Do 50 J=1,NCOL
        WORK(J)=WORK(J)-PROJ*Q(J,K)
   50  Continue
   30 Continue
      RN=DSQRT(ORCGDOT(NCOL,WORK,WORK))
      If(RN.le.TOL) Return
      IADD=1
      Return
      End

      Double Precision Function ORCGDOT(N,A,B)
      Implicit Real*8(A-H,O-Z)
      Integer N
      Dimension A(*),B(*)
C
      ORCGDOT=0.0D0
      Do 10 I=1,N
       ORCGDOT=ORCGDOT+A(I)*B(I)
   10 Continue
      Return
      End

      Subroutine ORCGLCB(NCOL,NTERM,IPRIM,COEF,BMAT,BROW,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NCOL,NTERM,IPRIM(*),IFAIL
      Dimension COEF(*),BMAT(NCOL,*),BROW(*)
C
C     Build the B row for a frozen ORACLE GIC represented as
C     COEFFS=primitive_id:coefficient.  This is the Fortran-side mirror
C     of Python build_gic_b_matrix for local symmetrized GICs.
C
      IFAIL=0
      Call ORCFCLR(NCOL,BROW)
      If(NCOL.le.0.or.NTERM.lt.0) then
       IFAIL=1
       Return
      EndIf
      Do 20 IT=1,NTERM
       IP=IPRIM(IT)
       If(IP.le.0) then
        IFAIL=2
        Return
       EndIf
       Do 10 J=1,NCOL
        BROW(J)=BROW(J)+COEF(IT)*BMAT(J,IP)
   10  Continue
   20 Continue
      Return
      End

      Subroutine ORCFTRN(NAT,NF,FAT,NR,RAT,MODE,BROW)
      Implicit Real*8(A-H,O-Z)
      Integer NAT,NF,NR,MODE,FAT(*),RAT(*)
      Dimension BROW(*)
C
      Call ORCFCLR(3*NAT,BROW)
      If(MODE.lt.1.or.MODE.gt.3) Return
      WF=1.0D0/DFloat(NF)
      WR=1.0D0/DFloat(NR)
      Do 10 I=1,NF
       IAT=FAT(I)
       IND=3*(IAT-1)+MODE
       BROW(IND)=BROW(IND)+WF
   10 Continue
      Do 20 I=1,NR
       IAT=RAT(I)
       IND=3*(IAT-1)+MODE
       BROW(IND)=BROW(IND)-WR
   20 Continue
      Return
      End

      Subroutine ORCFCDI(NAT,NF,FAT,NR,RAT,C,BROW,VALUE,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NAT,NF,NR,FAT(*),RAT(*)
      Dimension C(3,*),BROW(*),CF(3),CR(3),DEL(3),U(3)
C
      IFAIL=0
      Call ORCFCLR(3*NAT,BROW)
      Call ORCFCTR(NF,FAT,C,CF)
      Call ORCFCTR(NR,RAT,C,CR)
      Do 10 I=1,3
       DEL(I)=CF(I)-CR(I)
   10 Continue
      VALUE=DSQRT(ORCFDOT(DEL,DEL))
      If(VALUE.le.1.0D-10) then
       IFAIL=1
       Return
      EndIf
      Do 20 I=1,3
       U(I)=DEL(I)/VALUE
   20 Continue
      Do 30 I=1,NF
       IAT=FAT(I)
       Do 40 J=1,3
        BROW(3*(IAT-1)+J)=BROW(3*(IAT-1)+J)+U(J)/DFloat(NF)
   40  Continue
   30 Continue
      Do 50 I=1,NR
       IAT=RAT(I)
       Do 60 J=1,3
        BROW(3*(IAT-1)+J)=BROW(3*(IAT-1)+J)-U(J)/DFloat(NR)
   60  Continue
   50 Continue
      Return
      End

      Subroutine ORCFCAD(NAT,NF,FAT,IATOM,C,BROW,VALUE,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NAT,NF,IATOM,FAT(*)
      Dimension C(3,*),BROW(*),CF(3),DEL(3),U(3)
C
      IFAIL=0
      Call ORCFCLR(3*NAT,BROW)
      Call ORCFCTR(NF,FAT,C,CF)
      Do 10 I=1,3
       DEL(I)=CF(I)-C(I,IATOM)
   10 Continue
      VALUE=DSQRT(ORCFDOT(DEL,DEL))
      If(VALUE.le.1.0D-10) then
       IFAIL=1
       Return
      EndIf
      Do 20 I=1,3
       U(I)=DEL(I)/VALUE
   20 Continue
      Do 30 I=1,NF
       IAT=FAT(I)
       Do 40 J=1,3
        BROW(3*(IAT-1)+J)=BROW(3*(IAT-1)+J)+U(J)/DFloat(NF)
   40  Continue
   30 Continue
      Do 50 J=1,3
       BROW(3*(IATOM-1)+J)=BROW(3*(IATOM-1)+J)-U(J)
   50 Continue
      Return
      End

      Subroutine ORCFROT(NAT,NF,FAT,NR,RAT,PF,QF,PR,QR,MODE,C,
     $  BROW,VALUE,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NAT,NF,NR,PF,QF,PR,QR,MODE,FAT(*),RAT(*)
      Dimension C(3,*),BROW(*)
      Dimension FF(3,3),FR(3,3),DFF(3,3),DFR(3,3)
      Dimension R(3,3),DR(3,3)
C
      IFAIL=0
      Call ORCFCLR(3*NAT,BROW)
      If(MODE.lt.1.or.MODE.gt.3) then
       IFAIL=1
       Return
      EndIf
      Call ORCFRAM(NAT,NF,FAT,PF,QF,C,0,0,FF,DFF,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFRAM(NAT,NR,RAT,PR,QR,C,0,0,FR,DFR,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFRMAT(FF,FR,DFF,DFR,R,DR)
      Call ORCFQEX(R,DR,MODE,VALUE,DZERO,IFAIL)
      If(IFAIL.ne.0) Return
C
      Do 100 IAT=1,NAT
       Do 110 IAX=1,3
        Call ORCFRAM(NAT,NF,FAT,PF,QF,C,IAT,IAX,FF,DFF,IFAIL)
        If(IFAIL.ne.0) Return
        Call ORCFRAM(NAT,NR,RAT,PR,QR,C,IAT,IAX,FR,DFR,IFAIL)
        If(IFAIL.ne.0) Return
        Call ORCFRMAT(FF,FR,DFF,DFR,R,DR)
        Call ORCFQEX(R,DR,MODE,VTMP,DVAL,IFAIL)
        If(IFAIL.ne.0) Return
        BROW(3*(IAT-1)+IAX)=DVAL
  110  Continue
  100 Continue
      Return
      End

      Subroutine ORCFRMAT(FF,FR,DFF,DFR,R,DR)
      Implicit Real*8(A-H,O-Z)
      Dimension FF(3,3),FR(3,3),DFF(3,3),DFR(3,3)
      Dimension R(3,3),DR(3,3)
C
      Do 10 I=1,3
       Do 20 J=1,3
        R(I,J)=0.0D0
        DR(I,J)=0.0D0
        Do 30 K=1,3
         R(I,J)=R(I,J)+FF(K,I)*FR(K,J)
         DR(I,J)=DR(I,J)+DFF(K,I)*FR(K,J)
     $    +FF(K,I)*DFR(K,J)
   30   Continue
   20  Continue
   10 Continue
      Return
      End

      Subroutine ORCFRAM(NAT,NF,FAT,IP,IQ,C,IDAT,IDAX,FR,DFR,
     $  IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer NAT,NF,IP,IQ,IDAT,IDAX,FAT(*)
      Dimension C(3,*),FR(3,3),DFR(3,3)
      Dimension CEN(3),DCEN(3),PV(3),DPV(3),QV(3),DQV(3)
      Dimension PAX(3),DPAX(3),QRAW(3),DQRAW(3),QAX(3),DQAX(3)
      Dimension SRAW(3),DSRAW(3),SAX(3),DSAX(3),T1(3),T2(3)
C
      IFAIL=0
      Call ORCFCTR(NF,FAT,C,CEN)
      Call ORCFDCT(NF,FAT,IDAT,IDAX,DCEN)
      Do 10 I=1,3
       PV(I)=C(I,IP)-CEN(I)
       QV(I)=C(I,IQ)-CEN(I)
       DPV(I)=-DCEN(I)
       DQV(I)=-DCEN(I)
   10 Continue
      If(IP.eq.IDAT.and.IDAX.ge.1.and.IDAX.le.3) then
       DPV(IDAX)=DPV(IDAX)+1.0D0
      EndIf
      If(IQ.eq.IDAT.and.IDAX.ge.1.and.IDAX.le.3) then
       DQV(IDAX)=DQV(IDAX)+1.0D0
      EndIf
      Call ORCFUNI(PV,PAX,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFUDE(PV,DPV,PAX,DPAX,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFCRS(PAX,QV,QRAW)
      Call ORCFCRS(DPAX,QV,T1)
      Call ORCFCRS(PAX,DQV,T2)
      Do 20 I=1,3
       DQRAW(I)=T1(I)+T2(I)
   20 Continue
      Call ORCFUNI(QRAW,QAX,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFUDE(QRAW,DQRAW,QAX,DQAX,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFCRS(PAX,QAX,SRAW)
      Call ORCFCRS(DPAX,QAX,T1)
      Call ORCFCRS(PAX,DQAX,T2)
      Do 30 I=1,3
       DSRAW(I)=T1(I)+T2(I)
   30 Continue
      Call ORCFUNI(SRAW,SAX,IFAIL)
      If(IFAIL.ne.0) Return
      Call ORCFUDE(SRAW,DSRAW,SAX,DSAX,IFAIL)
      If(IFAIL.ne.0) Return
      Do 40 I=1,3
       FR(I,1)=PAX(I)
       FR(I,2)=QAX(I)
       FR(I,3)=SAX(I)
       DFR(I,1)=DPAX(I)
       DFR(I,2)=DQAX(I)
       DFR(I,3)=DSAX(I)
   40 Continue
      Return
      End

      Subroutine ORCFQEX(R,DR,MODE,VALUE,DVALUE,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Integer MODE
      Double Precision R,DR,VALUE,DVALUE,TR,DTR,KW,DKW
      Double Precision NUM,DNUM,K,DK,KN2,DKN2,KN,DKN
      Double Precision ANG,DEN,DANG,FAC,DFAC
      Dimension R(3,3),DR(3,3),K(3),DK(3)
C
      IFAIL=0
      TR=R(1,1)+R(2,2)+R(3,3)
      DTR=DR(1,1)+DR(2,2)+DR(3,3)
      If(TR.le.-1.0D0+1.0D-10) then
       IFAIL=1
       Return
      EndIf
      KW=0.5D0*DSQRT(TR+1.0D0)
      If(KW.le.1.0D-12) then
       IFAIL=1
       Return
      EndIf
      DKW=DTR/(8.0D0*KW)
C
      NUM=R(2,3)-R(3,2)
      DNUM=DR(2,3)-DR(3,2)
      K(1)=NUM/(4.0D0*KW)
      DK(1)=DNUM/(4.0D0*KW)-NUM*DKW/(4.0D0*KW*KW)
      NUM=R(3,1)-R(1,3)
      DNUM=DR(3,1)-DR(1,3)
      K(2)=NUM/(4.0D0*KW)
      DK(2)=DNUM/(4.0D0*KW)-NUM*DKW/(4.0D0*KW*KW)
      NUM=R(1,2)-R(2,1)
      DNUM=DR(1,2)-DR(2,1)
      K(3)=NUM/(4.0D0*KW)
      DK(3)=DNUM/(4.0D0*KW)-NUM*DKW/(4.0D0*KW*KW)
C
      KN2=K(1)*K(1)+K(2)*K(2)+K(3)*K(3)
      DKN2=2.0D0*(K(1)*DK(1)+K(2)*DK(2)+K(3)*DK(3))
      If(KN2.le.1.0D-20) then
       VALUE=2.0D0*K(MODE)
       DVALUE=2.0D0*DK(MODE)
       Return
      EndIf
      KN=DSQRT(KN2)
      DKN=DKN2/(2.0D0*KN)
      ANG=2.0D0*DATAN2(KN,KW)
      DEN=KN*KN+KW*KW
      DANG=2.0D0*(KW*DKN-KN*DKW)/DEN
      FAC=ANG/KN
      DFAC=DANG/KN-ANG*DKN/(KN*KN)
      VALUE=FAC*K(MODE)
      DVALUE=FAC*DK(MODE)+K(MODE)*DFAC
      Return
      End

      Subroutine ORCFCTR(NF,FAT,C,CEN)
      Implicit Real*8(A-H,O-Z)
      Integer NF,FAT(*)
      Dimension C(3,*),CEN(3)
C
      Call ORCFCLR(3,CEN)
      Do 10 I=1,NF
       IAT=FAT(I)
       Do 20 J=1,3
        CEN(J)=CEN(J)+C(J,IAT)
   20  Continue
   10 Continue
      Do 30 J=1,3
       CEN(J)=CEN(J)/DFloat(NF)
   30 Continue
      Return
      End

      Subroutine ORCFDCT(NF,FAT,IDAT,IDAX,DCEN)
      Implicit Real*8(A-H,O-Z)
      Integer NF,IDAT,IDAX,FAT(*)
      Dimension DCEN(3)
C
      Call ORCFCLR(3,DCEN)
      If(IDAT.le.0.or.IDAX.lt.1.or.IDAX.gt.3) Return
      Do 10 I=1,NF
       If(FAT(I).eq.IDAT) DCEN(IDAX)=DCEN(IDAX)+1.0D0/DFloat(NF)
   10 Continue
      Return
      End

      Subroutine ORCFUNI(V,U,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Dimension V(3),U(3)
C
      IFAIL=0
      DN=DSQRT(ORCFDOT(V,V))
      If(DN.le.1.0D-10) then
       IFAIL=1
       Return
      EndIf
      Do 10 I=1,3
       U(I)=V(I)/DN
   10 Continue
      Return
      End

      Subroutine ORCFUDE(V,DV,U,DU,IFAIL)
      Implicit Real*8(A-H,O-Z)
      Dimension V(3),DV(3),U(3),DU(3)
C
      IFAIL=0
      DN=DSQRT(ORCFDOT(V,V))
      If(DN.le.1.0D-10) then
       IFAIL=1
       Return
      EndIf
      PROJ=ORCFDOT(U,DV)
      Do 10 I=1,3
       DU(I)=(DV(I)-U(I)*PROJ)/DN
   10 Continue
      Return
      End

      Subroutine ORCFCRS(A,B,C)
      Implicit Real*8(A-H,O-Z)
      Dimension A(3),B(3),C(3)
C
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      Return
      End

      Double Precision Function ORCFDOT(A,B)
      Implicit Real*8(A-H,O-Z)
      Dimension A(3),B(3)
C
      ORCFDOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      Return
      End

      Subroutine ORCFCLR(N,A)
      Implicit Real*8(A-H,O-Z)
      Integer N
      Dimension A(*)
C
      Do 10 I=1,N
       A(I)=0.0D0
   10 Continue
      Return
      End
