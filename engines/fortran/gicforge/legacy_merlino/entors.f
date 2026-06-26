C=== CENTRAL / TERMINAL CLASSIFIERS ================================
      INTEGER FUNCTION CENTRAL_CLASS(IAt,JAt,IAn,EAn,IRing,C)
      INTEGER IAt,JAt, IAn(*), IRing(*)
      DOUBLE PRECISION EAn(*), C(3,*)
      DOUBLE PRECISION Distan, RCovCT
      DOUBLE PRECISION Value, Val0, BO
C---- BO stimata dal modello (clamp difensivo)
      Value = Distan(C,JAt,IAt,0)
      Val0  = RCovCT(IAn(JAt),IAn(IAt))
      IF (Value .LE. 0.0D0) Value = 1.0D-03
      IF (Val0  .LE. 0.0D0) Val0  = 1.0D-03
      BO = DEXP((Val0-Value)/3.0D-01)
      IF (BO .GT.  60.0D0) BO = 60.0D0
      IF (BO .LT. 1.0D-06) BO = 1.0D-06
C---- euristica semplice ma robusta
      IF (IRing(IAt).GT.0 .OR. IRing(JAt).GT.0) THEN
         CENTRAL_CLASS = 2
      ELSE IF (BO .GE. 1.6D0) THEN
         CENTRAL_CLASS = 3
      ELSE IF (BO .GE. 1.3D0) THEN
         CENTRAL_CLASS = 2
      ELSE
         CENTRAL_CLASS = 1
      ENDIF
      RETURN
      END

      INTEGER FUNCTION TERM_CLASS(Z)
      INTEGER Z
C---- mapping terminale compatto T1..T6
      IF (Z .EQ. 1) THEN
         TERM_CLASS = 1
      ELSE IF (Z .EQ. 6) THEN
         TERM_CLASS = 2
      ELSE IF (Z .EQ. 7 .OR. Z .EQ. 8 .OR. Z .EQ. 9) THEN
         TERM_CLASS = 4
      ELSE IF (Z .EQ. 16 .OR. Z .EQ. 17) THEN
         TERM_CLASS = 5
      ELSE IF (Z .EQ. 35 .OR. Z .EQ. 53) THEN
         TERM_CLASS = 6
      ELSE
C        default: trattalo come etero EN
         TERM_CLASS = 4
      ENDIF
      RETURN
      END
C=== BASE (CENTRAL) E MODULATORI (TERMINALI) =======================
      SUBROUTINE INIT_KNEW
      DOUBLE PRECISION B1C(6), B2C(6), B3C(6), B6C(6)
      DOUBLE PRECISION S1T(6,6), S2T(6,6)
      DOUBLE PRECISION S3T(6,6), S6T(6,6)
      COMMON /KBASEN/ B1C, B2C, B3C, B6C
      COMMON /KMODT/  S1T, S2T, S3T, S6T
      LOGICAL KNEW_INIT
      COMMON /KNEWFLG/ KNEW_INIT
      INTEGER i, j

C---- base per classi centrali (C1..C6; usiamo 1..5)
      DO 10 i=1,6
         B1C(i)=0.10D0
         B2C(i)=0.06D0
         B3C(i)=0.22D0
         B6C(i)=0.04D0
   10 CONTINUE
      B1C(2)=0.12D0
      B2C(2)=0.14D0
      B3C(2)=0.28D0
      B6C(2)=0.05D0
      B1C(3)=0.00D0
      B2C(3)=0.65D0
      B3C(3)=0.00D0
      B6C(3)=0.00D0
      B1C(4)=0.10D0
      B2C(4)=0.12D0
      B3C(4)=0.24D0
      B6C(4)=0.05D0
      B1C(5)=0.08D0
      B2C(5)=0.10D0
      B3C(5)=0.20D0
      B6C(5)=0.04D0

C---- inizializza modulatori con 1.0
      DO 30 i=1,6
         DO 20 j=1,6
            S1T(i,j)=1.0D0
            S2T(i,j)=1.0D0
            S3T(i,j)=1.0D0
            S6T(i,j)=1.0D0
   20    CONTINUE
   30 CONTINUE

C---- aggiusta S1 (sterico) secondo gerarchie
C     H(1), Csp3(2), Csp2(3 non usato qui), N/O/F(4), S/Cl(5), Br/I(6)
      S1T(1,1)=0.8D0
      S1T(1,2)=0.9D0
      S1T(2,2)=1.1D0
      S1T(2,4)=1.1D0
      S1T(2,5)=1.15D0
      S1T(4,4)=1.10D0
      S1T(5,5)=1.15D0
      S1T(6,6)=1.20D0
C     simmetrizza
      DO 50 i=1,6
         DO 40 j=i+1,6
            S1T(j,i)=S1T(i,j)
   40    CONTINUE
   50 CONTINUE

C---- S3 (periodicita' etano-like): leggero boost su Csp3-Csp3
      S3T(2,2)=1.10D0
      S3T(1,2)=1.05D0
      S3T(1,1)=1.00D0
      DO 70 i=1,6
         DO 60 j=i+1,6
            S3T(j,i)=S3T(i,j)
   60    CONTINUE
   70 CONTINUE

C---- S6: cresce con anisotropia / alogeni pesanti
      S6T(1,1)=0.60D0
      S6T(1,2)=0.70D0
      S6T(5,5)=1.10D0
      S6T(6,6)=1.20D0
      DO 90 i=1,6
         DO 80 j=i+1,6
            S6T(j,i)=S6T(i,j)
   80    CONTINUE
   90 CONTINUE

C---- S2: simile a S1 ma smorzato
      DO 110 i=1,6
         DO 100 j=1,6
            S2T(i,j)=1.0D0 + 0.8D0*(S1T(i,j)-1.0D0)
  100    CONTINUE
  110 CONTINUE

      KNEW_INIT = .TRUE.
      RETURN
      END

C=== SCALA ADDITIVA PER L'EXTRA DI K2 (BO) =========================
      DOUBLE PRECISION FUNCTION SCALE_T2(TX,TY)
      INTEGER TX, TY
      INTEGER a, b
      a = TX
      b = TY
      IF (a .GT. b) THEN
         a = TY
         b = TX
      ENDIF
C---- H-H:1.0 ; H-X:1.5 ; X-X:2.0 ; se presente Br/I -> 2.4
      IF (a.EQ.1 .AND. b.EQ.1) THEN
         SCALE_T2 = 1.0D0
      ELSE IF (a.EQ.1) THEN
         SCALE_T2 = 1.5D0
         IF (b.EQ.6) SCALE_T2 = 2.0D0
      ELSE
         SCALE_T2 = 2.0D0
         IF (a.EQ.6 .OR. b.EQ.6) SCALE_T2 = 2.4D0
      ENDIF
      RETURN
      END
*Deck TorsEn
      SUBROUTINE TORSEN(IAt,JAt,
     & MxBnd,NAtoms,IAn,EAn,IRing,NLen,IAtom,NBond,IBond,C,
     & V0,V1,V2,V3,V6,NSA,NSB,KREF,LREF,IFLAG)

C... argomenti
      INTEGER IAt, JAt, MxBnd, NAtoms, NLen
      INTEGER IAn(NAtoms), IRing(NAtoms)
      DOUBLE PRECISION EAn(NAtoms)
      INTEGER IAtom(2,*), NBond(NAtoms), IBond(MxBnd,NAtoms)
      DOUBLE PRECISION C(3,NAtoms)

C... output
      DOUBLE PRECISION V0, V1, V2, V3, V6
      INTEGER NSA, NSB, KREF, LREF, IFLAG

C... tabelle base/modulatori
      DOUBLE PRECISION B1C(6), B2C(6), B3C(6), B6C(6)
      DOUBLE PRECISION S1T(6,6), S2T(6,6), S3T(6,6), S6T(6,6)
      COMMON /KBASEN/ B1C, B2C, B3C, B6C
      COMMON /KMODT/  S1T, S2T, S3T, S6T
      LOGICAL KNEW_INIT
      COMMON /KNEWFLG/ KNEW_INIT

C... locali interi
      INTEGER i, j, IAUX, JAUX
      INTEGER ASUB(3), BSUB(3)
      INTEGER IREF_A, IREF_B
      INTEGER ISYM_A, ISYM_B, EQMAX_A, EQMAX_B
      INTEGER NGRID, k
      INTEGER TX, TY, ICENT

C... locali double
      DOUBLE PRECISION OFFA(3), OFFB(3)
      DOUBLE PRECISION K1, K2, K3, K6, K2_BASE
      DOUBLE PRECISION Distan, RCovCT
      DOUBLE PRECISION Value, Val0, BndOrd
      DOUBLE PRECISION S, Fmin, Eval, w, PI, EPS
      DOUBLE PRECISION SUMK1, SUMK2, SUMK3
      DOUBLE PRECISION BETH, BOETH, BBIP, BOBIP
      DOUBLE PRECISION XBO, PEXP, BARR

C... locali logici
      LOGICAL NEED6
      LOGICAL A_C3, B_C3, A_C2, B_C2, A_EVEN, B_EVEN

C... esterne
      INTEGER CENTRAL_CLASS, TERM_CLASS
      DOUBLE PRECISION SCALE_T2
      EXTERNAL ANALYZE_SIDE_EAN
      EXTERNAL CENTRAL_CLASS, TERM_CLASS, SCALE_T2
      EXTERNAL INIT_KNEW

C... init
      V0=0.0D0
      V1=0.0D0
      V2=0.0D0
      V3=0.0D0
      V6=0.0D0
      SUMK1=0.0D0
      SUMK2=0.0D0
      SUMK3=0.0D0
      NSA=0
      NSB=0
      KREF=0
      LREF=0
      IFLAG=0
      PI = DACOS(-1.0D0)
      EPS = 1.0D-10

C... esclusioni
      IF (NBond(IAt) .LE. 1 .OR. NBond(JAt) .LE. 1) THEN
         IFLAG = 1
         RETURN
      ENDIF
      IF (IRing(IAt).GT.0 .AND. IRing(JAt).GT.0) THEN
         IFLAG = 2
         RETURN
      ENDIF

C... sostituenti lato A
      DO 100 i=1, NBond(IAt)
         IAUX = IBond(i,IAt)
         IF (IAUX .EQ. JAt) GOTO 100
         IF (NSA .LT. 3) THEN
            NSA = NSA + 1
            ASUB(NSA) = IAUX
         ENDIF
  100 CONTINUE

C... sostituenti lato B
      DO 200 j=1, NBond(JAt)
         JAUX = IBond(j,JAt)
         IF (JAUX .EQ. IAt) GOTO 200
         IF (NSB .LT. 3) THEN
            NSB = NSB + 1
            BSUB(NSB) = JAUX
         ENDIF
  200 CONTINUE

      IF (NSA .LE. 0 .OR. NSB .LE. 0) THEN
         IFLAG = 3
         RETURN
      ENDIF

C... riferimenti e simmetrie da EAn
      CALL ANALYZE_SIDE_EAN(NSA,ASUB,EAn,IREF_A,EQMAX_A,ISYM_A)
      CALL ANALYZE_SIDE_EAN(NSB,BSUB,EAn,IREF_B,EQMAX_B,ISYM_B)
      KREF = ASUB(IREF_A)
      LREF = BSUB(IREF_B)

C... BO stimata e warning
      Value  = Distan(C,JAt,IAt,0)
      Val0   = RCovCT(IAn(JAt),IAn(IAt))
      IF (Value .LE. 0.0D0) Value = 1.0D-03
      IF (Val0  .LE. 0.0D0) Val0  = 1.0D-03
      BndOrd = DEXP((Val0-Value)/3.0D-01)
      IF (BndOrd .GT. 1.1D0) THEN
         WRITE(*,'(A,2I6,A,F6.3)') 'WARNING: high BO for I,J=',
     &     IAt, JAt, '  BO=', BndOrd
      ENDIF

C... offsets
      CALL ASSIGN_OFFSETS(NSA,IREF_A,OFFA)
      CALL ASSIGN_OFFSETS(NSB,IREF_B,OFFB)

C... serve n=6 per 3|2 o 2|3
      NEED6 = ( (NSA.EQ.3 .AND. NSB.EQ.2) .OR.
     &          (NSA.EQ.2 .AND. NSB.EQ.3) )

C... classe centrale + init nuove tabelle
      ICENT = CENTRAL_CLASS(IAt,JAt,IAn,EAn,IRing,C)
      IF (.NOT. KNEW_INIT) CALL INIT_KNEW

C... barriera extra da BO (clamp robusto)
      BETH  = 65.0D0
      BOETH =  2.0D0
      BBIP  =  2.5D0
      BOBIP =  1.3D0
      XBO   = BndOrd - 1.0D0
      IF (XBO .LT. 1.0D-12) XBO = 1.0D-12
      IF (XBO .GT. 50.0D0)  XBO = 50.0D0
      PEXP  = DLOG(BBIP/BETH) / DLOG(BOBIP-1.0D0)
      BARR  = BETH * ( XBO ** PEXP )

C... accumulo Vn (fattorizzato central/terminal)
      DO 500 i=1,NSA
         DO 490 j=1,NSB
            TX = TERM_CLASS( IAn(ASUB(i)) )
            TY = TERM_CLASS( IAn(BSUB(j)) )

C.......... K1,K3,K6
            K1 = B1C(ICENT) * S1T(TX,TY)
            K3 = B3C(ICENT) * S3T(TX,TY)
            IF (NEED6) THEN
               K6 = B6C(ICENT) * S6T(TX,TY)
            ELSE
               K6 = 0.0D0
            ENDIF

C.......... K2 = base + extra(BO)*scala(TX,TY)
            K2_BASE = B2C(ICENT) * S2T(TX,TY)
            K2      = K2_BASE
            IF (BndOrd .GT. 1.1D0) THEN
               K2 = K2 + 0.5D0 * BARR * SCALE_T2(TX,TY)
            ENDIF

C.......... proiezioni
            V1 = V1 + K1*DCOS( OFFA(i) - OFFB(j) )
            V2 = V2 + K2*DCOS( 2.0D0*(OFFA(i)-OFFB(j)) )
            V3 = V3 + K3*DCOS( 3.0D0*(OFFA(i)-OFFB(j)) )
            V6 = V6 + K6*DCOS( 6.0D0*(OFFA(i)-OFFB(j)) )

C.......... somme totali per eventuale fallback
            SUMK1 = SUMK1 + K1
            SUMK2 = SUMK2 + K2
            SUMK3 = SUMK3 + K3
  490    CONTINUE
  500 CONTINUE

C... sanita' numerica: elimina NaN
      IF (V1 .NE. V1) V1 = 0.0D0
      IF (V2 .NE. V2) V2 = 0.0D0
      IF (V3 .NE. V3) V3 = 0.0D0
      IF (V6 .NE. V6) V6 = 0.0D0

C... fallback V1 se lato terminale annulla proiezione
      IF ((NSA.EQ.1 .OR. NSB.EQ.1) .AND. DABS(V1).LT.EPS) THEN
         V1 = SUMK1
      ENDIF

C... SYMMETRY GATING (come prima)
      A_C3  = (NSA.EQ.3 .AND. EQMAX_A.EQ.3)
      B_C3  = (NSB.EQ.3 .AND. EQMAX_B.EQ.3)
      A_C2  = (NSA.EQ.2 .AND. EQMAX_A.EQ.2)
      B_C2  = (NSB.EQ.2 .AND. EQMAX_B.EQ.2)
      A_EVEN = A_C2 .OR. (NSA.EQ.3 .AND. EQMAX_A.EQ.2)
      B_EVEN = B_C2 .OR. (NSB.EQ.3 .AND. EQMAX_B.EQ.2)

      IF ((A_C3 .AND. B_C2) .OR. (B_C3 .AND. A_C2)) THEN
         V1 = 0.0D0
         V2 = 0.0D0
         V3 = 0.0D0
      ELSE IF (A_C3 .AND. B_C3) THEN
         V1 = 0.0D0
         V2 = 0.0D0
      ELSE IF (A_EVEN .OR. B_EVEN) THEN
         V1 = 0.0D0
         V3 = 0.0D0
      ELSE IF ((A_C3 .AND. .NOT. B_C2) .OR.
     &         (B_C3 .AND. .NOT. A_C2)) THEN
         V1 = 0.0D0
         V2 = 0.0D0
      ENDIF

C... V0 per min(E)=0 con griglia
      S    = V1 + V2 + V3 + V6
      NGRID = 720
      Fmin = 1.0D30
      DO 600 k=0, NGRID-1
         w = (2.0D0*PI*DBLE(k))/DBLE(NGRID)
         Eval = V1*DCOS(w)
     &        + V2*DCOS(2.0D0*w)
     &        + V3*DCOS(3.0D0*w)
     &        + V6*DCOS(6.0D0*w)
         IF (Eval .EQ. Eval) THEN
            IF (Eval .LT. Fmin) Fmin = Eval
         ENDIF
  600 CONTINUE
      IF (Fmin .GE. 1.0D29) Fmin = 0.0D0
      V0 = - ( S + Fmin )

      RETURN
      END
*Deck EAN_EQ
      LOGICAL FUNCTION EAN_EQ(A,B)
      DOUBLE PRECISION A, B, ATOL, RTOL, M
      ATOL = 1.0D-04
      RTOL = 1.0D-03
      M = DMAX1(DABS(A),DABS(B))
      EAN_EQ = (DABS(A-B) .LE. ATOL + RTOL*M)
      RETURN
      END
*Deck ANALYZE_SIDE_EAN
      SUBROUTINE ANALYZE_SIDE_EAN(ns,SUB,EAn,IREF,EQMAX,ISYM)
      INTEGER ns, SUB(3)
      DOUBLE PRECISION EAn(*)
      INTEGER IREF, EQMAX, ISYM
      DOUBLE PRECISION EA(3)
      LOGICAL E12, E13, E23, ALL_EQ, TWO_EQ
C     tipizza la funzione come LOGICAL
      LOGICAL EAN_EQ

      IF (ns .LE. 1) THEN
         IREF  = 1
         EQMAX = 0
         ISYM  = 0
         RETURN
      ENDIF

      EA(1) = EAn(SUB(1))
      EA(2) = EAn(SUB(2))
      IF (ns .EQ. 3) EA(3) = EAn(SUB(3))

      IF (ns .EQ. 3) THEN
         E12 = EAN_EQ(EA(1),EA(2))
         E13 = EAN_EQ(EA(1),EA(3))
         E23 = EAN_EQ(EA(2),EA(3))
         ALL_EQ = (E12 .AND. E13 .AND. E23)
         TWO_EQ = (E12 .OR.  E13 .OR.  E23)
         IF (ALL_EQ) THEN
            EQMAX = 3
            ISYM  = 30
            IREF  = 1
         ELSE IF (TWO_EQ) THEN
            EQMAX = 2
            ISYM  = 32
            IF (.NOT. E12 .AND. .NOT. E13) THEN
               IREF = 1
            ELSE IF (.NOT. E12 .AND. .NOT. E23) THEN
               IREF = 2
            ELSE
               IREF = 3
            ENDIF
         ELSE
            EQMAX = 1
            ISYM  = 33
            IREF = 1
            IF (EA(2) .GT. EA(IREF)) IREF = 2
            IF (EA(3) .GT. EA(IREF)) IREF = 3
         ENDIF
      ELSE IF (ns .EQ. 2) THEN
         E12 = EAN_EQ(EA(1),EA(2))
         IF (E12) THEN
            EQMAX = 2
            ISYM  = 20
            IREF  = 1
         ELSE
            EQMAX = 1
            ISYM  = 21
            IF (EA(1) .GE. EA(2)) THEN
               IREF = 1
            ELSE
               IREF = 2
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
*Deck ASSIGN_OFFSETS
      SUBROUTINE ASSIGN_OFFSETS(ns,iref,OFF)
      INTEGER ns, iref
      DOUBLE PRECISION OFF(3)
      DOUBLE PRECISION PI, a120, a90
      PI   = DACOS(-1.0D0)
      a120 = 2.0D0*PI/3.0D0
      a90  = 0.5D0*PI
      IF (ns .EQ. 3) THEN
         OFF(1)=0.0D0
         OFF(2)= a120
         OFF(3)=-a120
         IF (iref .EQ. 2) THEN
            OFF(2)=0.0D0
            OFF(1)= a120
            OFF(3)=-a120
         ELSE IF (iref .EQ. 3) THEN
            OFF(3)=0.0D0
            OFF(1)= a120
            OFF(2)=-a120
         ENDIF
      ELSE IF (ns .EQ. 2) THEN
         OFF(1)=-a90
         OFF(2)= a90
         IF (iref .EQ. 1) THEN
            OFF(1)= a90
            OFF(2)=-a90
         ENDIF
      ELSE
         OFF(1)=0.0D0
      ENDIF
      RETURN
      END
C===============================================================
C  SCAN_ALL_BONDS_RIGID — stampa V0,V1,V2,V3,V6 per ogni legame
C===============================================================
      SUBROUTINE SCAN_ALL_BONDS_RIGID(IOut, IPrint,
     &     MxBnd, NAtoms, IAn, EAn, IRing, NLen, IAtom,
     &     NBond, IBond, C)

      INTEGER IOut, IPrint, MxBnd, NAtoms, NLen
      INTEGER IAn(NAtoms), IRing(NAtoms)
      INTEGER IAtom(2,NLen), NBond(NAtoms), IBond(MxBnd,NAtoms)
      DOUBLE PRECISION EAN(NAtoms),C(3,NAtoms)

      INTEGER L, IAt, JAt
      DOUBLE PRECISION V0, V1, V2, V3, V6
      INTEGER NSA, NSB, KREF, LREF, IFLAG
      CHARACTER*2 IAnEl2,SK,SI,SJ,SL,SYI,SYJ
      INTEGER NHI,NYI,NHJ,NYJ,IAnYI,IAnYJ

      DO 100 L = 1, NLen
         IAt = IAtom(1, L)
         JAt = IAtom(2, L)

         CALL TORSEN(IAt, JAt, MxBnd, NAtoms, IAn, EAn, IRing, NLen, 
     &        IAtom, NBond, IBond, C, V0, V1, V2, V3, V6, NSA, NSB, 
     &        KREF, LREF, IFLAG)

         IF (IFLAG .NE. 0) THEN
            IF (IPrint .NE. 0) THEN
               IF (IFLAG .EQ. 1) THEN
                  WRITE(IOut,'(A,I6,2I6)') 'SKIP terminal bond  L,I,J=',
     &                 L, IAt, JAt
               ELSE IF (IFLAG .EQ. 2) THEN
                  WRITE(IOut,'(A,I6,2I6)') 'SKIP both in rings  L,I,J=',
     &                 L, IAt, JAt
               ELSE
                 WRITE(IOut,'(A,I6,2I6)') 'SKIP no substituents L,I,J=',
     &                 L, IAt, JAt
               ENDIF
            ENDIF
            GOTO 100
         ENDIF
        
         SK=IAnEl2(IAn(KRef))
         SI=IAnEl2(IAn(IAt))
         SJ=IAnEl2(IAn(JAt))
         SL=IAnEl2(IAn(LRef))
         SYI='  '
         SYJ='  '
         Call SymDih(IOut,IPrint,MxBnd,IAt,JAt,IAN,NBond,IBond,EAN,C,
     $    NHI,NYI,NHJ,NYJ,IAnYI,IAnYJ)
         If(NYI.gt.1) SYI=IAnEl2(IAnYI)
         If(NYJ.gt.1) SYJ=IAnEl2(IAnYJ)
         If(NHI.gt.1) then
          Write(IOut,'(3X,2I3,5X,2I3,5F8.2,6X,''H'',I1,A2,''-'',A2)')
     $     IAt, JAt, NSA, NSB, V0, V1, V2, V3, V6, NHI, SI, SJ
         ElseIf(NYI.gt.1) then
          Write(IOut,'(3X,2I3,5X,2I3,5F8.2,5X,A2,I1,A2,''-'',A2)')
     $     IAt, JAt, NSA, NSB, V0, V1, V2, V3, V6, SYI, NYI, SI, SJ
         ElseIf(NHJ.gt.1) then
          Write(IOut,'(3X,2I3,5X,2I3,5F8.2,6X,''H'',I1,A2,''-'',A2)')
     $     JAt, IAt, NSB, NSA, V0, V1, V2, V3, V6, NHJ, SJ, SI
         ElseIf(NYJ.gt.1) then
          Write(IOut,'(3X,2I3,5X,2I3,5F8.2,5X,A2,I1,A2,''-'',A2)')
     $     JAt, IAt, NSB, NSA, V0, V1, V2, V3, V6, SYJ,NYJ, SJ, SI
         Else
          Write(IOut,'(4I3,2X,2I3,5F8.2,3X,4A3)')
     $    KRef,IAt,JAt,LRef,NSA,NSB,V0,V1,V2,V3,V6,SK,SI,SJ,SL
         EndIf

  100 CONTINUE

      RETURN
      END
*Deck SymDih
      Subroutine SymDih(IOut,IPrint,MxBnd,JAt,KAt,IAN,NBond,IBond,EAN,C,
     $  NHJ,NYJ,NHK,NYK,IAnYI,IAnYJ)
      Implicit Real*8 (A-H,O-Z)
      Dimension IAN(*),NBond(*),IBond(MxBnd,*)
      Dimension C(3,*),EAN(*)
      NHJ=0
      NYJ=0
      NHK=0
      NYK=0
      ERef=0.0D0
      Do 10 ii=1,NBond(JAt)
       IAt=IBond(ii,JAt)
       If(IAt.eq.KAt) goto 10
       If(ERef.eq.0.0D0) ERef=EAn(IAt)
       If(IAn(IAt).eq.1) NHJ=NHJ+1
       If(EAN(IAt).eq.ERef) then
        NYJ=NYJ+1
        IAnYJ=IAn(IAt)
       EndIf
   10 Continue 
      Do 20 ii=1,NBond(KAt)
       IAt=IBond(ii,KAt)
       If(IAt.eq.JAt) goto 20
       If(ERef.eq.0.0D0) ERef=EAn(IAt)
       If(IAn(IAt).eq.1) NHK=NHK+1
       If(EAN(IAt).eq.ERef) then
        NYK=NYK+1
        IAnYI=IAn(IAt)
       EndIf
   20 Continue 
      Return
      End
