*Deck AAbs
       Subroutine AAbs(N,X,Y)
      Implicit Real*8(A-H,O-Z)
C
C     Get the absolute value for elements of array X and store in Y
C
      Dimension X(N), Y(N)
C
      Do 10 I = 1, N
       Y(I) = Abs(X(I))
   10 Continue
      Return
      End
*Deck AAdd
      Subroutine AAdd(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to do vector operation
C     C = A + B
C
      Dimension A(*), B(*), C(*)
      Do 10 I = 1, N
       C(I) = A(I) + B(I)
   10 continue
      Return
      End
*Deck ACASB
      Subroutine ACASB(N,A,B,C,S)
      Implicit Real*8(A-H,O-Z)
C
C     Perform the vector opeation:
C     C = A + S * B
C     where A, B and C are vectors of length N and S is a scalar.
C
      Dimension A(*), B(*), C(*)
      Do 10 I = 1, N
       C(I) = A(I) + S*B(I)
   10 continue
      Return
      End
*Deck Aclear
      Subroutine AClear(N,A)
      Implicit Real*8 (A-H,O-Z)
C
C     Clear N elements of A
C
      Dimension A(*)
      Zero =0.0D0
      Do 10 I=1,N
       A(I)=Zero
   10 continue
      Return
      End
*Deck AHpSrA
      Subroutine AHpSrA(N,Ind,A)
      Implicit Integer(A-Z)
C
C     Sort Real*8 array A using HeapSort, returning the new order in Ind.
C     This version sorts based on absolute value.
C
      Real*8 A(*), aK
      Integer Ind(*)
C
      L = N/2 + 1
      R = N
      Do 10 I = 1, N
       Ind(I) = I
   10 continue
      If(N.le.1) Return
   20 If(L.gt.1) then
        L = L - 1
        IndCur = Ind(L)
        aK = Abs(A(IndCur))
      else
        IndCur = Ind(R)
        aK = Abs(A(IndCur))
        Ind(R) = Ind(1)
        R = R - 1
        endIf
      If(L.gt.1.or.R.gt.1) then
        I = L
        J = 2*L
   30   If(J.le.R) then
          If(J.lt.R) then
            If(Abs(A(Ind(J))).lt.Abs(A(Ind(J+1)))) J = J + 1
            endIf
          If(aK.lt.Abs(A(Ind(J)))) then
            Ind(I) = Ind(J)
            I = J
            J = 2*J
          else
            J = R + 1
            endIf
          Goto 30
          endIf
        Ind(I) = IndCur
        Goto 20
        endIf
      Ind(1) = IndCur
      Return
      End
*Deck AHpSrS
      Subroutine AHpSrS(N,Ind,A)
      Implicit Real*8(A-H,O-Z)
C
C     Sort Real*8 array A using HeapSort, returning the new order in Ind.
C     This routine fudges A to make the sort stable.
C
      Parameter (Zero=0.0d0,FactR=1.d-15)
      Dimension Ind(*), A(*)
C
      Big = Zero
      Do 10 I = 1, N
       Big = Max(Big,Abs(A(I)))
   10 continue
      Fact = Big*FactR
      Do 20 I = 1, N
       II = N - I + 1
       A(I) = A(I) + Fact*II
   20 continue
      Call AHpSrt(N,Ind,A)
      Do 30 I = 1, N
       II = N - I + 1
       A(I) = A(I) - Fact*II
   30 continue
      Return
      End
*Deck AHpSrt
      Subroutine AHpSrt(N,Ind,A)
      Implicit Integer(A-Z)
C
C     Sort Real*8 array A using HeapSort, returning the new order in Ind.
C
      Real*8 A(N), aK
      Integer Ind(N)
C
      L = N/2 + 1
      R = N
      Do 10 I = 1, N
       Ind(I) = I
   10 continue
      If(N.le.1) Return
   20 If(L.gt.1) then
        L = L - 1
        IndCur = Ind(L)
        aK = A(IndCur)
      else
        IndCur = Ind(R)
        aK = A(IndCur)
        Ind(R) = Ind(1)
        R = R - 1
        endIf
      If(L.gt.1.or.R.gt.1) then
        I = L
        J = 2*L
   30   If(J.le.R) then
          If(J.lt.R) then
            If(A(Ind(J)).lt.A(Ind(J+1))) J = J + 1
            endIf
          If(aK.lt.A(Ind(J))) then
            Ind(I) = Ind(J)
            I = J
            J = 2*J
          else
            J = R + 1
            endIf
          Goto 30
          endIf
        Ind(I) = IndCur
        Goto 20
        endIf
      Ind(1) = IndCur
      Return
      End
*Deck AMove
      Subroutine AMove(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Move N words from A to B.  Do not use this routine for
C     overlapping arrays!
C
      Common/IO/in,iout,ipunch
      Dimension A(*), B(*)
C
        If(.not.(N.le.0.or.Loc(B(N)).lt.Loc(A(1)).or.
     $  Loc(B(1)).gt.Loc(A(N)).or.
     $  (Loc(B(1)).eq.Loc(A(1)).and.Loc(B(N)).eq.Loc(A(N))))) then
         write(IOut,'('' Overlap in AMove'')')
         Stop
        endif
      Do 10 I = 1, N
       B(I) = A(I)
   10 continue
      Return
      End
*Deck AMove1
      Subroutine AMove1(N,IndF,IndT,A)
      Implicit Real*8(A-H,O-Z)
C
C     Move N elements from A(IndF+1) to A(IndT+1), doing the right
C     thing regardless of whether IndF is larger or smaller than IndF.
C
      Dimension A(*)
C
      If(IndF.gt.IndT) then
        Do 10 I = 1, N
         A(IndT+I) = A(IndF+I)
   10   continue
      else if(IndF.lt.IndT) then
        Do 20 I = N, 1, -1
         A(IndT+I) = A(IndF+I)
   20   continue
        endIf
      Return
      End
*Deck ANeg
      Subroutine ANeg(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Vector B = -A.
C
      Dimension A(*), B(*)
      Do 10 I = 1, N
       B(I) = -A(I)
   10  continue
      Return
      End
*Deck AngDeg
      Function AngDeg(IFlag,N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Return the angle between vectors A and B.
C     Iflag=1 Angle is Returned in Radians.
C     Iflag=0 Angle is Returned in Degrees.
C
      Real*8 AngDeg, MDCutO
      Dimension A(1), B(1)
      Data Zero/0.0d0/, One/1.0d0/
C
      R1R2Sq = SProd(N,A,A)*SProd(N,B,B)
      RMin = MDCutO(0)**2
      CosA = Zero
      If(R1R2Sq.ge.RMin) CosA = SProd(N,A,B) / Sqrt(R1R2Sq)
      If(Abs(CosA).gt.One) CosA = Sign(One,CosA)
      If(IFlag.eq.0) then
        AngDeg = ACos(CosA) * Float(45) / ATan(One)
      else
        AngDeg = ACos(CosA)
        endIf
      Return
      End
*Deck ArMax1
      Function ArMax1(A,N,IfAbs,IMax)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns the maximum element of an array.  IfAbs
C     determines whether the absolute values are to be compared.
C     IMax is set to the index of the largest element.
C
      Dimension A(*)
      Logical IfAbs
      Save Zero
      Data Zero/0.0d0/
C
      If(N.lt.1) then
        ArMax1 = Zero
        IMax = 0
      else
        IMax = 1
        If(IfAbs) then
          AM = Abs(A(1))
          Do 10 I = 2, N
            AV = Abs(A(I))
            If(AV.gt.AM) then
              IMax = I
              AM = AV
              endIf
   10       Continue
          ArMax1 = AM
        else
          Do 20 I = 2, N
            If(A(I).gt.A(IMax)) IMax = I
   20       Continue
          ArMax1 = A(IMax)
          endIf
        endIf
      Return
      End
*Deck ArMin1
      Function ArMin1(A,N,IfAbs,IMin)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns the minimum element of an array.  IfAbs
C     determines whether the absolute values are to be compared.
C     IMin is set to the index of the largest element.
C
      Dimension A(*)
      Logical IfAbs
      Save Zero
      Data Zero/0.0d0/
C
      If(N.lt.1) then
        ArMin1 = Zero
        IMin = 0
      else
        IMin = 1
        If(IfAbs) then
          AM = Abs(A(1))
          Do 10 I = 2, N
            AV = Abs(A(I))
            If(AV.lt.AM) then
              IMin = I
              AM = AV
              endIf
   10       Continue
          ArMin1 = AM
        else
          Do 20 I = 2, N
            If(A(I).lt.A(IMin)) IMin = I
   20       Continue
          ArMin1 = A(IMin)
          endIf
        endIf
      Return
      End
*Deck AScale
      Subroutine AScale(N,S,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Vector B = S * A
C
      Dimension A(*),B(*)
      Do 10 I = 1, N
       B(I) = S*A(I)
   10 continue
      Return
      End
*Deck ASet
      Subroutine ASet(N,Val,X)
      Implicit Real*8(A-H,O-Z)
C
C     Set all elements of array X to Val.
C
C$Acc Routine Gang
      Dimension X(*)
C
C$Acc Loop Gang Vector Worker
      Do 10 I = 1, N
   10   X(I) = Val
      Return
      End
*Deck ASub
      Subroutine ASub(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to perform the vector operation
C     C = A - B
C     where C, A, and B are vectors of length N.
C
      Dimension C(*), A(*), B(*)
        Do 10 I = 1, N
         C(I) = A(I) - B(I)
   10   continue
      Return
      End
*Deck ASUNIT
      SUBROUTINE ASUNIT(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     THIS SUBROUTINE LOADS C WITH A UNIT VECTOR IN THE A-B DIRECTION.
C
      DIMENSION A(1), B(1), C(1)
      CALL ASUB(N,A,B,C)
      CALL AUNIT(N,C,C)
      RETURN
      END
*Deck AtomP1
      Subroutine AtomP1(IOut,PrintH,IUseRL,MapRow,IUseCL,MapCol,ColLab,
     $  NAtoms,IAn,E,MDim,NCols)
      Implicit Real*8(A-H,O-Z)
C
C     Print a matrix of NCols values for each atom.  IUse[R/C]L is
C     0 for raw numbers, 1 to use map array, 2 to use label array.
C
      Logical PrintH
      Parameter (MaxEl=200,NPerL=6)
      Character*(*) ColLab(*), Lab(NPerL)*11
      Dimension MapRow(*), MapCol(*), E(MDim,*), IAn(*), IEl(0:MaxEl)
 1000 Format(10X,11A)
 1010 Format(5X,11I11)
 1020 Format(I6,2X,A2,11F11.6)
C
      Call FillEl(0,MaxEl,IEl)
      Do 30 JOff = 0, (NCols-1), NPerL
        NDo = Min(NCols-JOff,NPerL)
        JStart = JOff + 1
        JEnd = JOff + NDo
        If(IUseCL.eq.2) then
          Do 10 J = 1, NDo
            Lab(J) = ' '
            LL = LinEnd(ColLab(JOff+J))
            ISt = Max((11-LL)/2+1,1)
            Lab(J)(ISt:) = ColLab(JOff+J)
   10     continue
          Write(IOut,1000) (Lab(J),J=1,NDo)
        else if(IUseCL.eq.1) then
          Write(IOut,1010) (MapCol(J),J=JStart,JEnd)
        else
          Write(IOut,1010) (J,J=JStart,JEnd)
          endIf
        Do 20 I = 1, NAtoms
          IP = I
          If(IUseRL.eq.1) IP = MapRow(IP)
          If(PrintH.or.IAn(I).ne.1) Write(IOut,1020) IP, IEl(IAn(I)),
     $      (E(I,J),J=JStart,JEnd)
   20     Continue
   30   Continue
      Return
      End
*Deck AtomPr
      Subroutine AtomPr(IOut,PrintH,NAtoms,IAn,E,MDim,NCols)
      Implicit Real*8(A-H,O-Z)
C
C     Print a matrix of NCols values for each atom.
C
      Logical PrintH
      Character*1 CDum(1)
      Dimension E(MDim,*), IAn(*), JJ(1)
      Save CDum, JJ
      Data CDum/' '/, JJ/0/
C
      Call AtomP1(IOut,PrintH,0,JJ,0,JJ,CDum,NAtoms,IAn,E,MDim,NCols)
      Return
      End
*Deck AUnit
      Subroutine AUnit(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Load B with a unit vector in the direction of A.
C     A and B can be the same.
C
      Dimension A(*), B(*)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      R = Sqrt(SProd(N,A,A))
      If(R.eq.Zero) then
        Call AClear(N,B)
      else
        Call AScale(N,One/R,A,B)
        endIf
      Return
      End
*Deck AUnitM
      Subroutine AUnitM(LT,NRI,ND,N,A)
      Implicit Real*8(A-H,O-Z)
C
C     Set matrix A to the unit matrix.
C
      Logical LT
      Dimension A(NRI,*)
      Save One
      Data One/1.0d0/
C
      If(LT) then
        NTT = NRI*(N*(N+1))/2
        Call AClear(NTT,A)
        Do 10 I = 1, N
          II = (I*(I+1))/2
   10     A(1,II) = One
      else
        Call AClear(NRI*ND*N,A)
        Do 20 I = 1, N
          II = ND*(I-1) + I
   20     A(1,II) = One
        endIf
      Return
      End
*Deck ChkAnM
      Subroutine ChkAnM(IOut,Natoms,AtMass,CReac,CProd,OK)
      Implicit Real*8 (A-H,O-Z)
C
C Check Angular Momentum
C     Checks that no total angular momentum is present along the linear
C     path connecting reactants  with products.
C     0 = Sum{ Atmass(I) * [CReac(I)/CProd(I)]}
C
C Input:
C     AtMass : (NAtoms) Atomic masses
C     CReac  : (3,NAtoms) Cartesian Coordinates of Reactants
C     CProd  : (3,NAtoms) Cartesian Coordinates of Products
C
C Output:
C     OK     : True if all Angular Moments components are small enough
C              (lower than Tol). False otherwise
C
C     Dimensions
      Integer NAtoms
C     Input
      Real*8 CProd(3,*), CReac(3,*), AtMass(*)
C     Output
      Logical OK
C     Local
      Integer ia, ix
      Real*8 XCross(3), XSum(3), Tol
      Save Tol
      Data Tol/1.0d-8/
 9000 Format(/,1X,78('-'),/,10X,'WARNING: SIGNIFICANT ANGULAR MOMENTUM',
     $ ' ALONG THE PATH',/,
     $ 11X,'Ax =',D12.5,'  Ay=',D12.5,'  Az=',D12.5,/,1X,78('-'))
C
      OK = .True.
      Call AClear(3,XSum)
      Do 100 ia = 1, NAtoms
        Call Aclear(3,XCross)
C       Call VProd(XCross,CReac(1,ia),CProd(1,ia))
        Do 110 ix = 1, 3
         XSum(ix) = XSum(ix) + XCross(ix)*AtMass(ia)
  110   continue
  100   Continue
      Do 200 ix = 1, 3
        If(Abs(XSum(ix)).gt.Tol) OK = .False.
  200   Continue
      If(.not.OK.and.IOut.gt.0) Write(IOut,9000) (XSum(ix), ix=1,3)
      Return
      End
*Deck CntMas
      Subroutine CntMas(ITrans,NCopy,NAtoms,AtMass,C,TotWt,COM)
      Implicit Real*8(A-H,O-Z)
C
C     Compute the coordinates of the center of mass (charge, etc.).
C     If ITrans is +/- 1, C is translated to the COM.
C
      Dimension AtMass(*), C(3,NAtoms,*), COM(3)
      Save Zero
      Data Zero/0.0d0/
C
      COM(1) = Zero
      COM(2) = Zero
      COM(3) = Zero
      TotWt = Zero
      Do 20 IC = 1, NCopy
        Do 10 IAt = 1, NAtoms
          TotWt = TotWt + AtMass(IAt)
          COM(1) = COM(1) + AtMass(IAt)*C(1,IAt,IC)
          COM(2) = COM(2) + AtMass(IAt)*C(2,IAt,IC)
   10     COM(3) = COM(3) + AtMass(IAt)*C(3,IAt,IC)
   20   Continue
      If(TotWt.gt.Zero) then
        COM(1) = COM(1) / TotWt
        COM(2) = COM(2) / TotWt
        COM(3) = COM(3) / TotWt
        endIf
      If(ITrans.ne.0) then
        Do 30 IC = 1, NCopy
          Call Transl(ITrans,3,NAtoms,COM,C(1,1,IC))
   30     Continue
        endIf
      Return
      End
*Deck CnvFct
      Real*8 Function CnvFct(Factor)
      Implicit Real*8(A-H,O-Z)
C
C     Conversion Factors for vibrational analysis
C
C     Description
C     -----------
C     Computes the required factor and give it in output
C
C     Input
C     -----
C     Factor :: character(len=*)
C         Identifier of the requested factor
C
C     Returns
C     -------
C     real*8
C         Conversion factor
C         au2amu  -> a.u. (m_e) to a.m.u ("proton" mass)
C         au2Deb  -> a.u. (e^-^.Bohr) to Debye conversion factor
C         au2cm1  -> a.u. (Hartree) to cm^-1^
C         au2kJM  -> a.u. (Hartree) to kJ/mol
C         au2kCaM -> a.u. (Hartree) to kCal/mol
C         hbar    -> h/(2.pi)                   in amu.Ang^2^.s^-1^
C         Fac0AU  -> 1/(h.c)                    in cm^-1^.Eh^-1^
C         Fac1AU  -> 1/(2.pi.h^1/2^.c^3/2^)     in cm^-1^.(abc).Eh^-1^
C         Fac2AU  -> 1/(4.pi^2^.c^2^)           in cm^-1^.(abc)^2^.Eh^-1^
C         Fac3AU  -> 1/(8.pi^3^.h^-1/2^.c^5/2^) in cm^-1^.(abc)^3^.Eh^-1^
C         Fac4AU  -> 1/(16.pi^4^.h^-1^.c^3^)    in cm^-1^.(abc)^4^.Eh^-1^
C         Fact1   -> 1/(2.pi.h^1/2^.c^3/2^)     in aAc.cm^-1^.aJ^-1^
C         Fact2   -> 1/(4.pi^2^.c^2^)           in aAc^2^.cm^-1^.aJ^-1^
C         Fact3   -> 1/(8.pi^3^.h^-1/2^.c^5/2^) in aAc^3^.cm^-1^.aJ^-1^
C         Fact4   -> 1/(16.pi^4^.h^-1^.c^3^)    in aAc^4^.cm^-1^.aJ^-1^
C         FactA   -> h/(8.pi^2^.c^2^)           in aAc^2^
C         FactB   -> h.c/(2.kB)                 in cm.K
C         FactG   -> 4.pi^2^.c/h                in cm.amu^-1^.Ang^-2^
C         FactC   -> kB/c^2^                    in abc^2^.cm^-1^.K^-1^
C         HC      -> h.c                        in aJ.cm == mdyn.Ang.cm
C         MWQ2q   -> Q->q conversionin          Bohr.amu^1/2^.cm^-1/2^
C         PiCH12 -> pi.(c/h)^1/2^               in aAc^-1^
C         ToUMA  -> amu to uma
C
C     Notes
C     -----
C     * `FactC' is replaced to correct relation in [CorAvr]
C       Old version: kB/(4.pi^2^)             in amu.K^-1^
C     * Conventions
C       `aAc':: amu^1/2^.Ang.cm^-1/2^
C       `abc':: amu^1/2^.Bohr.cm^-1/2^
C       `aJ':: attoJoule
C       `kB':: Boltzmann constant
C       `Na':: Avogadro constant:q
C       `Eh':: Hartree (energy unit)
C
CEND
C
C     Parameters
      Common /IO/ In, IOut, IPunch
      Common /PhyCon/ PhyCon(30)
C     Input
      Real*8 PhyCon
      Character*(*) Factor
C     Local
      Integer JJ(1)
      Real*8 hbar, F1, F2, F4, F10P3, F10P10, F10P18, FactG, HC, m2Ang,
     $  Pi
      Character StrUp*8
      Save JJ
      Save F1, F2, F4, F10P3, F10P10, F10P18
      Data JJ/0/
      Data F1/1.0D0/, F2/2.0D0/, F4/4.0D0/, F10P3/1.0D3/,
     $  F10P10/1.0D10/, F10P18/1.0D18/
C
      Call LinupC(Factor,StrUp)
C
      Pi    = F4*ATan(F1)
C     m2Ang : meter-to-Angstroms conversion factor
      m2Ang = F10P10
C     hbar  : Reduced Planck constant in amu.Ang^2.s^-1
      hbar  = PhyCon(4)*m2Ang**2/(F2*Pi*PhyCon(2))
C     FactG : 4.pi^2.c/h in in cm.amu^-1.Ang^-2
      FactG = F2*Pi*PhyCon(9)/hbar
C     HC    : h.c in attoJ.cm == mdyn.Ang.cm.
      HC    = PhyCon(4)*PhyCon(9)*F10P18
C
      If(StrUp.eq.'FAC0AU') then
        CnvFct = PhyCon(8)*F10P18/HC
      else if(StrUp.eq.'FAC1AU'.or.StrUp.eq.'FACT1') then
       CnvFct = F1/(Sqrt(FactG)*HC)
       If(StrUp.eq.'FAC1AU')
     $  CnvFct = CnvFct*PhyCon(8)*F10P18/PhyCon(1)
      else if(StrUp.eq.'FAC2AU'.or.StrUp.eq.'FACT2') then
       CnvFct = F1/(FactG*HC)
       If(StrUp.eq.'FAC2AU')
     $  CnvFct = CnvFct*PhyCon(8)*F10P18/PhyCon(1)**2
      else if(StrUp.eq.'FAC3AU'.or.StrUp.eq.'FACT3') then
        CnvFct = F1/(FactG*Sqrt(FactG)*HC)
        If(StrUp.eq.'FAC3AU')
     $    CnvFct = CnvFct*PhyCon(8)*F10P18/PhyCon(1)**3
      else if(StrUp.eq.'FAC4AU'.or.StrUp.eq.'FACT4') then
        CnvFct = F1/(FactG**2*HC)
        If(StrUp.eq.'FAC4AU')
     $    CnvFct = CnvFct*PhyCon(8)*F10P18/PhyCon(1)**4
      else if(StrUp.eq.'FACTG') then
        CnvFCt = FactG
      else if(StrUp.eq.'FACTA') then
        CnvFct = F1/(F2*FactG)
      else if(StrUp.eq.'FACTB') then
        CnvFct = PhyCon(4)*PhyCon(9)/(F2*PhyCon(10))
      else if(StrUp.eq.'FACTC') then
        CnvFct = PhyCon(10)*m2Ang**2
     $    /(PhyCon(9)**2*PhyCon(2)*PhyCon(1)**2)
      else if(StrUp.eq.'HBAR') then
        CnvFct = hbar
      else if(StrUp.eq.'TOUMA') then
        CnvFct = PhyCon(2)/PhyCon(12)
      else if(StrUp.eq.'PICH12') then
        CnvFct = Pi*Sqrt(PhyCon(9)*PhyCon(2)/PhyCon(4)/m2Ang**2)
      else if(StrUp.eq.'AU2DEB') then
        CnvFct = F10P10*PhyCon(3)*PhyCon(1)
      else if(StrUp.eq.'MWQ2Q') then
        CnvFCt = F1/(Sqrt(FactG)*PhyCon(1))
      else if(StrUp.eq.'AU2CM1') then
        CnvFct = PhyCon(8)/(PhyCon(4)*PhyCon(9))
      else if(StrUp.eq.'AU2AMU') then
        CnvFCt = PhyCon(12)/PhyCon(2)
      else if(StrUp.eq.'AU2KJM') then
        CnvFct = PhyCon(8)*PhyCon(5)/F10P3
      else if(StrUp.eq.'AU2KCAM') then
        CnvFct = PhyCon(8)*PhyCon(5)/(F10P3*PhyCon(6))
      else if(StrUp.eq.'FACG0AU') then
        CnvFct = F1
      else if(StrUp.eq.'FACG1AU') then
        CnvFct = F1/(PhyCon(1)*Sqrt(FactG))
      else if(StrUp.eq.'FACG2AU') then
        CnvFct = (hbar/(PhyCon(1)**2))
     $    *(F1/(F2*Pi*PhyCon(9)))
        CnvFCt = F1/(PhyCon(1)**2*FactG)
      else
       write(IOut,'(''Unrecognized keyword for CnvFct'')')
       Stop
      endIf
      Return
      End
*Deck DStJac
      Function DStJac(JSt,IEnd,N,NDimA,ATop,A)
      Implicit Real*8(A-H,O-Z)
C
C     Compute stopping factor for Jacobi.
C
      Dimension A(NDimA,*)
      Save Zero
      Data Zero/0.0d0/
C
      D = Zero
      Do 10 J = JSt, N
       Do 20 I = 1, Min(J-1,IEnd)
        S = A(I,J) / ATop
        D = D + S*S
   20  continue
   10 continue
      DStJac = D
      Return
      End
*Deck Eckart
      Subroutine Eckart(IOut,IPrint,Move,NAtoms,C,AtMass,COM,Linear,
     $  PMom,RotMat)
      Implicit Real*8 (A-H,O-Z)
C 
C    Ir Representation, Iz < Ix < Iy
C   IIr Representation, Iy < Iz < Ix
C  IIIr Representation, Ix < Iy < Iz
C    Il Representation, Iz < Iy < Ix
C   IIl Representation, Ix < Iz < Iy
C  IIIl Representation, Iy < Ix < Iz
C
      Logical Prol,Move,Linear,Spher
      Dimension AtMass(*),C(3,*)
      Dimension PMom(3),RotMat(3,3)
      Dimension COM(3),CS(3)
      Dimension Ini(3)
      Save Small,One,Three
      Data Small/1.0D-8/, One/1.0D0/, Three/3.0d0/
      CS(1)=ArMin1(PMom,3,.false.,Ini(1))
      CS(3)=ArMax1(PMom,3,.false.,Ini(3))
      Ini(2)=6-Ini(1)-Ini(3) 
      CS(2)=PMom(Ini(2))
      Dif13=CS(3)-CS(1) 
      Dif21=CS(2)-CS(1)
      Prol=(Dif13.lt.Dif21) 
      Linear=(CS(1).lt.Small)
      Spher=(.not.Linear.and.(CS(3)-CS(1)).lt.Small)
      If(Linear) then
       Write(IOut,'(/,'' This Molecule is a Linear Top'',2x)')
      ElseIf(Spher) then
       Write(IOut,'(/,'' This Molecule is a Spherical Top'',2x)')
      ElseIf(Dif21.lt.small) then
       Write(IOut,'(/,'' This Molecule is an Oblate Symmetric Top'',
     $    2x)')
      ElseIf(Dif13.lt.small) then
       Write(IOut,'(/,'' This Molecule is a Prolate Symmetric Top'',
     $    2X)') 
      ElseIf(Prol) then
       Write(IOut,'(/,'' This Molecule is a Nearly Prolate Asymmetric'',
     $   '' Top'',2X)')
      Else
       Write(IOut,'(/,'' This Molecule is a Nearly Oblate Asymmetric'', 
     $   '' Top'',2X)')
      EndIf
      If(move) then
       write(IOut,'(5X,''Ir Representation'')')
       Do 10 IAt=1,NAtoms
        Do 10 IXYZ=1,3
         C(IXYZ,IAT) =  C(IXYZ,IAt) - COM(IXYZ)
   10  continue
       If(Spher) Return
C        Call RigHnd(RotMat)
         Call RotF1(NAtoms,RotMat,C)
         Call AMove(3,CS,PMom)
         Do 30 IAt=1,NAtoms
        call AMove(3,C(1,IAt),CS)
        do 30 IXYZ=1,3
         C(Ini(IXYZ),IAt)=CS(IXYZ) 
   30  continue
      EndIf
      Return
      End
*Deck IAnEl2
      Character*2 Function IAnEl2(IAn)
      Implicit Integer(A-Z)
C
C     Returns the atomic symbol for a given atomic number.
C
C Input:
C     IAn    : Atomic Number
C     NoCase : Ignore case
C Output:
C     IAnEl2 : Atomic Symbol (A2)
C     Local
      Integer Pos
      Character AllSmb*238
C
      AllSmb='BqH HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV '
     $ //'CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSn'
     $ //'SbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIr'
     $ //'PtAuHgTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRf'
     $ //'DbSgBhHsMtDsRgCnNhFlMcLvTsOg'
      If(IAn.lt.0.or.IAn.gt.118) then
       IAnEl2='  '
      Else
       Pos=2*IAn+1
       IAnEl2=AllSmb(Pos:Pos+1)
      EndIf
      Return
      End
*Deck El2IAN
      Integer Function El2IAN(NoCase,AtSymb)
      Implicit Integer(A-Z)
C
C     Returns the atomic number for a given atomic symbol.
C
C Input:
C     NoCase : Ignore case while looking for atomic number
C     AtSymb : Atomic symbol (character)
C
C     Input
      Character AtSymb*2
      Logical NoCase
C     Local
      Integer i, Pos
      Character ASymb*2, TSymb*2
      Character AllSmb*238
C
      AllSmb='BqH HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV '
     $ //'CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSn'
     $ //'SbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIr'
     $ //'PtAuHgTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRf'
     $ //'DbSgBhHsMtDsRgCnNhFlMcLvTsOg'
      If(NoCase) then
        Call LinUpC(AtSymb,ASymb)
      else
        ASymb = AtSymb
      endIf
      If(ASymb(1:1).eq.' '.and.ASymb(2:2).ne.' ') then
       ASymb(1:1)=ASymb(2:2)
       ASymb(2:2)=' '
      EndIf
      El2IAN = -1
      Do 10 i=0,118
       Pos=2*i+1
       TSymb=AllSmb(Pos:Pos+1)
       If(NoCase) Call LinUpC(TSymb,TSymb)
       If(ASymb.eq.TSymb) then
        El2IAN=i
        Return
       EndIf
   10 Continue
      Return
      End
*Deck ElNeg
      Subroutine ElNeg(IOut,IPrint,NAtoms,IAn,ELnA)
      IMPLICIT NONE
      INTEGER IOut,IPrint,NAtoms,N,I,IAt
      Integer IAn(*)
      PARAMETER (N = 99)
      CHARACTER*3 SYMBOL(N)
      Real*8 EN(N),ElnA(*)
C     L'array SYMBOL contiene il simbolo dell'elemento (3 caratteri);
C     l'array EN contiene il valore di elettronegatività (4° parametro nel codice C++)
      DATA SYMBOL /'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ',
     +     'F  ','NE ','NA ','MG ','AL ','SI ','P  ','S  ',
     +     'CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ',
     +     'MN ','FE ','CO ','NI ','CU ','ZN ','GA ','GE ',
     +     'AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ',
     +     'NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ',
     +     'IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ',
     +     'LA ','CE ','PR ','ND ','PM ','SM ','EU ','GD ',
     +     'TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ',
     +     'TA ','W  ','RE ','OS ','IR ','PT ','AU ','HG ',
     +     'TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ',
     +     'AC ','TH ','PA ','U  ','NP ','PU ','AM ','CM ',
     +     'BK ','CF ','ES '/
      DATA EN /2.1D0,4.16D0,0.97D0,1.47D0,2.01D0,2.50D0,3.07D0,3.50D0,
     + 4.10D0,4.78D0,1.01D0,1.23D0,1.47D0,1.74D0,2.06D0,2.44D0,
     + 2.83D0,3.24D0,0.91D0,1.04D0,1.20D0,1.32D0,1.45D0,1.56D0,
     + 1.60D0,1.64D0,1.70D0,1.75D0,1.75D0,1.66D0,1.82D0,2.02D0,
     + 2.20D0,2.48D0,2.74D0,2.97D0,0.89D0,0.99D0,1.11D0,1.22D0,
     + 1.23D0,1.30D0,1.36D0,1.42D0,1.45D0,1.35D0,1.42D0,1.46D0,
     + 1.49D0,1.72D0,1.82D0,2.01D0,2.21D0,2.58D0,0.86D0,0.97D0,
     + 1.08D0,1.08D0,1.07D0,1.07D0,1.07D0,1.07D0,1.01D0,1.11D0,
     + 1.10D0,1.10D0,1.10D0,1.11D0,1.11D0,1.06D0,1.14D0,1.23D0,
     + 1.33D0,1.40D0,1.46D0,1.52D0,1.55D0,1.44D0,1.42D0,1.44D0,
     + 1.44D0,1.55D0,1.67D0,1.76D0,1.90D0,2.60D0,0.86D0,0.97D0,
     + 1.00D0,1.11D0,1.14D0,1.22D0,1.22D0,1.22D0,1.20D0,1.20D0,
     + 1.20D0,1.20D0,1.20D0/
      If(IPrint.gt.1) then
       write(IOut,'('' Electronegativity values:'')')
       write(IOut,'('' Z   Symbol   Electronegativity'')')
       Do 10 I = 1, N
        write(IOut,'(I3,2X,A3,3X,F5.2)') I, SYMBOL(I), EN(I)
  10   continue 
      EndIf
      Do 20 IAt=1,NAtoms
       ElnA(IAt)=EN(IAn(IAt))
  20  continue
      return
      end
*Deck Elipse
      Subroutine Elipse(NAtoms,IAn,C,Radius,Axis)
      Implicit Real*8(A-H,O-Z)
C
C     Estimate radii for an elipsoid of approximately the molecular
C     shape, using the maximum distance between atoms (including atomic
C     radius) as criterion.
C
      Dimension IAn(NAtoms), C(3,NAtoms), Radius(3), Axis(3,3), VIJ(3)
      Real*8 MDCutO
      Save Zero, Pt5, One
      Data Zero/0.0d0/, Pt5/0.5d0/, One/1.0d0/
C
C     First find the maximum distance.
C
      Call AClear(3,Radius)
      Call AClear(9,Axis)
      Small = MDCutO(0)
      RadMax = Zero
      Do 10 I = 1, NAtoms
C       RadI = SBondL(5,IAn(I),IAn(I))
        RadI = RCovCT(IAn(I),IAn(I)) 
        If(RadI.gt.RadMax) RadMax = RadI
   10   Continue
      Do 100 IAxis = 1, 3
        Do 50 I = 2, NAtoms
C          RadI = Pt5*SBondL(5,IAn(I),IAn(I))
           RadI = RCovCT(IAn(I),0)
          Do 50 J = 1, NAtoms
C           RadJ = Pt5*SBondL(5,IAn(J),IAn(J))
            RadJ = RCovCT(IAn(J),0)
            Call ASub(3,C(1,J),C(1,I),VIJ)
            Do 30 JAxis = 1, (IAxis-1)
              ADot1 = SProd(3,VIJ,Axis(1,JAxis))
              Call ACasB(3,VIJ,Axis(1,JAxis),VIJ,-ADot1)
   30         Continue
            RIJ = Sqrt(SProd(3,VIJ,VIJ))
            RIJE = RIJ + RadI + RadJ
            If(RIJE.gt.Radius(IAxis).and.RIJ.gt.Small) then
              Call AScale(3,(One/RIJ),VIJ,Axis(1,IAxis))
              Radius(IAxis) = RIJE
              endIf
   50       Continue
        If(Radius(IAxis).eq.Zero) then
          Radius(IAxis) = RadMax
          Do 70 JAxis = 1, 3
            Call AClear(3,VIJ)
            VIJ(JAxis) = One
            Do 60 KAxis = 1, (IAxis-1)
              Dot = SProd(3,Axis(1,KAxis),VIJ)
              Call ACasB(3,VIJ,Axis(1,KAxis),VIJ,-Dot)
   60         Continue
            RNew = Sqrt(SProd(3,VIJ,VIJ))
            If(RNew.gt.Small) then
              Call AScale(3,One/RNew,VIJ,VIJ)
              Goto 80
              endIf
   70       Continue
   80     Call AMove(3,VIJ,Axis(1,IAxis))
          endIf
  100   Continue
      Call AScale(3,Pt5,Radius,Radius)
      Return
      End
*Deck EpsEta
      Subroutine EpsEta(Eps,Eta)
      Implicit Real*8(A-H,O-Z)
C
C     Compute and return Eta, the smallest representable number,
C     and Eps, the smallest number for which 1+Eps.ne.1.
C
      Save Zero, One, Two
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/
C
      Eta = One
   10 T = Eta / Two
      If(T.ne.Zero) then
        Eta = T
        Goto 10
        endIf
      Eps = One
   20 T = One + (Eps/Two)
      If(T.ne.One) then
        Eps = Eps / Two
        Goto 20
        endIf
      Return
      End
*Deck ExpPE
      Function ExpPE(eps,rm,alph,r)
      Implicit Real*8 (A-H,O-Z)
C
C     computes the van der waals modified Morse function (Exp-PE) by
C     Yang, Sun and Deng, J.Phys.Chem. A 124, 2102-2107 (2020)
C     distances in Angstrom and energy in kJ/mol
C
      Data one,two,three /1.0d0,2.0d0,3.0d0/
      Save one,two,three
      term   = r/rm
      facexp = one - term
      cnexp1 = exp(alph*facexp)
      cnexp2 = exp(alph*facexp/two)
      cnpol2 = term**2
      cnpol4 = cnpol2*cnpol2
      cnpol  = cnpol4 - two*cnpol2 + three
      ExpPE = eps*(cnexp1-cnpol*cnexp2)
      return
      end
*Deck FilAMS
      Subroutine FilAMS(NAtoms,IAn,AMS)
      Implicit Real*8(A-H,O-Z)
C
C     Load the mean atomic masses.
C
      Parameter (MaxAn=83)
      Dimension IAn(*), AMS(*), AMSD(MaxAn)
      Save AMSD, Zero
      Data AMSD/   1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,
     $10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,
     $20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,
     $30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,
     $40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,
     $54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,
     $65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,
     $79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,
     $91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,
     $102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,
     $118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,
     $132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,
     $183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,
     $196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0/
      Data Zero/0.0d0/
C
      Do 10 I = 1, NAtoms
        If(IAn(I).le.0) then
          AMS(I) = Zero
        else if(IAn(I).le.MaxAn) then
          AMS(I) = AMSD(IAn(I))
        else
          write(IOut,*) 'Atomic number out of range in FilAMS.'
          Stop
        endif
   10   Continue
      Return
      End
*Deck FillEl
      Subroutine FillEl(ISt,IEnd,El)
      Implicit Integer(A-Z)
C
C     Load array El with the names of the elements from ISt to IEnd.
C     ISt can be zero, in which case El starts with Banquo, or -1, in
C     which case El starts with X
C
      Parameter (MinEl=-2,MaxEl=118)
      Dimension ElDat(MinEl:MaxEl), El(ISt:IEnd)
      Save ElDat, Quest
      Data ElDat/2hTV,2hX ,2HBq,2HH ,2HHe,2HLi,2HBe,2HB ,2HC ,
     $2HN ,2HO ,2HF ,2HNe,2HNa,2HMg,2HAl,2HSi,2HP ,2HS ,2HCl,2HAr,2HK ,
     $2HCa,2HSc,2HTi,2HV ,2HCr,2HMn,2HFe,2HCo,2HNi,2HCu,2HZn,2HGa,2HGe,
     $2HAs,2HSe,2HBr,2HKr,2HRb,2HSr,2HY ,2HZr,2HNb,2HMo,2HTc,2HRu,2HRh,
     $2HPd,2HAg,2HCd,2HIn,2HSn,2HSb,2HTe,2HI ,2HXe,2HCs,2HBa,2HLa,2HCe,
     $2HPr,2HNd,2HPm,2HSm,2HEu,2HGd,2HTb,2HDy,2HHo,2HEr,2HTm,2HYb,2HLu,
     $2HHf,2HTa,2HW ,2HRe,2HOs,2HIr,2HPt,2HAu,2HHg,2HTl,2HPb,2HBi,2HPo,
     $2HAt,2HRn,2HFr,2HRa,2HAc,2HTh,2HPa,2HU ,2hNp,2hPu,2hAm,2hCm,2hBk,
     $2hCf,2hEs,2hFm,2hMd,2hNo,2hLr,2hRf,2hDb,2hSg,2hBh,2hHs,2hMt,2hDs,
     $2hRg,2hCn,2hNh,2hFl,2hMc,2hLv,2hTs,2hOg/,Quest/2h??/
C
      ISt1 = Max(ISt,MinEl)
      IEnd1 = Min(IEnd,MaxEl)
      Do 10 I = ISt, (ISt1-1)
       El(I) = Quest
   10 continue
      Do 20 I = ISt1, IEnd1
       El(I) = ElDat(I)
   20 continue
      Do 30 I = (IEnd1+1), IEnd
       El(I) = Quest
   30 continue 
      Return
      End
*Deck FilIAn
      Subroutine FilIAn(El,IAn)
      Integer El2IAn
      Character*10 Num
      Character*2 El
      Save Num
      Data Num/'0123456789'/
  100 Format(A1)
      IAn=El2IAN(.True.,El)
      If(IAn.gt.-1) Return
      If(El(2:2).eq.' ') then
       do 10 i=1,10
        If(El(1:1).eq.Num(i:i)) IAn=IAn+i
   10  continue
      Else
       Do 20 i=1,10
        If(El(2:2).eq.Num(i:i)) IAn=Ian+i
        If(El(1:1).eq.Num(i:i)) IAn=IAn+(i-1)*10
   20  Continue
      EndIf
      Return
      End
*Deck FilMag 
      Subroutine FilMag(MNI,IAnI,JUse,MNO,RMass,ISpin,QMom,GFac)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine looks up an integer mass number in the tables and
C     returns its index and the nuclear properties. Negative values of
C     MNI fetch that element in the table. JUse is set to 0 if MNI is
C     out of range and to -1 if IAnI is out of range.
C
C     Masses are stored in AM: AM(I,J) is the mass for the Jth
C     most abundant isotope for atomic number I.  The corresponding mass
C     numbers (# of protons + # of neutrons) are stored in the array MN.
C     The masses are from A. H. Wapstra and K. Bos, Atomic and
C     Nuclear Data Tables, 1977, 19, 185.  The isotopic abundances from
C     the 1981-82 CRC tables were used.  Masses of atoms after Kr were
C     also taken from CRC.  For those atoms without tabulated natural
C     abundances, the isotopes of longest lifetime were selected.
C     ISpin holds 2*nuclear spin, GM the nuclear magnetic moment,
C     and QM the nuclear quadrupole moment.
C
C VARIABLES:
C   MNI    ... # of protons + # of neutrons
C   IAnI   ... Integer atomic number
C   JUse   ... Jth most abundant isotope
C   MNO    ...
C   RMass  ... Mass
C   ISpin  ... 2* the nuclear spin
C   QMom   ... Nuclear quadropole moment
C   GFac   ... Nuclear magnetic moment
C
      Parameter (MaxIso=4,MaxAn=109)
      Dimension AM(MaxIso,0:MaxAn), MN(MaxIso,0:MaxAn),
     $  IS(MaxIso,0:MaxAn), QM(MaxIso,0:MaxAn), GM(MaxIso,0:MaxAn)
      Save Zero, AM, MN, IS, GM, QM, Small
      Data Zero/0.0d0/, Small/1.0d-6/
      Data (MN(I,0),I=1,MaxIso)/MaxIso*0/,
     $  (AM(I,0),I=1,MaxIso)/MaxIso*0.0d0/,
     $  (IS(I,0),I=1,MaxIso)/MaxIso*0/,
     $  (GM(I,0),I=1,MaxIso)/MaxIso*0.0d0/,
     $  (QM(I,0),I=1,MaxIso)/MaxIso*0.0d0/
      Data MN(1,  1)/  1/, AM(1,  1)/   1.007825037d0/, IS(1,  1)/ 1/,
     $  GM(1,  1)/   2.792846000d0/, QM(1,  1)/   0.000000000d0/
      Data MN(2,  1)/  2/, AM(2,  1)/   2.014101787d0/, IS(2,  1)/ 2/,
     $  GM(2,  1)/   0.857438000d0/, QM(2,  1)/   0.286000000d0/
      Data MN(3,  1)/  3/, AM(3,  1)/   3.016049286d0/, IS(3,  1)/ 1/,
     $  GM(3,  1)/   2.978960000d0/, QM(3,  1)/   0.000000000d0/
      Data MN(4,  1)/  0/, AM(4,  1)/   0.000000000d0/, IS(4,  1)/ 0/,
     $  GM(4,  1)/   0.000000000d0/, QM(4,  1)/   0.000000000d0/
      Data MN(1,  2)/  4/, AM(1,  2)/   4.002603250d0/, IS(1,  2)/ 0/,
     $  GM(1,  2)/   0.000000000d0/, QM(1,  2)/   0.000000000d0/
      Data MN(2,  2)/  3/, AM(2,  2)/   3.016029297d0/, IS(2,  2)/ 1/,
     $  GM(2,  2)/  -2.127620000d0/, QM(2,  2)/   0.000000000d0/
      Data MN(3,  2)/  0/, AM(3,  2)/   0.000000000d0/, IS(3,  2)/ 0/,
     $  GM(3,  2)/   0.000000000d0/, QM(3,  2)/   0.000000000d0/
      Data MN(4,  2)/  0/, AM(4,  2)/   0.000000000d0/, IS(4,  2)/ 0/,
     $  GM(4,  2)/   0.000000000d0/, QM(4,  2)/   0.000000000d0/
      Data MN(1,  3)/  7/, AM(1,  3)/   7.016004500d0/, IS(1,  3)/ 3/,
     $  GM(1,  3)/   3.256424000d0/, QM(1,  3)/  -4.010000000d0/
      Data MN(2,  3)/  6/, AM(2,  3)/   6.015123200d0/, IS(2,  3)/ 2/,
     $  GM(2,  3)/   0.822047000d0/, QM(2,  3)/  -0.080800000d0/
      Data MN(3,  3)/  0/, AM(3,  3)/   0.000000000d0/, IS(3,  3)/ 0/,
     $  GM(3,  3)/   0.000000000d0/, QM(3,  3)/   0.000000000d0/
      Data MN(4,  3)/  0/, AM(4,  3)/   0.000000000d0/, IS(4,  3)/ 0/,
     $  GM(4,  3)/   0.000000000d0/, QM(4,  3)/   0.000000000d0/
      Data MN(1,  4)/  9/, AM(1,  4)/   9.012182500d0/, IS(1,  4)/ 3/,
     $  GM(1,  4)/  -1.177900000d0/, QM(1,  4)/   5.288000000d0/
      Data MN(2,  4)/  0/, AM(2,  4)/   0.000000000d0/, IS(2,  4)/ 0/,
     $  GM(2,  4)/   0.000000000d0/, QM(2,  4)/   0.000000000d0/
      Data MN(3,  4)/  0/, AM(3,  4)/   0.000000000d0/, IS(3,  4)/ 0/,
     $  GM(3,  4)/   0.000000000d0/, QM(3,  4)/   0.000000000d0/
      Data MN(4,  4)/  0/, AM(4,  4)/   0.000000000d0/, IS(4,  4)/ 0/,
     $  GM(4,  4)/   0.000000000d0/, QM(4,  4)/   0.000000000d0/
      Data MN(1,  5)/ 11/, AM(1,  5)/  11.009305300d0/, IS(1,  5)/ 3/,
     $  GM(1,  5)/   2.688637000d0/, QM(1,  5)/   4.059000000d0/
      Data MN(2,  5)/ 10/, AM(2,  5)/  10.012938000d0/, IS(2,  5)/ 6/,
     $  GM(2,  5)/   1.800650000d0/, QM(2,  5)/   8.459000000d0/
      Data MN(3,  5)/  0/, AM(3,  5)/   0.000000000d0/, IS(3,  5)/ 0/,
     $  GM(3,  5)/   0.000000000d0/, QM(3,  5)/   0.000000000d0/
      Data MN(4,  5)/  0/, AM(4,  5)/   0.000000000d0/, IS(4,  5)/ 0/,
     $  GM(4,  5)/   0.000000000d0/, QM(4,  5)/   0.000000000d0/
      Data MN(1,  6)/ 12/, AM(1,  6)/  12.000000000d0/, IS(1,  6)/ 0/,
     $  GM(1,  6)/   0.000000000d0/, QM(1,  6)/   0.000000000d0/
      Data MN(2,  6)/ 13/, AM(2,  6)/  13.003354839d0/, IS(2,  6)/ 1/,
     $  GM(2,  6)/   0.702411000d0/, QM(2,  6)/   0.000000000d0/
      Data MN(3,  6)/ 14/, AM(3,  6)/  14.003241000d0/, IS(3,  6)/ 0/,
     $  GM(3,  6)/   0.000000000d0/, QM(3,  6)/   0.000000000d0/
      Data MN(4,  6)/  0/, AM(4,  6)/   0.000000000d0/, IS(4,  6)/ 0/,
     $  GM(4,  6)/   0.000000000d0/, QM(4,  6)/   0.000000000d0/
      Data MN(1,  7)/ 14/, AM(1,  7)/  14.003074008d0/, IS(1,  7)/ 2/,
     $  GM(1,  7)/   0.403761000d0/, QM(1,  7)/   2.044000000d0/
      Data MN(2,  7)/ 15/, AM(2,  7)/  15.000108978d0/, IS(2,  7)/ 1/,
     $  GM(2,  7)/  -0.283190000d0/, QM(2,  7)/   0.000000000d0/
      Data MN(3,  7)/  0/, AM(3,  7)/   0.000000000d0/, IS(3,  7)/ 0/,
     $  GM(3,  7)/   0.000000000d0/, QM(3,  7)/   0.000000000d0/
      Data MN(4,  7)/  0/, AM(4,  7)/   0.000000000d0/, IS(4,  7)/ 0/,
     $  GM(4,  7)/   0.000000000d0/, QM(4,  7)/   0.000000000d0/
      Data MN(1,  8)/ 16/, AM(1,  8)/  15.994914640d0/, IS(1,  8)/ 0/,
     $  GM(1,  8)/   0.000000000d0/, QM(1,  8)/   0.000000000d0/
      Data MN(2,  8)/ 18/, AM(2,  8)/  17.999159390d0/, IS(2,  8)/ 0/,
     $  GM(2,  8)/   0.000000000d0/, QM(2,  8)/   0.000000000d0/
      Data MN(3,  8)/ 17/, AM(3,  8)/  16.999130600d0/, IS(3,  8)/ 5/,
     $  GM(3,  8)/  -1.893800000d0/, QM(3,  8)/  -2.558000000d0/
      Data MN(4,  8)/  0/, AM(4,  8)/   0.000000000d0/, IS(4,  8)/ 0/,
     $  GM(4,  8)/   0.000000000d0/, QM(4,  8)/   0.000000000d0/
      Data MN(1,  9)/ 19/, AM(1,  9)/  18.998403250d0/, IS(1,  9)/ 1/,
     $  GM(1,  9)/   2.628867000d0/, QM(1,  9)/   0.000000000d0/
      Data MN(2,  9)/  0/, AM(2,  9)/   0.000000000d0/, IS(2,  9)/ 0/,
     $  GM(2,  9)/   0.000000000d0/, QM(2,  9)/   0.000000000d0/
      Data MN(3,  9)/  0/, AM(3,  9)/   0.000000000d0/, IS(3,  9)/ 0/,
     $  GM(3,  9)/   0.000000000d0/, QM(3,  9)/   0.000000000d0/
      Data MN(4,  9)/  0/, AM(4,  9)/   0.000000000d0/, IS(4,  9)/ 0/,
     $  GM(4,  9)/   0.000000000d0/, QM(4,  9)/   0.000000000d0/
      Data MN(1, 10)/ 20/, AM(1, 10)/  19.992439100d0/, IS(1, 10)/ 0/,
     $  GM(1, 10)/   0.000000000d0/, QM(1, 10)/   0.000000000d0/
      Data MN(2, 10)/ 22/, AM(2, 10)/  21.991383700d0/, IS(2, 10)/ 0/,
     $  GM(2, 10)/   0.000000000d0/, QM(2, 10)/   0.000000000d0/
      Data MN(3, 10)/ 21/, AM(3, 10)/  20.993845300d0/, IS(3, 10)/ 3/,
     $  GM(3, 10)/  -0.661800000d0/, QM(3, 10)/  10.155000000d0/
      Data MN(4, 10)/  0/, AM(4, 10)/   0.000000000d0/, IS(4, 10)/ 0/,
     $  GM(4, 10)/   0.000000000d0/, QM(4, 10)/   0.000000000d0/
      Data MN(1, 11)/ 23/, AM(1, 11)/  22.989769700d0/, IS(1, 11)/ 3/,
     $  GM(1, 11)/   2.217520000d0/, QM(1, 11)/  10.400000000d0/
      Data MN(2, 11)/  0/, AM(2, 11)/   0.000000000d0/, IS(2, 11)/ 0/,
     $  GM(2, 11)/   0.000000000d0/, QM(2, 11)/   0.000000000d0/
      Data MN(3, 11)/  0/, AM(3, 11)/   0.000000000d0/, IS(3, 11)/ 0/,
     $  GM(3, 11)/   0.000000000d0/, QM(3, 11)/   0.000000000d0/
      Data MN(4, 11)/  0/, AM(4, 11)/   0.000000000d0/, IS(4, 11)/ 0/,
     $  GM(4, 11)/   0.000000000d0/, QM(4, 11)/   0.000000000d0/
      Data MN(1, 12)/ 24/, AM(1, 12)/  23.985045000d0/, IS(1, 12)/ 0/,
     $  GM(1, 12)/   0.000000000d0/, QM(1, 12)/   0.000000000d0/
      Data MN(2, 12)/ 26/, AM(2, 12)/  25.982595400d0/, IS(2, 12)/ 0/,
     $  GM(2, 12)/   0.000000000d0/, QM(2, 12)/   0.000000000d0/
      Data MN(3, 12)/ 25/, AM(3, 12)/  24.985839200d0/, IS(3, 12)/ 5/,
     $  GM(3, 12)/  -0.855460000d0/, QM(3, 12)/  19.940000000d0/
      Data MN(4, 12)/  0/, AM(4, 12)/   0.000000000d0/, IS(4, 12)/ 0/,
     $  GM(4, 12)/   0.000000000d0/, QM(4, 12)/   0.000000000d0/
      Data MN(1, 13)/ 27/, AM(1, 13)/  26.981541300d0/, IS(1, 13)/ 5/,
     $  GM(1, 13)/   3.641504000d0/, QM(1, 13)/  14.660000000d0/
      Data MN(2, 13)/  0/, AM(2, 13)/   0.000000000d0/, IS(2, 13)/ 0/,
     $  GM(2, 13)/   0.000000000d0/, QM(2, 13)/   0.000000000d0/
      Data MN(3, 13)/  0/, AM(3, 13)/   0.000000000d0/, IS(3, 13)/ 0/,
     $  GM(3, 13)/   0.000000000d0/, QM(3, 13)/   0.000000000d0/
      Data MN(4, 13)/  0/, AM(4, 13)/   0.000000000d0/, IS(4, 13)/ 0/,
     $  GM(4, 13)/   0.000000000d0/, QM(4, 13)/   0.000000000d0/
      Data MN(1, 14)/ 28/, AM(1, 14)/  27.976928400d0/, IS(1, 14)/ 0/,
     $  GM(1, 14)/   0.000000000d0/, QM(1, 14)/   0.000000000d0/
      Data MN(2, 14)/ 29/, AM(2, 14)/  28.976496400d0/, IS(2, 14)/ 1/,
     $  GM(2, 14)/  -0.555290000d0/, QM(2, 14)/   0.000000000d0/
      Data MN(3, 14)/ 30/, AM(3, 14)/  29.973771700d0/, IS(3, 14)/ 0/,
     $  GM(3, 14)/   0.000000000d0/, QM(3, 14)/   0.000000000d0/
      Data MN(4, 14)/  0/, AM(4, 14)/   0.000000000d0/, IS(4, 14)/ 0/,
     $  GM(4, 14)/   0.000000000d0/, QM(4, 14)/   0.000000000d0/
      Data MN(1, 15)/ 31/, AM(1, 15)/  30.973763400d0/, IS(1, 15)/ 1/,
     $  GM(1, 15)/   1.131600000d0/, QM(1, 15)/   0.000000000d0/
      Data MN(2, 15)/  0/, AM(2, 15)/   0.000000000d0/, IS(2, 15)/ 0/,
     $  GM(2, 15)/   0.000000000d0/, QM(2, 15)/   0.000000000d0/
      Data MN(3, 15)/  0/, AM(3, 15)/   0.000000000d0/, IS(3, 15)/ 0/,
     $  GM(3, 15)/   0.000000000d0/, QM(3, 15)/   0.000000000d0/
      Data MN(4, 15)/  0/, AM(4, 15)/   0.000000000d0/, IS(4, 15)/ 0/,
     $  GM(4, 15)/   0.000000000d0/, QM(4, 15)/   0.000000000d0/
      Data MN(1, 16)/ 32/, AM(1, 16)/  31.972071800d0/, IS(1, 16)/ 0/,
     $  GM(1, 16)/   0.000000000d0/, QM(1, 16)/   0.000000000d0/
      Data MN(2, 16)/ 34/, AM(2, 16)/  33.967867740d0/, IS(2, 16)/ 0/,
     $  GM(2, 16)/   0.000000000d0/, QM(2, 16)/   0.000000000d0/
      Data MN(3, 16)/ 33/, AM(3, 16)/  32.971459100d0/, IS(3, 16)/ 3/,
     $  GM(3, 16)/   0.643821000d0/, QM(3, 16)/  -6.780000000d0/
      Data MN(4, 16)/ 36/, AM(4, 16)/  35.967079000d0/, IS(4, 16)/ 0/,
     $  GM(4, 16)/   0.000000000d0/, QM(4, 16)/   0.000000000d0/
      Data MN(1, 17)/ 35/, AM(1, 17)/  34.968852729d0/, IS(1, 17)/ 3/,
     $  GM(1, 17)/   0.821874000d0/, QM(1, 17)/  -8.165000000d0/
      Data MN(2, 17)/ 37/, AM(2, 17)/  36.965902624d0/, IS(2, 17)/ 3/,
     $  GM(2, 17)/   0.684123000d0/, QM(2, 17)/  -6.435000000d0/
      Data MN(3, 17)/  0/, AM(3, 17)/   0.000000000d0/, IS(3, 17)/ 0/,
     $  GM(3, 17)/   0.000000000d0/, QM(3, 17)/   0.000000000d0/
      Data MN(4, 17)/  0/, AM(4, 17)/   0.000000000d0/, IS(4, 17)/ 0/,
     $  GM(4, 17)/   0.000000000d0/, QM(4, 17)/   0.000000000d0/
      Data MN(1, 18)/ 40/, AM(1, 18)/  39.962383100d0/, IS(1, 18)/ 0/,
     $  GM(1, 18)/   0.000000000d0/, QM(1, 18)/   0.000000000d0/
      Data MN(2, 18)/ 36/, AM(2, 18)/  35.967545605d0/, IS(2, 18)/ 0/,
     $  GM(2, 18)/   0.000000000d0/, QM(2, 18)/   0.000000000d0/
      Data MN(3, 18)/ 38/, AM(3, 18)/  37.962732200d0/, IS(3, 18)/ 0/,
     $  GM(3, 18)/   0.000000000d0/, QM(3, 18)/   0.000000000d0/
      Data MN(4, 18)/  0/, AM(4, 18)/   0.000000000d0/, IS(4, 18)/ 0/,
     $  GM(4, 18)/   0.000000000d0/, QM(4, 18)/   0.000000000d0/
      Data MN(1, 19)/ 39/, AM(1, 19)/  38.963707900d0/, IS(1, 19)/ 3/,
     $  GM(1, 19)/   0.391466000d0/, QM(1, 19)/   5.850000000d0/
      Data MN(2, 19)/ 41/, AM(2, 19)/  40.961825400d0/, IS(2, 19)/ 3/,
     $  GM(2, 19)/   0.214870000d0/, QM(2, 19)/   7.110000000d0/
      Data MN(3, 19)/ 40/, AM(3, 19)/  39.963998800d0/, IS(3, 19)/ 8/,
     $  GM(3, 19)/  -1.298100000d0/, QM(3, 19)/  -7.300000000d0/
      Data MN(4, 19)/  0/, AM(4, 19)/   0.000000000d0/, IS(4, 19)/ 0/,
     $  GM(4, 19)/   0.000000000d0/, QM(4, 19)/   0.000000000d0/
      Data MN(1, 20)/ 40/, AM(1, 20)/  39.962590700d0/, IS(1, 20)/ 0/,
     $  GM(1, 20)/   0.000000000d0/, QM(1, 20)/   0.000000000d0/
      Data MN(2, 20)/ 44/, AM(2, 20)/  43.955484800d0/, IS(2, 20)/ 0/,
     $  GM(2, 20)/   0.000000000d0/, QM(2, 20)/   0.000000000d0/
      Data MN(3, 20)/ 42/, AM(3, 20)/  41.958621800d0/, IS(3, 20)/ 0/,
     $  GM(3, 20)/   0.000000000d0/, QM(3, 20)/   0.000000000d0/
      Data MN(4, 20)/ 48/, AM(4, 20)/  47.952532000d0/, IS(4, 20)/ 0/,
     $  GM(4, 20)/   0.000000000d0/, QM(4, 20)/   0.000000000d0/
      Data MN(1, 21)/ 45/, AM(1, 21)/  44.955913600d0/, IS(1, 21)/ 7/,
     $  GM(1, 21)/   4.756483000d0/, QM(1, 21)/ -22.000000000d0/
      Data MN(2, 21)/  0/, AM(2, 21)/   0.000000000d0/, IS(2, 21)/ 0/,
     $  GM(2, 21)/   0.000000000d0/, QM(2, 21)/   0.000000000d0/
      Data MN(3, 21)/  0/, AM(3, 21)/   0.000000000d0/, IS(3, 21)/ 0/,
     $  GM(3, 21)/   0.000000000d0/, QM(3, 21)/   0.000000000d0/
      Data MN(4, 21)/  0/, AM(4, 21)/   0.000000000d0/, IS(4, 21)/ 0/,
     $  GM(4, 21)/   0.000000000d0/, QM(4, 21)/   0.000000000d0/
      Data MN(1, 22)/ 48/, AM(1, 22)/  47.947946700d0/, IS(1, 22)/ 0/,
     $  GM(1, 22)/   0.000000000d0/, QM(1, 22)/   0.000000000d0/
      Data MN(2, 22)/ 46/, AM(2, 22)/  45.952632700d0/, IS(2, 22)/ 0/,
     $  GM(2, 22)/   0.000000000d0/, QM(2, 22)/   0.000000000d0/
      Data MN(3, 22)/ 47/, AM(3, 22)/  46.951764900d0/, IS(3, 22)/ 5/,
     $  GM(3, 22)/  -0.788480000d0/, QM(3, 22)/  30.200000000d0/
      Data MN(4, 22)/ 49/, AM(4, 22)/  48.947870500d0/, IS(4, 22)/ 7/,
     $  GM(4, 22)/  -1.104170000d0/, QM(4, 22)/  24.700000000d0/
      Data MN(1, 23)/ 51/, AM(1, 23)/  50.943962500d0/, IS(1, 23)/ 7/,
     $  GM(1, 23)/   5.151400000d0/, QM(1, 23)/  -5.200000000d0/
      Data MN(2, 23)/ 50/, AM(2, 23)/  49.947161300d0/, IS(2, 23)/12/,
     $  GM(2, 23)/   3.347450000d0/, QM(2, 23)/  21.000000000d0/
      Data MN(3, 23)/  0/, AM(3, 23)/   0.000000000d0/, IS(3, 23)/ 0/,
     $  GM(3, 23)/   0.000000000d0/, QM(3, 23)/   0.000000000d0/
      Data MN(4, 23)/  0/, AM(4, 23)/   0.000000000d0/, IS(4, 23)/ 0/,
     $  GM(4, 23)/   0.000000000d0/, QM(4, 23)/   0.000000000d0/
      Data MN(1, 24)/ 52/, AM(1, 24)/  51.940509700d0/, IS(1, 24)/ 0/,
     $  GM(1, 24)/   0.000000000d0/, QM(1, 24)/   0.000000000d0/
      Data MN(2, 24)/ 53/, AM(2, 24)/  52.940651000d0/, IS(2, 24)/ 3/,
     $  GM(2, 24)/  -0.474540000d0/, QM(2, 24)/ -15.000000000d0/
      Data MN(3, 24)/ 50/, AM(3, 24)/  49.946046300d0/, IS(3, 24)/ 0/,
     $  GM(3, 24)/   0.000000000d0/, QM(3, 24)/   0.000000000d0/
      Data MN(4, 24)/ 54/, AM(4, 24)/  53.938882200d0/, IS(4, 24)/ 0/,
     $  GM(4, 24)/   0.000000000d0/, QM(4, 24)/   0.000000000d0/
      Data MN(1, 25)/ 55/, AM(1, 25)/  54.938046300d0/, IS(1, 25)/ 5/,
     $  GM(1, 25)/   3.453200000d0/, QM(1, 25)/  33.000000000d0/
      Data MN(2, 25)/  0/, AM(2, 25)/   0.000000000d0/, IS(2, 25)/ 0/,
     $  GM(2, 25)/   0.000000000d0/, QM(2, 25)/   0.000000000d0/
      Data MN(3, 25)/  0/, AM(3, 25)/   0.000000000d0/, IS(3, 25)/ 0/,
     $  GM(3, 25)/   0.000000000d0/, QM(3, 25)/   0.000000000d0/
      Data MN(4, 25)/  0/, AM(4, 25)/   0.000000000d0/, IS(4, 25)/ 0/,
     $  GM(4, 25)/   0.000000000d0/, QM(4, 25)/   0.000000000d0/
      Data MN(1, 26)/ 56/, AM(1, 26)/  55.934939300d0/, IS(1, 26)/ 0/,
     $  GM(1, 26)/   0.000000000d0/, QM(1, 26)/   0.000000000d0/
      Data MN(2, 26)/ 54/, AM(2, 26)/  53.939612100d0/, IS(2, 26)/ 0/,
     $  GM(2, 26)/   0.000000000d0/, QM(2, 26)/   0.000000000d0/
      Data MN(3, 26)/ 57/, AM(3, 26)/  56.935395700d0/, IS(3, 26)/ 1/,
     $  GM(3, 26)/   0.090623000d0/, QM(3, 26)/  16.000000000d0/
      Data MN(4, 26)/ 58/, AM(4, 26)/  57.933277800d0/, IS(4, 26)/ 0/,
     $  GM(4, 26)/   0.000000000d0/, QM(4, 26)/   0.000000000d0/
      Data MN(1, 27)/ 59/, AM(1, 27)/  58.933197800d0/, IS(1, 27)/ 7/,
     $  GM(1, 27)/   4.627000000d0/, QM(1, 27)/  42.000000000d0/
      Data MN(2, 27)/  0/, AM(2, 27)/   0.000000000d0/, IS(2, 27)/ 0/,
     $  GM(2, 27)/   0.000000000d0/, QM(2, 27)/   0.000000000d0/
      Data MN(3, 27)/  0/, AM(3, 27)/   0.000000000d0/, IS(3, 27)/ 0/,
     $  GM(3, 27)/   0.000000000d0/, QM(3, 27)/   0.000000000d0/
      Data MN(4, 27)/  0/, AM(4, 27)/   0.000000000d0/, IS(4, 27)/ 0/,
     $  GM(4, 27)/   0.000000000d0/, QM(4, 27)/   0.000000000d0/
      Data MN(1, 28)/ 58/, AM(1, 28)/  57.935347100d0/, IS(1, 28)/ 0/,
     $  GM(1, 28)/   0.000000000d0/, QM(1, 28)/   0.000000000d0/
      Data MN(2, 28)/ 60/, AM(2, 28)/  59.930789000d0/, IS(2, 28)/ 0/,
     $  GM(2, 28)/   0.000000000d0/, QM(2, 28)/   0.000000000d0/
      Data MN(3, 28)/ 62/, AM(3, 28)/  61.928346400d0/, IS(3, 28)/ 0/,
     $  GM(3, 28)/   0.000000000d0/, QM(3, 28)/   0.000000000d0/
      Data MN(4, 28)/ 61/, AM(4, 28)/  60.931058600d0/, IS(4, 28)/ 3/,
     $  GM(4, 28)/  -0.750020000d0/, QM(4, 28)/  16.200000000d0/
      Data MN(1, 29)/ 63/, AM(1, 29)/  62.929599200d0/, IS(1, 29)/ 3/,
     $  GM(1, 29)/   2.223300000d0/, QM(1, 29)/ -22.000000000d0/
      Data MN(2, 29)/ 65/, AM(2, 29)/  64.927792400d0/, IS(2, 29)/ 3/,
     $  GM(2, 29)/   2.381700000d0/, QM(2, 29)/ -20.400000000d0/
      Data MN(3, 29)/  0/, AM(3, 29)/   0.000000000d0/, IS(3, 29)/ 0/,
     $  GM(3, 29)/   0.000000000d0/, QM(3, 29)/   0.000000000d0/
      Data MN(4, 29)/  0/, AM(4, 29)/   0.000000000d0/, IS(4, 29)/ 0/,
     $  GM(4, 29)/   0.000000000d0/, QM(4, 29)/   0.000000000d0/
      Data MN(1, 30)/ 64/, AM(1, 30)/  63.929145400d0/, IS(1, 30)/ 0/,
     $  GM(1, 30)/   0.000000000d0/, QM(1, 30)/   0.000000000d0/
      Data MN(2, 30)/ 66/, AM(2, 30)/  65.926035200d0/, IS(2, 30)/ 0/,
     $  GM(2, 30)/   0.000000000d0/, QM(2, 30)/   0.000000000d0/
      Data MN(3, 30)/ 68/, AM(3, 30)/  67.924845800d0/, IS(3, 30)/ 0/,
     $  GM(3, 30)/   0.000000000d0/, QM(3, 30)/   0.000000000d0/
      Data MN(4, 30)/ 67/, AM(4, 30)/  66.927128900d0/, IS(4, 30)/ 5/,
     $  GM(4, 30)/   0.875479000d0/, QM(4, 30)/  15.000000000d0/
      Data MN(1, 31)/ 69/, AM(1, 31)/  68.925580900d0/, IS(1, 31)/ 3/,
     $  GM(1, 31)/   2.016590000d0/, QM(1, 31)/  17.100000000d0/
      Data MN(2, 31)/ 71/, AM(2, 31)/  70.924700600d0/, IS(2, 31)/ 3/,
     $  GM(2, 31)/   2.562270000d0/, QM(2, 31)/  10.700000000d0/
      Data MN(3, 31)/  0/, AM(3, 31)/   0.000000000d0/, IS(3, 31)/ 0/,
     $  GM(3, 31)/   0.000000000d0/, QM(3, 31)/   0.000000000d0/
      Data MN(4, 31)/  0/, AM(4, 31)/   0.000000000d0/, IS(4, 31)/ 0/,
     $  GM(4, 31)/   0.000000000d0/, QM(4, 31)/   0.000000000d0/
      Data MN(1, 32)/ 74/, AM(1, 32)/  73.921178800d0/, IS(1, 32)/ 0/,
     $  GM(1, 32)/   0.000000000d0/, QM(1, 32)/   0.000000000d0/
      Data MN(2, 32)/ 72/, AM(2, 32)/  71.922080000d0/, IS(2, 32)/ 0/,
     $  GM(2, 32)/   0.000000000d0/, QM(2, 32)/   0.000000000d0/
      Data MN(3, 32)/ 70/, AM(3, 32)/  69.924249800d0/, IS(3, 32)/ 0/,
     $  GM(3, 32)/   0.000000000d0/, QM(3, 32)/   0.000000000d0/
      Data MN(4, 32)/ 73/, AM(4, 32)/  72.923463900d0/, IS(4, 32)/ 9/,
     $  GM(4, 32)/  -0.879470000d0/, QM(4, 32)/ -19.600000000d0/
      Data MN(1, 33)/ 75/, AM(1, 33)/  74.921595500d0/, IS(1, 33)/ 3/,
     $  GM(1, 33)/   1.439470000d0/, QM(1, 33)/  31.400000000d0/
      Data MN(2, 33)/  0/, AM(2, 33)/   0.000000000d0/, IS(2, 33)/ 0/,
     $  GM(2, 33)/   0.000000000d0/, QM(2, 33)/   0.000000000d0/
      Data MN(3, 33)/  0/, AM(3, 33)/   0.000000000d0/, IS(3, 33)/ 0/,
     $  GM(3, 33)/   0.000000000d0/, QM(3, 33)/   0.000000000d0/
      Data MN(4, 33)/  0/, AM(4, 33)/   0.000000000d0/, IS(4, 33)/ 0/,
     $  GM(4, 33)/   0.000000000d0/, QM(4, 33)/   0.000000000d0/
      Data MN(1, 34)/ 80/, AM(1, 34)/  79.916520500d0/, IS(1, 34)/ 0/,
     $  GM(1, 34)/   0.000000000d0/, QM(1, 34)/   0.000000000d0/
      Data MN(2, 34)/ 78/, AM(2, 34)/  77.917304000d0/, IS(2, 34)/ 0/,
     $  GM(2, 34)/   0.000000000d0/, QM(2, 34)/   0.000000000d0/
      Data MN(3, 34)/ 82/, AM(3, 34)/  81.916709000d0/, IS(3, 34)/ 0/,
     $  GM(3, 34)/   0.000000000d0/, QM(3, 34)/   0.000000000d0/
      Data MN(4, 34)/ 76/, AM(4, 34)/  75.919206600d0/, IS(4, 34)/ 0/,
     $  GM(4, 34)/   0.000000000d0/, QM(4, 34)/   0.000000000d0/
      Data MN(1, 35)/ 79/, AM(1, 35)/  78.918336100d0/, IS(1, 35)/ 3/,
     $  GM(1, 35)/   2.106399000d0/, QM(1, 35)/  31.300000000d0/
      Data MN(2, 35)/ 81/, AM(2, 35)/  80.916290000d0/, IS(2, 35)/ 3/,
     $  GM(2, 35)/   2.270560000d0/, QM(2, 35)/  26.200000000d0/
      Data MN(3, 35)/  0/, AM(3, 35)/   0.000000000d0/, IS(3, 35)/ 0/,
     $  GM(3, 35)/   0.000000000d0/, QM(3, 35)/   0.000000000d0/
      Data MN(4, 35)/  0/, AM(4, 35)/   0.000000000d0/, IS(4, 35)/ 0/,
     $  GM(4, 35)/   0.000000000d0/, QM(4, 35)/   0.000000000d0/
      Data MN(1, 36)/ 84/, AM(1, 36)/  83.911506400d0/, IS(1, 36)/ 0/,
     $  GM(1, 36)/   0.000000000d0/, QM(1, 36)/   0.000000000d0/
      Data MN(2, 36)/ 86/, AM(2, 36)/  85.910614000d0/, IS(2, 36)/ 0/,
     $  GM(2, 36)/   0.000000000d0/, QM(2, 36)/   0.000000000d0/
      Data MN(3, 36)/ 82/, AM(3, 36)/  81.913483000d0/, IS(3, 36)/ 0/,
     $  GM(3, 36)/   0.000000000d0/, QM(3, 36)/   0.000000000d0/
      Data MN(4, 36)/ 83/, AM(4, 36)/  82.914134000d0/, IS(4, 36)/ 9/,
     $  GM(4, 36)/  -0.970670000d0/, QM(4, 36)/  25.900000000d0/
      Data MN(1, 37)/ 85/, AM(1, 37)/  84.911700000d0/, IS(1, 37)/ 5/,
     $  GM(1, 37)/   1.353030000d0/, QM(1, 37)/  27.600000000d0/
      Data MN(2, 37)/ 87/, AM(2, 37)/  86.909187000d0/, IS(2, 37)/ 3/,
     $  GM(2, 37)/   2.751240000d0/, QM(2, 37)/  13.350000000d0/
      Data MN(3, 37)/  0/, AM(3, 37)/   0.000000000d0/, IS(3, 37)/ 0/,
     $  GM(3, 37)/   0.000000000d0/, QM(3, 37)/   0.000000000d0/
      Data MN(4, 37)/  0/, AM(4, 37)/   0.000000000d0/, IS(4, 37)/ 0/,
     $  GM(4, 37)/   0.000000000d0/, QM(4, 37)/   0.000000000d0/
      Data MN(1, 38)/ 88/, AM(1, 38)/  87.905600000d0/, IS(1, 38)/ 0/,
     $  GM(1, 38)/   0.000000000d0/, QM(1, 38)/   0.000000000d0/
      Data MN(2, 38)/ 84/, AM(2, 38)/  83.913400000d0/, IS(2, 38)/ 0/,
     $  GM(2, 38)/   0.000000000d0/, QM(2, 38)/   0.000000000d0/
      Data MN(3, 38)/ 86/, AM(3, 38)/  85.909400000d0/, IS(3, 38)/ 0/,
     $  GM(3, 38)/   0.000000000d0/, QM(3, 38)/   0.000000000d0/
      Data MN(4, 38)/ 87/, AM(4, 38)/  86.908900000d0/, IS(4, 38)/ 9/,
     $  GM(4, 38)/  -1.092830000d0/, QM(4, 38)/  33.500000000d0/
      Data MN(1, 39)/ 89/, AM(1, 39)/  88.905400000d0/, IS(1, 39)/ 1/,
     $  GM(1, 39)/  -0.137420000d0/, QM(1, 39)/   0.000000000d0/
      Data MN(2, 39)/  0/, AM(2, 39)/   0.000000000d0/, IS(2, 39)/ 0/,
     $  GM(2, 39)/   0.000000000d0/, QM(2, 39)/   0.000000000d0/
      Data MN(3, 39)/  0/, AM(3, 39)/   0.000000000d0/, IS(3, 39)/ 0/,
     $  GM(3, 39)/   0.000000000d0/, QM(3, 39)/   0.000000000d0/
      Data MN(4, 39)/  0/, AM(4, 39)/   0.000000000d0/, IS(4, 39)/ 0/,
     $  GM(4, 39)/   0.000000000d0/, QM(4, 39)/   0.000000000d0/
      Data MN(1, 40)/ 90/, AM(1, 40)/  89.904300000d0/, IS(1, 40)/ 0/,
     $  GM(1, 40)/   0.000000000d0/, QM(1, 40)/   0.000000000d0/
      Data MN(2, 40)/ 91/, AM(2, 40)/  90.905300000d0/, IS(2, 40)/ 5/,
     $  GM(2, 40)/  -1.303620000d0/, QM(2, 40)/ -17.600000000d0/
      Data MN(3, 40)/ 92/, AM(3, 40)/  91.904600000d0/, IS(3, 40)/ 0/,
     $  GM(3, 40)/   0.000000000d0/, QM(3, 40)/   0.000000000d0/
      Data MN(4, 40)/ 94/, AM(4, 40)/  93.906100000d0/, IS(4, 40)/ 0/,
     $  GM(4, 40)/   0.000000000d0/, QM(4, 40)/   0.000000000d0/
      Data MN(1, 41)/ 93/, AM(1, 41)/  92.906000000d0/, IS(1, 41)/ 9/,
     $  GM(1, 41)/   6.170500000d0/, QM(1, 41)/ -32.000000000d0/
      Data MN(2, 41)/  0/, AM(2, 41)/   0.000000000d0/, IS(2, 41)/ 0/,
     $  GM(2, 41)/   0.000000000d0/, QM(2, 41)/   0.000000000d0/
      Data MN(3, 41)/  0/, AM(3, 41)/   0.000000000d0/, IS(3, 41)/ 0/,
     $  GM(3, 41)/   0.000000000d0/, QM(3, 41)/   0.000000000d0/
      Data MN(4, 41)/  0/, AM(4, 41)/   0.000000000d0/, IS(4, 41)/ 0/,
     $  GM(4, 41)/   0.000000000d0/, QM(4, 41)/   0.000000000d0/
      Data MN(1, 42)/ 98/, AM(1, 42)/  97.905500000d0/, IS(1, 42)/ 0/,
     $  GM(1, 42)/   0.000000000d0/, QM(1, 42)/   0.000000000d0/
      Data MN(2, 42)/ 92/, AM(2, 42)/  91.906300000d0/, IS(2, 42)/ 0/,
     $  GM(2, 42)/   0.000000000d0/, QM(2, 42)/   0.000000000d0/
      Data MN(3, 42)/ 95/, AM(3, 42)/  94.905840000d0/, IS(3, 42)/ 5/,
     $  GM(3, 42)/  -0.914200000d0/, QM(3, 42)/  -2.200000000d0/
      Data MN(4, 42)/ 96/, AM(4, 42)/  95.904600000d0/, IS(4, 42)/ 0/,
     $  GM(4, 42)/   0.000000000d0/, QM(4, 42)/   0.000000000d0/
      Data MN(1, 43)/ 99/, AM(1, 43)/  98.906300000d0/, IS(1, 43)/ 0/,
     $  GM(1, 43)/   0.000000000d0/, QM(1, 43)/   0.000000000d0/
      Data MN(2, 43)/  0/, AM(2, 43)/   0.000000000d0/, IS(2, 43)/ 0/,
     $  GM(2, 43)/   0.000000000d0/, QM(2, 43)/   0.000000000d0/
      Data MN(3, 43)/  0/, AM(3, 43)/   0.000000000d0/, IS(3, 43)/ 0/,
     $  GM(3, 43)/   0.000000000d0/, QM(3, 43)/   0.000000000d0/
      Data MN(4, 43)/  0/, AM(4, 43)/   0.000000000d0/, IS(4, 43)/ 0/,
     $  GM(4, 43)/   0.000000000d0/, QM(4, 43)/   0.000000000d0/
      Data MN(1, 44)/102/, AM(1, 44)/ 101.903700000d0/, IS(1, 44)/ 0/,
     $  GM(1, 44)/   0.000000000d0/, QM(1, 44)/   0.000000000d0/
      Data MN(2, 44)/ 99/, AM(2, 44)/  98.906100000d0/, IS(2, 44)/ 5/,
     $  GM(2, 44)/  -0.641300000d0/, QM(2, 44)/   7.900000000d0/
      Data MN(3, 44)/100/, AM(3, 44)/  99.903000000d0/, IS(3, 44)/ 0/,
     $  GM(3, 44)/   0.000000000d0/, QM(3, 44)/   0.000000000d0/
      Data MN(4, 44)/104/, AM(4, 44)/ 103.905500000d0/, IS(4, 44)/ 0/,
     $  GM(4, 44)/   0.000000000d0/, QM(4, 44)/   0.000000000d0/
      Data MN(1, 45)/103/, AM(1, 45)/ 102.904800000d0/, IS(1, 45)/ 1/,
     $  GM(1, 45)/  -0.088400000d0/, QM(1, 45)/   0.000000000d0/
      Data MN(2, 45)/  0/, AM(2, 45)/   0.000000000d0/, IS(2, 45)/ 0/,
     $  GM(2, 45)/   0.000000000d0/, QM(2, 45)/   0.000000000d0/
      Data MN(3, 45)/  0/, AM(3, 45)/   0.000000000d0/, IS(3, 45)/ 0/,
     $  GM(3, 45)/   0.000000000d0/, QM(3, 45)/   0.000000000d0/
      Data MN(4, 45)/  0/, AM(4, 45)/   0.000000000d0/, IS(4, 45)/ 0/,
     $  GM(4, 45)/   0.000000000d0/, QM(4, 45)/   0.000000000d0/
      Data MN(1, 46)/106/, AM(1, 46)/ 105.903200000d0/, IS(1, 46)/ 0/,
     $  GM(1, 46)/   0.000000000d0/, QM(1, 46)/   0.000000000d0/
      Data MN(2, 46)/104/, AM(2, 46)/ 103.903600000d0/, IS(2, 46)/ 0/,
     $  GM(2, 46)/   0.000000000d0/, QM(2, 46)/   0.000000000d0/
      Data MN(3, 46)/105/, AM(3, 46)/ 104.904600000d0/, IS(3, 46)/ 5/,
     $  GM(3, 46)/  -0.642000000d0/, QM(3, 46)/  66.000000000d0/
      Data MN(4, 46)/108/, AM(4, 46)/ 107.903890000d0/, IS(4, 46)/ 0/,
     $  GM(4, 46)/   0.000000000d0/, QM(4, 46)/   0.000000000d0/
      Data MN(1, 47)/107/, AM(1, 47)/ 106.905090000d0/, IS(1, 47)/ 1/,
     $  GM(1, 47)/  -0.113570000d0/, QM(1, 47)/   0.000000000d0/
      Data MN(2, 47)/109/, AM(2, 47)/ 108.904700000d0/, IS(2, 47)/ 1/,
     $  GM(2, 47)/  -0.130690000d0/, QM(2, 47)/   0.000000000d0/
      Data MN(3, 47)/  0/, AM(3, 47)/   0.000000000d0/, IS(3, 47)/ 0/,
     $  GM(3, 47)/   0.000000000d0/, QM(3, 47)/   0.000000000d0/
      Data MN(4, 47)/  0/, AM(4, 47)/   0.000000000d0/, IS(4, 47)/ 0/,
     $  GM(4, 47)/   0.000000000d0/, QM(4, 47)/   0.000000000d0/
      Data MN(1, 48)/114/, AM(1, 48)/ 113.903600000d0/, IS(1, 48)/ 0/,
     $  GM(1, 48)/   0.000000000d0/, QM(1, 48)/   0.000000000d0/
      Data MN(2, 48)/110/, AM(2, 48)/ 109.903000000d0/, IS(2, 48)/ 0/,
     $  GM(2, 48)/   0.000000000d0/, QM(2, 48)/   0.000000000d0/
      Data MN(3, 48)/111/, AM(3, 48)/ 110.904200000d0/, IS(3, 48)/ 1/,
     $  GM(3, 48)/  -0.594890000d0/, QM(3, 48)/   0.000000000d0/
      Data MN(4, 48)/112/, AM(4, 48)/ 111.902800000d0/, IS(4, 48)/ 0/,
     $  GM(4, 48)/   0.000000000d0/, QM(4, 48)/   0.000000000d0/
      Data MN(1, 49)/115/, AM(1, 49)/ 114.904100000d0/, IS(1, 49)/ 9/,
     $  GM(1, 49)/   5.540800000d0/, QM(1, 49)/  81.000000000d0/
      Data MN(2, 49)/113/, AM(2, 49)/ 112.904300000d0/, IS(2, 49)/ 9/,
     $  GM(2, 49)/   5.528900000d0/, QM(2, 49)/  79.900000000d0/
      Data MN(3, 49)/  0/, AM(3, 49)/   0.000000000d0/, IS(3, 49)/ 0/,
     $  GM(3, 49)/   0.000000000d0/, QM(3, 49)/   0.000000000d0/
      Data MN(4, 49)/  0/, AM(4, 49)/   0.000000000d0/, IS(4, 49)/ 0/,
     $  GM(4, 49)/   0.000000000d0/, QM(4, 49)/   0.000000000d0/
      Data MN(1, 50)/118/, AM(1, 50)/ 117.901800000d0/, IS(1, 50)/ 0/,
     $  GM(1, 50)/   0.000000000d0/, QM(1, 50)/   0.000000000d0/
      Data MN(2, 50)/116/, AM(2, 50)/ 115.902100000d0/, IS(2, 50)/ 0/,
     $  GM(2, 50)/   0.000000000d0/, QM(2, 50)/   0.000000000d0/
      Data MN(3, 50)/117/, AM(3, 50)/ 116.903100000d0/, IS(3, 50)/ 1/,
     $  GM(3, 50)/  -1.001050000d0/, QM(3, 50)/   0.000000000d0/
      Data MN(4, 50)/119/, AM(4, 50)/ 118.903400000d0/, IS(4, 50)/ 1/,
     $  GM(4, 50)/  -1.047290000d0/, QM(4, 50)/ -12.800000000d0/
      Data MN(1, 51)/121/, AM(1, 51)/ 120.903800000d0/, IS(1, 51)/ 5/,
     $  GM(1, 51)/   3.363400000d0/, QM(1, 51)/ -36.000000000d0/
      Data MN(2, 51)/123/, AM(2, 51)/ 122.904100000d0/, IS(2, 51)/ 7/,
     $  GM(2, 51)/   2.549800000d0/, QM(2, 51)/ -49.000000000d0/
      Data MN(3, 51)/  0/, AM(3, 51)/   0.000000000d0/, IS(3, 51)/ 0/,
     $  GM(3, 51)/   0.000000000d0/, QM(3, 51)/   0.000000000d0/
      Data MN(4, 51)/  0/, AM(4, 51)/   0.000000000d0/, IS(4, 51)/ 0/,
     $  GM(4, 51)/   0.000000000d0/, QM(4, 51)/   0.000000000d0/
      Data MN(1, 52)/130/, AM(1, 52)/ 129.906700000d0/, IS(1, 52)/ 0/,
     $  GM(1, 52)/   0.000000000d0/, QM(1, 52)/   0.000000000d0/
      Data MN(2, 52)/125/, AM(2, 52)/ 124.904400000d0/, IS(2, 52)/ 1/,
     $  GM(2, 52)/  -0.888280000d0/, QM(2, 52)/   0.000000000d0/
      Data MN(3, 52)/126/, AM(3, 52)/ 125.903200000d0/, IS(3, 52)/ 0/,
     $  GM(3, 52)/   0.000000000d0/, QM(3, 52)/   0.000000000d0/
      Data MN(4, 52)/128/, AM(4, 52)/ 127.904700000d0/, IS(4, 52)/ 0/,
     $  GM(4, 52)/   0.000000000d0/, QM(4, 52)/   0.000000000d0/
      Data MN(1, 53)/127/, AM(1, 53)/ 126.900400000d0/, IS(1, 53)/ 5/,
     $  GM(1, 53)/   2.813280000d0/, QM(1, 53)/ -71.000000000d0/
      Data MN(2, 53)/  0/, AM(2, 53)/   0.000000000d0/, IS(2, 53)/ 0/,
     $  GM(2, 53)/   0.000000000d0/, QM(2, 53)/   0.000000000d0/
      Data MN(3, 53)/  0/, AM(3, 53)/   0.000000000d0/, IS(3, 53)/ 0/,
     $  GM(3, 53)/   0.000000000d0/, QM(3, 53)/   0.000000000d0/
      Data MN(4, 53)/  0/, AM(4, 53)/   0.000000000d0/, IS(4, 53)/ 0/,
     $  GM(4, 53)/   0.000000000d0/, QM(4, 53)/   0.000000000d0/
      Data MN(1, 54)/132/, AM(1, 54)/ 131.904200000d0/, IS(1, 54)/ 0/,
     $  GM(1, 54)/   0.000000000d0/, QM(1, 54)/   0.000000000d0/
      Data MN(2, 54)/129/, AM(2, 54)/ 128.904800000d0/, IS(2, 54)/ 1/,
     $  GM(2, 54)/  -0.777980000d0/, QM(2, 54)/   0.000000000d0/
      Data MN(3, 54)/131/, AM(3, 54)/ 130.905100000d0/, IS(3, 54)/ 3/,
     $  GM(3, 54)/   0.691861000d0/, QM(3, 54)/ -11.400000000d0/
      Data MN(4, 54)/134/, AM(4, 54)/ 133.905400000d0/, IS(4, 54)/ 0/,
     $  GM(4, 54)/   0.000000000d0/, QM(4, 54)/   0.000000000d0/
      Data MN(1, 55)/133/, AM(1, 55)/ 132.905429000d0/, IS(1, 55)/ 0/,
     $  GM(1, 55)/   0.000000000d0/, QM(1, 55)/   0.000000000d0/
      Data MN(2, 55)/  0/, AM(2, 55)/   0.000000000d0/, IS(2, 55)/ 0/,
     $  GM(2, 55)/   0.000000000d0/, QM(2, 55)/   0.000000000d0/
      Data MN(3, 55)/  0/, AM(3, 55)/   0.000000000d0/, IS(3, 55)/ 0/,
     $  GM(3, 55)/   0.000000000d0/, QM(3, 55)/   0.000000000d0/
      Data MN(4, 55)/  0/, AM(4, 55)/   0.000000000d0/, IS(4, 55)/ 0/,
     $  GM(4, 55)/   0.000000000d0/, QM(4, 55)/   0.000000000d0/
      Data MN(1, 56)/138/, AM(1, 56)/ 137.905000000d0/, IS(1, 56)/ 0/,
     $  GM(1, 56)/   0.000000000d0/, QM(1, 56)/   0.000000000d0/
      Data MN(2, 56)/134/, AM(2, 56)/ 133.904300000d0/, IS(2, 56)/ 0/,
     $  GM(2, 56)/   0.000000000d0/, QM(2, 56)/   0.000000000d0/
      Data MN(3, 56)/135/, AM(3, 56)/ 134.905600000d0/, IS(3, 56)/ 3/,
     $  GM(3, 56)/   0.837943000d0/, QM(3, 56)/  16.000000000d0/
      Data MN(4, 56)/136/, AM(4, 56)/ 135.904400000d0/, IS(4, 56)/ 0/,
     $  GM(4, 56)/   0.000000000d0/, QM(4, 56)/   0.000000000d0/
      Data MN(1, 57)/139/, AM(1, 57)/ 138.906100000d0/, IS(1, 57)/ 7/,
     $  GM(1, 57)/   2.783200000d0/, QM(1, 57)/  20.000000000d0/
      Data MN(2, 57)/138/, AM(2, 57)/ 137.906800000d0/, IS(2, 57)/10/,
     $  GM(2, 57)/   3.713900000d0/, QM(2, 57)/  45.000000000d0/
      Data MN(3, 57)/  0/, AM(3, 57)/   0.000000000d0/, IS(3, 57)/ 0/,
     $  GM(3, 57)/   0.000000000d0/, QM(3, 57)/   0.000000000d0/
      Data MN(4, 57)/  0/, AM(4, 57)/   0.000000000d0/, IS(4, 57)/ 0/,
     $  GM(4, 57)/   0.000000000d0/, QM(4, 57)/   0.000000000d0/
      Data MN(1, 58)/140/, AM(1, 58)/ 139.905300000d0/, IS(1, 58)/ 0/,
     $  GM(1, 58)/   0.000000000d0/, QM(1, 58)/   0.000000000d0/
      Data MN(2, 58)/138/, AM(2, 58)/ 137.905700000d0/, IS(2, 58)/ 0/,
     $  GM(2, 58)/   0.000000000d0/, QM(2, 58)/   0.000000000d0/
      Data MN(3, 58)/142/, AM(3, 58)/ 141.909000000d0/, IS(3, 58)/ 0/,
     $  GM(3, 58)/   0.000000000d0/, QM(3, 58)/   0.000000000d0/
      Data MN(4, 58)/  0/, AM(4, 58)/   0.000000000d0/, IS(4, 58)/ 0/,
     $  GM(4, 58)/   0.000000000d0/, QM(4, 58)/   0.000000000d0/
      Data MN(1, 59)/141/, AM(1, 59)/ 140.907400000d0/, IS(1, 59)/ 5/,
     $  GM(1, 59)/   4.136000000d0/, QM(1, 59)/  -5.890000000d0/
      Data MN(2, 59)/  0/, AM(2, 59)/   0.000000000d0/, IS(2, 59)/ 0/,
     $  GM(2, 59)/   0.000000000d0/, QM(2, 59)/   0.000000000d0/
      Data MN(3, 59)/  0/, AM(3, 59)/   0.000000000d0/, IS(3, 59)/ 0/,
     $  GM(3, 59)/   0.000000000d0/, QM(3, 59)/   0.000000000d0/
      Data MN(4, 59)/  0/, AM(4, 59)/   0.000000000d0/, IS(4, 59)/ 0/,
     $  GM(4, 59)/   0.000000000d0/, QM(4, 59)/   0.000000000d0/
      Data MN(1, 60)/142/, AM(1, 60)/ 141.907500000d0/, IS(1, 60)/ 0/,
     $  GM(1, 60)/   0.000000000d0/, QM(1, 60)/   0.000000000d0/
      Data MN(2, 60)/143/, AM(2, 60)/ 142.909600000d0/, IS(2, 60)/ 7/,
     $  GM(2, 60)/  -1.065000000d0/, QM(2, 60)/ -63.000000000d0/
      Data MN(3, 60)/144/, AM(3, 60)/ 143.909900000d0/, IS(3, 60)/ 0/,
     $  GM(3, 60)/   0.000000000d0/, QM(3, 60)/   0.000000000d0/
      Data MN(4, 60)/146/, AM(4, 60)/ 145.912700000d0/, IS(4, 60)/ 0/,
     $  GM(4, 60)/   0.000000000d0/, QM(4, 60)/   0.000000000d0/
      Data MN(1, 61)/145/, AM(1, 61)/ 144.912700000d0/, IS(1, 61)/ 0/,
     $  GM(1, 61)/   0.000000000d0/, QM(1, 61)/   0.000000000d0/
      Data MN(2, 61)/147/, AM(2, 61)/ 146.915100000d0/, IS(2, 61)/ 0/,
     $  GM(2, 61)/   0.000000000d0/, QM(2, 61)/   0.000000000d0/
      Data MN(3, 61)/  0/, AM(3, 61)/   0.000000000d0/, IS(3, 61)/ 0/,
     $  GM(3, 61)/   0.000000000d0/, QM(3, 61)/   0.000000000d0/
      Data MN(4, 61)/  0/, AM(4, 61)/   0.000000000d0/, IS(4, 61)/ 0/,
     $  GM(4, 61)/   0.000000000d0/, QM(4, 61)/   0.000000000d0/
      Data MN(1, 62)/152/, AM(1, 62)/ 151.919500000d0/, IS(1, 62)/ 0/,
     $  GM(1, 62)/   0.000000000d0/, QM(1, 62)/   0.000000000d0/
      Data MN(2, 62)/147/, AM(2, 62)/ 146.914600000d0/, IS(2, 62)/ 7/,
     $  GM(2, 62)/  -0.814900000d0/, QM(2, 62)/ -25.900000000d0/
      Data MN(3, 62)/149/, AM(3, 62)/ 148.916900000d0/, IS(3, 62)/ 7/,
     $  GM(3, 62)/  -0.671800000d0/, QM(3, 62)/   7.400000000d0/
      Data MN(4, 62)/154/, AM(4, 62)/ 153.922000000d0/, IS(4, 62)/ 0/,
     $  GM(4, 62)/   0.000000000d0/, QM(4, 62)/   0.000000000d0/
      Data MN(1, 63)/153/, AM(1, 63)/ 152.920900000d0/, IS(1, 63)/ 5/,
     $  GM(1, 63)/   1.533100000d0/, QM(1, 63)/ 241.200000000d0/
      Data MN(2, 63)/151/, AM(2, 63)/ 150.919600000d0/, IS(2, 63)/ 5/,
     $  GM(2, 63)/   3.471800000d0/, QM(2, 63)/  90.300000000d0/
      Data MN(3, 63)/  0/, AM(3, 63)/   0.000000000d0/, IS(3, 63)/ 0/,
     $  GM(3, 63)/   0.000000000d0/, QM(3, 63)/   0.000000000d0/
      Data MN(4, 63)/  0/, AM(4, 63)/   0.000000000d0/, IS(4, 63)/ 0/,
     $  GM(4, 63)/   0.000000000d0/, QM(4, 63)/   0.000000000d0/
      Data MN(1, 64)/158/, AM(1, 64)/ 157.924100000d0/, IS(1, 64)/ 0/,
     $  GM(1, 64)/   0.000000000d0/, QM(1, 64)/   0.000000000d0/
      Data MN(2, 64)/156/, AM(2, 64)/ 155.922100000d0/, IS(2, 64)/ 0/,
     $  GM(2, 64)/   0.000000000d0/, QM(2, 64)/   0.000000000d0/
      Data MN(3, 64)/157/, AM(3, 64)/ 156.933900000d0/, IS(3, 64)/ 3/,
     $  GM(3, 64)/  -0.339900000d0/, QM(3, 64)/ 135.000000000d0/
      Data MN(4, 64)/160/, AM(4, 64)/ 159.927100000d0/, IS(4, 64)/ 0/,
     $  GM(4, 64)/   0.000000000d0/, QM(4, 64)/   0.000000000d0/
      Data MN(1, 65)/159/, AM(1, 65)/ 158.925000000d0/, IS(1, 65)/ 3/,
     $  GM(1, 65)/   2.014000000d0/, QM(1, 65)/ 143.200000000d0/
      Data MN(2, 65)/151/, AM(2, 65)/ 150.923000000d0/, IS(2, 65)/ 0/,
     $  GM(2, 65)/   0.000000000d0/, QM(2, 65)/   0.000000000d0/
      Data MN(3, 65)/  0/, AM(3, 65)/   0.000000000d0/, IS(3, 65)/ 0/,
     $  GM(3, 65)/   0.000000000d0/, QM(3, 65)/   0.000000000d0/
      Data MN(4, 65)/  0/, AM(4, 65)/   0.000000000d0/, IS(4, 65)/ 0/,
     $  GM(4, 65)/   0.000000000d0/, QM(4, 65)/   0.000000000d0/
      Data MN(1, 66)/164/, AM(1, 66)/ 163.928800000d0/, IS(1, 66)/ 0/,
     $  GM(1, 66)/   0.000000000d0/, QM(1, 66)/   0.000000000d0/
      Data MN(2, 66)/161/, AM(2, 66)/ 160.926600000d0/, IS(2, 66)/ 5/,
     $  GM(2, 66)/  -0.480600000d0/, QM(2, 66)/ 250.700000000d0/
      Data MN(3, 66)/162/, AM(3, 66)/ 161.926500000d0/, IS(3, 66)/ 0/,
     $  GM(3, 66)/   0.000000000d0/, QM(3, 66)/   0.000000000d0/
      Data MN(4, 66)/163/, AM(4, 66)/ 162.928400000d0/, IS(4, 66)/ 5/,
     $  GM(4, 66)/   0.672600000d0/, QM(4, 66)/ 264.800000000d0/
      Data MN(1, 67)/165/, AM(1, 67)/ 164.930300000d0/, IS(1, 67)/ 7/,
     $  GM(1, 67)/   4.173000000d0/, QM(1, 67)/ 358.000000000d0/
      Data MN(2, 67)/  0/, AM(2, 67)/   0.000000000d0/, IS(2, 67)/ 0/,
     $  GM(2, 67)/   0.000000000d0/, QM(2, 67)/   0.000000000d0/
      Data MN(3, 67)/  0/, AM(3, 67)/   0.000000000d0/, IS(3, 67)/ 0/,
     $  GM(3, 67)/   0.000000000d0/, QM(3, 67)/   0.000000000d0/
      Data MN(4, 67)/  0/, AM(4, 67)/   0.000000000d0/, IS(4, 67)/ 0/,
     $  GM(4, 67)/   0.000000000d0/, QM(4, 67)/   0.000000000d0/
      Data MN(1, 68)/166/, AM(1, 68)/ 165.930400000d0/, IS(1, 68)/ 0/,
     $  GM(1, 68)/   0.000000000d0/, QM(1, 68)/   0.000000000d0/
      Data MN(2, 68)/167/, AM(2, 68)/ 166.932000000d0/, IS(2, 68)/ 7/,
     $  GM(2, 68)/  -0.566500000d0/, QM(2, 68)/ 356.500000000d0/
      Data MN(3, 68)/168/, AM(3, 68)/ 167.932400000d0/, IS(3, 68)/ 0/,
     $  GM(3, 68)/   0.000000000d0/, QM(3, 68)/   0.000000000d0/
      Data MN(4, 68)/170/, AM(4, 68)/ 169.935500000d0/, IS(4, 68)/ 0/,
     $  GM(4, 68)/   0.000000000d0/, QM(4, 68)/   0.000000000d0/
      Data MN(1, 69)/169/, AM(1, 69)/ 168.934400000d0/, IS(1, 69)/ 1/,
     $  GM(1, 69)/  -0.231600000d0/, QM(1, 69)/   0.000000000d0/
      Data MN(2, 69)/  0/, AM(2, 69)/   0.000000000d0/, IS(2, 69)/ 0/,
     $  GM(2, 69)/   0.000000000d0/, QM(2, 69)/   0.000000000d0/
      Data MN(3, 69)/  0/, AM(3, 69)/   0.000000000d0/, IS(3, 69)/ 0/,
     $  GM(3, 69)/   0.000000000d0/, QM(3, 69)/   0.000000000d0/
      Data MN(4, 69)/  0/, AM(4, 69)/   0.000000000d0/, IS(4, 69)/ 0/,
     $  GM(4, 69)/   0.000000000d0/, QM(4, 69)/   0.000000000d0/
      Data MN(1, 70)/174/, AM(1, 70)/ 173.939000000d0/, IS(1, 70)/ 0/,
     $  GM(1, 70)/   0.000000000d0/, QM(1, 70)/   0.000000000d0/
      Data MN(2, 70)/171/, AM(2, 70)/ 170.936500000d0/, IS(2, 70)/ 1/,
     $  GM(2, 70)/   0.493670000d0/, QM(2, 70)/   0.000000000d0/
      Data MN(3, 70)/172/, AM(3, 70)/ 171.936600000d0/, IS(3, 70)/ 0/,
     $  GM(3, 70)/   0.000000000d0/, QM(3, 70)/   0.000000000d0/
      Data MN(4, 70)/173/, AM(4, 70)/ 172.938300000d0/, IS(4, 70)/ 5/,
     $  GM(4, 70)/  -0.679890000d0/, QM(4, 70)/ 280.000000000d0/
      Data MN(1, 71)/175/, AM(1, 71)/ 174.940900000d0/, IS(1, 71)/ 7/,
     $  GM(1, 71)/   2.232700000d0/, QM(1, 71)/ 349.000000000d0/
      Data MN(2, 71)/  0/, AM(2, 71)/   0.000000000d0/, IS(2, 71)/ 0/,
     $  GM(2, 71)/   0.000000000d0/, QM(2, 71)/   0.000000000d0/
      Data MN(3, 71)/  0/, AM(3, 71)/   0.000000000d0/, IS(3, 71)/ 0/,
     $  GM(3, 71)/   0.000000000d0/, QM(3, 71)/   0.000000000d0/
      Data MN(4, 71)/  0/, AM(4, 71)/   0.000000000d0/, IS(4, 71)/ 0/,
     $  GM(4, 71)/   0.000000000d0/, QM(4, 71)/   0.000000000d0/
      Data MN(1, 72)/180/, AM(1, 72)/ 179.946800000d0/, IS(1, 72)/ 0/,
     $  GM(1, 72)/   0.000000000d0/, QM(1, 72)/   0.000000000d0/
      Data MN(2, 72)/177/, AM(2, 72)/ 176.943500000d0/, IS(2, 72)/ 7/,
     $  GM(2, 72)/   0.793600000d0/, QM(2, 72)/ 336.500000000d0/
      Data MN(3, 72)/178/, AM(3, 72)/ 177.943900000d0/, IS(3, 72)/ 0/,
     $  GM(3, 72)/   0.000000000d0/, QM(3, 72)/   0.000000000d0/
      Data MN(4, 72)/179/, AM(4, 72)/ 178.946000000d0/, IS(4, 72)/ 9/,
     $  GM(4, 72)/  -0.640900000d0/, QM(4, 72)/ 379.300000000d0/
      Data MN(1, 73)/181/, AM(1, 73)/ 180.948000000d0/, IS(1, 73)/ 7/,
     $  GM(1, 73)/   2.371000000d0/, QM(1, 73)/ 317.000000000d0/
      Data MN(2, 73)/180/, AM(2, 73)/ 179.941500000d0/, IS(2, 73)/16/,
     $  GM(2, 73)/   0.000000000d0/, QM(2, 73)/   0.000000000d0/
      Data MN(3, 73)/  0/, AM(3, 73)/   0.000000000d0/, IS(3, 73)/ 0/,
     $  GM(3, 73)/   0.000000000d0/, QM(3, 73)/   0.000000000d0/
      Data MN(4, 73)/  0/, AM(4, 73)/   0.000000000d0/, IS(4, 73)/ 0/,
     $  GM(4, 73)/   0.000000000d0/, QM(4, 73)/   0.000000000d0/
      Data MN(1, 74)/184/, AM(1, 74)/ 183.951000000d0/, IS(1, 74)/ 0/,
     $  GM(1, 74)/   0.000000000d0/, QM(1, 74)/   0.000000000d0/
      Data MN(2, 74)/182/, AM(2, 74)/ 181.948300000d0/, IS(2, 74)/ 0/,
     $  GM(2, 74)/   0.000000000d0/, QM(2, 74)/   0.000000000d0/
      Data MN(3, 74)/183/, AM(3, 74)/ 182.950300000d0/, IS(3, 74)/ 1/,
     $  GM(3, 74)/   0.117785000d0/, QM(3, 74)/   0.000000000d0/
      Data MN(4, 74)/186/, AM(4, 74)/ 185.954300000d0/, IS(4, 74)/ 0/,
     $  GM(4, 74)/   0.000000000d0/, QM(4, 74)/   0.000000000d0/
      Data MN(1, 75)/187/, AM(1, 75)/ 186.956000000d0/, IS(1, 75)/ 5/,
     $  GM(1, 75)/   3.219700000d0/, QM(1, 75)/ 207.000000000d0/
      Data MN(2, 75)/185/, AM(2, 75)/ 184.953000000d0/, IS(2, 75)/ 5/,
     $  GM(2, 75)/   3.187100000d0/, QM(2, 75)/ 218.000000000d0/
      Data MN(3, 75)/  0/, AM(3, 75)/   0.000000000d0/, IS(3, 75)/ 0/,
     $  GM(3, 75)/   0.000000000d0/, QM(3, 75)/   0.000000000d0/
      Data MN(4, 75)/  0/, AM(4, 75)/   0.000000000d0/, IS(4, 75)/ 0/,
     $  GM(4, 75)/   0.000000000d0/, QM(4, 75)/   0.000000000d0/
      Data MN(1, 76)/190/, AM(1, 76)/ 189.958600000d0/, IS(1, 76)/ 0/,
     $  GM(1, 76)/   0.000000000d0/, QM(1, 76)/   0.000000000d0/
      Data MN(2, 76)/188/, AM(2, 76)/ 187.956000000d0/, IS(2, 76)/ 0/,
     $  GM(2, 76)/   0.000000000d0/, QM(2, 76)/   0.000000000d0/
      Data MN(3, 76)/189/, AM(3, 76)/ 188.958600000d0/, IS(3, 76)/ 3/,
     $  GM(3, 76)/   0.659933000d0/, QM(3, 76)/  85.600000000d0/
      Data MN(4, 76)/  0/, AM(4, 76)/   0.000000000d0/, IS(4, 76)/ 0/,
     $  GM(4, 76)/   0.000000000d0/, QM(4, 76)/   0.000000000d0/
      Data MN(1, 77)/193/, AM(1, 77)/ 192.963300000d0/, IS(1, 77)/ 3/,
     $  GM(1, 77)/   0.159200000d0/, QM(1, 77)/  75.100000000d0/
      Data MN(2, 77)/191/, AM(2, 77)/ 190.960900000d0/, IS(2, 77)/ 3/,
     $  GM(2, 77)/   0.146200000d0/, QM(2, 77)/  81.600000000d0/
      Data MN(3, 77)/  0/, AM(3, 77)/   0.000000000d0/, IS(3, 77)/ 0/,
     $  GM(3, 77)/   0.000000000d0/, QM(3, 77)/   0.000000000d0/
      Data MN(4, 77)/  0/, AM(4, 77)/   0.000000000d0/, IS(4, 77)/ 0/,
     $  GM(4, 77)/   0.000000000d0/, QM(4, 77)/   0.000000000d0/
      Data MN(1, 78)/195/, AM(1, 78)/ 194.964800000d0/, IS(1, 78)/ 1/,
     $  GM(1, 78)/   0.609500000d0/, QM(1, 78)/   0.000000000d0/
      Data MN(2, 78)/194/, AM(2, 78)/ 193.962800000d0/, IS(2, 78)/ 0/,
     $  GM(2, 78)/   0.000000000d0/, QM(2, 78)/   0.000000000d0/
      Data MN(3, 78)/196/, AM(3, 78)/ 195.965000000d0/, IS(3, 78)/ 0/,
     $  GM(3, 78)/   0.000000000d0/, QM(3, 78)/   0.000000000d0/
      Data MN(4, 78)/198/, AM(4, 78)/ 197.967500000d0/, IS(4, 78)/ 0/,
     $  GM(4, 78)/   0.000000000d0/, QM(4, 78)/   0.000000000d0/
      Data MN(1, 79)/197/, AM(1, 79)/ 196.966600000d0/, IS(1, 79)/ 3/,
     $  GM(1, 79)/   0.148159000d0/, QM(1, 79)/  54.700000000d0/
      Data MN(2, 79)/  0/, AM(2, 79)/   0.000000000d0/, IS(2, 79)/ 0/,
     $  GM(2, 79)/   0.000000000d0/, QM(2, 79)/   0.000000000d0/
      Data MN(3, 79)/  0/, AM(3, 79)/   0.000000000d0/, IS(3, 79)/ 0/,
     $  GM(3, 79)/   0.000000000d0/, QM(3, 79)/   0.000000000d0/
      Data MN(4, 79)/  0/, AM(4, 79)/   0.000000000d0/, IS(4, 79)/ 0/,
     $  GM(4, 79)/   0.000000000d0/, QM(4, 79)/   0.000000000d0/
      Data MN(1, 80)/202/, AM(1, 80)/ 201.970600000d0/, IS(1, 80)/ 0/,
     $  GM(1, 80)/   0.000000000d0/, QM(1, 80)/   0.000000000d0/
      Data MN(2, 80)/199/, AM(2, 80)/ 198.968300000d0/, IS(2, 80)/ 1/,
     $  GM(2, 80)/   0.505885000d0/, QM(2, 80)/   0.000000000d0/
      Data MN(3, 80)/200/, AM(3, 80)/ 199.968300000d0/, IS(3, 80)/ 0/,
     $  GM(3, 80)/   0.000000000d0/, QM(3, 80)/   0.000000000d0/
      Data MN(4, 80)/201/, AM(4, 80)/ 200.970300000d0/, IS(4, 80)/ 3/,
     $  GM(4, 80)/  -0.560220000d0/, QM(4, 80)/  38.600000000d0/
      Data MN(1, 81)/205/, AM(1, 81)/ 204.974500000d0/, IS(1, 81)/ 1/,
     $  GM(1, 81)/   1.638213000d0/, QM(1, 81)/   0.000000000d0/
      Data MN(2, 81)/203/, AM(2, 81)/ 202.972300000d0/, IS(2, 81)/ 1/,
     $  GM(2, 81)/   1.622257000d0/, QM(2, 81)/   0.000000000d0/
      Data MN(3, 81)/  0/, AM(3, 81)/   0.000000000d0/, IS(3, 81)/ 0/,
     $  GM(3, 81)/   0.000000000d0/, QM(3, 81)/   0.000000000d0/
      Data MN(4, 81)/  0/, AM(4, 81)/   0.000000000d0/, IS(4, 81)/ 0/,
     $  GM(4, 81)/   0.000000000d0/, QM(4, 81)/   0.000000000d0/
      Data MN(1, 82)/208/, AM(1, 82)/ 207.976600000d0/, IS(1, 82)/ 0/,
     $  GM(1, 82)/   0.000000000d0/, QM(1, 82)/   0.000000000d0/
      Data MN(2, 82)/204/, AM(2, 82)/ 203.973000000d0/, IS(2, 82)/ 0/,
     $  GM(2, 82)/   0.000000000d0/, QM(2, 82)/   0.000000000d0/
      Data MN(3, 82)/206/, AM(3, 82)/ 205.974500000d0/, IS(3, 82)/ 0/,
     $  GM(3, 82)/   0.000000000d0/, QM(3, 82)/   0.000000000d0/
      Data MN(4, 82)/207/, AM(4, 82)/ 206.975900000d0/, IS(4, 82)/ 1/,
     $  GM(4, 82)/   0.582190000d0/, QM(4, 82)/   0.000000000d0/
      Data MN(1, 83)/209/, AM(1, 83)/ 208.980400000d0/, IS(1, 83)/ 9/,
     $  GM(1, 83)/   4.110600000d0/, QM(1, 83)/ -51.600000000d0/
      Data MN(2, 83)/211/, AM(2, 83)/ 210.987300000d0/, IS(2, 83)/ 0/,
     $  GM(2, 83)/   0.000000000d0/, QM(2, 83)/   0.000000000d0/
      Data MN(3, 83)/  0/, AM(3, 83)/   0.000000000d0/, IS(3, 83)/ 0/,
     $  GM(3, 83)/   0.000000000d0/, QM(3, 83)/   0.000000000d0/
      Data MN(4, 83)/  0/, AM(4, 83)/   0.000000000d0/, IS(4, 83)/ 0/,
     $  GM(4, 83)/   0.000000000d0/, QM(4, 83)/   0.000000000d0/
      Data MN(1, 84)/209/, AM(1, 84)/ 208.982500000d0/, IS(1, 84)/ 1/,
     $  GM(1, 84)/   0.000000000d0/, QM(1, 84)/   0.000000000d0/
      Data MN(2, 84)/206/, AM(2, 84)/ 205.980500000d0/, IS(2, 84)/ 0/,
     $  GM(2, 84)/   0.000000000d0/, QM(2, 84)/   0.000000000d0/
      Data MN(3, 84)/207/, AM(3, 84)/ 206.981600000d0/, IS(3, 84)/ 0/,
     $  GM(3, 84)/   0.000000000d0/, QM(3, 84)/   0.000000000d0/
      Data MN(4, 84)/208/, AM(4, 84)/ 207.981300000d0/, IS(4, 84)/ 0/,
     $  GM(4, 84)/   0.000000000d0/, QM(4, 84)/   0.000000000d0/
      Data MN(1, 85)/211/, AM(1, 85)/ 210.987500000d0/, IS(1, 85)/ 0/,
     $  GM(1, 85)/   0.000000000d0/, QM(1, 85)/   0.000000000d0/
      Data MN(2, 85)/  0/, AM(2, 85)/   0.000000000d0/, IS(2, 85)/ 0/,
     $  GM(2, 85)/   0.000000000d0/, QM(2, 85)/   0.000000000d0/
      Data MN(3, 85)/  0/, AM(3, 85)/   0.000000000d0/, IS(3, 85)/ 0/,
     $  GM(3, 85)/   0.000000000d0/, QM(3, 85)/   0.000000000d0/
      Data MN(4, 85)/  0/, AM(4, 85)/   0.000000000d0/, IS(4, 85)/ 0/,
     $  GM(4, 85)/   0.000000000d0/, QM(4, 85)/   0.000000000d0/
      Data MN(1, 86)/222/, AM(1, 86)/ 222.017500000d0/, IS(1, 86)/ 0/,
     $  GM(1, 86)/   0.000000000d0/, QM(1, 86)/   0.000000000d0/
      Data MN(2, 86)/211/, AM(2, 86)/ 210.990600000d0/, IS(2, 86)/ 0/,
     $  GM(2, 86)/   0.000000000d0/, QM(2, 86)/   0.000000000d0/
      Data MN(3, 86)/210/, AM(3, 86)/ 209.989700000d0/, IS(3, 86)/ 0/,
     $  GM(3, 86)/   0.000000000d0/, QM(3, 86)/   0.000000000d0/
      Data MN(4, 86)/212/, AM(4, 86)/ 211.990700000d0/, IS(4, 86)/ 0/,
     $  GM(4, 86)/   0.000000000d0/, QM(4, 86)/   0.000000000d0/
      Data MN(1, 87)/223/, AM(1, 87)/ 223.019800000d0/, IS(1, 87)/ 3/,
     $  GM(1, 87)/   0.000000000d0/, QM(1, 87)/ 117.000000000d0/
      Data MN(2, 87)/212/, AM(2, 87)/ 211.996000000d0/, IS(2, 87)/ 0/,
     $  GM(2, 87)/   0.000000000d0/, QM(2, 87)/   0.000000000d0/
      Data MN(3, 87)/221/, AM(3, 87)/ 221.014200000d0/, IS(3, 87)/ 0/,
     $  GM(3, 87)/   0.000000000d0/, QM(3, 87)/   0.000000000d0/
      Data MN(4, 87)/  0/, AM(4, 87)/   0.000000000d0/, IS(4, 87)/ 0/,
     $  GM(4, 87)/   0.000000000d0/, QM(4, 87)/   0.000000000d0/
      Data MN(1, 88)/226/, AM(1, 88)/ 226.025400000d0/, IS(1, 88)/ 0/,
     $  GM(1, 88)/   0.000000000d0/, QM(1, 88)/   0.000000000d0/
      Data MN(2, 88)/  0/, AM(2, 88)/   0.000000000d0/, IS(2, 88)/ 0/,
     $  GM(2, 88)/   0.000000000d0/, QM(2, 88)/   0.000000000d0/
      Data MN(3, 88)/  0/, AM(3, 88)/   0.000000000d0/, IS(3, 88)/ 0/,
     $  GM(3, 88)/   0.000000000d0/, QM(3, 88)/   0.000000000d0/
      Data MN(4, 88)/  0/, AM(4, 88)/   0.000000000d0/, IS(4, 88)/ 0/,
     $  GM(4, 88)/   0.000000000d0/, QM(4, 88)/   0.000000000d0/
      Data MN(1, 89)/227/, AM(1, 89)/ 227.027800000d0/, IS(1, 89)/ 3/,
     $  GM(1, 89)/   1.100000000d0/, QM(1, 89)/ 170.000000000d0/
      Data MN(2, 89)/  0/, AM(2, 89)/   0.000000000d0/, IS(2, 89)/ 0/,
     $  GM(2, 89)/   0.000000000d0/, QM(2, 89)/   0.000000000d0/
      Data MN(3, 89)/  0/, AM(3, 89)/   0.000000000d0/, IS(3, 89)/ 0/,
     $  GM(3, 89)/   0.000000000d0/, QM(3, 89)/   0.000000000d0/
      Data MN(4, 89)/  0/, AM(4, 89)/   0.000000000d0/, IS(4, 89)/ 0/,
     $  GM(4, 89)/   0.000000000d0/, QM(4, 89)/   0.000000000d0/
      Data MN(1, 90)/232/, AM(1, 90)/ 232.038200000d0/, IS(1, 90)/ 0/,
     $  GM(1, 90)/   0.000000000d0/, QM(1, 90)/   0.000000000d0/
      Data MN(2, 90)/228/, AM(2, 90)/ 228.028700000d0/, IS(2, 90)/ 0/,
     $  GM(2, 90)/   0.000000000d0/, QM(2, 90)/   0.000000000d0/
      Data MN(3, 90)/229/, AM(3, 90)/ 229.031600000d0/, IS(3, 90)/ 0/,
     $  GM(3, 90)/   0.000000000d0/, QM(3, 90)/ 430.000000000d0/
      Data MN(4, 90)/230/, AM(4, 90)/ 230.033100000d0/, IS(4, 90)/ 0/,
     $  GM(4, 90)/   0.000000000d0/, QM(4, 90)/   0.000000000d0/
      Data MN(1, 91)/231/, AM(1, 91)/ 231.035900000d0/, IS(1, 91)/ 3/,
     $  GM(1, 91)/   2.010000000d0/, QM(1, 91)/-172.000000000d0/
      Data MN(2, 91)/234/, AM(2, 91)/ 234.043000000d0/, IS(2, 91)/ 0/,
     $  GM(2, 91)/   0.000000000d0/, QM(2, 91)/   0.000000000d0/
      Data MN(3, 91)/  0/, AM(3, 91)/   0.000000000d0/, IS(3, 91)/ 0/,
     $  GM(3, 91)/   0.000000000d0/, QM(3, 91)/   0.000000000d0/
      Data MN(4, 91)/  0/, AM(4, 91)/   0.000000000d0/, IS(4, 91)/ 0/,
     $  GM(4, 91)/   0.000000000d0/, QM(4, 91)/   0.000000000d0/
      Data MN(1, 92)/238/, AM(1, 92)/ 238.050800000d0/, IS(1, 92)/ 0/,
     $  GM(1, 92)/   0.000000000d0/, QM(1, 92)/   0.000000000d0/
      Data MN(2, 92)/234/, AM(2, 92)/ 234.040900000d0/, IS(2, 92)/ 0/,
     $  GM(2, 92)/   0.000000000d0/, QM(2, 92)/   0.000000000d0/
      Data MN(3, 92)/235/, AM(3, 92)/ 235.043900000d0/, IS(3, 92)/ 7/,
     $  GM(3, 92)/  -0.350000000d0/, QM(3, 92)/ 493.600000000d0/
      Data MN(4, 92)/236/, AM(4, 92)/ 236.045700000d0/, IS(4, 92)/ 0/,
     $  GM(4, 92)/   0.000000000d0/, QM(4, 92)/   0.000000000d0/
      Data MN(1, 93)/237/, AM(1, 93)/ 237.048000000d0/, IS(1, 93)/ 5/,
     $  GM(1, 93)/   3.140000000d0/, QM(1, 93)/ 388.600000000d0/
      Data MN(2, 93)/236/, AM(2, 93)/ 236.046600000d0/, IS(2, 93)/ 0/,
     $  GM(2, 93)/   0.000000000d0/, QM(2, 93)/   0.000000000d0/
      Data MN(3, 93)/  0/, AM(3, 93)/   0.000000000d0/, IS(3, 93)/ 0/,
     $  GM(3, 93)/   0.000000000d0/, QM(3, 93)/   0.000000000d0/
      Data MN(4, 93)/  0/, AM(4, 93)/   0.000000000d0/, IS(4, 93)/ 0/,
     $  GM(4, 93)/   0.000000000d0/, QM(4, 93)/   0.000000000d0/
      Data MN(1, 94)/242/, AM(1, 94)/ 242.058700000d0/, IS(1, 94)/ 0/,
     $  GM(1, 94)/   0.000000000d0/, QM(1, 94)/   0.000000000d0/
      Data MN(2, 94)/239/, AM(2, 94)/ 239.052200000d0/, IS(2, 94)/ 0/,
     $  GM(2, 94)/   0.000000000d0/, QM(2, 94)/   0.000000000d0/
      Data MN(3, 94)/240/, AM(3, 94)/ 240.054000000d0/, IS(3, 94)/ 0/,
     $  GM(3, 94)/   0.000000000d0/, QM(3, 94)/   0.000000000d0/
      Data MN(4, 94)/  0/, AM(4, 94)/   0.000000000d0/, IS(4, 94)/ 0/,
     $  GM(4, 94)/   0.000000000d0/, QM(4, 94)/   0.000000000d0/
      Data MN(1, 95)/243/, AM(1, 95)/ 243.061400000d0/, IS(1, 95)/ 5/,
     $  GM(1, 95)/   1.610000000d0/, QM(1, 95)/ 421.000000000d0/
      Data MN(2, 95)/241/, AM(2, 95)/ 241.056700000d0/, IS(2, 95)/ 0/,
     $  GM(2, 95)/   0.000000000d0/, QM(2, 95)/   0.000000000d0/
      Data MN(3, 95)/  0/, AM(3, 95)/   0.000000000d0/, IS(3, 95)/ 0/,
     $  GM(3, 95)/   0.000000000d0/, QM(3, 95)/   0.000000000d0/
      Data MN(4, 95)/  0/, AM(4, 95)/   0.000000000d0/, IS(4, 95)/ 0/,
     $  GM(4, 95)/   0.000000000d0/, QM(4, 95)/   0.000000000d0/
      Data MN(1, 96)/246/, AM(1, 96)/ 246.067400000d0/, IS(1, 96)/ 0/,
     $  GM(1, 96)/   0.000000000d0/, QM(1, 96)/   0.000000000d0/
      Data MN(2, 96)/245/, AM(2, 96)/ 245.065300000d0/, IS(2, 96)/ 0/,
     $  GM(2, 96)/   0.000000000d0/, QM(2, 96)/   0.000000000d0/
      Data MN(3, 96)/  0/, AM(3, 96)/   0.000000000d0/, IS(3, 96)/ 0/,
     $  GM(3, 96)/   0.000000000d0/, QM(3, 96)/   0.000000000d0/
      Data MN(4, 96)/  0/, AM(4, 96)/   0.000000000d0/, IS(4, 96)/ 0/,
     $  GM(4, 96)/   0.000000000d0/, QM(4, 96)/   0.000000000d0/
      Data MN(1, 97)/247/, AM(1, 97)/ 247.070200000d0/, IS(1, 97)/ 0/,
     $  GM(1, 97)/   0.000000000d0/, QM(1, 97)/   0.000000000d0/
      Data MN(2, 97)/  0/, AM(2, 97)/   0.000000000d0/, IS(2, 97)/ 0/,
     $  GM(2, 97)/   0.000000000d0/, QM(2, 97)/   0.000000000d0/
      Data MN(3, 97)/  0/, AM(3, 97)/   0.000000000d0/, IS(3, 97)/ 0/,
     $  GM(3, 97)/   0.000000000d0/, QM(3, 97)/   0.000000000d0/
      Data MN(4, 97)/  0/, AM(4, 97)/   0.000000000d0/, IS(4, 97)/ 0/,
     $  GM(4, 97)/   0.000000000d0/, QM(4, 97)/   0.000000000d0/
      Data MN(1, 98)/249/, AM(1, 98)/ 249.074800000d0/, IS(1, 98)/ 0/,
     $  GM(1, 98)/   0.000000000d0/, QM(1, 98)/   0.000000000d0/
      Data MN(2, 98)/250/, AM(2, 98)/ 250.076600000d0/, IS(2, 98)/ 0/,
     $  GM(2, 98)/   0.000000000d0/, QM(2, 98)/   0.000000000d0/
      Data MN(3, 98)/  0/, AM(3, 98)/   0.000000000d0/, IS(3, 98)/ 0/,
     $  GM(3, 98)/   0.000000000d0/, QM(3, 98)/   0.000000000d0/
      Data MN(4, 98)/  0/, AM(4, 98)/   0.000000000d0/, IS(4, 98)/ 0/,
     $  GM(4, 98)/   0.000000000d0/, QM(4, 98)/   0.000000000d0/
      Data MN(1, 99)/252/, AM(1, 99)/ 252.082900000d0/, IS(1, 99)/ 0/,
     $  GM(1, 99)/   0.000000000d0/, QM(1, 99)/   0.000000000d0/
      Data MN(2, 99)/253/, AM(2, 99)/ 253.084700000d0/, IS(2, 99)/ 0/,
     $  GM(2, 99)/   0.000000000d0/, QM(2, 99)/ 670.000000000d0/
      Data MN(3, 99)/254/, AM(3, 99)/ 254.088100000d0/, IS(3, 99)/ 0/,
     $  GM(3, 99)/   0.000000000d0/, QM(3, 99)/   0.000000000d0/
      Data MN(4, 99)/  0/, AM(4, 99)/   0.000000000d0/, IS(4, 99)/ 0/,
     $  GM(4, 99)/   0.000000000d0/, QM(4, 99)/   0.000000000d0/
      Data MN(1,100)/252/, AM(1,100)/ 252.082700000d0/, IS(1,100)/ 0/,
     $  GM(1,100)/   0.000000000d0/, QM(1,100)/   0.000000000d0/
      Data MN(2,100)/250/, AM(2,100)/ 250.079500000d0/, IS(2,100)/ 0/,
     $  GM(2,100)/   0.000000000d0/, QM(2,100)/   0.000000000d0/
      Data MN(3,100)/254/, AM(3,100)/ 254.087000000d0/, IS(3,100)/ 0/,
     $  GM(3,100)/   0.000000000d0/, QM(3,100)/   0.000000000d0/
      Data MN(4,100)/  0/, AM(4,100)/   0.000000000d0/, IS(4,100)/ 0/,
     $  GM(4,100)/   0.000000000d0/, QM(4,100)/   0.000000000d0/
      Data MN(1,101)/255/, AM(1,101)/ 255.090600000d0/, IS(1,101)/ 0/,
     $  GM(1,101)/   0.000000000d0/, QM(1,101)/   0.000000000d0/
      Data MN(2,101)/  0/, AM(2,101)/   0.000000000d0/, IS(2,101)/ 0/,
     $  GM(2,101)/   0.000000000d0/, QM(2,101)/   0.000000000d0/
      Data MN(3,101)/  0/, AM(3,101)/   0.000000000d0/, IS(3,101)/ 0/,
     $  GM(3,101)/   0.000000000d0/, QM(3,101)/   0.000000000d0/
      Data MN(4,101)/  0/, AM(4,101)/   0.000000000d0/, IS(4,101)/ 0/,
     $  GM(4,101)/   0.000000000d0/, QM(4,101)/   0.000000000d0/
      Data MN(1,102)/259/, AM(1,102)/ 259.101000000d0/, IS(1,102)/ 0/,
     $  GM(1,102)/   0.000000000d0/, QM(1,102)/   0.000000000d0/
      Data MN(2,102)/  0/, AM(2,102)/   0.000000000d0/, IS(2,102)/ 0/,
     $  GM(2,102)/   0.000000000d0/, QM(2,102)/   0.000000000d0/
      Data MN(3,102)/  0/, AM(3,102)/   0.000000000d0/, IS(3,102)/ 0/,
     $  GM(3,102)/   0.000000000d0/, QM(3,102)/   0.000000000d0/
      Data MN(4,102)/  0/, AM(4,102)/   0.000000000d0/, IS(4,102)/ 0/,
     $  GM(4,102)/   0.000000000d0/, QM(4,102)/   0.000000000d0/
      Data MN(1,103)/262/, AM(1,103)/ 262.109700000d0/, IS(1,103)/ 0/,
     $  GM(1,103)/   0.000000000d0/, QM(1,103)/   0.000000000d0/
      Data MN(2,103)/  0/, AM(2,103)/   0.000000000d0/, IS(2,103)/ 0/,
     $  GM(2,103)/   0.000000000d0/, QM(2,103)/   0.000000000d0/
      Data MN(3,103)/  0/, AM(3,103)/   0.000000000d0/, IS(3,103)/ 0/,
     $  GM(3,103)/   0.000000000d0/, QM(3,103)/   0.000000000d0/
      Data MN(4,103)/  0/, AM(4,103)/   0.000000000d0/, IS(4,103)/ 0/,
     $  GM(4,103)/   0.000000000d0/, QM(4,103)/   0.000000000d0/
      Data MN(1,104)/261/, AM(1,104)/ 261.108700000d0/, IS(1,104)/ 0/,
     $  GM(1,104)/   0.000000000d0/, QM(1,104)/   0.000000000d0/
      Data MN(2,104)/  0/, AM(2,104)/   0.000000000d0/, IS(2,104)/ 0/,
     $  GM(2,104)/   0.000000000d0/, QM(2,104)/   0.000000000d0/
      Data MN(3,104)/  0/, AM(3,104)/   0.000000000d0/, IS(3,104)/ 0/,
     $  GM(3,104)/   0.000000000d0/, QM(3,104)/   0.000000000d0/
      Data MN(4,104)/  0/, AM(4,104)/   0.000000000d0/, IS(4,104)/ 0/,
     $  GM(4,104)/   0.000000000d0/, QM(4,104)/   0.000000000d0/
      Data MN(1,105)/262/, AM(1,105)/ 262.114100000d0/, IS(1,105)/ 0/,
     $  GM(1,105)/   0.000000000d0/, QM(1,105)/   0.000000000d0/
      Data MN(2,105)/  0/, AM(2,105)/   0.000000000d0/, IS(2,105)/ 0/,
     $  GM(2,105)/   0.000000000d0/, QM(2,105)/   0.000000000d0/
      Data MN(3,105)/  0/, AM(3,105)/   0.000000000d0/, IS(3,105)/ 0/,
     $  GM(3,105)/   0.000000000d0/, QM(3,105)/   0.000000000d0/
      Data MN(4,105)/  0/, AM(4,105)/   0.000000000d0/, IS(4,105)/ 0/,
     $  GM(4,105)/   0.000000000d0/, QM(4,105)/   0.000000000d0/
      Data MN(1,106)/266/, AM(1,106)/ 266.121900000d0/, IS(1,106)/ 0/,
     $  GM(1,106)/   0.000000000d0/, QM(1,106)/   0.000000000d0/
      Data MN(2,106)/  0/, AM(2,106)/   0.000000000d0/, IS(2,106)/ 0/,
     $  GM(2,106)/   0.000000000d0/, QM(2,106)/   0.000000000d0/
      Data MN(3,106)/  0/, AM(3,106)/   0.000000000d0/, IS(3,106)/ 0/,
     $  GM(3,106)/   0.000000000d0/, QM(3,106)/   0.000000000d0/
      Data MN(4,106)/  0/, AM(4,106)/   0.000000000d0/, IS(4,106)/ 0/,
     $  GM(4,106)/   0.000000000d0/, QM(4,106)/   0.000000000d0/
      Data MN(1,107)/264/, AM(1,107)/ 264.124700000d0/, IS(1,107)/ 0/,
     $  GM(1,107)/   0.000000000d0/, QM(1,107)/   0.000000000d0/
      Data MN(2,107)/  0/, AM(2,107)/   0.000000000d0/, IS(2,107)/ 0/,
     $  GM(2,107)/   0.000000000d0/, QM(2,107)/   0.000000000d0/
      Data MN(3,107)/  0/, AM(3,107)/   0.000000000d0/, IS(3,107)/ 0/,
     $  GM(3,107)/   0.000000000d0/, QM(3,107)/   0.000000000d0/
      Data MN(4,107)/  0/, AM(4,107)/   0.000000000d0/, IS(4,107)/ 0/,
     $  GM(4,107)/   0.000000000d0/, QM(4,107)/   0.000000000d0/
      Data MN(1,108)/277/, AM(1,108)/   0.000000000d0/, IS(1,108)/ 0/,
     $  GM(1,108)/   0.000000000d0/, QM(1,108)/   0.000000000d0/
      Data MN(2,108)/  0/, AM(2,108)/   0.000000000d0/, IS(2,108)/ 0/,
     $  GM(2,108)/   0.000000000d0/, QM(2,108)/   0.000000000d0/
      Data MN(3,108)/  0/, AM(3,108)/   0.000000000d0/, IS(3,108)/ 0/,
     $  GM(3,108)/   0.000000000d0/, QM(3,108)/   0.000000000d0/
      Data MN(4,108)/  0/, AM(4,108)/   0.000000000d0/, IS(4,108)/ 0/,
     $  GM(4,108)/   0.000000000d0/, QM(4,108)/   0.000000000d0/
      Data MN(1,109)/268/, AM(1,109)/ 268.138800000d0/, IS(1,109)/ 0/,
     $  GM(1,109)/   0.000000000d0/, QM(1,109)/   0.000000000d0/
      Data MN(2,109)/  0/, AM(2,109)/   0.000000000d0/, IS(2,109)/ 0/,
     $  GM(2,109)/   0.000000000d0/, QM(2,109)/   0.000000000d0/
      Data MN(3,109)/  0/, AM(3,109)/   0.000000000d0/, IS(3,109)/ 0/,
     $  GM(3,109)/   0.000000000d0/, QM(3,109)/   0.000000000d0/
      Data MN(4,109)/  0/, AM(4,109)/   0.000000000d0/, IS(4,109)/ 0/,
     $  GM(4,109)/   0.000000000d0/, QM(4,109)/   0.000000000d0/
C
      If(IAnI.ge.0.and.IAnI.le.MaxAn) then
        If(MNI.le.0) then
          JUse = Max(-MNI,1)
          If(JUse.gt.MaxIso) then
            JUse = 0
          else if(MN(JUse,IAnI).eq.0) then
            JUse = 0
            endIf
        else
          JUse = 0
          Do 10 J = 1, MaxIso
            If(MN(J,IAnI).eq.MNI) JUse = J
   10       Continue
          endIf
      else
        JUse = -1
        endIf
      If(JUse.le.0) then
        If(MNI.ge.0) then
          RMassX = Float(MNI)
        else
          RMassX = Zero
          endIf
        ISpinX = 0
        QMomX = Zero
        GFacX = Zero
      else
        RMassX = AM(JUse,IAnI)
        ISpinX = IS(JUse,IAnI)
        QMomX = QM(JUse,IAnI)
        GFacX = GM(JUse,IAnI)
        endIf
      If(JUse.gt.0) then
        MNO = MN(JUse,IAnI)
      else if(MNI.gt.0) then
        MNO = MNI
      else
        MNO = 0
        endIf
      If(RMass.eq.Zero.or.Abs(RMass-RMassX).le.Small) RMass = RMassX
      If(ISpin.eq.-1) ISpin = ISpinX
      If(GFac.eq.Zero) GFac = GFacX
      If(QMom.eq.Zero) QMom = QMomX
      Return
      End
*Deck FixPh1
      Subroutine FixPh1(NBDim,NBasis,NC,NO,NV,NODim,NVDim,CMO,VO,VV)
      Implicit Real*8(A-H,O-Z)
C
C     Fix the phases of the NC+NO+NV orbitals in the specified rwf by
C     forcing the largest coefficient of each orbital to be positive.
C     If the dimension allow, VO and VV have transformation matrices
C     between occupieds and virtuals which are transformed in the same
C     way.  If NO=NV=0, then VO and VV are not used.
C
      Dimension CMO(NBDim,*), VO(NODim,*), VV(NVDim,*)
      Save Zero
      Data Zero/0.0d0/
C
      NOrbs = NO + NV + NC
      Do 20 I = 1, NOrbs
        JMax = IFixPh(1,NBasis,CMO(1,I),1)
        If(CMO(JMax,I).lt.Zero) then
          Call ANeg(NBasis,CMO(1,I),CMO(1,I))
          If(I.gt.NC.and.I.le.(NO+NC)) then
            Call ANeg(NO,VO(1,I-NC),VO(1,I-NC))
          else if(I.gt.(NC+NO)) then
            Call ANeg(NV,VV(1,I-NO-NC),VV(1,I-NO-NC))
            endIf
          endIf
   20   Continue
      Return
      End
*Deck ACos1
      Function ACos1(Arg)
      Implicit Real*8(A-H,O-Z)
C
C     This little function handles calling ACos when there
C     is concern about an argument out of range due to round
C     off error.  It truncates out of range values to +/- 1.
C
      Save One
      Data One/1.0d0/
C
      A = Arg
      If(Abs(A).gt.One) A = Sign(One,A)
      ACos1 = ACos(A)
      Return
      End
*Deck IArMx1
      Function IArMx1(A,N,IfAbs,IMax)
      Implicit Integer(A-Z)
C
C     This function returns the maximum element of an array.  IfAbs
C     determines whether the absolute values are to be compared.
C     IMax is set to the index of the largest element.
C
      Dimension A(*)
      Logical IfAbs
C
      If(N.lt.1) then
        IArMx1 = 0
        IMax = 0
      else
        IMax = 1
        If(IfAbs) then
          AM = Abs(A(1))
          Do 10 I = 2, N
            AV = Abs(A(I))
            If(AV.gt.AM) then
              IMax = I
              AM = AV
              endIf
   10       Continue
          IArMx1 = AM
        else
          Do 20 I = 2, N
            If(A(I).gt.A(IMax)) IMax = I
   20       Continue
          IArMx1 = A(IMax)
          endIf
        endIf
      Return
      End
*Deck IClear
      Subroutine IClear(N,IA)
      Integer N,IA,I
C     
C     Clear N elements of IA
C
      Dimension IA(*)
      Do 10 I=1,N
       IA(I)=0  
  10  Continue
      Return
      End
*Deck IrMax1
      Integer Function IrMax1(IA,N,IfAbs,IMax)
      Implicit Integer(A-Z)
C
C     This function returns the maximum element of an integer array.
C     IfAbs determines whether the absolute values are to be compared.
C     IMax is set to the index of the largest element.
C
      Dimension IA(*)
      Logical IfAbs
C
      If(N.lt.1) then
        IrMax1 = 0 
        IMax = 0
      else
        IMax = 1
        If(IfAbs) then
          IAM = Abs(IA(1))
          Do 10 I = 2, N
           IAV = Abs(IA(I))
            If(IAV.gt.IAM) then
              IMax = I
              IAM = IAV
              endIf
   10       Continue
          IrMax1 = IAM
        else
          Do 20 I = 2, N
            If(IA(I).gt.IA(IMax)) IMax = I
   20       Continue
          IrMax1 = IA(IMax)
          endIf
        endIf
      Return 
      End
*Deck IrMin1
      Integer Function IrMin1(IA,N,IfAbs,IMin)
      Implicit Integer(A-Z)
C       
C     This function returns the minimum element of an integer array.
C     IfAbs determines whether the absolute values are to be compared.
C     IMin is set to the index of the largest element.
C       
      Dimension IA(*)
      Logical IfAbs
C  
      If(N.lt.1) then
        IrMin1 = 0
        IMin = 0
      else
        IMin = 1
        If(IfAbs) then
          IAM = Abs(IA(1))
          Do 10 I = 2, N
            IAV = Abs(IA(I)) 
            If(IAV.lt.IAM) then
              IMin = I
              IAM = IAV
              endIf
   10       Continue
          IrMin1 = IAM
        else
          Do 20 I = 2, N
            If(IA(I).lt.IA(IMin)) IMin = I
   20       Continue
          IrMin1 = IA(IMin)
          endIf
        endIf
      Return
      End
*Deck IFixPh
      Integer Function IFixPh(NRI,N,X,LDX)
      Implicit Real*8(A-H,O-Z)
C
C     Return the index of an entry in a row or column of X to scale in
C     order to impose a phase convention.
C
      Dimension X(NRI,LDX,*)
      Save Thresh
      Data Thresh/1.d-6/
C
      If(NRI.eq.1) then
        CMax = Abs(X(1,1,1)) + Thresh
        JMax = 1
        Do 10 J = 2, N
          AC = Abs(X(1,1,J))
          If(AC.gt.CMax) then
            CMax = AC + Thresh
            JMax = J
            endIf
   10     Continue
      else
        CMax = Sqrt(X(1,1,1)**2+X(2,1,1)**2) + Thresh
        JMax = 1
        Do 20 J = 2, N
          AC = Sqrt(X(1,1,J)**2+X(2,1,J)**2)
          If(AC.gt.CMax) then
            CMax = AC + Thresh
            JMax = J
            endIf
   20     Continue
        endIf
      IFixPh = JMax
      Return
      End
*Deck IGet10
      Function IGet10(Num,Digit)
      Implicit Integer(A-Z)
C
C     Return one decimal digit from a number.
C
      IGet10 = Mod(Num/10**Digit,10)
      Return
      End
*Deck IPopVC
      Function IPopVc(N,A,Thresh,B)
      Implicit Real*8(A-H,O-Z)
C
C     Copy vector A to B, set all elements less than Thresh to 0
C     and return the number of zeroes.
C
      Dimension A(1), B(1)
      Save Zero
      Data Zero/0.0d0/
C
      ICount = 0
      Do 10 I = 1, N
        B(I) = A(I)
        If(Abs(B(I)).le.Thresh) then
          ICount = ICount + 1
          B(I) = Zero
          endIf
   10   Continue
      IPopVC = ICount
      Return
      End
*Deck Inv1
      Logical Function Inv1(A,N,IS,IAD1,IAD2,D,MDM,Det)
      Implicit Real*8(A-H,O-Z)
C
C     Inversion of square matrix a by means of the gauss-jordan
C     algorithm.  The value of the function is .true. if the
C     matrix is invertable.
C
C     April 72/RS9B, April 84 and February 85 Mike Frisch
C     Nov 94 Keith and Frisch
C
      Dimension A(MDM,MDM),IS(2,MDM),IAD1(MDM),IAD2(MDM),D(MDM)
      Save Zero, One, Small
      Data Zero/0.0d0/, One/1.0d0/, Small/1.0d-22/
C
      Call IClear(2*N,IS)
      ISt = Mod(N,4) + 1
      Do 500 IMA = 1, N
        B = Zero
        Do 20 L = 1, N
          If(IS(1,L).ne.1) then
            Do 10 M = 1, N
              If(IS(2,M).ne.1) then
                E = Abs(A(M,L))
                If(E.ge.B) then
                  B = E
                  I = L
                  K = M
                  endIf
                endIf
   10         Continue
            endIf
   20     Continue
        IS(1,I) = 1
        IS(2,K) = 1
        IAD1(K) = I
        IAD2(I) = K
        B = A(K,I)
        AbsB = Abs(B)
        If(AbsB.lt.Small) then
          Det = B
          Inv1 = .False.
          Return
          endIf
        T = -One / B
        Do 110 L = 1, N
  110     A(L,I) = T*A(L,I)
        T = -T
        A(K,I) = T
        Do 120 M = 1, N
  120     D(M) = A(M,I)
        Do 210 L = 1, (ISt-1)
          S = A(K,L)
          Do 200 M = 1, N
  200       A(M,L) = A(M,L) + S*D(M)
  210     A(K,L) = S*T
        Do 230 L = ISt, N, 4
          S = A(K,L)
          S1 = A(K,L+1)
          S2 = A(K,L+2)
          S3 = A(K,L+3)
          Do 220 M = 1, N
            A(M,L) = A(M,L) + S*D(M)
            A(M,L+1) = A(M,L+1) + S1*D(M)
            A(M,L+2) = A(M,L+2) + S2*D(M)
  220       A(M,L+3) = A(M,L+3) + S3*D(M)
          A(K,L) = S*T
          A(K,L+1) = S1*T
          A(K,L+2) = S2*T
  230     A(K,L+3) = S3*T
        Do 500 M = 1, N
  500     A(M,I) = D(M)
      Do 610 L = 1, N
        Do 600 J = 1, N
          K = IAD1(J)
  600     D(J) = A(L,K)
        Do 610 J = 1, N
  610     A(L,J) = D(J)
      Do 710 L = 1, N
        Do 700 J = 1, N
          K = IAD2(J)
  700     D(J) = A(K,L)
        Do 710 J = 1, N
  710     A(J,L) = D(J)
      Inv1 = .True.
      Det = B
      Return
      End
*Deck ISeq
      Subroutine ISeq(N,IOff,Inc,IX)
      Implicit Integer(A-Z)
C
C     Set N elements of IX to a sequence of integers
C
      Dimension IX(N)
C
      Do 10 I = 1, N
   10   IX(I) = IOff + I*Inc
      Return
      End
*Deck LClear
      Subroutine LClear(N,Log)
      Implicit Integer(A-Z)
C
C     Set all elements of logical array Log to .False.
C
      Logical Log(1)
C
      If(N.lt.1) Return
      Do 10 I = 1, N
   10   Log(I) = .False.
      Return
      End
*Deck LIdxTM
      Integer function LIdxTM(IOut,N,i,j,k,l)
      Implicit Integer(A-Z)
C
C Linear InDeX for Triangular Matrix
C     Give a unique index for a lower triangular matrix stored in linear
C     form
C
C Input:
C     N: dimension of the matrix (4 at most)
C     i,j,k,l: indexes of the element in the triangular matrix
C Output:
C     LIdxTM: index of the element in linear storage
C
C     Dimension
      Integer N, NMAX
      Parameter (NMAX=4)
C     Input
      Integer i, j, k, l
C     Local
      Integer LI(NMAX), a, La
C
      If(N.lt.0) then
C       Wrong dimension gives negative results
        LIdxTM = -1
      else if(N.eq.0) then
        LIdxTM = 1
      else if(N.eq.1) then
        LIdxTM = i
      else if(N.eq.2) then
C       Simple case of bidimensional matrix
        If(i.le.j) then
          LIdxTM = (j-1)*j/2+i
        else
          LIdxTM = (i-1)*i/2+j
          endIf
      else if(N.le.NMAX) then
        LI(1) = i
        LI(2) = j
        LI(3) = k
        If(N.ge.4) LI(4) = l
        a = 1
   10   If(LI(a).gt.LI(a+1)) then
          La      = LI(a)
          LI(a)   = LI(a+1)
          LI(a+1) = La
          a       = 0
          endIf
        If(a.lt.N-1) then
          a = a + 1
          Goto 10
          endIf
        LIdxTM = (LI(3)-1)*LI(3)*(LI(3)+1)/6+(LI(2)-1)*LI(2)/2+LI(1)
        If(N.ge.4)
     $    LIdxTM = LIdxTM + (LI(4)-1)*LI(4)*(LI(4)+1)*(LI(4)+2)/24
      else
C        Call GauErr('Dimension of the matrix is too high for LIdxTM')
       write(IOut,100)
  100  format('Dimension of the matrix is too high for LIdxTM')
       stop 
        endIf
      Return
      End
*Deck Linear
      Subroutine Linear(A,B,ND,N)
      Implicit Real*8(A-H,O-Z)
C
C     Places symmetric square array in linear form.
C
      Dimension A(ND,ND), B(*)
C
      Do 10 J = 1, N
       JJ = (J*(J-1))/2
       Do 20 I = 1, J
        B(JJ+I) = A(I,J)
   20  Continue
   10 Continue  
      Return
      End
*Deck LINTRP
      SUBROUTINE LINTRP(N,S,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     THIS ROUTINE LINEARLY INTERPOLATES BETWEEN ARRAYS A AND B,
C     FORMING
C
C     C(I) = S*A(I) + (1-S)*B(I)
C
      DIMENSION A(1), B(1), C(1)
      Save ONE
      DATA ONE/1.0D0/
C
      IF(N.LT.1) RETURN
      DO 10 I = 1, N
   10     C(I) = S*A(I) + (ONE-S)*B(I)
      RETURN
      END
*Deck LinUpC 
      Subroutine LinUpc(InStr,OutStr)
      Implicit Integer(A-Z)
C    
C     Translate a character string to upper case.
C    
      Character*(*) InStr, OutStr
C     
      IUA = IChar('A')
      ILA = IChar('a')
      ILZ = IChar('z')
      LenO = Min(Len(InStr),Len(OutStr))
      Do 10 I = 1, LenO
        ICI = IChar(InStr(I:I))
        If(ICI.ge.ILA.and.ICI.le.ILZ) ICI = ICI + IUA - ILA
        OutStr(I:I) = Char(ICI)
   10 Continue
      If(LenO.lt.Len(OutStr)) OutStr(LenO+1:) = ' '
      Return
      End
*Deck LLinBl
      Integer Function LLinBl(cline,find)
      Logical OK
C
C     subroutine for finding the first occurrence of character 'find' in a line
C     
      character*(*) cline, find 
      llinbl=-1
      OK = .False.
      do 10 i=1,len(cline),1
       if(OK) goto 10
       if(cline(i:i).eq.find) then
        llinbl=i
        OK = .true.
       EndIf
   10 continue
      return
      end  
*Deck LLinNB
      Integer Function LLinNB(cline,find)
      Logical OK
C
C     subroutine for finding the first occurrence of a character different from 'find' in a line
C
      character*(*) cline, find
      llinbl=-1
      OK = .False.
      do 10 i=1,len(cline),1
       if(OK) goto 10
       if(cline(i:i).ne.find) then
        llinNB=i
        OK = .true.
       EndIf
   10 continue
      return
      end
*Deck LlinCl 
      Subroutine LLinCl(cline)
      implicit real*8(a-h,o-z)
c
c     subroutine for clearing a string
c
      character*(*) cline, blank
      parameter (blank=' ')
c
      do 10 i = 1, len(cline), 1
       cline(i:i)=blank
   10  continue
      return
      end
*Deck LRmCom
      Subroutine LRmCom(cline)
      implicit real*8(a-h,o-z)
c     
c     subroutine for removing ',' from a string
c
      character*(*) cline, blank, coma
      parameter (blank=' ',coma=',')
c
      do 10 i = 1, len(cline), 1
       if(cline(i:i).eq.coma) cline(i:i)=blank
   10  continue
      return
      end
*Deck LRmEqu
      Subroutine LRmEqu(cline)
      implicit real*8(a-h,o-z)
c     
c     subroutine for removing '=' from a string
c
      character*(*) cline, blank, equal
      parameter (blank=' ',equal='=')
c
      do 10 i = 1, len(cline), 1
       if(cline(i:i).eq.equal) cline(i:i)=blank
   10  continue
      return
      end
*Deck LinEnd
      integer function linend(cline)
      implicit real*8(a-h,o-z)
c
c     function which returns the length of a character string,
c     excluding trailing blanks, tabs, nulls, and carriage returns.
c
      character*(*) cline, blank, null*1, tab*1, cr*1
      parameter (blank=' ')
c
      null = char(0)
      tab = char(9)
      cr = char(13)
      Do 5 i = len(cline), 1, -1
        if(cline(i:i).ne.blank.and.cline(i:i).ne.cr.and.
     $    cline(i:i).ne.null.and.
     $    cline(i:i).ne.tab) goto 10
    5     continue
      linend = 0
      return
   10 linend = i
      return
      end
*Deck LineSt
      integer function LineSt(cline)
      implicit real*8(a-h,o-z)
c
c     function which returns the start of a character string, skipping
c     leading blanks, tabs and nulls.
c
      character*(*) cline, blank
      character*1 null, tab
      parameter (blank=' ')
c
      null = char(0)
      tab = char(9)
      do 5 i = 1, len(cline), 1
          if(cline(i:i).ne.blank.and.
     $    cline(i:i).ne.null.and.
     $    cline(i:i).ne.tab) goto 10
    5     continue
      linest = 0
      return
   10 linest = i
      return
      end
*Deck LnK1e
      Subroutine Lnk1e(IDum)
      Integer IDum
      If(IDum.eq.0) then
       Write(*,'('' Job Killed'')')
       STOP
      EndIf
      Return
      End
*Deck LocC
      Function LocC(Char,Line)
      Implicit Integer(A-Z)
C
C     Return the position of char in line, or 0 if char is not in line.
C
      Character*(*) Char, Line
C
      LC = Len(Char)
      Do 10 I = 1, Len(Line)
        If(Line(I:I+LC-1).eq.Char) then
          LocC = I
          Return
          endIf
   10   Continue
      LocC = 0
      Return
      End
*Deck MatMP1
      Subroutine MatMP1(IThrsh,IOpt,ISign,LL,L,M,N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     This routine computes the matrix product
C
C     C(I,J) = C(I,J) +/- SUM(K) A(K,I) * B(K,J)
C
C     IOpt is 0/1 to initialize/add to C.
C     ISign is +1/-1 to add/subtract the product.
C
C     C must be distinct from both A and B.
C
      Real*8 MDCutO
      Dimension A(M,L), B(M,N), C(LL,N)
      Save Zero
      Data Zero/0.0d0/
C
      If(IOpt.eq.0) then
       Do 10 J = 1, N
        Do 10 I = 1, L
   10  C(I,J) = Zero
      endIf
      If(L.le.0.or.M.le.0.or.N.le.0) Return
      Sign = Float(ISign)
      Thresh = MDCutO(IThrsh)
      Do 50 J = 1, N
       Do 50 K = 1, M
        S = Sign * B(K,J)
        If(Abs(S).gt.Thresh) then
         Do 40 I = 1, L
   40    C(I,J) = C(I,J) + S*A(K,I)
        EndIf
   50  Continue
      Return
      End
*Deck MatMpy
      Subroutine MatMpy(L,M,N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     This routine computes the matrix product:
C
C     C(I,J) = Sum(K) A(I,K) * B(K,J)
C
C     C must be distinct from both A and B.
C
      Real*8 MDCutO
      Dimension A(L,M), B(M,N), C(L,N)
      Save Zero
      Data Zero/0.0d0/
C
      Do 10 J = 1, N
        Do 10 I = 1, L
   10     C(I,J) = Zero
        Do 20 J = 1, N
          Do 20 K = 1, M
            S = B(K,J)
            Do 20 I = 1, L
   20         C(I,J) = C(I,J) + A(I,K) * S
      Return
      End
*Deck MDCutO
      Function MDCutO(I)
      Implicit Real*8(A-H,O-Z)
      Real*8 MDCutO
C
C     Return a machine-dependent cutoff.
C     I = 0 ... A value which should be considered small compared to
C               1.0, taking into account round-off error.
C     I = 1 ... A value suitable for a general test for zero values.
C     I = 2 ... A value suitable for thresholding in rotation from
C               local atomic axes.
C
      If(I.eq.0) then
       MDCutO = 1.d-12
      else if(I.eq.1) then
       MDCutO = 1.d-30
      else if(I.eq.2) then
       MDCutO = 1.d-6
      endIf
      Return
      End
*Deck MkIntF
      Subroutine MkIntF(IOut,IPunch,T,TotWt,QRot,QElec,QNSp,dip,
     $  OutFil,Title)
      Implicit Real*8 (A-H,O-Z)
C INum = number (I2) for differentiating entries with the same molecular 
C weight)
      Character*(*) Title,OutFil
      Character*80 IntFil*80, LblDip*2
      Dimension dip(3)
      Save LblDip,Small
      Data LblDip /'00'/, Small/1.0d-5/
  100 Format(2I3,I1,A1,I1,F10.4,2I5,2F10.1,2I8)
      IEnd = LinEnd(OutFil)-3
      IntFil=OutFil
      Intfil(IEnd:IEnd+3)='.int'
      OPEN(IPunch,FILE=IntFil,STATUS='UNKNOWN')
      Rewind(IPunch)
      QTot = QRot*QElec*QNSp
      Mass = Int(TotWt)
      IFMin = 0
      IFMax = 150
      StLog1 = -8.0D0
      STLog2 = -8.0D0
      MxFrGH = 2000
      IFlag  = 1
      Index  = 5
      INum   = 1
      ITemp  = Int(T)
      Write(IPunch,'(A80)') Title
      Write(IPunch,100) IFlag,Mass,Index,LblDip(1:1),INum,QTot,IFMin,
     $  IFMax,STlog1,STlog2,MxFrGH,ITemp
      do 10 ixyz=1,3
       if(abs(dip(IXYZ)).gt.small) write(IPunch,'(3X,A2,I1,F10.3)')
     $   LblDip,IXYZ,Dip(IXYZ)
   10 Continue
      close(IPunch)
      Return
      End
*Deck MOfI 
      Subroutine MOfI(IOut,Iprint,Align,NCopy,NAtoms,C,AtMass,COM,TMom,
     $  PMom,EigVec)
      Implicit Real*8(A-H,O-Z)
C
C     Compute the principal moments of inertia (units are amu-bohr**2)
C     combined for NCopy sets of coordinates and return them along
C     with their eigenvalues and eigenvectors.  Align is true to flip
C     signs and permute eigenvectors to minimize the change from the
C     initial structure.
C
      Logical Align
      Dimension C(3,NAtoms,*), AtMass(*), COM(3), PMom(3), EigVec(3,3),
     $  T(6), E(9), E2(18), TMom(6), TRef(3,3), IPBest(3), TBest(3,3)
      Dimension IX(3),WA(18),Scr(6),Eig1(3,3),A(3,3)
C
C     Compute the position of the center of mass.
C
      Call CntMas(0,1,NAtoms,AtMass,C,TotWt,COM)
C
C     Translate the molecule and compute the principal moments.
C
      Call AClear(6,T)
      Do 30 IC = 1, NCopy
        Do 20 IAt = 1, NAtoms
          Wt = AtMass(IAt)
          X  = C(1,IAt,IC) - COM(1)
          Y  = C(2,IAt,IC) - COM(2)
          Z  = C(3,IAt,IC) - COM(3)
          T(1) = T(1) + Wt * (Y*Y+Z*Z)
          T(3) = T(3) + Wt * (X*X+Z*Z)
          T(6) = T(6) + Wt * (X*X+Y*Y)
          T(2) = T(2) - Wt * X * Y
          T(4) = T(4) - Wt * X * Z
          T(5) = T(5) - Wt * Y * Z
   20    Continue
   30   Continue
      Call AMove(6,T,TMom)
      PMom(1) = T(6)
      PMom(2) = T(3)
      PMom(3) = T(1)
      Call Aclear(9,EigVec)
      EigVec(1,1) = One
      EigVec(2,2) = One
      EigVec(3,3) = One
      Call HQRII1(IOut,3,1,3,0,T,PMom,3,EigVec,.true.,IErr,IX,WA,SCr,6)
C     Call DiagD(T,EigVec,PMom,3,E,E2,3,.False.)
      If(Align) then
        Call AUnitM(.False.,1,3,3,TRef)
        Call RotChk(IOut,0,0,EigVec,TRef,JPBest,IXBest,IPBest,TBest)
        If(JPBest.ne.0.or.IXBest.ne.0) then
          Call AMove(9,TBest,EigVec)
          Call AMove(3,PMom,E)
          Do 40 I = 1, 3
   40       PMom(I) = E(IPBest(I))
          endIf
        endIf
      Return
      End
*Deck MWeigV
      Subroutine MWeigV(Inv,NAtoms,Cin,RtMass,Cout)
      Implicit Real*8(A-H,O-Z)
C
C     Do/Undo mass-weighting of coordinates or forces
C     Inv = True     Divide by Root-Mass
C     Inv = False    Multiply by Root-Mass
C
      Dimension Cin(3,NAtoms),Cout(3,NAtoms),Rtmass(NAtoms)
      Logical Inv
      Save One
      Data One/1.d0/
C
      Do 100 IAtom = 1, NAtoms
        RtMasI = RtMass(IAtom)
        If(Inv) RtMasI = One/RtMass(IAtom)
        Call AScale(3,RtMasI,Cin(1,IAtom),Cout(1,IAtom))
  100   Continue
      Return
      End
*Deck NumChar
      Subroutine NumChar(Num,Str,ILen)
      Implicit Real*8 (A-H,O-Z)
C
C Form string Str with length ILen from number Num
C
      Character SNum*10,Str*4
      Data SNum/'0123456789'/
      IThous=Num/1000
      Num1=Num-IThous*1000
      IHund=Num1/100
      Num2=Num1-IHund*100
      Iten=Num2/10
      Num3=Num2-ITen*10
      IUn=Num3
      If(IThous.ne.0) then
       Str(1:1)=SNum(IThous+1:IThous+1)
       Str(2:2)=SNum(IHund+1:IHund+1)
       Str(3:3)=SNum(ITen+1:ITen+1)
       Str(4:4)=SNum(IUn+1:IUn+1)
       ILen=4
      ElseIf(IHund.ne.0) then
       Str(1:1)=SNum(IHund+1:IHund+1)
       Str(2:2)=SNum(ITen+1:ITen+1)
       Str(3:3)=SNum(IUn+1:IUn+1)
       ILen=3
      ElseIf(ITen.ne.0) then
       Str(1:1)=SNum(ITen+1:ITen+1)
       Str(2:2)=SNum(IUn+1:IUn+1)
       ILen=2
      Else
       Str(1:1)=SNum(IUn+1:IUn+1)
       ILen=1
      EndIf
      Return
      End
*Deck Numer
      FUNCTION NUMER(NGRP)
C
C     FORM A NUMBER FROM THE HOLLERITH DIGITS IN NGRP.
C
      DIMENSION NGRP(*), N(10)
      Save N
      DATA N/1H0, 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9/
C
      IDIGIT = NGRP(2)
      JDIGIT = NGRP(3)
      NUMER = 0
      N1 = -1
      N2 = -1
      DO 20 I=1,10
         IF (IDIGIT .EQ. N(I)) N1 = I - 1
         IF (JDIGIT .EQ. N(I)) N2 = I - 1
   20    CONTINUE
      NUMER = N1
      IF (N2 .NE. -1) NUMER = 10*N1 + N2
      IF (NUMER .LT. 0) NUMER = 0
C
      RETURN
      END
*Deck OutMat
      Subroutine OutMat(IOut,Key,X,M,N,MM,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Print matrix X.  ISmall is >0 to print all elements, <=0 to
C     print elements greater than 10**(Key-6)
C
      Dimension X(M,N)
C
      Call OutMtC(IOut,1,.False.,Key,1,1,1,X,1,1,M,N,MM,NN)
      Return
      End
*Deck OutMtC
      Subroutine OutMtC(IOut,IOpt,NSpinr,Key,NRI,IRI,NSpin,X,MSt,NSt,M,
     $  N,MM,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Print matrix X.  ISmall is >0 to print all elements, <=0 to
C     print elements greater than 10**(Key-6).  NRI is the number
C     of components (1 for real, 2 for complex), IRI which component
C     to print (-1=phase, 0=modulus, 1=real, 2=imaginary), NSpin is the
C     number of spinor terms (-1=just large, 1=just one spin, -2=l and
C     s, 2=a and b).  M,N,MM,NN refer to the number of elements allocated
C     and used (i.e., to print all of a matrix MSt:MEnd M=MM=MEnd-MSt+1).
C
C     IOpt .. Output format:
C             1 ... D14.6
C             2 ... F14.3
C             3 ... F20.8
C
      Parameter (NumCol=5)
      Logical NSpinr
      Character*1 LabX(NumCol)
      Dimension X(NRI,M,N), Y(NumCol), IRX(NumCol)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(5X,8(7X,I6,A1))
 1010 Format(I7,A,8D14.6)
 1020 Format(I7,A,8F14.3)
 1030 Format(I7,A,5F20.8)
 1040 Format(1X,8(13X,I6,A1))
C
      If(IOpt.lt.1.or.IOpt.gt.3) then
       write(IOut,'(''Illegal IOpt in OutMtC.'')')
       Stop
      EndIf 
      If(NRI.eq.1) then
        IRIU = 1
      else
        IRIU = IRI
        If(IRI.lt.-1.or.IRI.gt.2) then
         write(Iout,'('' Illegal IRI in OutMtC.'')')
         stop
        endif 
      endIf
      Call PrtThr(Key,Small)
      MS = MM*IAbs(NSpin)
      If(NSpinr) then
        NS = NN*IAbs(NSpin)
        NSpinX = NSpin
      else
        NS = NN
        NSpinX = 1
        endIf
      Do 100 ILower = 1, NS, NumCol
        IUpper = Min(ILower+NumCol-1,NS)
        Num = IUpper - ILower + 1
        Call LabSpn(NSpinX,ILower-2+NSt,Num,IRX,LabX)
        If(IOpt.le.2) then
          Write(IOut,1000) (IRX(I),LabX(I),I=1,Num)
        else
          Write(IOut,1040) (IRX(I),LabX(I),I=1,Num)
          endIf
        Do 90 I = 1, MS
          Call LabSpn(NSpin,I-2+MSt,1,IRX,LabX)
          If(IRI.eq.-1) then
            Do 10 J = ILower, IUpper
   10         Y(J-ILower+1) = ZATan2(X(2,I,J),X(1,I,J))
          else if(IRI.eq.0) then
            Do 20 J = ILower, IUpper
   20         Y(J-ILower+1) = Sqrt(X(1,I,J)**2+X(2,I,J)**2)
          else
            Do 30 J = ILower, IUpper
   30         Y(J-ILower+1) = X(IRI,I,J)
            endIf
          Do 40 J = 1, Num
            If(Key.ne.1.and.Abs(Y(J)).lt.Small) Y(J) = Zero
   40       Continue
          If(IOpt.eq.1) then
            Write(IOut,1010) IRX(1), LabX(1), (Y(J),J=1,Num)
          else if(IOpt.eq.2) then
            Write(IOut,1020) IRX(1), LabX(1), (Y(J),J=1,Num)
          else
            Write(IOut,1030) IRX(1), LabX(1), (Y(J),J=1,Num)
            endIf
   90     Continue
  100   Continue
      Return
      End
*Deck OutMtS
      Subroutine OutMtS(IOut,String,Num,ISmall,X,M,N,MM,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Print out a heading followed by a matrix.  If Num is
C     positive, it is also printed.
C
      Character*(*) String
      Dimension X(*)
C
      Call HedPrt(IOut,0,String,Num)
      Call OutMtC(IOut,1,.False.,ISmall,1,1,1,X,1,1,M,N,MM,NN)
      Return
      End
*Deck PhyFil
      Subroutine PhyFil(IOpt,LenPhy,PhyCon)
      Implicit None
C
C     Routine to fill the array PhyCon with values of useful physical
C     constants.  IOpt:
C
C     0 ... Default (2010).
C
C     1979 ... values used through Gaussian 86, taken from
C              Pure and applied chemistry, 51, 1 (1979).
C
C     1986 ... values used through Gaussian 98, taken from
C              E. R. Cohen, B. N. Taylor, "The 1986 Adjustment of the
C              Fundamental Physical Constants," report of the CODATA Task
C              group on Fundamental Constants, CODATA bulletin 63, Pergamon,
C              Elmsford, NY (1986).
C
C     1998 ... 1998 values, from P. J. Mohr and B. N. Taylor, Physics Today,
C              August 2000, Page BG6 which summarized the 1988 NIST publications.
C
C     2006 ... 2006 CoData values.
C
C     2010 ... 2010 CoData values.
C
      Integer IOpt,LenPhy,ICon,NConst,ToAng,ToKG,ToE,Planck,Avog,JPCal,
     $  MPerB,Hartre,SLight,Boltz,FineSC,EMKG,VolMol,EMM,PRM,LenDat,
     $  LenClr,I,J,GFree,Next
      Parameter (NConst=5,ToAng=1,ToKG=2,ToE=3,Planck=4,Avog=5,JPCal=6,
     $  MPerB=7,Hartre=8,SLight=9,Boltz=10,FineSC=11,EMKG=12,VolMol=13,
     $  EMM=14,PRM=15,GFree=16,Next=17,LenDat=30,
     $  LenClr=((LenDat-Next+1)*NConst))
      Real*8 One,Ten,Tp4,PhyCon(LenPhy), PhyDat(LenDat,NConst)
      Save PhyDat, One, Ten, Tp4
      Data One/1.0d0/, Ten/1.0d1/, Tp4/1.0d4/
C
C     Angstroms per Bohr
      Data (PhyDat(ToAng,I),I=1,NConst)/0.52917706D+00,0.529177249D+00,
     $  0.5291772083d0,0.52917720859d0,0.52917721092d0/
C
C     Kilograms per atomic mass unit
      Data (PhyDat(ToKG,I),I=1,NConst)/1.6605655D-27,1.6605402D-27,
     $  1.66053873D-27,1.660538782d-27,1.660538921d-27/
C
C     Coulombs per electron -- multiplied by speed of light/10 to produce
C     Electrostatic units (ESU) per electron.  Old values have lots of digits
C     to reproduce values previous in data statements.
      Data (PhyDat(ToE,I),I=1,NConst)/1.6021890717477622D-19,
     $  1.6021890717477622D-19,1.602176462d-19,1.602176487d-19,
     $  1.602176565d-19/
C
C     Planck constant, Joule-Seconds
      Data (PhyDat(Planck,I),I=1,NConst)/6.626176D-34,6.6260755D-34,
     $  6.62606876D-34,6.62606896d-34,6.62606957d-34/
C
C     Avogadro constant
      Data (PhyDat(Avog,I),I=1,NConst)/6.022045D+23,6.0221367D+23,
     $  6.02214199D+23,6.02214179d+23,6.02214129d+23/
C
C     Joules per calorie
      Data (PhyDat(JPCal,I),I=1,NConst)/NConst*4.184D+00/
C
C     Meters per Bohr, computed from other constants at run-time
      Data (PhyDat(MPerB,I),I=1,NConst)/NConst*0.0d0/
C
C     Joules per Hartre
      Data (PhyDat(Hartre,I),I=1,NConst)/4.359814D-18,4.3597482D-18,
     $  4.35974381D-18,4.35974394d-18,4.35974434d-18/
C
C     Speed of light, cm sec(-1)
      Data (PhyDat(SLight,I),I=1,NConst)/NConst*2.99792458d+10/
C
C     Boltzman constant, Joules per Kelvin
      Data (PhyDat(Boltz,I),I=1,NConst)/1.380662D-23,1.380658D-23,
     $  1.3806503d-23,1.3806504d-23,1.3806488d-23/
C
C     Inverse Fine structure constant, inverted at run-time.
      Data (PhyDat(FineSC,I),I=1,NConst)/137.03602D+00,137.0359895D+00,
     $  137.03599976d0,137.035999679d0,137.035999074d0/
C
C     Electron mass in KG, computed at run-time from other constants.
      Data (PhyDat(EMKG,I),I=1,NConst)/NConst*0.0d0/
C
C     Molar volume of ideal gas in m**3 at 273.15 K
      Data (PhyDat(VolMol,I),I=1,NConst)/22.41383d-3,22.41410d-3,
     $  22.413996d-3,22.413996d-3,22.4139679d-3/
C
C     Electron Magnetic Moment (J T-1) (sign flipped).  Was not present
C     in G86 or before, so 1988 value used.
C
      Data (PhyDat(EMM,I),I=1,NConst)/9.2847701D-24,9.2847701D-24,
     $  928.476362d-26,928.476377d-26,928.476430d-26/
C
C     Proton Rest Mass (Kg), not present in G86 or before so 1988 value used.
      Data (PhyDat(PRM,I),I=1,NConst)/1.672623D-27,1.672623D-27,
     $  1.67262158d-27,1.672621637d-27,1.672621777d-27/
C
C     Free-electron g-factor.
      Data (PhyDat(GFree,I),I=1,NConst)/3*2.002319304386d0,
     $  2.0023193043622d0,2.00231930436153d0/
      Data ((PhyDat(J,I),J=Next,LenDat),I=1,NConst)/LenClr*0.0d0/
C
      If(IOpt.eq.0) then
        ICon = 5
      else if(IOpt.eq.1979) then
        ICon = 1
      else if(IOpt.eq.1986) then
        ICon = 2
      else if(IOpt.eq.1998) then
        ICon = 3
      else if(IOpt.eq.2006) then
        ICon = 4
      else if(IOpt.eq.2010) then
        ICon = 5
      else
       write(*,*) 'Illegal IOpt in PhyFil.'
        Stop
      endIf
      If(LenPhy.lt.LenDat) then
        write(*,*) 'Output array is too short in PhyFil.'
        Stop
      EndIf
      Call AMove(LenDat,PhyDat(1,ICon),PhyCon)
      Call AClear(LenPhy-LenDat,PhyCon(LenDat+1))
C
C     Set some values derived from other constants
      PhyCon(ToE) = PhyCon(ToE)*PhyCon(SLight)/Ten
      PhyCon(MPerB) = PhyCon(ToAng) / Ten**10
      PhyCon(FineSC) = One / PhyCon(FineSC)
      PhyCon(EMKG) = PhyCon(Hartre)*Tp4 /
     $  (PhyCon(SLight)*PhyCon(FineSC))**2
      Return
      End
*Deck Place
      Subroutine Place(CutOff,ZH,CRef,Bl,Alpha,C,NonLin)
      Implicit Real*8(A-H,O-Z)
C
C     Place an atom when all previous atoms are collinear.
C     ZH is a unit vector along the collinear axis, CRef
C     is the atom referenced to the present atom, Bl and
C     Alpha are the bond lengths, and C is the output
C     vector.  The atom is placed in the plane defined
C     by ZH and the X-axis if possible; otherwise the Z-axis
C     is used.
C
      Logical NonLin
      Dimension ZH(3), CRef(3), C(3), CLoc(3), XH(3)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      CLoc(1) = Bl*Sin(Alpha)
      CLoc(2) = Zero
      CLoc(3) = -Bl*Cos(Alpha)
      NonLin = Abs(CLoc(1)).gt.CutOff
      XH(1) = One - ZH(1)*ZH(1)
      XH(2) = -ZH(1)*ZH(2)
      XH(3) = -ZH(1)*ZH(3)
      R = SProd(3,XH,XH)
      If(R.le.CutOff) then
        XH(1) = -ZH(3)*ZH(1)
        XH(2) = -ZH(3)*ZH(2)
        XH(3) = One - ZH(3)*ZH(3)
        R = SProd(3,XH,XH)
        endIf
      Do 10 I = 1, 3
   10   C(I) = CLoc(1)*XH(I) + CLoc(3)*ZH(I) + CRef(I)
      Return
      End
*Deck QuaIRC
      Subroutine QuaIrc(IOut,IPrint,NAtoms,AtMass,CReac,CProd,RotMat,
     $  NIter,Fail)
      Implicit Real*8 (A-H,O-Z)
C
C QUAternion-based procedure for IRC
C     Rotates the product coordinates until no total angular momentum

C     is present along the linear path connecting reactants with
C     products:
C     0 = Sum{ Atmass(I) * [CReac(I)/\CProd(I)]}
C
C Input:
C     AtMass : (NAtoms) Atomic masses
C     CReac  : (3,NAtoms) Atomic coordinates of the reactant(s)
C     CProd  : (3,NAtoms) Atomic coordinates of the product(s)
C     NIter  : Max. number of iterations
C
C Output:
C     RotMat : (3,3) Rotation matrix
C     CProd  : (3,NAtoms) Rotated atomic coord. of the product(s)
C     Fail   : Reports if superposition has failed or not
C
C     Dimension
      Integer NAtoms
C     Input
      Integer IOut, IPrint, NIter
      Real*8 CReac(3,*), AtMass(*)
C     Output
      Real*8 CProd(3,*), RotMat(3,3)
      Logical Fail
C     Local
      Integer IOut2, NCyc
      Real*8 A(3,3), CRCP(3,3), One
      Logical OK
      Save One
      Data One/1.0d0/
 7000 Format(' Superposition has been achieved after ',I3,' iterations')
 7001 Format(' Superposition was unsuccessful after ',I3,' iterations')
C
      If(IPrint.gt.0) then
        IOut2 = IOut
      else
        IOut2 = 0
        endIf
      Fail = .False.
      OK = .False.
      Call AClear(9,RotMat)
      RotMat(1,1) = One
      RotMat(2,2) = One
      RotMat(3,3) = One
      NCyc = 0
  100 If(.not.OK) then
        Call AClear(9,CRCP)
        Call QuaRot(NAtoms,AtMass,CReac,CProd,CRCP)
C       Call RotF1(NAtoms,CRCP,CProd)
        Call Aclear(9,A)
C       Call MatMpy(3,3,3,CRCP,RotMat,A)
        Call AMove(9,A,RotMat)
        NCyc = NCyc + 1
        If(NCyc.lt.NIter) then
          Call ChkAnm(IOut2,NAtoms,AtMass,CReac,CProd,OK)
          Goto 100
          endIf
        endIf
      If(IPrint.ge.1) then
        If(OK) then
          Write(IOut,7000) NCyc
        else
          Write(IOut,7001) NCyc
          endIf
        endIf
C     If(IPrint.ge.1.and.OK) Call MatOut(RotMat,3,3,3,3)
      Fail = .not.OK
      Return
      End
*Deck QuaRot
      Subroutine QuaRot(NAtoms,AtMass,C0,C,RotMat)
      Implicit Real*8 (A-H,O-Z)
C
C QUAternion-based ROTation matrix
C     Find the best superposition between 2 molecular structures using
C     quaternions.
C     See G.R. Kneller, Mol. Sim. 7, 113-119 (1991)
C               "       J. Chim. Phys. 88,2709-2715 (1991)
C Note:
C     QuaMat is the matrix M    of the above reference
C     RotMat is the matrix D(q) of the above reference
C
C Input:
C     AtMass : (NAtoms) Atomic masses
C     C0     : (3,NAtoms) Atomic coord. of the reference structure
C     C      : (3,NAtoms) Atomic coord. of the structure to superpose
C
C Output:
C     RotMat : (3,3) Rotation matrix
C
C     Dimension
      Integer NAtoms
C     Input
      Real*8 C(3,*), C0(3,*), AtMass(*)
C     Output
      Real*8 RotMat(3,3)
C     Local
      Integer ia
      Real*8 E(4),E2(8), EigVal(4), QuaMat(10), Quat(16), Diag0, Diag1,
     $  Diag2, Diag3, Q0Q0, Q0Q1, Q0Q2, Q0Q3, Q1, Q10, Q1Q1, Q1Q2, Q1Q3,
     $  Q2, Q2Q2, Q2Q3, Q3, Q3Q3, Q4, Q5, Q6, Q7, Q8, Q9, Two, X, X0,
     $  XY0, XZ0, Y, Y0, YX0, YZ0, Z, Z0, ZX0, ZY0
      Save Two
      Data Two/2.0D0/
C
      Call AClear(10,QuaMat)
      Do 100 ia = 1, NAtoms
        X   = C(1,ia)
        Y   = C(2,ia)
        Z   = C(3,ia)
        X0  = C0(1,ia)
        Y0  = C0(2,ia)
        Z0  = C0(3,ia)
        XY0 = X*Y0
        XZ0 = X*Z0
        YX0 = Y*X0
        YZ0 = Y*Z0
        ZX0 = Z*X0
        ZY0 = Z*Y0
        Diag0 = X*X + Y*Y + Z*Z + X0*X0 + Y0*Y0 +Z0*Z0
        Diag1 = Two*X*X0
        Diag2 = Two*Y*Y0
        Diag3 = Two*Z*Z0
        Q1 = Diag0 - Diag1 - Diag2 - Diag3
        Q3 = Diag0 - Diag1 + Diag2 + Diag3
        Q6 = Diag0 + Diag1 - Diag2 + Diag3
        Q10= Diag0 + Diag1 + Diag2 - Diag3
        Q2 = YZ0 - ZY0
        Q4 = ZX0 - XZ0
        Q7 = XY0 - YX0
        Q5 = -(XY0 + YX0)
        Q8 = -(XZ0 + ZX0)
        Q9 = -(YZ0 + ZY0)
        QuaMat(1)  = QuaMat(1)  + Q1 * AtMass(ia)
        QuaMat(3)  = QuaMat(3)  + Q3 * AtMass(ia)
        QuaMat(6)  = QuaMat(6)  + Q6 * AtMass(ia)
        QuaMat(10) = QuaMat(10) + Q10* AtMass(ia)
        QuaMat(2)  = QuaMat(2)  + Q2 * AtMass(ia) * Two
        QuaMat(4)  = QuaMat(4)  + Q4 * AtMass(ia) * Two
        QuaMat(7)  = QuaMat(7)  + Q7 * AtMass(ia) * Two
        QuaMat(5)  = QuaMat(5)  + Q5 * AtMass(ia) * Two
        QuaMat(8)  = QuaMat(8)  + Q8 * AtMass(ia) * Two
        QuaMat(8)  = QuaMat(8)  + Q8 * AtMass(ia) * Two
        QuaMat(9)  = QuaMat(9)  + Q9 * AtMass(ia) * Two
  100   Continue
C     +---------------+
C      DIAGONALIZATION
C     +---------------+
C     Call DiagD(QuaMat,Quat,EigVal,4,E,E2,4,.False.)
C     +------------------------+
C      BUILDING ROTATION MATRIX
C     +------------------------+
      Q0Q0 = Quat(1) * Quat(1)
      Q1Q1 = Quat(2) * Quat(2)
      Q2Q2 = Quat(3) * Quat(3)
      Q3Q3 = Quat(4) * Quat(4)
      RotMat(1,1) = Q0Q0 + Q1Q1 - Q2Q2 - Q3Q3
      RotMat(2,2) = Q0Q0 - Q1Q1 + Q2Q2 - Q3Q3
      RotMat(3,3) = Q0Q0 - Q1Q1 - Q2Q2 + Q3Q3
      Q0Q1 = Quat(1) * Quat(2)
      Q0Q2 = Quat(1) * Quat(3)
      Q0Q3 = Quat(1) * Quat(4)
      Q1Q2 = Quat(2) * Quat(3)
      Q1Q3 = Quat(2) * Quat(4)
      Q2Q3 = Quat(3) * Quat(4)
      RotMat(1,2) = Two * (Q1Q2 - Q0Q3)
      RotMat(2,1) = Two * (Q1Q2 + Q0Q3)
      RotMat(1,3) = Two * (Q1Q3 + Q0Q2)
      RotMat(3,1) = Two * (Q1Q3 - Q0Q2)
      RotMat(2,3) = Two * (Q2Q3 - Q0Q1)
      RotMat(3,2) = Two * (Q2Q3 + Q0Q1)
C     Call MatOut(RotMat,3,3,3,3)
      Return
      End
*Deck RCovA
      Function RCovA(IA,IB)
      Implicit Real*8(A-H,O-Z)
C
C     Return a covalent radius in Angstroms.
C
      Common /PhyCon/ PhyCon(30)
C
      RCovA = RCov(IA,IB) * PhyCon(1)
      Return
      End
*Deck RCov
      Function RCov(IA,IB)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns an estimated covalent bond distance between
C     atoms of atomic numbers IA and IB.  Since ghost atoms are assumed
C     to have 0 radius, setting IB to 0 returns the covalent radius of IA.
C
      Parameter (MaxAn=85)
      Dimension Rad(0:MaxAn)
      Save Rad
      Data Rad/0.0d0,
     $  0.643d0,0.643d0,2.457d0,1.909d0,1.587d0,1.436d0,1.209d0,
     $  1.096d0,1.020d0,0.945d0,2.986d0,2.646d0,2.400d0,2.192d0,
     $  2.060d0,1.890d0,1.795d0,1.701d0,3.836d0,3.288d0,2.721d0,
     $  2.494d0,2.305d0,2.230d0,2.211d0,2.211d0,2.192d0,2.173d0,
     $  2.211d0,2.362d0,2.381d0,2.305d0,2.268d0,2.192d0,2.154d0,
     $  2.116d0,4.082d0,3.609d0,3.061d0,2.740d0,2.532d0,2.457d0,
     $  2.400d0,2.362d0,2.362d0,2.419d0,2.532d0,2.797d0,2.721d0,
     $  2.665d0,2.646d0,2.570d0,2.513d0,2.476d0,4.441d0,3.742d0,
     $  3.194d0,3.118d0,3.118d0,3.099d0,3.080d0,3.061d0,3.496d0,
     $  3.042d0,3.005d0,3.005d0,2.986d0,2.967d0,2.948d0,2.948d0,
     $  2.948d0,2.721d0,2.532d0,2.457d0,2.419d0,2.381d0,2.400d0,
     $  2.457d0,2.532d0,2.816d0,2.797d0,2.778d0,2.759d0,2.759d0,
     $  2.740d0/
C
      RCov = Rad(Min(Max(IA,0),MaxAn)) + Rad(Min(Max(IB,0),MaxAn))
      Return
      End
*Deck RCovCT
      Function RCovCT(IA,JA)
      Implicit Real*8(A-H,O-Z)
C
C     Return the bond distance in Angstroms between elements IA and JA
C     using radii from M. Mantina, R. Valero, C. J. Cramer, and
C     D. G. Truhlar, listed in the CRC Handbook.
C
      Parameter (MaxAn=118)
      Dimension Rad(0:MaxAn)
      Save Rad
      Data Rad/0.0d0,0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,
     $  0.64D0,0.60D0,0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,
     $  1.00D0,1.01D0,2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,
     $  1.24D0,1.18D0,1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,
     $  1.17D0,1.16D0,2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,
     $  1.36D0,1.34D0,1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,
     $  1.36D0,1.36D0,2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,
     $  1.85D0,1.83D0,1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,
     $  1.74D0,1.64D0,1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,
     $  1.32D0,1.44D0,1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,
     $  2.01D0,1.90D0,1.84D0,1.83D0,1.80D0,1.80D0,1.73D0,1.68D0,1.68D0,
     $  1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0,1.57D0,1.49D0,1.43D0,
     $  1.41D0,1.34D0,1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,
     $  1.75D0,1.65D0,1.57D0/
C
      IA1 = Min(Max(IA,0),MaxAn)
      JA1 = Min(Max(JA,0),MaxAn)
      RCovCT = Rad(IA1) + Rad(JA1)
      Return
      End
*Deck RCov97
      Function RCov97(IA,IB)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns an estimated covalent bond distance (in Ang)
C     between atoms of atomic numbers IA and IB. Setting IB to 0 returns
C     the covalent radius of IA. Parameters for atoms heavier than At are
C     taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
C     W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).
C
      Parameter (MaxI=104)
      Dimension Radius(0:MaxI)
      Save Radius
      Data Radius/0.0d0,
     $  0.354D0,0.849D0,
     $  1.336D0,1.010D0,0.838D0,0.757D0,0.700D0,0.658D0,0.668D0,0.920D0,
     $  1.539D0,1.421D0,1.244D0,1.117D0,1.101D0,1.064D0,1.044D0,1.032D0,
     $  1.953D0,1.761D0,
     $  1.513D0,1.412D0,1.402D0,1.345D0,1.382D0,
     $  1.270D0,1.241D0,1.164D0,1.302D0,1.193D0,
     $                  1.260D0,1.197D0,1.211D0,1.190D0,1.192D0,1.147D0,
     $  2.260D0,2.052D0,
     $  1.698D0,1.564D0,1.473D0,1.467D0,1.322D0,
     $  1.478D0,1.332D0,1.338D0,1.386D0,1.403D0,
     $                  1.459D0,1.398D0,1.407D0,1.386D0,1.382D0,1.267D0,
     $  2.570D0,2.277D0,
     $  1.943D0,1.841D0,1.823D0,1.816D0,1.801D0,1.780D0,1.771D0,
     $  1.735D0,1.732D0,1.710D0,1.696D0,1.673D0,1.660D0,1.637D0,
     $  1.671D0,1.611D0,1.511D0,1.392D0,1.372D0,
     $  1.372D0,1.371D0,1.364D0,1.262D0,1.340D0,
     $                  1.518D0,1.459D0,1.512D0,1.500D0,1.545D0,1.420D0,
     $  2.880D0,2.512D0,1.983D0,1.721D0,1.711D0,1.684D0,1.666D0,
     $  1.657D0,1.660D0,1.801D0,1.761D0,1.750D0,1.724D0,1.712D0,
     $  1.689D0,1.679D0,1.698D0,
     $  1.850D0/
C
      RCov97 = Radius(Min(Max(IA,0),MaxI)) + Radius(Min(Max(IB,0),MaxI))
      Return
      End
*Deck RCv123
      Function RCv123(IA,JA,IBnd)
      Implicit Real*8(A-H,O-Z)
C
C     Return the bond distance in Angstroms between elements IA and JA
C     (IBnd is 1,2,3 for single, double or triple bond) using radii 
C     from P. Pyykko JPCA 119, 2326-2337 (2015).
C
      Parameter (MaxAn=118)
      Dimension Rad(3,MaxAn)
      Save Rad
      Data Rad/0.32D0,2*0.D0,0.46D0,2*0.d0,
     $  1.33D0,1.24D0,1.24d0,1.02d0,0.90d0,0.85d0,0.85d0,0.78d0,0.73d0,
     $  0.75d0,0.67d0,0.60d0,0.71d0,0.60d0,0.54d0,0.63d0,0.57d0,0.53d0,
     $  0.64d0,0.59d0,0.53d0,0.67d0,0.96d0,0.00d0,
     $  1.55d0,1.60d0,1.60d0,1.39d0,1.32d0,1.27d0,1.26d0,1.13d0,1.11d0,
     $  1.16d0,1.07d0,1.02d0,1.11d0,1.02d0,0.94d0,1.03d0,0.94d0,0.95d0,
     $  0.99d0,0.95d0,0.93d0,0.96d0,1.07d0,0.96d0,
     $  1.96d0,1.93d0,1.93d0,1.71d0,1.47d0,1.33d0,1.48d0,1.16d0,1.14d0,
     $  1.36d0,1.17d0,1.08d0,1.34d0,1.12d0,1.06d0,1.22d0,1.11d0,1.03d0,
     $  1.19d0,1.05d0,1.03d0,1.16d0,1.09d0,1.02d0,1.11d0,1.03d0,0.96d0,
     $  1.10d0,1.01d0,1.01d0,1.12d0,1.15d0,1.20d0,1.18d0,1.20d0,0.00d0,
     $  1.24d0,1.17d0,1.21d0,1.21d0,1.11d0,1.14d0,1.21d0,1.14d0,1.06d0,
     $  1.16d0,1.07d0,1.07d0,1.14d0,1.09d0,1.10d0,1.17d0,1.21d0,1.08d0,
     $  2.10d0,2.02d0,2.02d0,1.85d0,1.57d0,1.39d0,1.63d0,1.30d0,1.24d0,
     $  1.54d0,1.27d0,1.21d0,1.47d0,1.25d0,1.16d0,1.38d0,1.21d0,1.13d0,
     $  1.28d0,1.20d0,1.10d0,1.25d0,1.14d0,1.03d0,1.25d0,1.10d0,1.06d0,
     $  1.20d0,1.17d0,1.12d0,1.28d0,1.39d0,1.37d0,1.36d0,1.44d0,0.00d0,
     $  1.42d0,1.36d0,1.46d0,1.40d0,1.30d0,1.32d0,1.40d0,1.33d0,1.27d0,
     $  1.36d0,1.28d0,1.21d0,1.33d0,1.29d0,1.25d0,1.31d0,1.35d0,1.22d0,
     $  2.32d0,2.09d0,0.00d0,1.96d0,1.61d0,1.49d0,1.80d0,1.39d0,1.39d0,
     $  1.63d0,1.37d0,1.31d0,1.76d0,1.38d0,1.28d0,1.74d0,1.37d0,1.37d0,
     $  1.73d0,1.35d0,1.35d0,1.72d0,1.34d0,1.34d0,1.68d0,1.34d0,1.34d0,
     $  1.69d0,1.35d0,1.32d0,1.68d0,1.35d0,1.35d0,1.67d0,1.33d0,1.33d0,
     $  1.66d0,1.33d0,1.33d0,1.65d0,1.33d0,1.33d0,1.64d0,1.31d0,1.31d0,
     $  1.70d0,1.29d0,1.29d0,1.62d0,1.31d0,1.31d0,1.52d0,1.28d0,1.22d0,
     $  1.46d0,1.26d0,1.19d0,1.37d0,1.20d0,1.15d0,1.31d0,1.19d0,1.10d0,
     $  1.29d0,1.16d0,1.09d0,1.22d0,1.15d0,1.07d0,1.23d0,1.12d0,1.10d0,
     $  1.24d0,1.21d0,1.23d0,1.33d0,1.42d0,1.42d0,1.44d0,1.42d0,1.50d0,
     $  1.44d0,1.35d0,1.37d0,1.51d0,1.41d0,1.35d0,1.45d0,1.35d0,1.29d0,
     $  1.47d0,1.38d0,1.38d0,1.42d0,1.45d0,1.33d0,
     $  2.23d0,2.18d0,2.18d0,2.01d0,1.73d0,1.59d0,1.86d0,1.53d0,1.40d0,
     $  1.75d0,1.43d0,1.36d0,1.69d0,1.38d0,1.29d0,1.70d0,1.34d0,1.18d0,
     $  1.71d0,1.36d0,1.16d0,1.72d0,1.35d0,1.35d0,1.66d0,1.35d0,1.35d0,
     $  1.66d0,1.36d0,1.36d0,1.68d0,1.39d0,1.39d0,1.68d0,1.40d0,1.40d0,
     $  1.65d0,1.40d0,1.40d0,1.67d0,1.42d0,1.42d0,1.73d0,1.39d0,1.39d0,
     $  1.76d0,1.40d0,1.40d0,1.61d0,1.41d0,1.41d0,1.57d0,1.40d0,1.31d0,
     $  1.49d0,1.36d0,1.26d0,1.43d0,1.28d0,1.21d0,1.41d0,1.28d0,1.19d0,
     $  1.34d0,1.25d0,1.18d0,1.29d0,1.25d0,1.13d0,1.28d0,1.16d0,1.12d0,
     $  1.21d0,1.16d0,1.18d0,1.22d0,1.37d0,1.30d0,1.36d0,1.36d0,1.36d0,
     $  1.43d0,1.43d0,1.43d0,1.62d0,1.62d0,1.62d0,1.75d0,1.75d0,1.75d0,
     $  1.65d0,1.65d0,1.65d0,1.57d0,1.57d0,1.57d0/
C
      IA1 = Min(Max(IA,0),MaxAn)
      JA1 = Min(Max(JA,0),MaxAn)
      RCv123 = Rad(IBnd,IA1) + Rad(IBnd,JA1)
      Return
      End
*Deck RdData
      Subroutine RdData(In,Search,Error)
      Character*80 CLine 
      Character*4 Search,Found
      Logical Error
      Error=.True.
      Rewind(In)
   10 Read(In,'(A80)',end=20) CLine 
      If(CLine(1:1).ne.'#') goto 10
      Call LinUpc(CLine(2:5),Found(1:4))
      If(Found.ne.Search) goto 10
      Error=.False. 
   20 Continue
      Return
      End
*Deck RigHnd
      Subroutine RigHnd(RotVec) 
      Implicit Real*8(A-H,O-Z)
C     
C     Make sure RotVec is right hand frame
C     
      Dimension RotVec(3,3)
C     
      RotVec(1,3) = RotVec(2,1)*RotVec(3,2) - RotVec(3,1)*RotVec(2,2)
      RotVec(2,3) = RotVec(3,1)*RotVec(1,2) - RotVec(1,1)*RotVec(3,2)
      RotVec(3,3) = RotVec(1,1)*RotVec(2,2) - RotVec(2,1)*RotVec(1,2)
      Return
      End
*Deck RotChk
      Subroutine RotChk(IOut,IPrint,IOpt,T,TOld,JPBest,IXBest,IPBest,
     $  TBest)
      Implicit Real*8(A-H,O-Z)
C
C     Compare two rotation matrices and possibly apply 180 degree
C     rotations to T to better match TOld.
C     IOpt = 0 ... All permutations and sign flips considered.
C            1 ... Only rotations, no inversions.
C            2 ... Only sign flips of axes.
C
      Parameter (MinPrt=1,MinPrM=2)
      Character AxLab(3)*1, PxLab(3)*2, PxBest(3)*2
      Logical DoPrnt, DoPrMt, Pos0, PosT
      Dimension T(3,3), TOld(3,3), Tmp(3,3), TBest(3,3), S(3),
     $  IPermL(3,6), IPerm(3), IPBest(3)
      Save One, Small, IPerm, AxLab, PxLab
      Data One/1.0d0/, Small/1.d-6/, IPermL/1,2,3,2,3,1,3,1,2,
     $  2,1,3,1,3,2,3,2,1/, AxLab/'X','Y','Z'/
 1000 Format(' RotChk:  JP=',I1,' IX=',I2,1X,3A2,' Diff=',1PD9.2)
 1010 Format(' T for ',3A2,':')
 1020 Format(' RotChk:  IX=',I1,' Diff=',1PD9.2)
C
      DoPrnt = IPrint.ge.MinPrt
      DoPrMt = IPrint.ge.MinPrM
      If(DoPrMt) then
        Call OutMtS(IOut,'TOld',0,0,TOld,3,3,3,3)
        Call OutMtS(IOut,'T',0,0,T,3,3,3,3)
        endIf
      BestND = UDiff(3,T,TOld)
      Call AMove(9,T,TBest)
      If(DoPrnt) Write(IOut,1000) 0, 0, AxLab, BestND
      JPBest = 1
      IXBest = 0
      Call ISeq(3,0,1,IPBest)
      Do 5 I = 1, 3
    5   PxBest(I) = ' '//AxLab(1)
      Det0 = T(1,1)*(T(2,2)*T(3,3)-T(2,3)*T(3,2))
     $  -T(1,2)*(T(2,1)*T(3,3)-T(2,3)*T(3,1))
     $  +T(1,3)*(T(2,1)*T(3,2)-T(2,2)*T(3,1))
      Pos0 = Det0.ge.(-Small)
      If(IOpt.eq.2) then
        NPUse = 1
        LimIX = 0
      else
        NPUse = 6
        LimIX = 3
        endIf
      Do 40 JP = 1, NPUse
        Call IMove(3,IPermL(1,JP),IPerm)
        Do 30 IX = -3, LimIX
          If(IX.lt.0) then
C           Follow same order as old algorithm.
            PxLab(1) = '-'//AxLab(IPerm(1))
            PxLab(2) = '-'//AxLab(IPerm(2))
            PxLab(3) = '-'//AxLab(IPerm(3))
            Call ASet(3,-One,S)
            IX1 = 4 + IX
            PxLab(IX1) = ' '//AxLab(IPerm(IX1))
            S(IX1) = One
          else
            PxLab(1) = ' '//AxLab(IPerm(1))
            PxLab(2) = ' '//AxLab(IPerm(2))
            PxLab(3) = ' '//AxLab(IPerm(3))
            Call ASet(3,One,S)
            If(IX.gt.0) then
              PxLab(IX) = '-'//AxLab(IPerm(IX))
              S(IX) = -One
              endIf
            endIf
          Do 20 I = 1, 3
            Do 10 J = 1, 3
   10         Tmp(J,I) = S(J)*T(J,IPerm(I))
   20       Continue
          If(IOpt.eq.1) then
            DetT = Tmp(1,1)*(Tmp(2,2)*Tmp(3,3)-Tmp(2,3)*Tmp(3,2))
     $        -Tmp(1,2)*(Tmp(2,1)*Tmp(3,3)-Tmp(2,3)*Tmp(3,1))
     $        +Tmp(1,3)*(Tmp(2,1)*Tmp(3,2)-Tmp(2,2)*Tmp(3,1))
            PosT = DetT.ge.(-Small)
            If(PosT.neqv.Pos0) goto 30
            endIf
          If(DoPrMt) then
            Write(IOut,1010) (PxLab(J),J=1,3)
            Call OutMat(IOut,0,Tmp,3,3,3,3)
            endIf
          Diff = UDiff(3,Tmp,TOld)
          If(DoPrnt) Write(IOut,1000) JP, IX, PxLab, Diff
          If((Diff+Small).lt.BestND) then
            IXBest = IX
            JPBest = JP
            Call IMove(3,IPerm,IPBest)
            Do 25 J = 1, 3
   25         PxBest(J) = PxLab(J)
            BestND = Diff
            Call AMove(9,Tmp,TBest)
            endIf
   30     Continue
   40   Continue
      If(DoPrnt) then
        If(IOpt.eq.2) then
          IX1 = IXBest
          If(IX1.lt.0) IX1 = 4 + IX1
          Write(IOut,1020) IX1, BestND
        else
          Write(IOut,1000) JPBest, IXBest, PxBest, BestND
          endIf
        endIf
      If(DoPrMt) Call OutMtS(IOut,'TBest',0,0,TBest,3,3,3,3)
      Return
      End
*Deck RotF
      Subroutine RotF(NAtoms,T,A,B)
      Implicit Real*8(A-H,O-Z)
C     
C     Routine to rotate energy derivatives back to original cartesian
C     axes.  Note that this is a general routine, and is not resticted
C     to derivatives.
C     
C     Arguments:
C     
C     NAtoms ... Number of Atoms.
C     T ... (3 by 3) rotation matrix.
C     A ... Input vector of length 3*NAtoms,
C     B ... Output vector.
C     
      Dimension T(3,3), A(3,NAtoms), B(3,NAtoms)
C     
      Do 30 I = 1, NAtoms
       A1 = A(1,I)
       A2 = A(2,I)
       A3 = A(3,I)
       B(1,I) = T(1,1)*A1 + T(2,1)*A2 + T(3,1)*A3 
       B(2,I) = T(1,2)*A1 + T(2,2)*A2 + T(3,2)*A3
       B(3,I) = T(1,3)*A1 + T(2,3)*A2 + T(3,3)*A3
   30 Continue
      Return
      End
*Deck RotF1
      Subroutine RotF1(NAtoms,T,A)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to rotate energy derivatives back to original cartesian
C     axes in-place.  Note that this is a general routine, and is not
C     resticted to derivatives, although the transpose of T should be
C     provided for displacements.
C
C     Arguments:
C
C     NAtoms ... Number of Atoms.
C     T ... (3 by 3) rotation matrix.
C     A ... Input vector of length 3*NAtoms,
C
      Dimension T(3,3), A(3,NAtoms)
C
      Do 30 I = 1, NAtoms
       A1 = A(1,I)
       A2 = A(2,I)
       A3 = A(3,I)
       A(1,I) = T(1,1)*A1 + T(2,1)*A2 + T(3,1)*A3
       A(2,I) = T(1,2)*A1 + T(2,2)*A2 + T(3,2)*A3
       A(3,I) = T(1,3)*A1 + T(2,3)*A2 + T(3,3)*A3
   30 Continue
      Return
      End
*Deck RotFt
      Subroutine RotFt(NAtoms,T,A)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to rotate in place coordinates to a new orientation. 
C     It uses the transpose of  the rotational matrix      
C
C     Arguments:
C
C     NAtoms ... Number of Atoms.
C     T ... (3 by 3) rotation matrix.
C     A ... Input vector of length 3*NAtoms,
C
      Dimension T(3,3), A(3,NAtoms)
C
      Do 30 I = 1, NAtoms
        A1 = A(1,I)
        A2 = A(2,I)
        A3 = A(3,I)
        A(1,I) = T(1,1)*A1 + T(1,2)*A2 + T(1,3)*A3
        A(2,I) = T(2,1)*A1 + T(2,2)*A2 + T(2,3)*A3
   30   A(3,I) = T(3,1)*A1 + T(3,2)*A2 + T(3,3)*A3
      Return
      End
*Deck RotCon
      Subroutine RotCon(IOut,IPrint,Linear,NAtoms,PhyCon,PMom,
     $  RotGHz,RTemp)
      Implicit Real*8(A-H,O-Z)
C
C     Compute and print rotational constants.
C
      Logical Linear
      Dimension PhyCon(*), PMom(3),RotGHz(3),RTemp(3)
      Save Zero, One, Four, Eight, AGiga
      Data Zero/0.0d0/, One/1.0d0/, Four/4.0d0/, Eight/8.0d0/,
     $     AGiga/1.0d9/
 1000 Format(' Rotational constants (GHZ)   ',3F15.7)
 1100 Format(' Rotational temperatures (K)  ',3F15.7)
C
C     Get needed constants
C
      If(NAtoms.lt.2) Return
      Boltz  = Phycon(10)
      Planck = Phycon(4)
      ToMet  = Phycon(7)
      ToKG   = Phycon(2)
      Pi = Four * ATan(One)
      PiPi = Pi * Pi
      Con = Planck / (Boltz*Eight*PiPi)
      Con = (Con / ToKG)  *  (Planck / (ToMet*ToMet))
C
C     Compute rotational constants in GHZ.  Beware of linear molecules!
C
      Call AClear(3,RTemp)
      If(.not.Linear.and.PMom(1).ne.Zero) RTemp(1) = Con / PMom(1)
      If(PMom(2).ne.Zero) RTemp(2) = CON / PMom(2)
      If(PMom(3).ne.Zero) RTemp(3) = CON / PMom(3)
      RotGHz(1) = Boltz * RTemp(1) / AGiga / Planck
      RotGHz(2) = Boltz * RTemp(2) / AGiga / Planck
      RotGHz(3) = Boltz * RTemp(3) / AGiga / Planck
      If(IPrint.gt.0) then
       Write(IOut,1000) (RotGHz(i),i=1,3)
       Write(IOut,1100) (RTemp(i),i=1,3)
      EndIf 
      Return
      End
*Deck RVdWMK
      Function RVdWMK(IA)
      Implicit Real*8(A-H,O-Z)
C
C     Assignes van der Waals radii to atoms.
C
C     (Merz-Kollman values from G94 for I, II, III periods; Bondi
C      values for Ar and successive)
C
      Parameter (MaxIA=110)
      Real*8 Radius(MaxIA)
      Save Radius
      Data Radius/
C        H            He
     $ 1.20d+00,     1.20d+00,
C       Li            Be           B               C            N
     $ 1.37d+00,     1.45d+00,    1.45d+00,       1.50d+00,    1.50d+00,
C       O             F            Ne
     $ 1.40d+00,     1.35d+00,    1.30d+00,
C       Na            Mg           Al              Si           P
     $ 1.57d+00,     1.36d+00,    1.24d+00,       1.17d+00,    1.90d+00,
C       S             Cl           Ar
     $ 1.85d+00,     1.80d+00,    1.88d+00,
C       K             Ca (n.a.)
     $ 2.75d+00,     0.00d+00,
C       Sc(n.a.)      Ti (n.a)     V (n.a.)        Cr (n.a.)    Mn (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Fe            Co           Ni              Cu           Zn
     $ 0.00d+00,     0.00d+00,    1.63d+00,       1.40d+00,    1.39d+00,
C       Ga            Ge (n.a.)    As              Se           Br
     $ 1.87d+00,     1.86d+00,    2.00d+00,       2.00d+00,    1.95d+00,
C       Kr
     $ 2.02d+00,
C       Rb (n.a.)     Sr (n.a.)
     $ 0.00d+00,     0.00d+00,
C       Y  (n.a.)     Zr (n.a.)    Nb (n.a.)       Mo (n.a.)    Tc (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ru (n.a.)     Rh (n.a.)    Pd              Ag           Cd
     $ 0.00d+00,     0.00d+00,    1.63d+00,       1.72d+00,    1.58d+00,
C       In            Sn           Sb (n.a.)       Te           I
     $ 1.93d+00,     2.17d+00,    2.20d+00,       2.20d+00,    2.15d+00,
C       Xe
     $ 2.16d+00,
C       Cs (n.a.)     Ba (n.a.)
     $ 0.00d+00,     0.00d+00,
C       La (n.a.)     Ce (n.a.)    Pr (n.a.)       Nd (n.a.)    Pm (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Sm (n.a.)     Eu (n.a.)    Gd (n.a.)       Tb (n.a.)    Dy (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ho (n.a.)     Er (n.a.)    Tm (n.a.)        Yb (n.a.)    Lu (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Hf (n.a.)
     $ 0.00d+00,
C       Ta (n.a.)     W (n.a.)     Re (n.a.)       Os (n.a.)    Ir (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Pt            Au           Hg              Tl           Pb
     $ 1.72d+00,     1.66d+00,    1.55d+00,       1.96d+00,    1.02d+00,
C       Bi (n.a.)     Po (n.a.)    At (n.a.)       Rn (n.a.)    Fr (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ra (n.a.)     Ac (n.a.)    Th (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,
C       Pa (n.a.)     U            Np (n.a.)       Pu (n.a.)    Am (n.a.)
     $ 0.00d+00,     1.86d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Cm (n.a.)     Bk (n.a.)    Cf (n.a.)       Es (n.a.)    Fm (n.a.)
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Md (n.a.)
     $ 0.00d+00,   9*0.00d+00/
C
      RvdWMK = Radius(Min(Max(IA,0),MaxIA))
      Return
      End
*Deck RVdW97
      Function RVdW97(IA)
      Implicit Real*8(A-H,O-Z)
C
C     Assignes Van der Waals radii to atoms.
C     (Values from the UFF paper)
C
      Parameter (MaxIA=110)
      Real*8 Rii(0:MaxIA)
      Save Pt5,Rii
      Data Pt5/0.5d0/
C
C     Warning: the following are atomic diameters.
C
      Data Rii/0.000d+00,
C        H             He
     $ 2.886d+00,    2.362d+00,
C       Li            Be           B               C           N
     $ 2.451d+00,    2.745d+00,   4.083d+00,      3.851d+00,  3.660d+00,
C       O             F            Ne
     $ 3.500d+00,    3.364d+00,   3.243d+00,
C       Na            Mg           Al              Si          P
     $ 2.983d+00,    3.021d+00,   4.499d+00,      4.295d+00,  4.147d+00,
C       S             Cl           Ar
     $ 4.035d+00,    3.947d+00,   3.868d+00,
C       K             Ca
     $ 3.812d+00,    3.399d+00,
C       Sc            Ti           V               Cr          Mn
     $ 3.295d+00,    3.175d+00,   3.144d+00,      3.023d+00,  2.961d+00,
C       Fe            Co           Ni              Cu          Zn
     $ 2.912d+00,    2.872d+00,   2.834d+00,      3.495d+00,  2.763d+00,
C       Ga            Ge           As              Se          Br
     $ 4.383d+00,    4.280d+00,   4.230d+00,      4.205d+00,  4.189d+00,
C       Kr
     $ 4.141d+00,
C       Rb            Sr (+2)
     $ 4.114d+00,    3.641d+00,
C       Y  (+3)       Zr (+4)      Nb (+5)         Mo (+6)     Tc (+5)
     $ 3.345d+00,    3.124d+00,   3.165d+00,      3.052d+00,  2.998d+00,
C       Ru (+2)       Rh (+3)      Pd (+2)         Ag (+1)     Cd (+2)
     $ 2.963d+00,    2.929d+00,   2.899d+00,      3.148d+00,  2.848d+00,
C       In            Sn           Sb              Te          I
     $ 4.463d+00,    4.392d+00,   4.420d+00,      4.470d+00,  4.50d+00,
C       Xe
     $ 4.404d+00,
C       Cs            Ba (+2)
     $ 4.517d+00,    3.703d+00,
C       La (+3)       Ce (+3)      Pr (+3)         Nd (+3)      Pm (+3)
     $ 3.522d+00,    3.556d+00,   3.606d+00,      3.575d+00,  3.547d+00,
C       Sm (+3)       Eu (+3)      Gd (+3)         Tb (+3)      Dy (+3)
     $ 3.520d+00,    3.493d+00,   3.368d+00,      3.451d+00,  3.428d+00,
C       Ho (+3)       Er (+3)      Tm (+3) )        Yb (+3)     Lu (+3)
     $ 3.409d+00,    3.391d+00,   3.374d+00,      3.355d+00,  3.640d+00,
C       Hf (+4)
     $ 3.141d+00,
C       Ta (+5)       W (+4,+6)    Re (+5,+7)      Os (+6)      Ir (+3)
     $ 3.170d+00,    3.069d+00,   2.954d+00,      3.120d+00,  2.840d+00,
C       Pt            Au           Hg              Tl          Pb
     $ 2.754d+00,    3.293d+00,   2.705d+00,      4.347d+00,  4.297d+00,
C       Bi (+3)       Po (+2)      At              Rn (+4)      Fr
     $ 4.370d+00,    4.709d+00,   4.750d+00,      4.765d+00,  4.900d+00,
C       Ra (+2)       Ac (+3)      Th (+4)
     $ 3.677d+00,    3.478d+00,   3.396d+00,
C       Pa (+4)       U (+4)       Np (+4)         Pu (+4)      Am (+4)
     $ 3.424d+00,    3.395d+00,   3.424d+00,      3.424d+00,  3.381d+00,
C       Cm (+3)       Bk (+3)      Cf (+3)         Es (+3)      Fm (+3)
     $ 3.326d+00,    3.339d+00,   3.313d+00,      3.299d+00,  3.286d+00,
C       Md (+3)       No (+3)      Lw(+3)
     $ 3.274d+00,    3.248d+00,   3.236d+00,    7*3.500D+00/
C
      RVdW97 = Pt5*Rii(Min(Max(IA,0),MaxIA))
      Return
      End
*Deck SProd
      Function SProd(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns the scalar product of vectors A and B.
C     Same as SProd but a separate GPU kernel.
C
      Parameter (Zero=0.0D0)
      Dimension A(*), B(*)
C
      Prod = Zero
      Do 10 I = 1, N
   10   Prod = Prod + A(I) * B(I)
      SProd = Prod
      Return
      End
*Deck Square
      Subroutine Square(A,B,Max,N,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Places linear array in square form:
C
C     1  Input is square and just unpacked.
C     0  Input is symmetric and stored lower triangular.
C    -1  Input is antisymmetric and stored lower triangular
C        (including the 0 diagonal).
C     2  Input is lower triangular; expand A to upper triangular
C          matrix with 0's in lower triangular
C     3  Input is symmetric lower triangular and A and B are different.
C    -3  Input is antisymmetric lower triangular and A and B are different.
C
      Common/IO/In,IOut,IPunch
      Dimension A(*), B(Max,Max)
      Save Zero
      Data Zero/0.d0/
C
C     Key=1:  Expand square array.
C
      If(Key.eq.1) then
        Do 10 JX = N, 1, -1
          KX = N*(JX-1)
          Do 10 IX = N, 1, -1
            K = KX + IX
            B(IX,JX) = A(K)
   10       Continue
C
C     Key = 0:  Symmetric lower triangle.
C
      else if(Key.eq.0) then
        Do 20 JX = N, 1, -1
          KX = JX*(JX-1)/2
          Do 20 IX = JX, 1, -1
            K = KX + IX
            B(IX,JX) = A(K)
   20       Continue
        Do 30 J = 1, N
          Do 30 I = 1, (J-1)
   30       B(J,I) = B(I,J)
C
C     Key = -1:  Antisymmetric lower triangle.
C
      else if(Key.eq.-1) then
        Do 40 JX = N, 1, -1
          KX = JX*(JX-1)/2
          Do 40 IX = JX, 1, -1
            K = KX + IX
            B(IX,JX) = -A(K)
   40       Continue
        Do 50 J = 1, N
          Do 50 I = 1, J
   50       B(J,I) = -B(I,J)
C
C     Key=2:  Expand A to be an upper triangular matrix
C
      else if(Key.eq.2) then
        Do 60 JX = N, 1, -1
          KX = JX*(JX-1)/2
          Do 60 IX = JX, 1, -1
            K = KX + IX
            B(IX,JX) = A(K)
   60       Continue
        Do 70 JX = 1, (N-1)
          Do 70 IX = (JX+1), N
   70       B(IX,JX) = Zero
C
C     Key=3:  Expand lower triangular array efficiently
C
      else if(Key.eq.3) then
        ISt = Mod(N,4) + 1
        Do 110 I = 1, (ISt-1)
          II = (I*(I-1))/2
          Do 110 J = 1, I
            B(J,I) = A(II+J)
  110       B(I,J) = A(II+J)
        Do 130 I = ISt, N, 4
          II0 = (I*(I-1))/2
          II1 = II0 + I
          II2 = II1 + I + 1
          II3 = II2 + I + 2
          Do 120 J = 1, I
            AIJ0 = A(II0+J)
            AIJ1 = A(II1+J)
            AIJ2 = A(II2+J)
            AIJ3 = A(II3+J)
            B(J,I) = AIJ0
            B(J,I+1) = AIJ1
            B(J,I+2) = AIJ2
            B(J,I+3) = AIJ3
            B(I,J) = AIJ0
            B(I+1,J) = AIJ1
            B(I+2,J) = AIJ2
  120       B(I+3,J) = AIJ3
          AIJ1 = A(II1+I+1)
          AIJ2 = A(II2+I+1)
          AIJ3 = A(II3+I+1)
          B(I+1,I+1) = AIJ1
          B(I+1,I+2) = AIJ2
          B(I+1,I+3) = AIJ3
          B(I+2,I+1) = AIJ2
          B(I+3,I+1) = AIJ3
          AIJ2 = A(II2+I+2)
          AIJ3 = A(II3+I+2)
          B(I+2,I+2) = AIJ2
          B(I+2,I+3) = AIJ3
          B(I+3,I+2) = AIJ3
          AIJ3 = A(II3+I+3)
  130     B(I+3,I+3) = AIJ3
C
C     Key=-3:  Expand lower triangular array efficiently
C
      else if(Key.eq.-3) then
        ISt = Mod(N,4) + 1
        Do 210 I = 1, (ISt-1)
          II = (I*(I-1))/2
          Do 210 J = 1, I
            B(J,I) = -A(II+J)
  210       B(I,J) = A(II+J)
        Do 230 I = ISt, N, 4
          II0 = (I*(I-1))/2
          II1 = II0 + I
          II2 = II1 + I + 1
          II3 = II2 + I + 2
          Do 220 J = 1, I
            AIJ0 = A(II0+J)
            AIJ1 = A(II1+J)
            AIJ2 = A(II2+J)
            AIJ3 = A(II3+J)
            B(J,I) = -AIJ0
            B(J,I+1) = -AIJ1
            B(J,I+2) = -AIJ2
            B(J,I+3) = -AIJ3
            B(I,J) = AIJ0
            B(I+1,J) = AIJ1
            B(I+2,J) = AIJ2
  220       B(I+3,J) = AIJ3
          B(I+1,J) = AIJ1
          B(I+2,J) = AIJ2
          B(I+3,J) = AIJ3
          AIJ2 = A(II2+I+1)
          AIJ3 = A(II3+I+1)
          B(I+1,I+1) = Zero
          B(I+1,I+2) = -AIJ2
          B(I+1,I+3) = -AIJ3
          B(I+2,I+1) = AIJ2
          B(I+3,I+1) = AIJ3
          AIJ3 = A(II3+I+2)
          B(I+2,I+2) = Zero
          B(I+2,I+3) = -AIJ3
          B(I+3,I+2) = AIJ3
  230     B(I+3,I+3) = Zero
      else
       write(IOut,'('' Illegal value for Key in Square'')')
       Stop
      endIf
      Return
      End
*Deck SrtEig
      Subroutine SrtEig(IDir,Small,NDim,NV,N,Eig,V,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     Sort the eigenvectors in V by the eigenvalues in Eig.
C     IDir = +/-1 for increasing/decreasing order.  This version
C     does extra work in order to be careful to maintain the order
C     of (near)-degenerate eigenvectors.
C
      Dimension Eig(*), V(NDim,*), Scr(*)
C
      If(Abs(IDir).ne.1) Return
      Do 30 I = 1, (N-1)
        ISub = 1
        SX = Eig(1)
        Last = N - I + 1
        If(IDir.gt.0) then
          Do 10 J = 2, Last
            If(Eig(J).ge.(SX-Small)) then
              ISub = J
              SX = Eig(J)
              endIf
   10       Continue
        else
          Do 20 J = 2, Last
            If(Eig(J).le.(SX+Small)) then
              ISub = J
              SX = Eig(J)
              endIf
   20       Continue
          endIf
        If(ISub.ne.Last) then
          Call AMove1(Last-ISub,ISub,ISub-1,Eig)
          Eig(Last) = SX
          Call AMove(NV,V(1,ISub),Scr)
          Call AMove1(NDim*(Last-ISub),NDim*ISub,NDim*(ISub-1),V)
          endIf
   30   Continue
      Return
      End
*Deck St2Dat
      Subroutine St2Dat(String,IBDlim,IWDlim,ICur,IType,CType,IVal,
     $  RVal,CVal,LVal,SVal)
      Implicit Integer(A-Z)
C
C String to data
C     Convert a non-null string element to the proper Fortran Data type
C     (integer, real, complex, character string, logical) and increment
C     the reading cursor in String
C
C Input:
C     String : String to analyze
C     IBDlim : Block delimiters (quotes, parentheses) to be removed
C              -1: Ignore any block delimiter
C               0: Do not remove anything
C              .1: Remove parentheses
C              1.: Remove quotes
C     IWDlim : Word delimiters (spaces, commas) to separate each element
C                0: spaces (always included)
C               .1: tabs
C               1.: , ;
C              1..: =
C
C InOut:
C     ICur   : Position cursor in String, incremented on output
C
C Output:
C     IType  : Data type, as integer
C              -2: open delimiter (unclosed)
C              -1: end-of-line found
C               0: Nothing found (null)
C               1: Integer
C               2: Real
C               3: Complex
C               4: Logical
C               5: Character
C     CType  : Data type, as character
C              O: open delimiter (unclosed)
C              E: end-of-line found
C              N: Nothing found (null)
C              I: Integer
C              R: Real
C              C: Complex
C              L: Logical
C              S: Character
C     IVal   : Converted integer value
C     RVal   : Converted real value
C     CVal   : Converted complex value
C     LVal   : Converted logical value
C     SVal   : Converted character value
C
C Note:
C     FOR NOW, IBDlim AND IWDlim ARE NOT SUPPORTED
C
C     Dimensions
      Common/IO/in,iout,ipunch
      Integer MXWDLM, MXBDLM
      Parameter (MXWDLM=5, MXBDLM=6)
C     Input
      Integer IBDlim, IGet10, IWDlim
      Character*(*) String
C     InOut
      Integer ICur
C     Output
      Integer IType, IVal
      Real*8 RVal
      Complex*16 CVal
      Logical LVal
      Character CType*1, SVal*(*)
C     Local
      Integer ExpSgn, i, IBS, IComma, ICur1, IDot, IDot2, IExp, IType0,
     $  IType1, IWE, IWE1, IWS, IWS1, LBDlim, LChain, LenDec, LWDlim,
     $  NumDec, NumExp, NumInt, NumSgn
      Real*8 X, Zero
      Logical PotTyp(4), HasInt, HasDec, HasExp, InPar, NoBloc
      Character ABDlmE*(MXBDLM), ABDlmS*(MXBDLM), BDlimE*(MXBDLM),
     $  BDlimS*(MXBDLM), BlcTyp*4, CurS*1, UpBloc*8, WDLim*(MXWDLM)
      Save Zero
      Data Zero/0.0D0/
C
      ABDlmE = ')]}>''"'
      ABDlmS = '([{<''"'
      IWS = 0
      IWE = 0
      InPar = .False.
      LChain = len(String)
C     -- DEFINE ACCEPTED WORD DELIMITERS --
      WDlim  = ' '
      LWDlim = 1
      If(IWDlim.gt.0) then
        If(IGet10(IWDlim,0).eq.1) then
          WDlim  = WDlim(:LWDlim) // char(9)
          LWDlim = LWDlim + 1
        endIf
        If(IGet10(IWDlim,1).eq.1) then
          WDlim  = WDlim(:LWDlim) // ',;'
          LWDlim = LWDlim + 2
        endIf
        If(IGet10(IWDlim,2).eq.1) then
          WDlim  = WDlim(:LWDlim) // '='
          LWDlim = LWDlim + 1
        endIf
      endIf
C     -- DEFINE ACCEPTED BLOCK DELIMITERS --
      NoBloc = IBDlim.eq.-1
      BDlimS = ' '
      BDlimE = ' '
      LBDlim = 0
      If(IBDlim.gt.0) then
        If(IGet10(IBDlim,0).eq.1) then
          BDlimS = BDlimS(:LBDlim) // '([{<'
          BDlimE = BDlimE(:LBDlim) // ')]}>'
          LBDlim = LBDlim + 4
        endIf
        If(IGet10(IBDlim,0).eq.1) then
          BDlimS = BDlimS(:LBDlim) // '''"'
          BDlimE = BDlimE(:LBDlim) // '''"'
          LBDlim = LBDlim + 2
        endIf
      endIf
C     -- CASE: CURSOR OUTSIDE STRING --
      If(ICur.gt.LChain.or.ICur.lt.0) then
        IType = -1
        CType = 'E'
        Return
      endIf
C     -- FIND WORD --
C     ----- FIND BEGINNING -----
      ICur1 = ICur
  100 If(Index(WDlim(:LWDlim),String(ICur1:ICur1)).ge.1) then
        ICur1 = ICur1 + 1
        If(ICur1.le.LChain) then
          Goto 100
        else
          IType = 0
          CType = 'N'
          Return
        endIf
      else
        IWS = ICur1
        IBS = Index(ABDlmS,String(ICur1:ICur1))
        InPar = IBS.eq.1
      endIf
C     ----- FIND END -----
      If(IWS.gt.0) then
        If((NoBloc.and..not.InPar).or.IBS.eq.0) then
  200     If(Index(WDlim(:LWDlim),String(ICur1:ICur1)).eq.0) then
            ICur1 = ICur1 + 1
            If(ICur1.le.LChain) then
              Goto 200
            else
              IWE = LChain
            endIf
          else
            IWE = ICur1 - 1
          endIf
        else
  300     If(String(ICur1:ICur1).ne.ABDlmE(IBS:IBS)) then
            ICur1 = ICur1 + 1
            If(ICur1.le.LChain) then
              Goto 300
            else
              IType = -2
              CType = 'O'
              Return
            endIf
          else
            IWE = ICur1
            ICur1 = ICur1 + 1
          endIf
        endIf
      endIf
C     -- WORD ANALYSIS --
      If(IWS.gt.0.and.IWE.gt.0) then
        IComma = 0
        If(InPar) then
C         PARENTHESES ARE REMOVED AS BLOCK DELIMITERS FOR ANALYSIS
C         Look for possible (x,y), which can be a complex number
          i = Index(String(IWS:IWE),',')
          If(i.gt.0) then
            IWS1 = IWS + 1
            IComma = IWS + i - 1
            IWE1 = IComma - 1
          else
            IWS1 = IWS + 1
            IWE1 = IWE - 1
          endIf
        else if(IBS.gt.0) then
          If(NoBloc) then
            IWS1 = IWS + 1
            IWE1 = IWE
          else if(IBDlim.gt.0) then
            IWS1 = IWS + 1
            IWE1 = IWE - 1
          endIf
        else
          IWS1 = IWS
          IWE1 = IWE
        endIf
C       ---- ANALYSIS OF BLOCK/SUB-BLOCK ----
        IType1 = 0
C       PotTyp: Potential types
C          1: integer
C          2: real
C          3: logical
C          4: character
C       Note: Complexes are assumed to be to int/real combined  with ,
  400   PotTyp(1) = .True.
        PotTyp(2) = .True.
        PotTyp(3) = .True.
        PotTyp(4) = .True.
        BlcTyp = 'clri'
        IType0 = 0
        IDot   = 0
        IDot2  = 0
        IExp   = 0
        HasInt = .False.
        HasDec = .False.
        HasExp = .False.
        NumInt = 0
        NumDec = 0
        LenDec = 0
        NumSgn = 1
        NumExp = 0
        ExpSgn = 1
C       Remove spaces present after block delimiters
  420   If(String(IWS1:IWS1).eq.' ') then
          IWS1 = IWS1 + 1
          Goto 420
        endIf
        i = IWS1
  430   CurS = String(i:i)
        If(Index('1234567890',CurS).gt.0) then
          PotTyp(3) = .False.
          BlcTyp(2:2) = ' '
          If(PotTyp(1).or.PotTyp(2)) then
            If(IExp.gt.0) then
              NumExp = NumExp*10 + Mod(Index('1234567890',CurS),10)
              HasExp = .True.
            else if(IDot.gt.0) then
              NumDec = NumDec*10 + Mod(Index('1234567890',CurS),10)
              LenDec = LenDec + 1
              HasDec = .True.
            else
              NumInt = NumInt*10 + Mod(Index('1234567890',CurS),10)
              HasInt = .True.
            endIf
          endIf
        else if(Index('+-',CurS).gt.0) then
          PotTyp(3) = .False.
          BlcTyp(2:2) = ' '
          If(PotTyp(1).or.PotTyp(2)) then
            If(i.eq.IWS1) then
              If(CurS.eq.'-') NumSgn = -1
            else if(i.eq.IExp+1) then
              If(CurS.eq.'-') then
                ExpSgn = -1
                PotTyp(1) = .False.
                BlcTyp(4:4) = ' '
              endIf
            else
              PotTyp(1) = .False.
              PotTyp(2) = .False.
              BlcTyp(3:4) = '  '
            endIf
          endIf
        else if(CurS.eq.'.') then
          If(IDot.gt.0) then
C           Two points: Only string or logical
            PotTyp(1) = .False.
            PotTyp(2) = .False.
            BlcTyp(3:4) = '  '
            IDot2 = i
          else
            PotTyp(1) = .False.
            BlcTyp(4:4) = ' '
            IDot = i
          endIf
        else if(Index('DdEe',CurS).gt.0) then
          If(IExp.gt.0) then
            PotTyp(1) = .False.
            PotTyp(2) = .False.
            BlcTyp(3:4) = '  '
          else
            IExp = i
          endIf
        else if(CurS.eq.' ') then
          If(String(i+1:IWE1).eq.' ') then
            IWE1 = i
          else
            PotTyp(1) = .False.
            PotTyp(2) = .False.
            PotTyp(3) = .False.
            BlcTyp(2:4) = '   '
          endIf
        else
          PotTyp(1) = .False.
          PotTyp(2) = .False.
          BlcTyp(3:4) = '  '
        endIf
        If(BlcTyp.ne.'c   ') then
          i = i + 1
          If(i.le.IWE1) Goto 430
        endIf
C       Check imcomplete numbers
        If(PotTyp(1).or.PotTyp(2)) then
          If(.not.(HasInt.or.HasDec)) then
            PotTyp(1) = .False.
            PotTyp(2) = .False.
          endIf
          If(IExp.gt.0.and..not.HasExp) then
            PotTyp(1) = .False.
            PotTyp(2) = .False.
          endIf
        endIf
        If(PotTyp(1)) then
          IVal = NumSgn * NumInt * 10**NumExp
          IType0 = 1
        else if(PotTyp(2)) then
          RVal = Float(NumSgn)
     $      * (Float(NumInt) + Float(NumDec)*Float(10)**(-LenDec))
     $      * Float(10)**(ExpSgn*NumExp)
          IType0 = 2
        else
          If(PotTyp(3)) then
            Call LinUpC(String(IWS1:IWE1),UpBloc)
            If(IDot2.gt.0) then
              If(UpBloc.eq.'.TRUE.'.or.UpBloc.eq.'.T.') then
                LVal = .True.
              else if(UpBloc.eq.'.FALSE.'.or.UpBloc.eq.'.F.') then
                LVal = .False.
              else
                PotTyp(3) = .False.
              endIf
            else if(IDot.gt.0) then
              If(UpBloc.eq.'.TRUE'.or.UpBloc.eq.'.T') then
                LVal = .True.
              else if(UpBloc.eq.'.FALSE'.or.UpBloc.eq.'.F') then
                LVal = .False.
              else
                PotTyp(3) = .False.
              endIf
            else
              If(UpBloc.eq.'TRUE'.or.UpBloc.eq.'T') then
                LVal = .True.
              else if(UpBloc.eq.'FALSE'.or.UpBloc.eq.'F') then
                LVal = .False.
              else
                PotTyp(3) = .False.
              endIf
            endIf
          endIf
          If(.not.PotTyp(3)) then
            SVal = String(IWS1:IWE1)
            IType0 = 5
          else
            IType0 = 4
          endIf
        endIf
        If(IComma.gt.0.and.IType0.le.2) then
          If(IType0.eq.1) then
            CVal = Cmplx(Float(IVal), Zero)
          else
            CVal = Cmplx(RVal, Zero)
          endIf
          IType1 = IType0
          IWS1 = IComma + 1
          IWE1 = IWE - 1
          IComma = 0
          Goto 400
        endIf
C       FINAL CHECK OF TYPE
        If(IType0.gt.0.and.IType1.gt.0) then
          If(IType1.le.2) then
            X = Real(CVal)
            If(IType0.eq.1) then
              CVal = Cmplx( X, Float(IVal) )
              IType = 3
            else if(IType0.eq.2) then
              CVal = Cmplx( X, RVal )
              IType = 3
            else
              SVal = String(IWS:IWE)
              IType = 5
            endIf
          else
            SVal = String(IWS:IWE)
            IType = 5
          endIf
        else
          IType = IType0
        endIf
        If(IType.eq.1) then
          CType = 'I'
        else if(IType.eq.2) then
          CType = 'R'
        else if(IType.eq.3) then
          CType = 'C'
        else if(IType.eq.4) then
          CType = 'L'
        else if(IType.eq.5) then
          CType = 'S'
        endIf
        ICur = ICur1
      else
       write(IOut,'('' Wrong word analysis'')')
       Stop
      endIf
      Return
      End
*Deck St2Int
      Integer Function St2Int(Str,ErrVal)
      Implicit Integer(A-Z)
C
C     Convert string to an integer, returning ErrVal if the conversion
C     fails.
C
      Character*(*) Str, StrL*20
C
      St2Int = ErrVal
      LenS = LinEnd(Str)
      LL = Len(StrL)
      If(LenS.lt.1.or.LenS.gt.LL) Return
      StrL = ' '
      StrL(LL+1-LenS:) = Str(1:LenS)
      Read(StrL,'(I20)',Err=900) St2Int
      Return
  900 St2Int = ErrVal
      Return
      End
*Deck SUBSTR
      Subroutine SubStr(LINE,NMAX,ISTART,NVALUE)
      Implicit Real*8 (A-H,O-Z)
C
C   Finds the number(NVALUE) and starting position (ISTART) of substrings in
C   the string LINE. NMAX is the dimension of ISTART in the calling program
C
      character LINE*(*)
      Logical LEADSP
      Dimension ISTART(*)
      do 10 I = 1,NMAX
       ISTART(I) = 80
   10 Continue
      LMAX = LEN(LINE)
      LEADSP=.TRUE.
      NVALUE=0
      DO 20 I=1,LMAX
       IF (LEADSP.AND.LINE(I:I).NE.' ') THEN
        NVALUE=NVALUE+1
        IF(NVALUE.GT.NMAX) GO TO 20
        ISTART(NVALUE)=I
       END IF
       LEADSP=(LINE(I:I).EQ.' ')
   20 CONTINUE
      Return
      End
*Deck SymNum
      Function Symnum(Linear,PG)
      Implicit Real*8(A-H,O-Z)
C
C     Provide the calling routine with the rotational symmetry number
C     of the molecule.  A table look up is done based upon the point
C     group of the molecule.  See S. W. Benson, "Thermochemical
C     Kinetics, 2nd ed.", Wiley, New York, 1976, P49.  In addition, the
C     logical variable Linear is set true if the molecular point group
C     is D*H OR C*V.  PG is the point group as an unpacked character
C     string (the first part of the framework group string.
C
      Integer PG(*), Num(10)
      Logical Linear, Test
      Save One,Two,Twelve,F24,IHSt,Num,IHC,IHI,IHS,IHD,IHT,IHO
      Data One,Two,Twelve,F24/1.D0,2.D0,12.D0,24.D0/, IHSt/1h*/
      Data Num/1H0, 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9/
      Data IHC/1hC/, IHI/1hI/, IHS/1hS/, IHD/1hD/, IHT/1hT/, IHO/1hO/
C
      Symnum = One
      Linear = PG(2) .EQ. IHSt
      N = Max(Numer(PG),1)
      Test = .False.
      Do 10 I = 1, 10
   10   Test = Test.or.PG(3).eq.Num(I)
      If(.not.Test) N = N / 10
C
C     CI, CS, CN, CNH, CNV.
C
      If(PG(1).eq.IHC) then
        If(Linear.or.PG(2).eq.IHI.or.PG(2).eq.IHS) then
          SymNum = One
        else
          SymNum = Float(N)
          endIf
C
C     DN, DNH, DND.
C
      else if(PG(1).eq.IHD) then
        SymNum = Two * Float(N)
        If(Linear) SymNum = Two
C
C     SN.
C
      else if(PG(1).eq.IHS) then
        SymNum = Float(N) / Two
C
C     T, TD, O, OH.
C
      else if(PG(1).eq.IHT) then
        SymNum = Twelve
      else if(PG(1).eq.IHO) then
        SymNum = F24
      else
        SymNum = One
        endIf
      Return
      End
*Deck ThrEle
      Subroutine ThrEle(PhyCon,Multip,TUser,QEle,SEle,Eele,CEle)
      Implicit Real*8 (A-H,O-Z)
      Dimension PhyCon(*)
C
C     COMPUTE CONTRIBUTIONS DUE TO ELECTRONIC MOTION
C        IT IS ASSUMED THAT THE FIRST ELECTRONIC EXCITATION ENERGY
C        IS MUCH GREATER THAN KT.
C        QELEC-- PARTITION FUNCTION
C        QLELEC-- COMMON  LOGARITHM OF Q
C        QLNELE-- NATURAL LOGARITHM OF Q
C        EELEC-- INTERNAL ENERGY
C        CELEC-- HEAT CAPACITY
C        SELEC-- ENTROPY
C
      Save Zero,TStd,T0C
      Data Zero /0.0D0/
      Data TStd/298.15D0/,T0C/273.15D0/
      Avog   = PhyCon(5)
      Boltz  = PhyCon(10)
      Gas  = Avog * Boltz
      T = TUser
      IF(T.Eq.Zero) T = TStd 
      RT = Gas * T
      Degen = Float(Multip)
      QELec = Degen 
      QLElec = Log10(QElec)
      QLnEle = Log(QElec)
      EElec = Zero * RT
      CElec = Zero * Gas
      SElec = Log(Degen) * Gas
      Return
      End
*Deck ThrTra
      Subroutine ThrTra(PhyCon,TotWt,PUser,TUser,QTran,STran,ETran,
     $  CTran)
      Implicit Real*8 (A-H,O-Z)
      Dimension PhyCon(*)
C
C     COMPUTE CONTRIBUTIONS DUE TO TRANSLATION from
C     TotWt = mass in amu, PUser = Pressure in Atm, TUser = Temperature in K
C        ETRAN-- INTERNAL ENERGY
C        CTRAN-- CONSTANT V HEAT CAPACITY
C        STRAN-- ENTROPY
C        QTRAN-- PARTITION FUNCTION
C        QLTRAN-- COMMON LOG OF PARTITION FUNCTION
C        QLNTRA-- NATURAL LOG OF PARTITION FUNCTION
C
      Save Zero,One,OnePt5,Two,TwoPt5,Four,TStd,T0C
      Data Zero /0.0D0/
      Data One,OnePt5,Two,TwoPt5,Four/1.0D0,1.5D0,2.0D0,2.5D0,4.0D0/
      Data TStd/298.15D0/,T0C/273.15D0/
      TOKG   = PhyCon(2)
      BOLTZ  = PhyCon(10)
      PLANCK = PhyCon(4)
      AVOG   = PhyCon(5)
      VolMol = PhyCon(13)
C
C     Compute the GAS constant and Pi
C
      Gas  = Avog * Boltz
      Pi   = Four * ATan(One)
C
C     Compute translational contributions
C
      WEIGHT = TotWt * TOKG
      T = TUSER
      IF(T.EQ.ZERO) T = TSTD
      PATM = PUSER
      IF(PATM.EQ.ZERO) PATM = ONE
      P = PATM * T0C * Gas / VolMol
      RT = GAS * T
      DUM1 = BOLTZ * T
      DUM2 = (TWO*PI) ** ONEPT5
      ARG  = DUM1 ** ONEPT5  / PLANCK
      ARG  = (ARG/P) * (DUM1/PLANCK)
      ARG  = ARG * DUM2 * (WEIGHT/PLANCK)
      QTRAN = ARG * Sqrt(WEIGHT)
      QLTRAN = Log10(QTRAN)
      QLNTRA = Log(QTRAN)
      STRAN = GAS * (Log(QTRAN)+TWOPT5)
      ETRAN = ONEPT5 * RT
      CTRAN = ONEPT5 * GAS
      Return
      End
*Deck THrVib 
      Subroutine ThrVib(IOut,IWrite,NVib,PhyCon,PUser,TUser,Freq,VTemp,
     $  QVib,QZVib,SVib,EVib,CVib)
      Implicit Real*8 (A-H,O-Z)
C
C     COMPUTE CONTRIBUTIONS DUE TO VIBRATION.
C     COMPUTE VIBRATIONAL TEMPERATURES AND ZERO POINT VIBRATIONAL
C     ENERGY.  ONLY REAL FREQUENCIES ARE INCLUDED IN THE ANALYSIS.
C
      Logical Ing,DoJoul
      Real*8 JPCal,MdCutO,MxToVT
      Dimension PhyCon(*),Freq(*),VTemp(*)
      Common/LowVib/EVibN(20),CVibN(20),SVibN(20),QVibN(20),QZVibN(20)
      Save Zero,One,Four,Fract,Half,MxToVt,AKilo,T0C,TStd
      Data Zero/0.0D0/,One/1.0D0/,Four/4.0D0/,Fract/5.0D-2/,Half/5.0D-1/
      Data MxToVT/8.50D1/, AKilo/1.0D3/, T0C/2.7315D2/,  TStd/2.9815D2/
 1000 Format(1X,I4,' imaginary frequencies ignored.')
 1100 Format(' Frequencies are scaled by ',F10.5)
 1200 Format(' Temperature',F10.3,' Kelvin.  Pressure',F10.5,' Atm.')
 1300 Format(' Zero-point vibrational energy ',F12.1,' (Joules/Mol)',
     $      /,31X,F12.5,' (Kcal/Mol)')
 1400 Format(' Warning -- explicit consideration of',I4,
     $ ' degrees of freedom as',
     $  /,1X,'          vibrations may cause significant error')
 1500 Format(' Vibrational temperatures: ',5F9.2)
 1600 Format(10X,'(Kelvin)',9X,5F9.2)
 1700 Format(1X,26X,5F9.2)
C
C     TOKG     KILOGRAMS PER AMU.
C     BOLTZ    BOLTZMAN CONSTANT, IN JOULES PER KELVIN.
C     PLANCK   PLANCK CONSTANT, IN JOULE-SECONDS.
C     AVOG     AVOGADRO CONSTANT, IN MOL**(-1).
C     JPCAL    JOULES PER CALORIE.
C     TOMET    METRES PER BOHR.
C     HARTRE   JOULES PER HARTREE.
C
      Small  = MDCutO(0)
      TOKG   = PhyCon(2)
      PLANCK = PhyCon(4)
      AVOG   = PhyCon(5)
      JPCAL  = PhyCon(6)
      TOMET  = PhyCon(7)
      HARTRE = PhyCon(8)
      BOLTZ  = PhyCon(10)
      VolMol = PhyCon(13)
C
C     COMPUTE THE GAS CONSTANT, PI, PI**2, AND E.
C     COMPUTE THE CONVERSION FACTORS CAL PER JOULE AND KCAL PER JOULE.
C
      Gas  = Avog * Boltz
      Pi   = Four * ATan(One)
      PiPi = Pi * Pi
      ToCal  = One / JPCal
      ToKCal = ToCal / AKilo
C
C     Determine the temperature and pressure (default is 25 C and 1 atm).
C
      T = TUSER
      IF(T.EQ.ZERO) T = TSTD
      PATM = PUSER
      IF(PATM.EQ.ZERO) PATM = ONE
      P = PATM * T0C * Gas / VolMol
      If(IWrite.ge.1) Write(IOut,1200) T,PATM
C
C     Determine NImag and scale frequencies (in increasing order)
C
      NImag=0
      Do 10 IVib=1,NVib
   10  If(Freq(IVib).lt.Small) NImag = NImag + 1
      NDOF = NVib - NImag
      If(NImag.ne.0.and.IWrite.ge.0) Write(IOut,1000) NImag
C
      CON = PLANCK / BOLTZ
      EZPE = ZERO
      DO 20 I = 1, NDOF
        VTEMP(I) = FREQ(I+NImag) * PhyCon(9) * CON
   20   EZPE     = EZPE + FREQ(I+NImag) * PhyCon(9)
      EZPE = HALF * PLANCK * EZPE
      EZJ  = EZPE * AVOG
      EZKC = EZPE * TOKCAL * AVOG
      EZAU = EZPE / HARTRE
      If(IWrite.ge.2) Write(IOut,1300) EZJ, EZKC
C
C     COMPUTE THE NUMBER OF VIBRATIONS FOR WHICH MORE THAN 5  OF AN
C     ASSEMBLY OF MOLECULES WOULD EXIST IN VIBRATIONAL EXCITED STATES.
C     SPECIAL PRINTING FOR THESE MODES IS DONE TO ALLOW THE USER TO
C     EASILY TAKE INTERNAL ROTATIONS INTO ACCOUNT.  THE CRITERION
C     CORRESPONDS ROUGHLY TO A LOW FREQUENCY OF 1.9(10**13) HZ, OR
C     625 CM**(-1), OR A VIBRATIONAL TEMPERATURE OF 900 K.
C
      LOFREQ = 0
      THRESH = -T * Log(FRACT)
      DO 30 I=1,NDOF
        IF(VTEMP(I).LT.THRESH) LOFREQ = LOFREQ + 1
   30   CONTINUE
      If(IWrite.ge.2) then
        If(LoFreq.ne.0) Write(IOut,1400) LoFreq
        ITOP = Min(NDOF,5)
        Write(IOut,1500) (VTemp(I),I=1,ITop)
        If(NDOF.le.5) then
          Write(IOut,1600)
        else
          ITop = Min(NDOF,10)
          Write(IOut,1600) (VTemp(I),I=6,ITop)
          If(NDOF.gt.10) Write(IOut,1700) (VTemp(I),I=11,NDOF)
          endIf
        endIf
C
C     COMPUTE
C        EVIB--- THE VIBRATIONAL COMPONENT OF THE INTERNAL ENERGY.
C        CVIB--- THE VIBRATIONAL COMPONENT OF THE HEAT CAPACITY.
C        SVIB--- THE VIBRATIONAL COMPONENT OF THE ENTROPY.
C        QVIB--- THE VIBRATIONAL COMPONENT OF THE PARTITION FUNCTION.
C        QLVIB-- THE VIBRATIONAL COMPONENT OF LOG10 OF Q.
C        QLNVIB- THE VIBRATIONAL COMPONENT OF LN    OF Q.
C
      TWOT = T + T
      EVIB = ZERO
      CVIB = ZERO
      SVIB = ZERO
      QLVIB = ZERO
      QLNVIB = ZERO
      QZLVIB = ZERO
      QZLNVI = ZERO
      DO 40 I=1,NDOF
        TOVT  = VTEMP(I) / T
        TOV2T = VTEMP(I) / TWOT
C
C       COMPUTE CONTRIBUTIONS DUE TO THE I'TH VIBRATION.
C       FOR THE LOW FREQUENCY MODES THESE ARE STORED.
C       FOR ALL MODES THEY ARE ADDED INTO THE TOTAL.
C
C       WARNING!! THE FOLLOWING SLEAZE ONLY WORKS BECAUSE WE JUST DIVIDE
C       BY ETOVT AND (ETOVT-1).  IF THESE NUMBERS WERE USED FOR OTHER
C       PURPOSE MORE CARE WOULD BE NEEDED.  THE CUTOFF IS REQUIRED TO
C       AVOID OVERFLOW ON THE VAX.
        If(TOvT.gt.Small) then  
          ETOvT = Exp(Min(TOvT,MXTOVT))
          ONEME = ONE - ONE/ETOVT
          EM1   = ETOVT - ONE
          ECONT = TOVT  *  (HALF + ONE/EM1)
          CCONT = ETOVT *  (TOVT/EM1)**2
          SCONT = TOVT/EM1 - Log(ONEME)
          QCONT = EXP(-TOV2T) / ONEME
          QLCONT = Log10(QCONT)
          QLNCON = Log(QCONT)
          QZCONT = ONE / ONEME
          QZLCON = Log10(QZCONT)
          QZLNCO = Log(QZCONT)
        else
          ECont = Zero
          CCont = Zero
          SCont = One
          QCont = One
          QLCONT = Zero
          QLNCON = Zero
          QZCONT = One
          QZLCON = Zero
          QZLNCO = Zero
          endIf
        If(I.le.LoFreq) then
          EVIBN(I) = ECONT * RT
          CVIBN(I) = CCONT * GAS
          SVIBN(I) = SCONT * GAS
          QVIBN(I) = QCONT
          QZVIBN(I) = QZCONT
          endIf
        EVIB = EVIB + ECONT
        CVIB = CVIB + CCONT
        SVIB = SVIB + SCONT
        QLVIB = QLVIB + QLCONT
        QLNVIB = QLNVIB + QLNCON
        QZLVIB = QZLVIB + QZLCON
   40   QZLNVI = QZLNVI + QZLNCO
      EVIB = EVIB * RT
      CVIB = CVIB * GAS
      SVIB = SVIB * GAS
      QVib = Exp(QLnVib)
      QZVib = Exp(QZLnVi)
C
C     THE UNITS ARE NOW
C         E-- JOULES/MOL
C         C-- JOULES/MOL-KELVIN
C         S-- JOULES/MOL-KELVIN
C
      Return
      End
*Deck TMAtom
      Logical Function TMAtom(IAnI)
      Implicit Integer(A-Z)
C
C     Return whether the specified atom is an active transition metal
C     atom.
C
      TMAtom = (IAnI.ge.21.and.IAnI.le.30).or.
     $  (IAnI.ge.39.and.IAnI.le.48).or.
     $  (IAnI.ge.57.and.IAnI.le.80).or.(IAnI.ge.89.and.IAnI.le.112)
      Return
      End
*Deck TetRnd
      Logical Function TetRnd(Angle)
      Implicit Real*8(A-H,O-Z)
C
C     Test Angle for being near to tetrahedral and round it to the
C     exact value.  The value of this function is .true. if rounding
C     was performed.
C
      Save One, Three, F45, Tol
      Data One/1.0d0/, Three/3.0d0/, F45/45.0d0/, Tol/1.d-3/
C
      RefAng = ACos(-One/Three)
      Tol1 = Tol * ATan(One) / F45
      If(Abs(Angle-RefAng).lt.Tol1) then
        TetRnd = .True.
        Angle = RefAng
      else if(Abs(Angle+RefAng).lt.Tol1) then
        TetRnd = .True.
        Angle = -RefAng
      else
        TetRnd = .False.
        endIf
      Return
      End
*Deck Transl
      Subroutine Transl(ISign,NX,NVec,TrVect,X)
      Implicit Real*8(A-H,O-Z)
C
C     Translate the NVec vectors in X using TrVect.
C
      Dimension TrVect(NX), X(NX,NVec)
C
      If(ISign.ge.0) then
        Do 10 IX = 1, NX
          Do 10 IV = 1, NVec
   10       X(IX,IV) = X(IX,IV) + TrVect(IX)
      else
        Do 20 IX = 1, NX
          Do 20 IV = 1, NVec
   20       X(IX,IV) = X(IX,IV) - TrVect(IX)
        endIf
      Return
      End
*Deck Trspn
      Subroutine Trspn(NDim,N,T)
      Implicit Real*8(A-H,O-Z)
C
C     Transpose square matrix T.
C
      Dimension T(NDim,NDim)
C
      Do 10 I = 2, N
        Do 10 J = 1, (I-1)
          X = T(I,J)
          T(I,J) = T(J,I)
   10     T(J,I) = X
      Return
      End
*Deck UDiff
      Function UDiff(N,V,U)
      Implicit Real*8(A-H,O-Z)
C
C     Norm of the difference two unitary matrices -- |V(t).U-I|.
C
      Dimension V(N,N), U(N,N)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      Sum = Zero
      Do 30 I = 1, N
        Do 10 J = 1, (I-1)
   10     Sum = Sum + SProd(N,V(1,I),U(1,J))**2
        Sum = Sum + (SProd(N,V(1,I),U(1,I))-One)**2
        Do 20 J = (I+1), N
   20     Sum = Sum + SProd(N,V(1,I),U(1,J))**2
   30   Continue
      UDiff = Sqrt(Sum)
      Return
      End
*Deck VEC
      SUBROUTINE VEC(SMALL,OHOH,U,C,J,K)
      Implicit Real*8(A-H,O-Z)
      LOGICAL OHOH
      DIMENSION C(1),R(3),U(3)
      Save ZERO
      DATA ZERO/0.0D0/
C
      R2=ZERO
      JTEMP=(J-1)*3
      KTEMP=(K-1)*3
      DO 10 I=1,3
      R(I)=C(I+JTEMP)-C(I+KTEMP)
   10 R2=R2+R(I)*R(I)
      R2=Sqrt(R2)
      OHOH = R2 .LT. SMALL
      IF (OHOH) RETURN
      DO 20 I=1,3
   20 U(I)=R(I)/R2
      RETURN
      END
*Deck VibPrt
      Subroutine VibPrt(IOut,IPrint,NAtoms,NVib,NImag,Linear,
     $  FcScal,PUser,TUser,PhyCon,QVib,QZVib,SVib,EVib,CVib,
     $  Freq,Scr)
      Implicit Real*8 (A-H,O-Z)
      Logical Linear
      Dimension PhyCon(*),Freq(*),Scr(*)
      Save zero,one
      Data zero,one /0.0d0,1.0d0/
      QVib  = One
      QZVib = One
      EVib  = Zero
      CVib  = Zero
      SVib  = Zero
      FCScal= One
      If(NVib.eq.0) Return
      NImag = 0
      do 10 IV=1,NVib
       If(Freq(IV).lt.zero) NImag=NImag+1
       Freq(IV) = Freq(IV)*FCScal
   10 continue
      Write(IOut,'(5X,I5,'' Vibrations scaled by'',F10.5)') NVib,FCScal
      Write(IOut,'(6F12.2)') (Freq(I),I=1,NVib)
      call ThrVib(IOut,IPrint,NVib,PhyCon,PUser,TUser,Freq,Scr,QVib,
     $  QZVib,SVib,EVib,CVib)
      Return
      End   
*Deck VProd
      Subroutine VProd(VP,X,Y)
      Implicit Real*8(A-H,O-Z)
C
C     VP = X Cross Y
C
      Dimension VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      Return
      End
*Deck ZATan2
      Function ZATan2(Y,X)
      Implicit Real*8(A-H,O-Z)
C
C     This silly little function is necessary because the fortran ATan2
C     function is undefined for two zero arguments.  This routine returns
C     zero in that case, which is usually what a caller who wants to do
C     nothing in that case requires.
C
      Save Zero, One, Two
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/
C
      If(X.eq.Zero.and.Y.eq.Zero) then
        ZATan2 = Zero
      else if(X.eq.Zero) then
        ZATan2 = Sign(Two*ATan(One),Y)
      else if(Y.eq.Zero) then
        ZATan2 = Zero
      else
        ZATan2 = ATan2(Y,X)
        endIf
      Return
      End
