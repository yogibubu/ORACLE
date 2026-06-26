*Deck GenCB1
      Subroutine GenCB1(IOut,IAlg,Levels,ToAng,DLim,NBox,IBoxSt,IBox,
     $  IndTab,IP2Box,IAn,C,MxBond,NBond,IBond,IBType)
      Implicit Real*8(A-H,O-Z)
C
C     Build the basic, unsorted connectivity table for GenCBx.  The
C     boxes are big enough that there can only be bonds between adjacent
C     boxes.  So we loop over boxes and consider only atom pairs
C     involving that box and 27 others.
C
      Integer BinSeI
      Dimension IBoxSt(*), IBox(*), IndTab(*), IP2Box(4,*), IAn(*),
     $  C(3,*), NBond(*), IBond(MxBond,*), IBType(MxBond,*)
 1000 Format(' Maximum number of bonds=',I3,' exceeded for atom',I4,'.')
C
      DLimSq = DLim**2
      NBX = 2**Levels
      Do 60 IB = 1, NBox
        IBSt = IBoxSt(IB)
        IBEnd = IBoxSt(IB+1) - 1
        IBNum = IBox(IB)
        I1 = IndTab(IBSt)
        IBX = IP2Box(1,I1)
        IBY = IP2Box(2,I1)
        IBZ = IP2Box(3,I1)
        Do 50 IXI = -1, 1
          JBX = IBX + IXI
          If(JBX.lt.1.or.JBX.gt.NBX) goto 50
          Do 40 IYI = -1, 1
            JBY = IBY + IYI
            If(JBY.lt.1.or.JBY.gt.NBX) goto 40
            Do 30 IZI = -1, 1
              JBZ = IBZ + IZI
              If(JBZ.lt.1.or.JBZ.gt.NBX) goto 30
              JBNum = IndSrt(JBZ-1,JBY-1,JBX-1)
              If(JBNum.gt.IBNum) goto 30
              JB = BinSeI(0,JBNum,NBox,IBox)
              If(JB.eq.0) goto 30
              JBSt = IBoxSt(JB)
              JBEnd = IBoxSt(JB+1) - 1
              Do 20 IX = IBSt, IBEnd
                I = IndTab(IX)
                XI = C(1,I)
                YI = C(2,I)
                ZI = C(3,I)
                If(IB.eq.JB) then
                  LimJ = IX - 1
                else
                  LimJ = JBEnd
                  endIf
                Do 10 JX = JBSt, LimJ
                  J = IndTab(JX)
                  DX = C(1,J) - XI
                  DY = C(2,J) - YI
                  DZ = C(3,J) - ZI
                  RSq = DX*DX + DY*DY + DZ*DZ
                  If(RSq.le.DLimSq) then
                    RIJ = Sqrt(RSq)
                    IBondO = IRToBO(ToAng,IAn(I),IAn(J),RIJ,IAlg)
                    If(IBondO.ne.0) then
                      NBond(I) = NBond(I) + 1
                      If(NBond(I).gt.MxBond) then
                        Write(IOut,1000) MxBond, I
                        Call Lnk1E(0)
                        endIf
                      IBond(NBond(I),I) = J
                      IBType(NBond(I),I) = IBondO
                      NBond(J) = NBond(J) + 1
                      If(NBond(J).gt.MxBond) then
                        Write(IOut,1000) MxBond, J
                        Call Lnk1E(0)
                        endIf
                      IBond(NBond(J),J) = I
                      IBType(NBond(J),J) = IBondO
                      endIf
                    endIf
   10             Continue
   20           Continue
   30         Continue
   40       Continue
   50     Continue
   60   Continue
      Return
      End
*Deck GenCBx
      Subroutine GenCBx(IOut,IPrint,IAlg,ToAng,DLimI,MxBond,NAtoms,IAn,
     $  C,NBond,IBond,IBType,IP2Box,IndTab,IndBox,IBox,IBoxSt)
      Implicit Real*8(A-H,O-Z)
C
C     Generate connectivity based on bond distances alone.  The criteria
C     are contained in routine IRToBO.  IAlg=3 applies additional criteria
C     for UFF mechanics.  IAlg=-1 means nothing will be bonded.
C
      Dimension IAn(*), C(3,*), NBond(*), IBond(MxBond,*), XX(1),
     $  IBType(MxBond,*), IP2Box(4,*), IndTab(*), IndBox(*), IBox(*),
     $  IBoxSt(*), Shift(3)
      Save XX, DLimD, Zero
      Data XX/0.0d0/, DLimD/8.0d0/, Zero/0.0d0/
C
      IAlg0 = Mod(IAlg,10)
      Len = MxBond * NAtoms
      Call IClear(NAtoms,NBond)
      Call IClear(Len,IBond)
      Call IClear(Len,IBType)
      DLim = DLimI
      If(DLim.eq.Zero) DLim = DLimD
      Scaler = DLim
      Call SMxMin(NAtoms,C,SMaxX,Shift)
      Call SMx2Lv(0,SMaxX,Scaler,Levels)
      Call SrtPnt(NAtoms,C,Levels,Scaler,SMaxX,Shift,IP2Box,IndTab,
     $  IndBox)
      Call PntBox(NAtoms,IndTab,IndBox,NBox,IBox,IBoxSt)
      Call PrtBox(IOut,IPrint,Levels,Scaler,SMaxX,Shift,NBox,IBox,
     $  IBoxSt,IndTab,IP2Box)
      Call GenCB1(IOut,IAlg0,Levels,ToAng,DLim,NBox,IBoxSt,IBox,IndTab,
     $  IP2Box,IAn,C,MxBond,NBond,IBond,IBType)
      Call SrtCon(.False.,MxBond,NAtoms,NBond,IBond,IBType,XX,IndTab,
     $  IndBox,XX,IBox)
      If(IPrint.ge.2) Call PrtCon(IOut,'GenCBx connectivity:',.False.,
     $  NAtoms,MxBond,IAn,NBond,IBond,IBType,XX,0)
      Return
      End
*Deck IRToBO
      Function IRToBO(ToAng,IA,JA,RIJ,IAlg)
      Implicit Real*8(A-H,O-Z)
C
C     Return the order of the bond between atoms of atomic numbers IA
C     and JA, separated by RIJ, or 0 if there is no bond.
C
C     Generate connectivity based on bond distances alone.  The criteria
C     is whether the distances is no more than 30% longer than the 6-31g*
C     single bond length.  For the moment all bond types are determined
C     using Pauling bond orders.  For bondings involving atoms
C     from Ar - At, atomic radii are used to determine the single bond
C     length.  For atoms later than At, an arbitrary maximum bond length
C     of 3.5 Angstroms is used.
C
C     Attempt to identify resonance bonds for mechanics (IAlg=3)
C
      Save Zero, Half, TolA, RCC, RCN, RCO, RNN
      Data Zero/0.0d0/, Half/0.5d0/, TolA/0.03d0/, RCC/1.4d0/,
     $  RCN/1.34d0/, RCO/1.36d0/, RNN/1.34d0/
C
      BondOr = RToBO(ToAng,IA,JA,RIJ,IAlg)
      If(BondOr.gt.Zero) then
        IRToBO = Min(Max(IGFix(BondOr+Half),1),3)
        If(IAlg.eq.3) then
          RIJA = ToAng * RIJ
C         resonant C-C bond
          If(IA.eq.6.and.JA.eq.6.and.Abs(RIJA-RCC).lt.TolA) then
            IRToBO = -1
C         resonant N-N bond
          else if(IA.eq.7.and.JA.eq.7.and.Abs(RIJA-RNN).lt.TolA) then
            IRToBO = -1
C         resonant C-N bond
          else if(((IA.eq.6.and.JA.eq.7).or.(IA.eq.7.and.JA.eq.6)).and.
     $      Abs(RIJA-RCN).lt.TolA) then
            IRToBO = -1
C         resonant C-O bond
          else if(((IA.eq.6.and.JA.eq.8).or.(IA.eq.8.and.JA.eq.6)).and.
     $      Abs(RIJA-RCO).lt.TolA) then
            IRToBO = -1
            endIf
          endIf
      else
        IRToBO = 0
        endIf
      Return
      End
*Deck HpSort
      Subroutine HpSort(N,Ind,A)
      Implicit Integer(A-Z)
C
C     Sort integer array A using HeapSort, returning the new order in Ind.
C
      Integer Ind(N),A(N)
C
      L = N/2 + 1
      R = N
      Do 10 I = 1, N
      Ind(I) = I
   10 continue
      If(N.le.1) Return
C
   20 If(L.gt.1) then
        L = L - 1
        IndCur = Ind(L)
        K = A(IndCur)
      else
        IndCur = Ind(R)
        K = A(IndCur)
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
          If(K.lt.A(Ind(J))) then
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
*Deck PntBox
      Subroutine PntBox(N,IndTab,IndBox,NBox,IBox,IBoxSt)
      Implicit Integer(A-Z)
C
C Generate a table of where the points in each non-empty box can be
C     found in the sorted list.
C
      Dimension IndTab(*), IndBox(*), IBox(*), IBoxSt(*)
C
      NBox = 0
      Last = -1
      Do 10 I = 1, N
        IX = IndTab(I)
        IS = IndBox(IX)
        If(IS.ne.Last) then
          NBox = NBox + 1
          Last = IS
          IBox(NBox) = Last
          IBoxSt(NBox) = I
          endIf
   10   Continue
      IBox(NBox+1) = 0
      IBoxSt(NBox+1) = N + 1
      Return
      End
*Deck PrtBox
      Subroutine PrtBox(IOut,IPrint,Levels,Scaler,SMaxX,Shift,NBox,IBox,
     $  IBoxSt,IndTab,IP2Box)
      Implicit Real*8(A-H,O-Z)
C
C     Print out boxified center data.
C
      Parameter (MinPrt=1)
      Dimension Shift(3), IBox(*), IBoxSt(*), IndTab(*), IP2Box(4,*)
 1000 Format(' PrtBox:  NBox=',I6,' Levels=',I3,' BoxLen=',F8.2,
     $  ' SMaxX=',F12.2,/,'          Shift=',3F20.6)
 1010 Format(' Box',I7,' Number',I20,' centers from',I7,' to',I7,':')
 1020 Format('    ',7X,' I=',I6,' IX=',I6,' XYZ=',3I7,' Ind=',I20)
C
      If(IPrint.lt.MinPrt) Return
      Write(IOut,1000) NBox, Levels, Scaler, SMaxX, Shift
      Do 20 IB = 1, NBox
        IBSt = IBoxSt(IB)
        IBEnd = IBoxSt(IB+1) - 1
        Write(IOut,1010) IB, IBox(IB), IBSt, IBEnd
        If(IPrint.gt.MinPrt) then
          Do 10 I = IBSt, IBEnd
            IX = IndTab(I)
            IT = IndSrt(IP2Box(3,IX)-1,IP2Box(2,IX)-1,IP2Box(1,IX)-1)
            Write(IOut,1020) I, IX, (IP2Box(J,IX),J=1,3), IT
   10       Continue
          endIf
   20   Continue
      Return
      End
*Deck PrtCon
      Subroutine PrtCon(IOut,Label,RealCn,NAtoms,MxBond,IAn,NBond,IBond,
     $  IBType,RBType,NDo)
      Implicit Real*8(A-H,O-Z)
C
C     Print out connectivity information.
C     NDo > 1 : Print out connectivities 1 to NDo
C     NDo < 0 : Print out only connectivity -NDo
C
      Parameter (MaxL1=5)
      Character*(*) Label
      Logical RealCn
      Dimension IAn(*), NBond(NAtoms,*), IBond(MxBond,NAtoms,*),
     $  RBType(MxBond,NAtoms,*), IBType(MxBond,NAtoms,*)
 1000 Format(1X,A)
 1005 Format(1X,A,I3)
 1010 Format(' I=',I6,' IAn=',I3,' IBond, RBType=',5(I6,F6.2))
 1020 Format(('   ',6X,'     ',3X,'               ',5(I6,F6.2)))
 1030 Format(' I=',I6,' IAn=',I3,' IBond, IBType=',10I6)
 1040 Format(('   ',6X,'     ',3X,'               ',10I6))
C
      If(NDo.eq.0) then
        NStart = 1
        NEnd = 1
      else if(NDo.gt.0) then
        NStart = 1
        NEnd = NDo
      else
        NStart = -NDo
        NEnd = -NDo
        endif
      Do 30 II = NStart, NEnd
        L = LinEnd(Label)
        If(L.gt.0.and.NDo.eq.0) Write(IOut,1000) Label(1:L)
        If(L.gt.0.and.NDo.ne.0) Write(IOut,1005) Label(1:L), II
        If(RealCn) then
          Do 10 I = 1, NAtoms
            Write(IOut,1010) I, IAn(I),
     $        (IBond(J,I,II),RBType(J,I,II),J=1,Min(NBond(I,II),MaxL1))
            If(NBond(I,II).gt.MaxL1) Write(IOut,1020)
     $        (IBond(J,I,II),RBType(J,I,II),J=(MaxL1+1),NBond(I,II))
   10       Continue
        else
          Do 20 I = 1, NAtoms
            Write(IOut,1030) I, IAn(I),
     $        (IBond(J,I,II),IBType(J,I,II),J=1,Min(NBond(I,II),MaxL1))
            If(NBond(I,II).gt.MaxL1) Write(IOut,1040)
     $        (IBond(J,I,II),IBType(J,I,II),J=(MaxL1+1),NBond(I,II))
   20       Continue
          endIf
   30   continue
      Return
      End
*Deck RadixS
      Subroutine RadixS(N,SrtInd,BckInd,Npass)
      Implicit Integer(a-z)
C
C     The subroutine radix-sorts SrtInd.  BckInd will show what was the
C     original order.  Npass shows how many passes (1,2,3)
C
      Parameter (iBuck=2**8)
      Integer SrtInd(N),BckInd(N,2),Buck(0:iBuck-1,4)
C
      shft1 = iBuck
      shft2 = 1
      sw1 = 1
      sw2 = sw1 + 1
      sw3 = 3
      st1 = 1
      st2 = 2
C     Here some setup work
      Call Iclear(iBuck*2,Buck(0,sw3))
      Buck(0,sw3) = 1
      Do 10 jn = 1, N
       BckInd(jn,st2) = jn + 1
   10 continue
      BckInd(N,st2) = 0
C     This is the main part
      Do 40 I = 1, NPass
        Call IClear(iBuck*2,Buck(0,sw1))
        Do 30 lb = 0, (iBuck-1)
          J = Buck(lb,sw3)
          If(J.ne.0) then
   20       A = SrtInd(j)
C            B1 = (A-((A/shft1)*shft1))/shft2
            B1 = mod(A,shft1)/shft2
            B = B1
C           sw1 corresponds to first element
C           sw2 corresponds to last element
            ICur = Buck(B,sw2)
            Buck(B,sw2) = J
            BckInd(j,st1) = 0
            If(ICur.eq.0) then
              Buck(B,sw1) = J
            else
              BckInd(ICur,st1) = J
              endIf
            J = BckInd(j,st2)
            If(J.ne.0) goto 20
            endIf
   30     Continue
C       Here we just switch indices for next pass
        t = sw1
        sw1 = sw3
        sw2 = sw1 + 1
        sw3 = t
        t = st2
        st2 = st1
        st1 = t
        shft1 = shft1*iBuck
        shft2 = shft2*iBuck
        If(Shft1.le.Shft2) Shft1 = MaxInt(0)
   40  Continue
C
C     Now we can collect the results
      itot = 1
      Do 60 lb = 0, (iBuck-1)
        j = Buck(lb,sw3)
        If(J.ne.0) then
   50     BckInd(itot,st1) = J
          itot = itot + 1
          J = BckInd(j,st2)
          If(J.ne.0) goto 50
          endIf
   60   Continue
C     It would be nice to return the results in BckInd(1,1)
      If(st1.ne.1) Call IMove(N,BckInd(1,st1),BckInd(1,1))
      Return
      End
*Deck RToBO
      Function RToBO(ToAng,IA,JA,RIJ,IAlg)
      Implicit Real*8(A-H,O-Z)
C
C     Return the order of the bond between atoms of atomic numbers IA
C     and JA, separated by RIJ, or 0 if there is no bond.
C
C     IAlg = 0,1 ... Default choice of Radii (in order of precedence
C                    tabulated 6-31G* values, model B, pcm cov. radii, 3.5,
C                    with an upper limit 30% larger
C              2 ... Same as 1, but upper limit 20% larger.
C              3 ... Truhlar's radii and bond order.
C             -1 ... Bond order 0 (nothing bonded).
C
C     Generate connectivity based on bond distances alone.  The criteria
C     is whether the distances is no more than 30% longer than the 6-31g*
C     single bond length.  For the moment all bond types are determined
C     using Pauling bond orders.  For bondings involving atoms
C     from Ar - At, atomic radii are used to determine the single bond
C     length.  For atoms later than At, an arbitrary maximum bond length
C     of 3.5 Angstroms is used.
C
      Common/IO/In,IOut,Ipunch
      Save Zero, One, Extnd1, Extnd2, A, DefBnd, CTAlp
      Data Zero/0.0d0/, One/1.0d0/, Extnd1/1.3d0/, Extnd2/1.2d0/,
     $  A/0.3d0/, DefBnd/3.5d0/, CTAlp/2.4740d0/
C
      ToBohr = One / ToAng
      If(IAlg.eq.0.or.IAlg.eq.1.or.IAlg.eq.2) then
        RRef = SBondL(5,IA,JA)
        If(RRef.eq.Zero) RRef = BLModB(IA,JA)
        If(RRef.eq.Zero) RRef = RCovA(IA,JA)
        If(RRef.eq.Zero) RRef = DefBnd
        RRef = ToBohr*RRef
        Extend = Extnd1
        If(IAlg.eq.2) Extend = Extnd2
        RMax = Extend * RRef
        If(RIJ.le.RMax) then
          RToBO = Exp((RRef-RIJ)/(A*ToBohr))
        else
          RToBO = Zero
          endIf
      else if(IAlg.eq.3) then
        RRef = RCovCT(IA,JA)
        RIJA = RIJ*ToAng
        RToBO = Exp(-CTAlp*(RIJA-RRef))
      else if(IAlg.eq.-1) then
        RToBO = Zero
      else
        RToBO = Zero
        write(IOut,'('' Illegal IAlg in RToBO'')')
        Stop
        endIf
      Return
      End
*Deck SBondL
      Function SBondL(IType,IA,JA)
      Implicit Real*8(A-H,O-Z)
C
C     Return a standard single bond length between atoms types IA and
C     JA.  The appropriate saturated two-heavy atom or hydride molecule
C     is used (i.e., Li-H from LiH, Li-C from CH3Li).  IType is the
C     method to use:
C
C     IType = 1 ... Experimental (NYI).
C     IType = 2 ... MNDO/AM1 (NYI).
C     IType = 3 ... RHF/STO-3G (NYI).
C     IType = 4 ... RHF/3-21G (NYI).
C     IType = 5 ... RHF/6-31G*.
C
C     A value of 0.0 is returned if the distance is not available.
C
C     RRef is dimensioned for MaxAn pairs of atoms and MaxTyp types.
C     The DRefxx are just used to avoid exceessive continuation lines
C     in the data statement.
      Dimension RRef(153,5), DRef1A(70), DRef1B(83), DRef2A(70),
     $    DRef2B(83), DRef3A(70), DRef3B(83), DRef4A(70), DRef4B(83),
     $    DRef5A(70), DRef5B(83)
      Equivalence (RRef(1,1),DRef1A(1)), (RRef(71,1),DRef1B(1))
      Equivalence (RRef(1,2),DRef2A(1)), (RRef(71,2),DRef2B(1))
      Equivalence (RRef(1,3),DRef3A(1)), (RRef(71,3),DRef3B(1))
      Equivalence (RRef(1,4),DRef4A(1)), (RRef(71,4),DRef4B(1))
      Equivalence (RRef(1,5),DRef5A(1)), (RRef(71,5),DRef5B(1))
      Save MaxTyp, MaxAn, Zero, RRef
      Data MaxTyp/5/, MaxAn/17/, Zero/0.0d0/
      Data DRef1A/70*0.0d0/, DRef1B/83*0.0d0/, DRef2A/70*0.0d0/,
     $    DRef2B/83*0.0d0/, DRef3A/70*0.0d0/, DRef3B/83*0.0d0/,
     $    DRef4A/70*0.0d0/, DRef4B/83*0.0d0/
      Data DRef5A/
C
C                H-H    He-H    He-He   Li-H    Li-He   Li-Li   Be-H
     $          0.730d0,0.000d0,0.000d0,1.636d0,0.000d0,2.812d0,1.334d0,
C                Be-He   Be-Li   Be-Be    B-H    B-He    B-Li    B-Be
     $          0.000d0,2.469d0,2.123d0,1.188d0,0.000d0,2.234d0,1.903d0,
C                 B-B     C-H    C-He    C-Li    C-Be     C-B     C-C
     $          1.679d0,1.084d0,0.000d0,2.001d0,1.699d0,1.574d0,1.528d0,
C                 N-H    N-He    N-Li    N-Be     N-B     N-C     N-N
     $          1.002d0,0.000d0,1.750d0,1.503d0,1.389d0,1.453d0,1.413d0,
C                 O-H    O-He    O-Li    O-Be     O-B     O-C     O-N
     $          0.947d0,0.000d0,1.592d0,1.377d0,1.344d0,1.400d0,1.404d0,
C                 O-O     F-H    F-He    F-Li    F-Be     F-B     F-C
     $          1.396d0,0.911d0,0.000d0,1.555d0,1.366d0,1.313d0,1.365d0,
C                 F-N     F-O     F-F    Ne-H    Ne-He   Ne-Li   Ne-Be
     $          1.386d0,1.376d0,1.345d0,0.000d0,0.000d0,0.000d0,0.000d0,
C                Ne-B    Ne-C    Ne-N    Ne-O    Ne-F    Ne-Ne   Na-H
     $          0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,1.914d0,
C                Na-He   Na-Li   Na-Be   Na-B    Na-C    Na-N    Na-O
     $          0.000d0,2.999d0,2.742d0,2.537d0,2.324d0,2.080d0,1.921d0,
C                Na-F    Na-Ne   Na-Na   Mg-H    Mg-He   Mg-Li   Mg-Be
     $          1.885d0,0.000d0,3.189d0,1.718d0,0.000d0,2.847d0,2.531d0/
      Data DRef5B/
C                Mg-B    Mg-C    Mg-N    Mg-O    Mg-F    Mg-Ne   Mg-Na
     $          2.319d0,2.106d0,1.894d0,1.756d0,1.730d0,0.000d0,3.087d0,
C                Mg-Mg   Al-H    Al-He   Al-Li   Al-Be   Al-B    Al-C
     $          2.912d0,1.584d0,0.000d0,2.693d0,2.373d0,2.151d0,1.972d0,
C                Al-N    Al-O    Al-F    Al-Ne   Al-Na   Al-Mg   Al-Al
     $          1.771d0,1.697d0,1.640d0,0.000d0,2.963d0,2.771d0,2.613d0,
C                Si-H    Si-He   Si-Li   Si-Be   Si-B    Si-C    Si-N
     $          1.475d0,0.000d0,2.524d0,2.217d0,2.040d0,1.888d0,1.724d0,
C                Si-O    Si-F    Si-Ne   Si-Na   Si-Mg   Si-Al   Si-Si
     $          1.647d0,1.594d0,0.000d0,2.815d0,2.617d0,2.478d0,2.352d0,
C                 P-H    P-He    P-Li    P-Be     P-B     P-C     P-N
     $          1.403d0,0.000d0,2.375d0,2.075d0,1.902d0,1.860d0,1.706d0,
C                 P-O     P-F    P-Ne    P-Na    P-Mg    P-Al    P-Si
     $          1.650d0,1.599d0,0.000d0,2.683d0,2.479d0,2.341d0,2.266d0,
C                 P-P     S-H    S-He    S-Li    S-Be     S-B     S-C
     $          2.214d0,1.326d0,0.000d0,2.191d0,1.917d0,1.791d0,1.818d0,
C                 S-N     S-O     S-F    S-Ne    S-Na    S-Mg    S-Al
     $          1.695d0,1.654d0,1.612d0,0.000d0,2.515d0,2.317d0,2.196d0,
C                S-Si     S-P     S-S    Cl-H    Cl-He   Cl-Li   Cl-Be
     $          2.151d0,2.127d0,2.063d0,1.266d0,0.000d0,2.072d0,1.810d0,
C                Cl-B    Cl-C    Cl-N    Cl-O    Cl-F    Cl-Ne   Cl-Na
     $          1.754d0,1.785d0,1.732d0,1.670d0,1.613d0,0.000d0,2.397d0,
C                Cl-Mg   Cl-Al   Cl-Si   Cl-P    Cl-S    Cl-Cl
     $          2.211d0,2.111d0,2.068d0,2.072d0,2.034d0,1.990d0/
      LInd(I,J) = ((Max(I,J)*(Max(I,J)-1)/2)+Min(I,J))
C
      Dist = Zero
      If(IType.lt.1.or.IType.gt.MaxTyp.or.IA.lt.1.or.IA.gt.MaxAn.or.
     $   JA.lt.1.or.JA.gt.MaxAn) goto 999
      IJ = LInd(IA,JA)
      Dist = RRef(IJ,IType)
  999 SBondL = Dist
      Return
      End
*Deck SMxMin
      Subroutine SMxMin(N,C,SMaxX,Shift)
      Implicit Real*8(A-H,O-Z)
C
C     Generate the shift and max dimensions required to boxify the
C     points in C.
C
      Dimension C(3,*), Shift(3), SMin(3), SMax(3)
      Save OffSet
      Data OffSet/1.d-8/
C
      Do 20 IX = 1, 3
       SMax(IX) = C(IX,1)
       SMin(IX) = C(IX,1)
       Do 10 IAt = 2, N
        SMax(IX) = Max(SMax(IX),C(IX,IAt))
        SMin(IX) = Min(SMin(IX),C(IX,IAt))
   10  Continue
       SMin(IX) = SMin(IX) - OffSet
   20 Continue
      SMaxX = Max(SMax(1)-SMin(1),SMax(2)-SMin(2),SMax(3)-SMin(3))
      Call ANeg(3,SMin,Shift)
      Return
      End
*Deck SMx2Lv
      Subroutine SMx2Lv(LevIn,SMaxX,BoxLen,Levels)
      Implicit Real*8(A-H,O-Z)
C
C     Compute the number of levels given the box size and dimension of
C     the cell.
C
      Save Pt5, Two, SInc
      Data Pt5/0.5d0/, Two/2.0d0/, SInc/1.0d0/
C
      RNBox = (SMaxX+SInc) / BoxLen
      RLev = GLog(RNBox)/GLog(Two)
      Levels = Max(LevIn,IGFix(RLev+Pt5)+1)
      Return
      End
*Deck SrtCon
      Subroutine SrtCon(UseRB,MxBond,NAtoms,NBond,IBond,IBType,RBType,
     $  IBI,IBTI,RBTI,Ind)
      Implicit Real*8(A-H,O-Z)
C
C     Sort the connectivity table into increasing order in each row.
C
      Logical UseRB
      Dimension NBond(*), IBond(MxBond,*), IBType(MxBond,*),
     $  RBType(MxBond,*), IBI(*), IBTI(*), RBTI(*), Ind(*)
C
      If(UseRB) then
        Do 20 I = 1, NAtoms
         If(NBond(I).gt.1) then
          Call IMove(NBond(I),IBond(1,I),IBI)
          Call AMove(NBond(I),RBType(1,I),RBTI)
          Call HpSort(NBond(I),Ind,IBI)
          Do 10 J = 1, NBond(I)
           J1 = Ind(J)
           IBond(J,I) = IBI(J1)
           RBType(J,I) = RBTI(J1)
   10     Continue 
         endIf
   20   Continue
      else
        Do 40 I = 1, NAtoms
         If(NBond(I).gt.1) then
          Call IMove(NBond(I),IBond(1,I),IBI)
          Call IMove(NBond(I),IBType(1,I),IBTI)
          Call HpSort(NBond(I),Ind,IBI)
          Do 30 J = 1, NBond(I)
           J1 = Ind(J)
           IBond(J,I) = IBI(J1)
           IBType(J,I) = IBTI(J1)
   30     continue
         endIf
   40   Continue
      endIf
      Return
      End
*Deck SrtPnt
      Subroutine SrtPnt(NPnt,XYZ,Levels,Scaler,SMaxX,Shift,ip2Box,
     $  IndTab,SrtInd)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine sorts distributions and prepares
C     an index array for further manipulations
C
      Integer SrtInd(NPnt),ip2Box(4,NPnt),IndTab(NPnt,2),xBox,yBox,zBox
      Real*8 XYZ(3,NPnt),Shift(3)
      Save One,Pjunk
      Data Pjunk/1.234567898d-8/,One/1.0d0/
C
C     Determine if the system fits in the box with levels and scale1
C     If it doesn't, first levels will be sacrificed (changed), then
C     if the system is too large, scaler will be increased
      MxBit = MDNBPW(0)/3
C     Increase levels till the scaler is correct
   15 fNb = float(2**levels)
      If((fNb*scaler-Pjunk).lt.SMaxX) then
        Levels = Levels + 1
        Goto 15
        endIf
C     If there are too many levels, decrease scaler
      If(Levels.gt.MxBit) then
        Levels = MxBit
        fNb = float(2**levels)
        scaler = (SMaxX+Pjunk)/fNb
        endIf
      NPass = (3*Levels+7)/8
C
C     Boxify point charges
      Scale1 = One/Scaler
      Do 20 i=1,NPnt
        xBox = Int((XYZ(1,i)+Shift(1))*Scale1)
        yBox = Int((XYZ(2,i)+Shift(2))*Scale1)
        zBox = Int((XYZ(3,i)+Shift(3))*Scale1)
        ip2Box(1,i) = xBox + 1
        ip2Box(2,i) = yBox + 1
        ip2Box(3,i) = zBox + 1
        ip2Box(4,i) = 1
        SrtInd(i) = IndSrt(zBox,yBox,xBox)
   20   Continue
C
C     Now it is time to sort the distributions
C     Note: IndTab is dimensioned as IndTab(Nmax,2)
C
      Call RadixS(NPnt,SrtInd,IndTab,NPass)
      Return
      End
*Deck BinSeI
      Integer Function BinSeI(IRet,IVal,N,List)
      Implicit Real*8(A-H,O-Z)
C
C     Binary search for IVal in sorted List.  If IRet=1, the place where
C     Val should be inserted is returned if I is not in List; if IRet=0,
C     zero is returned on failure.
C
      Integer U
      Dimension List(*)
C
      L = 1
      U = N
   10 If(U.lt.L) then
        If(IRet.eq.0) then
          BinSeI = 0
        else
          BinSeI = L
          endIf
        Return
        endIf
      I = (U+L)/2
      If(List(I)-IVal) 20, 30, 40
   20 L = I + 1
      Goto 10
   30 BinSeI = I
      Return
   40 U = I - 1
      Goto 10
      End
*Deck BLModB
      Function BLModB(IAn,JAn)
      Implicit Real*8(A-H,O-Z)
C
C     Return the model B bond length between atoms with atomic
C     numbers IAn and JAn, or zero if the bond length is not
C     defined.  A few of these bond lengths were defined in the
C     original Gordon/Pople Model B.  The remaining values are
C     somewhat arbitrary guesses based on optimized HF/3-21G
C     structures.  All bonds to He are arbitrarily set to 1.0
C     and all bonds to Ne are set to 1.5.
C
C     Mike Frisch -- May, 1983.
C
      Dimension BL(153)
      Save MaxAn, Zero, Bl
      Data MaxAn/17/, Zero/0.0d0/
C              H-H    He-H   He-He  Li-H   Li-He  Li-Li  Be-H   Be-He
      Data BL/0.74d0,0.74d0,0.74d0,1.60d0,1.60d0,2.82d0,1.34d0,1.34d0,
C       Be-Li  Be-Be  B-H    B-He   B-Li   B-Be   B-B    C-H    C-He
     $ 2.55d0,4.96d0,1.19d0,1.19d0,2.26d0,1.85d0,1.60d0,1.08d0,1.08d0,
C       C-Li   C-Be   C-B    C-C    N-H    N-He   N-Li   N-Be   N-B
     $ 2.10d0,1.60d0,1.40d0,1.40d0,1.00d0,1.00d0,1.85d0,1.45d0,1.45d0,
C       N-C    N-N    O-H    O-He   O-Li   O-Be   O-B    O-C    O-N
     $ 1.37d0,1.35d0,0.96d0,0.96d0,1.55d0,1.40d0,1.30d0,1.36d0,1.30d0,
C       O-O    F-H    F-He   F-Li   F-Be   F-B    F-C    F-N    F-O
     $ 1.48d0,0.92d0,0.92d0,1.52d0,1.35d0,1.35d0,1.35d0,1.36d0,1.42d0,
C       F-F    Ne-H   Ne-He  Ne-Li  Ne-Be  Ne-B   Ne-C   Ne-N   Ne-O
     $ 1.42d0,0.92d0,0.92d0,1.52d0,1.35d0,1.35d0,1.35d0,1.36d0,1.42d0,
C       Ne-F   Ne-Ne  Na-H   Na-He  Na-Li  Na-Be  Na-B   Na-C   Na-N
     $ 1.42d0,1.42d0,1.93d0,1.93d0,3.04d0,2.76d0,2.56d0,2.32d0,2.04d0,
C       Na-O   Na-F   Na-Ne  Na-Na  Mg-H   Mg-He  Mg-Li  Mg-Be  Mg-B
     $ 1.90d0,1.86d0,1.86d0,3.23d0,17.3d0,1.73d0,2.85d0,2.53d0,2.33d0,
C       Mg-C   Mg-N   Mg-O   Mg-F   Mg-Ne  Mg-Na  Mg-Mg  Al-H   Al-He
     $ 2.00d0,1.88d0,1.75d0,1.73d0,1.73d0,3.10d0,2.65d0,1.60d0,1.60d0,
C       Al-Li  Al-Be  Al-B   Al-C   Al-N   Al-O   Al-F   Al-Ne  Al-Na
     $ 2.55d0,2.25d0,2.00d0,1.90d0,1.70d0,1.63d0,1.65d0,1.65d0,2.99d0,
C       Al-Mg  Al-Al  Si-H   Si-He  Si-Li  Si-Be  Si-B   Si-C   Si-N
     $ 2.60d0,2.34d0,1.49d0,1.49d0,2.45d0,2.20d0,1.90d0,1.74d0,1.65d0,
C       Si-O   Si-F   Si-Ne  Si-Na  Si-Mg  Si-Al  Si-Si  P-H    P-He
     $ 1.61d0,1.63d0,1.63d0,2.83d0,2.53d0,2.25d0,2.15d0,1.42d0,1.42d0,
C       P-Li   P-Be   P-B    P-C    P-N    P-O    P-F    P-Ne   P-Na
     $ 2.30d0,2.00d0,1.77d0,1.71d0,1.64d0,1.65d0,1.60d0,1.60d0,2.70d0,
C       P-Mg   P-Al   P-Si   P-P    S-H    S-He   S-Li   S-Be   S-B
     $ 2.48d0,2.22d0,2.15d0,2.10d0,1.35d0,1.35d0,2.22d0,1.90d0,1.72d0,
C       S-C    S-N    S-O    S-F    S-Ne   S-Na   S-Mg   S-Al   S-Si
     $ 1.75d0,1.65d0,1.68d0,1.66d0,1.66d0,2.47d0,1.95d0,2.17d0,2.12d0,
C       S-P    S-S    Cl-H   Cl-He  Cl-Li  Cl-Be  Cl-B   Cl-C   Cl-N
     $ 2.17d0,2.14d0,1.29d0,1.29d0,2.11d0,1.85d0,1.80d0,1.80d0,1.80d0,
C       Cl-O   Cl-F   Cl-Ne  Cl-Na  Cl-Mg  Cl-Al  Cl-Si  Cl-P   Cl-S
     $ 1.76d0,1.69d0,1.69d0,2.42d0,2.26d0,2.24d0,2.15d0,2.24d0,2.23d0,
C       Cl-Cl
     $ 2.19d0/
      LInd(I,J) = (((Max(I,J)*(Max(I,J)-1))/2)+Min(I,J))
C
      BlModB = Zero
      If(IAn.gt.0.and.JAn.gt.0.and.IAn.le.MaxAn.and.JAn.le.MaxAn)
     $    BlModB = Bl(LInd(IAn,JAn))
      Return
      End
*Deck IMove1
      Subroutine IMove1(N,IndF,IndT,A)
      Implicit Integer(A-Z)
C
C     Move N elements from A(IndF+1) to A(IndT+1), doing the right
C     thing regardless of whether IndF is larger or smaller than IndF.
C
      Dimension A(*)
C
      If(IndF.gt.IndT) then
        Do 10 I = 1, N
   10     A(IndT+I) = A(IndF+I)
      else if(IndF.lt.IndT) then
        Do 20 I = N, 1, -1
   20     A(IndT+I) = A(IndF+I)
        endIf
      Return
      End
*Deck IMove
      Subroutine IMove(N,A,B)
      Implicit Integer(A-Z)
C
C     Copy N elements from A to B.  Do not use this routine for
C     overlapping arrays!
C
      Dimension A(*), B(*)
       Do 10 I = 1, N
   10    B(I) = A(I)
      Return
      End
*Deck IndSrt
      Integer Function IndSrt(iz,iy,ix)
      Implicit Real*8(a-h,o-z)
C
C     This subroutine forms a single index for FMM purposes
C     It works in the following fasion
C     iBit=0
C     ix occupy iBit+1, iBit+4, iBit+7, ...
C     iy  iBit+2, iBit+5, iBit+8, ...
C     iz  iBit+3, iBit+6, iBit+9, ...
C     Such a scheme will allow us to sort all the levels in one pass
C     Also when we use Radix sort it will be true linear scaling
C
      Integer Arr(3)
C
      Arr(1) = iz
      Arr(2) = iy
      Arr(3) = ix
      Ind = 0
      ip = 1
      Do 20 i = 1, 3
        ib = arr(i)
        imv = ip
   10   ir = Mod(ib,2)
        ib = ib/2
C       Here we put the bit into its position
        If(ir.ne.0) Ind = Ind + imv
        imv = imv*8
C       If nothing left, we go to the next arr(i)
        If(ib.gt.0) goto 10
   20   ip = ip*2
      IndSrt = Ind
      Return
      End
*Deck GLog
      function glog (arg)
      real*8 glog
      real*8 arg
c
c     base e (natural) logarithm
c
      glog=dlog(arg)
      return
      end
*Deck IGFix
      function igfix(r)
      implicit real*8(a-h,o-z)
c
c     floating to integer conversion.
c
      igfix = r
      return
      end
*Deck MaxInt
      Function MaxInt(IType)
      Implicit Integer(A-Z)
C
C     Return the largest positive integer if IType=0, or the maximum
C     integer which can be stored exactly in a wp value if IType=1.
C
      NB = MDNBPW(0)
C     Following is only correct for IEEE fp.
      If(IType.eq.1) NB = Min(NB,54)
      MaxInt = -((-(2**(NB-2)))*2+1)
      Return
      End
*Deck MDNBPW
      Function MDNBPW(Base)
      Implicit Integer(A-Z)
C
C     Return the number of bits (Base=0 or 2) or digits (Base=10) in
C     an integer word.
C
      If(Base.eq.0.or.Base.eq.2) then
       MDNBPW = 32 
      else if(Base.eq.10) then
       MDNBPW = 9
      endIf
      Return
      End
*Deck Eigen
      SUBROUTINE Eigen(A,R,N,MV)
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN DESCENDING ORDER.
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C                   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C                   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C                       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C                       SEQUENCE)
C
      DIMENSION A(1),R(1)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
      Real*8  A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
     $         COSX2,SINCS,RANGE
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  DSQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DDSQRT.  DABS IN STATEMENT
C        62 MUST BE CHANGED TO DDABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
    5 RANGE=1.0D-12
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0D0
      IF(I-J) 20,15,20
   15 R(IJ)=1.0D0
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0D0
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM) 165,165,40
   40 ANORM=SQRT(ANORM*2.0d0)
      ANRMX=ANORM*RANGE/DFLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/DFLOAT(N)
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
   62 IF(DABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5D0*(A(LL)-A(MM))
   68 Y=-A(LM)/DSQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 SINX=Y/DSQRT(2.0d0*(1.0d0+(DSQRT(1.0D0-Y*Y))))
      SINX2=SINX*SINX
   78 COSX=DSQRT(1.0D0-SINX2)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0D0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(LL)-A(MM)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
      RETURN
      END
*Deck EckOld 
      Subroutine EckOld(IOut,IPrint,NAtoms,AtMass,C,Com,PMom,RotMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension AtMass(*),C(3,*),Com(3),PMom(3),RotMat(9)
      Dimension AOrd(3),REig(9)
C
C  PUT MOMENTS OF INERTIA IN ORDER SUCH THAT IC
C  IS THE 'MOST UNEQUAL' (SIC) ONE
C
C      AORD(1)=AEIG(1)
C      AORD(2)=AEIG(3)
C      AORD(3)=AEIG(6)
      call AMove(3,PMom,AOrd)
      call AMove(9,RotMat,REig)
      if(IPrint.gt.1) then
       write(IOut,'(3F12.5)') (AOrd(i),i=1,3)
       write(IOut,'(3F12.5)') (Reig(i),i=1,9)
       do 10 I=1,NAtoms
       write(IOut,'(4F12.5)') AtMAss(i),(C(j,i),j=1,3)
   10  continue
      EndIf     
      DO 142 I=1,3
      DO 43 J=I,3
      I3=3*I
      I2=I3-1
      I1=I3-2
      J3=3*J
      J2=J3-1
      J1=J3-2
      IF(AORD(I).LE.AORD(J)) GO TO 43
      TEMPA=AORD(I)
      TEMPR1=REIG(I1)
      TEMPR2=REIG(I2)
      TEMPR3=REIG(I3)
      AORD(I)=AORD(J)
      REIG(I1)=REIG(J1)
      REIG(I2)=REIG(J2)
      REIG(I3)=REIG(J3)
      AORD(J)=TEMPA
      REIG(J1)=TEMPR1
      REIG(J2)=TEMPR2
      REIG(J3)=TEMPR3
  43  CONTINUE
  142 CONTINUE
      R21=AORD(2)/AORD(1)
      R32=AORD(3)/AORD(2)
      IF(R32.GE.R21) GO TO 149
      TEMPA=AORD(3)
      TEMPR1=REIG(7)
      TEMPR2=REIG(8)
      TEMPR3=REIG(9)
      AORD(3)=AORD(1)
      REIG(7)=REIG(1)
      REIG(8)=REIG(2)
      REIG(9)=REIG(3)
      AORD(1)=TEMPA
      REIG(1)=TEMPR1
      REIG(2)=TEMPR2
      REIG(3)=TEMPR3
 149  CONTINUE
C
C     CALCULATE CARTESIAN COORDINATES IN PRINCIPAL AXES SYSTEM
C
      If(IPrint.lt.1) Return
      WRITE(IOut,345)
  345 FORMAT(//,2X,'CARTESIAN COORDINATES IN PRINCIPAL AXES SYSTEM: '
     *,//,' ATOM     MASS      X            Y          Z',/)
      DO 27 IAt=1,NAtoms
        XI=C(1,IAt)-Com(1) 
        YI=C(2,IAt)-Com(2)
        ZI=C(3,IAt)-Com(3)
        XPRINC=XI*REIG(1)+YI*REIG(2)+ZI*REIG(3)
        YPRINC=XI*REIG(4)+YI*REIG(5)+ZI*REIG(6)
        ZPRINC=XI*REIG(7)+YI*REIG(8)+ZI*REIG(9)
        WRITE(IOut,'(I5,4F12.5)') IAt,AtMass(IAt),XPRINC,YPRINC,ZPRINC
   27 CONTINUE
      Return
      End
*Deck HQRII1
      Subroutine HQRII1 (IOut,N,IEV1,IEVL,IOrd,AL,EVal,NVX,EVec,Schmid,
     $  IErr,IX,WA,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C Householder QR Inverse Interation method 1
C     Compute the eigenvalues (EVal) of a symmetric matrix stored in a
C     lower triangular form (AL)
C     Based on the original HQRII:
C       Y. Beppu and I. Ninomiya, Comput. Chem., 6, 87 (1982)
C     Improved by A.V. Bunge and C.F. Bunge, Comput. Chem., 10, 259 (1986)
C
C Input:
C     N      : Dimension of the symmetric matrix A (NxN)
C     IEV1   : Index of the first wanted eigenvector (default 1)
C     IEVL   : Index of the last wanted eigenvector
C       > If IEVL < IEV1, eigenvectors are not computed
C     IOrd   : Order of the eigenvalues
C              IOrd >= 0: eigenvalues in decreasing order
C              IOrd < 0:  eigenvalues in increasing order
C     AL     : (Nx(N+1)/2) lower triangular array of the symmetric
C              matrix form stored as a vector
C     NVX    : number of columns for EVec
C     Schmid : if False, the eigenvectors of degenerate or nearly-degenerate
C              are not orthogonalized among themselves. True is the normal
C              use
C
C Output:
C     EVal   : (N) eigenvalues of A, in the order chosen with IOrd
C     EVec   : (NVX,N) eigenvectors of A
C     IErr   : Error flag
C                0    normal succesful completion
C               33    Error: AL is a null matrix
C              129    Error: N, LV or NVX are outside permissible bounds
C
C Local:
C     IX     : (N) working array to store various indexes
C     WA     : (6,N) working array for the actual calculations
C     V      : Scratch array
C
C
C     Dimension
      Integer MDV, N, NTT, NVX
C     Input
      Integer IEV1, IEVL, IOrd, IOut
      Logical Schmid
C     Output
      Integer IErr
      Real*8 AL(*), EVal(*), EVec(*)
C     Local
      Integer IX(*), IV, j, LV, NV
      Real*8 WA(6,*), V(*), Boundr, Eps0, Eta, GoldRt, MDCutO, One, Pt5,
     $  Zero
C
      Save Boundr, GoldRt, One, Pt5, Zero
      Data Boundr/1.0D-6/, GoldRt/0.618033988749894D0/, One/1.0D0/,
     $  Pt5/0.5D0/, Zero/0.0D0/
C
 9100 Format(' ERROR: At least one dim. in HQRII1 is outside the ',
     $  'permissible bounds')
 9110 Format(' ERROR: Null matrix in subroutine HQRII1')
C
      NTT = N*(N+1)/2
C     Call TstCor(NTT,MDV,'HQRII1')
      If(MDV.lt.NTT) then
       write(IOut,'('' Too Small Scratch Vector'')')
      endif
      Call AClear(NTT,V)
C     Test0 to control null values
      Test0 = MDCutO(1)
C     Eps0 is smallest number so that 1.+Eps0 =/= 1.0
      Call EpsEta(Eps0,Eta)
      IV = IEV1
      LV = IEVL
      NV = NVX
      IErr = 0
      If(IV.le.0) IV = 1
c      If (N.lt.1 .or. N.gt.NX .or. (N.gt.NV .and. LV.ge.IV) .or.
c     $    LV.gt.N) then
C     Control that dimensions are acceptable
      If(N.lt.1 .or. (N.gt.NV .and. LV.ge.IV) .or. LV.gt.N) then
        Write(IOut, 9100)
        IERR = 129
        Return
      endIf
C     Obvious case, N = 1
      If(N.eq.1) then
         EVal(1) = AL(1)
         EVec(1) = One
         Return
      endIf
C
C     General Case
      IX(1) = 0
      Do 10 j = 2, N
       IX(j) = IX(j-1) + j - 1
   10 Continue
      NM1 = N - 1
      NVF = (LV-1)*NV
      If(NVF.lt.0) NVF = 0
      If(N.gt.2) Then
        NM2 = N - 2
        Do 100 k = 1, NM2
          KP1 = k + 1
          WA(2,K) = AL(k+IX(k))
          Scale = Zero
          Do 110 j = KP1, N
           Scale = Scale + Abs(AL(IX(j)+k))
  110     Continue 
          WA(1,K) = AL(IX(KP1)+k)
          If(Scale.gt.Zero) then
            ScaleI = One/Scale
            Sum = Zero
            Do 111 j = KP1, N
              WA(2,J) = AL(IX(j)+k)*ScaleI
              Sum  = Sum + WA(2,J)**2
  111       Continue
            S = Sign(Sqrt(Sum),WA(2,KP1))
            WA(1,K)   = -S*Scale
            WA(2,KP1) = WA(2,KP1) + S
            AL(IX(KP1)+k) = WA(2,KP1)*Scale
            H    = WA(2,KP1)*S
            HUNS = (H*Scale)*Scale
            HI   = One/H
            SumM = Zero
            Do 112 ii = KP1, N
             WA(5,II) = Zero
  112       Continue
            Do 113 i = KP1, N
              IM1   = i - 1
              Sum   = ZERO
              I0    = IX(i)
              W2I   = WA(2,i)
              NRest = Mod(IM1-KP1+1,6)
              Do 120 j = KP1, KP1+NRest-1
                Sum = Sum + WA(2,j)*AL(I0+j)
                WA(5,J) = WA(5,j) + W2I*AL(I0+j)
  120         Continue
              Do 121 j = KP1+NRest, IM1, 6
               Sum = Sum + WA(2,j)*AL(I0+j) + WA(2,j+1)*AL(I0+j+1)
     C                + WA(2,j+2)*AL(I0+j+2) + WA(2,j+3)*AL(I0+j+3)
     C                + WA(2,j+4)*AL(I0+j+4) + WA(2,j+5)*AL(I0+j+5)
  121         Continue
              Do 122 j = KP1+NRest, IM1, 6
                WA(5,j)   = WA(5,j)   + W2I*AL(I0+j)
                WA(5,j+1) = WA(5,j+1) + W2I*AL(I0+j+1)
                WA(5,j+2) = WA(5,j+2) + W2I*AL(I0+j+2)
                WA(5,j+3) = WA(5,j+3) + W2I*AL(I0+j+3)
                WA(5,j+4) = WA(5,j+4) + W2I*AL(I0+j+4)
                WA(5,j+5) = WA(5,j+5) + W2I*AL(I0+j+5)
  122         Continue
              WA(6,i) = W2I*AL(I0+i) + Sum
  113       Continue
            Do 114 i = KP1, N
              WA(1,i) = (WA(5,i)+WA(6,i))*HI
              SumM    = WA(1,i)*WA(2,i) + SumM
  114       Continue
            U = Pt5*SumM*HI
            Do 115 i = KP1,N
              I0      = IX(i)
              WA(1,i) = WA(2,i)*U - WA(1,i)
              W1I     = WA(1,i)
              W2I     = WA(2,i)
              NRest   = Mod(i-KP1+1,6)
              Do 123 j = KP1, KP1+NRest-1
               AL(I0+j) = W2I*WA(1,j) + W1I*WA(2,j) + AL(I0+j)
  123         Continue
              Do 124 j = KP1+NRest, i, 6
                AL(I0+j)   = W2I*WA(1,j)   + W1I*WA(2,j)   + AL(I0+j)
                AL(I0+j+1) = W2I*WA(1,j+1) + W1I*WA(2,j+1) + AL(I0+j+1)
                AL(I0+j+2) = W2I*WA(1,j+2) + W1I*WA(2,j+2) + AL(I0+j+2)
                AL(I0+j+3) = W2I*WA(1,j+3) + W1I*WA(2,j+3) + AL(I0+j+3)
                AL(I0+j+4) = W2I*WA(1,j+4) + W1I*WA(2,j+4) + AL(I0+j+4)
                AL(I0+j+5) = W2I*WA(1,j+5) + W1I*WA(2,j+5) + AL(I0+j+5)
  124         Continue
  115       Continue
          else
            HUnS = Zero
          endIf
          AL(IX(k)+k) = HUnS
  100   Continue
      endIf
      NM1NM1    = IX(NM1) + NM1
      NM1N      = IX(N)   + NM1
      NN        = NM1N    + 1
      WA(2,NM1) = AL(NM1NM1)
      WA(2,N)   = AL(NN)
      WA(1,NM1) = AL(NM1N)
      WA(1,N)   = Zero
      GERSCH    = Abs(WA(2,1)) + Abs(WA(1,1))
      Do 200 i = 1,NM1
       GERSCH = Max(Abs(WA(2,i+1))+Abs(WA(1,i))+Abs(WA(1,i+1)),GERSCH)
  200 Continue
C     Trap null matrix before it is too late.
      If(GERSCH.lt.Test0) then
        Write(IOut, 9110)
        IErr = 33
        Return
      endIf
      SumD   = Zero
      SumCOD = Zero
      Do 300 i = 1, N
        SumCOD =  SumCOD + Abs(WA(1,i))
        SumD   =  SumD   + Abs(WA(2,i))
  300 Continue
      Scale  = SumD + SumCOD
      ScaleI = One/Scale
      Do 400 i = 1, N
        WA(1,i)     = WA(1,i)*ScaleI
        WA(2,i)     = WA(2,i)*ScaleI
        WA(3,i)     = WA(1,i)
        EVal(i)     = WA(2,i)
        EVec(i+NVF) = EVal(i)
  400 Continue
      Eps    = Sqrt(Eps0)
      GERSCH = GERSCH*ScaleI
      Del    = GERSCH*Eps
      DelW5  = GERSCH*Eps0
      If(SumD/SumCOD.gt.Del) Then
C       QR method with origin shift.
        Do 500 k = N, 2, -1
  510     Continue
            l = k
  520       If (Abs(WA(3,l-1)).gt.Del) then
              l = l - 1
              If(l.gt.1) Goto 520
            endIf
            If(l.ne.k) then
              WW        = (EVal(k-1)+EVal(k))*Pt5
              R         = EVal(k) - WW
              Z         = WW - Sign(Sqrt(WA(3,k-1)**2 + R*R),WW)
              EE        = EVal(l) - Z
              EVal(l)   = EE
              FF        = WA(3,l)
              R         = Sqrt(EE*EE + FF*FF)
              RI        = One/R
              C         = EVal(l)*RI
              S         = WA(3,l)*RI
              WW        = EVal(l+1) - Z
              EVal(l)   = (FF*C + WW*S)*S + EE + Z
              EVal(l+1) = C*WW - S*FF
              Do 521 j = l+1, k-1
                R         = Sqrt(EVal(j)**2 + WA(3,j)**2)
                RI        = One/R
                WA(3,j-1) = S*R
                EE        = EVal(j)*C
                FF        = WA(3,j)*C
                C         = EVal(j)*RI
                S         = WA(3,j)*RI
                WW        = EVal(j+1) - Z
                EVal(j)   = (FF*C + WW*S)*S + EE +Z
                EVal(j+1) = C*WW - S*FF
  521         Continue
              WA(3,k-1) = EVal(k)*S
              EVal(k)   = EVal(k)*C + Z
              Goto 510
            endIf
  500   Continue
C       Straight selection sort of eigenvalues.
        Sorter = One
        If(IOrd.lt.0) Sorter = -One
        j = N
  600   Continue
          l  = 1
          ii = 1
          ll = 1
          Do 610 i = 2, j
            If((EVal(i)-EVal(l))*Sorter.le.Zero) Then
              l = i
            else
              ii = i
              ll = l
            endIf
  610     Continue
          If(ii.ne.ll) Then
            WW       = EVal(ll)
            EVal(ll) = EVal(ii)
            EVal(ii) = WW
          endIf
          j = ii - 1
        If(j.gt.1) GoTo 600
      endIf
      If(LV.ge.IV) Then
C       Inverse iteration for eigenvectors.
        FN   = FLOAT(N)
        Eps1 = Sqrt(FN)*Eps
        SEps = Sqrt(Eps)
        Eps2 = (GERSCH*Boundr)/(FN*SEps)
        RN   = Zero
        RA   = Eps*GoldRt
        i2   = (IV-2)*NV
        Do 700 i = IV, LV
          i2 = i2 + NV
          Do 710 j = 1, N
            WA(3,j) = Zero
            WA(4,j) = WA(1,j)
            WA(5,j) = EVec(NVF+j) - EVal(i)
            RN      = RN + RA
            If(RN.ge.Eps) RN = RN - Eps
            WA(6,j) = RN
  710     Continue
          Do 711 j = 1, NM1
            If(Abs(WA(5,j)).le.Abs(WA(1,j))) Then
              If(Abs(WA(1,j)).lt.Test0) WA(1,j) = Del
              WA(2,j)   = -WA(5,j)/WA(1,j)
              WA(5,j)   =  WA(1,j)
              T         =  WA(5,j+1)
              WA(5,j+1) =  WA(4,j)
              WA(4,j)   =  T
              WA(3,j)   =  WA(4,j+1)
              If(Abs(WA(3,j)).lt.Test0) WA(3,j) = Del
              WA(4,j+1)=  Zero
            else
              WA(2,j)   = -WA(1,j)/WA(5,j)
            endIf
            WA(4,j+1) = WA(3,j)*WA(2,j) + WA(4,j+1)
            WA(5,j+1) = WA(4,j)*WA(2,j) + WA(5,j+1)
  711     Continue
          If(Abs(WA(5,N)).lt.Test0) WA(5,N) = DelW5
          WNM15I = One/WA(5,NM1)
          WN5I   = One/WA(5,N)
          Do 712 Itere=1,2
            If(Itere.ne.1) Then
              Do 720 j = 1, NM1
                If(WA(5,j).eq.WA(1,j)) Then
                  T         = WA(6,j)
                  WA(6,j)   = WA(6,j+1)
                  WA(6,j+1) = T
                endIf
                WA(6,j+1) = WA(6,j)*WA(2,j) + WA(6,j+1)
  720         Continue
            endIf
            WA(6,N)   = WA(6,N)*WN5I
            WA(6,NM1) = (WA(6,NM1)-WA(6,N)*WA(4,NM1))*WNM15I
            VN        = Max(Abs(WA(6,N)),Abs(WA(6,NM1)))
            Do 721 k = NM2, 1, -1
              WA(6,k) = (WA(6,k)-WA(6,k+1)*WA(4,k)-WA(6,k+2)*WA(3,k))
     $          / WA(5,k)
              VN = Max(Abs(WA(6,k)),VN)
  721       Continue
            S     = Eps1/VN
            NRest = Mod(N,6)
            Do 722 j = 1, NRest
             WA(6,j) = S*WA(6,j)
  722       Continue 
            Do 723 j = 1+NRest, N, 6
              WA(6,j)   = S*WA(6,j)
              WA(6,j+1) = S*WA(6,j+1)
              WA(6,j+2) = S*WA(6,j+2)
              WA(6,j+3) = S*WA(6,j+3)
              WA(6,j+4) = S*WA(6,j+4)
              WA(6,j+5) = S*WA(6,j+5)
  723       Continue
  712     Continue
          Do 713 j = 1, N
           EVec(i2+j) = WA(6,j)
  713     Continue
  700   Continue
C Build indexing and upper triangular matrix
        IX(1) = 0
        Do 800 j = 2, N
         IX(j) = IX(j-1) - j + 1 + N
  800   Continue 
        ij = 0
        Do 900 j = 1, N
          Do 910 i = 1, j
            ij = ij + 1
            V(IX(i)+j) = AL(ij)
  910     Continue
  900   Continue
C Back transformation of eigenvectors.
        ig = 1
        i2 = (IV-2)*NV
        Do 1000 i = IV, LV
          i2 = i2 + NV
          Do 1010 j = 1, N
           WA(6,j) = EVec(i2+j)
 1010     Continue
          IM1 = i - 1
          If(N.gt.2) Then
            Do 1020 j = 1, NM2
              k  = N - j - 1
              k0 = IX(k)
              If(V(k0+k).ne.Zero) Then
                KP1   = k + 1
                Sum   = Zero
                NRest = Mod(N-KP1+1,6)
                Do 1030 kk = KP1, KP1+NRest-1
                 Sum = V(K0+kk)*WA(6,kk) + Sum
 1030           Continue
                Do 1031 kk = KP1+NRest, N, 6
                 Sum = V(k0+kk)*WA(6,kk) + V(k0+kk+1)*WA(6,KK+1)
     $                  + V(k0+kk+2)*WA(6,kk+2) + V(k0+kk+3)*WA(6,kk+3)
     $                  + V(k0+kk+4)*WA(6,kk+4) + V(k0+kk+5)*WA(6,kk+5)
     $                  + Sum
 1031           Continue
                S = -Sum/V(k0+k)
                Do 1032 kk = KP1, KP1+NRest-1
                 WA(6,kk) =  S*V(k0+kk) + WA(6,kk)
 1032           Continue 
                Do 1033 kk = KP1+NRest, N, 6
                  WA(6,kk)   =  S*V(k0+kk)   + WA(6,kk)
                  WA(6,kk+1) =  S*V(k0+kk+1) + WA(6,kk+1)
                  WA(6,kk+2) =  S*V(k0+kk+2) + WA(6,kk+2)
                  WA(6,kk+3) =  S*V(k0+kk+3) + WA(6,kk+3)
                  WA(6,kk+4) =  S*V(k0+kk+4) + WA(6,kk+4)
                  WA(6,kk+5) =  S*V(k0+kk+5) + WA(6,kk+5)
 1033           Continue
              endIf
 1020       Continue
          endIf
          j = ig
 1011     If (Abs(EVal(j)-Eval(i)).ge.Eps2) then
            j = j + 1
            If(j.le.i) GoTo 1011
          endIf
          ig = Min(j,i)
          NRest = Mod(N,6)
C
          If(ig.ne.i .and. Schmid) Then
C Degenerate eigenvalues.First,orthogonalize.
            KF = (ig-2)*NV
            Do 1012 k = ig, IM1
              KF  = KF + NV
              Sum = Zero
              Do 1021 j = 1, NRest
               Sum = EVec(KF+j)*WA(6,j) + Sum
 1021         Continue
              Do 1022 J = 1+ NRest, N, 6
                Sum = EVec(KF+j)*WA(6,j)     + EVec(KF+j+1)*WA(6,j+1)
     *              + EVec(KF+j+2)*WA(6,j+2) + EVec(KF+j+3)*WA(6,j+3)
     *              + EVec(KF+j+4)*WA(6,j+4) + EVec(KF+j+5)*WA(6,j+5)
     $              + Sum
 1022          Continue
              S = -Sum
              Do 1023 j = 1, NRest
               WA(6,j) = S*EVec(KF+j) + WA(6,j)
 1023         Continue
              Do 1024 j = 1+NRest, N, 6
                WA(6,j)   = S*EVec(KF+j)   + WA(6,j)
                WA(6,j+1) = S*EVec(KF+j+1) + WA(6,j+1)
                WA(6,j+2) = S*EVec(KF+j+2) + WA(6,j+2)
                WA(6,j+3) = S*EVec(KF+j+3) + WA(6,j+3)
                WA(6,j+4) = S*EVec(KF+j+4) + WA(6,j+4)
                WA(6,j+5) = S*EVec(KF+j+5) + WA(6,j+5)
 1024         Continue
 1012       Continue
          endIf
C Normalization
          Sum = Zero
          Do 1013 j = 1, NRest
           Sum = WA(6,j)**2 + Sum
 1013     Continue 
          Do 1014 j = 1+NRest, N, 6
           Sum = WA(6,j  )**2 + WA(6,j+1)**2 + WA(6,j+2)**2
     *           + WA(6,j+3)**2 + WA(6,j+4)**2 + WA(6,j+5)**2 + Sum
 1014     Continue
          S = One/Sqrt(Sum)
          Do 1015 j = 1, NRest
           EVec(i2+j) = S*WA(6,j)
 1015     Continue
          Do 1016 j = 1+NRest, N, 6
            EVec(i2+j  ) = S*WA(6,j  )
            EVec(i2+j+1) = S*WA(6,j+1)
            EVec(i2+j+2) = S*WA(6,j+2)
            EVec(i2+j+3) = S*WA(6,j+3)
            EVec(i2+j+4) = S*WA(6,j+4)
            EVec(i2+j+5) = S*WA(6,j+5)
 1016     Continue
 1000   Continue
      endIf
      Do 1100 i = 1, N
       EVal(i) = Scale*EVal(i)
 1100 Continue
      Return
      End

