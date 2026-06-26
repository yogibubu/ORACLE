*Deck Coord
      Subroutine Coord(In,IOut,IPunch,IPrint,MxAt,MaxNZ,MxBnd,MxBox,
     $ PhyCon,KWd,Multip,NAtoms,NFrag,NBond,NH,NZ,IAn,Isot,IFrag,IBond,
     $ IANZ,IZ,MapZAt,Linear,C,CZ,EAN,EANZ,AtMass,TotWT,PMom,RotGHz,
     $ RTemp,MultN,QMom,GFac,NSpec,IScr,Scr,Group)
      Implicit Real*8 (A-H,O-Z)
      Integer St2Int
C I/O
      Logical KWd(*),RdIsot
      Dimension IAn(*),IFrag(*),NBond(*),IBond(MxBnd,*),NH(*)
      Dimension ISot(*),MultN(*)
      Dimension IAnZ(*),IZ(4,*),MapZAt(*)
      Dimension QMom(*),GFac(*)
      Dimension PhyCon(*),C(3,*),EAN(*),PMom(*),RotGHz(*),RTemp(*)
      Dimension AtMass(*),CZ(3,*),EANZ(*)
      Dimension IScr(*),Scr(*)
C Local
      Logical DoBPCS,DoHBnd,DoEck,TstAng,TTest,Error,AllHB,JoinFr,Linear
      Logical Aver,OK,Clean
      Character*2 AtSymb(MxAt)
      Character*2 SymAt(MxAt),IAnEl2
      Character*20 StrInp,SL,SA,SD
      Character*6 Str1,Str2,Str3
      Character*16 Group, Group0
      Character*8 SymGroup
      Character*80 LinScr
      Dimension IStart(100),IZIZ(4)
      Dimension CInp(3),Values(MaxNZ)
      Dimension IstL(100),IStA(100),IStB(100)
      Dimension IDnHB(100),IHAtHB(100),ICcHB(100)
      Dimension TMom(6)
C Dimensions for Rotational Constants 
      Dimension XYZCM(3),Rotcm1(3),Dip(3),PMomB(3),RotMat(3,3)
      Dimension XSym(MxAt),YSym(MxAt),ZSym(MxAt)
C Dimensions for Z-Matrix
      Dimension LBl(MaxNZ),LAlpha(MaxNZ),LBeta(MaxNZ)
      Dimension BL(MaxNZ),Alpha(MaxNZ),Beta(MaxNZ)
      Character*4000 StBL,StAlph,StBet
      Character*4 SBVar,SAVar,SDVar
C Scratch for ZTOC
      Dimension A(MaxNZ),B(MaxNZ),D(MaxNZ),Alpha1(MaxNZ),Beta1(MaxNZ)
C Dimensions for Connectivity
      Dimension IBoxSt(MxBox),IBox(MxBox),IndTab(MxBox),IndBox(MxBox)
      Dimension IP2Box(4,MxBox),IBType(MxBnd,MxAt)
C Bonder are Bond Orders
C ACN and ACN1 are coordination numbers (from D3 and D4 of Grimme)
C SPIJ is the sum of bond orders
      Dimension Bonder(MxBnd,MxAt)
      Dimension ACN(MxAt),ACN1(MxAt),SPIJ(MxAt)
C Long Formats
 1000 Format(' Atomic Number =',I3,'; Isotope =',I3,
     $  '; Nuc.Spin Multip =',I5,/,' Mass =',F12.5,
     $  '; Quadrup.Mom. =',F10.5,'; Nuc.Magn.Mom. =',F10.5)
C Numerical Data
      pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
      ToAng=1.0D0
      Thresh=1.0d-10
C Set defaults and starting values
      Aver   = Kwd(5)
      DoEck  = KWd(6)
      DoBPCS = KWd(13)
      JoinFr = KWd(16)
      DoHBnd = KWd(17)
      RdIsot = .False.
      Clean  = Kwd(31)
      NAtoms = 0
      NFrag  = 0
      IAt    = 0
      IZAt   = 0
      IToAng = 1
      IAlg   = 1
      NVBL   = 0
      NFBL   = 0
      NVAlph = 0
      NFAlph = 0
      NVBeta = 0
      NFBeta = 0
      IEndL  = 0
      IEndA  = 0
      IEndB  = 0
      Linear = .False.
      call AClear(3,CInp)
      call IClear(4,IZIZ)
C Read Cartesian coordinates.  GICForge is deliberately xyz-only for
C molecular geometry: "provin" controls the run, but never carries geometry.
C Do not add silent fallbacks to FCHK, Z-matrix, or stdin coordinates here.
      OK = .False.
      INQUIRE( FILE='xyzin', EXIST=OK )
      If(.not.OK) then
       write(IOut,'('' File XYZIN does not exist'')')
       Stop
      EndIf
      InFil=8
      OPEN(InFil,FILE='xyzin',STATUS='OLD')
      Read(InFil,*) NAtoms
      If(NAtoms.le.1) then
       write(IOut,'(/'' Less than 2 Atoms: STOP'')')
       STOP
      EndIf
      write(IOut,'(/,'' Cartesian Coordinates from XYZ File for'',
     $ I5,'' Atoms'',/)') NAtoms
      Read(InFil,'(A80)') LinScr
      CALL FndGrp(LinScr, GROUP0, GROUP)
      IF (GROUP .EQ. ' ') THEN
       GROUP0='C1'
       GROUP='c1'
      END IF
      Do 10 IAt=1,NAtoms
       Call LlinCl(LinScr)
       Read(InFil,'(A80)') LinScr
       Call RdCart(IOut,IPrint,IAt-1,LinScr,NValue,IAnI,IsotI,IFragI,
     $  CInp)
       If(NValue.ne.0) then
        IAn(IAt)=IAnI
        Isot(IAt)=IsotI
        IFrag(IAt)=IFragI
        If(IFragI.gt.NFrag) NFrag=IFragI
        Call AMove(3,CInp,C(1,IAt))
       Else
        Write(IOut,'('' Error for Atom'',I5)') IAt
        STOP
       EndIf
   10 Continue
      Close(InFil)
C Enhanced connectivity generation
C      call GenCBx(IOut,IPrint,IAlg,IToAng,DLimI,MxBnd,NAtoms,IAn,C,
C    $   NBond,IBond,IBType,IP2Box,IndTab,IndBox,IBox,IBoxSt)   
C IBtype is the integer bond order
C Primitive construction of bonds
      call MkBnd(IOut,IPrint,MxBnd,NAtoms,IAn,NBond,IBond,C)
C the bonds are updated by MkCNA
      do 100 iat=1,NAtoms
       IFrag(IAt)=1
  100 continue    
C Set various atomic properties
      If(NFrag.eq.0) NFrag=1
      TotWT = 0.0d0
      Do 110 IAt=1,NAtoms
       RMass = 0.0d0
       ISpin = -1
       QMom1 = 0.0d0
       GFac1 = 0.0d0
       IAnI  = IAn(IAt)
       MNI   = ISot(IAt)
       Call FilMag(MNI,IAnI,JUse,MNO,RMass,ISpin,QMom1,GFac1)
       MultN(IAt) = ISpin+1
       QMom(IAt) = QMom1
       GFac(IAt) = GFac1
       ISot(IAt) = MNO 
       If(Aver) then
        Call FilAMS(1,IAn(IAt),AtMass(IAt))
       Else
        AtMass(IAt) = RMass
       EndIf
       TotWt = TotWt + AtMass(IAt)
       If(IPrint.gt.0) Write(IOut,1000) IAn(IAt),ISot(IAt),MultN(IAt),
     $  AtMass(IAt),QMom(IAt),GFac(IAt)
  110 continue  
C NSigma, NPi,NVal in IScr(INSig),IScr(INPi), IScr(INVal)
      INPI=1 
      INLP=INPI+NAtoms
C Updates bond topology (NBond,IBond)
      call MkCNA(IOut,IPrint,Multip,MxBnd,NAtoms,IAn,NBond,IBond,
     $  Iscr(INPI),IScr(INLP),C,ACN,ACN1,SPIJ)
C Make Bond Orders (Bonder) with  Del(i) in Scr(idel)
      IDel=1
      call BndOrd(IOut,IPrint,MxBnd,MxAt,NAtoms,IAn,NBond,IBond,
     $  IScr(INPI),IScr(INLP),Scr(IDel),Bonder,C)
C Find H-Bonds
      AllHB=.true.
      if(DoHBnd) then
       call FindHBnd(IOut,IPrint,MxBnd,AllHB,NAtoms,NFrag,NHB,
     $  IAn,NBond,IBond,IFrag,C,IDnHB,IHAtHB,ICcHB)
       If(NHB.gt.0) then
        write(IOut,'(I5,'' H-Bonds detected as non-covalent targets'',
     $   '' (not added to GIC topology)'')') NHB
       EndIf
      EndIf
C Make Coordination Numbers, Synthons, Sigma and Pi bonds, Lone Pairs and Unpair.El.
      call MKEAN(IOut,IPrint,MxBnd,NAtoms,IAN,NBond,IBond,C,Bonder,
     $ EAN,EANZ)
C Find Fragments    
      call FndFrg(MxBnd,NAtoms,IBond,NBond,IFrag,IScr)
      NFrag=IrMax1(IFrag,NAtoms,.True.,NAtFrM)
      If(NFrag.gt.1.and.JoinFr) then 
       AllHB=.false.
       call FindHBnd(IOut,IPrint,MxBnd,AllHB,NAtoms,NFrag,NHB,IAn,
     $   NBond,IBond,IFrag,C,IDnHB,IHAtHB,ICcHB)
       If(NHB.gt.0) then
        write(IOut,'(I5,'' Inter-Molecular H-Bonds detected but not'',
     $   '' added to GIC topology'')') NHB
       EndIf
      EndIf
      If(NFrag.gt.1) then
       If(JoinFr) then
        write(IOut,'('' STOP Because'', I3,'' Fragments Survive after'',
     $    '' JoinFr'')') NFrag
        Stop
       Else
        write(IOut,'(I3,'' Fragments'')') NFrag
       EndIf
      EndIf   
C Compute principal inertia axes and Eckart orientation 
      call MofI(IOut,Iprint,.false.,1,NAtoms,C,AtMass,XYZCM,TMom,PMom,
     $  RotMat)
      Do 130 IAt=1,NAtoms
       SymAt(IAt)=IAnEl2(IAn(IAt))
       DX=C(1,IAt)-XYZCM(1)
       DY=C(2,IAt)-XYZCM(2)
       DZ=C(3,IAt)-XYZCM(3)
       XSym(IAt)=RotMat(1,1)*DX+RotMat(2,1)*DY+RotMat(3,1)*DZ
       YSym(IAt)=RotMat(1,2)*DX+RotMat(2,2)*DY+RotMat(3,2)*DZ
       ZSym(IAt)=RotMat(1,3)*DX+RotMat(2,3)*DY+RotMat(3,3)*DZ
  130 Continue
      Call DETERMINE_POINT_GROUP(MxAt,NAtoms,SymAt, XSym,YSym,ZSym,
     $ 1.0D-4,SymGroup,ISymStat,DelSym)
      Write(IOut,'(/,'' Point Group from symm.f: '',A8)') SymGroup
      If(ISymStat.eq.0) then
       Write(IOut,'('' Symmetry quality: STRICT, max deviation ='',
     $  1PE10.3,'' Angstrom'')') DelSym
      ElseIf(ISymStat.eq.1) then
       Write(IOut,'('' Symmetry quality: QUASI, max deviation ='',
     $  1PE10.3,'' Angstrom'')') DelSym
      Else
       Write(IOut,'('' Symmetry quality: BROKEN, max deviation ='',
     $  1PE10.3,'' Angstrom'')') DelSym
      EndIf
      if(DoEck) then
       call Eckart(IOut,IPrint,DoEck,NAtoms,C,AtMass,XYZCM,Linear,PMom,
     $   RotMat)
       call CntMas(0,1,NAtoms,AtMass,C,TotM,Scr)
       if(Abs(SCr(1)+Scr(2)+Scr(3)).gt.1.0d-5) then
        write(IOut,'(/,'' Wrong Center of mass at'',3F10.5)')
     $   (Scr(i),i=1,3)
        STOP
       endif
CEnzo
C Set Equivalent Atoms
CEnzo
      endif
C Compute rotational parameters
      FactA = CnvFct('FactA')
      Fac1  = PhyCon(1)**2
      AvPMom=1.0d0
      AvPMoB=1.0d0
      do 120 I=1,3
       If(PMom(i).le.0.0d0) goto 120 
       AvPMom=AvPMom*PMom(i)
       PMomB(i)=PMom(i)/Fac1
       AvPMoB=AvPMoB*PMomB(i)
       Rotcm1(i)=FactA/PMom(i)
  120  Continue
      AvPMom=AvPMom**(1.0D0/3.0D0)
      AvPMoB=AvPMoB**(1.0d0/3.0d0)
      Call RotCon(IOut,IPrint,Linear,NAtoms,PhyCon,PMomB,RotGHz,RTemp)
      Return
      End
*Deck RdCart
      Subroutine RdCart(IOut,IPrint,IAt,CLine,NValue,IAn,ISotp,IFrag,C)
      Implicit Real*8 (A-H,O-Z)
      Integer St2Int,IType,IVal
      Character*(*) CLine
      Character CType*1,El*2,Sval*20
      Real*8 RVal
      Complex*16 CVal
      Logical LVal
      Dimension IStart(10),c(3)
  100 Format(' For Atom',I5)
      call SubStr(CLine,10,IStart,Nvalue)
      If(NValue.eq.0) Return
      Isotp=0
      IAt=IAt+1
C find atomic number or atomic symbol and isotope (separate by -)
      ICur=IStart(1) 
      If(CLine(ICur:ICur+1).eq.'D ') then
       IAn=1
       Isotp=2
       goto 10
      ElseIf(CLine(ICur:ICuri+1).eq.'T') then
       IAn=1
       Isotp=3
       goto 10
      EndIf
      IIso = 0
      If(CLine(ICur+2:ICur+2).eq.'-')  IIso = IStart(1)+3
      If(CLine(ICur+1:Icur+1).eq.'-')  IIso = IStart(1)+2
      If(IIso.eq.0) then
       El=CLine(IStart(1):IStart(1)+1)
      Else
       El=CLine(IStart(1):IIso-2)
      EndIf
      call FilIAn(El,Number)
      If(Number.gt.-1) then
       IAn=Number
       If(IIso.ne.0) then
        call St2Dat(CLine,0,0,IIso,IType,CType,IVal,RVal,CVal,LVal,
     $    SVal)
        Isotp = IVal
       EndIf
       goto 10
      EndIf
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      if(CType.eq.'S') then
       If(IIso.eq.0) then
        El=CLine(IStart(1):IStart(1)+1)
       Else
        El=CLine(IStart(1):IIso-2)
       EndIf
       call FilIAn(El,Number)
       if(Number.eq.-1) then
        write(IOut,'('' Wrong Atomic Number for Atom'',I5)') IAt
        STOP 
       Else
        IAn=Number
       EndIf
       If(IIso.ne.0) then
        call St2Dat(CLine,0,0,IIso,IType,CType,IVal,RVal,CVal,LVal,
     $    SVal)
        Isotp = IVal
       endif
      elseif(CType.eq.'I') then
       IAn = IVal
      else
       Write(IOut,'(/,''Wrong Atomic Number for Atom'',I5)') IAt
       STOP
      endif
   10 continue  
C find coordinates 
      Do 20 IXYZ=1,3 
       ICur=IStart(IXYZ+1)
       call St2Dat(CLine,0,0,ICur,IType,CType,IVal,C(ixyz),CVal,
     $  LVal,SVal)
       If(CType.ne.'R') then
        Write(IOut,'(/,''Wrong Coordinates for Atom'',I5)') IAt
        STOP
       EndIf
   20 continue
C find fragments
      If(NValue.eq.4) Return
      ICur=IStart(5)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(CType.ne.'I') then
       Write(IOut,'(/,''Wrong Fragment for Atom'',I5)') IAt
       STOP
      EndIf     
      IFrag=IVal 
      Return
      End
*Deck NewRdZ
      Subroutine NewRdZ(In,IOut,IPunch,IPrint,MaxNZ,NAtoms,NZ,Kwd,
     $  DoBPCS,LinScr,AtSymb,IAn,IANZ,Isot,MapZAt,IScr,C,CZ,Scr,NSpec)
      Implicit Real*8(A-H,O-Z)      
      Logical KWd(*),DoBPCS,RdIsot,GauCon
      Dimension IAn(*),IANZ(*),ISot(*),MapZAt(*),IScr(*)
      Dimension C(3,*),CZ(3,*),Scr(*)
      Character*2 AtSymb(*)
C Local
      Logical TstAng,TTest,Error,OK
      Character*20 StrInp,SL,SA,SD
      Character*6 Str1,Str2,Str3
      Character*80 LinScr
      Dimension IStart(100),IZIZ(4)
      Dimension CInp(3)
      Dimension IstL(100),IStA(100),IStB(100)
C Dimensions for Z-Matrix
      Dimension LBl(MaxNZ),LAlpha(MaxNZ),LBeta(MaxNZ),IZ(4,MaxNZ)
      Dimension A(MaxNZ),B(MaxNZ),D(MaxNZ),Alpha(MaxNZ),Beta(MaxNZ)
      Dimension BL(MaxNZ),Alpha1(MaxNZ),Beta1(MaxNZ)
      Character*4000 StBL,StAlph,StBet
      Character*4 SBVar,SAVar,SDVar
      Character*80 FilNam,OutFil
C Numerical parameters
      pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
      ToAng=1.0D0
      ToDeg=1.0d0/ToRad
      NAtoms = 0
      NFrag  = 0
      IAt    = 0
      IZAt   = 0
      IToAng = 1
      IAlg   = 1
      NVBL   = 0
      NFBL   = 0
      NVAlph = 0
      NFAlph = 0
      NVBeta = 0
      NFBeta = 0
      IEndL  = 0
      IEndA  = 0
      IEndB  = 0
      RdIsot = .False.
      GauCon = Kwd(4)       
      filnam='topo24'
C Read Z-Matrix
   10 if(IZAt.ne.0) then
       call LlinCl(LinScr)
       Read(In,'(A80)') LinScr
      EndIf
      call IClear(4,IZIZ)
      call AClear(3,CInp)
      call RdZMat(IOut,IPrint,IZAt,LinScr,NValue,IAnZI,IsotI,IZIZ,
     $  CInp,SL,SA,SD,NBLVar,NBLFix,NAlVar,NAlFix,NBeVar,NBeFix,LBlI,
     $  LAlphI,LBetaI,ILnSL,ILnSA,ILnSB)
      If(NValue.eq.0) goto 30
      IAnZ(IZAt)=IAnZI
      Isot(IZAt)=IsotI
      LBL(IZAt)=LBLI
      LAlpha(IZAt)=LAlphI
      LBeta(IZAt)=LBetaI
      BL(IZAt)=CInp(1)
      Alpha(IZAt)=CInp(2)*ToRad
      Beta(IZAt)=CInp(3)*ToRad
      IniL=IEndL+2
      IniA=IEndA+2
      IniB=IEndB+2
      IEndL=IniL+ILnSL
      IEndA=IniA+ILnSA
      IEndB=IniB+ILnSB
      StBL(IniL:IEndL)=SL(1:ILnSL)
      StAlph(IniA:IEndA)=SA(1:ILnSA)
      StBet(IniB:IEndB)=SD(1:ILnSB)
      do 20 i1=1,4
       IZ(i1,izat)=IZIZ(i1)
   20 Continue
      goto 10
   30 NZ=IZAt
      NZVar=NBLVar+NAlVar+NBeVar
      NZFix=NBLFix+NAlFix+NBeFix
      Call LinUpC(StBL,StBL)
      Call LinUpC(StAlph,StAlph)
      Call LinUpC(StBet,StBet)
      write(IOut,'(I5,'' Variable and'',I5,'' Fixed Parameters'',
     $  '' in Z-Matrix'')') NZVar,NZFix
      call SubStr(StBL,NZ,IStL,NvalL)
      call SubStr(StAlph,NZ,IStA,NValA)
      call SubStr(StBet,NZ,IStB,NValB)
      IStL(NZ+1)=IStL(NZ)+20
      IStA(NZ+1)=IStA(NZ)+20
      IStB(NZ+1)=IStB(NZ)+20
      NVZFnd=0
      NFZFnd=0
      call RdVarZ(In,IOut,IPrint,ToRad,NZ,NVar,BL,Alpha,Beta,StBL,
     $  StAlph,StBet,NVZFnd,LBl,LAlpha,LBeta,IStL,IStA,IStB,.true.)
      Write(IOut,'(I5,'' Different Variable Z-matrix Values Found'',
     $  ''  by RdVarZ'')') NVar
      If(NVZFnd.lt.NZVar) then
       call RdVarZ(In,IOut,IPrint,ToRad,NZ,NFix,BL,Alpha,Beta,StBl,
     $   StAlph,StBet,NFZFnd,LBl,LAlpha,LBeta,IStL,IStA,IStB,.false.)
       If(NFix.ne.0) Write(IOut,'(I5,'' Different Fixed Z-matrix'',
     $   '' Values Found by RdVarZ'')') NFix
      EndIf
C Possibly apply BDPCS3 corrections to Bond Lengths;
C Then set Kwd(13)=.false. in order to avoid double correction
      If(DoBPCS) then
       call MkBPCS(IOut,IPrint,NZ,IANZ,IZ,Bl)
       Kwd(13)=.false.
      EndIf
C Print Z-matrix
      ToAng=1.0d0
      ToDeg=1.0d0/ToRad
      call ZPrint(IOut,NZ,IANZ,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,ToAng,
     $  ToDeg)
C Write Z-matrix for CFour
      If(RdIsot) then
       read(In,*) NSpec
       OPEN(IPunch,FILE='msrin',STATUS='UNKNOWN')
       Rewind(IPunch)
       If(GauCon) then
        write(IPunch,'(''#m optim=(method=gaun,coord=zmat)'',
     $   '' geom=coord=zmat symmetry=isotopes vibr=gauconv'')')       
       Else        
        write(IPunch,'(''#m optim=(method=gaun,coord=zmat)'',
     $   '' geom=coord=zmat symmetry=isotopes'')')
       EndIf
       If(NSpec.lt.10) then
        write(IPunch,'(''#d niso='',I1,'' imax=20'')') NSpec
       Else
        write(IPunch,'(''#d niso='',I2,'' imax=20'')') NSpec
       EndIf
       write(IPunch,'('' '')')
       write(IOut,'(/,'' Z-matrix for CFOUR:'',I3,'' Lines'')') NZ
       call PrtZmt(IPunch,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $   LAlSh,LBeSh)
       call PrtVrZ(IPunch,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $   LAlSh,LBeSh)
       close(IPunch)
      EndIf  
C transform to Cartesian coordinates
C checking tetrahetral angles (TTest) and 0>teta<180 (TstAng)
      TstAng=.true.
      TTest=.true.
      call ZtoC(MaxNZ,NZ,IAnZ,IZ,Bl,Alpha,Beta,TTest,NAtoms,IAn,C,
     $  CZ,A,B,D,Alpha1,Beta1,IOut,Error,TstAng)
      call MkMpZA(IOut,NZ,IAnZ,NAtoms,MapZAt)
      return
      end
*Deck RdZMat
      Subroutine RdZMat(IOut,IPrint,IZAt,Cline,NValue,IAnZ,Isotp,IZIZ,
     $  C,SL,SA,SD,NBLVar,NBLFix,NAlVar,NAlFix,NBeVar,NBeFix,LBl,LAlpha,
     $  LBeta,ILnSL,ILnSA,ILnSD)
      Implicit Real*8 (A-H,O-Z)
      Logical TTest,Error,TsTAng,LVal
      Integer St2Int,IType,IVal
      Character*(*) CLine,SL,SA,SD
      Character CType*1,El*2,Sval*20,StrAdd*4
      Real*8 RVal
      Complex*16 CVal
      Dimension IStart(10),c(3),iziz(4)
  100 Format(' For Atom',I5)
      call LRmCom(CLine) 
      call SubStr(CLine,8,IStart,Nvalue)
      If(NValue.eq.0) Return
      IZAt=IZAt+1
      call NumChar(IZAt,StrAdd,LenStr)
      LBL=0
      LAlpha=0
      LBeta=0
      SL(1:4)='    '
      SA(1:4)='    '
      SD(1:4)='    '
      ILnSL=4
      ILnSA=4
      ILnSD=4
      Isotp=0
C
C find atomic number or atomic symbol and isotope (separate by -)
C
      ICur=IStart(1) 
      If(CLine(ICur:ICur+1).eq.'D ') then
       IAnZ=1
       Isotp=2
       goto 10
      ElseIf(CLine(ICur:ICuri+1).eq.'T') then
       IAnZ=1
       Isotp=3
       goto 10
      EndIf
      IIso = 0
      If(CLine(ICur+2:ICur+2).eq.'-')  IIso = IStart(1)+3
      If(CLine(ICur+1:Icur+1).eq.'-')  IIso = IStart(1)+2
      If(IIso.eq.0) then
       El=CLine(IStart(1):IStart(1)+1)
      Else
       El=CLine(IStart(1):IIso-2)
      EndIf
      call FilIAn(El,Number)
      If(Number.gt.-1) then
       IAnZ=Number
       If(IIso.ne.0) then
        call St2Dat(CLine,0,0,IIso,IType,CType,IVal,RVal,CVal,LVal,
     $    SVal)
        Isotp = IVal
       EndIf
       goto 10
      EndIf
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      if(CType.eq.'S') then
       If(IIso.eq.0) then
        El=CLine(IStart(1):IStart(1)+1)
       Else
        El=CLine(IStart(1):IIso-2)
       EndIf
       call FilIAn(El,Number)
        IAnZ = Number
       If(IIso.ne.0) then
        call St2Dat(CLine,0,0,IIso,IType,CType,IVal,RVal,CVal,LVal,
     $    SVal)
        Isotp = IVal
       endif
      elseif(CType.eq.'I') then
       IAnZ = IVal
      else
       Write(IOut,'('' Wrong Atomic Number for Atom'',I5)') IAt
       Stop
      endIf
   10 continue  
C
C read geometrical parameters
C
      If(IZAt.eq.1) goto 999 
      ICur = IStart(2)
      IVal = 0
      call LLinCl(SVal)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(Ctype.ne.'I'.or.ival.le.0.or.ival.gt.IZAt) then
       Write(IOut,'('' Wrong Connectivity for Atom'',I5)') IZAt
       Stop
      else
       IZIZ(1) = IVal
      endif
      ICur = IStart(3)  
      RVal = 0.0D0 
      call LLinCl(SVal)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(CType.eq.'S') then
       NBLVar=NBLVar+1
       ILnSL=Len(SVal)
       SL(1:ILnSL)=SVal(1:ILnSL)
       LBL=NBLVar
       C(1)=0.0d0
      ElseIf(Ctype.eq.'R') then 
       SL(1:1)='R'
       SL(2:Lenstr+1)=StrAdd(1:Lenstr)
       ILnSL=Lenstr+1
       LBL=0
       NBLFix=NBLFix+1
       C(1) = RVal
      Else
       Write(IOut,'('' Wrong Bond Length for Atom'',I5)') IZAt
       Stop
      EndIf
      IF(IZAt.eq.2) goto 999 
      ICur = IStart(4)
      IVal = 0
      call LLinCl(SVal)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(Ctype.ne.'I'.or.ival.le.0.or.ival.gt.IZAt) then
       Write(IOut,'('' Wrong Connectivity for Atom'',I5)') IZAt
       Stop
      else
       IZIZ(2) = IVal 
      endif
      ICur = IStart(5)
      RVal=0.0D0
      call LLinCl(SVal)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(Ctype.eq.'S') then
       NAlVar=NAlVar+1
       ILnSA=Len(SVal)
       SA(1:ILnSA)=SVal(1:ILnSA)
       LAlpha=NAlVar
       C(2)=0.0d0       
      elseif(CType.eq.'R') then
       SA(1:1)='A'
       call NumChar(IZAt-2,StrAdd,LenStr)
       SA(2:Lenstr+1)=StrAdd(1:LenStr)
       ILnSA=Lenstr+1
       LAlpha=0
       NAlFix=NAlFix+1
       C(2) = RVal
      else
       Write(IOut,'('' Wrong Valence Angle for Atom'',I5)') IZAt
       Stop
      endif
      If(IZAt.eq.3) goto 999 
      ICur = IStart(6)
      IVal = 0
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(Ctype.ne.'I'.or.ival.le.0.or.ival.gt.IZAt) then
       Write(IOut,'('' Wrong Connectivity for Atom'',I5)') IZAt
       Stop
      else
       IZIZ(3) =IVal
      EndIf
      ICur = IStart(7)
      RVal = 0.0d0
      call LLinCl(SVal)
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(Ctype.eq.'S') then
       NBeVar=NBeVar+1
       ILnSD=Len(SVal)
       SD(1:ILnSD)=SVal(1:ILnSD)
       LBeta=NBeVar
       C(3) = 0.0D0
       If(SVal(1:1).eq.'-') LBeta=-NZVar
      elseif(CType.eq.'R') then
       NBeFix=NBeFix+1
       C(3) = RVal
       SD(1:1)='D'
       call NumChar(IZAt-3,StrAdd,LenStr)
       SD(2:Lenstr+1)=StrAdd(1:Lenstr)
       ILnSD=LenStr+1
       LBeta=0
      else
       Write(IOut,'('' Wrong Dihedral Angle for Atom'',I5)') IZAt
       Stop
      endif
      ICur = IStart(8)
      IVal = 0
      call St2Dat(CLine,0,0,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      If(CType.eq.'I') IZIZ(4) = IVal 
  999 Continue
      Return
      End
*Deck RdVarZ
      Subroutine RdVarZ(In,IOut,IPrint,ToRad,NZ,NVar,BL,Alph,Bet,
     $  StBL,StAlph,StBet,NZFnd,LBL,LAlpha,LBeta,IStL,IStA,IStB,LZVar)
      Implicit Real*8 (A-H,O-Z)
      Integer St2Int,IType,IVal,IStart(10)
      Dimension BL(*),Alph(*),Bet(*)
      Character*4000 StBL,StAlph,StBet
      Character CType*1,Sval*20,StrInp*20,CLine*80,CLine1*80
      Dimension RValIn(1000)
      Dimension IStL(*),IStA(*),IStB(*)
      Real*8 RVal
      Complex*16 CVal
      Logical LVal
      Logical LZVar
      Dimension LBL(*),LAlpha(*),LBeta(*)
      NVar=0
   10 call LlinCl(CLine1)
      call LLinCl(CLine)
      Read(In,'(A80)') CLine1
      call LinUpC(CLine1,CLine)
      call LRmEqu(CLine)
      call SubStr(CLine,8,IStart,Nvalue)
      If(NValue.eq.0) goto 999
      NVar=NVar+1
      ICur=IStart(1)
      call St2Dat(CLine,0,100,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      if(CType.ne.'S') then
       write(IOut,'('' Ctype = '',A1,'' Error in Variable'',I5)') 
     $  CType,IFound
       Stop       
      else
       Ini=1
       Iend=Len(SVal)
       StrInp(Ini:IEnd)=SVal(Ini:IEnd)
      endif
      ICur=IStart(2)
      call St2Dat(CLine,0,100,ICur,IType,CType,IVal,RVal,CVal,LVal,SVal)
      if(CType.ne.'R') then
       write(IOut,'('' CType = '',A1,'' Error in Variable'',I5)')
     $   CType,IFound
       Stop
      EndIf
      RValIn(NVar)=RVal
      Do 20 II=1,NZ
       IniZL=IStL(II)+1
       InizA=IStA(II)+1
       InizB=IStB(II)+1
       IEndL=IStL(II+1)-1
       IEndA=IStA(II+1)-1
       IEndB=IStB(II+1)-1
       If(StrInp(Ini:IEnd).eq.StBL(InizL:IEndL)) then
        BL(II)=RValIn(NVar) 
        LBL(ii)=0
        If(LZVar) LBL(II)=NVar
        NZFnd=NZFnd+1
       ElseIf(StrInp(Ini:IEnd).eq.StAlph(InizA:IEndA)) then 
        Alph(II)=RValIn(NVar)*ToRad
        LAlpha(ii)=0
        If(LZVar) LAlpha(II)=NVar
        NZFnd=NZFnd+1
       ElseIf(StrInp(Ini:IEnd).eq.StBet(InizB:IEndB)) then
        Bet(II)=RValIn(NVar)*ToRad
        LBeta(ii)=0
        If(LZVar) LBeta(II)=NVar
        NZFnd=NZFnd+1
       ElseIf(STrInp(Ini:IEnd).eq.StBet(InizB+1:IEndB)) then
        If(StBet(InizB:InizB).eq.'-') then
         Bet(II)=-RValIn(NVar)*ToRad
         LBeta(ii)=0
         If(LZVar) LBeta(II)=-NVar
         NZFnd=NZFnd+1
        EndIf
       EndIf
   20 continue
      goto 10
  999 continue
C     write(IOut,'(/,'' Geometry in RdVarZ for'',I5,'' Atoms'')')NZ
C     do 30 II=1,NZ
C      write(IOut,'(3(F12.6,I3,2X))') BL(II),LBL(II),Alph(II)/ToRad,
C    $   LAlpha(II),Bet(II)/ToRad,LBeta(II)
C 30  Continue
      Return
      End
*Deck MkEAN
      Subroutine MKEAN(IOut,IPrint,MxBNd,NAtoms,IAN,NBond,IBond,C,
     $  Bonder,EAN,EANT)
      Implicit Real*8 (A-H,O-Z)
C The synthon is made by Delocalization (Del), Coordination (Crd) and 
C Rigidity (Rig)
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Dimension EAN(*),EANT(*),C(3,*),Bonder(MxBnd,*)
      Dimension teta0(3)
      Data zero /0.0d0/,one/1.0d0/,four/4.0d0/
      Data teta0/1.0947D+2,1.20d+2,1.80d+2/
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      if(IPrint.gt.1) write(IOut,'(/,''Bond Characteristics'')')
C Build EAN
      SynMax=zero
      SynMxT=zero
      Do 10 JAt=1,NAtoms
       NValJ=IAn(JAt)
       IF(IAn(JAt).gt.2) NValJ=NValJ-2
       If(IAn(JAt).gt.10) NValJ=NValJ-8
       NMaxJ=8-NValJ
       NEffJ=NMaxJ-NBond(JAt)+1
       if(NValJ.eq.1) then
        NEffJ=1
        t0j=teta0(3)
       else
        t0j=teta0(NEffJ)
       endif
       DelJ=one
       DelJT=one
       CrdJ=zero
       CrdJT=zero
       RigJ=zero
       NAngJ=0
       Do 20 IB=1,NBond(JAt)
        IAt=IBond(IB,JAt)
        Value=Distan(C,IAt,JAt,0)         
        Val0=RCovCT(IAn(IAt),IAn(JAt))
        BndOrd=exp((Val0-Value)/3.0d-01)
        BndOrT=Bonder(IB,JAt)
C for terminal atoms integer bond orders
        IBnd=Int(BndOrd+5.0d-1)
C       If(NBond(JAt).eq.1) BndOrd=Float(IBnd)
        if(IPrint.gt.1) write(IOut,'(''Atoms'',2I5,'' Distance'',f10.4,
     $    '' Bond Order'',f10.4)') IAt,JAt,Value,BndOrd
        DelJ=DelJ*BndOrd
        DelJT=DelJT*BndOrT 
        CrdJ=CrdJ+Float(IAn(IAt))*BndOrd
        CrdJT=CrdJT+Float(IAn(IAt))*BndOrT
        if(NBond(JAt).eq.1) go to 20
        Do 30 KB=IB,NBond(JAt)
         KAt=IBond(KB,JAt)
         if(IAt.eq.KAt) go to 30
         NAngJ=NAngJ+1
         Val=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
         Val1=val-t0j/todeg
         RigJ=RigJ+abs(sin(Val1))
   30   continue
   20  continue
       If(NAngJ.eq.0) NAngJ=1
       RigJ=RigJ/Float(NAngJ)
       DelJ=DelJ/Float(NBond(JAt))
       DelJT=DelJT/Float(NBond(JAt))
       CrdJ=CrdJ/Float(NBond(JAt))
       CrdJT=CrdJT/Float(NBond(JAt))
       SyntJ=CrdJ+RigJ+DelJ
       SyntJT=CrdJT+RigJ+DelJT 
       EAN(JAt)=SyntJ
       EANT(JAt)=SyntJT
       If(SyntJ.gt.SynMax) SynMax=SyntJ
       If(SyntJT.gt.SynMxT) SynMxT=SyntJT
       If(IPrint.gt.1) write(IOut,'('' IAN,NVal,NEff,NAng,T0,Teta,'',
     $   ''Coord,Del.,Rig.'',4I3,5F10.5)') IAn(JAt),NValJ,NEffJ,NAngJ,
     $   t0j,val*todeg,CrdJ,DelJ,RigJ
   10 continue
      if(IPrint.gt.0) write(IOut,'(''  Atom  IAN    EAN      EANT'')')
      do 40 IAt=1,NAtoms
       EAn(IAt)=Float(IAn(IAt))-0.495d0+EAN(IAT)/(SynMax+1.0d-1)
       EAnT(IAt)=Float(IAn(IAt))-0.495d0+EANT(IAT)/(SynMxT+1.0d-1)
       call ElNeg(IOut,IPrint,1,IAn(IAt),ELnA)
       if(IPrint.gt.0) write(IOut,'(2I5,2F10.5)')IAt,IAn(IAt),EAN(IAt),
     $  EAnT(IAt)
   40 continue
      return
      end
*Deck MkCNA
      Subroutine MkCNA(IOut,IPrint,Multip,MxBnd,NAtoms,IAn,NBond,IBond,
     $  NPi,NLP,C,ACN,ACN1,SPIJ)
      Implicit Real*8 (A-H,O-Z)
      Logical Error
      Parameter(MaxEl=200)
      Dimension IEl(0:MaxEl)
      Dimension C(3,*) 
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Dimension NPi(*),NLP(*)
      Dimension ACN(*),ACN1(*),SPIJ(*)
      Call FillEl(0,MaxEl,IEl)
      Call AClear(NAtoms,ACN)
      Call AClear(NAtoms,ACN1)
      Call AClear(NAtoms,SPIJ)
      Call IClear(NAtoms,NBond)
      Call IClear(NAtoms*MxBnd,IBond)
      Call IClear(NAtoms,NPi)
      Call IClear(NAtoms,NLP)
      ElnIJ=0.0D0
      Pt5=5.0D-1
      F1=1.0D0
      K1=1.6D1 
      Do 20 IAt=2,NAtoms
       IAI=IAn(IAt)
       Do 10 JAt=1,IAt-1
        IAJ=IAn(JAt)
        RIJ=Distan(C,IAt,JAt,0)
        R0IJ=RCovCT(IAI,IAJ)
C Grimme original
        REff=4.0D0*R0IJ/3.0D0
C Coherent with MkBnd
        REff = 1.15D0*R0IJ
C
        CNA=Pt5*(F1+DAbs(ElnIJ))*(F1+ERF((-K1*(RIJ-REff)/REff)))
        CNA1=F1/(F1+DExp(-K1*(REff/RIJ-F1)))
        PIJ=DExp((R0IJ-RIJ)/3.0D-1)
        if(PIJ.lt.CNA) PIJ=CNA
        ACN(IAt)=ACN(IAt)+CNA
        ACN(JAt)=ACN(JAt)+CNA
        ACN1(IAt)=ACN1(IAt)+CNA1
        ACN1(JAt)=ACN1(JAt)+CNA1 
        SPIJ(IAt)=SPIJ(IAt)+PIJ
        SPIJ(JAt)=SPIJ(JAt)+PIJ
        Test=CNA+0.2d0
        IF(Test.lt.-1.0d0) then
         write(IOut,'('' Distance of'',F10.5,'' between atoms'',2I5)') 
     $     RIJ,IAt,IAt
         stop
        ElseIf(Test.gt.1.0d0) then
         NBond(IAt)=NBond(IAt)+1
         NBond(JAt)=NBond(JAt)+1
         IBond(NBond(IAt),IAt)=JAt
         IBond(NBond(JAt),JAt)=IAt
        EndIf
   10  Continue
   20 Continue
       Write(IOut,'(''    Atom  Sigma  Pi  Lone Pairs  Unpaired'',
     $  ''  Bonded to'')')
      NTUp=0
      Do 30 IAt=1,NAtoms
       IAnI = IAn(IAt)
       NBI=NBond(IAt)
       ACNI = ACN(IAt)
       SPII = SPIJ(IAt)
       NSigma=Int(ACNI+1.0d-1)
       NSP=Int(SPII+4.0d-1)
       NPI(IAt)=NSP-NSigma
       NVal=IAnI
       If(IAnI.gt.2.and.IAnI.le.10) then
        NVal=NVal-2
       ElseIf(IAnI.gt.10.and.IAnI.le.18) then
        NVal=NVal-10
       EndIf
       If(IAnI.gt.1) then
        call SetLP(Multip,NVal,NSigma,NPi(IAt),NSP,NLP(IAt),NUp)
       Else
        NPI(IAt)=0
        NLP(IAt)=0
        NUP=0
       EndIf 
       NTUp=NTUp+NUP
       write(IOut,'(I5,2X,A2,I4,I5,I7,I12,2X,4I5)')IAt,IEl(IAn(IAt)),
     $   NSigma,NPI(IAt),NLP(IAt),NUP,(IBond(II,IAt),II=1,NBI)
   30 Continue
      If(Multip.ne.(NTup+1))Write(IOut,'('' Multiplicity: Expected ='',
     $  I2,'' Found ='',I2)') Multip,NTup+1
      Return
      End
*Deck SetLP
      Subroutine SetLP(Multip,NVal,NSigma,NPi,NSP,NLP,NUp) 
      NTest=Max(0,(NVal-NSIGMA-NPI))
      If(Multip.eq.1.and.NTest.eq.1) then
       NPI=NPI+1
       NTest=NTest-1
      EndIf
      NLP=Max(0,NTest)/2
      NUP=Max(0,(NTest-NLP*2))
      NTest1=NLP+NUP
      If(NTest1.gt.NTest) then
       NLP=NLP-1
       NUP=NUP+1
      EndIf 
      If(NSP.gt.NVal.and.NPI.gt.1) then
       NPI=NPI-2
       NUP=NUP+1
      EndIf 
      Return
      End
*Deck MkBnd
      Subroutine MkBnd(IOut,IPrint,MxBnd,NAtoms,IAn,NBond,IBond,C)
      Implicit Real*8 (A-H,O-Z)
      Implicit Integer (I-N)
      Dimension C(3,*) 
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Do 20 IAt=2,NAtoms
       IAI=IAn(IAt)
       Do 10 JAt=1,IAt-1
        IAJ=IAn(JAt)
        RIJ=Distan(C,IAt,JAt,0)
        RCvIJ=RCovCT(IAI,IAJ)
        Test=RIJ-RCvIJ
        IF(Test.lt.-1.0d0) then
         write(IOut,100) RIJ,IAt,IAt
         stop
        ElseIf(Test.lt.0.3d0) then
         NBond(IAt)=NBond(IAt)+1
         NBond(JAt)=NBond(JAt)+1
         IBond(NBond(IAt),IAt)=JAt
         IBond(NBond(JAt),JAt)=IAt
        EndIf
   10  Continue
   20 Continue
      If(IPrint.eq.0) Return
      write(IOut,200) NAtoms
      do 30 IAt=1,NAtoms
       write(IOut,300)IAt,NBond(IAt),(IBond(II,IAt),II=1,NBond(IAt))
   30 Continue
  100 Format('Distance of',F10.5,' between atoms',2I5)
  200 Format(' In MkBnd NAtoms =',I5)
  300 Format('Atom',I5,' forms',I2,' bonds with',6I5)
      Return
      End
*Deck FormZ
      Subroutine FormZ(IOut,IPrint,NAtoms,NZ,MaxNZ,IAnZ,C,IZ,Bl,Alpha,
     $  Beta,ZSymb,ZCheck,LBl,LAlpha,LBeta,Values,IUsed,ToAng,OK)
      Implicit Real*8(A-H,O-Z)
C
C     From the nuclear cordinates in C and the integers in the
C     Z-Matrix (IZ) build the Z-matrix.  If ZSymb is .false.,
C     all values are loaded into Bl, Alpha, and Beta.  If ZSymb
C     is .true., Values is loaded with the appropriate variables
C     instead.  If ZSymb is .true., ZCheck determines whether an
C     error message is printed if the constants in the symbolic
C     Z-Matrix are not compatible with the specified coordinates.
C     Note that the constants in a symbolic Z-Matrix are not
C     updated ... to have both constants and variables changed,
C     call this routine twice, with ZSymb .false. and then .true.
C     IUsed is an integer scratch vector 3*NZ long.  OK is returned
C     .true. if there are no problems.
C
C     Deal also with dummy atoms.
C
      Logical ZSymb, ZCheck, OK
      Real*8 MDCutO
      Dimension C(3,*), IZ(4,MaxNZ), Bl(*), Alpha(*), Beta(*), IAnZ(*),
     $  LBl(*), LAlpha(*), LBeta(*), Values(*), IUsed(*)
      Save One,Two,F45
      Data One,Two,F45/1.0D0,2.0D0,45.0d0/
 1000 Format(/,' Not coded to convert from cartesian to internal ',
     $          'coordinates when dummy atoms are used')
C
C     Define a statement function to calculate interatomic distances.
C
      Dist(I,J) = Sqrt((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     $                 +(C(3,I)-C(3,J))**2)
C
C     Get Z-matrix info, and check if we can handle it.
C
C      OK = NZ.eq.NAtoms
C      If(.not.OK) Write(IOut,1000)
C      If(.not.OK) Return
      Small = MDCutO(0)
      ToDeg = F45 / ATan(One)
C
C     Compute bond lengths.
C
      If(NZ.lt.2) Return
      Call IClear(3*NZ,IUsed)
      Do 140 IAt = 2, NZ
        BD = Dist(IAt,IZ(1,IAt))
        Call UpVar(ZSymb,ZCheck,IOut,1,IAt,LBl(IAt),BD,IUsed,Bl,Values,
     $    ToAng,OK)
        If(.not.OK) Return
  140   Continue
C
C     Compute bond angles (alpha) using the law of cosines.
C
      If(NZ.lt.3) Return
      Do 60 IAt = 3, NZ
        RIJ = Dist(IAt,IZ(1,IAt))
        RIK = Dist(IAt,IZ(2,IAt))
        RJK = Dist(IZ(1,IAt),IZ(2,IAt))
        Denom = Two * RIJ * RJK
        OK = Abs(Denom).ge.Small
        If(.not.OK) Return
        Arg = (RIJ*RIJ + RJK*RJK - RIK*RIK)  /  Denom
        Ang = ACos1(Arg)
        Call UpVar(ZSymb,ZCheck,IOut,2,IAt,LAlpha(IAt),Ang,IUsed,
     $    Alpha,Values,ToDeg,OK)
        If(.not.OK) Return
   60   Continue
C
C     Compute beta angles.  These may be dihedral (4-point) or normal
C     (3-point) angles.
C
      If(NZ.lt.4) Return
      Do 120 IAt = 4, NZ
        JAt = IZ(1,IAt)
        KAt = IZ(2,IAt)
        LAt = IZ(3,IAt)
C
C       Beta is a dihedral angle.  The dihedral angle is
C       IAt.JAt.KAt.LAt.
C
        If(IZ(4,IAt).eq.0) then
          Ang = Dihed(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))
C
C       Beta is a second bond angle.
C
        else
          RIJ = Dist(IAt,JAt)
          RIL = Dist(IAt,LAt)
          RJL = Dist(JAt,LAt)
          Denom = Two * RIJ * RJL
          OK = Abs(Denom).ge.Small
          If(.not.OK) Return
          Ang = ACos1((RIJ*RIJ+RJL*RJL-RIL*RIL)/Denom)
          endIf
        Call UpVar(ZSymb,ZCheck,IOut,3,IAt,LBeta(IAt),Ang,IUsed,
     $    Beta,Values,ToDeg,OK)
        If(.not.OK) Return
  120   Continue
      Return
      End
*Deck UpVar
      Subroutine UpVar(ZSymb,ZCheck,IOut,IType,IAt,IVar,X,IUsed,XArray,
     $  Values,Conv,OK)
      Implicit Real*8(A-H,O-Z)
C
C     This routine handles updating value X associated with atom
C     IAt and variable IVar.  It handles checking against constants
C     in the Z-matrix and for contradictory variable defininitions.
C
      Logical ZSymb, ZCheck, OK, OKX
      Dimension IUsed(*), XArray(*), Values(*)
      Character*2 Name(3)
      Save Name, Tol
      Data Name/'BL','AL','BE'/, Tol/1.D-6/
 1010 Format(' Atom ',I3,' needs constant  ',A,'=',F15.10,' but is ',
     $       F15.10)
 1020 Format(' Atom ',I3,' needs variable ',I3,'=',F15.10,' but is ',
     $       F15.10)
C
      If(.not.ZSymb) then
        XArray(IAt) = X
      else if(IVar.eq.0) then
        XW = X * Conv
        XAW = XArray(IAt) * Conv
        OKX = Abs(XArray(IAt)-X).lt.Tol
        OK = OK.and.OKX
        If(ZCheck.and..not.OKX)
     $    Write(IOut,1010) IAt, Name(IType), XAW, XW
      else if(IUsed(IAbs(IVar)).ne.0) then
        XTest = Values(IAbs(IVar))
        If(IVar.lt.0) XTest = -XTest
        OKX = Abs(XTest-X).lt.Tol
        OK = OK.and.OKX
        XW = X * Conv
        XAW = XTest * Conv
        If(ZCheck.and..not.OKX) Write(IOut,1020) IAt, IVar, XW, XAW
        XArray(IAt) = XTest
      else
        IUsed(IAbs(IVar)) = 1
        If(IVar.gt.0) Values(IVar) = X
        If(IVar.lt.0) Values(-IVar) = -X
        XArray(IAt) = X
        endIf
      Return
      End
*Deck MkMpZA
      Subroutine MkMpZA(IOut,NZ,IAnZ,NAtoms,MapZAt)
      Implicit Integer(A-Z)
C    
C     Generate a map from Z-matrix numbers to atom numbers.
C     
      Dimension IAnZ(*), MapZAt(*)
C     
      NAt = 0
      Do 10 IZ = 1, NZ
        If(IAnZ(IZ).ge.0)then
          NAt = NAt + 1
          MapZAt(IZ) = NAt
        else
          MapZAt(IZ) = 0
          endIf
   10   Continue
      If(NAt.ne.NAtoms) then
       write(IOut,'('' Consistency failure in MkMpZA'')')
       Stop
      EndIf
      Return
      End
*Deck InToCh
      Subroutine InToCh(Num,Str,ILen)
      Implicit Real*8 (A-H,O-Z)
C
C Form string Str with length 4 from number Num
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
      Str(1:1)=SNum(IThous+1:IThous+1) 
      Str(2:2)=SNum(IHund+1:IHund+1)
      Str(3:3)=SNum(ITen+1:ITen+1)
      Str(4:4)=SNum(IUn+1:IUn+1)
      ILen=4
      Return
      End
*Deck BldZmt
      Subroutine BldZmt(IOut,IPrint,NAtoms,MxBnd,IAn,NBond,IBond,C,IZ,
     $  LBL,LAlpha,LBeta,IANZ,CZ,EAN,EANZ,OK)
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*),CZ(3,*),EAN(*),EANZ(*)
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Dimension IAnZ(*),IZ(4,*),LBl(*),LAlpha(*),LBeta(*)
      Dimension IndZ(1000),IndI(1000) 
      Logical OK,FndZ(1000),FndI(1000)
      OK=.false.
      NItMax=NAtoms+1
C
C INdZ is the Z numbering in terms of Input numbering
C INdI is the input numbering in terms of Z numbering
C FndZ,FndI are the corresponding logical variables 
C
      Do 10 IAt=1,NAtoms
       IndZ(IAt)=IAt 
       IndI(IAt)=IAt
       FndZ(IAt)=.false.
       FndI(IAt)=.false.
       IAnZ(IAt)=IAn(IAt) 
       IZ(1,IAt)=0
       IZ(2,IAt)=0
       IZ(3,IAt)=0
       IZ(4,IAt)=0
   10 Continue
      Call AMove(3*NAtoms,C,CZ)
      Call AMove(NAtoms,EAN,EANZ)
      Call IniZmt(IOut,IPrint,MxBnd,NAtoms,NBond,IBond,IndI,IndZ,IZ,IAn,
     $  IAnZ,EAn,EAnZ,FndI,FndZ)
      NAtZ=Min0(4,NAtoms)
      NIter=0
      if(IPrint.gt.0) write(IOut,'('' Input Z-Mat. At.Numb. IZ(1,I)'',
     $  '' IZ(2,I) IZ(3,I)'')')
      do 20 IAtZ=1,NAtZ
       II=IndI(IAtZ)
       Call AMove(3,C(1,II),CZ(1,IAtZ))
       if(IAtZ.eq.1) then
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8)') II,IAtZ,IANZ(IAtZ)
       ElseIf(IAtz.eq.2) then
        LBL(IAtZ)=1
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8,I8)') II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ)
       ElseIf(IAtZ.eq.3) then
        LBL(IAtZ)=2
        LAlpha(IAtZ)=3
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8,2I8)') II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ),IZ(2,IAtZ)
       Else
        LBL(IAtZ)=(IAtZ-3)*3+1
        LAlpha(IAtZ)=LBL(IAtZ)+1
        LBeta(IAtZ)=LAlpha(IAtZ)+1
        If(IPrint.gt.0) write(IOut,'(I5,I6,I8,3I8)')II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ),IZ(2,IAtZ),IZ(3,IAtZ)
       EndIf
   20 continue          
      If(NAtoms.le.4) Return  
   30 OK=.false.
      NIter=NIter+1 
      call SetAtom(NAtoms,IAt,NAtZ,MxBnd,NBond,IBond,IndI,IndZ,FndI,
     $  FndZ,IAn,IAnZ,IZ,EAN,EANZ,OK)
      If(NAtZ.lt.NAtoms.and.NIter.lt.NItMax) goto 30  
      If(.not.OK) then 
       write(IOut,'(/,'' After'',I5,'' Iterations'')') NIter
       write(IOut,'(I5,'' Atoms Assigned Over'',I5)') NAtZ,NAtoms
       Stop
      EndIf 
C     write(IOut,'('' Input Z-Mat. At.Numb. IZ(1,I) IZ(2,I) IZ(3,I)'')')
      do 40 IAtZ=5,NAtZ  
       II=IndI(IAtZ)
       Call AMove(3,C(1,II),CZ(1,IAtZ))
       if(IAtZ.eq.1) then
        if(IPrint.gt.0)write(IOut,'(I5,I6,I8)') II,IAtZ,IANZ(IAtZ)
       ElseIf(IAtz.eq.2) then
        LBL(IAtZ)=1
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8,I8)') II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ) 
       ElseIf(IAtZ.eq.3) then
        LBL(IAtZ)=2
        LAlpha(IAtZ)=3
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8,2I8)') II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ),IZ(2,IAtZ)
       Else
        LBL(IAtZ)=(IAtZ-3)*3+1
        LAlpha(IAtZ)=LBL(IAtZ)+1
        LBeta(IAtZ)=LAlpha(IAtZ)+1 
        if(IPrint.gt.0) write(IOut,'(I5,I6,I8,3I8)')II,IAtZ,IANZ(IAtZ),
     $   IZ(1,IAtZ),IZ(2,IAtZ),IZ(3,IAtZ)
       EndIf
   40 continue   
      Write(IOut,'(/,13X,''Z-Matrix Built for'',I5,'' Atoms After'',I3,
     $  '' Iterations'')') NatZ,NIter
      Return
      End
*Deck SetAtom
      Subroutine SetAtom(NAtoms,IAt,NAtZ,MxBnd,NBond,IBond,IndI,IndZ,
     $  FndI,FndZ,IAn,IAnZ,IZ,EAN,EANZ,OK)
      Implicit Real*8 (A-H,O-Z)
      Logical OK,FndI(*),FndZ(*)
      Dimension NBond(*),IBond(MxBnd,*),IndI(*),IndZ(*),IAn(*),IAnZ(*)
      Dimension IZ(4,*)
      Dimension EAN(*),EANZ(*)
      OK=.False.
      Do 10 IAt=1,NAtoms
       If(FndI(IAt)) goto 10
       NBI=NBond(IAt)
       Do 20 J1=1,NBI
        If(OK) goto 20
        JAt=IBond(J1,IAt)
        If(.not.FndI(JAt)) goto 20
        NAtZ=NAtZ+1
        IndI(NAtZ)=IAt
        IndZ(IAt)=NAtZ
        FndZ(NAtZ)=.true.
        FndI(IAt)=.true.
        IAnZ(NAtZ)=IAn(IAt)
        EANZ(NAtZ)=EAN(IAt)
        JAtZ=IndZ(JAt)
        IZ(1,NAtZ)=JAtZ
        IZ(2,NAtZ)=IZ(1,JAtZ)
        IZ(3,NAtZ)=IZ(2,JAtZ)
        OK=.true.
   20  Continue
   10 Continue
      Return
      End
*Deck SymZMt
      Subroutine SymZMt(IOut,IPrint,NAtoms,NVar,IZ,LBL,LAlpha,LBeta,
     $  EANZ)
      Implicit Real*8 (A-H,O-Z)
      Logical EqBL,EqAlph
      Dimension IZ(4,*),LBL(*),LAlpha(*),LBeta(*),IEq(1000)
      Dimension EANZ(*)
      Thresh=1.0d-4
      If(NAtoms.eq.1) return
      NVarL=1
      LBL(2)=1
      If(NAtoms.eq.2) return
      If(Abs(EANZ(3)-EANZ(1)).lt.thresh) then
       LBL(3)=-1 
       LAlpha(3)=1
       NVARA=1
      Else
       LBL(3)=2
       LAlpha(3)=1
       NVArL=2
       NVARA=1
      EndIf
      If(NAtoms.eq.3) return
      if(IPrint.gt.0) write(IOut,'(/,'' Equivalent atoms'')')
      IEq(1)=1
      IAt=1
      if(IPrint.gt.0) write(IOut,'(I5,F10.5,I5)') IAt,EANZ(IAt),IEq(IAt)
      Do 100 IAt=2,NAtoms
       IEq(IAt)=0
       Do 110 JAt=1,IAt-1
        If(IEq(IAt).ne.0) goto 110
        If(ABs(EANZ(JAt)-EANZ(IAt)).lt.thresh) IEq(IAt)=JAt
  110  Continue
       If(IEq(IAt).eq.0) IEq(IAt)=IAt
       If(IPrint.gt.0) write(IOut,'(I5,F10.5,I5)')IAt,EANZ(IAt),IEq(IAt)
  100 Continue
      Do 10 IAt=4,Natoms
       JAt=IEq(IZ(1,IAt))
       KAt=IEq(IZ(2,IAt))
       LAt=IEq(IZ(3,IAt))
       LBl(IAt)=0
       LAlpha(IAt)=0
       LBeta(IAt)=0
       Do 20 ITst=1,IAt-1
        JTst=IEq(IZ(1,ITst))
        KTst=IEq(IZ(2,ITst))
        LTst=IEq(IZ(3,ITst)) 
        If(EqBL(IAt,ITst,IEq,IZ,EANZ)) LBL(IAT)=-Abs(LBL(ITst))
        If(EqAlph(IAt,ITst,IEq,IZ,EANZ)) LAlpha(IAt)=-Abs(LAlpha(ITst))
C       If(ABs(EANZ(LAt)-EANZ(LTst)).lt.thresh) LBeta(IAt)=
C    $    -Abs(LBeta(ITst))
   20  Continue
       If(LBL(IAt).eq.0) then
        NVarL=NVarL+1
        LBL(IAt)=NVarL
       EndIf
       If(LAlpha(IAt).eq.0) then
        NVarA=NVarA+1
        LAlpha(IAt)=NVarA
       EndIf
       If(LBeta(IAt).eq.0) then
        NVarB=NVarB+1
        LBeta(IAt)=NVarB
       EndIf
   10 Continue
      If(IPrint.gt.0)write(IOut,'(/,'' From SymZmt NVar ='',I5)')NVar 
      Return
      End
*Deck EqBL
      Logical Function EqBL(IAt,JAt,IEq,IZ,EANZ)
      Implicit Real*8 (A-H,O-Z)
      Logical OK1,OK2
      Dimension IEq(*),IZ(4,*)
      Dimension EANZ(*)
      thresh=1.0d-4
      EqBL=.false.
      IRef=IEq(IAt)
      IRef1=IEq(IZ(1,IAt))
      ITst=IEq(JAt) 
      ITst1=IEq(IZ(1,JAt))
      OK1=(ABS(EAnZ(IRef)-EANZ(ITst)).lt.thresh)
      if(OK1) then
       OK2=(ABS(EANZ(IRef1)-EANZ(ITst1)).lt.thresh)
      Else
       OK1=(ABS(EANZ(IRef)-EANZ(ITst1)).lt.thresh)
       OK2=(ABS(EANZ(IRef1)-EANZ(ITst)).lt.thresh)
      EndIf
      if(OK1.and.OK2) EqBL=.true.
      Return
      End
*Deck EqAlph
      Logical Function EqAlph(IAt,JAt,IEq,IZ,EANZ)
      Implicit Real*8 (A-H,O-Z) 
      Dimension IEq(*),IZ(4,*)
      Dimension EANZ(*)
      Thresh=1.0d-4
      EqAlph=.true.
      Tst1=ABS(EAnZ(IAT)-EANZ(JAt))+ABS(EANZ(IZ(2,IAt))
     $  -EANZ(IZ(2,JAt)))
      Tst2=ABs(EANZ(IAt)-EANZ(IZ(2,JAt)))+ABS(EANZ(IZ(2,IAt))
     $  -EANZ(JAt))
      Tst3=ABs(EANZ(IZ(1,IAt))-EANZ(IZ(1,JAt)))
      If(Tst1.gt.thresh.and.Tst2.gt.thresh) EqAlph=.false.
      If(Tst3.gt.thresh) EqAlph=.false.
      Return
      End
*Deck IniZmt
      Subroutine IniZmt(IOut,IPrint,MxBnd,NAtoms,NBond,IBond,IndI,IndZ,
     $  IZ,IAn,IAnZ,EAn,EAnZ,FndI,FndZ)
      Implicit Real*8 (A-H,O-Z)
      Dimension NBond(*),IBond(MxBnd,*),IndI(*),IndZ(*),IZ(4,*)
      Dimension IAn(*),IAnZ(*)
      Logical FndKAt,FndI(*),FndZ(*)
      Dimension EAN(*),EAnZ(*)      
C First Atom
      IAt=1
      IndI(1)=IAt
      IndZ(IAt)=1
      IAnZ(1)=IAn(IAt)
      EANZ(1)=EAN(IAt)      
      FndZ(1)=.true.
      FndI(IAt)=.true.
      If(NAtoms.eq.1) Return
C Second Atom    
      JAt=IBond(1,IAt)
      If(NBond(JAt).eq.1) then
       If(NBond(IAt).eq.1.and.Natoms.gt.2) then
        write(IOut,'('' Error in IniZmt for Atoms'',2I5)')IAt,JAt       
       EndIf
       II=JAt
       JAt=IAt
       IAt=II
      EndIf
      IndI(1)=IAt
      IndI(2)=JAt
      IndZ(IAt)=1
      IndZ(JAt)=2
      IZ(1,2)=1
      IAnZ(1)=IAn(IAt)
      IAnZ(2)=IAn(JAt)
      EANZ(1)=EAN(IAt)
      EANZ(2)=EAN(JAt)
      FndZ(1)=.true.
      FndZ(2)=.true.
      FndI(IAt)=.true.
      FndI(JAt)=.true.
      If(NAtoms.eq.2) Return      
C Third Atom
      If(NBond(IAt).gt.2) then
       JAt=IBond(1,IAt)
      ElseIf(NBond(IAt).eq.2) then
       JAt=IAt
       IAt=IBond(1,JAt)
      ElseIf(NBond(IAt).eq.1) then
       JAt=IBond(1,IAt)
       If(NBond(JAt).eq.1) then
        write(IOut,'(/,'' Error in  IniZmt for Atoms'',2I5)') IAt,JAt
        Stop
       EndIf
      EndIf
      FndKAt=.false.
      Do 10 JJ=1,NBond(JAt)
       If(FndKAt) goto 10
       If(NBond(JJ).eq.1) goto 10
       IF(IBond(JJ,JAt).eq.IAt) go to 10
       KAt=IBond(JJ,JAt)
       If(NBond(KK).eq.1) goto 10
       FndKAt=.true.
   10 Continue 
      If(.not.FndKAt.and.NAtoms.gt.3) then
       write(IOut,'('' Error in IniZmt for Atoms'',2I5)')IAt,JAt
       Stop
      EndIf 
      IndI(1)=IAt
      IndI(2)=JAt
      IndI(3)=KAt
      IndZ(IAt)=1
      IndZ(JAt)=2
      IndZ(KAt)=3
      IZ(1,2)=1
      IZ(1,3)=2
      IZ(2,3)=1
      IAnZ(1)=IAn(IAt)
      IAnZ(2)=IAn(JAt)
      IAnZ(3)=IAn(KAt)
      EANZ(1)=EAN(IAt)
      EANZ(2)=EAN(JAt)
      EANZ(3)=EAN(KAt)
      FndZ(1)=.true.
      FndZ(2)=.true.
      FndZ(3)=.true.
      FndI(IAt)=.true.
      FndI(JAt)=.true.
      FndI(KAt)=.true.
      If(NAtoms.eq.3) Return
C Fourth Atom
      If(NBond(KAt).eq.1) then
       If(NBond(JAt).gt.2) then
        LAt=IBond(3,JAt)
       ElseIf(NBond(IAt).gt.2) then
        LAt=IBond(3,IAt)
       Else
        write(IOut,'(/,'' Error in IniZmt for Atoms'',3I5)')IAt,JAt,KAt 
        Stop       
       EndIf
      Else
       LAt=IBond(2,KAt)
       If(LAt.eq.IAt.or.LAt.eq.JAt) then
        LAt=IBond(1,KAt)
        If(LAt.eq.IAt.or.LAt.eq.JAt) then
         IF(NBond(KAt).ge.3) then
          LAt=IBond(3,KAt)
         Else
          write(IOut,'(/,'' Error in IniZmt for Atoms'',3I5)') 
     $       IAt,JAt,KAt
         EndIf       
        EndIf 
       EndIf
      EndIf 
C Assign ZMat parameters
      IndI(1)=IAt
      IndI(2)=JAt
      IndI(3)=KAt
      IndI(4)=LAt
      IndZ(IAt)=1
      IndZ(JAt)=2
      IndZ(KAt)=3
      IndZ(LAt)=4
      IZ(1,1)=2
      IZ(2,1)=3
      IZ(1,2)=1
      IZ(2,2)=3
      IZ(1,3)=2
      IZ(2,3)=1
      IZ(1,4)=3
      IZ(2,4)=2
      IZ(3,4)=1
      IAnZ(1)=IAn(IAt)
      IAnZ(2)=IAn(JAt)
      IAnZ(3)=IAn(KAt)
      IAnZ(4)=IAn(LAt)
      EANZ(1)=EAN(IAt)
      EANZ(2)=EAN(JAt)
      EANZ(3)=EAN(KAt)
      EANZ(4)=EAN(LAt)
      FndZ(1)=.true.
      FndZ(2)=.true.
      FndZ(3)=.true.
      FndZ(4)=.true.
      FndI(IAt)=.true.
      FndI(JAt)=.true.
      FndI(KAt)=.true.
      FndI(LAt)=.true.
      Return
      End 
*Deck CtoZ
      Subroutine CtoZ(IOut,IPrint,MxBnd,NAtoms,IAn,NBond,IBond,IZ,
     $ LBL,LAlpha,LBeta,IAnZ,C,EAN,CZ,EANZ,BL,Alpha,Beta,Scr,IScr)
      Implicit Real*8 (A-H,O-Z)
      Logical OK
      Character Str1*6,Str2*6,Str3*6 
      Dimension NBond(*),IBond(MxBnd,*),IZ(4,*),IAN(*),IANZ(*)
      Dimension LBL(*),LAlpha(*),LBeta(*),IScr(*)
      Dimension C(3,*),CZ(3,*),EAN(*),EANZ(*),BL(*),Alpha(*),Beta(*)
      Dimension Scr(*)
C Conversion Factors
      ToAng=1.0d0
      Pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
C Build Z-matrix
      call BldZmt(IOut,IPrint,NAtoms,MxBnd,IAn,NBond,IBond,C,IZ,LBl,
     $  LAlpha,LBeta,IANZ,CZ,EAN,EANZ,OK)
      call SymZMt(IOut,IPrint,NAtoms,NVar,IZ,LBL,LAlpha,LBeta,EANZ)
      NZ=NAtoms
      if(IPrint.gt.0)write(IOut,'(/,'' From CtoZ: Cartesian'',
     $  '' Coordinates'')')
      do 88 IZAt=1,NZ
       if(IPrint.gt.0)write(IOut,'(I5,4F12.5)')IAn(IZAt),EAn(IZAt),
     $  (C(II,IZAt),II=1,3)
   88 continue
C Assign variables to Z-matrix
      call FormZ(IOut,IPrint,NAtoms,NZ,MaxNZ,IAnZ,CZ,IZ,Bl,Alpha,
     $  Beta,.False.,.False.,LBl,LAlpha,LBeta,Scr,IScr,ToAng,OK)
C write Z-matrix
      NZ=NAtoms 
      call PrtZmt(IOut,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $  LAlSh,LBeSh)
      call PrtVrZ(IOut,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $  LAlSh,LBeSh)
      Return
      End
*Deck PrtZmt
      Subroutine PrtZmt(IOut,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $  LAlSh,LBeSh)
      Implicit Real*8 (A-H,O-Z)
      Character Str1*6,Str2*6,Str3*6
      Logical FndA,FndB
      Dimension IANZ(*),IZ(4,*),LBL(*),LAlpha(*),LBeta(*)
      Dimension BL(*),Alpha(*),Beta(*)
      Dimension IEl(-2:204)
      call IClear(206,IEl)
      MaxEl=204
      call FillEl(-2,MaxEl,IEl) 
      LBlFix=0
      LAlFix=0
      LBeFix=0
      IAlSh=0
      IBeSh=0
      FndA=.false.
      FndB=.false.
      Do 5 II=2,NZ
       If(FndA) goto 5
       If(LAlpha(II).ne.0) then
        LAlSh=LAlpha(II)-1
        FndA=.true.
       EndIf
    5 Continue  
      Do 7 II=3,NZ
       If(FndB) goto 7
       If(LBeta(II).ne.0) then
        LBeSh=LBeta(II)-1
        FndB=.true.
       EndIf  
    7 Continue    
      Do 10 II=1,NZ
       IDx=IAnZ(II)
       If(II.gt.1) then
        If(LBL(II).ne.0) then
         ISI=Abs(LBl(II))
         Str1(1:1)=' '
         Str1(2:2)='R'
        Else
         LBlFix=LBlFix+1
         ISI=LBLFix
         Str1(1:1)='#'
         Str1(2:2)='F'
        EndIf
        Call InToCh(ISI,Str1(3:6),ILen1)
       EndIf
       If(II.gt.2) then
        If(LAlpha(II).ne.0) then
         ISI=Abs(LAlpha(II)-LAlSh)
         Str2(1:1)=' '
         Str2(2:2)='A'
        Else
         LAlFix=LAlFix+1
         ISI=LAlFix
         Str2(1:1)='#'
         Str2(2:2)='G'
        EndIf
        Call InToCh(ISI,Str2(3:6),ILen2)
       EndIf
       If(II.gt.3) then
        If(LBeta(II).ne.0) then
         ISI=Abs(LBeta(II)-LBeSh)
         Str3(1:1)=' '
         Str3(2:2)='D'
        Else
         LBeFix=LBeFix+1
         ISI=LBeFix
         Str3(1:1)='#'
         Str3(2:2)='H'
        EndIf
        Call InToCh(ISI,Str3(3:6),ILen3)
       EndIf
       If(II.eq.1) then
        write(IOut,'(1X,A2)') IEl(IDx)
       elseif(II.eq.2) then
        write(IOut,'(1X,A2,I5,1X,A6)') IEl(IDx),IZ(1,2),Str1
       elseif(II.eq.3) then
        write(IOut,'(1X,A2,2(I5,1X,A6))')IEl(IDX),IZ(1,3),Str1,
     $    IZ(2,3),Str2
       else
        write(IOut,'(1X,A2,3(I5,1X,A6),I5)')IEl(IDx),IZ(1,II),Str1,
     $   IZ(2,II),Str2,IZ(3,II),Str3,IZ(4,II)
       endif
  10  Continue
      write(IOut,'('' '')')
      Return
      End
*Deck PrtVrZ
      Subroutine PrtVrZ(IOut,NZ,IANZ,IZ,LBl,LAlpha,LBeta,BL,Alpha,Beta,
     $  LAlSh,LBeSh)
      Implicit Real*8 (A-H,O-Z)
      Dimension IANZ(*),IZ(4,*),LBL(*),LAlpha(*),LBeta(*)
      Dimension BL(*),Alpha(*),Beta(*)
      Character Str1*6,Str2*6,Str3*6
      Logical PrtBL,PrtAl,PrtBe
C Conversion Factors
      ToAng=1.0d0
      Pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
C
      LBlFix=0
      LAlFix=0
      LBeFix=0
      do 10 II=2,NZ
       If(II.gt.1) then
C       If(LBL(II).gt.0) then
        If(LBL(II).ne.0) then
         PrtBL=.true.
         do 15 JJ=1,II-1
          if(ABs(LBL(II)).eq.Abs(LBL(JJ))) PrtBL=.false.
   15    continue      
         ISI=Abs(LBl(II))
         Str1(1:1)=' '
         Str1(2:2)='R'
        Else
         PrtBL=.true.       
         LBlFix=LBlFix+1
         ISI=LBLFix
         Str1(1:1)=' '
         Str1(2:2)='F'
        EndIf
        Call InToCh(ISI,Str1(3:6),ILen1)
       EndIf
       If(II.gt.2) then
C       If(LAlpha(II).gt.0) then
        If(LAlpha(II).ne.0) then
         PrtAl=.true.
         do 25 JJ=1,II-1
          if(Abs(LAlpha(II)).eq.Abs(Lalpha(JJ))) PrtAl=.false.
   25    continue
         ISI=Abs(LAlpha(II)-LAlSh)
         Str2(1:1)=' '
         Str2(2:2)='A'
        Else
         PrtAl=.true.
         LAlFix=LAlFix+1
         ISI=LAlFix
         Str2(1:1)=' '
         Str2(2:2)='G'
        EndIf
        Call InToCh(ISI,Str2(3:6),ILen2)
       EndIf
       If(II.gt.3) then
        If(LBeta(II).ne.0) then
         PrtBe=.true.
         do 35 JJ=1,II-1
          if(LBeta(II).eq.LBeta(JJ)) PrtBe=.false.
   35    continue         
         ISI=Abs(LBeta(II)-LBeSh)
         Str3(1:1)=' '
         Str3(2:2)='D'
        Else
         PrtBe=.true.
         LBeFix=LBeFix+1
         ISI=LBeFix
         Str3(1:1)=' '
         Str3(2:2)='H'
        EndIf
        Call InToCh(ISI,Str3(3:6),ILen3)
       EndIf
       If(PrtBL) write(IOut,'(A6,'' ='',F10.5)') Str1,BL(II)
       If(II.gt.2.and.PrtAl) write(IOut,'(A6,'' ='',F10.5)') 
     $  Str2,Alpha(II)/ToRad
       If(II.gt.3.and.PrtBe) write(IOut,'(A6,'' ='',F10.5)') 
     $  Str3,Beta(II)/ToRad
  10  continue
C     write(IOut,'('' '')') 
      Return
      End
*Deck MkBPCS
      Subroutine MkBPCS(IOut,IPrint,NZ,IANZ,IZ,Bl)
      Implicit Real*8 (A-H,O-Z)
      Dimension IANZ(*),IZ(4,*),Bl(*)
      Save Zero,One,Two,Three,CVK,CHyp
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/
      Data CVK/1.1D-3/, HypK/2.5d-02/
      write(IOut,'(/,'' Bonded Atoms  rDSD Distance BDPCS3 Distance'')')
      Do 10 IAt=2,NZ
       BlOr=Bl(IAt)
       IAnI=IAnz(IAt)
       JAt=IZ(1,IAt)
       IAnJ=IAnZ(JAt)
       If(IANI.le.0.or.IAnJ.le.0) go to 10
       DltO=zero
       DltF=zero
       If(IAnI.eq.8.or.IAnJ.eq.8) DltO=one
       If(IAnI.eq.9.or.IAnJ.eq.9) DltF=one
       ZeffI=two
       ZeffJ=two
       If(IAnI.lt.3) ZeffI=one
       If(IAnJ.lt.3) ZeffJ=one
       If(IAnI.gt.10) ZeffI=three
       If(IAnJ.gt.10) ZeffJ=three
       Value=Blor
       Val0=RCovCT(IAnZ(IAt),IAnZ(JAt))
       BndOrd=exp((Val0-Value)/3.0d-01)
       Trm1=-CVK*sqrt(ZeffI*ZeffJ-one)*val0
       Trm2=SQrt(Abs(BndOrd-Two+DltO*(one-BndOrd)))
       Trm3=-HypK*DltF*(BndOrd-one)**2
       DltR=Trm1*Trm2+Trm3
C Do not apply the correction to H-Bonds
       If(BndOrd.lt.3.0d-01) DltR=0.0d0
       Bl(IAt)=BlOr+DltR
       write(IOut,'(2I5,7X,F8.5,7X,F8.5)') IAt,JAt,BlOr,Bl(IAt)
   10 continue
      write(IOut,'('' '')')
      return
      end
*Deck HedPrt
      Subroutine HedPrt(IOut,Width,String,Num)
      Implicit Integer(A-Z)
C
C     Print a header for a matrix.  If Width is >0, the
C     header is centered in a line of length Width.  If
C     Width=-1 or -2, then real/imaginary is appended to
C     the label.
C
      Character String*(*), Blank*80, LabRI*20
 1000 Format(A,A,A,I10)
 1010 Format(A,A,A)
 1020 Format(A,A,A,I10,A,':')
 1030 Format(A,A,A,A,':')
C
      Blank = ' '
      LStr = Max(LinEnd(String),1)
      If(Width.le.0) then
        LBl = 1
      else if(Num.gt.0) then
        LBl = (Width-LStr-10)/2 + 1
      else
        LBl = (Width-LStr)/2 + 1
        endIf
      If(Width.eq.-1) then
        LabRI = '(real)'
      else if(Width.eq.-2) then
        LabRI = '(imag)'
      else
        LabRI = ' '
        endIf
      LLabRI = LinEnd(LabRI)
      LBl = Max(LBl,1)
      If(Width.ge.0.and.Num.gt.0) then
        Write(IOut,1000) Blank(1:LBl),String(1:LStr),Blank(1:LBl),Num
      else if(Width.ge.0) then
        Write(IOut,1010) Blank(1:LBl), String(1:LStr), Blank(1:LBl)
      else if(Num.gt.0) then
        Write(IOut,1020) Blank(1:LBl),String(1:LStr),Blank(1:LBl),Num,
     $    LabRI(1:LLabRI)
      else
        Write(IOut,1030) Blank(1:LBl),String(1:LStr),Blank(1:LBl),
     $    LabRI(1:LLabRI)
        endIf
      Return
      End
*Deck LTOutC
      Subroutine LTOutC(IOut,IOpt,NRI,IRI,NSpin,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a matrix:
C
C     IOpt .. Output format:
C             1 ... D13.6
C             2 ... F13.4
C     NRI ... Number of components (1 for real, 2 for complex).
C     IRI ... Which component to print here (1=real, 2=imaginary).
C             Negative for anti-Hermetian matrices.
C     NSpin . -1 for just large component, 1 for just one spin component,
C             -2 for l and s components, 2 for a and b components.
C     N   ... Dimension of matrix.
C     A   ... Array to be printed.
C     Key ... <0: Suppress elements with absolute values less than
C                 10**(-6+Key)
C              1: Print complete matrix.
C
      Parameter (NumCol=5)
      Character*1 LabX(NumCol)
      Dimension A(NRI,*), S(NumCol), IRX(NumCol)
      Save Zero
      Data Zero/0.0d0/
 1010 Format(4X,10(7X,I6,A1))
 1020 Format(I7,A,D13.6,9D14.6)
 1030 Format(I7,A,F13.4,9F14.4)
C
      If(IOpt.lt.1.or.IOpt.gt.2) then
       write(IOut,'(/,''Illegal IOpt in LTOutC.'')')
       Stop
      EndIf
       Call PrtThr(Key,Thresh)
      NDo = IAbs(NSpin)*N
      IRIA = IAbs(IRI)
      Do 100 IStart = 1, NDo, NumCol
       IEnd = Min(IStart+NumCol-1,NDo)
       NCol = IEnd - IStart + 1
       Call LabSpn(NSpin,IStart-1,NCol,IRX,LabX)
       Write(IOut,1010) (IRX(I),LabX(I),I=1,NCol)
       Do 200 IRow = IStart, NDo
        ILim = Min(IRow-IStart+1,NumCol)
        L = (IRow*(IRow-1))/2 + IStart
        Do 10 I = 1, ILim
         S(I) = A(IRIA,L+I-1)
         If(Key.ne.1.and.Abs(S(I)).lt.Thresh) S(I) = Zero
   10   Continue
        If((IRI.eq.-1.or.IRI.eq.2)) then
         ILim1 = Min(ILim,IRow-IStart)
         Do 20 I = 1, ILim1
          S(I) = -S(I)
   20    Continue
        endIf
        Call LabSpn(NSpin,IRow-1,1,IRX,LabX)
        If(IOpt.eq.1) then
         Write(IOut,1020) IRX(1), LabX(1), (S(I),I=1,ILim)
        else
         Write(IOut,1030) IRX(1), LabX(1), (S(I),I=1,ILim)
        endIf
  200  Continue
  100 Continue
      Return
      End
*Deck LTOut
      Subroutine LTOut(IOut,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a symmetric matrix stored in compressed lower triangular form.
C
C        N         Dimension of matrix.
C        A         Array to be printed.
C        KEY    N<=0 ... Suppress elements with absolute values less
C                        than 10**(-6+N)
C                  1 ... Print complete matrix.
C
      Dimension A(*)
C
      Call LTOutC(IOut,1,1,1,1,N,A,Key)
      Return
      End
*Deck LTOutS
      Subroutine LTOutS(IOut,String,Num,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Print a heading, followed by a lower triangular matrix.
C     If Num is positive, it is also printed.
C
      Character*(*) String
      Dimension A(*)
C
      Call HedPrt(IOut,0,String,Num)
      Call LTOutC(IOut,1,1,1,1,N,A,Key)
      Return
      End
*Deck LabSpn
      Subroutine LabSpn(NSpin,IOff,N,IRX,LabX)
      Implicit Integer(A-Z)
C
C     Generate numbers and labels for a set of functions for printing.
C
      Dimension IRX(*)
      Character*1 LabX(*), LabLS(0:1), LabAB(0:1)
      Save LabLS, LabAB
      Data LabLS/'l','s'/, LabAB/'a','b'/
C
      If(NSpin.lt.0) then
        Do 10 I = 1, N
          IP = IOff + I
          IB = (IP-1) / N
          IRX(I) = IP - N*IB
         LabX(I) = LabLS(IB)
   10   Continue
      else if(NSpin.eq.1) then
        Do 20 I = 1, N
          IRX(I) = IOff + I
          LabX(I) = ' '
   20   Continue
      else
        Do 30 I = 1, N
          IP = IOff + I
          IB = 1 - Mod(IP,2)
          IRX(I) = (IP+1) / 2
          LabX(I) = LabAB(IB)
   30    Continue
        endIf
      Return
      End
*Deck PrtThr
      Subroutine PrtThr(Key,Thresh)
      Implicit Real*8(A-H,O-Z)
C
C     Decode the Key argument to print routines and return the threshold.
C
      Save Zero, Ten
      Data Zero/0.0d0/, Ten/10.0d0/
C
      If(Key.gt.0) then
        Thresh = Zero
      else
        Thresh = Ten**(Key-6)
        endIf
      Return
      End
*Deck DisMat
      Subroutine DisMat(NAtoms,IAn,C,IFlag,NColX,IOut,Error,ICrowd,
     $  Conver)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to print a lower triangular matrix of intermolecular
C     distances, given the coordinates.
C
C     NAtoms ... Number of atoms.
C     IAn    ... Integer vector containing the atomic numbers for
C                each of the centers.
C     C      ... Floating point vector containing the coordinates
C                of the NAtoms centers.
C     IFlag  ... Conversion flag.  This variable controls scaling
C                of the distances before printing.
C                1 ... Scale output distances by one.
C                2 ... Scale output distances by Conver.
C                3 ... Scale output distances by one/Conver
C     NColX  ... Number of columns per printed block.
C     IOut   ... Fortran output unit for printing.  If zero, prinitng is
C                suppressed.  If negative, only error messages are printed.
C     Error  ... Logical variable set to true if any pair of atoms is
C                less than 0.5 angstroms (TooClo).
C     ICrowd   ... 0,1 Flag an error on zero distances only
C                  2   Flag an error on distances < 0.5 A
C                  3   Do not flag errors.
C
C     Conver ... A unit conversion factor.
C
      Logical Error, ErrI
      Parameter (NEl=200, NCMax=11, MaxSav=20)
      Real*8 MDCutO
      Dimension C(3,*),S(NCMax),IEl(-2:NEl),IAn(NAtoms),ISav(MaxSav),
     $  JSav(MaxSav)
      Save TooClo, Thresh, Zero, One
      Data TooClo/0.5d0/,Thresh/1.0d-6/,Zero/0.0d0/,One/1.0d0/
 2001 Format(10X,11(5X,I6))
 2003 Format(I6,2X,A2,11F11.6)
 2004 Format(' Small interatomic distances encountered: ',10I6)
 2005 Format(' Small interatomic distances encountered:',
     $  (/,1X,2I5,1PD9.2))
C
      Error = .False.
      If(ICrowd.le.1) then
        Cut = MDCutO(0)
      else if(ICrowd.eq.2) then
        Cut = TooClo
      else
        If(IOut.le.0) Return
        Cut = -Float(1)
        endIf
      Scale = One
      If(IFlag.eq.2) Scale = Conver
      If(IFlag.eq.3) Scale = One / Conver
      If(IOut.le.0) then
        Cut = Cut/Scale
        Do 30 I = 2, NAtoms
          Error = .False.
          If(IAn(I).gt.0) then
            Do 10 J = 1, (I-1)
              R = Sqrt( (C(1,I)-C(1,J))**2 +
     $                  (C(2,I)-C(2,J))**2 +
     $                  (C(3,I)-C(3,J))**2 )
              Error = Error.or.(R.lt.Cut.and.IAn(J).gt.0)
   10         Continue
              If(Error) then
              Do 20 J = 1, (I-1)
                R = Sqrt( (C(1,I)-C(1,J))**2 +
     $                    (C(2,I)-C(2,J))**2 +
     $                    (C(3,I)-C(3,J))**2 )
                If(R.lt.Cut) then
                  If(IOut.ne.0) Write(-IOut,2005) I, J, R*Scale
                  STOP 
                  endIf
   20           Continue
              endIf
            endIf
   30     Continue
      else
        Call FillEl(-2,NEl,IEl)
        NCol = Min(NColX,NCMax)
        IStart = 1
        Kount = 0
        Do 60 IStart = 1, NAtoms, NCol
          M = 0
          IEnd = Min(IStart+NCol-1,NAtoms)
          Write(IOut,2001) (IR,IR=IStart,IEnd)
          Do 50 I = IStart, NAtoms
            M = M + 1
            IRange = Min(M,NCol)
            Do 40 IR = 1, IRange
              J = IStart + IR - 1
              Temp = (C(1,I)-C(1,J))**2
     $          + (C(2,I)-C(2,J))**2 +(C(3,I)-C(3,J))**2
              S(IR) = Zero
              If(Temp.gt.Thresh) S(IR) = Scale*Sqrt(Temp)
              If(I.ne.J) then
                ErrI = S(IR).lt.Cut.and.IAn(I).gt.0.and.IAn(J).gt.0
                Error = Error.or.ErrI
                If(ErrI.and.Kount.lt.MaxSav) then
                  Kount = Kount + 1
                  ISav(Kount) = I
                  JSav(Kount) = J
                  endIf
                endIf
   40         Continue
            Write(IOut,2003) I,IEl(IAn(I)),(S(IR),IR=1,IRange)
   50       Continue
   60     Continue
        If(Error) Write(IOut,2004) (ISav(I),JSav(I),I=1,Kount)
        endIf
      Return
      End
*Deck ZPrint
      Subroutine ZPrint(IOut,NZ,IANZ,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,
     $  ToAng,ToDeg)
      Implicit Real*8(A-H,O-Z)
C
C     Z-MATRIX PRINTING ROUTINE.
C     CONVERTS FROM INTERNAL (BOHR/RADIAN) UNITS TO EXTERNAL
C     (ANGSTROM/DEGREE) UNITS LOCALLY FOR PRINTING USING ToAng and ToDeg.
C
      Logical Coord
      Parameter (MaxEl=204)
      Dimension IANZ(*), LBL(*), LAlpha(*), LBeta(*)
      Dimension BL(*), ALPHA(*), BETA(*), IZ(4,*), IEl(-2:204)
      Save Zero, One, F45
      Data Zero/0.0d0/, One/1.0d0/, F45/45.0D0/
 1000 Format(1X,99('-'))
 1010 Format(1X,27X,'Z-MATRIX (ANGSTROMS AND DEGREES)')
 1020 Format('   CD    Cent   Atom    N1       Length/X',
     $  '        N2       Alpha/Y        N3        Beta/Z          J')
 2110 Format(1X,I6,1X,I6,2X,A2)
 2120 Format(1X,I6,9X,A2)
 2210 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,')')
 2220 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,')')
 2310 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,')')
 2320 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,')')
 2410 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,') ',I6,1X,F8.3,'(',I6,') ',
     $          I6)
 2420 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,') ',I6,1X,F8.3,'(',I6,') ',
     $          I6)
 2430 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,7X,F11.6,7X,F11.6)
 2440 Format(1X,I6,1X,6X,2X,A2,2X,I6,F11.6,7X,F11.6,7X,F11.6)
C
C     PRINT THE HEADING.
C
C     TODEG = F45 / ATan(ONE)
      Call FillEl(-2,MaxEl,IEl)
      WRITE(IOUT,1000)
      WRITE(IOUT,1010)
      WRITE(IOUT,1020)
      WRITE(IOUT,1000)
C
C     First card.
C
      If(NZ.lt.1) goto 900
      ICard = 1
      Idx   = IAnZ(1)
      Coord = IZ(1,1).lt.0.or.Bl(1).ne.Zero.or.Alpha(1).ne.Zero.or.
     $  Beta(1).ne.Zero
      PBl = Bl(1) * ToAng 
      PA = Alpha(1) * ToAng
      PB = Beta(1) * ToAng
      ICent = 0
      If(IAnZ(1).ge.0) ICent = 1
      If(IAnZ(1).ge.0.and..not.Coord)
     $  Write(IOut,2110) ICard,ICent,IEl(Idx)
      If(IAnZ(1).lt.0.and..not.Coord) Write(IOut,2120) ICard, IEl(Idx)
      If(IAnZ(1).ge.0.and.Coord)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,1), PBl, PA, PB
      If(IAnZ(1).lt.0.and.Coord)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,1), PBl, PA, PB
      If(NZ.eq.1) goto 900
C
C     Second card.
C
      NP1   = LBL(2)
      ICard = 2
      Idx   = IAnZ(2)
      PBl   = BL(2) * ToAng 
      PA = Alpha(2) * ToAng
      PB = Beta(2) * ToAng
      If(IAnZ(2).ge.0) ICent = ICent + 1
      If(IAnZ(2).ge.0.and.IZ(1,2).gt.0)
     $  Write(IOut,2210) ICard,ICent,IEl(Idx),IZ(1,2),PBl,NP1
      If(IAnZ(2).lt.0.and.IZ(1,2).gt.0)
     $  Write(IOut,2220) ICard,IEl(Idx),IZ(1,2),PBl,NP1
      If(IAnZ(2).ge.0.and.IZ(1,2).le.0)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,2), PBl, PA, PB
      If(IAnZ(2).lt.0.and.IZ(1,2).le.0)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,2), PBl, PA, PB
      If(NZ.eq.2) goto 900
C
C     Third card.
C
      NP1   = LBL(3)
      NP2   = LAlpha(3)
      ICard = 3
      Idx   = IAnZ(3)
      PBl   = Bl(3) * ToAng
      PA    = Alpha(3) * ToDeg
      If(IZ(1,3).le.0) PA = Alpha(3) * ToAng 
      PB = Beta(3) * Conver
      If(IAnZ(3).ge.0) ICent = ICent + 1
      If(IAnZ(3).ge.0.and.IZ(1,3).gt.0) Write(IOut,2310) ICard,ICent,
     $  IEl(Idx),IZ(1,3),PBl,NP1,IZ(2,3),PA,NP2
      If(IAnZ(3).lt.0.and.IZ(3,1).gt.0) Write(IOut,2320) ICard,
     $  IEl(Idx),IZ(1,3),PBl,NP1,IZ(2,3),PA,NP2
      If(IAnZ(3).ge.0.and.IZ(1,3).le.0)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,3), PBl, PA, PB
      If(IAnZ(3).lt.0.and.IZ(1,3).le.0)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,3), PBl, PA, PB
      If(NZ.eq.3) goto 900
C
C     Cards 4 through NZ.
C
      Do 500 ICard = 4, NZ
        NP1   = LBL(ICard)
        NP2   = LAlpha(ICard)  
        NP3   = LBeta(ICard)
        Idx   = IAnZ(ICard)
        PBl   = Bl(ICard) * ToAng
        PA    = Alpha(ICard) * ToDeg
        PB    = Beta(ICard)  * ToDeg
        If(IZ(1,ICard).le.0) PA = Alpha(ICard) * ToAng 
        If(IZ(1,ICard).le.0) PB = Beta(ICard) * ToAng
        If(IAnZ(ICard).ge.0) ICent = ICent + 1
        If(IAnZ(ICard).ge.0.and.IZ(1,ICard).gt.0) Write(IOut,2410)
     $    ICard,ICent,IEl(Idx),IZ(1,ICard),PBl,NP1,IZ(2,ICard),PA,
     $    NP2,IZ(3,ICard),PB,NP3,IZ(4,ICard)
        If(IAnZ(ICard).lt.0.and.IZ(1,ICard).gt.0) Write(IOut,2420)
     $    ICard,IEl(Idx),IZ(1,ICard),PBl,NP1,IZ(2,ICard),PA,
     $    NP2,IZ(3,ICard),PB,NP3,IZ(4,ICard)
        If(IAnZ(ICard).ge.0.and.IZ(1,ICard).le.0)
     $    Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,ICard),
     $    PBl, PA, PB
        If(IAnZ(ICard).lt.0.and.IZ(1,ICard).le.0)
     $    Write(IOut,2440) ICard, IEl(Idx), IZ(1,ICard), PBl, PA, PB
  500   Continue
  900 Write(IOut,1000)
      Return
      End
*Deck BndOrd
      Subroutine BndOrd(IOut,IPrint,MxBnd,MxAt,NAtoms,IAn,NBond,IBond,
     $  NPi,NLP,Del,Bonder,C)
      Implicit Real*8 (A-H,O-Z)
      Dimension IAn(*),NBond(*),IBond(MxBnd,*),Bonder(MxBnd,*)
      Dimension NPi(*),NLP(*)
      Dimension C(3,*)
C Local (Scratch in the calling program)
      Dimension Del(*)
      Thresh=1.0d-3
      Do 10 IAt=1,NAtoms
       IAnI=IAn(IAt)
       NBI=NBond(IAt)
       NPiI=NPI(IAt)
       NLPI=NLP(IAt)
       Del(IAt)=0.0d0
       Do 20 II=1,NBI
        Bonder(II,IAt)=1.0d0 
        JAt=IBond(II,IAt)
        IAnJ=IAn(JAt)
        NPJ=NPi(JAt)
        NLPJ=NLP(JAt)
        If(NPJ.gt.0) then
         Del(IAt)=Del(IAt)+1.0d0
        ElseIf(NPiI.gt.0.and.NLPJ.ne.0) then
         Del(IAt)=Del(IAt)+4.0d-01 
        EndIf
   20  Continue
C      write(IOut,'(''Atom,NBond,NPi,NLP,Del'',4I3,F8.3)') IAt,NBI,NPiI,
C    $  NLPI,Del(IAt)
   10 Continue
      Do 30 IAt=1,NAtoms
       IAnI=IAn(IAt)
       NBI=NBond(IAt)
       NPiI=NPi(IAt)
       NLPI=NLP(IAt)
       DelI=Del(IAt)
       BeffI=0.0d0
       If(NPiI.ne.0) then
        BeffI=Float(NPiI)/DelI
       ElseIf(NLPI.ne.0.and.DelI.gt.thresh) then
        BeffI=5.0D-02
       EndIf
       do 40 II=1,NBI
        JAt=IBond(II,IAt)
        IAnJ=IAn(JAt)
        NBJ=NBond(JAt)
        NPIJ=NPi(JAt)
        NLPJ=NLP(JAt)
        DelJ=Del(JAt)
        BeffJ=0.0d0
        If(NPiJ.ne.0) then
         BeffJ=Float(NPiJ)/DelJ
        ElseIf(NLPJ.ne.0.and.NPiI.ne.0) then
         BeffJ=5.0D-02
        EndIf
        If(BEffI.lt.thresh.or.BEffJ.lt.thresh) goto 40
        Bonder(II,IAt)=(BeffI+BeffJ)/2.0d0+1.0d0
C       Write(IOut,'(2I5,3F10.3)') IAt,JAt,BeffI,BeffJ,Bonder(II,IAt)
   40  continue
   30 continue
        Write(IOut,'(/,5X,''Bond'',4X,''Pauling Bond Order'',4X,
     $   ''Topol. Bond Order'')')
        do 50 IAt=1,NAtoms
         NBI=NBond(IAt)
         do 60 II=1,NBI
          JAt=IBond(II,IAt)  
          Value=Distan(C,IAt,JAt,0)
          Val0=RCovCT(IAn(IAt),IAn(JAt))
          PaulIJ=exp((Val0-Value)/3.0d-01)
          If(JAt.gt.IAt) Write(IOut,'(2I5,F15.3,F20.3)')IAt,JAt,
     $      PaulIJ,Bonder(II,IAt)
   60  continue
   50 continue     
      Return
      End
*Deck FChkRd
      Subroutine FChkRd(In,IOut,Thresh,NAtoms,IAn,PhyCon,C)
      Implicit Real*8 (A-H,O-Z)
      Character String*80
      Dimension C(*)
      Dimension IAn(*)
      Dimension PhyCon(*)
      ToAng=PhyCon(1)
      Read(In,*) String
      Read(In,*) String
      Read(In,'(A49,I12)') String,NAtoms  
      If(String(1:15).ne.'Number of atoms') then
       write(IOut,'('' Error in Input'')') 
       Stop
      EndIf 
   20 Read(In,'(A80)') String
      If(String(1:14).ne.'Atomic numbers') goto 20
      Read(In,'(6I12)') (IAn(IAt),IAt=1,NAtoms)
   30 Read(In,'(A80)') String
      If(String(1:7).ne.'Current') goto 30
      NAt3=3*NAtoms
      NRow=NAt3/5
      Do 40 IXYZ=1,NRow 
       Ini=5*(IXYZ-1)+1
       IEnd=Ini+4
       Read(In,'(5E16.8)') (C(JJ),JJ=Ini,IEnd)
   40 Continue
      Ini=NRow*5+1
      IEnd=NAt3-Ini
      If(IEnd.ne.0) Read(In,'(5E16.8)') (C(JJ),JJ=Ini,Ini+IEnd) 
      Do 50 I=1,NAt3
       C(I)=C(I)*ToAng
       If(Abs(C(I)).lt.Thresh) C(I)=0.0d0 
   50 Continue   
      Return
      End
*Deck ZToC
      Subroutine ZToC(MaxNZ,NZ,IAnZ,IZ,Bl,Alph,Bet,TTest,NAtoms,IAn,C,
     $  CZ,A,B,D,Alpha,Beta,IOut,Error,TstAng)
      Implicit Real*8(A-H,O-Z)
C1ZToC
C
C SYNOPSIS:
C   ZToC(MaxNZ,NZ,IAnZ,IZ,Bl,Alph,Bet,TTest,NAtoms,IAn,C,CZ,A,B,D,
C   Alpha,Beta,IOut,Error,TstAng)
C
C DESCRIPTION:
C   This subroutine computes the cartesian coordinates, given the
C   Z-matrix. This routine returns coordinates both with and without
C   the dummy atoms.	
C
C VARIABLES:
C   MaxNZ  ... Maximum number of lines in Z-matrix.
C   NZ     ... Number of lines in the Z-matrix.
C   IAnZ   ... The atomic numbers of the Z-matrix centers.
C   IZ     ... The integer components of the Z-matrix.
C   Bl     ... The bond-lengths from the Z-matrix.
C   Alph   ... The bond-angles from the Z-matrix.
C   Bet    ... The dihedral angles from the Z-matrix.
C   TTest  ... Logical flag to enable testing for tetrahedral angles.
C	       This feature is useful in obtaining exact tetrahedral
C	       angles.	If any are found and this flag is set, the exact
C	       value is used.
C   NAtoms ... Number of atoms (dummies removed), returned.
C   IAn    ... Atomic numbers (dummies removed), returned.
C   C      ... Cartesian Coordinates (dummies removed), returned.
C   CZ     ... Cartesian Coordinates with dummies, returned.
C   A      ... Scratch vector of length NZ.
C   B      ... Scratch vector of length NZ.
C   D      ... Scratch vector of length NZ.
C   Alpha  ... Scratch vector of length NZ.
C   Beta   ... Scratch vector of length NZ.
C   IOut   ... Output unit for error print.  If IOut=0, no print is done.
C   Error  ... A logical variable set to true if ZToC is unable to complete
C              its task.
C   TstAng ... Whether to require angles to be .ge.0.and.lt.180.
C
C2HISTORY
C   Updated error messages - J. L. Sonnenberg 7/2009
C?
      Logical TTest,Error,TstAng,TetRnd,NonLin
      Dimension IAnZ(*),IZ(4,*),Bl(*),Alpha(*),Beta(*),IAN(*),
     $  C(3,*),CZ(3,*),A(*),B(*),D(*),Alph(*),Bet(*),U1(3),U2(3),U3(3),
     $  U4(3),VJ(3),VP(3),V3(3)
      Save Zero,One,Two,CutOff,Four,RTo0
      Data Zero,One,Two,CutOff/0.D0,1.D0,2.D0,1.D-08/
      Data Four,RTo0/4.D0,1.D-12/
 1000 Format(1X,I7,' Z-matrix lines is greater than the maximum of',I7,
     $  ' in subroutine ZToC.')
 1010 Format(' Error on Z-matrix line number',I7,': ',
     $  ' invalid Beta angle type',I5,'.')
 1020 Format(' Error on Z-matrix line number',I7,': ',
     $  ' reference made to an undefined center.')
 1030 Format(' Error on Z-matrix line number',I7,': ',
     $  ' multiple references to a center on the same line.')
 1041 Format(' Error on Z-matrix line number',I7,': ',
     $      ' Atom',I7,' too close to atom',I7)
 1042 Format(' Error on Z-matrix line number',I7,': ',
     $  ' Normalization problem. Are the reference atoms co-linear?')
 1043 Format(' Error on Z-matrix line number',I7,': ',
     $       ' Atoms',3I7,' are colinear')
 1044 Format(' Error on Z-matrix line number',I7,': ',
     $       ' No solution for these values of Alpha and Beta')
 1050 Format(' Error on Z-matrix line number',I7,': ',
     $  ' angle Alpha is outside the valid range of 0 to 180.')
 1060 Format(' Bond length on Z-matrix line number',I7,
     $       ' is not positive.')
 1070 Format(' Error on Z-matrix line number',I7,': ',
     $  ' angle Beta is outside the valid range of 0 to 180.')
 1080 Format(1X,I5,' tetrahedral angles replaced.')
C
C     Check for potential overflow.
C
      Error = NZ.gt.MaxNZ
      If(Error) then
        If(IOut.ne.0) Write(IOut,1000) NZ, MaxNZ
        Return
        endIf
C
C     Check for nonsense in the connectivity.
C
      Do 13 I = 2, NZ
        If(IZ(1,I).le.0.and.IZ(2,I).eq.0.and.IZ(3,I).eq.0.and.
     $     IZ(4,I).eq.0) goto 13
        If(IZ(1,I).ge.I.or.IZ(2,I).ge.I.or.IZ(3,I).ge.I.or.
     $     IZ(1,I).le.0) then
          Error = .True.
          If(IOut.ne.0) Write(IOut,1020) I
          endIf
        If(I.eq.3) then
          If(IZ(2,I).le.0) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1020) I
            endIf
          If(IZ(1,I).eq.IZ(2,I)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1030) I
            endIf
        else if(I.ge.4) then
          If(IAbs(IZ(4,I)).gt.2) then
            Error = .true.
            If(IOut.ne.0) Write(IOut,1010) I, IZ(4,I)
            endIf
          If(IZ(3,I).le.0) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1020) I
            endIf
          If((IZ(1,I).eq.IZ(2,I)).or.(IZ(1,I).eq.IZ(3,I)).or.
     $       (IZ(2,I).eq.IZ(3,I))) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1030) I
            endIf
          endIf
   13   Continue
      If(Error) Return
      Pi = Four*ATan(One)
      Call AClear(3*NZ,CZ)
C
C     Move angles to local arrays and optionally test for tetrahedral
C     angles.  If TstAng, test Alpha for out of range 0 to 180 degrees.
C     Test for negative bond lengths.
C
      NumTet = 0
      Do 20 I = 1, NZ
        Alpha(I) = Alph(I)
        Beta(I)  = Bet(I)
        If(IZ(1,I).gt.0) then
          If(I.gt.1.and.Bl(I).le.Zero) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1060) I
            endIf
          If(I.gt.2.and.TstAng.and.
     $      (Alpha(I).le.Zero.or.Alpha(I).ge.Pi)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1050) I
            endIf
          If(TTest) then
            If(TetRnd(Alpha(I))) then
              Alph(I) = Alpha(I)
              NumTet = NumTet + 1
              endIf
            If(TetRnd(Beta(I))) then
              Bet(I)  = Beta(I)
              NumTet = NumTet + 1
              endIf
            endIf
          If(I.gt.3.and.IZ(4,I).ne.0.and.
     $      (Beta(I).le.Zero.or.Beta(I).ge.Pi)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1070) I
            endIf
          endIf
   20   Continue
      If(NumTet.ne.0.and.IOut.ne.0) Write(IOut,1080) NumTet
      If(Error) Return 
C
C     Atom 1.
C
      If(IZ(1,1).le.0) then
        CZ(1,1) = Bl(1)
        CZ(2,1) = Alpha(1)
        CZ(3,1) = Beta(1)
        endIf
C
C     Atom 2.
C
      If(NZ.ge.2) then
        If(IZ(1,2).le.0) then
          CZ(1,2) = Bl(2)
          CZ(2,2) = Alpha(2)
          CZ(3,2) = Beta(2)
        else
          CZ(1,2) = CZ(1,1)
          CZ(2,2) = CZ(2,1)
          CZ(3,2) = CZ(3,1) + Bl(2)
          endIf
        Call ASUnit(3,CZ(1,2),CZ(1,1),U2)
        endIf
C
C     Atom 3  Put in the plane defined by atoms 1 and 2 and
C     the X-axis.  If this is impossible, use the Z-axis.
C     Be careful of placing atoms until at least one is off-axis.
C
      Do 80 I = 3, NZ
        If(IZ(1,I).le.0) then
          CZ(1,I) = Bl(I)
          CZ(2,I) = Alpha(I)
          CZ(3,I) = Beta(I)
          Call ASUnit(3,CZ(1,I),CZ(1,1),U1)
          NonLin = (One-Abs(SProd(3,U1,U2))).gt.CutOff
        else
          If(IZ(1,I).eq.1) then
            Call ASUnit(3,CZ(1,IZ(1,I)),CZ(1,2),U1)
          else
            Call ASUnit(3,CZ(1,IZ(1,I)),CZ(1,1),U1)
            endIf
          Call Place(CutOff,U1,CZ(1,IZ(1,I)),Bl(I),Alpha(I),CZ(1,I),
     $      NonLin)
          endIf
        If(NonLin) goto 90
   80   Continue
      I = NZ
C
C     General Nth atom for non-linear molecule.
C
   90 K = I + 1
      Do 250 J = K, NZ
        DCAJ = Cos(Alpha(J))
        DSAJ = Sin(Alpha(J))
        DCBJ = Cos(Beta(J))
        DSBJ = Sin(Beta(J))
C
C       Cartesian coordinates
        If(IZ(1,J).le.0) then
          CZ(1,J) = Bl(J)
          CZ(2,J) = Alpha(J)
          CZ(3,J) = Beta(J)
C
C       Dihedral angle
        else if(IZ(4,J).eq.0) then
          Call Vec(CutOff,Error,U1,CZ,IZ(2,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(1,J),IZ(2,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(2,J)
            Return
            endIf
          Call VProd(VP,U1,U2)
          Arg = One - SProd(3,U1,U2)**2
          Error = Arg.lt.Zero
          If(Error) then
            If(IOut.ne.0) Write(IOut,1042) J
            Return
            endIf
          R = Sqrt(Arg)
          Error = R.lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
            Return
            endIf
          Call AScale(3,(One/R),VP,U3)
          Call VProd(U4,U3,U2)
          Do 130 I = 1, 3
            VJ(I) = (-U2(I)*DCAJ+U4(I)*DSAJ*DCBJ+U3(I)*DSAJ*DSBJ)*Bl(J)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  130     Continue
C
C       Second bond angle specified.
        else if(IAbs(IZ(4,J)).eq.1) then
          Call Vec(CutOff,Error,U1,CZ,IZ(1,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(2,J),IZ(1,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(1,J)
            Return
            endIf
          Zeta = -SProd(3,U1,U2)
          Denom = One - Zeta** 2
          Error = Denom.lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
            Return
            endIf
          A(J) = (-DCBJ+Zeta*DCAJ)/Denom
          B(J) = (DCAJ-Zeta*DCBJ)/Denom
          R = Zero
          Gamma = Pi / Two
          If(Abs(Zeta).ge.CutOff) then
            If(Zeta.lt.Zero) R = Pi
            Error = Denom.lt.Zero
            If(Error) then
              If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
              Return
              endIf
            Gamma = ATan(Sqrt(Denom)/Zeta) + R
            endIf
          D(J) = Zero
          If(Abs(Gamma+Alpha(J)+Beta(J)-Two*PI).ge.CutOff) then
            Arg = (One + A(J)*DCBJ - B(J)*DCAJ)  /  Denom
            Error = Arg.lt.Zero
            If(Error) then
              If(IOut.ne.0) Write(IOut,1044) J
              Return
              endIf
            D(J) = Float(IZ(4,J)) * Sqrt(Arg)
            endIf
          Call VProd(V3,U1,U2)
          Do 200 I = 1, 3
            U3(I) = A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
            VJ(I) = Bl(J)*U3(I)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  200     Continue
C
C       Angle to plane specified.
        else if(IAbs(IZ(4,J)).eq.2) then
          Call Vec(CutOff,Error,U1,CZ,IZ(1,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(2,J),IZ(1,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(1,J)
            Return
            endIf
          Zeta = -SProd(3,U1,U2)
          Call VProd(V3,U1,U2)
          V3MAG = Sqrt(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
          Denom = One - Zeta**2
          Error = Abs(Denom).lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1042) J
            Return
            endIf
          A(J) = V3MAG*DCBJ / Denom
          Arg = (One-DCAJ*DCAJ-A(J)*DCBJ*V3MAG) / Denom
          Error = Arg.lt.Zero
          If(Error) then
            If(IOut.ne.0) Write(IOut,1044) J
            Return
            endIf
          B(J) = Sqrt(Arg)
          If(IZ(4,J).ne.2) B(J) = -B(J)
          D(J) = B(J)*Zeta+DCAJ
          Do 240 I = 1, 3
            U3(I) = B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
            VJ(I) = Bl(J)*U3(I)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  240     Continue 
          endIf
  250   Continue
C
C     Eliminate dummy atoms, which are characterized by negative atomic
C     numbers.  Ghost atoms (zero atomic number) are retained.  Tidy up
C     the coordinates by zeroing elements less than RTo0.
C
      NAtoms = 0
      Do 280 I = 1, NZ
        If(IAnZ(I).ne.-1) then
          NAtoms = NAtoms + 1
          IAn(NAtoms) = IAnZ(I)
          C(1,NAtoms) = CZ(1,I)
          C(2,NAtoms) = CZ(2,I)
          C(3,NAtoms) = CZ(3,I)
          endIf
  280   Continue
      I = IPopVc(3*NAtoms,C,RTo0,C)
      Return
      End
*Deck FndGrp
      SUBROUTINE FndGrp(LINE, GROUP, GROUP_NORM)

      CHARACTER*(*) LINE, GROUP, GROUP_NORM
      CHARACTER*32 TOK(10), TWORK
      INTEGER NTOK, I
      LOGICAL IS_POINT_GROUP

      GROUP  = ' '
      GROUP_NORM = ' '

      CALL SPLIT_TOKENS(LINE, TOK, NTOK)

      DO 200 I = 1, NTOK

C        copia di lavoro per il gruppo
         TWORK = TOK(I)
         CALL TOLOWER(TWORK)

         IF (IS_POINT_GROUP(TWORK)) THEN
C           salva il gruppo ORIGINALE (non alterato)
            GROUP = TOK(I)
            GROUP_NORM = TWORK
         END IF

  200 CONTINUE

      RETURN
      END
