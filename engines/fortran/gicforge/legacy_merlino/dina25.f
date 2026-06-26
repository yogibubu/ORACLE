*Deck DiNa25
C
C New version started 2th November 2024
C This deck contains the driver and the routines
C Aclear,Amove,IClear,LinUpC,PhyFil
C and the function CnvFct
C
      Implicit Real*8 (A-H,O-Z)
C MaxEl  = Number of different elements
C MxScr  = Dimension of Scratch Array
C MxAt   = Maximum Number of Atoms
C MaxNZ  = Maximum Number of ZMatrix Elements
C MxBnd  = Maximum coordination
C MxBox  = Maximum number of boxes for new connectivity
C MxKwd  = Maximum Number of KeyWords
C LenPhy = Number of Physical Constants
C MxGNIC = Maximum Number of GNICS
C MxTrm  = Maximum Number of Terms in GNICs
C MxAtP  = Maximum Number of Atoms in a Primitive
C MxCyc  = Maximum Number of Cycles
C MxAtCy = Maximum Number of Atoms in a Cycle
C MxPot  = Max Terms in Potential Fitting
      Parameter(MaxEl=200)
      Parameter(MxScr=100000,MxAt=1000,MxBnd=10,MxKwd=50,LenPhy=30)
      Parameter(MxGNIC=1000,MxTrm=45,MxAtP=4,MxCyc=20,MxAtCy=10)
      Parameter(MxFrg=100,MxAtFr=30,MxAtB=100,MaxNZ=1000,MxBox=1000)
      Parameter(MxPot=20)
      Character*80 FilNam,InFil,OutFil,GauKwd,Title,React,LinScr
      Character*20 StrInp
CENZO
      Character*16 Group
      Common/IO/In,IOut,IPunch
      Common/PhyCon/PhyCon(30)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
      Integer St2Int,El2IAn,Istart(100)
      Logical Kwd(MxKwd),Loose
      Logical Linear,ImpDih,DoGDV,PrtPic,SyGNIC,DoBPCS,Clean
      Logical DoEck,Do1Dih,DoNorm,DVibRot,InvDst,LConn,DoSySt,Aver
      Logical DoBMat,DoMW,Inv1,DoScan,DoRig,RIgB,RigA,RigL,RigD,RigO
      Logical DoneC,DoB1,DoGNIC,RdIsot
      Logical DoG16,DoLocSVD
      Logical UsedPrim
      Logical TstAng,TTest,Error,PrtVal
      Dimension IEl(0:MaxEl)
      Dimension IScr(MxScr)
      Dimension Scr(MxScr)
C Dimensions for Atom Properties (including Coordinates)
      Dimension IAn(MxAt),IFrag(MxAt)
      Dimension ISot(MxAt),MultN(MxAt)
      Dimension IAnZ(MaxNZ),IZ(4,MaxNZ),MapZAt(MaxNZ)
      Dimension LBl(MaxNZ),LAlpha(MaxNZ),LBeta(MaxNZ)
      Dimension BL(MaxNZ),Alpha(MaxNZ),Beta(MaxNZ)
      Dimension QMom(MxAt),GFac(MxAt)
      Dimension PMom(3),RotGHz(3),RTemp(3),DBVib(3),AtMass(MxAt)
      Dimension XYZCM(3),RotMat(3),Rotcm1(3),PMomB(3)
      Dimension C(3,MxAt),CZ(3,MaxNZ),TMom(6)
C Scratch for ZTOC
      Dimension A(MaxNZ),B(MaxNZ),D(MaxNZ),Alpha1(MaxNZ),Beta1(MaxNZ)
C Dimensions for vibrations
      Dimension Frq(3*MxAt)
C Dimension for Fragments
      Dimension IFrsAt(MxFrg),ILstAt(MxFrg),NAtFr(MxFrg)
      Dimension IAtFr(MxAtFr,MxFrg),LConn(MxFrg)
C Dimensions for Connectivity
      Dimension NBond(MxAt),IBond(MxBnd,MxAt),NH(MxAt),IArom(MxAt)
      Dimension BndOrd(MxBnd,MxAt),InCyc(MxAt)
C Dimensions for Stretchings
      Dimension IAtomB(MxAtP,MxTrm,MxGNIC),NTermB(MxGNIC),ITVB(MxGNIC)
      Dimension NEqAt(2,MxGNIC),IFixB(MxGNIC)
      Dimension CoefB(MxTrm,MxGNIC),ValTB(MxGNIC)
C Dimensions for Bendings
      Dimension IAtomA(MxAtP,MxTrm,MxGNIC),NTermA(MxGNIC),ITVA(MxGNIC)
      Dimension IFixA(MxGNIC)
      Dimension CoefA(MxTrm,MxGNIC),ValTA(MxGNIC)
C Dimensions for Linear Bendings
      Dimension IAtomL(MxAtP,MxTrm,MxGNIC),NTermL(MxGNIC),ITVLA(MxGNIC)
      Dimension IFixL(MxGnic)
      Dimension CoefL(MxTrm,MxGNIC),ValTL(MxGNIC)
C Dimensions for Torsions
      Dimension IAtomD(MxAtP,MxTrm,MxGNIC),NTermD(MxGNIC),ITVD(MxGNIC)
      Dimension IPerD(MxGNIC),IFixD(MxGNIC)
      Dimension CoefD(MxTrm,MxGNIC),ValTD(MxGNIC)
C Dimensions for Out-of-Plane Bendings
      Dimension IAtomO(MxAtP,MxTrm,MxGNIC),NTermO(MxGNIC),ITVO(MxGNIC)
      Dimension IFixO(MxGNIC)
      Dimension CoefO(MxTrm,MxGNIC),ValTO(MxGNIC)
C Dimensions for Syntons
      Dimension EAN(MxAt),EANZ(MaxNZ)
C Dimensions for Cycles
      Dimension NatC(MxCyc),ICAt(MxAtCy,MxCyc)
      Dimension IAtCyc(MxAt),IBr(2,MxAt)
C Dimensions for primitives
      Dimension IAtmBR(MxAtP,MxTrm,MxGNIC),IPrimB(MxTrm,MxGNIC)
      Dimension IAtmAR(MxAtP,MxTrm,MxGNIC),IPrimA(MxTrm,MxGNIC)
      Dimension IAtmLR(MxAtP,MxTrm,MxGNIC),IPrimL(MxTrm,MxGNIC)
      Dimension IAtmDR(MxAtP,MxTrm,MxGNIC),IPrimD(MxTrm,MxGNIC)
      Dimension IAtmOR(MxAtP,MxTrm,MxGNIC),IPrimO(MxTrm,MxGNIC)
C Dimension for B and DB matrices
      Dimension BMat(3*MxAtB,MxGNIC)
C Dimensions for BDPCS3
      Dimension R0IJ(MxGNIC),RBIJ(MxGNIC)
C ACN and ACN1 are coordination numbers (from D3 and D4 of Grimme)
C SPIJ is the sum of bond orders
      Dimension ACN(MxAt),ACN1(MxAt),SPIJ(MxAt),Teta0(4)
C Dimensions for Potential Fitting
      Dimension IPot(MxPot),Cof(MxPot)
C Dimension for Gaussian Keywords
      Dimension IDeriv(4)
C Long Formats
 1000 Format(' Atomic Number =',I3,'; Isotope =',I3,
     $  '; Nuc.Spin Multip =',I5,/,' Mass =',F12.5,
     $  '; Quadrup.Mom. =',F10.5,'; Nuc.Magn.Mom. =',F10.5)
C Input/Output Files
      In=1
      IOut=2
      IPunch=7
C Input and Output Files are wayin and wayout
      filnam='topo24'
      ITpFil=0
      OPEN(In,FILE='provin',STATUS='OLD')
      OPEN(IOut,FILE='provout',STATUS='UNKNOWN')
      Rewind(IOut)
C which version of physical constants
      IOpt=2010
C Set physical constants and treshold for linear angles in radiants
C (171.9 degrees)
      Call PhyFil(IOpt,LenPhy,PhyCon)
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      TreshL=3.0d0
C Set Atomic Symbols
      call FillEl(0,MaxEl,IEl)
C Read keywords
      Call LlinCl(LinScr)
      Read(In,'(A80)') LinScr
      Write(IOut,'(A80)') LinScr
      Call FndKwd(IOut,IPrint,ModPCS,IDeriv,LinScr,Kwd,MxKwd)
      DoG16=Kwd(2)
      DoGDV=Kwd(3)
      Aver=Kwd(5)
      DoEck=Kwd(6)
      ImpDih=Kwd(7)
      DoGNIC=Kwd(8)
      InvDst=Kwd(9)
      DoSySt=Kwd(10)
      Do1Dih=.not.Kwd(23)
      DoNorm=Kwd(12)
      DoBPCS=Kwd(13)
      SyGNIC=Kwd(14).or.Kwd(24)
      DoBMat=Kwd(19)
      DoScan=Kwd(20)
      DoRig=Kwd(21)
      DoLocSVD=Kwd(22)
      RdIsot=.False.
      Clean=Kwd(31)
      DVIBRot=.False.
      Loose=Kwd(33)
      If(DoG16) ImpDih=.true.
      RigB=.False.
      RigA=.False.
      RigL=.False.
      RigD=.False.
      RigO=.False.
      If(DoBPCS) RigB=.True.
      If(DoScan.and.DoRig) then
       RigB=.True.
       RigA=.True.
       RigL=.True.
       RigO=.True.
      EndIf
C Read Title preceeded and followed by a blank card (Gaussian style)
      Read(In,'(A80)')    LinScr
      Read(In,'(A80)')    Title
      Write(IOut,'(/,A80)') Title
      Read(In,'(A80)')    LinScr
C Read charge and multiplicity
      call LlinCl(StrInp)
      Read(In,'(A20)') Strinp
      call SubStr(StrInp,2,IStart,NValue)
      If(NValue.ne.2) then
       write(IOut,'(A20)') Strinp
       write(IOut,'(''Wrong Charge or Multiplicity'')')
        Stop
      EndIf
      ICharg=St2Int(StrInp(IStart(1):IStart(2)-1),-999)
      If(ICharg.eq.-999) then
       write(IOut,'(''Wrong Charge'')')
       Stop
      EndIf
      Multip=St2Int(StrInp(IStart(2):len(StrInp)),-999)
      If(Multip.eq.-999) then
       write(IOut,'(''Wrong Multiplicity'')')
       Stop
      EndIf
C Read Cartesian geometry and set connectivity.  Coordinate generation from
C any non-XYZ source is owned by the Python layer.
      call Coord(In,IOut,IPunch,IPrint,MxAt,MaxNZ,MxBnd,MxBox,PhyCon,
     $  KWd,Multip,NAtoms,NFrag,NBond,NH,NZ,IAn,Isot,IFrag,IBond,IANZ,
     $  IZ,MapZAt,Linear,C,CZ,EAN,EANZ,AtMass,TotWt,PMom,RotGHz,RTemp,
     $  MultN,QMom,GFac,NSpec,IScr,Scr,Group)
      NTRot=6*NFrag
      If(Linear) NTRot=5*NFrag
C Print Coordinates and Rotational Constants
      if(DoEck) then
       Write(IOut,'(/,10X,''Cartesian Coords.in Eckart Orientation'')')
      else
       Write(IOut,'(/,10X,''Cartesian Coords. in Input Orientation'')')
      endif
      If(NFrag.eq.1) then
       Write(IOut,'(4X,''Atom'',4X,''Pauling EAN  Topolog.EAN'',6X,
     $   ''X'',10X,''Y'',12X,''Z'')')
      Else
       Write(IOut,'(4X,''Atom'',4X,''Pauling EAN  Topolog.EAN'',6X,
     $  ''X'',10X,''Y'',12X,''Z'',5X,''Fragment'')')
      EndIf
      Do 5 IAt=1,NAtoms
       If(NFrag.eq.1) then
        Write(IOut,'(I5,2X,A2,2(4X,F8.5),3F12.5)') IAt,IEl(IAn(IAt)),
     $   EAn(IAt),EANZ(IAt),(C(IXYZ,IAt),IXYZ=1,3)
       Else
        Write(IOut,'(I5,2X,A2,2(4X,F8.5),3F12.5,2X,I5)') IAt,
     $   IEl(IAn(IAt)),EAn(IAt),EANZ(IAt),(C(IXYZ,IAt),IXYZ=1,3),
     $   IFrag(IAt)
       EndIf
   5  Continue
      write(IOut,'(/,I5,'' Atoms with Total Mass ='',F10.3)') NAtoms,
     $  TotWt
      write(IOut,'('' Equilibrium Rotational Constants (MHz)'',
     $  3(1X,F12.5))') (1.0d+3*RotGHz(i),i=1,3)
      if(DVibRot) write(IOut,'('' Ground State Rotational Constants'',
     $  '' (MHz)  '',3F12.5)') (1.0d+3*RotGHz(i)+DBVib(i),i=1,3)
      write(IOut,'('' '')')
C Print interatomic distances
      Call LlinCl(LinScr)
      LinScr(1:22)=' Interatomic Distances (Angstrom)'
      Call HedPrt(IOut,0,LinScr,Num)
      ScalI=1.0d0
      Call DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
      write(IOut,'('' '')')
C Make partition functions
      call Therm1(In,IOut,IPrint,PhyCon,Multip,NAtoms,Linear,
     $ TotWt,ZPE,RTemp,DBVib,Frq)
C Make primitive bond lengths and valence angles including linear ones
      NLenR=0
      NAngR=0
      NLAngR=0
      NTotR=0
      call MkBAL(IOut,IPrint,MxBnd,MxAtP,MxTrm,NAtoms,NLenR,NAngR,
     $  NLAngR,IAn,NBond,IBond,NTermB,NTermA,NTermL,IAtmBR,IAtmAR,
     $  IAtmLR,CoefB,CoefA,CoefL,C,TreshL)
      NTotR=NLenR+NAngR+NLAngR
      if(NLenR.ne.0.and.DoBPCS) then
      write(IOut,'(/,'' Bonded Atoms  rDSD Distance BDPCS3 distance'')')
       do 10 il=1,NLenR
        IAt=IAtmBr(1,1,Il)
        JAt=IAtmBr(2,1,Il)
        R0IJ(IL)=Distan(C,IAt,JAt,0)
        RBIJ(IL)=RBPCS(IAt,JAt,IAn,R0IJ(IL))
        write(IOut,'(2I5,7X,F8.5,7X,F8.5)') IAt,JAt,R0IJ(IL),RBIJ(IL)
   10  continue
       write(IOut,'('' '')')
      endif
C Make primitive dihedrals
      NDihR=0
      IPrDih=0
      call MkPrmD(IOut,IPrDih,MxBnd,MxTrm,MxTrm,MxAtP,MxAtP,MxAtCy,
     $  MxCyc,NAtoms,NBond,NLenR,NDihR,NTotR,NCyc,IBond,NTermD,IAtmBR,
     $  IAtmDR,NAtC,ICAt,IAtCyc,CoefD,C,TreshL)
      NTotR=NTotR+NDihR
C Make primitive out-of-plane bendings (use improper dihedrals for G16)
      NOuplR=0
      call MkPrmO(IOut,IPrOup,MxBnd,MxTrm,MxAtP,NAtoms,NBond,
     $  NOuplR,IBond,NTermO,IAtmOR,CoefO,C,ImpDih,IAtCyc)
      NTotR=NTotR+NOuPlR
      if(NTotR.eq.0) NTotR=NLenR+NAngR+NLAngR+NDihR+NOuPlR
C Print summary of internal primitives
      If(.not.DoGNIC) then
       Write(IOut,'(/,I5,'' Atoms and'',I5,'' Internal Coordinates'')')
     $  NAtoms,3*NAtoms-NTRot
       Write(IOut,
     $ '(/,14X, '' Stretch.  Bend.  L. Bend. Tors.  Out-Pl. Total'')')
        Write(IOut,'('' Redundant  '',6I8)') NLenR,NAngR,NLAngR,NDihR,
     $  NOuPlR,NTotR
      EndIf
C Compute B matrix (not for G16)
      If(DoBPCS) then
       If(DoG16) then
        write(IOut,'('' BDPCS3 not available for G16'')')
        Stop
       EndIf
C      Call DrvBG(IOut,IPrint,Linear,DoBPCS,DoneC,MxAtP,MxTrm,NAtoms,
C    $   NLenR,NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,IAtmLR,
C    $   IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,C,Atmass,R0IJ,
C    $   RBIJ,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr)
      Call PCSGeo(IOut,IPrint,Linear,DoneC,MxAtP,MxTrm,NAtoms,NLenR,
     $  NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,IAtmLR,
     $  IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,C,Atmass,R0IJ,
     $  RBIJ,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr)
       If(.not.DoneC) then
        write(IOut,'('' Failed Conversion to Cartesian Coordinates'')')
        Stop
       EndIf
       Write(IOut,'(/,'' BDPCS3 Cartesian Coordinates'')')
       Write(IOut,'(3X,''Atom'',2X,''At.Numb. Eff.At.Numb.'',5X,''X'',
     $    10X,''Y'',12X,''Z'')')
       Write(IOut,'(I5,I7,7X,F8.5,3F12.5)')(IAt,IAn(IAt),EAn(IAt),
     $   (C(IXYZ,IAt),IXYZ=1,3),IAt=1,NAtoms)
      Call MofI(IOut,Iprint,.false.,1,NAtoms,C,AtMass,XYZCM,TMom,PMom,
     $  RotMat)
C Compute rotational parameters
       FactA = CnvFct('FactA')
       Fac1  = PhyCon(1)**2
       AvPMom=1.0D0
       AvPMoB=1.0D0
       do 85 I=1,3
        If(PMom(i).le.0.0d0) goto 85
        AvPMom=AvPMom*PMom(i)
        PMomB(i)=PMom(i)/Fac1
        AvPMoB=AvPMoB*PMomB(i)
        Rotcm1(i)=FactA/PMom(i)
   85  Continue
       AvPMom=AvPMom**(1.0D0/3.0D0)
       AvPMoB=AvPMoB**(1.0d0/3.0d0)
       Call RotCon(IOut,IPrint,Linear,NAtoms,PhyCon,PMomB,RotGHz,RTemp)
       Write(IOut,'(/,'' Equilibrium Rotational Constants (MHz)'',
     $   3(1X,F12.5))') (1.0d+3*RotGHz(i),i=1,3)
       if (DVibRot) then
        write (IOut,'('' Ground State Rotational Constants (MHz)'',
     $    3(1X,F12.5))') (1.0d+3*RotGHz(ii)+DBVib(ii),ii=1,3)
       EndIf
       Write(IOut,'('' '')')
C Print interatomic distances
       Call LlinCl(LinScr)
       LinScr(1:22)=' Interatomic Distances (Angstrom)'
       Call HedPrt(IOut,0,LinScr,Num)
       ScalI=1.0d0
       Call DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
       write(IOut,'('' '')')
      EndIf
      If(.not.DoGNIC) then
       goto 999
      Else
       write(IOut,'(/,'' Definition of Generalized Natural Internal'',
     $   '' Coordinates'')')
      EndIf
C Make cycles
      NExpCy=NLenR-NAtoms+NFrag
      IPrCyc=0
      call MkCyc(IOut,IPrCyc,MxBnd,MxAtP,MxTrm,
     $  MxAtCy,MxCyc,NBond,IBond,NLenR,NAngR,NDihR,NExpCy,
     $  IAtmBR,IAtmAR,IAtmDR,NCyc,NatC,ICAt,IAtCyc,IBr,IAn,EAN)
      If(NCyc.gt.0.or.NExpCy.gt.0) write(IOut,'(I3,
     $  '' Rings Found Over'',I3,'' Expected'',/)') NCyc,NExpCy
C Make bond GNICs.  Stretchings are final primitive bond coordinates:
C no local SALC or redundancy reduction is applied at this stage.
      NLen=0
      IPrBnd=0
      Call MkGNCB(IOut,IPrBnd,.False.,InvDst,MxBnd,MxTrm,MxAtP,NAtoms,
     $  IAn,NBond,NLen,IBond,NTermB,IAtomB,ITVB,IAtCyc,CoefB,C)
C     If(NCyc.gt.0.and.DoSySt) then
C      IPrtCB=0
C Cycles NYI for stretchings
C      Call CySalc(IOut,IprtCB,MxAtP,MxAtCy,MxTrm,NAtoms,NLen,
C    $   NLenR,NatC,ICAt,IAtCyc,IAtmBr,IAtomB,NTermB,ITVB,CoefB,EAn,C)
C     EndIf
C Make angle GNICs
C      NAng=0
      IPrAng=0
      call MkGNCA(IOut,IPrAng,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NCyc,MxAtCy,NAtC,ICAt,NBond,NAng,IBond,NTermA,IAtomA,IAn,
     $  IAtCyc,ITVA,CoefA,C,EAN,TreshL,DoLocSVD)
C Make ring coordinates for valence angles
      If(NCyc.gt.0) then
       Write(IOut,'('' Endocyclic Valence Angles'')')
       IPrtCA=0
       NAng00=NAng
       do 20 icyc=1,NCyc
        NAng0=NAng
        call CyGNA(IOut,IPrtCA,MxAtCy,MxAtP,MxTrm,NBond,NAtC,NAng,ICyc,
     $    ICAt,NTermA,IAtomA,ITVA,CoefA,DoLocSVD,NAtoms,C)
        If(NAtC(ICyc).eq.3) then
         write(IOut,'(1X,I2,A,3I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        ElseIf(NAtC(ICyc).eq.4) then
         write(IOut,'(1X,I2,A,4I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        ElseIf(NAtC(ICyc).eq.5) then
         write(IOut,'(1X,I2,A,5I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        ElseIf(NAtC(ICyc).eq.6) then
         write(IOut,'(1X,I2,A,6I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        ElseIf(NAtC(ICyc).eq.7) then
         write(IOut,'(1X,I2,A,7I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        ElseIf(NAtC(ICyc).eq.8) then
         write(IOut,'(1X,I2,A,8I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NAng-NAng0,' Valence angles'
        EndIf
   20  continue
      EndIf
C Make Linear Angle GNICs: make 2 angles for linear molecules
      NLAng=0
      IPrLAn=0
      call MkGNLA(IOut,IPrLAn,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NBond,NLAng,Linear,IBond,NTermL,IAtomL,IAn,CoefL,C,TreshL,
     $  DoLocSVD)
C Make dihedral GNICs
      NDih=0
      call MkGNCD(IOut,IPrint,MxBnd,MxTrm,MxTrm,MxAtP,MxAtP,MxAtCy,
     $  Do1Dih,NAtoms,IAn,NBond,NLenR,NDih,NTot,NCyc,IBond,NTermD,
     $  IAtmBR,IAtomD,IBr,NAtC,ICAt,IAtCyc,ITVD,IPerD,NEqAt,CoefD,C,
     $  EAN,TreshL,DoNorm)
      NButD=0
      Do 25 IDih=1,NDih
       If(ITVD(IDih).eq.2) NButD=NButD+1
   25 Continue
      NExoD=NDih-NButD
C Make ring coordinates for dihedra angles
      If(NCyc.gt.0) then
       Write(IOut,'(/,'' Exocyclic Dihedral Angles:'',I5)') NExoD
       Write(IOut,'('' Endocyclic Dihedral Angles'')')
       NDihCh=NDih
       IPrtCD=0
       do 30 ICyc=1,NCyc
        NDih0=NDih
        call CyGND(IOut,IPrtCD,MxAtCy,MxAtP,MxTrm,NAtC,NDih,ICyc,ICAt,
     $    NTermD,IAtomD,ITVD,CoefD,DoLocSVD,NAtoms,C)
        If(NAtC(ICyc).eq.3) then
         write(IOut,'(1X,I2,A,3I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        ElseIf(NAtC(ICyc).eq.4) then
         write(IOut,'(1X,I2,A,4I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        ElseIf(NAtC(ICyc).eq.5) then
         write(IOut,'(1X,I2,A,5I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        ElseIf(NAtC(ICyc).eq.6) then
         write(IOut,'(1X,I2,A,6I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        ElseIf(NAtC(ICyc).eq.7) then
         write(IOut,'(1X,I2,A,7I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        ElseIf(NAtC(ICyc).eq.8) then
         write(IOut,'(1X,I2,A,8I4,A,I2,A)')
     $     NAtC(ICyc),'-Membered Ring (',
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),'): ',
     $     NDih-NDih0,' Dihedral angles'
        EndIf
   30  continue
       If(NButD.gt.0) then
        Do 35 IBrid=1,NBrL
         write(IOut,'('' Butterfly GNIC Around Bond'',
     $    I5,''  -'',I5,'' Joining Rings'',2I3)')
     $    IBrL(1,IBrid),IBrL(2,IBrid),IBrL(3,IBrid),IBrL(4,IBrid)
   35   Continue
       EndIf
      endif
C Make Out-of-Plane GNICs
      NOupl=0
      call MkGNCO(IOut,IPrint,.True.,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NCyc,NBond,NOuPl,IBond,NTermO,IAtomO,IAn,IAtCyc,CoefO,C,ImpDih)
C
      Write(IOut,'(/,'' Out-of-plane angles:'',I5,/)') NOuPl
      NTot=NLen+NAng+NLAng+NDih+NOUPl
      NTarget=3*NAtoms-NTRot
      UsedPrim=.False.
      If(NTot.lt.NTarget) then
       Write(IOut,'(/,'' GNIC candidate count below vibrational rank;'',
     $ '' expanding to primitive candidates before pruning.'')')
       Write(IOut,'(''   Current GNIC candidates='',I5,
     $ '' target='',I5,'' primitive candidates='',I5)')
     $ NTot,NTarget,NTotR
       Call UsePrimitiveGICs(IOut,MxAtP,MxTrm,NLenR,NAngR,NLAngR,
     $ NDihR,NOuplR,NLen,NAng,NLAng,NDih,NOupl,IAtmBR,IAtmAR,
     $ IAtmLR,IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,
     $ IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,ITVB,ITVA,ITVLA,ITVD,
     $ ITVO,IFixB,IFixA,IFixL,IFixD,IFixO,CoefB,CoefA,CoefL,
     $ CoefD,CoefO)
       NTot=NLen+NAng+NLAng+NDih+NOupl
       If(NTot.lt.NTarget) then
        Write(IOut,'('' ERROR: primitive GIC candidates='',I5,
     $ '' below target vibrational rank='',I5)') NTot,NTarget
       Stop
       EndIf
       UsedPrim=.True.
      EndIf
CENZO Print Information on Torsions
      Indd=0
      do IPuf=1,NLenR
       IScr(indd+1)=IAtmBr(1,1,IPuf)
       Iscr(Indd+2)=IAtmBr(2,1,IPuf)
       Indd=Indd+2
      enddo
      Call DrvTrs(IOut,IPrint,MxBnd,NAtoms,NLenR,IAn,IAtCyc,NBond,IBond,
     $IScr,EAn,C)
CENZO
      NRed=NTot-NTarget
      IPrOrd=IPrint
      If(NRed.ne.0) IPrOrd=-1
      If(NRed.eq.0) then
       write(IOut,'('' All local redundancies have been '',
     $  ''Eliminated'',/)')
      else
       write(IOut,'('' Pre-pruning residual redundancies:'',I3,
     $  '' (handled by type-preserving rank check)'')') NRed
       write(IOut,'('' Coordinate definitions are printed after '',
     $  ''pruning only.'',/)')
      endif
      NTT=0
      IVlt=0
C Bond Lengths
      IType=1
      Ini=1
      NVar=NLen
      IniP=1
      NVarP=NLenR
      If(DoBPCS) write(IOut,'('' Bond Lengths corrected by BDPCS3'',/)')
      IFill=0
      If(RigB) IFill=1
      Do 40 Ir=1,NLen
       IFixB(Ir)=IFill
   40 Continue
C Stretchings remain primitive R(i,j) coordinates here.  Symmetry labels
C are assigned later by the global GICForge symmetry pass.
      call OrdRed(IOut,IVlt,IPrOrd,MxAtP,MxTrm,DoBPCS,IType,InvDst,NVar,
     $  Ini,IniP,NTermB,IAtomB,IPrimB,ITVB,IFixB,IAn,CoefB,ValTB,C,
     $  ImpDih,Clean)
C Valence Angles
      IType=2
      Ini=NLen+1
      NVar=NAng
      IniP=NLenR+1
      NVarP=NAngR
      IFill=0
      If(RigA) IFill=1
      Do 50 Ir=1,NAng
       IFixA(Ir)=IFill
   50 Continue
      If(SyGNIC) call SymOneGICBlock(IOut,'Bend',MxAtP,MxTrm,2,
     $ NAng,NTermA,IAtomA,ITVA,IFixA,IAn,CoefA)
      call OrdRed(IOut,IVlt,IPrOrd,MxAtP,MxTrm,DoBPCS,IType,.False.,
     $ NVar,Ini,IniP,NTermA,IAtomA,IPrimA,ITVA,IFixA,IAn,CoefA,ValTA,C,
     $ ImpDih,Clean)
C Linear Valence Angles
      IType=3
      Ini=NLen+NAng+1
      NVar=NLAng
      IniP=NLenR+NAngR+1
      NVarP=NLAngR
      IFill=0
      If(RigL) IFill=1
      Do 60 Ir=1,NLAng
       IFixL(Ir)=IFill
   60 Continue
      If(SyGNIC) call SymOneGICBlock(IOut,'Linear bend',MxAtP,MxTrm,3,
     $ NLang,NTermL,IAtomL,ITVLA,IFixL,IAn,CoefL)
      call OrdRed(IOut,IVlt,IPrOrd,MxAtP,MxTrm,DoBPCS,IType,.False.,
     $ NVar,Ini,IniP,NTermL,IAtomL,IPrimL,ITVLA,IFixL,IAn,CoefL,ValTL,C,
     $ ImpDih,Clean)
C Dihedrals
      Itype=4
      Ini=NLen+NAng+NLAng+1
      NVar=NDih
      IniP=NLenR+NAngR+NLAngR+1
      NVarP=NDihR
      IFill=0
      If(RigD) IFill=1
      Do 70 Ir=1,NDih
       IFixD(Ir)=IFill
   70 Continue
      If(SyGNIC) call SymOneGICBlock(IOut,'Torsion',MxAtP,MxTrm,4,
     $ NDih,NTermD,IAtomD,ITVD,IFixD,IAn,CoefD)
      call OrdRed(IOut,IVlt,IPrOrd,MxAtP,MxTrm,DoBPCS,Itype,.False.,
     $ NVar,Ini,IniP,NTermD,IAtomD,IPrimD,ITVD,IFixD,IAn,CoefD,ValTD,C,
     $ ImpDih,Clean)
C Out of Plane
      Itype=5
      Ini=NLen+NAng+NLAng+NDih+1
      NVar=NOuPl
      IniP=NLenR+NAngR+NLAngR+NDihR+1
      NVarP=NOuPlR
      IFill=0
      If(RigO) IFill=1
      Do 80 Ir=1,NOupl
       IFixO(Ir)=IFill
   80 Continue
C     OOP combinations are not printed by PrtOut yet
C     ("Combinations of Improper Dihedrals NYI"). Keep OOP primitives
C     unchanged here; residual redundancies are still pruned by type below.
      call OrdRed(IOut,IVlt,IPrOrd,MxAtP,MxTrm,DoBPCS,Itype,.False.,
     $ NVar,Ini,IniP,NTermO,IAtomO,IPrimO,ITVO,IFixO,IAn,CoefO,ValTO,C,
     $ ImpDih,Clean)
      NLenP=NLen
      NAngP=NAng
      NLangP=NLAng
      NDihP=NDih
      NOuplP=NOupl
      call PruneGICBlocks(IOut,IPrint,MxAtP,MxTrm,NAtoms,NTarget,NLen,
     $  NAng,
     $  NLAng,NOupl,NDih,NTermB,NTermA,NTermL,NTermD,NTermO,IAtomB,
     $  IAtomA,IAtomL,IAtomD,IAtomO,IPrimB,IPrimA,IPrimL,IPrimD,
     $  IPrimO,ITVB,ITVA,ITVLA,ITVD,ITVO,IFixB,IFixA,IFixL,IFixD,
     $  IFixO,CoefB,CoefA,CoefL,CoefD,CoefO,ValTB,ValTA,ValTL,ValTD,
     $  ValTO,C,DoBMat,BMat,Scr,ImpDih,DoLocSVD)
      NTot=NLen+NAng+NLAng+NDih+NOupl
      If(NTot.lt.NTarget.and..not.UsedPrim) then
       Write(IOut,'(/,'' Post-pruning GIC count below vibrational '',
     $ ''rank; retrying with primitive candidates.'')')
       Write(IOut,'(''   Current pruned GICs='',I5,
     $ '' target='',I5,'' primitive candidates='',I5)')
     $ NTot,NTarget,NTotR
       Call UsePrimitiveGICs(IOut,MxAtP,MxTrm,NLenR,NAngR,NLAngR,
     $ NDihR,NOuplR,NLen,NAng,NLAng,NDih,NOupl,IAtmBR,IAtmAR,
     $ IAtmLR,IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,
     $ IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,ITVB,ITVA,ITVLA,ITVD,
     $ ITVO,IFixB,IFixA,IFixL,IFixD,IFixO,CoefB,CoefA,CoefL,
     $ CoefD,CoefO)
       If(DoLocSVD) Call KeepTrueLinearCenters(IOut,MxAtP,MxTrm,
     $  NAtoms,NLAng,NTermL,IAtomL,IPrimL,ITVLA,IFixL,CoefL,ValTL)
       UsedPrim=.True.
       NLenP=NLen
       NAngP=NAng
       NLangP=NLAng
       NDihP=NDih
       NOuplP=NOupl
       call PruneGICBlocks(IOut,IPrint,MxAtP,MxTrm,NAtoms,NTarget,NLen,
     $  NAng,
     $  NLAng,NOupl,NDih,NTermB,NTermA,NTermL,NTermD,NTermO,IAtomB,
     $  IAtomA,IAtomL,IAtomD,IAtomO,IPrimB,IPrimA,IPrimL,IPrimD,
     $  IPrimO,ITVB,ITVA,ITVLA,ITVD,ITVO,IFixB,IFixA,IFixL,IFixD,
     $  IFixO,CoefB,CoefA,CoefL,CoefD,CoefO,ValTB,ValTA,ValTL,ValTD,
     $  ValTO,C,DoBMat,BMat,Scr,ImpDih,DoLocSVD)
      EndIf
      NTot=NLen+NAng+NLAng+NDih+NOupl
      If(NTot.ne.NTarget) then
       Write(IOut,'('' ERROR: final post-pruning non-redundant '',
     $ ''GIC count='',I5,'' differs from target vibrational rank='',
     $ I5)') NTot,NTarget
       Write(IOut,'('' ERROR: no Gaussian GIC coordinate block was '',
     $ ''written because the final basis is not rank complete.'')')
       Stop 1
      EndIf
      Write(IOut,
     $ '(/,28X,''Stretch.  Bend.  L. Bend. Tors.  Out-Pl. Total'')')
      Write(IOut,'('' Redundant  '',6I8)') NLenR,NAngR,NLAngR,NDihR,
     $  NOuplR,NTotR
      Write(IOut,'('' Pre pruning Non Redund.'',6I8)') NLenP,NAngP,
     $  NLangP,NDihP,NOuplP,NLenP+NAngP+NLangP+NDihP+NOuplP
      Write(IOut,'('' Final Non Redund.      '',6I8,/)') NLen,NAng,
     $  NLang,NDih,NOupl,NTot
      If(DoBMat) Write(IOut,'('' Machine-readable final B matrix: '',
     $ ''bmat.out'')')
      NTTsav=NTT
      If(SyGNIC) then
       Write(IOut,'(/,'' Final symmetrized GIC summary '',
     $ ''(Gaussian syntax)'')')
      Else
       Write(IOut,'(/,'' Final GIC summary (Gaussian syntax)'')')
      EndIf
      PrtVal=.True.
      NTTsum=0.0d0
      Call PrtBnd(IOut,MxAtP,MxTrm,InvDst,NLen,NTTsum,NTermB,IAtomB,
     $   ITVB,IFixB,IAn,CoefB,ValTB,C,PrtVal)
      Call PrtAng(IOut,MxAtP,MxTrm,NAng,NTTsum,NTermA,IAtomA,ITVA,
     $  IFixA,CoefA,ValTA,C,PrtVal)
      Call PrtLAn(IOut,MxAtP,MxTrm,NLAng,NTTsum,Linear,NTermL,IAtomL,
     $   ITVLA,IFixL,CoefL,ValTL,C,PrtVal)
      Call PrtDih(IOut,MxAtP,MxTrm,NDih,NTTsum,DoScan,NTermD,IAtomD,
     $   ITVD,IPerD,IFixD,CoefD,ValTD,C,Clean,PrtVal,.True.)
      Call PrtOut(IOut,MxAtP,MxTrm,NOupl,NTTsum,NTermO,IAtomO,ITVO,
     $   IFixO,CoefO,ValTO,C,ImpDih,PrtVal)
      NTT=NTTsav
C Gaussian Input
      if(DoG16.or.DoGDV) then
       OPEN(IPunch,FILE='gauin',STATUS='UNKNOWN')
       Rewind(IPunch)
       call SetGKw(IPunch,DoGNIC,SyGNIC,Loose,DoScan,ModPCS,IDeriv)
       Write(IPunch,'(/,A80,/)') Title
       Write(IPunch,'(2I3)') ICharg,Multip
       do 90 IAt=1,NAtoms
        If(Loose) then
         Write(IPunch,'(I3,3F18.10)') IAn(IAt),(C(ii,IAt),ii=1,3)
        Else
         Write(IPunch,'(I3,3F18.10)') IAn(IAt),(C(ii,IAt),ii=1,3)
        EndIf
   90  continue
       Write(IPunch,'('' '')')
C Bond Lengths
CENZO test
       PrtVal=.False.
CENZO
       Call PrtBnd(IPunch,MxAtP,MxTrm,InvDst,NLen,NTT,NTermB,IAtomB,
     $   ITVB,IFixB,IAn,CoefB,ValtB,C,PrtVal)
C Valence Angles
       Call PrtAng(IPunch,MxAtP,MxTrm,NAng,NTT,NTermA,IAtomA,ITVA,
     $  IFixA,CoefA,ValtA,C,PrtVal)
C Linear Angles
       Call PrtLAn(IPunch,MxAtP,MxTrm,NLAng,NTT,Linear,NTermL,IAtomL,
     $   ITVLA,IFixL,CoefL,ValtL,C,PrtVal)
C Dihedral Angles
       Call PrtDih(IPunch,MxAtP,MxTrm,NDih,NTT,DoScan,NTermD,IAtomD,
     $   ITVD,IPerD,IFixD,CoefD,ValtD,C,Clean,PrtVal,.True.)
C Out of Plane Angles(U) or Impr. Dihedrals (Only option in g16)
       Call PrtOut(IPunch,MxAtP,MxTrm,NOuPl,NTT,NTermO,IAtomO,ITVO,
     $   IFixO,CoefO,ValtO,C,ImpDih,PrtVal)
C Fragments
C     If(NFrag.gt.1) call OrdFrg(IPunch,IPrint,MxAtFr,NAtoms,NFrag,
C    $  IFrag,NAtFr,IFrsAt,IlstAt,IAtFr,LConn)
       If(NFrag.gt.1.and.DoGDV) write(IPunch,'('' Rotor(1)'')')
C rDSD
       If(ModPCS.eq.3) write(IPunch,'(/,''@3F12red.gbs'',/)')
       close(ipunch)
      EndIf
  999 Continue
      write(*,'('' Normal Termination of DiNa25'')')
      write(*,'('' DiNa25   output on file xxx.out'')')
      If(DoG16)  write(*,'('' G16-C01  input  on file xxx.gjf'')')
      If(DoGDV)  write(*,'('' GDV-J28  input  on file xxx.gjf'')')
      End
C===============================================================
C  DRIVER per il calcolo di energie torsionali
C===============================================================
*Deck DrvTrs
      Subroutine DrvTrs(IOut,IPrint,MxBnd,NAtoms,NLen,IAn,IRing,NBond,
     $  IBond,IAtom,EAn,C)
      INTEGER IOut, IPrint, MxBnd, NAtoms, NLen
      INTEGER IAn(NAtoms), IRing(NAtoms), IAtom(2,NAtoms)
      INTEGER NBond(NAtoms), IBond(MxBnd, NAtoms)
      DOUBLE PRECISION EAn(NAtoms),C(3, NAtoms)
C     Write(IOut,'(/,I5,'' Atoms'')') NAtoms
C     do 10 IAt=1,NAtoms
C      write(IOut,'('' Atom'',I5,'' Ian ='',I3,'' EAn ='',F8.4,
C    $   '' forms'',I2,'' bonds with atoms'',4I3)') IAt,IAn(IAt),
C    $   EAn(IAt),NBond(IAt),(IBond(IB,IAt),IB=1,NBond(IAt))
C 10  Continue
C     Write(IOut,'('' IRing'')')
C     Write(IOut,'(20I3)') (IRing(I),I=1,NAtoms)
C     Write(IOut,'('' Atoms Involved in'',I5,'' Bonds'')') NLen
C     Write(IOut,'(2I3)') (Iatom(1,ILen),IAtom(2,ILen),ILen=1,NLen)
      CALL SCAN_ALL_BONDS_RIGID(IOut, IPrint,
     &     MxBnd, NAtoms, IAn, EAn, IRing, NLen, IAtom,
     &     NBond, IBond, C)

      Return
      END
*Deck SetIO
      Subroutine SetIO(InFil,OutFil,Ext)
      Character*(*) InFil, OutFil
      Character*3 Ext
      IStart= LineSt(InFil)
      IEnd=IStart+LinEnd(InFil)
      OutFil=InFil(IStart:IEnd)
      OutFil(IEnd:IEnd)='.'
      OutFil(IEnd+1:IEnd+3)=Ext
      Return
      End
*Deck FndKwd
      Subroutine FndKwd(IOut,IPrint,ModPCS,IDeriv,CLine,Kwd,MxKwd)
      Implicit Real*8 (A-H,O-Z)
      Character*80  CLine
      Character*8   Test
      Logical Kwd(MxKwd),Error,FndPCS,FndDer
      Dimension Ini(50)
      Dimension IDeriv(4)
      call LClear(MxKwd,Kwd)
      FndPCS=.false.
      FndDer=.false.
      ModPCS=-1
      If(CLine(1:1).eq.'#') then
       CLine(1:1)=' '
      Else
       write(IOut,'('' No Input for DiNa25'')')
       Stop
      EndIf
      If(Cline(2:3).eq.'P ') then
       CLine(2:2)=' '
       IPrint=1
      EndIf
      call SubStr(CLine,MxKwd-1,Ini,NumKwd)
      Ini(NumKwd+1)=len(CLine)
      Error=.False.
      do 10 i=1,NumKwd
       I1=Ini(i)
       If(I.eq.NumKwd) then
        I2=I1+7
       Else
        I2=Min0(I1+7,Ini(i+1)-1)
       EndIf
C the default is to normalize GNICs
       Kwd(12)=.true.
       Call LinUpC(CLine(I1:I2),Test)
       If(Test(1:5).eq.'PRINT') then
        IPrint=1
       ElseIf(Test(1:3).eq.'G16') then
        Kwd(2)=.True.
        Kwd(7)=.True.
       ElseIF(Test(1:3).eq.'GDV') then
        Kwd(3)=.True.
       ElseIf(Test(1:5).eq.'CUBIC') then
        Kwd(4)=.True.
       ElseIf(Test(1:7).eq.'NATURAL') then
        Kwd(5)=.True.
       ElseIf(Test(1:6).eq.'ECKART') then
        Kwd(6)=.True.
       ElseIf(Test.eq.'IMPDIH') then
        Kwd(7)=.True.
       ElseIf(Test(1:4).eq.'GNIC') then
        Kwd(8)=.True.
       ElseIf(Test(1:7).eq.'INVDIST') then
        Kwd(9)=.True.
       ElseIf(Test(1:7).eq.'SYMMSTR') then
        Kwd(10)=.True.
       ElseIf(Test(1:6).eq.'ONEDIH') then
        Kwd(11)=.True.
       ElseIf(Test(1:8).eq.'NOONEDIH') then
        Kwd(23)=.True.
       ElseIf(Test(1:6).eq.'NONORM') then
        Kwd(12)=.False.
       ElseIf(Test(1:6).eq.'BDPCS3') then
        Kwd(13)=.True.
       ElseIf(Test(1:7).eq.'SYMMALL') then
        Kwd(14)=.True.
       ElseIf(Test(1:6).eq.'GICSYM') then
        Kwd(24)=.True.
       ElseIf(Test(1:6).eq.'SYCART') then
        Kwd(25)=.True.
       ElseIf(Test(1:6).eq.'FINDFR') then
        Kwd(15)=.True.
       ElseIf(Test(1:6).eq.'JOINFR') then
        Kwd(16)=.True.
       ElseIf(Test(1:5).eq.'HBOND') then
        Kwd(17)=.True.
       ElseIf(Test(1:4).eq.'BMAT') then
        Kwd(19)=.True.
       ElseIf(Test(1:4).eq.'SCAN') then
       Kwd(20)=.True.
      ElseIf(Test(1:5).eq.'RIGID') then
        Kwd(21)=.True.
       ElseIf(Test(1:6).eq.'LOCSVD') then
        Kwd(22)=.True.
      ElseIf(Test(1:5).eq.'CLEAN') then
        Kwd(31)=.true.
       ElseIf(Test(1:5).eq.'LOOSE') then
        Kwd(33)=.true.
       ElseIf(Test(1:3).eq.'OPT') then
        IDeriv(1)=1
       ElseIf(Test(1:4).eq.'FREQ') then
        IDeriv(2)=1
       ElseIf(Test(1:4).eq.'CUBIC') then
        IDeriv(3)=1
       ElseIf(Test(1:6).eq.'ANHARM') then
        IDeriv(4)=1
       ElseIf(.not.FndPCS) then
        call SetPCS(Test,ModPCS,FndPCS)
        If(ModPCS.lt.0) Error=.true.
       EndIf
C       If(.not.FndDer) Call SetDer(Test,IDeriv,FndDer)
C       ITDer=IDeriv(1)+IDeriv(2)+IDeriv(3)+IDeriv(4)
C       If(ModPCS.eq.-1.and.ITDer.eq.0) Error=.true.
C       EndIf
   10 continue
      If(Error) then
       write(IOut,'(/,'' Spurious Keywords in the Input '')')
C      write(IOut,'(A80)') CLine
       write(IOut,'('' The Following Keywords are Allowed'',/)')
       do 20 i=1,MxKwd
        Kwd(i)=.true.
   20  continue
      EndIf
      If(Kwd(2))  write(IOut,'('' G16       : Make G16 Input'')')
      If(Kwd(3))  write(IOut,'('' GDV       : Make GDV Input'')')
      If(Kwd(4))  write(IOut,'('' CUBIC     : Freq=Cubic(GDV) and'',
     $ '' Vibr=Gauconv MSR'')')
      If(Kwd(5))  write(IOut,'('' NATURAL   : Natural Isotopic'',
     $ '' Abundance '')')
      If(Kwd(6))  write(IOut,'('' ECKART    : Enforce Eckart'',
     $  '' Orientation'')')
      If(Kwd(7))  write(IOut,'('' IMPDIH    : Use Improper Dihedrals'',
     $   '' in place of Out-of-plane Bends'')')
      If(Kwd(8))  write(IOut,'('' GNIC      : Make Generalized'',
     $  ''  Natural Internal Coords.'')')
      If(Kwd(9))  write(IOut,'('' INVDIST   : Inv.Dist. for Stretching''
     $  )')
      If(Kwd(10)) write(IOut,'('' SYMMSTR   : Legacy keyword ignored'',
     $ '' for primitive stretchings'')')
      If(.not.Kwd(23)) write(IOut,'('' ONEDIH    : 1 Dihedral per Bond'',
     $ '' (default)'')')
      If(Kwd(23)) write(IOut,'('' NOONEDIH  : Use all non-ring '',
     $ ''dihedrals'')')
      If(.not.Kwd(12)) write(IOut,'('' NONORM    : Not Normalize'',
     $  '' Dihedral GNICS'')')
      If(Kwd(31)) write(IOut,'('' CLEAN     : Clean GNIC values'')')
      If(Kwd(33)) write(IOut,'('' LOOSE     : Loose Symmetry for'',
     $  '' Gaussian'')')
      If(Kwd(13)) write(IOut,'('' BDPCS3    : Make BDPCS3 Bond'',
     $  '' Lengths'')')
      If(Kwd(14)) write(IOut,'('' SYMMALL   : Symmetrize same-type '',
     $ ''GNIC blocks (legacy alias)'')')
      If(Kwd(24)) write(IOut,'('' GICSYM    : Symmetrize GIC blocks'',
     $ '' for downstream modules'')')
      If(Kwd(25)) write(IOut,'('' SYCART    : Write symmetrized '',
     $ ''Cartesian coordinates without changing frame'')')
      If(Kwd(15)) write(IOut,'('' FINDFR    : Find Fragments'')')
      If(Kwd(16)) write(IOut,'('' JOINFR    : Join Fragments'')')
      If(Kwd(17)) write(IOut,'('' HBOND     : Detect H-Bonds only'')')
      If(Kwd(19)) write(IOut,'('' BMAT      : Build B Matrix'')')
      If(Kwd(20)) write(IOut,'('' SCAN      : Scan for Soft DOF'')')
      If(Kwd(21)) write(IOut,'('' RIGID     : Freeze Hard Modes'',
     $  '' in Scan'')')
      If(Kwd(22)) write(IOut,'('' LOCSVD    : Local SVD GNIC blocks'',
     $  '' (experimental)'')')
      If(ModPCS.eq.0.or.error)write(IOut,'('' PCS0      : PCS0 Model'',
     $  '' (UFF)'')')
      If(ModPCS.eq.1.or.error)write(IOut,'('' PCS1      : PCS1 Model'',
     $  '' (HF3C)'')')
      If(ModPCS.eq.2.or.error)write(IOut,'('' HPCS2     : HPCS2 Model'',
     $  '' (B3LYP-D4/6-31G*)'')')
      If(ModPCS.eq.3.or.error)write(IOut,'('' DPCS3     : DPCS3 Model'',
     $  '' (rev-dsd-PBEP86D4/3F12-)'')')
      If(IDeriv(1).eq.1.or.error)write(IOut,'('' OPT       : Geometry'',
     $  ''  Optimization'')')
      If(IDeriv(2).eq.1.or.error)write(IOut,'('' FREQ      : Second'',
     $  '' Energy Derivatives'')')
      If(IDeriv(3).eq.1.or.error)write(IOut,'('' CUBIC     : Third'',
     $  '' Energy Derivatives'')')
      If(IDeriv(4).eq.1.or.error)write(IOut,'('' ANHARM    : Fourth'',
     $  '' Energy Derivatives (up to ijkk'')')
      If(Error) STOP
      Return
      End
*Deck Therm1
      Subroutine Therm1(In,IOut,IPrint,PhyCon,Multip,NAtoms,Linear,
     $ TotWt,ZPE,RTemp,DBVib,Frq)
      Implicit Real*8 (A-H,O-Z)
      Logical Linear
      Save TRef, PRef
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/
      Data TRef/2.9815d2/, PRef/1.0d0/
      Dimension PhyCon(*),RTemp(3),DBVib(3),Frq(*)
      P=PRef
      T=TRef
      pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
      NAt3 = 3*NAtoms
      NVib = NAt3-6
      If(Linear) NVib=NVib+1
C compute electronic partition function
      Boltz = PhyCon(10)
      Avog  = PhyCon(5)
      Gas   = Avog * Boltz
      Degen = Float(Multip)
      QElec = Degen
      SElec = Log(Degen) * Gas
      EElec = Zero
      CElec = Zero
C compute translational partition function
      Call ThrTra(PhyCon,TotWt,P,T,QTran,STran,ETran,CTran)
C compute rotational partition function
      if(NAtoms.eq.1) then
       QRot = one
       ERot = zero
       CRot = zero
       SRot = zero
      Else
       If(Linear) then
        QRot = T/RTemp(1)
        CRot = Gas
        ERot = Gas*T
        SRot = Gas*(Log(QRot)+one)
       Else
        QRot = T/RTemp(1)+T/RTemp(2)+T/RTemp(3)
        CRot = three*gas/two
        ERot = CRot*T
        SRot = GAS*(Log(QROT)+three/two)
       EndIf
      EndIf
      QVib=zero
      QZVib=zero
      QTot=QTran*QRot*QVib*QElec
      QZTot=QTran*QRot*QZVib*QElec
      write(IOut,'(/,5X,'' Partition Functions at '',F7.2,'' K'')') T
      Write(IOut,'(6X,32(''-''),2X)')
      write(IOut,'(5X,'' Translation         '',D12.5)') QTran
      write(IOut,'(5X,'' Rotation            '',D12.5)') QRot
      write(IOut,'(5X,'' Vibration (Bottom)  '',D12.5)') QVib
      write(IOut,'(5X,'' Vibration (from ZPE)'',D12.5)') QZVib
      write(IOut,'(5X,'' Electronic          '',D12.5)') QElec
      Write(IOut,'(6X,32(''-''),2X)')
      write(IOut,'(5X,'' Total (from Bottom) '',D12.5)') QTot
      write(IOut,'(5X,'' Total (from ZPE)    '',D12.5)') QZTot
      return
      end
*Deck ChrNum
      Subroutine ChrNum(String,IStart,IEnd,Number)
      Implicit Real*8 (A-H,O-Z)
C
C     FORM A NUMBER FROM THE Character DIGITS In String.
C
      Character String*80
      Character N*10
      Logical Found
      Save N
      DATA N/'0123456789'/
C
      NTot=IEnd-IStart
      Number=0
      NPow=-1
      do 10 ii=IEnd,IStart,-1
       if(String(ii:ii).eq.' ') goto 10
       found=.false.
       IBas=0
       do 20 j=1,10
        if(found) goto 20
        if(String(ii:ii).eq.N(j:j)) then
         IBas=j-1
         Npow=NPow+1
         found=.true.
        endif
   20  continue
       Number=Number+IBas*10**NPow
   10 continue
      return
      end
*Deck SetPCS
      Subroutine SetPCS(String,ModPCS,FndPCS)
      Character String*8
      Logical FndPCS
      ModPCS=-1
      FndPCS=.False.
      If(String(1:4).eq.'PCS0') then
       ModPCS=0
      ElseIf(String(1:4).eq.'PCS1') then
       ModPCS=1
      ElseIf(String(1:5).eq.'HPCS2') then
       ModPCS=2
      ElseIf(String(1:5).eq.'DPCS3') then
       ModPCS=3
      EndIf
      If(ModPCS.ge.0) FndPCS=.true.
      Return
      End
*Deck SetDer
      Subroutine SetDer(String,IDeriv,FndDer)
      Character String*8
      Logical FndDer
      Integer IDeriv(4)
      Call IClear(4,IDeriv)
      ITot=0
      FndDer=.false.
      If(String(1:3).eq.'OPT') then
       IDeriv(1)=1
       ITot=ITot+1
      EndIf
      If(String(1:4).eq.'FREQ') then
       IDeriv(2)=1
       ITot=ITot+1
      ElseIf(String(1:5).eq.'CUBIC') then
       IDERIV(3)=1
       ITot=ITot+1
      ElseIf(String(1:6).eq.'ANHARM') then
       IDeriv(4)=1
       ITot=ITot+1
      EndIf
C     If(ITot.gt.0) FndDer=.true.
      Return
      End
*Deck SetGKw
      Subroutine SetGKw(IOut,DoGNIC,SyGNIC,Loose,DoScan,ModPCS,IDeriv)
      Logical DoGNIC,SyGNIC,Loose,DoScan
      Integer IOut,ModPCS,IDeriv(4)
      Write(IOut,'(''%Nprocshared=8'')')
      Write(IOut,'(''%Mem=32GB'')')
      Write(IOut,'(''%chk=gicforge.chk'')')
      If(DoGNIC) then
       If(SyGNIC) then
        If(Ideriv(1).eq.1.and.IDeriv(2).eq.0) then
         Write(IOut,'(''#P geom=(readallgic,gicsymm) '')',
     $    advance='no')
        Else
         Write(IOut,'(''#P geom=(readallgic,gicallsymm) '')',
     $    advance='no')
        EndIf
       Else
        Write(IOut,'(''#P geom=readallgic '')',advance='no')
       EndIf
      Else
       Write(IOut,'(''#P'')',advance='no')
      EndIf
      If(Loose) Write(IOut,'(''Symm=loose '')',advance='no')
      If(ModPCS.eq.0) then
       Write(IOut,'(''UFF '')',advance='no')
      ElseIf(ModPCS.eq.1) then
       Write(IOut,'(''HF3C '')',advance='no')
      ElseIf(ModPCS.eq.2) then
       Write(IOut,'(''B3LYP EMPIRICALDISPERSION=GD4 6-31G* '')',
     $   advance='no')
      ElseIf(ModPCS.eq.3) then
       Write(IOut,'(''revDSDPBEP86D4 gen '')',advance='no')
      EndIf
      If(IDeriv(1).eq.1.or.DoScan) then
       If(ModPCS.eq.0) then
        Write(IOut,'(''OPT=nomicro '')',advance='no')
       ElseIf(ModPCS.gt.1.and..not.DoScan) then
        Write(IOut,'(''OPT=calcHFFC '')',advance='no')
       Else
        Write(IOut,'(''OPT '')',advance='no')
       EndIf
      EndIf
      If(IDeriv(2).eq.1) Write(IOut,'(''Freq=IntModes '')',advance='no')
      If(IDeriv(3).eq.1) Write(IOut,'(''Freq=Cubic '')',advance='no')
      If(IDeriv(4).eq.1) Write(IOut,'(''Freq=Anharm '')',advance='no')
      write(IOut,'(''Output=Pickett'')')
      Return
      End
*Deck UsePrimitiveGICs
      Subroutine UsePrimitiveGICs(IOut,MxAtP,MxTrm,NLenR,NAngR,
     $ NLangR,NDihR,NOuplR,NLen,NAng,NLang,NDih,NOupl,IAtmBR,
     $ IAtmAR,IAtmLR,IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,
     $ NTermO,IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,ITVB,ITVA,ITVLA,
     $ ITVD,ITVO,IFixB,IFixA,IFixL,IFixD,IFixO,CoefB,CoefA,CoefL,
     $ CoefD,CoefO)
      Implicit Real*8 (A-H,O-Z)
      Integer MxAtP,MxTrm,NLenR,NAngR,NLangR,NDihR,NOuplR
      Dimension IAtmBR(MxAtP,MxTrm,*),IAtmAR(MxAtP,MxTrm,*)
      Dimension IAtmLR(MxAtP,MxTrm,*),IAtmDR(MxAtP,MxTrm,*)
      Dimension IAtmOR(MxAtP,MxTrm,*)
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension IAtomB(MxAtP,MxTrm,*),IAtomA(MxAtP,MxTrm,*)
      Dimension IAtomL(MxAtP,MxTrm,*),IAtomD(MxAtP,MxTrm,*)
      Dimension IAtomO(MxAtP,MxTrm,*)
      Dimension ITVB(*),ITVA(*),ITVLA(*),ITVD(*),ITVO(*)
      Dimension IFixB(*),IFixA(*),IFixL(*),IFixD(*),IFixO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)

      NLen=NLenR
      NAng=NAngR
      NLang=NLangR
      NDih=NDihR
      NOupl=NOuplR

      Do 10 I=1,NLen
       NTermB(I)=1
       ITVB(I)=0
       IFixB(I)=0
       Do 11 K=1,MxTrm
        CoefB(K,I)=0.0D0
   11  Continue
       CoefB(1,I)=1.0D0
       Do 12 J=1,2
        IAtomB(J,1,I)=IAtmBR(J,1,I)
   12  Continue
   10 Continue

      Do 20 I=1,NAng
       NTermA(I)=1
       ITVA(I)=0
       IFixA(I)=0
       Do 21 K=1,MxTrm
        CoefA(K,I)=0.0D0
   21  Continue
       CoefA(1,I)=1.0D0
       Do 22 J=1,3
        IAtomA(J,1,I)=IAtmAR(J,1,I)
   22  Continue
   20 Continue

      Do 30 I=1,NLang
       NTermL(I)=1
       ITVLA(I)=0
       IFixL(I)=0
       Do 31 K=1,MxTrm
        CoefL(K,I)=0.0D0
   31  Continue
       CoefL(1,I)=1.0D0
       Do 32 J=1,4
        IAtomL(J,1,I)=IAtmLR(J,1,I)
   32  Continue
   30 Continue

      Do 40 I=1,NDih
       NTermD(I)=1
       ITVD(I)=0
       IFixD(I)=0
       Do 41 K=1,MxTrm
        CoefD(K,I)=0.0D0
   41  Continue
       CoefD(1,I)=1.0D0
       Do 42 J=1,4
        IAtomD(J,1,I)=IAtmDR(J,1,I)
   42  Continue
   40 Continue

      Do 50 I=1,NOupl
       NTermO(I)=1
       ITVO(I)=0
       IFixO(I)=0
       Do 51 K=1,MxTrm
        CoefO(K,I)=0.0D0
   51  Continue
       CoefO(1,I)=1.0D0
       Do 52 J=1,4
        IAtomO(J,1,I)=IAtmOR(J,1,I)
   52  Continue
   50 Continue

      Write(IOut,'(''   Primitive fallback active counts:'')')
      Write(IOut,'(''     Stretch='',I5,'' Bend='',I5,
     $ '' Linear='',I5,'' Torsion='',I5,'' Out-of-plane='',I5)')
     $ NLen,NAng,NLang,NDih,NOupl
      Return
      End
*Deck PMOMG
      Subroutine PMOMG(NAtoms,AtMass,C,PMom1)
      Implicit Real*8 (A-H,O-Z)
C Computes the first derivatives of the diagonal elements of
C the moment of inertia tensor w.r.t. Cartesian coordinates
C I/O
C NAtoms            = Number of Atoms
C C(3,NAtoms)       = Cartesian Coordinates
C AtMass(NAtoms)    = Atomic Masses
C PMom1(3,3*NAtoms) = Cartesian first derivatives of the diagonal
      Dimension AtMass(*),C(3,*),PMom1(3,*)
      Zero=0.0d0
      Two=2.0d0
      Do 10 IAt=1,NAtoms
       ix = 3*(i-1) + 1
       iy = 3*(i-1) + 2
       iz = 3*(i-1) + 3
C XX Component
       PMom1(1,ix) = zero
       PMom1(1,iy) = two*atmass(i)*c(2,i)
       PMom1(1,iz) = two*atmass(i)*c(3,i)
C YY Component
       PMom1(2,ix) = two*atmass(i)*c(1,i)
       PMom1(2,iy) = zero
       PMom1(2,iz) = two*atmass(i)*c(3,i)
C ZZ Component
       PMom1(3,ix) = two*atmass(i)*c(1,i)
       PMom1(3,iy) = two*atmass(i)*c(2,i)
       PMom1(3,iz) = zero
   10 Continue
      Return
      End
*Deck PMOMH
      Subroutine PMOMH(NAtoms,AtMass,C,PMom2)
      Implicit Real*8 (A-H,O-Z)
C Computes the second derivatives of the diagonal elements
C of the moment of inertia tensor w.r.t. Cartesian coordinates
C I/O
C NAtoms             = Number of Atoms (Nat3 = 3*NAtoms)
C C(3,NAtoms)        = Cartesian Coordinates
C AtMass(NAtoms)     = Atomic Masses
C PMom2(3,NAt3,NAt3) = Second derivatives of the diagonal elements of
C                      the Inertia Moment Tensor w.r.t. Cartesian Coordinates
      Dimension AtMass(*),C(3,*),PMom2(3,3*Natoms,*)
      Zero=0.0d0
      Two=2.0d0
      Do 10 IAt=1,NAtoms
       ix = 3*(i-1) + 1
       iy = 3*(i-1) + 2
       iz = 3*(i-1) + 3
C XX Component
       PMom2(1,ix,ix) = zero
       PMom2(1,iy,iy) = two*atmass(i)
       PMom2(1,iz,iz) = two*atmass(i)
C YY Component
       PMom2(2,ix,ix) = two*atmass(i)
       PMom2(2,iy,iy) = zero
       PMom2(2,iz,iz) = two*atmass(i)
C ZZ Component
       PMom2(3,ix,ix) = two*atmass(i)
       PMom2(3,iy,iy) = two*atmass(i)
       PMom2(3,iz,iz) = zero
   10 Continue
      Return
      End
