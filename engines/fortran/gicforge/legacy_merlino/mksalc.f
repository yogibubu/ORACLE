*Deck MkGNCB
      Subroutine MkGNCB(IOut,IPrint,DoSymm,InvDst,MxBnd,MxTrm,MxAtP,
     $  NAtoms,IAn,NBond,ILen,IBond,NTermB,IAtomB,ITVB,IAtCyc,CoefB,C)
      Implicit Real*8 (A-H,O-Z)
      Logical DoSymm,InvDst
      Dimension IAn(*),NBond(*),IBond(MxBnd,*),IAtCyc(*)
      Dimension IAtomB(MxAtP,MxTrm,*),NTermB(*),ITVB(*)
      Dimension CoefB(MxTrm,*),C(3,*)
C Local
      Dimension ITest(4),IAnT(4)
C Build Stretchings
      ILen=0
      Do 10 JAt=1,NAtoms
       NBJ=NBond(JAt)
       If(NBJ.eq.1.and.DoSymm) go to 10
       It=0
       call IClear(4,ITest)
       call IClear(4,IAnT)
       Do 20 IB=1,NBond(JAt)
        IAt=IBond(IB,JAt)
        NBI=NBond(IAt)
        If(DoSymm.and.NBI.eq.1) then
         It=It+1
         ITest(It)=IAt
         IAnT(It)=IAn(IAt)
         go to 20
        EndIf
        If(JAt.lt.IAt) then
C Cycles NYI
C         If(DoSymm) then
C          If(IAtCyc(JAt).ne.0.and.IAtCyc(IAt).ne.0) go to 20
C         EndIf
         ILen=ILen+1
         NTermB(ILen)=1
         If(InvDst) then
          ITVB(ILen)=9
         Else
          ITVB(ILen)=0
         EndIf
         IAtomB(1,1,ILen)=JAt
         IAtomB(2,1,ILen)=IAt
         CoefB(1,ILen)=1.0d0
        EndIf
   20  continue
       If(.not.DoSymm) go to 10
C     write(IOut,'(''Atom'',I3,'' Term.Atoms.'',I2)')JAt,IT
       If(It.eq.1) then
        ILen=ILen+1
        NTermB(ILen)=1
        If(InvDst) then
         ITVB(ILen)=9
        Else
         ITVB(ILen)=0
        EndIf
        IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
        IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
        CoefB(1,ILen)=1.0d0
       ElseIf(It.eq.2) then
        Call BLAX2(MxAtP,MxTrm,InvDst,JAt,ILen,ITest,IAnT,IAtomB,NTermB,
     $    ITVB,CoefB)
       ElseIf(It.eq.3) then
        Call BLAX3(MxAtP,MxTrm,InvDst,JAt,ILen,ITest,IAnT,IAtomB,NTermB,
     $    ITVB,CoefB)
       EndIf
   10 continue
      If(IPrint.gt.0) then
       NGICB=ILen
       write(IOut,'(/,I5,'' Bond GNICS'')') NGICB
       if(NGICB.eq.0) return
       Do 30 IGICB=1,NGICB
        NREDB=NTermB(IGICB)
        write(IOut,'('' BGNC('',I3,'')'')') IGICB
        Do 40 IREDB=1,NREDB
         I1=IAtomB(1,IRedB,IGICB)
         I2=IAtomB(2,IRedB,IGICB)
         Coef=CoefB(IRedB,IGICB)
         Value=Distan(C,I1,I2,0)
         write (IOUT,'(F8.4,'' B('',I3,'','',I3,'')'',2X,
     $     ''Value: '',F8.3)') Coef,I1,I2,Value
  40    continue
  30   continue
      endif
      return
      end
*Deck BLAX2
      Subroutine BLAX2(MxAtP,MxTrm,InvDst,JAt,ILen,ITest,IAnT,IAtomB,
     $  NTermB,ITVB,CoefB)
      Implicit Real*8 (A-H,O-Z)
      Logical InvDst
      Dimension ITest(4),IAnT(4)
      Dimension IAtomB(MxAtP,MxTrm,*),NTermB(*),ITVB(*)
      Dimension CoefB(MxTrm,*)
      If(IAnT(1).eq.IAnT(2)) then
       ILen=ILen+1
       NTermB(ILen)=2
       ITVB(ILen)=1
       IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
       IAtomB(1,2,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(2))
       CoefB(1,ILen)=1.0d0/sqrt(2.0d0)
       CoefB(2,ILen)=1.0d0/sqrt(2.0d0)
       ILen=ILen+1
       NTermB(ILen)=2
       ITVB(ILen)=2
       IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
       IAtomB(1,2,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(2))
       CoefB(1,ILen)=1.0d0/sqrt(2.0d0)
       CoefB(2,ILen)=-1.0d0/sqrt(2.0d0)
      Else
       ILen=ILen+1
       NTermB(ILen)=1
       If(InvDst) then
        ITVB(ILen)=9
       Else
        ITVB(ILen)=0
       EndIf
       CoefB(1,ILen)=1.0d0
       IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
       ILen=ILen+1
       NTermB(ILen)=1
       If(InvDst) then
        ITVB(ILen)=9
       Else
        ITVB(ILen)=0
       EndIf
       CoefB(1,ILen)=1.0d0
       IAtomB(1,2,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(2))
      EndIf
      Return
      End
*Deck BLAX3
      Subroutine BLAX3(MxAtP,MxTrm,InvDst,JAt,ILen,ITest,IAnT,IAtomB,
     $  NTermB,ITVB,CoefB)
      Implicit Real*8 (A-H,O-Z)
      Logical InvDst
      Dimension ITest(4),IAnT(4)
      Dimension IAtomB(MxAtP,MxTrm,*),NTermB(*),ITVB(*)
      Dimension CoefB(MxTrm,*)
      SD2=1.0d0/Sqrt(2.0d0)
      SD3=1.0d0/Sqrt(3.0d0)
      SD4=1.0d0/2.0d0
      If(IAnT(1).eq.IAnT(2).and.IAnT(1).eq.IAnT(3)) then
       ILen=ILen+1
       NTermB(ILen)=3
       ITVB(ILen)=3
       IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
       CoefB(1,ILen)=SD3
       IAtomB(1,2,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(2))
       CoefB(2,ILen)=SD3
       IAtomB(1,3,ILen)=Min0(JAt,ITest(3))
       IAtomB(2,3,ILen)=Max0(JAt,ITest(3))
       CoefB(3,ILen)=SD3
       ILen=ILen+1
       NTermB(ILen)=3
       ITVB(ILen)=4
       IAtomB(1,1,ILen)=Min0(JAt,ITest(1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(1))
       CoefB(1,ILen)=Sqrt(2.0d0)*SD4
       IAtomB(1,2,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(2))
       CoefB(2,ILen)=-SD4
       IAtomB(1,3,ILen)=Min0(JAt,ITest(3))
       IAtomB(2,3,ILen)=Max0(JAt,ITest(3))
       CoefB(3,ILen)=-SD4
       ILen=ILen+1
       NTermB(ILen)=2
       ITVB(ILen)=5
       IAtomB(1,1,ILen)=Min0(JAt,ITest(2))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(2))
       CoefB(1,ILen)=SD2
       IAtomB(1,2,ILen)=Min0(JAt,ITest(3))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(3))
       CoefB(2,ILen)=-SD2
       Return
      EndIf
      I1=0
      If(IAnT(1).eq.IAnT(2)) then
       I1=1
       I2=2
       I3=3
      ElseIf(IAnT(1).eq.IAnT(3)) then
       I1=1
       I2=3
       I3=2
      ElseIf(IAnT(2).eq.IAnT(3)) then
       I1=2
       I2=3
       I3=1
      EndIf
      If(I1.ne.0) then
       ILen=ILen+1
       NTermB(ILen)=2
       ITVB(ILen)=1
       IAtomB(1,1,ILen)=Min0(JAt,ITest(I1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(I1))
       IAtomB(1,2,ILen)=Min0(JAt,ITest(I2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(I2))
       CoefB(1,ILen)=1.0d0/sqrt(2.0d0)
       CoefB(2,ILen)=1.0d0/sqrt(2.0d0)
       ILen=ILen+1
       NTermB(ILen)=2
       ITVB(ILen)=2
       IAtomB(1,1,ILen)=Min0(JAt,ITest(I1))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(I1))
       IAtomB(1,2,ILen)=Min0(JAt,ITest(I2))
       IAtomB(2,2,ILen)=Max0(JAt,ITest(I2))
       CoefB(1,ILen)=1.0d0/sqrt(2.0d0)
       CoefB(2,ILen)=-1.0d0/sqrt(2.0d0)
       ILen=ILen+1
       NTermB(ILen)=1
       IAtomB(1,1,ILen)=Min0(JAt,ITest(I3))
       IAtomB(2,1,ILen)=Max0(JAt,ITest(I3))
       CoefB(1,ILen)=1.0d0
      Else
       do 20 ii=1,3
        Ilen=Ilen+ii
        NTermB(ILen)=1
        If(InvDst) then
         ITVB(ILen)=9
        Else
         ITVB(ILen)=0
        EndIf
        CoefB(1,ILen)=1.0d0
        IAtomB(1,1,ILen)=Min0(JAt,ITest(ii))
        IAtomB(2,1,ILen)=Max0(JAt,ITest(ii))
  20   continue
      EndIf
      Return
      End
*Deck MkGNLA
      Subroutine MkGNLA(IOut,IPrint,MxBond,MxGIcL,MxTerL,MaxAtL,
     $  NAtoms,NBond,NGICL,Linear,IBond,NTermL,IAtomL,IAn,CoefL,
     $  C,TreshL,DoLocSVD)
      Implicit Real*8 (A-H,O-Z)
      Dimension C(3,*)
      Dimension NBond(*),IBond(MxBond,*),IAn(*)
      Dimension NTermL(MxGICL),IAtomL(MaxAtL,MxTerL,MxGICL)
      Dimension CoefL(MxTerL,MxGICL)
      Logical Linear,DoLocSVD
      pi = dacos(-1.d0)
      ToDeg=1.80d+2/pi
      NGicL=0
C Build Valence Angles
      Do 30 JAt=1,NAtoms
       NBJ=NBond(JAt)
       NLPair=0
       Do 40 ii=1,NBJ-1
        IAt=IBond(ii,JAt)
        Do 50 kk=ii+1,NBJ
         KAt=IBond(kk,JAt)
         IF(KAt.lt.IAt) then
          LAt=IAt
          IAt=KAt
          KAt=LAt
         EndIf
         Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
         If(Value.lt.TreshL) go to 50
         If(DoLocSVD.and.NBJ.gt.2) go to 50
         If(DoLocSVD.and.NLPair.ge.3) go to 50
         NLPair=NLPair+1
         NGicL=NGicL+1
         NTermL(NGicL)=1
         IAtomL(1,1,NGicL)=IAt
         IAtomL(2,1,NGicL)=JAt
         IAtomL(3,1,NGicL)=KAt
         IAtomL(4,1,NGIcL)=-1
         CoefL(1,NGicL)=1.0d0
         NGicL=NGicL+1
         NTermL(NGicL)=1
         IAtomL(1,1,NGicL)=IAt
         IAtomL(2,1,NGicL)=JAt
         IAtomL(3,1,NGicL)=KAt
         IAtomL(4,1,NGIcL)=-2
         CoefL(1,NGicL)=1.0d0
   50   continue
   40  continue
   30 continue
      If(IPrint.gt.0) then
       write(IOut,'(/,I5,'' Linear Angle GNICS'')') NGICL
       if(NGICL.eq.0) return
       Do 60 IGICL=1,NGICL
        write(IOut,'('' LAngGNC('',I3,'')'')') IGICL
        I1=IAtomL(1,1,IGICL)
        I2=IAtomL(2,1,IGICL)
        I3=IAtomL(3,1,IGICL)
        I4=IAtomL(4,1,IGICL)
        Coef=CoefL(1,IGICL)
        Value=ValAng(C(1,I1),C(1,I2),C(1,I3))
        write (IOUT,'(F8.4,'' LA('',3(I3,'',''),I3,'')'',2X,
     $     ''Value: '',F8.3)') Coef,I1,I2,I3,I4,Value*ToDeg
  60   continue
      endif
      return
      End
*Deck MkGNCD
      Subroutine MkGNCD(IOut,IPrint,MxBnd,MxTrmB,MxTrmD,MxAtB,MxAtD,
     $  MxAtCy,Do1Dih,NAtoms,IAn,NBond,NLen,NDih,NTot,NCyc,IBond,NTermD,
     $  IAtomB,IAtomD,IBr,NAtC,ICAt,IAtCyc,ITVD,IPerD,NEqAt,CoefD,C,EAN,
     $  TreshL,DoNorm)
      Implicit Real*8 (A-H,O-Z)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Dimension NTermD(*),IAtomD(MxAtD,MxTrmD,*)
      Dimension IAtomB(MxAtB,MxTrmB,*),IBr(2,*),ITVD(*),IPerD(*)
      Dimension NatC(*),ICAt(MxAtCy,*),IAtCyc(*),NEqAt(2,*)
      Dimension CoefD(MxTrmD,*),C(3,*),EAN(*)
      Logical join2c,koin2c,join3c,koin3c,Do1Dih,DoNorm
      pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
C Build Dihedral GNICs
      ISoft=0
      Do 10 ILen=1,NLen
       Join2C=.False.
       Koin2C=.False.
       Join3C=.False.
       Koin3c=.False.
       JAt=IAtomB(1,1,ILen)
       KAt=IAtomB(2,1,ILen)
       If(N3Cyc.gt.0) then
        do 20 ii=1,N3Cyc
         ITst=IAt3C(1,ii)
         IF(ITSt.eq.JAt) Join3C=.true.
         If(ITst.eq.KAt) Koin3C=.true.
   20   continue
       EndIf
       if(Join3C.or.Koin3C) go to 10
       IBut=0
       Do 30 jk=1,NBrL
        If(IBut.eq.2) go to 30
        If(JAT.eq.IBrL(1,jk).and.KAt.eq.IBrL(2,jk))IBut=2
        If(JAT.eq.IBrL(2,jk).and.KAt.eq.IBrL(1,jk))IBut=2
   30  Continue
       If(IBut.eq.0.and.IAtCyc(JAt).ne.0.and.IAtCyc(KAt).ne.0) go to 10
       NBJ=NBond(JAt)
       NBK=NBond(KAt)
       if(NBJ.eq.1.or.NBK.eq.1) go to 10
       ITerm=0
       If(IBut.ne.0) then
        call BtFly(IOut,IPrint,MxBnd,MxAtCy,MxAtD,MxTrmD,NDih,JAt,KAt,
     $    NCyc,NBond,IBond,NAtC,ICAt,IAtCyc,ITVD,NTermD,IAtomD,CoefD,C,
     $    TreshL,DoNorm)
       ElseIf(Do1Dih) then
C Deterministic non-ring torsion: one dihedral per central bond.  The
C substituents are selected by local priority, with linear-angle and
C small-ring exclusions applied before the final choice.
        call PickDih(IOut,IPrint,MxBnd,JAt,KAt,NBJ,NBK,IBond,IAtCyc,
     $   IAt,LAt,NEqAt(1,ILen),NEqAt(2,ILen),EAn,C,TreshL,IFLAG)
        If(IFLAG.ne.0) goto 10
        NDih=NDih+1
        Call OrbitDih(IOut,IPrint,MxBnd,MxTrmD,NDih,JAt,KAt,NBJ,NBK,
     $    IBond,IAtCyc,IAn,NBond,IAt,LAt,NTermD,IAtomD,CoefD,C,TreshL,
     $    DoNorm,ITerm)
        NTermD(NDih)=ITerm
CENZO
        NSmx=NEqAt(1,ILen)
        If(NEqAt(2,ILen).gt.NSmx) NSmx=NEqAt(2,ILen)
        If(NSMx.eq.3) then
         ITVD(NDih)=-3
         IPerd(NDih)=3
        ElseIf(NSmx.eq.2) then
         ITVD(NDih)=-2
         IPerD(NDih)=2
        Else
         ITVD(NDIH)=-1
         IPerd(NDih)=1
        EndIf
CENZO
       Else
        Call AllDih(IOut,IPrint,MxBnd,MxAtD,MxTrmD,NDih,JAt,KAt,
     $  NBond,IBond,ITVD,NTermD,IAtomD,CoefD,C,TreshL,DoNorm)
       EndIf
   10 continue
      return
      End
*Deck PickDih
      Subroutine PickDih(IOut,IPrint,MxBnd,JAt,KAt,NBJ,NBK,IBond,
     $ IAtCyc,IAt,LAt,NEqJ,NEqK,EAn,C,TreshL,IFLAG)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,IPrint,MxBnd,JAt,KAt,NBJ,NBK,IAt,LAt,NEqJ,NEqK
      Integer IBond(MxBnd,*),IAtCyc(*),IFLAG
      Integer JJ,KK,ICand,LCand,NBJB,NBKB,IC
      Dimension EAn(*),C(3,*)
      Logical Join2C,Better
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Real*8 BJ,BK,BestJ,BestK,Tresh

      Tresh=5.0d-4
      IFLAG=1
      IAt=0
      LAt=0
      NEqJ=1
      NEqK=1
      BestJ=-1.0d30
      BestK=-1.0d30

      Do 10 JJ=1,NBJ
       ICand=IBond(JJ,JAt)
       If(ICand.eq.KAt) go to 10
       Value=ValAng(C(1,ICand),C(1,JAt),C(1,KAt))
       If(Value.gt.TreshL) go to 10
       Join2C=.False.
       If(N2Cyc.gt.0) then
        Do 20 IC=1,N2Cyc
         If(IAt2C(1,IC).eq.ICand) Join2C=.True.
   20   Continue
       EndIf
       If(Join2C) go to 10
       Do 30 KK=1,NBK
        LCand=IBond(KK,KAt)
        If(LCand.eq.JAt.or.LCand.eq.ICand) go to 30
        Value=ValAng(C(1,JAt),C(1,KAt),C(1,LCand))
        If(Value.gt.TreshL) go to 30
        Join2C=.False.
        If(N2Cyc.gt.0) then
         Do 40 IC=1,N2Cyc
          If(IAt2C(1,IC).eq.LCand) Join2C=.True.
   40    Continue
        EndIf
        If(Join2C) go to 30
        BJ=EAn(ICand)
        BK=EAn(LCand)
        Better=.False.
        If(IFLAG.ne.0) then
         Better=.True.
        ElseIf(BJ.gt.BestJ+Tresh) then
         Better=.True.
        ElseIf(DAbs(BJ-BestJ).lt.Tresh) then
         If(BK.gt.BestK+Tresh) then
          Better=.True.
         ElseIf(DAbs(BK-BestK).lt.Tresh) then
C Stable tie-breakers: prefer more substituted ends, then lower atom labels.
          If(NBJ.gt.NBJB) then
           Better=.True.
          ElseIf(NBJ.eq.NBJB.and.NBK.gt.NBKB) then
           Better=.True.
          ElseIf(NBJ.eq.NBJB.and.NBK.eq.NBKB) then
           If(ICand.lt.IAt) then
            Better=.True.
           ElseIf(ICand.eq.IAt.and.LCand.lt.LAt) then
            Better=.True.
           EndIf
          EndIf
         EndIf
        EndIf
        If(Better) then
         IFLAG=0
         IAt=ICand
         LAt=LCand
         BestJ=BJ
         BestK=BK
         NBJB=NBJ
         NBKB=NBK
        EndIf
   30  Continue
   10 Continue

      If(IFLAG.ne.0) then
       If(IPrint.gt.0) write(IOut,'('' No valid priority torsion '',
     $  ''around bond'',2I5)') JAt,KAt
       Return
      EndIf
      If(IPrint.gt.0) then
       If(NEqJ.gt.1.or.NEqK.gt.1) write(IOut,'('' ONEDIH orbit-closed '',
     $  ''representative around bond'',2I5,'' eq classes'',2I3)')
     $  JAt,KAt,NEqJ,NEqK
      EndIf
      Return
      End
*Deck OrbitDih
      Subroutine OrbitDih(IOut,IPrint,MxBnd,MxTrmD,NDih,JAt,KAt,NBJ,NBK,
     $  IBond,IAtCyc,IAn,NBond,SIAt,SLAt,NTermD,IAtomD,CoefD,C,TreshL,
     $  DoNorm,ITerm)
      Implicit Real*8 (A-H,O-Z)
      Logical DoNorm
      Integer IOut,IPrint,MxBnd,MxTrmD,NDih,JAt,KAt,NBJ,NBK,SIAt,SLAt
      Integer IBond(MxBnd,*),IAtCyc(*),IAn(*),NBond(*),NTermD(*),ITerm
      Integer IAtomD(4,MxTrmD,*)
      Real*8 CoefD(MxTrmD,*),C(3,*),TreshL
      Integer JJ,KK,IAt,LAt
      Real*8 Value

      ITerm=0
      Do 10 JJ=1,NBJ
       IAt=IBond(JJ,JAt)
       If(IAt.eq.KAt) go to 10
       Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
       If(Value.gt.TreshL) go to 10
       If(IAtCyc(IAt).ne.IAtCyc(SIAt)) go to 10
       If(IAn(IAt).ne.IAn(SIAt)) go to 10
       If(NBond(IAt).ne.NBond(SIAt)) go to 10
       Do 20 KK=1,NBK
        LAt=IBond(KK,KAt)
        If(LAt.eq.JAt) go to 20
        Value=ValAng(C(1,JAt),C(1,KAt),C(1,LAt))
        If(Value.gt.TreshL) go to 20
        If(LAt.eq.IAt) go to 20
        If(IAtCyc(LAt).ne.IAtCyc(SLAt)) go to 20
        If(IAn(LAt).ne.IAn(SLAt)) go to 20
        If(NBond(LAt).ne.NBond(SLAt)) go to 20
        If(ITerm.ge.MxTrmD) go to 10
        ITerm=ITerm+1
        IAtomD(1,ITerm,NDih)=IAt
        IAtomD(2,ITerm,NDih)=JAt
        IAtomD(3,ITerm,NDih)=KAt
        IAtomD(4,ITerm,NDih)=LAt
        CoefD(ITerm,NDih)=1.0d0
   20   Continue
   10 Continue
      If(ITerm.eq.0) then
       ITerm=1
       IAtomD(1,1,NDih)=SIAt
       IAtomD(2,1,NDih)=JAt
       IAtomD(3,1,NDih)=KAt
       IAtomD(4,1,NDih)=SLAt
       CoefD(1,NDih)=1.0d0
      EndIf
      If(.not.DoNorm) Return
      Do 30 II=1,ITerm
       CoefD(II,NDih)=CoefD(II,NDih)/SQrt(Float(ITerm))
   30 Continue
      Return
      End
*Deck BtFly
      Subroutine BtFly(IOut,IPrint,MxBnd,MxAtCy,MxAtD,MxTrmD,NDih,JAt,
     $  KAt,NCyc,NBond,IBond,NAtC,ICAt,IAtCyc,ITVD,NTermD,IAtomD,CoefD,
     $  C,TreshL,DoNorm)
      Implicit Real*8 (A-H,O-Z)
      Logical DoNorm,Join2C
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Dimension NBond(*),IBond(MxBnd,*),NTermD(*),IAtomD(MxAtD,MxTrmD,*)
      Dimension ITVD(*),NatC(*),ICAt(MxAtCy,*),IAtCyc(*)
      Dimension CoefD(MxTrmD,*),C(3,*)
      Pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
      ITerm=0
      NBJ=NBond(JAt)
      NBK=NBond(KAt)
      Do 10 jj=1,NBJ
       IAt=IBond(jj,JAt)
       if(IAt.eq.KAt) go to 10
       Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
       If(Value.gt.TreshL) go to 10
       If(IAtCyc(IAt).eq.0) go to 10
       Do 20 kk=1,NBK
        LAt=IBond(kk,KAt)
        If(LAt.eq.JAt) go to 20
        Value=ValAng(C(1,JAt),C(1,KAt),C(1,LAt))
        If(Value.gt.TreshL) go to 20
        If(LAt.eq.IAt) go to 20
        If(IAtCyc(LAt).eq.0) go to 20
        call TstCyc(IAt,LAt,NCyc,MxAtCy,NAtC,ICAt,ISC)
        If(ISC.ne.0) go to 20
        If(ITerm.eq.0) NDih=NDih+1
        ITerm=ITerm+1
        IAtomD(1,ITerm,NDih)=IAt
        IAtomD(2,ITerm,NDih)=JAt
        IAtomD(3,ITerm,NDih)=KAt
        IAtomD(4,ITerm,NDih)=LAt
   20  Continue
   10 Continue
      CoefD(1,NDih)=1.0d0
      CoefD(2,NDih)=-1.0d0
      NTermD(NDih)=ITerm
      ITVD(NDih)=2
      If(.not.DoNorm) Return
      Do 30 ii=1,ITerm
       CoefD(ii,NDih)=CoefD(ii,NDih)/SQrt(Float(ITerm))
   30 Continue
      Return
      End
*Deck AllDih
      Subroutine AllDih(IOut,IPrint,MxBnd,MxAtD,MxTrmD,NDih,JAt,KAt,
     $  NBond,IBond,ITVD,NTermD,IAtomD,CoefD,C,TreshL,DoNorm)
      Implicit Real*8 (A-H,O-Z)
      Logical DoNorm,Join2C
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Dimension NBond(*),IBond(MxBnd,*),NTermD(*),IAtomD(MxAtD,MxTrmD,*)
      Dimension ITVD(*)
      Dimension CoefD(MxTrmD,*),C(3,*)
      Pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
      ITerm=0
      NBJ=NBond(JAt)
      NBK=NBond(KAt)
      Do 10 jj=1,NBJ
       IAt=IBond(jj,JAt)
       if(IAt.eq.KAt) go to 10
       Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
       If(Value.gt.TreshL) go to 10
       Join2C=.false.
       If(N2Cyc.gt.0) then
        do 20 jl=1,N2Cyc
         IF(IAt2C(1,jl).eq.IAt) Join2C=.true.
   20   continue
       EndIf
       if(join2c) go to 10
       Do 30 kk=1,NBK
        LAt=IBond(kk,KAt)
        If(LAt.eq.JAt) go to 30
        Value=ValAng(C(1,JAt),C(1,KAt),C(1,LAt))
        If(Value.gt.TreshL) go to 30
        If(LAt.eq.IAt) go to 30
        Join2C=.false.
        If(N2Cyc.gt.0) then
         do 40 jm=1,N2Cyc
          IF(IAt2C(1,jm).eq.LAt) Join2C=.true.
   40    continue
        EndIf
        if(Join2c) go to 30
        If(ITerm.eq.0) NDih=NDih+1
        ITerm=ITerm+1
        IAtomD(1,ITerm,NDih)=IAt
        IAtomD(2,ITerm,NDih)=JAt
        IAtomD(3,ITerm,NDih)=KAt
        IAtomD(4,ITerm,NDih)=LAt
        CoefD(ITerm,NDih)=1.0d0
   30  Continue
   10 Continue
      NTermD(NDih)=ITerm
      ITVD(NDih)=0
      If(.not.DoNorm) Return
      Do 50 ii=1,ITerm
       CoefD(ii,NDih)=1.0d0/SQrt(Float(ITerm))
   50 Continue
      Return
      End
*Deck MkGNCO
      Subroutine MkGNCO(IOut,IPrint,DoGNIC,MxBond,MxGICO,MxTrmO,MaxAtO,
     $  NAtoms,NCyc,NBond,NGICO,IBond,NTermO,IAtomO,IAn,IAtCyc,CoefO,C,
     $  ImpDih)
      Implicit Real*8 (A-H,O-Z)
      Logical DoGNIC,ImpDih
      Dimension IAn(*),IAtCyc(*),NBond(*),IBond(MxBond,*)
      Dimension NTermO(MxGICO),IAtomO(MaxAtO,MxTrmO,MxGICO)
      Dimension CoefO(MxTrmO,MxGICO),C(3,*)
C
C NBond(NAtoms) = number of bonds for each atom
C IBond(MxBond,NAtoms) = atoms bonded to each atom
C MxGICO = max. number of OUPL GICS;
C MaxAtO = Maximum Number of Atoms in primitives for OUPL GICs
C MxTrmO = Maximum number of primitives in OUPL GICS
C NTermO(MxGICO) = number of terms in OUPL GIC
C NGICO = number of OUPL GICs
C NTermO(NGICO) = number of components for each OUPL GIC
C CoefO(MxTrmO,NGICO) = coefficients of components for each OUPL GIC
C IAtomO(MaxAtO,MxTrmO,NGICO) = atoms defining the components of each OUPL GIC
C
      pi = dacos(-1.d0)
      ToDeg=1.80d+2/pi
      NGICO=0
C Build Out-of-Plane Angles
      Do 10 IAt=1,NAtoms
       If(NBond(IAt).ne.3) go to 10
       IC1=IAtCyc(IAt)
       JAt=IBond(1,IAt)
       IC2=IAtCyc(JAt)
       KAt=IBond(2,IAt)
       IC3=IAtCyc(KAt)
       LAt=IBond(3,IAt)
       IC4=IAtCyc(LAt)
       If(IC1.gt.0.and.IC2.gt.0.and.IC3.gt.0.and.IC4.gt.0) then
        If(IPrint.gt.0) write(IOut,'('' Around'',I5,
     $    '' No out-of-plane for'',3I5)')IAt,JAt,KAt,LAt
        If(DoGNIC) go to 10
       endif
       NGICO=NGICO+1
       NTermO(NGicO)=1
       IAtomO(1,1,NGicO)=IAt
       IAtomO(2,1,NGicO)=JAt
       IAtomO(3,1,NGicO)=KAt
       IAtomO(4,1,NGicO)=LAt
       CoefO(1,NGICO)=1.0d0
   10 continue
      if(IPrint.gt.0) write(IOut,'(/,I5,'' Out-of-Plane GNICS'')')
     $  NGICO
      if(NGICO.eq.0) return
      Do 20 IGICO=1,NGICO
       NREDO=NTermO(IGICO)
       Do 30 IREDO=1,NREDO
        I1=IAtomO(1,IRedO,IGICO)
        I2=IAtomO(2,IRedO,IGICO)
        I3=IAtomO(3,IRedO,IGICO)
        I4=IAtomO(4,IRedO,IGICO)
        Coef=CoefO(IRedO,IGICO)
        If(ImpDih) then
         ValG16=Dihed(C(1,I2),C(1,I1),C(1,I4),C(1,I3))*ToDeg
         If(IPrint.gt.0) write (IOUT,'(I5,''  D('',3(I3,'',''),I3,
     $    '')'',2X,''Value: '',F9.3)')NOupl,JAt,IAt,LAt,KAt,ValG16
        Else
         Value=Outang(C(1,I1),C(1,I2),C(1,I3),C(1,I4))*ToDeg
         If(IPrint.gt.0) write (IOUT,'(I5,''  O('',3(I3,'',''),I3,
     $    '')'',2X,''Value: '',F9.3)')NOupl,JAt,IAt,LAt,KAt,Value
        EndIf
  30   continue
  20  continue
      return
      End
*Deck RedDih
      Subroutine RedDih(IOut,IPrint,MxBnd,ILen,JAt,KAt,NBJ,NBK,IBond,
     $  IAt,LAt,NeqAt,BndOrd,EAn)
      Implicit None
C IO
      Integer IOut,IPrint,MxBnd,ILen,JAt,KAt,NBJ,NBK,IAt,LAt
      Integer IBond(MxBnd,*),NEqAt(2,*)
      Real*8 BndOrd,EAn(*)
C Local
      Integer I1,J1,K1,L1,JI,JJ,KK,KL,IMnI,IAvI,IMxI,IMnL,IAvL,IMxL
      Integer IBJ(3),IBL(3)
      Real*8 ArMax1,ArMin1,EMnI,EAvI,EMxI,EMnL,EAvL,EMxL,Tresh
      Real*8 ETJ(3),ETK(3)
C IAt is the atom bonded to JAt with the largest EAn
C LAt is the atom bonded to KAt with the largest EAn
C If JAt (KAt) has three bonds and two bonded atoms are equal IAt (LAt) is the third one
      Tresh=5.0d-4
      J1=JAt
      K1=KAt
      NEqAt(1,ILen)=1
      NEqAt(2,ILen)=1
      Call IClear(3,IBJ)
      Call IClear(3,IBL)
      Call Aclear(3,ETJ)
      Call AClear(3,ETK)
      I1=IBond(1,J1)
      If(I1.eq.K1) I1=IBond(2,J1)
      If(NBJ.gt.2) then
       ji=0
       do 10 jj=1,NBJ
        I1=IBond(jj,J1)
        if(I1.eq.K1) go to 10
        ji=ji+1
        ETJ(ji)=EAn(I1)
        IBJ(ji)=I1
   10  continue
       If(Abs(ETJ(1)-ETJ(2)).lt.tresh) NEqAt(1,ILen)=2
       EmxI=ArMax1(ETJ,3,.true.,IMxI)
       I1=IBJ(IMxI)
       If(NBJ.eq.4) then
        EmnI=ArMin1(ETJ,3,.true.,IMnI)
        IAvI=6-IMxI-IMnI
        EAvI=ETJ(IAvI)
        If(Abs(EMxI-EMnI).lt.tresh) then
         NEqAt(1,ILen)=3
        ElseIf(Abs(EMnI-EAvI).lt.tresh) then
         NEqAt(1,ILen)=2
        ElseIf(Abs(EMxI-EAvI).lt.tresh) then
         NEqAt(1,ILen)=2
         I1=IBJ(IMnI)
        EndIf
       EndIf
      EndIf
      IAT=I1
      L1=IBond(1,K1)
      If(L1.eq.J1.or.L1.eq.I1) L1=IBond(2,K1)
      NEqAt(2,ILen)=1
      If(NBK.gt.2) then
       kl=0
       do 20 kk=1,NBK
        L1=IBond(kk,K1)
        if(L1.eq.J1.or.L1.eq.I1) go to 20
        kl=kl+1
        ETK(kl)=EAn(L1)
        IBL(kl)=L1
   20  continue
       If(Abs(ETK(1)-ETK(2)).lt.tresh) NEqAt(2,ILen)=2
       EMxL=ArMax1(ETK,3,.true.,IMxL)
       L1=IBL(IMxL)
       If(NBK.eq.4) then
        EMnL=ArMin1(ETK,3,.true.,IMnL)
        IAvL=6-IMxL-IMnL
        EAvL=ETK(IAvL)
        If(Abs(EMxL-EMnL).lt.tresh) then
         NEqAt(2,ILen)=3
        ElseIf(Abs(EMnL-EAvL).lt.tresh) then
         NEqAt(2,ILen)=2
        ElseIf(Abs(EMxL-EAvL).lt.tresh) then
         NEqAt(2,ILen)=2
         L1=IBL(IMnL)
        EndIf
       EndIf
      EndIf
      LAt=L1
      Return
      End
*Deck Gen3At
      Subroutine Gen3At(Iout,IPrint,MxBond,MaxAtA,MxTrmA,ICoord,IAt,
     $  IAn,NBond,IBond,IatCyc,NTermA,IAtomA,CoefA)
      Implicit None
C
C Dimensions
      Integer MxBond, MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,NAng,ICoord,IAt,IAn(*),NBond(*)
      Integer IBond(MxBond,*),IAtCyc(*)
C Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer JAt,KAt,LAt,ICycI,ICycJ,ICycK,ICycL,IProd,IAB,IBC
      Integer I1,I2,I3,I4,IMom1,IMom2,IMom3,IRMax1,IRMin1
      Integer Mom(3),IM(3)
      Real*8 Den1,Den2
C
 9999 Format(/,'GNIC3A: Invalid Atomic Number')
      JAt=IBond(1,IAt)
      KAt=IBond(2,IAt)
      LAt=IBond(3,IAt)
      ICyCI=IAtCyc(IAt)
      ICycJ=IAtCyc(JAt)
      ICycK=IAtCyc(KAt)
      ICycL=IAtCyc(LAt)
      If(ICycI.ne.0.and.ICycJ.ne.0.and.ICycK.ne.0.and.ICycL.ne.0) return
      Mom(1)=IAn(JAt)
      Mom(2)=IAn(KAt)
      Mom(3)=IAn(LAt)
      IProd = Mom(1)*Mom(2)*Mom(3)
      IMom1  = IrMin1(Mom,3,.True.,IM(1))
      IMom3  = IrMax1(Mom,3,.True.,IM(3))
      IM(2)  = 6 - IM(1) - IM(3)
      IMom2  = Mom(IM(2))
      IAB = Abs(IMom1-IMom2)
      IBC = Abs(IMom2-IMom3)
      If(IProd.eq.0) then
        Write(IOut,9999)
        Stop
      elseIf((IMom3-IMom1).eq.0) then
C 3 equal substituents: use 1 different from 2 and 3
       I1=JAt
       I2=IAt
       I3=KAt
       I4=LAt
C 2 equal substituents (for 3 different consider equal the most similar)
      elseif(IAB.lt.IBC) then
       I1=im(3)
       I2=IAt
       I3=im(1)
       I4=im(2)
      elseif(IAB.gt.IBC) then
       I1=im(1)
       I2=IAt
       I3=im(2)
       I4=im(3)
      endif
      If(IAtCyc(I2).eq.0.and.IAtCyc(I3).eq.0) then
       ICoord=ICoord+1
       NTerma(ICoord)=3
       Den1=Sqrt(6.0d0)
       CoefA(1,ICoord)=2.0d0/Den1
       CoefA(2,ICoord)=-1.0d0/Den1
       CoefA(3,ICoord)=-1.0d0/Den1
       IAtomA(1,1,ICoord)=JAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=KAt
       IAtomA(1,2,ICoord)=LAt
       IAtomA(2,2,ICoord)=IAt
       IAtomA(3,2,ICoord)=JAt
       IAtomA(1,3,ICoord)=LAt
       IAtomA(2,3,ICoord)=IAt
       IAtomA(3,3,ICoord)=KAt
      endIf
      Icoord=ICoord+1
      NTermA(ICoord)=2
      Den2=Sqrt(2.0d0)
      CoefA(1,ICoord)=1.0D0/Den2
      CoefA(2,ICoord)=-1.0D0/Den2
      IAtomA(1,1,ICoord)=LAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=JAt
      IAtomA(1,2,ICoord)=LAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=KAt
      return
      end
*Deck OrdRed
      Subroutine OrdRed(IOut,IPunch,IPrint,MaxAtG,MaxTer,DoBPCS,Itp,
     $  InvDst,NVar,Ini,IniP,NTerm,IAtom,IPrim,ITPV,IFixG,IAn,Coef,
     $  ValTot,C,ImpDih,Clean)
C the variables with a final P refer to primitives
      Implicit Real*8 (A-H,O-Z)
      Logical DoBPCS,InvDst,ImpDih,Clean,PrtOrd
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),IPrim(MaxTer,*),ITPV(*)
      Dimension IAn(*),IFixG(*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
      Dimension NAtG(5)
      Character CT*1,CType*5,LbVb*52,Lbl*9,Lbl1*6
      Data CType/'BALDO'/
      DAta LbVb/'SymDRockScisSciLWaggTwstAsyDEEeeT2xxT2yyT2zz B1GEEUU'/
      Data NatG/2,3,4,4,4/
      pi = dacos(-1.d0)
      ToDeg=1.80d+2/pi
      If(NVar.lt.1) return
      PrtOrd=IPrint.ge.0
      IOutW=IOut
      If(.not.PrtOrd) Then
       IOutW=98
       Open(IOutW,File='ordred.tmp',Status='Unknown')
       Rewind(IOutW)
      EndIf
      NAT=NAtG(Itp)
      CT=CType(ITp:ITp)
      IValP=0
      IJKL=0
      I3=0
      I4=0
      Do 10 IGIC=1,NVar
       ValTot(IGic)=0.0d0
       Call LLinCl(Lbl)
       IValG=IGIC+Ini-1
       NRED=NTerm(IGIC)
       ITLab=ITPV(IGIC)
       LBL1=' Free '
       If(IFixG(IGIC).ne.0) LBL1='Frozen'
       If(ITp.eq.1) then
        If(ITLab.eq.0.or.ITLab.gt.10) then
         Lbl(1:7)='Stretch'
        ElseIf(ITLab.eq.1) then
         Lbl(1:9)='XY2 S.Str'
        ElseIf(ITLab.eq.2) then
         Lbl(1:9)='XY2 A.Str'
        ElseIf(ITLab.eq.3) then
         Lbl(1:9)='XY3 S.Str'
        ElseIf(ITLab.eq.4) then
         Lbl(1:9)='XY3 A.Str'
        ElseIf(ITLab.eq.5) then
         Lbl(1:7)='XY3 Def'
        ElseIf(ITLab.eq.6) then
         Lbl(1:7)='XY3 Def'
        ElseIf(ITLab.eq.7) then
         Lbl(1:9)='Ring Brth'
        ElseIf(ITLab.eq.8) then
         Lbl(1:9)='Ring BDef'
        ElseIf(ITLab.eq.9) then
         Lbl(1:)=' 1/R'
        ElseIf(ITLab.eq.10) then
         Lbl(1:6)='H-bond'
        EndIf
       EndIf
       If(ITp.eq.2) then
        If(ITLab.eq.0.or.ITLab.gt.16) then
         Lbl(1:4)='Bend'
        ElseIf(ITLab.eq.1) then
         Lbl(1:9)='Symm Def '
        ElseIf(ITLab.eq.2) then
         Lbl(1:9)='Rocking  '
        ElseIf(ITLab.eq.3) then
         Lbl(1:9)='Scissor  '
        ElseIf(ITLab.eq.4) then
         Lbl(1:9)='Wagging  '
        ElseIf(ITLab.eq.5) then
         Lbl(1:9)='Twisting '
        ElseIf(ITLab.eq.6) then
         Lbl(1:9)='Asymm Def'
        ElseIf(ITLab.eq.7) then
         Lbl(1:9)='Asymm Def'
        ElseIf(ITLab.eq.14) then
         Lbl(1:9)='Ring ADef'
        ElseIf(ITLab.eq.15) then
         Lbl(1:9)='Bic Brid '
        ElseIf(ITLab.eq.16) then
         Lbl(1:9)='Spir Bend'
        Else
         InLb=(ITLab-1)*4+1
         Lbl(1:4)=LbVb(InLb:InLb+3)
        EndIf
       EndIf
       If(ITp.eq.4) then
        If(ITLab.eq.0.or.ITLab.gt.2) then
         Lbl(1:9)='Torsion  '
        ElseIf(ITLab.eq.1) then
         Lbl(1:9)='Ring Puck'
        ElseIf(ITLab.eq.2) then
         Lbl(1:9)='Butterfly'
        EndIf
       EndIf
C       write(IOut,'('' '')')
       Do 20 IRED=1,NRED
        I1=IAtom(1,IRed,IGIC)
        I2=IAtom(2,IRed,IGIC)
        If(NAt.gt.2) I3=IAtom(3,IRed,IGIC)
        If(NAt.gt.3) I4=IAtom(4,IRed,IGIC)
        If(IGIC.eq.1) then
         IValP=IRed
         IJKL=IValP
        else
         call FndRed(MaxAtG,MaxTer,I1,I2,I3,I4,IGIC,IJKL,NTerm,IPrim,
     $    IAtom)
        endif
C       Write(IOut,'(''IJKL='',I3)') IJKL
        If(IJKL.eq.0) then
         IValP=IValP+1
         IJKL=IValP
        EndIf
        IPrim(IRed,IGIC)=IJKL
        CoefP=Coef(IRed,IGIC)
        IValR=IJKL+IniP-1
        if(IPunch.ne.0) Write(IPunch,'(2I4,F10.5)') IValG,IValR,CoefP
        If(Itp.eq.1) then
         Value=Distan(C,I1,I2,0)
         ValTot(IGIC)=ValTot(IGIC)+Value
         if(IRed.eq.1) then
          write (IOutW,'(I3,3X,A9,3X,'' Prim:'',I3,4X,''Atoms'',
     $     2I4,12X,''Value='',F8.4,4X,'' Coeff.='',F8.3,3X,A6)')IValG,
     $     Lbl(1:9),IValR,I1,I2,Value,CoefP,Lbl1
         Else
          write (IOutW,'(18X,'' Prim:'',I3,4X,''Atoms'',2I4,12X,
     $      ''Value='',F8.4,4X,'' Coeff.='',F8.3)')IValR,I1,I2,
     $      Value,CoefP
         EndIf
        ElseIf(ITp.eq.2) then
         ValRad=ValAng(C(1,I1),C(1,I2),C(1,I3))
         ValTot(IGIC)=ValTot(IGIC)+CoefP*ValRad
         Value=ValRad*ToDeg
         if(IRed.eq.1) then
          write (IOutW,'(I3,3X,A9,3X,'' Prim:'',I3,4X,''Atoms'',
     $      3I4,8X,''Value='',F8.3,4X,'' Coeff.='',F8.3,3X,A6)')
     $      IValG,Lbl(1:9),IValR,I1,I2,I3,Value,CoefP,LBL1
         Else
          write (IOutW,'(18X,'' Prim:'',I3,4X,''Atoms'',
     $      3I4,8X,''Value='',F8.3,4X,'' Coeff.='',F8.3)') IValR,I1,
     $      I2,I3,Value,CoefP
         EndIf
        ElseIf(ITp.eq.3) then
         ValRad=ValAng(C(1,I1),C(1,I2),C(1,I3))
         ValTot(IGIC)=ValTot(IGIC)+CoefP*ValRad
         Value=ValRad*ToDeg
         if(IRed.eq.1) then
          write (IOutW,'(I3,'' Lin.Bending   '','' Prim:'',I3,4X,
     $     ''Atoms'',4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3,
     $     3X,A6)') IValG,IValR,I1,I2,I3,I4,Value,CoefP,LBL1
        Else
          write (IOutW,'(18X,'' Prim:'',I3,4X,''Atoms'',
     $      4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3)')IValR,
     $      I1,I2,I3,I4,Value,CoefP
         EndIf
        ElseIf(ITp.eq.4) then
         ValRad=Dihed(C(1,I1),C(1,I2),C(1,I3),C(1,I4))
         ValTot(IGIC)=ValTot(IGIC)+CoefP*ValRad
         Value=ValRad*ToDeg
         if(IRed.eq.1) then
          write (IOutW,'(I3,3X,A9,3X,'' Prim:'',I3,4X,''Atoms'',
     $      4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3,3X,A6)')
     $      IValG,Lbl(1:9),IValR,I1,I2,I3,I4,Value,CoefP,LBL1
         else
          write (IOutW,'(18X,'' Prim:'',I3,4X,''Atoms'',
     $     4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3)')IValR,
     $     I1,I2,I3,I4,Value,CoefP
         EndIf
        ElseIf(Itp.eq.5) then
         ValRad=OutAng(C(1,I1),C(1,I2),C(1,I3),C(1,I4))
         If(ImpDih) ValRad=Dihed(C(1,I2),C(1,I1),C(1,I4),C(1,I3))
         ValTot(IGIC)=ValTot(IGIC)+CoefP*ValRad
         Value=ValRad*ToDeg
         if(IREd.eq.1) then
          If(ImpDih) then
           write (IOutW,'(I3,'' Improper Dih. '','' Prim:'',I3,4X,
     $      ''Atoms'',4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3,
     $      3X,A6)') IValG,IValR,I2,I1,I4,I3,Value,CoefP,LBL1
          Else
           write (IOutW,'(I3,'' Out-of-Plane  '','' Prim:'',I3,4X,
     $      ''Atoms'',4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3,
     $      3X,A6)') IValG,IValR,I1,I2,I3,I4,Value,CoefP,LBL1
          EndIf
         Else
          If(ImpDih) then
           write(IOut,'('' Combinations of Improper Dihedrals NYI'')')
           Stop
          Else
           write (IOutW,'(18X,'' Prim:'',I3,4X,
     $      ''Atoms'',4I4,4X,''Value='',F8.3,4X,'' Coeff.='',F8.3)')
     $      IValR,I1,I2,I3,I4,Value,CoefP
          EndIf
         EndIf
        Else
         write(IOut,'('' From OrdRed: Only IType=1,2,3,4,5 allowed'')')
         stop
        EndIf
  20   continue
C Clean values close to n*pi
       If(ITp.ne.1.and.Clean) then
        ThrAng=2.0d-3
        ValT1=Abs(ValTot(IGIC))
        Test0=ValT1
        Test1=Abs(ValT1-pi)
        Test2=Abs(ValT1-2.0*pi)
        If(Test0.lt.ThrAng) then
         Valt0=0.0d0
        ElseIf(Test1.lt.ThrAng) then
         ValT0=Pi
        ElseIf(Test2.lt.ThrAng) then
         ValT0=0.0d0
        Else
         ValT0=ValTot(IGIC)
        EndIf
        ValTot(IGic)=ValT0
       EndIf
       If(NRed.gt.1)Write(IOutW,'(6X,''Normalized GIC Value ='',F8.3)')
     $  ValTot(IGic)
       write(IOutW,'(100(''-''))')
  10  continue
      If(.not.PrtOrd) Close(IOutW,Status='Delete')
      return
      end

*Deck PrtPckVal
      Subroutine PrtPckVal(IOut,NVar,NTerm,IAtom,ITPV,ValTot)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NVar,ITPV(*),NTerm(*)
      Dimension IAtom(4,15,*)
      Dimension ValTot(*)
      Character S1*4,S2*4
      Integer IVar,IPair,JVar,IV,JV,IV1,JV1,L1,L2
      IPair=0
      IVar=1
   10 If(IVar.gt.NVar) Return
      If(ITPV(IVar).ne.1) Then
       IVar=IVar+1
       Go To 10
      EndIf
      If(IVar.eq.NVar) Return
      If(ITPV(IVar+1).ne.1) Then
       IVar=IVar+1
       Go To 10
      EndIf
      JVar=IVar+1
      IPair=IPair+1
      Call IntoCh(IVar,S1,L1)
      Call IntoCh(JVar,S2,L2)
      IV=IAtom(1,1,IVar)
      JV=IAtom(2,1,IVar)
      IV1=IAtom(1,1,JVar)
      JV1=IAtom(2,1,JVar)
      Write(IOut,'('' RPck'',A4,''(Value='',F12.6,'') = D('',I3,
     $  3('','',I3),'')'')') S1,ValTot(IVar),IAtom(1,1,IVar),
     $  IAtom(2,1,IVar),IAtom(3,1,IVar),IAtom(4,1,IVar)
      Write(IOut,'('' RPck'',A4,''(Value='',F12.6,'') = D('',I3,
     $  3('','',I3),'')'')') S2,ValTot(JVar),IAtom(1,1,JVar),
     $  IAtom(2,1,JVar),IAtom(3,1,JVar),IAtom(4,1,JVar)
      IVar=IVar+2
      Go To 10
      End
*Deck PickPckAtoms
      Subroutine PickPckAtoms(IAtom,NTerm,IVar,I1,I2,I3,I4)
      Implicit Integer (A-Z)
      Integer IAtom(4,15,*),NTerm(*)
      I1=0
      I2=0
      I3=0
      I4=0
      NTrm=NTerm(IVar)
      Do 10 ITr=1,NTrm
       If(IAtom(1,ITr,IVar).ne.0.or.IAtom(2,ITr,IVar).ne.0.or.
     $    IAtom(3,ITr,IVar).ne.0.or.IAtom(4,ITr,IVar).ne.0) then
        I1=IAtom(1,ITr,IVar)
        I2=IAtom(2,ITr,IVar)
        I3=IAtom(3,ITr,IVar)
        I4=IAtom(4,ITr,IVar)
        Return
       EndIf
   10 Continue
      If(NTrm.gt.0) then
       I1=IAtom(1,NTrm,IVar)
       I2=IAtom(2,NTrm,IVar)
       I3=IAtom(3,NTrm,IVar)
       I4=IAtom(4,NTrm,IVar)
      EndIf
      Return
      End
*Deck FndRed
      Subroutine FndRed(MxIAt,MaxTer,IAt,JAt,KAt,LAt,NVar,IJKL,NTerm,
     $  IPrm,IAtP)
      Implicit Integer (A-H,O-Z)
      Dimension NTerm(*),IPrm(MaxTer,*),IAtP(MxIAt,MaxTer,*)
      IJKL=0
      do 10 IVar=1,NVar
       NTrI=NTerm(IVar)
       do 20 ITr=1,NTrI
        if(IAt.ne.IAtP(1,ITr,IVar)) go to 20
        if(JAt.ne.IAtP(2,ITr,IVar)) go to 20
        if(KAt.eq.0) then
         IJKL=IPrm(ITr,IVar)
         return
        else
         if(KAt.ne.IAtP(3,ITr,IVar)) go to 20
         if(LAt.eq.0) then
          IJKL=IPrm(ITr,IVar)
          return
         endif
         if(LAt.ne.IAtP(4,ITr,IVar)) go to 20
         IJKL=IPrm(ITr,IVar)
         return
        endif
   20  continue
   10 continue
      return
      end
*Deck CySalc
      Subroutine CySalc(IOut,IPrint,MxAt,MxAtCy,MxTrm,NAtoms,NLen,
     $  NLenR,NatC,ICAt,IAtCyc,IAtmBr,IAtomB,NTermB,ITVB,CoefB,EAn,C)
      Implicit Real*8 (A-H,O-Z)
      Logical DoJac
      Dimension NatC(*),ICAt(MxAtCy,*)
      Dimension IAtCyc(*)
      Dimension IAtmBR(MxAt,MxTrm,*)
      Dimension IAtomB(MxAt,MxTrm,*),NTermB(*),ITVB(*)
      Dimension CoefB(MxTrm,*)
      Dimension EAn(*),C(3,*)
C Local
      Dimension IX(200),IBb(100),JBb(100)
      Dimension V(1000)
      Tresh=1.0d-5
      IBCyc=0
      NLen0=NLenR
CTemp
      DoJac=.false.
C Determine number of Bonds (and the involved atoms) belonging to Cycles
      Do 10 ILen=1,NLen0
       IAt=IAtmBr(1,1,ILen)
       JAt=IAtmBr(2,1,ILen)
       If(IAtCyc(IAt).eq.0.or.IAtCyc(JAt).eq.0) go to 10
       IBCyc=IBCyc+1
       IBb(IBCyc)=IAt
       JBb(IBCyc)=JAt
   10 Continue
      NBCyc=IBCyc
C Split scratch for diagonalization
C (Jacobi uses a square matrix, whereas HQRII1 uses a low-triangular matrix)
      I1=1
C V(I1) Huckel matrix (lower triangular or square)
      If(DoJac) then
C V(I2) Eigenvalues
       I2=I1+NBCyc*NBCyc
      else
       I2=I1+NBCyc*(NBCyc+1)/2
      EndIf
C V(I3)=EVec
      I3=I2+NBCyc
C V(I4)=WA or Scratch
      I4=I3+NBCyc*NBCyc
      If(DoJac) then
       I6=I4
      Else
C V(I5)=Scr
       I5=I4+6*NBCyc
       I6=I5+NBCyc*(NBCyc+1)/2+1
      EndIf
      If(I6.gt.1000) then
       write(IOut,'('' Outside Memory in CySalc: I6='',I10)') I6
       Stop
      Else
       Call AClear(I6,V)
      EndIf
C Build nearest neighbour matrix (Square for Jacobi and Lower Triangle for HQRII1)
      V(1)=1.0d0
      Do 20 INN=1,NBCyc
       I1At=IBb(INN)
       I2At=JBb(INN)
       Do 30 JNN=1,INN
        J1At=IBb(JNN)
        J2At=JBb(JNN)
        LInd1=(INN-1)*NBCyc+JNN
        Lind2=(JNN-1)*NBCyc+INN
        If(.not.DoJac)LInd1=(INN-1)*INN/2+JNN
        If(INN.eq.JNN) then
         V(LInd1)=1.0d0
         go to 30
        EndIf
        If(I1At.ne.J1At.and.I1At.ne.J2At) then
         If(I2At.ne.J1At.and.I2At.ne.J2At) go to 30
        EndIf
        V(LInd1)=2.0d-01
        If(DoJac) V(LInd2)=2.0d-01
   30  continue
   20 continue
      If(IPrint.gt.1) then
       Write(IOut,'(/,''Nearest Neighbour Stretching Matrix'')')
       do 110 I=1,NBCyc
        If(DoJac) then
         Ini=(I-1)*NBCyc+1
         Iend=Ini+NBCyc-1
        Else
         Ini=(I-1)*I/2+1
         Iend=ini+i-1
        EndIf
        Write(IOut,'(6F10.5)') (V(J),J=ini,iend)
  110  continue
      EndIf
      If(DoJac) then
       ThrDgn=0.0d0
C      call Jacobi(IOut,IPrint,.true.,.true.,0,-1,ThrDgn,NBCyc,0,NBCyc,
C    $   V(I1),V(I2),NBCyc,NBCyc,V(I3))
C      call Eigen(V(I1),V(I3),NBCyc,0)
       do 120 ii=1,NBCyc
        Lind=ii*NBCyc
        V(I2+ii-1)=V(I1+Lind-1)
  120  continue
      Else
       call HQRII1(IOut,NBCyc,1,NBCyc,0,V(I1),V(I2),NBCyc,V(I3),.true.,
     $  IErr,IX,V(I4),V(I5),I6-I5)
       if(IErr.ne.0) then
        write(IOut,'('' HQRII1: IErr='',I5)') IErr
        Stop
       EndIf
      EndIf
      If(IPrint.gt.1) then
       do 100 I=1,NBCyc
        Write(IOut,'(/,'' EigenValue:'',F12.5)') V(I2+I-1)
        Ini=(I-1)*NBCyc
        IEnd=Ini+NBCyc-1
        Write(IOut,'(6F10.5)') (V(I3+ii),ii=Ini,IEnd)
  100  Continue
      EndIf
      do 40 IH=1,NBCyc
       ILen=IH+NLen
       NTermB(ILen)=NBCyc
       ITVB(ILen)=7
       Ini=(IH-1)*NBCyc
       If(IPrint.gt.0) write(IOut,'(/,''Ring Breath'',I3,''('',I3,
     $   '')'')') IH,ILen
       Do 50 JH=1,NBCyc
        IEVec=I3+Ini+JH-1
        CoefB(JH,ILen)=V(IEVec)
        IAtomB(1,JH,ILen)=IBb(JH)
        IAtomB(2,JH,ILen)=JBb(JH)
        If(IPrint.gt.0) write(IOut,'('' Term'',I2,'' Atoms'',2I3,
     $    '' Coeff.'',F10.5)') JH,IBB(JH),JBB(JH),V(IEVec)
   50  Continue
   40 Continue
      NLen=NLen+NBCyc
      Return
      End
*Deck Jacobi
      Subroutine Jacobi(IOut,IPrint,InitV,FulMix,IConv,ISrtE,ThrDgn,N,
     $  NO,NDimA,A,Eig,NDimV,NV,V)
      Implicit Real*8(A-H,O-Z)
C
C     Diagonalize real symmetric matrix a by Jacobi rotations:
C     InitV  ... Whether V is to be initialized to the unit matrix, so
C                that eigenvectors of A are returned, or whether
C                rotations are to performed on the existing elements of V.
C     FulMix ... .True. to diagonalize completely; .False. to avoid
C                mixing vectors which are degenerate to within ThrDgn.
C     IConv  ... 0 full convergence.
C                N remove only off-diagonals within 10^-N of the largest.
C     ISrtE  ... -2/-1/0/1 don't return eigenvalues/return decreasing/
C                return unordered/return increasing.
C     N      ... Dimension of problem.
C     NO     ... If non-zero, only the OV block of A is swept.
C     NDimA  ... Allocated dimension of A.
C     A      ... Input matrix in square form.
C     Eig    ... Eigenvalues.
C     NDimV  ... Allocated dimension of V.
C     NV     ... Leading dimension used in V.  Need not be N if InitV
C                is false.
C     V      ... Eigenvectors.
C
      Parameter (MinPrt=0)
      Logical InitV, FulMix
      Dimension A(NDimA,*), V(NDimV,*), Eig(*)
      Real*8 MDCutO
      Save Zero, Pt25, Pt5, Pt99, One, Two, Five, StpFac
      Data Zero/0.0d0/, Pt25/0.25d0/, Pt5/0.5d0/, Pt99/0.99d0/,
     $  One/1.0d0/, Two/2.0d0/, Five/5.0d0/, StpFac/1.0d-7/
 1000 Format(' Jacobi:  N=',I6,' NSweep=',I6,' Done=',1PD9.2,' Thrsh=',
     $  1PD9.2,' NDid=',I12,0PF6.1,'%')
 1010 Format(' Jacobi:  N=',I6,' NO=',I6,' IConv=',I2,' OffTop=',1PD9.2,
     $  ' ATop=',1PD9.2,' Done=',1PD9.2,' NSweep=',I6,' NTot=',I10)
C
      If(N.eq.1) then
       If(InitV) V(1,1) = One
       If(ISrtE.ge.-1) Eig(1) = A(1,1)
       Return
      endIf
      Call EpsEta(Eps,Eta)
      Small = MDCutO(0)
      If(InitV) then
       If(NV.lt.N) then
        write(IOut,'(/,'' InitV and NV<N in Jacobi.'')')
        STOP
       EndIf
       Do 30 J = 1, N
        Do 20 I = 1, NV
         V(I,J) = Zero
   20   Continue
        V(J,J) = One
   30  Continue
      EndIf
      If(NO.gt.0) then
        AvgF = Max(Float(NO*(N-NO)),One)
        JSt = NO + 1
        IEnd = NO
      else
        AvgF = Float(N*(N-1)/2)
        JSt = 2
        IEnd = N - 1
        endIf
      AvgF = One/AvgF
      NTot = 0
C
C     Find the absolutely largest element of a.
C
C     First check the off-diagonal elements:
      ATop = Zero
      Do 50 J = JSt, N
       Do 55 I = 1, Min(J-1,IEnd)
        ATop = Max(ATop,Abs(A(I,J)))
   55  Continue
   50 Continue
      OffTop = ATop
      Done = Small / Float(10)
      If(IConv.gt.0) Done = Max(Done,OffTop/Float(10)**IConv)
      NSweep = 0
C     Now check the diagonal elements:
      Do 60 J = 1, N
       ATop = Max(ATop,Abs(A(J,J)))
   60 Continue
C
C     If matrix is already effectively diagonal, put diagonal elements
C     in Eig and return.
C
      If(OffTop.lt.Done) goto 800
C
C     Calculate the stopping criterion -- dstop.
C
      D = DStJac(JSt,IEnd,N,NDimA,ATop,A)
      DStop = D*StpFac
C
C     Calculate the threshold.  To make thrsh different than any matrix
C     element of A, multiply by 0.99
C
      Thrsh = Max(Sqrt(D*AvgF)*ATop*Pt99,Done)
C
C     Start a sweep
C
   90 NSweep = NSweep + 1
      NDid = 0
      Do 260 JCol = JSt, N
        Do 250 IRow = 1, Min(JCol-1,IEnd)
          AIJ = A(IRow,JCol)
C
C         Compare the off-diagonal element with Thrsh.
C
          AbsAIJ = Abs(AIJ)
          If(AbsAIJ.ge.Thrsh) then
            AII = A(IRow,IRow)
            AJJ = A(JCol,JCol)
            S = AJJ - AII
            AbsS = Abs(S)
C           Don't rotate IRow and JCol if they would still be degenerate.
            If(FulMix.or.AbsS.ge.ThrDgn.or.AbsAIJ.ge.ThrDgn) then
C
C             Check to see if the chosen rotation is less than the
C             rounding error.
C
              If(AbsAIJ.ge.(Eps*AbsS)) then
                NDid = NDid + 1
C
C               Round if the angle is close to 45 degrees.
C
                Test = Eps*AbsAIJ
                If(AbsS.le.Test) then
                  S = Sqrt(Two)
                  C = S
                else
                  T = AIJ/S
                  S = Pt25/Sqrt(Pt25+T*T)
                  C = Sqrt(Pt5+S)
                  S = Two*T*S/C
                  endIf
C
C               Calculation of the new elements of matrix A.
C
                Do 150 I = 1, IRow
                 T = A(I,IRow)
                 U = A(I,JCol)
                 A(I,IRow) = C*T-S*U
                 A(I,JCol) = S*T+C*U
  150           Continue
                Do 170 I = (IRow+2), JCol
                 T = A(I-1,JCol)
                 U = A(IRow,I-1)
                 A(I-1,JCol) = S*U+C*T
                 A(IRow,I-1) = C*U-S*T
  170           Continue
                A(JCol,JCol) = S*AIJ+C*AJJ
                A(IRow,IRow) = C*A(IRow,IRow)-S*(C*AIJ-S*AJJ)
                Do 190 J = JCol, N
                 T = A(IRow,J)
                 U = A(JCol,J)
                 A(IRow,J) = C*T-S*U
                 A(JCol,J) = S*T+C*U
  190           Continue
C
C               Rotation completed.
C
                Do 210 I = 1, NV
                 T = V(I,IRow)
                 V(I,IRow) = C*T-V(I,JCol)*S
                 V(I,JCol) = S*T+V(I,JCol)*C
  210           Continue
C
C               Calculate the new norm d and compare with DStop.
C               Recalculate DStop and Thrsh to discard rounding errors.
C
                S = AIJ / ATop
                D = D - S*S
                If(D.le.DStop) then
                  D = DStJac(JSt,IEnd,N,NDimA,ATop,A)
                  DStop = D*StpFac
                  endIf
                Thrsh = Min(Max(Sqrt(D*AvgF)*ATop,Done)*Pt99,Thrsh)
                endIf
              endIf
            endIf
  250     Continue
  260   Continue
      NTot = NTot + NDid
      Frac = Float(100)*Float(NDid)*AvgF
      If(IPrint.ge.MinPrt)
     $  Write(IOut,1000) N, NSweep, Done, Thrsh, NDid, Frac
      If(Thrsh.ge.Done) then
        If(NDid.eq.0) then
          If(NDid.eq.0.or.NO.gt.0) Thrsh = Thrsh/Five
        else
          D = DStJac(JSt,IEnd,N,NDimA,ATop,A)
          DStop = D*StpFac
          Thrsh = Min(Max(Sqrt(D*AvgF)*ATop,Done)*Pt99,Thrsh)
          endIf
        Goto 90
        endIf
C
C     Place eigenvalues in Eig and possibly sort.
C
  800 If(ISrtE.ge.-1) then
        Do 810 J = 1, N
         Eig(J) = A(J,J)
  810   Continue
        Call SrtEig(ISrtE,Small,NDimV,NV,N,Eig,V,A)
        endIf
      If(IPrint.ge.MinPrt)
     $  Write(IOut,1010) N, NO, IConv, OffTop, ATop, Done, NSweep, NTot
      Return
      End
*Deck TstHuk
      Subroutine TstHuk(IOut)
      Implicit Real*8 (A-H,O-Z)
      Dimension H(6,6),EVec(6,6),EVal(6)
      Zero=0.0d0
      do 10 I=1,6
       H(I,I)=1.0d0
       If(I.gt.1) then
        H(I,I-1)=0.2d0
        H(I-1,I)=0.2d0
       else
        H(I,6)=0.2d0
        H(6,I)=0.2d0
       endif
   10 continue
      write(IOut,'(/''Puffo: H matrix'')')
      do 15 i=1,6
       write(IOut,'(6F10.5)') (H(I,J),J=1,6)
   15 continue
        IPrint = -1
        IConv = 8
        Zero = 0.0d0
        Call Jacobi(IOut,IPrint,.True.,.True.,IConv,0,Zero,6,0,6,H,
     $    EVal,6,6,EVec)
      do 20 I=1,6
       write(IOut,'(''Eigenvalue'',F10.5)') EVal(i)
       write(IOut,'(6F10.5)') (EVec(J,I),J=1,6)
   20 continue
      return
      end
*Deck PrtOut
      Subroutine PrtOut(IOut,MaxAtG,MaxTer,NVar,NTT,NTerm,IAtom,ITPV,
     $  IFixO,Coef,ValTot,C,ImpDih,PrtVal)
      Implicit Real*8 (A-H,O-Z)
      Logical ImpDih,DoG16,PrtVal
      Character StrVar*4
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),ITPV(*),IFixO(*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      If(NVar.eq.0) return
      do 100 IVar=1,NVar
       Call IntoCh(IVar,StrVar(1:4),Len4)
       NTT=NTT+1
       NTrmI=NTerm(IVar)
       If(NTrmI.eq.1) then
        IAt1=IAtom(1,1,IVar)
        IAt2=IAtom(2,1,IVar)
        IAt3=IAtom(3,1,IVar)
        IAt4=IAtom(4,1,IVar)
C       Value1=OutAngOLd(C(1,IAt2),C(1,IAt1),C(1,IAt3),C(1,IAt4))*ToDeg
        If(ImpDih) then
         Value=Dihed(C(1,IAt2),C(1,IAt1),C(1,IAt4),C(1,IAt3))*ToDeg
         If(IFixO(IVar).eq.0) then
          If(PrtVal) then
            write(IOut,'('' ImpD'',A4,''(Value='',F10.5,'') = D('',I3,
     $       3('','',I3),'')'')')StrVar(1:4),Value,IAt2,IAt1,IAt4,IAt3
          Else
            write(IOut,'('' ImpD'',A4,'' = D('',I3,
     $       3('','',I3),'')'')')StrVar(1:4),IAt2,IAt1,IAt4,IAt3
          EndIf
         Else
          If(PrtVal) then
          write(IOut,'('' ImpD'',A4,''(Frozen,Value='',F10.5,'') = D('',
     $      I3,3('','',I3),'')'')')StrVar(1:4),Value,IAt2,IAt1,IAt4,IAt3
          Else
           write(IOut,'('' ImpD'',A4,''(Frozen) = D('',
     $      I3,3('','',I3),'')'')')StrVar(1:4),IAt2,IAt1,IAt4,IAt3

          EndIf
         EndIf
        Else
         Value=Outang(C(1,IAt1),C(1,IAt2),C(1,IAt3),C(1,IAt4))*ToDeg
         If(IFixO(IVar).eq.0) then
          If(PrtVal) then
          write(IOut,'('' OuPl'',A4,''(Value='',F10.5,'') = U('',I3,
     $     3('','',I3),'')'')') StrVar(1:4),Value,IAt1,IAt2,IAt3,IAt4
          Else
          write(IOut,'('' OuPl'',A4,'' = U('',I3,
     $     3('','',I3),'')'')') StrVar(1:4),IAt1,IAt2,IAt3,IAt4
          EndIf
         Else
          if(PrtVal) then
          write(IOut,'('' OuPl'',A4,''(Frozen,Value='',F10.5,'') = U('',
     $     I3,3('','',I3),'')'')')StrVar(1:4),Value,IAt1,IAt2,IAt3,IAt4
          Else
           write(IOut,'('' OuPl'',A4,''(Frozen) = U('',
     $      I3,3('','',I3),'')'')')StrVar(1:4),IAt1,IAt2,IAt3,IAt4
          EndIf
         EndIf
        EndIf
       EndIf
       go to 100
       If(ImpDih) then
        Write(IOut,'('' Combinations of Improper Dihedrals NYI'')')
        Stop
       EndIf
       write(IOut,'('' UGNIC'',A4,''=['',F12.8,''*U('',3(I3,'',''),I3,
     $  '')'')',advance='no') StrVar(1:4),Coef(1,IVar),
     $  (IAtom(ii,1,IVar),ii=1,4)
       if(NTrmI.gt.2) then
        do 110 i4=2,NTrmI-1
         if(Coef(i4,IVar).gt.0.d0) then
         write(IOut,'(''+'',F12.8,''*U('',3(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,4)
         else
           write(IOut,'(''-'',F12.8,''*U('',3(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,4)
         endif
  110   continue
       endif
  100 continue
      return
      end

*Deck PrtPckQP
      Subroutine PrtPckQP(IOut,NVar,ITPV,ValTot,PrtVal)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NVar,ITPV(*),IPair,IVar,JVar
      Logical PrtVal
      Dimension ValTot(*)
      Character S1*4,S2*4,SP*4
C
C     Derive Gaussian functional GICs from ring-puckering components
C     already printed by PrtDih.  CyGND generates RPck coordinates in
C     consecutive pairs.  Each complete pair defines:
C       Q   = SQRT(RPck_i*RPck_i + RPck_j*RPck_j)
C       Phi = ATAN2(RPck_j, RPck_i)
C
      IPair=0
      IVar=1
   10 If(IVar.gt.NVar) Return
      If(ITPV(IVar).ne.1) Then
       IVar=IVar+1
       Go To 10
      EndIf
      If(IVar.eq.NVar) Return
      If(ITPV(IVar+1).ne.1) Then
       IVar=IVar+1
       Go To 10
      EndIf
      JVar=IVar+1
      IPair=IPair+1
      Call IntoCh(IVar,S1,L1)
      Call IntoCh(JVar,S2,L2)
      Call IntoCh(IPair,SP,LP)
      QVal=DSqrt(ValTot(IVar)*ValTot(IVar)+
     $ ValTot(JVar)*ValTot(JVar))
      PhiVal=DAtan2(ValTot(JVar),ValTot(IVar))
      If(PrtVal) Then
       Write(IOut,1000) SP,QVal,S1,S1,S2,S2
       Write(IOut,1010) SP,PhiVal,S2,S1
      Else
       Write(IOut,1020) SP,S1,S1,S2,S2
       Write(IOut,1030) SP,S2,S1
      EndIf
      IVar=IVar+2
      Go To 10
 1000 Format(' QPck',A4,'(Value=',F10.5,')=SQRT(RPck',A4,
     $ '*RPck',A4,'+RPck',A4,'*RPck',A4,')')
 1010 Format(' PhiP',A4,'(Value=',F10.5,')=ATAN2(RPck',A4,
     $ ',RPck',A4,')')
 1020 Format(' QPck',A4,'=SQRT(RPck',A4,'*RPck',A4,
     $ '+RPck',A4,'*RPck',A4,')')
 1030 Format(' PhiP',A4,'=ATAN2(RPck',A4,',RPck',A4,
     $ ')')
      End
*Deck PrtBnd
      Subroutine PrtBnd(IOut,MaxAtG,MaxTer,InvDst,NVar,NTT,NTerm,IAtom,
     $  ITPV,IFixB,IAn,Coef,Valtot,C,PrtVal)
      Implicit Real*8 (A-H,O-Z)
      Logical InvDst,PrtVal
      Character StrVar*4
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),ITPV(*),IAn(*),IFixB(*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
C Set for MxVar=999
      If(NVar.eq.0) return
      do 100 IVar=1,NVar
       Call IntoCh(IVar,StrVar(1:4),Len4)
       Value=ValTot(IVar)
       NTT=NTT+1
       IAt=IAtom(1,1,IVar)
       JAt=IAtom(2,1,IVar)
       NTrmI=NTerm(IVar)
       If(NTrmI.eq.1) then
        If(IFixB(IVar).eq.0) then
         If(.not.InvDst) then
          If(PrtVal) then
          write(IOut,'('' Stre'',A4,''(Value='',F7.4,'')=R('',I3,'','',
     $      I3,'')'')') StrVar(1:4),Value,IAt,JAt
          Else
           write(IOut,'('' Stre'',A4,'' =R('',I3,'','',
     $       I3,'')'')') StrVar(1:4),IAt,JAt
          EndIf
         Else
          If(PrtVal) then
          write(IOut,'('' SPF'',A4,''(Value='',F7.4,''=1/R('',I3,'','',
     $      I3,'')'')') StrVar(1:4),1.0D0/Value,IAt,JAt
          Else
           write(IOut,'('' SPF'',A4,''=1/R('',I3,'','',
     $       I3,'')'')') StrVar(1:4),IAt,JAt
          EndIf
         EndIf
        Else
         RIJ=Distan(C,IAt,JAt,0)
         If(PrtVal) then
         write(IOut,'('' Stre'',A4,''(Frozen,Value='',F8.5,'') = R('',
     $     I3,'','',I3,'')'')') StrVar(1:4),RIJ,IAt,JAt
         Else
          write(IOut,'('' Stre'',A4,''(Frozen) = R('',
     $      I3,'','',I3,'')'')') StrVar(1:4),IAt,JAt
         EndIf
        endif
        go to 100
       endif
C Inverse distances not used for bond combinations
C      If(InvDst) then
C       write(IOut,'('' Inverse distances not allowed for bond'',
C    $    '' combinations'')')
C       STOP
C      EndIf
       write(IOut,'('' Stre'',A4,''=['',F7.4,''*R('',I3,'','',I3,
     $    '')'')',advance='no') StrVar(1:4),Coef(1,IVar),
     $    (IAtom(ii,1,IVar),ii=1,2)
       if(NTrmI.gt.2) then
        do 110 i4=2,NTrmI-1
         if(Coef(i4,IVar).gt.0.d0) then
          write(IOut,'(''+'',F6.4,''*R('',I3,'','',I3,'')'')',
     $    advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,2)
         else
          write(IOut,'(''-'',F6.4,''*R('',I3,'','',I3,'')'')',
     $    advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,2)
         endif
  110   continue
       endif
       if(Coef(NTrmI,IVar).gt.0.d0) then
         write(IOut,'(''+'',F6.4,''*R('',I3,'','',I3,'')]'')')
     $    DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,2)
       else
        write(IOut,'(''-'',F6.4,''*R('',I3,'','',I3,'')]'')')
     $     DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,2)
       endif
  100 continue
      return
      end
*Deck PrtAng
      Subroutine PrtAng(IOut,MaxAtG,MaxTer,NVar,NTT,NTerm,IAtom,
     $  ITPV,IFixA,Coef,ValTot,C,PrtVal)
      Implicit Real*8 (A-H,O-Z)
      Character LbVb*52,Lbl*4,StrVar*4
      Logical PrtVal
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),ITPV(*),IFixA(*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
      DAta LbVb/'SymDRockScisSciLWaggTwstAsyDEEeeT2xxT2yyT2zz B1GEEUU'/
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
C Set for MxVar=999
      If(NVar.eq.0) return
      do 100 IVar=1,NVar
       Call IntoCh(IVar,StrVar(1:4),Len4)
       If(ITPV(IVar).eq.0) then
        Lbl(1:4)='Bend'
       ElseIf(ITPV(IVar).gt.16) then
        Lbl(1:4)='HCAn'
       ElseIf(ITPV(IVAr).eq.14) then
        Lbl(1:4)='RDef'
       ElseIf(ITPV(IVAr).eq.15) then
        Lbl(1:4)='ByBr'
       ElseIf(ITPV(IVAr).eq.16) then
        Lbl(1:4)='Spir'
       Else
        InLb=(ITPV(IVAR)-1)*4+1
        Lbl(1:4)=LbVb(InLb:InLb+3)
       EndIf
       NTT=NTT+1
       NTrmI=NTerm(IVar)
       Value=ValTot(IVar)
       If(NTrmI.eq.1) then
        Value=ValTot(IVar)*ToDeg
        IAt1=IAtom(1,1,IVar)
       IAt2=IAtom(2,1,IVar)
       IAt3=IAtom(3,1,IVar)
C       Value=ValAng(C(1,IAt1),C(1,IAt2),C(1,IAt3))*ToDeg
        If(ITPV(IVar).gt.16) then
         Lbl(1:4)='HCAn'
        Else
         Lbl(1:4)='Bend'
        EndIf
        If(IFixA(IVar).eq.0) then
         If(PrtVal) then
           write(IOut,'(1X,A4,A4,''(Value='',F10.5,'') = A('',I3,
     $     '','',I3,'','',I3,'')'')') Lbl(1:4),StrVar(1:4),Value,
     $     IAt1,IAt2,IAt3
         Else
           write(IOut,'(1X,A4,A4,'' = A('',I3,
     $     '','',I3,'','',I3,'')'')') Lbl(1:4),StrVar(1:4),
     $     IAt1,IAt2,IAt3
         EndIf
        Else
         If(PrtVal) then
          write(IOut,'(1X,A4,A4,''(Frozen,Value='',F10.5,'') = A('',
     $    I3,'','',I3,'','',I3,'')'')') Lbl(1:4),StrVar(1:4),Value,
     $    IAt1,IAt2,IAt3
         Else
          write(IOut,'(1X,A4,A4,''(Frozen) = A('',
     $     I3,'','',I3,'','',I3,'')'')') Lbl(1:4),StrVar(1:4),
     $     IAt1,IAt2,IAt3
         EndIf
        EndIf
       go to 100
       EndIf
       If(IFixA(IVar).eq.0) then
        If(PrtVal) then
         write(IOut,'(1X,A4,A4,''(Value='',F8.5,'')=['',F12.8,''*A('',
     $    2(I3,'',''),I3,'')'')',advance='no')LBl(1:4),StrVar(1:4),
     $    Value,Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3)
        Else
         write(IOut,'(1X,A4,A4,'' =['',F12.8,''*A('',
     $    2(I3,'',''),I3,'')'')',advance='no')LBl(1:4),StrVar(1:4),
     $    Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3)
        EndIf
       Else
        If(PrtVal) then
         write(IOut,'(1X,A4,A4,''(Frozen,Value='',F8.5,'')=['',F8.5,
     $    ''*A('',2(I3,'',''),I3,'')'')',advance='no')LBl(1:4),
     $    StrVar(1:4),Value,Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3)
        Else
         write(IOut,'(1X,A4,A4,''(Frozen)=['',F8.5,
     $    ''*A('',2(I3,'',''),I3,'')'')',advance='no')LBl(1:4),
     $    StrVar(1:4),Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3)
        EndIf
       EndIf
       if(NTrmI.gt.2) then
        do 110 i4=2,NTrmI-1
         if(Coef(i4,IVar).gt.0.d0) then
          write(IOut,'(''+'',F12.8,''*A('',2(I3,'',''),I3,'')'')',
     $    advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,3)
         else
          write(IOut,'(''-'',F12.8,''*A('',2(I3,'',''),I3,'')'')',
     $    advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,3)
         endif
  110   continue
       endif
       if(Coef(NTrmI,IVar).gt.0.d0) then
        write(IOut,'(''+'',F12.8,''*A('',2(I3,'',''),I3,'')]'')')
     $    DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,3)
       else
        write(IOut,'(''-'',F12.8,''*A('',2(I3,'',''),I3,'')]'')')
     $     DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,3)
       endif
  100 continue
      return
      end
*Deck PrtLAn
      Subroutine PrtLAn(IOut,MaxAtG,MaxTer,NVar,NTT,Linear,NTerm,
     $  IAtom,ITPV,IFixL,Coef,ValTot,C,PrtVal)
      Implicit Real*8 (A-H,O-Z)
      Character StrVar*4
      Logical Linear,PrtVal
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),ITPV(*),IFixL(*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
C Set for MxVar=999
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      If(NVar.eq.0) return
      do 100 IVar=1,NVar
       Call IntoCh(IVar,StrVar(1:4),Len4)
       NTT=NTT+1
       NTrmI=NTerm(IVar)
       I1=IAtom(1,1,IVar)
       I2=IAtom(2,1,IVar)
       I3=IAtom(3,1,IVar)
       I4=-4
       I5=IAtom(4,1,IVar)
       If(Linear) I4=0
       If(NTrmI.eq.1) then
        IAt1=IAtom(1,1,IVar)
        IAt2=IAtom(2,1,IVar)
        IAt3=IAtom(3,1,IVar)
        Value=ValAng(C(1,IAt1),C(1,IAt2),C(1,IAt3))*ToDeg
        If(IFixL(IVar).eq.0) then
         If(PrtVal) then
          write(IOut,'('' LAng'',A4,''(Value='',F8.3,'') = L('',I3,
     $     4('','',I3),'')'')') StrVar(1:4),Value,IAt1,IAt2,IAt3,I4,I5
         Else
          write(IOut,'('' LAng'',A4,'' = L('',I3,
     $     4('','',I3),'')'')') StrVar(1:4),IAt1,IAt2,IAt3,I4,I5
         EndIf
        Else
         If(PrtVal) then
          write(IOut,'('' LAng'',A4,''(Frozen,Value='',F8.3,'') = L('',
     $     I3,4('','',I3),'')'')')StrVar(1:4),Value,IAt1,IAt2,IAt3,I4,I5
         Else
          write(IOut,'('' LAng'',A4,''(Frozen) = L('',
     $     I3,4('','',I3),'')'')')StrVar(1:4),IAt1,IAt2,IAt3,I4,I5
         EndIf
        EndIf
        go to 100
       EndIf
       If(PrtVal) then
        write(IOut,'('' LGIC'',A4,''(Value='',F7.4,'') =['',F12.8,
     $   ''*L('',4(I3,'',''),I3,'')'')',advance='no')StrVar(1:4),
     $   Value,Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3),I4,I5
       Else
        write(IOut,'('' LGIC'',A4,'' =['',F12.8,
     $   ''*L('',4(I3,'',''),I3,'')'')',advance='no')StrVar(1:4),
     $   Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,3),I4,I5
       EndIf
       if(NTrmI.gt.2) then
        do 110 i44=2,NTrmI-1
         if(Coef(i44,IVar).gt.0.d0) then
          write(IOut,'(''+'',F12.8,''*L('',4(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i44,IVar)),(IAtom(ii,i44,IVar),
     $     ii=1,3),Abs(I4),I5
         else
          write(IOut,'(''-'',F12.8,''*L('',4(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i44,IVar)),(IAtom(ii,i44,IVar),
     $     ii=1,3),I4,I5
         endif
  110   continue
       endif
       I5=IAtom(4,NTrmI,IVar)
       if(Coef(NTrmI,IVar).gt.0.d0) then
        write(IOut,'(''+'',F12.8,''*L('',4(I3,'',''),I3,'')]'')')
     $    DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,3),
     $    I4,I5
       else
        write(IOut,'(''-'',F12.8,''*L('',4(I3,'',''),I3,'')]'')')
     $    DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,3),
     $    I4,I5
       endif
  100 continue
      return
      end
*Deck PrtDih
      Subroutine PrtDih(IOut,MaxAtG,MaxTer,NVar,NTT,DoScan,NTerm,IAtom,
     $  ITPV,IPerD,IFixD,Coef,ValTot,C,Clean,PrtVal,PrtRingInt)
      Implicit Real*8 (A-H,O-Z)
      Logical DoScan,Clean,PrtVal,PrtRingInt
      Character Lbl*4,StrVar*4
      Dimension NTerm(*),ITPV(*),IPerD(*),IFixD(*)
      Dimension IAtom(MaxAtG,MaxTer,*)
      Dimension Coef(MaxTer,*),ValTot(*),C(3,*)
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      If(NVar.eq.0) return
      IScan=0
      do 100 IVar=1,NVar
       Call IntoCh(IVar,StrVar(1:4),Len4)
       Value=ValTot(IVar)*ToDeg
CENZO
       NStep=0
       Step=180.0D0
       If(ITpV(IVar).eq.-1) then
        NStep=1
       ElseIf(ITpV(IVar).eq.-2) then
        If(IPerD(IVar).eq.1) then
         NStep=1
        EndIf
       ElseIf(ITpV(IVar).eq.-3) then
        Step=120.0D0
        If(IPerD(IVar).eq.2) then
         NStep=1
        ElseIf(IPerD(IVar).eq.1) then
         NStep=2
        EndIf
       EndIf
       If(.not.Clean) ValRef=Value
       If(ITpV(IVar).eq.1) then
        If(.not.PrtRingInt) go to 100
        Lbl(1:4)='RPck'
       ElseIf(ITpV(IVar).eq.2) then
        Lbl(1:4)='BtFl'
       ElseIf(ITpV(IVar).eq.3) then
        Lbl(1:4)='Libr'
       Else
        Lbl(1:4)='Tors'
       EndIf
       NTT=NTT+1
       NTrmI=NTerm(IVar)
       If(NTrmI.eq.1) then
        Value=ValTot(IVar)*ToDeg
        IAt1=IAtom(1,1,IVar)
        IAt2=IAtom(2,1,IVar)
        IAt3=IAtom(3,1,IVar)
        IAt4=IAtom(4,1,IVar)
C       Value=Dihed(C(1,IAt1),C(1,IAt2),C(1,IAt3),C(1,IAt4))*ToDeg
        If(ITpV(IVar).gt.0) then
         If(PrtVal) then
          write(IOut,'('' Dihe'',A4,''(Value='',F10.5,'') = D('',I3,
     $     3('','',I3),'')'')')StrVar(1:4),Value,IAt1,IAt2,IAt3,IAt4
         Else
          write(IOut,'('' Dihe'',A4,'' = D('',I3,
     $     3('','',I3),'')'')')StrVar(1:4),IAt1,IAt2,IAt3,IAt4
         EndIf
        Else
         If(DoScan.and.NStep.gt.0) then
          IScan=IScan+1
          write(IOut,'('' Scan'',A4,''(Value='',F10.5,'', NSteps='',
     $      I2,'', StepSize='',F5.1,'') = D('',I3,3('','',I3),'')'')')
     $      StrVar(1:4),ValRef,NStep,Step,IAt1,IAt2,IAt3,IAt4
         Else
          If(PrtVal) then
            write(IOut,'('' Dihe'',A4,''(Value='',F10.5,'') = D('',
     $      I3,3('','',I3),'')'')')StrVar(1:4),Value,IAt1,IAt2,IAt3,IAt4
          Else
            write(IOut,'('' Dihe'',A4,'' = D('',
     $      I3,3('','',I3),'')'')')StrVar(1:4),IAt1,IAt2,IAt3,IAt4
          EndIf
         EndIf
        EndIf
        go to 100
       EndIf
       Value=ValTot(IVar)
       If(ITPV(IVar).eq.1) Then
        If(PrtVal) then
         write(IOut,'(1X,A4,A4,''(Inactive,Value='',F10.5,
     $    '')=['',F12.8,''*D('',3(I3,'',''),I3,'')'')',
     $    advance='no') Lbl(1:4),StrVar(1:4),Value,Coef(1,IVar),
     $    (IAtom(ii,1,IVar),ii=1,4)
        Else
         write(IOut,'(1X,A4,A4,''(Inactive)=['',F12.8,
     $    ''*D('',3(I3,'',''),I3,'')'')',advance='no')
     $    Lbl(1:4),StrVar(1:4),Coef(1,IVar),
     $    (IAtom(ii,1,IVar),ii=1,4)
        EndIf
       Else
        If(PrtVal) then
         write(IOut,'(1X,A4,A4,''(Value='',F10.5,'')=['',F12.8,
     $    ''*D('',3(I3,'',''),I3,'')'')',advance='no')
     $    Lbl(1:4),StrVar(1:4),Value,Coef(1,IVar),
     $    (IAtom(ii,1,IVar),ii=1,4)
        Else
         write(IOut,'(1X,A4,A4,'' =['',F12.8,''*D('',
     $    3(I3,'',''),I3,'')'')',advance='no') Lbl(1:4),
     $    StrVar(1:4),Coef(1,IVar),(IAtom(ii,1,IVar),ii=1,4)
        EndIf
       EndIf
       if(NTrmI.gt.2) then
        do 110 i4=2,NTrmI-1
         if(Coef(i4,IVar).gt.0.d0) then
          write(IOut,'(''+'',F12.8,''*D('',3(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,4)
         else
          write(IOut,'(''-'',F12.8,''*D('',3(I3,'',''),I3,'')'')',
     $     advance='no') DAbs(Coef(i4,IVar)),(IAtom(ii,i4,IVar),ii=1,4)
         endif
  110   continue
       endif
       if(Coef(NTrmI,IVar).gt.0.d0) then
        write(IOut,'(''+'',F12.8,''*D('',3(I3,'',''),I3,'')]'')')
     $   DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,4)
       else
        write(IOut,'(''-'',F12.8,''*D('',3(I3,'',''),I3,'')]'')')
     $   DAbs(Coef(NTrmI,IVar)),(IAtom(ii,NTrmI,IVar),ii=1,4)
       endif
  100 continue
      Call PrtPckQP(IOut,NVar,ITPV,ValTot,PrtVal)
      return
      end
