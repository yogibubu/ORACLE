*Deck GICPrune
C=======================================================================
C  GICPRUNE
C
C  Final type-local redundancy pruning for GICForge.
C
C  The reduction is intentionally block-local: bends, linear bends,
C  torsions, and out-of-plane coordinates are pruned independently.
C  Stretchings are kept as primitive bond coordinates and are not
C  pruned here.
C
C  Method: modified Gram-Schmidt over the B rows of one block.  A coordinate is
C  retained if its B row increases the numerical rank of its own type block.
C=======================================================================

      Subroutine PruneGICBlocks(IOut,IPrint,MxAtP,MxTrm,NAtoms,NTarget,
     $ NLen,
     $ NAng,NLAng,NOupl,NDih,NTermB,NTermA,NTermL,NTermD,NTermO,
     $ IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,IPrimB,IPrimA,IPrimL,
     $ IPrimD,IPrimO,ITVB,ITVA,ITVLA,ITVD,ITVO,IFixB,IFixA,IFixL,
     $ IFixD,IFixO,CoefB,CoefA,CoefL,CoefD,CoefO,ValTB,ValTA,ValTL,
     $ ValTD,ValTO,C,DoBMat,BMat,Scr,ImpDih,DoLocSVD)
      Implicit Real*8 (A-H,O-Z)
      Integer MxAtP,MxTrm,NAtoms,NTarget,NLen,NAng,NLAng,NOupl,NDih
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension IAtomB(MxAtP,MxTrm,*),IAtomA(MxAtP,MxTrm,*)
      Dimension IAtomL(MxAtP,MxTrm,*),IAtomD(MxAtP,MxTrm,*)
      Dimension IAtomO(MxAtP,MxTrm,*)
      Dimension IPrimB(MxTrm,*),IPrimA(MxTrm,*),IPrimL(MxTrm,*)
      Dimension IPrimD(MxTrm,*),IPrimO(MxTrm,*)
      Dimension ITVB(*),ITVA(*),ITVLA(*),ITVD(*),ITVO(*)
      Dimension IFixB(*),IFixA(*),IFixL(*),IFixD(*),IFixO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)
      Dimension ValTB(*),ValTA(*),ValTL(*),ValTD(*),ValTO(*)
      Dimension C(3,*),BMat(3*NAtoms,*),Scr(*)
      Logical KeepA(1000),KeepL(1000),KeepD(1000),KeepO(1000)
      Logical DoB1,DoBMat,ImpDih,DoLocSVD
      Integer CountKeep
      Character*16 Label

      If(NAtoms.le.0) return
      NTot=NLen+NAng+NLAng+NDih+NOupl
      If(NTot.le.0) return
      NVib=NTarget

      DoB1=.False.
      Call MkBNew(IOut,0,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,NLAng,
     $ NOupl,NDih,IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,NTermB,NTermA,
     $ NTermL,NTermD,NTermO,CoefB,CoefA,CoefL,CoefD,CoefO,C,BMat,
     $ ImpDih)

      Write(IOut,'(/,'' Type-local residual GIC redundancy pruning'')')
      Write(IOut,'(''   Method: modified Gram-Schmidt on B rows;'',
     $ '' coordinate types are not mixed.'')')

      Write(IOut,'(''   Stretch: kept all '',I5,
     $ '' primitive coordinates.'')') NLen

      If(NTot.lt.NVib) then
       Write(IOut,'(''   Current GIC count='',I5,
     $ '' target vibrational rank='',I5)') NTot,NVib
       Write(IOut,'(''   Candidate counts: Stretch='',I5,'' Bend='',I5,
     $ '' Linear='',I5,'' Torsion='',I5,'' Out-of-plane='',I5)')
     $ NLen,NAng,NLAng,NDih,NOupl
       Write(IOut,'('' ERROR: GIC set is below the vibrational rank '',
     $ ''before pruning.'')')
       Write(IOut,'('' ERROR: LOCSVD/full primitive fallback must '',
     $ ''generate at least target-rank candidates before pruning.'')')
       Stop 1
      EndIf

      If(NTot.eq.NVib) then
       Write(IOut,'(''   Current GIC count='',I5,
     $ '' target vibrational rank='',I5)') NTot,NVib
       NBasis=0
       Call SeedGICBasis(IOut,NAtoms,NTot,0,BMat,Scr,NBasis)
       If(NBasis.eq.NVib) Then
        Write(IOut,'(''   No pruning performed: the set is already '',
     $  ''minimal and full rank.'')')
        Write(IOut,'(''   Final active GIC counts:'')')
        Write(IOut,'(''     Stretch='',I5,'' Bend='',I5,
     $  '' Linear='',I5,'' Torsion='',I5,'' Out-of-plane='',I5)')
     $  NLen,NAng,NLAng,NDih,NOupl
        If(DoBMat) then
         Call WriteGICBMat(NAtoms,NTot,BMat)
        EndIf
        Return
       EndIf
       Write(IOut,'(''   Minimal-count set has B-rank '',I5,
     $ '' / '',I5,''; pruning dependent rows before fallback.'')')
     $ NBasis,NVib
      EndIf

      NAng0=NAng
      NLang0=NLAng
      NDih0=NDih

      If(.not.DoLocSVD) then
       NBasis=0
       Call SeedGICBasis(IOut,NAtoms,NLen,0,BMat,Scr,NBasis)

       Label='Bend'
       IOff=NLen
       Call PruneOneBlockAgainst(IOut,Label,NAtoms,NAng,IOff,BMat,
     $  KeepA,Scr,NBasis,NVib,NKeep)
       If(NKeep.lt.NAng) Call PackGICBlock(MxAtP,MxTrm,NAng,KeepA,
     $  NTermA,IAtomA,IPrimA,ITVA,IFixA,CoefA,ValTA)

       Label='Linear bend'
       IOff=NLen+NAng0
       Call PruneOneBlockAgainst(IOut,Label,NAtoms,NLAng,IOff,BMat,
     $  KeepL,Scr,NBasis,NVib,NKeep)
       If(NKeep.lt.NLAng) Call PackGICBlock(MxAtP,MxTrm,NLAng,KeepL,
     $  NTermL,IAtomL,IPrimL,ITVLA,IFixL,CoefL,ValTL)

       Label='Torsion'
       IOff=NLen+NAng0+NLang0
       Call PruneOneBlockAgainst(IOut,Label,NAtoms,NDih,IOff,BMat,
     $  KeepD,Scr,NBasis,NVib,NKeep)
       If(NKeep.lt.NDih) Call PackGICBlock(MxAtP,MxTrm,NDih,KeepD,
     $  NTermD,IAtomD,IPrimD,ITVD,IFixD,CoefD,ValTD)

       Label='Out-of-plane'
       IOff=NLen+NAng0+NLang0+NDih0
       Call PruneOneBlockAgainst(IOut,Label,NAtoms,NOupl,IOff,BMat,
     $  KeepO,Scr,NBasis,NVib,NKeep)
       If(NKeep.lt.NOupl) Call PackGICBlock(MxAtP,MxTrm,NOupl,KeepO,
     $  NTermO,IAtomO,IPrimO,ITVO,IFixO,CoefO,ValTO)
       GoTo 500
      EndIf

      Call InitKeep(NAng,KeepA)
      Call InitKeep(NLAng,KeepL)
      Call InitKeep(NDih,KeepD)
      Call InitKeep(NOupl,KeepO)
      NBasis=0
      Call SeedGICBasis(IOut,NAtoms,NLen,0,BMat,Scr,NBasis)

      Write(IOut,'(''   Priority order: Stretch, Exocyclic bend, '',
     $ ''Linear, Exocyclic torsion, Butterfly, Ring bend, Ring '',
     $ ''torsion, Out-of-plane.'')')

      Label='Exocyclic bend'
      IOff=NLen
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NAng,IOff,
     $ BMat,ITVA,4,KeepA,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Linear bend'
      IOff=NLen+NAng0
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NLAng,IOff,
     $ BMat,ITVLA,1,KeepL,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Exocyclic torsion'
      IOff=NLen+NAng0+NLang0
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NDih,IOff,
     $ BMat,ITVD,2,KeepD,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Butterfly torsion'
      IOff=NLen+NAng0+NLang0
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NDih,IOff,
     $ BMat,ITVD,3,KeepD,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Ring bend'
      IOff=NLen
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NAng,IOff,
     $ BMat,ITVA,5,KeepA,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Ring torsion'
      IOff=NLen+NAng0+NLang0
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NDih,IOff,
     $ BMat,ITVD,6,KeepD,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      Label='Out-of-plane'
      IOff=NLen+NAng0+NLang0+NDih0
      Call PruneSelectedBlockAgainst(IOut,Label,NAtoms,NOupl,IOff,
     $ BMat,ITVO,7,KeepO,Scr,NBasis,NVib,NKeep,NSelect)
      Call ReportPruneProgress(IOut,Label,NBasis,NVib)

      If(CountKeep(NAng,KeepA).lt.NAng) Call PackGICBlock(MxAtP,MxTrm,
     $ NAng,KeepA,NTermA,IAtomA,IPrimA,ITVA,IFixA,CoefA,ValTA)
      If(CountKeep(NLAng,KeepL).lt.NLAng) Call PackGICBlock(MxAtP,
     $ MxTrm,NLAng,KeepL,NTermL,IAtomL,IPrimL,ITVLA,IFixL,CoefL,ValTL)
      If(CountKeep(NDih,KeepD).lt.NDih) Call PackGICBlock(MxAtP,MxTrm,
     $ NDih,KeepD,NTermD,IAtomD,IPrimD,ITVD,IFixD,CoefD,ValTD)
      If(CountKeep(NOupl,KeepO).lt.NOupl) Call PackGICBlock(MxAtP,
     $ MxTrm,NOupl,KeepO,NTermO,IAtomO,IPrimO,ITVO,IFixO,CoefO,ValTO)

  500 Continue
      Write(IOut,'(''   Final active GIC counts:'')')
      Write(IOut,'(''     Stretch='',I5,'' Bend='',I5,'' Linear='',I5,
     $ '' Torsion='',I5,'' Out-of-plane='',I5)') NLen,NAng,NLAng,NDih,
     $ NOupl
      NTot=NLen+NAng+NLAng+NDih+NOupl
      If(NTot.ne.NVib) then
       Write(IOut,'('' ERROR: post-pruning GIC count='',I5,
     $ '' differs from target vibrational rank='',I5)') NTot,NVib
       Write(IOut,'('' ERROR: interim block counts: Stretch='',I5,
     $ '' Bend='',I5,'' Linear='',I5,'' Torsion='',I5,
     $ '' Out-of-plane='',I5)') NLen,NAng,NLAng,NDih,NOupl
       Return
      EndIf
      If(DoBMat) then
       Call MkBNew(IOut,0,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,NLAng,
     $ NOupl,NDih,IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,NTermB,NTermA,
     $ NTermL,NTermD,NTermO,CoefB,CoefA,CoefL,CoefD,CoefO,C,BMat,
     $ ImpDih)
       Call WriteGICBMat(NAtoms,NTot,BMat)
      EndIf
      Return
      End

*Deck ReportPruneProgress
      Subroutine ReportPruneProgress(IOut,Label,NBasis,NTarget)
      Integer IOut,NBasis,NTarget
      Character*(*) Label
      Write(IOut,'(''     Rank after '',A,'': '',I5,'' / '',I5)')
     $ Label,NBasis,NTarget
      If(NBasis.ge.NTarget) Write(IOut,'(''     Target reached after '',
     $ A)') Label
      Return
      End

*Deck InitKeep
      Subroutine InitKeep(NVar,Keep)
      Integer NVar,I
      Logical Keep(*)
      Do 10 I=1,NVar
       Keep(I)=.False.
   10 Continue
      Return
      End

*Deck CountKeep
      Integer Function CountKeep(NVar,Keep)
      Integer NVar,I
      Logical Keep(*)
      CountKeep=0
      Do 10 I=1,NVar
       If(Keep(I)) CountKeep=CountKeep+1
   10 Continue
      Return
      End

*Deck KeepTrueLinearCenters
      Subroutine KeepTrueLinearCenters(IOut,MxAtP,MxTrm,NAtoms,NLAng,
     $ NTermL,IAtomL,IPrimL,ITVLA,IFixL,CoefL,ValTL)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,MxAtP,MxTrm,NAtoms,NLAng
      Integer CountKeep
      Dimension NTermL(*),IAtomL(MxAtP,MxTrm,*),IPrimL(MxTrm,*)
      Dimension ITVLA(*),IFixL(*),CoefL(MxTrm,*),ValTL(*)
      Dimension NPair(1000)
      Logical Keep(1000)
      If(NAtoms.gt.1000) Return
      Do 10 I=1,NAtoms
       NPair(I)=0
   10 Continue
      Do 20 I=1,NLAng
       Keep(I)=.False.
       JAt=IAtomL(2,1,I)
       Mode=IAtomL(4,1,I)
       If(Mode.eq.-1.and.JAt.ge.1.and.JAt.le.NAtoms)
     $  NPair(JAt)=NPair(JAt)+1
   20 Continue
      Do 30 I=1,NLAng
       JAt=IAtomL(2,1,I)
       If(JAt.ge.1.and.JAt.le.NAtoms) then
        If(NPair(JAt).le.1) Keep(I)=.True.
       Else
        Keep(I)=.True.
       EndIf
   30 Continue
      NOld=NLAng
      If(CountKeep(NLAng,Keep).lt.NLAng) then
       Call PackGICBlock(MxAtP,MxTrm,NLAng,Keep,NTermL,IAtomL,
     $ IPrimL,ITVLA,IFixL,CoefL,ValTL)
       Write(IOut,'(''   LOCSVD primitive fallback true-linear '',
     $ ''filter: kept '',I5,'' of '',I5,'' linear coordinates.'')')
     $ NLAng,NOld
      EndIf
      Return
      End

*Deck SelectPruneMode
      Logical Function SelectPruneMode(ITPV,Mode)
      Integer ITPV,Mode
      Logical IsExoTorsion,IsButterflyTorsion,IsRingBend,IsRingTorsion
      SelectPruneMode=.False.
      If(Mode.eq.1) then
       SelectPruneMode=.True.
      ElseIf(Mode.eq.2) then
       If(IsExoTorsion(ITPV)) SelectPruneMode=.True.
      ElseIf(Mode.eq.3) then
       If(IsButterflyTorsion(ITPV)) SelectPruneMode=.True.
      ElseIf(Mode.eq.4) then
       If(.not.IsRingBend(ITPV)) SelectPruneMode=.True.
      ElseIf(Mode.eq.5) then
       If(IsRingBend(ITPV)) SelectPruneMode=.True.
      ElseIf(Mode.eq.6) then
       If(IsRingTorsion(ITPV)) SelectPruneMode=.True.
      ElseIf(Mode.eq.7) then
       SelectPruneMode=.True.
      EndIf
      Return
      End

*Deck IsRingBend
      Logical Function IsRingBend(ITPV)
      Integer ITPV
      IsRingBend=(ITPV.eq.14)
      Return
      End

*Deck IsRingTorsion
      Logical Function IsRingTorsion(ITPV)
      Integer ITPV
      IsRingTorsion=(ITPV.eq.1)
      Return
      End

*Deck IsButterflyTorsion
      Logical Function IsButterflyTorsion(ITPV)
      Integer ITPV
      IsButterflyTorsion=(ITPV.eq.2)
      Return
      End

*Deck IsExoTorsion
      Logical Function IsExoTorsion(ITPV)
      Integer ITPV
      IsExoTorsion=(ITPV.le.0)
      Return
      End

*Deck PruneSelectedBlockAgainst
      Subroutine PruneSelectedBlockAgainst(IOut,Label,NAtoms,NVar,IOff,
     $ BMat,ITPV,Mode,Keep,Scr,NBasis,NTarget,NKeep,NSelect)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NAtoms,NVar,IOff,Mode,NBasis,NTarget,NKeep,NSelect
      Integer ITPV(*)
      Character*(*) Label
      Dimension BMat(3*NAtoms,*),Scr(3*NAtoms,*)
      Logical Keep(*),SelectPruneMode
      Logical AnyRem,Selected
      Real*8 Norm,Norm0,Dot,TAbs,TRel
      Data TAbs/1.0D-10/, TRel/1.0D-08/

      NCart=3*NAtoms
      NKeep=0
      NSelect=0
      If(NVar.le.0) return

      Do 100 I=1,NVar
       Selected=SelectPruneMode(ITPV(I),Mode)
       If(.not.Selected) GoTo 100
       NSelect=NSelect+1
       If(NBasis.ge.NTarget) GoTo 100
       Row=IOff+I
       Norm0=0.0D0
       Do 110 K=1,NCart
        Scr(K,NBasis+1)=BMat(K,Row)
        Norm0=Norm0+Scr(K,NBasis+1)*Scr(K,NBasis+1)
  110  Continue
       Norm0=DSqrt(Norm0)
       If(Norm0.le.TAbs) GoTo 100

       Do 130 J=1,NBasis
        Dot=0.0D0
        Do 120 K=1,NCart
         Dot=Dot+Scr(K,NBasis+1)*Scr(K,J)
  120   Continue
        Do 125 K=1,NCart
         Scr(K,NBasis+1)=Scr(K,NBasis+1)-Dot*Scr(K,J)
  125   Continue
  130  Continue

       Norm=0.0D0
       Do 140 K=1,NCart
        Norm=Norm+Scr(K,NBasis+1)*Scr(K,NBasis+1)
  140  Continue
       Norm=DSqrt(Norm)
       If(Norm.gt.TAbs.and.Norm.gt.TRel*Norm0) then
        NKeep=NKeep+1
        NBasis=NBasis+1
        Keep(I)=.True.
        Do 150 K=1,NCart
         Scr(K,NBasis)=Scr(K,NBasis)/Norm
  150   Continue
       EndIf
  100 Continue

      If(NSelect.le.0) Return
      If(NKeep.lt.NSelect) then
       Write(IOut,'(''   '',A,'': kept '',I5,'' of '',I5,
     $ '' coordinates; removed '',I5,'' residual redundant.'')')
     $ Label,NKeep,NSelect,NSelect-NKeep
       Write(IOut,'(''     Removed local indices:'')')
       AnyRem=.False.
       Do 210 I=1,NVar
        If(SelectPruneMode(ITPV(I),Mode).and..not.Keep(I)) then
         Write(IOut,'(I6)',advance='no') I
         AnyRem=.True.
        EndIf
  210  Continue
       If(AnyRem) Write(IOut,'('' '')')
      Else
       Write(IOut,'(''   '',A,'': kept all '',I5,
     $ '' coordinates.'')') Label,NSelect
      EndIf
      Return
      End

*Deck SeedGICBasis
      Subroutine SeedGICBasis(IOut,NAtoms,NVar,IOff,BMat,Scr,NBasis)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NAtoms,NVar,IOff,NBasis
      Dimension BMat(3*NAtoms,*),Scr(3*NAtoms,*)
      Real*8 Norm,Norm0,Dot,TAbs,TRel
      Data TAbs/1.0D-10/, TRel/1.0D-08/

      NCart=3*NAtoms
      Do 100 I=1,NVar
       Row=IOff+I
       Norm0=0.0D0
       Do 10 K=1,NCart
        Scr(K,NBasis+1)=BMat(K,Row)
        Norm0=Norm0+Scr(K,NBasis+1)*Scr(K,NBasis+1)
   10  Continue
       Norm0=DSqrt(Norm0)
       If(Norm0.le.TAbs) go to 100
       Do 30 J=1,NBasis
        Dot=0.0D0
        Do 20 K=1,NCart
         Dot=Dot+Scr(K,NBasis+1)*Scr(K,J)
   20   Continue
        Do 25 K=1,NCart
         Scr(K,NBasis+1)=Scr(K,NBasis+1)-Dot*Scr(K,J)
   25   Continue
   30  Continue
       Norm=0.0D0
       Do 40 K=1,NCart
        Norm=Norm+Scr(K,NBasis+1)*Scr(K,NBasis+1)
   40  Continue
       Norm=DSqrt(Norm)
       If(Norm.gt.TAbs.and.Norm.gt.TRel*Norm0) then
        NBasis=NBasis+1
        Do 50 K=1,NCart
         Scr(K,NBasis)=Scr(K,NBasis)/Norm
   50   Continue
       EndIf
  100 Continue
      Return
      End

*Deck PruneOneBlockAgainst
      Subroutine PruneOneBlockAgainst(IOut,Label,NAtoms,NVar,IOff,BMat,
     $ Keep,Scr,NBasis,NTarget,NKeep)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NAtoms,NVar,IOff,NBasis,NTarget,NKeep
      Character*(*) Label
      Dimension BMat(3*NAtoms,*),Scr(3*NAtoms,*)
      Logical Keep(*)
      Logical AnyRem
      Real*8 Norm,Norm0,Dot,TAbs,TRel
      Data TAbs/1.0D-10/, TRel/1.0D-08/

      NCart=3*NAtoms
      NKeep=0
      If(NVar.le.0) return

      Do 10 I=1,NVar
       Keep(I)=.False.
   10 Continue

      Do 100 I=1,NVar
       If(NBasis.ge.NTarget) go to 100
       Row=IOff+I
       Norm0=0.0D0
       Do 110 K=1,NCart
        Scr(K,NBasis+1)=BMat(K,Row)
        Norm0=Norm0+Scr(K,NBasis+1)*Scr(K,NBasis+1)
  110  Continue
       Norm0=DSqrt(Norm0)
       If(Norm0.le.TAbs) go to 100

       Do 130 J=1,NBasis
        Dot=0.0D0
        Do 120 K=1,NCart
         Dot=Dot+Scr(K,NBasis+1)*Scr(K,J)
  120   Continue
        Do 125 K=1,NCart
         Scr(K,NBasis+1)=Scr(K,NBasis+1)-Dot*Scr(K,J)
  125   Continue
  130  Continue

       Norm=0.0D0
       Do 140 K=1,NCart
        Norm=Norm+Scr(K,NBasis+1)*Scr(K,NBasis+1)
  140  Continue
       Norm=DSqrt(Norm)
       If(Norm.gt.TAbs.and.Norm.gt.TRel*Norm0) then
        NKeep=NKeep+1
        NBasis=NBasis+1
        Keep(I)=.True.
        Do 150 K=1,NCart
         Scr(K,NBasis)=Scr(K,NBasis)/Norm
  150   Continue
       EndIf
  100 Continue

      If(NKeep.lt.NVar) then
       Write(IOut,'(''   '',A,'': kept '',I5,'' of '',I5,
     $ '' coordinates; removed '',I5,'' residual redundant.'')')
     $ Label,NKeep,NVar,NVar-NKeep
       Write(IOut,'(''     Removed local indices:'')')
       AnyRem=.False.
       Do 210 I=1,NVar
        If(.not.Keep(I)) then
         Write(IOut,'(I6)',advance='no') I
         AnyRem=.True.
        EndIf
  210  Continue
       If(AnyRem) Write(IOut,'('' '')')
      Else
       Write(IOut,'(''   '',A,'': kept all '',I5,
     $ '' coordinates.'')') Label,NVar
      EndIf
      Return
      End

*Deck WriteGICBMat
      Subroutine WriteGICBMat(NAtoms,NInt,BMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension BMat(3*NAtoms,*)
      NCart=3*NAtoms
      Open(77,File='bmat.out',Status='Unknown')
      Rewind(77)
      Write(77,'(A)') '# merlino.gicforge.bmatrix.v1'
      Write(77,'(A)') '# row col value; rows are final GICs'
      Write(77,'(2I8)') NInt,NCart
      Do 20 IInt=1,NInt
       Do 10 ICart=1,NCart
        Write(77,'(2I8,1X,D24.16)') IInt,ICart,BMat(ICart,IInt)
   10  Continue
   20 Continue
      Close(77)
      Return
      End

*Deck PruneOneBlock
      Subroutine PruneOneBlock(IOut,Label,NAtoms,NVar,IOff,BMat,Keep,
     $ Scr,NKeep)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NAtoms,NVar,IOff,NKeep
      Character*(*) Label
      Dimension BMat(3*NAtoms,*),Scr(3*NAtoms,*)
      Logical Keep(*)
      Logical AnyRem
      Real*8 Norm,Norm0,Dot,TAbs,TRel
      Data TAbs/1.0D-10/, TRel/1.0D-08/

      NCart=3*NAtoms
      NKeep=0
      If(NVar.le.0) return

      Do 10 I=1,NVar
       Keep(I)=.False.
   10 Continue

      Do 100 I=1,NVar
       Row=IOff+I
       Norm0=0.0D0
       Do 110 K=1,NCart
        Scr(K,NKeep+1)=BMat(K,Row)
        Norm0=Norm0+Scr(K,NKeep+1)*Scr(K,NKeep+1)
  110  Continue
       Norm0=DSqrt(Norm0)
       If(Norm0.le.TAbs) go to 100

       Do 130 J=1,NKeep
        Dot=0.0D0
        Do 120 K=1,NCart
         Dot=Dot+Scr(K,NKeep+1)*Scr(K,J)
  120   Continue
        Do 125 K=1,NCart
         Scr(K,NKeep+1)=Scr(K,NKeep+1)-Dot*Scr(K,J)
  125   Continue
  130  Continue

       Norm=0.0D0
       Do 140 K=1,NCart
        Norm=Norm+Scr(K,NKeep+1)*Scr(K,NKeep+1)
  140  Continue
       Norm=DSqrt(Norm)
       If(Norm.gt.TAbs.and.Norm.gt.TRel*Norm0) then
        NKeep=NKeep+1
        Keep(I)=.True.
        Do 150 K=1,NCart
         Scr(K,NKeep)=Scr(K,NKeep)/Norm
  150   Continue
       EndIf
  100 Continue

      If(NKeep.lt.NVar) then
       Write(IOut,'(''   '',A,'': kept '',I5,'' of '',I5,
     $ '' coordinates; removed '',I5,'' residual redundant.'')')
     $ Label,NKeep,NVar,NVar-NKeep
       Write(IOut,'(''     Removed local indices:'')')
       AnyRem=.False.
       Do 210 I=1,NVar
        If(.not.Keep(I)) then
         Write(IOut,'(I6)',advance='no') I
         AnyRem=.True.
        EndIf
  210  Continue
       If(AnyRem) Write(IOut,'('' '')')
      Else
       Write(IOut,'(''   '',A,'': kept all '',I5,
     $ '' coordinates.'')') Label,NVar
      EndIf
      Return
      End

*Deck PackGICBlock
      Subroutine PackGICBlock(MxAtP,MxTrm,NVar,Keep,NTerm,IAtom,IPrim,
     $ ITPV,IFixG,Coef,ValTot)
      Implicit Real*8 (A-H,O-Z)
      Integer MxAtP,MxTrm,NVar
      Integer I,J,T,A
      Dimension NTerm(*),IAtom(MxAtP,MxTrm,*),IPrim(MxTrm,*),ITPV(*)
      Dimension IFixG(*),Coef(MxTrm,*),ValTot(*)
      Logical Keep(*)

      J=0
      Do 100 I=1,NVar
       If(.not.Keep(I)) go to 100
       J=J+1
       If(J.ne.I) then
        NTerm(J)=NTerm(I)
        ITPV(J)=ITPV(I)
        IFixG(J)=IFixG(I)
        ValTot(J)=ValTot(I)
        Do 30 T=1,MxTrm
         IPrim(T,J)=IPrim(T,I)
         Coef(T,J)=Coef(T,I)
         Do 20 A=1,MxAtP
          IAtom(A,T,J)=IAtom(A,T,I)
   20    Continue
   30   Continue
       EndIf
  100 Continue
      NVar=J
      Return
      End
