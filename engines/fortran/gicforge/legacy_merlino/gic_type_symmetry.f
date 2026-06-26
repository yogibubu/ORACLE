*Deck GICTypeSym
C=======================================================================
C  GIC TYPE-LOCAL SYMMETRIZATION
C
C  Generic symmetry combinations for GICForge.
C
C  Scope:
C   - operates only inside one coordinate family at a time
C   - never mixes stretches with bends, torsions, linear bends, or OOPs
C   - acts only on generic one-term coordinates (ITPV = 0)
C   - leaves specialized ring/puckering/butterfly coordinates untouched
C
C  The driver currently calls this for bends, linear bends and torsions.
C  Stretchings are kept as primitive R(i,j) coordinates.  OOP signatures
C  are implemented here for completeness, but the GICForge driver keeps
C  OOP primitives unchanged until the Gaussian writer supports OOP
C  linear-combination output.
C
C  This mirrors the Python policy: group comparable primitives first, build
C  sum/difference combinations inside each homogeneous block, then let the
C  B-rank pruning remove any remaining dependencies.
C=======================================================================

      Subroutine SymOneGICBlock(IOut,Label,MaxAtG,MaxTer,Itp,NVar,
     $ NTerm,IAtom,ITPV,IFixG,IAn,Coef)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,MaxAtG,MaxTer,Itp,NVar
      Character*(*) Label
      Dimension NTerm(*),IAtom(MaxAtG,MaxTer,*),ITPV(*),IFixG(*)
      Dimension IAn(*),Coef(MaxTer,*)
      Logical Used(1000),SameGICSig,Skip
      Integer Group(100)
      Integer ONTerm(100)
      Integer I,J,NG,IG,IX,T,A
      Integer OAtom(4,15,100)
      Real*8 OCoef(15,100)

      If(NVar.le.1) return
      Do 5 I=1,NVar
       Used(I)=.False.
    5 Continue

      Do 100 I=1,NVar
       If(Used(I)) go to 100
       Used(I)=.True.
       NG=1
       Group(1)=I

C      Do not combine coordinates that are already semantic GICs.
C      ITPV != 0 labels ring/puckering/butterfly/etc. coordinates.
       If(ITPV(I).ne.0.or.NTerm(I).ne.1) go to 100

       Do 20 J=I+1,NVar
        If(Used(J)) go to 20
        If(ITPV(J).ne.0.or.NTerm(J).ne.1) go to 20
C       Frozen and active coordinates must not enter the same symmetry block.
        If(IFixG(J).ne.IFixG(I)) go to 20
        If(SameGICSig(MaxAtG,MaxTer,Itp,I,J,IAtom,IAn)) then
         NG=NG+1
         If(NG.le.100) Group(NG)=J
         Used(J)=.True.
        EndIf
   20  Continue

       If(NG.le.1) go to 100
C      A symmetric sum over NG primitives needs NG terms in one GIC.
       If(NG.gt.MaxTer) then
        Write(IOut,'(''   '',A,'': same-type symmetry group of size '',
     $ I5,'' skipped; MaxTer too small.'')') Label,NG
        go to 100
       EndIf

       Skip=.False.
       Do 30 IG=1,NG
        IX=Group(IG)
        ONTerm(IG)=NTerm(IX)
        If(ONTerm(IG).gt.MaxTer) Skip=.True.
        Do 32 T=1,ONTerm(IG)
         OCoef(T,IG)=Coef(T,IX)
         Do 31 A=1,MaxAtG
          OAtom(A,T,IG)=IAtom(A,T,IX)
   31    Continue
   32   Continue
   30  Continue
       If(Skip) go to 100

C      First coordinate: totally symmetric normalized sum.
       IX=Group(1)
       NTerm(IX)=NG
       Den=DSqrt(DBLE(NG))
       Do 40 IG=1,NG
        Coef(IG,IX)=OCoef(1,IG)/Den
        Do 41 A=1,MaxAtG
         IAtom(A,IG,IX)=OAtom(A,1,IG)
   41   Continue
   40  Continue

C      Remaining coordinates: local orthonormal differences.  OrdRed and the
C      final B-rank pruning decide which of these differences are independent.
       Do 60 IG=2,NG
        IX=Group(IG)
        NTerm(IX)=2
        Coef(1,IX)=OCoef(1,IG-1)/DSqrt(2.0D0)
        Coef(2,IX)=-OCoef(1,IG)/DSqrt(2.0D0)
        Do 50 A=1,MaxAtG
         IAtom(A,1,IX)=OAtom(A,1,IG-1)
         IAtom(A,2,IX)=OAtom(A,1,IG)
   50   Continue
   60  Continue

       Write(IOut,'(''   '',A,'': symmetrized same-type group starting'',
     $ '' at '',I5,'' size '',I5)') Label,I,NG
  100 Continue
      Return
      End

*Deck SameGICSig
      Logical Function SameGICSig(MaxAtG,MaxTer,Itp,I,J,IAtom,IAn)
      Implicit Real*8 (A-H,O-Z)
      Integer MaxAtG,MaxTer,Itp,I,J
      Dimension IAtom(MaxAtG,MaxTer,*),IAn(*)
      Integer AI1,AI2,AI3,AI4,AJ1,AJ2,AJ3,AJ4
      Logical SamePairZ,SameTripleZ
      SameGICSig=.False.

      AI1=IAtom(1,1,I)
      AI2=IAtom(2,1,I)
      AI3=IAtom(3,1,I)
      AI4=IAtom(4,1,I)
      AJ1=IAtom(1,1,J)
      AJ2=IAtom(2,1,J)
      AJ3=IAtom(3,1,J)
      AJ4=IAtom(4,1,J)

      If(Itp.eq.1) then
       If(SamePairZ(IAn(AI1),IAn(AI2),IAn(AJ1),IAn(AJ2)))
     $  SameGICSig=.True.
      ElseIf(Itp.eq.2.or.Itp.eq.3) then
       If(IAn(AI2).ne.IAn(AJ2)) return
       If(SamePairZ(IAn(AI1),IAn(AI3),IAn(AJ1),IAn(AJ3)))
     $  SameGICSig=.True.
      ElseIf(Itp.eq.4) then
       If(IAn(AI2).ne.IAn(AJ2)) return
       If(IAn(AI3).ne.IAn(AJ3)) return
       If(IAn(AI1).ne.IAn(AJ1)) return
       If(IAn(AI4).ne.IAn(AJ4)) return
       SameGICSig=.True.
      ElseIf(Itp.eq.5) then
       If(IAn(AI2).ne.IAn(AJ2)) return
       If(SameTripleZ(IAn(AI1),IAn(AI3),IAn(AI4),
     $  IAn(AJ1),IAn(AJ3),IAn(AJ4))) SameGICSig=.True.
      EndIf
      Return
      End

*Deck SamePairZ
      Logical Function SamePairZ(A,B,C,D)
      Implicit None
      Integer A,B,C,D
      SamePairZ=.False.
      If(A.eq.C.and.B.eq.D) SamePairZ=.True.
      If(A.eq.D.and.B.eq.C) SamePairZ=.True.
      Return
      End

*Deck SameTripleZ
      Logical Function SameTripleZ(A,B,C,D,E,F)
      Implicit None
      Integer A,B,C,D,E,F
      Integer X1(3),X2(3),I,J,T
      X1(1)=A
      X1(2)=B
      X1(3)=C
      X2(1)=D
      X2(2)=E
      X2(3)=F
      Do 10 I=1,2
       Do 20 J=I+1,3
        If(X1(J).lt.X1(I)) then
         T=X1(I)
         X1(I)=X1(J)
         X1(J)=T
        EndIf
        If(X2(J).lt.X2(I)) then
         T=X2(I)
         X2(I)=X2(J)
         X2(J)=T
        EndIf
   20  Continue
   10 Continue
      SameTripleZ=.False.
      If(X1(1).eq.X2(1).and.X1(2).eq.X2(2).and.
     $ X1(3).eq.X2(3)) SameTripleZ=.True.
      Return
      End
