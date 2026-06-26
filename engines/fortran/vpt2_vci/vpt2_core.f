C     ORACLE independent small-space VPT2 helper routines.
C     Fixed-form Fortran77. Uses the same dimensionless normal-coordinate
C     convention as the Python VPT2/VCI core.

      Double Precision Function M4XElem(NL,NR,IPow)
      Integer NL,NR,IPow,I,J,K
      Double Precision V(0:128),W(0:128),X
      Do 10 I=0,128
         V(I)=0.0D0
         W(I)=0.0D0
10    Continue
      If(NR.lt.0 .or. NR.gt.128 .or. NL.lt.0 .or. NL.gt.128) Then
         M4XElem=0.0D0
         Return
      End If
      V(NR)=1.0D0
      If(IPow.eq.0) Then
         M4XElem=V(NL)
         Return
      End If
      Do 50 K=1,IPow
         Do 20 I=0,128
            W(I)=0.0D0
20       Continue
         Do 30 I=0,127
            If(Abs(V(I)).le.1.0D-30) Go To 30
            X=Sqrt(Dble(I+1))/Sqrt(2.0D0)
            W(I+1)=W(I+1)+X*V(I)
30       Continue
         Do 40 I=1,128
            If(Abs(V(I)).le.1.0D-30) Go To 40
            X=Sqrt(Dble(I))/Sqrt(2.0D0)
            W(I-1)=W(I-1)+X*V(I)
40       Continue
         Do 45 J=0,128
            V(J)=W(J)
45       Continue
50    Continue
      M4XElem=V(NL)
      Return
      End

      Double Precision Function M4TermEl(NMode,Left,Right,NIdx,Idx)
      Integer NMode,Left(NMode),Right(NMode),NIdx,Idx(NIdx)
      Integer Pow(64),I,Mode
      Double Precision Val,M4XElem
      If(NMode.gt.64) Then
         M4TermEl=0.0D0
         Return
      End If
      Do 10 I=1,NMode
         Pow(I)=0
10    Continue
      Do 20 I=1,NIdx
         Mode=Idx(I)
         If(Mode.lt.1 .or. Mode.gt.NMode) Then
            M4TermEl=0.0D0
            Return
         End If
         Pow(Mode)=Pow(Mode)+1
20    Continue
      Val=1.0D0
      Do 30 I=1,NMode
         Val=Val*M4XElem(Left(I),Right(I),Pow(I))
         If(Abs(Val).le.1.0D-30) Go To 40
30    Continue
40    Continue
      M4TermEl=Val
      Return
      End

      Double Precision Function M4HarmE(NMode,State,Freq)
      Integer NMode,State(NMode),I
      Double Precision Freq(NMode),E
      E=0.0D0
      Do 10 I=1,NMode
         E=E+Freq(I)*(Dble(State(I))+0.5D0)
10    Continue
      M4HarmE=E
      Return
      End

      Subroutine M4VPT2Basis(NMode,NState,Basis,Freq,NC3,C3Idx,C3Val,
     $                       NC4,C4Idx,C4Val,E0,E1,E2,ETot)
      Integer NMode,NState,NC3,NC4
      Integer Basis(NState,NMode),C3Idx(3,NC3),C4Idx(4,NC4)
      Double Precision Freq(NMode),C3Val(NC3),C4Val(NC4)
      Double Precision E0(NState),E1(NState),E2(NState),ETot(NState)
      Integer I,J,K,M,Left(64),Right(64)
      Double Precision Den,Coup,El,M4TermEl,M4HarmE
      Do 100 I=1,NState
         Do 10 M=1,NMode
            Left(M)=Basis(I,M)
10       Continue
         E0(I)=M4HarmE(NMode,Left,Freq)
         E1(I)=0.0D0
         E2(I)=0.0D0
         Do 20 K=1,NC4
            El=M4TermEl(NMode,Left,Left,4,C4Idx(1,K))
            E1(I)=E1(I)+C4Val(K)*El
20       Continue
         Do 60 J=1,NState
            If(J.eq.I) Go To 60
            Do 30 M=1,NMode
               Right(M)=Basis(J,M)
30          Continue
            Coup=0.0D0
            Do 40 K=1,NC3
               El=M4TermEl(NMode,Right,Left,3,C3Idx(1,K))
               Coup=Coup+C3Val(K)*El
40          Continue
            If(Abs(Coup).le.1.0D-30) Go To 60
            Den=E0(I)-M4HarmE(NMode,Right,Freq)
            If(Abs(Den).gt.1.0D-8) E2(I)=E2(I)+Coup*Coup/Den
60       Continue
         ETot(I)=E0(I)+E1(I)+E2(I)
100   Continue
      Return
      End
