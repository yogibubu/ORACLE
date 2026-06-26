C     ORACLE independent dense small-space VCI helpers.
C     Fixed-form Fortran77. Large spaces will use a Davidson driver later.

      Subroutine M4VCIBasis(NMode,MaxQ,MaxSt,NState,Basis,Info)
      Integer NMode,MaxQ,MaxSt,NState,Basis(MaxSt,NMode),Info
      Integer QMin(64),QMax(64),ClsMin(4),ClsMax(4),I
      If(NMode.gt.64) Then
         Info=2
         Return
      End If
      Do 10 I=1,NMode
         QMin(I)=0
         QMax(I)=MaxQ
10    Continue
      Do 20 I=1,4
         ClsMin(I)=-1
         ClsMax(I)=-1
20    Continue
      Call M4VCIBasisCtl(NMode,MaxQ,MaxSt,QMin,QMax,ClsMin,ClsMax,
     $                   NState,Basis,Info)
      Return
      End

      Subroutine M4VCIBasisCtl(NMode,MaxQ,MaxSt,QMin,QMax,ClsMin,
     $                         ClsMax,NState,Basis,Info)
      Integer NMode,MaxQ,MaxSt,NState,Basis(MaxSt,NMode),Info
      Integer QMin(NMode),QMax(NMode),ClsMin(4),ClsMax(4)
      Integer State(64),I,K,Sum,Done,OK,M4BasisAccept
      If(NMode.gt.64) Then
         Info=2
         Return
      End If
      NState=0
      Info=0
      Do 30 I=1,NMode
         If(QMin(I).lt.0 .or. QMax(I).lt.QMin(I)) Then
            Info=3
            Return
         End If
         State(I)=QMin(I)
30    Continue
      Done=0
40    Continue
      Sum=0
      Do 50 I=1,NMode
         Sum=Sum+State(I)
50    Continue
      OK=M4BasisAccept(NMode,State,MaxQ,ClsMin,ClsMax)
      If(OK.eq.1) Then
         NState=NState+1
         If(NState.gt.MaxSt) Then
            Info=1
            Return
         End If
         Do 60 I=1,NMode
            Basis(NState,I)=State(I)
60       Continue
      End If
      K=NMode
70    Continue
      If(K.le.0) Then
         Done=1
      Else If(State(K).lt.QMax(K)) Then
         State(K)=State(K)+1
         Do 80 I=K+1,NMode
            State(I)=QMin(I)
80       Continue
      Else
         K=K-1
         Go To 70
      End If
      If(Done.eq.0) Go To 40
      Return
      End

      Integer Function M4BasisAccept(NMode,State,MaxQ,ClsMin,ClsMax)
      Integer NMode,State(NMode),MaxQ,ClsMin(4),ClsMax(4)
      Integer I,Tot,NExc
      Tot=0
      NExc=0
      Do 10 I=1,NMode
         Tot=Tot+State(I)
         If(State(I).gt.0) NExc=NExc+1
10    Continue
      If(Tot.gt.MaxQ) Then
         M4BasisAccept=0
         Return
      End If
      If(NExc.ge.1 .and. NExc.le.4) Then
         If(ClsMin(NExc).ge.0 .and. Tot.lt.ClsMin(NExc)) Then
            M4BasisAccept=0
            Return
         End If
         If(ClsMax(NExc).ge.0 .and. Tot.gt.ClsMax(NExc)) Then
            M4BasisAccept=0
            Return
         End If
      End If
      M4BasisAccept=1
      Return
      End

      Subroutine M4VCIHarm(NMode,NState,Basis,Freq,H)
      Integer NMode,NState,Basis(NState,NMode)
      Double Precision Freq(NMode),H(NState,NState),E
      Integer I,J,K
      Do 20 I=1,NState
         Do 10 J=1,NState
            H(I,J)=0.0D0
10       Continue
20    Continue
      Do 40 I=1,NState
         E=0.0D0
         Do 30 K=1,NMode
            E=E+Freq(K)*(Dble(Basis(I,K))+0.5D0)
30       Continue
         H(I,I)=E
40    Continue
      Return
      End

      Double Precision Function M4X4Diag(N)
      Integer N
C     <n|x**4|n> for x=(a+a+)/sqrt(2)
      M4X4Diag=0.75D0*(2.0D0*Dble(N)*Dble(N)+
     $         2.0D0*Dble(N)+1.0D0)
      Return
      End
