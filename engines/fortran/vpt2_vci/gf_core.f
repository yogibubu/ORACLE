C     ORACLE independent Wilson-GF helper routines.
C     Fixed-form Fortran77.  The caller supplies non-redundant F and G.

      Subroutine M4GF(N,F,G,Eval,Work,Info)
      Integer N,Info
      Double Precision F(N,N),G(N,N),Eval(N),Work(*)
      Integer I,J,K,LD
      Double Precision Sum

      LD=N
      Info=0
C     Work layout: 1 G eigenvectors, 1+N*N G eigenvalues,
C     1+N*N+N G-half, 1+2*N*N+N transformed matrix.
      Call M4Copy(N,G,Work(1))
      Call M4Jacobi(N,Work(1),Work(1+N*N),Work(1+N*N+N),Info)
      If(Info.ne.0) Return
      Do 20 I=1,N
         If(Work(N*N+I).le.0.0D0) Then
            Info=2
            Return
         End If
20    Continue
      Call M4Zero(N,Work(1+N*N+N))
      Do 50 I=1,N
         Do 40 J=1,N
            Sum=0.0D0
            Do 30 K=1,N
               Sum=Sum+Work(I+(K-1)*LD)*Sqrt(Work(N*N+K))*
     $             Work(J+(K-1)*LD)
30          Continue
            Work(1+N*N+N+(J-1)*LD+I-1)=Sum
40       Continue
50    Continue
      Call M4Triple(N,Work(1+N*N+N),F,Work(1+2*N*N+N))
      Call M4Jacobi(N,Work(1+2*N*N+N),Eval,Work(1),Info)
      Return
      End

      Subroutine M4Triple(N,A,B,C)
      Integer N,I,J,K,L
      Double Precision A(N,N),B(N,N),C(N,N),Sum
      Do 40 I=1,N
         Do 30 J=1,N
            Sum=0.0D0
            Do 20 K=1,N
               Do 10 L=1,N
                  Sum=Sum+A(I,K)*B(K,L)*A(J,L)
10             Continue
20          Continue
            C(I,J)=Sum
30       Continue
40    Continue
      Return
      End

      Subroutine M4Copy(N,A,B)
      Integer N,I,J
      Double Precision A(N,N),B(*)
      Do 20 I=1,N
         Do 10 J=1,N
            B(I+(J-1)*N)=A(I,J)
10       Continue
20    Continue
      Return
      End

      Subroutine M4Zero(N,A)
      Integer N,I,J
      Double Precision A(N,N)
      Do 20 I=1,N
         Do 10 J=1,N
            A(I,J)=0.0D0
10       Continue
20    Continue
      Return
      End

      Subroutine M4Jacobi(N,A,Eval,Evec,Info)
      Integer N,Info
      Double Precision A(N,N),Eval(N),Evec(N,N)
      Integer I,J,P,Q,Iter,MaxIt
      Double Precision App,Aqq,Apq,Phi,C,S,Aip,Aiq,Tol,Off
      Info=0
      MaxIt=100*N*N
      Do 20 I=1,N
         Do 10 J=1,N
            Evec(I,J)=0.0D0
10       Continue
         Evec(I,I)=1.0D0
20    Continue
      Do 100 Iter=1,MaxIt
         Off=0.0D0
         P=1
         Q=1
         Do 40 I=1,N-1
            Do 30 J=I+1,N
               If(Abs(A(I,J)).gt.Off) Then
                  Off=Abs(A(I,J))
                  P=I
                  Q=J
               End If
30          Continue
40       Continue
         Tol=1.0D-12
         If(Off.lt.Tol) Go To 120
         App=A(P,P)
         Aqq=A(Q,Q)
         Apq=A(P,Q)
         Phi=0.5D0*Atan2(2.0D0*Apq,Aqq-App)
         C=Cos(Phi)
         S=Sin(Phi)
         Do 50 I=1,N
            Aip=A(I,P)
            Aiq=A(I,Q)
            A(I,P)=C*Aip-S*Aiq
            A(I,Q)=S*Aip+C*Aiq
50       Continue
         Do 60 J=1,N
            Aip=A(P,J)
            Aiq=A(Q,J)
            A(P,J)=C*Aip-S*Aiq
            A(Q,J)=S*Aip+C*Aiq
60       Continue
         A(P,Q)=0.0D0
         A(Q,P)=0.0D0
         Do 70 I=1,N
            Aip=Evec(I,P)
            Aiq=Evec(I,Q)
            Evec(I,P)=C*Aip-S*Aiq
            Evec(I,Q)=S*Aip+C*Aiq
70       Continue
100   Continue
      Info=1
      Return
120   Continue
      Do 130 I=1,N
         Eval(I)=A(I,I)
130   Continue
      Return
      End
