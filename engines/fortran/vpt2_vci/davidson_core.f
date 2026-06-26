C     ORACLE independent Davidson helper routines.
C     These routines know only vectors, diagonals and residuals.

      Subroutine M4DavCorr(N,Theta,Diag,Resid,Corr,Info)
      Integer N,Info,I
      Double Precision Theta,Diag(N),Resid(N),Corr(N),Den
      Info=0
      Do 10 I=1,N
         Den=Theta-Diag(I)
         If(Abs(Den).lt.1.0D-10) Then
            If(Den.lt.0.0D0) Then
               Den=-1.0D-10
            Else
               Den=1.0D-10
            End If
         End If
         Corr(I)=Resid(I)/Den
10    Continue
      Return
      End

      Double Precision Function M4Dnrm2(N,X)
      Integer N,I
      Double Precision X(N),Sum
      Sum=0.0D0
      Do 10 I=1,N
         Sum=Sum+X(I)*X(I)
10    Continue
      M4Dnrm2=Sqrt(Sum)
      Return
      End

      Subroutine M4Orth(N,M,Q,V,Norm)
      Integer N,M,I,J
      Double Precision Q(N,M),V(N),Norm,Dot
      Double Precision M4Dnrm2
      If(M.gt.0) Then
         Do 30 J=1,M
            Dot=0.0D0
            Do 10 I=1,N
               Dot=Dot+Q(I,J)*V(I)
10          Continue
            Do 20 I=1,N
               V(I)=V(I)-Dot*Q(I,J)
20          Continue
30       Continue
      End If
      Norm=M4Dnrm2(N,V)
      If(Norm.gt.1.0D-12) Then
         Do 40 I=1,N
            V(I)=V(I)/Norm
40       Continue
      End If
      Return
      End
