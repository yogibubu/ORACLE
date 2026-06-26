*Deck LocSVDJacobi
C=======================================================================
C  Local SVD helper for GICForge.
C
C  The caller supplies a small symmetric Gram matrix G = B B^T in A.
C  Jacobi diagonalization returns eigenvalues in descending order and
C  eigenvectors in EVec(:,mode).  These eigenvectors are the local GNIC
C  coefficients for the primitive rows used to build G.
C
C  This file is intentionally self-contained and does not require LAPACK.
C=======================================================================
      Subroutine LocSVDJacobi(NMax,N,A,EVal,EVec,Rank)
      Implicit Real*8 (A-H,O-Z)
      Integer NMax,N,Rank
      Dimension A(NMax,*),EVal(*),EVec(NMax,*)
      Real*8 TAbs,TRel
      Data TAbs/1.0D-10/,TRel/1.0D-08/

      If(N.le.0) then
       Rank=0
       Return
      EndIf

      Do 20 I=1,N
       Do 10 J=1,N
        EVec(I,J)=0.0D0
   10  Continue
       EVec(I,I)=1.0D0
   20 Continue

      Do 200 ISweep=1,100
       P=1
       Q=1
       AMx=0.0D0
       Do 40 I=1,N-1
        Do 30 J=I+1,N
         If(DAbs(A(I,J)).gt.AMx) then
          AMx=DAbs(A(I,J))
          P=I
          Q=J
         EndIf
   30   Continue
   40  Continue
       If(AMx.le.1.0D-12) go to 210

       App=A(P,P)
       Aqq=A(Q,Q)
       Apq=A(P,Q)
       Tau=(Aqq-App)/(2.0D0*Apq)
       If(Tau.ge.0.0D0) then
        T=1.0D0/(Tau+DSqrt(1.0D0+Tau*Tau))
       Else
        T=-1.0D0/(-Tau+DSqrt(1.0D0+Tau*Tau))
       EndIf
       C=1.0D0/DSqrt(1.0D0+T*T)
       S=T*C

       A(P,P)=App-T*Apq
       A(Q,Q)=Aqq+T*Apq
       A(P,Q)=0.0D0
       A(Q,P)=0.0D0

       Do 60 K=1,N
        If(K.ne.P.and.K.ne.Q) then
         Akp=A(K,P)
         Akq=A(K,Q)
         A(K,P)=C*Akp-S*Akq
         A(P,K)=A(K,P)
         A(K,Q)=S*Akp+C*Akq
         A(Q,K)=A(K,Q)
        EndIf
   60  Continue

       Do 70 K=1,N
        Vkp=EVec(K,P)
        Vkq=EVec(K,Q)
        EVec(K,P)=C*Vkp-S*Vkq
        EVec(K,Q)=S*Vkp+C*Vkq
   70  Continue
  200 Continue

  210 Continue
      Do 220 I=1,N
       EVal(I)=A(I,I)
       If(EVal(I).lt.0.0D0.and.DAbs(EVal(I)).lt.1.0D-12)
     $  EVal(I)=0.0D0
  220 Continue

      Do 260 I=1,N-1
       IBest=I
       VBest=EVal(I)
       Do 230 J=I+1,N
        If(EVal(J).gt.VBest) then
         IBest=J
         VBest=EVal(J)
        EndIf
  230  Continue
       If(IBest.ne.I) then
        Tmp=EVal(I)
        EVal(I)=EVal(IBest)
        EVal(IBest)=Tmp
        Do 240 K=1,N
         Tmp=EVec(K,I)
         EVec(K,I)=EVec(K,IBest)
         EVec(K,IBest)=Tmp
  240   Continue
       EndIf
  260 Continue

      Rank=0
      Ref=EVal(1)
      If(Ref.lt.0.0D0) Ref=0.0D0
      Do 280 I=1,N
       If(EVal(I).gt.TAbs*TAbs.and.EVal(I).gt.TRel*TRel*Ref) then
        Rank=Rank+1
        IMax=1
        CMax=DAbs(EVec(1,I))
        Do 270 K=2,N
         If(DAbs(EVec(K,I)).gt.CMax) then
          IMax=K
          CMax=DAbs(EVec(K,I))
         EndIf
  270   Continue
        If(EVec(IMax,I).lt.0.0D0) then
         Do 275 K=1,N
          EVec(K,I)=-EVec(K,I)
  275    Continue
        EndIf
       EndIf
  280 Continue

      Return
      End
