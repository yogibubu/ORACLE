*Deck PathDVR
C Fortran77 DVR kernels for ORACLE.
C Python prepares grids from Gaussian; this program solves Hamiltonians.
C
C dvrin format:
C MODE
C
C MODE=1  One-dimensional grid DVR
C   NPTS NLEVELS IBOUND
C   q_i V_i                     i=1,NPTS
C   IBOUND=0 sinc/nonperiodic, IBOUND=1 Fourier/periodic
C
C MODE=2  One-dimensional Gaussian basis
C   NGRID NBASIS NLEVELS
C   q_i V_i                     i=1,NGRID
C   Gaussian centers and widths are generated from the interval.
C
C MODE=3  Two-dimensional product-grid DVR
C   N1 N2 NLEVELS IBOUND1 IBOUND2 G11 G22 G12
C   q1_i                        i=1,N1
C   q2_j                        j=1,N2
C   V_ij                        i=1,N1, j=1,N2
      Program PathDVR
      Implicit Real*8 (A-H,O-Z)
      Parameter (MaxPts=2000,MaxGrid=2000,MaxBas=400)
      Parameter (MxTri=MaxPts*(MaxPts+1)/2)
      Parameter (MxVec=MaxPts*MaxPts)
      Dimension Q(MaxGrid),Vpot(MaxGrid),Q2(MaxGrid),V2D(MaxGrid)
      Dimension Ham(MxTri),Eval(MaxPts),EVec(MxVec)
      Dimension WA(6*MaxPts),Scr(MxTri+MaxPts),IScr(MaxPts)
      Common /DVRIO/ IOut
      IOut=6
      Open(10,File='dvrin',Status='old')
      Read(10,*) Mode
      If(Mode.eq.1) then
       Call RunGrid1D(10,IOut,MaxPts,MxTri,MaxGrid,MxVec,Q,Vpot,Ham,
     $  Eval,EVec,WA,Scr,IScr)
      ElseIf(Mode.eq.2) then
       Call RunGauss1D(10,IOut,MaxPts,MaxGrid,MaxBas,MxTri,MxVec,Q,
     $  Vpot,Ham,Eval,EVec,WA,Scr,IScr)
      ElseIf(Mode.eq.3) then
       Call RunGrid2D(10,IOut,MaxPts,MaxGrid,MxTri,MxVec,Q,Q2,V2D,
     $  Ham,Eval,EVec,WA,Scr,IScr)
      Else
       Write(IOut,'(''Unsupported DVR mode:'',I8)') Mode
       Stop
      EndIf
      Close(10)
      End

*Deck RunGrid1D
      Subroutine RunGrid1D(In,IOut,MaxPts,MxTri,MaxGrid,MxVec,Q,Vpot,
     $ Ham,Eval,EVec,WA,Scr,IScr)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q(MaxGrid),Vpot(MaxGrid),Ham(MxTri),Eval(MaxPts)
      Dimension EVec(MxVec),WA(6,*),Scr(*),IScr(*)
      Logical Schmid
      Schmid=.True.
      Read(In,*) NPts,NLevels,IBound
      If(NPts.lt.2.or.NPts.gt.MaxPts) Stop 'Invalid NPTS'
      If(NLevels.lt.1) NLevels=1
      If(NLevels.gt.NPts) NLevels=NPts
      Do 10 I=1,NPts
       Read(In,*) Q(I),Vpot(I)
   10 Continue
      Call CheckGrid(IOut,NPts,Q,DX)
      Call BuildHam1D(IOut,NPts,IBound,Q,Vpot,DX,Ham)
      MDV=MxTri+MaxPts
      Call DVRHQRII(IOut,NPts,1,NLevels,-1,Ham,Eval,NPts,EVec,Schmid,
     $ IErr,IScr,WA,Scr,MDV)
      If(IErr.ne.0) Stop 'DVRHQRII failed'
      Call WriteReport('grid-1d',NPts,NLevels,DX,Eval)
      Call WriteLevels(NLevels,Eval)
      Call WriteVectors(NPts,NLevels,Q,Vpot,EVec)
      Return
      End

*Deck RunGauss1D
      Subroutine RunGauss1D(In,IOut,MaxPts,MaxGrid,MaxBas,MxTri,MxVec,Q,
     $ Vpot,Ham,Eval,EVec,WA,Scr,IScr)
      Implicit Real*8 (A-H,O-Z)
      Parameter (MxBas=400)
      Dimension Q(MaxGrid),Vpot(MaxGrid),Ham(MxTri),Eval(MaxPts)
      Dimension EVec(MxVec),WA(6,*),Scr(*),IScr(*)
      Dimension S(MxBas*(MxBas+1)/2),SEval(MxBas),SVec(MxBas*MxBas)
      Dimension HOr(MxBas*(MxBas+1)/2),Phi(MxBas)
      Dimension Cen(MxBas),Alp(MxBas),XMat(MxBas*MxBas)
      Logical Schmid
      Schmid=.True.
      Read(In,*) NGrid,NBas,NLevels
      If(NGrid.lt.2.or.NGrid.gt.MaxGrid) Stop 'Invalid NGRID'
      If(NBas.lt.2.or.NBas.gt.MaxBas) Stop 'Invalid NBAS'
      If(NLevels.lt.1) NLevels=1
      If(NLevels.gt.NBas) NLevels=NBas
      Do 10 I=1,NGrid
       Read(In,*) Q(I),Vpot(I)
   10 Continue
      Call CheckGrid(IOut,NGrid,Q,DX)
      Call InitGauss(NBas,Q(1),Q(NGrid),Cen,Alp)
      Call BuildGaussSH(NGrid,NBas,Q,Vpot,DX,Cen,Alp,S,Ham)
      MDV=MxTri+MaxPts
      Call DVRHQRII(IOut,NBas,1,NBas,-1,S,SEval,NBas,SVec,Schmid,
     $ IErr,IScr,WA,Scr,MDV)
      If(IErr.ne.0) Stop 'Overlap diagonalization failed'
      Call BuildOrtho(NBas,SEval,SVec,XMat)
      Call TransformHam(NBas,Ham,XMat,HOr)
      Call DVRHQRII(IOut,NBas,1,NLevels,-1,HOr,Eval,NBas,EVec,Schmid,
     $ IErr,IScr,WA,Scr,MDV)
      If(IErr.ne.0) Stop 'Gaussian Hamiltonian failed'
      Call WriteReport('gaussian-1d',NGrid,NLevels,DX,Eval)
      Call WriteLevels(NLevels,Eval)
      Call WriteGaussVectors(NGrid,NBas,NLevels,Q,Vpot,Cen,Alp,XMat,
     $ EVec,Phi)
      Return
      End

*Deck RunGrid2D
      Subroutine RunGrid2D(In,IOut,MaxPts,MaxGrid,MxTri,MxVec,Q1,Q2,V2D,
     $ Ham,Eval,EVec,WA,Scr,IScr)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q1(MaxGrid),Q2(MaxGrid),V2D(MaxGrid),Ham(MxTri)
      Dimension Eval(MaxPts),EVec(MxVec),WA(6,*),Scr(*),IScr(*)
      Logical Schmid
      Schmid=.True.
      Read(In,*) N1,N2,NLevels,IB1,IB2,G11,G22,G12
      NTot=N1*N2
      If(N1.lt.2.or.N2.lt.2.or.NTot.gt.MaxPts) Stop 'Invalid 2D grid'
      If(NLevels.lt.1) NLevels=1
      If(NLevels.gt.NTot) NLevels=NTot
      Do 10 I=1,N1
       Read(In,*) Q1(I)
   10 Continue
      Do 20 J=1,N2
       Read(In,*) Q2(J)
   20 Continue
      Do 40 I=1,N1
       Do 30 J=1,N2
        K=(I-1)*N2+J
        Read(In,*) V2D(K)
   30  Continue
   40 Continue
      Call CheckGrid(IOut,N1,Q1,DX1)
      Call CheckGrid(IOut,N2,Q2,DX2)
      Call BuildHam2D(IOut,N1,N2,IB1,IB2,Q1,Q2,DX1,DX2,V2D,G11,G22,
     $ G12,Ham)
      MDV=MxTri+MaxPts
      Call DVRHQRII(IOut,NTot,1,NLevels,-1,Ham,Eval,NTot,EVec,Schmid,
     $ IErr,IScr,WA,Scr,MDV)
      If(IErr.ne.0) Stop '2D Hamiltonian failed'
      Call WriteReport('grid-2d',NTot,NLevels,DX1,Eval)
      Call WriteLevels(NLevels,Eval)
      Call WriteVectors2D(N1,N2,NLevels,Q1,Q2,V2D,EVec)
      Return
      End

*Deck CheckGrid
      Subroutine CheckGrid(IOut,NPts,Q,DX)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q(*)
      DX=Q(2)-Q(1)
      If(DX.eq.0.0D0) Stop 'DVR grid has zero spacing'
      Tol=1.0D-6*DAbs(DX)
      Do 10 I=2,NPts-1
       If(DAbs((Q(I+1)-Q(I))-DX).gt.Tol) Stop 'Unequal DVR grid'
   10 Continue
      Return
      End

*Deck BuildHam1D
      Subroutine BuildHam1D(IOut,NPts,IBound,Q,Vpot,DX,Ham)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q(*),Vpot(*),Ham(*)
      Parameter (H2CM=219474.6313705D0)
      Real*8 Length
      Pi=4.0D0*DAtan(1.0D0)
      If(IBound.eq.0) then
       Call BuildSinc(NPts,Vpot,DX,Ham,H2CM,Pi)
      ElseIf(IBound.eq.1) then
       Length=Q(NPts)+DX-Q(1)
       Call BuildFourier(NPts,Vpot,Length,Ham,H2CM,Pi)
      Else
       Stop 'Unsupported boundary'
      EndIf
      Return
      End

*Deck BuildSinc
      Subroutine BuildSinc(NPts,Vpot,DX,Ham,H2CM,Pi)
      Implicit Real*8 (A-H,O-Z)
      Dimension Vpot(*),Ham(*)
      IJ=0
      Do 20 I=1,NPts
       Do 10 J=1,I
        IJ=IJ+1
        If(I.eq.J) then
         Ham(IJ)=H2CM*Pi*Pi/(6.0D0*DX*DX)+Vpot(I)
        Else
         Idiff=I-J
         Sign=1.0D0
         If(Mod(Idiff,2).ne.0) Sign=-1.0D0
         Ham(IJ)=H2CM*Sign/(DX*DX*DFloat(Idiff*Idiff))
        EndIf
   10  Continue
   20 Continue
      Return
      End

*Deck BuildFourier
      Subroutine BuildFourier(NPts,Vpot,Length,Ham,H2CM,Pi)
      Implicit Real*8 (A-H,O-Z)
      Real*8 Length
      Dimension Vpot(*),Ham(*)
      IJ=0
      Do 30 I=1,NPts
       Do 20 J=1,I
        IJ=IJ+1
        Sum=0.0D0
        Do 10 M=0,NPts-1
         KMode=M
         If(M.gt.NPts/2) KMode=M-NPts
         Wave=2.0D0*Pi*DFloat(KMode)/Length
         Ang=2.0D0*Pi*DFloat(KMode*(I-J))/DFloat(NPts)
         Sum=Sum+0.5D0*Wave*Wave*DCos(Ang)
   10   Continue
        Ham(IJ)=H2CM*Sum/DFloat(NPts)
        If(I.eq.J) Ham(IJ)=Ham(IJ)+Vpot(I)
   20  Continue
   30 Continue
      Return
      End

*Deck BuildHam2D
      Subroutine BuildHam2D(IOut,N1,N2,IB1,IB2,Q1,Q2,DX1,DX2,V2D,G11,
     $ G22,G12,Ham)
      Implicit Real*8 (A-H,O-Z)
      Parameter (MaxAx=300,MaxTri=45150,H2CM=219474.6313705D0)
      Dimension Q1(*),Q2(*),V2D(*),Ham(*)
      Dimension T1(MaxTri),T2(MaxTri),Zero1(MaxAx),Zero2(MaxAx)
      If(N1.gt.MaxAx.or.N2.gt.MaxAx) Stop '2D axis too large'
      Do 5 I=1,MaxAx
       Zero1(I)=0.0D0
       Zero2(I)=0.0D0
    5 Continue
      Pi=4.0D0*DAtan(1.0D0)
      If(IB1.eq.0) then
       Call BuildSinc(N1,Zero1,DX1,T1,H2CM,Pi)
      Else
       Call BuildFourier(N1,Zero1,Q1(N1)+DX1-Q1(1),T1,H2CM,Pi)
      EndIf
      If(IB2.eq.0) then
       Call BuildSinc(N2,Zero2,DX2,T2,H2CM,Pi)
      Else
       Call BuildFourier(N2,Zero2,Q2(N2)+DX2-Q2(1),T2,H2CM,Pi)
      EndIf
      NTot=N1*N2
      IJ=0
      Do 30 A=1,NTot
       I1=(A-1)/N2+1
       J1=A-(I1-1)*N2
       Do 20 B=1,A
        I2=(B-1)/N2+1
        J2=B-(I2-1)*N2
        IJ=IJ+1
        Val=0.0D0
        If(J1.eq.J2) Val=Val+G11*TriGet(T1,I1,I2)
        If(I1.eq.I2) Val=Val+G22*TriGet(T2,J1,J2)
        If(A.eq.B) Val=Val+V2D(A)
C       Cross metric is handled in Python for now unless it is zero.
        If(DAbs(G12).gt.0.0D0) then
         If(A.eq.B) Val=Val+0.0D0
        EndIf
        Ham(IJ)=Val
   20  Continue
   30 Continue
      Return
      End

*Deck TriGet
      Function TriGet(A,I,J)
      Implicit Real*8 (A-H,O-Z)
      Dimension A(*)
      II=I
      JJ=J
      If(JJ.gt.II) then
       IT=II
       II=JJ
       JJ=IT
      EndIf
      K=II*(II-1)/2+JJ
      TriGet=A(K)
      Return
      End

*Deck InitGauss
      Subroutine InitGauss(NBas,XMin,XMax,Cen,Alp)
      Implicit Real*8 (A-H,O-Z)
      Dimension Cen(*),Alp(*)
      DX=(XMax-XMin)/DFloat(Max0(NBas-1,1))
      If(DX.le.0.0D0) Stop 'Bad Gaussian interval'
      A0=1.0D0/(DX*DX)
      Do 10 I=1,NBas
       Cen(I)=XMin+DFloat(I-1)*DX
       Alp(I)=A0
   10 Continue
      Return
      End

*Deck GaussPhi
      Function GaussPhi(X,C,A)
      Implicit Real*8 (A-H,O-Z)
      Pi=4.0D0*DAtan(1.0D0)
      GaussPhi=(2.0D0*A/Pi)**0.25D0*DExp(-A*(X-C)*(X-C))
      Return
      End

*Deck BuildGaussSH
      Subroutine BuildGaussSH(NGrid,NBas,Q,Vpot,DX,Cen,Alp,S,H)
      Implicit Real*8 (A-H,O-Z)
      Parameter (H2CM=219474.6313705D0)
      Dimension Q(*),Vpot(*),Cen(*),Alp(*),S(*),H(*)
      IJ=0
      Do 30 I=1,NBas
       Do 20 J=1,I
        IJ=IJ+1
        Ai=Alp(I)
        Aj=Alp(J)
        As=Ai+Aj
        D2=(Cen(I)-Cen(J))**2
        Ovl=DSqrt(2.0D0)*(Ai*Aj)**0.25D0/DSqrt(As)
     $   *DExp(-Ai*Aj*D2/As)
        T=H2CM*Ovl*(Ai*Aj/As)*(1.0D0-2.0D0*Ai*Aj*D2/As)
        VInt=0.0D0
        Do 10 K=1,NGrid
         W=DX
         If(K.eq.1.or.K.eq.NGrid) W=0.5D0*DX
         VInt=VInt+W*GaussPhi(Q(K),Cen(I),Ai)*Vpot(K)
     $    *GaussPhi(Q(K),Cen(J),Aj)
   10   Continue
        S(IJ)=Ovl
        H(IJ)=T+VInt
   20  Continue
   30 Continue
      Return
      End

*Deck BuildOrtho
      Subroutine BuildOrtho(N,SEval,SVec,XMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension SEval(*),SVec(N,*),XMat(N,*)
      Thr=1.0D-10
      Do 20 I=1,N
       Den=DSqrt(Max(SEval(I),Thr))
       Do 10 J=1,N
        XMat(J,I)=SVec(J,I)/Den
   10  Continue
   20 Continue
      Return
      End

*Deck TransformHam
      Subroutine TransformHam(N,H,X,HOr)
      Implicit Real*8 (A-H,O-Z)
      Integer A,B
      Dimension H(*),X(N,*),HOr(*)
      IJ=0
      Do 30 I=1,N
       Do 20 J=1,I
        IJ=IJ+1
        Sum=0.0D0
        Do 11 A=1,N
         Do 10 B=1,N
          Sum=Sum+X(A,I)*TriGet(H,A,B)*X(B,J)
   10    Continue
   11   Continue
        HOr(IJ)=Sum
   20  Continue
   30 Continue
      Return
      End

*Deck WriteReport
      Subroutine WriteReport(Name,NPts,NLevels,DX,Eval)
      Implicit Real*8 (A-H,O-Z)
      Character*(*) Name
      Dimension Eval(*)
      Open(20,File='dvrout',Status='unknown')
      Write(20,'(''ORACLE Fortran77 Path DVR'')')
      Write(20,'(''Solver  : '',A)') Name
      Write(20,'(''Points  :'',I8)') NPts
      Write(20,'(''Levels  :'',I8)') NLevels
      Write(20,'(''DX      :'',D20.10)') DX
      Write(20,'(/,''Lowest levels (cm-1)'')')
      Do 10 I=1,NLevels
       Write(20,'(I6,2D20.10)') I-1,Eval(I),Eval(I)-Eval(1)
   10 Continue
      Close(20)
      Return
      End

*Deck WriteLevels
      Subroutine WriteLevels(NLevels,Eval)
      Implicit Real*8 (A-H,O-Z)
      Dimension Eval(*)
      Open(30,File='dvr_levels.csv',Status='unknown')
      Write(30,'(A)') 'state,energy_cm-1,energy_above_ground_cm-1'
      Do 10 I=1,NLevels
       Write(30,'(I8,'','',D24.14,'','',D24.14)')
     $  I-1,Eval(I),Eval(I)-Eval(1)
   10 Continue
      Close(30)
      Return
      End

*Deck WriteVectors
      Subroutine WriteVectors(NPts,NLevels,Q,Vpot,EVec)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q(*),Vpot(*),EVec(NPts,*)
      Open(40,File='dvr_vectors.csv',Status='unknown')
      Write(40,'(A)') 'point,q_au,potential_cm-1,state,amplitude'
      Do 20 I=1,NPts
       Do 10 K=1,NLevels
        Write(40,'(I8,'','',D24.14,'','',D24.14,'','',I8,'','',
     $   D24.14)') I,Q(I),Vpot(I),K-1,EVec(I,K)
   10  Continue
   20 Continue
      Close(40)
      Return
      End

*Deck WriteVectors2D
      Subroutine WriteVectors2D(N1,N2,NLevels,Q1,Q2,V2D,EVec)
      Implicit Real*8 (A-H,O-Z)
      Integer P
      Dimension Q1(*),Q2(*),V2D(*),EVec(N1*N2,*)
      Open(40,File='dvr_vectors.csv',Status='unknown')
      Write(40,'(A)')
     $ 'point,q1_au,q2_au,potential_cm-1,state,amplitude'
      Do 30 I=1,N1
       Do 20 J=1,N2
        P=(I-1)*N2+J
        Do 10 K=1,NLevels
         Write(40,'(I8,'','',D24.14,'','',D24.14,'','',D24.14,
     $    '','',I8,'','',D24.14)')
     $    P,Q1(I),Q2(J),V2D(P),K-1,EVec(P,K)
   10   Continue
   20  Continue
   30 Continue
      Close(40)
      Return
      End

*Deck WriteGaussVectors
      Subroutine WriteGaussVectors(NGrid,NBas,NLevels,Q,Vpot,Cen,Alp,
     $ XMat,EVec,Phi)
      Implicit Real*8 (A-H,O-Z)
      Dimension Q(*),Vpot(*),Cen(*),Alp(*),XMat(NBas,*),EVec(NBas,*)
      Dimension Phi(*)
      Open(40,File='dvr_vectors.csv',Status='unknown')
      Write(40,'(A)') 'point,q_au,potential_cm-1,state,amplitude'
      Do 30 I=1,NGrid
       Do 5 B=1,NBas
        Phi(B)=GaussPhi(Q(I),Cen(B),Alp(B))
    5  Continue
       Do 20 K=1,NLevels
        Amp=0.0D0
        Do 10 B=1,NBas
         C=0.0D0
         Do 9 J=1,NBas
          C=C+XMat(B,J)*EVec(J,K)
    9    Continue
         Amp=Amp+Phi(B)*C
   10   Continue
        Write(40,'(I8,'','',D24.14,'','',D24.14,'','',I8,'','',
     $   D24.14)') I,Q(I),Vpot(I),K-1,Amp
   20  Continue
   30 Continue
      Close(40)
      Return
      End
