*Deck DvrAClear
      Subroutine DvrAClear(N,A)
      Implicit Real*8 (A-H,O-Z)
      Dimension A(*)
      Do 10 I=1,N
       A(I)=0.0D0
   10 Continue
      Return
      End
*Deck DvrEpsEta
      Subroutine DvrEpsEta(Eps,Eta)
      Implicit Real*8(A-H,O-Z)
      Save Zero, One, Two
      Data Zero/0.0D0/, One/1.0D0/, Two/2.0D0/
      Eta = One
   10 T = Eta / Two
      If(T.ne.Zero) then
       Eta = T
       Goto 10
      EndIf
      Eps = One
   20 T = One + (Eps/Two)
      If(T.ne.One) then
       Eps = Eps / Two
       Goto 20
      EndIf
      Return
      End
*Deck DvrMDCutO
      Function DvrMDCutO(I)
      Implicit Real*8 (A-H,O-Z)
      If(I.eq.0) then
       DvrMDCutO = 1.0D-12
      ElseIf(I.eq.1) then
       DvrMDCutO = 1.0D-30
      Else
       DvrMDCutO = 1.0D-6
      EndIf
      Return
      End
*Deck DVRHQRII
      Subroutine DVRHQRII (IOut,N,IEV1,IEVL,IOrd,AL,EVal,NVX,EVec,
     $  Schmid,IErr,IX,WA,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C Householder QR Inverse Interation method 1
C     Compute the eigenvalues (EVal) of a symmetric matrix stored in a
C     lower triangular form (AL)
C     Based on the original HQRII:
C       Y. Beppu and I. Ninomiya, Comput. Chem., 6, 87 (1982)
C     Improved by A.V. Bunge and C.F. Bunge, Comput. Chem., 10, 259 (1986)
C
C Input:
C     N      : Dimension of the symmetric matrix A (NxN)
C     IEV1   : Index of the first wanted eigenvector (default 1)
C     IEVL   : Index of the last wanted eigenvector
C       > If IEVL < IEV1, eigenvectors are not computed
C     IOrd   : Order of the eigenvalues
C              IOrd >= 0: eigenvalues in decreasing order
C              IOrd < 0:  eigenvalues in increasing order
C     AL     : (Nx(N+1)/2) lower triangular array of the symmetric
C              matrix form stored as a vector
C     NVX    : number of columns for EVec
C     Schmid : if False, the eigenvectors of degenerate or nearly-degenerate
C              are not orthogonalized among themselves. True is the normal
C              use
C
C Output:
C     EVal   : (N) eigenvalues of A, in the order chosen with IOrd
C     EVec   : (NVX,N) eigenvectors of A
C     IErr   : Error flag
C                0    normal succesful completion
C               33    Error: AL is a null matrix
C              129    Error: N, LV or NVX are outside permissible bounds
C
C Local:
C     IX     : (N) working array to store various indexes
C     WA     : (6,N) working array for the actual calculations
C     V      : Scratch array
C
C
C     Dimension
      Integer MDV, N, NTT, NVX
C     Input
      Integer IEV1, IEVL, IOrd, IOut
      Logical Schmid
C     Output
      Integer IErr
      Real*8 AL(*), EVal(*), EVec(*)
C     Local
      Integer IX(*), IV, j, LV, NV
      Real*8 WA(6,*), V(*), Boundr, Eps0, Eta, GoldRt, DvrMDCutO
      Real*8 One, Pt5, Zero
C
      Save Boundr, GoldRt, One, Pt5, Zero
      Data Boundr/1.0D-6/, GoldRt/0.618033988749894D0/, One/1.0D0/,
     $  Pt5/0.5D0/, Zero/0.0D0/
C
 9100 Format(' ERROR: At least one dim. in DVRHQRII is outside the ',
     $  'permissible bounds')
 9110 Format(' ERROR: Null matrix in subroutine DVRHQRII')
C
      NTT = N*(N+1)/2
C     Call TstCor(NTT,MDV,'DVRHQRII')
      If(MDV.lt.NTT) then
       write(IOut,'('' Too Small Scratch Vector'')')
      endif
      Call DvrAClear(NTT,V)
C     Test0 to control null values
      Test0 = DvrMDCutO(1)
C     Eps0 is smallest number so that 1.+Eps0 =/= 1.0
      Call DvrEpsEta(Eps0,Eta)
      IV = IEV1
      LV = IEVL
      NV = NVX
      IErr = 0
      If(IV.le.0) IV = 1
c      If (N.lt.1 .or. N.gt.NX .or. (N.gt.NV .and. LV.ge.IV) .or.
c     $    LV.gt.N) then
C     Control that dimensions are acceptable
      If(N.lt.1 .or. (N.gt.NV .and. LV.ge.IV) .or. LV.gt.N) then
        Write(IOut, 9100)
        IERR = 129
        Return
      endIf
C     Obvious case, N = 1
      If(N.eq.1) then
         EVal(1) = AL(1)
         EVec(1) = One
         Return
      endIf
C
C     General Case
      IX(1) = 0
      Do 10 j = 2, N
       IX(j) = IX(j-1) + j - 1
   10 Continue
      NM1 = N - 1
      NVF = (LV-1)*NV
      If(NVF.lt.0) NVF = 0
      If(N.gt.2) Then
        NM2 = N - 2
        Do 100 k = 1, NM2
          KP1 = k + 1
          WA(2,K) = AL(k+IX(k))
          Scale = Zero
          Do 110 j = KP1, N
           Scale = Scale + Abs(AL(IX(j)+k))
  110     Continue
          WA(1,K) = AL(IX(KP1)+k)
          If(Scale.gt.Zero) then
            ScaleI = One/Scale
            Sum = Zero
            Do 111 j = KP1, N
              WA(2,J) = AL(IX(j)+k)*ScaleI
              Sum  = Sum + WA(2,J)**2
  111       Continue
            S = Sign(Sqrt(Sum),WA(2,KP1))
            WA(1,K)   = -S*Scale
            WA(2,KP1) = WA(2,KP1) + S
            AL(IX(KP1)+k) = WA(2,KP1)*Scale
            H    = WA(2,KP1)*S
            HUNS = (H*Scale)*Scale
            HI   = One/H
            SumM = Zero
            Do 112 ii = KP1, N
             WA(5,II) = Zero
  112       Continue
            Do 113 i = KP1, N
              IM1   = i - 1
              Sum   = ZERO
              I0    = IX(i)
              W2I   = WA(2,i)
              NRest = Mod(IM1-KP1+1,6)
              Do 120 j = KP1, KP1+NRest-1
                Sum = Sum + WA(2,j)*AL(I0+j)
                WA(5,J) = WA(5,j) + W2I*AL(I0+j)
  120         Continue
              Do 121 j = KP1+NRest, IM1, 6
               Sum = Sum + WA(2,j)*AL(I0+j) + WA(2,j+1)*AL(I0+j+1)
     C                + WA(2,j+2)*AL(I0+j+2) + WA(2,j+3)*AL(I0+j+3)
     C                + WA(2,j+4)*AL(I0+j+4) + WA(2,j+5)*AL(I0+j+5)
  121         Continue
              Do 122 j = KP1+NRest, IM1, 6
                WA(5,j)   = WA(5,j)   + W2I*AL(I0+j)
                WA(5,j+1) = WA(5,j+1) + W2I*AL(I0+j+1)
                WA(5,j+2) = WA(5,j+2) + W2I*AL(I0+j+2)
                WA(5,j+3) = WA(5,j+3) + W2I*AL(I0+j+3)
                WA(5,j+4) = WA(5,j+4) + W2I*AL(I0+j+4)
                WA(5,j+5) = WA(5,j+5) + W2I*AL(I0+j+5)
  122         Continue
              WA(6,i) = W2I*AL(I0+i) + Sum
  113       Continue
            Do 114 i = KP1, N
              WA(1,i) = (WA(5,i)+WA(6,i))*HI
              SumM    = WA(1,i)*WA(2,i) + SumM
  114       Continue
            U = Pt5*SumM*HI
            Do 115 i = KP1,N
              I0      = IX(i)
              WA(1,i) = WA(2,i)*U - WA(1,i)
              W1I     = WA(1,i)
              W2I     = WA(2,i)
              NRest   = Mod(i-KP1+1,6)
              Do 123 j = KP1, KP1+NRest-1
               AL(I0+j) = W2I*WA(1,j) + W1I*WA(2,j) + AL(I0+j)
  123         Continue
              Do 124 j = KP1+NRest, i, 6
                AL(I0+j)   = W2I*WA(1,j)   + W1I*WA(2,j)   + AL(I0+j)
                AL(I0+j+1) = W2I*WA(1,j+1) + W1I*WA(2,j+1) + AL(I0+j+1)
                AL(I0+j+2) = W2I*WA(1,j+2) + W1I*WA(2,j+2) + AL(I0+j+2)
                AL(I0+j+3) = W2I*WA(1,j+3) + W1I*WA(2,j+3) + AL(I0+j+3)
                AL(I0+j+4) = W2I*WA(1,j+4) + W1I*WA(2,j+4) + AL(I0+j+4)
                AL(I0+j+5) = W2I*WA(1,j+5) + W1I*WA(2,j+5) + AL(I0+j+5)
  124         Continue
  115       Continue
          else
            HUnS = Zero
          endIf
          AL(IX(k)+k) = HUnS
  100   Continue
      endIf
      NM1NM1    = IX(NM1) + NM1
      NM1N      = IX(N)   + NM1
      NN        = NM1N    + 1
      WA(2,NM1) = AL(NM1NM1)
      WA(2,N)   = AL(NN)
      WA(1,NM1) = AL(NM1N)
      WA(1,N)   = Zero
      GERSCH    = Abs(WA(2,1)) + Abs(WA(1,1))
      Do 200 i = 1,NM1
       GERSCH = Max(Abs(WA(2,i+1))+Abs(WA(1,i))+Abs(WA(1,i+1)),GERSCH)
  200 Continue
C     Trap null matrix before it is too late.
      If(GERSCH.lt.Test0) then
        Write(IOut, 9110)
        IErr = 33
        Return
      endIf
      SumD   = Zero
      SumCOD = Zero
      Do 300 i = 1, N
        SumCOD =  SumCOD + Abs(WA(1,i))
        SumD   =  SumD   + Abs(WA(2,i))
  300 Continue
      Scale  = SumD + SumCOD
      ScaleI = One/Scale
      Do 400 i = 1, N
        WA(1,i)     = WA(1,i)*ScaleI
        WA(2,i)     = WA(2,i)*ScaleI
        WA(3,i)     = WA(1,i)
        EVal(i)     = WA(2,i)
        EVec(i+NVF) = EVal(i)
  400 Continue
      Eps    = Sqrt(Eps0)
      GERSCH = GERSCH*ScaleI
      Del    = GERSCH*Eps
      DelW5  = GERSCH*Eps0
      If(SumD/SumCOD.gt.Del) Then
C       QR method with origin shift.
        Do 500 k = N, 2, -1
  510     Continue
            l = k
  520       If (Abs(WA(3,l-1)).gt.Del) then
              l = l - 1
              If(l.gt.1) Goto 520
            endIf
            If(l.ne.k) then
              WW        = (EVal(k-1)+EVal(k))*Pt5
              R         = EVal(k) - WW
              Z         = WW - Sign(Sqrt(WA(3,k-1)**2 + R*R),WW)
              EE        = EVal(l) - Z
              EVal(l)   = EE
              FF        = WA(3,l)
              R         = Sqrt(EE*EE + FF*FF)
              RI        = One/R
              C         = EVal(l)*RI
              S         = WA(3,l)*RI
              WW        = EVal(l+1) - Z
              EVal(l)   = (FF*C + WW*S)*S + EE + Z
              EVal(l+1) = C*WW - S*FF
              Do 521 j = l+1, k-1
                R         = Sqrt(EVal(j)**2 + WA(3,j)**2)
                RI        = One/R
                WA(3,j-1) = S*R
                EE        = EVal(j)*C
                FF        = WA(3,j)*C
                C         = EVal(j)*RI
                S         = WA(3,j)*RI
                WW        = EVal(j+1) - Z
                EVal(j)   = (FF*C + WW*S)*S + EE +Z
                EVal(j+1) = C*WW - S*FF
  521         Continue
              WA(3,k-1) = EVal(k)*S
              EVal(k)   = EVal(k)*C + Z
              Goto 510
            endIf
  500   Continue
C       Straight selection sort of eigenvalues.
        Sorter = One
        If(IOrd.lt.0) Sorter = -One
        j = N
  600   Continue
          l  = 1
          ii = 1
          ll = 1
          Do 610 i = 2, j
            If((EVal(i)-EVal(l))*Sorter.le.Zero) Then
              l = i
            else
              ii = i
              ll = l
            endIf
  610     Continue
          If(ii.ne.ll) Then
            WW       = EVal(ll)
            EVal(ll) = EVal(ii)
            EVal(ii) = WW
          endIf
          j = ii - 1
        If(j.gt.1) GoTo 600
      endIf
      If(LV.ge.IV) Then
C       Inverse iteration for eigenvectors.
        FN   = FLOAT(N)
        Eps1 = Sqrt(FN)*Eps
        SEps = Sqrt(Eps)
        Eps2 = (GERSCH*Boundr)/(FN*SEps)
        RN   = Zero
        RA   = Eps*GoldRt
        i2   = (IV-2)*NV
        Do 700 i = IV, LV
          i2 = i2 + NV
          Do 710 j = 1, N
            WA(3,j) = Zero
            WA(4,j) = WA(1,j)
            WA(5,j) = EVec(NVF+j) - EVal(i)
            RN      = RN + RA
            If(RN.ge.Eps) RN = RN - Eps
            WA(6,j) = RN
  710     Continue
          Do 711 j = 1, NM1
            If(Abs(WA(5,j)).le.Abs(WA(1,j))) Then
              If(Abs(WA(1,j)).lt.Test0) WA(1,j) = Del
              WA(2,j)   = -WA(5,j)/WA(1,j)
              WA(5,j)   =  WA(1,j)
              T         =  WA(5,j+1)
              WA(5,j+1) =  WA(4,j)
              WA(4,j)   =  T
              WA(3,j)   =  WA(4,j+1)
              If(Abs(WA(3,j)).lt.Test0) WA(3,j) = Del
              WA(4,j+1)=  Zero
            else
              WA(2,j)   = -WA(1,j)/WA(5,j)
            endIf
            WA(4,j+1) = WA(3,j)*WA(2,j) + WA(4,j+1)
            WA(5,j+1) = WA(4,j)*WA(2,j) + WA(5,j+1)
  711     Continue
          If(Abs(WA(5,N)).lt.Test0) WA(5,N) = DelW5
          WNM15I = One/WA(5,NM1)
          WN5I   = One/WA(5,N)
          Do 712 Itere=1,2
            If(Itere.ne.1) Then
              Do 720 j = 1, NM1
                If(WA(5,j).eq.WA(1,j)) Then
                  T         = WA(6,j)
                  WA(6,j)   = WA(6,j+1)
                  WA(6,j+1) = T
                endIf
                WA(6,j+1) = WA(6,j)*WA(2,j) + WA(6,j+1)
  720         Continue
            endIf
            WA(6,N)   = WA(6,N)*WN5I
            WA(6,NM1) = (WA(6,NM1)-WA(6,N)*WA(4,NM1))*WNM15I
            VN        = Max(Abs(WA(6,N)),Abs(WA(6,NM1)))
            Do 721 k = NM2, 1, -1
              WA(6,k) = (WA(6,k)-WA(6,k+1)*WA(4,k)-WA(6,k+2)*WA(3,k))
     $          / WA(5,k)
              VN = Max(Abs(WA(6,k)),VN)
  721       Continue
            S     = Eps1/VN
            NRest = Mod(N,6)
            Do 722 j = 1, NRest
             WA(6,j) = S*WA(6,j)
  722       Continue
            Do 723 j = 1+NRest, N, 6
              WA(6,j)   = S*WA(6,j)
              WA(6,j+1) = S*WA(6,j+1)
              WA(6,j+2) = S*WA(6,j+2)
              WA(6,j+3) = S*WA(6,j+3)
              WA(6,j+4) = S*WA(6,j+4)
              WA(6,j+5) = S*WA(6,j+5)
  723       Continue
  712     Continue
          Do 713 j = 1, N
           EVec(i2+j) = WA(6,j)
  713     Continue
  700   Continue
C Build indexing and upper triangular matrix
        IX(1) = 0
        Do 800 j = 2, N
         IX(j) = IX(j-1) - j + 1 + N
  800   Continue
        ij = 0
        Do 900 j = 1, N
          Do 910 i = 1, j
            ij = ij + 1
            V(IX(i)+j) = AL(ij)
  910     Continue
  900   Continue
C Back transformation of eigenvectors.
        ig = 1
        i2 = (IV-2)*NV
        Do 1000 i = IV, LV
          i2 = i2 + NV
          Do 1010 j = 1, N
           WA(6,j) = EVec(i2+j)
 1010     Continue
          IM1 = i - 1
          If(N.gt.2) Then
            Do 1020 j = 1, NM2
              k  = N - j - 1
              k0 = IX(k)
              If(V(k0+k).ne.Zero) Then
                KP1   = k + 1
                Sum   = Zero
                NRest = Mod(N-KP1+1,6)
                Do 1030 kk = KP1, KP1+NRest-1
                 Sum = V(K0+kk)*WA(6,kk) + Sum
 1030           Continue
                Do 1031 kk = KP1+NRest, N, 6
                 Sum = V(k0+kk)*WA(6,kk) + V(k0+kk+1)*WA(6,KK+1)
     $                  + V(k0+kk+2)*WA(6,kk+2) + V(k0+kk+3)*WA(6,kk+3)
     $                  + V(k0+kk+4)*WA(6,kk+4) + V(k0+kk+5)*WA(6,kk+5)
     $                  + Sum
 1031           Continue
                S = -Sum/V(k0+k)
                Do 1032 kk = KP1, KP1+NRest-1
                 WA(6,kk) =  S*V(k0+kk) + WA(6,kk)
 1032           Continue
                Do 1033 kk = KP1+NRest, N, 6
                  WA(6,kk)   =  S*V(k0+kk)   + WA(6,kk)
                  WA(6,kk+1) =  S*V(k0+kk+1) + WA(6,kk+1)
                  WA(6,kk+2) =  S*V(k0+kk+2) + WA(6,kk+2)
                  WA(6,kk+3) =  S*V(k0+kk+3) + WA(6,kk+3)
                  WA(6,kk+4) =  S*V(k0+kk+4) + WA(6,kk+4)
                  WA(6,kk+5) =  S*V(k0+kk+5) + WA(6,kk+5)
 1033           Continue
              endIf
 1020       Continue
          endIf
          j = ig
 1011     If (Abs(EVal(j)-Eval(i)).ge.Eps2) then
            j = j + 1
            If(j.le.i) GoTo 1011
          endIf
          ig = Min(j,i)
          NRest = Mod(N,6)
C
          If(ig.ne.i .and. Schmid) Then
C Degenerate eigenvalues.First,orthogonalize.
            KF = (ig-2)*NV
            Do 1012 k = ig, IM1
              KF  = KF + NV
              Sum = Zero
              Do 1021 j = 1, NRest
               Sum = EVec(KF+j)*WA(6,j) + Sum
 1021         Continue
              Do 1022 J = 1+ NRest, N, 6
                Sum = EVec(KF+j)*WA(6,j)     + EVec(KF+j+1)*WA(6,j+1)
     *              + EVec(KF+j+2)*WA(6,j+2) + EVec(KF+j+3)*WA(6,j+3)
     *              + EVec(KF+j+4)*WA(6,j+4) + EVec(KF+j+5)*WA(6,j+5)
     $              + Sum
 1022          Continue
              S = -Sum
              Do 1023 j = 1, NRest
               WA(6,j) = S*EVec(KF+j) + WA(6,j)
 1023         Continue
              Do 1024 j = 1+NRest, N, 6
                WA(6,j)   = S*EVec(KF+j)   + WA(6,j)
                WA(6,j+1) = S*EVec(KF+j+1) + WA(6,j+1)
                WA(6,j+2) = S*EVec(KF+j+2) + WA(6,j+2)
                WA(6,j+3) = S*EVec(KF+j+3) + WA(6,j+3)
                WA(6,j+4) = S*EVec(KF+j+4) + WA(6,j+4)
                WA(6,j+5) = S*EVec(KF+j+5) + WA(6,j+5)
 1024         Continue
 1012       Continue
          endIf
C Normalization
          Sum = Zero
          Do 1013 j = 1, NRest
           Sum = WA(6,j)**2 + Sum
 1013     Continue
          Do 1014 j = 1+NRest, N, 6
           Sum = WA(6,j  )**2 + WA(6,j+1)**2 + WA(6,j+2)**2
     *           + WA(6,j+3)**2 + WA(6,j+4)**2 + WA(6,j+5)**2 + Sum
 1014     Continue
          S = One/Sqrt(Sum)
          Do 1015 j = 1, NRest
           EVec(i2+j) = S*WA(6,j)
 1015     Continue
          Do 1016 j = 1+NRest, N, 6
            EVec(i2+j  ) = S*WA(6,j  )
            EVec(i2+j+1) = S*WA(6,j+1)
            EVec(i2+j+2) = S*WA(6,j+2)
            EVec(i2+j+3) = S*WA(6,j+3)
            EVec(i2+j+4) = S*WA(6,j+4)
            EVec(i2+j+5) = S*WA(6,j+5)
 1016     Continue
 1000   Continue
      endIf
      Do 1100 i = 1, N
       EVal(i) = Scale*EVal(i)
 1100 Continue
      Return
      End
