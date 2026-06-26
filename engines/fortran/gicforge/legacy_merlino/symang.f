*Deck MkGNCA 
      Subroutine MkGNCA(IOut,IPrint,MxBnd,MxGIcA,MxTerA,MaxAtA,
     $  NAtoms,NCyc,MxAtCy,NAtC,ICAt,NBond,NGICA,IBond,NTermA,IAtomA,
     $  IAn,IAtCyc,ITVA,CoefA,C,EAN,TreshL,DoLocSVD)
      Implicit Real*8 (A-H,O-Z)
      Logical DoLocSVD
      Logical IsEndoAngleInCycles
      Dimension NAtC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*),IAn(*),IAtCyc(*)
      Dimension NTermA(*),IAtomA(MaxAtA,MxTerA,*),ITVA(*) 
      Dimension C(3,*),CoefA(MxTerA,*),EAN(*)
      pi = dacos(-1.d0)
      ToDeg=1.80d+2/pi
C     NGicA=0
C Build Valence Angles 
      write(IOut,'(/,'' Exocyclic Valence Angles'')')
      write(IOut,'('' Center  Symmetry  Valence Angles'')')
      If(DoLocSVD) then
       Call MkGNCALocSVD(IOut,IPrint,MxBnd,MxGIcA,MxTerA,MaxAtA,
     $  NAtoms,NCyc,MxAtCy,NAtC,ICAt,NBond,NGICA,IBond,NTermA,IAtomA,
     $  IAtCyc,ITVA,CoefA,C,TreshL)
       GoTo 80
      EndIf
      Do 30 JAt=1,NAtoms
       NBJ=NBond(JAt)
       if(NBJ.eq.1) go to 30
       if(NBJ.eq.2) then
        if(IAtCyc(JAt).gt.0) go to 30
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
          If(IsEndoAngleInCycles(NCyc,MxAtCy,NAtC,ICAt,IAt,JAt,KAt))
     $     go to 50
          If(Value.gt.TreshL) go to 50 
          NGicA=NGicA+1
          NTermA(NGicA)=1
          IAtomA(1,1,NGicA)=IAt
          IAtomA(2,1,NGicA)=JAt
          IAtomA(3,1,NGicA)=KAt 
          CoefA(1,NGicA)=1.0d0
          NFree=1
          write(IOut,'(I4,5X,2I4,14X,I2)') JAt,IAt,KAt,NFree
   50    continue
   40   continue
        go to 30
       endif
       if(NBJ.eq.3) then
C       call Gen3At(IOut,IPrint,MxBnd,MaxAtA,MxTerA,NGICA,JAt,IAn,
C    $    NBond,IBond,IAtCyc,NTermA,IAtomA,CoefA)
        call C2V3At(Iout,IPrint,MxBnd,MaxAtA,MxTerA,NGICA,JAt,
     $    IAn,NBond,IBond,IatCyc,NTermA,IAtomA,ITVA,CoefA,EAN)
        go to 30
       endif 
       if(NBJ.eq.4) then
        call FourAt(IOut,IPrint,MxBnd,MaxAtA,MxTerA,NGICA,JAt,IAn,
     $   NBond,IBond,IAtCyc,NTermA,IAtomA,ITVA,CoefA,C,EAN)
        go to 30
       endif 
       if(NBJ.gt.4) then
        call HighCoordAt(IOut,IPrint,MxBnd,MaxAtA,MxTerA,NGICA,JAt,
     $   NBond,IBond,NTermA,IAtomA,ITVA,CoefA,C)
        go to 30
       endif
   30 continue
   80 continue
      If(IPrint.gt.0) then
       write(IOut,'(/,I5,'' Angle GNICS'')') NGICA
       if(NGICA.eq.0) return
       Do 60 IGICA=1,NGICA
        NREDA=NTermA(IGICA)
        write(IOut,'('' AngGNC('',I3,'')'')') IGICA
        Do 70 IREDA=1,NREDA
         I1=IAtomA(1,IRedA,IGICA)
         I2=IAtomA(2,IRedA,IGICA)
         I3=IAtomA(3,IRedA,IGICA)
         Coef=CoefA(IRedA,IGICA)
         Value=ValAng(C(1,I1),C(1,I2),C(1,I3))*ToDeg
         write (IOUT,'(F8.4,'' A('',2(I3,'',''),I3,'')'',2X,
     $     ''Value: '',F8.3)') Coef,I1,I2,I3,Value 
  70    continue
  60   continue
      endif
      return
      End
*Deck MkGNCALocSVD
      Subroutine MkGNCALocSVD(IOut,IPrint,MxBnd,MxGIcA,MxTerA,
     $ MaxAtA,NAtoms,NCyc,MxAtCy,NAtC,ICAt,NBond,NGICA,IBond,NTermA,
     $ IAtomA,IAtCyc,ITVA,CoefA,C,TreshL)
      Implicit Real*8 (A-H,O-Z)
      Parameter(MxLoc=45,MxCart=3000)
      Dimension NAtC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*),IAtCyc(*)
      Dimension NTermA(*),IAtomA(MaxAtA,MxTerA,*),ITVA(*)
      Dimension CoefA(MxTerA,*),C(3,*)
      Dimension LAt1(MxLoc),LAt2(MxLoc),BLoc(MxCart,MxLoc)
      Dimension G(MxLoc,MxLoc),EVal(MxLoc),EVec(MxLoc,MxLoc)
      Dimension B(3,4),DB(3,4,3,4),IB(4)
      Integer Rank
      Logical Endo,IsEndoAngleInCycles

      NCart=3*NAtoms
      If(NCart.gt.MxCart) then
       Write(IOut,'('' LOCSVD angle block skipped: too many '',
     $ ''atoms'',I6)') NAtoms
       Return
      EndIf

      Do 200 JAt=1,NAtoms
       NBJ=NBond(JAt)
       If(NBJ.le.1) GoTo 200
       NPrim=0
       Do 220 II=1,NBJ-1
        IAt=IBond(II,JAt)
        Do 230 KK=II+1,NBJ
         KAt=IBond(KK,JAt)
         Endo=IsEndoAngleInCycles(NCyc,MxAtCy,NAtC,ICAt,IAt,JAt,KAt)
         If(Endo) GoTo 230
         Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
         If(Value.gt.TreshL) GoTo 230
         If(KAt.lt.IAt) then
          LAt=IAt
          IAt=KAt
          KAt=LAt
         EndIf
         NPrim=NPrim+1
         If(NPrim.gt.MxLoc) then
          Write(IOut,'('' LOCSVD angle block too large at center'',I6)')
     $     JAt
          NPrim=MxLoc
          GoTo 240
         EndIf
         LAt1(NPrim)=IAt
         LAt2(NPrim)=KAt
  230   Continue
  220  Continue
  240  Continue
       If(NPrim.eq.0) GoTo 200

       Call AClear(MxCart*MxLoc,BLoc)
       Do 260 IP=1,NPrim
        IAt=LAt1(IP)
        KAt=LAt2(IP)
        Call AClear(12,B)
        Call AClear(144,DB)
        Call DBBend(1,IAt,JAt,KAt,B,IB,C,DB)
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        KXYZ=3*(KAt-1)
        Do 250 IC=1,3
         BLoc(IXYZ+IC,IP)=B(IC,1)
         BLoc(JXYZ+IC,IP)=B(IC,2)
         BLoc(KXYZ+IC,IP)=B(IC,3)
  250   Continue
  260  Continue

       Do 290 IP=1,NPrim
        Do 280 JP=1,NPrim
         Sum=0.0D0
         Do 270 IC=1,NCart
          Sum=Sum+BLoc(IC,IP)*BLoc(IC,JP)
  270    Continue
         G(IP,JP)=Sum
  280   Continue
  290  Continue

       Call LocSVDJacobi(MxLoc,NPrim,G,EVal,EVec,Rank)
       If(Rank.gt.0) write(IOut,'(I4,6X,A6,5X,I2)')
     $  JAt,'LOCSVD',Rank
       Do 330 IM=1,Rank
        If(NGICA.ge.MxGIcA) then
         Write(IOut,'('' Too many angle GNICs in LOCSVD'')')
         Return
        EndIf
        NGICA=NGICA+1
        ITVA(NGICA)=0
        NTerm=0
        Do 310 IP=1,NPrim
         If(DAbs(EVec(IP,IM)).le.1.0D-12) GoTo 310
         NTerm=NTerm+1
         If(NTerm.gt.MxTerA) then
          Write(IOut,'('' Too many LOCSVD angle terms at center'',I6)')
     $    JAt
          NTerm=MxTerA
          GoTo 320
         EndIf
         CoefA(NTerm,NGICA)=EVec(IP,IM)
         IAtomA(1,NTerm,NGICA)=LAt1(IP)
         IAtomA(2,NTerm,NGICA)=JAt
         IAtomA(3,NTerm,NGICA)=LAt2(IP)
  310   Continue
  320   Continue
        NTermA(NGICA)=NTerm
  330  Continue
  200 Continue
      Return
      End

*Deck IsEndoAngleInCycles
      Logical Function IsEndoAngleInCycles(NCyc,MxAtCy,NAtC,ICAt,
     $ IAt,JAt,KAt)
      Implicit None
      Integer NCyc,MxAtCy,IAt,JAt,KAt
      Integer NAtC(*),ICAt(MxAtCy,*)
      Integer ICyc,N,I,Prev,Next
      IsEndoAngleInCycles=.False.
      If(NCyc.le.0) Return
      Do 10 ICyc=1,NCyc
       N=NAtC(ICyc)
       If(N.le.2) GoTo 10
       Do 20 I=1,N
        If(ICAt(I,ICyc).ne.JAt) GoTo 20
        Prev=I-1
        If(Prev.le.0) Prev=N
        Next=I+1
        If(Next.gt.N) Next=1
        If((ICAt(Prev,ICyc).eq.IAt.and.ICAt(Next,ICyc).eq.KAt).or.
     $     (ICAt(Prev,ICyc).eq.KAt.and.ICAt(Next,ICyc).eq.IAt)) then
         IsEndoAngleInCycles=.True.
         Return
        EndIf
   20  Continue
   10 Continue
      Return
      End
*Deck C2V3At
      Subroutine C2V3At(Iout,IPrint,MxBnd,MaxAtA,MxTrmA,ICoord,IAt,
     $  IAn,NBond,IBond,IatCyc,NTermA,IAtomA,ITVA,CoefA,EAN)
      Implicit None
C               
C Dimensions
      Integer MxBnd, MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,IAn(*),NBond(*) 
      Integer IBond(MxBnd,*),IAtCyc(*),ITVA(*) 
      Real*8 EAN(*)
C Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer JAt,KAt,LAt,I1,J1,K1,L1,IDiff,ICycI,ICycJ,ICycK,ICycL,NAng 
      Logical EqJK,EqJL,EqKl
      Real*8 DenSD,DenRK,Tresh 
C
      tresh=5.0d-4
      I1=IAt
      J1=IBond(1,IAt)
      K1=IBond(2,IAt)
      L1=IBond(3,IAt)
      EqJK=Abs(EAN(J1)-EAN(K1)).lt.tresh 
      EqJL=Abs(EAN(J1)-EAN(L1)).lt.tresh
      EqKL=Abs(EAN(K1)-EAN(L1)).lt.tresh  
C IDiff= J1 for 3 equal atoms and for K1=L1
C Idiff= hydrogen or terminal atom for 3 different atoms
      If(.not.EqJK.and..not.EqJL.and..not.EqKl) then
       IDiff=J1
       If(IAn(K1).eq.1) then 
        IDiff=K1
       ElseIf(IAn(L1).eq.1) then
        IDiff=L1
       ElseIf(NBond(K1).eq.1) then
        IDiff=K1
       ElseIf(NBond(L1).eq.1) then
        IDiff=L1
       EndIf
      ElseIf(EqJK.and.EqJL) then
       IDiff=J1
      ElseIf(EqJK.and..not.EqJL) then
       IDiff=L1
      ElseIF(EqJL.and..not.EqJK) then
       IDiff=K1
      ElseIf(EqKL.and..not.EqJK) then
       IDiff=J1
      EndIf
      ICyCI=IAtCyc(I1)
      ICycJ=IAtCyc(J1)
      ICycK=IAtCyc(K1)
      ICycL=IAtCyc(L1)
C IDiff= atom not involved in cycles if I1 is involved in a cycle
      If(ICycI.ne.0) then
       If(ICycI.ne.0) IDiff=-J1
       If(ICycJ.ne.0.and.ICycK.ne.0) IDiff=-L1
       If(ICycJ.ne.0.and.ICycL.ne.0) IDiff=-K1
       If(ICycK.ne.0.and.ICycL.ne.0) IDiff=-J1
      endif
      If(Abs(IDiff).eq.J1) then
       JAt=J1
       KAt=K1
       LAt=L1
      ElseIf(Abs(IDiff).eq.K1) then
       JAt=K1
       KAt=J1
       LAt=L1
      ElseIf(Abs(IDiff).eq.L1) then
       JAt=L1
       KAt=K1
       LAt=J1
      Endif
      If(ICycI.ne.0) then
       If(IAtCyc(JAt).ne.0) then
        return
       Else
        NAng=1
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       EndIf
      Else
       NAng=2
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
C SD (Symmetric Deformation) Symmetry Coordinate
       ICoord=ICoord+1
       NTerma(ICoord)=3
       ITVA(ICoord)=1
       DenSD=Sqrt(6.0d0)
       CoefA(1,ICoord)=2.0d0/DenSD
       CoefA(2,ICoord)=-1.0d0/DenSD
       CoefA(3,ICoord)=-1.0d0/DenSD
       IAtomA(1,1,ICoord)=KAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=LAt
       IAtomA(1,2,ICoord)=JAt
       IAtomA(2,2,ICoord)=IAt
       IAtomA(3,2,ICoord)=KAt
       IAtomA(1,3,ICoord)=JAt
       IAtomA(2,3,ICoord)=IAt
       IAtomA(3,3,ICoord)=LAt
      endif 
C RK (Rocking) Symmetry Coordinate
      Icoord=ICoord+1 
      NTermA(ICoord)=2
      ITVA(ICoord)=2
      DenRK=SQrt(2.0d0)
      CoefA(1,ICoord)=1.0D0/DenRK
      CoefA(2,ICoord)=-1.0D0/DenRK
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=KAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=LAt
      return
      end
*Deck FourAt
      Subroutine FourAt(IOut,IPrint,MxBnd,MaxAtA,MxTrmA,ICoord,IAt,
     $  IAn,NBond,IBond,IAtCyc,NTermA,IAtomA,ITVA,CoefA,C,EAN)
      Implicit real*8 (A-H,O-Z)
C IAt    is the central atom 
C ICOORD is the number of already built GNICs
C Dimensions
      Integer MxBnd, MaxAtA, MxTrmA
C Input
      Integer IAn(*),NBond(*),ITVA(*)
      Integer IBond(MxBnd,*),IAtCyc(*)
      Real*8 C(3,*),EAN(*)
C Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Logical FrozJ,FrozK,FrozL,FrozM
      Logical TstK,TstL,TstM
C
      pi = dacos(-1.d0)
      ToDeg=1.80d+2/pi 
      treshL=5.0d+0
      tresh=5.0d-4
      JAt=IBond(1,IAt)
      KAt=IBond(2,IAt)
      LAt=IBond(3,IAt)
      MAt=IBond(4,IAt)
      NPivT=0 
C if IAT do not belong to cycles there are no restrictions
      if(IAtCyc(IAt).ne.0) then
       do 10 j1=1,4
        jj=IBond(j1,IAt)
        if(IAtCyc(jj).ne.0) NPivT=NPivT+1
   10  continue 
      endif
      if(NPivT.eq.4) then
       If(IAtCyc(IAt).eq.2) then
        If(IPrint.eq.1) write(IOut,'('' Spiro Compound Around Atom'',
     $    I5)')IAt
        call SpiAng(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,
     $    MAt,NTermA,IAtomA,ITVA,CoefA)
        return
       Else
        write(IOut,'('' Connectivity Error for Atom'',I5)') IAt
       EndIf
      EndIf
C     write(IOut,'(/,'' Origin. Atoms'',4I3)') (IBond(ii,iat),ii=1,4)
      call Ord4At(IOut,JAt,KAt,LAt,MAt,NPivT,NEq,IAtCyc,Tresh,EAn)
C if IAT do not belong to cycles there are no restrictions
      FrozJ=(IAtCyc(JAt).ne.0.and.IAtCyc(IAt).ne.0)
      FrozK=(IAtCyc(KAt).ne.0.and.IAtCyc(IAt).ne.0)
      FrozL=(IAtCyc(LAt).ne.0.and.IAtCyc(IAt).ne.0)
      FrozM=(IAtCyc(MAt).ne.0.and.IAtCyc(IAt).ne.0)
C     write(IOut,'(  '' Reord.  Atoms'',4I3)') JAt,KAt,LAt,MAt
C     write(IOut,'(  ''Reord.Cycles  '',4L3)') FrozJ,FrozK,FrozL,FrozM
C     write(IOut,'(  '' N. of Pivots '',I3)') NPivT
C     write(IOut,'(  '' N. of Eq.At. '',I3)') NEq
      AngJK=ValAng(C(1,JAt),C(1,IAt),C(1,KAt))*ToDeg
      tstk=Abs(AngJK-1.80d+2).lt.treshL
      AngJL=ValAng(C(1,JAt),C(1,IAt),C(1,LAt))*ToDeg
      tstl=Abs(AngJL-1.80d+2).lt.treshL
      AngJM=ValAng(C(1,JAt),C(1,IAt),C(1,MAt))*ToDeg
      tstm=Abs(AngJM-1.80d+2).lt.treshL
C D4h symmetry
      If(tstk.or.tstl.or.tstm)then 
       call D4h4At(Iout,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
     $   FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA,C) 
       return
      endif
C Td symmetry
      if(neq.eq.4) then
       call Td4At(Iout,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,
     $   MAt,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
       return
      endif
C C3v symmetry
      if(neq.eq.3.or.npivt.eq.1.or.npivt.eq.3) then
       call WXY3(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
     $  NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
c      call C3V4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
c    $   NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
       return
      endif 
C C2v or Cs symmetry
      if(neq.le.2) then
C       call C2V4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
C    $   NEq,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
        call W2XY2(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
     $    NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
       return
      endif
C Cs symmetry
c      if(neq.eq.1) then
c       call Cs4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
c     $   FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
c       return
c      endif
C C1 symmetry uses nearly equivalent atoms in cycles,or H or terminal atoms
      if(neq.eq.0) then
       call LooseS(IOut,IPrint,MxBnd,NAtoms,IAt,IAn,NBond,IBond,JAt,KAt,
     $   LAt,MAt)
       if(NPivT.eq.1.or.npivt.eq.3) then
C pseudo c3v symmetry 
        call WXY3(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
     $   NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
C       call C3V4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,
C    $   MAt,NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
       else
        NEq1=1
C pseudo Cs symmetry
C       call C2v4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,
C    $     MAt,NEq1,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
        call W2XY2(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,LAt,MAt,
     $    NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
       endif
      endif
      return
      end
*Deck HighCoordAt
      Subroutine HighCoordAt(IOut,IPrint,MxBnd,MaxAtA,MxTrmA,ICoord,IAt,
     $  NBond,IBond,NTermA,IAtomA,ITVA,CoefA,C)
      Implicit None
C
C High-coordination fallback: explicit pairwise valence angles for centers
C with coordination > 4.  These are not symmetry-adapted here; the global
C rank pruning keeps the final set non-redundant.
C
C Dimensions
      Integer IOut,IPrint,MxBnd,MaxAtA,MxTrmA,ICoord,IAt
      Integer NBond(*),IBond(MxBnd,*)
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*),C(3,*)
      Double Precision ValAng
C Local
      Integer II,KK,I1,I2,ILeft,IRight,NAng
      Real*8 Value,ToDeg,Pi
      Pi = Dacos(-1.D0)
      ToDeg = 1.80D+2 / Pi
      NAng = 0
      Do 20 II=1,NBond(IAt)-1
       I1 = IBond(II,IAt)
       Do 10 KK=II+1,NBond(IAt)
        I2 = IBond(KK,IAt)
        ILeft = I1
        IRight = I2
        If(ILeft.gt.IRight) then
         ILeft = I2
         IRight = I1
        EndIf
        ICoord = ICoord + 1
        NTermA(ICoord) = 1
        ITVA(ICoord) = 17
        CoefA(1,ICoord) = 1.0D0
        IAtomA(1,1,ICoord) = ILeft
        IAtomA(2,1,ICoord) = IAt
        IAtomA(3,1,ICoord) = IRight
        NAng = NAng + 1
        If(IPrint.gt.0) then
         Value = ValAng(C(1,ILeft),C(1,IAt),C(1,IRight))*ToDeg
         Write(IOut,'(I4,6X,''HCAn'',1X,2I4,14X,I2,4X,''Value='',F8.3)')
     $     IAt,ILeft,IRight,NAng,Value
        EndIf
   10   Continue
   20 Continue
      Return
      End
*Deck C2V4At
      Subroutine C2V4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,NEq,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
      Implicit None
C
C ICoord is the starting number of the coordinates to be built
C IAt is always the central atom
C JAt and KAt are equivalent; if negative they are frozen
C               
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
      Logical FrozJ,FrozK,FrozL,FrozM
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NEq,NAng
      Real*8 DenSci,DenRk
C
      DenSci=Sqrt(6.0d0)
      DenRk=2.0d0 
      If(FrozJ.and.FrozK) then
       NAng=4
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
      ElseIf(FrozL.and.FrozM) then
       NAng=4
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
      ElseIf(FrozJ.and.FrozL) then
       NAng=4
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
      ElseIf(FrozJ.and.FrozM) then
       NAng=4
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
      Else
       NAng=5
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
      EndIf
C Rocking symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=-1.0d0/DenRk
      CoefA(3,ICoord)=1.0d0/DenRk
      CoefA(4,ICoord)=-1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C Wagging symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=5
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=1.0d0/DenRk
      CoefA(3,ICoord)=-1.0d0/DenRk
      CoefA(4,ICoord)=-1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C Twisting symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=6
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=-1.0d0/DenRk
      CoefA(3,ICoord)=-1.0d0/DenRk
      CoefA(4,ICoord)=1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C Sciss1 symmetry coordinate
      If(.not.FrozJ.and..not.FrozK) then
       If(.not.FrozL.and..not.FrozM) then
        ICoord=ICoord+1
        NTerma(ICoord)=2
        ITVA(ICoord)=3
        CoefA(1,ICoord)=Sqrt(5.0d0)/DenSci
        CoefA(2,ICoord)=1.0d0/DenSci
        IAtomA(1,1,ICoord)=JAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=KAt
        IAtomA(1,2,ICoord)=LAt
        IAtomA(2,2,ICoord)=IAt
        IAtomA(3,2,ICoord)=MAt
C Sciss2 symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=2
        ITVA(ICoord)=3
        CoefA(1,ICoord)=1.0d0/DenSci
        CoefA(2,ICoord)=Sqrt(5.0d0)/DenSci
        IAtomA(1,1,ICoord)=JAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=KAt
        IAtomA(1,2,ICoord)=LAt
        IAtomA(2,2,ICoord)=IAt
        IAtomA(3,2,ICoord)=MAt
       else
C Sciss JIK symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=JAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=KAt
       endif
      EndIf
      If(FrozJ.and.FrozK) then 
       If(.not.FrozL.and..not.FrozM) then
C Sciss LIM symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=LAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=MAt
       EndIf
      EndIf
      If(FrozJ.and.FrozL) then
       If(.not.FrozK.and..not.FrozM) then
C Sciss KIM symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=KAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=MAt
       EndIf
      EndIf
      If(FrozJ.and.FrozM) then
       If(.not.FrozK.and..not.FrozL) then
C Sciss KIL symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=KAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=LAt
       EndIf
      EndIf
      If(FrozK.and.FrozL) then
       If(.not.FrozJ.and..not.FrozM) then
C Sciss JIM symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=JAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=MAt
       EndIf
      EndIf
      If(FrozK.and.FrozM) then
       If(.not.FrozJ.and..not.FrozL) then
C Sciss JIL symmetry coordinate
        ICoord=ICoord+1
        NTerma(ICoord)=1
        ITVA(ICoord)=4
        CoefA(1,ICoord)=1.0d0
        IAtomA(1,1,ICoord)=JAt
        IAtomA(2,1,ICoord)=IAt
        IAtomA(3,1,ICoord)=LAt
       EndIf
      EndIf
      return
      end 
*Deck C3V4At
      Subroutine C3V4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,
     $  CoefA)
      Implicit None
C               
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
      Logical FrozJ,FrozK,FrozL,FrozM 
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NAng,NEq,NPivT
      Real*8 DenSD,DenAD,DenAD1,DenRK,DenRK1
C
      If(NPivT.eq.0) then
       NAng=5
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      Elseif(NPivT.eq.1) then
       If(.not.FrozJ) then
        write(IOut,'('' Wrong Free Atom'',I5,''  With Frozen'',3I5,
     $    '' Atoms Around'',I5)') JAt,KAt,LAt,MAt,IAt
        Stop
       EndIf
       NAng=5
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      Elseif(NPivT.eq.2) then
       NAng=3
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      Elseif(NPivT.eq.3) then 
       If(FrozJ) then
        write(IOut,'('' Wrong Frozen Atom'',I5,'' Around'',I5)')JAt,IAt
        Stop
       EndIf
       NAng=2
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      Else
       NAng=0
       write(IOut,'('' No Free Angles Around Center'',I5)')IAt
       Stop
      EndIf
      DenSD=Sqrt(6.0d0)
      DenAD=Sqrt(6.0d0)
      DenAD1=Sqrt(2.0d0)
      DenRk=Sqrt(6.0d0)
      DenRk1=Sqrt(2.0d0)
      If(NPivT.lt.2) then 
C SD (Symmetric Deformation) symmetry coordinate
       ICoord=ICoord+1
       NTerma(ICoord)=6
       ITVA(ICoord)=1
       CoefA(1,ICoord)=1.0d0/DenSD
       CoefA(2,ICoord)=1.0d0/DenSD
       CoefA(3,ICoord)=1.0d0/DenSD
       CoefA(4,ICoord)=-1.0d0/DenSD
       CoefA(5,ICoord)=-1.0d0/DenSD
       CoefA(6,ICoord)=-1.0d0/DenSD 
       IAtomA(1,1,ICoord)=JAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=KAt
       IAtomA(1,2,ICoord)=JAt
       IAtomA(2,2,ICoord)=IAt
       IAtomA(3,2,ICoord)=LAt
       IAtomA(1,3,ICoord)=JAt
       IAtomA(2,3,ICoord)=IAt
       IAtomA(3,3,ICoord)=MAt
       IAtomA(1,4,ICoord)=KAt
       IAtomA(2,4,ICoord)=IAt
       IAtomA(3,4,ICoord)=LAt
       IAtomA(1,5,ICoord)=KAt
       IAtomA(2,5,ICoord)=IAt
       IAtomA(3,5,ICoord)=MAt
       IAtomA(1,6,ICoord)=LAt
       IAtomA(2,6,ICoord)=IAt
       IAtomA(3,6,ICoord)=MAt
C Rocking (RK) symmetry coordinate
       ICoord=ICoord+1
       NtermA(ICoord)=3
       ITVA(ICoord)=2
       CoefA(1,ICoord)=2.0d0/DenRK
       CoefA(2,ICoord)=-1.0d0/DenRK
       CoefA(3,ICoord)=-1.0d0/DenRK
       IAtomA(1,1,ICoord)=KAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=LAt
       IAtomA(1,2,ICoord)=KAt
       IAtomA(2,2,ICoord)=IAt
       IAtomA(3,2,ICoord)=MAt
       IAtomA(1,3,ICoord)=LAt
       IAtomA(2,3,ICoord)=IAt
       IAtomA(3,3,ICoord)=MAt
      EndIf
C Rocking (RK1) symmetry coordinate
      If(NPivT.lt.3) then
       ICoord=ICoord+1
       NtermA(ICoord)=2
       ITVA(ICoord)=2
       CoefA(1,ICoord)=1.0d0/DenRK1
       CoefA(2,ICoord)=-1.0d0/DenRK1
       IAtomA(1,1,ICoord)=KAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=MAt
       IAtomA(1,2,ICoord)=LAt
       IAtomA(2,2,ICoord)=IAt
       IAtomA(3,2,ICoord)=MAt
      endif
C AD (Asymmetric Deformation) Symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=3
      ITVA(ICoord)=7
      CoefA(1,ICoord)=2.0d0/DenAD
      CoefA(2,ICoord)=-1.0d0/DenAD
      CoefA(3,ICoord)=-1.0d0/DenAD
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=KAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=LAt
      IAtomA(1,3,ICoord)=JAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=MAt
C Asymmetric Deformation (AD1) symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=7
      CoefA(1,ICoord)=1.0d0/DenAD1
      CoefA(2,ICoord)=-1.0d0/DenAD1
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      return
      end
*Deck Cs4At
      Subroutine Cs4At(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
      Implicit None  
C
C ICoord is the starting number of the coordinates to be built
C IAt is always the central atom
C JAt and KAt are equivalent and can be frozen
C LAt, MAt can be different and frozen
C              
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
      Logical FrozJ,FrozK,FrozL,FrozM
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NAng
      Real*8 DenSci,DenRk
C
      DenSci=Sqrt(6.0d0)
      DenRk=2.0d0 
      If(FrozJ.and.FrozK) then
       NAng=4
        write(IOut,'(I4,6X,''Cs  '',1X,2I4,11X,2I4,10X,2I4,14X,I2)')
     $    IAt,JAt,KAt,JAt,KAt,LAt,MAt,NAng
      ElseIf(FrozL.and.FrozM) then
       NAng=4
        write(IOut,'(I4,6X,''Cs  '',1X,2I4,11X,2I4,10X,2I4,14X,I2)')
     $    IAt,JAt,KAt,LAt,MAt,JAt,KAt,NAng
      Else
       NAng=5
        write(IOut,'(I4,6X,''Cs  '',1X,2I4,29X,4I4,5X,I2)')
     $    IAt,JAt,KAt,JAt,KAt,LAt,MAt,NAng
      EndIf
C Sciss1 symmetry coordinate
      If(.not.FrozJ.and..not.FrozK) then
       ICoord=ICoord+1
       NTerma(ICoord)=1
       ITVA(ICoord)=3
       CoefA(1,ICoord)=1.0d0
       IAtomA(1,1,ICoord)=JAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=KAt
      endIf
C Sciss2 symmetry coordinate
      if(.not.FrozL.and..not.FrozM) then
       ICoord=ICoord+1
       NTerma(ICoord)=1
       ITVA(ICoord)=3
       CoefA(1,ICoord)=1.0d0
       IAtomA(1,1,ICoord)=LAt
       IAtomA(2,1,ICoord)=IAt
       IAtomA(3,1,ICoord)=MAt
      endif
C Rocking symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=-1.0d0/DenRk
      CoefA(3,ICoord)=1.0d0/DenRk
      CoefA(4,ICoord)=-1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C Wagging symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=5
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=1.0d0/DenRk
      CoefA(3,ICoord)=-1.0d0/DenRk
      CoefA(4,ICoord)=-1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C Twisting symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=6
      CoefA(1,ICoord)=1.0d0/DenRk
      CoefA(2,ICoord)=-1.0d0/DenRk
      CoefA(3,ICoord)=-1.0d0/DenRk
      CoefA(4,ICoord)=1.0d0/DenRk
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
      return
      end 
*Deck D4h4At
      Subroutine D4h4At(Iout,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA,C)
      Implicit None
C IAt is the central atom              
C the angles J1,I1,K1 and L1,I1,M1 are of 90 degrees
C the angles J1,I1,M1 and K1,I1,L1 are of 90 degrees 
C the out-of-plane angles are managed in the GNICLA routine
C Dimensions
      Integer MxBnd, MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*),C(3,*)
C Local
      Integer I1,J1,K1,L1,M1
      Logical OKJK,OKJL,OKJM 
      Logical FrozJ,FrozK,FrozL,FrozM
      Real*8 Den,AJK,AJL,AJM,ValAng
C
      Den=2.0d0
      I1=IAt
      J1=JAt
      K1=KAt
      L1=LAt
      M1=MAt
      AJK=valang(C(1,J1),C(1,I1),C(1,K1))
      AJL=valang(C(1,J1),C(1,I1),C(1,L1))
      AJM=valang(C(1,J1),C(1,I1),C(1,M1))
      OKJK=AJK.eq.0.9d+2
      OKJL=AJL.eq.0.9d+2
      OKJM=AJM.eq.0.9d+2
      If(OKJK.and.OKJM) go to 10
      If(OKJK.and.OKJL) then
       L1=MAt
       M1=LAt
       go to 10
      endif
      If(OKJL.and.OKJM) then
       K1=LAt
       L1=KAt
      endif
   10 continue
      write(IOut,'('' Td4At:  Central ='',I5,'' Free   ='',4I5,
     $  30X,'' : Three angles in the Plane'')') IAt,JAt,KAt,
     $  LAt,MAt
C B1g symmetry coordinate
      ICoord=ICoord+1
      NTerma(ICoord)=4
      ITVA(ICoord)=12
      CoefA(1,ICoord)=1.0d0/Den
      CoefA(2,ICoord)=-1.0d0/Den
      CoefA(3,ICoord)=1.0d0/Den
      CoefA(4,ICoord)=-1.0d0/Den
      IAtomA(1,1,ICoord)=J1
      IAtomA(2,1,ICoord)=I1
      IAtomA(3,1,ICoord)=K1
      IAtomA(1,2,ICoord)=K1
      IAtomA(2,2,ICoord)=I1
      IAtomA(3,2,ICoord)=L1
      IAtomA(1,3,ICoord)=L1
      IAtomA(2,3,ICoord)=I1
      IAtomA(3,3,ICoord)=M1
      IAtomA(1,4,ICoord)=J1
      IAtomA(2,4,ICoord)=I1
      IAtomA(3,4,ICoord)=M1
C Eu1 symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=13
      CoefA(1,ICoord)=1.0d0/Den
      CoefA(2,ICoord)=-1.0d0/Den
      CoefA(3,ICoord)=-1.0d0/Den
      CoefA(4,ICoord)=1.0d0/Den
      IAtomA(1,1,ICoord)=J1
      IAtomA(2,1,ICoord)=I1
      IAtomA(3,1,ICoord)=K1
      IAtomA(1,2,ICoord)=K1
      IAtomA(2,2,ICoord)=I1
      IAtomA(3,2,ICoord)=L1
      IAtomA(1,3,ICoord)=L1
      IAtomA(2,3,ICoord)=I1
      IAtomA(3,3,ICoord)=M1
      IAtomA(1,4,ICoord)=J1
      IAtomA(2,4,ICoord)=I1
      IAtomA(3,4,ICoord)=M1
C Eu2 symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=13
      CoefA(1,ICoord)=1.0d0/Den
      CoefA(2,ICoord)=1.0d0/Den
      CoefA(3,ICoord)=-1.0d0/Den
      CoefA(4,ICoord)=-1.0d0/Den
      IAtomA(1,1,ICoord)=J1
      IAtomA(2,1,ICoord)=I1
      IAtomA(3,1,ICoord)=K1
      IAtomA(1,2,ICoord)=K1
      IAtomA(2,2,ICoord)=I1
      IAtomA(3,2,ICoord)=L1
      IAtomA(1,3,ICoord)=L1
      IAtomA(2,3,ICoord)=I1
      IAtomA(3,3,ICoord)=M1
      IAtomA(1,4,ICoord)=J1
      IAtomA(2,4,ICoord)=I1
      IAtomA(3,4,ICoord)=M1
      return
      end
*Deck Td4At
      Subroutine Td4At(Iout,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,CoefA)
      Implicit None
C               
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
C Input/Output
      Logical FrozJ,FrozK,FrozL,FrozM
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Real*8 DenEa,DenEb,DenT2
      write(IOut,'('' Td4At:  Central ='',I5,'' Free   ='',4I5,
     $  30X,'' : Five angles'')') IAt,JAt,KAt,LAt,MAt
C
      DenEa=Sqrt(12.0d0)
      DenEb=Sqrt(4.0d0)
      DenT2=Sqrt(2.0d0)
C Ea symmetry coordinate
      ICoord=ICoord+1
      NTerma(ICoord)=6
      ITVA(ICoord)=8
      CoefA(1,ICoord)=2.0d0/DenEa
      CoefA(2,ICoord)=-1.0d0/DenEa
      CoefA(3,ICoord)=-1.0d0/DenEa
      CoefA(4,ICoord)=-1.0d0/DenEa
      CoefA(5,ICoord)=-1.0d0/DenEa
      CoefA(6,ICoord)=2.0d0/DenEa 
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=KAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=LAt
      IAtomA(1,3,ICoord)=JAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=MAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=LAt
      IAtomA(1,5,ICoord)=KAt
      IAtomA(2,5,ICoord)=IAt
      IAtomA(3,5,ICoord)=MAt
      IAtomA(1,6,ICoord)=LAt
      IAtomA(2,6,ICoord)=IAt
      IAtomA(3,6,ICoord)=MAt
C Eb symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=4
      ITVA(ICoord)=8
      CoefA(1,ICoord)=1.0d0/DenEb
      CoefA(2,ICoord)=-1.0d0/DenEb
      CoefA(3,ICoord)=-1.0d0/DenEb
      CoefA(4,ICoord)=1.0d0/DenEb
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      IAtomA(1,3,ICoord)=KAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=LAt
      IAtomA(1,4,ICoord)=KAt
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=MAt
C T2x symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=9
      CoefA(1,ICoord)=1.0d0/DenT2
      CoefA(2,ICoord)=-1.0d0/DenT2
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=KAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
C T2y symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=10
      CoefA(1,ICoord)=1.0d0/DenT2
      CoefA(2,ICoord)=-1.0d0/DenT2
      IAtomA(1,1,ICoord)=KAt  
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt  
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
C T2z symmetry coordinate
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=11
      CoefA(1,ICoord)=1.0d0/DenT2
      CoefA(2,ICoord)=-1.0d0/DenT2
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=KAt
      IAtomA(1,2,ICoord)=LAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      return
      end
*Deck LooseS 
      Subroutine LooseS(IOut,IPrint,MxBnd,NAtoms,IAt,IAn,NBond,IBond,
     $  JAt,KAt,LAt,MAt)
      Implicit Integer (A-Z)
      Dimension IAN(*),NBond(*),IBond(MxBnd,*)
      Dimension ITst(MxBnd),ISV(MxBnd)
      do 10 i1=1,NAtoms
       ITst(i1)=IAn(IBond(i1,IAt))
       if(NBond(IBond(I1,IAt)).eq.1) itrm=i1
   10 continue
      ISmall=IrMin1(ITSt,NAtoms,.true.,ism)
      if(ISmall.eq.1) then
       ifnd=ism
      elseif(itrm.ne.0) then
       ifnd=itrm
      else
       JAt=0
       return
      endif
      do 20 i1=1,NAtoms 
       inew=i1+ifnd-1
       if(inew.gt.natoms) inew=inew-natoms
       ISv(i1)=IBond(inew,IAt)
   20 continue
       JAt=isv(1)
       KAt=isv(2)
       LAt=isv(3)
       If(NAtoms.eq.3) then
        MAt=0
       else
        MAt=isv(4)
       endif
       return
       end
*Deck Ord4At
      Subroutine Ord4At(IOut,JAt,KAt,LAt,MAt,NPivT,NEq,IAtCyc,Tresh,EAn)
      Implicit None
      Real*8 Tresh,EAn(*)
      Integer IOut,JAt,KAt,LAt,MAt,J1,K1,L1,M1,NPivT,NEq
      Integer IAtCyc(*)
      J1=JAt
      K1=KAt
      L1=LAt
      M1=MAt
      NEq=0
      If(Abs(EAN(J1)-EAN(K1)).lt.tresh) then
       neq=neq+1
       If(Abs(EAN(J1)-EAN(L1)).lt.tresh) then 
        neq=neq+2
        If(Abs(EAN(J1)-EAN(M1)).lt.tresh) neq=neq+1
       ElseIf(Abs(EAN(J1)-EAN(M1)).lt.tresh) then
        neq=neq+2
        LAt=M1
        MAt=L1
       ElseIf(Abs(EAN(L1)-EAN(M1)).lt.tresh) then
        neq=neq+1
       endif
      ElseIf(Abs(EAN(J1)-EAN(L1)).lt.tresh) then
        neq=neq+1
        JAt=J1
        KAt=L1
        LAt=K1
        MAt=M1 
       If(Abs(EAN(J1)-EAN(M1)).lt.tresh) then
        neq=neq+2
        JAt=J1
        KAt=L1
        LAt=M1
        MAt=K1
       ElseIf(Abs(EAN(K1)-EAN(M1)).lt.tresh) then
        neq=neq+1
        JAt=J1
        KAt=L1
        LAt=K1
        MAt=M1
       endif
      ElseIf(Abs(EAN(J1)-EAN(M1)).lt.tresh) then
       neq=neq+1
       JAt=J1
       KAt=M1
       LAt=K1
       MAt=L1
       If(Abs(EAN(K1)-EAN(L1)).lt.tresh) neq=neq+1
      ElseIf(Abs(EAN(K1)-EAN(L1)).lt.tresh) then
       neq=neq+1
       JAt=K1
       KAt=L1
       LAt=J1
       MAt=M1
       If(Abs(EAN(L1)-EAN(M1)).lt.tresh) then
        neq=neq+2
        JAt=K1
        KAt=L1
        LAt=M1
        MAt=J1
       endif
      ElseIf(Abs(EAN(K1)-EAN(M1)).lt.tresh) then
       neq=neq+1
       JAt=K1
       KAt=M1
       LAt=J1
       MAt=L1
      ElseIf(Abs(EAN(L1)-EAN(M1)).lt.tresh) then
       neq=neq+1
       JAt=L1
       KAt=M1
       LAt=J1
       MAt=K1
      endif
      if(NPivT.eq.0) return
      if(NPivT.eq.1) then
       If(IAtCyc(JAt).ne.0) return
       If(IAtCyc(KAt).ne.0) then
        J1=JAt
        K1=KAt
        JAt=K1
        KAt=J1
        Return
       ElseIf(IAtCyc(LAt).ne.0) then
        J1=JAt
        L1=LAt
        JAt=L1
        LAt=J1
        Return
       ElseIf(IAtCyc(MAt).ne.0) then 
        J1=JAt
        M1=MAt
        JAt=M1
        MAt=J1
        Return 
       EndIf
      Elseif(NPivT.eq.3) then
       If(IAtCyc(JAt).eq.0) return
       If(IAtCyc(KAt).eq.0) then
        J1=JAt
        K1=KAt
        JAt=K1
        KAt=J1
       ElseIf(IAtCyc(LAt).eq.0) then
        J1=JAt
        L1=LAt
        JAt=L1
        LAt=J1
       ElseIf(IAtCyc(MAt).eq.0) then
        J1=JAt
        M1=MAt
        JAt=M1
        MAt=J1
       EndIf
       If(Abs(EAn(KAt)-EAn(LAt)).lt.tresh) then
        If(Abs(EAn(LAt)-EAn(MAt)).gt.tresh) then
         K1=KAt
         M1=MAt
         KAt=M1
         MAt=K1
        EndIf 
       ElseIf(Abs(EAn(LAt)-EAn(MAt)).gt.tresh) then
        K1=KAt
        L1=LAt
        KAt=L1
        LAt=K1
       EndIf
      EndIf
      Return
      End
*Deck WXY3
      Subroutine WXY3(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,
     $  CoefA)
      Implicit None
C               
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
      Logical FrozJ,FrozK,FrozL,FrozM 
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NAng,NEq,NPivT
C All frozen Atoms
      If(NPivT.eq.4) then
       NAng=0
       write(IOut,'('' No Free Angles Around Center'',I5)')IAt
       Return
      EndIf
C At Most Three Frozen Atoms
C bending 2WXY1-WXY2-WXY3
      ICoord=ICoord+1
      NtermA(ICoord)=3
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/2.0d0
      CoefA(2,ICoord)=-1.0d0/4.0d0
      CoefA(3,ICoord)=CoefA(2,ICoord)
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=KAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=LAt
      IAtomA(1,3,ICoord)=JAt
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=MAt
C bending WXY2-WXY3
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/Sqrt(2.0d0)
      CoefA(2,ICoord)=-CoefA(1,ICoord)
      IAtomA(1,1,ICoord)=JAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
      IAtomA(1,2,ICoord)=JAt
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=MAt
      If(NPivT.eq.3) then
       If(FrozJ) then
        write(IOut,'('' Wrong Frozen Atom'',I5,'' Around'',I5)')JAt,IAt
        Stop
       EndIf
       NAng=2
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
       Return
      EndIf
C bending Y2XY3
      ICoord=ICoord+1
      ITVA(ICoord)=0
      NTermA(ICoord)=1
      CoefA(1,ICoord)=1.0D0
      IAtomA(1,1,ICoord)=LAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=MAt
      if(NPivT.eq.2) then
       NAng=3
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
       return
      endif
C At most One frozen atom
C bending Y1XY2 
      ICoord=ICoord+1
      ITVA(ICoord)=0
      NTermA(ICoord)=1
      CoefA(1,ICoord)=1.0D0
      IAtomA(1,1,ICoord)=KAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=LAt
C bending Y1XY3
      ICoord=ICoord+1 
      ITVA(ICoord)=0
      NTermA(ICoord)=1
      CoefA(1,ICoord)=1.0D0
      IAtomA(1,1,ICoord)=KAt
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=MAt
      NAng=5
      If(NPivT.eq.1) then
       If(.not.FrozJ) then
        write(IOut,'('' Wrong Free Atom'',I5,''  With Frozen'',3I5,
     $    '' Atoms Around'',I5)') JAt,KAt,LAt,MAt,IAt
        Stop
       EndIf
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      ElseIf(NPivT.eq.0) then
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C3v',NAng
      endif
      return
      end
*Deck W2XY2
      Subroutine W2XY2(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,NEq,NPivT,FrozJ,FrozK,FrozL,FrozM,NTermA,IAtomA,ITVA,
     $  CoefA)
      Implicit None
C               
C Dimensions
      Integer MaxAtA, MxTrmA
C Input
      Integer IOut,IPrint,ICoord,IAt,JAt,KAt,LAt,MAt
      Logical FrozJ,FrozK,FrozL,FrozM 
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NAng,NEq,NPivT,INot1,INot2,IYes1,IYes2
C All frozen Atoms
      If(NPivT.eq.4) then
       NAng=0
       write(IOut,'('' No Free Angles Around Center'',I5)')IAt
       Return
      EndIf
      NAng=4
C At Most Two Frozen Atoms
      Inot1=JAt
      Inot2=KAt
      Iyes1=LAt
      Iyes2=MAt
      If(FrozJ.and.FrozL) then
       Inot1=JAt
       Inot2=LAt
       Iyes1=KAt
       Iyes2=MAt
      elseIf(FrozJ.and.FrozM) then
       Inot1=JAt
       Inot2=MAt
       Iyes1=KAt
       Iyes2=LAt
      elseIf(FrozK.and.FrozL) then
       Inot1=KAt
       Inot2=LAt
       Iyes1=JAt
       Iyes2=MAt
      elseIf(FrozK.and.FrozM) then
       Inot1=KAt
       Inot2=MAt
       Iyes1=JAt
       Iyes2=LAt
      elseIf(FrozL.and.FrozM) then  
       Inot1=LAt
       Inot2=MAt
       Iyes1=JAt
       Iyes2=KAt
      endif
C bending 2Y1IY2-Y1IN1-Y1IN2
      ICoord=ICoord+1
      NtermA(ICoord)=3
      ITVA(ICoord)=1
      CoefA(1,ICoord)=2.0d0/SQrt(6.0d0)
      CoefA(2,ICoord)=-1.0d0/Sqrt(6.0d0)
      CoefA(3,ICoord)=CoefA(2,ICoord) 
      IAtomA(1,1,ICoord)=IYes1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IYes2
      IAtomA(1,2,ICoord)=INot1
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IYes1
      IAtomA(1,3,ICoord)=INot1
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=IYes2
C bending Y1IN1-Y1IN2
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/SQrt(2.0d0)
      CoefA(2,ICoord)=-1.0d0/SQrt(2.0d0)
      IAtomA(1,1,ICoord)=INot1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IYes1
      IAtomA(1,2,ICoord)=INot1
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IYes2
C bending 2Y1IY2-Y2IN1-Y2IN2
      ICoord=ICoord+1
      NtermA(ICoord)=3
      ITVA(ICoord)=1
      CoefA(1,ICoord)=2.0d0/SQrt(6.0d0)
      CoefA(2,ICoord)=-1.0d0/SQrt(6.0d0)
      CoefA(3,ICoord)=CoefA(2,ICoord)
      IAtomA(1,1,ICoord)=IYes1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IYes2
      IAtomA(1,2,ICoord)=INot2
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IYes1
      IAtomA(1,3,ICoord)=INot2
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=IYes2
C bending Y2IN1-Y2IN2
      ICoord=ICoord+1
      NtermA(ICoord)=2
      ITVA(ICoord)=2
      CoefA(1,ICoord)=1.0d0/SQrt(2.0d0)
      CoefA(2,ICoord)=-1.0d0/SQrt(2.0d0)
      IAtomA(1,1,ICoord)=INot2
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IYes1
      IAtomA(1,2,ICoord)=INot2
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IYes2
      If(FrozJ.or.FrozK.or.FrozL.or.FrozM) then
       NAng=4
       If(NEq.eq.2) then
        write(IOut,'(I4,6X,''C2v '',1X,2I4,'' +'',2I4,1X,2I4,10X,2I4,
     $    14X,I2)') IAt,JAt,KAt,LAt,MAt,INot1,INot2,IYes1,IYes2,NAng
       Else
        write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
       EndIf
       Return
      EndIf
C N1IN2 bending
      ICoord=ICoord+1
      NtermA(ICoord)=1
      ITVA(ICoord)=0
      CoefA(1,ICoord)=1.0d0
      IAtomA(1,1,ICoord)=INot1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=INot2
      NAng=5
      If(NEq.eq.2) then
       write(IOut,'(I4,6X,A3,8X,I2)') IAt,'C2v',NAng
      Else
       write(IOut,'(I4,6X,A2,9X,I2)') IAt,'Cs',NAng
      EndIf
      Return
      End
*Deck SpiAng
      Subroutine SpiAng(IOut,IPrint,MaxAtA,MxTrmA,ICoord,IAt,JAt,KAt,
     $  LAt,MAt,NTermA,IAtomA,ITVA,CoefA)
      Implicit Real*8 (A-H,O-Z)
      Common/onec/N1Cyc,IAt1C(2,100)
      Dimension NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Dimension CoefA(MxTrmA,*)
      Dimension ITst(4)
      Do 10 ic=1,N1Cyc
       If(IAt1C(1,IC).eq.JAt) ITst(1)=IAt1C(2,IC)
       If(IAt1C(1,IC).eq.KAt) ITst(2)=IAt1C(2,IC)
       If(IAt1C(1,IC).eq.LAt) ITst(3)=IAt1C(2,IC)
       If(IAt1C(1,IC).eq.MAt) ITst(4)=IAt1C(2,IC)
   10 Continue   
      IAng1=JAt
      If(ITst(2).eq.ITst(1)) then
       IAng2=KAt
       IAng3=LAt
       IAng4=MAt
      Else
       If(ITst(3).eq.ITst(1)) then
        IAng2=LAt
        IAng3=KAt
        IAng4=MAt 
       Else
        IAng2=MAt
        IAng3=KAt
        IAng4=LAt
       EndIf
      EndIf
      NAng=3
      write(IOut,'(I4,6X,''Spiro '',58X,I2)') IAt, NAng
      Cff=1.0/Sqrt(2.0d0)
C + + - - symmetry coordinate
      ICoord=ICoord+1
      ITVA(ICoord)=16
      NTerma(ICoord)=4
      CoefA(1,ICoord)=Cff
      CoefA(2,ICoord)=Cff 
      CoefA(3,ICoord)=-Cff
      CoefA(4,ICoord)=-Cff
      IAtomA(1,1,ICoord)=IAng1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IAng3
      IAtomA(1,2,ICoord)=IAng1
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IAng4
      IAtomA(1,3,ICoord)=IAng2
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=IAng3
      IAtomA(1,4,ICoord)=IAng2
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=IAng4
C + - - + symmetry coordinate
      ICoord=ICoord+1
      ITVA(ICoord)=16
      NTerma(ICoord)=4 
      CoefA(1,ICoord)=Cff
      CoefA(2,ICoord)=-Cff 
      CoefA(3,ICoord)=-Cff
      CoefA(4,ICoord)=Cff 
      IAtomA(1,1,ICoord)=IAng1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IAng3
      IAtomA(1,2,ICoord)=IAng1
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IAng4
      IAtomA(1,3,ICoord)=IAng2
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=IAng3
      IAtomA(1,4,ICoord)=IAng2
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=IAng4
C + - + - symmetry coordinate 
      ICoord=ICoord+1
      ITVA(ICoord)=16
      NTerma(ICoord)=4 
      CoefA(1,ICoord)=Cff
      CoefA(2,ICoord)=-Cff
      CoefA(3,ICoord)=Cff
      CoefA(4,ICoord)=-Cff  
      IAtomA(1,1,ICoord)=IAng1
      IAtomA(2,1,ICoord)=IAt
      IAtomA(3,1,ICoord)=IAng3
      IAtomA(1,2,ICoord)=IAng1
      IAtomA(2,2,ICoord)=IAt
      IAtomA(3,2,ICoord)=IAng4
      IAtomA(1,3,ICoord)=IAng2
      IAtomA(2,3,ICoord)=IAt
      IAtomA(3,3,ICoord)=IAng3
      IAtomA(1,4,ICoord)=IAng2
      IAtomA(2,4,ICoord)=IAt
      IAtomA(3,4,ICoord)=IAng4
      Return
      End
