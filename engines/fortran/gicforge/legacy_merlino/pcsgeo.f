*Deck PCSGeo
      Subroutine PCSGeo(IOut,IPrint,Linear,DoneC,MxAtP,MxTrm,NAtoms,
     $  NLenR,NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,IAtmLR,
     $  IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,C,Atmass,R0IJ,
     $  RBIJ,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr)
      Implicit Real*8 (A-H,O-Z)
      Logical NoG,DoB1,DoMW,Linear,DoneC
      Common/PhyCon/ToAng,PhyCon(29)
      Dimension C(3,*),AtMass(*),R0IJ(*),RBIJ(*),Scr(*)
      Dimension IAtmBR(MxAtP,MxTrm,*),IAtmAR(MxAtP,MxTrm,*)
      Dimension IAtmLR(MxAtP,MxTrm,*),IAtmDR(MxAtP,MxTrm,*)
      Dimension IAtmOR(MxAtP,MxTrm,*),IAn(*),IScr(*)
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)
      NoG=.False.
      DoB1=.False.
      DoMW=.False.
      IQ0=1
C Scr(IQ0)=Q0(NTotR)
      ISQ=IQ0+NTotR
C SCr(ISQ)=Sq=qfin-q0
      IQK=ISQ+NTotR 
C SCr(IQK)=Qk(NTotR) int. coords. at iteration K
      IDltQ=IQK+NTotR
C SCr(IDltQ)=Deltaq(NTotR) at iteration k 
      ICOld=IDltQ+NTotR 
C Scr(ICOld)=Cartesian Coordinates at iteration k (3*NAtoms)
      ICNew=ICOld+3*NAtoms 
C SCr(ICNew)=Cartesian Coordinates at iteration k+1 
      IniB=ICNew+3*NAtoms
      IniG=IniB+3*NAtoms*NTotR
      IniGm1=IniG+NTotR*NTotR 
      InBtm1=IniGm1+NTotR*NTotR
      IEnd=InBTm1+3*NAtoms*NTotR
C Scr(IniB) = B Matrix, Scr(IniG) = G Matrix, Scr(IniGm1) = G-1 Matrix
C Scr(InBtm1) = (B+)-1 Matrix all made in MkGB routine
C Initialize Cartesian Coordinates
      Call AMove(3*NAtoms,C,Scr(ICOld))
      Call AMove(3*NAtoms,Scr(ICOld),Scr(ICNew))
C Initialize Internal Coordinates 
      Call MkInt(IOut,IPrint,NAtoms,MxAtP,MxTrm,NTotR,NLenR,NAngR,
     $  NLAngR,NDihR,NOuplR,IAtmBR,IAtmAR,IAtmLR,IAtmDR,IAtmOR,
     $  Scr(ICOld),Scr(IQ0))
C Set Sq vector for BPCS
      call AClear(NTotR,Scr(ISQ))
      do 10 ir=1,NLenR
       Scr(ISQ+ir-1)=RBIJ(ir)-R0IJ(ir)
   10 Continue
C Initialize DeltaQ
      Call Aclear(NTotR,Scr(IDltQ)) 
      NCMax=20
      NCyc=0
   20 Continue
      If(NCyc.gt.NCMax) then
       write(IOut,'(/,'' Conversion from Internal to Cartesian Not''
     $   '' Converged After'',I3,'' Iterations'')') NCMax
       Stop
      EndIf
      NCyc=NCyc+1
      Call AMove(3*NAtoms,Scr(ICNew),Scr(ICOld)) 
      Call MkBG(IOut,IPrint,Linear,NoG,DoB1,DoMW,MxAtP,MxTrm,
     $  NAtoms,NLenR,NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,
     $  IAtmLR,IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,
     $  Scr(ICNew),Atmass,CoefB,CoefA,CoefL,CoefD,CoefO,Scr(IniB),Iscr)
C Compute new Cartesians: x(K+1)=x(k)+[(B+)-1(k)[sq-(qk-1-q0)]
      Call MkCart(IOut,IPrint,NAtoms,NTotR,DoneC,DltCMx,DltCAv,IAn,
     $  Scr(ICOld),Scr(ICNew),Scr(ISQ),Scr(IDltQ),Scr(InBtm1))
C Compute New Internal Coordinates
      Call MkInt(IOut,IPrint,NAtoms,MxAtP,MxTrm,NTotR,NLenR,NAngR,
     $  NLAngR,NDihR,NOuplR,IAtmBR,IAtmAR,IAtmLR,IAtmDR,IAtmOR,
     $  Scr(ICNew),Scr(IQK))
      Do 30 Iint=1,NTotR 
       Scr(IDltQ+iint-1)=Scr(IQK+iint-1)-Scr(IQ0+iint-1)
   30 continue
      If(.not.DoneC) go to 20
      Write(IOut,'('' MkCart converged after'',I3,'' Cycles'')') NCyc
      Write(IOut,'('' Max.Diff. ='',D10.3,''; Aver.Diff. ='',D10.3)')
     $  DltCMx,DltCAv
      Call AMove(3*NAtoms,Scr(ICNew),C)
      Return
      End

*Deck MkInt
      Subroutine MkInt(IOut,IPrint,NAtoms,MxAtP,MxTrm,NTot,NLen,NAng,
     $  NLAng,NDih,NOuPl,IAtmB,IAtmA,IAtmL,IAtmD,IAtmO,C,R)
      Implicit Real*8 (A-H,O-Z)
      Dimension C(3,*),R(*)
      Dimension IAtmB(MxAtP,MxTrm,*),IAtmA(MxAtP,MxTrm,*)
      Dimension IAtmL(MxAtP,MxTrm,*),IAtmD(MxAtP,MxTrm,*)
      Dimension IAtmO(MxAtP,MxTrm,*)
C Angles internally in radiants
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      TrAv=1.0D-06
      TrMax=1.0D-05
      IntMax=0
      DltIMx=0.0D0
      DlTot=0.0D0
      Do 10 ILen=1,NLen
       IRed=ILen
       IAt=IAtmB(1,1,ILen)
       JAt=IAtmB(2,1,ILen)
       R(IRed)=Distan(C,IAt,JAt,0)  
       If(IPrint.gt.0) Write(IOut,'(I3,'' Bond Length   '',2I5,10X,
     $   F12.5)') IRed,IAt,JAt,R(IRed)
   10 Continue
      If(NAng.gt.0) then
       Do 20 IAng=1,NAng
        IRed=NLen+IAng
        IAt=IAtmA(1,1,IAng)
        JAt=IAtmA(2,1,IAng)
        KAt=IAtmA(3,1,IAng)
        R(IRed)=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        If(IPrint.gt.0) Write(IOut,'(I3,'' Valence Angle '',3I5,5X,
     $    F12.5)') IRed,IAt,JAt,KAt,R(IRed)*ToDeg
   20  Continue
      EndIf
      If(NLAng.gt.0) then
       Do 30 ILAng=1,NLAng
        IRed=NLen+NAng+ILAng
        IAt=IAtmL(1,1,ILAng)
        JAt=IAtmL(2,1,ILAng)
        KAt=IAtmL(3,1,ILAng)
        R(IRed)=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        If(IPrint.gt.0) Write(IOut,'(I3,'' Linear Angle  '',3I5,5X,
     $    F12.5)') IRed,IAt,JAt,KAt,R(IRed)*ToDeg
   30  Continue
      EndIf
      If(NDih.gt.0) then
       Do 40 IDih=1,NDih
        IRed=NLen+NAng+NLAng+IDih
        IAt=IAtmD(1,1,IDih)
        JAt=IAtmD(2,1,IDih)
        KAt=IAtmD(3,1,IDih)
        LAt=IAtmD(4,1,IDih)
        R(IRed)=Dihed(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))
        If(IPrint.gt.0) Write(IOut,'(I3,'' Dihedral Angle'',4I5,
     $    F12.5)') IRed,IAt,JAt,KAt,LAt,R(IRed)*ToDeg
   40  Continue
      EndIf
      If(NOupl.gt.0) then
       Do 50 IOupl=1,NOupl
        IRed=NLen+NAng+NLAng+NDih+IOupl
        IAt=IAtmO(1,1,IOupl)
        JAt=IAtmO(2,1,IOupl)
        KAt=IAtmO(3,1,IOupl)
        LAt=IAtmO(4,1,IOupl)
        R(IRed)=OutAng(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))
        If(IPrint.gt.0) Write(IOut,'(I3,'' Out Pl. Angle '',4I5,
     $    F12.5)') IRed,IAt,JAt,KAt,LAt,R(IRed)*ToDeg
   50  Continue
      EndIf
      Return
      End 
*Deck MkCart 
      Subroutine MkCart(IOut,IPrint,NAtoms,NGIC,DoneC,DiffMx,DiffAv,IAn,
     $  COld,CNew,SQ,DQk,G1BMat)
      Implicit Real*8 (A-H,O-Z)
      Logical DoneC
      Dimension IAn(*)
      Dimension COld(3,*),CNew(3,*),SQ(*),DQk(*),G1BMat(3*NAtoms,*)
C Angles in radiants
      DoneC=.False.
      DiffT=0.0d0
      DiffMx=0.0d0
      TrAv=1.0d-9
      TrMx=1.0d-8
      Do 10 IAt=1,NAtoms
       Do 20 ICoord=1,3
        CNew(ICoord,IAt)=COld(ICoord,IAt)
        IXYZ=3*(IAt-1)+ICoord
        Do 30 Iint=1,NGIC
         CNew(ICoord,IAt)=CNew(ICoord,IAt)+G1BMat(IXYZ,Iint)*
     $     (SQ(Iint)-DQk(Iint))
   30   continue
        Diff=Abs(CNew(ICoord,IAt)-COld(ICoord,IAt))
        If(Diff.gt.DiffMx) DiffMx=Diff
        DiffT=DiffT+Diff
   20  continue
   10 continue
      DiffAv=DiffT/Float(3*NAtoms)
      If(DiffAv.lt.TrAv.and.DiffMx.lt.TrMx) DoneC=.True.  
      If(IPrint.eq.0) Return 
      Write(IOut,'('' Maximum Difference:'',D10.5,5X,
     $  '' Average Difference:'',D10.5)') DiffMx,DiffAv
      If(DiffAv.lt.TrAv.and.DiffMx.lt.TrMx) then
       Write(IOut,'('' Conversion Converged'')')
       Write(IOut,'(/,'' New Cartesian Coordinates'')')
       Write(IOut,'(3X,''Atom'',2X,''At.Numb.'',5X, ''X'',10X,''Y'',
     $  12X,''Z'')')
       Write(IOut,'(I5,I7,2X,3F12.5)')(IAt,IAn(IAt),
     $    (CNew(IXYZ,IAt),IXYZ=1,3),IAt=1,NAtoms)
       Write(IOut,'('' '')')
      EndIf
      Return
      End
*Deck MkBG 
      Subroutine MkBG(IOut,IPrint,Linear,NoG,DoB1,DoMW,MxAtP,MxTrm,
     $  NAtoms,NLen,NAng,NLAng,NOupl,NDih,NTot,IAn,IAtmB,IAtmA,
     $  IAtmL,IAtmD,IAtmO,NTermB,NTermA,NTermL,NTermD,NTermO,C,
     $  Atmass,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr)
      Implicit Real*8 (A-H,O-Z)
      Logical Linear,NoG,DoB1,DoMW
      Common/PhyCon/ToAng,PhyCon(29)
      Dimension C(3,*),AtMass(*),Scr(*)
      Dimension IAtmB(MxAtP,MxTrm,*),IAtmA(MxAtP,MxTrm,*)
      Dimension IAtmL(MxAtP,MxTrm,*),IAtmD(MxAtP,MxTrm,*)
      Dimension IAtmO(MxAtP,MxTrm,*),IAn(*),IScr(*)
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)
      ToBohr = 1.0d0/ToAng
      IniB=1 
C SCr(IniB)=BMat(3*NAtoms,NTot)
      IniG=IniB+3*NAtoms*NTot
C Scr(IniG)=GMat(NTot,NTot)
      IniGm1=IniG+NTot*NTot
C Scr(IniGm1)=G-1(NTot,NTot)
      InBTm1=IniGm1+NTot*NTot
C Scr(InBTm1)=(BB+)-1B(3*NAtoms,NTot)
      IEnd=InBTm1+3*NAtoms*NTot
C
C Here both B and G are in Angstroms whereas they are in Bohr in Gaussian 
C There is no difference for stretchings, but there is a factor 0.529177 for 
C angles
C Compute B Matrix      
      DoB1=.true.
      Call MkBNew(IOut,IPrint,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,
     $  NLAng,NOupl,NDih,IAtmB,IAtmA,IAtmL,IAtmD,IAtmO,NTermB,
     $  NTermA,NTermL,NTermD,NTermO,CoefB,CoefA,CoefL,CoefD,CoefO,
     $  C,Scr(IniB),.False.)
      If(NoG) Return
C Compute BB+ (DoMW=.false.) or G (DoMW=.true.) Matrix
      Call MakeG(IOut,IPrint,DoMW,NAtoms,NTot,AtMass,Scr(IniB),
     $  Scr(IniG))
C Compute (BB+)-1 or G-1
      NTInt=3*NAtoms-6
      If(Linear) NTInt=NTInt+1
      NRed=NTot-NTInt
      Call Amove(NTot*NTot,Scr(IniG),Scr(IniGM1))
      Call MakGm1(IOut,IPrint,DoMW,NTot,NRed,Scr(IniGm1),Scr(InBTm1),
     $  IScr)
C Compute (B+)-1 = (BB+)-1B (3*NAtoms columns)
      Call AClear(3*NAtoms*NTot,Scr(InBTm1))
      Call MkGm1B(IOut,IPrint,NAtoms,NTot,Scr(IniB),Scr(IniGm1),
     $  Scr(InBTm1))
      Return
      End
