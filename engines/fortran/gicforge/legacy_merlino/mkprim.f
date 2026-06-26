*Deck FndFrg 
      Subroutine FndFrg(MxBnd,NAtoms,IBond,NBond,IFrg,IV)
      Implicit Integer(A-Z)
C
C     Internal Coordinates Fragments Builder
C     ====================================================
C
C     Description
C     -----------
C     Uses the information regarding the connectivity to find out if
C     separate fragments are present in the molecule
C
C     Input
C     -----
C     MxBnd  :: Maximum atomic valence of the molecule
C     IBond  :: (MxBnd,NAtoms) Indexes of the atoms bonded to each atom
C     NBond  :: (NAtoms) Number of atoms bonded to each atom
C
C     Output
C     ------
C     IFrg   :: (NAtoms) Index of the fragment containing each atom
C
CEND
C
C     Dimensions
      Integer MDV, MxBnd, NAtoms
C     Input
      Integer IBond(MxBnd,*), NBond(*), InToWP
C     Output
      Integer IFrg(*)
C     Local
      Integer IV(*), i, ia, ja, ka
      Logical Check
C
      Call IClear(NAtoms,IFrg)
      Call IClear(NAtoms,IV)
C
C     Finds the first atom of the fragment
      Do 100 ia = 1, NAtoms
C       First loop over all the atoms, used to identify the first atom
C       of the fragment
        Check = .False.
        ja = 1
  110   If(IFrg(ja).eq.0) then
          Check    = .True.
          IFrg(ja) = ia
          IV(ja)   = ia
        else
          ja = ja + 1
          If(ja.le.NAtoms) Goto 110
          endIf
        If(.not.Check) Goto 10
C       The first atom of the fragment has been found, assign the flags
C       to the atoms bonded to this first atom
        Do 111 ka = 1, NBond(ja)
         IFrg(IBond(ka,ja)) = ia
  111   Continue
C       Starts the general loop over all the atoms. This loop is
C       repeated NAtoms time, to be sure that all the atoms have been
C       tested.
        Do 112 ja = 1, NAtoms
C         "Real" loop.
          Do 120 ka = 1, NAtoms
            If(IFrg(ka).eq.ia.and.IV(ka).eq.0) then
              IV(ka) = ia
              Do 130 i = 1, NBond(ka)
               IFrg(IBond(i,ka)) = ia
  130         Continue
              endIf
  120       Continue
  112     Continue
  100   Continue
   10 Continue
      call IMove(NAtoms,IV,IFrg) 
      Return
      End
*Deck MkBAL
      Subroutine MkBAL(IOut,IPrint,MxBnd,MxAtP,MxTrm,NAtoms,ILen,IAng,
     $ ILAng,IAn,NBond,IBond,NTermB,NTermA,NTermL,IAtomB,IAtomA,IAtomL,
     $ CoefB,CoefA,CoefL,C,TreshL)
      Implicit Real*8 (A-H,O-Z)
C The argument of cos is in radiants
C The initial values of ILen,IAng,ILAng are 0
      Dimension IAn(*),NBond(*),IBond(MxBnd,*)
      Dimension IAtomB(MxAtP,MxTrm,*),NTermB(*)
      Dimension IAtomA(MxAtP,MxTrm,*),NTermA(*)
      Dimension IAtomL(MxAtP,MxTrm,*),NTermL(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*)
      Dimension CoefL(MxTrm,*),C(3,*)
C Build Stretchings
      Do 10 JAt=1,NAtoms
       Do 20 IB=1,NBond(JAt)
        IAt=IBond(IB,JAt)
        if(NBond(JAt).eq.1) go to 40 
C Builds valence angles 
C (including linear ones for values larger than treshl in radiants)
        Do 30 KB=IB,NBond(JAt)
         KAt=IBond(KB,JAt)
         if(IAt.eq.KAt) go to 30
C Val1 is in radiants
         Val1=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
         If(Val1.lt.treshl) then
          IAng=IAng+1
          NTermA(IAng)=1
          CoefA(1,IAng)=1.0d0 
          IAtomA(2,1,IAng)=JAt 
          IAtomA(1,1,IAng)=Min0(IAt,KAt)
          IAtomA(3,1,IAng)=Max0(IAt,KAt)
         else
          ILAng=ILAng+1
          NTermL(ILAng)=1
          CoefL(1,ILAng)=1.0d0
          IAtomL(2,1,ILAng)=JAt
          IAtomL(1,1,ILAng)=Min0(IAt,KAt)
          IAtomL(3,1,ILAng)=Max0(IAt,KAt)
          IAtomL(4,1,ILAng)=-1
          ILAng=ILAng+1
          NTermL(ILAng)=1
          CoefL(1,ILAng)=1.0d0
          IAtomL(2,1,ILAng)=JAt
          IAtomL(1,1,ILAng)=Min0(IAt,KAt)
          IAtomL(3,1,ILAng)=Max0(IAt,KAt)
          IAtomL(4,1,ILAng)=-2
         endif
   30   continue
   40   If(IAt.lt.JAt) go to 20
        ILen=ILen+1
        NTermB(ILen)=1
        IAtomB(1,1,ILen)=JAt
        IAtomB(2,1,ILen)=IAt
        CoefB(1,ILen)=1.0d0
   20  continue
   10 continue
      return
      end
*Deck Distan
      Function Distan(C,iatom,jatom,IScale)
      Implicit Real*8 (A-H,O-Z)
      Dimension C(3,*),Ctemp(3)
      Common/PHYCON/ToAng,PhyCon(29)
C
C  Compute the distance between iatom and jatom
C  IScale = 0  Don't scale
C         = 1  Convert Bohr to Angs
C         = 2  Convert Angs to Bohr
C
      Call ASub(3,C(1,iatom),C(1,jatom),Ctemp)
      Distan = Sqrt(sprod(3,Ctemp,Ctemp))
      If(IScale.eq.1) then
        Distan = ToAng*Distan
      else if(IScale.eq.2) then
        Distan = Distan/ToAng
        endIf
      Return
      End
*Deck ValAng 
      Function ValAng(CI,CJ,CK)
      Implicit Real*8(A-H,O-Z)
C     
C     Calculation of the value of the valence angle in RADIANTS
C     
      Dimension CI(3),CJ(3),CK(3)
      Save Zero, One, Small
      Data Zero/0.0d0/, One/1.d0/, Small/1.0d-20/
C     
      DJI = Zero
      DJK = Zero
      DotJ = Zero
      Do 10 M = 1, 3
        DJIM = CI(M) - CJ(M)
        DJKM = CK(M) - CJ(M)
        DotJ = DotJ + DJIM*DJKM
        DJI = DJI + DJIM**2
        DJK = DJK + DJKM**2
   10 continue
      If(Abs(DotJ).ge.Small)
     $  DotJ = DotJ / Sqrt(DJI*DJK)
      DotJ = Max(Min(DotJ,One),-One)
      ValAng = ACos(DotJ)
      Return
      End
*Deck MkPrmD 
      Subroutine MkPrmD(IOut,IPrint,MxBnd,MxTrmB,MxTrmD,MxAtB,MxAtD,
     $  MxAtCy,MxCyc,NAtoms,NBond,NLen,NDih,NTot,NCyc,IBond,NTermD,
     $  IAtomB,IAtomD,NAtC,ICAt,IAtCyc,CoefD,C,TreshL)
      Implicit Real*8 (A-H,O-Z)
      Logical Fnd3C,Fnd4C
      Dimension NBond(*),IBond(MxBnd,*)
      Dimension NTermD(*),IAtomD(MxAtD,MxTrmD,*)
      Dimension IAtomB(MxAtB,MxTrmB,*)
      Dimension NatC(*),ICAt(MxAtCy,*),IAtCyc(*)
      Dimension CoefD(MxTrmD,*),C(3,*)
      pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
      NExpCy=NLen-NAtoms+1
      NDih=0
C Build Primitive Dihedrals
      Do 10 ILen=1,NLen 
       JAt=IAtomB(1,1,ILen)
       KAt=IAtomB(2,1,ILen)
       NBJ=NBond(JAt)
       NBK=NBond(KAt)
       If(NBJ.eq.1.or.NBK.eq.1) go to 10
       Do 20 jj=1,NBJ
        IAt=IBond(jj,JAt)
        NBI=NBond(IAt)
        if(IAt.eq.KAt) go to 20
        Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        If(Value.gt.TreshL) go to 20
        Do 30 kk=1,NBK
         LAt=IBond(kk,KAt)
         If(LAt.eq.JAt) go to 30
         Value=ValAng(C(1,JAt),C(1,KAt),C(1,LAt))
         If(Value.gt.TreshL) go to 30 
         If(LAt.eq.IAt) then
C Build 3-membered cycles
          Fnd3C=.false.
          if(NCyc.gt.0) call Cyc3At(MxAtCy,IAt,JAt,KAt,NCyc,NAtC,ICAt,
     $      Fnd3C)
          if(.not.Fnd3c) then
           NCyc=NCyc+1
           NAtC(NCyc)=3
           ICAt(1,NCyc)=IAt
           ICAt(2,NCyc)=JAt
           ICAt(3,NCyc)=KAt
           IAtCyc(IAt)=NCyc
           IAtCyc(JAt)=NCyc
           IAtCyc(KAt)=NCyc
          endif
          go to 30
         EndIf
         NDih=NDih+1 
         NTermD(NDih)=1
         IAtomD(1,1,NDih)=IAt
         IAtomD(2,1,NDih)=JAt
         IAtomD(3,1,NDih)=KAt
         IAtomD(4,1,NDih)=LAt
         CoefD(1,NDih)=1.0d0
C Build 4-membered cycles
C        If(NCyc.lt.NExpCy) then
          Fnd4C=.false.
          if(NCyc.gt.0) call Cyc4At(MxAtCy,IAt,JAt,KAt,LAt,NCyc,NAtC,
     $      ICAt,Fnd4C)
          If(.not.Fnd4C) then 
           do 40 ii=1,NBI          
            MAt=IBond(ii,IAt)
            If(MAt.eq.LAt) then
             NCyc=NCyc+1
             NAtC(NCyc)=4
             ICAt(1,NCyc)=IAt
             ICAt(2,NCyc)=JAt
             ICAt(3,NCyc)=KAt
             ICAt(4,NCyc)=LAt
             IAtCyc(IAt)=NCyc
             IAtCyc(JAt)=NCyc
             IAtCyc(KAt)=NCyc
             IAtCyc(LAt)=NCyc
            endif
   40      continue
          EndIf
C        EndIf
         If(IPrint.gt.0) then
          If(NDih.eq.1) write(IOut,'(/,'' Dihedrals'')')
          Value=Dihed(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))*ToDeg
          write (IOUT,'(I5,''  D('',3(I3,'',''),I3,'')'',2X,''Value: '',
     $      F9.3)') NTot+NDih,IAt,JAt,KAt,LAt,Value
         endif
   30   continue
   20  continue
   10 continue 
      if(NCyc.gt.0.and.Iprint.gt.0) then
       do 50 icyc=1,NCyc
        write(IOut,'('' Ring Number'',I3,'' Contains'',I2,'' Atoms:'',
     $    4I4)') ICyc,NAtC(ICyc),(ICAt(ii,ICyc),ii=1,NAtc(ICyc))
   50  continue
      endif
      return
      End
*Deck MkPrmDOld 
      Subroutine MkPrmDOld(IOut,IPrint,MxBnd,MxTrmB,MxTrmD,MxAtB,MxAtD,
     $  NAtoms,NBond,NLen,NDih,NTot,IBond,NTermD,IAtomB,IAtomD,CoefD,
     $  C,TreshL)
      Implicit Real*8 (A-H,O-Z)
      Dimension NBond(*),IBond(MxBnd,*)
      Dimension NTermD(*),IAtomD(MxAtD,MxTrmD,*)
      Dimension IAtomB(MxAtB,MxTrmB,*)
      Dimension CoefD(MxTrmD,*),C(3,*)
      pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
      NDih=0
C Build Primitive Dihedrals
      Do 10 ILen=1,NLen 
       JAt=IAtomB(1,1,ILen)
       KAt=IAtomB(2,1,ILen)
       NBJ=NBond(JAt)
       NBK=NBond(KAt)
       If(NBJ.eq.1.or.NBK.eq.1) go to 10
       Do 20 jj=1,NBJ
        IAt=IBond(jj,JAt)
        if(IAt.eq.KAt) go to 20
        Value=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        If(Value.gt.TreshL) go to 20
        Do 30 kk=1,NBK
         LAt=IBond(kk,KAt)
         If(LAt.eq.JAt) go to 30
         Value=ValAng(C(1,JAt),C(1,KAt),C(1,LAt))
         If(Value.gt.TreshL) go to 30 
         NDih=NDih+1 
         NTermD(NDih)=1
         IAtomD(1,1,NDih)=IAt
         IAtomD(2,1,NDih)=JAt
         IAtomD(3,1,NDih)=KAt
         IAtomD(4,1,NDih)=LAt
         CoefD(1,NDih)=1.0d0
         If(IPrint.gt.0) then
          If(NDih.eq.1) write(IOut,'(/,'' Dihedrals'')')
          Value=Dihed(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))*ToDeg
          write (IOUT,'(I5,''  D('',3(I3,'',''),I3,'')'',2X,''Value: '',
     $      F9.3)') NTot+NDih,IAt,JAt,KAt,LAt,Value
         endif
   30   continue
   20  continue
   10 continue 
      return
      End
*Deck MkPrmO
      Subroutine MkPrmO(IOut,IPrint,MxBnd,MxTrmO,MxAtO,NAtoms,NBond,
     $  NOupl,IBond,NTrmO,IAtomO,CoefO,C,DoG16,IAtCyc)
      Implicit Real*8 (A-H,O-Z)
      Logical DoG16
      Dimension NBond(*),IBond(MxBnd,*),IAtCyc(*)
      Dimension NTrmO(*),IAtomO(MxAtO,MxTrmO,*) 
      Dimension CoefO(MxTrmO,*),C(3,*)
      pi = acos(-1.d0)
      ToDeg=1.80d+2/pi
      NOuPl=0
C Build Primitive Out-of-Plane Angles
      Do 10 IAt=1,NAtoms
       If(NBond(IAt).ne.3) go to 10
       JAt=IBond(1,IAt)
       KAt=IBond(2,IAt)
       LAt=IBond(3,IAt)
       If(IAtCyc(IAt).gt.0.and.IAtCyc(JAt).gt.0.and.
     $    IAtCyc(KAt).gt.0.and.IAtCyc(LAt).gt.0) go to 10
       NOupl=NOupl+1
       NTrmO(NOupl)=1
       IAtomO(1,1,NOupl)=IAt
       IAtomO(2,1,NOupl)=JAt
       IAtomO(3,1,NOupl)=KAt
       IAtomO(4,1,NOupl)=LAt
       CoefO(1,NOupl)=1.0d0
       If(IPrint.gt.0) then 
        If(NOuPl.eq.1) write(IOut,'(/,'' Out of Plane Angles'')')
        Value0=OutAngOld(C(1,JAt),C(1,IAt),C(1,KAt),C(1,LAt))*ToDeg
        Value=Outang(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))*ToDeg
        If(DoG16) then 
         ValG16=Dihed(C(1,JAt),C(1,IAt),C(1,LAt),C(1,KAt))*ToDeg
         write (IOUT,'(I5,''  D('',3(I3,'',''),I3,'')'',2X,
     $    ''Value: '',F9.3)')NOupl,JAt,IAt,LAt,KAt,ValG16 
        Else
         Value=Outang(C(1,IAt),C(1,JAt),C(1,KAt),C(1,LAt))*ToDeg
         write (IOUT,'(I5,''  D('',3(I3,'',''),I3,'')'',2X,
     $    ''Value: '',F9.3)')NOupl,IAt,JAt,KAt,LAt,Value  
        EndIf
       EndIf
   10 continue  
      return
      End
*Deck OutAngOld
      Real*8 function OutAngOld(CI,CJ,CK,CL)
C
C Computes the value of an out-of-plane angle
C i:1st, j:2nd, k:3rd, l:4th)
C the fourth atom out of the plane of the first 3
C the second atom is the central atom
C
      implicit none
      integer natoms, iat, jat, kat, lat, ixyz
      real*8  rij(3), dij, eij(3), rkj(3), dkj, ekj(3)
      real*8  rlj(3), dlj, elj(3), rmj(3), dmj, emj(3), sint, sprod
      real*8  ci(3),  cj(3), ck(3), cl(3)
C Compute distance vectors:
      do 10 ixyz=1,3
       rij(ixyz)=ci(ixyz)-cj(ixyz)
       rkj(ixyz)=ck(ixyz)-cj(ixyz)
       rlj(ixyz)=cl(ixyz)-cj(ixyz)
 10   continue
C Normalize the distance vectors
      dij = sqrt(sprod(3,rij,rij))
      eij = rij/dij
      dkj = sqrt(sprod(3,rkj,rkj))
      ekj = rkj/dkj
      dlj = sqrt(sprod(3,rlj,rlj))
      elj = rlj/dlj
C Compute the out-of-plane coordinate
      call vprod(rmj,rij,rkj)
      dmj = sqrt(sprod(3,rmj,rmj))
      emj = rmj/dmj
      sint = sprod(3,emj,elj)
      OutAngOld = asin(sint)
      return
      end
*Deck Dihed
      Function Dihed(CI,CJ,CK,CL)
      Implicit Real*8(A-H,O-Z)
C        
C     Return the dihedral angle between the planes defined by I-J and
C     K-L.  Compute the normalized vectors A=LAt-KAt, B=IAt-JAt, and
C     D=KAt-JAt.  Then compute the normals PIJK (perpendicular to the
C     IAt.JAt.KAt plane) and PJKL (perpendicular to the JAt.KAt.LAt
C     plane).  The angle between PIJK and PJKL is the dihedral angle.
C     If BxA.D is negative, the dihedral angle is also.  If the bond
C     angle is 180, the dihedral angle is undefined and is set to zero.
C        
      Dimension CI(3), CJ(3), CK(3), CL(3), A(3), B(3), D(3), PIJK(3),
     $  PJKL(3), BXA(3)
      Save Zero, One, Four, Tol
      Data Zero/0.0d0/, One/1.0d0/, Four/4.0d0/, Tol/1.0d-6/
C        
      Pi = Four*ATan(One)
      Call ASUnit(3,CL,CK,A)
      Call ASUnit(3,CI,CJ,B)
      Call ASUnit(3,CK,CJ,D)
      Call VProd(PIJK,B,D)
      Call VProd(PJKL,A,D)
      Call VProd(BXA,B,A)
      PDotP = SProd(3,PIJK,PJKL)
      ALen = Sqrt(SProd(3,PIJK,PIJK))
      BLen = Sqrt(SProd(3,PJKL,PJKL))
      ABLen = ALen*BLen
      If(Abs(ABLen).lt.Tol) then
        Ang = Zero
      else
        BADotD = SProd(3,BXA,D)
        Ang = Sign(ACos1(PDotP/(ALen*BLen)),BADotD)
        If(Abs(Ang).lt.Tol) Ang = Zero
        If(Abs(Abs(Ang)-Pi).lt.Tol) Ang = Pi
        endIf
      Dihed = Ang
      Return
      End
*Deck PrtPrm
      Subroutine PrtPrm(IOut,MaxAtG,MaxTer,NGic,Itype,IAtomG,NPV)
      Implicit Integer (A-H,O-Z)
      Character CT*1,CType*5
      Dimension NAtG(5),IAtomG(MaxAtG,MaxTer,*)
      Data CType/'BALDO'/
      Data NatG/2,3,4,4,4/
      If(NGic.lt.1) return
      If(IType.lt.1.or.IType.gt.5) then
       write(IOut,'('' From PrtPrm: wrong type of primitive'')')
       stop
      endif
      NAT=NAtG(Itype)
      CT=CType(IType:IType)
      do 10 IGic=1,NGic
       NPV=NPV+1
       I1=IAtomG(1,1,IGic)
       J1=IAtomG(2,1,IGic)
       If(Itype.eq.1) then
        Write(IOut,'(1X,A1,2I5,20X,''\'',I5)') CT,I1,J1,NPV 
        go to 10
       Endif
       K1=IAtomG(3,1,IGic)
       IF(Itype.eq.2) then
        Write(IOut,'(1X,A1,3I5,15X,''\'',I5)') CT,I1,J1,K1,NPV
        go to 10
       Endif
       L1=IAtomG(4,1,IGic)
       Write(IOut,'(1X,A1,4I5,10X,''\'',I5)') CT,I1,J1,K1,L1,NPV 
   10 continue
      return
      end
*Deck RBPCS
      Function RBPCS(IAt,JAt,IAn,RIJ)
      Implicit Real*8 (A-H,O-Z)
      Dimension IAN(*)
      Save Zero,One,Two,Three,CVK,CHyp
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/
      Data CVK/1.1D-3/, HypK/2.5d-02/
      BlOr=RIJ 
      IAnI=IAn(IAt)
      IAnJ=IAn(JAt)
      DltO=zero
      DltF=zero
      If(IAnI.eq.8.or.IAnJ.eq.8) DltO=one
      If(IAnI.eq.9.or.IAnJ.eq.9) DltF=one
      ZeffI=two
      ZeffJ=two
      If(IAnI.lt.3) ZeffI=one
      If(IAnJ.lt.3) ZeffJ=one
      If(IAnI.gt.10) ZeffI=three
      If(IAnJ.gt.10) ZeffJ=three
      Value=Blor
      Val0=RCovCT(IAnI,IAnJ)
      BndOrd=exp((Val0-Value)/3.0d-01)
      Trm1=-CVK*sqrt(ZeffI*ZeffJ-one)*val0
      Trm2=SQrt(Abs(BndOrd-Two+DltO*(one-BndOrd)))
      Trm3=-HypK*DltF*(BndOrd-one)**2
      DltR=Trm1*Trm2+Trm3
C Do not apply the correction to H-Bonds
      If(BndOrd.lt.3.0d-01) DltR=0.0d0 
      RBPCS=BlOr+DltR
      return
      end
*Deck Cyc3At 
      Subroutine Cyc3At(MxAtCy,IAt,JAt,KAt,NCyc,NAtC,ICAt,Fnd3C)
      Dimension NatC(*),ICAt(MxAtCy,*)
      Logical Fnd3C
      Fnd3c=.false.
      if(NCyc.eq.0) return
      do 10 ICyc=1,NCyc
       If(NAtC(ICyc).ne.3) go to 10
       is=0 
       do 20 ii=1,3
        LAt=ICAt(ii,ICyc)
        if(LAt.eq.IAt.or.LAt.eq.JAt.or.LAt.eq.KAt)is=is+1
   20  continue
       if(is.eq.3) then
        Fnd3C=.true.
        return
       endif
   10 continue
      return
      end
*Deck Cyc4At 
      Subroutine Cyc4At(MxAtCy,IAt,JAt,KAt,LAt,NCyc,NAtC,ICAt,Fnd4C)
      Dimension NatC(*),ICAt(MxAtCy,*)
      Logical Fnd4C
      Fnd4c=.false.
      if(NCyc.eq.0) return 
      do 10 ICyc=1,NCyc
       If(NAtC(ICyc).ne.4) go to 10
       is=0 
       do 20 ii=1,4
        MAt=ICAt(ii,ICyc)
        if(MAt.eq.IAt.or.MAt.eq.JAt.or.MAt.eq.KAt.or.MAt.eq.LAt)is=is+1
   20  continue
       if(is.eq.4) then
        Fnd4C=.true.
        return
       endif
   10 continue
      return
      end
*Deck FindHBnd
      Subroutine FindHBnd(IOut,IPrint,MxBnd,AllHB,NAtoms,NFrag,NHB,
     $  IAn,NBond,IBond,IFrag,C,IDn,IHAt,ICc)
      Implicit Real*8 (A-H,O-Z)
      Logical FndHB,AllHB,GemIJ,VicIJ
      Dimension C(3,*)
      Dimension IAn(*),IFrag(*),NBond(*),IBond(MxBnd,*)
      Dimension IDn(*),ICc(*),IHAt(*)
      INCLUDE 'bdpcs3_hbond_params.inc'
C Find H-bonds without modifying the covalent topology.
C GIC construction must use only the covalent graph; these contacts are
C non-covalent targets for geometry correction/constraints.
      NHB=0
      ToDeg=45.0d0/ATan(1.0d0)
      If(.Not.AllHB.and.NFrag.eq.1) return
      Do 10 IAt=1,NAtoms
       RIJMin=BDPCS3_HB_SEARCH_CUTOFF
       FndHB=.false.
       IAI=IAn(IAt)
       If(IAI.ne.1) go to 10
       If(NBond(IAt).lt.1) go to 10
       IDon=IBond(1,IAt)
       IADon=IAn(IDon)
C N-H, O-H and S-H donors.
       If(IADon.ne.7.and.IADon.ne.8.and.IADon.ne.16) go to 10
       Do 20 IAcc=1,NAtoms
        If(IAcc.eq.IAt.or.IAcc.eq.IDon) go to 20
        IAAcc=IAn(IAcc)
        If(IAAcc.ne.7.and.IAAcc.ne.8) go to 20
        If(IFrag(IDon).eq.IFrag(IAcc)) then
         If(.not.AllHB) go to 20
        EndIf
        RIJ=Distan(C,IAt,IAcc,0)
        if(RIJ.ge.BDPCS3_HB_SEARCH_CUTOFF) go to 20
        AngXHY=ValAng(C(1,IDon),C(1,IAt),C(1,IAcc))*ToDeg
        If(AngXHY.lt.BDPCS3_HB_ANGLE_MIN) go to 20
C do not consider geminal atoms
        GemIJ=.false.
        Do 30 IGJ=1,NBond(IAcc)
         VicIJ=.false.
         IGAt=IBond(IGJ,IAcc)
         If(IGAt.eq.IAt.or.IGAt.eq.IDon) go to 30
         If(IGAt.eq.IDon) then
           Write(IOut,'('' No H-Bond between Geminal Atoms:'',2I5)')
     $       IAcc,IDon
           GemIJ=.true.
          go to 20
         EndIf
C do not consider vicinal atoms
         Do 40 IVJ=1,NBond(IDon)
          IVAt=IBond(IVJ,IDon)
          If(IVAt.eq.IAt.or.IVAt.eq.IAcc) go to 40
          If(IVAt.eq.IGAt) then
           Write(IOut,'('' No H-Bond between Vicinal Atoms:'',2I5)')
     $       IAcc,IDon
           VicIJ=.true.
           go to 20
          EndIf
   40    Continue
   30   Continue
C only 1 H-bond for each Hydrogen atom
        If(GemIJ.or.VicIJ) go to 20
        If(RIJ.ge.RIJMin) go to 20
        FndHB=.true.
        RIJMin=RIJ
        JHBI=IAcc
   20  Continue
       If(.not.FndHB) go to 10
       NHB=NHB+1
       IDn(NHB)=IDon
       ICc(NHB)=JHBI
       IHAt(NHB)=IAt
   10 Continue
      Return
      End
*Deck MkHBnd
      Subroutine MkHBnd(IOut,IPrint,MxBnd,AllHB,NAtoms,NFrag,NHB,IAn,
     $  NBond,IBond,IFrag,C)
      Implicit Real*8 (A-H,O-Z)
      Logical AllHB
      Dimension C(3,*)
      Dimension IAn(*),IFrag(*),NBond(*),IBond(MxBnd,*)
C Local
      Dimension IDn(100),ICc(100),IHAt(100)
      Call FindHBnd(IOut,IPrint,MxBnd,AllHB,NAtoms,NFrag,NHB,IAn,
     $  NBond,IBond,IFrag,C,IDn,IHAt,ICc)
      if(NHB.gt.0) then
       NHBr=0
       Do 50 i1=1,NHB
        IAt=IDn(i1)
        JAt=IHAt(i1)
        KAt=ICc(i1)
C only 1 H-bond for each Donor-Acceptor pair
        Do 60 i2=1,i1
         if(i2.eq.i1) go to 60
         IAt2=IDn(i2)
         KAt2=ICc(i2)
         If(IAt.eq.IAt2.and.KAt.eq.KAt2) go to 50
   60   Continue
        NHBr=NHBr+1
        NBond(JAt)=NBond(JAt)+1
        NBond(KAt)=NBond(KAt)+1
        IBond(NBond(JAt),JAt)=KAt
        IBond(NBond(KAt),KAt)=JAt
        RIJ=Distan(C,IAt,KAt,0)
        Write(IOut,'('' Atoms:'',2I5,'' at Dist.='',F10.5,
     $    '' are bridged by H'',I5)') IAt,KAt,RIJ,JAt
   50  Continue
       NHB=NHBr
       if(AllHB) write(IOut,'(I5,'' H-Bonds added to topology'')') NHB
       if(.not.allHB) write(IOut,'(I5,'' Inter-Molecular H-Bonds'',
     $   '' Added to Topology'')') NHB
      EndIf
      Return
      End
*Deck OrdFrg
      Subroutine OrdFrg(IOut,IPrint,MxAtFr,NAtoms,NFrag,IFrag,NAtFr,
     $  IFrsAt,IlstAt,IAtFr,LConn)
      Implicit Real*8 (A-H,O-Z)
      Dimension IFrag(*),IFrsAt(*),ILstAt(*)
      Dimension NAtFr(*),IAtFr(MxAtFr,*)
      Logical LConn(*)
C IFrsAt(IFrag) = first atom of fragment IFrag
C ILstAt(IFrag) = last  atom of fragment Ifrag
      Call IClear(NFrag,NAtFr)
      Do 10 I=1,NFrag
       IFrsAt(I)=100000
       ILstAt(I)=0
   10 Continue 
      Do 20 IAt=1,NAtoms
       IFr=IFrag(IAt)
       NAtFr(Ifr)=NAtFr(IFr)+1
       IAtFr(NAtFr(IFr),IFr)=IAt
       if(IAt.lt.IFrsAt(IFr)) IFrsAt(IFr)=IAt       
       if(IAt.gt.ILstAt(IFr)) ILstAt(IFr)=IAt
   20 Continue
      Do 30 IFr=1,NFrag
       LConn(IFr)=.True.
       If(NAtFr(IFr).eq.1) go to 30
       Do 40 II=2,NAtFr(IFr)
        if(IAtFr(II,IFr)-IAtFr(II-1,IFr).ne.1) LConn(IFr)=.False.
   40  Continue
   30 Continue 
      Do 50 Ifr=1,NFrag
       If(LConn(IFr)) then
        write(IOut,'(''Mol'',I1,''=Fragment('',I3,''-'',I3,'')'')')IFr,
     $   IFrsAt(IFr),ILstAt(IFr)
       Else
        write(IOut,'(''Mol'',I1,''=Fragment('',I3)',advance='no')IFr,
     $    IAtFr(1,IFr)
        do 60 IAt=2,NAtFr(IFr)
         write(IOut,'('','',I3)',advance='no')IAtFr(IAt,IFr)
   60   continue
        write(IOut,'('')'')') 
       EndIf
   50 Continue
      write(IOut,'(''RMol1(Frozen)=Rotor(Mol1)'')')
      do 70 I=2,NFrag
       write(IOut,'(''Rotor(Mol'',I1,'')'')')I
   70 Continue
      Return
      End
*Deck MkBNew
      Subroutine MkBNew(IOut,IPrint,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,
     $  NLAng,NOupl,NDih,IAtmB,IAtmA,IAtmL,IAtmD,IAtmO,NTermB,NTermA,
     $  NTermL,NTermD,NtermO,CoefB,CoefA,CoefL,CoefD,CoefO,C,BMat,
     $  ImpDih)
      Implicit Real*8 (A-H,O-Z)
      Logical DoB1,ImpDih
      Dimension IAtmB(MxAtP,MxTrm,*),IAtmA(MxAtP,MxTrm,*)
      Dimension IAtmL(MxAtP,MxTrm,*),IAtmD(MxAtP,MxTrm,*)
      Dimension IAtmO(MxAtP,MxTrm,*)
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)
      Dimension C(3,*),BMat(*)
C Local
      Dimension B(3,4),DB(3,4,3,4)
      Dimension IB(4)
C
C Compute B and DB matrices IB is the equivalent of IAtmB, etc.
C Here distances are in Angstroms, whereas in Gaussian they are in Bohrs
C Note that there are no differences for stretchings, whereas there is a factor 
C 0.529177 for angles
C
      NTot=NLen+NAng+NLAng+NDih+NOupl 
      If(NTot.eq.0.or.NAtoms.eq.0) Return
      Call AClear(3*NAtoms*NTot,BMat)
      If(NLen.gt.0) then
       Do 100 ILen=1,NLen
        NTerm=NTermB(ILen)
        Do 110 ITerm=1,NTerm
         IAt=IAtmB(1,ITerm,ILen)
         JAt=IAtmB(2,ITerm,ILen)
         call AClear(12,B)
         If(DoB1) call AClear(144,DB) 
         call DBStr(1,IAt,JAt,B,IB,C,DB)    
         IXYZ=3*(IAt-1)
         JXYZ=3*(JAt-1)
         Do 120 IP=1,3
          Ind0=3*NAtoms*(ILen-1)+IP
          Ind1=Ind0+IXYZ
          Ind2=Ind0+JXYZ
          BMat(Ind1)=BMat(Ind1)+CoefB(ITerm,ILen)*B(IP,1)
          BMat(Ind2)=BMat(Ind2)+CoefB(ITerm,ILen)*B(IP,2)
  120    Continue    
  110   Continue 
  100  Continue
      EndIf
      IA0=NLen 
      If(NAng.gt.0) then
       Do 200 IAng=1,NAng
        NTerm=NTermA(IAng)
        Do 210 ITerm=1,NTerm
         IAt=IAtmA(1,ITerm,IAng)
         JAt=IAtmA(2,ITerm,IAng)
         KAt=IAtmA(3,ITerm,IAng)
         call AClear(12,B)
         If(DoB1) call Aclear(144,DB)
         call DBBend(1,IAt,JAt,KAt,B,IB,C,DB)
         IXYZ=3*(IAt-1)
         JXYZ=3*(JAt-1)
         KXYZ=3*(KAt-1)
         Do 220 IP=1,3
          Ind0=3*NAtoms*(IA0+IAng-1)+IP
          Ind1=Ind0+IXYZ
          Ind2=Ind0+JXYZ 
          Ind3=Ind0+KXYZ
          BMat(Ind1)=BMat(Ind1)+CoefA(ITerm,IAng)*B(IP,1)
          BMat(Ind2)=BMat(Ind2)+CoefA(ITerm,IAng)*B(IP,2)
          BMat(Ind3)=BMat(Ind3)+CoefA(ITerm,IAng)*B(IP,3)
  220    Continue
  210   Continue
  200  Continue
      EndIf
      IA0=IA0+NAng
      If(NLang.gt.0) then
       Do 300 ILAng=1,NLAng
        NTerm=NTermL(ILAng)    
        Do 310 ITerm=1,NTerm
         IAt=IAtmL(1,ITerm,ILAng)
         JAt=IAtmL(2,ITerm,ILAng)
         KAt=IAtmL(3,ITerm,ILAng)
         LL=IAtmL(4,ITerm,ILAng)
         Call AClear(12,B)
         If(DoB1) Call AClear(144,DB)
         If(LL.eq.-1.or.LL.eq.-2) then
          Call DBLinPy(IAt,JAt,KAt,LL,B,C)
         Else
          call DBLBnd(1,IAt,JAt,KAt,LL,B,IB,C,DB)
         EndIf
         IXYZ=3*(IAt-1)
         JXYZ=3*(JAt-1)
         KXYZ=3*(KAt-1)
         Do 320 IP=1,3
          Ind0=3*NAtoms*(IA0+ILAng-1)+IP
          Ind1=Ind0+IXYZ
          Ind2=Ind0+JXYZ
          Ind3=Ind0+KXYZ
          BMat(Ind1)=BMat(Ind1)+CoefL(ITerm,ILAng)*B(IP,1)
          BMat(Ind2)=BMat(Ind2)+CoefL(ITerm,ILang)*B(IP,2)
          BMat(Ind3)=BMat(Ind3)+CoefL(ITerm,ILAng)*B(IP,3)  
  320    Continue
  310   Continue
  300  Continue
      EndIf
      IA0=IA0+NLAng
      If(NDih.gt.0) then 
       Do 400 IDih=1,NDih
        NTerm=NTermD(IDih)
        Do 410 ITerm=1,NTerm
         IAt=IAtmD(1,ITerm,IDih)
         JAt=IAtmD(2,ITerm,IDih)
         KAt=IAtmD(3,ITerm,IDih)
         LAt=IAtmD(4,ITerm,IDih)
         Call Aclear(12,B)
         If(DoB1) Call AClear(144,DB)
         call DBTors(1,IAt,JAt,KAt,LAt,B,IB,C,DB)
         IXYZ=3*(IAt-1)
         JXYZ=3*(JAt-1)
         KXYZ=3*(KAt-1)
         LXYZ=3*(LAt-1) 
         Do 420 IP=1,3
          Ind0=3*NAtoms*(IA0+IDih-1)+IP
          Ind1=Ind0+IXYZ
          Ind2=Ind0+JXYZ
          Ind3=Ind0+KXYZ
          Ind4=Ind0+LXYZ 
          BMat(Ind1)=BMat(Ind1)+CoefD(ITerm,IDih)*B(IP,1)
          BMat(Ind2)=BMat(Ind2)+CoefD(ITerm,IDih)*B(IP,2)
          BMat(Ind3)=BMat(Ind3)+CoefD(ITerm,IDih)*B(IP,3)
          BMat(Ind4)=BMat(Ind4)+CoefD(ITerm,IDih)*B(IP,4)
  420    Continue
  410   Continue
  400  Continue
      EndIf
      IA0=IA0+NDih
      If(NOupl.gt.0) then
       Do 500 IOupl=1,NOupl
        NTerm=NTermO(IOupl)
        Do 510 ITerm=1,NTerm
         IAt=IAtmO(1,ITerm,IOupl)
         JAt=IAtmO(2,ITerm,IOupl)
         KAt=IAtmO(3,ITerm,IOupl)
         LAt=IAtmO(4,ITerm,IOupl)
         call Aclear(12,B)
         If(DoB1) call Aclear(144,DB)
         If(ImpDih) then
          call DBTors(1,JAt,IAt,LAt,KAt,B,IB,C,DB)
          IXYZ=3*(JAt-1)
          JXYZ=3*(IAt-1)
          KXYZ=3*(LAt-1)
          LXYZ=3*(KAt-1)
         Else
          call DBOOPl(IAt,JAt,KAt,LAt,B,IB,C,DB,COST,SINT)
          IXYZ=3*(IAt-1)
          JXYZ=3*(JAt-1)
          KXYZ=3*(KAt-1)
          LXYZ=3*(LAt-1)
         EndIf
         Do 520 IP=1,3
          Ind0=3*NAtoms*(IA0+IOupl-1)+IP
          Ind1=Ind0+IXYZ
          Ind2=Ind0+JXYZ
          Ind3=Ind0+KXYZ
          Ind4=Ind0+LXYZ
          BMat(Ind1)=BMat(Ind1)+CoefO(ITerm,IOupl)*B(IP,1)
          BMat(Ind2)=BMat(Ind2)+CoefO(ITerm,IOupl)*B(IP,2)
          BMat(Ind3)=BMat(Ind3)+CoefO(ITerm,IOupl)*B(IP,3)
          BMat(Ind4)=BMat(Ind4)+CoefO(ITerm,IOupl)*B(IP,4)
  520    Continue 
  510   Continue
  500  Continue
      EndIf 
      NTot=IA0+NOupl
C Print B matrix
      If(IPrint.eq.0) Return  
      Write(IOut,'('' B Matrix'')')
      Do 600 IInt=1,NTot
       Write(IOut,'('' Curvilinear Internal:'',I5)') IInt
       Ini=3*NAtoms*(IInt-1)+1
       IEnd=Ini+3*NAtoms-1
       Write(IOut,'(6F10.5)') (BMat(ii),ii=ini,iend)
  600 Continue
      Return
      End

*Deck DBLinPy
      Subroutine DBLinPy(I,J,K,Mode,B,C)
      Implicit Real*8(A-H,O-Z)
      Dimension B(3,4),C(3,*),CLoc(3,3)
      Integer At(3)
      Data H/1.0D-4/
      At(1)=I
      At(2)=J
      At(3)=K
      Do 30 IA=1,3
       Do 20 IC=1,3
        Do 10 JA=1,3
         Do 5 JC=1,3
          CLoc(JC,JA)=C(JC,At(JA))
    5    Continue
   10   Continue
        CLoc(IC,IA)=CLoc(IC,IA)+H
        Call LinPyVal(CLoc,Mode,VP)
        CLoc(IC,IA)=CLoc(IC,IA)-2.0D0*H
        Call LinPyVal(CLoc,Mode,VM)
        B(IC,IA)=(VP-VM)/(2.0D0*H)
   20  Continue
   30 Continue
      Return
      End

*Deck LinPyVal
      Subroutine LinPyVal(CLoc,Mode,Val)
      Implicit Real*8(A-H,O-Z)
      Dimension CLoc(3,3),U(3),V(3),Axis(3),E1(3),E2(3),BVec(3)
      RN1=0.0D0
      RN2=0.0D0
      Do 10 IC=1,3
       U(IC)=CLoc(IC,1)-CLoc(IC,2)
       V(IC)=CLoc(IC,3)-CLoc(IC,2)
       RN1=RN1+U(IC)*U(IC)
       RN2=RN2+V(IC)*V(IC)
   10 Continue
      RN1=DSqrt(RN1)
      RN2=DSqrt(RN2)
      If(RN1.lt.1.0D-12.or.RN2.lt.1.0D-12) then
       Val=0.0D0
       Return
      EndIf
      Do 20 IC=1,3
       U(IC)=U(IC)/RN1
       V(IC)=V(IC)/RN2
   20 Continue
      Axis(1)=1.0D0
      Axis(2)=0.0D0
      Axis(3)=0.0D0
      If(DAbs(U(1)).gt.0.9D0) then
       Axis(1)=0.0D0
       Axis(2)=1.0D0
      EndIf
      E1(1)=U(2)*Axis(3)-U(3)*Axis(2)
      E1(2)=U(3)*Axis(1)-U(1)*Axis(3)
      E1(3)=U(1)*Axis(2)-U(2)*Axis(1)
      RN=DSqrt(E1(1)*E1(1)+E1(2)*E1(2)+E1(3)*E1(3))
      If(RN.lt.1.0D-12) then
       Val=0.0D0
       Return
      EndIf
      Do 30 IC=1,3
       E1(IC)=E1(IC)/RN
   30 Continue
      E2(1)=U(2)*E1(3)-U(3)*E1(2)
      E2(2)=U(3)*E1(1)-U(1)*E1(3)
      E2(3)=U(1)*E1(2)-U(2)*E1(1)
      Do 40 IC=1,3
       BVec(IC)=V(IC)+U(IC)
   40 Continue
      Val=0.0D0
      If(Mode.eq.-1) then
       Do 50 IC=1,3
        Val=Val+BVec(IC)*E1(IC)
   50  Continue
      Else
       Do 60 IC=1,3
        Val=Val+BVec(IC)*E2(IC)
   60  Continue
      EndIf
      Return
      End
*Deck MakeB
      Subroutine MakeB(IOut,IPrint,MxAtP,MxTrm,NAtoms,NLenR,NAngR,
     $  NLAngR,NOuplR,NDihR,IAtmBR,IAtmAR,IAtmLR,IAtmDR,IAtmOR,C,BMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension IAtmBR(MxAtP,MxTrm,*),IAtmAR(MxAtP,MxTrm,*)
      Dimension IAtmLR(MxAtP,MxTrm,*),IAtmDR(MxAtP,MxTrm,*)
      Dimension IAtmOR(MxAtP,MxTrm,*)
      Dimension C(3,*),BMat(*)
C Local
      Dimension B(3,4),DB(3,4,3,4)
      Dimension IB(4)
C
C Compute B and DB matrices IB is the equivalent of IAtmBr, etc.
C Here distances are in Angstroms, whereas in Gaussian they are in Bohrs
C Note that there are no differences for stretchings, whereas there is a factor 
C 0.529177 for angles
C
      Call AClear(9*NAtoms*NAtoms,BMat)
      If(NLenR.gt.0) then
       Do 100 ILen=1,NLenR
        IAt=IAtmBR(1,1,ILen)
        JAt=IAtMBR(2,1,ILen)
        call AClear(12,B)
        call AClear(144,DB) 
        call DBStr(1,IAt,JAt,B,IB,C,DB)    
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        Do 102 IP=1,3
         Ind0=3*NAtoms*(ILen-1)+IP
         Ind1=Ind0+IXYZ
         Ind2=Ind0+JXYZ
         BMat(Ind1)=B(IP,1)
         BMat(Ind2)=B(IP,2)
  102   Continue    
  100  Continue
      EndIf
      IA0=NLenR 
      If(NAngR.gt.0) then
       Do 110 IAng=1,NAngR
        IAt=IAtmAR(1,1,IAng)
        JAt=IAtmAR(2,1,IAng)
        KAt=IAtmAr(3,1,IAng)
        call AClear(12,B)
        call Aclear(144,DB)
        call DBBend(1,IAt,JAt,KAt,B,IB,C,DB)
CENZO
C       Write(IOut,'(''Bend'',I3)') IAng
C       write(IOut,'(I3,3F12.5)') IAt,(B(IPPO,1,1),IPPO=1,3)
C       write(IOut,'(I3,3F12.5)') JAt,(B(IPPO,2,1),IPPO=1,3)
C       write(IOut,'(I3,3F12.5)') KAt,(B(IPPO,3,1),IPPO=1,3)
CENZO
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        KXYZ=3*(KAt-1)
        Do 112 IP=1,3
         Ind0=3*NAtoms*(IA0+IAng-1)+IP
         Ind1=Ind0+IXYZ
         Ind2=Ind0+JXYZ 
         Ind3=Ind0+KXYZ
         BMat(Ind1)=B(IP,1)
         BMat(Ind2)=B(IP,2)
         BMat(Ind3)=B(IP,3)
  112   Continue
  110  Continue
      EndIf
      IA0=IA0+NAngR
      If(NLangR.gt.0) then
       Do 120 ILAng=1,NLAngR
        IAt=IAtmLR(1,1,ILAng)
        JAt=IAtmLR(2,1,ILAng)
        KAt=IAtmLR(3,1,ILAng)
        LL=0
        Call AClear(12,B)
        Call AClear(144,DB)
        call DBLBnd(1,IAt,JAt,KAt,LL,B,IB,C,DB)
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        KXYZ=3*(KAt-1)
        Do 122 IP=1,3
         Ind0=3*NAtoms*(IA0+ILAng-1)+IP
         Ind1=Ind0+IXYZ
         Ind2=Ind0+JXYZ
         Ind3=Ind0+KXYZ
         BMat(Ind1)=B(IP,1)
         BMat(Ind2)=B(IP,2)
         BMat(Ind3)=B(IP,3)  
  122   Continue
  120  Continue
      EndIf
      IA0=IA0+NLAngR
      If(NDihR.gt.0) then 
       Do 130 IDih=1,NDihR
        IAt=IAtmDR(1,1,IDih)
        JAt=IAtmDR(2,1,IDih)
        KAt=IAtmDR(3,1,IDih)
        LAt=IAtmDR(4,1,IDih)
        Call Aclear(12,B)
        Call AClear(144,DB)
        call DBTors(1,IAt,JAt,KAt,LAt,B,IB,C,DB)
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        KXYZ=3*(KAt-1)
        LXYZ=3*(LAt-1) 
        Do 132 IP=1,3
         Ind0=3*NAtoms*(IA0+IDih-1)+IP
         Ind1=Ind0+IXYZ
         Ind2=Ind0+JXYZ
         Ind3=Ind0+KXYZ
         Ind4=Ind0+LXYZ 
         BMat(Ind1)=B(IP,1)
         BMat(Ind2)=B(IP,2)
         BMat(Ind3)=B(IP,3)
         BMat(Ind4)=B(IP,4)
  132   Continue
  130  Continue
      EndIf
      IA0=IA0+NDihR
      If(NOuplr.gt.0) then
       Do 140 IOupl=1,NOuplR
        IAt=IAtmOR(1,1,IOupl)
        JAt=IAtmOR(2,1,IOupl)
        KAt=IAtmOR(3,1,IOupl)
        LAt=IAtmOR(4,1,IOupl)
        call Aclear(12,B)
        call Aclear(144,DB)
        call DBOOPl(IAt,JAt,KAt,LAt,B,IB,C,DB,COST,SINT)
        IXYZ=3*(IAt-1)
        JXYZ=3*(JAt-1)
        KXYZ=3*(KAt-1)
        LXYZ=3*(LAt-1)
        Do 142 IP=1,3
         Ind0=3*NAtoms*(IA0+IOupl-1)+IP
         Ind1=Ind0+IXYZ
         Ind2=Ind0+JXYZ
         Ind3=Ind0+KXYZ
         Ind4=Ind0+LXYZ
         BMat(Ind1)=B(IP,1)
         BMat(Ind2)=B(IP,2)
         BMat(Ind3)=B(IP,3)
         BMat(Ind4)=B(IP,4)
  142   Continue
  140  Continue
      EndIf 
      NTotR=IA0+NOuplR
C Print B matrix
      If(IPrint.eq.0) Return  
      Write(IOut,'('' B Matrix'')')
      Do 150 IInt=1,NTotR
       Write(IOut,'('' Redundant Internal:'',I5)') IInt
       Ini=3*NAtoms*(IInt-1)+1
       IEnd=Ini+3*NAtoms-1
       Write(IOut,'(6F10.5)') (BMat(ii),ii=ini,iend)
  150 Continue
       Return
       End
*Deck MkNewR
      Subroutine MkNewR(IOut,IPrint,NAtoms,MxAtP,MxTrm,NTot,NLen,NAng,
     $  NLAng,NDih,NOuPl,IAtmB,IAtmA,IAtmL,IAtmD,IAtmO,DoneI,IntMax,
     $  DltIMx,DltIAv,C,R,RRef)
      Implicit Real*8 (A-H,O-Z)
      Logical DoneI
      Dimension C(3,*),R(*),RRef(*)
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
      DoneI=.False.
      Do 10 ILen=1,NLen
       IRed=ILen
       IAt=IAtmB(1,1,ILen)
       JAt=IAtmB(2,1,ILen)
       R(IRed)=Distan(C,IAt,JAt,0)  
       Delta=Abs(R(IRed)-RRef(IRed))
       If(Delta.gt.DltIMx) then
        IntMax=IRed 
        DltIMx=Delta
       EndIf
       DlTot=DlTot+Delta
       If(IPrint.gt.0) Write(IOut,'(I3,'' Bond Length   '',2I5,10X,
     $  '' Expected'', F12.5,'' Actual'',F12.5)') IRed,IAt,JAt,
     $  RRef(IRed),R(IRed)
   10 Continue
      If(NAng.gt.0) then
       Do 20 IAng=1,NAng
        IRed=NLen+IAng
        IAt=IAtmA(1,1,IAng)
        JAt=IAtmA(2,1,IAng)
        KAt=IAtmA(3,1,IAng)
        R(IRed)=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        Delta=Abs(R(IRed)-RRef(IRed))
        If(Delta.gt.DltIMx) then
         IntMax=IRed
         DltIMx=Delta
        EndIf
        DlTot=DlTot+Delta
        If(IPrint.gt.0) Write(IOut,'(I3,'' Valence Angle '',3I5,5X,
     $  '' Expected'', F12.5,'' Actual'',F12.5)') IRed,IAt,JAt,KAt,
     $  RRef(IRed)*ToDeg,R(IRed)*ToDeg
   20  Continue
      EndIf
      If(NLAng.gt.0) then
       Do 30 ILAng=1,NLAng
        IRed=NLen+NAng+ILAng
        IAt=IAtmL(1,1,ILAng)
        JAt=IAtmL(2,1,ILAng)
        KAt=IAtmL(3,1,ILAng)
        R(IRed)=ValAng(C(1,IAt),C(1,JAt),C(1,KAt))
        Delta=Abs(R(IRed)-RRef(IRed))
        If(Delta.gt.DltIMx) then
         IntMax=IRed
         DltIMx=Delta
        EndIf
        DlTot=DlTot+Delta
        If(IPrint.gt.0) Write(IOut,'(I3,'' Linear Angle  '',3I5,5X,
     $  '' Expected'', F12.5,'' Actual'',F12.5)') IRed,IAt,JAt,KAt,
     $  RRef(IRed)*ToDeg,R(IRed)*ToDeg
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
        Delta=Abs(R(IRed)-RRef(IRed))
        If(Delta.gt.DltIMx) then
         IntMax=IRed
         DltIMx=Delta
        EndIf
        DlTot=DlTot+Delta
        If(IPrint.gt.0) Write(IOut,'(I3,'' Dihedral Angle'',4I5,
     $  '' Expected'', F12.5,'' Actual'',F12.5)') IRed,IAt,JAt,KAt,
     $  LAt,RRef(IRed)*ToDeg,R(IRed)*ToDeg
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
        Delta=Abs(R(IRed)-RRef(IRed))
        If(Delta.gt.DltIMx) then
         IntMax=IRed
         DltIMx=Delta
        EndIf
        DlTot=DlTot+Delta
        If(IPrint.gt.0) Write(IOut,'(I3,'' Out Pl. Angle '',4I5,
     $  '' Expected'', F12.5,'' Actual'',F12.5)') IRed,IAt,JAt,KAt,
     $  LAt,RRef(IRed)*ToDeg,R(IRed)*ToDeg
   50  Continue
      EndIf
      DltIAv=DlTot/Float(NTot)
      If(DltIAv.lt.TrAv.and.DltIMx.lt.TrMax) DOneI=.True.
      Return
      End 
*Deck MkNewC
      Subroutine MkNewC(IOut,IPrint,NAtoms,NGIC,DoneC,DiffMx,DiffAv,IAn,
     $  COld,CNew,QOld,QNew,G1BMat)
      Implicit Real*8 (A-H,O-Z)
      Logical DoneC
      Dimension IAn(*)
      Dimension COld(3,*),CNew(3,*),QOld(*),QNew(*),G1BMat(3*NAtoms,*)
C Angles in radiants
      DoneC=.False.
      DiffT=0.0d0
      DiffMx=0.0d0
      TrAv=1.0d-6
      TrMx=1.0d-5
      Do 10 IAt=1,NAtoms
       Do 20 ICoord=1,3
        CNew(ICoord,IAt)=COld(ICoord,IAt)
        IXYZ=3*(IAt-1)+ICoord
        Do 30 Iint=1,NGIC
         CNew(ICoord,IAt)=CNew(ICoord,IAt)+G1BMat(IXYZ,Iint)*
     $     (QNew(Iint)-QOld(Iint))
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
*Deck MakeG
      Subroutine MakeG(IOut,IPrint,DoMW,NAtoms,NInt,AtMass,BMAt,GMAt)
      Implicit Real*8 (A-H,O-Z)
      Logical DoMW
      Dimension AtMass(*),BMat(3*NAtoms,*),GMat(NInt,*)
      Call AClear(NInt*NInt,GMat)
      IJ=0
      Do 10 I=1,NInt
       Do 20 J=1,NInt
        GMAt(I,J)=0.0d0
        Do 30 K=1,3*NAtoms
         If(DoMW) then
          IAt=Float((K-1)/3+1)
          FactK=1.0d0/AtMass(IAt)
         Else
          FactK=1.0d0
         EndIf
         GMat(I,J)=GMat(I,J)+BMAt(K,I)*BMAt(K,J)*FactK
   30   Continue
   20  Continue
   10 Continue
C Print G matrix
      If(IPrint.eq.0) Return
      If(DoMW) then
       Write(IOut,'(/,'' G Matrix'')')
      Else
       Write(IOut,'(/,'' BB+ Matrix'')')
      EndIf
      Do 150 IInt=1,NInt
       Write(IOut,'('' Redundant Internal:'',I5)') IInt
       Write(IOut,'(6F10.5)') (GMat(iint,ii),ii=1,NInt)
  150 Continue
      Return
      End
*Deck MakeGW
      Subroutine MakeGW(IOut,IPrint,DoMW,NAtoms,NInt,AtMass,WInt,
     $ BMat,GMat)
      Implicit Real*8 (A-H,O-Z)
      Logical DoMW
      Dimension AtMass(*),WInt(*),BMat(3*NAtoms,*),GMat(NInt,*)
C Weighted internal-coordinate metric:
C   G(I,J)=sqrt(WI)*sqrt(WJ)*B(I) M^-1 B(J)^T.
C WInt are physical least-squares weights, not row scale factors.
      Call AClear(NInt*NInt,GMat)
      Do 10 I=1,NInt
       WI=DSqrt(DMax1(0.0d0,WInt(I)))
       Do 20 J=1,NInt
        WJ=DSqrt(DMax1(0.0d0,WInt(J)))
        GMAt(I,J)=0.0d0
        Do 30 K=1,3*NAtoms
         If(DoMW) then
          IAt=Float((K-1)/3+1)
          FactK=1.0d0/AtMass(IAt)
         Else
          FactK=1.0d0
         EndIf
         GMat(I,J)=GMat(I,J)+WI*WJ*BMAt(K,I)*BMAt(K,J)*FactK
   30   Continue
   20  Continue
   10 Continue
      If(IPrint.eq.0) Return
      Write(IOut,'(/,'' Weighted internal-coordinate metric'')')
      Do 150 IInt=1,NInt
       Write(IOut,'('' Internal:'',I5,'' Weight:'',D12.5)')
     $  IInt,WInt(IInt)
       Write(IOut,'(6F10.5)') (GMat(iint,ii),ii=1,NInt)
  150 Continue
      Return
      End
*Deck MakGm1
      Subroutine MakGm1(IOut,IPrint,DoMW,NTotR,NRed,GMat,D,IScr)
      Implicit Real*8 (A-H,O-Z)
      Logical Inv1,DoMW
      Dimension GMat(*),D(*),IScr(*)
      If(NRed.gt.0) then
       Call GenInv(IOut,IPrint,NTotR,NRed,GMat,IScr,D)
      Else
       InIS=1
       InIAD1=INIS+2*NTotR
       InIAd2=InIAd1+NTotR
       If(.not.Inv1(GMat,NTotR,IScr(InIS),IScr(InIAD1),IScr(InIAD2),D,
     $   NTotR,Det)) then
        Write(IOut,'('' Inversion of G Matrix Failed'')')
        Stop
       EndIf
      EndIf
      If(IPrint.eq.0) Return
      If(DoMW) then
       Write(IOut,'(/,'' G-1 Matrix'')')
      Else
       Write(IOut,'(/,'' (BB+)-1 Matrix'')')
      EndIf
      Do 150 IR=1,NTotR
       Write(IOut,'('' Redundant Internal:'',I5)') IR
       Ini=NTotR*(IR-1)+1
       IEnd=Ini+NTotR-1
       Write(IOut,'(6F10.5)') (GMat(ii),ii=ini,iend) 
  150 Continue
      Return
      End
*Deck GenInv
      Subroutine GenInv(IOut,IPrint,NInt,NRed,GMat,IScr,V)
      Implicit Real*8 (A-H,O-Z)
      Dimension IScr(*) 
      Dimension GMat(Nint,*),V(*)
      Tresh=1.0D-05
C Split scratch for diagonalization
C V(I1) G  matrix (lower triangular)
      I1=1 
C V(I2) Eigenvalues
      I2=I1+NInt*(NInt+1)/2
C V(I3)=EVec
      I3=I2+NInt
C V(I4)=WA or Scratch
      I4=I3+NInt*NInt
C V(I5)=Scr
      I5=I4+6*NInt
      I6=I5+NInt*(NInt+1)/2+1
C
      IJ=0
      Do 10 I=1,NInt
       Do 20 J=1,I
        IJ=IJ+1
        V(IJ)=GMat(I,J)
   20  Continue
   10 Continue
      Call HQRII1(IOut,NInt,1,NInt,0,V(I1),V(I2),NInt,V(I3),.true.,
     $  IErr,IScr,V(I4),V(I5),I6-I5)
      If(IErr.ne.0)write(IOut,'(/,'' After HQRII1: IErr ='',I3,
     $  '' NInt ='',I3)') IErr,NInt
      Do 30 I=1,NInt
       EVal=Abs(V(I2+I-1))
       IF(EVal.lt.tresh) then
        V(I2+I-1)=0.0d0
        IniVec=(I-1)*NInt+I3
        Call AClear(NInt,V(IniVec))
       EndIf
       If(IPrint.gt.0) then
        write(IOut,'(''EigenValue'',D12.5)') V(I2+I-1)
        Ini=(I-1)*NInt+I3-1
        write(IOut,'(7D12.5)') (V(Ini+J),J=1,NInt)
       EndIf
   30 Continue
      Call MkGGm1(IOut,IPrint,NInt,NRed,V(I2),V(I3),V(I4))
      Call AMove(NInt*NInt,V(I4),GMat)
      Return
      End
*Deck MkGGm1 
      Subroutine MkGGm1(IOut,IPrint,N,NRed,EVal,EVec,GM1) 
      Implicit Real*8 (A-H,O-Z)
      Dimension EVal(*),EVec(N,*),GM1(N,*)
      Tresh=1.0d-5
      NRed=0
      Do 10 I=1,N
       If(ABs(EVal(I)).gt.tresh) then
        EVal(I)=1.0d0/EVal(I)
       Else
        EVal(I)=0.0d0
        NRed=NRed+1
       EndIf
   10 Continue
      Do 20 I=1,N
       Do 30 J=1,N
        GM1(I,J)=0.0d0
        Do 40 K=1,N
         GM1(I,J)=GM1(I,J)+EVec(I,K)*EVal(K)*EVec(J,K)
   40   Continue
   30  Continue
   20 Continue
      If(IPrint.eq.0) Return   
      Write(IOut,'(/,'' From MkGGm1:'',I3,'' null eigenvalues'')') NRed
      Write(IOut,'(/,''GM1 Matrix'')')
      Do 50 I=1,N
       Write(IOut,'('' Redundant Internal:'',I5)') I
       Write(IOut,'(7D11.5)') (GM1(I,JJ),JJ=1,N)
   50 Continue
      Return
      End
*Deck MkGm1B
      Subroutine MkGm1B(IOut,IPrint,NAtoms,NInt,BMAt,G1MAt,G1BMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension BMat(3*NAtoms,*),G1Mat(NInt,*),G1BMat(3*NAtoms,*)
      Call AClear(NInt*3*NAtoms,G1BMat)
      Do 10 I=1,3*NAtoms
       Do 20 J=1,NInt
        G1BMAt(I,J)=0.0d0
        Do 30 K=1,NInt
         G1BMat(I,J)=G1BMat(I,J)+BMat(I,K)*G1MAt(K,J)
   30   Continue
   20  Continue
   10 Continue
C Print Gm1B matrix
      If(IPrint.eq.0) Return
      Write(IOut,'(/,'' (BB+)-1B Matrix'')')
      Do 40 IXYZ=1,3*NAtoms
       Write(IOut,'('' Cartesian:'',I5)') IXYZ
       Write(IOut,'(6F10.5)') (G1BMat(IXYZ,Int),Int=1,NInt)
   40 Continue
      Return
      End
*Deck MkGm1BW
      Subroutine MkGm1BW(IOut,IPrint,NAtoms,NInt,WInt,BMAt,G1MAt,
     $ G1BMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension WInt(*),BMat(3*NAtoms,*),G1Mat(NInt,*)
      Dimension G1BMat(3*NAtoms,*)
C Effective Cartesian update matrix for raw internal displacements dQ:
C   dX = B'^T (B'B'^T)^-1 S dQ, with B'=S B.
      Call AClear(NInt*3*NAtoms,G1BMat)
      Do 10 I=1,3*NAtoms
       Do 20 J=1,NInt
        WJ=DSqrt(DMax1(0.0d0,WInt(J)))
        G1BMAt(I,J)=0.0d0
        Do 30 K=1,NInt
         WK=DSqrt(DMax1(0.0d0,WInt(K)))
         G1BMat(I,J)=G1BMat(I,J)+BMat(I,K)*WK*G1MAt(K,J)*WJ
   30   Continue
   20  Continue
   10 Continue
      If(IPrint.eq.0) Return
      Write(IOut,'(/,'' Weighted (BB+)-1B Matrix'')')
      Do 40 IXYZ=1,3*NAtoms
       Write(IOut,'('' Cartesian:'',I5)') IXYZ
       Write(IOut,'(6F10.5)') (G1BMat(IXYZ,Int),Int=1,NInt)
   40 Continue
      Return
      End
*Deck DrPMom
      Subroutine DrPMom(NAtoms,PMom1,C,AtMass)
      Implicit Real*8 (A-H,O-Z)
C Computes the first-order Cartesian derivatives of the diagonal
C entries of the moment of inertia tensor.
      Dimension PMom1(3,*),AtMass(*),C(3,*)
      Data zero/0.0d0/,two/2.0d0/
C PMom1(3,3*NAtoms) contains first derivatives of PMom
      do 10 i = 1, Natoms
C XX component:
       PMom1(1,3*(i-1)+1) = zero
       PMom1(1,3*(i-1)+2) = two*AtMass(i)*C(i,2)
       PMom1(1,3*(i-1)+3) = two*AtMass(i)*C(i,3)
C YY component:
       PMom1(2,3*(i-1)+1) = two*AtMass(i)*C(i,1)
       PMom1(2,3*(i-1)+2) = zero
       PMom1(2,3*(i-1)+3) = two*AtMass(i)*C(i,3)
C ZZ component:
       PMom1(3,3*(i-1)+1) = two*AtMass(i)*C(i,1)
       PMom1(3,3*(i-1)+2) = two*AtMass(i)*C(i,2)
       PMom1(3,3*(i-1)+3) = zero
   10 Continue
      Return
      End
*Deck TRBM
      Subroutine TRBM(Natoms,FlagM,AtMass,Crd,Bas)
      Implicit Real*8(A-H,O-Z)
C
C Translation-Rotational B Matrix
C     Compute the six basis vectors for rotations and translation in
C     mass weighted (if FlagM=.True.) Cartesian coordinates.
C     See M.Page and J.W. McIver : J.Chem.Phys. 88, 922-935 (1988)
C
C Input:
C     FlagM  : Apply mass weighting if True
C O   AtMass : (NAtoms) Atomic masses
C     Crd    : (3,NAtoms) Atomic coordinates
C
C Output:
C     Bas    : (NAt3,6) Rotational and translational basis array
C
C     Dimensions
      Integer NAtoms
C     Input
      Real*8 Crd(3,NAtoms), AtMass(*)
      Logical FlagM
C     Output
      Real*8 Bas(3*Natoms,6)
C     Local
      Integer i, ia
      Real*8 One, Wt, X, Y, Z
      Save One
      Data One/1.0D0/
C
      Call AClear(3*Natoms*6,Bas)
      i = 1
      Do 100 ia = 1, NAtoms
        Wt = Sqrt(AtMass(ia))
        If(.not.FlagM) Wt = One
        X          =  Crd(1,ia)
        Y          =  Crd(2,ia)
        Z          =  Crd(3,ia)
        Bas(i  ,1) =  Wt
        Bas(i+1,2) =  Wt
        Bas(i+2,3) =  Wt
        Bas(i+1,4) =  Wt*Z
        Bas(i+2,4) = -Wt*Y
        Bas(i  ,5) = -Wt*Z
        Bas(i+2,5) =  Wt*X
        Bas(i  ,6) =  Wt*Y
        Bas(i+1,6) = -Wt*X
        i          = i+3
  100   Continue
      Return
      End
*Deck Bend
      Subroutine Bend(NoInt,I,J,K,B,IB,C)
      Implicit Real*8(A-H,O-Z)
C
C     Adapted from the normal coordinate analysis program of Schachtschneider.
C
      Common/IO/ in,iout,ipunch
      Dimension B(3,4,*),IB(4,*),C(3,*)
C
      Call BendX(NoInt,I,J,K,B,IB,C,IFail)
      If(IFail.eq.1) then
       write(IOut,'('' Near-zero distance in Bend'')')
       Stop 
      else if(IFail.eq.2) then
       write(IOut,'('' Linear angle in Bend'')')
       Stop 
      else if(IFail.ne.0) then
       write(IOut,'('' Unknown failure in Bend'')')
       Stop 
      endIf
      Return
      End
*Deck BendX
      Subroutine BendX(NoInt,I,J,K,B,IB,C,IFail)
      Implicit Real*8(A-H,O-Z)
C
C     Adapted from the normal coordinate analysis program of Schachtschneider.
C
      Real*8 MDCutO
      Dimension B(3,4,*),IB(4,*),C(3,*),RJI(3),RJK(3),EJI(3),EJK(3)
      Save Zero, One
      Data Zero,One/0.D0,1.D0/
C
      Small = MDCutO(0)
      IB(1,NoInt) = I
      IB(2,NoInt) = J
      IB(3,NoInt) = K
      IB(4,NoInt) = 0
      DJISQ = Zero
      DJKSQ = Zero
      Do 10 M = 1, 3
        RJI(M) = C(M,I) - C(M,J)
        RJK(M) = C(M,K) - C(M,J)
        DJISQ = DJISQ + RJI(M)**2
        DJKSQ = DJKSQ + RJK(M)**2
   10 Continue
      DJI = Sqrt(DJISQ)
      DJK = Sqrt(DJKSQ)
      If(Abs(DJI).lt.Small.or.Abs(DJK).lt.Small) then
        IFail = 1
        Return
        endIf
      DotJ = Zero
      Do 20 M = 1, 3
        EJI(M) = RJI(M) / DJI
        EJK(M) = RJK(M) / DJK
        DotJ = DotJ + EJI(M)*EJK(M)
   20 Continue
      SinJ  = Sqrt(One-Min(DotJ**2,One))
      If(SinJ.lt.Small) then
        IFail = 2
        Return
        endIf
      Do 30 M = 1, 3
        B(M,3,NoInt) = ((DotJ*EJK(M)-EJI(M)))/(DJK*SinJ)
        B(M,1,NoInt) = ((DotJ*EJI(M)-EJK(M)))/(DJI*SinJ)
        B(M,2,NoInt) = -B(M,1,NoInt)-B(M,3,NoInt)
   30 Continue
      IFail = 0
      Return
      End
*Deck DBBend
      Subroutine DBBend(NoInt,I,J,K,B,IB,C,DB)
      Implicit Real*8(A-H,O-Z)
C
C     Compute bend matrix elements for Wilson B-matrix and their
C     cartesian first derivative.
C
C     Philippe Y. Ayala /Oct. 93
C
      Dimension B(3,4,*),DB(3,4,3,4,*),IB(4,*),C(*),RJI(3),RJK(3),
     $  EJI(3),EJK(3),BB(3,4,2),IBB(4,2),DBB(3,4,3,4,2)
      Save Zero, One
      Data Zero/0.d0/,One/1.d0/
C
      IAIND=3*(I-1)
      JAIND=3*(J-1)
      KAIND=3*(K-1)
      IB(1,NoInt)=I
      IB(2,NoInt)=J
      IB(3,NoInt)=K
      IB(4,NoInt)=0
      DJISQ=Zero
      DJKSQ=Zero
      Do 120 M=1,3
        RJI(M)=C(M+IAIND)-C(M+JAIND)
        RJK(M)=C(M+KAIND)-C(M+JAIND)
        DJISQ = DJISQ + RJI(M)**2
        DJKSQ = DJKSQ + RJK(M)**2
  120 Continue
      DJI =Sqrt(DJISQ)
      DJK =Sqrt(DJKSQ)
      DOTJ = Zero
      Do 132 M=1,3
        EJI(M)=RJI(M)/DJI
        EJK(M)=RJK(M)/DJK
        DOTJ=DOTJ+EJI(M)*EJK(M)
  132 Continue
      SINJ=Sqrt(ONE-DOTJ**2)
      Do 144 M=1,3
        B(M,3,NoInt)=((DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
        B(M,1,NoInt)=((DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
        B(M,2,NoInt)=-B(M,1,NoInt)-B(M,3,NoInt)
 144  Continue
C
C     Compute matrix elements needed for partial first derivatives
C
      Call AClear(24,BB)
      Call AClear(288,DBB)
      Call DBStr(1,I,J,BB,IBB,C,DBB)
      Call DBStr(2,J,K,BB,IBB,C,DBB)
      NoIntIJ=1
      NoIntJK=2
      Do 200 M=1,3
        Do 190 N=M,3
          DB(M,3,N,3,NoInt)=( -SINJ*B(N,3,NoInt)*EJK(M)
     $      + DOTJ*DBB(M,2,N,2,NoIntJK)
     $      - B(M,3,NoInt)*(EJK(N)*SINJ
     $      + DOTJ*B(N,3,NoInt)*DJK))/(DJK*SINJ)
          DB(M,3,N,1,NoInt)=(-SINJ*B(N,1,NoInt)*EJK(M)
     $      - DBB(M,1,N,1,NoIntIJ)
     $      - B(M,3,NoInt)*DJK*DOTJ*B(N,1,NoInt))/(DJK*SINJ)
          DB(M,1,N,1,NoInt)=(-SINJ*B(N,1,NoInt)*EJI(M)
     $      + DOTJ*DBB(M,1,N,1,NoIntIJ)
     $      - B(M,1,NoInt)*(DOTJ*B(N,1,NoInt)*DJI
     $      + EJI(N)*SINJ))/(DJI*SINJ)
          DB(M,1,N,3,NoInt)=(-SINJ*B(N,3,NoInt)*EJI(M)
     $      - DBB(M,2,N,2,NoIntJK)
     $      - B(M,1,NoInt)*DJI*DOTJ*B(N,3,NoInt))/(DJI*SINJ)
C
C         Make use of symmetry and linear dependency relationship
C
          DB(N,3,M,3,NoInt)=DB(M,3,N,3,NoInt)
          DB(N,1,M,1,NoInt)=DB(M,1,N,1,NoInt)
          DB(N,1,M,3,NoInt)=DB(M,3,N,1,NoInt)
          DB(N,3,M,1,NoInt)=DB(M,1,N,3,NoInt)
          DB(M,3,N,2,NoInt)=-DB(M,3,N,1,NoInt)-DB(M,3,N,3,NoInt)
          DB(N,3,M,2,NoInt)=-DB(N,3,M,1,NoInt)-DB(N,3,M,3,NoInt)
          DB(M,1,N,2,NoInt)=-DB(M,1,N,1,NoInt)-DB(M,1,N,3,NoInt)
          DB(N,1,M,2,NoInt)=-DB(N,1,M,1,NoInt)-DB(N,1,M,3,NoInt)
          DB(N,2,M,3,NoInt)=DB(M,3,N,2,NoInt)
          DB(M,2,N,3,NoInt)=DB(N,3,M,2,NoInt)
          DB(N,2,M,1,NoInt)=DB(M,1,N,2,NoInt)
          DB(M,2,N,1,NoInt)=DB(N,1,M,2,NoInt)
  190    Continue
  200   Continue
      Do 290 M=1,3
       Do 300 N=1,3
        DB(M,2,N,2,NoInt)=-DB(M,1,N,2,NoInt)-DB(M,3,N,2,NoInt)
  300  Continue
  290 Continue
      Return
      End
*Deck DBLBnd
      Subroutine DBLBnd(ICor,I,J,K,LL,B,IB,C,DB)
      Implicit Real*8(A-H,O-Z)
C
C     Calculate B-matrix for linear angle bend type coordinates
C     ICor  : index of the coordinate in the arrays
C     I,J,K : atom numbers
C     LL    : type flag and may contain the fourth atom number
C             LL > -10        : 3 atomic linear bend
C               -LL = (ITyp-1)*3+2+IAx, where ITyp 1 or 2 for the two
C               linear bend and IAx shows the axel for the projection.
C               If LL is -1 or -2 set up an axel and give it back in LL.
C             LL < -10
C             -1 or -2 - 10*L : use atom L in the linear 4 atom
C                               bend definition
C     B     : B-matrix array
C     IB    : atom indexes
C     C     : Cartesian coordinates
C
C     O. Farkas, 1997, WSU.
C
      Common /IO/ In, IOut, IPunch
      Dimension B(3,4,*),IB(4,*),C(3,*),R(3),BB(3,4),IBB(4),
     $  RJI(3),RJK(3),DB(3,4,3,4,*),DBB(3,4,3,4),IP(4), VP(3)
      Save Zero, One, Two, Eps, Half
      Data Zero/0.D0/, One/1.D0/, Two/2.D0/, Eps/1.D-13/, Half/0.5D0/
 1000 Format(' DBLBnd:  LL=',I10)
C
      L = 0
      If(LL.lt.0.and.LL.ge.-2) then
C
C       The axes are not decided yet, which is an error.
C
        Write(IOut,'('' Error in DBLBnd: LL='',I10)') LL
        Stop
      else if(LL.lt.-10) then
C
C     This is a new type linear bend with 4 atoms, get the 4th number
C     and the type
C
        LLL = IAbs(LL)
        L = LLL/10
        LLL = Mod(LLL,10)
C
C     Get the axel for the 3 atomic projected linear bend
C
      else
        L = 0
        LLL = Mod(-LL,3)+1
        endIf
C
      IB(1,ICor) = I
      IB(2,ICor) = J
      IB(3,ICor) = K
      IB(4,ICor) = L
      If(L.eq.0) then
C
C       This is the linear bend with 3 atoms
C
        Call AClear(3,R)
        R(LLL) = One
        Call ASub(3,C(1,I),C(1,J),RJI)
        RJI(LLL) = Zero
        Call ASub(3,C(1,K),C(1,J),RJK)
        RJK(LLL) = Zero
        Call VProd(VP,R,RJI)
        SJI = SProd(3,VP,VP)
        If(SJI.gt.Eps) then
          SJI = One/SJI
        else
          SJI = Zero
          endIf
        Call AScale(3,SJI,VP,VP)
        DB(1,1,1,1,ICor) = + Two*VP(1)*(R(2)*VP(3)-R(3)*VP(2))
        DB(1,2,1,1,ICor) = - DB(1,1,1,1,ICor)
        DB(2,1,1,1,ICor) = + Two*VP(2)*(R(2)*VP(3)-R(3)*VP(2))+R(3)*SJI
        DB(2,2,1,1,ICor) = - DB(2,1,1,1,ICor)
        DB(3,1,1,1,ICor) = + Two*VP(3)*(R(2)*VP(3)-R(3)*VP(2))-R(2)*SJI
        DB(3,2,1,1,ICor) = - DB(3,1,1,1,ICor)
        DB(1,1,2,1,ICor) = + Two*VP(1)*(R(3)*VP(1)-R(1)*VP(3))-R(3)*SJI
        DB(1,2,2,1,ICor) = - DB(1,1,2,1,ICor)
        DB(2,1,2,1,ICor) = + Two*VP(2)*(R(3)*VP(1)-R(1)*VP(3))
        DB(2,2,2,1,ICor) = - DB(2,1,2,1,ICor)
        DB(3,1,2,1,ICor) = + Two*VP(3)*(R(3)*VP(1)-R(1)*VP(3))+R(1)*SJI
        DB(3,2,2,1,ICor) = - DB(3,1,2,1,ICor)
        DB(1,1,3,1,ICor) = + Two*VP(1)*(R(1)*VP(2)-R(2)*VP(1))+R(2)*SJI
        DB(1,2,3,1,ICor) = - DB(1,1,3,1,ICor)
        DB(2,1,3,1,ICor) = + Two*VP(2)*(R(1)*VP(2)-R(2)*VP(1))-R(1)*SJI
        DB(2,2,3,1,ICor) = - DB(2,1,3,1,ICor)
        DB(3,1,3,1,ICor) = + Two*VP(3)*(R(1)*VP(2)-R(2)*VP(1))
        DB(3,2,3,1,ICor) = - DB(3,1,3,1,ICor)
        Call AMove(3,VP,B(1,1,ICor))
        Call VProd(VP,RJK,R)
        SJK = SProd(3,VP,VP)
        If(SJK.gt.Eps) then
          SJK = One/SJK
        else
          SJK = Zero
          endIf
        Call AScale(3,SJK,VP,VP)
        DB(1,3,1,3,ICor) = - Two*VP(1)*(R(2)*VP(3)-R(3)*VP(2))
        DB(1,2,1,3,ICor) = - DB(1,3,1,3,ICor)
        DB(2,3,1,3,ICor) = - Two*VP(2)*(R(2)*VP(3)-R(3)*VP(2))-R(3)*SJK
        DB(2,2,1,3,ICor) = - DB(2,3,1,3,ICor)
        DB(3,3,1,3,ICor) = - Two*VP(3)*(R(2)*VP(3)-R(3)*VP(2))+R(2)*SJK
        DB(3,2,1,3,ICor) = - DB(3,3,1,3,ICor)
        DB(1,3,2,3,ICor) = - Two*VP(1)*(R(3)*VP(1)-R(1)*VP(3))+R(3)*SJK
        DB(1,2,2,3,ICor) = - DB(1,3,2,3,ICor)
        DB(2,3,2,3,ICor) = - Two*VP(2)*(R(3)*VP(1)-R(1)*VP(3))
        DB(2,2,2,3,ICor) = - DB(2,3,2,3,ICor)
        DB(3,3,2,3,ICor) = - Two*VP(3)*(R(3)*VP(1)-R(1)*VP(3))-R(1)*SJK
        DB(3,2,2,3,ICor) = - DB(3,3,2,3,ICor)
        DB(1,3,3,3,ICor) = - Two*VP(1)*(R(1)*VP(2)-R(2)*VP(1))-R(2)*SJK
        DB(1,2,3,3,ICor) = - DB(1,3,3,3,ICor)
        DB(2,3,3,3,ICor) = - Two*VP(2)*(R(1)*VP(2)-R(2)*VP(1))+R(1)*SJK
        DB(2,2,3,3,ICor) = - DB(2,3,3,3,ICor)
        DB(3,3,3,3,ICor) = - Two*VP(3)*(R(1)*VP(2)-R(2)*VP(1))
        DB(3,2,3,3,ICor) = - DB(3,3,3,3,ICor)
        Call AMove(3,VP,B(1,3,ICor))
        Call AAdd(3,B(1,1,ICor),B(1,3,ICor),B(1,2,ICor))
        Call ANeg(3,B(1,2,ICor),B(1,2,ICor))
        Do 10 II = 1, 3
          Call ANeg(3,DB(1,1,II,1,ICor),DB(1,1,II,2,ICor))
          Call AMove(3,DB(1,1,II,1,ICor),DB(1,2,II,2,ICor))
          Call ANeg(3,DB(1,3,II,3,ICor),DB(1,3,II,2,ICor))
          Call AAdd(3,DB(1,3,II,3,ICor),DB(1,2,II,2,ICor),
     $                                  DB(1,2,II,2,ICor))
          Call AClear(3,DB(1,1,II,3,ICor))
          Call AClear(3,DB(1,3,II,1,ICor))
   10     Continue
      else if(LLL.eq.1) then
C
C     This is the sum of two bond angles
C
        IP(1) = 1
        IP(2) = 2
        IP(3) = 4
        Call DBBend(1,I,J,L,BB,IBB,C,DBB)
        Call AMove(6,BB,B(1,1,ICor))
        Call AMove(3,BB(1,3),B(1,4,ICor))
        Do 20 II = 1, 3
          Do 20 IC = 1, 3
            Do 20 JJ = 1, 3
              Call AMove(3,DBB(1,JJ,IC,II),DB(1,IP(JJ),IC,IP(II),ICor))
   20         Continue
        IP(1) = 3
        Call DBBend(1,K,J,L,BB,IBB,C,DBB)
        Call AMove(3,BB,B(1,3,ICor))
        Call AAdd(3,BB(1,2),B(1,2,ICor),B(1,2,ICor))
        Call AAdd(3,BB(1,3),B(1,4,ICor),B(1,4,ICor))
        Do 30 II = 1, 3
          Do 30 IC = 1, 3
            Do 30 JJ = 1, 3
              Call AAdd(3,DBB(1,JJ,IC,II),DB(1,IP(JJ),IC,IP(II),ICor),
     $                                    DB(1,IP(JJ),IC,IP(II),ICor))
   30         Continue
      else if(LLL.eq.2) then
C
C     This is a dihedral angle combination
C
        IP(1) = 1
        IP(2) = 2
        IP(3) = 4
        IP(4) = 3
        Call DBTors(1,I,J,L,K,BB,IBB,C,DBB)
        Do 40 II = 1, 4
          Call AMove(3,BB(1,II),B(1,IP(II),ICor))
          Do 40 IC = 1, 3
            Do 40 JJ = 1, 4
              Call AMove(3,DBB(1,JJ,IC,II),DB(1,IP(JJ),IC,IP(II),ICor))
   40         Continue
        IP(2) = 4
        IP(3) = 2
        Call DBTors(1,I,L,J,K,BB,IBB,C,DBB)
        Do 50 II = 1, 4
          Call ASub(3,B(1,IP(II),ICor),BB(1,II),B(1,IP(II),ICor))
          Do 50 IC = 1, 3
            Do 50 JJ = 1, 4
              Call ASub(3,DB(1,IP(JJ),IC,IP(II),ICor),DBB(1,JJ,IC,II),
     $          DB(1,IP(JJ),IC,IP(II),ICor))
   50         Continue
        Call AScale(12,Half,B(1,1,ICor),B(1,1,ICor))
        Call AScale(144,Half,DB(1,1,1,1,ICor),DB(1,1,1,1,ICor))
        endIf
      Return
      End
*Deck DBOOp2
      Subroutine DBOOp2(I,J,K,L,B,IB,C,DB,CosT,SinT)
C
C     Out-of-plane B matrix
C     Wilson, Decius, and Cross, pp. 58-60
C
C     I is the central atom, J and K the in-plane atoms,
C     and L the out-of-plane atom.
C
      Implicit Real*8(A-H,O-Z)
      Real*8 NProd,NrmRij,NrmRik,NrmRil,NrmVP1,mOne
      Dimension B(3,4), C(3,*), IB(*), Rij(3), Rik(3), Ril(3),
     $  VP1(3), DB(3,4,3,4), dnvp1(3,4),dnril(3,4),
     $  DDNVP(3,4,3,4),DDNIL(3,4,3,4), dVP1(3,4,3),
     $  dRij(3,4,3), dRik(3,4,3), dRil(3,4,3), ddVP1(3,4,3,4,3),
     $  dfdX(3,4), TMPij(3),Tmpik(3),Tmp1(3),Tmp2(3), TmpijX(3),
     $  TmpijY(3), TmpikX(3), TmpikY(3)
      Save mOne, Zero, One, Three, Tol, dRij, dRik, dRil
      Data mOne/-1.0d0/, Zero/0.d0/, One/1.d0/, Three/3.d0/,Tol/1.d-10/
      Data dRij/-1.0d0,2*0.0d0,1.0d0,9*0.0d0,-1.0d0,2*0.d0,1.0d0,
     $  9*0.0d0,-1.0d0,2*0.0d0,1.0d0,6*0.0d0/
      Data dRik/-1.0d0,5*0.0d0,1.0d0,6*0.0d0,-1.0d0,5*0.d0,1.0d0,
     $  6*0.0d0,-1.0d0,5*0.0d0,1.0d0,3*0.0d0/
      Data dRil/-1.0d0,8*0.0d0,1.0d0,3*0.0d0,-1.0d0,8*0.d0,1.0d0,
     $  3*0.0d0,-1.0d0,8*0.0d0,1.0d0/
C
C     ==================================================================
C       T = | F(x) | = sqrt[ (F(x)**2) ]
C
C       F(x) = Arcos( f(X1,X2,X3,X4) - ASin(1.0)
C
C                    1         1        -->   -->
C         with f =  ----   *  ----   * (Ril . VP1)
C                  ||Ril||   ||Vp1||
C
C              -->   -->    -->
C         and  VP1 = Rij /\ Rik
C     ==================================================================
C       d        1        /  d      \     /      F(x)         \
C      --  T  =  - . 2 . |  -- F(x)  | . |  -----------------  |
C    dX        2        \ dX      /     \ sqrt[ (F(x)**2) ] /
C     ==================================================================
C       d  d       /  d  d       \     /      F(x)         \
C      -- --  T = |  -- --  F(X)  | . |  -----------------  |
C      dY dX       \ dY dX       /     \ sqrt[ (F(x)**2) ] /
C
C                   /  d      \     / /        1           d      \     /  d             \     /      1       F(x)          d      \
C               +  |  -- F(x)  | . | |  --------------- . --  F(x) | + |  -- F(x) . F(x) |  . | 2 . - - .  ------------  . -- F(x)  |
C                   \ dX      /     \ \ sqrt[(F(x)**2)]   dY      /     \ dX             /     \      2    [F(x)**2]**3/2  dY      /
C     ==================================================================
C       d                 /        -1        \   d
C       --   F(X)       = | ---------------- | * --  f(X)
C       dX                \ sqrt ( 1 - f**2) /   dx,y,z,(i,j,k,l)
C     ==================================================================
C       df     -->   -->      1    d    1            1    d     1
C       --  = (Ril . VP1) [ ---- . -- ----     +   ---- . --  ----   ]
C       dX                 ||Ril|| dX ||Vp1||     ||Vp1|| dx  ||Ril||
C
C                1         1      -->   d  -->   -->   d   -->
C           +  ----   *  ----   [ Ril . -- VP1 + VP1 . --  Ril ]
C             ||Ril||   ||Vp1||         dX             dX
C     ==================================================================
C    d    d           d  /  /        -1        \   df               \
C    --  --  (F(x) = --  |  | ---------------- | * --               |
C    dY  dX          dY  \  \ sqrt ( 1 - f**2) /   dx,y,z,(i,j,k,l) /
C
C                   /        1         /   /   d   d       \
C    =              | ---------------- | * |  --  --  f(X) |
C                   \ sqrt ( 1 - f**2) /   \  dY   dX      /
C
C       / df  \   /  1                                   df   \
C    +  | --  | . |  - . (( 1 - f**2)**-3/2) . (-2 . f . --   |
C       \ dX  /   \  2                                   dY  /
C     ==================================================================
C     Alternatively B and DB can be calculated using the OOPl and
C     DBOOPl routines defined for Wilson angles with the following
C     alteration ;
C     if(sint.le.zero)
C      B(ifox,jfox)=-B(ifox,jfox)
C     DB(ICoor1,IAt1,ICoor2,IAt2) = -DB(ICoor1,IAt1,ICoor2,IAt2)
C     ==================================================================
C     DBOOP2 and OOPL2 could be a bit compacted using intermediate
C     arrays. See the role of Q1 in e.g. DBDI, DBDIP.
C     ==================================================================
C
      IB(1) = I
      IB(2) = J
      IB(3) = K
      IB(4) = L
      Do 10 KK = 1,3
        Rij(KK) = C(KK,J) - C(KK,I)
        Rik(KK) = C(KK,K) - C(KK,I)
        Ril(KK) = C(KK,L) - C(KK,I)
   10 Continue
      NrmRij = Sqrt(SProd(3,Rij,Rij))
      NrmRil = Sqrt(SProd(3,Ril,Ril))
      NrmRik = Sqrt(SProd(3,Rik,Rik))
      If(NrmRij.ge.Tol.and.NrmRil.ge.Tol.and.NrmRik.ge.Tol) then
        Call VProd(VP1,Rij,Rik)
        NrmVP1 = Sqrt(SProd(3,VP1,VP1))
        NProd = SProd(3,Ril,VP1)/(NrmVP1*NrmRil)
        ff = (ACos(NPRod) - Asin(One))
        T = Abs(ff)
        Signe = ff/T
        CosT = Cos(T)
        SinT = Sin(T)
        Fac1 = mOne/(Sqrt(One -(NProd*NProd)))
        Fac2 = Fac1*NPRod/(One -(NProd*NProd))
        Do 40 ICoor=1,3
          Do 40 IAt=1,4
            Do 20 IFox =1,3
              tmpij(IFox) = dRij(ICoor,IAt,IFox)
              tmpik(IFox) = dRik(ICoor,IAt,IFox)
   20       Continue
            Call VProd(Tmp1,Rij,Tmpik)
            Call VProd(Tmp2,Tmpij,Rik)
            Do 30 IFox=1,3
             dVP1(ICoor,IAt,IFox) = Tmp1(IFox) + Tmp2(Ifox)
   30       Continue
   40       Continue
        Do 70 ICoor1=1,3
          Do 70 IAt1=1,4
            Do 70 ICoor2=1,3
              Do 70 IAt2=1,4
                Do 50 IFox =1,3
                  tmpijY(IFox) = dRij(ICoor1,IAt1,IFox)
                  tmpikY(IFox) = dRik(ICoor1,IAt1,IFox)
                  tmpijX(IFox) = dRij(ICoor2,IAt2,IFox)
                  tmpikX(IFox) = dRik(ICoor2,IAt2,IFox)
   50           Continue 
                Call VProd(Tmp1,TmpijY,TmpikX)
                Call VProd(Tmp2,TmpijX,tmpikY)
                Do 60 IFox=1,3
                  ddVP1(ICoor1,IAt1,ICoor2,IAt2,IFox) = Tmp1(IFox) +
     $              Tmp2(Ifox)
   60           Continue
   70           Continue
        Do 90 ICoor=1,3
          Do 90 IAt=1,4
            t = Zero
            Do 80 IFox=1,3
              t = t + dRil(ICoor,IAt,IFox) * Ril(Ifox)
   80       Continue
            dnril(ICoor,IAt) = -t / (NrmRil*NrmRil*NrmRil)
   90       Continue
        Do 110 ICoor1=1,3
          Do 110 IAt1=1,4
            Do 110 ICoor2=1,3
              Do 110 IAt2=1,4
                t1 = Zero
                t2 = Zero
                Do 100 IFox=1,3
                  t1 = t1 + Ril(IFox) * dRil(ICoor2,IAt2,IFox)
                  t2 = t2 + dRil(ICoor1,IAt1,IFox) *
     $              dRil(ICoor2,IAt2,IFox)
  100           Continue 
                DDNIL(ICoor1,IAt1,ICoor2,IAt2) =
     $            -t2 / (NrmRil*NrmRil*NrmRil) -
     $            t1  * Three * dnril(ICoor1,IAt1) / (NrmRil*NrmRil)
  110  Continue
        Do 130 ICoor1=1,3
         Do 130 IAt1=1,4
          t = Zero
          Do 120 ifox=1,3
           t = t + ( dVP1(ICoor1,IAt1,IFox) * VP1(IFox) )
  120     continue
          dnvp1(ICoor1,IAt1)= -t / (NrmVP1*NrmVP1*NrmVP1)
  130   continue 
        Do 150 ICoor1=1,3
          Do 150 IAt1=1,4
            Do 150 ICoor2=1,3
              Do 150 IAt2=1,4
                t1 = Zero
                t2 = Zero
                Do 140 IFox=1,3
                  t1 = t1 + ( VP1(IFox) * dVP1(ICoor2,IAt2,IFox) )
  140             t2 = t2 +
     $              dVP1(ICoor1,IAt1,IFox)*dVP1(ICoor2,IAt2,IFox) +
     $              VP1(IFox) * ddVP1(ICoor1,IAt1,ICoor2,Iat2,IFox)
                t =  ( Three * dnvp1(ICoor1,IAt1) ) / (NrmVP1*NrmVP1)
  150           DDNVP(ICoor1,IAt1,ICoor2,IAt2) =
     $            -t2 / (NrmVP1*NrmVP1*NrmVP1) -
     $            t1 * Three * dnvp1(ICoor1,IAt1) / (NrmVP1*NrmVP1)
C
        Do 170 ICoor=1,3
          Do 170 IAt =1,4
            t1 = SProd(3,Ril,VP1) * (dnvp1(ICoor,Iat)/NrmRil +
     $        dnril(ICoor,Iat)/NrmVP1)
            t2 = Zero
            Do 160 IFox = 1,3
  160         t2 = t2 + (1/(NrmRil*NrmVP1)) *
     $        (Ril(IFox)*dVP1(ICoor,Iat,IFox) +
     $        VP1(IFox)*dRil(ICoor,Iat,IFox) )
            dfdx(ICoor,Iat) = t1 + t2
  170       B(ICoor,IAt) = Signe * Fac1 * dfdx(ICoor,IAt)
C
        Do 190 ICoor1=1,3
          Do 190 IAt1=1,4
            Do 190 ICoor2=1,3
              Do 190 IAt2 =1,4
                t1 = Zero
                t2 = Zero
                t3 = Zero
                t4 = Zero
                Do 180 ifox=1,3
                  t1 = t1 + Ril(IFox) * dVP1(ICoor1,Iat1,IFox) +
     $              VP1(IFox) * dRil(ICoor1,Iat1,IFox)
                  t2 = t2 + Ril(IFox) * VP1(IFox)
                  t3 = t3 + Ril(IFox) * dVP1(ICoor2,Iat2,IFox) +
     $              VP1(IFox) * dRil(ICoor2,Iat2,IFox)
  180             t4 = t4 +
     $              dRil(ICoor1,Iat1,IFox) * dVP1(ICoor2,Iat2,IFox) +
     $              Ril(IFox) * ddVP1(Icoor1,IAt1,ICoor2,IAt2,Ifox) +
     $              dVP1(ICoor1,Iat1,IFox) * dRil(ICoor2,IAt2,IFox)
                tt1 = dnvp1(ICoor2,Iat2) / NrmRil +
     $            dnril(ICoor2,Iat2) / NrmVP1
                tt2 = dnril(ICoor1,Iat1) * dnvp1(ICoor2,Iat2) +
     $            DDNVP(ICoor1,IAt1,ICoor2,IAt2) / NrmRil +
     $            dnvp1(ICoor1,Iat1) * dnril(ICoor2,IAt2) +
     $            DDNIL(Icoor1,Iat1,Icoor2,Iat2) / NrmVP1
                tt3 = dnvp1(ICoor1,Iat1) / NrmRil +
     $            dnril(ICoor1,Iat1) / NrmVP1
                tt4 = NrmVP1 * NrmRil
                t = t1 * tt1 + t2 * tt2 + t3 * tt3 + t4 / tt4
  190           DB(icoor1,iat1,icoor2,iat2) = signe * ( Fac1 * t +
     $            Fac2 * dfdx(ICoor1,IAt1) * dfdx(ICoor2,IAt2) )
        endIf
      Return
      End
*Deck DBOOPl
      Subroutine DBOOPl(I,J,K,L,B,IB,C,DB,CosT,SinT)
      Implicit Real*8(A-H,O-Z)
C
C     Out-of-plane first and second order B matrices
C     Wilson, Decius, and Cross, pp. 58-60
C     See also Tuzun et al: JCC 18 1997 pp. 1804
C
C     I is the central atom, J and K the in-plane atoms,
C     and L the out-of-plane atom.
C
      Dimension B(3,4), C(3,*), IB(*), VP(3,4),
     $ DB(3,4,3,4), BP(3,4), N(3,3), NRC(3,3),
     $ IBP(4), R(3,4), D2(4), D(4), DE(3,4)
      Save Zero, One, Tol, N, NRC
      Data Zero/0.d0/, One/1.d0/, Tol/1.d-10/
      Data ((N(N1,N2),N1=1,3),N2=1,3)/0,1,-1,-1,0,1,1,-1,0/
      Data ((NRC(N1,N2),N1=1,3),N2=1,3)/0,3,2,3,0,1,2,1,0/
C
      IB(1) = I
      IB(2) = J
      IB(3) = K
      IB(4) = L
      Do 15 KK = 2,4
        D2(KK) = Zero
        Do 10 LL = 1,3
          R(LL,KK) = C(LL,IB(KK)) - C(LL,I)
   10     D2(KK) = D2(KK) + R(LL,KK)*R(LL,KK)
   15   D(KK) = Sqrt(D2(KK))
      If(Min(D(2),D(3),D(4)).ge.Tol) then
        Call VProd(VP(1,2),R(1,3),r(1,4))
        Call VProd(VP(1,3),R(1,4),r(1,2))
        Call VProd(VP(1,4),R(1,2),r(1,3))
        RNorm = D(2)*D(3)*D(4)
        E = SProd(3,R(1,4),VP(1,4)) / RNorm
        CosP = SProd(3,R(1,2),R(1,3)) / (D(2)*D(3))
        SinP = Sqrt(One - CosP*CosP)
        If(Abs(SinP).ge.Tol) then
          SinT = E / SinP
          CosT = Sqrt(One - SinT*SinT)
          If(Abs(CosT).ge.Tol) then
            Call AClear(12,BP)
            Call AClear(144,DB)
            Call DBBend(1,J,I,K,BP,IBP,C,DB)
            Do 16 KK = 1,3
              ATemp = BP(KK,1)
              BP(KK,1) = BP(KK,2)
              BP(KK,2) = ATemp
              Do 16 LL = 1,3
                Do 16 MM = 1,3
                  ATemp = DB(MM,1,LL,KK)
                  DB(MM,1,LL,KK) = DB(MM,2,LL,KK)
                  DB(MM,2,LL,KK) = ATemp
   16             Continue
            Do 17 KK = 1,3
              Do 17 LL = 1,3
                Do 17 MM = 1,3
                  ATemp = DB(MM,KK,LL,1)
                  DB(MM,KK,LL,1) = DB(MM,KK,LL,2)
                  DB(MM,KK,LL,2) = ATemp
   17             Continue
            T1 = COSP*SINT
            T2 = COST*SINP
            T3 = SINT*SINP
            T4 = COST*COSP
            Do 25 KK = 1,3
              Do 20 LL = 2,4
                DE(KK,LL) = VP(KK,LL)/RNorm - R(KK,LL)*E/D2(LL)
   20           B(KK,LL) = (DE(KK,LL) - T1*BP(KK,LL)) / T2
              DE(KK,1) = -(DE(KK,2) + DE(KK,3) + DE(KK,4))
   25         B(KK,1) = -(B(KK,2) + B(KK,3) + B(KK,4))
           Call AScale(144,-T1,DB,DB)
           Do 30 KK = 1,3
             Do 30 LL = 2,4
               DB(KK,LL,KK,LL) = DB(KK,LL,KK,LL) - E/D2(LL)
               Do 30 MM = 1,3
                 DB(KK,LL,MM,LL) = DB(KK,LL,MM,LL) +
     $             2*R(KK,LL)*R(MM,LL)*E/(D2(LL)*D2(LL))
                 Do 30 NN = 2,4
                   DB(KK,LL,MM,NN) = DB(KK,LL,MM,NN) -
     $               DE(MM,NN)*R(KK,LL)/D2(LL) -
     $               VP(KK,LL)*R(MM,NN)/(RNorm*D2(NN)) +
     $               T3*(B(KK,LL)*B(MM,NN)+BP(KK,LL)*BP(MM,NN)) -
     $               T4*(B(KK,LL)*BP(MM,NN)+BP(KK,LL)*B(MM,NN))
                   If(LL.ne.NN.and.KK.ne.MM)
     $               DB(KK,LL,MM,NN) = DB(KK,LL,MM,NN) +
     $               N(KK,MM)*N(LL-1,NN-1)*
     $               R(NRC(KK,MM),1+NRC(LL-1,NN-1))/RNorm
   30              Continue
            Do 40 LL = 2,4
              Do 40 KK = 1,3
                Do 40 MM = 1,3
   40             DB(KK,1,MM,LL) = -(DB(KK,2,MM,LL) +
     $              DB(KK,3,MM,LL) + DB(KK,4,MM,LL))
            Do 50 LL = 1,4
              Do 50 KK = 1,3
                Do 50 MM = 1,3
   50             DB(KK,LL,MM,1) = -(DB(KK,LL,MM,2) +
     $              DB(KK,LL,MM,3) + DB(KK,LL,MM,4))
            Call AScale(144,One/T2,DB,DB)
            endIf
          endIf
        endIf
      Return
      End
*Deck DBStr1
      Subroutine DBStr1(I,J,B,IB,C,DB)
      Implicit Real*8(A-H,O-Z)
C
C     Compute stretch matrix elements for Wilson B-matrix and their
C     cartesian first derivatives.
C
C     Philippe Y. Ayala /Oct. 93
C
      Dimension B(3,2), DB(3,2,3,2), IB(2), C(3,*)
      Save One
      Data One/1.0d0/
C
      IB(1) = I
      IB(2) = J
      DX = C(1,J) - C(1,I)
      DY = C(2,J) - C(2,I)
      DZ = C(3,J) - C(3,I)
      DIJI = One / Sqrt(DX*DX+DY*DY+DZ*DZ)
      B(1,2) = DX*DIJI
      B(1,1) = -B(1,2)
      B(2,2) = DY*DIJI
      B(2,1) = -B(2,2)
      B(3,2) = DZ*DIJI
      B(3,1) = -B(3,2)
      DB(1,1,1,1) = (One-B(1,1)*B(1,1))*DIJI
      DB(1,1,1,2) = -DB(1,1,1,1)
      DB(2,1,2,1) = (One-B(2,1)*B(2,1))*DIJI
      DB(2,1,2,2) = -DB(2,1,2,1)
      DB(3,1,3,1) = (One-B(3,1)*B(3,1))*DIJI
      DB(3,1,3,2) = -DB(3,1,3,1)
      DB(2,1,1,2) = B(2,1)*B(1,1)*DIJI
      DB(1,1,2,2) = DB(2,1,1,2)
      DB(1,1,2,1) = -DB(2,1,1,2)
      DB(3,1,1,2) = B(3,1)*B(1,1)*DIJI
      DB(1,1,3,2) = DB(3,1,1,2)
      DB(1,1,3,1) = -DB(3,1,1,2)
      DB(3,1,2,2) = B(3,1)*B(2,1)*DIJI
      DB(2,1,3,2) = DB(3,1,2,2)
      DB(2,1,3,1) = -DB(3,1,2,2)
      DB(1,2,1,2)=-DB(1,1,1,2)
      DB(1,2,2,2)=-DB(1,1,2,2)
      DB(2,2,2,2)=-DB(2,1,2,2)
      DB(1,2,3,2)=-DB(1,1,3,2)
      DB(2,2,3,2)=-DB(2,1,3,2)
      DB(3,2,3,2)=-DB(3,1,3,2)
      Return
      End
*Deck DBStr
      Subroutine DBStr(NoInt,I,J,B,IB,C,DB)
      Implicit Real*8(A-H,O-Z)
C
C     Compute stretch matrix elements for Wilson B-matrix and their
C     cartesian first derivatives.
C
C     Philippe Y. Ayala /Oct. 93
C
      Dimension B(3,4,*), DB(3,4,3,4,*), IB(4,*), C(3,*)
      Save One
      Data One/1.0d0/
C
      IB(1,NoInt) = I
      IB(2,NoInt) = J
      IB(3,NoInt) = 0
      IB(4,NoInt) = 0
      DX = C(1,J) - C(1,I)
      DY = C(2,J) - C(2,I)
      DZ = C(3,J) - C(3,I)
      DIJI = One / Sqrt(DX*DX+DY*DY+DZ*DZ)
      B(1,2,NoInt) = DX*DIJI
      B(1,1,NoInt) = -B(1,2,NoInt)
      B(2,2,NoInt) = DY*DIJI
      B(2,1,NoInt) = -B(2,2,NoInt)
      B(3,2,NoInt) = DZ*DIJI
      B(3,1,NoInt) = -B(3,2,NoInt)
      DB(1,1,1,1,NoInt) = (One-B(1,1,NoInt)*B(1,1,NoInt))*DIJI
      DB(1,1,1,2,NoInt) = -DB(1,1,1,1,NoInt)
      DB(2,1,2,1,NoInt) = (One-B(2,1,NoInt)*B(2,1,NoInt))*DIJI
      DB(2,1,2,2,NoInt) = -DB(2,1,2,1,NoInt)
      DB(3,1,3,1,NoInt) = (One-B(3,1,NoInt)*B(3,1,NoInt))*DIJI
      DB(3,1,3,2,NoInt) = -DB(3,1,3,1,NoInt)
      DB(2,1,1,2,NoInt) = B(2,1,NoInt)*B(1,1,NoInt)*DIJI
      DB(1,1,2,2,NoInt) = DB(2,1,1,2,NoInt)
      DB(2,1,1,1,NoInt) = -DB(2,1,1,2,NoInt)
      DB(1,1,2,1,NoInt) = -DB(2,1,1,2,NoInt)
      DB(3,1,1,2,NoInt) = B(3,1,NoInt)*B(1,1,NoInt)*DIJI
      DB(1,1,3,2,NoInt) = DB(3,1,1,2,NoInt)
      DB(3,1,1,1,NoInt) = -DB(3,1,1,2,NoInt)
      DB(1,1,3,1,NoInt) = -DB(3,1,1,2,NoInt)
      DB(3,1,2,2,NoInt) = B(3,1,NoInt)*B(2,1,NoInt)*DIJI
      DB(2,1,3,2,NoInt) = DB(3,1,2,2,NoInt)
      DB(3,1,2,1,NoInt) = -DB(3,1,2,2,NoInt)
      DB(2,1,3,1,NoInt) = -DB(3,1,2,2,NoInt)
      DB(1,2,1,1,NoInt)=-DB(1,1,1,1,NoInt)
      DB(2,2,1,1,NoInt)=-DB(2,1,1,1,NoInt)
      DB(3,2,1,1,NoInt)=-DB(3,1,1,1,NoInt)
      DB(1,2,2,1,NoInt)=-DB(1,1,2,1,NoInt)
      DB(2,2,2,1,NoInt)=-DB(2,1,2,1,NoInt)
      DB(3,2,2,1,NoInt)=-DB(3,1,2,1,NoInt)
      DB(1,2,3,1,NoInt)=-DB(1,1,3,1,NoInt)
      DB(2,2,3,1,NoInt)=-DB(2,1,3,1,NoInt)
      DB(3,2,3,1,NoInt)=-DB(3,1,3,1,NoInt)
      DB(1,2,1,2,NoInt)=-DB(1,1,1,2,NoInt)
      DB(2,2,1,2,NoInt)=-DB(2,1,1,2,NoInt)
      DB(3,2,1,2,NoInt)=-DB(3,1,1,2,NoInt)
      DB(1,2,2,2,NoInt)=-DB(1,1,2,2,NoInt)
      DB(2,2,2,2,NoInt)=-DB(2,1,2,2,NoInt)
      DB(3,2,2,2,NoInt)=-DB(3,1,2,2,NoInt)
      DB(1,2,3,2,NoInt)=-DB(1,1,3,2,NoInt)
      DB(2,2,3,2,NoInt)=-DB(2,1,3,2,NoInt)
      DB(3,2,3,2,NoInt)=-DB(3,1,3,2,NoInt)
      Return
      End
*Deck DBTors
      Subroutine DBTors(NoInt,I,J,K,L,B,IB,C,DB)
      Implicit Real*8(A-H,O-Z)
C
C     Compute torsion matrix elements for Wilson B-matrix and their
C     cartesian first derivatives.
C
C     Philippe Y. Ayala /Oct. 93
C
      Dimension B(3,4,*),IB(4,*),C(*),DB(3,4,3,4,*),DCR1(3,3,4),
     $  DCR2(3,3,4),RIJ(3),RJK(3),RKL(3),EIJ(3),EJK(3),EKL(3),CR1(3),
     $  CR2(3),BB(3,4,5),IBB(4,5),DBB(3,4,3,4,3)
      Integer Alpha,Beta
      Save Zero,One,Two
      Data Zero/0.d0/,One/1.d0/,Two/2.0d0/
C
      IAIND=3*(I-1)
      JAIND=3*(J-1)
      KAIND=3*(K-1)
      LAIND=3*(L-1)
      IB(1,NoInt)=I
      IB(2,NoInt)=J
      IB(3,NoInt)=K
      IB(4,NoInt)=L
      DIJSQ=Zero
      DJKSQ=Zero
      DKLSQ=Zero
      Do 124 M=1,3
        RIJ(M)=C(M+JAIND)-C(M+IAIND)
        DIJSQ=DIJSQ+RIJ(M)**2
        RJK(M)=C(M+KAIND)-C(M+JAIND)
        DJKSQ=DJKSQ+RJK(M)**2
        RKL(M)=C(M+LAIND)-C(M+KAIND)
  124   DKLSQ=DKLSQ+RKL(M)**2
      DIJ =Sqrt(DIJSQ)
      DJK =Sqrt(DJKSQ)
      DKL =Sqrt(DKLSQ)
      Do 136 M=1,3
        EIJ(M)=RIJ(M)/DIJ
        EJK(M)=RJK(M)/DJK
  136   EKL(M)=RKL(M)/DKL
C
C     Compute matrix elements needed for partial first derivatives
C
      Call AClear(60,BB)
      Call AClear(432,DBB)
      Call DBStr(1,I,J,BB,IBB,C,DBB)
      Call DBStr(2,J,K,BB,IBB,C,DBB)
      Call DBStr(3,K,L,BB,IBB,C,DBB)
      Call Bend(4,I,J,K,BB,IBB,C)
      Call Bend(5,J,K,L,BB,IBB,C)
      NoIntIJ=1
      NoIntJK=2
      NoIntKL=3
      Alpha=4
      Beta=5
C
      CR1(1)=EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2)=EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3)=EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
      CR2(1)=EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2)=EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3)=EJK(1)*EKL(2)-EJK(2)*EKL(1)
      DOTPJ = -(EIJ(1)*EJK(1)+EIJ(2)*EJK(2)+EIJ(3)*EJK(3))
      DOTPK = -(EJK(1)*EKL(1)+EJK(2)*EKL(2)+EJK(3)*EKL(3))
      SINPJ =Sqrt(ONE-DOTPJ**2)
      SINPK =Sqrt(ONE-DOTPK**2)
C
C     Compute B-matrix elements
C
      Do 164 M=1,3
        SMI=-CR1(M)/(DIJ*SINPJ*SINPJ)
        B(M,1,NoInt)=SMI
        F1=(CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
        F2=(DOTPK*CR2(M))/(DJK*SINPK*SINPK)
        SMJ=F1-F2
        B(M,2,NoInt)=SMJ
        SML= CR2(M)/(DKL*SINPK*SINPK)
        B(M,4,NoInt)=SML
  164   B(M,3,NoInt)=(-SMI-SMJ-SML)
C
C     Compute partial derivatives of cross-products
C
      Do 200 M=1,3
        Do 200 LL=1,4
          If(LL.eq.1) then
            LIJ=1
            LJK=3
            LKL=3
          else if(LL.eq.2) then
            LIJ=2
            LJK=1
            LKL=3
          else if(LL.eq.3) then
            LIJ=3
            LJK=2
            LKL=1
          else if(LL.eq.4) then
            LIJ=3
            LJK=3
            LKL=2
            endIf
          DCR1(1,M,LL)=(DBB(2,2,M,LIJ,NoIntIJ)*BB(3,2,NoIntJK)
     $      + BB(2,2,NoIntIJ)*DBB(3,2,M,LJK,NoIntJK)
     $      - DBB(3,2,M,LIJ,NoIntIJ)*BB(2,2,NoIntJK)
     $      - BB(3,2,NoIntIJ)*DBB(2,2,M,LJK,NoIntJK))
          DCR2(1,M,LL)=(DBB(2,2,M,LJK,NoIntJK)*BB(3,2,NoIntKL)
     $      + BB(2,2,NoIntJK)*DBB(3,2,M,LKL,NoIntKL)
     $      - DBB(3,2,M,LJK,NoIntJK)*BB(2,2,NoIntKL)
     $      - BB(3,2,NoIntJK)*DBB(2,2,M,LKL,NoIntKL))
          DCR1(2,M,LL)=(DBB(3,2,M,LIJ,NoIntIJ)*BB(1,2,NoIntJK)
     $      + BB(3,2,NoIntIJ)*DBB(1,2,M,LJK,NoIntJK)
     $      - DBB(1,2,M,LIJ,NoIntIJ)*BB(3,2,NoIntJK)
     $      - BB(1,2,NoIntIJ)*DBB(3,2,M,LJK,NoIntJK))
          DCR2(2,M,LL)=(DBB(3,2,M,LJK,NoIntJK)*BB(1,2,NoIntKL)
     $      + BB(3,2,NoIntJK)*DBB(1,2,M,LKL,NoIntKL)
     $      - DBB(1,2,M,LJK,NoIntJK)*BB(3,2,NoIntKL)
     $      - BB(1,2,NoIntJK)*DBB(3,2,M,LKL,NoIntKL))
          DCR1(3,M,LL)=(DBB(1,2,M,LIJ,NoIntIJ)*BB(2,2,NoIntJK)
     $      + BB(1,2,NoIntIJ)*DBB(2,2,M,LJK,NoIntJK)
     $      - BB(2,2,NoIntIJ)*DBB(1,2,M,LJK,NoIntJK)
     $      - DBB(2,2,M,LIJ,NoIntIJ)*BB(1,2,NoIntJK))
  200     DCR2(3,M,LL)=(DBB(1,2,M,LJK,NoIntJK)*BB(2,2,NoIntKL)
     $      + BB(1,2,NoIntJK)*DBB(2,2,M,LKL,NoIntKL)
     $      - BB(2,2,NoIntJK)*DBB(1,2,M,LKL,NoIntKL)
     $      - DBB(2,2,M,LJK,NoIntJK)*BB(1,2,NoIntKL))
C
C     Compute partial derivatives of B-matrix elements
C
      Do 250 M=1,3
        Do 250 N=1,3
          Do 250 LL=1,4
            If(LL.eq.1) then
              LIJ=1
              LJK=3
              LKL=3
              LAlpha=1
              LBeta=4
            else if(LL.eq.2) then
              LIJ=2
              LJK=1
              LKL=3
              LAlpha=2
              LBeta=1
            else if(LL.eq.3) then
              LIJ=3
              LJK=2
              LKL=1
              LAlpha=3
              LBeta=2
            else if(LL.eq.4) then
              LIJ=3
              LJK=3
              LKL=2
              LAlpha=4
              LBeta=3
              endIf
            DB(M,1,N,LL,NoInt)=-(DCR1(M,N,LL)
     $        -CR1(M)*( BB(N,LIJ,NoIntIJ)/DIJ
     $        +TWO*BB(N,LAlpha,Alpha)*DOTPJ/SINPJ ))/(DIJ*SINPJ*SINPJ)
            DB(M,4,N,LL,NoInt)=(DCR2(M,N,LL)
     $        -CR2(M)*( BB(N,LKL,NoIntKL)/DKL
     $        + TWO*BB(N,LBeta,Beta)*DOTPK/SINPK ))/(DKL*SINPK*SINPK)
  250       DB(M,2,N,LL,NoInt)=( ( DCR1(M,N,LL)*(DJK-DIJ*DOTPJ)
     $        + CR1(M)*(BB(N,LJK,NoIntJK)-DOTPJ*BB(N,LIJ,NoIntIJ)
     $        + SINPJ*DIJ*BB(N,LAlpha,Alpha) )
     $        - CR1(M)*(DJK-DIJ*DOTPJ)*(BB(N,LJK,NoIntJK)/DJK
     $        + BB(N,LIJ,NoIntIJ)/DIJ
     $        +TWO*DOTPJ*BB(N,LAlpha,Alpha)/SINPJ) )
     $        / (DJK*DIJ*SINPJ*SINPJ))
     $        - (( DCR2(M,N,LL)*DOTPK - SINPK*BB(N,LBeta,Beta)*CR2(M)
     $        - DOTPK*CR2(M)*( BB(N,LJK,NoIntJK)/DJK
     $        + TWO*DOTPK*BB(N,LBeta,Beta)/SINPK ) )
     $        / (DJK*SINPK*SINPK))
      Do 350 M=1,3
        Do 350 N=1,3
          Do 350 KK=1,4
  350       DB(M,3,N,KK,NoInt)=-(DB(M,1,N,KK,NoInt)
     $        + DB(M,2,N,KK,NoInt) + DB(M,4,N,KK,NoInt) )
      Return
      End
*Deck OutAng
      Real*8 Function OutAng(CI,CJ,CK,CL)
      Implicit Real*8 (A-H,O-Z)
      Dimension CI(3),CJ(3),CK(3),CL(3) 
      Dimension RIJ(3),RIK(3),RIL(3),VP1(3)
C
C Computes the value of an out-of-plane angle
C i:1st, j:2nd, k:3rd, l:4th)
C the fourth atom out of the plane of the first 3
C the first atom is the central atom
C
      Tol=1.0d-3
      Do 10 KK = 1,3
        Rij(KK) = CJ(KK) - CI(KK)
        Rik(KK) = CK(KK) - CI(KK)
   10   Ril(KK) = CL(KK) - CI(KK)
      SNrmRij = Sqrt(SProd(3,Rij,Rij))
      SNrmRil = Sqrt(SProd(3,Ril,Ril))
      SNrmRik = Sqrt(SProd(3,Rik,Rik))
      If(SNrmRij.ge.Tol.and.SNrmRil.ge.Tol.and.SNrmRik.ge.Tol) then
        Call VProd(VP1,Rij,Rik)
        SNrmVP1 = Sqrt(SProd(3,VP1,VP1))
        SNProd = SProd(3,Ril,VP1)/(SNrmVP1*SNrmRil)
        ff = (ACos(SNPRod) - Asin(1.0d0))
        T = Abs(ff)
        Signe = ff/T
        CosT = Cos(T)
        SinT = Sqrt(1.0d0-CosT*CosT)
        OutAng=ff
       Else
        OutAng=0.0d0
       EndIf
       Return
       End
