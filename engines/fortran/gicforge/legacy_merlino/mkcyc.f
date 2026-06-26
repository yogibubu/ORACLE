*Deck MkCyc
      Subroutine MkCyc(IOut,IPrint,MxBnd,MxAtG,MxTrm,MxAtCy,MxCyc,
     $  NBond,IBond,NLen,NAng,NDih,NExpCy,IAtomB,IAtomA,IAtomD,
     $  NCyc,NatC,ICAt,IAtCyc,IBr,IAn,EAN)
      Implicit Real*8 (A-H,O-Z)
      Common/onec/N1Cyc,IAt1C(2,100)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
      Dimension NBond(*),IBond(MxBnd,*)
      Dimension IAtomB(MxAtG,MxTrm,*),IAtomA(MxAtG,MxTrm,*)
      Dimension IAtomD(MxAtG,MxTrm,*)
      Dimension NatC(*),ICAt(MxAtCy,*),IAtCyc(*),IAN(*)
      Dimension EAN(*)
      Logical ReNumb
C find cycles
      if(NCyc.lt.NExpCy) call Cy5(IOut,IPrint,MxBnd,MxAtG,MxTrm,NBond,
     $  IBond,NDih,IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,NAtC,ICAt,
     $  IAtCyc)
      if(NCyc.lt.NExpCy) call Cy6(IOut,IPrint,MxBnd,MxAtG,MxTrm,NBond,
     $  IBond,NDih,IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,NAtC,ICAt,
     $  IAtCyc)
      if(NCyc.lt.NExpCy) call Cy7(IOut,IPrint,MxAtG,MxAtG,MxTrm,MxTrm,
     $  NAng,NDih,NBond,IAtomA,IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,
     $  NAtC,ICAt,IAtCyc)
      if(NCyc.lt.NExpCy) call Cy8(IOut,IPrint,MxBnd,MxAtCy,MxAtG,MxTrm,
     $  NDih,NBond,IBond,IAtomD,MxCyc,NCyc,NExpCy,NAtC,ICAt,
     $  IAtCyc)
      If(IPrint.gt.0) write(IOut,'(/,'' The Molecule has'',I3,
     $  '' Cycles'')')NCyc
C Set canonical atom numbering.
C Merlino convention: cyclic order starts from the Prelog/CIP priority
C implied by atom type and local connectivity; the input atom index is
C only the final tie-break. Ring-puckering GICs generated later by
C CyGND use this order.
C Do not call SymCyc here: a symmetry-based shift changes the phase
C origin and would make the Fortran and Python QPck/PhiP conventions
C diverge for otherwise identical rings.
      do 10 icyc = 1, NCyc
       if(IPrint.gt.0) write(IOut,'(/,I2,''-Membered Cycle'')')
     $   NAtC(ICyc)
       if(IPrint.gt.0) write(IOut,'('' Input Numbering    '',10I5)')
     $   (ICAt(i,ICyc),i=1,NAtC(ICyc))
       ReNumb=.true.
       call CanCyc(IOut,IPrint,ReNumb,MxBnd,MxAtCy,NAtC(ICyc),
     $   ICAt(1,ICyc),NBond,IBond,EAN)
       if(ReNumb.and.IPrint.gt.0) write(IOut,'('' Canonical'',
     $   '' Numbering'',10I5)') (ICAt(i,ICyc),i=1,NAtC(ICyc))
   10 continue
C find atoms common to 2 or 3 cycles
      if(NCyc.ge.3) call At3Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      if(NCyc.ge.2) call At2Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      if(NCyc.ge.1) call At1Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      if(N1Cyc.gt.0) write(Iout,'(I3,'' Atoms in Only 1 Cycle '')')N1Cyc
      if(N2Cyc.gt.0) write(IOut,'(I3,'' Atoms Linking 2 Cycles'')')N2Cyc
      if(N3Cyc.gt.0) write(IOut,'(I3,'' Atoms Linking 3 Cycles'')')N3Cyc
C find bonds, angles and dihedrals common to 2 Cycles
      if(NCyc.ge.2) then
       Do 20 ICyc=1,NCyc-1
        Do 30 JCyc=ICyc+1,NCyc
         Call BrCyc(IOut,IPrint,MxAtCy,NCyc,NAtC,ICAt,ICyc,JCyc,NAtIdC)
   30   Continue
   20  Continue      
      EndIf
      if(IPrint.gt.0) write(IOut,'('' '')')
      return 
      end
*Deck Cy5 
      Subroutine Cy5(IOut,IPrint,MxBnd,MaxAtD,MxTrmD,NBond,IBond,NDih,
     $  IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,NAtC,ICAt,IAtCyc)
      Implicit None
C
C     Dimensions
      Integer MxBnd,MaxAtD,MxTrmD,MxAtCy,MxCyc
C     Input
      Integer IOut,IPrint,NDih,NExpCy,NCyc
      Integer IBond(MxBnd,*), NBond(*)
      Integer IAtomD(MaxAtD,MxTrmD,*)
      Logical EqCyc 
C     Output
      Integer IAtCyc(*), ICAt(MxAtCy,*), NAtC(*)
C     Local
      Integer I,II,LL,IAt,JAt,KAt,LAt,MAt,NAt,N5
      Integer IC5(5)
C Find five term cycles using the informations on the dihedral angles
      Do 10 i = 1, NDih
       IAt = IAtomD(1,1,i)
       JAt = IAtomD(2,1,i)
       KAt = IAtomD(3,1,i)
       LAt = IAtomD(4,1,i)
       if(LAt.eq.IAt) go to 10 
       Do 20 ii = 1, NBond(IAt)
        MAt = IBond(ii,IAt)
        If(MAt.ne.JAt.and.MAt.ne.KAt.and.MAt.ne.LAt) then
         Do 30 ll = 1, NBond(LAt)
          NAt = IBond(ll,LAt)
          If(MAt.eq.NAt) then
           N5 = 5
           IC5(1) = IAt
           IC5(2) = JAt
           IC5(3) = KAt
           IC5(4) = LAt
           IC5(5) = MAt
           If(.not.EqCyc(MxAtCy,NCyc,NatC,ICAt,N5,IC5)) then
            If(NCyc.ge.NExpCy) Return
            If(NCyc.ge.MxCyc) then
             write(IOut,'('' WARNING: maximum number of cycles ('',I3,
     $       '') reached in Cy5; remaining cycles ignored'')') MxCyc
             Return
            EndIf
            NCyc=NCyc+1
            NAtC(NCyc) = N5
            Call IMove(5,IC5,ICAt(1,NCyc))
            IAtCyc(IAt)  = NCyc
            IAtCyc(JAt)  = NCyc
            IAtCyc(KAt)  = NCyc
            IAtCyc(LAt)  = NCyc
            IAtCyc(MAt)  = NCyc
            If(NCyc.eq.NExpCy) Return
           endif
          endif
   30    Continue
        endIf
   20  Continue
   10 Continue
      Return
      End
*Deck Cy6 
      Subroutine Cy6(IOut,IPrint,MxBnd,MaxAtD,MxTrmD,NBond,IBond,NDih,
     $  IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,NAtC,ICAt,IAtCyc)
      Implicit None
C
C     Dimensions
      Integer MxBnd,MaxAtD,MxTrmD,MxAtCy,MxCyc
C     Input
      Integer IOut,IPrint,NDih,NExpCy,NCyc
      Integer IBond(MxBnd,*),NBond(*)
      Integer IAtomD(MaxAtD,MxTrmD,*)
      Logical EqCyc 
C     Output
      Integer IAtCyc(*),ICAt(MxAtCy,*),NAtC(*)
C     Local
      Integer i,j,n,nn,k,jj,l,kk,m,ll,mm,N6
      Integer IC6(6)
C Find Six term cycles using the informations on the dihedral angles
      Do 300 i = 1, NDih
       j  = IAtomD(1,1,i)
       n  = IAtomD(2,1,i)
       nn = IAtomD(3,1,i)
       k  = IAtomD(4,1,i)
       if(k.eq.j) go to 300
       Do 310 jj = 1, NBond(j)
        L = IBond(jj,j)
        If(l.ne.n.and.l.ne.nn.and.l.ne.k) then
         Do 320 kk = 1, NBond(K)
          M = IBond(kk,K)
          If(L.ne.M.and.M.ne.N.and.M.ne.NN.and.M.ne.j) then
           Do 330 ll = 1,NBond(L)
            mm = IBond(ll,l)
            If(MM.eq.M) then
             N6 = 6
             IC6(1) = j
             IC6(2) = n
             IC6(3) = nn
             IC6(4) = k
             IC6(5) = l
             IC6(6) = m  
             If(.not.EqCyc(MxAtCy,NCyc,NatC,ICAt,N6,IC6)) then
              If(NCyc.ge.NExpCy) return
              If(NCyc.ge.MxCyc) then
               write(IOut,'('' WARNING: maximum number of cycles ('',
     $         I3,'') reached in Cy6; remaining cycles ignored'')')
     $         MxCyc
               return
              EndIf
              NCyc=NCyc+1
              NAtC(NCyc)=N6
              Call IMove(6,IC6,ICAt(1,NCyc))
              IAtCyc(j)  = NCyc
              IAtCyc(k)  = NCyc
              IAtCyc(l)  = NCyc
              IAtCyc(n)  = NCyc
              IAtCyc(nn) = NCyc
              IAtCyc(m)  = NCyc
              If(NCyc.eq.NExpCy) return
             EndIf
            EndIf
  330      Continue
          EndIf
  320    Continue
        EndIf
  310  Continue
  300 Continue
      Return
      End
*Deck Cy7        
      Subroutine Cy7(IOut,IPrint,MaxAtA,MaxAtD,MxTrmA,MxTrmD,
     $  NAng,NDih,NBond,IAtomA,IAtomD,MxAtCy,MxCyc,NCyc,NExpCy,
     $  NAtC,ICAt,IAtCyc)
      Implicit None 
C               
C     Dimensions
      Integer MaxAtA, MaxAtD, MxTrmA, MxTrmD, MxAtCy,MxCyc
C     Input
      Integer IOut, IPrint, NAng, NDih, NCyc,NExpCy
      Integer NBond(*),IAtomA(MAxAtA,MxTrmA,*),IAtomD(MaxAtD,MxTrmD,*)
      Logical EqCyc
C     Output
      Integer IAtCyc(*), ICAt(MxAtCy,*), NAtC(*)
C     Local 
      Integer i, i1, i2, i3, i4, j, j1, j2, j3, j4, k, l, inew
      Integer i1at, i7at, i17at, ic7(7),N7
C
      do 10 i=1,NAng
       if(NCyc.ge.NExpCy) return
       i1at=IAtomA(3,1,i)
       i7at=IAtomA(1,1,i)
       i17at=IAtomA(2,1,i)
       if(NBond(i1at).eq.1.or.NBond(i7at).eq.1) go to 10
c      if(IPrint.eq.3) write(IOut,'('' Angle ('',I5,'')'',3I5)')
c    $   i,i7at,i17at,i1at
       do 20 j=1,NDih
        i1=IAtomD(1,1,j)
        i2=IAtomD(2,1,j) 
        i3=IAtomD(3,1,j)
        i4=IAtomD(4,1,j)
        if(NBond(I1).eq.1.or.NBond(I4).eq.1) go to 20
        if(i1.ne.i1at) go to 20
c        if(IPrint.eq.3) write(IOut,'('' First dihedral ('',I5,'')'',
c    $    4I5)')  j,i1,i2,i3,i4
        do 30 k=1,NDih
         j1=IAtomD(1,1,k)
         j2=IAtomD(2,1,k)
         j3=IAtomD(3,1,k) 
         j4=IAtomD(4,1,k)
         if(j.eq.k) go to 30
         if(NBond(j1).eq.1.or.NBond(j4).eq.1) go to 30 
         if(j1.ne.i4) go to 30
         if(j4.ne.i17at) go to 30 
         if(j2.eq.i1.or.j2.eq.i2.or.j2.eq.i3) go to 30 
         if(j3.eq.i1.or.j3.eq.i2.or.j3.eq.i3) go to 30
         if(j4.eq.i1.or.j4.eq.i2.or.j4.eq.i3) go to 30
c        if(IPrint.eq.3)  write(IOut,'('' Second dihedral ('',I5,'')'',
c    $     4I5)') k,j1,j2,j3,j4
         IC7(1)=i1 
         IC7(2)=i2 
         IC7(3)=i3  
         IC7(4)=i4 
         IC7(5)=j2 
         IC7(6)=j3  
         IC7(7)=j4  
         N7=7
         If(EqCyc(MxAtCy,NCyc,NAtC,ICAt,N7,IC7)) go to 30
         If(NCyc.ge.NExpCy) return
         If(NCyc.ge.MxCyc) then
          write(IOut,'('' WARNING: maximum number of cycles ('',I3,
     $    '') reached in Cy7; remaining cycles ignored'')') MxCyc
          return
         EndIf
         NCyc=NCyc+1
         NAtC(NCyc)=7
         Call IMove(7,IC7,ICAt(1,NCyc))
         do 40 inew=1,7
          IatCyc(IC7(inew))=IAtCyc(IC7(inew))+1
   40    continue 
         If(IPrint.gt.1) then
          write(IOut,'(/,'' In Cy7'')')
          write(IOut,'('' 7-membered Ring'',7I5)') (IC7(l),l=1,7)
         EndIf
   30   continue
   20  continue 
   10 continue
      return
      end
*Deck Cy8        
      Subroutine Cy8(IOut,IPrint,MxBnd,MxAtCy,MxAtD,MxTrmD,NDih,NBond,
     $  IBond,IAtomD,MxCyc,NCyc,NExpCy,NAtC,ICAt,IAtCyc)
      Implicit None 
C               
C     Dimensions
      Integer MxBnd,MxAtCy,MxAtD,MxTrmD,MxCyc
C     Input
      Integer IOut,IPrint,NDih,NCyc,NExpCy
      Integer NBond(*),IBond(MxBnd,*),IAtomD(MxAtD,MxTrmD,*)
      Logical EqCyc
C     Output
      Integer IAtCyc(*),ICAt(MxAtCy,*),NAtC(*)
C     Local 
      Integer i,i1,i2,i3,i4,j,j1,j2,j3,j4,k,l,kk,ll,k1,l1,inew
      Integer i1at, i8at, i18at, i28at, ic8(8)
      Integer N8
      Logical Found1,Found2,Revers
C
      if(NDih.le.2) return
      do 20 j=2,NDih
       if(NCyc.ge.NExpCy) return
       i1=IAtomD(1,1,j)
       i2=IAtomD(2,1,j) 
       i3=IAtomD(3,1,j)
       i4=IAtomD(4,1,j)
       if(NBond(I1).eq.1.or.NBond(I4).eq.1) go to 20
       do 30 k=1,j-1 
        j1=IAtomD(1,1,k)
        j2=IAtomD(2,1,k)
        j3=IAtomD(3,1,k) 
        j4=IAtomD(4,1,k)
        if(NBond(j1).eq.1.or.NBond(J4).eq.1) go to 30
        if(i1.eq.j1.or.i1.eq.j2.or.i1.eq.j3.or.i1.eq.j4) go to 30
        if(i2.eq.j1.or.i2.eq.j2.or.i2.eq.j3.or.i2.eq.j4) go to 30
        if(i3.eq.j1.or.i3.eq.j2.or.i3.eq.j3.or.i3.eq.j4) go to 30
        if(i4.eq.j1.or.i4.eq.j2.or.i4.eq.j3.or.i4.eq.j4) go to 30
        revers=.false.
        found1=.false.
        do 40 kk=1,NBond(j1)
         k1=IBond(kk,j1)
         if(found1) go to 40
         if(k1.eq.i1) then
          found1=.true.
         elseif(k1.eq.i4) then
          found1=.true.
          revers=.true.
         endif
   40   continue
        if(.not.Found1) go to 30
        found2=.false.
        do 45 ll=1,NBond(j4)
         l1=IBond(ll,j4)
         if(found2) go to 45
         if(l1.eq.i1.and.revers) then
          found2=.true.
         elseif(l1.eq.i4.and..not.revers) then
          found2=.true.
         endif
   45   continue 
        if(.not.Found2) go to 30
        IC8(1)=i1 
        IC8(2)=i2 
        IC8(3)=i3  
        IC8(4)=i4 
        if(revers) then
         IC8(5)=j4
         IC8(6)=j3
         IC8(7)=j2
         IC8(8)=j1
        else
         IC8(5)=j1
         IC8(6)=j2 
         IC8(7)=j3  
         IC8(8)=j4
        endif  
        N8=8
        If(EqCyc(MxAtCy,NCyc,NAtC,ICAt,N8,IC8)) go to 30
        If(NCyc.ge.NExpCy) return
        If(NCyc.ge.MxCyc) then
         write(IOut,'('' WARNING: maximum number of cycles ('',I3,
     $   '') reached in Cy8; remaining cycles ignored'')') MxCyc
         return
        EndIf
        NCyc=NCyc+1
        NAtC(NCyc)=8
        Call IMove(8,IC8,ICAt(1,NCyc))
        do 50 inew=1,8
         IatCyc(IC8(inew))=IAtCyc(IC8(inew))+1
   50   continue 
        If(IPrint.gt.1) then
         write(IOut,'(/,'' In Cy8'')')
         write(IOut,'('' 8-membered Ring'',8I5)') (IC8(l),l=1,8)
        EndIf
   30  continue 
   20 continue
      return
      end
*Deck BrCyc
      Subroutine BrCyc(IOut,IPrint,MxAtCy,NCyc,NAtC,ICAt,ICyc,JCyc,
     *  NAtIdC)
      Implicit None
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
      Integer NAtC(*),ICAt(MxAtCy,*)
      Integer MxAtCy,NCyc,ICyc,JCyc,I,J,IAt,JAt,NCI,NCJ,IOut,IPrint
      Integer NBrL,NBrA,NBrD,IBrL,IBrA,IBrD,NAtIdC,ii
      Integer IBic(10)
      If(NCyc.lt.2) then
       write(IOut,'('' Only'',I2,'' Cycles'')') NCyc
       STOP
      ElseIf(ICyc.gt.NCyc) then
       write(IOut,'('' ICyc='',I3,'' But Only'',I3,
     $   '' Cycles Are Present'')') ICyc,NCyc
       STOP
      ElseIf(JCyc.gt.NCyc) then 
       write(IOut,'('' JCyc='',I3,'' But Only'',I3,
     $   '' Cycles Are Present'')') JCyc,NCyc
       STOP
      EndIf
      NCI=NAtC(ICyc)
      NCJ=NAtC(JCyc)
      NAtIdC=0
      Call IClear(10,IBic)
      Do 10 I=1,NCI
       IAt=ICAt(i,ICyc)
       Do 20 J=1,NCJ 
        JAt=ICAt(j,JCyc)
        if(IAt.eq.JAt) then
         NAtIdC=NAtIdC+1
         IBic(NAtIdC)=IAt
        EndIf
   20  Continue
   10 Continue
      If(NAtIdc.eq.2) then
       NBrL=NBrL+1
       IBrL(1,NBrL)=IBic(1)
       IBrL(2,NBrL)=IBic(2)
       IBrL(3,NBrL)=ICyc
       IBrL(4,NBrL)=JCyc
       if(IPrint.gt.0) write(IOut,'(''  BrCyc bond :'',2I5,
     $   '' Between Cycles'',2I3)')(IBrL(ii,NBrL),ii=1,4)
      ElseIf(NAtIdc.eq.3) then
       NBrA=NBrA+1
       IBrA(1,NBrA)=IBic(1)
       IBrA(2,NBrA)=IBic(2)
       IBrA(3,NBrA)=IBic(3)
       IBrA(4,NBrA)=ICyc
       IBrA(5,NBrA)=JCyc
       write(IOut,'(''  BrCyc angle:'',3I5,'' Between Cycles'',2I3)')
     $   (IBrA(ii,NBrA),ii=1,5)
      ElseIf(NAtIDc.eq.4) then
       NBrD=NBrD+1
       IBrD(1,NBrD)=IBic(1)
       IbrD(2,NBrD)=IBic(2)
       IBrD(3,NBrD)=IBic(3)
       IBrD(4,NBrD)=IBic(4)
       IBrD(5,NBrD)=ICyc
       IBrD(6,NBrD)=JCyc
       write(IOut,'(''  BrCyc dihed:'',4I5,'' Between Cycles'',2I3)')
     $   (IBrD(ii,NBrD),ii=1,6)
      EndIf 
      Return
      End
*Deck EqCyc
      Logical Function EqCyc(MxAtCy,NCyc,NAtC,ICAt,NX,ICX)
      Implicit None
      Integer NAtC(*),ICAt(MxAtCy,*),ICX(*)
      Integer MxAtCy,NCyc,NX,ICyc,IRef,JRef,NId
      EqCyc = .False.
      If(NCyc.lt.1) return
      Do 10 ICyc = 1,NCyc
       If(NAtC(ICyc).eq.NX) then
        NId = 0
        Do 20 IRef = 1,NX
         Do 30 JRef = 1,NX
          If(ICAt(JRef,ICyc).eq.ICX(IRef)) NId = NId + 1
   30    Continue
   20   Continue
        If(NId.eq.NX) then
         EqCyc = .True.
         Return
        EndIf
       EndIf
  10  Continue
      Return
      End
*Deck CanCyc
      Subroutine CanCyc(IOut,IPrint,ReNumb,MxBnd,MxAtCy,NAtC,ICAt,
     $  NBond,IBond,EAN)
      Implicit None
      Integer IOut,IPrint,MxBnd,MxAtCy,NAtC
      Integer ICAt(*),NBond(*),IBond(MxBnd,*)
      Real*8 EAN(*)
      Logical ReNumb
C Local
      Integer ICAtOK(MxAtCy),IBest(MxAtCy),ICand(MxAtCy)
      Integer I,II,KK,J,KAt,JAt,Kkii,ICK,IC3K
      Integer IAt,I1,I2,IM,Ini,IEnd,IStart,IDir,IP,Idx
      Integer NBest,IBetter
      Logical BetterCyc
      Integer IrMin1
C
C First rebuild a connected cyclic traversal from the detected cycle.
      call IMove(NAtC,ICAt,ICAtOK)
      IAt=IrMin1(ICAtOK,NAtC,.True.,IM)
      I1=ICAtOK(1)
      I2=ICAtOK(IM)
      ICAtOK(1)=I2
      ICAtOK(IM)=I1
      Ini=1
      IEnd=NAtC-1
      do 40 II=Ini,IEnd
       IAt=ICAtOK(II)
       Kkii=II+1
       do 50 KK=Ini+1,IEnd+1
        KAt=ICAtOK(KK)
        do 60 J=1,NBond(KAt)
         JAt=IBond(J,KAt)
         If(JAt.eq.IAt) Kkii=KK
   60   continue
   50  continue
       ICK=ICAtOK(Kkii)
       IC3K=ICAtOK(II+1)
       ICAtOK(II+1)=ICK
       ICAtOK(Kkii)=IC3K
   40 continue
C
C Then select the Prelog-first rotation and direction.  Atomic number
C is primary, local degree is secondary, and input atom index is only
C the final deterministic tie-break.  This mirrors the Python generator.
      NBest=0
      do 100 IStart=1,NAtC
       do 110 IDir=1,2
        do 120 IP=1,NAtC
         if(IDir.eq.1) then
          Idx=IStart+IP-1
          if(Idx.gt.NAtC) Idx=Idx-NAtC
         else
          Idx=IStart-IP+1
          if(Idx.lt.1) Idx=Idx+NAtC
         endif
         ICand(IP)=ICAtOK(Idx)
  120   continue
        if(NBest.eq.0) then
         call IMove(NAtC,ICand,IBest)
         NBest=1
        else
         if(BetterCyc(NAtC,ICand,IBest,NBond,IBond,MxBnd,EAN))
     $    then
          call IMove(NAtC,ICand,IBest)
         endif
        endif
  110  continue
  100 continue
      if(ReNumb) then
       call IMove(NAtC,IBest,ICAt)
      else
       write(IOut,'('' CanCyc: Prelog Numbering'',10I5)')
     $   (IBest(I),I=1,NAtC)
      endif
      return
      end
*Deck BetterCyc
      Logical Function BetterCyc(NAtC,ICand,IBest,NBond,IBond,
     $ MxBnd,EAN)
      Implicit None
      Integer NAtC,ICand(*),IBest(*),NBond(*),MxBnd
      Integer IBond(MxBnd,*)
      Real*8 EAN(*)
      Integer I,IA,IB,K,KMax
      Real*8 Diff,ZA,ZB,ExoZ
      BetterCyc=.False.
      do 10 I=1,NAtC
       IA=ICand(I)
       IB=IBest(I)
       Diff=EAN(IA)-EAN(IB)
       if(Diff.gt.1.0D-8) then
        BetterCyc=.True.
        return
       elseif(Diff.lt.-1.0D-8) then
        return
       endif
       if(NBond(IA).gt.NBond(IB)) then
        BetterCyc=.True.
        return
       elseif(NBond(IA).lt.NBond(IB)) then
        return
       endif
       KMax=NBond(IA)
       if(NBond(IB).gt.KMax) KMax=NBond(IB)
       do 20 K=1,KMax
        ZA=ExoZ(IA,K,NAtC,ICand,NBond,IBond,MxBnd,EAN)
        ZB=ExoZ(IB,K,NAtC,IBest,NBond,IBond,MxBnd,EAN)
        Diff=ZA-ZB
        if(Diff.gt.1.0D-8) then
         BetterCyc=.True.
         return
        elseif(Diff.lt.-1.0D-8) then
         return
        endif
   20  continue
       if(IA.lt.IB) then
        BetterCyc=.True.
        return
       elseif(IA.gt.IB) then
        return
       endif
   10 continue
      return
      end
*Deck ExoZ
      Real*8 Function ExoZ(IAt,KOrd,NAtC,IRing,NBond,IBond,
     $ MxBnd,EAN)
      Implicit None
      Integer IAt,KOrd,NAtC,IRing(*),NBond(*),MxBnd
      Integer IBond(MxBnd,*)
      Real*8 EAN(*),ZL(20),ZT
      Integer I,J,N,KAt
      Logical InRing
      ExoZ=-1.0D0
      N=0
      do 10 I=1,NBond(IAt)
       KAt=IBond(I,IAt)
       if(.not.InRing(KAt,NAtC,IRing)) then
        if(N.lt.20) then
         N=N+1
         ZL(N)=EAN(KAt)
        endif
       endif
   10 continue
      if(N.lt.KOrd) return
      do 20 I=1,N-1
       do 30 J=I+1,N
        if(ZL(J).gt.ZL(I)) then
         ZT=ZL(I)
         ZL(I)=ZL(J)
         ZL(J)=ZT
        endif
   30  continue
   20 continue
      ExoZ=ZL(KOrd)
      return
      end
*Deck InRing
      Logical Function InRing(IAt,NAtC,IRing)
      Implicit None
      Integer IAt,NAtC,IRing(*),I
      InRing=.False.
      do 10 I=1,NAtC
       if(IRing(I).eq.IAt) then
        InRing=.True.
        return
       endif
   10 continue
      return
      end
*Deck TstCyc
      Subroutine TstCyc(IAt,JAt,NCyc,MxAtCy,NAtC,ICAt,ISmCyc)
      Logical LICyc,LJCyc
      Dimension NatC(*),ICAt(MxAtCy,*)
      ISmCyc=0 
      do 10 icyc=1,NCyc
       N1=NAtC(icyc)
       if(ISmCyc.ne.0) go to 10
       LICyc=.false.
       LJCyc=.false.
       do 20 i1=1,N1
        KAt=ICAt(i1,icyc)
        If(KAt.eq.IAt) LICyc=.true.
        If(KAt.eq.JAt) LJCyc=.true.
   20  continue
       If(LICyc.and.LJCyc) ISmCyc=ICyc
   10 continue
      return
      end
*Deck CyGNA  
      Subroutine CyGNA(IOut,IPrint,MxAtCy,MaxAtA,MxTrmA,NBond,
     $  NAtC,NAng,ICyc,IAtomC,NTermA,IAtomA,ITVA,CoefA,DoLocSVD,
     $  NAtoms,C)
      Implicit None 
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
C Input
      Integer NAtC(*),IAtomC(MxAtCy,*),NBond(*)
      Integer IOut,IPrint,MxAtCy,MaxAtA,MxTrmA,NAng,ICyc,NAtoms
      Integer N2Cyc,N3Cyc,NBrL,NBrA,NBrD
      Integer IAt2C,IAt3C,IBrL,IBrA,IBrD
      Logical DoLocSVD
      Real*8 C(3,*)
C Input/Output
      Integer NTermA(*),IAtomA(MaxAtA,MxTrmA,*),ITVA(*)
      Real*8 CoefA(MxTrmA,*)
C Local
      Integer NAtPrm,NAng0,NAtCyc,Ii,Ip,I4,IAll3C,jj
      Logical ValAng
C
C builds ring coordinates for valence angles
C for symmetry reasons the first angle is: 341,345,234,456 
C for 4-,5-,6-,and 7-membered rings 
C
      if(N3Cyc.ge.NAtC(ICyc)) then
       iall3c=0
       do 10 ii=1,N3Cyc
        do 20 jj=1,NAtC(ICyc)
         if(NBond(IAtomC(jj,ICyc)).gt.3) go to 20
         if(IAt3c(1,ii).eq.IAtomC(jj,ICyc)) iall3c=iall3c+1
   20   continue
   10  continue
       if(IAll3c.ge.NAtC(ICyc)) return  
      endif
      NAtPrm=3
      ValAng=.true.
      NAng0=NAng
      NAtCyc=NAtC(ICyc)
      If(DoLocSVD) then
       call CyGNSVD(IOut,IPrint,MaxAtA,MxTrmA,MxAtCy,ValAng,NAng,
     $  ICyc,NAtC,IAtomC,NTermA,ITVA,IAtomA,CoefA,NAtoms,C)
      Else
       call CycAng(IOut,Iprint,MaxAtA,MxTrmA,MxAtCy,ValAng,NAng,
     $  ICyc,NAtC,IAtomC,NTermA,ITVA,IAtomA,CoefA)
      EndIf
      if(IPrint.gt.0) then
       write(IOut,'(/,I2,''-membered Cycle'')') NAtCyc
       do 30 ip=1,NAtCyc-3
        write(IOut,'(/,'' Angle'',I5)') NAng0+ip
        do 40 i4=1,NtermA(NAng0+ip)
         write(IOut,'(F8.4,''*A('',2(I3,'',''),I3,'')'')')
     $     CoefA(i4,NAng0+ip),(IAtomA(ii,i4,NAng0+ip),ii=1,NAtPrm)
   40   continue
   30  continue
      endif 
      return
      end
*Deck CyGND   
      Subroutine CyGND(IOut,IPrint,MxAtCy,MaxAtD,MxTrmD,NAtC,
     $  NDih,ICyc,IAtomC,NTermD,IAtomD,ITVD,CoefD,DoLocSVD,
     $  NAtoms,C)
      Implicit none 
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
C Input
      Integer NAtC(*),IAtomC(MxAtCy,*)
      Integer IOut,IPrint,MxAtCy,MaxAtD,MxTrmD,NDih,ICyc,NAtoms
      Integer N2Cyc,N3Cyc,NBrL,NBrA,NBrD,IAt2C,IAt3C,IBrL,IBrA,IBrD
      Logical DoLocSVD
      Real*8 C(3,*)
C Input/Output
      Integer NTermD(*),IAtomD(MaxAtD,MxTrmD,*),ITVD(*)
      Real*8 CoefD(MxTrmD,*)
C Local
      Logical ValAng
      Integer NAtPrm,NDih0,NAtCyc,ii,jj,Ip,I4,iall3c,iall2c
C            
C builds ring coordinates for dihedral angles
C for symmetry reasons the first dihedral is always: n123  
C where n is the last atom in the ring
C if all the atoms of the cycle join 3 cycles there are no free dihedrals
      if(N3Cyc.ge.NAtC(ICyc)) then
       iall3c=0        
       do 10 ii=1,N3Cyc 
        do 20 jj=1,NAtC(ICyc)
         if(IAt3c(1,ii).eq.IAtomC(jj,ICyc)) iall3c=iall3c+1
   20   continue
   10  continue
       if(IAll3c.ge.NAtC(ICyc)) then
        write(IOut,'('' No Free Dihedrals for Cycle'',I3)') ICyc
        return
       endif 
      endif
      NAtPrm=4
      ValAng=.false.
      NDih0=NDih
      If(DoLocSVD) then
       call CyGNSVD(IOut,IPrint,MaxAtD,MxTrmD,MxAtCy,ValAng,NDih,
     $  ICyc,NAtC,IAtomC,NTermD,ITVD,IAtomD,CoefD,NAtoms,C)
      Else
       call CycAng(IOut,Iprint,MaxAtD,MxTrmD,MxAtCy,ValAng,NDih,
     $  ICyc,NAtC,IAtomC,NTermD,ITVD,IAtomD,CoefD)
      EndIf
      if(IPrint.gt.0) then
       write(IOut,'(/,I2,''-membered Cycle'')') NAtCyc
       do 30 ip=1,NAtCyc-3
        write(IOut,'(/,'' Dihedral'',I5)') NDih0+ip
        do 40 i4=1,NTermD(NDih0+ip) 
         write(IOut,'(F8.4,''*D('',3(I3,'',''),I3,'')'')')
     $     CoefD(i4,NDih0+ip),(IAtomD(ii,i4,NDih0+ip),ii=1,NAtPrm)
   40   continue
   30  continue
      endif 
      return
      end
*Deck CyGNSVD
      Subroutine CyGNSVD(IOut,IPrint,MxAt,MxTrm,MxAtCy,ValAng,IGnic,
     $ ICyc,NAtC,IAtomC,NTerm,ITV,IAtomG,Coeff,NAtoms,C)
      Implicit Real*8 (A-H,O-Z)
      Parameter(MxLoc=20,MxCart=3000)
      Logical ValAng
      Integer IOut,IPrint,MxAt,MxTrm,MxAtCy,IGnic,ICyc,NAtoms
      Integer NAtC(*),IAtomC(MxAtCy,*),NTerm(*),ITV(*)
      Integer IAtomG(MxAt,MxTrm,*)
      Dimension Coeff(MxTrm,*),C(3,*)
      Dimension BLoc(MxCart,MxLoc),G(MxLoc,MxLoc),EVal(MxLoc)
      Dimension EVec(MxLoc,MxLoc),B(3,4),DB(3,4,3,4),IB(4)
      Dimension Ref(MxLoc,MxLoc),EKeep(MxLoc,MxLoc)
      Dimension Used(MxLoc)
      Integer Rank,RankUse,Used

      NCyc=NAtC(ICyc)
      If(NCyc.le.3) Return
      If(NCyc.gt.MxLoc) then
       Write(IOut,'('' LOCSVD ring block too large for cycle'',I5)')
     $  ICyc
       Return
      EndIf
      NCart=3*NAtoms
      If(NCart.gt.MxCart) then
       Write(IOut,'('' LOCSVD ring block skipped: too many atoms'',I6)')
     $  NAtoms
       Return
      EndIf
      If(ValAng) then
       IStart=3
       If(NCyc.eq.6) IStart=2
       If(NCyc.eq.7) IStart=4
       If(NCyc.eq.8) IStart=2
      Else
       IStart=NCyc
      EndIf

      Call AClear(MxCart*MxLoc,BLoc)
      Do 30 ITerm=1,NCyc
       I1=ITerm+IStart-1
       If(I1.gt.NCyc) I1=I1-NCyc
       I2=I1+1
       If(I2.gt.NCyc) I2=I2-NCyc
       I3=I2+1
       If(I3.gt.NCyc) I3=I3-NCyc
       I4=I3+1
       If(I4.gt.NCyc) I4=I4-NCyc
       IAt=IAtomC(I1,ICyc)
       JAt=IAtomC(I2,ICyc)
       KAt=IAtomC(I3,ICyc)
       LAt=IAtomC(I4,ICyc)
       Call AClear(12,B)
       Call AClear(144,DB)
       If(ValAng) then
        Call DBBend(1,IAt,JAt,KAt,B,IB,C,DB)
       Else
        Call DBTors(1,IAt,JAt,KAt,LAt,B,IB,C,DB)
       EndIf
       Do 20 IC=1,3
        BLoc(3*(IAt-1)+IC,ITerm)=B(IC,1)
        BLoc(3*(JAt-1)+IC,ITerm)=B(IC,2)
        BLoc(3*(KAt-1)+IC,ITerm)=B(IC,3)
        If(.not.ValAng) BLoc(3*(LAt-1)+IC,ITerm)=B(IC,4)
   20  Continue
   30 Continue

      Call AClear(MxLoc*MxLoc,Ref)
      Pi=4.0D0*DATan(1.0D0)
      DNC=DFloat(NCyc)
      VNorm=DSqrt(2.0D0/DNC)
      VNorm1=DSqrt(1.0D0/DNC)
      Do 45 IVar=1,NCyc-3
       IVar2=IVar/2
       IVar1=IVar
       If(IVar.eq.1) IVar1=IVar+1
       If(IVar.eq.4) IVar1=IVar-1
       If(IVar.eq.5) IVar1=IVar-1
       If(IVar.eq.6) IVar1=IVar-2
       Do 44 ITerm=1,NCyc
        SNum=DFloat(2*IVar1*(ITerm-1))
        Val=Pi*SNum/DNC
        If(IVar.eq.2*IVar2) then
         Ref(ITerm,IVar)=VNorm*DSin(Val)
        ElseIf(IVar.lt.NCyc-3) then
         Ref(ITerm,IVar)=VNorm*DCos(Val)
        Else
         Ref(ITerm,IVar)=VNorm1*DCos(DFloat(ITerm-1)*Pi)
        EndIf
   44  Continue
   45 Continue

      Do 60 I=1,NCyc
       Do 50 J=1,NCyc
        Sum=0.0D0
        Do 40 IC=1,NCart
         Sum=Sum+BLoc(IC,I)*BLoc(IC,J)
   40   Continue
        G(I,J)=Sum
   50  Continue
   60 Continue

      Call LocSVDJacobi(MxLoc,NCyc,G,EVal,EVec,Rank)
      RankUse=Rank
      If(RankUse.gt.NCyc-3) RankUse=NCyc-3
      Call IClear(MxLoc,Used)
      Call AClear(MxLoc*MxLoc,EKeep)
      Do 80 IM=1,RankUse
       Best=-1.0D0
       DotBest=0.0D0
       JBest=1
       Do 72 JM=1,Rank
        If(Used(JM).ne.0) GoTo 72
        Dot=0.0D0
        Do 70 ITerm=1,NCyc
         Dot=Dot+Ref(ITerm,IM)*EVec(ITerm,JM)
   70   Continue
        Score=DAbs(Dot)
        If(Score.gt.Best) then
         Best=Score
         DotBest=Dot
         JBest=JM
        EndIf
   72  Continue
       Used(JBest)=1
       Sgn=1.0D0
       If(DotBest.lt.0.0D0) Sgn=-1.0D0
       Do 74 ITerm=1,NCyc
        EKeep(ITerm,IM)=Sgn*EVec(ITerm,JBest)
   74  Continue
       If(Best.lt.5.0D-1) Write(IOut,'('' WARNING: LOCSVD ring '',
     $ ''mode weakly matches CycAng'',3I5,F10.5)') ICyc,IM,JBest,Best
   80 Continue
      Do 84 IM=1,RankUse
       Do 82 ITerm=1,NCyc
        EVec(ITerm,IM)=EKeep(ITerm,IM)
   82  Continue
   84 Continue
      If(IPrint.gt.0) then
       If(ValAng) then
        Write(IOut,'('' LOCSVD ring angle modes'',2I5)') ICyc,RankUse
       Else
        Write(IOut,'('' LOCSVD ring dihedral modes'',2I5)')
     $   ICyc,RankUse
       EndIf
      EndIf

      Do 100 IM=1,RankUse
       IGnic=IGnic+1
       NTerm(IGnic)=NCyc
       ITV(IGnic)=14
       If(.not.ValAng) ITV(IGnic)=1
       Do 90 ITerm=1,NCyc
        I1=ITerm+IStart-1
        If(I1.gt.NCyc) I1=I1-NCyc
        I2=I1+1
        If(I2.gt.NCyc) I2=I2-NCyc
        I3=I2+1
        If(I3.gt.NCyc) I3=I3-NCyc
        I4=I3+1
        If(I4.gt.NCyc) I4=I4-NCyc
        IAtomG(1,ITerm,IGnic)=IAtomC(I1,ICyc)
        IAtomG(2,ITerm,IGnic)=IAtomC(I2,ICyc)
        IAtomG(3,ITerm,IGnic)=IAtomC(I3,ICyc)
        If(.not.ValAng) IAtomG(4,ITerm,IGnic)=IAtomC(I4,ICyc)
        Coeff(ITerm,IGnic)=EVec(ITerm,IM)
   90  Continue
  100 Continue
      Return
      End
*Deck SymCyc
      Subroutine SymCyc(IOut,IPrint,ReNumb,MxBnd,MxAtCy,ICyc,NAtC,ICAt,
     $  NBond,IBond,EAn)
      Implicit Real*8 (A-H,O-Z)
C IO  
      Dimension NatC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*)
      Dimension EAN(*)
C Local
      Dimension ICAtOK(MxAtCy),ICAtSy(MxAtCy)
      Dimension Neq(MxAtCy)
      Logical Found,ReNumb,Even
C
      Even=.false.
      Found=.false.
      tresh=3.0d-5
      NAtCI=NAtC(ICyc)
      IHalf=NAtCI/2
      if(2*IHalf.eq.NAtCI) Even=.true.
      istold=1
      jstold=2
      istok=1
      jstok=2
      Aver=0.0d0
      call IClear(NAtCI,ICAtOK)
      call IClear(NAtCI,ICAtSy)
      call IClear(NAtCI,NEq)
      Do 10 I1=1,NAtCI
       IAt=ICAt(I1,ICyc)
       EAI=EAn(IAt)
       Aver=Aver+EAI
       If(Found) go to 10
       Do 20 J1=1,NAtCI
        If(Found) go to 20
        JAt=ICAt(J1,ICyc)
        EAJ=EAn(JAt) 
        If(IAt.eq.JAt) go to 20
        If(Abs(EAI-EAJ).gt.tresh) go to 20
        Neq(I1)=NEq(I1)+1
        do 30 k1=1,NBond(IAt)
         if(Found) go to 30
         KAt=IBond(k1,IAt)
         If(KAt.ne.JAt) go to 30
         istold=i1 
         jstold=j1
         Found=.true.
   30   continue
   20  continue
   10 continue
      AVer=AVer/Float(NAtCI)
      If(Found) then
       if(jstold.ne.(istold+1)) then
        isv=istold
        istold=jstold
        jstold=isv 
       endif
       if(IPrint.gt.0) write(IOut,'('' Equal and Bonded Atoms'',I5,
     $   '' and'',I5,'' in positions'',I3,'' and'',I3,'' of cycle'',
     $   I3)')ICAt(IStOld,ICyc),ICAt(JStOld,ICyc),IStOld,JStOld,ICyc 
      Else
       If(Even) then
        IStOK=3
        JStOK=IStOK+1
        IF(JStOK.gt.NAtCI) JStOK=JStOK-NAtCI
       Else
        IstOK=2+NAtCI/2
        JStOK=IstOK+1
        if(JStOK.gt.NAtCI) JStOK=JStOK-NAtCI 
       EndIf
       diff=0.0d0
       do 40 i1=1,NAtCI
        If(NEq(i1).gt.0) go to 40
        IAt=ICAt(I1,ICyc)
        diffI=Abs(EAN(IAt)-Aver)
        if(diffI.gt.diff) then
         diff=diffI
         IStOld=i1
         JStOld=IStOld+1
         If(JStOld.gt.NAtCI) JStOld=JStOld-NAtCI
        endif
   40  continue 
       if(IPrint.gt.0) write(IOut,'('' Different Atom'',I5,'' from '',
     $  ''position'',I3,'' to position'',I3,'' of cycle'',I3)') 
     $  ICAt(IStOld,ICyc),ISTOld,IStOK,ICyc
      EndIf
      ICAtOK(IStOK)=ICAt(IStOld,ICyc)
      ICAtOK(JStOK)=ICAt(JStOld,ICyc)
      do 50 I1=3,NAtCI
       INew=JStOK+I1-2
       IF(INew.gt.NAtCI) INew=INew-NAtCI
       IOld=JStold+I1-2
       IF(IOld.gt.NAtcI) IOld=IOld-NAtCI
       ICAtOK(INew)=ICAt(IOld,ICyc)
   50 continue
      if(Found) then
       IAt=ICAtOK(1)
       JAt=ICAtOK(2)
       If(JAt.lt.IAt) then
        call IMove(NAtCI,ICAtOK,ICAtSy)
        ISv=ICAtOK(1)
        ICAtOK(1)=ICAtOK(2)
        ICAtOK(2)=ISV
        do 60 ii=3,NAtCI
         ICAtOK(ii)=ICAtSy(NAtCI+3-ii)
   60   continue
       endif
      endif
      if(Renumb) then
       call IMove(NAtCI,ICAtOK,ICAt(1,ICyc))
      else 
       write(IOut,'('' NewCyc: Symmetric Numbering'',10I5)')(ICAtOK(i),
     $   i=1,NAtCI)
      endif
      return
      end
*Deck At3Cyc
      Subroutine At3Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      Implicit Real*8 (A-H,O-Z)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Dimension NatC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*)
      N3Cyc=0
C     write(IOut,'(''Entering At3Cyc'')')
      do 10 ICyc=3,NCyc
       NAtCI=NAtC(ICyc)
       do 20 I1=1,NAtCI
        IAt=ICAt(I1,ICyc)
        do 30 JCyc=2,ICyc-1
         NAtCJ=NAtC(JCyc)   
         do 40 J1=1,NAtCJ
          JAt=ICAt(J1,JCyc)
          If(JAt.ne.IAt) go to 40
          do 50 KCyc=1,JCyc-1
           NAtCK=NAtC(KCyc)
           do 60 K1=1,NAtCK
            KAt=ICAt(K1,KCyc)
            If(KAt.ne.JAt) go to 60
            N3Cyc=N3Cyc+1
            iat3c(1,n3cyc)=IAt
            iat3c(2,n3cyc)=KCyc
            iat3c(3,n3cyc)=JCyc
            iat3c(4,n3cyc)=ICyc
   60      continue
   50     continue
   40    continue
   30   continue
   20  continue
   10 continue
      return
      end
*Deck At2Cyc
      Subroutine At2Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      Implicit Real*8 (A-H,O-Z)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Dimension NatC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*)
      Logical Also3
      N2Cyc=0
      If(NCyc.lt.2) return
      do 10 ICyc=2,NCyc
       NAtCI=NAtC(ICyc)
       do 20 I1=1,NAtCI
        IAt=ICAt(I1,ICyc)
        do 30 JCyc=1,ICyc-1
         NAtCJ=NAtC(JCyc)   
         do 40 J1=1,NAtCJ
          JAt=ICAt(J1,JCyc)
          If(JAt.ne.IAt) go to 40
          Also3=.false.
          If(N3Cyc.ne.0) then
           do 50 kk=1,N3Cyc
            If(IAt.eq.IAt3C(1,kk)) Also3=.true.
   50      continue
          EndIf
          If(.not.also3) then
           N2Cyc=N2Cyc+1
           iat2c(1,n2cyc)=IAt
           iat2c(2,n2cyc)=JCyc
           iat2c(3,n2cyc)=ICyc
          EndIf
   40    continue
   30   continue
   20  continue
   10 continue
      return
      end
*Deck At1Cyc
      Subroutine At1Cyc(IOut,IPrint,MxAtCy,MxBnd,NCyc,NAtC,ICAt,
     $  NBond,IBond)
      Implicit Real*8 (A-H,O-Z)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/onec/N1Cyc,IAt1C(2,100)
      Dimension NatC(*),ICAt(MxAtCy,*)
      Dimension NBond(*),IBond(MxBnd,*)
      Logical Also2,Also3
      N1Cyc=0
      do 10 ICyc=1,NCyc
       NAtCI=NAtC(ICyc)
       do 20 I1=1,NAtCI
        also2=.false.
        also3=.false.
        IAt=ICAt(I1,ICyc)
        do 30 I2Cyc=1,N2Cyc
         JAt=IAt2C(1,I2Cyc)
         If(IAt.eq.JAt) Also2=.true.
   30   continue
        do 40 I3Cyc=1,N3Cyc 
         JAt=IAt3C(1,I3Cyc)
         If(IAt.eq.JAt) Also3=.true.
   40   continue 
        if(also2.or.also3) go to 20
        N1Cyc=N1Cyc+1
        iat1c(1,n1cyc)=IAt
        iat1c(2,n1cyc)=ICyc
   20  continue
   10 continue
      return
      end
*Deck CycAng
      Subroutine CycAng(IOut,Iprint,MxAt,MxTrm,MxAtCy,ValAng,IGnic,ICyc,
     $  NAtC,IAtomC,NTerm,ITV,IAtomG,Coeff)
      Implicit none
C Input
      Logical ValAng,Even
      Integer NAtC(*),IAtomC(MxAtCy,*)
      Integer IOut,IPrint,MxAtCy,MxAt,MxTrm,IGnic,ICyc,IGnic0
C Input/Output
      Integer NTerm(*),IAtomG(MxAt,MxTrm,*),ITV(*)
      Real*8 Coeff(MxTrm,*)
C Local
      Real*8 Pi,Val,VNorm,VNorm1,SNum,DNC
      Integer NCyc,IStart,IVar,ITerm,IAng1,IAng2,IAng3,IAng4
      Integer ISin,MMin,MMax,MEven,ii,igg,ip,iv2,ivar1
C valang=.true. for valence angles, =.false. for dihedrals
C ignic=number of already defined GNICs
C icyc=number of examined cycle
C natc(icyc) = number of atoms in cycle icyc
C
      pi=4.0d0*ATan(1.0d0)
      Ncyc=NAtC(ICyc)
      ignic0=ignic
      if(NCyc.eq.3) return
C the first atom is determined by symmetry considerations
      if(ValAng) then
       IStart=3
       if(ncyc.eq.6) IStart=2
       if(ncyc.eq.7) IStart=4
       if(ncyc.eq.8) IStart=2 
      else
       IStart=NCyc
      endif 
      DNC=Float(NCyc)
      VNorm=SQrt(2.0d0/DNC)
      VNorm1=SQrt(1.0d0/DNC) 
      do 10 ivar=1,NCyc-3  
       IGnic=IGnic+1
       iv2=ivar/2
       Even=(ivar.eq.(iv2*2))
       do 20 ITerm=1,Ncyc
        iang1=ITerm+istart-1
        if(iang1.gt.ncyc) iang1=iang1-ncyc
        iang2=iang1+1
        if(iang2.gt.ncyc) iang2=iang2-ncyc
        iang3=iang2+1
        if(iang3.gt.ncyc) iang3=iang3-ncyc
        iang4=iang3+1
        if(iang4.gt.ncyc) iang4=iang4-ncyc 
        IAtomG(1,ITerm,IGnic)=IAtomC(IAng1,ICyc)
        IAtomG(2,ITerm,IGnic)=IAtomC(IAng2,ICyc)
        IAtomG(3,ITerm,IGnic)=IAtomC(IAng3,ICyc)
        NTerm(IGnic)=NCyc
        ITV(IGnic)=14
        If(.not.ValAng) then
         IAtomG(4,ITerm,IGnic)=IAtomC(IAng4,ICyc)
         ITV(IGnic)=1
        endif
        ivar1=ivar
        if(ivar.eq.1) ivar1=ivar+1
        if(ivar.eq.4) ivar1=ivar-1
        if(ivar.eq.5) ivar1=ivar-1
        if(ivar.eq.6) ivar1=ivar-2
        SNum=float(2*(ivar1)*(iterm-1))
        Val=pi*SNum/DNC
        If(Even) then
         Coeff(ITerm,IGnic)=VNorm*sin(Val)
        Elseif(IVar.lt.NCyc-3) then
         Coeff(ITerm,IGnic)=VNorm*cos(Val)
        Else
C Different normalization
         Coeff(ITerm,IGnic)=VNorm1*cos(float(iterm-1)*pi)
        EndIf
   20  continue
   10 continue 
      ignic=ignic0+NCyc-3
c     write(IOut,'(/,'' In CycAng'')')
c     do 30 ip=1,NCyc-3
c      igg=ignic0+ip
c      write(IOut,'('' Angle GNIC'',I3)') IGG
c      do 40 iterm=1,NCyc
c       write(IOut,'('' Term ='',I3)') Iterm
c       if(ValAng) then
c        write(IOut,'('' Angle  '',3I5,'' Coeff.'',f8.5)') 
c    $   (IAtomG(ii,ITerm,igg),ii=1,3),coeff(iterm,igg)
c       else
c        write(IOut,'('' Torsion'',4I5,'' Coeff.'',f8.5)')
c    $   (IAtomG(ii,ITerm,igg),ii=1,4),coeff(iterm,igg)
c       endif
c  40  continue
c  30 continue
      return
      end
