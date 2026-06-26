*Deck SymDih
      Subroutine SymDih(IOut,IPrint,MxBnd,JAt,KAt,IAN,NBond,IBond,EAN,C,
     $  Found)
      Implicit Real*8 (A-H,O-Z)
      Dimension IAN(*),NBond(*),IBond(MxBnd,*)
      Dimension C(3,*),EAN(*)
      Logical Found
      NHJ=0
      NYJ=0
      NHK=0
      NYK=0
      ERef=0.0D0
      If(NBond(JAt).lt.2) go to 20
      Do 10 ii=1,NBond(JAt)
       IAt=IBond(ii,JAt)
       If(IAt.eq.KAt) goto 10
       If(ERef.eq.0.0D0) ERef=EAn(IAt)
       If(IAn(IAt).eq.1) NHJ=NHJ+1
       If(EAN(IAt).eq.ERef) NYJ=NYJ+1
   10 Continue 
      If(NHJ.eq.(NBond(JAt)-1)) then 
       Write(IOut,'(''   W-XH'',I1,'' rotor around Bond:'',I3,'' -'',
     $   I3)') NHJ, KAt, JAt
       Found=.true.        
       Return
      ElseIf(NYJ.eq.(NBond(JAt)-1)) then 
       Write(IOut,'(''   W-XY'',I1,'' rotor around Bond:'',I3,'' -'',
     $   I3)') NYJ, KAt, JAt
       Found=.true.
       Return
      EndIf 
   20 If(NBond(KAt).lt.2) go to 40
      Do 30 ii=1,NBond(KAt)
       IAt=IBond(ii,KAt)
       If(IAt.eq.JAt) goto 30
       If(ERef.eq.0.0D0) ERef=EAn(IAt)
       If(IAn(IAt).eq.1) NHK=NHK+1
       If(EAN(IAt).eq.ERef) NYK=NYK+1
   30 Continue 
      If(NHK.eq.(NBond(KAt)-1)) then 
       Write(IOut,'(''   W-XH'',I1,'' rotor around Bond:'',I3,'' -'',
     $   I3)') NHK, JAt, KAt
       Found=.true.
       Return 
      ElseIf(NYK.eq.(NBond(KAt)-1)) then 
       Write(IOut,'(''   W-XY'',I1,'' rotor around Bond:'',I3,'' -'',
     $   I3)') NYK, JAt, KAt
       Found=.true.
       Return
      EndIf
   40 continue
      Return
      End

