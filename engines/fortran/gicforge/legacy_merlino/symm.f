*Deck EqvAtm
      Subroutine EqvAtm(MAXAT, NAt, SYMBOL, GROUP, GROUP_NORM, FAM, 
     $  N, OK, X, Y, Z, CSIZE, CLASS)
      IMPLICIT NONE

      INTEGER MAXAT, NAT, I, J, IAT
      CHARACTER*2 SYMBOL(MAXAT)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION XP, YP, ZP
      DOUBLE PRECISION TOL, DELTA_MAX

C     ---- INFO EXTRA ----
      CHARACTER*16  GROUP, GROUP_NORM
      LOGICAL OK

C     ---- SIMMETRIA ----
      CHARACTER*8 FAM
      INTEGER N
      INTEGER NOPS
      DOUBLE PRECISION R(3,3,200)

C     ---- EQUIVALENZE ----
      INTEGER NCLASS
      INTEGER CLASS(MAXAT,MAXAT)
      INTEGER CSIZE(MAXAT)

C     ---- CLASSIFICAZIONE DEL GRUPPO ----
      CALL CLASSIFY_GROUP(GROUP_NORM, FAM, N, OK)

      IF (.NOT. OK) THEN
         WRITE(*,*) 'ERROR: UNKNOWN POINT GROUP ', GROUP
         STOP
      END IF

C     ---- COSTRUZIONE OPERAZIONI DI SIMMETRIA ----
      CALL BUILD_GROUP_OPS(GROUP_NORM, FAM, N, R, NOPS)

      WRITE(*,*) 'GROUP  = ', GROUP
      WRITE(*,*) 'NUMBER OF SYMMETRY OPERATIONS = ', NOPS
      WRITE(*,'('' NATOMS ='',I4)')  NAT
      DO IAT=1,NAT
       Write(*,*) SYMBOL(IAT),X(IAT),Y(IAT),Z(IAT)
      ENDDO

C     --- FIND_EQUIVALENT_ATOMS ---
C       - assumes Cartesian coordinates
C       - center of mass at origin
C       - exact point-group operations in R(:,:,K)
C       - strict tolerance = 1.0D-4 (exact symmetry)
C       - loose tolerance  = 1.0D-2 (quasi-symmetry)

      TOL = 1.0D-4
      NCLASS = 0 

      CALL FIND_EQUIVALENT_ATOMS
     &   (MAXAT, NAT, SYMBOL, X, Y, Z,
     &    R, NOPS, TOL,
     &    NCLASS, CLASS, CSIZE, DELTA_MAX)

      WRITE(*,*) 'NUMBER OF EQUIVALENCE CLASSES = ', NCLASS

      DO 500 I = 1, NCLASS
       WRITE(*,'(I4,2X,I4,2X,*(I4))')
     &   I, CSIZE(I), (CLASS(I,J), J=1,CSIZE(I))
  500 CONTINUE

      IF (DELTA_MAX .LT. TOL) THEN
       WRITE(*,*) 'SYMMETRY: STRICT (max deviation = ', DELTA_MAX,' Å)'
      ELSE IF (DELTA_MAX .LT. TOL*1.0D2) THEN
       WRITE(*,*) 'SYMMETRY: QUASI (max deviation = ', DELTA_MAX,' Å)'
      ELSE
       WRITE(*,*) 'SYMMETRY: BROKEN (max deviation = ', DELTA_MAX,' Å)'
      END IF

      RETURN
      END

*Deck READ_XYZ_WITH_INFO
      SUBROUTINE READ_XYZ_WITH_INFO
     &     (NAT, SYMBOL, X, Y, Z, GROUP, GROUP_NORM)

      INTEGER NAT
      CHARACTER*2 SYMBOL(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      CHARACTER*(*) GROUP, GROUP_NORM

      INTEGER I, IOS
      CHARACTER*256 LINE

      GROUP  = ' '
      GROUP_NORM = ' '

      OPEN(10, FILE='xyzin', STATUS='OLD')

      READ(10,*) NAT
      READ(10,'(A)') LINE

      CALL PARSE_SECOND_LINE(LINE, GROUP, GROUP_NORM)

      DO 100 I = 1, NAT
         READ(10,*,IOSTAT=IOS) SYMBOL(I), X(I), Y(I), Z(I)
         IF (IOS .NE. 0) THEN
            WRITE(*,*) 'ERROR READING ATOM ', I
            STOP
         END IF
  100 CONTINUE

      CLOSE(10)

      RETURN
      END
*Deck PARSE_SECOND_LINE
      SUBROUTINE PARSE_SECOND_LINE(LINE, GROUP, GROUP_NORM)

      CHARACTER*(*) LINE, GROUP, GROUP_NORM
      CHARACTER*32 TOK(10), TWORK
      INTEGER NTOK, I
      LOGICAL IS_POINT_GROUP

      GROUP  = ' '
      GROUP_NORM = ' '

      CALL SPLIT_TOKENS(LINE, TOK, NTOK)

      DO 200 I = 1, NTOK

C        copia di lavoro per il gruppo
         TWORK = TOK(I)
         CALL TOLOWER(TWORK)

         IF (IS_POINT_GROUP(TWORK)) THEN
C           salva il gruppo ORIGINALE (non alterato)
            GROUP = TOK(I)
            GROUP_NORM = TWORK
         END IF

  200 CONTINUE

      RETURN
      END
*Deck SPLIT_TOKENS
      SUBROUTINE SPLIT_TOKENS(LINE, TOK, NTOK)

      CHARACTER*(*) LINE
      CHARACTER*(*) TOK(*)
      INTEGER NTOK

      INTEGER I, START, LENL

      LENL = LEN(LINE)
      NTOK = 0
      START = 1

      DO 300 I = 1, LENL+1
         IF (I .GT. LENL .OR. LINE(I:I) .EQ. ' ') THEN
            IF (I .GT. START) THEN
               NTOK = NTOK + 1
               TOK(NTOK) = LINE(START:I-1)
            END IF
            START = I + 1
         END IF
  300 CONTINUE

      RETURN
      END
*Deck IS_POINT_GROUP
      LOGICAL FUNCTION IS_POINT_GROUP(TOK)

      CHARACTER*(*) TOK
      CHARACTER*8 FAM
      INTEGER N
      LOGICAL OK

      CALL CLASSIFY_GROUP(TOK, FAM, N, OK)

      IS_POINT_GROUP = OK
      RETURN
      END
*Deck CLASSIFY_GROUP
      SUBROUTINE CLASSIFY_GROUP(GROUP, FAM, N, OK)

      CHARACTER*(*) GROUP, FAM
      INTEGER N
      LOGICAL OK

      INTEGER L
      LOGICAL GOTN

      OK  = .FALSE.
      FAM = ' '
      N   = 0
      L   = LEN_TRIM(GROUP)

C     ---- casi banali ----
      IF (GROUP .EQ. 'c1') THEN
         FAM = 'C1'
         OK  = .TRUE.
         RETURN
      END IF

      IF (GROUP .EQ. 'ci') THEN
         FAM = 'Ci'
         OK  = .TRUE.
         RETURN
      END IF

      IF (GROUP .EQ. 'cs') THEN
         FAM = 'Cs'
         OK  = .TRUE.
         RETURN
      END IF

C     ---- gruppi cubici ----
      IF (GROUP .EQ. 't'  .OR. GROUP .EQ. 'td' .OR.
     &    GROUP .EQ. 'th' .OR. GROUP .EQ. 'o'  .OR.
     &    GROUP .EQ. 'oh' .OR. GROUP .EQ. 'i'  .OR.
     &    GROUP .EQ. 'ih') THEN
         FAM = GROUP
         OK  = .TRUE.
         RETURN
      END IF

C     ---- gruppi Cn, Cnv, Cnh ----
      IF (GROUP(1:1) .EQ. 'c') THEN
         CALL GET_N_FROM_GROUP(GROUP, N, GOTN)
         IF (GOTN) THEN
            IF (INDEX(GROUP,'v') .GT. 0) THEN
               FAM = 'Cnv'
            ELSE IF (INDEX(GROUP,'h') .GT. 0) THEN
               FAM = 'Cnh'
            ELSE
               FAM = 'Cn'
            END IF
            OK = .TRUE.
            RETURN
         END IF
      END IF

C     ---- gruppi Dn, Dnd, Dnh ----
      IF (GROUP(1:1) .EQ. 'd') THEN
         CALL GET_N_FROM_GROUP(GROUP, N, GOTN)
         IF (GOTN) THEN
            IF (INDEX(GROUP,'d') .GT. 0) THEN
               FAM = 'Dnd'
            ELSE IF (INDEX(GROUP,'h') .GT. 0) THEN
               FAM = 'Dnh'
            ELSE
               FAM = 'Dn'
            END IF
            OK = .TRUE.
            RETURN
         END IF
      END IF

      RETURN
      END
*Deck TOLOWER
      SUBROUTINE TOLOWER(STR)

      CHARACTER*(*) STR
      INTEGER I, C

      DO I = 1, LEN(STR)
         C = ICHAR(STR(I:I))
         IF (C .GE. 65 .AND. C .LE. 90) THEN
            STR(I:I) = CHAR(C + 32)
         END IF
      END DO

      RETURN
      END
*Deck GET_N_FROM_GROUP
      SUBROUTINE GET_N_FROM_GROUP(GROUP, N, OK)

      CHARACTER*(*) GROUP
      INTEGER N
      LOGICAL OK

      CHARACTER*16 NUM
      INTEGER I, J, C

      OK = .FALSE.
      N  = 0
      NUM = ' '

      J = 0
      DO 100 I = 2, LEN_TRIM(GROUP)
         C = ICHAR(GROUP(I:I))
         IF (C .GE. 48 .AND. C .LE. 57) THEN
            J = J + 1
            NUM(J:J) = GROUP(I:I)
         ELSE
            IF (J .GT. 0) GOTO 200
         END IF
  100 CONTINUE

  200 CONTINUE
      IF (J .GT. 0) THEN
         READ(NUM(1:J),*) N
         OK = .TRUE.
      END IF

      RETURN
      END

*Deck BUILD_GROUP_OPS
      SUBROUTINE BUILD_GROUP_OPS(GROUP, FAM, N, R, NOPS)

      CHARACTER*(*) GROUP, FAM
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)

      NOPS = 0

C     ---- C groups ----
      IF (FAM .EQ. 'C') THEN
         CALL OPS_CN(N, R, NOPS)
         RETURN
      END IF

      IF (FAM .EQ. 'Cv') THEN
         CALL OPS_CNV(N, R, NOPS)
         RETURN
      END IF

      IF (FAM .EQ. 'Ch') THEN
         CALL OPS_CNH(N, R, NOPS)
         RETURN
      END IF

C     ---- D groups ----
      IF (FAM .EQ. 'D') THEN
         CALL OPS_DN(N, R, NOPS)
         RETURN
      END IF

      IF (FAM .EQ. 'Dnh') THEN
         CALL OPS_DNH(N, R, NOPS)
         RETURN
      END IF

      IF (FAM .EQ. 'Dnd') THEN
         CALL OPS_DND(N, R, NOPS)
         RETURN
      END IF

C     ---- Cubic groups ----
      IF (GROUP .EQ. 'Td') THEN
         CALL OPS_TD(R, NOPS)
         RETURN
      END IF

      IF (GROUP .EQ. 'Oh') THEN
         CALL OPS_OH(R, NOPS)
         RETURN
      END IF

      WRITE(*,*) 'ERROR: GROUP NOT IMPLEMENTED ', GROUP
      STOP
      END

*Deck ADD_OP
      SUBROUTINE ADD_OP(R, NOPS, A)
      IMPLICIT NONE
      INTEGER NOPS
      DOUBLE PRECISION R(3,3,*), A(3,3)
      INTEGER I, J

      NOPS = NOPS + 1

      DO I = 1, 3
       DO J = 1, 3
         R(I,J,NOPS) = A(I,J)
       END DO
      END DO

      RETURN
      END

*Deck ADD_IDENTITY
      SUBROUTINE ADD_IDENTITY(R, NOPS)

      DOUBLE PRECISION R(3,3,*)
      INTEGER NOPS

      NOPS = NOPS + 1

      R(1,1,NOPS) = 1.0D0
      R(1,2,NOPS) = 0.0D0
      R(1,3,NOPS) = 0.0D0
      R(2,1,NOPS) = 0.0D0
      R(2,2,NOPS) = 1.0D0
      R(2,3,NOPS) = 0.0D0
      R(3,1,NOPS) = 0.0D0
      R(3,2,NOPS) = 0.0D0
      R(3,3,NOPS) = 1.0D0

      RETURN
      END
*Deck ADD_INVERSION
      SUBROUTINE ADD_INVERSION(R, NOPS)

      DOUBLE PRECISION R(3,3,*)
      INTEGER NOPS

      NOPS = NOPS + 1

      R(1,1,NOPS) = -1.0D0
      R(1,2,NOPS) =  0.0D0
      R(1,3,NOPS) =  0.0D0
      R(2,1,NOPS) =  0.0D0
      R(2,2,NOPS) = -1.0D0
      R(2,3,NOPS) =  0.0D0
      R(3,1,NOPS) =  0.0D0
      R(3,2,NOPS) =  0.0D0
      R(3,3,NOPS) = -1.0D0

      RETURN
      END
*Deck ADD_MIRROR_XY
      SUBROUTINE ADD_MIRROR_XY(R, NOPS)

      DOUBLE PRECISION R(3,3,*)
      INTEGER NOPS

      NOPS = NOPS + 1

      R(1,1,NOPS) =  1.0D0
      R(1,2,NOPS) =  0.0D0
      R(1,3,NOPS) =  0.0D0
      R(2,1,NOPS) =  0.0D0
      R(2,2,NOPS) =  1.0D0
      R(2,3,NOPS) =  0.0D0
      R(3,1,NOPS) =  0.0D0
      R(3,2,NOPS) =  0.0D0
      R(3,3,NOPS) = -1.0D0

      RETURN
      END

*Deck ROTZ
      SUBROUTINE ROTZ(ANG, A)
      DOUBLE PRECISION ANG, A(3,3)
      DOUBLE PRECISION C, S

      C = DCOS(ANG)
      S = DSIN(ANG)

      A(1,1)= C
      A(1,2)=-S
      A(1,3)= 0.0D0
      A(2,1)= S
      A(2,2)= C
      A(2,3)= 0.0D0
      A(3,1)= 0.0D0
      A(3,2)= 0.0D0
      A(3,3)= 1.0D0
      RETURN
      END

*Deck ROTX
      SUBROUTINE ROTX(ANG, A)
      DOUBLE PRECISION ANG, A(3,3)
      DOUBLE PRECISION C, S

      C = DCOS(ANG)
      S = DSIN(ANG)

      A(1,1)= 1.0D0
      A(1,2)= 0.0D0
      A(1,3)= 0.0D0
      A(2,1)= 0.0D0
      A(2,2)= C
      A(2,3)=-S
      A(3,1)= 0.0D0
      A(3,2)= S
      A(3,3)= C
      RETURN
      END

*Deck ROTY
      SUBROUTINE ROTY(ANG, A)
      DOUBLE PRECISION ANG, A(3,3)
      DOUBLE PRECISION C, S

      C = DCOS(ANG)
      S = DSIN(ANG)

      A(1,1)=  C
      A(1,2)=  0.0D0
      A(1,3)=  S
      A(2,1)=  0.0D0
      A(2,2)=  1.0D0
      A(2,3)=  0.0D0
      A(3,1)= -S
      A(3,2)=  0.0D0
      A(3,3)=  C
      RETURN
      END

*Deck ROT_AXIS_PI
      SUBROUTINE ROT_AXIS_PI(U, R)
      DOUBLE PRECISION U(3), R(3,3)

C  Rotazione di PI attorno al versore U (|U|=1)

      R(1,1) =  2*U(1)*U(1) - 1
      R(2,2) =  2*U(2)*U(2) - 1
      R(3,3) =  2*U(3)*U(3) - 1

      R(1,2) =  2*U(1)*U(2)
      R(2,1) =  R(1,2)

      R(1,3) =  2*U(1)*U(3)
      R(3,1) =  R(1,3)

      R(2,3) =  2*U(2)*U(3)
      R(3,2) =  R(2,3)

      RETURN
      END

*Deck ROT_AXIS
      SUBROUTINE ROT_AXIS(U, ANG, R)
      DOUBLE PRECISION U(3), ANG, R(3,3)
      DOUBLE PRECISION C, S, OC
      DOUBLE PRECISION UX, UY, UZ

C     Rotazione di ANG attorno al versore U (|U|=1)

      UX = U(1)
      UY = U(2)
      UZ = U(3)

      C  = DCOS(ANG)
      S  = DSIN(ANG)
      OC = 1.0D0 - C

      R(1,1) = C + OC*UX*UX
      R(1,2) = OC*UX*UY - S*UZ
      R(1,3) = OC*UX*UZ + S*UY

      R(2,1) = OC*UY*UX + S*UZ
      R(2,2) = C + OC*UY*UY
      R(2,3) = OC*UY*UZ - S*UX

      R(3,1) = OC*UZ*UX - S*UY
      R(3,2) = OC*UZ*UY + S*UX
      R(3,3) = C + OC*UZ*UZ

      RETURN
      END

*Deck REF_XZ
      SUBROUTINE REF_XZ(A)
      DOUBLE PRECISION A(3,3)
      A(1,1)= 1.0D0
      A(1,2)= 0.0D0
      A(1,3)= 0.0D0
      A(2,1)= 0.0D0
      A(2,2)=-1.0D0
      A(2,3)= 0.0D0
      A(3,1)= 0.0D0
      A(3,2)= 0.0D0
      A(3,3)= 1.0D0
      RETURN
      END

*Deck REF_XY
      SUBROUTINE REF_XY(A)
      DOUBLE PRECISION A(3,3)
      A(1,1)= 1.0D0
      A(1,2)= 0.0D0
      A(1,3)= 0.0D0
      A(2,1)= 0.0D0
      A(2,2)= 1.0D0
      A(2,3)= 0.0D0
      A(3,1)= 0.0D0
      A(3,2)= 0.0D0
      A(3,3)=-1.0D0
      RETURN
      END

*Deck REF_YZ
      SUBROUTINE REF_YZ(A)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3)
      A(1,1) = -1.0D0
      A(1,2) =  0.0D0
      A(1,3) =  0.0D0
      A(2,1) =  0.0D0
      A(2,2) =  1.0D0
      A(2,3) =  0.0D0
      A(3,1) =  0.0D0
      A(3,2) =  0.0D0
      A(3,3) =  1.0D0
      RETURN
      END

*Deck INV 
      SUBROUTINE INV(A)
      DOUBLE PRECISION A(3,3)
      A(1,1)=-1.0D0
      A(1,2)= 0.0D0
      A(1,3)= 0.0D0
      A(2,1)= 0.0D0
      A(2,2)=-1.0D0
      A(2,3)= 0.0D0
      A(3,1)= 0.0D0
      A(3,2)= 0.0D0
      A(3,3)=-1.0D0
      RETURN
      END

*Deck MATMUL3
      SUBROUTINE MATMUL3(A, C, B)
      DOUBLE PRECISION A(3,3), C(3,3), B(3,3)
      INTEGER I, J, K
      DOUBLE PRECISION S

      DO 30 I=1,3
         DO 20 J=1,3
            S = 0.0D0
            DO 10 K=1,3
               S = S + A(I,K)*C(K,J)
   10       CONTINUE
            B(I,J) = S
   20    CONTINUE
   30 CONTINUE
      RETURN
      END

*Deck OPS_CN
      SUBROUTINE OPS_CN(N, R, NOPS)
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION PI, ANG
      INTEGER K

      PI = 4.0D0*DATAN(1.0D0)

      DO 100 K = 0, N-1
         ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
         CALL ROTZ(ANG, A)
         CALL ADD_OP(R, NOPS, A)
  100 CONTINUE

      RETURN
      END

*Deck OPS_CNV
      SUBROUTINE OPS_CNV(N, R, NOPS)
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3), S0(3,3), B(3,3)
      DOUBLE PRECISION PI, ANG
      INTEGER K

      PI = 4.0D0*DATAN(1.0D0)

C     rotazioni Cn^k
      DO 100 K = 0, N-1
         ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
         CALL ROTZ(ANG, A)
         CALL ADD_OP(R, NOPS, A)
  100 CONTINUE

C     riflessioni verticali: C^k * sigma(xz)
      CALL REF_XZ(S0)
      DO 200 K = 0, N-1
         ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
         CALL ROTZ(ANG, A)
         CALL MATMUL3(A, S0, B)
         CALL ADD_OP(R, NOPS, B)
  200 CONTINUE

      RETURN
      END

*Deck OPS_CHN
      SUBROUTINE OPS_CNH(N, R, NOPS)
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3), SH(3,3), B(3,3)
      DOUBLE PRECISION PI, ANG
      INTEGER K

      PI = 4.0D0*DATAN(1.0D0)

C     Cn^k
      DO 100 K = 0, N-1
         ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
         CALL ROTZ(ANG, A)
         CALL ADD_OP(R, NOPS, A)
  100 CONTINUE

C     sigma_h * Cn^k
      CALL REF_XY(SH)
      DO 200 K = 0, N-1
         ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
         CALL ROTZ(ANG, A)
         CALL MATMUL3(SH, A, B)
         CALL ADD_OP(R, NOPS, B)
  200 CONTINUE

      RETURN
      END

*Deck OPS_DN
      SUBROUTINE OPS_DN(N, R, NOPS)
      IMPLICIT NONE
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)

      DOUBLE PRECISION A(3,3), U(3)
      DOUBLE PRECISION PI, ANG
      INTEGER K

      PI = 4.0D0*DATAN(1.0D0)

      NOPS = 0

C === Cn^k ===
      DO K = 0, N-1
       ANG = 2.0D0*PI*DBLE(K)/DBLE(N)
       CALL ROTZ(ANG, A)
       CALL ADD_OP(R, NOPS, A)
      END DO

C === N rotazioni C2 attorno ad assi nel piano XY ===
      DO K = 0, N-1
       ANG = PI*DBLE(K)/DBLE(N)

       U(1) = DCOS(ANG)
       U(2) = DSIN(ANG)
       U(3) = 0.0D0

       CALL ROT_AXIS_PI(U, A)
       CALL ADD_OP(R, NOPS, A)
      END DO

      RETURN
      END

*Deck OPS_DNH
      SUBROUTINE OPS_DNH(N, R, NOPS)
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3), SH(3,3), B(3,3)
      INTEGER I0, I

C     prima costruisci Dn dentro R
      I0 = NOPS
      CALL OPS_DN(N, R, NOPS)

C     poi aggiungi sigma_h * ogni operazione già in lista
      CALL REF_XY(SH)

      DO 100 I = I0+1, NOPS
         CALL MATMUL3(SH, R(1,1,I), B)
         CALL ADD_OP(R, NOPS, B)
  100 CONTINUE

      RETURN
      END

*Deck REF_XEQY
      SUBROUTINE REF_XEQY(A)
      DOUBLE PRECISION A(3,3)
      A(1,1)=0.0D0
      A(1,2)=1.0D0
      A(1,3)=0.0D0
      A(2,1)=1.0D0
      A(2,2)=0.0D0
      A(2,3)=0.0D0
      A(3,1)=0.0D0
      A(3,2)=0.0D0
      A(3,3)=1.0D0
      RETURN
      END

*Deck OPS_DND
      SUBROUTINE OPS_DND(N, R, NOPS)
      INTEGER N, NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION SD(3,3), B(3,3)
      INTEGER I0, I

      I0 = NOPS
      CALL OPS_DN(N, R, NOPS)

      CALL REF_XEQY(SD)

      DO 100 I = I0+1, NOPS
         CALL MATMUL3(SD, R(1,1,I), B)
         CALL ADD_OP(R, NOPS, B)
  100 CONTINUE

      RETURN
      END

*Deck OPS_TD
      SUBROUTINE OPS_TD(R, NOPS)

      INTEGER NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3), B(3,3)
      INTEGER P(3,6)
      INTEGER IP, I
      INTEGER S1, S2, S3
      DOUBLE PRECISION D
      DOUBLE PRECISION DET3

C     permutazioni degli assi
      DATA P / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1 /

      NOPS = 0

C     ---- ROTAZIONI PROPRIE (12) ----
      DO 100 IP = 1, 6

C        solo permutazioni PARI
         IF (IP .EQ. 2 .OR. IP .EQ. 3 .OR. IP .EQ. 6) GOTO 100

         DO 90 S1 = -1, 1, 2
         DO 80 S2 = -1, 1, 2
         DO 70 S3 = -1, 1, 2

            CALL MAT_ZERO(A)

            A(1,P(1,IP)) = DBLE(S1)
            A(2,P(2,IP)) = DBLE(S2)
            A(3,P(3,IP)) = DBLE(S3)

            D = DET3(A)

C           rotazioni proprie: det = +1
            IF (D .GT. 0.5D0) THEN
               CALL ADD_OP(R, NOPS, A)
            END IF

   70    CONTINUE
   80    CONTINUE
   90    CONTINUE

  100 CONTINUE

C     ---- IMPROPRIE (12) = inversione * rotazioni ----
      DO 200 I = 1, 12
         CALL MAT_COPY(R(1,1,I), B)
         CALL MAT_INV(B)
         CALL ADD_OP(R, NOPS, B)
  200 CONTINUE

C     check finale
      IF (NOPS .NE. 24) THEN
         WRITE(*,*) 'WARNING: OPS_TD generated ', NOPS,
     &              ' operations (expected 24)'
      END IF

      RETURN
      END
*Deck MAT_COPY
      SUBROUTINE MAT_COPY(A, B)
      DOUBLE PRECISION A(3,3), B(3,3)
      INTEGER I,J
      DO 20 I=1,3
         DO 10 J=1,3
            B(I,J) = A(I,J)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END

*Deck MAT_INV
      SUBROUTINE MAT_INV(A)
      DOUBLE PRECISION A(3,3)
      INTEGER I,J
      DO 20 I=1,3
         DO 10 J=1,3
            A(I,J) = -A(I,J)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END

*Deck OPS_OH
      SUBROUTINE OPS_OH(R, NOPS)

      INTEGER NOPS
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION A(3,3)
      INTEGER P(3,6)
      INTEGER IP, I
      INTEGER S1, S2, S3
      DOUBLE PRECISION D
      DOUBLE PRECISION DET3

C     permutazioni degli assi
      DATA P / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1 /

      NOPS = 0

      DO 100 IP = 1, 6

         DO 90 S1 = -1, 1, 2
         DO 80 S2 = -1, 1, 2
         DO 70 S3 = -1, 1, 2

C           costruisci matrice con una sola entry per riga
            CALL MAT_ZERO(A)

            A(1, P(1,IP)) = DBLE(S1)
            A(2, P(2,IP)) = DBLE(S2)
            A(3, P(3,IP)) = DBLE(S3)

C           accetta solo se determinante = +/-1 (matrici ortogonali)
            D = DET3(A)
            IF (D .GT.  0.5D0 .OR. D .LT. -0.5D0) THEN
               CALL ADD_OP(R, NOPS, A)
            END IF

   70    CONTINUE
   80    CONTINUE
   90    CONTINUE

  100 CONTINUE

C     check (deve essere 48)
      IF (NOPS .NE. 48) THEN
         WRITE(*,*) 'WARNING: OPS_OH generated ', NOPS,
     &              ' operations (expected 48)'
      END IF

      RETURN
      END

*Deck MAT0
      SUBROUTINE MAT_ZERO(A)
      DOUBLE PRECISION A(3,3)
      INTEGER I, J

      DO 20 I = 1, 3
         DO 10 J = 1, 3
            A(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE

      RETURN
      END

*Deck DET3
      DOUBLE PRECISION FUNCTION DET3(A)
      DOUBLE PRECISION A(3,3)

      DET3 =
     &  A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
     &- A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
     &+ A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      RETURN
      END

*Deck MATCH_ATOM
      SUBROUTINE MATCH_ATOM
     &   (XP, YP, ZP, SYM, NAT,
     &    SYMBOL, X, Y, Z,
     &    TOL_LOOSE, IFOUND, OK, DELTA)
      IMPLICIT NONE

      INTEGER NAT
      CHARACTER*2 SYM, SYMBOL(*)
      DOUBLE PRECISION XP, YP, ZP
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION TOL_LOOSE

      INTEGER IFOUND
      LOGICAL OK
      DOUBLE PRECISION DELTA

      INTEGER I
      DOUBLE PRECISION DX, DY, DZ
      DOUBLE PRECISION D2, D2BEST, D2MAX

      OK     = .FALSE.
      IFOUND = 0
      DELTA  = 1.0D99

      D2MAX  = TOL_LOOSE*TOL_LOOSE
      D2BEST = D2MAX

      DO I = 1, NAT

         IF (SYMBOL(I) .NE. SYM) CYCLE

         DX = XP - X(I)
         DY = YP - Y(I)
         DZ = ZP - Z(I)

         D2 = DX*DX + DY*DY + DZ*DZ

         IF (D2 .LE. D2MAX) THEN
            IF (.NOT. OK .OR. D2 .LT. D2BEST) THEN
               D2BEST = D2
               IFOUND = I
               OK     = .TRUE.
            END IF
         END IF

      END DO

      IF (OK) THEN
         DELTA = DSQRT(D2BEST)
      END IF

      RETURN
      END

*Deck FIND_EQUIVALENT_ATOMS
      SUBROUTINE FIND_EQUIVALENT_ATOMS
     &   (MAXAT, NAT, SYMBOL, X, Y, Z,
     &    R, NOPS, TOL,
     &    NCLASS, CLASS, CSIZE, DELTA_MAX)
      IMPLICIT NONE

C     ---- INPUT ----
      INTEGER MAXAT, NAT, NOPS
      CHARACTER*2 SYMBOL(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION R(3,3,*)
      DOUBLE PRECISION TOL

C     ---- OUTPUT ----
      INTEGER NCLASS
      INTEGER CLASS(MAXAT,MAXAT)
      INTEGER CSIZE(MAXAT)
      DOUBLE PRECISION DELTA_MAX

C     ---- LOCALI ----
      LOGICAL USED(NAT)
      INTEGER I, J, K
      INTEGER IA, IB
      INTEGER OLD_SIZE
      DOUBLE PRECISION XP, YP, ZP
      CHARACTER*2 SYM
      LOGICAL OK
      DOUBLE PRECISION TOL_LOOSE, DELTA

C     Loose Tolerance
      TOL_LOOSE = TOL*1.0D2   
      DELTA_MAX = 0.0D0
C     ---- INIZIALIZZAZIONE ----
      DO I = 1, NAT
         USED(I) = .FALSE.
      END DO

      NCLASS = 0

C     ---- LOOP PRINCIPALE SU ATOMI ----
      DO I = 1, NAT

         IF (USED(I)) CYCLE

C        nuova classe
         NCLASS = NCLASS + 1
         CSIZE(NCLASS) = 1
         CLASS(NCLASS,1) = I
         USED(I) = .TRUE.

C        chiusura transitiva dell’orbita
  200    CONTINUE
         OLD_SIZE = CSIZE(NCLASS)

         DO J = 1, OLD_SIZE
            IA  = CLASS(NCLASS,J)
            SYM = SYMBOL(IA)

            DO K = 1, NOPS

C              applica operazione di simmetria
               XP = R(1,1,K)*X(IA) + R(1,2,K)*Y(IA)
     &              + R(1,3,K)*Z(IA)
               YP = R(2,1,K)*X(IA) + R(2,2,K)*Y(IA)
     &              + R(2,3,K)*Z(IA)
               ZP = R(3,1,K)*X(IA) + R(3,2,K)*Y(IA)
     &              + R(3,3,K)*Z(IA)

C              trova atomo equivalente
               CALL MATCH_ATOM(XP,YP,ZP,SYM,NAT,
     &                          SYMBOL,X,Y,Z,
     &                          TOL_LOOSE,IB,OK,DELTA)

C              controllo ROBUSTO: OK + IB valido
               IF (OK .AND. IB.GE.1 .AND. IB.LE.NAT) THEN
                  DELTA_MAX = MAX(DELTA_MAX, DELTA)
                  IF (.NOT. USED(IB)) THEN
                     CSIZE(NCLASS) = CSIZE(NCLASS) + 1
                     CLASS(NCLASS,CSIZE(NCLASS)) = IB
                     USED(IB) = .TRUE.
                  END IF
               END IF

            END DO
         END DO

C        ripeti se la classe è cresciuta
         IF (CSIZE(NCLASS) .GT. OLD_SIZE) GOTO 200

      END DO

C  DELTA_MAX = maximum atomic deviation under all symmetry operations
C              (measure of quasi-symmetry quality)

      RETURN
      END

*Deck DETERMINE_POINT_GROUP
      SUBROUTINE DETERMINE_POINT_GROUP
     & (MAXAT, NAT, SYMBOL, X, Y, Z,
     &  TOL_STRICT,
     &  GROUP_NAME, SYM_STATUS, DELTA_MAX)
      IMPLICIT NONE

C=================================================================
C  High-level determination of molecular point group (Fortran77)
C  Orchestrates symmetry detection, analysis and classification.
C  Also enumerates and prints all subgroups of the final group.
C=================================================================

C ---- INPUT ----
      INTEGER MAXAT, NAT
      CHARACTER*2 SYMBOL(MAXAT)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION TOL_STRICT

C ---- OUTPUT ----
      CHARACTER*8 GROUP_NAME
      INTEGER SYM_STATUS
      DOUBLE PRECISION DELTA_MAX

C ---- LOCALS ----
      INTEGER NOPS, STATUS
      DOUBLE PRECISION R_GROUP(3,3,200)

      INTEGER MAX_CN, NC2, NC3, NC4, NC5
      LOGICAL HAS_I, HAS_SIGMA, HAS_SN
      INTEGER FAMILY

      LOGICAL ISLIN
      DOUBLE PRECISION RINV(3,3)
      DOUBLE PRECISION TOL_LOOSE, DELINV
      LOGICAL OKINV

C ---- Subgroup data ----
      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)
      INTEGER NFOUND
      INTEGER NSUB_ALL(MAXSUB)
      DOUBLE PRECISION R_SUB_ALL(3,3,MAXSUB,MAXOPS)

      INTEGER NUNIQ
      INTEGER NSUB_UNIQ(MAXSUB)
      DOUBLE PRECISION R_SUB_UNIQ(3,3,MAXSUB,MAXOPS)

      CHARACTER*8 SUB_NAME
      DOUBLE PRECISION R_TMP(3,3,MAXOPS)

      CHARACTER*8 NAME_LIST(50)
      INTEGER NNAME

      LOGICAL FOUND
      INTEGER I, J

C ---- EXTERNAL ----
      LOGICAL IS_LINEAR

C=================================================================
C 0. LINEAR MOLECULE CHECK (EARLY EXIT)
C=================================================================
      ISLIN = IS_LINEAR(NAT, X, Y, Z, TOL_STRICT)

      IF (ISLIN) THEN

C        inversion matrix
         RINV(1,1) = -1.0D0
         RINV(1,2) =  0.0D0
         RINV(1,3) =  0.0D0
         RINV(2,1) =  0.0D0
         RINV(2,2) = -1.0D0
         RINV(2,3) =  0.0D0
         RINV(3,1) =  0.0D0
         RINV(3,2) =  0.0D0
         RINV(3,3) = -1.0D0

         TOL_LOOSE = 1.0D2 * TOL_STRICT

         CALL TEST_OPERATION(NAT, SYMBOL, X, Y, Z,
     &                      RINV, TOL_LOOSE,
     &                      OKINV, DELINV)

         IF (OKINV) THEN
            GROUP_NAME = 'Dinfh'
         ELSE
            GROUP_NAME = 'Cinfv'
         END IF

         SYM_STATUS = 0
         DELTA_MAX  = DELINV
         RETURN
      END IF

C=================================================================
C 1. BUILD FULL SYMMETRY GROUP
C=================================================================
      CALL BUILD_SYMMETRY_GROUP
     & (NAT, SYMBOL, X, Y, Z,
     &  TOL_STRICT,
     &  R_GROUP, NOPS,
     &  DELTA_MAX, STATUS)

C=================================================================
C 2. ANALYZE GROUP (INVARIANTS)
C=================================================================
      CALL ANALYZE_SYMMETRY_GROUP
     & (R_GROUP, NOPS,
     &  MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN)

C=================================================================
C 3. IDENTIFY GROUP FAMILY
C=================================================================
      CALL IDENTIFY_GROUP_FAMILY
     & (MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  FAMILY)

C=================================================================
C 4. ASSIGN FINAL POINT GROUP NAME
C=================================================================
      CALL ASSIGN_POINT_GROUP_NAME
     & (FAMILY, MAX_CN,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  GROUP_NAME)

C=================================================================
C 5. SYMMETRY QUALITY FLAG
C=================================================================
      IF (DELTA_MAX .LT. TOL_STRICT) THEN
         SYM_STATUS = 0
      ELSE IF (DELTA_MAX .LT. 1.0D2*TOL_STRICT) THEN
         SYM_STATUS = 1
      ELSE
         SYM_STATUS = 2
      END IF

C     GICForge uses this routine as a library call.  Return the point group
C     and quality flag here; subgroup enumeration is kept below for the
C     standalone symmetry workflow but is intentionally not run by GICForge.
      RETURN

C=================================================================
C 5b. ENUMERATE SUBGROUPS
C=================================================================
      CALL ENUMERATE_SUBGROUPS
     & (R_GROUP, NOPS,
     &  R_SUB_ALL, NSUB_ALL,
     &  NFOUND,
     &  1.0D-6)

C=================================================================
C 5c. FILTER AND PRINT SUBGROUPS
C=================================================================
      CALL FILTER_CONJUGATE_SUBGROUPS
     & (R_SUB_ALL, NSUB_ALL, NFOUND,
     &  R_GROUP, NOPS,
     &  R_SUB_UNIQ, NSUB_UNIQ, NUNIQ,
     &  1.0D-6)

      WRITE(*,*) ' '
      WRITE(*,*) 'UNIQUE SUBGROUPS (up to conjugacy):'

C      CALL PRINT_SUBGROUPS
C    & (R_SUB_UNIQ, NSUB_UNIQ, NUNIQ)

      NNAME = 0

      IF (GROUP_NAME .NE. 'C1') THEN
         WRITE(*,'(A,A8)') ' Full Group: ', GROUP_NAME
      END IF

      DO I = 1, NUNIQ

         CALL COPY_SUBGROUP_FROM_4D
     &        (R_SUB_UNIQ, I, R_TMP, NSUB_UNIQ(I))

         CALL CLASSIFY_SUBGROUP
     &        (R_TMP, NSUB_UNIQ(I), SUB_NAME)

         FOUND = .FALSE.
         DO J = 1, NNAME
            IF (SUB_NAME .EQ. NAME_LIST(J)) THEN
               FOUND = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. FOUND) THEN
            NNAME = NNAME + 1
            NAME_LIST(NNAME) = SUB_NAME

            WRITE(*,'(A,A8)') ' Representative subgroup: ', SUB_NAME
         END IF

      END DO

      RETURN
      END

*Deck BUILD_SYMMETRY_GROUP
      SUBROUTINE BUILD_SYMMETRY_GROUP
     & (NAT, SYMBOL, X, Y, Z,
     &  TOL_STRICT,
     &  R_GROUP, NOPS,
     &  DELTA_MAX, STATUS)
      IMPLICIT NONE

C     ---- INPUT ----
      INTEGER NAT
      CHARACTER*2 SYMBOL(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION TOL_STRICT

C     ---- OUTPUT ----
      INTEGER NOPS, STATUS
      DOUBLE PRECISION R_GROUP(3,3,*)
      DOUBLE PRECISION DELTA_MAX

C     ---- LOCALS ----
      INTEGER NCAND, I, J
      DOUBLE PRECISION R_CAND(3,3,200)
      DOUBLE PRECISION DELTA_OP
      DOUBLE PRECISION TOL_LOOSE
      LOGICAL OK, IS_NEW

C     ---- EXTERNAL ----
      LOGICAL SAME_OP

      TOL_LOOSE = 1.0D2 * TOL_STRICT
      DELTA_MAX = 0.0D0
      NOPS      = 0

C     1. candidate generators
      CALL GENERATE_CANDIDATE_OPS(R_CAND, NCAND)

C     2. test + dedup
      DO I = 1, NCAND

         CALL TEST_OPERATION(NAT, SYMBOL, X, Y, Z,
     &                      R_CAND(:,:,I), TOL_LOOSE,
     &                      OK, DELTA_OP)

         IF (.NOT. OK) CYCLE

         IS_NEW = .TRUE.
         DO J = 1, NOPS
            IF (SAME_OP(R_CAND(:,:,I), R_GROUP(:,:,J), 1.0D-6)) THEN
               IS_NEW = .FALSE.
               EXIT
            END IF
         END DO

         IF (IS_NEW) THEN
            NOPS = NOPS + 1
            R_GROUP(:,:,NOPS) = R_CAND(:,:,I)
            DELTA_MAX = MAX(DELTA_MAX, DELTA_OP)
         END IF

      END DO

C     3. closure
      CALL COMPLETE_GROUP(R_GROUP, NOPS, 1.0D-6)

C     4. check closure (warning only)
      CALL CHECK_GROUP_CLOSURE(R_GROUP, NOPS, OK)
      IF (.NOT. OK) THEN
         STATUS = 1
      ELSE
         STATUS = 0
      END IF

      RETURN
      END

*Deck ANALYZE_SYMMETRY_GROUP
      SUBROUTINE ANALYZE_SYMMETRY_GROUP
     & (R_GROUP, NOPS,
     &  MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN)
      IMPLICIT NONE

C=================================================================
C  Analyze symmetry operations and extract group invariants.
C  Correct for all molecular point groups (C, D, T, O, I).
C=================================================================

C ---- INPUT ----
      INTEGER NOPS
      DOUBLE PRECISION R_GROUP(3,3,*)

C ---- OUTPUT ----
      INTEGER MAX_CN, NC2, NC3, NC4, NC5
      LOGICAL HAS_I, HAS_SIGMA, HAS_SN

C ---- LOCALS ----
      INTEGER I, NROT
      DOUBLE PRECISION DET, TRACE
      DOUBLE PRECISION TOL

C ---- EXTERNAL ----
      DOUBLE PRECISION DET3

      TOL = 1.0D-6

C ---- INITIALIZE ----
      MAX_CN    = 0
      NC2 = 0
      NC3 = 0
      NC4 = 0
      NC5 = 0
      HAS_I     = .FALSE.
      HAS_SIGMA = .FALSE.
      HAS_SN    = .FALSE.

C=================================================================
C  Loop over symmetry operations
C=================================================================
      DO I = 1, NOPS

         DET = DET3(R_GROUP(:,:,I))
         TRACE = R_GROUP(1,1,I)
     &         + R_GROUP(2,2,I)
     &         + R_GROUP(3,3,I)

C---------------------------------------------------------------
C  Improper operations (det = -1)
C---------------------------------------------------------------
         IF (DET .LT. -0.5D0) THEN

C           Inversion: trace = -3
            IF (DABS(TRACE + 3.0D0) .LT. TOL) THEN
               HAS_I = .TRUE.

C           Proper reflection plane: trace = +1
            ELSE IF (DABS(TRACE - 1.0D0) .LT. TOL) THEN
               HAS_SIGMA = .TRUE.

C           Other improper rotations (Sn)
            ELSE
               HAS_SN = .TRUE.
            END IF

C---------------------------------------------------------------
C  Proper rotations (det = +1)
C---------------------------------------------------------------
         ELSE

            CALL GET_ROTATION_ORDER(R_GROUP(:,:,I), NROT)

            IF (NROT .GT. 0) THEN
               MAX_CN = MAX(MAX_CN, NROT)
               IF (NROT .EQ. 2) NC2 = NC2 + 1
               IF (NROT .EQ. 3) NC3 = NC3 + 1
               IF (NROT .EQ. 4) NC4 = NC4 + 1
               IF (NROT .EQ. 5) NC5 = NC5 + 1
            END IF

         END IF

      END DO

      RETURN
      END

*Deck IDENTIFY_GROUP_FAMILY
      SUBROUTINE IDENTIFY_GROUP_FAMILY
     & (MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  FAMILY)
      IMPLICIT NONE

C=================================================================
C  Identify abstract point-group family:
C     1 = C (cyclic)
C     2 = D (dihedral)
C     3 = T (tetrahedral: T, Td, Th)
C     4 = O (octahedral: O, Oh)
C     5 = I (icosahedral: I, Ih)
C=================================================================

C ---- INPUT ----
      INTEGER MAX_CN, NC2, NC3, NC4, NC5
      LOGICAL HAS_I, HAS_SIGMA, HAS_SN

C ---- OUTPUT ----
      INTEGER FAMILY

C ---- DEFAULT ----
      FAMILY = 1

C=================================================================
C  Icosahedral family
C=================================================================
C  Presence of any genuine C5 axis
      IF (NC5 .GT. 0) THEN
         FAMILY = 5
         RETURN
      END IF

C=================================================================
C  Octahedral family
C=================================================================
C  Requires C4 and C3 axes
      IF (NC4 .GT. 0 .AND. NC3 .GE. 4) THEN
         FAMILY = 4
         RETURN
      END IF

C=================================================================
C  Tetrahedral family
C=================================================================
C  Requires multiple NON-COLLINEAR C3 axes:
C    - not a single principal C3 with powers
C    - NC3 must be large
C    - MAX_CN <= 3 (no higher principal axis)
      IF (NC3 .GE. 4 .AND. MAX_CN .LE. 3) THEN
         FAMILY = 3
         RETURN
      END IF

C=================================================================
C  Dihedral family
C=================================================================
C  Requires:
C    - principal axis n >= 2
C    - multiple C2 perpendicular axes
      IF (MAX_CN .GE. 2 .AND. NC2 .GE. MAX_CN) THEN
         FAMILY = 2
         RETURN
      END IF

C=================================================================
C  Cyclic family (fallback)
C=================================================================
      FAMILY = 1
      RETURN
      END

*Deck ASSIGN_POINT_GROUP_NAME
      SUBROUTINE ASSIGN_POINT_GROUP_NAME
     & (FAMILY, MAX_CN,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  GROUP_NAME)
      IMPLICIT NONE

      INTEGER FAMILY, MAX_CN
      LOGICAL HAS_I, HAS_SIGMA, HAS_SN
      CHARACTER*8 GROUP_NAME

C     ---- DEFAULT ----
      GROUP_NAME = 'C1'

C=================================================================
C  Cubic families MUST exit immediately
C=================================================================

C     Tetrahedral family
      IF (FAMILY .EQ. 3) THEN
         IF (HAS_I) THEN
            GROUP_NAME = 'Th'
         ELSE
            GROUP_NAME = 'Td'
         END IF
         RETURN
      END IF

C     Octahedral family
      IF (FAMILY .EQ. 4) THEN
         IF (HAS_I) THEN
            GROUP_NAME = 'Oh'
         ELSE
            GROUP_NAME = 'O'
         END IF
         RETURN
      END IF

C     Icosahedral family
      IF (FAMILY .EQ. 5) THEN
         IF (HAS_I) THEN
            GROUP_NAME = 'Ih'
         ELSE
            GROUP_NAME = 'I'
         END IF
         RETURN
      END IF

C=================================================================
C  Dihedral family (n >= 2 only)
C=================================================================

      IF (FAMILY .EQ. 2) THEN
         IF (MAX_CN .LT. 2) THEN
            GROUP_NAME = 'Th'
            RETURN
         END IF

         IF (HAS_SIGMA) THEN
            WRITE(GROUP_NAME,'(A,I1,A)') 'D', MAX_CN, 'h'
         ELSE IF (HAS_SN) THEN
            WRITE(GROUP_NAME,'(A,I1,A)') 'D', MAX_CN, 'd'
         ELSE
            WRITE(GROUP_NAME,'(A,I1)') 'D', MAX_CN
         END IF
         RETURN
      END IF

C=================================================================
C  Cyclic family (LAST, cannot override cubic)
C=================================================================

      IF (FAMILY .EQ. 1) THEN
         IF (MAX_CN .GE. 2) THEN
            IF (HAS_SIGMA) THEN
               WRITE(GROUP_NAME,'(A,I1,A)') 'C', MAX_CN, 'v'
            ELSE IF (HAS_I) THEN
               WRITE(GROUP_NAME,'(A,I1,A)') 'C', MAX_CN, 'h'
            ELSE
               WRITE(GROUP_NAME,'(A,I1)') 'C', MAX_CN
            END IF
         ELSE
            IF (HAS_I) THEN
               GROUP_NAME = 'Ci'
            ELSE IF (HAS_SIGMA) THEN
               GROUP_NAME = 'Cs'
            ELSE
               GROUP_NAME = 'C1'
            END IF
         END IF
         RETURN
      END IF

      RETURN
      END

*Deck GENERATE_CANDIDATE_OPS
      SUBROUTINE GENERATE_CANDIDATE_OPS(R_CAND, NCAND)
      IMPLICIT NONE

C     ============================================================
C     Candidate symmetry operation generators for ALL
C     molecular point groups:
C     Cn, Dn, Sn, Cnv, Dnh, Dnd,
C     Td, Th, O, Oh, I, Ih
C     ============================================================

      INTEGER NCAND
      DOUBLE PRECISION R_CAND(3,3,*)

C     ---- locals ----
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
      DOUBLE PRECISION U(3)
      DOUBLE PRECISION PI, ANG, NORM
      INTEGER I, K

C     ---- allowed principal rotation orders ----
      INTEGER NORD
      PARAMETER (NORD = 5)
      INTEGER ORDERS(NORD)
      DATA ORDERS /2,3,4,5,6/

      PI = 4.0D0 * DATAN(1.0D0)
      NCAND = 0

C     ============================================================
C     Identity
C     ============================================================
      CALL IDMAT(A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     ============================================================
C     Inversion
C     ============================================================
      A = 0.0D0
      A(1,1) = -1.0D0
      A(2,2) = -1.0D0
      A(3,3) = -1.0D0
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     ============================================================
C     Principal-axis rotations Cn (Z axis)
C     ============================================================
      DO I = 1, NORD
         CALL ROTZ(2.0D0 * PI / DBLE(ORDERS(I)), A)
         NCAND = NCAND + 1
         R_CAND(:,:,NCAND) = A
      END DO

C     ============================================================
C     Perpendicular C2 axes (Dn)
C     ============================================================
      DO K = 0, 5
         ANG = DBLE(K) * 30.0D0 * PI / 180.0D0
         U(1) = DCOS(ANG)
         U(2) = DSIN(ANG)
         U(3) = 0.0D0
         CALL ROT_AXIS(U, PI, A)
         NCAND = NCAND + 1
         R_CAND(:,:,NCAND) = A
      END DO

C     ============================================================
C     Horizontal reflection sigma_h
C     ============================================================
      CALL REF_XY(A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     ============================================================
C     Vertical reflections sigma_v (planes containing Z)
C     ============================================================
C     XZ plane
      CALL REF_XZ(A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     YZ plane
      CALL REF_YZ(A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     Rotated sigma_v (45 degrees)
      DO K = 1, 3
         ANG = DBLE(K) * 45.0D0 * PI / 180.0D0
         U(1) = DCOS(ANG)
         U(2) = DSIN(ANG)
         U(3) = 0.0D0
         CALL REF_PLANE(U, A)
         NCAND = NCAND + 1
         R_CAND(:,:,NCAND) = A
      END DO

C     ============================================================
C     Improper rotations Sn
C     ============================================================
      CALL REF_XY(B)
      DO I = 1, NORD
         CALL ROTZ(2.0D0 * PI / DBLE(ORDERS(I)), A)
         CALL MATMUL3(B, A, C)
         NCAND = NCAND + 1
         R_CAND(:,:,NCAND) = C
      END DO

C     ============================================================
C     Cubic group generators
C     ============================================================

C     ---- C4 about X and Y (for O, Oh) ----
      U(1) = 1.0D0
      U(2) = 0.0D0
      U(3) = 0.0D0
      CALL ROT_AXIS(U, PI/2.0D0, A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

      U(1) = 0.0D0
      U(2) = 1.0D0
      U(3) = 0.0D0
      CALL ROT_AXIS(U, PI/2.0D0, A)
      NCAND = NCAND + 1
      R_CAND(:,:,NCAND) = A

C     ---- C3 along body diagonals (for T, Td, O, Oh) ----
      DO I = 1, 4
         U(1) =  1.0D0
         U(2) =  1.0D0
         U(3) =  1.0D0
         IF (I.EQ.2) U(3) = -1.0D0
         IF (I.EQ.3) U(2) = -1.0D0
         IF (I.EQ.4) U(1) = -1.0D0

         NORM = DSQRT(U(1)**2 + U(2)**2 + U(3)**2)
         U(1) = U(1)/NORM
         U(2) = U(2)/NORM
         U(3) = U(3)/NORM

         CALL ROT_AXIS(U, 2.0D0*PI/3.0D0, A)
         NCAND = NCAND + 1
         R_CAND(:,:,NCAND) = A
      END DO

      RETURN
      END

*Deck TEST_OPERATION
      SUBROUTINE TEST_OPERATION
     & (NAT, SYMBOL, X, Y, Z,
     &  ROP, TOL_LOOSE,
     &  OK, DELTA_OP)
      IMPLICIT NONE

C     ============================================================
C     Test a symmetry operation ROP for a molecule already
C     centered at the center of mass and Eckart-oriented.
C
C     Robust for ALL molecular point groups:
C     Cn, Dn, Sn, Cnv, Dnh, Td, Oh, Ih
C     ============================================================

C     ---- INPUT ----
      INTEGER NAT
      CHARACTER*2 SYMBOL(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION ROP(3,3)
      DOUBLE PRECISION TOL_LOOSE

C     ---- OUTPUT ----
      LOGICAL OK
      DOUBLE PRECISION DELTA_OP

C     ---- LOCALS ----
      INTEGER I, IB
      DOUBLE PRECISION XP, YP, ZP
      DOUBLE PRECISION DELTA
      LOGICAL FOUND
      LOGICAL USED(1000)
      DOUBLE PRECISION DET
      DOUBLE PRECISION ORTHTOL

C     ---- EXTERNAL ----
      DOUBLE PRECISION DET3

      ORTHTOL = 1.0D-6

C     ---- inizializzazione ----
      OK        = .TRUE.
      DELTA_OP = 0.0D0

      DO I = 1, NAT
         USED(I) = .FALSE.
      END DO

C     ============================================================
C     Check orthogonality and determinant
C     ============================================================
      CALL CHECK_ORTHOGONAL(ROP, ORTHTOL, OK)
      IF (.NOT. OK) RETURN

      DET = DET3(ROP)
      IF (DABS(DABS(DET) - 1.0D0) .GT. ORTHTOL) THEN
         OK = .FALSE.
         RETURN
      END IF

C     ============================================================
C     Loop over atoms
C     ============================================================
      DO I = 1, NAT

C        apply operation
         XP = ROP(1,1)*X(I) + ROP(1,2)*Y(I) + ROP(1,3)*Z(I)
         YP = ROP(2,1)*X(I) + ROP(2,2)*Y(I) + ROP(2,3)*Z(I)
         ZP = ROP(3,1)*X(I) + ROP(3,2)*Y(I) + ROP(3,3)*Z(I)

C        iniettività obbligatoria (Td / Oh / Ih)
         CALL MATCH_ATOM_UNUSED
     &      (XP, YP, ZP, SYMBOL(I), NAT,
     &       SYMBOL, X, Y, Z,
     &       USED, TOL_LOOSE,
     &       IB, FOUND, DELTA)

         IF (.NOT. FOUND) THEN
            OK = .FALSE.
            RETURN
         END IF

         USED(IB) = .TRUE.
         DELTA_OP = MAX(DELTA_OP, DELTA)

      END DO

      RETURN
      END

*Deck CLASSIFY_OP
      SUBROUTINE CLASSIFY_OP
     & (R,
     &  HAS_I,
     &  HAS_SIGMA,
     &  HAS_SN,
     &  HAS_CN,
     &  HAS_C2,
     &  MAX_CN)
      IMPLICIT NONE

C     ============================================================
C     General classification of a symmetry operation R.
C     Works for ALL molecular point groups.
C     ============================================================

C     ---- INPUT ----
      DOUBLE PRECISION R(3,3)
      INTEGER MAX_CN

C     ---- OUTPUT ----
      LOGICAL HAS_I
      LOGICAL HAS_SIGMA
      LOGICAL HAS_SN
      LOGICAL HAS_C2
      LOGICAL HAS_CN(MAX_CN)

C     ---- LOCALS ----
      DOUBLE PRECISION TRACE, COSTH, DET
      DOUBLE PRECISION TOL
      INTEGER NROT

C     ---- EXTERNAL ----
      DOUBLE PRECISION DET3

      TOL = 1.0D-6

C     ---- initialize ----
      NROT = 0

C     ============================================================
C     Determinant
C     ============================================================
      DET = DET3(R)

C     ============================================================
C     Inversion
C     ============================================================
      IF (DET .LT. 0.0D0) THEN
         TRACE = R(1,1) + R(2,2) + R(3,3)
         IF (DABS(TRACE + 3.0D0) .LT. TOL) THEN
            HAS_I = .TRUE.
            RETURN
         END IF
      END IF

C     ============================================================
C     Reflection (sigma)
C     det = -1, trace = +1
C     ============================================================
      IF (DET .LT. 0.0D0) THEN
         TRACE = R(1,1) + R(2,2) + R(3,3)
         IF (DABS(TRACE - 1.0D0) .LT. TOL) THEN
            HAS_SIGMA = .TRUE.
            RETURN
         END IF
      END IF

C     ============================================================
C     Proper rotations Cn
C     ============================================================
      IF (DET .GT. 0.0D0) THEN
         CALL GET_ROTATION_ORDER(R, NROT)
         IF (NROT .GE. 2 .AND. NROT .LE. MAX_CN) THEN
            HAS_CN(NROT) = .TRUE.
            IF (NROT .EQ. 2) HAS_C2 = .TRUE.
            RETURN
         END IF
      END IF

C     ============================================================
C     Improper rotations Sn
C     det = -1, trace = 1 + 2 cos(2π/n)
C     ============================================================
      IF (DET .LT. 0.0D0) THEN
         TRACE = R(1,1) + R(2,2) + R(3,3)
         COSTH = 0.5D0 * (TRACE - 1.0D0)
         IF (DABS(COSTH) .LE. 1.0D0) THEN
            HAS_SN = .TRUE.
         END IF
      END IF

      RETURN
      END

*Deck IDMAT
      SUBROUTINE IDMAT(A)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3)
      INTEGER I, J

      DO I = 1,3
         DO J = 1,3
            A(I,J) = 0.0D0
         END DO
      END DO

      DO I = 1,3
         A(I,I) = 1.0D0
      END DO

      RETURN
      END

*Deck GET_ROTATION_ORDER
      SUBROUTINE GET_ROTATION_ORDER(R, NROT)
      IMPLICIT NONE

C=================================================================
C  Determine the order n of a proper rotation matrix R (C_n).
C  Returns NROT = 0 if R is identity, improper, or non-integer.
C  Robust for n up to NMAX.
C=================================================================

C ---- INPUT ----
      DOUBLE PRECISION R(3,3)

C ---- OUTPUT ----
      INTEGER NROT

C ---- PARAMETERS ----
      INTEGER NMAX
      PARAMETER (NMAX = 12)

C ---- LOCALS ----
      DOUBLE PRECISION TRACE, COSTH
      DOUBLE PRECISION COSN, ERR
      DOUBLE PRECISION TOL
      DOUBLE PRECISION PI
      INTEGER N, K
      DOUBLE PRECISION DET3
      DOUBLE PRECISION P(3,3), TMP(3,3)

      PI  = 3.141592653589793D0
      TOL = 1.0D-6
      NROT = 0

C-----------------------------------------------------------------
C  Exclude improper operations
C-----------------------------------------------------------------
      IF (DET3(R) .LT. 0.0D0) RETURN

C-----------------------------------------------------------------
C  Compute trace
C-----------------------------------------------------------------
      TRACE = R(1,1) + R(2,2) + R(3,3)

C-----------------------------------------------------------------
C  Exclude identity explicitly
C-----------------------------------------------------------------
      IF (DABS(TRACE - 3.0D0) .LT. TOL) RETURN

C-----------------------------------------------------------------
C  Compute cos(theta)
C-----------------------------------------------------------------
      COSTH = 0.5D0 * (TRACE - 1.0D0)

C     Clamp for numerical safety
      IF (COSTH .GT.  1.0D0) COSTH =  1.0D0
      IF (COSTH .LT. -1.0D0) COSTH = -1.0D0

C=================================================================
C  Try integer orders
C=================================================================
      DO 20 N = 2, NMAX

         COSN = DCOS(2.0D0 * PI / DBLE(N))
         ERR  = DABS(COSTH - COSN)

C        Quick reject
         IF (ERR .GT. TOL) GOTO 20

C        Verify R^N ≈ I
         CALL IDMAT(P)

         DO 10 K = 1, N
            CALL MATMUL3(P, R, TMP)

C           Inline copy TMP -> P
            P(1,1) = TMP(1,1)
            P(1,2) = TMP(1,2)
            P(1,3) = TMP(1,3)
            P(2,1) = TMP(2,1)
            P(2,2) = TMP(2,2)
            P(2,3) = TMP(2,3)
            P(3,1) = TMP(3,1)
            P(3,2) = TMP(3,2)
            P(3,3) = TMP(3,3)
 10      CONTINUE

         IF (DABS(P(1,1)-1.0D0).LT.TOL .AND.
     &       DABS(P(2,2)-1.0D0).LT.TOL .AND.
     &       DABS(P(3,3)-1.0D0).LT.TOL .AND.
     &       DABS(P(1,2)).LT.TOL .AND.
     &       DABS(P(1,3)).LT.TOL .AND.
     &       DABS(P(2,1)).LT.TOL .AND.
     &       DABS(P(2,3)).LT.TOL .AND.
     &       DABS(P(3,1)).LT.TOL .AND.
     &       DABS(P(3,2)).LT.TOL) THEN

            NROT = N
            RETURN
         END IF

 20   CONTINUE

      RETURN
      END

*Deck IS_IDMAT_TOL
      LOGICAL FUNCTION IS_IDMAT_TOL(A, TOL)
      IMPLICIT NONE

C     Check whether A is identity matrix within tolerance TOL

      DOUBLE PRECISION A(3,3), TOL
      INTEGER I, J

      IS_IDMAT_TOL = .TRUE.

      DO I = 1, 3
         DO J = 1, 3
            IF (I .EQ. J) THEN
               IF (DABS(A(I,J) - 1.0D0) .GT. TOL) THEN
                  IS_IDMAT_TOL = .FALSE.
                  RETURN
               END IF
            ELSE
               IF (DABS(A(I,J)) .GT. TOL) THEN
                  IS_IDMAT_TOL = .FALSE.
                  RETURN
               END IF
            END IF
         END DO
      END DO

      RETURN
      END

*Deck CHECK_GROUP_CLOSURE
      SUBROUTINE CHECK_GROUP_CLOSURE(R_VALID, NOPS, OK)
      IMPLICIT NONE

C     ============================================================
C     Robust closure test for a set of symmetry operations.
C     Valid for ALL molecular point groups (Td, Oh, Ih included).
C     ============================================================

C     ---- INPUT ----
      INTEGER NOPS
      DOUBLE PRECISION R_VALID(3,3,*)

C     ---- OUTPUT ----
      LOGICAL OK

C     ---- LOCALS ----
      INTEGER I, J, K
      DOUBLE PRECISION P(3,3)
      LOGICAL FOUND
      DOUBLE PRECISION TOL

C     ---- EXTERNAL ----
      LOGICAL SAME_OP

      TOL = 1.0D-6
      OK  = .TRUE.

C     ============================================================
C     Loop over operation pairs
C     ============================================================
      DO I = 1, NOPS
         DO J = 1, NOPS

C           P = Ri * Rj
            CALL MATMUL3(R_VALID(:,:,I), R_VALID(:,:,J), P)

C           Clean numerical noise (important!)
            CALL CLEAN_ROTATION(P, TOL)

C           Search P in the valid set
            FOUND = .FALSE.
            DO K = 1, NOPS
               IF (SAME_OP(P, R_VALID(:,:,K), TOL)) THEN
                  FOUND = .TRUE.
                  GOTO 20
               END IF
            END DO
 20         CONTINUE

            IF (.NOT. FOUND) THEN
               OK = .FALSE.
               RETURN
            END IF

         END DO
      END DO

      RETURN
      END

*Deck CLEAN_ROTATION
      SUBROUTINE CLEAN_ROTATION(R, TOL)
      IMPLICIT NONE
      DOUBLE PRECISION R(3,3), TOL
      INTEGER I, J

      DO I = 1, 3
         DO J = 1, 3
            IF (DABS(R(I,J)) .LT. TOL) R(I,J) = 0.0D0
            IF (DABS(R(I,J) - 1.0D0) .LT. TOL) R(I,J) = 1.0D0
            IF (DABS(R(I,J) + 1.0D0) .LT. TOL) R(I,J) = -1.0D0
         END DO
      END DO
      RETURN
      END

*Deck SAME_OP
      LOGICAL FUNCTION SAME_OP(A, B, TOL)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), B(3,3), TOL
      DOUBLE PRECISION D
      INTEGER I, J

      D = 0.0D0
      DO I = 1, 3
         DO J = 1, 3
            D = D + (A(I,J) - B(I,J))**2
         END DO
      END DO

      D = DSQRT(D)
      SAME_OP = (D .LT. TOL)
      RETURN
      END

*Deck COMPLETE_GROUP
      SUBROUTINE COMPLETE_GROUP(R, NOPS, TOL)
      IMPLICIT NONE

C     ============================================================
C     Complete a symmetry group by closure.
C     Robust for ALL molecular point groups (Td, Oh, Ih).
C     ============================================================

      INTEGER NOPS
      DOUBLE PRECISION R(3,3,*), TOL

C     ---- LOCALS ----
      INTEGER I, J, K
      DOUBLE PRECISION P(3,3)
      LOGICAL FOUND, CHANGED, OK
      INTEGER NOPS_OLD
      INTEGER MAXOPS
      PARAMETER (MAXOPS = 200)

C     ---- EXTERNAL ----
      LOGICAL SAME_OP
      DOUBLE PRECISION DET3

 100  CONTINUE
      CHANGED  = .FALSE.
      NOPS_OLD = NOPS

      DO I = 1, NOPS_OLD
         DO J = 1, NOPS_OLD

C           Product
            CALL MATMUL3(R(:,:,I), R(:,:,J), P)

C           Clean numerical noise
            CALL CLEAN_ROTATION(P, TOL)

C           Orthogonality check
            CALL CHECK_ORTHOGONAL(P, TOL, OK)
            IF (.NOT. OK) GOTO 10

C           Determinant check
            IF (DABS(DABS(DET3(P)) - 1.0D0) .GT. TOL) GOTO 10

C           Search in current set
            FOUND = .FALSE.
            DO K = 1, NOPS
               IF (SAME_OP(P, R(:,:,K), TOL)) THEN
                  FOUND = .TRUE.
                  EXIT
               END IF
            END DO

C           Add new operation
            IF (.NOT. FOUND) THEN
               IF (NOPS .GE. MAXOPS) THEN
                  WRITE(*,*) 'ERROR: group too large'
                  STOP
               END IF
               NOPS = NOPS + 1
               R(:,:,NOPS) = P
               CHANGED = .TRUE.
            END IF

 10         CONTINUE
         END DO
      END DO

      IF (CHANGED) GOTO 100

      RETURN
      END

*Deck IS_LINEAR
      LOGICAL FUNCTION IS_LINEAR(NAT, X, Y, Z, TOL)
      IMPLICIT NONE
      INTEGER NAT, I
      DOUBLE PRECISION X(*), Y(*), Z(*), TOL
      DOUBLE PRECISION U(3), R(3), CROSS(3), NORMU, NORMC

C     trova atomo più lontano
      NORMU = 0.0D0
      DO I = 1, NAT
         R(1) = X(I)
         R(2) = Y(I)
         R(3) = Z(I)
         IF (R(1)*R(1)+R(2)*R(2)+R(3)*R(3) .GT. NORMU) THEN
            U = R
            NORMU = R(1)*R(1)+R(2)*R(2)+R(3)*R(3)
         END IF
      END DO

      NORMU = DSQRT(NORMU)
      IF (NORMU .LT. TOL) THEN
         IS_LINEAR = .TRUE.
         RETURN
      END IF

      U = U / NORMU

C     verifica collinearità
      DO I = 1, NAT
         R(1) = X(I)
         R(2) = Y(I)
         R(3) = Z(I)

C        cross product U x R
         CROSS(1) = U(2)*R(3) - U(3)*R(2)
         CROSS(2) = U(3)*R(1) - U(1)*R(3)
         CROSS(3) = U(1)*R(2) - U(2)*R(1)

         NORMC = DSQRT(CROSS(1)**2 + CROSS(2)**2 + CROSS(3)**2)
         IF (NORMC .GT. TOL) THEN
            IS_LINEAR = .FALSE.
            RETURN
         END IF
      END DO

      IS_LINEAR = .TRUE.
      RETURN
      END

*Deck CHECK_ORTHOGONAL
      SUBROUTINE CHECK_ORTHOGONAL(R, TOL, OK)
      IMPLICIT NONE
      DOUBLE PRECISION R(3,3), TOL
      LOGICAL OK
      DOUBLE PRECISION S
      INTEGER I, J, K

      OK = .TRUE.

      DO I = 1, 3
         DO J = 1, 3
            S = 0.0D0
            DO K = 1, 3
               S = S + R(K,I)*R(K,J)
            END DO
            IF (I .EQ. J) THEN
               IF (DABS(S - 1.0D0) .GT. TOL) OK = .FALSE.
            ELSE
               IF (DABS(S) .GT. TOL) OK = .FALSE.
            END IF
         END DO
      END DO

      RETURN
      END

*Deck MATCH_ATOM_UNUSED
      SUBROUTINE MATCH_ATOM_UNUSED
     & (XP, YP, ZP, SYM, NAT,
     &  SYMBOL, X, Y, Z,
     &  USED, TOL,
     &  IBEST, FOUND, DELTA)
      IMPLICIT NONE

      INTEGER NAT, IBEST, I
      CHARACTER*2 SYM, SYMBOL(*)
      DOUBLE PRECISION XP, YP, ZP
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION TOL, DELTA, DMIN, D
      LOGICAL USED(*), FOUND

      FOUND = .FALSE.
      IBEST = -1
      DMIN  = 1.0D30

      DO I = 1, NAT
         IF (.NOT. USED(I)) THEN
            IF (SYMBOL(I) .EQ. SYM) THEN
               D = (XP-X(I))**2 + (YP-Y(I))**2 + (ZP-Z(I))**2
               D = DSQRT(D)
               IF (D .LT. TOL .AND. D .LT. DMIN) THEN
                  DMIN  = D
                  IBEST = I
                  FOUND = .TRUE.
               END IF
            END IF
         END IF
      END DO

      IF (FOUND) THEN
         DELTA = DMIN
      ELSE
         DELTA = 1.0D30
      END IF

      RETURN
      END

*Deck REF_PLANE
      SUBROUTINE REF_PLANE(U, R)
      IMPLICIT NONE
      DOUBLE PRECISION U(3), R(3,3)
      INTEGER I, J

C     R = I - 2 * U * U^T

      DO I = 1, 3
         DO J = 1, 3
            R(I,J) = -2.0D0 * U(I) * U(J)
         END DO
      END DO

      R(1,1) = R(1,1) + 1.0D0
      R(2,2) = R(2,2) + 1.0D0
      R(3,3) = R(3,3) + 1.0D0

      RETURN
      END

C----------------------------------------------------------------
C----------------------------------------------------------------
*Deck BUILD_SUBGROUP
      SUBROUTINE BUILD_SUBGROUP
     & (R_GEN, NGEN,
     &  R_GROUP, NOPS,
     &  R_SUB, NSUB,
     &  TOL)
      IMPLICIT NONE

C=================================================================
C  Build the subgroup generated by R_GEN inside the full group G.
C=================================================================

C ---- INPUT ----
      INTEGER NGEN, NOPS
      DOUBLE PRECISION R_GEN(3,3,NGEN)
      DOUBLE PRECISION R_GROUP(3,3,NOPS)
      DOUBLE PRECISION TOL

C ---- OUTPUT ----
      INTEGER NSUB
      DOUBLE PRECISION R_SUB(3,3,*)

C ---- LOCALS ----
      INTEGER I, J, K
      DOUBLE PRECISION P(3,3)
      LOGICAL FOUND, CHANGED
      LOGICAL SAME_OP

C-----------------------------------------------------------------
C  Initialize subgroup with generators
C-----------------------------------------------------------------
      NSUB = 0

      DO I = 1, NGEN
         FOUND = .FALSE.
         DO J = 1, NSUB
            IF (SAME_OP(R_GEN(:,:,I), R_SUB(:,:,J), TOL)) THEN
               FOUND = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. FOUND) THEN
            NSUB = NSUB + 1
            R_SUB(:,:,NSUB) = R_GEN(:,:,I)
         END IF
      END DO

C-----------------------------------------------------------------
C  Closure inside G
C-----------------------------------------------------------------
 100  CONTINUE
      CHANGED = .FALSE.

      DO I = 1, NSUB
         DO J = 1, NSUB

            CALL MATMUL3(R_SUB(:,:,I), R_SUB(:,:,J), P)

C           check P belongs to full group G
            FOUND = .FALSE.
            DO K = 1, NOPS
               IF (SAME_OP(P, R_GROUP(:,:,K), TOL)) THEN
                  FOUND = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. FOUND) GOTO 200

C           check P already in subgroup
            FOUND = .FALSE.
            DO K = 1, NSUB
               IF (SAME_OP(P, R_SUB(:,:,K), TOL)) THEN
                  FOUND = .TRUE.
                  EXIT
               END IF
            END DO

            IF (.NOT. FOUND) THEN
               NSUB = NSUB + 1
               R_SUB(:,:,NSUB) = P
               CHANGED = .TRUE.
            END IF

         END DO
      END DO

      IF (CHANGED) GOTO 100
      RETURN

C-----------------------------------------------------------------
C  Product escaped G → generators do not define a subgroup of G
C-----------------------------------------------------------------
 200  CONTINUE
      NSUB = 0
      RETURN
      END

*Deck ENUMERATE_SUBGROUPS
      SUBROUTINE ENUMERATE_SUBGROUPS
     & (R_GROUP, NOPS,
     &  R_SUB_ALL, NSUB_ALL,
     &  NFOUND,
     &  TOL)
      IMPLICIT NONE

      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)

C ---- INPUT ----
      INTEGER NOPS
      DOUBLE PRECISION R_GROUP(3,3,MAXOPS)
      DOUBLE PRECISION TOL

C ---- OUTPUT ----
      INTEGER NFOUND
      INTEGER NSUB_ALL(MAXSUB)
      DOUBLE PRECISION R_SUB_ALL(3,3,MAXSUB,MAXOPS)

C ---- LOCALS ----
      INTEGER I, J
      INTEGER NGEN, NSUB
      DOUBLE PRECISION R_GEN(3,3,3)
      DOUBLE PRECISION R_SUB(3,3,MAXOPS)

      NFOUND = 0

C ---- one-generator subgroups ----
      DO I = 1, NOPS

         R_GEN(:,:,1) = R_GROUP(:,:,1)
         R_GEN(:,:,2) = R_GROUP(:,:,I)
         NGEN = 2

         CALL BUILD_SUBGROUP
     &        (R_GEN, NGEN,
     &         R_GROUP, NOPS,
     &         R_SUB, NSUB,
     &         TOL)

         IF (NSUB .GT. 0) THEN
            CALL ADD_SUBGROUP_IF_NEW
     &           (R_SUB, NSUB,
     &            R_SUB_ALL, NSUB_ALL,
     &            NFOUND, TOL)
         END IF

      END DO

C ---- two-generator subgroups ----
      DO I = 1, NOPS
         DO J = I+1, NOPS

            R_GEN(:,:,1) = R_GROUP(:,:,1)
            R_GEN(:,:,2) = R_GROUP(:,:,I)
            R_GEN(:,:,3) = R_GROUP(:,:,J)
            NGEN = 3

            CALL BUILD_SUBGROUP
     &           (R_GEN, NGEN,
     &            R_GROUP, NOPS,
     &            R_SUB, NSUB,
     &            TOL)

            IF (NSUB .GT. 0) THEN
               CALL ADD_SUBGROUP_IF_NEW
     &              (R_SUB, NSUB,
     &               R_SUB_ALL, NSUB_ALL,
     &               NFOUND, TOL)
            END IF

         END DO
      END DO

      RETURN
      END

*Deck ADD_SUBGROUP_IF_NEW
      SUBROUTINE ADD_SUBGROUP_IF_NEW
     & (R_SUB, NSUB,
     &  R_SUB_ALL, NSUB_ALL,
     &  NFOUND,
     &  TOL)
      IMPLICIT NONE

      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)

      INTEGER NSUB, NFOUND
      INTEGER NSUB_ALL(MAXSUB)
      DOUBLE PRECISION R_SUB(3,3,MAXOPS)
      DOUBLE PRECISION R_SUB_ALL(3,3,MAXSUB,MAXOPS)
      DOUBLE PRECISION TOL

      INTEGER I, J
      LOGICAL SAME_OP
      LOGICAL FOUND

      DO I = 1, NFOUND

         IF (NSUB_ALL(I) .NE. NSUB) GOTO 10

         FOUND = .TRUE.
         DO J = 1, NSUB
            IF (.NOT. SAME_OP(R_SUB(:,:,J),
     &                         R_SUB_ALL(:,:,I,J), TOL)) THEN
               FOUND = .FALSE.
               EXIT
            END IF
         END DO

         IF (FOUND) RETURN
 10      CONTINUE
      END DO

      NFOUND = NFOUND + 1
      NSUB_ALL(NFOUND) = NSUB

      DO J = 1, NSUB
         R_SUB_ALL(:,:,NFOUND,J) = R_SUB(:,:,J)
      END DO

      RETURN
      END

*Deck PRINT_SUBGROUPS
      SUBROUTINE PRINT_SUBGROUPS
     & (R_SUB_ALL, NSUB_ALL, NFOUND)
      IMPLICIT NONE

      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)

      INTEGER NFOUND
      INTEGER NSUB_ALL(MAXSUB)
      DOUBLE PRECISION R_SUB_ALL(3,3,MAXSUB,MAXOPS)

      INTEGER I, J
      DOUBLE PRECISION RLOC(3,3)

      WRITE(*,*) ' '
      WRITE(*,*) 'SUBGROUPS FOUND: ', NFOUND
      WRITE(*,*) '---------------------------'

      DO I = 1, NFOUND
         WRITE(*,'(A,I3,A,I4)') ' Subgroup ', I,
     &        '  order = ', NSUB_ALL(I)

         DO J = 1, NSUB_ALL(I)

            RLOC(1,1) = R_SUB_ALL(1,1,I,J)
            RLOC(1,2) = R_SUB_ALL(1,2,I,J)
            RLOC(1,3) = R_SUB_ALL(1,3,I,J)
            RLOC(2,1) = R_SUB_ALL(2,1,I,J)
            RLOC(2,2) = R_SUB_ALL(2,2,I,J)
            RLOC(2,3) = R_SUB_ALL(2,3,I,J)
            RLOC(3,1) = R_SUB_ALL(3,1,I,J)
            RLOC(3,2) = R_SUB_ALL(3,2,I,J)
            RLOC(3,3) = R_SUB_ALL(3,3,I,J)

            CALL PRINT_OP(RLOC)
         END DO

         WRITE(*,*) '---------------------------'
      END DO

      RETURN
      END

*Deck PRINT_OP
      SUBROUTINE PRINT_OP(R)
      IMPLICIT NONE
      DOUBLE PRECISION R(3,3)

      WRITE(*,'(3F12.6)') R(1,1), R(1,2), R(1,3)
      WRITE(*,'(3F12.6)') R(2,1), R(2,2), R(2,3)
      WRITE(*,'(3F12.6)') R(3,1), R(3,2), R(3,3)
      WRITE(*,*) ' '
      RETURN
      END

*Deck CONJUGATE SUBGROUP
      LOGICAL FUNCTION CONJUGATE_SUBGROUP
     & (R_SUB1, NS1,
     &  R_SUB2, NS2,
     &  R_GROUP, NG,
     &  TOL)
      IMPLICIT NONE

C=================================================================
C  Check whether two subgroups are conjugate in the full group G.
C  H2 = g * H1 * g^{-1} for some g in G
C=================================================================

      INTEGER MAXOPS
      PARAMETER (MAXOPS = 200)

C ---- INPUT ----
      INTEGER NS1, NS2, NG
      DOUBLE PRECISION R_SUB1(3,3,MAXOPS)
      DOUBLE PRECISION R_SUB2(3,3,MAXOPS)
      DOUBLE PRECISION R_GROUP(3,3,MAXOPS)
      DOUBLE PRECISION TOL

C ---- LOCALS ----
      INTEGER I, J, K
      DOUBLE PRECISION G(3,3), GINV(3,3)
      DOUBLE PRECISION TMP1(3,3), TMP2(3,3)
      LOGICAL FOUND
      LOGICAL SAME_OP

C ---- Quick reject ----
      IF (NS1 .NE. NS2) THEN
         CONJUGATE_SUBGROUP = .FALSE.
         RETURN
      END IF

C=================================================================
C  Loop over all g in G
C=================================================================
      DO I = 1, NG

         G(1,1) = R_GROUP(1,1,I)
         G(1,2) = R_GROUP(1,2,I)
         G(1,3) = R_GROUP(1,3,I)
         G(2,1) = R_GROUP(2,1,I)
         G(2,2) = R_GROUP(2,2,I)
         G(2,3) = R_GROUP(2,3,I)
         G(3,1) = R_GROUP(3,1,I)
         G(3,2) = R_GROUP(3,2,I)
         G(3,3) = R_GROUP(3,3,I)

C        Inverse of orthogonal matrix: transpose
         GINV(1,1) = G(1,1)
         GINV(1,2) = G(2,1)
         GINV(1,3) = G(3,1)
         GINV(2,1) = G(1,2)
         GINV(2,2) = G(2,2)
         GINV(2,3) = G(3,2)
         GINV(3,1) = G(1,3)
         GINV(3,2) = G(2,3)
         GINV(3,3) = G(3,3)

C---------------------------------------------------------------
C  Check if g * H1 * g^{-1} = H2
C---------------------------------------------------------------
         FOUND = .TRUE.

         DO J = 1, NS1

C           TMP1 = g * h
            CALL MATMUL3(G, R_SUB1(:,:,J), TMP1)

C           TMP2 = (g*h) * g^{-1}
            CALL MATMUL3(TMP1, GINV, TMP2)

C           check TMP2 in H2
            FOUND = .FALSE.
            DO K = 1, NS2
               IF (SAME_OP(TMP2, R_SUB2(:,:,K), TOL)) THEN
                  FOUND = .TRUE.
                  EXIT
               END IF
            END DO

            IF (.NOT. FOUND) EXIT
         END DO

         IF (FOUND) THEN
            CONJUGATE_SUBGROUP = .TRUE.
            RETURN
         END IF

      END DO

      CONJUGATE_SUBGROUP = .FALSE.
      RETURN
      END

*Deck FILTER_CONJUGATE_SUBGROUPS
      SUBROUTINE FILTER_CONJUGATE_SUBGROUPS
     & (R_SUB_ALL, NSUB_ALL, NFOUND,
     &  R_GROUP, NOPS,
     &  R_SUB_UNIQ, NSUB_UNIQ, NUNIQ,
     &  TOL)
      IMPLICIT NONE

C=================================================================
C  Filter subgroups up to conjugacy.
C  Keeps one representative per conjugacy class.
C  Fortran77 strict implementation.
C=================================================================

      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)

C ---- INPUT ----
      INTEGER NFOUND, NOPS
      INTEGER NSUB_ALL(MAXSUB)
      DOUBLE PRECISION R_SUB_ALL(3,3,MAXSUB,MAXOPS)
      DOUBLE PRECISION R_GROUP(3,3,MAXOPS)
      DOUBLE PRECISION TOL

C ---- OUTPUT ----
      INTEGER NUNIQ
      INTEGER NSUB_UNIQ(MAXSUB)
      DOUBLE PRECISION R_SUB_UNIQ(3,3,MAXSUB,MAXOPS)

C ---- LOCALS ----
      INTEGER I, J, K
      LOGICAL IS_CONJ
      LOGICAL CONJUGATE_SUBGROUP
      DOUBLE PRECISION R_TMP1(3,3,MAXOPS)
      DOUBLE PRECISION R_TMP2(3,3,MAXOPS)

      NUNIQ = 0

C=================================================================
C  Loop over all enumerated subgroups
C=================================================================
      DO I = 1, NFOUND

C        Copy subgroup I from 4D storage into 3D temporary array
         CALL COPY_SUBGROUP_FROM_4D
     &        (R_SUB_ALL, I, R_TMP1, NSUB_ALL(I))

         IS_CONJ = .FALSE.

C        Compare with already accepted unique subgroups
         DO J = 1, NUNIQ

            CALL COPY_SUBGROUP_FROM_4D
     &           (R_SUB_UNIQ, J, R_TMP2, NSUB_UNIQ(J))

            IF (CONJUGATE_SUBGROUP
     &          (R_TMP1, NSUB_ALL(I),
     &           R_TMP2, NSUB_UNIQ(J),
     &           R_GROUP, NOPS,
     &           TOL)) THEN

               IS_CONJ = .TRUE.
               EXIT
            END IF

         END DO

C        New conjugacy class found
         IF (.NOT. IS_CONJ) THEN
            NUNIQ = NUNIQ + 1
            NSUB_UNIQ(NUNIQ) = NSUB_ALL(I)

C           Store representative subgroup explicitly into 4D array
            DO K = 1, NSUB_ALL(I)
               R_SUB_UNIQ(1,1,NUNIQ,K) = R_TMP1(1,1,K)
               R_SUB_UNIQ(1,2,NUNIQ,K) = R_TMP1(1,2,K)
               R_SUB_UNIQ(1,3,NUNIQ,K) = R_TMP1(1,3,K)
               R_SUB_UNIQ(2,1,NUNIQ,K) = R_TMP1(2,1,K)
               R_SUB_UNIQ(2,2,NUNIQ,K) = R_TMP1(2,2,K)
               R_SUB_UNIQ(2,3,NUNIQ,K) = R_TMP1(2,3,K)
               R_SUB_UNIQ(3,1,NUNIQ,K) = R_TMP1(3,1,K)
               R_SUB_UNIQ(3,2,NUNIQ,K) = R_TMP1(3,2,K)
               R_SUB_UNIQ(3,3,NUNIQ,K) = R_TMP1(3,3,K)
            END DO

         END IF

      END DO

      RETURN
      END

*Deck COPY_SUBGROUP_FROM_4D
      SUBROUTINE COPY_SUBGROUP_FROM_4D
     & (RALL, ISUB, ROUT, NSUB)
      IMPLICIT NONE

      INTEGER MAXOPS, MAXSUB
      PARAMETER (MAXOPS = 200, MAXSUB = 200)

      INTEGER ISUB, NSUB
      DOUBLE PRECISION RALL(3,3,MAXSUB,MAXOPS)
      DOUBLE PRECISION ROUT(3,3,MAXOPS)

      INTEGER I

      DO I = 1, NSUB
         ROUT(1,1,I) = RALL(1,1,ISUB,I)
         ROUT(1,2,I) = RALL(1,2,ISUB,I)
         ROUT(1,3,I) = RALL(1,3,ISUB,I)
         ROUT(2,1,I) = RALL(2,1,ISUB,I)
         ROUT(2,2,I) = RALL(2,2,ISUB,I)
         ROUT(2,3,I) = RALL(2,3,ISUB,I)
         ROUT(3,1,I) = RALL(3,1,ISUB,I)
         ROUT(3,2,I) = RALL(3,2,ISUB,I)
         ROUT(3,3,I) = RALL(3,3,ISUB,I)
      END DO

      RETURN
      END

*Deck COPY_SUBGROUP
      SUBROUTINE COPY_SUBGROUP
     & (RIN, ROUT, NSUB)
      IMPLICIT NONE

      INTEGER MAXOPS
      PARAMETER (MAXOPS = 200)

      INTEGER NSUB
      DOUBLE PRECISION RIN(3,3,MAXOPS)
      DOUBLE PRECISION ROUT(3,3,MAXOPS)

      INTEGER I

      DO I = 1, NSUB
         ROUT(1,1,I) = RIN(1,1,I)
         ROUT(1,2,I) = RIN(1,2,I)
         ROUT(1,3,I) = RIN(1,3,I)
         ROUT(2,1,I) = RIN(2,1,I)
         ROUT(2,2,I) = RIN(2,2,I)
         ROUT(2,3,I) = RIN(2,3,I)
         ROUT(3,1,I) = RIN(3,1,I)
         ROUT(3,2,I) = RIN(3,2,I)
         ROUT(3,3,I) = RIN(3,3,I)
      END DO

      RETURN
      END

*Deck CLASSIFY_SUBGROUP
      SUBROUTINE CLASSIFY_SUBGROUP
     & (R_SUB, NSUB,
     &  NAME)
      IMPLICIT NONE

C=================================================================
C  Classify a subgroup given its symmetry operations.
C  Uses the same logic as the main point-group classifier.
C=================================================================

      INTEGER MAXOPS
      PARAMETER (MAXOPS = 200)

C ---- INPUT ----
      INTEGER NSUB
      DOUBLE PRECISION R_SUB(3,3,MAXOPS)

C ---- OUTPUT ----
      CHARACTER*8 NAME

C ---- LOCALS ----
      INTEGER MAX_CN, NC2, NC3, NC4, NC5
      LOGICAL HAS_I, HAS_SIGMA, HAS_SN
      INTEGER FAMILY

C=================================================================
C  Analyze subgroup invariants
C=================================================================
      CALL ANALYZE_SYMMETRY_GROUP
     & (R_SUB, NSUB,
     &  MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN)

C=================================================================
C  Identify family
C=================================================================
      CALL IDENTIFY_GROUP_FAMILY
     & (MAX_CN, NC2, NC3, NC4, NC5,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  FAMILY)

C=================================================================
C  Assign point group name
C=================================================================
      CALL ASSIGN_POINT_GROUP_NAME
     & (FAMILY, MAX_CN,
     &  HAS_I, HAS_SIGMA, HAS_SN,
     &  NAME)

      RETURN
      END
