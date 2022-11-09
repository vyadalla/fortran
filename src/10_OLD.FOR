C
C******************************************************************************
C
C      VERSION FEVRIER 1986                                             TEQ00070
C                                                                       TEQ00080
C  LINKEE POUR LE VAX 3800 LE 10/10/89
C
      INTEGER TITRE(74),ESP(8),PIED(9),KART(80)                         TEQ00110
      INTEGER DICO(215)                                                 TEQ00120
      INTEGER JDIN(14),JNUAN(2)                                         TEQ00130
      INTEGER TEQ0,TEQ1,TEQ2,TEQ3,TEQ4,TEQ5,TEQ6,TEQ7,TEQ8,TEQ9         TEQ00140
      INTEGER TEQA,TEQB,TEQC,TEQD,TEQE,TEQF,TEQG,TEQH,TEQI,TEQJ,TEQK    TEQ00150
      INTEGER TEQW,TEQZ,TEQL,TEQM                                       TEQ00160
      INTEGER ITOUT(11)                                                 TEQ00170
      INTEGER ILABEL(10)                                                TEQ00180
C                                                                       TEQ00190
C                                                                       TEQ00240
      DIMENSION IBUF(1804),Z(500),LABEL(4,200)                          TEQ00250
      DIMENSION LLOAD(32)                                               TEQ00260
      COMMON /POTBIA/ WWWW(1500)                                        TEQ00280
      COMMON /TITRES/ITITR(21,50)                                       TEQ00290
      COMMON /UNIT/IUNIT(20)                                            TEQ00300
      COMMON /RECORD/PLA33(173)                                         TEQ00310
      COMMON /ESPION/ IESP,LLNOM,ILIST,INORD,NCREP                      TEQ00320
      COMMON /DESSIN/ IBUF                                              TEQ00330
      COMMON /OPTCD/ IFEET,ITEST,INUAN(2),IRAIDS                        TEQ00340
      COMMON /MARQ/ PLA21(6000)                                         TEQ00350
      COMMON /BUFFER/ PLA13(2400)                                       TEQ00360
      COMMON /TECK/IC14,IC13,IECHW,INOM,IATA,IBART,IQUAL                TEQ00370
      COMMON /PLACE/ KART,PLA1(9)                                       TEQ00380
      COMMON E,N,X(500),Y(500),IX(4,500),PLA(6000)                      TEQ00390
      COMMON /WORK/ IY(4,500),JBJ(2400)                                 TEQ00400
      COMMON /CANAD/ KX(500),KY(500),PLA20(624) ,ICANAD(200)            TEQ00410
C      DORBYL 6000---->300000                                           TEQ00420
      COMMON /PER/ NPER,PLA2(6000)                                      TEQ00430
      COMMON /IPERM/ NPERM,PLA3(1000)                                   TEQ00440
      COMMON /DIVER/ ISUP,ILANG,IW28,PLA4(6),IDIN,KLAS(4),JPROJ         TEQ00450
      COMMON /SORTI/ PLA6(80)                                           TEQ00460
      COMMON /ENTREE/IPLA7(42)                                          TEQ00470
      COMMON /SOMEL/ PLA8(2)                                            TEQ00480
      COMMON /SUITW/ PLA10(8),MARG                                      TEQ00490
      COMMON /GOUSD/ PLA11(6)                                           TEQ00500
      COMMON /CARTO/ ICAS,XCAR(3,8),TITRE                               TEQ00510
      COMMON /DORBYL/NSTUD,ISTUD(2000)                                  TEQ00520
      COMMON /DOUBL/NDUP,IRED1(2,80),IRED2(2,80),NDUPL(80),IOPD(80)     TEQ00530
      COMMON /RAIDIS/NBRD(400)                                          TEQ00540
      COMMON /REJCOM/LABEL,NLABEL                                       TEQ00550
C                                                                       TEQ00560
C                                                                       TEQ00570
      EQUIVALENCE (PLA13(1101),Z(1))                                    TEQ00580
      EQUIVALENCE (IUNIT(1),INPUT)                                      TEQ00590
      EQUIVALENCE (IUNIT(3),IUNI3)                                      TEQ00600
      EQUIVALENCE (IUNIT(6),IUNI6)                                      TEQ00610
      EQUIVALENCE (IUNIT(7),IUNI7)                                      TEQ00620
      EQUIVALENCE (IUNIT(13),IUNI13)                                    TEQ00630
      EQUIVALENCE (IUNIT(17),IUNI17)                                    TEQ00640
C                                                                       TEQ00650
C                                                                       TEQ00730
C==========================================================I            TEQ00740
C  OVERLAY          6                 7             8      I            TEQ00750
C==========================================================I            TEQ00760
C  BARRES       TEQ1 DIMEN     TEQ3 BARRE     TEQ4 IMETA   I            TEQ00770
C               TEQ5 TABAR                                 I            TEQ00780
C==========================================================I            TEQ00790
C  NOEUDS       TEQ1 DIMEN     TEQ2 TOPOL                  I            TEQ00800
C                                   TOFF                   I            TEQ00810
C==========================================================I            TEQ00820
C  FLECHE       TEQ1 DIMEN                    TEQZ KFLCH   I            TEQ00830
C==========================================================I            TEQ00840
C  BOULON/STUDS TEQW BOUDEZ                                I            TEQ00850
C               DUPLN STUDS                                I            TEQ00860
C==========================================================I            TEQ00870
C  MONTAGE      TEQ1 DIMEN     TEQ9 MONTA                  I            TEQ00880
C  EPURE        TEQ1 DIMEN     TEQ8 EPURE                  I            TEQ00890
C==========================================================I            TEQ00900
C  ARCHI        TEQI ARCHI                                 I            TEQ00910
C==========================================================I            TEQ00920
C  DIVERS                                  TEQ0 AMETA      I            TEQ00930
C                                          TREMY COTAT..   I            TEQ00940
C==========================================================I            TEQ00950
C  DETAIL        TEQ6 NABLA    TEQJ PROCU   TEQ7 PROCT     I            TEQ00960
C==========================================================I            TEQ00970
C  NOMENCLATURE  TEQ6 REWRX    TEQK CBOUL   TEQH MAYEN     I            TEQ00980
C                     TEQTI    TBOUL BULET  TEQE MAYAN     I            TEQ00990
C                                           TEQG BULON     I            TEQ01000
C==========================================================I            TEQ01010
C  GRUGE         TEQB GRUG4    TEQF TROUA   TEQM GRUG3/SED.I            TEQ01020
C  FICHE         TEQC GRUGA    TEQD GRUF/2  TEQL DUPLI/GU  I            TEQ01030
C==========================================================I            TEQ01040
C                                                                       TEQ01050
      DATA TEQ0/'TEQ0'/,TEQ1/'TEQ1'/,TEQ2/'TEQ2'/,TEQ3/'TEQ3'/          TEQ01060
      DATA TEQ4/'TEQ4'/,TEQ5/'TEQ5'/,TEQ6/'TEQ6'/,TEQ7/'TEQ7'/          TEQ01070
      DATA TEQ8/'TEQ8'/,TEQ9/'TEQ9'/,TEQA/'TEQA'/,TEQB/'TEQB'/          TEQ01080
      DATA TEQC/'TEQC'/,TEQD/'TEQD'/,TEQE/'TEQE'/,TEQF/'TEQF'/          TEQ01090
      DATA TEQG/'TEQG'/,TEQH/'TEQH'/,TEQI/'TEQI'/,TEQJ/'TEQJ'/          TEQ01100
      DATA TEQK/'TEQK'/,TEQW/'TEQW'/,TEQZ/'TEQZ'/,TEQL/'TEQL'/          TEQ01110
      DATA TEQM/'TEQM'/                                                 TEQ01120
C                                                                       TEQ01130
C    0 RIEN          1 NOEUDS/POTEAU    2 BARRE/EPURE/MONTAGE           TEQ01140
C    3 DIVERS        4 ARCHI            5 DETAIL                        TEQ01150
C    6 BOULON/GOUJON/OFF                                                TEQ01160
C                                                                       TEQ01170
      DATA LLOAD/1,1,2,8*3,0,0,3,4,2*2,3,5,3,0,3,1,6,6,0,7,4*0,6/       TEQ01180
      DATA ESP/6,'E','S','P','I','O','N',0/                             TEQ01190
      DATA PIED/4,'P','I','E','D',2,'C','M',0/                          TEQ01200
      DATA JDIN/8,'D','I','N',' ','7','9','9','0',3,'E','D','F',0/      TEQ01210
      DATA DICO/5,'N','O','E','U','D',4,'P','O','T','E',3,'B','A','R',4 TEQ01220
     *,'T','R','E','M',3,'P','E','R', 3,'C','O','M',3,'C','O','T',3,'R',TEQ01230
     *'E','P',4,'N','O','T','A',3,'A','X','E',4,'T','R','O','U',3,'F','ITEQ01240
     *','N',5,'N','O','M','E','N',5,'C','O','N','T','O',5,'A','R','C','HTEQ01250
     *','I',5,'E','P','U','R','E',5,'M','O','N','T','A',4,'A','J','O','UTEQ01260
     *',6,'D','E','T','A','I','L',3,'A','R','C',4,'C','L','O','U',      TEQ01270
     *4,'C','A','R','T',6,'F','L','E','C','H','E',8,'*','*','B','O','U',TEQ01280
     *'L','O','N',4,'G','O','U','J',4,'T','E','S','T',5,'*','*','O','F',TEQ01290
     &'F',4,'D','A','T','A',4,'N','O','R','D',                          TEQ01300
     *9,'E','N','T','R','E',' ','A','X','E',6,'N','U','A','N','C','E',  TEQ01310
     *11,'D','U','P','L','I','C','A','T','I','O','N',                   TEQ01320
     *6,'C','L','A','S','S','E',7,'P','R','O','J','E','C','T',          TEQ01330
     *5,'M','A','R','G','E',5,'S','T','A','R','T',5,'R','E','J','E','T',TEQ01340
     *0/                                                                TEQ01350
      DATA ITOUT /4,'T','O','U','T','OU',3,'A','L','L',0/               TEQ01360
      DATA IBL/'    '/                                                  TEQ01370
      DATA JNUAN/'E242','    '/                                         TEQ01380
      DATA ILABEL/'DRAW','N BY',' TEQ','UILA',' SYS','TEM ','ALBI',     TEQ01390
     *' 63.','42.0','7.55'/                                             TEQ01400
C
C      VAX                                                              TEQ01560
C                                                                       TEQ01570
      OPEN(FILE='AKJ7.KIR',UNIT=1,STATUS='OLD')
C      DEFINE FILE 10(810,200,U,IC10)                                    TEQ01580
      OPEN (FILE='CODEBAR.DAT',UNIT=10,RECL=400,FORM='UNFORMATTED',
     * ACCESS='DIRECT',STATUS='OLD')
C      DEFINE FILE 11(10000,346,U,IUUU)                                  TEQ01590
C      DEFINE FILE 11(30000,346,U,IUUU)                                  TEQ01590
      OPEN(UNIT=11,FILE='ZNOMEN.DAT',ACCESS='DIRECT',
     * RECL=692,STATUS='OLD',FORM='UNFORMATTED')
C      DEFINE FILE 12(400,1584,U,IC12)                                   TEQ01600
      OPEN(UNIT=12,FILE='MAT.SOM',ACCESS='DIRECT',
     * RECL=3168,STATUS='OLD',FORM='UNFORMATTED')
C      DEFINE FILE 13(4000,52,U,IC13)                                    TEQ01610
      OPEN(UNIT=13,FILE='ATA.SOM',ACCESS='DIRECT',
     * RECL=104,STATUS='OLD',FORM='UNFORMATTED')
C      DEFINE FILE 15(1100,1600,U,IC15)                                  TEQ01620
      OPEN(UNIT=15,FILE='PRO.SOM',ACCESS='DIRECT',
     * STATUS='OLD',FORM='UNFORMATTED',RECL=3200)
C      DEFINE FILE 16(2000,40,U,IC16)                                    TEQ01630
      OPEN (FILE='FER.SOM',UNIT=16,RECL=80,FORM='UNFORMATTED',
     * ACCESS='DIRECT',STATUS='OLD')
      DEFINE FILE 17(600,1240,U,IC17)                                   TEQ01640
C      OPEN (FILE='RONDELLE.DAT',UNIT=17,RECL=2480,FORM='UNFORMATTED',
C     * ACCESS='DIRECT',STATUS='OLD')
                                                                        TEQ01650
 2000 FORMAT(80A1)                                                      TEQ02300
 3002 FORMAT(10X,80A1)                                                  TEQ02310
 3000 FORMAT('1C8:  NOMBRE D UNITES TEQUILA ',I4,' NOEUDS',I4,' BARRES')TEQ02320
 3001 FORMAT(' C9:',2F10.2,F10.6,2F10.2,'HORS DESSIN',F10.2,' TOTAL ',  TEQ02330
     *F10.2)                                                            TEQ02340

C      CALL DAREAD(IBUF,80,200,15)
C      write(*,*) IBUF(1:80)
C      write(*,*) "PEDRO"
C                                                                       TEQ02350
      DO 9999 I=1,20                                                    TEQ02360
      IUNIT(I)=I                                                        TEQ02370
 9999 CONTINUE                                                          TEQ02380
      IESP=0                                                            TEQ02430
      IQUAL=0                                                           TEQ02440
      IARCH=0                                                           TEQ02450
      ICAS=0                                                            TEQ02460
      ISUP=0                                                            TEQ02470
      INORD=0                                                           TEQ02480
      IDIN=0                                                            TEQ02490
      INTAX=0                                                           TEQ02500
      IBART=2                                                           TEQ02510
      IBUF(1)=3                                                         TEQ02520
      ILANG=0                                                           TEQ02530
      ITEST=0                                                           TEQ02540
      MPER=0                                                            TEQ02550
      IFEET=0                                                           TEQ02560
      NSTUD=0                                                           TEQ02570
      NPERM=0                                                           TEQ02580
      IATA=1                                                            TEQ02590
      JPROJ=0                                                           TEQ02600
      KLAS(1)=0                                                         TEQ02610
      MARG=0                                                            TEQ02620
      NLABEL=0                                                          TEQ02630
      INUAN(1)=JNUAN(1)                                                 TEQ02640
      INUAN(2)=JNUAN(2)                                                 TEQ02650
C                                                                       TEQ02660
      ILIST=0                                                           TEQ02670
      NDUP=0                                                            TEQ02680
      CALL TAB(NBRD,1,400,0,0)                                          TEQ02690
      IRAIDS=1                                                          TEQ02700
      CALL OVLY(TEQ1,IRC)                                               TEQ02710
      CALL OVLY(TEQ2,IRC)                                               TEQ02720
      CALL OVLY(TEQZ,IRC)                                               TEQ02730
C                                                                       TEQ02740
C                                                                       TEQ02830
      CALL IBENA(IBUF,270,IUNI7)                                        TEQ02840
      CALL TAB(IBUF,301,1804,0,0)                                       TEQ02850
      PO=60.                                                            TEQ02860
      PA=0.                                                             TEQ02870
      P2=0.0005                                                         TEQ02880
      M=0                                                               TEQ02890
      INOM=0                                                            TEQ02900
      CALL PNUMA(0.,0.,M,-0.5,-0.5)                                     TEQ02910
      CALL PCARA(-.1,1.8,0,ILABEL,                                      TEQ02920
     *40,0.2,0.3,0.,1.)                                                 TEQ02930
      IU17=-IUNI17                                                      TEQ02940
      CALL DAWRIT(IBUF,2480,IC17,IU17)                                  TEQ02950
C      write(*,*) INPUT, IUNIT(:)
      CALL CARTOU                                                       TEQ02960
CDH
      DO I=1,7
C         CALL DROOG2
      END DO
      CALL TAB(KX,1,500,0.,0)                                           TEQ02970
      N=0                                                               TEQ02980
      KLOAD=1                                                           TEQ02990
    1 ILOAD=KLOAD                                                       TEQ03000
      READ(INPUT,2000,END=115)KART                                      TEQ03010
      WRITE(IUNI3,3002)KART                                             TEQ03030
      KI=1                                                              TEQ03040
      KF=80                                                             TEQ03050
      CALL TEXTZ(KART,KI,80,ESP,JOB)                                    TEQ03060
      IF(JOB    .EQ.0)GO TO  777                                        TEQ03070
      IF(JOB.EQ.1)IESP=1                                                TEQ03080
  777 CONTINUE                                                          TEQ03090
      CALL TEXTZ(KART,KI,80,PIED,JOB)                                   TEQ03100
      IF(JOB    .EQ.0)GO TO  1                                          TEQ03110
      IF(JOB.LT.3)IFEET=JOB                                             TEQ03120
      CALL TEXTZ(KART,KI,80,JDIN,KOD0)                                  TEQ03130
      IF(KOD0.EQ.1)IDIN=1                                               TEQ03140
      IF(KOD0.EQ.2)IDIN=2                                               TEQ03150
      CALL TEXTZ(KART,KI,80,DICO,JOB)                                   TEQ03160
      IF(JOB    .EQ.0)GO TO  1                                          TEQ03170
C                                                                       TEQ03180
C        1 NOEUD   2 POTEAU  3 BARRE                                    TEQ03190
C        4 TREMY   5 PERCA   6 COMMENTAIRES                             TEQ03200
C        7 COTAT   8 REPERAGE9 NOTA                                     TEQ03210
C       10 AXE    11 TROU   12 FIN                                      TEQ03220
C       13 NOMEN  14 CONTOUR15 ARCHI                                    TEQ03230
C       16 EPURE  17 MONTA  18 AJOUTER                                  TEQ03240
C       19 DETAIL 20 ARC    21 CLOU                                     TEQ03250
C       22 CARTOU 23 FLECHE 24 **BOULON                                 TEQ03260
C       25 GOUJON 26 TEST   27 **OFF                                    TEQ03270
C       28 LISTER 29 NORD  30 ENTRE AXE                                 TEQ03280
C       31 NUANCE 32 DUPLICATION                                        TEQ03290
C       33 CLASSE 34 PROJET 35 MARGE                                    TEQ03300
C       36 START  37 REJET                                              TEQ03310
C                                                                       TEQ03320
C      NOMBRE MAXI DE MOTS CLEFS                                        TEQ03330
C                                                                       TEQ03340
      IF(JOB.EQ.35)MARG=1                                               TEQ03350
      IF(JOB.GT.37)GO TO 1                                              TEQ03360
      KLOAD=LLOAD(JOB)                                                  TEQ03370
      IF(JOB.LE.3)GO TO(2,14,3),JOB                                     TEQ03380
      IF(JOB.EQ.12)GO TO 15                                             TEQ03390
      IF(JOB.EQ.13)GO TO 16                                             TEQ03400
      IF(JOB.EQ.15)GO TO 22                                             TEQ03410
      IF(JOB.EQ.16)GO TO 18                                             TEQ03420
      IF(JOB.EQ.17)GO TO 20                                             TEQ03430
      IF(JOB.EQ.19)GO TO 109                                            TEQ03440
      IF(JOB.EQ.21)GO TO 51                                             TEQ03450
      IF(JOB.GE.23)GO TO 100                                            TEQ03460
      CALL OVLY(TEQ0,IRC)                                               TEQ03470
      CALL OVLY(TEQJ,IRC)                                               TEQ03480
  100 CONTINUE                                                          TEQ03490
      GO TO(2,14,3,4,5,6,7,8,9,12,13,15,16,17,22,18,19,21,109,19,51,52, TEQ03500
     *53,54,55,56,57,58,59,60,61,62,63,64,1,66,37,1),JOB                TEQ03510
   66 CALL LYS(KART,KI,80,R,IRAIDS,KOD)                                 TEQ03520
      GO TO 1                                                           TEQ03530
C                                                                       TEQ03540
C      DUPLICATION                                                      TEQ03550
C                                                                       TEQ03560
   62 CALL OVLY(TEQW,IRC)                                               TEQ03570
      CALL DUPLC(0)                                                     TEQ03580
       GO TO 1                                                          TEQ03590
   61 IF(KART(KI).NE.IBL)GO TO 161                                      TEQ03600
      KI=KI+1                                                           TEQ03610
      GO TO 61                                                          TEQ03620
  161 CALL PACKN(INUAN,4,KART(KI),1,8)                                  TEQ03630
      GO TO 1                                                           TEQ03640
   60 INTAX=1                                                           TEQ03650
      GO TO 1                                                           TEQ03660
   59 CALL LYS(KART,KI,80,RR,INORD,KOD)                                 TEQ03670
      IF(INORD.LT.0.OR.INORD.GT.4)INORD=0                               TEQ03680
      GO TO 1                                                           TEQ03690
   58 ILIST=1                                                           TEQ03700
      CALL TEXTZ(KART,KI,80,ITOUT,KOD)                                  TEQ03710
      IF(KOD.EQ.1)ILIST=2                                               TEQ03720
      GO TO 1                                                           TEQ03730
   63 DO 635 I=1,4                                                      TEQ03740
      KLAS(I)=0                                                         TEQ03750
      CALL LYS(KART,KI,80,RR,KLAS(I),KOD)                               TEQ03760
      IF(KOD.EQ.3)KLAS(I)=RR                                            TEQ03770
  635 CONTINUE                                                          TEQ03780
      GO TO 1                                                           TEQ03790
   64 JPROJ=1                                                           TEQ03800
      GO TO 1                                                           TEQ03810
C                                                                       TEQ03820
C      NOEUD                                                            TEQ03830
C                                                                       TEQ03840
    2 CONTINUE                                                          TEQ03850
      IF(KLOAD.EQ.0)KLOAD=ILOAD                                         TEQ03860
      IF(ILOAD.NE.1)CALL OVLY(TEQ2,IRC)                                 TEQ03870
      CALL OVLY(TEQZ,IRC)                                               TEQ03880
      KENTRY=1                                                          TEQ03890
      CALL TOPOL (KENTRY)                                               TEQ03900
      GO TO 1                                                           TEQ03910
C                                                                       TEQ03920
C      POTEAU                                                           TEQ03930
C                                                                       TEQ03940
   14 CONTINUE                                                          TEQ03950
      IF(ILOAD.GT.2)CALL OVLY(TEQ1,IRC)                                 TEQ03960
      IF(ILOAD.NE.1)CALL OVLY(TEQ2,IRC)                                 TEQ03970
      CALL OVLY(TEQZ,IRC)                                               TEQ03980
      KENTRY=2                                                          TEQ03990
      CALL TOPOL (KENTRY)                                               TEQ04000
      GO TO 1                                                           TEQ04010
C                                                                       TEQ04020
C      FLECHE                                                           TEQ04030
C                                                                       TEQ04040
   53 CONTINUE                                                          TEQ04050
      IF(ILOAD.GT.2)CALL OVLY(TEQ1,IRC)                                 TEQ04060
      IF(ILOAD.NE.1)CALL OVLY(TEQ2,IRC)                                 TEQ04070
      CALL OVLY(TEQZ,IRC)                                               TEQ04080
      CALL KFLCH(0,N,X,Y,Z)                                             TEQ04090
      GO TO 1                                                           TEQ04100
C                                                                       TEQ04110
C      **OFFSET (COMPATIBILITE KIR)                                     TEQ04120
C                                                                       TEQ04130
   57 CONTINUE                                                          TEQ04140
      IF(ILOAD.NE.1)CALL OVLY(TEQ2,IRC)                                 TEQ04150
      CALL OVLY(TEQZ,IRC)                                               TEQ04160
      CALL TOFF                                                         TEQ04170
      GO TO 1                                                           TEQ04180
C                                                                       TEQ04190
C      BARRE                                                            TEQ04200
C                                                                       TEQ04210
    3 CONTINUE                                                          TEQ04220
      IF(ILOAD.GT.2)CALL OVLY(TEQ1,IRC)                                 TEQ04230
      CALL OVLY(TEQ3,IRC)                                               TEQ04240
      CALL OVLY(TEQ4,IRC)                                               TEQ04250
      CALL BARRE                                                        TEQ04260
      P1=0.38                                                           TEQ04270
      P3=1.                                                             TEQ04280
      KP=2                                                              TEQ04290
      KMETR=1                                                           TEQ04300
      CALL SISSY(1,IUNI13)                                              TEQ04310
      CALL OVLY(TEQ5,IRC)                                               TEQ04320
      KLOAD=8                                                           TEQ04330
      CALL TABAR                                                        TEQ04340
  555 CONTINUE                                                          TEQ04350
      GO TO 1                                                           TEQ04360
C                                                                       TEQ04370
C      EPURE                                                            TEQ04380
C                                                                       TEQ04390
   18 CONTINUE                                                          TEQ04400
      IF(ILOAD.GT.2)CALL OVLY(TEQ1,IRC)                                 TEQ04410
      CALL OVLY(TEQ8,IRC)                                               TEQ04420
      CALL EPURE                                                        TEQ04430
      CALL SISSY(1,IUNI13)                                              TEQ04440
      P1=0.38                                                           TEQ04450
      P3=1.                                                             TEQ04460
  107 KMETR=2                                                           TEQ04470
      GO TO 555                                                         TEQ04480
C                                                                       TEQ04490
C      MONTAGE                                                          TEQ04500
C                                                                       TEQ04510
   20 CONTINUE                                                          TEQ04520
      IF(ILOAD.GT.2)CALL OVLY(TEQ1,IRC)                                 TEQ04530
      CALL OVLY(TEQ9,IRC)                                               TEQ04540
      CALL MONTA                                                        TEQ04550
      P1=0.22                                                           TEQ04560
      CALL SISSY(1,IUNI13)                                              TEQ04570
      KMETR=2                                                           TEQ04580
      P2=0.                                                             TEQ04590
      P3=0.65                                                           TEQ04600
      GO TO 107                                                         TEQ04610
C                                                                       TEQ04620
C      DETAIL (PUNCH)                                                   TEQ04630
C                                                                       TEQ04640
  109 CONTINUE                                                          TEQ04650
C                                                                       TEQ04660
C      VERIFICATION PUBLIC PRIVE                                        TEQ04670
C                                                                       TEQ04680
C                                                                       TEQ04740
      CALL OVLY(TEQ6,IRC)                                               TEQ04750
      CALL OVLY(TEQJ,IRC)                                               TEQ04760
      CALL OVLY(TEQ7,IRC)                                               TEQ04770
      CALL DTATQ                                                        TEQ04780
      GO TO 1                                                           TEQ04800
C                                                                       TEQ04810
C      DIVERS .......                                                   TEQ04820
C                                                                       TEQ04830
    4 KENTRY=1                                                          TEQ04840
      CALL TREMY (KENTRY)                                               TEQ04850
      PO=PO+15.                                                         TEQ04860
      GO TO 1                                                           TEQ04870
    5 CALL PERCA                                                        TEQ04880
      PO=PO+3.                                                          TEQ04890
      GO TO 1                                                           TEQ04900
    6 CALL TEXTZ(KART,KI,80,DICO,JOB)                                   TEQ04910
      IF(JOB    .EQ.0)GO TO  1                                          TEQ04920
      IF(JOB.GT.3)GO TO 7666                                            TEQ04930
      GO TO(10,6,11),JOB                                                TEQ04940
 7666 CONTINUE                                                          TEQ04950
      KI=KI+1                                                           TEQ04960
      GO TO 6                                                           TEQ04970
   10 CALL CONOE                                                        TEQ04980
      PO=PO+3.                                                          TEQ04990
      GO TO 1                                                           TEQ05000
   11 CALL COBAR                                                        TEQ05010
      PO=PO+3.                                                          TEQ05020
      GO TO 1                                                           TEQ05030
    7 CALL COTAT                                                        TEQ05040
      PO=PO+N*0.05*400.*E                                               TEQ05050
      GO TO 1                                                           TEQ05060
    8 CONTINUE                                                          TEQ05070
C                                                                       TEQ05080
C      ELIMINER MODULE REPPORT EVENTUEL                                 TEQ05090
C                                                                       TEQ05100
      IF(KART(KI).EQ.DICO(3))GO TO 1                                    TEQ05110
      CALL FILES                                                        TEQ05120
      PO=PO+N*0.02*400.*E                                               TEQ05130
      GO TO 1                                                           TEQ05140
    9 CALL NOTA                                                         TEQ05150
      PO=PO+3.                                                          TEQ05160
      GO TO 1                                                           TEQ05170
   12 CALL AXES(KART,KI,KF)                                             TEQ05180
      PO=PO+10.                                                         TEQ05190
      GO TO 1                                                           TEQ05200
   13 KENTRY=2                                                          TEQ05210
      CALL TREMY (KENTRY)                                               TEQ05220
      PO=PO+10.                                                         TEQ05230
      GO TO 1                                                           TEQ05240
C                                                                       TEQ05250
C      CLOU                                                             TEQ05260
C                                                                       TEQ05270
   51 CALL LYS(KART,KI,80,R,IBART,KOD)                                  TEQ05280
      KI=KI+1                                                           TEQ05290
      IF(KOD.EQ.4)GO TO 51                                              TEQ05300
      IPH5=0                                                            TEQ05310
      KI=KI-1                                                           TEQ05320
  151 CALL LYS(KART,KI,80,R,IPH5,KOD)                                   TEQ05330
      IF(KOD.NE.2)GO TO 1                                               TEQ05340
      IBART=IBART+10*IPH5                                               TEQ05350
      CALL LYS(KART,KI,80,R,IQUAL,KOD)                                  TEQ05360
      IF(KOD.NE.2)IQUAL=0                                               TEQ05370
      GO TO 1                                                           TEQ05380
C                                                                       TEQ05390
C      CARTOUCHE SPECIAL                                                TEQ05400
C                                                                       TEQ05410
   52 CALL CARLC(0)                                                     TEQ05420
      GO TO 1                                                           TEQ05430
C                                                                       TEQ05440
C      NOMENCLATURE                                                     TEQ05450
C                                                                       TEQ05460
   16 INOM=1                                                            TEQ05470
   50 CALL LYS(KART,KI,80,R,INOM,KOD)                                   TEQ05480
      KI=KI+1                                                           TEQ05490
      IF(KOD-2)1,1,50                                                   TEQ05500
   17 CALL BMETA                                                        TEQ05510
      PO=PO+30.                                                         TEQ05520
      GO TO 1                                                           TEQ05530
   19 CALL AMETA                                                        TEQ05540
      PO=PO+30.                                                         TEQ05550
      GO TO 1                                                           TEQ05560
   21 CONTINUE                                                          TEQ05570
      NPER=MPER                                                         TEQ05580
      CALL PERRA                                                        TEQ05590
      MPER=NPER                                                         TEQ05600
      GO TO 1                                                           TEQ05610
   22 IARCH=1                                                           TEQ05620
      GO TO 1                                                           TEQ05630
C                                                                       TEQ05640
C     **BOULON   (COMPATIBILITE KIR)                                    TEQ05650
C                                                                       TEQ05660
   54 CONTINUE                                                          TEQ05670
      CALL OVLY(TEQW,IRC)                                               TEQ05680
      CALL BOUDES(KART,KI)                                              TEQ05690
      GOTO 1                                                            TEQ05700
C                                                                       TEQ05710
C      GOUJON                                                           TEQ05720
C                                                                       TEQ05730
   55 CONTINUE                                                          TEQ05740
      CALL OVLY(TEQW,IRC)                                               TEQ05750
      CALL STUDS(NSTUD,ISTUD)                                           TEQ05760
      GO TO 1                                                           TEQ05770
   56 ITEST=1                                                           TEQ05780
      GO TO 1                                                           TEQ05790
C                                                                       TEQ05800
C      REJET                                                            TEQ05810
C                                                                       TEQ05820
   37  CALL REJET                                                       TEQ05830
       write(*,*) "PEDRO: REJET fini"
       GO TO 1                                                          TEQ05840
C                                                                       TEQ05850
C                                                                       TEQ05860
  115 IATA=2                                                            TEQ05870
C                                                                       TEQ05880
C      FIN DU DESSIN,DEBUT DU TRAITEMENT DES FICHES ET DES NOMENCLATURESTEQ05890
C                                                                       TEQ05900
   15 CALL POSA(U,V)                                                    TEQ05910
      IN=N                                                              TEQ05920
      IP=IC14-6                                                         TEQ05930
      CALL OVLY(TEQ0,IRC)                                               TEQ05940
      CALL OVLY(TEQ2,IRC)                                               TEQ05950
      CALL TPOS                                                         TEQ05960
      N99=999                                                           TEQ05970
      CALL PNUMA(0.,0.,N99,0.,0.)                                       TEQ05980
C                                                                       TEQ05990
C      FIN SI PASSAGE DE TEST SANS NOMENCLATURE OU FICHE                TEQ06000
C                                                                       TEQ06010
C                                                                       TEQ06060
      IF(ITEST.EQ.1)STOP                                                TEQ06070
      CALL OVLY(TEQI,IRC)                                               TEQ06080
      IF(IARCH.EQ.1)CALL ARCHI                                          TEQ06100
      INOM1=INOM                                                        TEQ06120
      IF(INOM.EQ.0)INOM=1                                               TEQ06130
      IF(INOM)35,35,34                                                  TEQ06140
   34 CONTINUE                                                          TEQ06150
      IF(KMETR.EQ.2)GO TO 134                                           TEQ06160
      N=IN                                                              TEQ06170
      CALL OVLY(TEQB,IRC)                                               TEQ06180
      CALL OVLY(TEQF,IRC)                                               TEQ06190
      CALL OVLY(TEQM,IRC)                                               TEQ06200
      NPER=MPER                                                         TEQ06210
      ITEST=INTAX                                                       TEQ06220
      CALL GRUGE                                                        TEQ06230
      CALL OVLY(TEQC,IRC)                                               TEQ06240
      CALL OVLY(TEQD,IRC)                                               TEQ06250
      CALL OVLY(TEQL,IRC)                                               TEQ06260
      KARM=0                                                            TEQ06270
      CALL GRUF(KARM)                                                   TEQ06280
      PA=1.0+0.7*IW28                                                   TEQ06290
      N=IN                                                              TEQ06300
      IECHW=2-KMETR                                                     TEQ06310
      INOM=INOM1                                                        TEQ06320
  134 CONTINUE                                                          TEQ06330
      IF(INOM1.EQ.0)GOTO 33                                             TEQ06340
      CALL OVLY(TEQA,IRC)                                               TEQ06350
      CALL OVLY(TEQK,IRC)                                               TEQ06360
      CALL OVLY(TEQE,IRC)                                               TEQ06370
C                                                                       TEQ06400
      IBUF(1)=0                                                         TEQ06410
      CALL POUT(IBUF)                                                   TEQ06420
      CALL MAYAN                                                        TEQ06430
      PA=PA+0.05                                                        TEQ06440
      IF(KMETR.EQ.2)GO TO 33                                            TEQ06450
      IF((NPER+NPERM).GT.0)GO TO 133                                    TEQ06460
      IF(IC13)33,33,133                                                 TEQ06470
  133 CONTINUE                                                          TEQ06480
      CALL OVLY(TEQG,IRC)                                               TEQ06490
      CALL BULON                                                        TEQ06500
      CALL OVLY(TEQH,IRC)                                               TEQ06510
      PA=PA+0.25                                                        TEQ06520
      CALL MAYEN                                                        TEQ06530
      GOTO 33                                                           TEQ06550
 35   IF(IW28.NE.0)GOTO 34                                              TEQ06560
 33   N=IN                                                              TEQ06570
      STOP                                                              TEQ06680
      END                                                               TEQ06690
      SUBROUTINE TPOS                                                   TEQ06700
C                                                                       TEQ06710
C      VERSION FEVRIER 1986                                             TEQ06720
C        TEQ.FT0                                                        TEQ06730
C                                                                       TEQ06740
      COMMON /SUITW/ XW,YW,XNIV,INIVO,IOP,IQT,XBAR,ICADR                TEQ06750
      COMMON /DIVER/ ISUP,ILANG,IW28,W1,IKOT(5)                         TEQ06760
      COMMON /RAIDIS/NBRE(400)                                                  
      COMMON /OPTCD/IDUM,ITEST,VBD(2),IRAIDS                            TEQ06780
      COMMON /CARTO/ICAS,XCAR(3,8),ITTRR(74)                            TEQ06790
      COMMON /PLACE/IPPAG,ICLIEN(56),IDAT(2),IAFF(2),NUMRO(7)           TEQ06800
C                                                                       TEQ06810
      DIM=80.+2.*IKOT(4)                                                TEQ06820
      NCASE=3                                                           TEQ06830
      IF(DIM.LT.90.)NCASE=2                                             TEQ06840
      XECR=138.+21.*ICADR-XW                                            TEQ06850
      CALL CARTC(0,XECR-20.,-YW)                                        TEQ06860
      CALL PLUME(0)                                                     TEQ06870
      NB=99                                                             TEQ06880
      IDXX=0                                                            TEQ06890
C     IF(NRAID.GT.0)NB=NB-1                                             TEQ06900
      CALL PNUMA(XECR,-YW,NB,0.,0.)                                     TEQ06910
C     WRITE(6,7786)ITEST                                                        
C7786 FORMAT(' ITEST=',I12)                                                     
      IF(ITEST.EQ.1)GO TO 1                                             TEQ06930
      CALL PACKN(IAFF,4,ITTRR(1),2,4)                                   TEQ06940
      CALL PACKN(IDAT,4,ITTRR(64),2,4)                                  TEQ06950
      DO 2 I=1,56                                                       TEQ06960
    2 ICLIEN(I)=ITTRR(I+4)                                              TEQ06970
      DO 3 I=1,7                                                        TEQ06980
    3 NUMRO(I)=ITTRR(I+67)                                              TEQ06990
      CALL RAIDF(-10.,0.,NCASE,NRAID)                                         TE
      IDXX=(NRAID+NCASE-1)/NCASE                                        TEQ07010
    1 CONTINUE                                                          TEQ07020
      DXX=IDXX*21.+IDXX/5*1.                                            TEQ07030
      CALL PNUMA(DXX,0.,NB,0.,0.)                                       TEQ07040
      RETURN                                                            TEQ07050
      END                                                               TEQ07060
      SUBROUTINE REJET                                                  TEQ07070
C                                                                       TEQ07080
C                                                                       TEQ07090
C     COMMANDE REJET                                                    TEQ07100
C                                                                       TEQ07110
C      VERSION JUIN 1984 PARIS                                          TEQ07120
C                                                                       TEQ07130
      INTEGER DICO(25),KART(80),IREP(2),LABEL(4,200),OUT                TEQ07140
      DIMENSION RLABEL(4,200)                                           TEQ07150
      COMMON /REJCOM/LABEL,NLABEL                                       TEQ07160
      EQUIVALENCE (LABEL(1,1),RLABEL(1,1))                              TEQ07170
      COMMON/UNIT/INPUT,IDUMMY,OUT                                      TEQ07180
      DATA DICO /5,'R','E','J','E','T',3,'F','I','N',                   TEQ07190
     *2,'D','X',2,'D','Y',7,'P','A','R','T','I','E','L',0/              TEQ07200
      DATA IBL/'    '/,ICOT/''''/,IEGAL/'='/                            TEQ07210
      INUM=0                                                            TEQ07220
C                                                                       TEQ07230
C    DECODAGE CARTE REJET                                               TEQ07240
C                                                                       TEQ07250
      INUM=0                                                            TEQ07260
      KNUM=1                                                            TEQ07270
  20  READ(INPUT,2000)KART                                              TEQ07280
      IREJET=0                                                          TEQ07290
      DX=9999.                                                          TEQ07300
      DY=9999.                                                          TEQ07310
 2000 FORMAT(80A1)                                                      TEQ07320
      LDEB=1                                                            TEQ07330
  40  CALL TEXTZ(KART,LDEB,80,DICO,KOD)                                 TEQ07340
      KI=LDEB                                                           TEQ07350
      IF(KOD.EQ.2.AND.IREJET.EQ.0)GO TO 500                             TEQ07360
      IF(KOD.EQ.0)GO TO 20                                              TEQ07370
      IF(KOD.LE.5)GO TO 10                                              TEQ07380
      DO 1 I=LDEB,80                                                    TEQ07390
      IF(KART(I)-ICOT)1,2,1                                             TEQ07400
    1 CONTINUE                                                          TEQ07410
    2 I=I+1                                                             TEQ07420
      IF(I.GT.80)GO TO 20                                               TEQ07430
      DO 3 J=I,80                                                       TEQ07440
      IF(KART(J)-ICOT)3,4,3                                             TEQ07450
    3 CONTINUE                                                          TEQ07460
      RETURN                                                            TEQ07470
    4 DO 6 L=1,2                                                        TEQ07480
    6 IREP(L)=IBL                                                       TEQ07490
      KK=J-I                                                            TEQ07500
      CALL PACKN(IREP,4,KART(I),1,KK)                                   TEQ07510
      INUM=INUM+1                                                       TEQ07520
      IF(INUM.GT.200)GO TO 333                                          TEQ07530
      LABEL(1,INUM)=IREP(1)                                             TEQ07540
      LABEL(2,INUM)=IREP(2)                                             TEQ07550
C                                                                       TEQ07560
C                                                                       TEQ07570
      LDEB=J+1                                                          TEQ07580
      IF(LDEB.GE.80)GO TO 20                                            TEQ07590
 700  DX=9999.                                                          TEQ07600
      DY=9999.                                                          TEQ07610
      GO TO 40                                                          TEQ07620
C                                                                       TEQ07630
C                                                                       TEQ07640
C                                                                       TEQ07650
   10 GO TO (100,200,300,400,110),KOD                                   TEQ07660
100   CONTINUE                                                          TEQ07670
      DX=9999.                                                          TEQ07680
      DY=9999.                                                          TEQ07690
      CALL TEXTZ(KART,KI,80,DICO,KOD)                                   TEQ07700
      IF(KOD.NE.5)GO TO 130                                             TEQ07710
 110  DX=-9999.                                                         TEQ07720
      DY=-9999.                                                         TEQ07730
 130  IREJET=1                                                          TEQ07740
      LDEB=KI                                                           TEQ07750
      GO TO 500                                                         TEQ07760
200   CONTINUE                                                          TEQ07770
      NLABEL=INUM                                                       TEQ07780
      RETURN                                                            TEQ07790
300   CONTINUE                                                          TEQ07800
340   CALL LYS(KART,KI,80,R,I,KODE)                                     TEQ07810
      IF(KODE.EQ.1)GO TO 700                                            TEQ07820
      IF(KODE.EQ.4.AND.KART(KI).EQ.IEGAL)GO TO 320                      TEQ07830
      DX=R                                                              TEQ07840
      IF(KODE.EQ.2)DX=I                                                 TEQ07850
      LDEB=KI                                                           TEQ07860
      GO TO 40                                                          TEQ07870
320   KI=KI+1                                                           TEQ07880
      GO TO 340                                                         TEQ07890
C                                                                       TEQ07900
400   CONTINUE                                                          TEQ07910
440   CALL LYS(KART,KI,80,R,I,KODE)                                     TEQ07920
      IF(KODE.EQ.4.AND.KART(KI).EQ.IEGAL)GO TO 420                      TEQ07930
      DY=R                                                              TEQ07940
      IF(KODE.EQ.2)DY=I                                                 TEQ07950
      LDEB=KI                                                           TEQ07960
      GO TO 500                                                         TEQ07970
420   KI=KI+1                                                           TEQ07980
      GO TO 440                                                         TEQ07990
C                                                                       TEQ08000
C   MISE A JOUR TABLEAU LABEL                                           TEQ08010
C                                                                       TEQ08020
500    DO 501 LN=KNUM,INUM                                              TEQ08030
       RLABEL(3,LN)=DX                                                  TEQ08040
       RLABEL(4,LN)=DY                                                  TEQ08050
501    CONTINUE                                                         TEQ08060
       KNUM=INUM+1                                                      TEQ08070
       IF(KOD.EQ.2)GO TO 200                                            TEQ08080
       GO TO 40                                                         TEQ08090
 333   WRITE(OUT,633)                                                   TEQ08100
 633   FORMAT('  E : PLUS DE 200 REPERES ')                             TEQ08110
       WRITE(2,634)                                                     TEQ08100
 634   FORMAT('  E : PLUS DE 200 REPERES ')                             TEQ08110
       RETURN                                                           TEQ08120
      END                                                               TEQ08130
