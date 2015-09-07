C***********************************************************************
C*   RESCAL:                                                   *
C*      PROGRAM  TO CALCULATE  COOPER-NATHANS  RESOLUTION      *
C*      MATRIX AND/OR RIS0  METHOD  CONTRIBUTIONS TO PHONON    *
C*      WIDTH. THE WIDTHS  OF BRAGG  SCANS ARE  CALCULATED,    *
C*      AND THE VANADIUM  WIDTH, FOR A TRIPLE-AXIS SPECTRO-    *
C*      METER. THE PARAMETERS ARE ENTERED INTERACTIVELY.       *
C*                                                             *
C*    WRITTEN BY:   M. HARGREAVE & P. HULLAH,                  *
C*                    DEPT. OF PHISICS,                        *
C*                      P. C. L.,                              *
C*                        115, NEW CAVENDISH ST.,              *
C*                          LONDON. W1M 8JS.                   *
C*                             (TEL 01-486 5811 EXT 6384)      *
C*                                                JULY 1979    *
C*                                       UPDATED APRIL 1980    *
C*  Modernised code (GFORTRAN) by E. FARHI , ILL/CS (2015)     *
C*  Version $Date$
C***************************************************************
C
C Compile with: ifort    -o rescal rescal.for
C               gfortran -o rescal rescal.for
C               g77      -o rescal rescal.for
C
C List of routines: PROGRAM RESCAL
C
C Computational routines:
C	TRNVCT
C	BRAG(AGM)
C	PHON(AGM)
C	RESOL(AGM)
C	RECLAT
C	AFILL
C	VIVF(VI,VF)		called by RESOL
C	CONVRT(A,AZ)	called by RESOL
C	DIAG(A,ADA,B)	called by RESOL
C	TRANS			called by TRNVCT
C	RESCON		called by PHON

	BLOCK DATA
C
	CHARACTER*8 NAM(42)
	COMMON/CHNM/NAM
	DATA NAM/
     1  'DM','DA','ETAM','ETAA','ETAS','SM','SS','SA','KFIX','FX',
     2	'ALF1','ALF2','ALF3','ALF4','BET1','BET2','BET3','BET4',
     3	'AS','BS','CS','AA','BB','CC','AX','AY','AZ','BX','BY','BZ',
     4	'QH','QK','QL','EN','DH','DK','DL','DE','GH','GK',
     5	'GL','GMOD' /
	DOUBLE PRECISION DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1,ALF2,
     1	ALF3,ALF4,BET1,BET2,BET3,BET4,AS,BS,CS,AA,BB,CC,AX,AY,AZ,BX,
     2	BY,BZ,QH,QK,QL,EN,DH,DK,DL,DE,GH,GK,GL,GMOD
	COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1,ALF2,
     1	ALF3,ALF4,BET1,BET2,BET3,BET4,AS,BS,CS,AA,BB,CC,AX,AY,AZ,BX,
     2	BY,BZ,QH,QK,QL,EN,DH,DK,DL,DE,GH,GK,GL,GMOD
	END

C***********************************************************************
C*   RESCAL:                                                           *
C*     MAIN PROGRAM using an Interactive prompt with commands          *
C***********************************************************************

 	PROGRAM RESCAL
 	
	PARAMETER(NVARS=42)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	
	INTEGER AGM
	
	DIMENSION A(4,4) ! this is the resolution function
	COMMON/MATRIX/A	
	
	COMMON/ERROR/IERR
	
	CHARACTER CUNIT*5
	COMMON/UNIT/CUNIT
	
	COMMON/EXTRA/XXC,YYC,F,W,WI,WF
	
	! variables which are used
	CHARACTER*1024 LINE
	CHARACTER*1024 KEYWORD,VALUE
	CHARACTER*8   NAM(42)
	COMMON/CHNM/NAM
	
	DOUBLE PRECISION PARS
c set default values for parameters (PARS)
	
	COMMON PARS(42)
	DATA PARS/
     1	3.35500, 3.35500, 30.00000, 35.00000, 10.00000, 1.00000, 
     1	1.00000, -1.00000, 2.662, 2.00000, 60.0, 60.00000,
     1	60.0, 60.0, 60.0, 60.0, 60.0, 
     1	60.0, 6.28320, 6.28320, 6.28,90.00000, 90.00000, 
     1	90.00000, 1.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
     1	1.00000, 1.00000, 0.00000, 0.00000, 0.00000,  0.00000, 
     1	0.00000, 0.00000, 1.00000, 0.00000, 0.00000, 1.00000, 
     1	0.00000 /

C INITIALISE ===================================================================
	write (*,*)'                                                    '
      write (*,*)'****************************************************'
      write (*,*)'                                                    '
      write (*,*)'                    RESCAL                          '
      write (*,*)'PROGRAMME  TO CALCULATE  COOPER-NATHANS  RESOLUTION '
	write (*,*)'MATRIX AND/OR RIS0  METHOD  CONTRIBUTIONS TO PHONON '
	write (*,*)'WIDTH. THE WIDTHS  OF BRAGG  SCANS ARE  CALCULATED, '
	write (*,*)'AND THE VANADIUM  WIDTH, FOR A TRIPLE-AXIS SPECTRO- '
	write (*,*)'METER. THE PARAMETERS ARE ENTERED INTERACTIVELY.'
      write (*,*)''
      write (*,*)'(c) M. HARGREAVE & P. HULLAH (1980). E. FARHI (2015)'
      write (*,*)'    Version $Date$'
      write (*,*)'****************************************************'
	WRITE(*,*) 'Type HELP at prompt for help. EXIT to kindly end.'

	F     = .4826
	CUNIT = '[meV]'
	WRITE(*,*) 'Using energy units as: ',CUNIT
	
C-----------------------------------------------------------------------
C   INTERACTIVE COMMAND AND PARAMETER INTERPRETER.
C
C
1     IERR=0
C we read a command
8	FORMAT('*Rescal>',$)
	WRITE(*,8)
	READ(*,'(A)') LINE
	AGM=index(LINE,' ')
	KEYWORD=trim(LINE(:AGM-1))	! the command/variable
	VALUE = trim(LINE(AGM+1:))	! the optional argument
	IF (len_trim(VALUE) .EQ. 0) VALUE=char(0)
	AGM=1

c we treat RESCAL commands
c	HELP
	IF ( KEYWORD(1:1) == 'H' .OR. KEYWORD(1:1) == 'h') THEN
	WRITE(*,*)'RESCAL: Computation of the Cooper-Nathans resolution'
	WRITE(*,*)'(c) M. HARGREAVE & P. HULLAH, (1980). E. FARHI (2015)'
	WRITE(*,*)'Valid commands at prompt:'
	WRITE(*,*)'LIST    - list all parameters <NAME>=<VALUE>'
	WRITE(*,*)'LOAD FN - load parameter file from <FN>'
	WRITE(*,*)'SAVE FN - save parameter file to <FN>'
	WRITE(*,*)'EXIT    - exit (or Ctrl-C)'
	WRITE(*,*)'BRAG 1  - radial, tangential, vertical bragg width, V'
	WRITE(*,*)'          width and bragg energy width'
	WRITE(*,*)'BRAG 2  - bragg width along scan direction'
	WRITE(*,*)'PHON 1  - phonon width by Cooper-Nathans method'
	WRITE(*,*)'PHON 2  - phonon width by RIS0 method'
	WRITE(*,*)'RESOL 1 - resolution volumes'
	WRITE(*,*)'RESOL 2 - res. matrix with x-axis along the Q-vector'
	WRITE(*,*)'RESOL 3 - res. matrix w.r.t. rlu'
	WRITE(*,*)'RESOL 4 - diagonalized res. matrix w.r.t. rlu and'
	WRITE(*,*)'        direction cos(phi) principal axes w.r.t. rlu'
	WRITE(*,*)'<NAME>  - display parameter <NAME> value.'
	WRITE(*,*)'<N> <V> - set parameter <N> value to <V>, e.g. DM 3.3'
	WRITE(*,*)'          Variable names are case sensitive (UPPER)'
	WRITE(*,*)''
	WRITE(*,*)'Commands can be shortened (LIST->LI), case insensitive'
	
	END IF 
c	EXIT
	IF ( KEYWORD(1:2) == 'EX' .OR. KEYWORD(1:2) == 'ex') STOP
c	SAVE <FILE>
	IF ( KEYWORD(1:3) == 'SAV' .OR. KEYWORD(1:3) == 'sav') THEN
	  IF (VALUE .EQ. char(0)) THEN
	    status = getcwd( line )
      if ( status .eq. 0 ) then
        write(*,*) 'Current directory is: '
        write(*,*) trim(line)
      end if
	    WRITE(*,*) 'Enter filename to write (space to abort): '
	    READ(*,'(A)') VALUE
	    END IF
	  IF (len_trim(VALUE) .GT. 0) THEN
	    OPEN(UNIT=24,FILE=trim(value),ERR=601)
	    DO I=1,NVARS
	      WRITE(24,'(F10.5, A12, A8)') PARS(I),'          % ', NAM(I)
	    END DO
	    WRITE(24,*) '% RESCAL Cooper-Nathans parameters (42*F10.5)'
	    CLOSE(UNIT=24)
	    WRITE(*,*) 'Saved parameters to: ', trim(value)
	  END IF
	END IF
c	LOAD <FILE>
	IF ( KEYWORD(1:2) == 'FI' .OR. KEYWORD(1:2) == 'LO' 
     1.OR. KEYWORD(1:2) == 'fi' .OR. KEYWORD(1:2) == 'lo') THEN
	  IF (VALUE .EQ. char(0)) THEN
	    WRITE(*,*) 'The file should contains 42 values (F10.5),'
	    WRITE(*,*) '  one per line, possibly followed by comments.'
	    status = getcwd( line )
      if ( status .eq. 0 ) then
        write(*,*) 'Current directory is: '
        write(*,*) trim(line)
      end if
	    WRITE(*,*) 'Enter filename to read (space to abort): '
	    READ(*,'(A)') VALUE
	  END IF
	  IF (len_trim(VALUE) .GT. 0) THEN
	    OPEN(UNIT=24,FILE=trim(VALUE),ERR=601)
	    DO I=1,NVARS
	      READ(24,'(F10.5)') PARS(I)
	    END DO
	    CLOSE(UNIT=24)
	    WRITE(*,*) 'Loaded parameters from: ', trim(value)
	  END IF
	END IF
c	LIST
	IF ( KEYWORD(1:2) == 'LI' .OR. KEYWORD(1:2) == 'li') THEN
	  DO I=1,NVARS
	    WRITE(*,*) NAM(I),'=',PARS(I) ! display current or new assignment
	  END DO
	END IF
c we treat variable display/assignment (there are 42 parameters with CN)
c	<VAR> or <VAR> <VALUE>
601	DO I=1,NVARS
	  IF ( trim(NAM(i)) == trim(KEYWORD) ) THEN
C     Read next Token, which is char(0) if nothing follows
	    IF (VALUE .NE. char(0)) READ(VALUE, '(E20.10)') PARS(I)
	    WRITE(*,*) NAM(I),'=',PARS(I) ! display current assignment
	  END IF
	END DO
c compute stuff before looping back to prompt

	CALL RECLAT
	IF (IERR.NE.0) THEN
	  WRITE(*,901) IERR, 'RECLAT'
	  GOTO 1
	ENDIF
C   RECIPROCAL LATTICE.
	CALL TRNVCT
	IF (IERR.NE.0) THEN
	  WRITE(*,901) IERR, 'TRNVCT'
	  GOTO 1
	ENDIF
C   TRANSFORM VECTORS Q, D & G.
	CALL AFILL
	IF (IERR.NE.0) THEN
	  WRITE(*,901) IERR, 'AFILL'
	  GOTO 1
	ENDIF
901	FORMAT( ' RESCAL: ERROR ',I4,' IN ',A)

c	Computation calls
      AGM=1
c	BRAG
      IF ( KEYWORD(1:2) == 'BR' .OR. KEYWORD(1:2) == 'br') THEN
	  IF (VALUE .NE. char(0)) READ(VALUE, *) AGM
	  CALL BRAG(AGM)
	END IF
c	PHON
	IF ( KEYWORD(1:2) == 'PH' .OR. KEYWORD(1:2) == 'ph') THEN
	  IF (VALUE .NE. char(0)) READ(VALUE, *) AGM
	  CALL PHON(AGM)
	END IF
c	RESOL
	IF ( KEYWORD(1:2) == 'RE' .OR. KEYWORD(1:2) == 're')  THEN
	  IF (VALUE .NE. char(0)) READ(VALUE, *) AGM
	  CALL RESOL(AGM)
	END IF

	GO TO 1

	END PROGRAM RESCAL

C =====================================================================
C                        COMPUTATIONAL ROUTINES
C =====================================================================

C
C
C =====================================================================
	SUBROUTINE TRNVCT
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION KFIX,EPS1
	PARAMETER(EPS1=1.D-5)
	COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1
     1	,ALF2,ALF3,ALF4,BET(4),A(3),DUMM(9),Q(3),EN,D(3),DE,
     2	G(3),GMOD
	COMMON /EXTRA/XXC,YYC,F,W,WI,WF
	COMMON /S/S(3,3),B(3)
	COMMON/TRSG/SX,SY,SXX,SYY,GX,GY,GXX,GYY
	COMMON /ERROR/IERR
	CHARACTER*40 MATER(5)
	DATA MATER/
	1	' TRNVCT: Q Not In Scattering Plane',
	2	' TRNVCT: Scan Not In Scattering Plane',
	3	' TRNVCT: check  scattering triangle',
	4	' TRNVCT: Gradient too small (GH,GK)',
	5	' TRNVCT: Gradient Has Component In Z-DIR'/
C---------------------------------------------------------------------------
C   TRANSFORM AND CHECK Q, D & G VECTORS.
C
	IERR=0
	XXC=TRANS(Q,1)
	YYC=TRANS(Q,2)
	IF(ABS(TRANS(Q,3)).GT.EPS1)THEN
	IERR=1
	GOTO 99
	ENDIF
C---------------------------------------------------------------------------
C
	SX=TRANS(D,1)
	SY=TRANS(D,2)
	IF(ABS(TRANS(D,3)).GT.EPS1) THEN
	IERR=2
	GOTO 99
	ENDIF
C
	DEN=SQRT(XXC*XXC+YYC*YYC)
	IF (DEN.LT.EPS1) THEN
	IERR=3
	GOTO 99
	ENDIF
	SXX=(SX*XXC+SY*YYC)/DEN
	SYY=(SY*XXC-SX*YYC)/DEN
      GX=TRANS(G,1)
      GY=TRANS(G,2)
      GG=SQRT(GX*GX+GY*GY)
	IF (GG.LT.EPS1) THEN
	IERR=4
	GOTO 99
	ENDIF
      GX=GMOD*GX/GG
      GY=GMOD*GY/GG
	IF(ABS(TRANS(G,3)).GT.EPS1)THEN
	IERR=5
	GOTO 99
	ENDIF
	 GXX= (GX*XXC+GY*YYC)/(DEN )
      GYY= (GY*XXC-GX*YYC)/(DEN )
C---------------------------------------------------------------------------
	W=EN
	WDF=W*F/.4826
	WS=KFIX*KFIX/.4826
	IF(FX.EQ.2) GOTO 10
	WI=WS
	WF=WI-WDF
	GOTO 11
10	WF=WS
	WI=WF+WDF
11	RETURN
C
C
99	CONTINUE
	WRITE(*,501) MATER(IERR)
501	FORMAT(A)
	RETURN
	END SUBROUTINE TRNVCT
C
C
C =====================================================================
	SUBROUTINE BRAG(AGM)
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(4,4)
	CHARACTER*5 CUNIT
	INTEGER AGM
	DOUBLE PRECISION KFIX,EPS1
	PARAMETER(EPS1=1.D-5)
	COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1
     1	,ALF2,ALF3,ALF4,BET(4),DUMMA(3),DUMM(9),Q(3),EN,D(3),DE,
     2	G(3),GMOD
C	COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,
C     1	ALF1,ALF2,ALF3,ALF4,BET1,BET2,BET3,BET4

	COMMON /EXTRA/ XXC,YYC,F,W,WI,WF
	COMMON /UNIT/CUNIT
	COMMON /MATRIX/A
	COMMON /S/S(3,3),B(3)
	COMMON/TRSG/SX,SY,SXX,SYY,GX,GY,GXX,GYY
	LOGICAL LZER
C
C
C
	HUITLOG2=8.D0*LOG(2.D0)
	LZER=.FALSE.
	DO 11 I=1,4
	IF (ABS(A(I,I)).LT.EPS1) LZER=.TRUE.
11	CONTINUE
C
	IF (LZER) THEN
	WRITE(*,*) ' BRAG: Sorry  Pb with matrix -> I am unable to
	1	divide by zero '
	RETURN
	ENDIF
c
	IF(AGM.GT.1)GO TO 4
C---------------------------------------------------------------------------
C   WIDTHS IN RECIPROCRAL ANGSTROMS
	DQX=SQRT(HUITLOG2/A(1,1))
	DQY=SQRT(HUITLOG2/A(2,2))
	DQZ=SQRT(HUITLOG2/A(3,3))
C   ENERGY AND VANADIUM WIDTHS.
	DEE=SQRT(HUITLOG2/A(4,4))
	DVN=(A(2,4)*A(1,1)-A(1,4)*A(1,2))**2/(A(1,1)*(A(2,2)*A(1,1)
     1	-A(1,2)**2))
	DVN=SQRT(HUITLOG2/(A(4,4)-A(1,4)**2/A(1,1)-DVN))
	WRITE(*,5)DQX,DQY,DQZ,CUNIT,DVN,CUNIT,DEE,CUNIT
    5	FORMAT(' BRAGG Widths, Radial,tangential, Vertical (HWHM) [ANG-1]'/
	1 ' DQR=',F9.5,' DQT=',F9.5,' DQV=',F9.5/
	2 ' Energy Widths (HWHM) ',A5/
	3 ' DVN=',F9.5,' ',A5,' DEE=',F9.6,' ',A5/)
	RETURN
C---------------------------------------------------------------------------
C   WIDTH IN TERMS OF NUMBER OF STEPS.
4	SXY=A(1,1)*SXX*SXX+2.*A(1,2)*SXX*SYY+A(2,2)*SYY*SYY
	AB=SQRT(SXY)
	IF(AB.GT.EPS1) THEN
	W=2.3548/AB
	WRITE(*,13)W
13	FORMAT(' FWHM Along (DH,DK,DL): No Of Steps = ',F11.5/)
	ELSE
	WRITE(*,601) D(1),D(2),D(3)
601	FORMAT(' BRAG: CHECK DH DK DL :',3(X,F12.3))
	ENDIF
	RETURN
	END SUBROUTINE BRAG
C
C
C =====================================================================
	SUBROUTINE PHON(AGM)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION KFIX
	DOUBLE PRECISION G1,G2,G3,BM(4,4),BXX,BYY
        DIMENSION AM(4,4)
	CHARACTER*5 CUNIT
	COMMON DM,DA,ETAM,ETAA,ETAS,DUMMY(32),DE
        COMMON/UNIT/CUNIT
	COMMON/EXTRA/XXC,YYC,F,W,WI,WF
	COMMON /MATRIX/AM
	COMMON/TRSG/SX,SY,SXX,SYY,GX,GY,GXX,GYY
	INTEGER AGM
	HUITLOG2=8.D0*LOG(2.D0)
	DO 10 I=1,4
	DO 10 J=1,4
10      BM(I,J)=AM(I,J)
	BXX=GXX
	BYY=GYY
	IF(AGM.GT.1)GO TO 2
C---------------------------------------------------------------------------
C   PHONON WIDTH BY COOPER-NATHANS METHOD.
1	G1=1.0D0/(BM(1,1)+BXX*BXX*BM(4,4)+2.0D0*BXX*BM(1,4))
	G2=BM(1,2)+BXX*BM(2,4)-BYY*(-BM(1,4)-BXX*BM(4,4))
	G2=1.0D0/(BM(2,2)+BYY*BYY*BM(4,4)+2.0D0*BYY*BM(2,4)-G1*G2*G2)
	G3=BM(1,2)+BXX*BM(2,4)-BYY*(-BM(1,4)-BXX*BM(4,4))
	G3=-BYY*BM(4,4)-BM(2,4)-G1*(-BM(1,4)-BXX*BM(4,4))*G3
	G3=BM(4,4)-G1*(-BM(1,4)-BXX*BM(4,4))**2-G2*G3*G3
	WIDTH=HUITLOG2*(SXX**2+SYY**2+DE**2)/G3
	WIDTH=SQRT(WIDTH/(DE-SXX*BXX-SYY*BYY)**2)
        WRITE(*,3)CUNIT,WIDTH
3	FORMAT(' Cooper - Nathans Method; Units Are ANGS-1 ',A5,
     1 ' Or Mixture'/' Phonon Width =',F11.5)
	RETURN
C---------------------------------------------------------------------------
C   PHONON WIDTH BY RISO METHOD.
2	IF(ABS(ETAS).GE.0.001) THEN
	WRITE(*,5)ETAS
5	FORMAT('WARNING RIS0 method does not use ETAS ',F10.5)
4	W=SQRT((SXX*SXX+SYY*SYY+DE*DE)/(DE-SXX*BXX-SYY*BYY)**2)
	END IF
	WRITE(*,6)
6	FORMAT(' RIS0 Method'/' Phonon Width Contributions')
	CALL RESCON
	RETURN
	END SUBROUTINE PHON
C
C
C =====================================================================
	SUBROUTINE RESOL(AGM)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(4,4),ADA(4,4),B(4,4),AZ(4,4)
	CHARACTER*5 CUNIT
	INTEGER AGM
	COMMON /UNIT/CUNIT
	COMMON /MATRIX/A
	CHARACTER*2 IAX(4),FWHM*4
	CHARACTER*2 IAXP(4),STRVI*4,STRVF*4
	CHARACTER*54 MESOR(5)
	DATA IAX/'X','Y','Z','W'/,FWHM/'FWHM'/
	DATA IAXP/'X''','Y''','Z''','W'''/,STRVI/'VI ='/,STRVF/'VF ='/
	DATA MESOR/
C                123456789012345678901234567890123456789012345678901234
	1	' ********  Resolution Volumes,   Units Are   [ANGS-3]',
	2	' ******  Resolution Matrix, X-AXIS Along Q [ANGS-1] &',
	3	' Resolution Matrix, Axes WRT Recip. Lattice R.l.u. & ',
	4	' ******** Diagonalised Resolution Matrix In R.l.u. & ',
	5	'  Direction Cosines Of Axes W.r.t. Reciprocal Lattice'/
C
	GO TO (11,12,13,14),AGM
	IF(I.EQ.1)GO TO 11
	RETURN
C --------------------------------------------------------------------
C	CAS RESOL 1 ON CALCULE LES VOLUMES VI VF
C
11	CALL VIVF(VI,VF)
	WRITE(*,411) MESOR(1),STRVI,VI,STRVF,VF
411	FORMAT(A//2(5X,A,G12.5))
	RETURN
CC
C --------------------------------------------------------------------
C	CAS RESOL 2 ON CALCULE LA MATRICE SUIVANT Q
C
C
12	WRITE(*,413) MESOR(2),CUNIT,IAX
	WRITE(*,424) (IAX(I),(A(I,J),J=1,4),I=1,4)
	RETURN
C --------------------------------------------------------------------
C	CAS RESOL 3 ON CALCULE LA MATRICE SUIVANT THE Reciprocal Lattice
C
C
13	CALL CONVRT(A,AZ)
	WRITE(*,413)MESOR(3),CUNIT,IAX
413	FORMAT(2A//5(10X,A))
412	FORMAT(A//5(10X,A))
	WRITE(*,424) (IAX(I),(AZ(J,I),J=1,4),I=1,4)
	RETURN
C --------------------------------------------------------------------
C	CAS RESOL 4 ON DIAGONALISE
C
14	CALL DIAG(A,ADA,B)
	WRITE(*,413) MESOR(4),CUNIT,IAXP
	DO 31 I=1,4
	WRITE(*,424),IAXP(I),(ADA(I,K),K=1,4)
31	CONTINUE
424	FORMAT(2X,A,4E12.4)
	TMP=LOG(2.D0)
	DO 33 I=1,4
33	ADA(I,I)= 2.D0*SQRT(TMP/ABS(ADA(I,I)))
   	WRITE(*,412) MESOR(5),IAX,FWHM
	DO 35 I=1,4
	WRITE(*,444) IAXP(I),(B(J,I),J=1,4),ADA(I,I)
35	CONTINUE
444	FORMAT (4(2X,A,5F12.5))
 	RETURN
	END SUBROUTINE RESOL
C
C
C =====================================================================
	SUBROUTINE VIVF(VI,VF)
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION KFIX,PI
	PARAMETER(PI=3.141592653589793239D0)
	COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,
     1	ALF1,ALF2,ALF3,ALF4,BET1,BET2,BET3,BET4
	COMMON /EXTRA/ XXC,YYC,F,W,WI,WF
	AKI=SQRT(WI*F/.4826)
	AKF=SQRT(WF*F/.4826)
	PI2=PI*PI
	STHM=PI2/(DM*AKI)**2
	COTM=1./STHM-1.
	STHA=PI2/(DA*AKF)**2
	COTA=1./STHA-1.
	VI=AKI**3*BET1*BET2*ETAM*ALF1*ALF2
	B=STHM*(2.*ETAM)**2+BET1**2+BET2**2
	B=B*(ALF1**2+ALF2**2+4.*ETAM**2)
	VI=VI*SQRT(COTM/B)
	B=STHA*(2.*ETAA)**2+BET3**2+BET4**2
	B=B*(ALF3**2+ALF4**2+4.*ETAA**2)
	VF=AKF**3*BET3*BET4*ETAA*ALF3*ALF4*SQRT(COTA/B)
	RETURN
	END SUBROUTINE VIVF 
C
C
C =====================================================================
	SUBROUTINE CONVRT(A,AZ)
C	CREATE MATRIX T WHICH CONVERTS R.L.U. AXES TO COOPER-NATHANS AXES
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON /S/S(3,3),B(3)
	COMMON /EXTRA/XXC,YYC,F,W,WI,WF
	DIMENSION A(4,4),AZ(4,4)
	DOUBLE PRECISION  CC(4),T(4,4),COS,SIN,DEN
	DEN=SQRT(XXC*XXC+YYC*YYC)
	COS=XXC/DEN
	SIN=YYC/DEN
	DO 1 I=1,4
	DO 1 J=1,4
1	AZ(I,J)=A(I,J)
	DO 11 I=1,3
	T(1,I)=S(1,I)*COS+S(2,I)*SIN
	T(2,I)=S(2,I)*COS-S(1,I)*SIN
	T(3,I)=S(3,I)
	T(4,I)=0.
11	T(I,4)=0.
	T(4,4)=1.
C	POSTMULTIPLY C-N MATRIX BY S MATRIX, RESULT CALLED AZ
	DO 21 I=1,4
	DO 31 J=1,4
	CC(J)=0.
	DO 31 K=1,4
31	CC(J)=CC(J)+AZ(I,K)*T(K,J)
	DO 21 J=1,4
21	AZ(I,J)=CC(J)
C	PREMULTIPLY A BY S TRANSPOSE, RESULT CALLED AZ
	DO 41 I=1,4
	DO 51 J=1,4
	CC(J)=0.
	DO 51 K=1,4
51	CC(J)=CC(J)+AZ(K,I)*T(K,J)
	DO 41 J=1,4
41	AZ(J,I)=CC(J)
C	NOW RESOLUTION MATRIX IS W.R.T. RECIPROCAL CRYSTAL LATTICE.
	RETURN
	END SUBROUTINE CONVRT
C
C
C =====================================================================
	SUBROUTINE DIAG(A,ADA,B)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(16),ADA(16),B(16),ARMAX(16),JRMAX(16)
	DATA N,ND,E/4,4,1.D-29/
	NDN=ND*N
	DO 1 K=1,NDN
	ADA(K)=A(K)
	B(K)=0.
    1	CONTINUE
	DO 2 K=1,N
	KK=K*(ND+1)-ND
	ARMAX(K)=0.
	B(KK)=1.
	DO 3 L=K,N
	IF(L-K)4,3,4
    4	KL=K+ND*(L-1)
	Y=ABS(ADA(KL))
	IF(ARMAX(K)-Y)5,3,3
    5	ARMAX(K)=Y
	JRMAX(K)=L
    3	CONTINUE
    2	CONTINUE
   11	AMAX=0.
	DO 6 K=1,N
	Y=ABS(ARMAX(K))
	IF(AMAX-Y)7,6,6
    7	AMAX=Y
	I=K
    6	CONTINUE
	J=JRMAX(I)
	IF(E-AMAX)8,9,9
    8	NDI=ND*(I-1)
	NDJ=ND*(J-1)
	II=I+NDI
	JJ=J+NDJ
	IJ=I+NDJ
	JI=J+NDI
	AII=ADA(II)
	AJJ=ADA(JJ)
	AIJ=ADA(IJ)
	Y=2.*AIJ
	X=AII-AJJ
	T=DSIGN(1.D0,X)*Y/(ABS(X)+SQRT(X**2+Y**2))
	TSQ=T**2
	C=1./SQRT(ABS(1.+TSQ))
	TY=T*Y
	S=T*C
	CSQ=C**2
	ADA(II)=CSQ*(AII+TY+AJJ*TSQ)
	ADA(JJ)=CSQ*(AJJ-TY+AII*TSQ)
	ADA(IJ)=0.
	ADA(JI)=0.
	DO 10 K=1,N
	JTES=(K-I)*(K-J)
	NDK=ND*(K-1)
	KI=K+NDI
	KJ=K+NDJ
	IF(JTES)13,12,13
   13	JK=J+NDK
	IK=I+NDK
	ADA(KI)=C*ADA(IK)+S*ADA(JK)
	ADA(KJ)=-S*ADA(IK)+C*ADA(JK)
	ADA(JK)=ADA(KJ)
	ADA(IK)=ADA(KI)
   12	X=B(KI)
	B(KI)=C*X+S*B(KJ)
	B(KJ)=-S*X+C*B(KJ)
   10	CONTINUE
	ARMAX(I)=0.
	DO 14 K=1,N
	IF(K-I)15,14,15
   15	IK=I+ND*(K-1)
	Y=ABS(ADA(IK))
	IF(ARMAX(I)-Y)16,14,14
   16	ARMAX(I)=Y
	JRMAX(I)=K
   14	CONTINUE
	ARMAX(J)=0.
	DO 17 K=1,N
	IF(K-J)18,17,18
   18	JK=J+ND*(K-1)
	Y=ABS(ADA(JK))
	IF(ARMAX(J)-Y)19,17,17
   19	ARMAX(J)=Y
	JRMAX(J)=K
   17	CONTINUE
	DO 20 K=1,N
	ITES=(K-I)*(K-J)
	KI=K+NDI
	KJ=K+NDJ
	IF(ITES)21,20,21
   21	X=ABS(ADA(KI))
	Y=ABS(ADA(KJ))
	JR=J
	IF(X-Y)22,22,23
   23	Y=X
	JR=I
   22	IF(ARMAX(K)-Y)24,20,20
   24	ARMAX(K)=Y
	JRMAX(K)=JR
   20	CONTINUE
	GOTO 11
9	CONTINUE
	RETURN
	END SUBROUTINE DIAG
	
C
C
C =====================================================================
	SUBROUTINE RECLAT
C	11-05-78
C	25-06-79 UPDATED
C	OVERLAY #11
C	COMPUTATION OF THE TRANSFORMATION MATRIX B
C	AND ORIENTATION MATRIX S
C	REQUIRED: VECTOR A=UNIT CELL SIDES, VECTOR ALFA=CELL ANGLES
C	A1=VECTOR IN DIRECTION ACO, A2=VECTOR IN
C	X=SUM(S(1,I)*Q(I)) AND Y=SUM(S(2,I)*Q(I)) ARE CARTESIAN
C	COORDINATES IN SCATTERING PLANE, TRANSFORMED FROM VECTOR
C	Q W.R.T. CRYSTAL LATTICE.  SUM(S(3,I)*Q(I)) SHOULD=0
     
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION PARS(42),PI
	PARAMETER(PI=3.141592653589793239D0)
	DIMENSION A(3),BB(3,3),V1(3),V2(3),V3(3),U(3,3),COSA(3)
     1	,SINA(3),B(3),COSB(3),SINB(3)
	COMMON DUMMY(18),AQ(3),ALFA(3),A1(3),A2(3)
	COMMON /S/S(3,3),BLAH(3)
	COMMON /ERROR/IERR

	ZD=2.D0*PI
	RD=2*PI/360.D0
	CC=0.
	DO 1 I=1,3
	A(I)=AQ(I)/ZD
	IF(ABS(A(I)).GT..00001) GOTO 6
	IERR=1
	GO TO 100
6	COSA(I)=COS(ALFA(I)*RD)
	SINA(I)=SIN(ALFA(I)*RD)
    1	CC=CC+COSA(I)*COSA(I)
	CC=1.+2.*COSA(1)*COSA(2)*COSA(3)-CC
	IF(CC.LE.0.) GOTO 23
	CC=SQRT(CC)
	J=2
	K=3
	DO 2 I=1,3
	B(I)=SINA(I)/(A(I)*CC)
	COSB(I)=(COSA(J)*COSA(K)-COSA(I))/(SINA(J)*SINA(K))
	SINB(I)=SQRT(1.-COSB(I)*COSB(I))
	J=K
    2	K=I
	BB(1,1)=B(1)
	BB(2,1)=0.
	BB(3,1)=0.
	BB(1,2)=B(2)*COSB(3)
	BB(2,2)=B(2)*SINB(3)
	BB(3,2)=0.
	BB(1,3)=B(3)*COSB(2)
	BB(2,3)=-B(3)*SINB(2)*COSA(1)
	BB(3,3)=1/A(3)
C	GENERATION OF ORIENTATION MATRIX REC. LATTICE TO SCATTERING PLANE
	DO 3 I=1,3
	C1=0.
	C2=0.
	DO 4 J=1,3
	C1=C1+BB(I,J)*A1(J)
    4	C2=C2+BB(I,J)*A2(J)
	V1(I)=C1
    3	V2(I)=C2
	V3(1)=V1(2)*V2(3)-V1(3)*V2(2)
	V3(2)=V1(3)*V2(1)-V1(1)*V2(3)
	V3(3)=V1(1)*V2(2)-V1(2)*V2(1)
	V2(1)=V3(2)*V1(3)-V3(3)*V1(2)
	V2(2)=V3(3)*V1(1)-V3(1)*V1(3)
	V2(3)=V3(1)*V1(2)-V3(2)*V1(1)
	C1=(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
	C2=(V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
	C3=V3(1)**2+V3(2)**2+V3(3)**2
	IF(ABS(C1)-EPS)14,14,15
   14	IERR=5
	GOTO 100
   15	IF(ABS(C2)-EPS)16,16,17
   16	IERR=6
	GOTO 100
   17	IF(ABS(C3)-EPS)18,18,19
   18	IERR=7
	GOTO 100
   19	CONTINUE
	C1=SQRT(C1)
	C2=SQRT(C2)
	C3=SQRT(C3)
	DO 5 I=1,3
	U(1,I)=V1(I)/C1
	U(2,I)=V2(I)/C2
    5	U(3,I)=V3(I)/C3
	DO 7 K=1,3
	DO 7 M=1,3
	SS=0.
	DO 8 L=1,3
    8	SS=SS+U(K,L)*BB(L,M)
	S(K,M)=SS
    7	CONTINUE
  100	IF(IERR)200,200,20
   20	GOTO(22,22,22,23,24,24,24),IERR
   22	WRITE(*,30)
	GOTO 200
   23	WRITE(*,31)
	IERR=1
	GOTO 200
   24	WRITE(*,32)
  200	CONTINUE
	RETURN
   30	FORMAT(' RECLAT: Check Lattice Spacings (AS,BS,CS)'/)
   31	FORMAT(' RECLAT: Check Cell Angles (AA,BB,CC)'/)
   32	FORMAT(' RECLAT: Check Scattering Plane (AX....BZ)'/)
	END SUBROUTINE RECLAT
C
C
C =====================================================================
	SUBROUTINE RESCON
C	RISO PROGRAM FOR FINDING CONTRIBUTIONS TO WIDTH OF A CONST-Q
C	SCAN OF TRIPLE AXIS SPECTROMETER.
C	A* AND B* - ORTHOGONAL RECIPROCAL LATTICE DIMENSIONS IN SCAT PLANE
C	TAUM & TAUA - 2*PI/(PLANE SPACING) FOR MONO. AND ANAL.
C	R- REACTOR. M - MONOCHROMATOR. S- SAMPLE. A - ANALYSER.
C	D - DETECTOR.
C	** FOR CONVERSION NOTE 4.1377MEV =1THZ **
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (PI=3.141592653589793239D0)
	COMMON AA(5),SM,SS,SA,AAA(6)
C
	CHARACTER*5 CUNIT
C
	COMMON /UNIT/CUNIT
	COMMON /EXTRA/XXC,YYC,FF,W,WI,WF
	COMMON/TRSG/SS1,SS2,SXX,SYY,GG1,GG2,GXX,GYY
	DIMENSION A(6),B(4),C(6),D(7),E(4),F(4)
C
C
C	C2=.0002908882
C	C1=6.2831853
C
	C1=2.*PI
	C2= C1/(360.D0*60.D0)
	A(1)=XXC
	A(2)=YYC
	A(3)=WI
	A(4)=WF
	DO 120 I=1,2
	C(I)=AAA(2+I)*C2
  120	C(I+3)=AAA(7-I)*C2
	C(3)=AA(3)*C2
	C(6)=AA(4)*C2
	T1=C1/AA(1)
	T2=C1/AA(2)
	IF(ABS(A(1))-.000001)163,163,170
163	A(1)=1E-06
170	P1=SS*(-SM)
	P2=-SA*(-SM)
	A(5)=0.695*SQRT(A(3))
	A(6)=0.695*SQRT(A(4))
	A7=SQRT(A(1)**2+A(2)**2)
	A9=ATAN(ABS(A(2)/A(1)))
	A9=SIGN(1.D0,A(1)*A(2))*A9+0.5*PI*(1.-SIGN(1.D0,A(1)))
	C1=(A(5)**2+A7**2-A(6)**2)/(2*A7*A(5))
	B1=ATAN(SQRT(1.-C1**2)/C1)
	IF(B1)300,310,310
  300	B1=B1+PI
  310	C1=(A(6)**2+A7**2-A(5)**2)/(2*A(6)*A7)
	B2=ATAN(SQRT(1.-C1**2)/C1)
	IF(B2)350,360,360
  350	B2=B2+PI
  360	S1=T1/(2*A(5))
	B(1)=ATAN(S1/SQRT(1.-S1**2))
	S1=T2/(2.*A(6))
	B(3)=ATAN(S1/SQRT(1.-S1**2))
	DO 630 J1=1,2
	J2=(J1-1)*3
	I1=1+J2
	I2=2+J2
	I3=3+J2
	I4=4+J2
	D(I1)=C(I1)*C(I2)/SQRT(C(I1)**2+C(I2)**2)
	C1=SQRT(4./(C(I1)**2+C(I2)**2)+1/C(I3)**2)
	D(I2)=(C(I2)**2-C(I1)**2)/(C(I2)**2+C(I1)**2)/C1
	D(I3)=2.*C(I2)**2/(C(I1)**2+C(I2)**2)/C1
	IF(C(I2)-C(I1))510,530,510
  510	D(I4)=2.*C(I2)**2/(C(I2)**2-C(I1)**2)
	GOTO 540
  530	D(I4)=1.E37
  540	J3=1+2*(J1-1)
	E(J3)=2.*A(2+J1)*D(I1)/TAN(B(J3))
	E(J3+1)=2.*A(2+J1)*D(I2)/TAN(B(J3))
	F(J3)=A(4+J1)*D(I1)/SIN(B(J3))
	B(J3+1)=ATAN(D(I4)*TAN(B(J3)))
	IF(B(J3+1))610,620,620
  610	B(J3+1)=B(J3+1)+3.14159265
  620	F(J3+1)=A(4+J1)*D(I3)/SIN(B(J3+1))
  630 	CONTINUE
	S1=GG1*FF/.4826
	S2=GG2*FF/.4826
	S1P=S1+(1+SM)*(S2*COS(A9)-S1*SIN(A9))*SIN(A9)
	S2=S2+(1+SM)*(S1*SIN(A9)-S2*COS(A9))*COS(A9)
	S1=S1P
	X=A9+B1*P1+B(1)
	W1=-S1*F(1)*COS(X)-S2*F(1)*SIN(X)+E(1)
	X=A9+B1*P1+B(2)
	W2=-S1*F(2)*COS(X)-S2*F(2)*SIN(X)+E(2)
	X=A9+(3.14159265-B2*P1)-B(3)*P2
	W3=-S1*F(3)*COS(X)-S2*F(3)*SIN(X)+E(3)
	X=A9+(3.14159265-B2*P1)-B(4)*P2
	W4=-S1*F(4)*COS(X)-S2*F(4)*SIN(X)+E(4)
	FACT=W*.4826/FF
	W1=SQRT(W1*W1+W2*W2)*FACT
	W2=SQRT(W3*W3+W4*W4)*FACT
	WW=SQRT(W1*W1+W2*W2)
	WRITE(*,101)CUNIT
  101	FORMAT(' RESCON: Units Are ANG-1, ',A,' OR Mixture')
	WRITE(*,100)W1,W2,WW
  100	FORMAT(' W( Mono.)= ',F10.5/
     1 ' W( Anal.)= ',F10.5/' W( Total)= ',F10.5)
	RETURN
	END SUBROUTINE RESCON
C
C
C =====================================================================
	SUBROUTINE AFILL
C	THIS GENERATES COOPER-NATHANS RESOLUTION MATRIX.
C	AGREES WITH MATRIX GENERATED BY INDEPENDENT METHOD BY R. PYNN.
C	COORDINATE AXES REVERSED W.R.T. COOPER-NATHANS CONVENTION.
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION ALP(3),BET(3),GAM(3),A(4,4),AD(4,4),ALG(8)
	DOUBLE PRECISION MOD,KFIX,PI
	PARAMETER (PI=3.141592653589793239D0)
	COMMON/EXTRA/XXC,YYC,F,W,WI,WF
	COMMON DM,DA,ETM,ETA,ETS,SM,SS,SA,KFIX,FX,ALF(8)
	COMMON/MATRIX/A
   	COMMON /ERROR/IERR
	EQUIVALENCE (ALG(1),ALZ),(ALG(2),ALM),(ALG(3),ALA)
	1 ,(ALG(4),AL3),(ALG(5),BET0),(ALG(6),BET1),(ALG(7),BET2)
	2 ,(ALG(8),BET3)

	FOC=-SS*(-SM)
	MOD=SA*(-SM)
	PIT=2.D0*PI/(360.D0*60.D0)
	DO 9 I=1,8
    9	ALG(I)=ALF(I)*PIT
	ETAM=ETM*PIT
	ETAA=ETA*PIT
	ETAS=ETS*PIT
	Q=SQRT(XXC*XXC+YYC*YYC)
	AOM=W*F
	IF(FX.EQ.2.)GOTO 5
	AKI=KFIX
	AKF=AKI*AKI-AOM
	IF(AKF.LT.0.) GOTO 6
	AKF=SQRT(AKF)
	GOTO 7
    5	AKF=KFIX
	AKI=AKF*AKF+AOM
	IF(AKI.LT.0.) GOTO 6
	AKI=SQRT(AKI)
    7	ALAM=AKI/AKF
	BE=-(Q*Q-2.*AKI*AKI+AOM)/(2.*AKI*AKF)
	IF(ABS(BE).GT.1.0)GO TO 10
	AL=SQRT(1.-BE*BE)*FOC
	B =-1.*(Q*Q-AOM)/(2.*Q*AKF)
	IF(ABS(B).GT.1.0)GO TO 10
	AA=SQRT(1.-B*B)*FOC
	ALP(1)=-B/AL
	ALP(2)=AA/AL
	ALP(3)=1./(AL*2.*AKF)
	SB=(Q*Q+AOM)/(2.*Q*AKI)
	IF(ABS(SB).GT.1.0)GO TO 10
	SASA=SQRT(1.-SB*SB)*FOC
	BET(1)=-SB/AL
	BET(2)=SASA/AL
	BET(3)=BE/(AL*2.*AKF)
	GAM(1)=0.
	GAM(2)=0.
	GAM(3)=-1./(2.*AKF)
	IF(AKI.LT.(PI/DM).OR.AKF.LT.(PI/DA)) GO TO 12
	TOM=-PI/(DM*SQRT(AKI*AKI-(PI/DM)**2))
	TOA=(MOD*PI/DA)/SQRT(AKF*AKF-(PI/DA)**2)
	A1=TOM/(AKI*ETAM)
	A2=1./(AKI*ETAM)
	A3=1./(AKI*ALM)
	A4=1./(AKF*ALA)
	A5=TOA/(AKF*ETAA)
	A6=-1./(AKF*ETAA)
	A7=2.*TOM/(AKI*ALZ)
	A8=1./(AKI*ALZ)
	A9=2.*TOA/(AL3*AKF)
	A10=-1./(AL3*AKF)
	B0=A1*A2+A7*A8
	B1=A2*A2+A3*A3+A8*A8
	B2=A4*A4+A6*A6+A10*A10
	B3=A5*A5+A9*A9
	B4=A5*A6+A9*A10
	C=-1.*(ALAM-BE)/AL
	E=-1.*(BE*ALAM-1.)/AL
	AP=A1*A1+2.*B0*C+B1*C*C
	1+B2*E*E+B3*ALAM*ALAM+2.*B4*ALAM*E+A7*A7
	D0=B1-1./AP*(B0+B1*C)**2
	D1=B2-1./AP*(B2*E+B4*ALAM)**2
	D2=B3-1./AP*(B3*ALAM+B4*E)**2
	D3=2.*B4-2./AP*(B2*E+B4*ALAM)*(B3*ALAM+B4*E)
	D4=-2./AP*(B0+B1*C)*(B2*E+B4*ALAM)
	D5=-2./AP*(B0+B1*C)*(B3*ALAM+B4*E)
	DO 3 I = 1,3
	DO 3 J = 1,3
	II=I
	JJ=J
	IF(I.EQ.3)  II=4
	IF(J.EQ.3)  JJ=4
    3   A(II,JJ)=D0*ALP(I)*ALP(J)
	1+D1*BET(I)*BET(J)+D2*GAM(I)*GAM(J)
	2+.5*D3*(BET(I)*GAM(J)+BET(J)*GAM(I))
	3+.5*D4*(ALP(I)*BET(J)+ALP(J)*BET(I))
	4+.5*D5*(ALP(I)*GAM(J)+ALP(J)*GAM(I))
	I1=1
	I2=2
	I3=3
	SOM=-PI/(DM*AKI)
	SOA=MOD*PI/(DA*AKF)
	DUM1=(4*SOM**2*ETAM**2+BET0**2)*AKI**2
	A112=1./DUM1+1./(BET1*AKI)**2
	DUM1=(4*SOA**2*ETAA**2+BET3**2)*AKF**2
		A122=1./DUM1+1./(BET2*AKF)**2
	A(3,3)=A112*A122/(A112+A122)
	DETR=0.
	DETB=1.
	DETC=1.
	IF(ETAS.LT.0.00005) GO TO 4
 	DETS=1./(Q*ETAS)**2
	DETR=1./(DETS+A(2,2))
	DETB=DETS/(DETS+A(2,2))
	DETC=DETS/(DETS+A(3,3))
    4	FAC=SQRT(DETB*DETC)
	AD(1,1)=A(1,1)-A(1,2)**2*DETR
	AD(1,2)=A(1,2)*DETB
	AD(1,3)=A(1,3)
	AD(1,4)=(A(1,4)-A(1,2)*A(2,4)*DETR)*F
	AD(2,1)=AD(1,2)
	AD(2,2)=A(2,2)*DETB
	AD(2,3)=A(2,3)
	AD(2,4)=A(2,4)*DETB*F
	AD(3,1)=AD(1,3)
	AD(3,2)=AD(2,3)
	AD(3,3)=A(3,3)*DETC
	AD(3,4)=A(3,4)
	AD(4,1)=AD(1,4)
	AD(4,2)=AD(2,4)
	AD(4,3)=AD(3,4)
	AD(4,4)=(A(4,4)-A(4,2)**2*DETR)*F**2
	DO 22 I=1,4
	DO 22 J=1,4
22	A(I,J)=5.545*AD(I,J)
	DO 21 I=1,4
	A(I,2)=A(I,2)*(-SM)
21	A(2,I)=A(2,I)*(-SM)
	RETURN
6	WRITE(*,8)KFIX,W,FX
8	FORMAT(' *** AFILL: Error: KFIX,EN,FX= ',3F12.6)
	IERR=1
	RETURN
10	WRITE(*,11)AKI,AKF,Q
11	FORMAT(' *** AFILL: KI,KF,Q Triangle Will Not Close'/' KI,KF,Q =',
     1 3G14.7)
	IERR=1
	RETURN
   12 WRITE(*,13)
   13 FORMAT(' *** AFILL: KI OR KF Cannot Be Obtained')
	IERR=1
	RETURN
	END SUBROUTINE AFILL
C
C
C =====================================================================
	FUNCTION TRANS(A,I)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(3)
	COMMON/S/S(3,3),B(3)
	TRANS=S(I,1)*A(1)+S(I,2)*A(2)+S(I,3)*A(3)
	RETURN
	END FUNCTION TRANS
	

