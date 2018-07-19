  !     ****************************************************************

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2016-12-20  Time: 18:18:16 (by J. Ollivier, ILL)

  !     PROGRAMM MUPHOCOR/5       REICHARDT 5.11.1991
  !     ****************************************************************


  ! Entry parameters:
  !
  !      E0:INC. ENERGY(MEV)             FP:FLIGHTPATH(CM)
  !      XNEL:CHANNEL OF ELASTIC LINE    CW:CHANNEL WIDTH(MYS)
  !      FIMI:MIN.SCATT.ANGLE            FIMA:MAX.SCATT.ANGLE
  !      UNT:CONST.BACKGROUND
  !      NSPEC: NUMBER OF DIFFERENT ATOMIC SPECIES
  !      HOX:   UPPER LIMIT OF DENSITY OF STATES
  !      TEMPO: TEMPERATURE IN KELVIN
  !      DW:    MEAN DEBY WALLER COEFFICIENT
  !      AMASI: ATOMIC MASS
  !      SIGI : SIGMA
  !      CONC:  CONCENTRATION
  !      ALFI:  SCATTERING POWER
  !      NPHO:     NUMBER OF MULTI PHONON TERMS,ITM:NUMBER OF ITERATIONS
  !      ITM:      TOTAL NUMBER OF ITERATIONS
  !      IVIT=0(1):ITERATION BY DIFFERENCE(QUOTIENT) METHOD
  !      IDW=0:    DW COEFF. KEPT CONST.=INPUT VALUE,DW=1:DW COEFF.ITERATED
  !      IEMP=1(0):CORRECTIONS FOR COUNTER EFFICIENCY WILL BE(NOT BE) DONE
  !      IRES=1(0):DATA WILL BE(NOT BE) CORRECTED FOR SPECTR. RESOLUTION
  !
  !COMPILATION:
  !
  !	make (with the appropriate Makefile)
  !
  !USAGE:
  !
  ! ./mupho2017 < input_md.txt > mupho.log
  !
  !

  DIMENSION fnu(1024),quo(1024),gnu0(1024),gnu(1024),gav(1024)
  DIMENSION gges(1024),dsig(1024)
  DIMENSION fin(1024),fout(1024)
  DIMENSION cmo(10),gnum(1024),gnup(1024)
  DIMENSION tetad(100),tempur(100),cv(100),cvr(100)
  DIMENSION quo1(1024),fnu1(1024),fpi(2,1024),quofg(1024),gns(1024)
  DIMENSION dwi(2),rmai(2,20),gai(2,1024),gesi(2,1024)
  DIMENSION gom1(2,20,1024),rm1(20),c(20),rm2(20)
  DIMENSION eres(3000)

  COMMON/reso/eres,ifold
  COMMON/multi/gnd(1024),lm,lm1,lm2,dele
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)
  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv, &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss
  COMMON/bb2/hoe(1024),efin(1024),exe(1024),exa(1024),com,coma,comi,  &
       qc(1024),qma(1024),qmi(1024),dqup(1024,20)
  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff
  COMMON/par2/parq(1024),iquo,lpar

  !     NUMBER OF MULTI-PHONON TERMS RESTRICTED TO 20 DUE TO DIMENSION

12 FORMAT(16I5)
20 FORMAT(8F10.4)
26 FORMAT(8F10.1)
31 FORMAT(10F12.3)
222 FORMAT(//,' ****************************************************')
45 FORMAT(//,' ******************************************************  &
       ****************************************************************')


  ifold=0
  DO i=1,1024
     dsig(i)=0.
  END DO

  DO i=1,3000
     eres(i)=0.
  END DO

  !     Input of data
  CALL indata

  !     Input of LTEM temperature values for which specific heat should be
  !     calculated

  READ 12,ltem
  IF(ltem > 0) READ 20,(tempur(l),l=1,ltem)

  !     END OF THE INPUT DATA  IF IRES=0 AND SECOND PART IS NOT USED

  igr=itm+1
  IF(ivit > 1)igr=ivit

  !     reminder: LM is the number of descret energy intervals

  lm2=2*lm
  lm1=lm+1
  lm5=lm2

  !     In the case of energy gain spectra (IE0<0) LM5 is at maximum
  !     equal to the incident energy divided by the width of the energy
  !     interval

  IF(ie0 > 0) GO TO 225
  lm5=e0/dele
  IF(lm5 > lm2) lm5=lm2
225 CONTINUE

  !     Reminder: NM = Number of TOF input data
  npho1=npho-1

  ngu=hot(1)/dele+.5
  ngo=hot(nm)/dele+.5

  ngu1=ngu+1
  ngo1=ngo+1

  PRINT 240,ngu,ngo
240 FORMAT(2I5)

  !     In subroutine PREP1 the temperature factor plus several other
  !     functions of Q or omega are calculated and stored in the
  !     respective data arrays

  CALL prep1(lm2,lm5,e0,dele)

  ngoo=MIN0(lm5,ngo)

  PRINT 241,ngoo
241 FORMAT(//,' NGOO:',i3)

  !     * Extrapolation of the energy spectrum from the time spectrum ***
  !     * This point is very important in the case of delta(w)< delta(t)*
  !     L: energy channels, N: time channels


  !loop257:  DO  l=ngu1,ngoo
  !  DO  n=1,nm
  !    dif=hoe(l)-hot(n)
  !    IF(dif <= 0.0) THEN
  !      z1(l)=z0(n-1)+(z0(n)-z0(n-1))*(hoe(l)-hot(n-1))/(hot(n)-hot(n-1))
  !      CYCLE loop257
  !    END IF
  !  END DO
  !END DO loop257

  DO  l=ngu1,ngoo
     DO  n=1,nm
        dif=hoe(l)-hot(n)
        IF(dif <= 0.0) THEN
           z1(l)=z0(n-1)+(z0(n)-z0(n-1))*(hoe(l)-hot(n-1))/(hot(n)-hot(n-1))
           EXIT
        END IF
     END DO
  END DO



  !     *****************************************************************

  DO  l=ngo1,lm5
     z1(l)=0.
  END DO

  DO  l=ngu1,lm5
     ef=efin(l)
     ho=ie0*hoe(l)
     ex=EXP(.5*ho/temp)
     !       Reminder: delta(w) is proportional to (E_f*SQRT(E_f))^(-1)
     z1(l)  = z1(l)*(e0/ef)*SQRT(e0/ef)
     gnu0(l)= z1(l)*exe(l)*ex/dqup(l,1)
  END DO

  !     ***** w^2 extrapolation for small w ****
  afa=0.

  DO  i=1,3
     afa=afa+gnu0(ngu+i)/(hoe(ngu+i))**2
  END DO
  afa=afa/3.
  DO  l=1,ngu
     gnu0(l)=afa*hoe(l)**2
     z1(l)=0.
  END DO

  !     ***************************************

  IF(ivit /= 0) THEN
     DO  l=1,lm5
        IF(gnu0(l) < 0.)gnu0(l)=0.
     END DO
  END IF

  !     ******** Calculation of M^(-n) ********
  !     ***** in units of the neutron mass ****
  DO  n=1,npho
     DO  i=1,nspec
        rmai(i,n)=(1.0087/amasi(i))**n
     END DO
  END DO
  !     ***************************************

  DO  l=lm1,lm2
     gnu(l)=0.
  END DO

  !     ********* First guess for g(w) *********
  !     Reichardt: L=1, LM5 !!!!
  DO  l=1,lm
     gnu(l)=gnu0(l)*EXP(.625*qc(l)*dw)
  END DO

  DO  l=lm-20,lm
     gnu(l)=(gnu(lm-20)/400.)*(lm-l)**2
  END DO

  sum=0.
  DO  l=1,lm
     sum=sum+gnu(l)
  END DO
  sum=sum*dele


  DO  l=1,lm2
     gnu(l)=gnu(l)/sum
     gnum(l)=0.
     fnu1(l)=gnu(l)
  END DO
  gnum(1)=gnu(1)*dele
  DO  l=2,lm
     gnum(l)=gnum(l-1)+gnu(l)*dele
  END DO

  CALL corrfa(lm ,gnu,dele,quofg,nspec)
  CALL separ(lm,gnu,dele,quofg,fpi,nspec)
  PRINT 45
  WRITE(6,308)
308 FORMAT(//,'            OMEGA  TOF DISTR.   INPUT.DISTR.  GNU-0.APPR.  &
       INT(OM)     FNU/GNU      FNU1          FNU2         ')
  WRITE(6,309)(l,hoe(l),z1(l),gnu0(l),gnu(l),gnum(l),quofg(l),fpi(1,l),fpi(2,l),l=1,lm5)
309 FORMAT(i5,2F12.3,2E14.5,2F12.4,2E14.5)
  ilau=0

  !**********************************************************************
  !     BEGINNING OF THE ITERATIONS *************************************
  !**********************************************************************

320 CONTINUE
  ilau=ilau+1

  PRINT 330,ilau
330 FORMAT(//,' ITERATION NUMBER:',i3)

  IF(ilau > igr)ivit=0
  IF(ilau == 1) GO TO 345

  !     ********* Normalization of GNU *********
  sum=0.
  DO  l=1,lm
     sum=sum+gnu(l)
  END DO
  sum=sum*dele

  DO  l=1,lm
     gnu(l)=gnu(l)/sum
  END DO

  CALL corrfa(lm ,gnu,dele,quofg,nspec)
  CALL separ(lm,gnu,dele,quofg,fpi,nspec)

  !     ****************************************

345 CONTINUE

  !     ********* Calulation of T_1 = F(w)/(2hw*sinh(hw/2kT)) ********

  !     Reichardt: L=1, LM2 !!!!
  DO  l=1,lm
     rexe=1./exe(l)
     gns(l)=gnu(l)*rexe
     DO  i=1,nspec
        fpi(i,l)=fpi(i,l)*rexe
     END DO
  END DO

  !     **************************************************************

  !     ******** Calculation of the Deby-Waller factor **************
  !     *** DW = Int(h^2/2M * F(w)/w * cosh(hw/2kT))/sinh(hw/2kT) ***
  !     Attention: the integral runs from 1 to LM, i.e. the cut-off *
  !                of the spectrum as defined by the user input *****

  DO  i=1,nspec
     dwi(i)=0.
     DO  l=1,lm
        dwi(i)=dwi(i)+fpi(i,l)*exa(l)
     END DO
     dwi(i)=dwi(i)*dele*1.0087/amasi(i)
  END DO

  dw1=0.

  DO  l=1,lm
     dw1=dw1+gns(l)*exa(l)
  END DO

  dw1=dw1*dele*rma

  PRINT 360,dw1,dwi(1),dwi(2)
360 FORMAT(/,' DEBYE-WALLER COEFF. FROM ITERATED DISTR.(GNU,FNU1,FNU2)  &
       :',3E14.5)

  !     **************************************************************

  DO  i=1,nspec

     !     **********Storage of T_1(-LM,+LM) in the array GND*************

     DO  l=1,lm
        gnd(l)=fpi(i,lm+1-l)
     END DO
     DO  l=lm1,lm2
        gnd(l)=fpi(i,l-lm)
     END DO

     !     ***************************************************************

     DO  l=1,lm2
        gom1(i,1,l)=fpi(i,l)
     END DO

     !       Here is the loop over the multiphonon terms.
     !       GND is the so far determined density of states function devided
     !       by h*omega and sinh(h*omega/2kT), i.e. T_1
     !       GOM1(n) contains the different T_n terms.
     rfc=1.

     DO  n=1,npho1
        !         RFC = 1/(n+1)!
        rfc=rfc/(n+1.)
        CALL vipho(gnum,gnup,n)
        DO  l=1,lm2
           gom1(i,n+1,l)=gnup(l)*rfc
           gnum(l)=gnup(l)
        END DO
     END DO

     IF(idw == 1) dw=dwi(i)

     IF(ipr /= 0 .AND. ilau == 1) THEN
        PRINT 390
390     FORMAT(//,'        OMEGA      F(OM)        G1(OM)       G2(OM)  &
             G3(OM)       G4(OM)       G5(OM)       G6(OM)       G7(OM) G8(OM)')
        WRITE(6,392)(l, hoe(l),fnu(l),(gom1(i,n,l),n=1,8),l=1,lm2)
392     FORMAT(i5,f9.3,9E13.4)
        WRITE(6,222)
     END IF

     !       Construction of the scattering law using eq.2.4

     DO  l=1,lm5
        gav(l) = 0.
        gges(l) = 0.
        ho = ie0*hoe(l)
        ex=EXP(.5*ho/temp) !  EX= detailed balance factor

        !         Factors related to the integration over Q. Attention: w dependent!
        qrma=dw*qma(l)
        qrmi=dw*qmi(l)
        exma=EXP(-qrma)
        exmi=EXP(-qrmi)
        cof0=exmi-exma    !  Here COFO = I_0
        DO  n=1,npho     !  Calculation of (I_N)(x) = -x^N*exp(-x) - N*I_(N-1)
           cof0 = n*cof0+exmi*qrmi**n-exma*qrma**n
           gav(l) = gav(l) + rmai(i,n) * dqup(l,n) * gom1(i,n,l)
           gom1(i,n,l) = (1/ex)*cof0*rmai(i,n)* gom1(i,n,l)/(2.*dw**(n+1))
           gges(l) = gges(l) + gom1(i,n,l)
        END DO
        gav(l)=gav(l)/(dqup(l,1)*rmai(i,1))
     END DO

     IF(ifold == 1) THEN
        CALL fold(gav,fout,lm5,dele)
        DO l=1,lm5
           gav(l) = fout(l)
        END DO
        CALL fold(gges,fout,lm5,dele)
        DO l=1,lm5
           gges(l) = fout(l)
        END DO
     END IF

     DO l=1,lm5
        ho=ie0*hoe(l)
        ex=EXP(.5*ho/temp)
        gges(l)=gges(l)/(dqup(l,1)*rmai(i,1))
        gai(i,l)=gav(l)*exe(l)
        gesi(i,l)=gges(l)*exe(l)*ex
     END DO

     !     ****** Creation of plot file for multi-phonon contributions *****

     IF(ilau >= itm) THEN
        IF(ifold == 1) THEN
           DO jj=1,npho
              DO l=1,lm2
                 fin(l) = gom1(i,jj,l)
              END DO
              CALL fold(fin,fout,lm2,dele)
              DO l=1,lm2
                 gom1(i,jj,l) = fout(l)
              END DO
           END DO
        END IF

        sumz = 0.
        sumc = 0.
        DO jj = 1, lm2
           sumz= sumz + z1(jj)
           DO jjj = 1,npho
              sumc = sumc + gom1(i,jjj,jj)
           END DO
        END DO
        sumz=sumz*dele
        sumc=sumc*dele

        DO l=1,lm5
           DO n=1,npho
              dsig(l)=dsig(l)+alfi(i)*gom1(i,n,l)
           END DO
        END DO

        CLOSE(12)
        IF(i == 1) OPEN(12,FILE= 'multi1.plot')
        IF(i == 2) OPEN(12,FILE=  'multi2.plot')
        REWIND 12
        WRITE(12,403)(l, hoe(l),z1(l)/sumz,(gom1(i,1,l)+gom1(i,2,l)+  &
             gom1(i,3,l)+gom1(i,4,l))/sumc, (gom1(i,n,l)/sumc,n=1,4),l=1,lm2)
403     FORMAT(i5,f9.3,6E13.4)
        CLOSE(12)
     END IF
  END DO  ! *** End of loop over I **


  DO  l=1,lm5
     quo(l)=1.
     gav(l)=gai(1,l)*alfi(1)+gai(2,l)*alfi(2)
     gges(l)=gesi(1,l)*alfi(1)+gesi(2,l)*alfi(2)
     IF(gges(l) /= 0.) quo(l)=gnu0(l)/gges(l)
     !        print*,hoe(l),gnu0(l),gges(l),quo(l)
  END DO

  ! PRINT*,'in it'
  !       stop
  sum=0.
  sum2=0.
  quav=0.


  !     Reichardt only normalized up to LM !!!!!!!!!!!
  DO  l=1,ngoo
     sum=sum+gges(l)
     sum2=sum2+gnu0(l)
  END DO
  DO  l=1,lm
     quav=quav+quo(l)
     !        print*,l,hoe(l),quo(l),quav
  END DO

  !      stop
  rsum=sum/sum2
  quav=rsum*quav/lm

  !     Reichardt: L=1, LM5 !!!!

  DO  l=1,lm
     gnum(l)=0.
     gnup(l)=0.
     quo1(l)=0.

     !       ******** New guess for f(w) using the quotient method ********
     !  IF(ivit == 0)GO TO 431
     ! gnu(l)=gnu(l)*quo(l)*rsum
     ! GO TO 432

     !       ******* New guess for f(w) using the difference method *******
     !  431   CONTINUE
     !  gnu(l)=gnu(l)+gav(l)*(rsum*quo(l)-1.)
     !  IF(gnu0(l) == 0.)gnu(l)=0.
     !  432   CONTINUE

     ! This above replaced by:

     IF(ivit /= 0) THEN    ! New guess for f(w) using the quotient method
        gnu(l)=gnu(l)*quo(l)*rsum
     ELSE                  !  New guess for f(w) using the difference method *******
        gnu(l)=gnu(l)+gav(l)*(rsum*quo(l)-1.)
        IF(gnu0(l) == 0.)gnu(l)=0.
     END IF

     IF(gnu(l) <= 0.) CYCLE
     gnupc  =gnu(l)/(3.*(hoe(l)**2))
     gnup(l)=1./gnupc**.3333333

     !       **************************************************************

     !       Reminder: FNU1 is the first guess for f(w)
     !       Attention GNU not yet normalized

     IF(fnu1(l) /= 0.) quo1(l)=gnu(l)/fnu1(l)
  END DO

  !     Reichardt: These loops are missing !!!!
  !     for a smooth cut off
  DO  l=lm1,lm2
     gnu(l)=0.
  END DO
  DO  l=lm-20,lm
     gnu(l)=(gnu(lm-20)/400.)*(lm-l)**2
  END DO



  !     ********** Calculation of first moment of f(w) **********
  !     s.a. 248
  DO  l=1,lm5
     gnum(l)=0.
  END DO
  gnum(1)=gnu(1)*dele
  DO  l=2,lm
     gnum(l)=gnum(l-1)+gnu(l)*dele
  END DO
  !     *********************************************************

  sum1=0.
  DO   l=1,lm
     sum1=sum1+ABS(rsum*quo(l)-quav)
  END DO
  sum1=sum1/(lm*quav)

  WRITE(6,444),rsum
444 FORMAT(/,' Correction factor=',e14.5)

  WRITE(6,445),quav,sum1
445 FORMAT(/,' AVERAGE OF (INP./CALC.)=',e14.5,'   MEAN RELATIVE  &
       DEVIATION FROM AVERAGE=',e14.5)

  IF(ipr /= 0) GO TO 450

  !     **** ILAU is the number of the actual iteration; ****
  !     **** ITM is the maximum number of iterations;    ****
  !     **** Next iteration if ILAU less than ITM ****
  IF(ilau < itm) GO TO 320


450 CONTINUE

  CALL corrfa(lm,gnu,dele,quofg,nspec)
  CALL separ (lm,gnu,dele,quofg,fpi,nspec)

  DO  l=1,lm5
     fnu(l)=conci(1)*fpi(1,l)+conci(2)*fpi(2,l)
  END DO

  WRITE(6,455)
455 FORMAT(//,'         OMEGA   INP. DIST.  CALC. DIST.   INP./CALC.  &
       GNU(OM)       INT(OM)   DEBYE-OM  GNU-N/GNU-0  FNU/GNU      FNU')
  WRITE(6,460)(l, hoe(l),gnu0(l)*rsum,(gges(l)/sum),(rsum*quo(l)),  &
       (gnu(l)/rsum),gnum(l),gnup(l),quo1(l),quofg(l),fnu(l),l=1,lm5)
460 FORMAT(i5,f9.3,3E13.4, 2X,e13.4,f12.4,f10.2,2F11.3,e14.4)

  !     **** Next iteration if ILAU less than ITM ****
  IF(ilau < itm) GO TO 320


  sdsig=0.
  DO l=1,lm5
     sdsig = sdsig + dsig(l)
  END DO
  sdsig=sdsig*dele

  DO l=1,lm2
     dsig(l) = dsig(l)/sdsig
  END DO

  OPEN(14,FILE = 'DSDO_f90.plot')
  REWIND 14

  WRITE(14,408)(l, hoe(l),z1(l)/sumz, dsig(l), l=1,lm2)
408 FORMAT(i5,f9.3,2E13.4)
  CLOSE(14)

  OPEN(10,FILE = 'GDOS_plot_f90.dat')
  REWIND 10
  ! WRITE(10,461) lm5
  ! 461  FORMAT(i5)
  WRITE(10,462)(l, hoe(l), gnu(l), (gnu0(l)/(dele*sum2)),  &
       (gges(l)/(dele*sum)), l=1,lm5)
462 FORMAT(i5,f9.3,3E13.4)
  CLOSE(10)


  ie=1.3*lm
  !      CALL SPLOT(1,IE,HOE,GNU,0.,GNUM,1.)
  DO  l=1,ie
     gnup(l)=0.
  END DO
  !     CALL LINPLT(
  !    1' PHONON DENSITY OF STATES FOR EQUIDISTENT ENERGY INTERVALS $',
  !    21,IE,3,IE,0,HOE,GNU,GNUP)
  CALL moment(gnu,hoe,lm,dele,cmo)
  WRITE(6,480)
480 FORMAT(//,'DEBYE CUTOFF FREQUENCIES (IN MEV AND THZ)= (.33333*(  &
       N+3)*(OM**N)AV)**1/N                                             ')
  WRITE(6,485)
485 FORMAT(/,'ARR.     -2        -1         0        +1        +2  &
       +3        +4        +5        +6        +7                ')
  WRITE(6,490)(cmo(m),m=1,10)
490 FORMAT(4X,10F10.3)
  DO  m=1,10
     cmo(m)=cmo(m)/4.1355
  END DO
  WRITE(6,490)(cmo(m),m=1,10)

  IF(ltem == 0) GO TO 505
  DO  l=1,ltem
     te=tempur(l)
     CALL heat(lm,hoe,dele,fnu,te,clv,tet)
     cv(l)=clv
     tetad(l)=tet
     cvr(l)=clv/(te**3)
  END DO
  PRINT 45
  PRINT 495
495 FORMAT(//,' SPECIFIC HEAT IN MJ/(K*G-ATOM) : ')
  PRINT 496
496 FORMAT(//,'       T(K)       T**2        CV(T)        CV(T)/T  &
       CV(T)/T**3     THETA-D(K)  ')
  DO  l=1,ltem
     ttt2=tempur(l)**2
     cvinvt=cv(l)/tempur(l)
     PRINT 497,tempur(l),ttt2,cv(l),cvinvt,cvr(l),tetad(l)
497  FORMAT(2F12.2,4E14.4)
  END DO
505 CONTINUE
  !     ------------------------------------------------------------------
  !     CALCULATION OF G(OM) FOR EQUIDISTANT TIME INTERVALS
  DO  l=1,lm5
     quo(l)=gges(l)-gnu(l)
     PRINT*,l,hoe(l),gnu0(l)
  END DO
  CALL equit(lm5,quo,quav/rsum)
  !     ------------------------------------------------------------------
  !     ANALYSIS OF INCOMPLETE ENERGY LOSS SPECTRUM
  IF(iloss == 1) CALL incomp(quav,gnu,gom1,rmai,dwi,gges,gnum,gnup)
1000 CONTINUE
  STOP
END PROGRAM


!***********************************************************************

SUBROUTINE fold(fin,fout,l,dele)
  !***********************************************************************


  REAL, INTENT(IN)                         :: fin(1024)
  REAL, INTENT(OUT)                        :: fout(1024)
  INTEGER, INTENT(IN)                      :: l
  REAL, INTENT(IN)                         :: dele

  DIMENSION eres(3000)

  COMMON/reso/eres,ifold

  DO j=1,l
     en=j*dele
     ien = 10*en
     res=eres(ien)
     IF(res == 0) GO TO 200
     res=res/dele
     m=3*res
     IF((j-m) < 1) THEN
        ml=1
     ELSE
        ml=j-m
     END IF
     IF((j+m) > l) THEN
        mu=l
     ELSE
        mu=j+m
     END IF
     fout(j)=0.
     fnorm=0.
     DO k=ml,mu
        fout(j)=fout(j)+fin(k)*EXP(-((j-k)/res)**2)
        fnorm=fnorm+EXP(-((j-k)/res)**2)
     END DO
     fout(j)=fout(j)/fnorm
200  CONTINUE
  END DO

END SUBROUTINE fold


!***********************************************************************

SUBROUTINE indata
  !***********************************************************************

  DIMENSION title(20)

  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv,  &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss
  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/multi/gnd(1024),lm,lm1,lm2,dele
  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff
  COMMON/par2/parq(1024),iquo,lpar
10 FORMAT('OUTPUT MUPHOCOR/5 *** EXTENDED VERSION  SEPT. 84      ')
15 FORMAT(20A4)
20 FORMAT(//,' PROBLEM:',20A4)
30 FORMAT(16I5)
40 FORMAT(8F10.4)
45 FORMAT(//,' **************************************************************')

  pi=3.1415927
  pir=pi/180.
  ndim=1024

  DO  n=1,ndim
     z0(n)=0.
     z1(n)=0.
     un(n)=0.
     parq(n)=0.
  END DO

  alfi(1)=1.
  alfi(2)=0.
  PRINT 10
  READ  15,(title(i),i=1,20)
  PRINT 20,(title(i),i=1,20)

  READ  40,e0,fp,cw,xnel,fimi,fima,abk,aemp
  !     E0:INC. ENERGY(MEV)             FP:FLIGHTPATH(CM)
  !     XNEL:CHANNEL OF ELASTIC LINE    CW:CHANNEL WIDTH(MYS)
  !     FIMI:MIN.SCATT.ANGLE            FIMA:MAX.SCATT.ANGLE
  !     UNT:CONST.BACKGROUND

  IF(aemp == 0.)aemp=5.60

  READ 30,lm,nspec,iquo,lpar
  !     NSPEC: NUMBER OF DIFFERENT ATOMIC SPECIES

  READ 40,hox,temp0,dw
  !     HOX:   UPPER LIMIT OF DENSITY OF STATES
  !     TEMPO: TEMPERATURE IN KELVIN
  !     DW:    MEAN DEBY WALLER COEFFICIENT

  dele=hox/lm
  !     DELE:  DELTA ENERGY FOR THE CALCULATION

  temp=temp0/11.605
  !     TEMP:  TEMPERATURE IN MEV

  t0=SQRT(522.724/e0)*fp
  !     TO:  TIME OF ARRIVAL OF ELASTICALLY SCATTERED NEUTRONS IN MYS

  t0r=t0/cw
  !     TOR: T0 IN MULTIPLES OF CHANNELS

  DO  i=1,nspec
     READ 40,amasi(i),sigi(i),conci(i)
  END DO

  !     AMASI: ATOMIC MASS
  !     SIGI : SIGMA
  !     CONC:  CONCENTRATION
  !     ALFI:  SCATTERING POWER

  conci(2)=1.-conci(1)
  IF(nspec == 2) READ 40,(parq(i),i=1,lpar)

  READ 30,  npho,itm,ivit,idw,iemp,ires,ipr,iloss

  !     NPHO:     NUMBER OF MULTI PHONON TERMS,ITM:NUMBER OF ITERATIONS
  !     ITM:      TOTAL NUMBER OF ITERATIONS
  !     IVIT=0(1):ITERATION BY DIFFERENCE(QUOTIENT) METHOD
  !     IDW=0:    DW COEFF. KEPT CONST.=INPUT VALUE,DW=1:DW COEFF.ITERATED
  !     IEMP=1(0):CORRECTIONS FOR COUNTER EFFICIENCY WILL BE(NOT BE) DONE
  !     IRES=1(0):DATA WILL BE(NOT BE) CORRECTED FOR SPECTR. RESOLUTION

  !     ------------------------------------------------------------------
  IF(iemp == 0) PRINT 110
110 FORMAT(//,'INPUT DISTRIBUTION HAS ALREADY BEEN CORRECTED FOR COUNTER EFFICIENCY')
  IF(iemp == 1) PRINT 115
115 FORMAT(//,'CORRECTION FOR COUNTER EFFICIENCY IS DONE IN THIS PROGRAM ')
  PRINT 130, e0,fp,cw,abk
130 FORMAT(/ ,'INC. ENERGY=',f10.5,'   FLIGHT PATH=',f10.5,'   CHANNEL WIDTH=',f10.5,'   ABS.COEFF.=',f10.5)
  PRINT 135, fimi,fima
135 FORMAT(/ ,'RANGE OF SCATTERING ANGLES:',f10.3,'TO ',f10.3)
  PRINT 140,hox
140 FORMAT(/,'UPPER LIMIT OF PHONON DENSITY OF STATES:',f10.4)
  PRINT 145,temp0,dw
145 FORMAT(/,' TEMPERATURE:',f10.3,'   MEAN DEBYE-WALLER COEFF.:',e14. 5)
  PRINT 150
150 FORMAT(' IN THIS PROGRAM VALUES FOR DW-COEFF. ARE IN (1/MEV), TO OBTAIN VALUES IN (ANG**2) MULTIPLY BY 2.0717 ')
  IF(idw == 0) PRINT 155
155 FORMAT(/,'DEBYE-WALLER COEF. FIXED=INPUT VALUE                ')
  IF(idw == 1) PRINT 160
160 FORMAT(/,'DEBYE-WALLER COEF. DETERM. FROM ITERATED DISTR.     ')
  IF(nspec == 1) GO TO 180
  PRINT 170
170 FORMAT(//,' MULTI-PHONON CONTRIBUTIONS CALCULATED SEPARATELY FOR BOTH COMPONENTS ')
  CALL const
180 CONTINUE
  DO  i=1,nspec
     PRINT 182,amasi(i),sigi(i),conci(i),alfi(i)
  END DO
182 FORMAT(/,' ATOMIC MASS:',f10.3,'   SIGMA:',f8.4,'   CONC:',f8.4,' SCATT.POWER:',f8.4)
  rma=1.0087/amasi(1)
  IF(nspec == 1) GO TO 200
  PRINT 195
195 FORMAT(/,' PARAMETERS FOR FNU/GNU:')
  PRINT 40,(parq(i),i=1,8)
  rma=1.0087/ameff
200 CONTINUE
  PRINT 205,npho
205 FORMAT(/,'NUMBER OF MULTI-PHONON-TERMS:',i3)
  PRINT 208,itm
208 FORMAT(/,' TOTAL NUMBER OF ITERATIONS :',i3)
  IF(ivit == 0) PRINT 210
210 FORMAT('     ITERATION BY DIFFERENCE METHOD')
  IF(ivit > 0) PRINT 215,ivit
215 FORMAT('     ITERATION BY QUOTIENT METHOD UP TO IT. NR. :',i3)

  CALL intof

  RETURN
END SUBROUTINE indata

!***********************************************************************
!***********************************************************************

SUBROUTINE intof
  !***********************************************************************
  !***********************************************************************

  DIMENSION z00(1024),un0(1024)

  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv,  &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss
  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/multi/gnd(1024),lm,lm1,lm2,dele
  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff
  COMMON/par2/parq(1024),iquo,lpar

30 FORMAT(16I5)
40 FORMAT(8F10.4)
45 FORMAT(//,' ******************************************************  &
       ****************************************************************')

  READ 30,isour,nuu,noo,nu,no,iglu
  !     NU:  First channel of spektrum NO: Last channel of spektrum
  !     NUU: First channel of TOF distribution used for calculation
  !     NOO: Last channel  of TOF distribution used for calculation

  nm=no-nu+1
  !     NM: Number of TOF input data

  READ 40,unt,fun
  !     UNT: CONSTANT BACKGROUND
  !     FUN: MULTIPLICATION FACTOR FOR TIME DEPENDENT BACKGROUND

  IF(isour /= 1) THEN
     READ 40,(z00(n),n=nuu,noo)    !     Z00: TIME OF FLIGHT DISTRIBUTION
     IF(fun /= 0.) READ 40,(un0(n),n=nuu,noo)     !     UN0: TIME DEPENDENT BACKGROUND
  END IF

  !     Shift data to the left
  DO  n=1,nm
     z0(n)=z00(nu+n-1)
     un(n)=un0(nu+n-1)
  END DO

  PRINT 160,xnel,nu,no
160 FORMAT(//,'CHAN.OF EL.LINE:',f10.2,'   FIRST CHANNEL:',i5,'   LAST CHANNEL:',i5)
  PRINT 235,unt
235 FORMAT(/,'CONSTANT BACKGROUND:',f10.3)
  PRINT 240,nuu,noo
240 FORMAT(//,' INPUT TIME OF FLIGH DISTRIBUTION FROM CHANNEL :',i4,'  &
       TO CHANNEL :',i4)
  PRINT 245, (z00(n),n=nuu,noo)
245 FORMAT(10F12.3)
  IF(fun == 0.)GO TO 295
  PRINT 260,fun
260 FORMAT(//,'TIME DEPENDENT BACKGROUND -- HAS TO BE MULTIPLIED BY  &
       A FACTOR:',f10.5)
  PRINT 245, (un0(n),n=nuu,noo)

  !     SMOOTHING OF TIME DEPENDENT BACKGROUND
  IF(iglu /= 0) THEN
     DO  i=1,iglu
        CALL glat0(un,z1,nm)
        DO  n=1,nm
           un(n)=z1(n)
           z1(n)=0.
        END DO
     END DO
     PRINT 280,iglu
280  FORMAT(//,'   TIME DEPENDENT BACKGROUND AFTER SMOOTHING -- NUMBER  &
          OF PROCESSES:',i2)
     PRINT 245, (un(n),n=1,nm)
  END IF


  !     Subtraction of background
  DO  n=1,nm
     z0(n)=z0(n)-fun*un(n)
  END DO

295 CONTINUE

  bac=unt
  DO  n=1,nm
     delk=nu-1+n-xnel
     !       Reminder: TOR IS THE TIME OF ARRIVAL OF ELASTICALLY SCATTERED
     !       NEUTRONS IN MULTIPLES OF CHANNELS
     !       HOT: Energy gain or loss respectively
     hot(n)=e0*(1.-(t0r/(t0r+delk))**2)
     !       Correction of constant backgroud for counter efficiency
     IF(iemp == 0) bac=unt/(1.0-EXP(-aemp/SQRT(e0-hot(n))))
     z0(n)=z0(n)-bac
  END DO

  !     Resolution corrections
  IF(ires /= 0) THEN
     PRINT 45
     PRINT 305
305  FORMAT(//,'INPUT DISTRIBUTION WILL BE CORRECTED FOR EXPERIMENTAL RESOLUTION')
     CALL rescor(z0,z1)

     DO  n=1,nm
        z0(n)=z1(n)
     END DO
  END IF

  PRINT 45
  ie0=-1
  vorz=.5*(nu+no)-xnel

  IF (vorz /= 0.0) THEN
     IF(vorz < 0.0) THEN      !     Energy gain spektrum
        ie0=+1 ! 320
        nmh=nm/2

        !     Exchange first with last number and so forth
        DO  n=1,nmh
           nn1=nm+1-n
           zint=z0(n)
           z0(n)=z0(nn1)
           z0(nn1)=zint
           zint=ABS(hot(n))
           hot(n)=ABS(hot(nn1))
           hot(nn1)=zint
   	END DO

   	hot(nmh+1)=ABS(hot(nmh+1))
   	PRINT 330
330 	FORMAT(//,50H Energy gain mode                                       )
     ELSE                      !     Energy loss spektrum

335 	PRINT 340
340 	FORMAT(//,50H Energy loss  mode                                      )
     END IF  ! if(vorz<0)
     emp=1.
     DO  n=1,nm
        ef=e0+ie0*hot(n)
        sqef=SQRT(ef)
        IF(iemp == 1) emp=1.-EXP(-aemp/sqef)
        z0(n)=z0(n)*EXP(5.025*abk/sqef)/emp
        un(n)=z0(n)
     END DO

     CALL mean(dele)

     PRINT 360
360  FORMAT(/50H Absolute energy transfer                             )
     PRINT 245, (hot(n),n=1,nm)
     PRINT 365
365  FORMAT(/,' INPUT DISTRIBUTION AFTER CORRECTIONS FOR BACKGROUND,RESOLUTION,COUNTER EFFICIENCY AND ABSORPTION  ')
     PRINT 245, (z0(n),n=1,nm)
  END IF ! 1000 CONTINUE



  RETURN
END SUBROUTINE intof


!***********************************************************************
!***********************************************************************

SUBROUTINE prep1(lm2,lm5,e0,dele)
  !***********************************************************************
  !***********************************************************************

  !     In subroutine PREP1 the temperature factor plus several other
  !     functions of Q or omega are calculated and stored in the
  !     respective data arrays


  INTEGER, INTENT(IN)                      :: lm2
  INTEGER, INTENT(IN)                      :: lm5
  REAL, INTENT(IN)                         :: e0
  REAL, INTENT(IN)                         :: dele
  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv,  &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss

  COMMON/bb2/hoe(1024),efin(1024),exe(1024),exa(1024),com,coma,comi,  &
       qc(1024),qma(1024),qmi(1024),dqup(1024,20)

  DO  l=1,lm2

     !       HOE = energy in [meV]
     en=(l-.5)*dele
     hoe(l)=en

     !       EX = e^(energy/2kT)
     ex=EXP(.5*en/temp)

     !       EXE = energy * sinh(energy/2kT)
     exe(l)=en*(ex-1./ex)

     !       EXA = energy * cosh(energy/2kT)
     exa(l)=ex+1./ex
  END DO

  !     Reminder: PIR = pi in radiant, FIMA and FIMI: maximum and minimum
  !               scattering angles

  fi0=.5*(fima+fimi)*pir
  dfi=.5*(fima-fimi)*pir

  com=COS(fi0)*COS(dfi)
  coma=COS(pir*fima)
  comi=COS(pir*fimi)

  DO  l=1,lm5
     !       EF = energy of scattered neutrons
     ef = e0+ie0*hoe(l)

     !       QC = h^2/2m * (ki^2 + kf^2 -2ki*kf*cos(theta)) =
     !       momentum transfer squared times h^2/2m, m denoting
     !       the neutron mass

     qc(l)  = e0+ef-2.*SQRT(e0*ef)*com
     qma(l) = e0+ef-2.*SQRT(e0*ef)*coma
     qmi(l) = e0+ef-2.*SQRT(e0*ef)*comi
     efin(l)= ef

     DO  n=1,npho
        n1=n+1
        dqup(l,n)=(qma(l)**n1-qmi(l)**n1)/(2.*n1)
     END DO
  END DO
  RETURN
END SUBROUTINE prep1




!***********************************************************************
!***********************************************************************

SUBROUTINE equit(lm5,quo,quav)
  !***********************************************************************
  !***********************************************************************


  INTEGER, INTENT(IN)                      :: lm5
  REAL, INTENT(IN)                         :: quo(1024)
  REAL, INTENT(IN)                         :: quav

  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/multi/gnd(1024),lm,lm1,lm2,dele
  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv,  &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss
  COMMON/bb2/hoe(1024),efin(1024),exe(1024),exa(1024),com,coma,comi,  &
       qc(1024),qma(1024),qmi(1024),dqup(1024,20)
45 FORMAT(//,' ******************************************************  &
       ****************************************************************')

  nma=nm                                                             !   last channel - first channel
  homax=hoe(lm5)

  IF(hot(nm) <= homax)GO TO 530                                      !   lm = 2*no of energy intervals
  DO  n=1,nm
     IF(hot(n) > homax)EXIT
  END DO
  nma=n-1
530 CONTINUE







  !loop537:  DO  n=1,nma
  !  DO  l=2,lm5
  !    dif=hot(n)-hoe(l)
  !    IF(dif <= 0.0) THEN
  !      GO TO   533
  !    END IF
  !    533     z1(n)=quo(l-1)+(quo(l)-quo(l-1))*(hot(n)-hoe(l-1))/dele
  !    CYCLE loop537
  !  END DO
  !END DO loop537

  ! This above replaced by:

  DO  n=1,nma
     DO  l=2,lm5
        dif=hot(n)-hoe(l)
        IF(dif <= 0.0) THEN
          z1(n)=quo(l-1)+(quo(l)-quo(l-1))*(hot(n)-hoe(l-1))/dele
          EXIT
        END IF
     END DO
  END DO



  zsum = 0.

  DO  n=1,nma

     ho=ie0*hot(n)
     ef=e0+ho
     ex=EXP(.5*ho/temp)
     afa=ho*(ex-1./ex)
     en=(comi-coma)*(e0+ef-2.*SQRT(e0*ef)*com)
     z0(n)=un(n)*afa*ex*e0/(en*ef*ef)

     PRINT*,hot(n),z0(n),z1(n),un(n)
     z1(n)=z0(n)/quav-z1(n)

     IF(z1(n) > 0.) GO TO 539
     z1(n)=0.
     !          print*,'scheiss, jetzt hammas'
539  CONTINUE
     !        print*,'ich bin bei 539 vorbei'

     IF((n > 1).AND.(n <= lm5/2)) zsum=zsum+(hot(n)-hot(n-1))*z1(n)

     qma(n)=hot(n)/4.1355
  END DO
  !      print*,'aus'
  !      stop

  DO  n=1,nma
     z1(n) = z1(n)/zsum
     qmi(n)= z1(n)*4.1355
  END DO


  PRINT 45
  PRINT 560
560 FORMAT(//,'PHONON DENSITY OF STATES FOR EQUIDISTANT TIME-INTERV  &
       ALS                                                             ')
  WRITE(6,565)
565 FORMAT(//,'            OMEGA    TOF DISTR.   INPUT DISTR.  &
       G(OM)             NU(THZ)     G(NU)                           ')
  WRITE(6,570)(n,hot(n),un(n),z0(n),z1(n),qma(n),qmi(n),n=1,nma)

  OPEN(11,FILE='DOS.dat')
  REWIND 11
  WRITE(11,570)(n,hot(n),un(n),z0(n),z1(n),qma(n),qmi(n),n=1,nma)
  CLOSE(11)
570 FORMAT(i5,f12.3,f12.2,e16.5,e18.5,f16.3,e16.5)
  RETURN
END SUBROUTINE equit

!***********************************************************************

SUBROUTINE incomp(quav,gnu,gom1,rmai,dwi,gges,gnum,gnup)
  !     THIS PART OF THE PROGRAM ALLOWS TO ANALYSE SPECTRA IN ENERGY LOSS
  !     IF ONE CANNOT DETERMINE THE TOTAL SPECTRUM DUE TO A TOO LOW PRIM.
  !     ENERGY
  !     THE E-GAIN SP. IS USED FOR NORMALISATION AND MULTI-PHONON CORR.

  REAL, INTENT(IN OUT)                     :: quav
  REAL, INTENT(OUT)                        :: gnu(1024)
  REAL, INTENT(IN)                         :: gom1(2,20,1024)
  REAL, INTENT(IN)                         :: rmai(2,20)
  REAL, INTENT(IN)                         :: dwi(2)
  REAL, INTENT(OUT)                        :: gges(1024)
  REAL, INTENT(OUT)                        :: gnum(1024)
  REAL, INTENT(OUT)                        :: gnup(1024)
  DIMENSION gnu0(1024), gesi(2,1024)

  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/multi/gnd(1024),lm,lm1,lm2,dele
  COMMON/bb1/pi,pir,fimi,fima,abk,aemp,npho,itm,idw,iemp,ires,ipr,iv,  &
       it,iplot,hox,temp,dw,unt,fun,nspec,rma,t0,t0r,ie0,iloss
  COMMON/bb2/hoe(1024),efin(1024),exe(1024),exa(1024),com,coma,comi,  &
       qc(1024),qma(1024),qmi(1024),dqup(1024,20)
  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff
30 FORMAT(16I5)
40 FORMAT(8F10.4)
45 FORMAT(//,' ******************************************************  &
       ****************************************************************')
  PRINT 45
  PRINT 60
60 FORMAT(//,50H analysis of the energy loss spectrum               )
  READ 30,ires

  CALL intof
  lm5=e0/dele
  PRINT 45
  DO  l=1,lm5
     en=(l-.5)*dele
     ef=e0-en
     qc(l)=e0+ef-2.*SQRT(e0*ef)*com
     qma(l)=e0+ef-2.*SQRT(e0*ef)*coma
     qmi(l)=e0+ef-2.*SQRT(e0*ef)*comi
     efin(l)=ef
     DO  n=1,npho
        n1=n+1
        dqup(l,n)=(qma(l)**n1-qmi(l)**n1)/(2.*n1)
     END DO
  END DO
  ngu=hot(1)/dele+.5
  ngo=hot(nm)/dele+.5
  ngu1=ngu+1
  ngo1=ngo+1
  ngoo=MIN0(lm5,ngo)

!
!   loop170:  DO  l=ngu1,ngoo
!      DO  n=1,nm
!         dif=hoe(l)-hot(n)
!         IF(dif <= 0.0) THEN
!            GO TO   160
!         END IF
! 160     z1(l)=z0(n-1)+(z0(n)-z0(n-1))*(hoe(l)-hot(n-1))/(hot(n)-hot(n-1))
!         CYCLE loop170
!      END DO
!   END DO loop170

  !  This above replaced by (same assomewhere above):
  DO  l=ngu1,ngoo
     DO  n=1,nm
        dif=hoe(l)-hot(n)
        IF(dif <= 0.0) THEN
           z1(l)=z0(n-1)+(z0(n)-z0(n-1))*(hoe(l)-hot(n-1))/(hot(n)-hot(n-1))
           EXIT
        END IF
     END DO
  END DO
  DO  l=ngo1,lm5
     z1(l)=0.
  END DO
  DO  l=ngu1,lm5
     ef=efin(l)
     ex=EXP(-.5*hoe(l)/temp)
     gnu0(l)=z1(l)*exe(l)*ex*(e0/ef)*SQRT(e0/ef)/(quav*dqup(l,1))
  END DO
  afa=0.
  DO  i=1,3
     afa=afa+gnu0(ngu+i)/(hoe(ngu+i))**2
  END DO
  afa=afa/3.
  DO  l=1,ngu
     gnu0(l)=afa*hoe(l)**2
     z1(l)=0.
  END DO
  DO  l=1,lm5
     ef=efin(l)
     hot(l)=t0r*(SQRT(e0/ef)-1.)
     gnum(l)=0.
  END DO
  DO  i=1,nspec
     dw=dwi(i)
     DO  l=1,lm5
        ex=EXP(-.5*hoe(l)/temp)
        gges(l)=0.
        qrma=dw*qma(l)
        qrmi=dw*qmi(l)
        exma=EXP(-qrma)
        exmi=EXP(-qrmi)
        cof0=exmi-exma
        DO  n=1,npho
           cof0=n*cof0+exmi*qrmi**n-exma*qrma**n
           IF(n == 1)cof=cof0
           gges(l)=gges(l)+cof0*rmai(i,n)*gom1(i,n,l)/(2.*dw**(n+1))
        END DO
        gges(l)=gges(l)/(dqup(l,1)*rmai(i,1))
        gesi(i,l)=gges(l)*exe(l)
     END DO
  END DO
  DO  l=1,lm5
     gges(l)=gesi(1,l)*alfi(1)+gesi(2,l)*alfi(2)
     !     GNU(L)=2.*DW**2*(GNU0(L)-GGES(L))*DQUP(L,1)/COF+GNU(L)
     gnup(l)=gnu0(l)/gges(l)
     gnu(l)=gnup(l)*gnu(l)
  END DO
  gnum(1)=gnu(1)*dele
  DO  l=2,lm5
     gnum(l)=gnum(l-1)+gnu(l)*dele
  END DO
  lm=hox/dele
  rn=1./gnum(lm)
  !     DO 240 L=1,LM5
  !     GNUP(L)=GNU(L)*RN
  ! 240 CONTINUE
  PRINT 260
260 FORMAT(//,'          OMEGA   INPUT DISTR. CALC.(E-GAIN)   INP/CALC  &
       GNU(L)       INT(GNU)      CHAN.SHIFT  ')
  PRINT 270,  (l,hoe(l),gnu0(l),gges(l),gnup(l),gnu(l),gnum(l),hot(l ),l=1,lm5)
270 FORMAT(i5,f10.3,2E14.5,f11.4,e16.5,f12.5,f12.2)
  RETURN
END SUBROUTINE incomp
!***********************************************************************

SUBROUTINE const
  !***********************************************************************

  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff

  !     Calculation of alpha(i)

  DO  i=1,2
     alfi(i)=conci(i)*sigi(i)/amasi(i)
  END DO

  alfn=alfi(1)+alfi(2)
  alfi(1)=alfi(1)/alfn
  alfi(2)=alfi(2)/alfn

  rnor=1./(alfi(1)*conci(2)-alfi(2)*conci(1))

  ac1=amasi(1)*conci(1)
  ac2=amasi(2)*conci(2)

  !     QFG0=S_0 : ratio of mean mass over effective mass
  qfg0=(ac1+ac2)/(alfi(1)*amasi(1)+alfi(2)*amasi(2))

  !     QFG1=S_1 : c(1) over alpha(1)
  qfg1=conci(1)/alfi(1)
  qfg2=conci(2)/alfi(2)

  !     AMEFF: effective mass
  ameff=(sigi(1)*conci(1)+sigi(2)*conci(2))/alfn

  qomav=(amasi(2)*conci(1)+amasi(1)*conci(2))/(amasi(2)*alfi(1)+amasi(1)*alfi(2))
  PRINT 50,qfg0,qfg1,qfg2,quomav
50 FORMAT(//,' LIMITS OF FNU/GNU(0,1,2):',3F12.4,'   RATIO OF <OM**2>  &
       :',f10.3)
  PRINT 60,ameff
60 FORMAT(/,' EFFECTIVE MASS:',f10.3)
  RETURN
END SUBROUTINE const


!***********************************************************************
!***********************************************************************

SUBROUTINE corrfa(lm,gnu,dele,quofg,nspec)
  !***********************************************************************
  !***********************************************************************


  INTEGER, INTENT(IN)                      :: lm
  REAL, INTENT(IN)                         :: gnu(1024)
  REAL, INTENT(IN)                         :: dele
  REAL, INTENT(OUT)                        :: quofg(1024)
  INTEGER, INTENT(IN OUT)                  :: nspec

  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff
  COMMON/par2/parq(1024),iquo,lpar

  lm5=2*lm
  IF(nspec == 2) GO TO 55

  DO  l=1,lm5
     quofg(l)=1.
  END DO

  !     Goto end
  GO TO 700

55 CONTINUE

  !     If IQUO is unequal 0 an external model for the ratio of F/G has
  !     to be provided. The model is given through the array PARQ
  !
  !  IF(iquo == 0) GO TO 70
  !   DO  l=1,lpar
  !     quofg(l)=parq(l)
  !   END DO
  !
  !   GO TO 700
  !   70 CONTINUE

  ! ----------------------------------------
  !  This above replaced by
  ! ----------------------------------------
  IF (iquo /= 0) THEN
     DO l=1,lpar
        quofg(l) = parq(l)
     END DO
     GO TO 700  ! End of corrfa
  END IF


  !     LG is the cut-off energy in units of DELE
  lg=parq(1)/dele+1
  xlg=FLOAT(lg)
  IF(parq(2) /= 0.) qfg1=parq(2)

  fint=0.
  fint1=0.
  fint2=0.
  fint4=0.

  DO  l=1,lg
     xlr=FLOAT(l)/xlg
     xl2=(l-.5)
     xl2=xl2*xl2

     !       QUOFG(L) = S_0 +(S_1 - S_0) * L^2
     quofg(l)=qfg0+(qfg1-qfg0)*xlr*xlr

     fint=fint+gnu(l)
     fint1=fint1+gnu(l)*quofg(l)
     fint2=fint2+gnu(l)*xl2
     fint4=fint4+gnu(l)*xl2*xl2
  END DO

  !     FINT: Integrated G(omega)
  fint=fint*dele

  !     FINT1: Integrated S(omega)*G(omega)
  fint1=fint1*dele

  !     QUO0: S_@
  quoo=(1.-fint1)/(1.-fint)

  !     OMAV2: second moment of G(omega)
  omav2=0.
  DO  l=1,lm
     xl2=(l-.5)
     xl2=xl2*xl2
     omav2=omav2+xl2*gnu(l)
  END DO

  lg1=lg+1
  DO  l=lg1,lm5
     quofg(l)=quoo
  END DO

700 CONTINUE
  RETURN
END SUBROUTINE corrfa


!***********************************************************************
!***********************************************************************

SUBROUTINE separ(lm,gnu,dele,quofg,fpi,nspec)
  !***********************************************************************
  !***********************************************************************

  INTEGER, INTENT(IN)                      :: lm
  REAL, INTENT(IN)                         :: gnu(1024)
  REAL, INTENT(IN)                         :: dele
  REAL, INTENT(IN)                         :: quofg(1024)
  REAL, INTENT(OUT)                        :: fpi(2,1024)
  INTEGER, INTENT(IN OUT)                  :: nspec

  COMMON/par1/amasi(2),sigi(2),conci(2),alfi(2),rnor,qfg0,qfg1,qfg2, ameff

  lm5=2*lm
  IF(nspec == 2) GO TO 55
  DO  l=1,lm5
     fpi(1,l)=gnu(l)
     fpi(2,l)=0.
  END DO

  GO TO 150

55 CONTINUE
  fint=0.
  DO  l=1,lm
     fint=fint+gnu(l)*quofg(l)
  END DO

  !     FINT: Integrated S(omega)*G(omega)
  fint=fint*dele

  !     FPI: partial densities of state
  DO  l=1,lm5
     fpi(1,l)= rnor*(conci(2)-alfi(2)*quofg(l)/fint)*gnu(l)
     fpi(2,l)=-rnor*(conci(1)-alfi(1)*quofg(l)/fint)*gnu(l)
  END DO

  sum1=0.
  sum2=0.
  DO  l=1,lm
     sum1=sum1+fpi(1,l)
     sum2=sum2+fpi(2,l)
  END DO
  sum1=sum1*dele
  sum2=sum2*dele
  PRINT 110,fint,sum1,sum2
110 FORMAT(//,' INTEGRALS (FNU,FNU1,FNU2):',3F12.5)
150 CONTINUE
  RETURN
END SUBROUTINE separ


!***********************************************************************
!***********************************************************************

SUBROUTINE vipho(gnum,gnup,ncase)
  !***********************************************************************
  !***********************************************************************

  REAL, INTENT(IN)                         :: gnum(1024)
  REAL, INTENT(OUT)                        :: gnup(1024)
  INTEGER, INTENT(IN OUT)                  :: ncase


  COMMON/multi/gnd(1024),lm,lm1,lm2,dele

  !     Folding of GNUM with GND. Output is stored in GNUP.
  !     Attention: The folding takes place in the range from 0 to 2 E_max

  !     NCASE = 1 for first run

  IF(ncase == 1)GO TO 180

  DO  i=1,lm
     gnup(i)=0.
     DO  k=1,lm2
        k1=lm2+1-k
        k2=i+k-1-lm
        IF(k2 > 0) THEN
           GO TO    95
        END IF
94      k2=-k2+1
95      gnup(i)=gnup(i)+gnd(k1)*gnum(k2)
     END DO
     gnup(i)=gnup(i)*dele
  END DO

  DO  i=lm1,lm2
     gnup(i)=0.
     lo=3*lm-i+1
     DO  k=1,lo
        gnup(i)=gnup(i)+gnd(lm2+1-k)*gnum(i+k-1-lm)
     END DO
     gnup(i)=gnup(i)*dele
  END DO

  GO TO 200

180 CONTINUE
  DO  i=1,lm2
     gnup(i)=0.
     lo=lm2-i+1
     DO  k=1,lo
        gnup(i)=gnup(i)+gnd(lm2+1-k)*gnd(k+i-1)
     END DO
     gnup(i)=gnup(i)*dele
  END DO

200 RETURN
END SUBROUTINE vipho

!***********************************************************************
!***********************************************************************

SUBROUTINE mean(dele)
  !***********************************************************************
  !***********************************************************************

  !     Subroutine MEAN averages the time spectrum in the region where
  !     there are several time channels per corresponding energy channel



  REAL, INTENT(IN)                         :: dele
  COMMON/avar/z0(1024),z1(1024),hot(1024),nm,ndum,iglu,un(1024)

  z0(nm+1)=z0(nm)
  z0(nm+2)=z0(nm)
  z0(nm+3)=z0(nm)
  DO  n=4,nm
     rdh=(hot(n)-hot(n-1))/dele
     rdh=1./ABS(rdh)
     fa1=0.
     fa2=0.
     fa3=0.
     IF(rdh <= 1.)GO TO 108
     fa1=.5*(rdh-1.)
     IF(rdh <= 3.)GO TO 108
     fa1=1.
     fa2=.5*(rdh-3.)
     IF(rdh <= 5.)GO TO 108
     fa2=1.
     IF(n >= 4)fa3=1.
108  CONTINUE
     fnorm=1./(1.+2.*(fa1+fa2+fa3))
     z1(n)=fnorm*(z0(n)+fa1*(z0(n-1)+z0(n+1))+fa2*(z0(n-2)+z0(n+2))+  &
          fa3*(z0(n-3)+z0(n+3)))
  END DO
  DO  n=4,nm
     z0(n)=z1(n)
  END DO
  RETURN
END SUBROUTINE mean

!***********************************************************************
!***********************************************************************

SUBROUTINE moment(gnu,ho,lmax,dho,cmo)
  !***********************************************************************
  !***********************************************************************



  REAL, INTENT(IN OUT)                     :: gnu(1024)
  REAL, INTENT(IN)                         :: ho(1024)
  INTEGER, INTENT(IN)                      :: lmax
  REAL, INTENT(IN)                         :: dho
  REAL, INTENT(OUT)                        :: cmo(10)

  sum=0.
  DO  l=1,lmax
     sum=sum+gnu(l)
  END DO
  sum=sum*dho
  DO  l=1,lmax
     gnu(l)=gnu(l)/sum
  END DO
  DO  m=1,10
     mex=m-3
     cmo(m)=0.
     DO  l=1,lmax
        cmo(m)=cmo(m)+gnu(l)*ho(l)**mex
     END DO
     cmo(m)=cmo(m)*dho
  END DO
  DO  m=1,10
     cmo(m)=(cmo(m)*m)/3.
     ex=m-3.
     IF(m == 3)ex=1.
     ex=1./ex
     cmo(m)=cmo(m)**ex
  END DO
  cmo(3)=0.
  DO  l=1,lmax
     cmo(3)=cmo(3)+gnu(l)*ALOG(ho(l))
  END DO
  cmo(3)=EXP(.3333333+cmo(3)*dho)
  RETURN
END SUBROUTINE moment
!     ******************************************************************

SUBROUTINE heat(lc,w,dva,cta,t,cv,tetad)


  INTEGER, INTENT(IN)                      :: lc
  REAL, INTENT(IN)                         :: w(1024)
  REAL, INTENT(IN)                         :: dva
  REAL, INTENT(IN)                         :: cta(1024)
  REAL, INTENT(IN)                         :: t
  REAL, INTENT(OUT)                        :: cv
  REAL, INTENT(OUT)                        :: tetad

  !     FREQUENCY SCALE IN MEV
  cv=0.
  DO  l=1,lc
     wr=11.605*w(l)/t
     IF(wr > 100.) GO TO 30
     ex=EXP(wr)
     IF(wr > 50.) GO TO 26
     cv=cv+ex*wr*wr*cta(l)/((ex-1.)**2)
     CYCLE
26   cv=cv+wr*wr*cta(l)/ex
30   CONTINUE
  END DO
  cvr=cv*dva
  cv=cvr*24.9429
  IF(cvr > .25) GO TO 50
  cvri=cvr**(-.33333333)
  tetad=4.27133*t*cvri*(1.-9.679*EXP(-2.9576*cvri))
  GO TO 100
50 CONTINUE
  IF(cvr > .6) GO TO 55
  dcvr=cvr-.26982
  x2=84.+158.238*dcvr+323.91*dcvr**3
  tetad=500.*t/x2
  GO TO 100
55 CONTINUE
  dcvr=20.*(1.-cvr)
  x2=dcvr*(1.+.03571429*dcvr+1.448728E-3*dcvr**2+6.24974E-5*dcvr**3  &
       +3.8E-6*dcvr**4)
  tetad=t*SQRT(x2)
100 CONTINUE
  RETURN
END SUBROUTINE heat

!***********************************************************************
!***********************************************************************

SUBROUTINE rescor(z0,z00)
  !***********************************************************************
  !***********************************************************************


  REAL, INTENT(IN)                         :: z0(1024)
  REAL, INTENT(OUT)                        :: z00(1024)
  COMMON/ros/nu,no,e0,fp,cw,xnel
  COMMON/reso/eres,ifold
  DIMENSION ho(1024),quo1(1024),quo2(1024),ERR(1024),res(100)
  DIMENSION eres(3000),den(1024)
  DIMENSION  z1(1024),z2(1024),tau(1024), z01(1024)
  DIMENSION rimp(100),nit(1024),den2(1024)

  xex=2.3
  !     DELL = ratio of the distances Fermi-sample over sample detector
  dell=70./300.

  !     NUS and NOS define the range of time channels where the correction
  !     is to be applied.
  !     IGL: number of smoothing processes
  !     ITM: number of iterations for the partial deconvolution

  READ(5,30)nus,nos,igl,itm
30 FORMAT(4I5)

  itm1=itm+1

  !     TAUE: primary resolution; TAUT: secondary resolution

  READ(5,35)taue,taut,dfp,unt,ferr
  IF(ferr == 0.) ferr=1.
35 FORMAT(8F10.3)

  nm=no-nu+1
  nms=nos-nus+1
  nud=nus-nu

  WRITE(6,37)taue,taut,unt,ferr
37 FORMAT(/,' DELTA-T(PRIM.)=',f8.2,'   DELTA-T(SEC.)=',f8.2,'   EXP.  &
       BACKGROUND=',f10.2,'   FACTOR FOR DELTA-T=',f9.3)

  !     T0: time of arrival of elastically scattered neutrons, flight path
  !     FP in cm
  t0=SQRT(522.724/e0)*fp
  t0=t0/cw

  DO  n=1,nm
     !       DELK: time channel with respect to elastic line
     delk=nu-1+n-xnel
     !       RT3: conversion factor from delta(w) to delta(t)
     !C      Reichardt:
     !C      RT3=(1.+DELK/T0)**3+DELL
     !C      TAU(N)=SQRT(TAUT**2+(TAUE*RT3)**2)
     rt3=(1.+delk/t0)**3
     tau(n)=SQRT(taut**2+(taue*(rt3+dell))**2)
     den(n)=2. * (e0/(t0*cw)) * (tau(n)/rt3)
     !       H0: hw
     ho(n)=e0*(1.-(t0/(t0+delk))**2)

     ihonew = ABS(10*ho(n))
     IF(n == 1) THEN
        ihoold=ihonew
        eres(ihoold)=den(n)
     ELSE
        IF(ihonew > ihoold) THEN
           f = ABS(ihonew-ihoold)
           DO i=1,(ihonew-ihoold)
              eres(ihoold+i)=eres(ihoold) + i/f*(den(n)-eres(ihoold))
           END DO
           ihoold=ihonew
        END IF
        IF(ihonew < ihoold) THEN
           f = ABS(ihonew-ihoold)
           DO i=1,ihoold-ihonew
              eres(ihoold-i)=eres(ihoold) + i/f*(den(n)-eres(ihoold))
           END DO
           ihoold=ihonew
        END IF
     END IF


     zab=ABS(z0(n))
     !       ZO: measured spectrum
     z01(n)=z0(n)
     IF(zab == 0.)zab=1.
     !       ERR: statistical error of measured spectrum
     ERR(n)=ferr*SQRT(zab+2.*unt)/zab
  END DO

  OPEN(13,FILE=  &
       'resolution.plot')
  REWIND 13

  WRITE(13,41)(n, ho(n),den(n),n=1,nm)
41 FORMAT(i5,2F9.3)
  CLOSE(13)

  IF(itm == 0) THEN
     DO n=1,nm
        z00(n)=z0(n)
     END DO
     ifold=1
     RETURN
  END IF

  WRITE(6,42)igl
42 FORMAT(//,'NUMBER OF SMOOTHING PROCESSES:',i3)

  WRITE(6,43)itm
43 FORMAT(/,'MAXIMUM NUMBER OF ITERATIONS:',i3)

  !     *******RIMP(n): Rel. resolution obtained after n iterations *******

  !     RIMP(n)=sqrt((RIMP(n-1)^4+RIMP(n-1)^8)/(1+RIMP(n-1)^4+RIMP(n-1)^8))
  rimp(1)=1.
  DO  n=2,itm1
     rim=rimp(n-1)
     rimp(n)=rim*(1.+rim)/(1.+rim+rim*rim)
  END DO
  DO  n=1,itm1
     rimp(n)=SQRT(rimp(n))
  END DO
  WRITE(6,314)
314 FORMAT(//,50H resolution tau(n)/tau(0) n=0,1,2,--,itm            )
  WRITE(6,315)(rimp(n),n=1,itm1)
315 FORMAT(10F12.4)

  !     ******************************************************************

  !     ********* Smoothing of input spectrum ****************************
  IF(igl == 0)GO TO 58
  DO  n=1,nms
     z00(n)=z0(nud+n)
  END DO
  DO  i=1,igl
     CALL glat0(z00,z01,nms)
     DO  n=1,nms
        z00(n)=z01(n)
     END DO
  END DO

  DO  n=1,nm
     z01(n)=z0(n)
  END DO
  DO  n=1,nms
     z01(n+nud)=z00(n)
  END DO

  WRITE(6,53)nus,nos
53 FORMAT(/,'SMOOTHING FROM CHANNEL NUMBER',i5,'TO ',i5)
58 CONTINUE

  !     ******************************************************************

  !     ***** Special treatment of the limits of the spectrum ************

  DO  n=1,nm
     nug=tau(n)/cw
     IF(nug-n < 0) THEN
        GO TO    61
     END IF
  END DO

61 nug=nug+1
  DO  n=1,nm
     nog=tau(nm+1-n)/cw
     IF(nog-n <= 0) THEN
        GO TO    63
     END IF
  END DO

63 nog=nm-nog-1

  WRITE(6,72)nm,nug,nog
72 FORMAT(//,3I10)
110 FORMAT(10F12.1)
  nug1=nug-1

  DO  n=1,nug1
     tau(n)=(n-1)*cw+.00001
  END DO

  !     ******************************************************************

  DO  n=1,nm
     nit(n)=0
     z1(n)=z01(n)
     z2(n)=0.
     quo1(n)=0.
     quo2(n)=0.
  END DO
  !     ------------------------------------------------------------------
502 FORMAT(10F12.2)

  !     *********** beginning of iterations ******************************

  it=0
82 it=it+1
  idev=0
  DO  n=2,nog

     !       RES: resolution function normalized to 1
     jm=tau(n)/cw
     DO  j=1,jm
        ex=(2.*j*cw/tau(n))**xex
        res(j)=EXP(-.69315*ex)
     END DO
     sur=1.
     DO  j=1,jm
        sur=sur+2.*res(j)
     END DO
     rsur=1./sur
     DO  j=1,jm
        res(j)=res(j)*rsur
     END DO

     !       ***** folding of Z1 and RES **********
     z2(n)=z1(n)*rsur
     DO  j=1,jm
        z2(n)=z2(n)+(z1(n-j)+z1(n+j))*res(j)
     END DO
     !       **************************************
     quo1(n)=z2(n)/z01(n)
     dev=ABS(quo1(n)-1.)
     !       Abbruchbedingung eq. (4.2)
     IF(dev > ERR(n))GO TO 98
     idev=idev+1
     GO TO 99
98   CONTINUE
     !       *** new partially deconvoluted spectrum **
     z1(n)=z1(n)/quo1(n)
     !       ******************************************
     nit(n)=nit(n)+1
99   CONTINUE
     quo2(n)=z1(n)/z01(n)
  END DO
  IF(it-itm < 0) THEN
     GO TO    82
  END IF

  !     ************** End of iteration loop *****************************

200 CONTINUE

  DO  n=1,nm
     z00(n)=z0(n)+z1(n)-z01(n)
     nl=nit(n)+1
     den2(n)=den(n)*rimp(nl)
  END DO

  WRITE(6,145)
145 FORMAT(//,'   N    OMEGA    TAU    DEL-OM   Z-IMP   Z-IMP-SM  &
       Z-CALC   III/III   LIMIT     Z-CORR-SM  CORR.F.   Z-CORR  N-IT  D  &
       EL-OM                            ')
  WRITE(6,146)
146 FORMAT(//,'                                  I         II       I  &
       II                             IV      IV/II                     ')
  WRITE(6,150)(n,ho(n),tau(n),den(n),z0(n),z01(n),z2(n),quo1(n),ERR(  &
       n),z1(n),quo2(n),z00(n),nit(n),den2(n),n=1,nm)
150 FORMAT(i5,3F8.2,     3F10.0,2F9.3,f13.0,f9.3,f10.0,i5,f8.2)
  RETURN
END SUBROUTINE rescor


!***********************************************************************
!***********************************************************************

SUBROUTINE glat0(f0,f1,MAX)
!***********************************************************************
!***********************************************************************




  REAL, INTENT(IN)                         :: f0(1024)
  REAL, INTENT(OUT)                        :: f1(1024)
  INTEGER, INTENT(IN)                      :: MAX

  MAX1=MAX-1
  max2=MAX-2
  DO  m=3,max2
     f1(m)=-6.*f0(m-2)+24.*f0(m-1)+34.*f0(m)+24.*f0(m+1)-6.*f0(m+2)
     f1(m)=f1(m)/70.
  END DO
  f1(1)=62.*f0(1)+18.*f0(2)-6.*f0(3)-10.*f0(4)+6.*f0(5)
  f1(2)=18.*f0(1)+26.*f0(2)+24.*f0(3)+12.*f0(4)-10.*f0(5)
  f1(MAX1)=-10.*f0(MAX-4)+12.*f0(MAX-3)+24.*f0(max2)+26.*f0(MAX1)+18.*f0(MAX)
  f1(MAX)=6.*f0(MAX-4)-10.*f0(MAX-3)-6.*f0(max2)+18.*f0(MAX1)+62.*f0 (MAX)
  f1(1)=f1(1)/70.
  f1(2)=f1(2)/70.
  f1(MAX1)=f1(MAX1)/70.
  f1(MAX)=f1(MAX)/70.
  RETURN
END SUBROUTINE glat0
!     *****************************************************************

SUBROUTINE splot(nu,no,xax,y1,y1m,y2,y2m)

  INTEGER, INTENT(IN)                      :: nu
  INTEGER, INTENT(IN)                      :: no
  REAL, INTENT(IN OUT)                     :: xax(1024)
  REAL, INTENT(IN)                         :: y1(1024)
  REAL, INTENT(OUT)                        :: y1m
  REAL, INTENT(IN)                         :: y2(1024)
  REAL, INTENT(IN)                         :: y2m
  DIMENSION  ysc1(6),ysc2(6)
  !     CHARACTER*121 POS
  CHARACTER (LEN=1) :: pos(121)*1

  y1m=0.
  DO  i=nu,no
     IF(y1(i) > y1m) y1m=y1(i)
  END DO
  DO  i=1,6
     ysc1(i)=.2*y1m*i
     ysc2(i)=.2*y2m*i
  END DO
  PRINT 75
75 FORMAT(//,' G(OM)  (EQUIDISTANT ENERGY INTERVALS)' )
  PRINT 80,(ysc1(i),i=1,6)
80 FORMAT(  22X,f8.3,5(12X,f8.3))
  DO  i=1,121
     pos(i)=' '
  END DO
  DO  i=1,120
     pos(i)='-'
  END DO
  DO  i=1,7
     pos(20*(i-1)+1)='I'
  END DO
  PRINT 95,pos
95 FORMAT( 9X,121A1)
  DO  i=nu,no
     DO  j=1,121
        pos(j)=' '
     END DO
     pos(1)='I'
     DO  j=1,6
        pos(20*j+1)=']'
     END DO
     !     POS(1:1)='I'
     ipos=100*y2(i)/y2m+.5
     IF(ipos > 121) ipos=121
     pos(ipos+1)='*'
     ipos=100*y1(i)/y1m+.5
     IF (ipos < 0) ipos=0
     pos(ipos+1)='+'
     !     POS(IPOS+1:IPOS+1)='*'
     PRINT 190,xax(i),pos
190  FORMAT(f8.2,1X,121A1)
     ! 190 FORMAT(A121)
  END DO
  DO  i=1,120
     pos(i)='-'
  END DO
  DO  i=1,7
     pos(20*(i-1)+1)='I'
  END DO
  PRINT 95,pos
  PRINT 80,(ysc2(i),i=1,6)
  PRINT 215
215 FORMAT(' INT ( G(OM) )    ' )
  RETURN
END SUBROUTINE splot
