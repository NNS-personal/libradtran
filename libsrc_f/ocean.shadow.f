c      IMPLICIT REAL*8 (A-B,D-H,O-Z), COMPLEX*16 (C)
c      REAL*8 R(4,4)
c      CN1=(1D0,0D0)
c      CN2=(1.33D0, -0.0D0)
C    5 READ (6,*) DMUI,DMUR,PHII,PHIR
c      DMUI=0.5
c      DMUR=0.5
c      PHII=0
c      PHIR=0
c      PHII=PHII*DACOS(-1D0)/180D0
c      PHIR=PHIR*DACOS(-1D0)/180D0
c      SIGMA2=0.04D0
c      CALL RMATR (CN1,CN2,SIGMA2,DMUI,PHII,DMUR,PHIR,R)
c      DO 1 I=1,4
c         R(I,1)=R(I,1)*DMUI
c         R(I,2)=R(I,2)*DMUI
c         R(I,3)=R(I,3)*DMUI
c         R(I,4)=R(I,4)*DMUI
c         WRITE (6,1000) R(I,1),R(I,2),R(I,3),R(I,4)
c    1 CONTINUE
C      GO TO 5
c1000  FORMAT(4D14.6)
c      STOP
c      END
 
C**************************************************************
 
C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR
C   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS INCLUDED.

C   FOR ALL FORMULAS AND DEFINITIONS, SEE THE PAPER
C   M. I. Mishchenko and L. D. Travis, Satellite retrieval
C   of aerosol properties over the ocean using polarization as well as
C   intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-
C   17013 (1997).

C   INPUT INFORMATION:
 
C   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER)
C   DMUI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   PHII = INCIDENT AZIMUTH ANGLE (IN RADIANS)
C   DMUR = ABS(COSINE OF THE REFLECTION ZENITH ANGLE)
C   PHIR = REFLECTION AZIMUTH ANGLE IN RADIANS

C   OUTPUT PARAMETERS:

C   R - (4X4) REFLECTION MATRIX [R(1,1)=REFLECTION FUNCTION]
 
 
      SUBROUTINE RMATR (CN1,CN2,SIGMA2,DMUI,PHII,DMUR,PHIR,R)
c      SUBROUTINE RMATR (CN1, CN2, SIGMA2,DMUI,PHII,DMUR,PHIR)
      IMPLICIT REAL*8 (A-B,D-H,O-Z), COMPLEX*16 (C)
      REAL*8 R(4,4)
 
C   CARTEZIAN COMPONENTS OF THE UNIT VECTORS OF THE INCIDENT AND
C   SCATTERED BEAMS
 
 
c      WRITE (*,*) "Test", SIGMA2, CN1, CN2, DMUI, PHII, DMUR, PHIR, R

c      C1=(REC1, IMC1)
c      C2=(REC2, IMC2)
c      CN2=(1.33D0, -0.0D0)
      IF (DABS(DMUI-1D0).LT.1D-10) DMUI=0.999999999999D0
      IF (DABS(DMUR-1D0).LT.1D-10) DMUR=0.999999999999D0
      DCOSI=DCOS(PHII)
      DSINI=DSIN(PHII)
      DCOSR=DCOS(PHIR)
      DSINR=DSIN(PHIR)
      DSI=DSQRT(1D0-DMUI*DMUI)
      DSR=DSQRT(1D0-DMUR*DMUR)
      VI1=DSI*DCOSI
      VI2=DSI*DSINI
      VI3=-DMUI
      VR1=DSR*DCOSR
      VR2=DSR*DSINR
      VR3=DMUR
 
C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION
 
      UNIT1=VI1-VR1
      UNIT2=VI2-VR2
      UNIT3=VI3-VR3
      FACT1=UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR=DSQRT(1D0/FACT1)
 
C    FRESNEL REFLECTION COEFFICIENTS
 
      XI1=FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2=1D0 - (1D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2=CDSQRT(CXI2)
      C1=CN1*XI1
      C2=CN2*CXI2
      CRPER=(C1-C2)/(C1+C2)
      C1=CN2*XI1
      C2=CN1*CXI2
      CRPAR=(C1-C2)/(C1+C2)
 
C    CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
 
      TI1=-DMUI*DCOSI
      TI2=-DMUI*DSINI
      TI3=-DSI
 
      TR1=DMUR*DCOSR
      TR2=DMUR*DSINR
      TR3=-DSR
 
      PI1=-DSINI
      PI2=DCOSI
      PI3=0D0
 
      PR1=-DSINR
      PR2=DCOSR
      PR3=0D0
 
      PIKR=PI1*VR1+PI2*VR2+PI3*VR3
      PRKI=PR1*VI1+PR2*VI2+PR3*VI3
      TIKR=TI1*VR1+TI2*VR2+TI3*VR3
      TRKI=TR1*VI1+TR2*VI2+TR3*VI3
 
      E1=PIKR*PRKI
      E2=TIKR*TRKI
      E3=TIKR*PRKI
      E4=PIKR*TRKI
 
      CF11=E1*CRPER+E2*CRPAR
      CF12=-E3*CRPER+E4*CRPAR
      CF21=-E4*CRPER+E3*CRPAR
      CF22=E2*CRPER+E1*CRPAR
 
 1000 FORMAT(4D15.6)
 
C   CALCULATION OF THE STOKES REFLECTION MATRIX
 
      VP1=VI2*VR3-VI3*VR2
      VP2=VI3*VR1-VI1*VR3
      VP3=VI1*VR2-VI2*VR1
      DMOD=VP1*VP1+VP2*VP2+VP3*VP3
      DMOD=DMOD*DMOD
 
      RDZ2=UNIT3*UNIT3
      RDZ4=RDZ2*RDZ2
 
      DCOEFF=1D0/(4D0*DMUI*DMUR*DMOD*RDZ4*2D0*SIGMA2)
      DEX= -(UNIT1*UNIT1 + UNIT2*UNIT2)/(2D0*SIGMA2*RDZ2)
      DEX=DEXP(DEX)
      DCOEFF=DCOEFF*FACT1*FACT1*DEX
C     CE: To obtain pure Fresnel reflection, set DCOEFF=1
C      DCOEFF=1D0
      
      AF=0.5D0*DCOEFF
      AF11=CDABS(CF11)
      AF12=CDABS(CF12)
      AF21=CDABS(CF21)
      AF22=CDABS(CF22)
      AF11=AF11*AF11
      AF12=AF12*AF12
      AF21=AF21*AF21
      AF22=AF22*AF22
 
      R(1,1)=(AF11+AF12+AF21+AF22)*AF
      R(1,2)=(AF11-AF12+AF21-AF22)*AF
      R(2,1)=(AF11-AF22+AF12-AF21)*AF
      R(2,2)=(AF11-AF12-AF21+AF22)*AF
 
      CI=(0D0, -1D0)
 
      C21=DCONJG(CF21)
      C22=DCONJG(CF22)
      CTTTP=CF11*DCONJG(CF12)
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22
 
      R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
      R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
      R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
      R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
      R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
      R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
      R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
      R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
      R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
      R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
      R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
      R(4,4)=    ( CTTPP-CTPPT)*DCOEFF
 
C  SHADOWING
      P=DACOS(-1D0)
      S1=DSQRT(2D0*SIGMA2/P)
      S3=1D0/(DSQRT(2D0*SIGMA2))
      S2=S3*S3
      XI=DMUI
      XXI=XI*XI
      DCOT=XI/DSQRT(1D0-XXI)
      T1=DEXP(-DCOT*DCOT*S2)
      T2=DERFC(DCOT*S3)
      SHADOWI=0.5D0*(S1*T1/DCOT-T2)
      XI=DMUR
      XXI=XI*XI
      DCOT=XI/DSQRT(1D0-XXI)
      T1=DEXP(-DCOT*DCOT*S2)
      T2=DERFC(DCOT*S3)
      SHADOWR=0.5D0*(S1*T1/DCOT-T2)
      SHADOW=1D0/(1D0+SHADOWI+SHADOWR)
c CE Uncomment this to remove shadowing and have pure Cox and Munk
c      SHADOW=1D0
      DO 1 I=1,4
         DO 1 J=1,4
            R(I,J)=R(I,J)*SHADOW

 1    CONTINUE
      RETURN
      END