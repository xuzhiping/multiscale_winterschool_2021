
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 N_ELEM=3,NSTVTO=2,NSTVTT=14,NSTV=18)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
       INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY

       REAL*8 AINTW(4),XII(4,2),XI(2),dNdxi(NNODE,2),
     1 VJACOB(2,2),dNdx(NNODE,2),VJABOBINV(2,2),AN(4),BP(2,NDOFEL),
     2 DP(2),SDV(NSTV),BB(3,NDOFEL),CMAT(3,3),EPS(3),STRESS(3),
     3 VNI(2,NDOFEL),ULOC(2),PHASENOD(NNODE)
     
       REAL*8 DTM,THCK,HIST,CLPAR,GCPAR,EMOD,ENU,PARK,ENG
 
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
C
   ******************************************************************
C     ==================================================================
       IF ((JTYPE.EQ.ONE).OR.(JTYPE.EQ.TWO)) THEN

       TIMEZ=USRVAR(JELEM,17,1)
       IF (TIMEZ.LT.TIME(2)) THEN
        USRVAR(JELEM,17,1)=TIME(2)
        USRVAR(JELEM,18,1)=ZERO
       ELSE
        USRVAR(JELEM,18,1)=USRVAR(JELEM,18,1)+ONE
       ENDIF
       STEPITER=USRVAR(JELEM,18,1)

       CLPAR=PROPS(1)
       GCPAR =PROPS(2)
       THCK = PROPS(3)

       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO

       IF (JTYPE.EQ.ONE) THEN
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
       ELSEIF (JTYPE.EQ.TWO) THEN
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF

       DO INPT=1,INNODE
C     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTO
          SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
C
C     Local coordinates of the integration point
        XI(1) = XII(INPT,1)
        XI(2) = XII(INPT,2) 
C     Shape functions and local derivatives
        IF (JTYPE.EQ.ONE) THEN
         CALL SHAPEFUNT(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.TWO) THEN
         CALL SHAPEFUN(AN,dNdxi,XI)
        ENDIF
C     Jacobian
        DO I = 1,2
         DO J = 1,2
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
C        
        DTM = ZERO
        DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        END IF
C     Inverse of Jacobian
        VJABOBINV(1,1)=VJACOB(2,2)/DTM
        VJABOBINV(1,2)=-VJACOB(1,2)/DTM
        VJABOBINV(2,1)=-VJACOB(2,1)/DTM
        VJABOBINV(2,2)=VJACOB(1,1)/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,2
          dNdx(K,I) = ZERO
          DO J = 1,2
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
C
C     Calculating B matrix (B=LN)
       DO INODE=1,NNODE
        BP(1,INODE)=dNdx(INODE,1)
        BP(2,INODE)=dNdx(INODE,2)
       END DO
C
        PHASE=ZERO
        DPHASE=ZERO
        DO I=1,NDOFEL
         PHASE=PHASE+AN(I)*U(I)
        END DO
        DO I=1,NDOFEL
         DPHASE=DPHASE+AN(I)*DU(I,1)
        END DO
C        
        IF (STEPITER.EQ.ZERO) THEN
          SDV(1)=PHASE-DPHASE
        ELSE
          SDV(1)=PHASE
        ENDIF
C
C     Gradient
        DO I=1,2
         DP(I)=ZERO
        END DO
        DO I=1,2
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*U(J)
         END DO
        END DO
C
        IF (STEPITER.EQ.ZERO) THEN
         ENGN=USRVAR(JELEM,13,INPT)
        ELSE
         ENGN=USRVAR(JELEM,16,INPT)
        ENDIF
C        
        HISTN=USRVAR(JELEM,16,INPT)
        IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
        ELSE
         HIST=HISTN
        ENDIF
        SDV(2)=HIST

        DO I=1,NNODE
         DO K=1,NNODE
          DO J=1,2
           AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1      THCK*GCPAR*CLPAR*AINTW(INPT)
          END DO
          AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*THCK*
     1     AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
         END DO
        END DO
C        
        DO I=1,NDOFEL
         DO J=1,2
           RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1      AINTW(INPT)*DTM*THCK
         END DO
         RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*THCK*
     1    ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
        END DO
C
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         USRVAR(JELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
       END DO
       
C     ==================================================================
      ELSEIF ((JTYPE.EQ.THREE).OR.(JTYPE.EQ.FOUR)) THEN
      STEPITER=USRVAR(JELEM-N_ELEM,18,1)

       EMOD = PROPS(1)
       ENU = PROPS(2)
       THCK = PROPS(3)
       PARK = PROPS(4)

       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO

       IF (JTYPE.EQ.THREE) THEN
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = ONE
        AINTW(1) = HALF
       ELSEIF (JTYPE.EQ.FOUR) THEN
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = FOUR
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF

       DO INPT=1,INNODE
        DO I=1,NSTVTT
          SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
        END DO

        XI(1) = XII(INPT,1)
        XI(2) = XII(INPT,2) 

        IF (JTYPE.EQ.THREE) THEN
         CALL SHAPEFUNT(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.FOUR) THEN
         CALL SHAPEFUN(AN,dNdxi,XI)
        ENDIF

        IY=ZERO
        DO I = 1,NNODE
         IX=IY+1
         IY=IX+1
         VNI(1,IX)=AN(I)
         VNI(1,IY)=ZERO
         VNI(2,IX)=ZERO
         VNI(2,IY)=AN(I)
        END DO

        DO I = 1,2
         DO J = 1,2
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO

        DTM = ZERO
        DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT
        ENDIF

        VJABOBINV(1,1)=VJACOB(2,2)/DTM
        VJABOBINV(1,2)=-VJACOB(1,2)/DTM
        VJABOBINV(2,1)=-VJACOB(2,1)/DTM
        VJABOBINV(2,2)=VJACOB(1,1)/DTM

        DO K = 1,NNODE
         DO I = 1,2
          dNdx(K,I) = ZERO
          DO J = 1,2
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO

       IY=0
       DO INODE=1,NNODE
        IX=IY+1
        IY=IX+1
        BB(1,IX)= dNdx(INODE,1)
        BB(1,IY)= ZERO
        BB(2,IX)= ZERO
        BB(2,IY)= dNdx(INODE,2)
        BB(3,IX)= dNdx(INODE,2)
        BB(3,IY)= dNdx(INODE,1)
       END DO

        DO I=1,3
         DO J=1,3
          CMAT(I,J)=ZERO
         END DO
        END DO
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
C
        DO J=1,2
         ULOC(J)=ZERO
        END DO
        DO J=1,2
         DO I=1,NDOFEL
          ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
        END DO  
        DO J=1,2
         SDV(J)=ULOC(J)
        END DO
C   
        IF (STEPITER.EQ.ZERO) THEN
         PHASE=USRVAR(JELEM-N_ELEM,15,INPT)
        ELSE
         PHASE=USRVAR(JELEM-N_ELEM,14,INPT)
        ENDIF
C
        SDV(14)=PHASE

        DO J=1,3
         EPS(J)=ZERO
        END DO
        DO I=1,3
         DO J=1,NDOFEL
          EPS(I)=EPS(I)+BB(I,J)*U(J)    
         END DO
        END DO
        DO J=1,3
         SDV(J+2)=EPS(J)
        END DO

        DO K1=1,3
         STRESS(k1)=ZERO
        END DO
        DO K1=1,3
         DO K2=1,3
          STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
         END DO
        END DO
        DO J=1,3
         SDV(J+5)=STRESS(J)*((ONE-PHASE)**TWO+PARK)
        END DO
        DO J=1,3
         SDV(J+8)=STRESS(J)
        END DO

        ENG=ZERO
        DO K2=1,3
         ENG=ENG+STRESS(K2)*EPS(K2)*HALF
        END DO
        SDV(12)=ENG*((ONE-PHASE)**TWO+PARK)
        SDV(13)=ENG
        ENERGY(2)=ENG

        DO K=1,NDOFEL
         DO L=1,NDOFEL
          DO I=1,3
           DO J=1,3
            AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1       BB(J,L)*DTM*THCK*((ONE-PHASE)**TWO+PARK)
           END DO
          END DO
         END DO
        END DO

        DO K1=1,NDOFEL
         DO K4=1,3
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1      THCK*((ONE-PHASE)**TWO+PARK)
         END DO
        END DO
C       
        DO I=1,NSTVTT
         SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
         USRVAR(JELEM-N_ELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
        END DO
       END DO
      ENDIF
C      
      RETURN
      END

      SUBROUTINE SHAPEFUN(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))

      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END

      SUBROUTINE SHAPEFUNT(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)


      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)

      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
      RETURN
      END      

       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=3,NSTV=18) 
       DATA NEWTON,TOLER/40,1.D-6/ 
      
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)

       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)

       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO

       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO 

       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO

       DO K1=1, NTENS
        DO K2=1, NTENS
         STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
       END DO 

       NELEMAN=NOEL-TWO*N_ELEM
       IF (NPT.EQ.3) THEN
        NPT=4
       ELSEIF (NPT.EQ.4) THEN
        NPT=3
       ENDIF

       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO

       RETURN
       END
