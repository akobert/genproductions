C ************************************************************
C Source for the library implementing a bias function that 
C populates the large pt tale of the leading jet. 
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) ptj_bias_target_ptj : target ptj value
C    > (double precision) ptj_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [ptj(evt)/mean_ptj]^enhancement_power
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the 
C relevant parameters for this bias module.
C
C  parameters = {}
C

C ************************************************************
C Helpers
C ************************************************************

      DOUBLE PRECISION FUNCTION mass_Z(im)
        implicit none
        double precision im
        double precision ialphaL
        double precision ialphaR
        double precision ixbar
        double precision isigmaL
        double precision isigmaR
        double precision ia2
        double precision iaI
        double precision inL
        double precision inR
        double precision massnorm
        double precision ix
        double precision var
        double precision iA
        double precision iB

C      Fit mass spectrum with normal errors. Misses the low end by a
C      bit.
C        ialpha =  0.428533
C        ixbar  = -209.945
C        isigma =  51.7438
C        iquadr =  0.000160548

C      Fit mass spectrum with reduced errors at the low end. Better low
C      end, worse shoulder.
C        ialpha =  0.368106
C        ixbar = -231.695
C        isigma =  59.7796
C        iquadr = -0.000328664

C       New Fit

        ialphaL = 1.4032
        ialphaR = 0.74733
        ixbar = 0.17873
        isigmaL =  .070031
        isigmaR =  .083945
        ia2 = -.083845
        iaI = -.0093637
        inL = 11.971
        inR = 4.3462

C       Normalized Mass
        massnorm = im / 400.

C       Transformed Mass
        ix = massnorm + ia2*massnorm*massnorm + iaI/(massnorm+0.001)

        IF (ix - ixbar .LE. 0.0) THEN
          var = (ix - ixbar)/isigmaL
          iA = (inL/abs(ialphaL))**inL * exp((-abs(ialphaL)**2) / 2.)
          iB = inL/abs(ialphaL) - abs(ialphaL)
        ELSE
          var = (ix - ixbar)/isigmaR
          iA = (inR/abs(ialphaR))**inR * exp((-abs(ialphaR)**2) / 2.)
          iB = inR/abs(ialphaR) - abs(ialphaR)
        ENDIF

        IF (var .LT. -ialphaL) THEN
          mass_Z = iA * (iB - var)**(-inL)
        ELSE IF (var .LE. ialphaR) THEN
          mass_Z = exp(-1./2. * (var**2))
        ELSE
          mass_Z = iA * (iB + var)**(-inR)
        ENDIF

      END

      DOUBLE PRECISION FUNCTION pt(ipt)
        implicit none
        double precision ipt
        double precision ialphaL
        double precision ialphaR
        double precision ixbar
        double precision isigmaL
        double precision isigmaR
        double precision ia2
        double precision iaI
        double precision inL
        double precision inR
        double precision ptnorm
        double precision ix
        double precision var
        double precision iA
        double precision iB


C       New Fit

        ialphaL = 29.346
        ialphaR = .49846
        ixbar = 0.087651
        isigmaL =  .014843
        isigmaR =  .028275
        ia2 = 4.99
        iaI = .000042703
        inL = .12488
        inR = 3.2167

C       Normalized pt
        ptnorm = (ipt - 20.) / (1000. - 20.)

C       Transformed pt
        ix = ptnorm + ia2*ptnorm*ptnorm + iaI/(ptnorm+0.001)

        IF (ix - ixbar .LE. 0.0) THEN
          var = (ix - ixbar)/isigmaL
          iA = (inL/abs(ialphaL))**inL * exp((-abs(ialphaL)**2) / 2.)
          iB = inL/abs(ialphaL) - abs(ialphaL)
        ELSE
          var = (ix - ixbar)/isigmaR
          iA = (inR/abs(ialphaR))**inR * exp((-abs(ialphaR)**2) / 2.)
          iB = inR/abs(ialphaR) - abs(ialphaR)
        ENDIF

        IF (var .LT. -ialphaL) THEN
          pt = iA * (iB - var)**(-inL)
        ELSE IF (var .LE. ialphaR) THEN
          pt = exp(-1./2. * (var**2))
        ELSE
          pt = iA * (iB + var)**(-inR)
        ENDIF

      END

      subroutine bias_wgt(p, original_weight, bias_weight)
          implicit none
C
C Parameters
C
          include '../../maxparticles.inc'         
          include '../../nexternal.inc'
C
C Accessing the details of the event
C
          include '../../run_config.inc'
          include '../../lhe_event_infos.inc'
C
C Event kinematics
C
          double precision p(0:3,nexternal)
          double precision original_weight, bias_weight
          double precision rho
          double precision mZp
          double precision pTZp
          double precision pTy

C Cut variables, from run card          
          double precision p4(4)
          integer i
          integer j
          integer k
          double precision mass_Z
          double precision mass_y
          double precision mass_weight
          double precision pt
          double precision pt_weight
C
C Global variables
C
C
C Mandatory common block to be defined in bias modules
C
          double precision stored_bias_weight
          data stored_bias_weight/1.0d0/          
          logical impact_xsec, requires_full_event_info
C         Don't unweight the bias: we really want a flat distribution
          data impact_xsec/.True./
C         Of course this module does not require the full event
C         information (color, resonances, helicities, etc..)
          data requires_full_event_info/.True./ 
          common/bias/stored_bias_weight,impact_xsec,
     &                requires_full_event_info
C
C Accessingt the details of the event
C
C --------------------
C BEGIN IMPLEMENTATION
C --------------------
          include '../../run.inc'
          include '../../cuts.inc'

          include '../bias.inc'

          p4(1) = 0.
          p4(2) = 0.
          p4(3) = 0.
          p4(4) = 0.
          mZp = -1.
          pTZp = -1.

          DO i=1,npart
            IF (jpart(1,i) .eq. 55) THEN
              p4(1) = pb(1,i)
              p4(2) = pb(2,i)
              p4(3) = pb(3,i)
              p4(4) = pb(0,i)
              mZp = sqrt(p4(4)**2 - p4(1)**2 - p4(2)**2 - p4(3)**2)
              pTZp = sqrt(p4(1)**2 + p4(2)**2)
            ELSE IF (jpart(1,i) .EQ. 22) THEN
              p4(1) = pb(1,i)
              p4(2) = pb(2,i)
              p4(3) = pb(3,i)
              p4(4) = pb(0,i)
              pTy = MAX(sqrt(p4(1)**2 + p4(2)**2), pTy)
            ENDIF
          ENDDO

C Rough cuts on mass, rho, and pT
          IF ((mZp .GT. 0.) .AND. (pTZp .GT. 0.)) THEN
            rho = LOG(mZp**2 / pTZp**2)
            IF (
     &         ((-8.0 .LT. rho) .AND. (rho .LT. -0.5)) .AND.
     &         ((20 .LT. pTZp) .AND. (pTZp .LT. 1000.)) .AND.
     &         ((5 .LT. mZp) .AND. (mZp .LT. 400.))) THEN

              mass_weight = mass_Z(mZp)
              pt_weight = pt(pTy)
              IF ((mass_weight .GT. 0) .AND. (pt_weight .GT. 0))
THEN
                bias_weight = 1. / (mass_weight * pt_weight)
              ELSE
                bias_weight = 0.00000001
              ENDIF
            ELSE
              bias_weight = 0.00000001
            ENDIF
          ELSE
C            print*,"jkl;"
            bias_weight = 0.00000001
          ENDIF
          RETURN
      END subroutine bias_wgt

