      SUBROUTINE INIAGS

      INCLUDE "paramveg.h"
      INCLUDE "comjac.h"
      
C LG- 05-2001, added the vegetation type number 3 to represent a special
C     class for coniferous trees in terms of some applied constants for the
C     physiological model

      CGAMMA25(1)=45.
      CGAMMA25(2)=2.8
      CGAMMA25(3)=45.
      CGAMMAQ10(1)=1.5
      CGAMMAQ10(2)=1.5
      CGAMMAQ10(3)=1.5
      C25=25.
      C10=10.
      CGM25(1)=7.0      ! normally 7.0, reset to smaller value to check the Rstom
                        ! Reinder Ronda has proposed value of 3.8 for tropical 
			! forest !
      CGM25(2)=17.5
      CGM25(3)=7.0
      CGMQ10(1)=2.0
      CGMQ10(2)=2.0
      CGMQ10(3)=2.0
      CDOT3=.3
      CGMT1(1)=5.
      CGMT1(2)=13.
      CGMT1(3)=5.
      CGMT2(1)=28.
      CGMT2(2)=36.
      CGMT2(3)=28.
      CGC(1)=.25	! gmin,c, see report by van de Kassteele
      CGC(2)=.25
      CGC(3)=.25	! 1., ESS_lg_20110325+ the value for coniferous forest was one but this resulted in way too
                        !   large resistance and small CO2 flux for Hyytiala ! a value of 0 gives too small uptake rates for
      CF0(1)=.89	! tropical forest canopies !!!
      CF0(2)=.85
      CF0(3)=.4       ! 0.4 selected for coniferous trees (Kassteele)
      CRICO(1)=-0.07
      CRICO(2)=-0.15
      CRICO(3)=-0.15	 !different CRICO for vegetation type 3
      CCS=340.
      CEPSILON0(1)=.017
      CEPSILON0(2)=.014
      CEPSILON0(3)=.017
      CAMAX25(1)=2.2	! normally 2.2, reset to smaller value to check the Rstom
                        ! Reinder Ronda has proposed value of 1 for tropical 
			! forest !

      CAMAX25(2)=1.7
      CAMAX25(3)=0.45	! normally 0.45, different CAMAX25 for vegetation type 3
      CAMAXQ10(1)=2.
      CAMAXQ10(2)=2.
      CAMAXQ10(3)=2.
      CAMAXT1(1)=8.
      CAMAXT1(2)=13.
      CAMAXT1(3)=8.
      CAMAXT2(1)=38.
      CAMAXT2(2)=38.
      CAMAXT2(3)=38.
      C9=9.
      C1DOT6=1.6
      CMA=28.9
      CMC02=44.
      C1000=1000.
      CMV=18.
      CEXPAR=.7
      RETURN
      END
