C ---------------PARAMETER FILE CHEMISTRY--------------------------------
C
C LG- This parameter file contains the different parameter and its values
C     of the chemistry scheme for ECHAM and RACMO, and contains the 
C     listing of the trace gases distinguished within the hydrocarbon
C     chemistry code and its assigned numbers 
C -----------------------------------------------------------------------------
      
C LG- NLEVV is the number of levels within the vegetation and
C     NLEVT is the total no. of levels, NLEVATM is the no. of levels
C     in the atmosphere and NLONT30 and NLATT30 are the no. of longitudes/
C     latitudes for the T30 ECHAM resolution, NTRAC is the number of 
C     long-lived species which are transported, if the wrong number of 
C     tracers is declared then the user is urged to reset the number of 
C     species to the right value (10 long lived and 12 short lived species
C     and radon for the ECHAM basic chemistry code, 23 long lived, 16 
C     short lived species and radon and three sulfur species for the 
C     ECHAM CBM4 chemistry scheme and 24 long lived and 12 short lived 
C     species and radon for TM3 CBM4 chemistry scheme, March 1998, 
C     see "prechem.f")

      INTEGER NLEVV,NLEVV_ML,NLEVS_ML,NLEVATM,NLEVT,NTRAC,NG3X,NTRACT,NSPG3X,
     +        NNLON,NLONT30,NLATT30,NLONT21,NLATT21,NTRAC_ECH,NTRAC_CBM4_ECH,
     +        NTRAC_CBM4_TM3,NTRAC_TM3,NLEV_TM3,NTRAC_CHEM_MIM,NTRAC_MECCA1,    ! ESS_lg_20100225+ NTRAC_MECCA1 
     +        NTRAC_MECCA1v16,NTRAC_MECCA1v19,NTR,NSPEC_EMISVOC,NRESUM,NCCTILES,! ESS_lg_20100225+ NTRAC_MECCA1 (v1.6 and v1.9c)
     +        NTRAC_AER,NTRAC_GAS                                               ! ESS_lg_20100726+ NTRAC_AER and NTRAC_GAS

C LG- in the basic version of the this chemistry code, being implemented in
C     ECHAM only the long lived species are being transported, however in
C     The 1-D model for all the species the dynamical processes are being
C     considered and thus resembles NTRAC the total number of species. 
C     NG3X, which represents the number of short lived species is and set to 
C     0. The parameter NTR, which resembles the number of tracers for
C     which the dry deposition is calculated should at least be larger 
C     than NTRAC (NTRACT) + the extra number of tracers being considered
C     within the dry deposition scheme but not in the chemistry scheme
C     (see EC4_VDIFF.f)

C LG- 08-2002, normally 50 tracers with NTRAC=34 in the CBM4_ECH scheme, 
C LG- 2002,    added the monoterpens (alpha + beta pinene), NTRACT=52 with NTRAC=36 
C LG- 08-2002, added CH3CL AND CHCL3, NTRACT=54 with NTRAC=38
C LG- 03-2003, added the sesquiterpenes as a group, NTRACT=55 with NTRAC=39
C LG_ 12-2003, added HONO, NTRACT=56 with NTRAC=40
C LG_ 08-2004, added CH2OHO2H (HMHP) and RCHOHO2H (HAHP), NTRACT=58 with NTRAC=42
C LG_ 12-2005, added CH3OH and CH3CN (acetonitrile) NTRACT=60 with NTRAC=44
C LG_ 03-2006, added 2-double and 1-double bond sesquiterpenes and extra
C              monoterpene (terpinolene) NTRACT=63 with NTRAC=47
C LG_ 06-2010, added hexane, 1,3butadiene and tetramethylbenzene to test 
C              anthropogenic alkenes, NTRACT=66 with NTRAC=50

! ESS_lg_2010, added the option to use the MECCA1 chemistry scheme

      PARAMETER (NLEVV=4,NLEVV_ML=2,NLEVS_ML=NLEVV_ML,NLEVATM=19,
     +           NLEVT=NLEVATM+NLEVV,
     +           NNLON=1,NLONT30=96,NLATT30=48,NLONT21=64,NLATT21=32,
     +           NLEV_TM3=19,NTR=200,NSPEC_EMISVOC=9,  ! ESS_lg_20070802+ 9 VOC emission species ! ESS_lg_20100326+ updated NTR, 20100707+ should be larger then NTRAC
     +           NRESUM=0,NCCTILES=1, ! mz_lg_20060628+
     +           NTRAC_ECH=24,NTRAC_CBM4_ECH=66,NTRAC_CBM4_TM3=38,NTRAC_TM3=30,  
     +           NTRAC_CHEM_MIM=42,
     +           NTRAC_MECCA1v16=67,  ! ESS_lg_20100707+ added NTRAC_MECCA1 for MECCA1, v1.6 
     +           NTRAC_MECCA1v19=189, ! ESS_lg_20100707+ added NTRAC_MECCA1 for MECCA1, version MIM3mam including sulphur, #27 of MECCA1 preprocessing 
     +           NTRAC_MECCA1=189,    ! ESS_lg_20100225+ added NTRAC_MECCA1, the number should resemble the parameter NSPEC in messy_mecca_kpp_g_mem 
     +           NTRAC=66,NG3X=0,NTRACT=NTRAC+NG3X)    ! ESS_lg_20100301+ moved to here the final assignment of the actual no. of tracers, 
                                                        ! NTRAC should resemble the number of gas and aerosol tracers with NTRAC_AER     
                                                        ! being assigned in oned.f

C LG- declaration of term for recalculation of
C     concentrations of chemistry scheme in molecules cm-3 to
C     ppb's (see PRECHEM.f), definition of specific mass of air and avogadro

      REAL ZMAIR,AVO,ZGASC

      PARAMETER (ZMAIR=28.9644,AVO=6.022E23,ZGASC=8.3144)

C LG- definition of maximum canopy height and minimum thickness of the lowest
C     canopy layer which are being used in different routines. Added 07-2002,
C     the maximum LAI to cut-off the LAI inferred from the NDVI data

      REAL MAX_CANHEIGHT, MINTHICK, LAIMAX

      PARAMETER (MAX_CANHEIGHT=40.,MINTHICK=0.25, LAIMAX=7.)

C LG- definition of parameters used to define the dimension of the
C     budget arrays, NHOR is number of horizontal compartiments,
C     NVER is no. of vertical compartiments and NFIELD is no. of fields
C     NDIM=NHOR*NVER, NREAC is the no. of reactions which are considered
C     for the budgets 

      INTEGER NHORIZ,NVERT,NFIELD,NDIM,NREAC

      PARAMETER (NHORIZ=1,NVERT=4,NFIELD=9,NDIM=NHORIZ*NVERT,
     +           NREAC=60) ! ESS_lg_20100609+ added to deal with check bounds compilation

C LG- parameters needed for budget calculations, IBWDL: large scale wet
C     depositon, IBWDC: convective wet deposition, IBNOYCORR: budgets after
C     correction of NOy (see TRANSP.f) and IBNEG: budgets after correction
C     for negative values (see TRANSP.f)

      INTEGER IBEMIS,IBDDEP,IBCHEM,IBWDL,IBWDC,IBNOYCORR,
     *        IBNEG,IBTRANS

      PARAMETER(IBEMIS=1,IBDDEP=3,IBCHEM=2,IBWDL=4,IBWDC=5,
     *          IBNOYCORR=7,IBNEG=8,IBTRANS=9)

