	  Ðm  è   k820309    h          14.0        thU                                                                                                           
       messy_airsea.f90 MESSY_AIRSEA       C       DP SP AIRSEA_READ_NML_CTRL AIRSEA_ORO AIRSEA_LAYERTHICKNESS AIRSEA_CALC_DENSITY AIRSEA_CALC_CAIR AIRSEA_CALC_SCHMIDT_AIR AIRSEA_CALC_SCHMIDT_SEA AIRSEA_CALC_OSTWALD AIRSEA_CALC_KL_SEA AIRSEA_CALC_KL_SEA_WC AIRSEA_CALC_KL_AIR AIRSEA_CALC_KL_TOT AIRSEA_CALC_WC AIRSEA_CALC_HENRY AIRSEA_WATER_CONC AIRSEA_DELTA_CONC AIRSEA_FLUX CALC_R1 CALC_R2 O3_DEPOSITION BESSK0_S BESSK0_V BESSK1_V BESSK1_S BESSI0_S BESSI0_V BESSI1_S BESSI1_V MODSTR MODVER L_CLIM L_TENDENCY L_WHITECAP L_CH4 L_HCHO L_CH3OH L_C2H6 L_C2H4 L_CH3CHO L_C3H8 L_C3H6 L_CH3COCH3 NASI L_REQUEST_ASI ASINAME IDASI_CH4 IDASI_HCHO IDASI_CH3OH IDASI_C2H6 IDASI_C2H4 IDASI_CH3CHO IDASI_C3H8 IDASI_C3H6 IDASI_CH3COCH3 HENRY_A HENRY_B MOLVOL ALPHA L_ASI L_NASI TRAC_ID NUM_TRAC_USED ZCIM ZCC2H4 ZCC3H6                                                     
       DP SP                                                                                                                                                                                                                                                                                                                 Cairsea                                                                                                                    C1.0                          D@                                                       D@                                                       D@                                                       D@                                 	                      D@                                 
                      D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                                                                                                 	               9           @                                     	                    p          p 	           p 	                         +                                               	       	             p          p 	           p 	                                                                                                                                    1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                    	               9                                               	              
      p          p 	           p 	                                                                        	              
      p          p 	           p 	                                                                         	              
      p          p 	           p 	                                                                   !     	              
      p          p 	           p 	                                                                     "     	                    p          p 	           p 	                                                                     #                                                        $     	                    p          p 	           p 	                                                                     %                                                        &     
                                                   '     
                                                   (     
       #         @                                   )                    #STATUS *   #IOU +             D                                 *                      
  @                               +           #         @                                   ,                    #LOSEA -   #SEALANDMASK .   #SEAICE /             D                                 -                    7              &                                                     
                                 .                   
 8             &                                                     
                                 /                   
 9             &                                           %         H                               0                    
       #GEOPOT 1                                          
                                 1     
      #         @                                   2                    #TEMP 3   #SPHUM 4   #PRESS 5   #DENSITY 6                                                                                                          
                                 3                   
 ,             &                                                     
                                 4                   
 -             &                                                     
                                 5                   
 .             &                                                     D                                6                   
 /              &                                           #         @                                   7                    #TEMP 8   #SPHUM 9   #PRESS :   #CAIR ;                                                                                                 
                                 8                   
 3             &                                                     
                                 9                   
 4             &                                                     
                                 :                   
 5             &                                                     D                                ;                   
 6              &                                           #         @                                   <                   #AIRSEA_CALC_SCHMIDT_AIR%SQRT =   #AIRSEA_CALC_SCHMIDT_AIR%SIZE >   #TEMP ?   #SC_AIR @   #NUM_TRAC A   #MOLW B   #PRESS C                @                             =     SQRT              @                             >     SIZE        0  
 @                              ?                   
 A             &                                                     D                                @                   
 C              &                                                     
                                  A                     
                                 B     
                
                                 C                   
 B             &                                           #         @                                   D                   #AIRSEA_CALC_SCHMIDT_SEA%SIZE E   #SST F   #SC_SEA G   #NUM_TRAC H                @                             E     SIZE        0  
 @                              F                   
 ?             &                                                     D                                G                   
 >              &                                                     
                                  H           #         @                                  I                    #SST J   #HENRY K   #OST L             
                                 J                   
 F             &                                                     
                                 K                   
 G             &                                                     D                                L                   
 H              &                                           #         @                                   M                   #AIRSEA_CALC_KL_SEA%SQRT N   #WIND O   #SC_SEA P   #K_VEL Q   #PROMA R   #LOSEA S                @                             N     SQRT           
                                 O                   
              &                                                     
                                 P                   
               &                                                     D                                Q                   
 !              &                                                     
                                  R                     
                                  S                    "             &                                           #         @                                   T                   #AIRSEA_CALC_KL_SEA_WC%SQRT U   #WIND V   #SC_SEA W   #K_VEL X   #PROMA Y   #LOSEA Z   #WC [   #HENRY \   #SST ]                @                             U     SQRT           
                                 V                   
 #             &                                                     
                                 W                   
 $             &                                                     D                                X                   
 %              &                                                     
                                  Y                     
                                  Z                    &             &                                                     
                                 [                   
 '             &                                                     
  @                              \                   
 (             &                                                     
  @                              ]                   
 )             &                                           #         @                                   ^                   #AIRSEA_CALC_KL_AIR%SQRT _   #AIRSEA_CALC_KL_AIR%SIZE `   #KG_VEL a   #ZUST b   #U c   #V d   #SC_AIR e   #PROMA f   #LOSEA g                @                             _     SQRT              @                             `     SIZE           D                                a                   
               &                                                  0  
 @                              b                   
              &                                                  0  
 @                              c                   
              &                                                     
                                 d                   
              &                                                     
  @                              e                   
              &                                                     
                                  f                     
                                  g                                 &                                           #         @                                   h                   #AIRSEA_CALC_KL_TOT%SIZE i   #KL_SEA j   #KL_AIR k   #KL_TOT l   #ALPHA m   #HENRY n   #TEMP o   #PROMA p   #LOSEA q                @                             i     SIZE           
                                 j                   
 	             &                                                     
                                 k                   
 
             &                                                     D                                l                   
               &                                                     
                                 m     
             0  
 @                              n                   
              &                                                     
                                 o                   
              &                                                     
                                  p                     
                                  q                                 &                                           #         @                                   r                    #WIND s   #WC t   #LOSEA u   #PROMA v             
                                 s                   
 0             &                                                     D                                t                   
 1              &                                                     
                                  u                    2             &                                                     
                                  v           #         @                                   w                   #AIRSEA_CALC_HENRY%EXP x   #PA y   #PB z   #TEMP {   #HENRY |                @                             x     EXP           
                                 y     
                
                                 z     
                
                                 {                   
 O             &                                                     D                                |                   
 P              &                                           #         @                                   }                    #CW ~   #NUM_TRAC              D                                ~                   
 =              &                                                     
                                             #         @                                                       #CX    #KH    #CW    #PRESS    #DELTA_C    #LOSEA    #PROMA    #NUM_TRAC              
                                                    
 I             &                                                     
                                                    
 K             &                                                     
                                                    
 J             &                                                     
                                                    
 L             &                                                     D                                                   
 M              &                                                     
                                                      N             &                                                     
                                                       
                                             #         @                                                       #DELTA_C    #KL    #FLUX              
                                                    
 :             &                                                     
                                                    
 ;             &                                                     D                                                   
 <              &                                           #         @                                                      #R1    #WIND    #ZUST              D                                                   
               &                                                     
                                                    
              &                                                     
                                                    
              &                                           #         @                                                      #R2    #ZUST    #SC_AIR              D                                                   
               &                                                     
                                                    
              &                                                     
                                                    
              &                                           #         @                                                      #O3_DEPOSITION%MAX    #O3_DEPOSITION%LOG    #O3_DEPOSITION%EXP    #O3_DEPOSITION%SQRT    #O3_DEPOSITION%SIZE    #LSEA    #SST    #WIND10    #UM1SL    #VM1SL    #GEOPOTSL     #PRESS ¡   #TEMP ¢   #ZUST £   #SRFL ¤   #SPHUM ¥   #AHFS ¦   #AHFL §   #GYRE ¨   #FRAC_COAST ©   #CHLOROFYLL ª   #PBLH «   #IM ¬   #NO3M ­   #DMS ®   #C2H4 ¯   #C3H6 °   #KORG ±   #PROMA ²   #ZOUTPUT ³   #VDO3_OCEAN ´                 @                                 MAX               @                                 LOG               @                                 EXP               @                                 SQRT              @                                   SIZE           
                                                                   &                                                     
                                                    
 w             &                                                     
                                                    
 x             &                                                     
                                                    
 y             &                                                     
                                                    
 z             &                                                     
                                                     
 {             &                                                     
                                 ¡                   
 |             &                                                     
                                 ¢                   
 }             &                                                     
                                 £                   
 ~             &                                                     
                                 ¤                   
              &                                                     
                                 ¥                   
              &                                                     
                                 ¦                   
              &                                                     
                                 §                   
              &                                                     
                                 ¨                   
              &                                                     
                                 ©                   
              &                                                     
                                 ª                   
              &                                                     
                                 «                   
              &                                                     
                                 ¬                   
              &                                                     
                                 ­                   
              &                                                     
                                 ®                   
              &                                                     
                                 ¯                   
              &                                                     
                                 °                   
              &                                                     
                                 ±     
                
                                  ²                     D                                ³                   
               &                   &                                                     D                                ´                   
               &                                           %         @                               µ                   	       #BESSK0_S%LOG ¶   #BESSK0_S%EXP ·   #BESSK0_S%SQRT ¸   #X ¹                                                                     @                            ¶     LOG               @                            ·     EXP               @                            ¸     SQRT           
  @                              ¹     	      (        `                                º                    c               	    #BESSK0_V%LOG »   #BESSK0_V%EXP ¼   #BESSK0_V%SQRT ½   #BESSK0_V%SIZE ¾   #BESSK0_V%ALL ¿   #X À   p          H r ¾     7	S	O p        j            j                                      H r ¾     7	S	O p        j            j                                                                                                                      @                            »     LOG               @                            ¼     EXP               @                            ½     SQRT               @                            ¾     SIZE               @                            ¿     ALL        0  
 @                              À                   	 ]             &                                           (        `                                Á                    l               	    #BESSK1_V%LOG Â   #BESSK1_V%EXP Ã   #BESSK1_V%SQRT Ä   #BESSK1_V%SIZE Å   #BESSK1_V%ALL Æ   #X Ç   p          H r Å     7	S	O p        j            j                                      H r Å     7	S	O p        j            j                                                                                                                      @                            Â     LOG               @                            Ã     EXP               @                            Ä     SQRT               @                            Å     SIZE               @                            Æ     ALL        0  
 @                              Ç                   	 f             &                                           %         @                               È                   	       #BESSK1_S%LOG É   #BESSK1_S%EXP Ê   #BESSK1_S%SQRT Ë   #X Ì                                                                     @                            É     LOG               @                            Ê     EXP               @                            Ë     SQRT           
  @                              Ì     	      %         @                               Í                   	       #BESSI0_S%ABS Î   #BESSI0_S%EXP Ï   #BESSI0_S%SQRT Ð   #BESSI0_S%REAL Ñ   #X Ò                                                                     @                            Î     ABS               @                            Ï     EXP               @                            Ð     SQRT               @             @              Ñ     REAL           
  @                              Ò     	      (        `                               Ó                    Z               	    #BESSI0_V%ABS Ô   #BESSI0_V%EXP Õ   #BESSI0_V%SQRT Ö   #BESSI0_V%SIZE ×   #BESSI0_V%REAL Ø   #X Ù   p          H r ×     7	S	O p        j            j                                      H r ×     7	S	O p        j            j                                                                                                                      @                            Ô     ABS               @                            Õ     EXP               @                            Ö     SQRT               @                            ×     SIZE               @             @              Ø     REAL        0  
 @                              Ù                   	 S             &                                           %         @                               Ú                   	       #BESSI1_S%ABS Û   #BESSI1_S%EXP Ü   #BESSI1_S%SQRT Ý   #BESSI1_S%REAL Þ   #X ß                                                                     @                            Û     ABS               @                            Ü     EXP               @                            Ý     SQRT               @             @              Þ     REAL           
  @                              ß     	      (        `                               à                    v               	    #BESSI1_V%ABS á   #BESSI1_V%EXP â   #BESSI1_V%SQRT ã   #BESSI1_V%SIZE ä   #BESSI1_V%REAL å   #X æ   p          H r ä     7	S	O p        j            j                                      H r ä     7	S	O p        j            j                                                                                                                      @                            á     ABS               @                            â     EXP               @                            ã     SQRT               @                            ä     SIZE               @             @              å     REAL        0  
 @                              æ                   	 o             &                                                  &      fn#fn "   Æ     b   uapp(MESSY_AIRSEA )   È  F   J  MESSY_MAIN_CONSTANTS_MEM ,     p       DP+MESSY_MAIN_CONSTANTS_MEM ,   ~  p       SP+MESSY_MAIN_CONSTANTS_MEM    î         MODSTR    u         MODVER    ù  @       L_CLIM    9  @       L_TENDENCY    y  @       L_WHITECAP    ¹  @       L_CH4    ù  @       L_HCHO    9  @       L_CH3OH    y  @       L_C2H6    ¹  @       L_C2H4    ù  @       L_CH3CHO    9  @       L_C3H8    y  @       L_C3H6    ¹  @       L_CH3COCH3    ù  q       NASI    j	         L_REQUEST_ASI    þ	         ASINAME    
  q       IDASI_CH4      q       IDASI_HCHO    |  q       IDASI_CH3OH    í  q       IDASI_C2H6    ^  q       IDASI_C2H4    Ï  q       IDASI_CH3CHO    @  q       IDASI_C3H8    ±  q       IDASI_C3H6    "  q       IDASI_CH3COCH3             HENRY_A    '         HENRY_B    »         MOLVOL    O         ALPHA    ã         L_ASI    w  @       L_NASI    ·         TRAC_ID    K  @       NUM_TRAC_USED      @       ZCIM    Ë  @       ZCC2H4      @       ZCC3H6 %   K  ]       AIRSEA_READ_NML_CTRL ,   ¨  @   a   AIRSEA_READ_NML_CTRL%STATUS )   è  @   a   AIRSEA_READ_NML_CTRL%IOU    (  p       AIRSEA_ORO !        a   AIRSEA_ORO%LOSEA '   $     a   AIRSEA_ORO%SEALANDMASK "   °     a   AIRSEA_ORO%SEAICE &   <  y       AIRSEA_LAYERTHICKNESS -   µ  @   a   AIRSEA_LAYERTHICKNESS%GEOPOT $   õ  Ò       AIRSEA_CALC_DENSITY )   Ç     a   AIRSEA_CALC_DENSITY%TEMP *   S     a   AIRSEA_CALC_DENSITY%SPHUM *   ß     a   AIRSEA_CALC_DENSITY%PRESS ,   k     a   AIRSEA_CALC_DENSITY%DENSITY !   ÷  Æ       AIRSEA_CALC_CAIR &   ½     a   AIRSEA_CALC_CAIR%TEMP '   I     a   AIRSEA_CALC_CAIR%SPHUM '   Õ     a   AIRSEA_CALC_CAIR%PRESS &   a     a   AIRSEA_CALC_CAIR%CAIR (   í  Å       AIRSEA_CALC_SCHMIDT_AIR -   ²  =      AIRSEA_CALC_SCHMIDT_AIR%SQRT -   ï  =      AIRSEA_CALC_SCHMIDT_AIR%SIZE -   ,     a   AIRSEA_CALC_SCHMIDT_AIR%TEMP /   ¸     a   AIRSEA_CALC_SCHMIDT_AIR%SC_AIR 1   D  @   a   AIRSEA_CALC_SCHMIDT_AIR%NUM_TRAC -     @   a   AIRSEA_CALC_SCHMIDT_AIR%MOLW .   Ä     a   AIRSEA_CALC_SCHMIDT_AIR%PRESS (   P          AIRSEA_CALC_SCHMIDT_SEA -   Ý   =      AIRSEA_CALC_SCHMIDT_SEA%SIZE ,   !     a   AIRSEA_CALC_SCHMIDT_SEA%SST /   ¦!     a   AIRSEA_CALC_SCHMIDT_SEA%SC_SEA 1   2"  @   a   AIRSEA_CALC_SCHMIDT_SEA%NUM_TRAC $   r"  e       AIRSEA_CALC_OSTWALD (   ×"     a   AIRSEA_CALC_OSTWALD%SST *   c#     a   AIRSEA_CALC_OSTWALD%HENRY (   ï#     a   AIRSEA_CALC_OSTWALD%OST #   {$         AIRSEA_CALC_KL_SEA (   %  =      AIRSEA_CALC_KL_SEA%SQRT (   T%     a   AIRSEA_CALC_KL_SEA%WIND *   à%     a   AIRSEA_CALC_KL_SEA%SC_SEA )   l&     a   AIRSEA_CALC_KL_SEA%K_VEL )   ø&  @   a   AIRSEA_CALC_KL_SEA%PROMA )   8'     a   AIRSEA_CALC_KL_SEA%LOSEA &   Ä'  »       AIRSEA_CALC_KL_SEA_WC +   (  =      AIRSEA_CALC_KL_SEA_WC%SQRT +   ¼(     a   AIRSEA_CALC_KL_SEA_WC%WIND -   H)     a   AIRSEA_CALC_KL_SEA_WC%SC_SEA ,   Ô)     a   AIRSEA_CALC_KL_SEA_WC%K_VEL ,   `*  @   a   AIRSEA_CALC_KL_SEA_WC%PROMA ,    *     a   AIRSEA_CALC_KL_SEA_WC%LOSEA )   ,+     a   AIRSEA_CALC_KL_SEA_WC%WC ,   ¸+     a   AIRSEA_CALC_KL_SEA_WC%HENRY *   D,     a   AIRSEA_CALC_KL_SEA_WC%SST #   Ð,  È       AIRSEA_CALC_KL_AIR (   -  =      AIRSEA_CALC_KL_AIR%SQRT (   Õ-  =      AIRSEA_CALC_KL_AIR%SIZE *   .     a   AIRSEA_CALC_KL_AIR%KG_VEL (   .     a   AIRSEA_CALC_KL_AIR%ZUST %   */     a   AIRSEA_CALC_KL_AIR%U %   ¶/     a   AIRSEA_CALC_KL_AIR%V *   B0     a   AIRSEA_CALC_KL_AIR%SC_AIR )   Î0  @   a   AIRSEA_CALC_KL_AIR%PROMA )   1     a   AIRSEA_CALC_KL_AIR%LOSEA #   1  ¿       AIRSEA_CALC_KL_TOT (   Y2  =      AIRSEA_CALC_KL_TOT%SIZE *   2     a   AIRSEA_CALC_KL_TOT%KL_SEA *   "3     a   AIRSEA_CALC_KL_TOT%KL_AIR *   ®3     a   AIRSEA_CALC_KL_TOT%KL_TOT )   :4  @   a   AIRSEA_CALC_KL_TOT%ALPHA )   z4     a   AIRSEA_CALC_KL_TOT%HENRY (   5     a   AIRSEA_CALC_KL_TOT%TEMP )   5  @   a   AIRSEA_CALC_KL_TOT%PROMA )   Ò5     a   AIRSEA_CALC_KL_TOT%LOSEA    ^6  p       AIRSEA_CALC_WC $   Î6     a   AIRSEA_CALC_WC%WIND "   Z7     a   AIRSEA_CALC_WC%WC %   æ7     a   AIRSEA_CALC_WC%LOSEA %   r8  @   a   AIRSEA_CALC_WC%PROMA "   ²8         AIRSEA_CALC_HENRY &   :9  <      AIRSEA_CALC_HENRY%EXP %   v9  @   a   AIRSEA_CALC_HENRY%PA %   ¶9  @   a   AIRSEA_CALC_HENRY%PB '   ö9     a   AIRSEA_CALC_HENRY%TEMP (   :     a   AIRSEA_CALC_HENRY%HENRY "   ;  ^       AIRSEA_WATER_CONC %   l;     a   AIRSEA_WATER_CONC%CW +   ø;  @   a   AIRSEA_WATER_CONC%NUM_TRAC "   8<         AIRSEA_DELTA_CONC %   Ô<     a   AIRSEA_DELTA_CONC%CX %   `=     a   AIRSEA_DELTA_CONC%KH %   ì=     a   AIRSEA_DELTA_CONC%CW (   x>     a   AIRSEA_DELTA_CONC%PRESS *   ?     a   AIRSEA_DELTA_CONC%DELTA_C (   ?     a   AIRSEA_DELTA_CONC%LOSEA (   @  @   a   AIRSEA_DELTA_CONC%PROMA +   \@  @   a   AIRSEA_DELTA_CONC%NUM_TRAC    @  g       AIRSEA_FLUX $   A     a   AIRSEA_FLUX%DELTA_C    A     a   AIRSEA_FLUX%KL !   B     a   AIRSEA_FLUX%FLUX    §B  d       CALC_R1    C     a   CALC_R1%R1    C     a   CALC_R1%WIND    #D     a   CALC_R1%ZUST    ¯D  f       CALC_R2    E     a   CALC_R2%R2    ¡E     a   CALC_R2%ZUST    -F     a   CALC_R2%SC_AIR    ¹F  Ý      O3_DEPOSITION "   H  <      O3_DEPOSITION%MAX "   ÒH  <      O3_DEPOSITION%LOG "   I  <      O3_DEPOSITION%EXP #   JI  =      O3_DEPOSITION%SQRT #   I  =      O3_DEPOSITION%SIZE #   ÄI     a   O3_DEPOSITION%LSEA "   PJ     a   O3_DEPOSITION%SST %   ÜJ     a   O3_DEPOSITION%WIND10 $   hK     a   O3_DEPOSITION%UM1SL $   ôK     a   O3_DEPOSITION%VM1SL '   L     a   O3_DEPOSITION%GEOPOTSL $   M     a   O3_DEPOSITION%PRESS #   M     a   O3_DEPOSITION%TEMP #   $N     a   O3_DEPOSITION%ZUST #   °N     a   O3_DEPOSITION%SRFL $   <O     a   O3_DEPOSITION%SPHUM #   ÈO     a   O3_DEPOSITION%AHFS #   TP     a   O3_DEPOSITION%AHFL #   àP     a   O3_DEPOSITION%GYRE )   lQ     a   O3_DEPOSITION%FRAC_COAST )   øQ     a   O3_DEPOSITION%CHLOROFYLL #   R     a   O3_DEPOSITION%PBLH !   S     a   O3_DEPOSITION%IM #   S     a   O3_DEPOSITION%NO3M "   (T     a   O3_DEPOSITION%DMS #   ´T     a   O3_DEPOSITION%C2H4 #   @U     a   O3_DEPOSITION%C3H6 #   ÌU  @   a   O3_DEPOSITION%KORG $   V  @   a   O3_DEPOSITION%PROMA &   LV  ¤   a   O3_DEPOSITION%ZOUTPUT )   ðV     a   O3_DEPOSITION%VDO3_OCEAN    |W  Â       BESSK0_S    >X  <      BESSK0_S%LOG    zX  <      BESSK0_S%EXP    ¶X  =      BESSK0_S%SQRT    óX  @   a   BESSK0_S%X    3Y  û      BESSK0_V    .[  <      BESSK0_V%LOG    j[  <      BESSK0_V%EXP    ¦[  =      BESSK0_V%SQRT    ã[  =      BESSK0_V%SIZE     \  <      BESSK0_V%ALL    \\     a   BESSK0_V%X    è\  û      BESSK1_V    ã^  <      BESSK1_V%LOG    _  <      BESSK1_V%EXP    [_  =      BESSK1_V%SQRT    _  =      BESSK1_V%SIZE    Õ_  <      BESSK1_V%ALL    `     a   BESSK1_V%X    `  Â       BESSK1_S    _a  <      BESSK1_S%LOG    a  <      BESSK1_S%EXP    ×a  =      BESSK1_S%SQRT    b  @   a   BESSK1_S%X    Tb  Õ       BESSI0_S    )c  <      BESSI0_S%ABS    ec  <      BESSI0_S%EXP    ¡c  =      BESSI0_S%SQRT    Þc  =      BESSI0_S%REAL    d  @   a   BESSI0_S%X    [d  ü      BESSI0_V    Wf  <      BESSI0_V%ABS    f  <      BESSI0_V%EXP    Ïf  =      BESSI0_V%SQRT    g  =      BESSI0_V%SIZE    Ig  =      BESSI0_V%REAL    g     a   BESSI0_V%X    h  Õ       BESSI1_S    çh  <      BESSI1_S%ABS    #i  <      BESSI1_S%EXP    _i  =      BESSI1_S%SQRT    i  =      BESSI1_S%REAL    Ùi  @   a   BESSI1_S%X    j  ü      BESSI1_V    l  <      BESSI1_V%ABS    Ql  <      BESSI1_V%EXP    l  =      BESSI1_V%SQRT    Êl  =      BESSI1_V%SIZE    m  =      BESSI1_V%REAL    Dm     a   BESSI1_V%X 