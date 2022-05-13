      subroutine gridcheck(oronew,psnew,error,indexf)
***********************************************************************
*                                                                     * 
*             TRAJECTORY MODEL SUBROUTINE GRIDCHECK                   *
*                                                                     *
***********************************************************************
*                                                                     * 
*             AUTHOR:      G. WOTAWA                                  *
*             DATE:        1997-08-06                                 *
*                                                                     * 
*             Update:      1998-12, global fields allowed, A. Stohl   * 
*                                                                     * 
***********************************************************************
*                                                                     *
* DESCRIPTION:                                                        *
*                                                                     *
* THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
* LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
* ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
* GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
* WITHIN ONE FLEXTRA RUN) IS CHECKED IN THE ROUTINE "READWIND" AT ANY *
* CALL.                                                               *
*                                                                     *
* OUTPUT       error .true.   - can not read grid specifications      *
*              error .false.  - normal                                *
*              oronew .true.  - Terrain heights given in grib files   *
*              oronew .false. - Terrain heights not specified in the  *
*                               grib files (old file standard)        *
*                                                                     * 
* XLON0                geographical longitude of lower left gridpoint *
* XLAT0                geographical latitude of lower left gridpoint  *
* NX                   number of grid points x-direction              *
* NY                   number of grid points y-direction              *
* DX                   grid distance x-direction                      *
* DY                   grid distance y-direction                      *
* NUVZ                 number of grid points for horizontal wind      *
*                      components in z direction                      *
* NWZ                  number of grid points for vertical wind        *
* sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
*                      points of the polar stereographic grid):       *
*                      used to check the CFL criterion                *
*                      component in z direction                       *
* UVHEIGHT(1)-         heights of gridpoints where u and v are        *
* UVHEIGHT(NUVZ)       given                                          *
* WHEIGHT(1)-          heights of gridpoints where w is given         *
* WHEIGHT(NWZ)                                                        *
*                                                                     *
***********************************************************************
*
      include 'parecmwf.h'
      include 'comecmwf.h'

      integer i,ifn,ifield,j,k,iumax,iwmax,numskip,LOOPI
      
C LG- added 
      integer indexf      
      
      real xaux1,xaux2,yaux1,yaux2,sizesouth,sizenorth,xauxa
      logical error,oronew,psnew

* VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

C dimension of isec2 at least (22+n), where n is the number of parallels or
C meridians in a quasi-regular (reduced) Gaussian or lat/long grid

C dimension of zsec2 at least (10+nn), where nn is the number of vertical
C coordinate parameters

      integer isec0(2),isec1(56),isec2(22+nxmax+nymax),isec3(2)
      integer isec4(64),inbuff(jpack),ilen,iswap,ierr,iword,lunit
      real zsec2(60+2*nuvzmax),zsec3(2),zsec4(jpunp)
      character*1 yoper,opt
      data yoper/'D'/

      error=.false.
      oronew=.false.

C LG- added the switch for using the ecmwf surface pressure to reinitialize
C     the vertical column (see also the file containing the initial profiles).
C     Also added is the definition of a default surface pressure. This parameter
C     was originally initialized in parecmwf.h.

      psnew=.false.
      p0=101325.

C LG- end

      iumax=0
      iwmax=0
*
      if(ideltas.gt.0) then
        ifn=1
      else
        ifn=numbwf
      endif
	
*
* OPENING OF DATA FILE (GRIB CODE)
*
C LG- slightly modified to deal with reading in files twice

5     call pbopen(lunit,path(3)(1:len(3))//wfname(indexf),'r',ierr)
      if(ierr.lt.0) goto 999

      ifield=0   
10    ifield=ifield+1
*
* GET NEXT FIELDS
*
      call pbgrib(lunit,inbuff,jpack,ilen,ierr)

      if(ierr.eq.-1) goto 30    ! EOF DETECTED
      if(ierr.lt.-1) goto 999   ! ERROR DETECTED
 
      ierr=1

C Check whether we are on a little endian or on a big endian computer
*********************************************************************

c      if (inbuff(1).eq.1112101447) then         ! little endian, swap bytes
c        iswap=1+ilen/4
c        call swap32(inbuff,iswap)
c      else if (inbuff(1).ne.1196575042) then    ! big endian
c        stop 'subroutine gridcheck: corrupt GRIB data'
c      endif

!      print *,'start gribex'

      call gribex(isec0,isec1,isec2,zsec2,isec3,zsec3,isec4,
     &            zsec4,jpunp,inbuff,jpack,iword,yoper,ierr)

!      print *,'end gribex',jpunp,inbuff,jpack,iword,yoper,ierr
!      read (*,*)
      
! ESS_lg_20071007+ when an error occurs with gribex, check carefully if the 
!     path to the grib features is listed in the .bashrc file

      if (ierr.ne.0) goto 999   ! ERROR DETECTED

!      print *,'ifield',ifield

      if(ifield.eq.1) then
        nxfield=isec2(2)
        ny=isec2(3)

!        print *,'isec2',isec2(2),isec2(3),isec2(5),isec2(8)
!        read (*,*)

        xaux1=float(isec2(5))/1000.
        xaux2=float(isec2(8))/1000.
        if(xaux1.gt.180) xaux1=xaux1-360.0
        if(xaux2.gt.180) xaux2=xaux2-360.0
        if(xaux1.lt.-180) xaux1=xaux1+360.0
        if(xaux2.lt.-180) xaux2=xaux2+360.0
        if (xaux2.lt.xaux1) xaux2=xaux2+360.
        yaux1=float(isec2(7))/1000.
        yaux2=float(isec2(4))/1000.
        xlon0=xaux1
        ylat0=yaux1
        dx=(xaux2-xaux1)/float(nxfield-1)
        dy=(yaux2-yaux1)/float(ny-1)
        nlev_ec=isec2(12)/2-1

C Check whether fields are global
C If they contain the poles, specify polar stereographic map 
C projections using the stlmbr- and stcm2p-calls
************************************************************

        xauxa=abs(xaux2+dx-360.-xaux1)
        if (xauxa.lt.0.001) then
          nx=nxfield+1                 ! field is cyclic
          xglobal=.true.
        else
          nx=nxfield
          xglobal=.false.
        endif
        xauxa=abs(yaux1+90.)

        if (xglobal.and.xauxa.lt.0.001) then
          sglobal=.true.               ! field contains south pole
C Enhance the map scale by factor 3 (*2=6) compared to north-south
C map scale
          sizesouth=6.*(switchsouth+90.)/dy

C LG-     call outcommented (only relevant for global domain or
C         grids covering the poles

C           call stlmbr(southpolemap,-90.,0.)
C           call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth,
C      +    sizesouth,switchsouth,180.)

C LG-     end

          switchsouthg=(switchsouth-ylat0)/dy
        else
          sglobal=.false.
          switchsouthg=999999.
        endif
        xauxa=abs(yaux2-90.)
        if (xglobal.and.xauxa.lt.0.001) then
          nglobal=.true.               ! field contains north pole
C Enhance the map scale by factor 3 (*2=6) compared to north-south
C map scale
          sizenorth=6.*(90.-switchnorth)/dy

C LG-     call outcommented (only relevant for global domain or
C         grids covering the poles

C           call stlmbr(northpolemap,90.,0.)
C           call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth,
C      +    sizenorth,switchnorth,180.)

C LG-     end

          switchnorthg=(switchnorth-ylat0)/dy
        else
          nglobal=.false.
          switchnorthg=999999.
        endif
      endif

!      print *,'isec(6)',isec1(6)

      if(isec1(6).eq.129) oronew=.true.
      k=isec1(8)
      if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
      if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1) 

      if(isec1(6).eq.129) then
        do 20 j=0,ny-1

!          print *,'129',j,ny,nxfield,ga,xglobal

          do 21 i=0,nxfield-1
21          oro(i,j)=zsec4(nxfield*(ny-j-1)+i+1)/ga
          if (xglobal) oro(nx-1,j)=oro(0,j)
20        continue

!          print *,'129, endif'
c	  read (*,*)
      endif

C LG- assigning the surface pressure to the parameter p0

      if(isec1(6).eq.134) then
        psnew=.true.
        do 23 j=0,ny-1
          do 22 i=0,nxfield-1
22          ps(i,j,1,1)=zsec4(nxfield*(ny-j-1)+i+1)
23        continue
      endif

C LG- end

!      print *,'end loop 10'

      goto 10                      !! READ NEXT LEVEL OR PARAMETER
*
* CLOSING OF INPUT DATA FILE
*

!      print *,'start pblcose'
c      read (*,*)

30    call pbclose(lunit,ierr)     !! FINNISHED READING / CLOSING GRIB FILE

C LG- writing some messages concerning the orography and surface pressure
C     of the ECMWF data

      IF (oronew .AND. indexf .EQ. 1) THEN
       WRITE(*,'(1a,f7.2,1a)')
     &   ' gridcheck.f, the orography (height) of the ECMWF grid square is: ',
     &     oro(1,1),' m'
       WRITE(*,'(1a)')
     &   ' ENTER to continue'
       READ (*,*)
      ENDIF

      IF (psnew .AND. indexf .EQ. 1) THEN
       WRITE(*,'(1a,f10.2,1a)')
     &   ' gridcheck.f, psnew=.TRUE., p0 is set to a value: ',
     &     ps(1,1,1,1),' Pa'
       p0=ps(1,1,1,1)
       WRITE(*,'(1a)')
     &   ' ENTER to continue'
       READ (*,*)
      ENDIF

C LG- end

      nuvz=iumax
      nwz =iwmax
      if(nuvz.eq.nlev_ec) nwz=nlev_ec+1

      if (nx.gt.nxmax) then                         
	write(*,*) nx
      write(*,*) 'Error (gridcheck.f): Too many grid points in x direction.'
      write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC).'
      error=.true.
      return
      endif

      if (ny.gt.nymax) then                         
      write(*,*) 'Error (gridcheck.f): Too many grid points in y direction.'
      write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC).'
      error=.true.
      return
      endif

      if (nuvz.gt.nuvzmax) then                         
      write(*,*) 'Error (gridcheck.f): Too many u,v grid points in z '//
     +'direction.'
      write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC).'
      error=.true.
      return
      endif

      if (nwz.gt.nwzmax) then                         
      write(*,*) 'Error (gridcheck.f): Too many w grid points in z '//
     +'direction.'
      write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC).'
      error=.true.
      return
      endif

C LG- some additional warnings, introduced based on own experiences

      if (dx.eq.0.) then                         
      write(*,*) 'Error (gridcheck.f): input files not available'
      write(*,*) 'Check input directory (.Z ext) or defined pathway'
      write(*,*) '(see ecmwf_pathnames in case/input/ECMWF/)'
      error=.true.
      return
      endif

C Output of grid info
*********************

      if (indexf .EQ. 1) then
        write(*,*)
        write(*,*)
        write(*,'(a,2i7)') '# of vertical levels: ',nuvz,nwz
        write(*,*)
        write(*,'(a)') ' Domain of ECMWF input data:'
        write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ',
     +  xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
        write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range:  ',
     +  ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy

C LG- enter to continue

        print *,' ENTER to continue'
        read (*,*)
      endif

C LG- end

C Compute often used aux variables to convert geografical into grid coord.
***************************************************************************

      xthelp=180./pi/r_earth/dx
      ythelp=180./pi/r_earth/dy

C Compute grid distances in radians
***********************************

      dxradinv = 1./( dx * pi / 180. )
      dyradinv = 1./( dy * pi / 180. )

* CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
* PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM
* wheight(i) IS THE HEIGHT OF THE i-th MODEL HALF LEVEL (=INTERFACE BETWEEN
* 2 MODEL LEVELS) IN THE "ETA" SYSTEM

      numskip=nlev_ec-nuvz  ! number of ecmwf model layers not used
                            ! by trajectory model
      do 40 i=1,nwz
        j=numskip+i
        k=nlev_ec+1+numskip+i
        akm(nwz-i+1)=zsec2(10+j)
        bkm(nwz-i+1)=zsec2(10+k)

C ==========================================================================
C LG-   implementing this program in the SCM, it turned out that an error
C       occured here due to the fact that the parameter zsec2 had weird 
C       values, resulting in unrealistic values of akm and bkm, which are
C       used to determine the vertical coordinate system of the ECMWF data.
C       It turned out that this was due to the use of the compilitation
C       option -r8 (double precision data), whereas the option -r4 should
C       be used or the libary emosR64. Laurens Ganzeveld, December, 2001)
C =========================================================================

        if (akm(nwz-i+1).gt.1.e20) then                         
         write(*,*) 'Error reading grib data, check compilation options'
         write(*,*) 'in Makefile the option -r4 instead of -r8 must be used'
         write(*,*) 'Stop called in gridcheck.f'
         STOP
	endif

C LG-   end

C LG-   the height is given here in a relative value between 0-1, in contrast
C       to the height of the SCM. This is done by dividing by the pressure
C       instead of multiplying by the pressure (normally, the b term is
C       multiplied with the surface pressure, dividing by p0 gives the 
C       equation below)

40      wheight(nwz-i+1)=akm(nwz-i+1)/p0+bkm(nwz-i+1)

* CALCULATION OF uvheight, akz, bkz
* akz,bkz ARE THE DISCRETIZATION PARAMETERS FOR THE MODEL LEVELS
* uvheight(i) IS THE HEIGHT OF THE i-th MODEL LEVEL IN THE "ETA" SYSTEM

      do 45 i=1,nuvz
        uvheight(i)=0.5*(wheight(i+1)+wheight(i))
        akz(i)=0.5*(akm(i+1)+akm(i))
        bkz(i)=0.5*(bkm(i+1)+bkm(i))
45      continue

C If vertical coordinates decrease with increasing altitude, multiply by -1.
C This means that also the vertical velocities have to be multiplied by -1.
****************************************************************************

      if (uvheight(1).lt.uvheight(nuvz)) then
        zdirect=1.
      else
        zdirect=-1.
        do 55 i=1,nuvz
55        uvheight(i)=zdirect*uvheight(i)
        do 65 i=1,nwz
65        wheight(i)=zdirect*wheight(i)
      endif

C Compute minimum and maximum height of modelling domain
********************************************************

      heightmin=max(uvheight(1),wheight(1))
      heightmax=min(uvheight(nuvz),wheight(nwz))

      return

999   write(*,*)  
      write(*,*) ' ###########################################'//
     &           '###### '
      write(*,*) ' MODEL SUBROUTINE GRIDCHECK CAN NOT OPEN'
      write(*,*) ' INPUT DATA FILE: '//wfname(ifn)
      write(*,*) ' ###########################################'//
     &           '###### '
      write(*,*)
      write(*,'(a)') '!!! PLEASE INSERT A NEW CD-ROM AND   !!!'
      write(*,'(a)') '!!! PRESS ANY KEY TO CONTINUE...     !!!'
      write(*,'(a)') '!!! ...OR TERMINATE FLEXTRA PRESSING !!!'
      write(*,'(a)') '!!! THE "X" KEY...                   !!!'
      write(*,*)
      read(*,'(a)') opt
      if(opt.eq.'X') then
        error=.true.
      else
        goto 5
      endif

      return
      end
