      SUBROUTINE reademis_netcdf(imonth1,imonth2,nlev,fname,varname,emis)

  !    ifort -L/usr/local/netCDF/lib -lnetcdf
  
      INCLUDE 'netcdf.inc'
 
      INTEGER :: imonth1, imonth2, nlev
      CHARACTER(len=70) :: fname,varname

      INTEGER :: i, ii, j, ilon, ilat, ilev, nlon, ngl, ncid, emisid, status
      PARAMETER (ilon=720, ilat=360, ilev=6)

      REAL :: emis(12,ilev,ilon,ilat)
      
      INTEGER :: nlondata(ilon),ngldata(ilat),nlevdata(ilev)

      ! LG- end

      ! define the grid size
      nlon = ilon   ! LG- normally 128
      ngl  = ilat   ! LG- normally 64
      nlev = ilev   ! #6 levels

      ! define the grid
      nlondata = (/ ((360. * REAL(i-1) / REAL(nlon)),i=1,nlon) /)
      ngldata  = (/ (90. - (180. * REAL(J-1) / REAL(ngl)),J=1,ngl) /) ! N -> S
      nlevdata = (/(i,i=1,nlev)/) 

      write(*,'(1a)')' Start input_netcdf, opening and reading netCDF file'
      
      ! LG- opening the input file to read, the open statement is different
      !     compared to that use to write a file!

      CALL nf(nf_open(fname,NF_NOWRITE,ncid))
      CALL nf(nf_inq_varid(ncid,varname,emisid))

      DO ii = imonth1,imonth2 ! begin time loop

         ! set the time (using the unit as defined above with nf90_put_att)
         time = REAL(ii-1)

         ! read data

         WRITE(*,'(1a,i2,1a,1a)'),' Reading emissions for month: ', ii,' param: ',varname
         CALL read_nc_file(ii,ncid,emisid,emis)

      END DO ! end time loop

      CALL close_nc_file(ncid)

      write(*,'(1a)')' End input_netcdf, closing netCDF file'

      RETURN
      END

      !*****************************************************************************

      SUBROUTINE nf(status) ! turns nf_* function into SUBROUTINE + checks status
        INTEGER :: status
        IF (status /= nf_noerr) THEN
           WRITE (*,*), 'netcdf error: ', nf_strerror(status)
           STOP 'stopped'
        ENDIF
      END SUBROUTINE nf

      !*****************************************************************************

      SUBROUTINE read_nc_file(timestep,ncid,emisid,emis)

        integer :: timestep,ncid,emisid,start2d(3),cnt2d(3)
        real    :: emis(12,6,720,360)

        ! syntax: nf90_get_var(ncid, varid, values, start, cnt)
        ! start:  start in netcdf variable
        ! cnt:    number of netcdf variable points
        ! values: starting point of the fortran variable
        start2d = (/ 1,    1,         timestep /)

        ! read 2d data

        CALL nf(nf_get_var_real(ncid,emisid,emis)) ! emis.

      END SUBROUTINE read_nc_file

      !*****************************************************************************

      SUBROUTINE close_nc_file(ncid)

        integer ncid

        CALL nf(nf_close(ncid))

      END SUBROUTINE close_nc_file
