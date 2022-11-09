module bloc_commun
  implicit double precision (a-h,o-z)
  integer, parameter :: imax=122
  integer, parameter :: jmax=65
  integer, parameter :: jsepar=50
  integer, parameter :: jeq=28
  real, parameter :: spv=-1.e+32
  real, parameter :: spvMin=-1.5e+32
  real, parameter :: spvMax=-0.5e+32
  integer, parameter :: jmtt=60
  integer, parameter :: imtt=120
  real*4 wdata3D(imax,jmax,2) !2 parce qu'au pire c'est un vecteur (2 dimensions)
  real*4 wdatai3D(imtt,jmtt,2)
  REAL, PARAMETER :: xi1 = 28.500, dxi = 3.00
  REAL, PARAMETER :: yj1 =-79.500, dyj = 3.00
  integer, PARAMETER :: iberp =56 , ibera = 103
  real xlon1, ylat1, dlat, dlong

  !--DEFINITION OF CONSTANTS.
  !  pi    : pi
  !  radian: value of one radian in degrees.
  !  degre : value of one degre in radian .
  !  separ : facteur de separation entre spv / normal value.
  real*4, parameter :: pi     = 4.0d0 * atan(1.0d0)
  real*4, parameter :: radian = pi/180.0
  real*4, parameter :: degre  = 180.0/pi
  real*4, parameter :: untour = 360.0d0
  real*4, parameter :: epsil = 0.1d-9
  real*4, parameter :: zero = 0.0
  real*4, parameter ::  one  = 1.0
  real*4, parameter :: separ = 0.5 + epsil

  !--DEFINITION OF ROTATION ANGLES

  real*4, parameter :: alpha =    0.0
  real*4, parameter :: beta  = -111.0

  save
end module bloc_commun

subroutine check(status)
  USE NETCDF
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: status
  if(status /= nf90_noerr) then 
    write(*,*)"Error : ", trim(nf90_strerror(status))
    stop 
  end if
end subroutine check  


program OceanGrideChange
  USE NETCDF
  USE bloc_commun
  IMPLICIT NONE

  TYPE variable
    CHARACTER(nf90_max_name) :: name
    integer :: itype
    integer :: netcdfId
    integer :: OutnetcdfId
    integer, dimension(:), allocatable :: dimIndex
    integer, dimension(:), allocatable :: OutdimIndex
    integer, dimension(:), allocatable :: dimSize
    integer, dimension(:), allocatable :: dimStart
    integer :: nbdim
  END TYPE variable

  integer, parameter :: mx=120
  integer, parameter :: my=65
  integer, parameter :: mz=20
  real*4 wdatx(imax,jmax), wdaty(imax,jmax)
  real*4 valgu(imtt,jmtt), valgv(imtt,jmtt)
  real*4 :: ttlon(imtt)
  real*4 :: ttlat(jmtt)
  integer t, ii, n, k, l
  integer :: j, i, jmin
  integer :: returnval
  character(nf90_max_name) :: inputfile, outputfile, ifnocompress_char
  real ylon1,dylon,xlat1,dxlat

  integer nio0p, njo0p

  type(variable), dimension(:), allocatable :: listVariable, listVariable2
  double precision, dimension(:,:,:), allocatable :: Value3D
  double precision, dimension(:,:,:,:), allocatable :: Value4Du
  double precision, dimension(:,:,:,:), allocatable :: Value4Dv
  double precision, dimension(:,:,:,:), allocatable :: Value4Dalbq
  double precision, dimension(:,:), allocatable :: Valueh
  double precision, dimension(:), allocatable :: ValueVar

  !variable pour l'ouverture ecriture du netcdf
  integer :: intputID, outputID, outdimid, RecordDimID
  integer :: unlimDimID, nbDim, nbVar, nbAtt
  integer :: nbvarDim, dimSize, varID
  integer :: variableType, outvarid
  integer, dimension(nf90_max_var_dims) :: varDimID
  character(nf90_max_name) :: varName, dimName, attName
  integer, dimension(nf90_max_var_dims) :: CorrespTabDimID, InverseCorrespTabDimID
  integer :: nbOutDim, nbOutVar
  double precision, dimension(:), allocatable :: valueDbl
  integer :: totaltime, nbexistvariable
  integer :: deflate_level, ifnocompress_int
  logical :: ifnocompress
  deflate_level=1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!              VARIABLE LIST            !!!!
  !!!! Becarefull albq need to be at the top !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(listVariable(29))
  !averaged lead fraction
  listVariable(1)%name="albq"!(time, ptlat, ptlon) ;
  listVariable(1)%itype=1
  !averaged salinity
  listVariable(2)%name="salt"!(time, tdepth, ptlat, ptlon)
  listVariable(2)%itype=1
  !averaged zonal velocity component
  listVariable(3)%name="u"!(time, tdepth, pulat, pulon)
  listVariable(3)%itype=2
  !averaged meridional velocity component
  listVariable(4)%name="v"!(time, tdepth, pulat, pulon) ;
  listVariable(4)%itype=2
  !averaged vertical velocity component
  listVariable(5)%name="w"!(time, wdepth, ptlat, ptlon) ;
  listVariable(5)%itype=1
  !averaged zonal barotropic momentum
  listVariable(6)%name="ubar"!(time, pulat, pulon) ;
  listVariable(6)%itype=2
  !averaged meridional barotropic momentum
  listVariable(7)%name="vbar"!(time, pulat, pulon) ;
  listVariable(7)%itype=2
  !averaged sea surface height
  listVariable(8)%name="ssh"!(time, ptlat, ptlon) ;
  listVariable(8)%itype=1
  !averaged SST
  listVariable(9)%name="sst"!(time, ptlat, ptlon) ;
  listVariable(9)%itype=1
  !averaged sea surface salinity
  listVariable(10)%name="sss"!(time, ptlat, ptlon) ;
  listVariable(10)%itype=1
  !averaged surface heat flux
  listVariable(11)%name="shflx"!(time, ptlat, ptlon) ;
  listVariable(11)%itype=1
  !averaged surface freshwater flux
  listVariable(12)%name="sfflx"!(time, ptlat, ptlon) ;
  listVariable(12)%itype=1
  !averaged depth of ocean surface mixed layer
  listVariable(13)%name="zmix"!(time, ptlat, ptlon) ;
  listVariable(13)%itype=1
  !averaged depth of convection
  listVariable(14)%name="zcnv"!(time, ptlat, ptlon) ;
  listVariable(14)%itype=1
  !averaged G-M slope
  listVariable(15)%name="msl"!(time, ptlat, ptlon) ;
  listVariable(15)%itype=1
  !averaged ice thickness
  listVariable(16)%name="hice"!(time, ptlat, ptlon) ;
  listVariable(16)%itype=1
  !averaged ice production
  listVariable(17)%name="hicp"!(time, ptlat, ptlon) ;
  listVariable(17)%itype=1
  !averaged snow thickness
  listVariable(18)%name="hsn"!(time, ptlat, ptlon) ;
  listVariable(18)%itype=1
  !averaged snow precipitation
  listVariable(19)%name="snow"!(time, ptlat, ptlon) ;
  listVariable(19)%itype=1
  !averaged ice temperature
  listVariable(20)%name="tice"!(time, ptlat, ptlon) ;
  listVariable(20)%itype=1
  !averaged heat flux at ice base
  listVariable(21)%name="fb"!(time, ptlat, ptlon) ;
  listVariable(21)%itype=1
  !averaged zonal ice velocity
  listVariable(22)%name="uice"!(time, pulat, pulon) ;
  listVariable(22)%itype=2
  !averaged meridional ice velocity
  listVariable(23)%name="vice"!(time, pulat, pulon) ;
  listVariable(23)%itype=2
  !averaged zonal wind stress
  listVariable(24)%name="wsx"!(time, pulat, pulon) ;
  listVariable(24)%itype=1
  !averaged meridional wind stress
  listVariable(25)%name="wsy"!(time, pulat, pulon) ;
  listVariable(25)%itype=1
  !meridional overturning streamfunction
  listVariable(26)%name="moc"!(time, sfdepth, sflat, basidx) ;
  listVariable(26)%itype=0
  !meridional heat transport
  listVariable(27)%name="mht"!(time, sflat, basidx) ;
  listVariable(27)%itype=0
  !meridional salt transport
  listVariable(28)%name="mst"!(time, sflat, basidx) ;
  listVariable(28)%itype=0
  !averaged potential temperature
  listVariable(29)%name="temp"!(time, tdepth, ptlat, ptlon)
  listVariable(29)%itype=1

  !get argument for filename
  ifnocompress=.FALSE.
  call getarg(1,inputfile)
  call getarg(2,outputfile)
  call getarg(3,ifnocompress_char)
  read (ifnocompress_char,'(I1)') ifnocompress_int
  if(ifnocompress_int.eq.1) then
    ifnocompress=.TRUE.
  endif
  write(*,'(A,L)') "No compression = ",ifnocompress




  dlat=dyj
  dlong=dxi
  xlon1=xi1
  ylat1=yj1
  nio0p=imtt
  njo0p=jmtt
  call gridtt(ttlon,ttlat,imtt,jmtt)
  xlat1=ttlat(1)
  dxlat=ttlat(2)-ttlat(1)
  ylon1=ttlon(1)
  dylon=ttlon(2)-ttlon(1)





  !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! OPEN INPUT FILE !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!
  call check(nf90_open(inputfile, nf90_nowrite, intputID))
  call check(nf90_inquire(intputID, nbDim, nbVar,  unlimitedDimId = unlimDimID))

  !get totaltime
  call check(nf90_inquire(intputID, unlimitedDimId = RecordDimID))
  call check(nf90_inquire_dimension(intputID, RecordDimID, len = dimSize))
  totaltime=dimSize
  write(*,'(A,I2)') "Total time in the file = ", totaltime

  !compte combien de variable de la la liste existe dans le fichier input
  nbexistvariable=0
  do i=1, size(listVariable)
    if(nf90_inq_varid(intputID, listVariable(i)%name, varID).eq.nf90_noerr) nbexistvariable=nbexistvariable+1
  enddo

  !load variable architecture
  write(*,'(A)') "Variable find in the input file ( idx) name ):"
  allocate(listVariable2(nbexistvariable))
  j=1
  do i=1, size(listVariable)
    if(nf90_inq_varid(intputID, listVariable(i)%name, varID).eq.nf90_noerr) then
      listVariable2(j)%name=listVariable(i)%name
      listVariable2(j)%itype=listVariable(i)%itype
      call check(nf90_inquire_variable(intputID, varID, varName, ndims = nbvarDim, dimids = varDimID))
      if(varName.eq.listVariable2(j)%name) then
        listVariable2(j)%netcdfId=varID
        allocate(listVariable2(j)%dimIndex(nbvarDim))
        listVariable2(j)%dimIndex(:)=varDimID(1:nbvarDim)
        allocate(listVariable2(j)%dimSize(nbvarDim))
        allocate(listVariable2(j)%dimStart(nbvarDim))
        do k=1, nbvarDim
          call check(nf90_inquire_dimension(intputID, varDimID(k), len = dimSize))
          listVariable2(j)%dimSize(k)=dimSize
          listVariable2(j)%dimStart(k)=1
        enddo
        listVariable2(j)%nbdim=nbvarDim
        write(*,'(A3,I3,A2,A)') '   ', listVariable2(j)%netcdfId, ') ', trim(listVariable2(j)%name)
        !write(*,'(4(I3))')listVariable2(j)%dimIndex
        j=j+1
      endif
    endif
  enddo


  !constuire le tableau de correspondance d'index de dimension entre le fichier d'entree et de sortie
  write(*,'(A)',ADVANCE='NO') "Construct index table..."
  CorrespTabDimID=-1
  InverseCorrespTabDimID=0
  l=3 !commence à trois parce que 1 et 2 sont le nouveau lon et lat
  do i=1,size(listVariable2)
    if(listVariable(i)%itype.ne.0) then
      jmin=3
      InverseCorrespTabDimID(listVariable2(i)%dimIndex(1))=1
      InverseCorrespTabDimID(listVariable2(i)%dimIndex(2))=2
    else
      jmin=1
    endif
    do j=jmin, size(listVariable2(i)%dimIndex)
      k=1
      do while(CorrespTabDimID(k).ne.listVariable2(i)%dimIndex(j))
        k=k+1
        if(k.gt.nf90_max_var_dims) exit
      enddo
      if(k.eq.nf90_max_var_dims+1) then
        CorrespTabDimID(l)=listVariable2(i)%dimIndex(j)
        InverseCorrespTabDimID(listVariable2(i)%dimIndex(j))=l
        l=l+1
      endif
    enddo
  enddo
  nbOutDim=l-1
  nbOutVar=0
  write(*,*) "OK"


  write(*,'(A)') "Create output file header :"
  !write target file
  if(ifnocompress) then
    call check(nf90_create(trim(outputfile), NF90_CLASSIC_MODEL, outputID))
  else
    call check(nf90_create(trim(outputfile), nf90_hdf5, outputID))
  endif

  !define dimension, variable associated and attribut
  write(*,'(A)',ADVANCE='NO') "   Define dimensions and copy attributs..."
  call check(nf90_inquire(intputID, unlimitedDimId = RecordDimID))
  do i=1, nbOutDim
    !get name and size of dimension
    if(i.eq.1) then
      dimName="lon"
      dimSize=120
    elseif(i.eq.2) then
      dimName="lat"
      dimSize=60
    else
      call check(nf90_inquire_dimension(intputID, CorrespTabDimID(i), name = dimName, len = dimSize))
    endif

    !define dimension with separation between unlimited and normal dimension
    if(CorrespTabDimID(i).eq.RecordDimID) then
      call check(nf90_def_dim(outputID, dimName, NF90_UNLIMITED, outdimid))
    else
      call check(nf90_def_dim(outputID, dimName, dimSize, outdimid))
    endif

    !define variable associate to dimension and copy/create attribus
    if(i.eq.1) then
      if(ifnocompress) then
        call check(nf90_def_var(outputID, "lon", NF90_FLOAT, (/ 1 /), outvarid))
      else
        call check(nf90_def_var(outputID, "lon", NF90_FLOAT, (/ 1 /), outvarid, shuffle = .TRUE., deflate_level=deflate_level))
      endif
      call check(nf90_put_att(outputID, outvarid, "long_name", "longitude coordinate"))
      call check(nf90_put_att(outputID, outvarid, "standard_name", "longitude"))
      call check(nf90_put_att(outputID, outvarid, "units", "degrees_east"))
      call check(nf90_put_att(outputID, outvarid, "axis", "X"))
      nbOutVar=nbOutVar+1
    elseif(i.eq.2) then
      if(ifnocompress) then
        call check(nf90_def_var(outputID, "lat", NF90_FLOAT, (/ 2 /), outvarid))
      else
        call check(nf90_def_var(outputID, "lat", NF90_FLOAT, (/ 2 /), outvarid, shuffle = .TRUE., deflate_level=deflate_level))
      endif
      call check(nf90_put_att(outputID, outvarid, "long_name", "latitude coordinate"))
      call check(nf90_put_att(outputID, outvarid, "standard_name", "latitude"))
      call check(nf90_put_att(outputID, outvarid, "units", "degrees_north"))
      call check(nf90_put_att(outputID, outvarid, "axis", "Y"))
      nbOutVar=nbOutVar+1
    else
      call check(nf90_inq_varid(intputID, dimName, varID))
      call check(nf90_inquire_variable(intputID, varID, xtype = variableType, nAtts =nbAtt))
      if(ifnocompress) then
        call check(nf90_def_var(outputID, dimName, variableType, (/ i /), outvarid))
      else
        call check(nf90_def_var(outputID, dimName, variableType, (/ i /), outvarid, shuffle = .TRUE., deflate_level=deflate_level))
      endif
      do j=1, nbAtt
        call check(nf90_inq_attname(intputID, varID, j, attName))
        call check(nf90_copy_att(intputID, varID, attName, outputID, outvarid))
      enddo
      nbOutVar=nbOutVar+1
    endif
  enddo
  write(*,*) "OK"


  !define variable from the list
  write(*,'(A)',ADVANCE='NO') "   Define variables and copy attributs..."
  do i=1, size(listVariable2)
    allocate(listVariable2(i)%OutdimIndex(size(listVariable2(i)%dimIndex)))
    do j=1, size(listVariable2(i)%dimIndex)
      listVariable2(i)%OutdimIndex(j)=InverseCorrespTabDimID(listVariable2(i)%dimIndex(j))
    enddo

    call check(nf90_inquire_variable(intputID, listVariable2(i)%netcdfId, nAtts =nbAtt))
    if(ifnocompress) then
      call check(nf90_def_var(outputID, listVariable2(i)%name, NF90_DOUBLE, listVariable2(i)%OutdimIndex, outvarid))
    else
      call check(nf90_def_var(outputID, listVariable2(i)%name, NF90_DOUBLE, listVariable2(i)%OutdimIndex, outvarid, shuffle = .TRUE., deflate_level=deflate_level))
    endif
    listVariable2(i)%OutnetcdfId=outvarid
    !and copy attribus
    do j=1, nbAtt
      call check(nf90_inq_attname(intputID, listVariable2(i)%netcdfId, j, attName))
      call check(nf90_copy_att(intputID, listVariable2(i)%netcdfId, attName, outputID, listVariable2(i)%OutnetcdfId))
    enddo
  enddo
  nbOutVar=nbOutVar+size(listVariable2)
  write(*,*) "OK"


  !finish the configuration of the output file and starting put data
  write(*,'(A)',ADVANCE='NO') "   Close definition mode..."
  call check(nf90_enddef(outputID))
  write(*,*) "OK"

  !copy dimension data
  write(*,'(A)',ADVANCE='NO') "Copy dimensions in output file..."
  call check(nf90_inquire(intputID, unlimitedDimId = RecordDimID))
  do i=1, nbOutDim
    if(i.eq.1) then
      call check(nf90_put_var(outputID, i, ttlon, (/ 1 /), (/ 120 /)))
    elseif(i.eq.2) then
      call check(nf90_put_var(outputID, i, ttlat, (/ 1 /), (/ 60 /)))
    else
      call check(nf90_inquire_dimension(intputID, CorrespTabDimID(i), name = dimName, len = dimSize))
      call check(nf90_inq_varid(intputID, dimName, varID))
      call check(nf90_inquire_variable(intputID, varID, xtype = variableType, nAtts =nbAtt))
      do j=1, nbOutVar
        returnval=nf90_inq_varid(outputID, dimName, outvarid)
        if(outvarid.ne.-1) exit
      enddo
      allocate(valueDbl(dimSize))
      call check(nf90_get_var(intputID, varID, valueDbl, (/ 1 /), (/ dimSize /)))
      call check(nf90_put_var(outputID, outvarid, valueDbl, (/ 1 /), (/ dimSize /)))
      deallocate(valueDbl)
    endif
  enddo
  write(*,*) "OK"


  !copy data already interpolate
  write(*,'(A)',ADVANCE='NO') "Copy variable already interpolate..."
  do n=1, size(listVariable2)
    if(listVariable2(n)%itype.eq.0) then
      call check(nf90_inq_varid(intputID, listVariable2(n)%name, varID))
      call check(nf90_inquire_variable(intputID, listVariable2(n)%netcdfId, xtype = variableType, nAtts =nbAtt))
      listVariable2(n)%dimSize(listVariable2(n)%nbdim)=1
      allocate(valueDbl(product(listVariable2(n)%dimSize)))
      do t=1, totaltime
        listVariable2(n)%dimStart(listVariable2(n)%nbdim)=t
        call check(nf90_get_var(intputID, varID, valueDbl, listVariable2(n)%dimStart, listVariable2(n)%dimSize))
        call check(nf90_put_var(outputID, listVariable2(n)%OutnetcdfId, valueDbl, listVariable2(n)%dimStart, listVariable2(n)%dimSize))
      enddo
      deallocate(valueDbl)
    endif
  enddo
  write(*,*) "OK"


  !get h value for undef verification
  write(*,'(A)') "Get undef zone..."
  call check(nf90_inq_varid(intputID, "h", varID))
  allocate(ValueVar(120*65))
  call check(nf90_get_var(intputID, varID, ValueVar, start = (/1,1/), count = (/120,65/)))
  allocate(Valueh(120,65))
  Valueh=reshape(ValueVar, (/120,65/))
  deallocate(ValueVar)
  write(*,*) "OK"

  !Open data time by time and process
  write(*,'(A)') "Processing..."
  allocate(Value4Dalbq(120, 65, 1, 1))
  allocate(Value3D(120, 65, 1))
  do t=1,totaltime
    write(*,'(A,I5,A,I5,A)',ADVANCE='NO') "   t = ", t, '/', totaltime, ' '
    n=1
    do while(n.le.size(listVariable2))

      !affichage de progression
      if(listVariable2(n)%itype.eq.0) then
        write(*,'(A)',ADVANCE='NO') '*'
      elseif(listVariable2(n)%itype.eq.1) then
        write(*,'(A)',ADVANCE='NO') '.'
      else
        write(*,'(A)',ADVANCE='NO') '^'
      endif


      if (listVariable2(n)%itype.ne.0) then !n'interpole pas le type 0
        varName=listVariable2(n)%name

        !call CF_READ2D(TRIM(name, varName, tk, imax-2, jmax, 1, w1)
        !limite la dimension temporelle pour lire pas de temps par pas de temps
        listVariable2(n)%dimStart(listVariable2(n)%nbdim)=t
        listVariable2(n)%dimSize(listVariable2(n)%nbdim)=1

        !initialise le pointeur de donnee monodimensionnel
        allocate(ValueVar(product(listVariable2(n)%dimSize)))

        !lit la variable
        call check(nf90_get_var(intputID, listVariable2(n)%netcdfId, ValueVar, start = listVariable2(n)%dimStart, count = listVariable2(n)%dimSize))

        !verifie si 3D ou 4D et reshape en conséquense pour stoquer dans value4D
        if(listVariable2(n)%nbdim.eq.3) then
          allocate(Value4Du(120, 65, 1, 1))
          Value3D=reshape(ValueVar, (/120,65,1/))
          Value4Du(:,:,1,:)=Value3D(:,:,:)
        else
          allocate(Value4Du(120, 65, listVariable2(n)%dimSize(3), 1))
          Value4Du=reshape(ValueVar, (/120,65,listVariable2(n)%dimSize(3),1/))
        endif
        deallocate(ValueVar)

        !storage albq like reference variable for other
        if(varName.eq."albq") then
          Value4Dalbq=Value4Du
        endif

        !open also the next variable if it's a vector type
        if (listVariable2(n)%itype.eq.2) then
          varName=listVariable2(n+1)%name
          listVariable2(n+1)%dimStart(listVariable2(n+1)%nbdim)=t
          listVariable2(n+1)%dimSize(listVariable2(n+1)%nbdim)=1
          allocate(ValueVar(product(listVariable2(n+1)%dimSize)))
          call check(nf90_get_var(intputID, listVariable2(n+1)%netcdfId, ValueVar, start = listVariable2(n+1)%dimStart, count = listVariable2(n+1)%dimSize))

          if(listVariable2(n+1)%nbdim.eq.3) then
            allocate(Value4Dv(120, 65, 1, 1))
            Value3D=reshape(ValueVar, (/120,65,1/))
            Value4Dv(:,:,1,:)=Value3D(:,:,:)
          else
            allocate(Value4Dv(120, 65, listVariable2(n)%dimSize(3), 1))
            Value4Dv=reshape(ValueVar, (/120,65,listVariable2(n)%dimSize(3),1/))
          endif
          deallocate(ValueVar)
        endif

        !chaque profondeur est traité indépendament
        do k=1, size(Value4Du,3)

          !si pas assez de glace met à zero et verifie les valeurs undef
          do i=1,120
            do j=1,65
              if(Valueh(i,j).lt.-0.9e+32) then
                Value4Du(i,j,k,1)=spv
                if(allocated(Value4Dv)) Value4Dv(i,j,k,1)=spv
              elseif(Value4Dalbq(i,j,1,1).lt.0.05) then !! si pas assez glace on met à zero
                if( (varName.eq."hice").or.(varName.eq."hicp").or.(varName.eq."hsn").or.(varName.eq."snow").or.(varName.eq."tice").or.(varName.eq."uice").or.(varName.eq."vice") ) then
                  Value4Du(i,j,k,1)=0.0
                  if(allocated(Value4Dv)) Value4Dv(i,j,k,1)=0.0
                endif
              endif
            enddo
          enddo

          wdata3D(:,:,1)=Value4Du(:,:,k,1)
          wdata3D(121,:,1)=Value4Du(1,:,k,1)
          wdata3D(122,:,1)=Value4Du(2,:,k,1)
          if(allocated(Value4Dv)) then
            wdata3D(:,:,2)=Value4Dv(:,:,k,1)
            wdata3D(121,:,2)=Value4Dv(1,:,k,1)
            wdata3D(122,:,2)=Value4Dv(2,:,k,1)
          else
            wdata3D(:,:,2)=0.0
          endif


          !cyclic correspondance
          do j=2,jeq
            wdata3D(1,j,:) = wdata3D(imax-1,j,:) !1<-121 (grille 120 65)
            wdata3D(imax,j,:) = wdata3D(2,j,:)
            do ii=ibera-5,ibera+5
              wdata3D(ii,jmax,:) = spv
            enddo
            do ii=iberp-5,iberp+5
              wdata3D(ii,jsepar,:) = spv
            enddo
          enddo


          !
          ! Interpolation
          !
          if (listVariable2(n)%itype.eq.2) then
            do i=1,imax
              do j=1,jmax
                wdatx(i,j)=wdata3D(i,j,1) !composante 1 vecteur
                wdaty(i,j)=wdata3D(i,j,2) !composante 2 vecteur
              enddo
            enddo
            call mercatv(ttlon,ttlat,wdatx,wdaty,valgu,valgv) 
            do i=1,imtt
              do j=1,jmtt
                wdatai3D(i,j,1)=valgu(i,j)
                wdatai3D(i,j,2)=valgv(i,j)
              enddo
            enddo

            !Put interpolate data in output file
            if(listVariable2(n)%nbdim.eq.3) then
              call check(nf90_put_var(outputID, listVariable2(n)%OutnetcdfId, wdatai3D(:,:,1), (/1,1,t/),  (/imtt,jmtt,1/)))
              call check(nf90_put_var(outputID, listVariable2(n+1)%OutnetcdfId, wdatai3D(:,:,2), (/1,1,t/),  (/imtt,jmtt,1/)))
            else
              call check(nf90_put_var(outputID, listVariable2(n)%OutnetcdfId, wdatai3D(:,:,1), (/1,1,k,t/),  (/imtt,jmtt,1,1/)))
              call check(nf90_put_var(outputID, listVariable2(n+1)%OutnetcdfId, wdatai3D(:,:,2), (/1,1,k,t/),  (/imtt,jmtt,1,1/)))
            endif
          else
            call mercat(ttlon,ttlat)
            !Put interpolate data in output file
            if(listVariable2(n)%nbdim.eq.3) then
              call check(nf90_put_var(outputID, listVariable2(n)%OutnetcdfId, wdatai3D(:,:,1), (/1,1,t/),  (/imtt,jmtt,1/)))
            else
              call check(nf90_put_var(outputID, listVariable2(n)%OutnetcdfId, wdatai3D(:,:,1), (/1,1,k,t/),  (/imtt,jmtt,1,1/)))
            endif
          endif



        enddo
        deallocate(Value4Du)
        if(allocated(Value4Dv)) deallocate(Value4Dv)
      endif
      if (listVariable2(n)%itype.eq.2) then
        n=n+2
      else
        n=n+1
      endif
    enddo
    write(*,*)
  enddo


  !close output file
  call check(nf90_close(outputID))

  write(*,*)'End of OceanGrideChange'
  write(*,*)'--------------------------------------'
end program OceanGrideChange

subroutine mercat(ttlon,ttlat)
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !--INTERPOLATION OF SCALAR DATA PLACED AT THE CENTER OF THE
  !--ELEMENTS OF A TWO RECTANGULAR GRID ONTO ONE GRID.
  !--LONGITUDE-LATITUDE COORDINATES FOR BOTH GRIDS.
  !--DATA ARE OUTPUTS OF THE PROGRAM OTI.
  !--This subroutine is identical to mercat.f except for input and ouputs
  !
  !--M.A.Morales Maqueda, 11-IV-1994.
  !--modified by H.GOOSSE 15-IV-1994
  !--modified by H.GOOSSE + M.A.M. Maqueda 15-IV-1994
  !--modified by H.GOOSSE 16-V-1994
  !  modif : 14/10/94

  !--(alpha,beta): (latitude,longitude) of the north pole of the new grid.
  !
  USE bloc_commun

  integer, parameter :: nsmax = 2
  integer :: jm, i, j

  real*4 :: ttlon(imtt)
  real*4 :: ttlat(jmtt)

  integer :: gxw, iwp, jwp, nprt
  real*4 ::  xaj1, yai1, dxaj, dyai, dxw, dyw, xxx, yyy, du, dd, dr, dl
  real*4 :: dsxw, dsyw, dcxw, dcyw, dxa, dya, nn0, nn1, xw, yw, rd, ru, rr, rl, unsdtx, unsdty
  real*4 :: gwlon(0:imax), galat(0:imax)
  real*4 :: gwlat(0:jmax+1), galon(0:jmax+1)
  real*4 :: valad, valcd, valau, valc, valcu, valcl, valcr
  real*4 :: val(0:imax,0:jmax+1)
  real*4 :: whigri(imtt,jmtt)



  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

  1001  format(A32,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
  1000  format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
  1111  format(3(F7.2,1X,F7.3,1X),I3,A)

  jm=jmax


  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !  3 ) definition de la nouvelle grille .                              |
  !-----------------------------------------------------------------------

  !-----
  !--DEFINE INTERPOLATING GRID WW.
  !
  !     call gridtt(ttlon,ttlat,imtt,jmtt)

  !--DEFINE ORIGINAL GRIDS.
  !----
  xaj1 =  90. + yj1
  yai1 =  90. + beta + untour - xi1
  dxaj = dyj
  dyai = -dxi
  do i=0,imax
    gwlon(i) = xi1 + dxi * DFLOAT(i-1)
    galat(i) = 90. + beta + untour - gwlon(i)
  enddo
  do j=0,jmax+1
    gwlat(j) = yj1 + dyj * DFLOAT(j-1)
    galon(j) = 90. + gwlat(j)
  enddo

  !        write(6,*) 'galon :'
  !        write(6,'(20F6.1)') (galon(j),j=0,jmax+1)
  !        write(6,*) 'galat :'
  !        write(6,'(20F6.1)') (galat(i),i=0,imax)

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

  !--COOMPUTE DE CORRESPONDANCE BETWEEN GRIDS

  call choiceg(ttlat,ttlon,imtt,jmtt,whigri)

  !       open (15,file='choiceg.dat')
  !       do 350 j=1,jmtt
  ! !       write(15,'(122(F8.3))') (whigri(i,j),i=1,imtt)
  !         write(15,'(122(i1))') (int(whigri(i,j)),i=1,imtt)
  !  350  continue
  !       close (15)

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !  4 ) Traitement de la nouvelle grille colonne par colonne .          |
  !-----------------------------------------------------------------------

  !--MAIN DO-LOOP.

  do j=1,jmtt
    do i=1,imtt
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !--debut du traitement de la colonne (i,j) :
      if (whigri(i,j).eq.1.) then
        dxw = ttlon(i)
        dyw = ttlat(j)
        gxw = xi1 + mod(dxw-xi1+untour, untour)
        xxx = ( gxw - xi1 ) / dxi + 0.5
        iwp = nint(xxx) 
        iwp = max(0,min(imax-1,iwp))
        dr = gwlon(iwp+1) - gxw
        dl = gxw - gwlon(iwp)
        yyy = ( dyw - yj1 ) / dyj + 0.5
        jwp = nint(yyy) 
        jwp = max(1,min(jmax,jwp))
        du = gwlat(jwp+1) - dyw
        dd = dyw - gwlat(jwp)
      else
        !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        yw   = ttlat(j)
        dyw  = yw*radian
        dsyw = sin(dyw)
        dcyw = cos(dyw)
        xw   = mod(ttlon(i)-beta, untour)
        dxw  = xw*radian
        dsxw = sin(dxw)
        dcxw = cos(dxw)
        !--COMPUTE COORDINATES ON THE SS GRID OF A POINT OF THE AA GRID.
        dya = asin(dcyw*dcxw) * degre
        dxa = atan2(dcyw*dsxw,-dsyw) * degre
        dxa = mod(dxa+untour, untour)
        !---
        yyy = ( dxa - xaj1 ) / dxaj + 0.5
        jwp = nint(yyy) 
        jwp = max(0,min(jmax,jwp))
        du = galon(jwp+1) - dxa
        dd = dxa - galon(jwp)
        xxx = ( dya - yai1 ) / dyai + 0.5
        iwp = nint(xxx) 
        iwp = max(0,min(imax-1,iwp))
        dr = galat(iwp+1) - dya
        dl = dya - galat(iwp)

      endif
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      !--Pour verification :
      !        goto 550
      nn0=999999.99
      nn1=999999.99
      if (ttlon(i).ge.369. .or. ttlon(i).le.-1. ) then
        nprt = nprt + 1
        if (nprt.ge.nn0 .and.  nprt.le.nn1 ) then
          write(99,*) 'nprt, whigri(i,j) :', nprt, whigri(i,j)
          write(99,*) 'i,j, iwp, jwp :'
          write(99,*)  i,j, iwp, jwp
          write(99,*) 'dl, dr, dd, du :'
          write(99,*)  dl, dr, dd, du
          if ( whigri(i,j).eq.1.0d0) then
            write(99,*) 'dxw, dyw :'
            write(99,*)  dxw, dyw
          else
            write(99,*) 'dxa, dya :'
            write(99,*)  dxa, dya 
          endif
          !            write(99,*) ' ttlon, ttlat :', ttlon(i), ttlat(j)
          !            write(99,*) ' gwlon, gwlat :', gwlon(iwp), gwlat(jwp)
          !            write(99,*) ' galon, galat :', galon(jwp), galat(iwp)
        endif
      endif

      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !   5 ) Interpolation a partir des 4 voisins iwp/iwp+1,jwp/jwp+1 .     |
      !-----------------------------------------------------------------------

      !--POINT (i,j). A POINT IS CONSIDERED TO BE A LAND POINT 
      !--IF THE NEAREST DATA POINT IS A LAND POINT.

      unsdtx = 1.0 / (dl + dr)
      rl = dl * unsdtx
      rr = dr * unsdtx
      unsdty = 1.0 / (dd + du)
      rd = dd * unsdty
      ru = du * unsdty

      !--debut du traitement du point (i,j) :
      valcl = wdata3D(iwp,jwp,1)
      valcr = wdata3D(iwp+1,jwp,1)
      !if (valcl.eq.spv) then
      if ((valcl.gt.spvMin).and.(valcl.lt.spvMax).or.(valcl.eq.0.0)) then
        if (rl.le.0.5) then
          valad = valcl
        else
          valad = valcr
          valcd = valcr
        endif
      else
        if ((valcr.gt.spvMin).and.(valcr.lt.spvMax).or.(valcr.eq.0.0)) then
          !if (valcr.eq.spv) then
          if (rr.le.0.5) then
            valad = valcr
          else
            valad = valcl
            valcd = valcl
          endif
        else
          valad = valcl * rr + valcr * rl
          valcd = valad
        endif
      endif

      valcl = wdata3D(iwp,jwp+1,1)
      valcr = wdata3D(iwp+1,jwp+1,1)
      !if (valcl.eq.spv) then
      if ((valcl.gt.spvMin).and.(valcl.lt.spvMax).or.(valcl.eq.0.0)) then
        if (rl.le.0.5) then
          valau = valcl
        else
          valau = valcr
          valcu = valcr
        endif
      else
        if ((valcr.gt.spvMin).and.(valcr.lt.spvMax).or.(valcr.eq.0.0)) then
          !if (valcr.eq.spv) then
          if (rr.le.0.5) then
            valau = valcr
          else
            valau = valcl
            valcu = valcl
          endif
        else
          valau = valcl * rr + valcr * rl
          valcu = valau
        endif
      endif

      if ((valad.gt.spvMin).and.(valad.lt.spvMax).or.(valad.eq.0.0)) then
        !if (valad.eq.spv) then
        if (rd.le.0.5) then
          wdatai3D(i,j,1)  = spv
        else
          if ((valau.gt.spvMin).and.(valau.lt.spvMax).or.(valau.eq.0.0)) then
            !if (valau.eq.spv) then
            wdatai3D(i,j,1) = spv
          else
            valc        = valcu
            wdatai3D(i,j,1) = valcu
          endif
        endif
      else
        if ((valau.gt.spvMin).and.(valau.lt.spvMax).or.(valau.eq.0.0)) then
          !if (valau.eq.spv) then
          if (ru.le.0.5) then
            wdatai3D(i,j,1) = spv
          else
            valc        = valcd
            wdatai3D(i,j,1) = valcd
          endif
        else
          valc        = valcd * ru + valcu * rd
          wdatai3D(i,j,1) = valc
        endif
      endif


      !--Pour verification : 
      nn0=999999.99
      nn1=999999.99
      if (ttlon(i).ge.369. .or. ttlon(i).le.-1. ) then
        if (nprt.ge.nn0 .and.  nprt.le.nn1 ) then
          write(99,*) 'val(i,i+1,/j,j+1) ='
          write(99,'(4F10.4)') val(iwp,jwp), val(iwp+1,jwp), val(iwp,jwp+1), val(iwp+1,jwp+1)
          !            write(99,*) 'vala(i,j,1) =', vala(i,j,1)
          write(99,*) 'wdatai3D(i,j,1) =', wdatai3D(i,j,1)
        endif
      endif

    enddo
  enddo

  return

end

subroutine gridtt(ttlon,ttlat,imtt,jmtt)
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !
  !--DEFINE INTERPOLATING GRID TT.
  !--LATITUDES AND LONGITUDES CORRESPOND TO THE CENTERS OF THE GRID ELEMENTS.
  !
  implicit double precision (a-h,o-z)
  !
  integer :: imtt, jmtt,i, j
  real*4 :: xlong1, delx, dely
  real*4 :: ttlon(imtt),ttlat(jmtt)
  !     xlong1 = 23
  xlong1 = 0.0
  delx=360.0/real(imtt)
  dely=180.0/real(jmtt)
  do i=1,imtt
    !       ttlon(i)=xlong1+real(i-1)*delx+0.5*delx
    ttlon(i)=xlong1+real(i-1)*delx
    ttlon(i)=mod(ttlon(i),360.0d0)
  enddo
  do j=1,jmtt
    ttlat(j)=-90.+real(j-1)*dely+0.5*dely
  enddo
  return
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
end

subroutine choiceg(ttlat,ttlong,imtt,jmtt,whigri)
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !
  !               THIS ROUTINE DETERMINE AT WHICH ORIGINAL GRID (AA OR WW)
  !               EACH POINT OF TT CORRESPOND.
  !
  implicit double precision (a-h,o-z)
  !
  integer :: imtt, jmtt, i,j
  real*4 :: whigri(imtt,jmtt)
  real*4 :: ttlat(jmtt),ttlong(imtt)
  real*4 :: xoncri
  !     
  !     write(6,*) 'begining of choiceg'
  do i=1,imtt
    do j=1,jmtt
      whigri(i,j)=1.
    enddo
  enddo
  !
  !     write(6,*) 'after 10'
  do i=1,imtt
    do j=1,jmtt
      if (ttlat(j).gt.0.0.and.ttlat(j).le.8.0) then
        if ((ttlong(i).ge.290.).or.(ttlong(i).lt.30.)) then
          whigri(i,j)=2.
          !             write(10,*) ttlat(j),ttlong(i)
        endif
      endif
      if ((ttlat(j).gt.8.).and.(ttlat(j).le.10.)) then
        xoncri=281.+(ttlat(j)-8.)/(10.-8.)*(276.-281.)
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.30.) whigri(i,j)=2.
      endif
      if ((ttlat(j).gt.10.).and.(ttlat(j).le.15.)) then
        xoncri=276.+(ttlat(j)-10.)/(15.-10.)*(270.-276.)
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.30.) whigri(i,j)=2.
      endif
      if ((ttlat(j).gt.15.).and.(ttlat(j).le.20.)) then
        xoncri=270.+(ttlat(j)-15.)/(20.-15.)*(260.-270.)
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.30) whigri(i,j)=2.
      endif
      if ((ttlat(j).gt.20.).and.(ttlat(j).le.30.)) then
        xoncri=260.+(ttlat(j)-20.)/(30.-20.)*(260.-260.)
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.30.) whigri(i,j)=2.
      endif
      if ((ttlat(j).gt.30.).and.(ttlat(j).le.68.)) then
        xoncri=260.+(ttlat(j)-30.)/(65.-30.)*(260.-260.)
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.50.) whigri(i,j)=2.
      endif
      if ((ttlat(j).gt.67.).and.(ttlat(j).le.90.)) then
        xoncri=0
        if (ttlong(i).ge.xoncri.or.ttlong(i).lt.360.) whigri(i,j)=2.
      endif
      !
    enddo
  enddo
  !
  !     write(6,*) 'end of choiceg'
  return
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
end
  !
subroutine mercatv(ttlon,ttlat,wdatx,wdaty,valgu,valgv)
  !
  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !--INTERPOLATION OF SCALAR DATA PLACED AT THE CENTER OF THE
  !--ELEMENTS OF A TWO RECTANGULAR GRID ONTO ONE GRID.
  !--LONGITUDE-LATITUDE COORDINATES FOR BOTH GRIDS.
  !--DATA ARE OUTPUTS OF THE PROGRAM OTI.
  !
  !--M.A.Morales Maqueda, 11-IV-1994.
  !--modified by H.GOOSSE 15-IV-1994
  !--modified by H.GOOSSE + M.A.M. Maqueda 15-IV-1994
  !--modified by H.GOOSSE 16-V-1994
  !--modified by JMC 19/09/95, adapted for vector, (derived from "provect").
  !  modif : 21/09/95

  !--(alpha,beta): (latitude,longitude) of the north pole of the new grid.
  !
  USE bloc_commun

  integer, parameter :: nsmax = 2

  real*4 :: ttlon(imtt)
  real*4 :: ttlat(jmtt)

  real*4 :: galat(0:imax), gwlon(0:imax)
  real*4 :: gwlat(0:imax), galon(0:jmax+1)
  real*4 :: cxw(0:imax), sxw(0:imax), cyw(0:jmax+1), syw(0:jmax+1)
  real*4 :: cya(0:imax), sya(0:imax), cxa(0:jmax+1), sxa(0:jmax+1)

  real*4 :: wdatx(imax,jmax), wdaty(imax,jmax)
  real*4 :: valx(0:imax,0:jmax+1), valy(0:imax,0:jmax+1), valz(0:imax,0:jmax+1)

  real*4 :: valgu(imtt,jmtt), valgv(imtt,jmtt)
  real*4 :: cxt(imtt), sxt(imtt)
  real*4 :: cyt(jmtt), syt(jmtt)
  real*4 :: whigri(imtt,jmtt)

  integer :: im, jm, i, j, gxw, iwp, jwp, nprt, nncrv
  real*4 ::  xaj1, yai1, dxaj, dyai, dxw, dyw, xxx, yyy, du, dd, dr, dl, unsdtx, unsdty, valdw, valxd, valzd, valup, valxu, valyu, valyd, valzu
  real*4 :: dxa, dya, nn0, nn1, rd, ru, rr, rl, ylim, ylim1, ylim2, valg, vvx, vvy, vvz




  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

  1001 format(A32,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
  1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
  1111 format(3(F7.2,1X,F7.3,1X),I3,A)

  !--READ DATA.
  im=imax
  jm=jmax

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

  !--Initialisation :
  do j=0,jmax+1
    do i=0,imax
      !        valu(i,j) = spv
      !        valv(i,j) = spv
      valx(i,j) = 0.
      valy(i,j) = 0.
      valz(i,j) = 0.
    enddo
  enddo

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !  3 ) definition de la nouvelle grille .                              |
  !-----------------------------------------------------------------------

  !-----
  !--DEFINE INTERPOLATING GRID WW.
  !
  call gridtt(ttlon,ttlat,imtt,jmtt)
  do i=1,imtt
    cxt(i) = cos(radian*(ttlon(i)-beta))
    sxt(i) = sin(radian*(ttlon(i)-beta))
  enddo
  do j=1,jmtt
    cyt(j) = cos(radian*ttlat(j))
    syt(j) = sin(radian*ttlat(j))
  enddo

  !--DEFINE ORIGINAL GRIDS.
  !----
  !     xi1=xlon1
  !     dxi=dlong
  !     yj1=ylat1
  !     dyj=dlat

  xaj1 =  90. + yj1
  yai1 =  90. + beta + untour - xi1
  dxaj = dyj
  dyai = -dxi
  do i=0,imax
    gwlon(i) = xi1 + dxi * DFLOAT(i-1)
    galat(i) = 90. + beta + untour - gwlon(i)
    cxw(i) = cos(radian*(gwlon(i)-beta))
    sxw(i) = sin(radian*(gwlon(i)-beta))
    cya(i) = cos(radian*galat(i))
    sya(i) = sin(radian*galat(i))
  enddo
  do j=0,jmax+1
    gwlat(j) = yj1 + dyj * DFLOAT(j-1)
    galon(j) = 90. + gwlat(j)
    cyw(j) = cos(radian*gwlat(j))
    syw(j) = sin(radian*gwlat(j))
    cxa(j) = cos(radian*galon(j))
    sxa(j) = sin(radian*galon(j))
  enddo

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

  !--COOMPUTE DE CORRESPONDANCE BETWEEN GRIDS

  call choiceg(ttlat,ttlon,imtt,jmtt,whigri)

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !--calcul des 3 composantes dans le repere fixe .

  ylim1 = 69.0
  ylim2 = 294.0
  !- critere d'apartenace a Grille AA : y > min[ 69., max(0., 294. - x) ]
  do j=1,jm
    do i=1,im
      ylim = min(ylim1, max(zero, ylim2-gwlon(i)) )
      if (gwlat(j).le.ylim) then
        !- WW :
        valx(i,j) = -sxw(i)*wdatx(i,j)-syw(j)*cxw(i)*wdaty(i,j)
        valy(i,j) =  cxw(i)*wdatx(i,j)-syw(j)*sxw(i)*wdaty(i,j)
        valz(i,j) =  cyw(j)*wdaty(i,j)

      else
        !- AA :
        valz(i,j) =  sxa(j)*wdaty(i,j)-sya(i)*cxa(j)*wdatx(i,j)
        valy(i,j) =  cxa(j)*wdaty(i,j)+sya(i)*sxa(j)*wdatx(i,j)
        valx(i,j) = -cya(i)*wdatx(i,j)
      endif
    enddo
  enddo

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
  !  4 ) Traitement de la nouvelle grille colonne par colonne .          |
  !-----------------------------------------------------------------------

  !--MAIN DO-LOOP.

  do j=1,jmtt
    do i=1,imtt
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !--debut du traitement de la colonne (i,j) :
      if (whigri(i,j).eq.1.) then
        dxw = ttlon(i)
        dyw = ttlat(j)
        gxw = xi1 + mod(dxw-xi1+untour, untour)
        xxx = ( gxw - xi1 ) / dxi + 0.5
        iwp = nint(xxx)
        iwp = max(0,min(imax-1,iwp))
        dr = gwlon(iwp+1) - gxw
        dl = gxw - gwlon(iwp)
        yyy = ( dyw - yj1 ) / dyj + 0.5
        jwp = nint(yyy)
        jwp = max(0,min(jmax,jwp))
        du = gwlat(jwp+1) - dyw
        dd = dyw - gwlat(jwp)
      else
        !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        !--COMPUTE COORDINATES ON THE SS GRID OF A POINT OF THE AA GRID.
        dya = asin(cyt(j)*cxt(i)) * degre
        dxa = atan2(cyt(j)*sxt(i),-syt(j)) * degre
        dxa = mod(dxa+untour, untour)
        !---
        yyy = ( dxa - xaj1 ) / dxaj + 0.5
        jwp = nint(yyy)
        jwp = max(0,min(jmax,jwp))
        du = galon(jwp+1) - dxa
        dd = dxa - galon(jwp)
        xxx = ( dya - yai1 ) / dyai + 0.5
        iwp = nint(xxx)
        iwp = max(0,min(imax-1,iwp))
        dr = galat(iwp+1) - dya
        dl = dya - galat(iwp)

      endif
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      !--Pour verification :
      !        goto 550
      nn0=999999.99
      nn1=999999.99
      if (ttlon(i).ge.369. .or. ttlon(i).le.-1. ) then
        nprt = nprt + 1
        if (nprt.ge.nn0 .and.  nprt.le.nn1 ) then
          write(99,*) 'nprt, whigri(i,j) :', nprt, whigri(i,j)
          write(99,*) 'i,j, iwp, jwp :'
          write(99,*)  i,j, iwp, jwp
          write(99,*) 'dl, dr, dd, du :'
          write(99,*)  dl, dr, dd, du
          if ( whigri(i,j).eq.1.0d0) then
            write(99,*) 'dxw, dyw :'
            write(99,*)  dxw, dyw
          else
            write(99,*) 'dxa, dya :'
            write(99,*)  dxa, dya
          endif
          !            write(99,*) ' ttlon, ttlat :', ttlon(i), ttlat(j)
          !            write(99,*) ' gwlon, gwlat :', gwlon(iwp), gwlat(jwp)
          !            write(99,*) ' galon, galat :', galon(jwp), galat(iwp)
        endif
      endif

      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !   5 ) Interpolation a partir des 4 voisins iwp/iwp+1,jwp/jwp+1 .     |
      !-----------------------------------------------------------------------

      !--POINT (i,j). A POINT IS CONSIDERED TO BE A LAND POINT
      !--IF THE NEAREST DATA POINT IS A LAND POINT.

      unsdtx = 1.0 / (dl + dr)
      rl = dl * unsdtx
      rr = dr * unsdtx
      unsdty = 1.0 / (dd + du)
      rd = dd * unsdty
      ru = du * unsdty

      !--debut du traitement du point (i,j,k) :

      if ((wdatx(iwp,jwp).gt.spvMin).and.(wdatx(iwp,jwp).lt.spvMax).or.(wdatx(iwp,jwp).eq.0.0)) then
        !if (wdatx(iwp,jwp).eq.spv) then
        if (rl.le.separ) then
          valdw = spv
        else
          valdw = wdatx(iwp+1,jwp)
          valxd = valx(iwp+1,jwp)
          valyd = valy(iwp+1,jwp)
          valzd = valz(iwp+1,jwp)
        endif
      else
        if ((wdatx(iwp+1,jwp).gt.spvMin).and.(wdatx(iwp+1,jwp).lt.spvMax).or.(wdatx(iwp+1,jwp).eq.0.0)) then
          !if (wdatx(iwp+1,jwp).eq.spv) then
          if (rr.le.separ) then
            valdw = spv
          else
            valdw = wdatx(iwp,jwp)
            valxd = valx(iwp,jwp)
            valyd = valy(iwp,jwp)
            valzd = valz(iwp,jwp)
          endif
        else
          valdw = rr*wdatx(iwp,jwp)+rl*wdatx(iwp+1,jwp)
          valxd = rr*valx(iwp,jwp) + rl*valx(iwp+1,jwp)
          valyd = rr*valy(iwp,jwp) + rl*valy(iwp+1,jwp)
          valzd = rr*valz(iwp,jwp) + rl*valz(iwp+1,jwp)
        endif
      endif
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !
      if ((wdatx(iwp,jwp+1).gt.spvMin).and.(wdatx(iwp,jwp+1).lt.spvMax).or.(wdatx(iwp,jwp+1).eq.0.0)) then
        !if (wdatx(iwp,jwp+1).eq.spv) then
        if (rl.le.separ) then
          valup = spv
        else
          valup = wdatx(iwp+1,jwp+1)
          valxu = valx(iwp+1,jwp+1)
          valyu = valy(iwp+1,jwp+1)
          valzu = valz(iwp+1,jwp+1)
        endif
      else
        if ((wdatx(iwp+1,jwp+1).gt.spvMin).and.(wdatx(iwp+1,jwp+1).lt.spvMax).or.(wdatx(iwp+1,jwp+1).eq.0.0)) then
          !if (wdatx(iwp+1,jwp+1).eq.spv) then
          if (rr.le.separ) then
            valup = spv
          else
            valup = wdatx(iwp,jwp+1)
            valxu = valx(iwp,jwp+1)
            valyu = valy(iwp,jwp+1)
            valzu = valz(iwp,jwp+1)
          endif
        else
          valup = rr*wdatx(iwp,jwp+1)+rl*wdatx(iwp+1,jwp+1)
          valxu = rr*valx(iwp,jwp+1) + rl*valx(iwp+1,jwp+1)
          valyu = rr*valy(iwp,jwp+1) + rl*valy(iwp+1,jwp+1)
          valzu = rr*valz(iwp,jwp+1) + rl*valz(iwp+1,jwp+1)
        endif
      endif
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      !
      if ((valdw.gt.spvMin).and.(valdw.lt.spvMax).or.(valdw.eq.0.0)) then
        !if (valdw.eq.spv) then
        if (rd.le.separ) then
          valg  = spv
          valgu(i,j)  = spv
          valgv(i,j)  = spv
        else
          if ((valup.gt.spvMin).and.(valup.lt.spvMax).or.(valup.eq.0.0)) then
            !if (valup.eq.spv) then
            valg  = spv
            valgu(i,j)  = spv
            valgv(i,j)  = spv
          else
            valg = valup
            valgu(i,j) = -sxt(i)*valxu + cxt(i)*valyu
            valgv(i,j) =  cyt(j)*valzu - syt(j) * ( cxt(i)*valxu + sxt(i)*valyu )
            !               valgv(i,j,k) =  valzu / cyt(j)
          endif
        endif
      else
        !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if ((valup.gt.spvMin).and.(valup.lt.spvMax).or.(valup.eq.0.0)) then
          !if (valup.eq.spv) then
          if (ru.le.separ) then
            valg  = spv
            valgu(i,j)  = spv
            valgv(i,j)  = spv
          else
            valg = valdw
            valgu(i,j) = -sxt(i)*valxd + cxt(i)*valyd
            valgv(i,j) =  cyt(j)*valzd - syt(j) * ( cxt(i)*valxd + sxt(i)*valyd )
            !               valgv(i,j,k) =  valzd / cyt(j)
          endif
        else
          valg  = rd*valup + ru*valdw
          vvx = rd*valxu + ru*valxd
          vvy = rd*valyu + ru*valyd
          vvz = rd*valzu + ru*valzd
          valgu(i,j) = -sxt(i)*vvx + cxt(i)*vvy
          valgv(i,j) =  cyt(j)*vvz - syt(j) * ( cxt(i)*vvx + sxt(i)*vvy )
          !               valgv(i,j,k) =  vvz / cyt(j)
        endif
      endif
      !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      nncrv=1
      if (nncrv.eq.0) valgu(i,j) = valg

      !         if (whigri(i,j).eq.2) write(6,*) valgu(i,j,k)
      !--Pour verification :
      if (ttlon(i).ge.369. .or. ttlon(i).le.-1. ) then
        if (nprt.ge.nn0 .and.  nprt.le.nn1 ) then
          !            write(99,*) 'valu(i,i+1,/j,j+1) ='
          !            write(99,'(4F10.4)') valu(iwp,jwp,1), valu(iwp+1,jwp,1),
          !    &                          valu(iwp,jwp+1,1), valu(iwp+1,jwp+1,1)
          write(99,*) 'valgu(i,j) =', valgu(i,j)
        endif
      endif

      !--fin du traitement de la colonne (i,j) .
    enddo
  enddo
  !
  return

  !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
end
