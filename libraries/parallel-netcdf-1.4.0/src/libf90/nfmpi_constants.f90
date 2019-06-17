!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: nfmpi_constants.f90 1468 2013-10-26 16:53:18Z wkliao $
!
! This file lists the Fortran constants used by PnetCDF
!
!

!
! external netcdf data types:
!
      integer, parameter, public :: &
      nf_byte = 1, &
      nf_int1 = nf_byte, &
      nf_char = 2, &
      nf_short = 3, &
      nf_int2 = nf_short, &
      nf_int = 4, &
      nf_float = 5, &
      nf_real = nf_float, &
      nf_double = 6, &
      nf_ubyte = 7, &
      nf_ushort = 8, &
      nf_uint = 9, &
      nf_int64 = 10, &
      nf_uint64 = 11, &
      nf_string = 12

!
! default fill values:
!
      integer, parameter, public :: &
      nf_fill_byte = -127, &
      nf_fill_int1 = nf_fill_byte, &
      nf_fill_char = 0, &
      nf_fill_short = -32767, &
      nf_fill_int2 = nf_fill_short, &
      nf_fill_int = -2147483647

      real, parameter, public :: &
      nf_fill_float = 9.9692099683868690e+36, &
      nf_fill_real = nf_fill_float

      double precision, parameter, public :: &
      nf_fill_double = 9.9692099683868690e+36

      integer*8, parameter, public :: &
      nf_fill_int64 = -9223372036854775806_8

!
! mode flags for opening and creating a netcdf dataset:
!
      integer, parameter, public :: &
      nf_nowrite = 0, &
      nf_write = 1, &
      nf_clobber = 0, &
      nf_noclobber = 4, &
      nf_fill = 0, &
      nf_nofill = 256, &
      nf_lock = 1024, &
      nf_share = 2048, &
      nf_64bit_offset = 512, &
      nf_64bit_data = 16, &
      nf_32bit = 16777216, &
      nf_sizehint_default = 0, &
      nf_align_chunk = -1, &
      nf_format_classic = 1, &
      nf_format_64bit = 2, &
      nf_format_64bit_data = 5

!
! size argument for defining an unlimited dimension:
!
      integer, parameter, public :: &
      nf_unlimited = 0

!
! global attribute id:
!
      integer, parameter, public :: &
      nf_global = 0

!
! implementation limits:
!
      integer, parameter, public :: &
      nf_max_dims = 1024, &
      nf_max_attrs = 8192, &
      nf_max_vars = 8192, &
      nf_max_name = 256, &
      nf_max_var_dims = nf_max_dims

!
! error codes:
!
      integer, parameter, public :: &
      NF_NOERR        =  0, &  ! No Error
      NF2_ERR         = -1, &  ! Returned for all errors in the v2 API
      NF_EBADID       = -33, & ! Not a netcdf id
      NF_ENFILE       = -34, & ! Too many netcdfs open
      NF_EEXIST       = -35, & ! netcdf file exists and NF_NOCLOBBER
      NF_EINVAL       = -36, & ! Invalid Argument
      NF_EPERM        = -37, & ! Write to read only
      NF_ENOTINDEFINE = -38, & ! Operation not allowed in data mode
      NF_EINDEFINE    = -39, & ! Operation not allowed in define mode
      NF_EINVALCOORDS = -40, & ! Index exceeds dimension bound
      NF_EMAXDIMS     = -41, & ! NF_MAX_DIMS exceeded
      NF_ENAMEINUSE   = -42, & ! String match to name in use
      NF_ENOTATT      = -43, & ! Attribute not found
      NF_EMAXATTS     = -44, & ! NF_MAX_ATTRS exceeded
      NF_EBADTYPE     = -45, & ! Not a netcdf data type
      NF_EBADDIM      = -46, & ! Invalid dimension id or name
      NF_EUNLIMPOS    = -47, & ! NF_UNLIMITED in the wrong index
      NF_EMAXVARS     = -48, & ! NF_MAX_VARS exceeded
      NF_ENOTVAR      = -49, & ! Variable not found
      NF_EGLOBAL      = -50, & ! Action prohibited on NF_GLOBAL varid
      NF_ENOTNC       = -51, & ! Not a netcdf file
      NF_ESTS         = -52, & ! In Fortran, string too short
      NF_EMAXNAME     = -53, & ! NF_MAX_NAME exceeded
      NF_EUNLIMIT     = -54, & ! NF_UNLIMITED size already in use
      NF_ENORECVARS   = -55, & ! nc_rec op when there are no record vars
      NF_ECHAR        = -56, & ! Attempt to convert between text & numbers
      NF_EEDGE        = -57, & ! Edge+start exceeds dimension bound
      NF_ESTRIDE      = -58, & ! Illegal stride
      NF_EBADNAME     = -59, & ! Attribute or variable name contains illegal characters
      NF_ERANGE       = -60, & ! Math result not representable
      NF_ENOMEM       = -61, & ! Memory allocation (malloc) failure
      NF_EVARSIZE     = -62, & ! One or more variable sizes violate format constraints
      NF_EDIMSIZE     = -63, & ! Invalid dimension size
      NF_ETRUNC       = -64    ! File likely truncated or possibly corrupted

! netCDF-4 error codes (copied from netCDF release)
      integer, parameter, public :: &
      NF_EHDFERR      = -101, & ! Error at HDF5 layer. 
      NF_ECANTREAD    = -102, & ! Can't read. 
      NF_ECANTWRITE   = -103, & ! Can't write. 
      NF_ECANTCREATE  = -104, & ! Can't create. 
      NF_EFILEMETA    = -105, & ! Problem with file metadata. 
      NF_EDIMMETA     = -106, & ! Problem with dimension metadata. 
      NF_EATTMETA     = -107, & ! Problem with attribute metadata. 
      NF_EVARMETA     = -108, & ! Problem with variable metadata. 
      NF_ENOCOMPOUND  = -109, & ! Not a compound type. 
      NF_EATTEXISTS   = -110, & ! Attribute already exists. 
      NF_ENOTNC4      = -111, & ! Attempting netcdf-4 operation on netcdf-3 file.   
      NF_ESTRICTNC3   = -112, & ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      NF_ENOTNC3      = -113, & ! Attempting netcdf-3 operation on netcdf-4 file.   
      NF_ENOPAR       = -114, & ! Parallel operation on file opened for non-parallel access.   
      NF_EPARINIT     = -115, & ! Error initializing for parallel access.   
      NF_EBADGRPID    = -116, & ! Bad group ID.   
      NF_EBADTYPID    = -117, & ! Bad type ID.   
      NF_ETYPDEFINED  = -118, & ! Type has already been defined and may not be edited. 
      NF_EBADFIELD    = -119, & ! Bad field ID.   
      NF_EBADCLASS    = -120, & ! Bad class.   
      NF_EMAPTYPE     = -121, & ! Mapped access for atomic types only.   
      NF_ELATEFILL    = -122, & ! Attempt to define fill value when data already exists. 
      NF_ELATEDEF     = -123, & ! Attempt to define var properties, like deflate, after enddef.
      NF_EDIMSCALE    = -124, & ! Probem with HDF5 dimscales.
      NF_ENOGRP       = -125, & ! No group found.
      NF_ESTORAGE     = -126, & ! Can't specify both contiguous and chunking.
      NF_EBADCHUNK    = -127, & ! Bad chunksize.
      NF_ENOTBUILT    = -128, & ! Attempt to use feature that was not turned on when netCDF was built.
      NF_EDISKLESS    = -129    ! Error in using diskless  access.


! PnetCDF error codes start here
      integer, parameter, public :: &
      NF_ESMALL                   = -201, & ! size of off_t too small for format
      NF_ENOTINDEP                = -202, & ! Operation not allowed in collective data mode
      NF_EINDEP                   = -203, & ! Operation not allowed in independent data mode
      NF_EFILE                    = -204, & ! Unknown error in file operation
      NF_EREAD                    = -205, & ! Unknown error in reading file
      NF_EWRITE                   = -206, & ! Unknown error in writting to file
      NF_EOFILE                   = -207, & ! file open/creation failed
      NF_EMULTITYPES              = -208, & ! Multiple types used in memory data
      NF_EIOMISMATCH              = -209, & ! Input/Output data amount mismatch
      NF_ENEGATIVECNT             = -210, & ! Negative count is specified
      NF_EUNSPTETYPE              = -211, & ! Unsupported etype in memory MPI datatype
      NF_EINVAL_REQUEST           = -212, & ! invalid nonblocking request ID
      NF_EAINT_TOO_SMALL          = -213, & ! MPI_Aint not large enough to hold requested value
      NF_ENOTSUPPORT              = -214, & ! feature is not yet supported
      NF_ENULLBUF                 = -215, & ! trying to attach a NULL buffer
      NF_EPREVATTACHBUF           = -216, & ! previous attached buffer is found
      NF_ENULLABUF                = -217, & ! no attached buffer is found
      NF_EPENDINGBPUT             = -218, & ! pending bput is found, cannot detach buffer
      NF_EINSUFFBUF               = -219, & ! attached buffer is too small
      NF_ENOENT                   = -220, & ! File does not exist when calling nfmpi_open()
      NF_EINTOVERFLOW             = -221    ! Overflow when type cast to 4-byte integer

! header inconsistency error starts in -250
      integer, parameter, public :: &
      NF_EMULTIDEFINE             = -250, & ! NC definitions on multiprocesses conflict
      NF_EMULTIDEFINE_OMODE       = -251, & ! file create/open modes are inconsistent among processes
      NF_EMULTIDEFINE_DIM_NUM     = -252, & ! inconsistent number of dimensions
      NF_EMULTIDEFINE_DIM_SIZE    = -253, & ! inconsistent size of dimension
      NF_EMULTIDEFINE_DIM_NAME    = -254, & ! inconsistent dimension names
      NF_EMULTIDEFINE_VAR_NUM     = -255, & ! inconsistent number of variables
      NF_EMULTIDEFINE_VAR_NAME    = -256, & ! inconsistent variable name
      NF_EMULTIDEFINE_VAR_NDIMS   = -257, & ! inconsistent variable's number of dimensions
      NF_EMULTIDEFINE_VAR_DIMIDS  = -258, & ! inconsistent variable's dimid
      NF_EMULTIDEFINE_VAR_TYPE    = -259, & ! inconsistent variable's data type
      NF_EMULTIDEFINE_VAR_LEN     = -260, & ! inconsistent variable's size
      NF_EMULTIDEFINE_NUMRECS     = -261, & ! inconsistent number of records
      NF_EMULTIDEFINE_VAR_BEGIN   = -262, & ! inconsistent variable file begin offset (internal use)
      NF_EMULTIDEFINE_ATTR_NUM    = -263, & ! inconsistent number of attributes */
      NF_EMULTIDEFINE_ATTR_SIZE   = -264, & ! inconsistent memory space used by attribute (internal use)
      NF_EMULTIDEFINE_ATTR_NAME   = -265, & ! inconsistent attribute name
      NF_EMULTIDEFINE_ATTR_TYPE   = -266, & ! inconsistent attribute type
      NF_EMULTIDEFINE_ATTR_LEN    = -267, & ! inconsistent attribute length
      NF_EMULTIDEFINE_ATTR_VAL    = -268, & ! inconsistent attribute value
      NF_ECMODE                   = NF_EMULTIDEFINE_OMODE,      &
      NF_EDIMS_NELEMS_MULTIDEFINE = NF_EMULTIDEFINE_DIM_NUM,    &
      NF_EDIMS_SIZE_MULTIDEFINE   = NF_EMULTIDEFINE_DIM_SIZE,   &
      NF_EDIMS_NAME_MULTIDEFINE   = NF_EMULTIDEFINE_DIM_NAME,   &
      NF_EVARS_NELEMS_MULTIDEFINE = NF_EMULTIDEFINE_VAR_NUM,    &
      NF_EVARS_NAME_MULTIDEFINE   = NF_EMULTIDEFINE_VAR_NAME,   &
      NF_EVARS_NDIMS_MULTIDEFINE  = NF_EMULTIDEFINE_VAR_NDIMS,  &
      NF_EVARS_DIMIDS_MULTIDEFINE = NF_EMULTIDEFINE_VAR_DIMIDS, &
      NF_EVARS_TYPE_MULTIDEFINE   = NF_EMULTIDEFINE_VAR_TYPE,   &
      NF_EVARS_LEN_MULTIDEFINE    = NF_EMULTIDEFINE_VAR_LEN,    &
      NF_ENUMRECS_MULTIDEFINE     = NF_EMULTIDEFINE_NUMRECS,    &
      NF_EVARS_BEGIN_MULTIDEFINE  = NF_EMULTIDEFINE_VAR_BEGIN

! error handling modes:
!
      integer, parameter, public :: &
      nf_fatal = 1, &
      nf_verbose = 2

! now we can define nfmpi_unlimited properly
      integer(KIND=MPI_OFFSET_KIND), parameter, public :: &
      nfmpi_unlimited = 0

! NULL request for non-blocking I/O APIs
      integer, parameter, public :: &
      NF_REQ_NULL = -1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!
! functions in the fortran interface
!

!
! netcdf data types:
!
      integer, parameter, public :: &
      ncbyte = 1, &
      ncchar = 2, &
      ncshort = 3, &
      nclong = 4, &
      ncfloat = 5, &
      ncdouble = 6

!
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!
      integer, parameter, public :: &
      ncrdwr = 1, &     ! read/write, 0 => readonly
      nccreat = 2, &    ! in create phase, cleared by ncendef
      ncexcl = 4, &     ! on create destroy existing file
      ncindef = 8, &    ! in define mode, cleared by ncendef
      ncnsync = 16, &   ! synchronise numrecs on change (x'10')
      nchsync = 32, &   ! synchronise whole header on change (x'20')
      ncndirty = 64, &  ! numrecs has changed (x'40')
      nchdirty = 128, & ! header info has changed (x'80')
      ncfill = 0, &     ! prefill vars on endef and increase of record, the default behavior
      ncnofill = 256, & ! do not fill vars on endef and increase of record (x'100')
      nclink = 32768    ! isa link (x'8000')

!
!     'mode' arguments for nccreate and ncopen
!
      integer, parameter, public :: &
      ncnowrit = 0, &
      ncwrite = ncrdwr, &
      ncclob = nf_clobber, &
      ncnoclob = nf_noclobber

!
!     'size' argument to ncdimdef for an unlimited dimension
!
      integer, parameter, public :: &
      ncunlim = 0

!
!     attribute id to put/get a global attribute
!
      integer, parameter, public :: &
      ncglobal  = 0

!
!     advisory maximums:
!
      integer, parameter, public :: &
      maxncop = 32, &
      maxncdim = 100, &
      maxncatt = 2000, &
      maxncvar = 2000

!     not enforced
      integer, parameter, public :: &
      maxncnam = 128, &
      maxvdims = maxncdim

!
!     global netcdf error status variable
!     initialized in error.c
!
      integer, parameter, public :: &
      ncnoerr = nf_noerr, &         ! no error
      ncebadid = nf_ebadid, &       ! not a netcdf id
      ncenfile = -31, &             ! nc_syserr too many netcdfs open
      nceexist = nf_eexist, &       ! netcdf file exists && ncnoclob
      nceinval = nf_einval, &       ! invalid argument
      nceperm = nf_eperm, &         ! write to read only
      ncenotin = nf_enotindefine, & ! operation not allowed in data mode
      nceindef = nf_eindefine, &    ! operation not allowed in define mode
      ncecoord = nf_einvalcoords, & ! coordinates out of domain
      ncemaxds = nf_emaxdims, &     ! maxncdims exceeded
      ncename = nf_enameinuse, &    ! string match to name in use
      ncenoatt = nf_enotatt, &      ! attribute not found
      ncemaxat = nf_emaxatts, &     ! maxncattrs exceeded
      ncebadty = nf_ebadtype, &     ! not a netcdf data type
      ncebadd = nf_ebaddim, &       ! invalid dimension id
      nceunlim = nf_eunlimpos, &    ! ncunlimited in the wrong index
      ncemaxvs = nf_emaxvars, &     ! maxncvars exceeded
      ncenotvr = nf_enotvar, &      ! variable not found
      nceglob = nf_eglobal, &       ! action prohibited on ncglobal varid
      ncenotnc = nf_enotnc, &       ! not a netcdf file
      ncests = nf_ests, &
      ncentool = nf_emaxname, &
      ncfoobar = 32, &
      ncsyserr = -31

!
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!
      integer, parameter, public :: &
      ncfatal = 1, &
      ncverbos = 2

!
!     default fill values.  these must be the same as in the c interface.
!
      integer, parameter, public :: &
      filbyte = -127, &
      filchar = 0, &
      filshort = -32767, &
      fillong = -2147483647

      real, parameter, public :: &
      filfloat = 9.9692099683868690e+36

      double precision, parameter, public :: &
      fildoub = 9.9692099683868690e+36

