dnl
dnl  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
dnl  See COPYRIGHT notice in top-level directory.
dnl
dnl $Id: test_iget.m4 1484 2013-11-02 14:57:28Z wkliao $
dnl

divert(-1)

dnl This is m4 source.
dnl Process using m4 to produce FORTRAN language file.

changequote([,])

undefine([index])dnl

dnl Macros

dnl Upcase(str)
dnl
define([Upcase],[dnl
translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl NFT_ITYPE(type)
dnl
define([NFT_ITYPE], [NFT_[]Upcase($1)])

dnl ARITH_VAR1(itype, value)
dnl
define([ARITH_VAR1], [ifelse($1, text, ichar($2), $2)])

dnl ARITH3(itype, value)
dnl
define([ARITH3], [ifelse($1, text, ichar($2($3:$3)), $2($3))])

dnl  DATATYPE(funf_suffix)
dnl
define([DATATYPE], [dnl
ifelse($1, text, character(len=MAX_NELS) $2,
ifelse($1, int1, NF_INT1_T $2$3,
ifelse($1, int2, NF_INT2_T $2$3,
ifelse($1, int, integer $2$3,
ifelse($1, int8, NF_INT8_T $2$3,
ifelse($1, real, real $2$3,
ifelse($1, double, doubleprecision $2$3)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
])

dnl  DATATYPE_VAR1(funf_suffix)
dnl
define([DATATYPE_VAR1], [dnl
ifelse($1, text, character $2,
ifelse($1, int1, NF_INT1_T $2,
ifelse($1, int2, NF_INT2_T $2,
ifelse($1, int, integer $2,
ifelse($1, int8, NF_INT8_T $2,
ifelse($1, real, real $2,
ifelse($1, double, doubleprecision $2)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
)[]dnl
])

dnl TEST_NFMPI_IGET_VAR1(TYPE)
dnl
define([TEST_NFMPI_IGET_VAR1],[dnl
        subroutine test_nf90mpi_iget_var1_$1()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer index2indexes
        double precision hash4
        logical equal, inRange3, in_internal_range

        integer ncid
        integer i
        integer j
        integer err
        integer nok      
        integer(kind=MPI_OFFSET_KIND) index(MAX_RANK)
        doubleprecision expect
        logical canConvert     
        DATATYPE_VAR1($1, value)
        doubleprecision val
        integer err_w, reqid(1), st(1)

        nok = 0

        err = nf90mpi_open(comm, testfile, NF90_NOWRITE, info, & 
                         ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_open: ', err)
        do 1, i = 1, numVars
            canConvert = (var_type(i) .eq. NF90_CHAR) .eqv. &
                         (NFT_ITYPE($1) .eq. NFT_TEXT)
            do 2, j = 1, var_rank(i)
                index(j) = 1
2           continue
            err = nf90mpi_iget_var(BAD_ID,i,value, reqid(1),index)
            if (err .ne. NF90_EBADID) &
                call errore('bad ncid: ', err)
            err = nf90mpi_iget_var(ncid,BAD_VARID, value, reqid(1), index)
            if (err .ne. NF90_ENOTVAR) &
                call errore('bad var id: ', err)
            do 3, j = 1, var_rank(i)
                index(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid,i,value, reqid(1),index)
                if (err .ne. NF90_EINVALCOORDS) &
                    call errore('bad index: ', err)
                index(j) = 1
3           continue
            do 4, j = 1, var_nels(i)
                err = index2indexes(j, var_rank(i), var_shape(1,i),  &
                                    index)
                if (err .ne. NF90_NOERR) &
                    call error('error in index2indexes 1')
                expect = hash4( var_type(i), var_rank(i), index,  &
                                NFT_ITYPE($1) )
                err = nf90mpi_iget_var(ncid,i,value, reqid(1),index)
                if (err .eq. NF90_NOERR) &
                    err_w = nf90mpi_wait_all(ncid, 1, reqid, st)
                if (canConvert) then
                    if (inRange3(expect,var_type(i),  &
                                 NFT_ITYPE($1))) then
                        if (in_internal_range(NFT_ITYPE($1), &
                                              expect)) then
                            if (st(1) .ne. 0) then
                                call errore('nf90mpi_iget_var: ',st(1))
                            else
                                val = ARITH_VAR1($1, value)
                                if (.not. equal(val, expect,  &
                                                var_type(i),  &
                                                NFT_ITYPE($1))) then
                                    call errord('unexpected: ', val)
                                else
                                    nok = nok + 1
                                end if
                            end if
                        else
                            if (st(1) .ne. NF90_ERANGE) &
                                call errore('Range error: ', st(1))
                        end if
                    else
                        if (st(1) .ne. 0  .and. st(1) .ne. NF90_ERANGE) &
                            call errore('OK or Range error: ', st(1))
                    end if
                else
                    if (err .ne. NF90_ECHAR) &
                        call errore('wrong type: ', err)
                end if
4           continue
1       continue
        err = nf90mpi_close(ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_close: ',  err)
        call print_nok(nok)
        end
])

dnl TEST_NFMPI_IGET_VAR(TYPE)
dnl
define([TEST_NFMPI_IGET_VAR],[dnl
        subroutine test_nf90mpi_iget_var_$1()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer index2indexes
        double precision hash4
        logical equal, inRange3, in_internal_range

        integer ncid
        integer i
        integer j
        integer err
        logical allInExtRange   
        logical allInIntRange   
        integer nels
        integer nok      
        integer(kind=MPI_OFFSET_KIND) index(MAX_RANK)
        doubleprecision expect(MAX_NELS)
        logical canConvert     
        DATATYPE($1, value, (MAX_NELS))
        doubleprecision val
        integer err_w, reqid(1), st(1)      

        nok = 0

        err = nf90mpi_open(comm, testfile, NF90_NOWRITE, info,  &
                         ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_open: ', err)
        do 1, i = 1, numVars
            canConvert = (var_type(i) .eq. NF90_CHAR) .eqv. &
                         (NFT_ITYPE($1) .eq. NFT_TEXT)
            err = nf90mpi_iget_var(BAD_ID, i, value,reqid(1))
            if (err .ne. NF90_EBADID) &
                call errore('bad ncid: ', err)
            err = nf90mpi_iget_var(ncid, BAD_VARID, value,reqid(1))
            if (err .ne. NF90_ENOTVAR) &
                call errore('bad var id: ', err)
            nels = 1
            do 3, j = 1, var_rank(i)
                nels = nels * var_shape(j,i)
3           continue
            allInExtRange = .true.
            allInIntRange = .true.
            do 4, j = 1, var_nels(i)
                err = index2indexes(j, var_rank(i), var_shape(1,i),  &
                                    index)
                if (err .ne. NF90_NOERR) &
                    call error('error in index2indexes 1')
                expect(j) = hash4( var_type(i), var_rank(i), index,  &
                                NFT_ITYPE($1) )
                if (inRange3(expect(j),var_type(i), NFT_ITYPE($1))) then
                    allInIntRange = allInIntRange .and. &
                        in_internal_range(NFT_ITYPE($1), expect(j))
                else
                    allInExtRange = .false.
                end if
4           continue
            err = nf90mpi_iget_var(ncid, i, value,reqid(1), count=var_shape(:,i))
            if (err .eq. NF90_NOERR) &
                err_w = nf90mpi_wait_all(ncid,1,reqid,st)
            if (canConvert) then
                if (allInExtRange) then
                    if (allInIntRange) then
                        if (st(1) .ne. 0)  &
                            call errore('nf90mpi_iget_var: ', st(1))
                    else
                        if (st(1) .ne. NF90_ERANGE) &
                            call errore('Range error: ', st(1))
                    endif
                else
                    if (st(1) .ne. 0  .and. st(1) .ne. NF90_ERANGE) &
                        call errore('Range error: ', st(1))
                endif
                do 5, j = 1, var_nels(i)
                    if (inRange3(expect(j),var_type(i), &
                                 NFT_ITYPE($1)) .and. &
                        in_internal_range(NFT_ITYPE($1), &
                                                expect(j))) then
                        val = ARITH3($1, value, j)
                        if (.not. equal(val, expect(j),  &
                                        var_type(i),  &
                                        NFT_ITYPE($1))) then
                            call errord('unexpected: ', val)
                        else
                            nok = nok + 1
                        end if
                    endif
5               continue
            else
                if (err .ne. NF90_ECHAR) &
                        call errore('wrong type: ', err)
            end if
1       continue
        err = nf90mpi_close(ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_close: ',  err)
        call print_nok(nok)
        end
])


dnl TEST_NFMPI_IGET_VARA(TYPE)
dnl
define([TEST_NFMPI_IGET_VARA],[dnl
        subroutine test_nf90mpi_iget_vara_$1()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer roll, index2indexes
        double precision hash4
        logical equal, inRange3, in_internal_range

        integer ncid
        integer d
        integer i
        integer j
        integer k
        integer err
        logical allInExtRange   
        logical allInIntRange   
        integer nels
        integer nslabs
        integer nok      
        integer(kind=MPI_OFFSET_KIND) start(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) edge(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) index(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) mid(MAX_RANK)
        logical canConvert     
        DATATYPE($1, value, (MAX_NELS))
        doubleprecision expect(MAX_NELS)
        doubleprecision val
        integer ud_shift
        integer err_w, reqid(1), st(1)      

        nok = 0

        err = nf90mpi_open(comm, testfile, NF90_NOWRITE, info, &
                         ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_open: ', err)
        do 1, i = 1, numVars
            canConvert = (var_type(i) .eq. NF90_CHAR) .eqv.  &
                         (NFT_ITYPE($1) .eq. NFT_TEXT)
            if (.not.(var_rank(i) .le. MAX_RANK)) stop 'assert'
            if (.not.(var_nels(i) .le. MAX_NELS)) stop 'assert'
            do 2, j = 1, var_rank(i)
                start(j) = 1
                edge(j) = 1
2           continue
            err = nf90mpi_iget_var(BAD_ID, i, value,reqid(1), start, edge)
            if (err .ne. NF90_EBADID) &
                call errore('bad ncid: ', err)
            err = nf90mpi_iget_var(ncid, BAD_VARID, value,reqid(1), start, edge)
            if (err .ne. NF90_ENOTVAR) &
                call errore('bad var id: ', err)
            do 3, j = 1, var_rank(i)
                start(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge)
                if (err .ne. NF90_EINVALCOORDS) &
                    call errore('bad index: ', err)

                start(j) = 1
                edge(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge)
                if (err == NF90_NOERR) then
                    err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    if (st(1) .ne. NF90_EINVALCOORDS) then
                        call errore('bad index/edge: ', st(1))
                    endif
                else
                    if (canConvert .and.  err .ne. NF90_EEDGE) &
                        call errore('bad index/edge: ', err)
                endif
                edge(j) = 1
3           continue

!           /* Check non-scalars for correct error returned even when */
!           /* there is nothing to get (edge(j).eq.0) */
            if (var_rank(i) .gt. 0) then
                do 10, j = 1, var_rank(i)
                    edge(j) = 0
10              continue
                err = nf90mpi_iget_var(BAD_ID, i, value,reqid(1), start, edge)
                if (err .ne. NF90_EBADID)  &
                    call errore('bad ncid: ', err)
                err = nf90mpi_iget_var(ncid, BAD_VARID, value,reqid(1), start, edge)
                if (err .ne. NF90_ENOTVAR)  &
                    call errore('bad var id: ', err)
                do 11, j = 1, var_rank(i)
                    if (var_dimid(j,i) .gt. 1) then     !/* skip record dim */
                        start(j) = var_shape(j,i) + 1
                        err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge)
                        if (err .ne. NF90_EINVALCOORDS) &
                            call errore( &
                            'Error:nf90mpi_iget_var: ', err)
                        start(j) = 1
                    endif
11              continue
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge)
                if (canConvert) then
                    if (err .ne. NF90_NOERR) then
                        call error(nf90mpi_strerror(err))
                    else
                        err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    endif
                else
                    if (err .ne. NF90_ECHAR) &
                        call errore('wrong type: ', err)
                endif
                do 12, j = 1, var_rank(i)
                    edge(j) = 1
12              continue
            endif

!           Choose a random point dividing each dim into 2 parts
!           get 2^rank (nslabs) slabs so defined
            nslabs = 1
            do 4, j = 1, var_rank(i)
                mid(j) = roll( var_shape(j,i) )
                nslabs = nslabs * 2
4           continue
!           bits of k determine whether to get lower or upper part of dim 
            do 5, k = 1, nslabs
                nels = 1
                do 6, j = 1, var_rank(i)
                    if (mod(ud_shift((k-1), -(j-1)), 2) .eq. 1) then
                        start(j) = 1
                        edge(j) = mid(j)
                    else
                        start(j) = 1 + mid(j)
                        edge(j) = var_shape(j,i) - mid(j)
                    end if
                    nels = nels * edge(j)
6               continue
                allInIntRange = .true.
                allInExtRange = .true.
                do 7, j = 1, nels
                    err = index2indexes(j, var_rank(i), edge, index)
                    if (err .ne. NF90_NOERR) &
                        call error('error in index2indexes 1')
                    do 8, d = 1, var_rank(i)
                        index(d) = index(d) + start(d) - 1
8                   continue
                    expect(j) = hash4(var_type(i), var_rank(i), index,  &
                                      NFT_ITYPE($1))
                    if (inRange3(expect(j),var_type(i),  &
                                 NFT_ITYPE($1))) then
                        allInIntRange =  &
                            allInIntRange .and. &
                            in_internal_range(NFT_ITYPE($1), expect(j))
                    else
                        allInExtRange = .false.
                    end if
7               continue
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge)
                if (err == NF90_NOERR) &
                    err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                if (canConvert) then
                    if (allInExtRange) then
                        if (allInIntRange) then
                            if (st(1) .ne. 0) &
                                call errore( &
                                    'nf90mpi_iget_var:',st(1))
                        else
                            if (st(1) .ne. NF90_ERANGE) &
                                call errore('Range error: ', st(1))
                        end if
                    else
                        if (st(1) .ne. 0 .and. st(1) .ne. NF90_ERANGE) &
                            call errore('OK or Range error: ', st(1))
                    end if
                    do 9, j = 1, nels
                        if (inRange3(expect(j),var_type(i), &
                                     NFT_ITYPE($1)) .and. &
                            in_internal_range(NFT_ITYPE($1), expect(j))) &
                                then
                            val = ARITH3($1, value, j)
                            if (.not.equal(val,expect(j), &
                                           var_type(i),NFT_ITYPE($1)))  &
                                    then
                                call error( &
                                    'value read not that expected')
                                if (verbose) then
                                    call error(' ')
                                    call errori('varid: ', i)
                                    call errorc('var_name: ', &
                                        var_name(i))
                                    call errori('element number: %d ',  &
                                                j)
                                    call errord('expect: ', expect(j))
                                    call errord('got: ', val)
                                end if
                            else
                                nok = nok + 1
                            end if
                        end if
9                   continue
                else
                    if (nels .gt. 0  .and. err .ne. NF90_ECHAR) &
                        call errore('wrong type: ', err)
                end if
5           continue
1       continue
        err = nf90mpi_close(ncid)
        if (err .ne. NF90_NOERR) &
            call errorc('nf90mpi_close: ', nf90mpi_strerror(err))
        call print_nok(nok)
        end
])dnl


dnl TEST_NFMPI_IGET_VARS(TYPE)
dnl
define([TEST_NFMPI_IGET_VARS],dnl
[dnl
        subroutine test_nf90mpi_iget_vars_$1()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer roll, index2indexes
        double precision hash4
        logical equal, inRange3, in_internal_range

        integer ncid
        integer d
        integer i
        integer j
        integer k
        integer m
        integer err
        logical allInExtRange   
        logical allInIntRange   
        integer nels
        integer nslabs
        integer nstarts         
        integer nok             
        integer(kind=MPI_OFFSET_KIND) start(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) edge(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) index(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) index2(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) mid(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) count(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) sstride(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) stride(MAX_RANK)
        logical canConvert     
        DATATYPE($1, value, (MAX_NELS))
        doubleprecision expect(MAX_NELS)
        doubleprecision val
        integer ud_shift
        integer err_w, reqid(1), st(1)      

        nok = 0

        err = nf90mpi_open(comm, testfile, NF90_NOWRITE, info, &
                         ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_open: ', err)
        do 1, i = 1, numVars
            canConvert = (var_type(i) .eq. NF90_CHAR) .eqv.  &
                         (NFT_ITYPE($1) .eq. NFT_TEXT)
            if (.not.(var_rank(i) .le. MAX_RANK)) stop 'assert'
            if (.not.(var_nels(i) .le. MAX_NELS)) stop 'assert'
            do 2, j = 1, var_rank(i)
                start(j) = 1
                edge(j) = 1
                stride(j) = 1
2           continue
            err = nf90mpi_iget_var(BAD_ID, i, value,reqid(1), start, edge, stride)
            if (err .ne. NF90_EBADID) &
                call errore('bad ncid: ', err)
            err = nf90mpi_iget_var(ncid, BAD_VARID, value,reqid(1), start, edge, stride)
            if (err .ne. NF90_ENOTVAR) &
                call errore('bad var id: ', err)
            do 3, j = 1, var_rank(i)
                start(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge, stride)
                if (err .ne. NF90_EINVALCOORDS) &
                    call errore('bad index: ', err)

                start(j) = 1
                edge(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge, stride)
                if (err == NF90_NOERR) then
                    err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    if (st(1) .ne. NF90_EINVALCOORDS) &
                        call errore('bad index: ', st(1))
                else
                    if (canConvert .and. err .ne. NF90_EEDGE) &
                        call errore('bad edge: ', err)
                endif
                edge(j) = 1
                stride(j) = 0
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge, stride)
                if (err .ne. NF90_ESTRIDE) &
                    call errore('bad stride: ', err)
                stride(j) = 1
3           continue
!               Choose a random point dividing each dim into 2 parts
!               get 2^rank (nslabs) slabs so defined
            nslabs = 1
            do 4, j = 1, var_rank(i)
                mid(j) = roll( var_shape(j,i) )
                nslabs = nslabs * 2
4           continue
!           bits of k determine whether to get lower or upper part of dim
!           choose random stride from 1 to edge
            do 5, k = 1, nslabs
                nstarts = 1
                do 6, j = 1, var_rank(i)
                    if (mod(ud_shift(k-1, j-1), 2) .eq. 1) then
                        start(j) = 1
                        edge(j) = mid(j)
                    else
                        start(j) = 1 + mid(j)
                        edge(j) = var_shape(j,i) - mid(j)
                    end if
                    if (edge(j) .gt. 0) then
                        sstride(j) = 1 + roll(edge(j))
                    else
                        sstride(j) = 1
                    end if
                    nstarts = nstarts * stride(j)
6               continue
                do 7, m = 1, nstarts
                    err = index2indexes(m, var_rank(i), sstride,  &
                                        index)
                    if (err .ne. NF90_NOERR) &
                        call error('error in index2indexes')
                    nels = 1
                    do 8, j = 1, var_rank(i)
                        count(j) = 1 + (edge(j) - index(j)) /  &
                                        stride(j)
                        nels = nels * count(j)
                        index(j) = index(j) + start(j) - 1
8                   continue
!                   Random choice of forward or backward 
!    /* TODO
!                   if ( roll(2) ) then
!                       for (j = 0 j < var_rank(i) j++) {
!                           index(j) += (count(j) - 1) * stride(j)
!                           stride(j) = -stride(j)
!                       }
!                   end if
!    */
                    allInIntRange = .true.
                    allInExtRange = .true.
                    do 9, j = 1, nels
                        err = index2indexes(j, var_rank(i), count,  &
                                            index2)
                        if (err .ne. NF90_NOERR) &
                            call error('error in index2indexes() 1')
                        do 10, d = 1, var_rank(i)
                            index2(d) = index(d) + (index2(d)-1) *  &
                                        stride(d)
10                      continue
                        expect(j) = hash4(var_type(i), var_rank(i),  &
                                          index2, NFT_ITYPE($1))
                        if (inRange3(expect(j),var_type(i), &
                                     NFT_ITYPE($1))) then
                            allInIntRange =  &
                                allInIntRange .and. &
                                in_internal_range(NFT_ITYPE($1),  &
                                                  expect(j))
                        else
                            allInExtRange = .false.
                        end if
9                   continue
                    err = nf90mpi_iget_var(ncid, i,value,reqid(1), index, count,stride)
                    if (err == NF90_NOERR) &
                        err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    if (canConvert) then
                        if (allInExtRange) then
                            if (allInIntRange) then
                                if (st(1) .ne. 0) &
                                    call error(nf90mpi_strerror(st(1)))
                            else
                                if (st(1) .ne. NF90_ERANGE) &
                                    call errore('Range error: ', st(1))
                            end if
                        else
                            if (st(1) .ne. 0 .and. st(1) .ne. NF90_ERANGE) &
                                call errore('OK or Range error: ',st(1))
                        end if
                        do 11, j = 1, nels
                            if (inRange3(expect(j),var_type(i), &
                                NFT_ITYPE($1)) .and. &
                                in_internal_range(NFT_ITYPE($1),  &
                                                  expect(j))) then
                                val = ARITH3($1, value, j)
                                if (.not.equal(val, expect(j), &
                                    var_type(i), NFT_ITYPE($1))) then
                                    call error( &
                                        'value read not that expected')
                                    if (verbose) then
                                        call error(' ')
                                        call errori('varid: ', i)
                                        call errorc('var_name: ',  &
                                                    var_name(i))
                                        call errori('element number: ', &
                                                    j)
                                        call errord('expect: ',  &
                                                    expect(j))
                                        call errord('got: ', val)
                                    end if
                                else
                                    nok = nok + 1
                                end if
                            end if
11                      continue
                    else
                        if (nels .gt. 0 .and. err .ne. NF90_ECHAR) &
                            call errore('wrong type: ', err)
                    end if
7               continue
5           continue

1       continue
        err = nf90mpi_close(ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_close: ', err)
        call print_nok(nok)
        end
])dnl


dnl TEST_NFMPI_IGET_VARM(TYPE)
dnl
define([TEST_NFMPI_IGET_VARM],dnl
[dnl
        subroutine test_nf90mpi_iget_varm_$1()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer roll, index2indexes
        double precision hash4
        logical equal, inRange3, in_internal_range

        integer ncid
        integer d
        integer i
        integer j
        integer k
        integer m
        integer err
        logical allInExtRange   
        logical allInIntRange   
        integer nels
        integer nslabs
        integer nstarts         
        integer nok             
        integer(kind=MPI_OFFSET_KIND) start(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) edge(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) index(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) index2(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) mid(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) count(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) sstride(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) stride(MAX_RANK)
        integer(kind=MPI_OFFSET_KIND) imap(MAX_RANK)
        logical canConvert     
        DATATYPE($1, value, (MAX_NELS))
        doubleprecision expect(MAX_NELS)
        doubleprecision val
        integer ud_shift
        integer err_w, reqid(1), st(1)      

        nok = 0

        err = nf90mpi_open(comm, testfile, NF90_NOWRITE, info, &
                         ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_open: ', err)
        do 1, i = 1, numVars
            canConvert = (var_type(i) .eq. NF90_CHAR) .eqv.  &
                         (NFT_ITYPE($1) .eq. NFT_TEXT)
            if (.not.(var_rank(i) .le. MAX_RANK)) stop 'assertion'
            if (.not.(var_nels(i) .le. MAX_NELS)) stop 'assertion'
            do 2, j = 1, var_rank(i)
                start(j) = 1
                edge(j) = 1
                stride(j) = 1
                imap(j) = 1
2           continue
            err = nf90mpi_iget_var(BAD_ID, i, value,reqid(1), start, edge, stride, imap)
            if (err .ne. NF90_EBADID) &
                call errore('bad ncid: ', err)
            err = nf90mpi_iget_var(ncid, BAD_VARID, value,reqid(1), start, edge, stride, imap)
            if (err .ne. NF90_ENOTVAR) &
                call errore('bad var id: ', err)
            do 3, j = 1, var_rank(i)
                start(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge, stride, imap)
                if (err .ne. NF90_EINVALCOORDS) &
                    call errore('bad index: ', err)

                start(j) = 1
                edge(j) = var_shape(j,i) + 1
                err = nf90mpi_iget_var(ncid, i, value,reqid(1), start, edge, stride, imap)
                if (err == NF90_NOERR) then
                    err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    if (st(1) .ne. NF90_EINVALCOORDS) &
                        call errore('bad index: ', st(1))
                else
                    if (canConvert .and. err .ne. NF90_EEDGE) &
                        call errore('bad edge: ', err)
                endif

                edge(j) = 1
                stride(j) = 0
                err = nf90mpi_iget_var(ncid, i, value, reqid(1), start, edge, stride, imap)
                if (err .ne. NF90_ESTRIDE) &
                    call errore('bad stride: ', err)
                stride(j) = 1
3           continue
!           Choose a random point dividing each dim into 2 parts 
!           get 2^rank (nslabs) slabs so defined 
            nslabs = 1
            do 4, j = 1, var_rank(i)
                mid(j) = roll( var_shape(j,i) )
                nslabs = nslabs * 2
4           continue
!           /* bits of k determine whether to get lower or upper part 
!            * of dim
!            * choose random stride from 1 to edge */
            do 5, k = 1, nslabs
                nstarts = 1
                do 6, j = 1, var_rank(i)
                    if (mod(ud_shift((k-1), -(j-1)), 2) .ne. 0) then
                        start(j) = 1
                        edge(j) = mid(j)
                    else
                        start(j) = 1 + mid(j)
                        edge(j) = var_shape(j,i) - mid(j)
                    end if
                    if (edge(j) .gt. 0) then
                        stride(j) = 1+roll(edge(j))
                    else
                        stride(j) = 1
                    end if
                    sstride(j) = stride(j)
                    nstarts = nstarts * stride(j)
6               continue
                do 7, m = 1, nstarts
                    err = index2indexes(m, var_rank(i), sstride, index)
                    if (err .ne. NF90_NOERR) &
                        call error('error in index2indexes')
                    nels = 1
                    do 8, j = 1, var_rank(i)
                        count(j) = 1 + (edge(j) - index(j)) /  &
                                        stride(j)
                        nels = nels * count(j)
                        index(j) = index(j) + start(j) - 1
8                   continue
!                   Random choice of forward or backward 
!    /* TODO
!                   if ( roll(2) ) then
!                       for (j = 0 j < var_rank(i) j++) {
!                           index(j) += (count(j) - 1) * stride(j)
!                           stride(j) = -stride(j)
!                       }
!                   end if
!     */
                    if (var_rank(i) .gt. 0) then
                        imap(1) = 1
                        do 9, j = 2, var_rank(i)
                            imap(j) = imap(j-1) * count(j-1)
9                       continue
                    end if
                    allInIntRange = .true.
                    allInExtRange = .true.
                    do 10, j = 1, nels
                        err = index2indexes(j, var_rank(i), count,  &
                                            index2)
                        if (err .ne. NF90_NOERR) &
                            call error('error in index2indexes 1')
                        do 11, d = 1, var_rank(i)
                            index2(d) = index(d) + (index2(d)-1) *  &
                                        stride(d)
11                      continue
                        expect(j) = hash4(var_type(i), var_rank(i),  &
                                          index2, NFT_ITYPE($1))
                        if (inRange3(expect(j),var_type(i), &
                                     NFT_ITYPE($1))) then
                            allInIntRange =  &
                                allInIntRange .and. &
                                in_internal_range(NFT_ITYPE($1), &
                                                  expect(j))
                        else
                            allInExtRange = .false.
                        end if
10                  continue
                    err = nf90mpi_iget_var(ncid,i, value, reqid(1), index, count, stride, imap)
                    if (err == NF90_NOERR) &
                        err_w = nf90mpi_wait_all(ncid,1,reqid,st)
                    if (canConvert) then
                        if (allInExtRange) then
                            if (allInIntRange) then
                                if (st(1) .ne. 0) &
                                    call error(nf90mpi_strerror(st(1)))
                            else
                                if (st(1) .ne. NF90_ERANGE) &
                                    call errore('Range error: ', st(1))
                            end if
                        else
                            if (st(1) .ne. 0 .and. st(1) .ne. NF90_ERANGE) &
                                call errore('OK or Range error: ',st(1))
                        end if
                        do 12, j = 1, nels
                            if (inRange3(expect(j),var_type(i), &
                                         NFT_ITYPE($1)) .and. &
                                in_internal_range(NFT_ITYPE($1), &
                                                  expect(j))) then
                                val = ARITH3($1, value, j)
                                if (.not.equal(val, expect(j), &
                                               var_type(i),  &
                                               NFT_ITYPE($1))) then
                                    call error( &
                                        'value read not that expected')
                                    if (verbose) then
                                        call error(' ')
                                        call errori('varid: ', i)
                                        call errorc('var_name: ',  &
                                                var_name(i))
                                        call errori('element number: ', &
                                                    j)
                                        call errord('expect: ',  &
                                                    expect(j))
                                        call errord('got: ', val)
                                    end if
                                else
                                    nok = nok + 1
                                end if
                            end if
12                      continue
                    else
                        if (nels .gt. 0 .and. err .ne. NF90_ECHAR) &
                            call errore('wrong type: ', err)
                    end if
7               continue
5           continue
1       continue
        err = nf90mpi_close(ncid)
        if (err .ne. NF90_NOERR) &
            call errore('nf90mpi_close: ',  err)
        call print_nok(nok)
        end
])dnl

divert(0)dnl
dnl If you see this line, you can ignore the next one.
! Do not edit this file. It is produced from the corresponding .m4 source */

!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!

TEST_NFMPI_IGET_VAR1(text)
#ifdef NF_INT1_T
TEST_NFMPI_IGET_VAR1(int1)
#endif
#ifdef NF_INT2_T
TEST_NFMPI_IGET_VAR1(int2)
#endif
TEST_NFMPI_IGET_VAR1(int)
TEST_NFMPI_IGET_VAR1(int8)
TEST_NFMPI_IGET_VAR1(real)
TEST_NFMPI_IGET_VAR1(double)

TEST_NFMPI_IGET_VAR(text)
#ifdef NF_INT1_T
TEST_NFMPI_IGET_VAR(int1)
#endif
#ifdef NF_INT2_T
TEST_NFMPI_IGET_VAR(int2)
#endif
TEST_NFMPI_IGET_VAR(int)
TEST_NFMPI_IGET_VAR(int8)
TEST_NFMPI_IGET_VAR(real)
TEST_NFMPI_IGET_VAR(double)

TEST_NFMPI_IGET_VARA(text)
#ifdef NF_INT1_T
TEST_NFMPI_IGET_VARA(int1)
#endif
#ifdef NF_INT2_T
TEST_NFMPI_IGET_VARA(int2)
#endif
TEST_NFMPI_IGET_VARA(int)
TEST_NFMPI_IGET_VARA(int8)
TEST_NFMPI_IGET_VARA(real)
TEST_NFMPI_IGET_VARA(double)

TEST_NFMPI_IGET_VARS(text)
#ifdef NF_INT1_T
TEST_NFMPI_IGET_VARS(int1)
#endif
#ifdef NF_INT2_T
TEST_NFMPI_IGET_VARS(int2)
#endif
TEST_NFMPI_IGET_VARS(int)
TEST_NFMPI_IGET_VARS(int8)
TEST_NFMPI_IGET_VARS(real)
TEST_NFMPI_IGET_VARS(double)

TEST_NFMPI_IGET_VARM(text)
#ifdef NF_INT1_T
TEST_NFMPI_IGET_VARM(int1)
#endif
#ifdef NF_INT2_T
TEST_NFMPI_IGET_VARM(int2)
#endif
TEST_NFMPI_IGET_VARM(int)
TEST_NFMPI_IGET_VARM(int8)
TEST_NFMPI_IGET_VARM(real)
TEST_NFMPI_IGET_VARM(double)

