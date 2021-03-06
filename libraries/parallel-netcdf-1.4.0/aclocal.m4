dnl $Id: aclocal.m4 1468 2013-10-26 16:53:18Z wkliao $
dnl UD macros for netcdf configure


dnl Convert a string to all uppercase.
dnl
define([uppercase],
[translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl
dnl Check for an m4(1) preprocessor utility.
dnl
AC_DEFUN(UD_PROG_M4,
[
    AC_CHECKING(for m4 preprocessor)
    case "${M4-unset}" in
	unset) AC_CHECK_PROGS(M4, m4 gm4, m4) ;;
	*) AC_CHECK_PROGS(M4, $M4 m4 gm4, m4) ;;
    esac
    AC_MSG_CHECKING(m4 flags)
    case "${M4FLAGS-unset}" in
	unset) M4FLAGS=-B10000 ;;
    esac
    AC_MSG_RESULT($M4FLAGS)
    AC_SUBST(M4FLAGS)
])

dnl
dnl Check for an ar(1) utility.
dnl
AC_DEFUN(UD_PROG_AR,
[
    AC_CHECKING(for ar utility)
    case "${AR-unset}" in
	unset) AC_CHECK_PROGS(AR, ar, ar) ;;
	*) AC_CHECK_PROGS(AR, $AR ar, ar) ;;
    esac
    AC_MSG_CHECKING(ar flags)
    case "${ARFLAGS-unset}" in
	unset) ARFLAGS=cru ;;
    esac
    AC_MSG_RESULT($ARFLAGS)
    AC_SUBST(ARFLAGS)
])

dnl
dnl Check for an nm(1) utility.
dnl
AC_DEFUN(UD_PROG_NM,
[
    AC_CHECKING(for nm utility)
    case "${NM-unset}" in
	unset) AC_CHECK_PROGS(NM, nm, nm) ;;
	*) AC_CHECK_PROGS(NM, $NM nm, nm) ;;
    esac
    AC_MSG_CHECKING(nm flags)
    case "${NMFLAGS-unset}" in
	unset) NMFLAGS= ;;
    esac
    AC_MSG_RESULT($NMFLAGS)
    AC_SUBST(NMFLAGS)
])

dnl
dnl Set the top-level source-directory.
dnl
AC_DEFUN(UD_SRCDIR,
[
    AC_MSG_CHECKING(for top-level source-directory)
    SRCDIR=`(cd $srcdir && pwd)`
    AC_MSG_RESULT($SRCDIR)
    AC_SUBST(SRCDIR)
])

dnl
dnl Check for a Standard C compiler.  Prefer a native one over the
dnl GNU one to reduce the chance that the environment variable LIBS
dnl will have to be set to reference the GNU C runtime library.
dnl
AC_DEFUN(UD_PROG_CC,
[
    # Because we must have a C compiler, we treat an unset CC
    # the same as an empty CC.
    case "${CC}" in
	'')
	    case `uname` in
		ULTRIX)
		    # The native ULTRIX C compiler isn't standard.
		    ccs='gcc cc'
		    ;;
		*)
		    # xlc is before c89 because AIX's sizeof(long long)
		    # differs between the two.
		    #
		    ccs='xlc c89 acc cc gcc'
		    ;;
	    esac
	    for cc in $ccs; do
		AC_CHECK_PROG(CC, $cc, $cc)
		case "$CC" in
		    '') ;;
		    *)  break
			;;
		esac
	    done
	    case "${CC}" in
		'')
		    AC_MSG_ERROR("Could not find C compiler")
		    ;;
	    esac
	    ;;
    esac
    #
    # On some systems, a discovered compiler nevertheless won't
    # work (due to licensing, for example); thus, we check the
    # compiler with a test program.
    # 
    AC_MSG_CHECKING(C compiler \"$CC\")
    AC_TRY_COMPILE(, ,
	AC_MSG_RESULT(works),
	AC_MSG_RESULT(failed to compile test program))
    AC_SUBST(CC)
    case "$CC" in
	*gcc*)
	    GCC=yes		# Expected by autoconf(1) macros
	    ;;
    esac
    case `uname -sr` in
	'HP-UX A.09'*)
	    AC_DEFINE(_HPUX_SOURCE)
	    ;;
    esac
])

dnl
dnl Check for a C++ compiler.  Prefer a native one over the
dnl GNU one to reduce the chance that the environment variable LIBS
dnl will have to be set to reference the GNU C runtime library.
dnl
AC_DEFUN(UD_PROG_CXX,
[
    case "${CXX-unset}" in
	unset)
	    case `uname` in
		AIX)
		    preferred_cxx='xlC'
		    ;;
	    esac
	    possible_cxxs="${preferred_cxx} CC cxx c++ g++ gcc"
	    ;;
	'') AC_MSG_WARN("Empty CXX variable")
	    possible_cxxs=
	    ;;
	*)  possible_cxxs=$CXX
	    ;;
    esac
    case "${possible_cxxs}" in
	'') CXX=
	    ;;
	*)  AC_LANG_SAVE()
	    AC_LANG_CPLUSPLUS()
	    for cxx in $possible_cxxs; do
		AC_CHECK_PROG(CXX, $cxx, $cxx)
		case "$CXX" in
		    '') ;;
		    *)  # On some systems, a discovered compiler nevertheless
			# won't work (because it's a script to a non-existant
			# executable, for example); thus, we check the
			# compiler with a test program.  We also test
			# for <iostream.h> and the standard C++ library
			# because we need these to work.
			# 
			AC_MSG_CHECKING(C++ compiler \"$CXX\")
			AC_TRY_RUN(
 			    [
				#include <iostream.h>
				int main() {
				    cout << "";
				    return 0;
				}
			    ],
			    [
				AC_MSG_RESULT(works)
				break
			    ],
			    [
				AC_MSG_WARN($CXX failed on test program)
				CXX=
				unset ac_cv_prog_CXX
			    ])
			;;
		esac
	    done
	    AC_LANG_RESTORE()
	    case "${CXX}" in
		'') AC_MSG_WARN("Could not find working C++ compiler")
		    AC_MSG_WARN(Setting CXX to the empty string)
		    ;;
	    esac
	    ;;
    esac
    case "${CXX}" in
	'') AC_MSG_WARN(The C++ interface will not be built)
	    ;;
    esac
    AC_SUBST(CXX)
    case `uname` in
	'HP-UX A.09'*)
	    AC_DEFINE(_HPUX_SOURCE)
	    ;;
    esac
])


dnl
dnl like AC_LONG_DOUBLE, except checks for 'long long'
dnl
AC_DEFUN(UD_C_LONG_LONG,
[AC_MSG_CHECKING(for long long)
AC_CACHE_VAL(ac_cv_c_long_long,
[if test "$GCC" = yes; then
  ac_cv_c_long_long=yes
else
AC_TRY_RUN([int main() {
long long foo = 0;
return(sizeof(long long) < sizeof(long)); }],
ac_cv_c_long_long=yes, ac_cv_c_long_long=no, :)
fi])dnl
AC_MSG_RESULT($ac_cv_c_long_long)
if test $ac_cv_c_long_long = yes; then
  AC_DEFINE(HAVE_LONG_LONG)
fi
])

dnl 
dnl UD_CHECK_IEEE
dnl If the 'double' is not an IEEE double
dnl or the 'float' is not and IEEE single,
dnl define NO_IEEE_FLOAT
dnl
AC_DEFUN(UD_CHECK_IEEE,
[
AC_MSG_CHECKING(for IEEE floating point format)
AC_TRY_RUN([#ifndef NO_FLOAT_H
#include <float.h>
#endif

#define EXIT_NOTIEEE	1
#define EXIT_MAYBEIEEE	0

int
main()
{
#if	defined(FLT_RADIX)	&& FLT_RADIX != 2
		return EXIT_NOTIEEE;
#elif	defined(DBL_MAX_EXP)	&& DBL_MAX_EXP != 1024
		return EXIT_NOTIEEE;
#elif	defined(DBL_MANT_DIG)	&& DBL_MANT_DIG != 53
		return EXIT_NOTIEEE;
#elif 	defined(FLT_MAX_EXP)	&& !(FLT_MAX_EXP == 1024 || FLT_MAX_EXP == 128)
		return EXIT_NOTIEEE;
#elif	defined(FLT_MANT_DIG)	&& !(FLT_MANT_DIG == 53 || FLT_MANT_DIG == 24)
		return EXIT_NOTIEEE;
#else
	/* (assuming eight bit char) */
	if(sizeof(double) != 8)
		return EXIT_NOTIEEE;
	if(!(sizeof(float) == 4 || sizeof(float) == 8))
		return EXIT_NOTIEEE;

	return EXIT_MAYBEIEEE;
#endif
}],ac_cv_c_ieeefloat=yes, ac_cv_c_ieeefloat=no, :)
AC_MSG_RESULT($ac_cv_c_ieeefloat)
if test x$ac_cv_c_ieeefloat = xno; then
  AC_DEFINE(NO_IEEE_FLOAT)
fi
])

dnl Check for utility for generating makefile dependencies.
dnl Should only be used at the UPC.
dnl
AC_DEFUN(UD_PROG_CC_MAKEDEPEND,
[
    AC_MSG_CHECKING(how to make dependencies)
    case `uname -s` in
	IRIX*|OSF1)
	    CC_MAKEDEPEND='cc -M'
	    ;;
	SunOS)
	    case `uname -r` in
		4*)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		5*|*)
		    CC_MAKEDEPEND='cc -xM'
		    ;;
	    esac
	    ;;
	ULTRIX)
	    case `uname -m` in
		RISC)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		VAX)	# Can't handle prototypes in netcdf.h
		    ;;
	    esac
	    ;;
	AIX)	# Writes to .u files rather than standard out
	    ;;
	HP-UX)	# Writes escaped newlines to standard error
	    ;;
    esac
    case "${CC_MAKEDEPEND}" in
	'')
	    CC_MAKEDEPEND=false
	    ;;
    esac
    AC_MSG_RESULT($CC_MAKEDEPEND)
    AC_SUBST(CC_MAKEDEPEND)
])


dnl Check for Fortran-90 compiler.
dnl
AC_DEFUN(UD_PROG_F90,
[
    case "${F90+set}" in
	set)
	    AC_MSG_CHECKING(user-defined Fortran-90 compiler \"$F90\")
            AC_LANG_PUSH([Fortran])
            AC_COMPILE_IFELSE(
               [AC_LANG_SOURCE([
		   subroutine foo(bar)
		   integer, intent(in) :: bar
		   end subroutine foo
               ])],
               [AC_MSG_RESULT(works)],
               [AC_MSG_RESULT(failed to compile test program)
		unset F90
               ]
            )
            AC_LANG_POP([Fortran])
	    ;;
	*)
	    case "${FC+set}" in
		set)
		    F90=$FC
		    F90FLAGS="${F90FLAGS-${FFLAGS--O}}"
		    F90LIBS="${F90LIBS-${FLIBS}}"
		    cat <<EOF >conftest.f90
			program foo
			call bar(1)
			end program foo
			subroutine bar(bof)
			integer, intent(in) :: bof
			end subroutine bar
EOF
		    AC_MSG_CHECKING(\"$F90\" as Fortran-90 compiler)
		    doit='$F90 -o conftest ${F90FLAGS} conftest.f90 ${F90LIBS}'
		    if AC_TRY_EVAL(doit); then
			doit=./conftest
			if AC_TRY_EVAL(doit); then
			    AC_MSG_RESULT(works)
			else
			    AC_MSG_RESULT(failed to build executable program)
			    unset F90
			fi
		    else
			AC_MSG_RESULT(failed to build test program)
			unset F90
		    fi
		    ${RM} -f conftest*
		    ;;
	    esac
	    case "${F90-unset}" in
		unset)
		    cat <<EOF >conftest.f90
			program foo
			call bar(1)
			end program foo
			subroutine bar(bof)
			integer, intent(in) :: bof
			end subroutine bar
EOF
		    for f90 in xlf90 pgf90 f90; do
			AC_CHECK_PROG(F90, $f90, $f90)
			case "${F90}" in
			    '')
				;;
			    *)
				AC_MSG_CHECKING(Fortran-90 compiler \"$F90\")
				doit='$F90 -o conftest ${F90FLAGS} conftest.f90 ${F90LIBS}'
				if AC_TRY_EVAL(doit); then
				    doit=./conftest
				    if AC_TRY_EVAL(doit); then
					AC_MSG_RESULT(works)
					break
				    else
					AC_MSG_RESULT(
					    failed to build executable program)
					unset F90
					unset ac_cv_prog_F90
				    fi
				else
				    AC_MSG_RESULT(failed to build test program)
				    unset F90
				    unset ac_cv_prog_F90
				fi
				;;
			esac
		    done
		    ${RM} -f conftest*
		    case "${F90}" in
			'') AC_MSG_WARN(
			    "Could not find working Fortran-90 compiler")
			    ;;
		    esac
		    ;;
	    esac
	    ;;
    esac
    case "${F90}" in
	'')
	    AC_MSG_WARN("The Fortran-90 interface will not be built")
	    ;;
	*f90*)
	    case `uname -s` in
		IRIX*)
		    NETCDF_MOD=NETCDF.mod
		    ;;
		*)
		    NETCDF_MOD=netcdf.mod
		    ;;
	    esac
	    AC_SUBST(NETCDF_MOD)
	    ;;
    esac
    AC_SUBST(F90)
    AC_SUBST(F90FLAGS)
    AC_SUBST(F90LIBS)
])

dnl Check for Fortran-77 compiler.
dnl
AC_DEFUN(UD_PROG_FC,
[
    AC_BEFORE([UD_FORTRAN_TYPES])
    case "${FC+set}" in
	set)
	    case "$FC" in
		'')
		    AC_MSG_WARN(Fortran-77 compiler is explicitly set to null)
		    ;;
		*)
		    AC_MSG_CHECKING(user-defined Fortran-77 compiler \"$FC\")
                    AC_LANG_PUSH([Fortran 77])
                    AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                       [AC_MSG_RESULT(works)],
                       [AC_MSG_RESULT(failed to compile test program)
                        FC=
                       ]
                    )
                    AC_LANG_POP([Fortran 77])
		    ;;
	    esac
	    ;;
	*)
            dnl FC is not set, let's try F90 as FC if F90 is set
	    case "${F90+set}" in
		set)
		    FC=$F90
		    FFLAGS="${FFLAGS-${F90FLAGS--O}}"
		    FLIBS="${FLIBS-${F90LIBS-}}"
		    AC_MSG_CHECKING(\"$FC\" as Fortran-77 compiler)
                    AC_LANG_PUSH([Fortran])
                    AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                       [AC_MSG_RESULT(works)],
                       [AC_MSG_RESULT(failed to compile test program)
			unset FC
                       ]
                    )
                    AC_LANG_POP([Fortran])
		    ;;
		*)
		    ;;
	    esac
	    case "${FC-unset}" in
		unset)
		    case `uname -sr` in
			AIX*)
			    # xlf90(1) thinks fortran/ftest.F has bad syntax.
			    forts="xlf f77 gfortran"
			    ;;
			BSD/OS*|FreeBSD*)
			    forts="f77 fort77 g77 gfortran"
			    ;;
			HP-UX*)
			    # f77(1) doesn't have the -L option.
			    forts=fort77
			    FLIBS=-lU77
			    ;;
			IRIX*)
			    # f90(1) can't link with c89(1)-compiled objects
			    forts=f77
			    ;;
			IRIX64*)
			    forts='f77 g77 gfortran fort77 '
			    ;;
			Linux*)
			    forts="pgf90 f77 fort77 g77 gfortran"
			    ;;
			OSF1*)
			    # The use of f90(1) results in the following for
			    # an unknown reason (`make' works in the fortran/
			    # directory):
			    # f90 -c -I../libsrc ftest.F 
			    # Last chance handler: pc = 0xa971b8, 
			    # sp = 0x3fece0, ra = 0xa971b8
			    # Last chance handler: internal exception: unwinding
			    forts="f77 gfortran"
			    ;;
			'SunOS 4'*)
			    forts='f77 g77 gfortran fort77'
			    ;;
			'SunOS 5'*)
			    # SunOS's f90(1) has problems passing a C `char'
			    # as a Fortran `integer*1' => use f77(1)
			    forts="pgf90 f77"
			    ;;
			sn*|UNICOS*|unicos*)
			    forts="fort77 cf77 f77 g77 gfortran f90"
			    ;;
			*)
			    forts="xlf fort77 ghf77 f77 cf77 g77 gfortran xlf90 f90"
			    ;;
		    esac
		    for fc in $forts; do
			AC_CHECK_PROG(FC, $fc, $fc)
			case "${FC}" in
			    '')
				;;
			    *)
				#
				# On some systems, a discovered compiler
				# nevertheless won't work (due to licensing,
				# for example); thus, we check the compiler
				# with a test program.
				# 
                                AC_LANG_PUSH([Fortran])
                                AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                                   [AC_MSG_RESULT(works)],
                                   [AC_MSG_RESULT(failed to compile test program)
			            unset FC
				    unset ac_cv_prog_FC
                                   ]
                                )
                                AC_LANG_POP([Fortran])
				;;
			esac
		    done
		    ${RM} -f conftest.*
		    case "${FC}" in
			'') AC_MSG_WARN(
				"Could not find working Fortran-77 compiler")
			    ;;
		    esac
		    ;;
	    esac
	    ;;
    esac
    case "${FC}" in
	'') AC_MSG_WARN("The Fortran-77 interface will not be built")
	    ;;
    esac
    AC_SUBST(FC)
    AC_SUBST(FFLAGS)
    AC_SUBST(FLIBS)
    #
    # Set the make(1) macro for compiling a .F file.
    #
    case "${FPP-}" in
    '')
	AC_MSG_CHECKING(for Fortran .F compiler)
	AC_MSG_RESULT($COMPILE_F)
	case "${COMPILE_F-unset}" in
	unset)
	    case "${FC}" in
	    '')
		COMPILE_F=
		;;
	    *)
		cat >conftest.h <<\EOF
#define J 1
EOF
                AC_LANG_PUSH([Fortran])
                AC_FC_SRCEXT([F])
		AC_MSG_CHECKING(if Fortran-77 compiler handles *.F files)
                AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#include "conftest.h"
#define N 5
		  real r(J,N)
                   ])],
                   [AC_MSG_RESULT(yes)
                    COMPILE_F='$(COMPILE.f)'],
                   [AC_MSG_RESULT(no)
		    COMPILE_F=
                   ]
                )
                AC_FC_SRCEXT([f])
                AC_LANG_POP([Fortran])
		${RM} -f conftest.h
		;;
	    esac
	    ;;
	esac
	;;
    *)
	unset COMPILE_F
	;;
    esac
    case "${COMPILE_F-}" in
	'') UD_PROG_FPP;;
    esac
    AC_SUBST(COMPILE_F)
    FPPFLAGS=${FPPFLAGS-}
    AC_SUBST(FPPFLAGS)
])


dnl Check for Fortran preprocessor.
dnl
AC_DEFUN(UD_PROG_FPP,
[
    AC_MSG_CHECKING(for Fortran preprocessor)
    case "$FPP" in
    '')
	AC_REQUIRE([AC_PROG_CPP])
	FPP="$CPP"
	;;
    esac
    AC_MSG_RESULT($FPP)
    AC_SUBST(FPP)
])


dnl Check for a Fortran type equivalent to a netCDF type.
dnl
dnl UD_CHECK_FORTRAN_NCTYPE(forttype, possibs, nctype)
dnl
AC_DEFUN(UD_CHECK_FORTRAN_NCTYPE,
[
    AC_MSG_CHECKING(for Fortran-equivalent to netCDF \"$3\")
dnl     for type in $2; do
dnl         cat >conftest.f <<EOF
dnl                $type foo
dnl                end
dnl EOF
dnl         doit='$FC -c ${FFLAGS} conftest.f'
dnl         if AC_TRY_EVAL(doit); then
dnl             break
dnl         fi
dnl     done
dnl     ${RM} -f conftest.f conftest.o

    AC_LANG_PUSH([Fortran])
    for type in $2; do
        AC_COMPILE_IFELSE(
           [AC_LANG_SOURCE([
               $type foo
               end
           ])],
           [break]
        )
    done
    AC_LANG_POP([Fortran])
    AC_DEFINE_UNQUOTED($1, $type)
    AC_MSG_RESULT($type)
    $1=$type
])


dnl Check for a Fortran type equivalent to a C type.
dnl
dnl UD_CHECK_FORTRAN_CTYPE(v3forttype, v2forttype, ctype, min, max)
dnl
AC_DEFUN(UD_CHECK_FORTRAN_CTYPE,
[
    AC_MSG_CHECKING(for Fortran-equivalent to C \"$3\")
    AC_LANG_PUSH([Fortran])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
        subroutine sub(values, minval, maxval)
        implicit        none
        $2              values(5), minval, maxval
        minval = values(2)
        maxval = values(4)
        if (values(2) .ge. values(4)) then
            minval = values(4)
            maxval = values(2)
        endif
        end
       ])],
       [found_f2c=yes], [found_f2c=no]
    )
    AC_LANG_POP([Fortran])
    if test "x${found_f2c}" = xyes ; then
dnl     cat >conftest.f <<EOF
dnl         subroutine sub(values, minval, maxval)
dnl         implicit        none
dnl         $2              values(5), minval, maxval
dnl         minval = values(2)
dnl         maxval = values(4)
dnl         if (values(2) .ge. values(4)) then
dnl             minval = values(4)
dnl             maxval = values(2)
dnl         endif
dnl         end
dnl EOF
dnl     doit='$FC -c ${FFLAGS} conftest.f'
dnl     if AC_TRY_EVAL(doit); then
dnl         mv conftest.o conftestf.o
	cat >conftest.c <<EOF
#include <limits.h>
#include <float.h>
void main()
{
$3		values[[]] = {0, $4, 0, $5, 0};
$3		minval, maxval;
int	$FCALLSCSUB($3*, $3*, $3*);
$FCALLSCSUB(values, &minval, &maxval);
return(!(minval == $4 && maxval == $5));
}
EOF
	doit='$CC -o conftest ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} conftest.c conftestf.o ${LIBS}'
	if AC_TRY_EVAL(doit); then
	    doit=./conftest
	    if AC_TRY_EVAL(doit); then
		AC_MSG_RESULT($2)
		$1=$2
		AC_DEFINE_UNQUOTED($1,$2)
	    else
		AC_MSG_RESULT(no equivalent type)
		unset $1
	    fi
	else
	    AC_MSG_ERROR(Could not compile-and-link conftest.c and conftestf.o)
	fi
    else
	AC_MSG_ERROR(Could not compile conftest.f)
    fi
    ${RM} -f conftest*
    unset found_f2c
])


dnl Check for a Fortran data type.
dnl
dnl UD_CHECK_FORTRAN_TYPE(varname, ftypes)
dnl
AC_DEFUN(UD_CHECK_FORTRAN_TYPE,
[
    AC_LANG_PUSH([Fortran])
    for ftype in $2; do
	AC_MSG_CHECKING(for Fortran \"$ftype\")
        AC_COMPILE_IFELSE(
           [AC_LANG_SOURCE([
               subroutine sub(value)
               $ftype value
               end
           ])],
           [AC_MSG_RESULT(yes)
	    $1=$ftype
	    AC_DEFINE_UNQUOTED($1, $ftype)
            break],
           [AC_MSG_RESULT(no)]
        )
    done
    AC_LANG_POP([Fortran])
])


dnl Check for the name format of a Fortran-callable C routine.
dnl
dnl UD_CHECK_FCALLSCSUB
AC_DEFUN([UD_CHECK_FCALLSCSUB],
[
    dnl AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
	'') ;;
	*)
	    AC_REQUIRE([UD_PROG_NM])
	    AC_BEFORE([UD_CHECK_FORTRAN_CTYPE])
	    AC_BEFORE([UD_CHECK_CTYPE_FORTRAN])
	    AC_MSG_CHECKING(for C-equivalent to Fortran routine \"SUB\")
            AC_FC_FUNC(sub, [FCALLSCSUB])
            AC_MSG_RESULT($FCALLSCSUB)
	    ;;
    esac
])


dnl Check for a C type equivalent to a Fortran type.
dnl
dnl UD_CHECK_CTYPE_FORTRAN(ftype, ctypes, fmacro_root)
dnl
AC_DEFUN(UD_CHECK_CTYPE_FORTRAN,
[
    cat >conftestf.f <<EOF
           $1 values(4)
           integer status, sub
           data values /-1, -2, -3, -4/
           status = sub(values)
           call exit(status)
           end
EOF
    ac_cv_ctype_fortran=no
    AC_MSG_CHECKING(if Fortran \"$1\" is )
    for ctype in $2; do
	dnl AC_MSG_CHECKING(if Fortran \"$1\" is C \"$ctype\")
	cat >conftest.c <<EOF
	    int $FCALLSCSUB($ctype values[[4]])
	    {
		return(values[[1]] != -2 || values[[2]] != -3);
	    }
EOF
	doit='$CC -c ${CPPFLAGS} ${CFLAGS} conftest.c'
	if AC_TRY_EVAL(doit); then
	    doit='$FC ${FFLAGS} -c conftestf.f'
	    if AC_TRY_EVAL(doit); then
	        doit='$FC -o conftest ${FFLAGS} ${FLDFLAGS} conftestf.o conftest.o ${LIBS}'
	        if AC_TRY_EVAL(doit); then
		    doit=./conftest
		    if AC_TRY_EVAL(doit); then
		        dnl AC_MSG_RESULT(yes)
		        AC_MSG_RESULT(\"$ctype\" in C)
		        cname=`echo $ctype | tr ' abcdefghijklmnopqrstuvwxyz' \
			    _ABCDEFGHIJKLMNOPQRSTUVWXYZ`
		        AC_DEFINE_UNQUOTED(NF_$3[]_IS_C_$cname)
                        ac_cv_ctype_fortran=yes
		        break
		    fi
	        else
		    AC_MSG_ERROR(Could not link conftestf.o and conftest.o)
	        fi
	    else
		AC_MSG_ERROR(Could not compile conftestf.f)
	    fi
	else
	    AC_MSG_ERROR(Could not compile conftest.c)
	fi
    done
    ${RM} -f conftest*

    if test "$ac_cv_ctype_fortran" = no ; then
        AC_MSG_RESULT(no correspond data type in C)
    fi
    unset ac_cv_ctype_fortran
])


dnl Get information about Fortran data types.
dnl
AC_DEFUN([UD_FORTRAN_TYPES],
[
    dnl AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
    '')
	;;
    *)
	AC_REQUIRE([UD_CHECK_FCALLSCSUB])
	UD_CHECK_FORTRAN_TYPE(NF_INT1_T, integer*1 byte "integer(kind=1)")
	UD_CHECK_FORTRAN_TYPE(NF_INT2_T, integer*2 "integer(kind=2)")
	UD_CHECK_FORTRAN_TYPE(NF_INT8_T, integer*8 "integer(kind=8)")

	case "${NF_INT1_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "signed char" short int long, INT1)
		;;
	esac
	case "${NF_INT2_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT2_T, short int long, INT2)
		;;
	esac
	case "${NF_INT8_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT8_T, int long "long long", INT8)
		;;
	esac
	UD_CHECK_CTYPE_FORTRAN(integer, int long, INT)
	UD_CHECK_CTYPE_FORTRAN(real, float double, REAL)
	UD_CHECK_CTYPE_FORTRAN(doubleprecision, double float, DOUBLEPRECISION)

	UD_CHECK_FORTRAN_NCTYPE(NCBYTE_T, byte integer*1 integer, byte)

	UD_CHECK_FORTRAN_NCTYPE(NCSHORT_T, integer*2 integer, short)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_SHORT_T, $NCSHORT_T, short, SHRT_MIN, SHRT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCLONG_T, integer*4 integer, long)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_INT_T, integer, int, INT_MIN, INT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCFLOAT_T, real*4 real, float)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_FLOAT_T, $NCFLOAT_T, float, FLT_MIN, FLT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCDOUBLE_T, real*8 doubleprecision real, double)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_DOUBLE_T, $NCDOUBLE_T, double, DBL_MIN, DBL_MAX)
	;;
    esac
])


dnl Setup for making a manual-page database.
dnl
AC_DEFUN(UD_MAKEWHATIS,
[
    #
    # NB: We always want to define WHATIS to prevent the
    # $(MANDIR)/$(WHATIS) make(1) target from being just $(MANDIR)/ and
    # conflicting with the (directory creation) target with the same name.
    #
    WHATIS=whatis
    case `uname -sr` in
	BSD/OS*|FreeBSD*)
	    # Can't generate a user-database -- only /usr/share/man/whatis.db.
	    MAKEWHATIS_CMD=
	    ;;
	'IRIX64 6.5'|'IRIX 6.5')
	    MAKEWHATIS_CMD='/usr/lib/makewhatis -M $(MANDIR) $(MANDIR)/whatis'
	    ;;
	'IRIX 6'*)
	    # Can't generate a user-database.
	    MAKEWHATIS_CMD=
	    ;;
	HP-UX*)
	    # Can't generate a user-database -- only /usr/lib/whatis.
	    MAKEWHATIS_CMD=
	    ;;
	'Linux '*)
	    # /usr/sbin/makewhatis doesn't work
	    MAKEWHATIS_CMD=
	    ;;
	ULTRIX*)
	    # Can't generate a user-database -- only /usr/lib/whatis.
	    MAKEWHATIS_CMD=
	    ;;
	*)
	    if test -r /usr/man/windex; then
		WHATIS=windex
	    fi
	    AC_CHECK_PROGS(prog, catman makewhatis /usr/lib/makewhatis)
	    case "$prog" in
		*catman*)
		    MAKEWHATIS_CMD=$prog' -w -M $(MANDIR)'
		    ;;
		*makewhatis*)
		    MAKEWHATIS_CMD=$prog' $(MANDIR)'
		    ;;
	    esac
	    ;;
    esac
    AC_SUBST(WHATIS)
    AC_SUBST(MAKEWHATIS_CMD)
    AC_MSG_CHECKING(for manual-page index command)
    AC_MSG_RESULT($MAKEWHATIS_CMD)
])


dnl Check for the math library.
dnl
AC_DEFUN(UD_CHECK_LIB_MATH,
[
    AC_CHECKING(for math library)
    case "${MATHLIB}" in
	'')
	    AC_CHECK_LIB(c, tanh, MATHLIB=,
		    AC_CHECK_LIB(m, tanh, MATHLIB=-lm, MATHLIB=))
	    ;;
	*)
	    AC_MSG_RESULT($MATHLIB (user defined))
	    ;;
    esac
    AC_SUBST(MATHLIB)
])

dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_SRCEXT],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile preprocessed .$1 files],
                ac_cv_fc_pp_srcext_$1,
[ac_ext=$1
ac_fcflags_pp_srcext_save=$ac_fcflags_srcext
ac_fcflags_srcext=
ac_cv_fc_pp_srcext_$1=unknown
case $ac_ext in #(
  [[fF]]77) ac_try=f77-cpp-input;; #(
  *) ac_try=f95-cpp-input;;
esac
for ac_flag in none -ftpp -fpp -Tf "-fpp -Tf" -xpp=fpp -Mpreprocess "-e Z" \
               -cpp -xpp=cpp -qsuffix=cpp=$1 "-x $ac_try" +cpp -Cpp; do
  test "x$ac_flag" != xnone && ac_fcflags_srcext="$ac_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 0
#include <ac_nonexistent.h>
      choke me
#endif]])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 1
#include <ac_nonexistent.h>
      choke me
#endif]])],
       [],
       [ac_cv_fc_pp_srcext_$1=$ac_flag; break])])
done
rm -f conftest.$ac_objext conftest.$1
ac_fcflags_srcext=$ac_fcflags_pp_srcext_save
])
if test "x$ac_cv_fc_pp_srcext_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran could not compile preprocessed .$1 files])])
else
  ac_fc_srcext=$1
  if test "x$ac_cv_fc_pp_srcext_$1" = xnone; then
    ac_fcflags_srcext=""
    FCFLAGS_[]$1[]=""
  else
    ac_fcflags_srcext=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_[]$1[]=$ac_cv_fc_pp_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
AC_LANG_POP(Fortran)dnl
])# UD_FC_PP_SRCEXT

dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_DEFINE],
[AC_LANG_PUSH([Fortran])dnl
ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_ext_saved=$ac_ext
UD_FC_PP_SRCEXT([F])
ac_ext=F
AC_CACHE_CHECK([how to define symbols for preprocessed Fortran],
  [ac_cv_fc_pp_define],
[ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_cv_fc_pp_define=unknown
ac_fc_pp_define_FCFLAGS_save=$FCFLAGS
for ac_flag in -D -WF,-D -Wp,-D -Wc,-D
do
  FCFLAGS="$ac_fc_pp_define_FCFLAGS_save ${ac_flag}FOOBAR ${ac_flag}ZORK=42"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#ifndef FOOBAR
      choke me
#endif
#if ZORK != 42
      choke me
#endif]])],
    [ac_cv_fc_pp_define=$ac_flag])
  test x"$ac_cv_fc_pp_define" != xunknown && break
done
FCFLAGS=$ac_fc_pp_define_FCFLAGS_save
])
ac_fc_srcext=$ac_fc_pp_define_srcext_save
if test "x$ac_cv_fc_pp_define" = xunknown; then
  FC_DEFINE=
  m4_default([$2],
             [AC_MSG_ERROR([Fortran does not allow to define preprocessor symbols], 77)])
else
  FC_DEFINE=$ac_cv_fc_pp_define
  $1
fi
ac_ext=$ac_ext_saved
AC_SUBST([FC_DEFINE])dnl
AC_LANG_POP([Fortran])dnl
])

# AC_PROG_FC_MOD
# ---------------
dnl Note that Mac OSX file system is case-insensitive, so this function does
dnl not work precisely on Mac. Hence, we check whether the file system can
dnl find the file with name in lowercase. If not, we say the mod file is in
dnl uppercase.
AC_DEFUN([UD_PROG_FC_UPPERCASE_MOD],
[
AC_REQUIRE([UD_FC_MODULE_EXTENSION])
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([if Fortran 90 compiler capitalizes .mod filenames])
AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([
        module conftest
        end module conftest
    ])]
)
dnl ac_try='$F90 ${F90FLAGS} conftest.f90 ${F90LIBS}>&AS_MESSAGE_LOG_FD'
dnl AC_TRY_EVAL(ac_try)
if test -f conftest.${FC_MODEXT} ; then
   ac_cv_prog_f90_uppercase_mod=no
else
   ac_cv_prog_f90_uppercase_mod=yes
   ${RM} -f CONFTEST.${FC_MODEXT}
fi
AC_MSG_RESULT($ac_cv_prog_f90_uppercase_mod)
${RM} -f conftest*
AC_LANG_POP(Fortran)
])

dnl steal from autoconf 2.69
# UD_FC_MODULE_EXTENSION
# ----------------------
# Find the Fortran 90 module file extension.  The module extension is stored
# in the variable FC_MODEXT and empty if it cannot be determined.  The result
# or "unknown" is cached in the cache variable ac_cv_fc_module_ext.
AC_DEFUN([UD_FC_MODULE_EXTENSION],
[AC_CACHE_CHECK([Fortran 90 module extension], [ac_cv_fc_module_ext],
[AC_LANG_PUSH(Fortran)
mkdir conftest.dir
cd conftest.dir
ac_cv_fc_module_ext=unknown
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [ac_cv_fc_module_ext=`ls | sed -n 's,conftest_module\.,,p'`
   if test x$ac_cv_fc_module_ext = x; then
dnl Some F90 compilers use upper case characters for the module file name.
     ac_cv_fc_module_ext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
   fi])
cd ..
rm -rf conftest.dir
AC_LANG_POP(Fortran)
])
FC_MODEXT=$ac_cv_fc_module_ext
if test "$FC_MODEXT" = unknown; then
  FC_MODEXT=
fi
AC_SUBST([FC_MODEXT])dnl
])

dnl steal from autoconf 2.69
# UD_FC_MODULE_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ---------------------------------------------------------------------
# Find a flag to include Fortran 90 modules from another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_flag variable.
# It may contain significant trailing whitespace.
#
# Known flags:
# gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated), -Jdir for writing)
# g95: -I dir (-fmod=dir for writing)
# SUN: -Mdir, -M dir (-moddir=dir for writing;
#                     -Idir for includes is also searched)
# HP: -Idir, -I dir (+moddir=dir for writing)
# IBM: -Idir (-qmoddir=dir for writing)
# Intel: -Idir -I dir (-mod dir for writing)
# Absoft: -pdir
# Lahey: -mod dir
# Cray: -module dir, -p dir (-J dir for writing)
#       -e m is needed to enable writing .mod files at all
# Compaq: -Idir
# NAGWare: -I dir
# PathScale: -I dir  (but -module dir is looked at first)
# Portland: -module dir (first -module also names dir for writing)
# Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.', then -I)
#                    (-Am indicates how module information is saved)
AC_DEFUN([UD_FC_MODULE_FLAG],[
AC_CACHE_CHECK([Fortran 90 module inclusion flag], [ac_cv_fc_module_flag],
[AC_LANG_PUSH([Fortran])
ac_cv_fc_module_flag=unknown
mkdir conftest.dir
cd conftest.dir
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [cd ..
   ac_fc_module_flag_FCFLAGS_save=$FCFLAGS
   # Flag ordering is significant for gfortran and Sun.
   for ac_flag in -M -I '-I ' '-M ' -p '-mod ' '-module ' '-Am -I'; do
     # Add the flag twice to prevent matching an output flag.
     FCFLAGS="$ac_fc_module_flag_FCFLAGS_save ${ac_flag}conftest.dir ${ac_flag}conftest.dir"
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ac_cv_fc_module_flag="$ac_flag"])
     if test "$ac_cv_fc_module_flag" != unknown; then
       break
     fi
   done
   FCFLAGS=$ac_fc_module_flag_FCFLAGS_save
])
rm -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_flag" != unknown; then
  FC_MODINC=$ac_cv_fc_module_flag
  $1
else
  FC_MODINC=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag for module search path])])
fi
AC_SUBST([FC_MODINC])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODINC in #(
  *\ ) FC_MODINC=$FC_MODINC'${ac_empty}' ;;
esac])dnl
])


dnl steal from autoconf 2.69
# UD_FC_MODULE_OUTPUT_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------------
# Find a flag to write Fortran 90 module information to another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_output_flag variable.
# It may contain significant trailing whitespace.
#
# For known flags, see the documentation of AC_FC_MODULE_FLAG above.
AC_DEFUN([UD_FC_MODULE_OUTPUT_FLAG],[
AC_CACHE_CHECK([Fortran 90 module output flag], [ac_cv_fc_module_output_flag],
[AC_LANG_PUSH([Fortran])
mkdir conftest.dir conftest.dir/sub
cd conftest.dir
ac_cv_fc_module_output_flag=unknown
ac_fc_module_output_flag_FCFLAGS_save=$FCFLAGS
# Flag ordering is significant: put flags late which some compilers use
# for the search path.
for ac_flag in -J '-J ' -fmod= -moddir= +moddir= -qmoddir= '-mod ' \
	      '-module ' -M '-Am -M' '-e m -J '; do
  FCFLAGS="$ac_fc_module_output_flag_FCFLAGS_save ${ac_flag}sub"
  AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
    [cd sub
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ac_cv_fc_module_output_flag="$ac_flag"])
     cd ..
     if test "$ac_cv_fc_module_output_flag" != unknown; then
       break
     fi])
done
FCFLAGS=$ac_fc_module_output_flag_FCFLAGS_save
cd ..
rm -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_output_flag" != unknown; then
  FC_MODOUT=$ac_cv_fc_module_output_flag
  $1
else
  FC_MODOUT=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag to write module information to])])
fi
AC_SUBST([FC_MODOUT])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODOUT in #(
  *\ ) FC_MODOUT=$FC_MODOUT'${ac_empty}' ;;
esac])dnl
])
