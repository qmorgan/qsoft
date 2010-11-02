AC_INIT(Astro,1.0)
AC_COPYRIGHT(Copyright Edward Rosten (2009))


orig_CXXFLAGS="$CXXFLAGS"

AC_PROG_AWK
AC_PROG_CXX

AC_LANG(C++)

################################################################################
#
# Useful macros
#

dnl APPEND(var, value)
dnl This appends vale to a shell variable var
define(APPEND, [$1="$$1 $2"])



dnl TEST_AND_SET_CXXFLAG(flag, [program], [run])
dnl
dnl This attempts to compile a and run program with a certain compiler flag.
dnl If no program is given, then the minimal C++ program is compiled, and 
dnl this tests just the validity of the compiler flag. 
dnl
define([TEST_AND_SET_CXXFLAG],[
	if test "$3" == ""
	then
		AC_MSG_CHECKING([if compiler flag $1 works])	
	else
		AC_MSG_CHECKING([$3])	
	fi
	save_CXXFLAGS="$CXXFLAGS"
	APPEND(CXXFLAGS, [$1])

	m4_if([$2],[],[pushdef(prog, [int main(){}])], [pushdef(prog, [$2])])

	m4_if([$4],[run],
		[AC_RUN_IFELSE([prog], [cvd_conf_test=1],[cvd_conf_test=0], [cvd_conf_test=0])],
		[AC_COMPILE_IFELSE([prog], [cvd_conf_test=1],[cvd_conf_test=0])]
	)
	
	
	popdef([prog])

	if test $cvd_conf_test = 1
	then
		AC_MSG_RESULT(yes)
		ts_success=yes
	else
		AC_MSG_RESULT(no)
		CXXFLAGS="$save_CXXFLAGS"
		ts_success=no
	fi
])


TEST_AND_SET_CXXFLAG(-Wall)
TEST_AND_SET_CXXFLAG(-Wextra)
TEST_AND_SET_CXXFLAG(-W)

################################################################################
#
# Check for libcvd
#
APPEND(LDFLAGS, -lcvd)
AC_MSG_CHECKING(for libCVD)
comp=no
AC_COMPILE_IFELSE([
	#include <cvd/fast_corner.h>
	#include <cvd/image_io.h>
	#include <cvd/cpu_hacks.h>
	int main(){ CVD::Image<CVD::byte> i = CVD::img_load(""); img_save(i, ""); return 0;}
], [comp=yes], [comp=no])

AC_MSG_RESULT([$comp])

if test $comp = no 
then
	AC_MSG_CHECKING(If libcvd needs -msse2)
	APPEND(CXXFLAGS,[ -msse2])

	AC_COMPILE_IFELSE([
		#include <cvd/fast_corner.h>
		#include <cvd/image_io.h>
		int main(){ CVD::Image<CVD::byte> i = CVD::img_load(""); return 0;}
	], [comp=yes], [comp=no])
	AC_MSG_RESULT($comp)
fi

if test $comp = no
then
	AC_MSG_ERROR([A working copy of libCVD is required for this program.
Get libCVD here: http://savannah.nongnu.org/projects/libcvd
])
fi

################################################################################
#
# Check for TooN
#
toon=yes
AC_CHECK_HEADERS(TooN/TooN.h, [], [toon=no])

comp=no
AC_MSG_CHECKING([TooN])
AC_COMPILE_IFELSE([
	#include <TooN/TooN.h>
	int main(){TooN::Vector<2> v;return 0;}
], [comp=yes], [toon=no])

AC_MSG_RESULT($comp)

if test $toon = no
then
	AC_MSG_ERROR([A working copy of TooN is required for this program.
Get TooN here: http://savannah.nongnu.org/projects/TooN
])
fi


################################################################################
#
# Check for GVars
#
gvars=yes
AC_CHECK_HEADERS(gvars3/instances.h, [], [gvars=no])
AC_CHECK_LIB(GVars3, main, [], [gvars=no])

AC_ARG_ENABLE(headless, [AS_HELP_STRING([--enable-headless-gvars],[Use headless GVars (useful for batch jobs)])])

if test "$enable_headless" = yes
then
	LDFLAGS="$LDFLAGS -lGVars3_headless"
else
	LDFLAGS="$LDFLAGS -lGVars3"
fi


comp=no
AC_MSG_CHECKING([gvars])
AC_COMPILE_IFELSE([
	#include <gvars3/instances.h>
	int main(){GVars3::GUI.LoadFile("");}
], [comp=yes], [gvars=no])

AC_MSG_RESULT($comp)

if test $gvars = no
then
	AC_MSG_ERROR([A working copy of GVars3 is required for this program.
Get GVars3 here: http://savannah.nongnu.org/projects/libcvd
])
fi

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],[compile with debugging enabled])])
if test "$enable_debug" = yes
then
	CXXFLAGS="$CXXFLAGS -DSAFE_VECTOR -DTOON_DEBUG -DCVD_IMAGE_DEBUG"
else
	TEST_AND_SET_CXXFLAG(-O3)
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
fi


AC_OUTPUT(Makefile)
