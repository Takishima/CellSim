/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

#ifdef __clang__
#  define COMPILER_VERSION (__clang_major__ * 10000		\
			 + __clang_minor__ * 100		\
			 + __clang_patchlevel__)
#  if COMPILER_VERSION < 30100
#     error This code requires at least clang-3.1
#  endif 
#elif (defined __GNUC__)
#  define COMPILER_VERSION (__GNUC__ * 10000		\
		       + __GNUC_MINOR__ * 100	\
		       + __GNUC_PATCHLEVEL__)
#  if COMPILER_VERSION < 40700
#      error This code requires at least GCC-4.7.x
#  endif
#endif

#if (defined __GNUC__) || (defined __clang__)
# define GCC_DIAG_STR(s) #s 
# define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y)
# define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x) 
# define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
# define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push)		\
     GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x)) 
# define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop)
#else
# define GCC_DIAG_OFF(x) 
# define GCC_DIAG_ON(x) 
#endif

#ifdef __clang__
# define CLANG_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(clang diagnostic x)
# define CLANG_DIAG_OFF(x) CLANG_DIAG_PRAGMA(push)		\
     GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x)) 
# define CLANG_DIAG_ON(x) CLANG_DIAG_PRAGMA(pop)
#else
# define CLANG_DIAG_OFF(x)
# define CLANG_DIAG_ON(x)
#endif


#if (defined HALIDI_2009) and (defined HALIDI_2012)
#  error Cannot define HALIDI_2009 *and* HALIDI_2012
#elif (not defined HALIDI_2009) and (not defined HALIDI_2012)
#  error Must define one of HALIDI_2009 *or* HALIDI_2012
#endif

#ifdef PLC_ACTIVATION
#  define _PLC_CONST 1
#else
#  define _PLC_CONST 0
#endif

#ifdef ELECTRIC_ACTIVATION
#  define _ELEC_CONST 1
#else
#  define _ELEC_CONST 0
#endif

#ifdef POTASSIUM_CHORIDE_ACTIVATION
#  define _K_CL_CONST 1
#else
#  define _K_CL_CONST 0
#endif

#if (_PLC_CONST + _ELEC_CONST + _K_CL_CONST) > 1
#  error Must define only one of PLC_ACTIVATION, ELECTRIC_ACTIVATION and POTASSIUM_CHORIDE_ACTIVATION
#elif (_PLC_CONST + _ELEC_CONST + _K_CL_CONST) == 0
#  error Must define at least one of PLC_ACTIVATION, ELECTRIC_ACTIVATION and POTASSIUM_CHORIDE_ACTIVATION
#endif

#endif //CONFIG_HPP_INCLUDED
