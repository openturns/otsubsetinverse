
#ifndef OTSUBSETINVERSE_PRIVATE_HXX
#define OTSUBSETINVERSE_PRIVATE_HXX

/* From http://gcc.gnu.org/wiki/Visibility */
/* Generic helper definitions for shared library support */
#if defined _WIN32 || defined __CYGWIN__
#define OTSUBSETINVERSE_HELPER_DLL_IMPORT __declspec(dllimport)
#define OTSUBSETINVERSE_HELPER_DLL_EXPORT __declspec(dllexport)
#define OTSUBSETINVERSE_HELPER_DLL_LOCAL
#else
#if __GNUC__ >= 4
#define OTSUBSETINVERSE_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
#define OTSUBSETINVERSE_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
#define OTSUBSETINVERSE_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define OTSUBSETINVERSE_HELPER_DLL_IMPORT
#define OTSUBSETINVERSE_HELPER_DLL_EXPORT
#define OTSUBSETINVERSE_HELPER_DLL_LOCAL
#endif
#endif

/* Now we use the generic helper definitions above to define OTSUBSETINVERSE_API and OTSUBSETINVERSE_LOCAL.
 * OTSUBSETINVERSE_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
 * OTSUBSETINVERSE_LOCAL is used for non-api symbols. */

#ifndef OTSUBSETINVERSE_STATIC /* defined if OT is compiled as a DLL */
#ifdef OTSUBSETINVERSE_DLL_EXPORTS /* defined if we are building the OT DLL (instead of using it) */
#define OTSUBSETINVERSE_API OTSUBSETINVERSE_HELPER_DLL_EXPORT
#else
#define OTSUBSETINVERSE_API OTSUBSETINVERSE_HELPER_DLL_IMPORT
#endif /* OTSUBSETINVERSE_DLL_EXPORTS */
#define OTSUBSETINVERSE_LOCAL OTSUBSETINVERSE_HELPER_DLL_LOCAL
#else /* OTSUBSETINVERSE_STATIC is defined: this means OT is a static lib. */
#define OTSUBSETINVERSE_API
#define OTSUBSETINVERSE_LOCAL
#endif /* !OTSUBSETINVERSE_STATIC */


#endif // OTSUBSETINVERSE_PRIVATE_HXX

