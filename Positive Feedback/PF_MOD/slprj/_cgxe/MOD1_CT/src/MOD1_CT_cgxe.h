#ifndef __MOD1_CT_cgxe_h__
#define __MOD1_CT_cgxe_h__

/* Include files */
#include "simstruc.h"
#include "rtwtypes.h"
#include "multiword_types.h"
#include "emlrt.h"
#include "covrt.h"
#include "cgxert.h"
#include "C:\Program Files\MATLAB\R2017a\toolbox\shared\dsp\vision\matlab\include\HostLib_FFT.h"
#include "C:\Program Files\MATLAB\R2017a\toolbox\shared\dspblks\extern\include\HostLib_rtw.h"
#include "blas.h"
#define rtInf                          (mxGetInf())
#define rtMinusInf                     (-(mxGetInf()))
#define rtNaN                          (mxGetNaN())
#define rtIsNaN(X)                     ((int)mxIsNaN(X))
#define rtIsInf(X)                     ((int)mxIsInf(X))

extern unsigned int cgxe_MOD1_CT_method_dispatcher(SimStruct* S, int_T method,
  void* data);

#endif
