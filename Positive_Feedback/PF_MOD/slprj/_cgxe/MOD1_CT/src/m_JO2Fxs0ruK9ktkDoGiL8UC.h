#ifndef __JO2Fxs0ruK9ktkDoGiL8UC_h__
#define __JO2Fxs0ruK9ktkDoGiL8UC_h__

/* Include files */
#include "simstruc.h"
#include "rtwtypes.h"
#include "multiword_types.h"
#include "slexec_vm_zc_functions.h"

/* Type Definitions */
#ifndef struct_md7278973d290b7e345807af41d4a7c9c9
#define struct_md7278973d290b7e345807af41d4a7c9c9

struct md7278973d290b7e345807af41d4a7c9c9
{
  int32_T S0_isInitialized;
  real_T W0_HostLib[137];
};

#endif                                 /*struct_md7278973d290b7e345807af41d4a7c9c9*/

#ifndef typedef_dsp_FFT_0
#define typedef_dsp_FFT_0

typedef struct md7278973d290b7e345807af41d4a7c9c9 dsp_FFT_0;

#endif                                 /*typedef_dsp_FFT_0*/

#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  real_T f1[2];
  real_T f2[2];
} struct_T;

#endif                                 /*typedef_struct_T*/

#ifndef typedef_b_struct_T
#define typedef_b_struct_T

typedef struct {
  char_T f1[6];
} b_struct_T;

#endif                                 /*typedef_b_struct_T*/

#ifndef typedef_c_struct_T
#define typedef_c_struct_T

typedef struct {
  char_T f1[4];
  char_T f2[8];
  char_T f3[6];
  char_T f4[7];
} c_struct_T;

#endif                                 /*typedef_c_struct_T*/

#ifndef typedef_d_struct_T
#define typedef_d_struct_T

typedef struct {
  char_T f1[6];
  char_T f2[6];
} d_struct_T;

#endif                                 /*typedef_d_struct_T*/

#ifndef typedef_e_struct_T
#define typedef_e_struct_T

typedef struct {
  char_T f1[2];
  char_T f2[9];
} e_struct_T;

#endif                                 /*typedef_e_struct_T*/

#ifndef struct_mdCaRKN9DEuZYBgsloVJSc4F
#define struct_mdCaRKN9DEuZYBgsloVJSc4F

struct mdCaRKN9DEuZYBgsloVJSc4F
{
  int32_T isInitialized;
  real_T pSampleRateInherit;
  creal_T pPeriodogramMatrix[51200];
  real_T pWindowData[51200];
  real_T pWindowPower;
  real_T pW[25601];
  real_T pNumAvgsCounter;
  real_T pNewPeriodogramIdx;
  dsp_FFT_0 pFFT;
  real_T pFrameCounter;
  real_T pFrameDelay;
  real_T pSR;
};

#endif                                 /*struct_mdCaRKN9DEuZYBgsloVJSc4F*/

#ifndef typedef_dsp_simulink_CrossSpectrumEstimator
#define typedef_dsp_simulink_CrossSpectrumEstimator

typedef struct mdCaRKN9DEuZYBgsloVJSc4F dsp_simulink_CrossSpectrumEstimator;

#endif                                 /*typedef_dsp_simulink_CrossSpectrumEstimator*/

#ifndef typedef_InstanceStruct_JO2Fxs0ruK9ktkDoGiL8UC
#define typedef_InstanceStruct_JO2Fxs0ruK9ktkDoGiL8UC

typedef struct {
  SimStruct *S;
  dsp_simulink_CrossSpectrumEstimator sysobj;
  boolean_T sysobj_not_empty;
  boolean_T isInitialized;
  creal_T Z[102400];
  creal_T P[51200];
  creal_T b_P[51200];
  real_T y[102400];
  real_T U0[102400];
  creal_T varargout_1[25601];
  real_T varargin_6[51200];
  real_T varargin_7[51200];
  real_T b_y[51200];
  real_T c_y[51200];
  real_T unusedExpr[25601];
  real_T F[25601];
  real_T a[51200];
  real_T b[51200];
  real_T b_unusedExpr[51200];
  real_T b_F[25601];
  creal_T x[51200];
  creal_T Pos_unscaled[25601];
  real_T w0[51200];
  void *emlrtRootTLSGlobal;
  real_T (*u0)[51200];
  real_T (*u1)[51200];
  creal_T (*b_y0)[25601];
} InstanceStruct_JO2Fxs0ruK9ktkDoGiL8UC;

#endif                                 /*typedef_InstanceStruct_JO2Fxs0ruK9ktkDoGiL8UC*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
extern void method_dispatcher_JO2Fxs0ruK9ktkDoGiL8UC(SimStruct *S, int_T method,
  void* data);

#endif
