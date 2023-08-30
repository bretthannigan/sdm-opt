/* Include files */

#include "MOD1_CT_cgxe.h"
#include "m_JO2Fxs0ruK9ktkDoGiL8UC.h"

unsigned int cgxe_MOD1_CT_method_dispatcher(SimStruct* S, int_T method, void
  * data)
{
  if (ssGetChecksum0(S) == 926501644 &&
      ssGetChecksum1(S) == 2230881906 &&
      ssGetChecksum2(S) == 368894390 &&
      ssGetChecksum3(S) == 3689493567) {
    method_dispatcher_JO2Fxs0ruK9ktkDoGiL8UC(S, method, data);
    return 1;
  }

  return 0;
}
