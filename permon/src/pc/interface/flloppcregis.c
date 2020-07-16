
#include <permonpc.h>

FLLOP_EXTERN PetscErrorCode PCCreate_Dual(PC);
  
#undef __FUNCT__  
#define __FUNCT__ "FllopPCRegisterAll"
PetscErrorCode  FllopPCRegisterAll()
{
  PetscFunctionBegin;
  TRY( PCRegister(PCDUAL, PCCreate_Dual) );
  PetscFunctionReturn(0);
}
