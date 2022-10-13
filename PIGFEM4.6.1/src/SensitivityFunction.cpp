#include "SensitivityFunction.h"




SensitivityFunction::~SensitivityFunction()
{

}


SensitivityFunction* SensitivityFunction::allocate_and_copy()
{
	SensitivityFunction* new_func = allocate();
	copy(new_func);
	return new_func;
}