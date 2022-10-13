#ifndef _SUBSCALE_NONLINEAR_H_
#define _SUBSCALE_NONLINEAR_H_
#include "SubscaleModel.h"



class SubscaleModelNonlinear : public SubscaleModel
{
	public:

		SubscaleModelNonlinear(int *argc, char***argv);

		virtual void set_macro_strain(std::vector<double> strain);
};


#endif