#ifndef _SUBSCALE_LINEAR_H_
#define _SUBSCALE_LINEAR_H_
#include "SubscaleModel.h"



class SubscaleModelLinear : public SubscaleModel
{
	public:

		SubscaleModelLinear(int *argc, char***argv);
		virtual ~SubscaleModelLinear();

		virtual void set_macro_strain(std::vector<double> strain);
};


#endif