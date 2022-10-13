#include "SubscaleModelNonlinear.h"
#include "Mesh.h"
#include "ProblemNonlinearStructural_Multiscale.h"


SubscaleModelNonlinear::SubscaleModelNonlinear(int *argc, char***argv)
{
	_mesh = new Mesh(argc, argv);
	_mesh->set_periodic(true);
	_prob = new ProblemNonlinearStructural_Multiscale;
	_prob->attach_mesh(_mesh);
	_prob->attach_subscale_model(this);
}



void SubscaleModelNonlinear::set_macro_strain(std::vector<double> strain)
{
	_prob->set_vec_parameter("macro strain", strain);
	_strain_set = true;
}