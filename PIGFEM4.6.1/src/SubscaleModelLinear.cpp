#include "SubscaleModelLinear.h"
#include "Mesh.h"
#include "ProblemLinearElasticity.h"
#include "BodyLoadMacroStrainLinear.h"


SubscaleModelLinear::SubscaleModelLinear(int *argc, char***argv)
{
	_mesh = new Mesh(argc, argv);
	_mesh->set_periodic(true);
	_prob = new ProblemLinearElasticity;
	_prob->attach_mesh(_mesh);
	_prob->attach_subscale_model(this);

	// Set up the external body load
	BodyLoad* load = new BodyLoadMacroStrainLinear;
	_prob->add_body_load(load);

}
SubscaleModelLinear::~SubscaleModelLinear()
{
	delete *(_prob->body_loads_begin());
}



void SubscaleModelLinear::set_macro_strain(std::vector<double> strain)
{
	BodyLoad* load = *(_prob->body_loads_begin());
	load->set_vec_parameter("macro strain", strain);
	_strain_set = true;
}