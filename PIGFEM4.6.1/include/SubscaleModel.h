#ifndef _SUBSCALE_H_
#define _SUBSCALE_H_
#include <vector>
#include "material.h"

//Predeclarations
class Problem;
class Mesh;


class SubscaleModel
{
	protected:

		Mesh* _mesh;
		Problem* _prob;
		bool _strain_set;
		std::vector<double> _origin;

	public:

		SubscaleModel();
		SubscaleModel(int *argc, char***argv);
		virtual ~SubscaleModel();

		Mesh* get_mesh() {return _mesh;};
		Problem* get_problem() {return _prob;};

		virtual void set_macro_strain(std::vector<double> strain) = 0;

		void homogenize_stress_strain(const std::vector<double>& curr_macro_strain, std::vector<double>& avg_stress, std::vector<double>& perturbation_strain);

		void write_stress_strain(std::string filename, id_type step_iter, const std::vector<double>& current_macro_strain,
								 const std::vector<double>& stress, const std::vector<double>& perturbation_strain);

		std::vector<double> compute_macro_displacement(const std::vector<double>& curr_macro_strain, const std::vector<double> coords);

		void init();
		void solve();

	protected:

		void stress_strain_kernel(std::vector<double>& strain, std::vector<double>& stress,
								  const std::vector<std::vector<double> >& shape_grad,
								  Material* mat, Material::input_params& input,
								  const std::vector<double>& elem_U_curr,
								  const std::vector<double>& curr_macro_strain);

		void determine_origin();
};


#endif