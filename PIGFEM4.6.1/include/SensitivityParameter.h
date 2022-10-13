/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _SENSITIVITY_PARAMETER_H_
#define _SENSITIVITY_PARAMETER_H_
#include "common.h"
#include "petscvec.h"
#include "petscmat.h"
#include "material.h"
#include "DenseMatrix.h"

class Elem;
class Problem;
class NodalData;

class SensitivityParameter
{
	protected:
	
		// A unique parameter id for each parmeter being analyzed
		id_type _id;
		Problem* _prob;
		std::vector<bool> _assemble_elem;
		std::vector<SensitivityParameter*> _parameters; // Note: This is a special case for when we are assembling an entire matrix of rhs

		PetscErrorCode set_dP(pvector& dPdd, std::vector<double>& dP_el, std::vector<id_type>& elem_dofs, id_type ngfd);


		// The functions that need to be implemented to update internal variable senitivities
		virtual void KernelVolumetricUpdate(const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											Material* mat, Material::sensitivity_input_params& input,
									 		const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr) {err_message("Must create the function KernelVolumetricUpdate to have a valid SensitivityParameter");};
		// B-matrix variant
		virtual void KernelVolumetricUpdate(const std::vector<double>& shape, const DenseMatrix<double>& B,
											Material* mat, Material::sensitivity_input_params& input,
									 		const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr) {err_message("Must create the function KernelVolumetricUpdate to have a valid SensitivityParameter");};
		virtual void KernelCohesiveUpdate(const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
										  Material* mat, Material::sensitivity_input_params& input,
										  const std::vector<double>& coh_U_curr, std::vector<double>& dcoh_U_curr) {err_message("Must create the function KernelCohesiveUpdate to have a valid SensitivityParameter");};



		// Handle the elemental assembly
		virtual void assemble_elem_dPdd(Elem* el, std::vector<double>& dP_el, std::vector<double>& elem_sol, Problem* prob) {err_message("Must create the function assemble_elem_dPdd to have a valid SensitivityParameter");};;
		void assemble_elem_dPdd_mat(Elem* el, DenseMatrix<double>& dP_el, std::vector<double>& elem_sol, Problem* prob);
		void updateSensitivityElemISVs(Elem* el, std::vector<double>& elem_sol, std::vector<double>& delem_sol, Problem* prob);
		void updateSensitivityElemISVs_mat(Elem* el, std::vector<double>& elem_sol, std::vector<std::vector<double> >& delem_sol, Problem* prob);

	public:

		SensitivityParameter();
		SensitivityParameter(std::vector<SensitivityParameter*>& param) {_parameters = param;};
		virtual ~SensitivityParameter() {};


		void set_id(id_type id) {_id = id;};
		id_type get_id() {return _id;};



		virtual sensitivity_parameter_type get_type() {return MATRIX_ASSEMBLER;};
		virtual void attachProblem(Problem* prob) {_prob = prob;};
		virtual void init();
		virtual void setup() {};
		void set_material_element_assemblies(Problem* prob, std::vector<std::string>& names);
		void set_shape_element_assemblies(Problem* prob, std::vector<id_type>& ids);
		bool assembleElem(id_type l_elem) {return _assemble_elem[l_elem];};



		virtual PetscErrorCode assemble_dPdd(pvector& dPdd, Problem* prob);
		virtual PetscErrorCode assemble_load_derivatives(Vec* dUpdd, Vec* dFfdd, Problem* prob) {return 0;}; // Load derivatives are different than material or shape I think
		virtual PetscErrorCode assemble_dPdd_mat(Mat* dPfdd, Mat* dPpdd, std::vector<SensitivityParameter*> parameters, Problem* prob);
		virtual PetscErrorCode assemble_load_derivatives_mat(Mat* dUpdd, Mat* dFfdd, std::vector<SensitivityParameter*> parameters, Problem* prob) {return 0;}; // Load derivatives are different than material or shape I think

		void updateSensitivityISVs(Problem* prob, NodalData* dU_dd);
		void updateSensitivityISVs_mat(Problem* prob, Mat* dU_dd, Mat* dUp_dd = NULL); // Can exclude Up portion
};





#endif