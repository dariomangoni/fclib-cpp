/* FCLIB Copyright (C) 2011--2016 FClib project
*
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
* Contact: fclib-project@lists.gforge.inria.fr
*/
/*
* fctst.c
* ----------------------------------------------
* frictional contact test
*/

#include "fclib_class.hpp"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <algorithm>


/* generate random sparse matrix */
void random_matrix(std::shared_ptr<fclib::fclib_matrix_CPP> mat, int m, int n, fclib::fclib_matrix_CPP::mat_format_t mat_format, bool create_info)
{
	int nz_temp;
	auto nzmax_temp = std::min((m + n) + m*n / (10 + rand() % 10), m*n);
	switch (mat_format)
	{
	case fclib::fclib_matrix_CPP::mat_format_t::COO:
		nz_temp = nzmax_temp;
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSR:
		nz_temp = -2;
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSC:
		nz_temp = -1;
		break;
	default:
		nz_temp = 0;
	}

	mat->Resize(m, n, nz_temp, nzmax_temp);

	switch (mat_format)
	{
	case fclib::fclib_matrix_CPP::mat_format_t::COO:
		for (auto j = 0; j < mat->GetNZmax(); j++)
		{
			mat->SetElementP(j,rand() % mat->GetRows());
			mat->SetElementI(j, rand() % mat->GetColumns());
		}
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSC:
	{
		auto l = mat->GetNZmax() / mat->GetColumns();
		mat->SetElementP(0, 0);
		for (auto j = 1; j < mat->GetColumns(); j++) mat->SetElementP(j, mat->GetElementP(j-1) + l);
		for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementI(j, rand() % mat->GetColumns());
	}
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSR:
	{
		auto l = mat->GetNZmax() / mat->GetRows();
		mat->SetElementP(0, 0);
		for (auto j = 1; j < mat->GetRows(); j++) mat->SetElementP(j, mat->GetElementP(j-1) + l);
		for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementI(j, rand() % mat->GetRows());
	}
		break;
	default: ;
	}



	for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementX(j, static_cast<double>(rand()) / static_cast<double>(RAND_MAX));

	if (create_info)
		mat->SetComment("A random matrix");

}

/* generate random vector */
void random_vector(std::shared_ptr<std::vector<double>> vect, int n)
{
	vect->resize(n);
	for (auto row_sel = 0; row_sel<vect->size(); ++ row_sel)
		(*vect)[row_sel] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

}


/* generate random global problem */
void random_global_problem(fclib::fclib_global_problem_CPP& problem_out,
						int global_dofs,
						int contact_points,
						int neq,
						int spacedim,
						fclib::fclib_matrix_CPP::mat_format_t mat_format,
						bool create_G,
						bool create_info)
{
	problem_out.SetSpaceDim(spacedim);

	random_matrix(problem_out.matM, global_dofs, global_dofs, mat_format, create_info);
	random_matrix(problem_out.matH, global_dofs, spacedim*contact_points, mat_format, create_info);
	random_vector(problem_out.mu, contact_points);
	random_vector(problem_out.f, global_dofs);
	random_vector(problem_out.w, spacedim*contact_points);

	if (neq == 0 || !create_G)
	{
		problem_out.matG.reset();
		problem_out.b.reset();
	}
	else
	{
		random_matrix(problem_out.matG, global_dofs, neq, mat_format, create_info);
		random_vector(problem_out.b, problem_out.matG->GetColumns());
	}

	if (create_info)
	{
		problem_out.SetTitle(std::string("A random global problem"));
		problem_out.SetDescription(std::string("With random matrices"));
		problem_out.SetMathInfo(std::string("And fake math"));
	}


}

void random_global_solutions(fclib::fclib_global_problem_CPP& global_prob, fclib::fclib_solution_CPP& sol_vect)
{

		random_vector(sol_vect.v, global_prob.matM->GetColumns());
		random_vector(sol_vect.u, global_prob.matH->GetColumns());
		random_vector(sol_vect.r, global_prob.matH->GetColumns());

		if (global_prob.matG)
			random_vector(sol_vect.l, global_prob.matG->GetColumns());

		else
			sol_vect.l.reset();
}

/* generate random global solutions */
void random_global_solutions(fclib::fclib_global_problem_CPP& global_prob, std::vector<fclib::fclib_solution_CPP>& sol_vect, int num_sol)
{
	sol_vect.resize(num_sol);
	for (auto i = 0; i < sol_vect.size(); i++)
	{
		random_global_solutions(global_prob, sol_vect[i]);
	}
}

/* generate random local problem */
void random_local_problem(fclib::fclib_local_problem_CPP& local_prob, int contact_points, int neq, int spacedim, fclib::fclib_matrix_CPP::mat_format_t mat_format, bool create_info)
{

	local_prob.SetSpaceDim(spacedim);

	random_matrix(local_prob.matW, spacedim*contact_points, spacedim*contact_points,mat_format, create_info);

	if (neq && rand() % 2)
	{
		random_matrix(local_prob.matV, spacedim*contact_points, neq, mat_format, create_info);
		random_matrix(local_prob.matR, neq, neq,mat_format, create_info);
		random_vector(local_prob.s, neq);
	}
	else
	{
		local_prob.matV.reset();
		local_prob.matR.reset();
		local_prob.s.reset();
	}

	random_vector(local_prob.mu, contact_points);
	random_vector(local_prob.q, spacedim*contact_points);

	if (create_info)
	{
		local_prob.SetTitle(std::string("A random global problem"));
		local_prob.SetDescription(std::string("With random matrices"));
		local_prob.SetMathInfo(std::string("And fake math"));
	}

}


void random_local_solutions(fclib::fclib_local_problem_CPP& local_prob, fclib::fclib_solution_CPP& sol_vect)
{

		sol_vect.v.reset();
		random_vector(sol_vect.u, local_prob.matW->GetColumns());
		random_vector(sol_vect.r, local_prob.matW->GetColumns());

		if (local_prob.matR)
			random_vector(sol_vect.l, local_prob.matR->GetColumns());
		else
			sol_vect.l.reset();
}

/* generate random local solutions */
void random_local_solutions(fclib::fclib_local_problem_CPP& local_prob, std::vector<fclib::fclib_solution_CPP>& sol_vect, int num_sol)
{
	sol_vect.resize(num_sol);
	for (auto i = 0; i < sol_vect.size(); i++)
	{
		random_local_solutions(local_prob, sol_vect[i]);
	}
}

///* compare matrix infos */
static int compare_matrix_infos(std::shared_ptr<fclib::fclib_matrix_CPP> a, std::shared_ptr<fclib::fclib_matrix_CPP> b)
{
	if ((a == nullptr || a->GetRows() < 0) && (b == nullptr || b->GetRows() < 0))
		return 1;
	if ((a == nullptr || a->GetRows() < 0) != (b == nullptr || b->GetRows() < 0))
		return 0;

	if (a->GetComment() != b->GetComment() ||
		a->GetConditioning() != b->GetConditioning() ||
		a->GetDeterminant() != b->GetDeterminant() ||
		a->GetRank() != b->GetRank())
		return 0;

	return 1;
}

/* compare two matrices */
static int compare_matrices(std::shared_ptr<fclib::fclib_matrix_CPP> a, std::shared_ptr<fclib::fclib_matrix_CPP> b)
{

	if ((a == nullptr || a->GetRows() < 0) && (b == nullptr || b->GetRows() < 0))
		return 1;
	if ((a == nullptr || a->GetRows() < 0) != (b == nullptr || b->GetRows() < 0))
		return 0;

	if (a->GetNZmax() != b->GetNZmax() ||
		a->GetColumns() != b->GetColumns() ||
		a->GetRows() != b->GetRows() ||
		a->GetNZ() != b->GetNZ())
	{
		fprintf(stderr,
			"ERROR: dimensions of matrix differ:\n"
			"a->GetNZmax() = %d, b->GetNZmax() = %d\n"
			"a->GetColumns() = %d, b->GetColumns() = %d\n"
			"a->GetRows() = %d, b->GetRows() = %d\n"
			"a->GetNZ() = %d, b->GetNZ() = %d\n",
			a->GetNZmax(), b->GetNZmax(),
			a->GetColumns(), b->GetColumns(),
			a->GetRows(), b->GetRows(),
			a->GetNZ(), b->GetNZ());
		return 0;
	}

	if (a->GetNZ() >= 0)
	{
		for (auto i = 0; i < a->GetNZmax(); i++)
		{
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For matrix in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					 i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For this matrix in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					 i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
		}
	}
	else if (a->GetNZ() == -1)
	{
		for (auto i = 0; i < a->GetColumns() + 1; i++)
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For this matrix in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					 i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}

		for (auto i = 0; i < a->GetNZmax(); i++)
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For this matrix in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					 i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
	}
	else if (a->GetNZ() == -2)
	{
		for (auto i = 0; i < a->GetRows() + 1; i++)
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For this matrix in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					 i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}

		for (auto i = 0; i < a->GetNZmax(); i++)
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For this matrix in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					 i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
	}

	for (auto i = 0; i < a->GetNZmax(); i++)
		if ((*(a->x))[i] != (*(b->x))[i])
		{
			fprintf(stderr,
				"ERROR: For this matrix in {a, b} (*(a->x)) [%d] != (*(b->x)) [%d] => %g != %g\n",
				 i, i, (*(a->x))[i], (*(b->x))[i]);
			return 0;
		}

	if (a->HasInfo() && b->HasInfo() && a->GetComment() != b->GetComment())
	{
		fprintf(stderr, "ERROR: matrix infos differ\n");
		return 0;
	}

	return 1;
}

/* compare two vectors */
static int compare_vectors(std::shared_ptr<std::vector<double>> a, std::shared_ptr<std::vector<double>> b)
{
	bool is_different = true;
	if ((a == nullptr || a->size() == 0) && (a == nullptr || a->size() == 0))
		return 1;
	if ((a == nullptr || a->size() == 0) != (a == nullptr || a->size() == 0))
		return 0;

	for (auto row_sel = 0; row_sel<a->size(); ++row_sel)
	{
		if ((*a)[row_sel] != (*b)[row_sel])
			return false;
	}
	return is_different;
}


/* compare global problems */
static int compare_global_problems(fclib::fclib_global_problem_CPP& a, fclib::fclib_global_problem_CPP& b)
{
	if (!compare_matrices(a.matM, b.matM) ||
		!compare_matrices(a.matH, b.matH) ||
		!compare_matrices(a.matG, b.matG) ||
		!compare_vectors(a.mu, b.mu) ||
		!compare_vectors(a.f, b.f) ||
		!compare_vectors(a.b, b.b) ||
		!compare_vectors(a.w, b.w) ||
		a.GetSpaceDim() != b.GetSpaceDim()
		//|| !compare_infos((*(a->i))nfo, (*(b->i))nfo)
		) return 0;

	return 1;
}

/* compare local problems */
static int compare_local_problems(fclib::fclib_local_problem_CPP& a, fclib::fclib_local_problem_CPP& b)
{
	if (!compare_matrices(a.matW, b.matW) ||
		!compare_matrices(a.matV, b.matV) ||
		!compare_matrices(a.matR, b.matR) ||
		!compare_vectors(a.mu, b.mu) ||
		!compare_vectors(a.q, b.q) ||
		!compare_vectors(a.s, b.s) ||
		a.GetSpaceDim() != b.GetSpaceDim()
		//||!compare_infos((*(a->i))nfo, (*(b->i))nfo)
		) return 0;

	return 1;
}

/* compare solutions */
static int compare_solutions(fclib::fclib_solution_CPP& a, fclib::fclib_solution_CPP& b)
{
	if (!compare_vectors(a.v, b.v) ||
		!compare_vectors(a.u, b.u) ||
		!compare_vectors(a.r, b.r) ||
		!compare_vectors(a.l, b.l)) return 0;

	return 1;
}

int main(int argc, char **argv)
{
	srand(static_cast<unsigned int>(time(nullptr)));

	int space_dimension = 2;
	int contact_points = 100;
	int neq = 100;
	auto mat_format = fclib::fclib_matrix_CPP::mat_format_t::CSR;
	int num_guesses_out = 3;
	bool create_info = true;

	// only global
	int global_dofs = 100;
	bool create_G = true;

	{


		fclib::fclib_global_problem_CPP problem_out, problem_in;
		fclib::fclib_solution_CPP solution_out, solution_in;
		std::vector<fclib::fclib_solution_CPP> guesses_out, guesses_in;
		guesses_out.resize(num_guesses_out);


		random_global_problem(problem_out, global_dofs, contact_points,	neq, space_dimension, mat_format, create_G, create_info);
		random_global_solutions(problem_out, solution_out);
		random_global_solutions(problem_out, guesses_out, num_guesses_out);


		problem_out.write_problem("output_file.hdf5");
		solution_out.write_solution("output_file.hdf5");
		fclib::fclib_solution_CPP::write_guesses(guesses_out, "output_file.hdf5", num_guesses_out);

		problem_in.read_problem("output_file.hdf5");
		solution_in.read_solution("output_file.hdf5");
		int num_guesses_in = fclib::fclib_solution_CPP::read_guesses(guesses_in, "output_file.hdf5");

		printf("Comparing written and read global problem data ...\n");
		assert(compare_matrices(problem_out.matM, problem_in.matM));
		assert(compare_matrices(problem_out.matH, problem_in.matH));
		assert(compare_matrices(problem_out.matG, problem_in.matG));
		assert(compare_vectors(problem_out.mu, problem_in.mu));
		assert(compare_vectors(problem_out.f, problem_in.f));
		assert(compare_vectors(problem_out.b, problem_in.b));
		assert(compare_vectors(problem_out.w, problem_in.w));
		assert(problem_out.GetSpaceDim() == problem_in.GetSpaceDim());
		assert(compare_global_problems(problem_out, problem_in) && "ERROR: written/read problem comparison failed");
		assert(compare_solutions(solution_out, solution_out) && "ERROR: written/read solution comparison failed");
		assert(num_guesses_out == num_guesses_in && "ERROR: numbers of written and read guesses differ");
		for (auto i = 0; i < num_guesses_in; i++)
		{
			assert(compare_solutions(guesses_in[i], guesses_out[i]) && "ERROR: written/read guess comparison failed");
		}

		printf("All comparisons PASSED\n");
		remove("output_file.hdf5");

	}

	{



		fclib::fclib_local_problem_CPP problem_out, problem_in;
		fclib::fclib_solution_CPP solution_out, solution_in;
		std::vector<fclib::fclib_solution_CPP> guesses_out, guesses_in;
		guesses_out.resize(num_guesses_out);

		random_local_problem(problem_out, contact_points, neq, space_dimension, mat_format, create_info);
		random_local_solutions(problem_out, solution_out);
		random_local_solutions(problem_out, guesses_out, num_guesses_out);


		problem_out.write_problem("output_file.hdf5");
		solution_out.write_solution("output_file.hdf5");
		fclib::fclib_solution_CPP::write_guesses(guesses_out, "output_file.hdf5", num_guesses_out);

		problem_in.read_problem("output_file.hdf5");

		solution_in.read_solution("output_file.hdf5");
		int num_guesses_in = fclib::fclib_solution_CPP::read_guesses(guesses_in, "output_file.hdf5");

		printf("Comparing written and read local problem data ...\n");

		//printf("Computing merit function ...\n");
		//double error1 = fclib_merit_local(problem, MERIT_1, solution);
		//double error2 = fclib_merit_local(p, MERIT_1, s);
		//printf("Error for initial problem = %12.8e\n", error1);
		//printf("Error for read problem = %12.8e\n", error2);

		assert(compare_matrices(problem_out.matW, problem_in.matW));
		assert(compare_matrices(problem_out.matV, problem_in.matV));
		assert(compare_matrices(problem_out.matR, problem_in.matR));
		assert(compare_vectors(problem_out.mu, problem_in.mu));
		assert(compare_vectors(problem_out.q, problem_in.q));
		assert(compare_vectors(problem_out.s, problem_in.s));
		assert(problem_out.GetSpaceDim() == problem_in.GetSpaceDim());

		assert(compare_solutions(solution_out, solution_out) && "ERROR: written/read solution comparison failed");
		assert(num_guesses_out == num_guesses_in && "ERROR: numbers of written and read guesses differ");
		for (auto i = 0; i < num_guesses_in; i++)
		{
			assert(compare_solutions(guesses_in[i], guesses_out[i]) && "ERROR: written/read guess comparison failed");
		}

		printf("All comparisons PASSED\n");

		remove("output_file.hdf5");

	}


	return 0;
}
