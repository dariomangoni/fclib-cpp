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


/* allocate matrix info */
//static struct fclib_matrix_info* matrix_info(struct fclib_matrix *mat, char *comment)
//{
//	struct fclib_matrix_info *info;
//
//	MM(info = malloc(sizeof(struct fclib_matrix_info)));
//	MM(info->comment = malloc(strlen(comment) + 1));
//	strcpy(info->comment, comment);
//	info->conditioning = rand();
//	info->determinant = rand();
//	info->rank = mat->GetRows();
//
//	return info;
//}

/* generate random sparse matrix */
void random_matrix(std::shared_ptr<fclib::fclib_matrix_CPP> mat, int m, int n, fclib::fclib_matrix_CPP::mat_format_t mat_format)
{
	int nz_temp;
	auto nzmax_temp = std::min((m + n) + m*n / (10 + rand() % 10), m*n);
	switch (mat_format)
	{
	case fclib::fclib_matrix_CPP::mat_format_t::COO:
		nz_temp = nzmax_temp;
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSR:
		nz_temp = -1;
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSC:
		nz_temp = -2;
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
	case fclib::fclib_matrix_CPP::mat_format_t::CSR:
	{
		auto l = mat->GetNZmax() / mat->GetColumns();
		mat->SetElementP(0, 0);
		for (auto j = 0; j < mat->GetColumns(); j++) mat->SetElementP(j + 1, mat->GetElementP(j) + l);
		for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementI(j, rand() % mat->GetColumns());
	}
		break;
	case fclib::fclib_matrix_CPP::mat_format_t::CSC:
	{
		auto l = mat->GetNZmax() / mat->GetRows();
		mat->SetElementP(0, 0);
		for (auto j = 0; j < mat->GetRows(); j++) mat->SetElementP(j + 1, mat->GetElementP(j) + l);
		for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementI(j, rand() % mat->GetRows());
	}
		break;
	default: ;
	}



	for (auto j = 0; j < mat->GetNZmax(); j++) mat->SetElementX(j, static_cast<double>(rand()) / static_cast<double>(RAND_MAX));

	//TODO: info
	//if (rand()) mat->info = matrix_info(mat, "A random matrix");
	//else mat->info = NULL;
}

/* generate random vector */
void random_vector(std::shared_ptr<std::vector<double>> vect, int n)
{
	vect->resize(n);
	for (auto row_sel = 0; row_sel<vect->size(); ++ row_sel)
		(*vect)[row_sel] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

}

///* allocate problem info */
//static struct fclib_info* problem_info(char *title, char *desc, char *math)
//{
//	struct fclib_info *info;
//
//	MM(info = malloc(sizeof(struct fclib_info)));
//	MM(info->title = malloc(strlen(title) + 1));
//	strcpy(info->title, title);
//	MM(info->description = malloc(strlen(desc) + 1));
//	strcpy(info->description, desc);
//	MM(info->math_info = malloc(strlen(math) + 1));
//	strcpy(info->math_info, math);
//
//	return info;
//}


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

	random_matrix(problem_out.matM, global_dofs, global_dofs, mat_format);
	random_matrix(problem_out.matH, global_dofs, spacedim*contact_points, mat_format);
	random_vector(problem_out.mu, contact_points);
	random_vector(problem_out.f, global_dofs);
	random_vector(problem_out.w, spacedim*contact_points);

	if (neq == 0 || !create_G)
	{
		problem_out.matG = nullptr;
		problem_out.b = nullptr;
	}
	else
	{
		random_matrix(problem_out.matG, global_dofs, neq, mat_format);
		random_vector(problem_out.b, problem_out.matG->GetColumns());
	}

	if (create_info)
		problem->info = problem_info("A random global problem", "With random matrices", "And fake math");
	else problem->info = NULL;

}

void random_global_solutions(fclib::fclib_global_problem_CPP& global_prob, fclib::fclib_solution_CPP& sol_vect)
{

		random_vector(sol_vect.v, global_prob.matM->GetColumns());
		random_vector(sol_vect.u, global_prob.matH->GetColumns());
		random_vector(sol_vect.r, global_prob.matH->GetColumns());

		if (global_prob.matG)
			random_vector(sol_vect.l, global_prob.matG->GetColumns());

		else
			sol_vect.l = nullptr;
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
void random_local_problem(fclib::fclib_local_problem_CPP& local_prob, int contact_points, int neq, int spacedim, fclib::fclib_matrix_CPP::mat_format_t mat_format)
{

	local_prob.SetSpaceDim(spacedim);

	random_matrix(local_prob.matW, spacedim*contact_points, spacedim*contact_points,mat_format);
	if (neq && rand() % 2)
	{
		random_matrix(local_prob.matV, spacedim*contact_points, neq, mat_format);
		random_matrix(local_prob.matR, neq, neq,mat_format);
		random_vector(local_prob.s, neq);
	}
	else
	{
		local_prob.matV = nullptr;
		local_prob.matR = nullptr;
		local_prob.s = nullptr;
	}
	random_vector(local_prob.mu, contact_points);
	random_vector(local_prob.q, spacedim*contact_points);

	//if (rand() % 2) problem->info = problem_info("A random local problem", "With random matrices", "And fake math");
	//else problem->info = NULL;

}


void random_local_solutions(fclib::fclib_local_problem_CPP& local_prob, fclib::fclib_solution_CPP& sol_vect)
{

		sol_vect.v = nullptr;
		random_vector(sol_vect.u, local_prob.matW->GetColumns());
		random_vector(sol_vect.r, local_prob.matW->GetColumns());

		if (local_prob.matR)
			random_vector(sol_vect.l, local_prob.matR->GetColumns());

		else
			sol_vect.l = nullptr;
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
//static int compare_matrix_infos(struct fclib_matrix_info *a, struct fclib_matrix_info *b)
//{
//	if (!a && !b) return 1;
//	else if ((a && !b) || (!a && b)) return 0;
//	else if (strcmp(a->comment, b->comment) != 0 ||
//		a->conditioning != b->conditioning ||
//		a->determinant != b->determinant ||
//		a->rank != b->rank) return 0;
//
//	return 1;
//}

/* compare two matrices */
static int compare_matrices(char *name, std::shared_ptr<fclib::fclib_matrix_CPP> a, std::shared_ptr<fclib::fclib_matrix_CPP> b)
{
	int i;

	if (!a && !b) return 1;
	else if ((a && !b) || (!a && b)) return 0;

	if (a->GetNZmax() != b->GetNZmax() ||
		a->GetColumns() != b->GetColumns() ||
		a->GetRows() != b->GetRows() ||
		a->GetNZ() != b->GetNZ())
	{
		fprintf(stderr,
			"ERROR: dimensions of matrix %s differ:\n"
			"a->GetNZmax() = %d, b->GetNZmax() = %d\n"
			"a->GetColumns() = %d, b->GetColumns() = %d\n"
			"a->GetRows() = %d, b->GetRows() = %d\n"
			"a->GetNZ() = %d, b->GetNZ() = %d\n", name,
			a->GetNZmax(), b->GetNZmax(),
			a->GetColumns(), b->GetColumns(),
			a->GetRows(), b->GetRows(),
			a->GetNZ(), b->GetNZ());
		return 0;
	}

	if (a->GetNZ() >= 0)
	{
		for (i = 0; i < a->GetNZmax(); i++)
		{
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					name, i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					name, i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
		}
	}
	else if (a->GetNZ() == -1)
	{
		for (i = 0; i < a->GetColumns() + 1; i++)
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					name, i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}

		for (i = 0; i < a->GetNZmax(); i++)
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					name, i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
	}
	else if (a->GetNZ() == -2)
	{
		for (i = 0; i < a->GetRows() + 1; i++)
			if ((*(a->p))[i] != (*(b->p))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->p)) [%d] != (*(b->p)) [%d] => %d != %d\n",
					name, i, i, (*(a->p))[i], (*(b->p))[i]);
				return 0;
			}

		for (i = 0; i < a->GetNZmax(); i++)
			if ((*(a->i))[i] != (*(b->i))[i])
			{
				fprintf(stderr,
					"ERROR: For %s in {a,b} (*(a->i)) [%d] != (*(b->i)) [%d] => %d != %d\n",
					name, i, i, (*(a->i))[i], (*(b->i))[i]);
				return 0;
			}
	}

	for (i = 0; i < a->GetNZmax(); i++)
		if ((*(a->x))[i] != (*(b->x))[i])
		{
			fprintf(stderr,
				"ERROR: For %s in {a, b} (*(a->x)) [%d] != (*(b->x)) [%d] => %g != %g\n",
				name, i, i, (*(a->x))[i], (*(b->x))[i]);
			return 0;
		}

	//if (!compare_matrix_infos((*(a->i))nfo, (*(b->i))nfo))
	//{
	//	fprintf(stderr, "ERROR: matrix %s infos differ\n", name);
	//	return 0;
	//}

	return 1;
}

/* compare two vectors */
static int compare_vectors(std::shared_ptr<std::vector<double>> a, std::shared_ptr<std::vector<double>> b)
{
	bool is_different = true;
	if (a->size() != b->size())
		return 1;

	for (auto row_sel = 0; row_sel<a->size(); ++row_sel)
	{
		if ((*a)[row_sel] != (*b)[row_sel])
			return false;
	}
	return is_different;
}

///* compare problem infos */
//static int compare_infos(struct fclib_info *a, struct fclib_info *b)
//{
//	if (!a && !b) return 1;
//	else if ((a && !b) || (!a && b)) return 0;
//	else if (strcmp(a->title, b->title) != 0 ||
//		strcmp(a->description, b->description) != 0 ||
//		strcmp(a->math_info, b->math_info) != 0) return 0;
//
//	return 1;
//}

/* compare global problems */
static int compare_global_problems(fclib::fclib_global_problem_CPP& a, fclib::fclib_global_problem_CPP& b)
{
	if (!compare_matrices("M", a.matM, b.matM) ||
		!compare_matrices("H", a.matH, b.matH) ||
		!compare_matrices("G", a.matG, b.matG) ||
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
	if (!compare_matrices("W", a.matW, b.matW) ||
		!compare_matrices("V", a.matV, b.matV) ||
		!compare_matrices("R", a.matR, b.matR) ||
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

	bool global = true;
	if (global)
	{
		int space_dimension = 2;
		int global_dofs = 100;
		int contact_points = 100;
		int neq = 100;
		int allfine = 0;
		bool create_G = true;
		auto mat_format = fclib::fclib_matrix_CPP::mat_format_t::COO;
		int num_guesses_out = 3;
		bool create_info = false; //TODO

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

		assert(compare_global_problems(problem_out, problem_in) && "ERROR: written/read problem comparison failed");
		assert(compare_solutions(solution_out, solution_out) && "ERROR: written/read solution comparison failed");
		assert(num_guesses_out == num_guesses_in && "ERROR: numbers of written and read guesses differ");
		for (auto i = 0; i < num_guesses_in; i++)
		{
			assert(compare_solutions(guesses_in[i], guesses_out[i]) && "ERROR: written/read guess comparison failed");
		}

		printf("All comparisons PASSED\n");
	}
	else
	{

		int space_dimension = 2;
		//int local_dofs = 100;
		int contact_points = 100;
		int neq = 100;
		bool create_G = true;
		auto mat_format = fclib::fclib_matrix_CPP::mat_format_t::COO;
		int num_guesses_out = 3;
		bool create_info = false; //TODO

		fclib::fclib_local_problem_CPP problem_out, problem_in;
		fclib::fclib_solution_CPP solution_out, solution_in;
		std::vector<fclib::fclib_solution_CPP> guesses_out, guesses_in;
		guesses_out.resize(num_guesses_out);



		random_local_problem(problem_out, contact_points, neq, space_dimension, mat_format);
		random_local_solutions(problem_out, solution_out);
		random_local_solutions(problem_out, guesses_out, num_guesses_out);


		problem_out.write_problem("output_file.hdf5");
		solution_out.write_solution("output_file.hdf5");
		fclib::fclib_solution_CPP::write_guesses(guesses_out, "output_file.hdf5", num_guesses_out);

		problem_in.read_problem("output_file.hdf5");
		solution_in.read_solution("output_file.hdf5");
		int num_guesses_in = fclib::fclib_solution_CPP::read_guesses(guesses_in, "output_file.hdf5");

		printf("Comparing written and read local problem data ...\n");

		printf("Computing merit function ...\n");

		//double error1 = fclib_merit_local(problem, MERIT_1, solution);
		//double error2 = fclib_merit_local(p, MERIT_1, s);
		//printf("Error for initial problem = %12.8e\n", error1);
		//printf("Error for read problem = %12.8e\n", error2);

		assert(compare_local_problems(problem_out, problem_in) && "ERROR: written/read problem comparison failed");
		assert(compare_solutions(solution_out, solution_out) && "ERROR: written/read solution comparison failed");
		assert(num_guesses_out == num_guesses_in && "ERROR: numbers of written and read guesses differ");
		for (auto i = 0; i < num_guesses_in; i++)
		{
			assert(compare_solutions(guesses_in[i], guesses_out[i]) && "ERROR: written/read guess comparison failed");
		}

		printf("All comparisons PASSED\n");
	}

	remove("output_file.hdf5");

	return 0;
}
