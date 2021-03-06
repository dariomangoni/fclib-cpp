#ifndef _fclib_hpp_
#define _fclib_hpp_

#include <string>
#include <vector>
#include <memory>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <cassert>
#include <fstream>

/* choose api version */
#define H5Gcreate_vers 2
#define H5Gopen_vers 2

/* useful macros */
#define ASSERT(Test, ...)\
  do {\
  if (! (Test)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
    fprintf (stderr, __VA_ARGS__);\
    fprintf (stderr, "\n"); exit (1); } } while (0)

#define IO(Call) ASSERT ((Call) >= 0, "ERROR: HDF5 call failed")
#define MM(Call) ASSERT ((Call), "ERROR: out of memory")


#ifndef FCLIB_APICOMPILE
#define FCLIB_APICOMPILE
#endif

namespace fclib
{

	class FCLIB_APICOMPILE HDF5handler
	{
	private:
		hid_t file_id;
		unsigned flags_used = 0;


	public:
		HDF5handler(std::string path, unsigned flags, bool create_if_not_exist = false) : flags_used(flags)
		{
			bool file_good = false;

			// OPEN statements
			if (flags_used == H5F_ACC_RDONLY ||   /*absence of rdwr => rd-only        */
				flags_used == H5F_ACC_RDWR  /*open for read and write           */
				)    
			{
				bool file_exist;
				{ file_exist = std::ifstream(path).good(); }

				if (file_exist) /* HDF5 outputs lots of warnings when file does not exist */
				{
					file_id = H5Fopen(path.c_str(), flags, H5P_DEFAULT);
					if (file_id < 0)
						fprintf(stderr, "ERROR: opening file failed\n");
					else
						file_good = true;
				}
				else
					if (create_if_not_exist)
						flags_used = H5F_ACC_EXCL;
			}


			// CREATE statements
			if (flags_used == H5F_ACC_EXCL ||    /*fail if file already exists*/
				flags_used == H5F_ACC_TRUNC ||    /*overwrite existing files*/
				flags_used == H5F_ACC_CREAT)       /*create non-existing files*/
			{
				file_id = H5Fcreate(path.c_str(), flags_used, H5P_DEFAULT, H5P_DEFAULT);
				if (file_id < 0)
					fprintf(stderr, "ERROR: creating file failed\n");
				else
					file_good = true;
			}

			assert(file_good);
		}

		hid_t MakeGroup(std::string name) const { return MakeGroup(file_id, name); }

		static hid_t MakeGroup(hid_t loc_id, std::string name)
		{
			if (H5Lexists(loc_id, name.c_str(), H5P_DEFAULT))
			{
				return H5Gopen(loc_id, name.c_str(), H5P_DEFAULT);
			}
			else
				return H5Gcreate(loc_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		}

		~HDF5handler()	{ IO(H5Fclose(file_id)); }

		int AddIntAttribute(std::string parent_path, std::string attr_name, int attr_value) const
		{
			hid_t group_id, dataspace_id, attr_id;
			hsize_t dims = 1;

			assert(flags_used == H5F_ACC_RDWR || flags_used == H5F_ACC_TRUNC);


			if (H5Lexists(file_id, parent_path.c_str(), H5P_DEFAULT))
			{
				IO(group_id = H5Gopen(file_id, parent_path.c_str(), H5P_DEFAULT));
				dataspace_id = H5Screate_simple(1, &dims, nullptr);
				attr_id = H5Acreate(group_id, attr_name.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
				IO(H5Awrite(attr_id, H5T_NATIVE_INT, &attr_value));
				IO(H5Aclose(attr_id));
				IO(H5Gclose(group_id));
			}
			else
			{
				IO(group_id = MakeGroup(file_id, parent_path));
				dataspace_id = H5Screate_simple(1, &dims, nullptr);
				attr_id = H5Acreate(group_id, attr_name.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
				IO(H5Awrite(attr_id, H5T_NATIVE_INT, &attr_value));
				IO(H5Aclose(attr_id));
				IO(H5Gclose(group_id));
			}

			return 1;
		}

		hid_t GetID() const { return file_id; }
	};


	inline bool getStringFromField(hid_t id, std::string& destination, std::string field)
	{
		if (H5LTfind_dataset(id, field.c_str()))
		{
			H5T_class_t class_id;
			hsize_t dim;
			size_t size;

			IO(H5LTget_dataset_info(id, field.c_str(), &dim, &class_id, &size));
			std::vector<char> char_temp(size);
			IO(H5LTread_dataset_string(id, field.c_str(), char_temp.data()));
			destination.assign(char_temp.data());
			return true;
		}

		destination.clear();
		return false;
	}
	
	class FCLIB_APICOMPILE fclib_matrix_CPP
	{
	public:
		std::shared_ptr<std::vector<int>> p = nullptr;
		std::shared_ptr<std::vector<int>> i = nullptr;
		std::shared_ptr<std::vector<double>> x = nullptr;

		int* p_ptr = nullptr;
		int* i_ptr = nullptr;
		double* x_ptr = nullptr;

	protected:

		int nzmax = 0;
		int m = -1;
		int n = -1;
		int nz = 0;

		// optional
		bool has_info = false;
		std::string comment;
		double conditioning = 0.0;
		double determinant = 0.0;
		int rank = 0;


	public:

		enum class mat_format_t
		{
			COO = 0,
			CSC = -1,
			CSR = -2
		};

		fclib_matrix_CPP()
		{
			p = std::make_shared<std::vector<int>>();
			i = std::make_shared<std::vector<int>>();
			x = std::make_shared<std::vector<double>>();
			p_ptr = p->data();
			i_ptr = i->data();
			x_ptr = x->data();
		}

		fclib_matrix_CPP(std::shared_ptr<std::vector<int>> p_in,
			std::shared_ptr<std::vector<int>> i_in,
			std::shared_ptr<std::vector<double>> x_in):
		p(p_in), i(i_in), x(x_in), p_ptr(p_in->data()), i_ptr(i_in->data()), x_ptr(x_in->data())
		{		}

		virtual ~fclib_matrix_CPP(){}

		int GetRows() const { return m; }
		int GetColumns() const { return n; }
		int GetNZmax() const { return nzmax; }
		int GetNZ() const { return nz; }

		void SetComment(std::string comment_in) { has_info = true; comment = comment_in; }
		void SetConditioning(double determinant_in) { has_info = true; determinant = determinant_in; }
		void SetDeterminant(double determinant_in) { has_info = true; determinant = determinant_in; }
		void SetRank(int rank_in) { has_info = true; rank = rank_in; }

		std::string GetComment() const { return comment; }
		double GetConditioning() const { return determinant; }
		double GetDeterminant() const { return determinant; }
		int GetRank() const { return rank; }

		bool HasInfo() const { return has_info; }
		
		void SetElementP(int index, int value) const { p_ptr[index] = value; }
		void SetElementI(int index, int value) const { i_ptr[index] = value; }
		void SetElementX(int index, double value) const { x_ptr[index] = value; }

		int GetElementP(int index) const { return p_ptr[index]; }
		int GetElementI(int index) const { return i_ptr[index]; }
		double GetElementX(int index) const		{ return x_ptr[index]; }

		void Resize(int m_in, int n_in, int nz_in, int nzmax_in = 0)
		{
			assert(p && i && x && "Resizing allowed only for internal memory handling");

			switch (nz_in)
			{
			case -2: // CSR
				p->resize(m_in + 1);
				i->resize(nzmax_in);
				break;
			case -1: // CSC
				p->resize(n_in + 1);
				i->resize(nzmax_in);
				break;
			default:
				if (nz_in < -2)
					throw std::bad_typeid();

				// triplets|COO
				p->resize(nz_in);
				i->resize(nz_in);
				break;
			}

			x->resize(nzmax_in);
			n = n_in;
			m = m_in;
			nz = nz_in;
			nzmax = nzmax_in;

		}

		void read_matrix(hid_t id)
		{
			int nzmax_temp, m_temp, n_temp, nz_temp;
			IO(H5LTread_dataset_int(id, "nzmax", &nzmax_temp));
			IO(H5LTread_dataset_int(id, "m", &m_temp));
			IO(H5LTread_dataset_int(id, "n", &n_temp));
			IO(H5LTread_dataset_int(id, "nz", &nz_temp));

			Resize(m_temp, n_temp, nz_temp, nzmax_temp);

			IO(H5LTread_dataset_int(id, "p", p_ptr ));
			IO(H5LTread_dataset_int(id, "i", i_ptr));
			IO(H5LTread_dataset_double(id, "x", x_ptr ));

			has_info = false;
			if (H5LTfind_dataset(id, "conditioning"))
			{
				has_info = true;
				IO(H5LTread_dataset_double(id, "conditioning", &conditioning));
			}
			if (H5LTfind_dataset(id, "determinant"))
			{
				has_info = true;
				IO(H5LTread_dataset_double(id, "determinant", &determinant));
			}
			if (H5LTfind_dataset(id, "rank"))
			{
				has_info = true;
				IO(H5LTread_dataset_int(id, "rank", &rank));
			}
			if (H5LTfind_dataset(id, "comment"))
			{
				has_info = true;
				getStringFromField(id, comment, std::string("comment"));
			}

		}

		void write_matrix(hid_t id) const
		{
			hsize_t dim = 1;

			IO(H5LTmake_dataset_int(id, "nzmax", 1, &dim, &nzmax));
			IO(H5LTmake_dataset_int(id, "m", 1, &dim, &m));
			IO(H5LTmake_dataset_int(id, "n", 1, &dim, &n));
			IO(H5LTmake_dataset_int(id, "nz", 1, &dim, &nz));

			hsize_t dim_p;
			hsize_t dim_i;
			hsize_t dim_x = nzmax;
			switch (nz)
			{
			case -2: // CSR
				dim_p = m + 1;
				dim_i = nzmax;
				break;
			case -1: // CSC
				dim_p = n + 1;
				dim_i = nzmax;
				break;
			default: // triplets|COO
				dim_p = nz;
				dim_i = nz;
				break;
			}

			IO(H5LTmake_dataset_int(id, "p", 1, &dim_p, p_ptr));
			IO(H5LTmake_dataset_int(id, "i", 1, &dim_i, i_ptr));
			IO(H5LTmake_dataset_double(id, "x", 1, &dim_x, x_ptr));


			if (has_info)
			{
				dim = 1;
				if (comment.c_str()) IO(H5LTmake_dataset_string(id, "comment", comment.c_str()));
				IO(H5LTmake_dataset_double(id, "conditioning", 1, &dim, &conditioning));
				IO(H5LTmake_dataset_double(id, "determinant", 1, &dim, &determinant));
				IO(H5LTmake_dataset_int(id, "rank", 1, &dim, &rank));
			}
		}

	};

	class FCLIB_APICOMPILE fclib_problem_CPP
	{
	protected:
		// optional
		bool has_info = false;
		std::string title;
		std::string description;
		std::string math_info;

		int spacedim = 0;

	public:
		fclib_problem_CPP(){}
		virtual ~fclib_problem_CPP() {}

		void write_problem_info(hid_t id) const
		{
			if (has_info)
			{
				if (!title.empty()) IO(H5LTmake_dataset_string(id, "title", title.c_str()));
				if (!description.empty()) IO(H5LTmake_dataset_string(id, "description", description.c_str()));
				if (!math_info.empty()) IO(H5LTmake_dataset_string(id, "math_info", math_info.c_str()));
			}
		}

		void read_problem_info(hid_t id)
		{
			has_info = getStringFromField(id, title, std::string("title")) | getStringFromField(id, description, std::string("description")) | getStringFromField(id, math_info, std::string("math_info"));
		}

		void SetSpaceDim(int space_dimension) { spacedim = space_dimension; }
		int GetSpaceDim() const { return spacedim; }

		void SetTitle(std::string title_in) { has_info = true; title = title_in; }
		void SetDescription(std::string description_in) { has_info = true; description = description_in; }
		void SetMathInfo(std::string math_info_in) { has_info = true; math_info = math_info_in; }

		std::string GetTitle() const {return title;}
		std::string GetDescription() const {return description;}
		std::string GetMathInfo() const {return math_info;}

		bool HasInfo() const { return has_info; }

		virtual void read_problem(std::string path) = 0;
		virtual void write_problem(std::string path) const = 0;
		virtual void read_vectors(hid_t group_id) const = 0;
		virtual void write_vectors(hid_t group_id) const = 0;

	};

	class FCLIB_APICOMPILE fclib_local_problem_CPP : public fclib_problem_CPP
	{
	public:
		std::shared_ptr<fclib_matrix_CPP> matW, matV, matR;
		std::shared_ptr<std::vector<double>> mu, q, s;

		double* mu_ptr = nullptr;
		double* q_ptr = nullptr;
		double* s_ptr = nullptr;

	protected:

		void initializeVectors()
		{
			mu = std::make_shared<std::vector<double>>();
			q = std::make_shared<std::vector<double>>();
			s = std::make_shared<std::vector<double>>();
			mu_ptr = mu->data();
			q_ptr = q->data();
			s_ptr = s->data();
		}

		void initializeMatrices()
		{
			matW = std::make_shared<fclib_matrix_CPP>();
			matV = std::make_shared<fclib_matrix_CPP>();
			matR = std::make_shared<fclib_matrix_CPP>();
		}

	public:
		virtual ~fclib_local_problem_CPP() {}

		fclib_local_problem_CPP() {
			initializeVectors();
			initializeMatrices();
		}

		fclib_local_problem_CPP(std::shared_ptr<std::vector<double>> mu_in,
			std::shared_ptr<std::vector<double>> q_in,
			std::shared_ptr<std::vector<double>> s_in)
			: mu(mu_in), q(q_in), s(s_in)
		{
			initializeMatrices();
		}

		fclib_local_problem_CPP(std::shared_ptr<fclib_matrix_CPP> matW_in,
			std::shared_ptr<fclib_matrix_CPP> matV_in,
			std::shared_ptr<fclib_matrix_CPP> matR_in)
			: matW(matW_in), matV(matV_in), matR(matR_in)
		{
			initializeVectors();
		}

		fclib_local_problem_CPP(std::shared_ptr<std::vector<double>> mu_in,
			std::shared_ptr<std::vector<double>> q_in,
			std::shared_ptr<std::vector<double>> s_in,
			std::shared_ptr<fclib_matrix_CPP> matW_in,
			std::shared_ptr<fclib_matrix_CPP> matV_in,
			std::shared_ptr<fclib_matrix_CPP> matR_in)
			: matW(matW_in), matV(matV_in), matR(matR_in), mu(mu_in), q(q_in), s(s_in)
		{
		}


		void read_vectors(hid_t group_id) const override
		{

			q->resize(matW->GetRows());
			IO(H5LTread_dataset_double(group_id, "q", q->data()));


			assert(matW->GetRows() % spacedim == 0 && "ERROR: number of W rows is not divisble by the spatial dimension");
			mu->resize(matW->GetRows() / spacedim);
			IO(H5LTread_dataset_double(group_id, "mu", mu->data()));


			if (matR && matR->GetRows()>=0)
			{
				s->resize(matR->GetRows());
				IO(H5LTread_dataset_double(group_id, "s", s->data()));
			}

		}

		void write_vectors(hid_t group_id) const override
		{
			hsize_t dim;

			dim = matW->GetRows();
			assert(q->data() && "ERROR: q must be given");
			IO(H5LTmake_dataset_double(group_id, "q", 1, &dim, q->data()));

			ASSERT(dim % spacedim == 0, "ERROR: number of W rows is not divisble by the spatial dimension");
			dim = matW->GetRows() / spacedim;
			IO(H5LTmake_dataset_double(group_id, "mu", 1, &dim, mu->data()));

			if (matV && matV->GetRows() >= 0)
			{
				dim = matR->GetRows();
				assert(s->data() && "ERROR: s must be given if R is present");
				IO(H5LTmake_dataset_double(group_id, "s", 1, &dim, s->data()));
			}
		}

		void read_problem(std::string path) override
		{
			HDF5handler file_handler(path, H5F_ACC_RDONLY);
			hid_t  file_id, main_id, id;
			file_id = file_handler.GetID();

			if (!H5Lexists(file_id, "/fclib_local", H5P_DEFAULT))
			{
				fprintf(stderr, "ERROR: spurious input file :: fclib_local group does not exists");
			}

			IO(main_id = H5Gopen(file_id, "/fclib_local", H5P_DEFAULT));
			IO(H5LTread_dataset_int(file_id, "/fclib_local/spacedim", &spacedim));

			IO(id = H5Gopen(file_id, "/fclib_local/W", H5P_DEFAULT));
			matW->read_matrix(id);
			IO(H5Gclose(id));



			if (H5Lexists(file_id, "/fclib_local/V", H5P_DEFAULT))
			{
				IO(id = H5Gopen(file_id, "/fclib_local/V", H5P_DEFAULT));
				matV->read_matrix(id);
				IO(H5Gclose(id));

				IO(id = H5Gopen(file_id, "/fclib_local/R", H5P_DEFAULT));
				matR->read_matrix(id);
				IO(H5Gclose(id));
			}


			IO(id = H5Gopen(file_id, "/fclib_local/vectors", H5P_DEFAULT));
			read_vectors(id);
			IO(H5Gclose(id));


			if (H5Lexists(file_id, "/fclib_local/info", H5P_DEFAULT))
			{
				IO(id = H5Gopen(file_id, "/fclib_local/info", H5P_DEFAULT));
				read_problem_info(id);
				IO(H5Gclose(id));
			}

			IO(H5Gclose(main_id));
		}

		void write_problem(std::string path) const override
		{
			HDF5handler file_handler(path, H5F_ACC_TRUNC);
			hid_t  file_id, parent_group_id, group_id;
			hsize_t dim = 1;
			file_id = file_handler.GetID();

			IO(parent_group_id = HDF5handler::MakeGroup(file_id, "/fclib_local"));

			assert((spacedim == 2 || spacedim == 3) && "ERROR: space dimension must be 2 or 3");
			IO(H5LTmake_dataset_int(file_id, "/fclib_local/spacedim", 1, &dim, &spacedim));

			assert(matW->GetRows()>=0 && "ERROR: W must be given");
			IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_local/W"));
			matW->write_matrix(group_id);
			IO(H5Gclose(group_id));

			if (matV && matR && matV->GetRows() >= 0 && matR->GetRows() >= 0)
			{
				IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_local/V"));
				matV->write_matrix(group_id);
				IO(H5Gclose(group_id));

				IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_local/R"));
				matR->write_matrix(group_id);
				IO(H5Gclose(group_id));
			}
			else 
				assert(((matR.get() == nullptr || matR->GetRows() >= 0) == (matV.get() == nullptr || matV->GetRows() >= 0 )) && "ERROR: V and R must be defined at the same time");

			IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_local/vectors"));
			write_vectors(group_id);
			IO(H5Gclose(group_id));

			if (has_info)
			{
				IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_local/info"));
				write_problem_info(group_id);
				IO(H5Gclose(group_id));
			}

			IO(H5Gclose(parent_group_id));
		}
		

	};


	class FCLIB_APICOMPILE fclib_global_problem_CPP : public fclib_problem_CPP
	{
	public:
		std::shared_ptr<fclib_matrix_CPP> matM, matH, matG;

		std::shared_ptr<std::vector<double>> mu, f, b, w;
	protected:
		//TODO
		int m = 0;
		int p = 0;
		int nc = 0;

		void initializeVectors()
		{
			mu = std::make_shared<std::vector<double>>();
			f = std::make_shared<std::vector<double>>();
			b = std::make_shared<std::vector<double>>();
			w = std::make_shared<std::vector<double>>();
		}

		void initializeMatrices()
		{
			matM = std::make_shared<fclib_matrix_CPP>();
			matH = std::make_shared<fclib_matrix_CPP>();
			matG = std::make_shared<fclib_matrix_CPP>();
		}

	public:
		virtual ~fclib_global_problem_CPP(){}

		fclib_global_problem_CPP(){
			initializeVectors();
			initializeMatrices();
		}

		fclib_global_problem_CPP(std::shared_ptr<std::vector<double>> mu_in,
								std::shared_ptr<std::vector<double>> f_in,
								std::shared_ptr<std::vector<double>> b_in,
								std::shared_ptr<std::vector<double>> w_in)
		: mu(mu_in), f(f_in), b(b_in), w(w_in)
		{
			initializeMatrices();
		}

		fclib_global_problem_CPP(std::shared_ptr<fclib_matrix_CPP> matM_in,
								std::shared_ptr<fclib_matrix_CPP> matH_in,
								std::shared_ptr<fclib_matrix_CPP> matG_in)
		: matM(matM_in), matH(matH_in), matG(matG_in)
		{
			initializeVectors();
		}

		fclib_global_problem_CPP(std::shared_ptr<std::vector<double>> mu_in,
			                    std::shared_ptr<std::vector<double>> f_in,
			                    std::shared_ptr<std::vector<double>> b_in,
			                    std::shared_ptr<std::vector<double>> w_in,
			                    std::shared_ptr<fclib_matrix_CPP> matM_in,
			                    std::shared_ptr<fclib_matrix_CPP> matH_in,
			                    std::shared_ptr<fclib_matrix_CPP> matG_in)
		: matM(matM_in), matH(matH_in), matG(matG_in), mu(mu_in), f(f_in), b(b_in), w(w_in)
		{
		}


		void read_problem(std::string path) override
		{
			HDF5handler file_handler(path, H5F_ACC_RDONLY);
			hid_t  file_id, main_id, id;
			file_id = file_handler.GetID();

			if (!H5Lexists(file_id, "/fclib_global", H5P_DEFAULT))
			{
				fprintf(stderr, "ERROR: spurious input file :: fclib_global group does not exists");
			}

			IO(main_id = H5Gopen(file_id, "/fclib_global", H5P_DEFAULT));
			IO(H5LTread_dataset_int(file_id, "/fclib_global/spacedim", &spacedim));

			IO(id = H5Gopen(file_id, "/fclib_global/M", H5P_DEFAULT));
			matM->read_matrix(id);
			IO(H5Gclose(id));

			IO(id = H5Gopen(file_id, "/fclib_global/H", H5P_DEFAULT));
			matH->read_matrix(id);
			IO(H5Gclose(id));

			if (H5Lexists(file_id, "/fclib_global/G", H5P_DEFAULT))
			{
				IO(id = H5Gopen(file_id, "/fclib_global/G", H5P_DEFAULT));
				matG->read_matrix(id);
				IO(H5Gclose(id));
			}

			IO(id = H5Gopen(file_id, "/fclib_global/vectors", H5P_DEFAULT));
			read_vectors(id);
			IO(H5Gclose(id));

			if (H5Lexists(file_id, "/fclib_global/info", H5P_DEFAULT))
			{
				IO(id = H5Gopen(file_id, "/fclib_global/info", H5P_DEFAULT));
				read_problem_info(id);
				IO(H5Gclose(id));
			}

			IO(H5Gclose(main_id));
		}

		void write_problem(std::string path) const override
		{
			HDF5handler file_handler(path, H5F_ACC_TRUNC);
			hid_t file_id, group_id, dataset_id;
			hsize_t dim = 1;

			file_id = file_handler.GetID();

			IO(group_id = HDF5handler::MakeGroup(file_id, "/fclib_global"));

			assert((spacedim == 2 || spacedim == 3) && "ERROR: space dimension must be 2 or 3");
			IO(H5LTmake_dataset_int(file_id, "/fclib_global/spacedim", 1, &dim, &spacedim));

			assert(matM && "ERROR: M must be given");
			IO(dataset_id = HDF5handler::MakeGroup(file_id, "/fclib_global/M"));
			matM->write_matrix(dataset_id);
			IO(H5Gclose(dataset_id));

			assert(matH && "ERROR: H must be given");
			IO(dataset_id = HDF5handler::MakeGroup(file_id, "/fclib_global/H"));
			matH->write_matrix(dataset_id);
			IO(H5Gclose(dataset_id));

			if (matG && matG->GetRows() >= 0)
			{
				IO(dataset_id = HDF5handler::MakeGroup(file_id, "/fclib_global/G"));
				matG->write_matrix(dataset_id);
				IO(H5Gclose(dataset_id));
			}

			IO(dataset_id = HDF5handler::MakeGroup(file_id, "/fclib_global/vectors"));
			write_vectors(dataset_id);
			IO(H5Gclose(dataset_id));

			if (has_info)
			{
				IO(dataset_id = HDF5handler::MakeGroup(file_id, "/fclib_global/info"));
				write_problem_info(dataset_id);
				IO(H5Gclose(dataset_id));
			}

			IO(H5Gclose(group_id));
		}

		void read_vectors(hid_t group_id) const override
		{
			f->resize(matM->GetRows());
			IO(H5LTread_dataset_double(group_id, "f", f->data()));

			assert(matH->GetColumns() % spacedim == 0 && "ERROR: number of H columns is not divisble by the spatial dimension");
			w->resize(matH->GetColumns());
		    mu->resize(matH->GetColumns()/spacedim);

			IO(H5LTread_dataset_double(group_id, "w", w->data()));
			IO(H5LTread_dataset_double(group_id, "mu", mu->data()));

			if (matG && matG->GetRows() >= 0)
			{
				b->resize(matG->GetColumns());
				IO(H5LTread_dataset_double(group_id, "b", b->data()));
			}
		}

		void write_vectors(hid_t group_id) const override
		{
			hsize_t dim;

			dim = matM->GetRows();
			assert(f->data() && "ERROR: f must be given");
			IO(H5LTmake_dataset_double(group_id, "f", 1, &dim, f->data()));

			dim = matH->GetColumns();
			assert(w->data() && mu->data() && "ERROR: w and mu must be given");
			IO(H5LTmake_dataset_double(group_id, "w", 1, &dim, w->data()));
			assert(dim % spacedim == 0 && "ERROR: number of H columns is not divisble by the spatial dimension");
			dim = matH->GetColumns() / spacedim;
			IO(H5LTmake_dataset_double(group_id, "mu", 1, &dim, mu->data()));

			if (matG && matG->GetRows() >= 0)
			{
				dim = matG->GetColumns();
				assert(b->data() && "ERROR: b must be given if G is present");
				IO(H5LTmake_dataset_double(group_id, "b", 1, &dim, b->data()));
			}

		}

		

		};



		class FCLIB_APICOMPILE fclib_solution_CPP
		{
		public:
			std::shared_ptr<std::vector<double>> v, u, r, l;

		public:
			explicit fclib_solution_CPP()
			{
				v = std::make_shared<std::vector<double>>();
				u = std::make_shared<std::vector<double>>();
				r = std::make_shared<std::vector<double>>();
				l = std::make_shared<std::vector<double>>();
			}

			fclib_solution_CPP(std::shared_ptr<std::vector<double>> v_in,
								std::shared_ptr<std::vector<double>> u_in,
								std::shared_ptr<std::vector<double>> r_in,
								std::shared_ptr<std::vector<double>> l_in)
			: v(v_in), u(u_in), r(r_in), l(l_in) {}

			~fclib_solution_CPP(){}

			void write_solution(std::string path) const
			{
				HDF5handler file_handler(path, H5F_ACC_RDWR);
				
				auto file_id = file_handler.GetID();

				hid_t  group_id;
				int nv, nr, nl;


				if (H5Lexists(file_id, "/solution", H5P_DEFAULT)) /* cannot overwrite existing datasets */
				{
					fprintf(stderr, "ERROR: a solution has already been written to this file\n");
					assert(0);
				}


				if (!read_solution_sizes(file_id, nv, nr, nl)) assert(0);

				IO(group_id = HDF5handler::MakeGroup(file_id, "/solution"));
				write_solution(group_id);
				IO(H5Gclose(group_id));

			}

			void write_solution(hid_t file_id) const
			{
				if (l && l->size()>0)
				{
					auto nl = static_cast<hsize_t>(l->size());
					if (nl) IO(H5LTmake_dataset_double(file_id, "l", 1, &nl, l->data()));
				}

				if (v && v->size()>0)
				{
					auto nv = static_cast<hsize_t>(v->size());
					IO(H5LTmake_dataset_double(file_id, "v", 1, &nv, v->data()));
				}

				auto nr = static_cast<hsize_t>(r->size());
				assert(nr && "ERROR: contact constraints must be present");
				IO(H5LTmake_dataset_double(file_id, "u", 1, &nr, u->data()));
				IO(H5LTmake_dataset_double(file_id, "r", 1, &nr, r->data()));
			}

			void read_solution(std::string path)
			{
				HDF5handler file_handler(path, H5F_ACC_RDONLY);
				auto file_id = file_handler.GetID();

				hid_t  group_id;
				int nv, nr, nl;

				if (!read_solution_sizes(file_id, nv, nr, nl)) assert(0);

				IO(group_id = H5Gopen(file_id, "/solution", H5P_DEFAULT));
				read_solution(group_id, nv, nr, nl);
				IO(H5Gclose(group_id));

			}


			void read_solution(hid_t file_id, hsize_t nv, hsize_t nr, hsize_t nl)
			{

				if (nv)
				{
					v->resize(nv);
					IO(H5LTread_dataset_double(file_id, "v", v->data()));
				}
				else
					v.reset();

				if (nl)
				{
					l->resize(nl);
					IO(H5LTread_dataset_double(file_id, "l", l->data()));
				}
				else
					l.reset();

				assert(nr && "ERROR: contact constraints must be present");
				u->resize(nr);
				IO(H5LTread_dataset_double(file_id, "u", u->data()));

				r->resize(nr);
				IO(H5LTread_dataset_double(file_id, "r", r->data()));
			}


			static int read_solution_sizes(hid_t file_id, int& nv, int& nr, int& nl) //TODO: nv, nl, nr may should be data members?
			{

				if (H5Lexists(file_id, "/fclib_global", H5P_DEFAULT))
				{
					IO(H5LTread_dataset_int(file_id, "/fclib_global/M/n", &nv));
					IO(H5LTread_dataset_int(file_id, "/fclib_global/H/n", &nr));
					if (H5Lexists(file_id, "/fclib_global/G", H5P_DEFAULT))
					{
						IO(H5LTread_dataset_int(file_id, "/fclib_global/G/n", &nl));
					}
					else nl = 0;
				}
				else if (H5Lexists(file_id, "/fclib_local", H5P_DEFAULT))
				{
					nv = 0;
					IO(H5LTread_dataset_int(file_id, "/fclib_local/W/n", &nr));
					if (H5Lexists(file_id, "/fclib_local/R", H5P_DEFAULT))
					{
						IO(H5LTread_dataset_int(file_id, "/fclib_local/R/n", &nl));
					}
					else nl = 0;
				}
				else
				{
					fprintf(stderr, "ERROR: neither global nor local problem has been stored. Global or local have to be stored before solutions or guesses\n");
					return 0;
				}

				return 1;
			}

			static int write_guesses(std::vector<fclib_solution_CPP>& guesses, std::string path, int number_of_guesses)
			{
				hid_t  parent_group_id, group_id;
				int nv, nr, nl;
				hsize_t dim = 1;
				char num[128];
				HDF5handler file_handler(path, H5F_ACC_RDWR, true);
				auto file_id = file_handler.GetID();


				if (H5Lexists(file_id, "/guesses", H5P_DEFAULT)) /* cannot overwrite existing datasets */
				{
					fprintf(stderr, "ERROR: some guesses have already been written to this file\n");
					return 0;
				}


				if (!read_solution_sizes(file_id, nv, nr, nl)) return 0;

				IO(parent_group_id = HDF5handler::MakeGroup(file_id, "/guesses"));
				IO(H5LTmake_dataset_int(file_id, "/guesses/number_of_guesses", 1, &dim, &number_of_guesses));

				for (auto i = 0; i < number_of_guesses; i++)
				{
					snprintf(num, 128, "%d", i + 1);
					IO(group_id = HDF5handler::MakeGroup(parent_group_id, num));
					guesses[i].write_solution(group_id);
					IO(H5Gclose(group_id));
				}

				IO(H5Gclose(parent_group_id));

				return 1;
			}

			static int read_guesses(std::vector<fclib_solution_CPP>& guesses, std::string path)
			{
				hid_t  main_id, group_id;
				int nv, nr, nl;
				char num[128];
				int number_of_guesses;

				HDF5handler file_handler(path, H5F_ACC_RDWR);
				auto file_id = file_handler.GetID();

				if (!read_solution_sizes(file_id, nv, nr, nl)) return 0;

				if (H5Lexists(file_id, "/guesses", H5P_DEFAULT))
				{
					IO(main_id = H5Gopen(file_id, "/guesses", H5P_DEFAULT));

					IO(H5LTread_dataset_int(file_id, "/guesses/number_of_guesses", &number_of_guesses));

					guesses.resize(number_of_guesses);
					for (auto i = 0; i < guesses.size(); ++i)
					{
						snprintf(num, 128, "%d", i + 1);
						IO(group_id = H5Gopen(main_id, num, H5P_DEFAULT));
						guesses[i].read_solution(group_id, nv, nr, nl);
						IO(H5Gclose(group_id));
					}

					IO(H5Gclose(main_id));
				}

				return number_of_guesses;

			}
		};



	inline int fclib_create_int_attributes_in_info(std::string path, std::string attr_name, int attr_value)
	{
		HDF5handler file_handler(path, H5F_ACC_RDWR);
		return file_handler.AddIntAttribute(std::string("/fclib_local/info"), attr_name, attr_value);
	}


} //end namespace fclib




#endif