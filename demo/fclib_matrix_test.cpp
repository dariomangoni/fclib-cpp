#include "fclib_class.hpp"

int main()
{

	//std::string path = "C:/workspace/fclib-library/Capsules/Capsules-i94-48.hdf5";

	fclib::HDF5handler file_handler(std::string("C:/workspace/fclib-library/Capsules/Capsules-i94-48.hdf5"), H5F_ACC_RDONLY);

	hid_t  file_id, main_id, id;

	file_id = file_handler.GetID();

	fclib::fclib_matrix_CPP mat;

	IO(main_id = H5Gopen(file_id, "/fclib_local", H5P_DEFAULT));
	//IO(H5LTread_dataset_int(file_id, "/fclib_global/spacedim", &problem->spacedim));

	IO(id = H5Gopen(file_id, "/fclib_local/W", H5P_DEFAULT));
	mat.read_matrix(id);
	IO(H5Gclose(id));

	IO(H5Gclose(main_id));

	mat.DisplayMatrix();
	getchar();

}