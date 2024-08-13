#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <cstring>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstdint>

extern "C" {
void to_c_string(const char* fortran_str, char* c_str, int fortran_str_len) {
    int length = fortran_str_len;
    while (length > 0 && fortran_str[length - 1] == ' ') {
        --length; // Decrement length to ignore trailing spaces
    }
    strncpy(c_str, fortran_str, length);
    c_str[length] = '\0'; // Null-terminate the C string
}
// A helper function to check the return values and print errors
void check_hdf5_error(herr_t status, const char* message) {
    if (status < 0) {
        std::cerr << "HDF5 error: " << message << std::endl;
        exit(EXIT_FAILURE);
    }
}
// Function to read a 1D dataset from HDF5
void read_dataset_1D(hid_t file_id, const char* name, void* buffer, const hsize_t* dims, size_t elem_size_bytes, const hsize_t* chunk_dims) {

    //std::cerr << "Reading 1D dataset: " << name << std::endl;

    // Open dataset
    hid_t dataset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "Error opening dataset" << std::endl;
        return;
    }

    // Get dataspace
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        std::cerr << "Error getting dataspace for dataset" << std::endl;
        H5Dclose(dataset_id);
        return;
    }

    // Verify dimensions match
    hsize_t actual_dims[1] = {0};
    int rank = H5Sget_simple_extent_dims(dataspace_id, actual_dims, NULL);
    if (rank != 1) {
        std::cerr << "Dimension mismatch: expected 1, got " << rank << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return;
    }

    //std::cerr << "Expected dimensions: ";
    //std::cerr << dims[0] << std::endl;

    //std::cerr << "Actual dimensions: ";
    //std::cerr << actual_dims[0] << std::endl;

    // Initialize total bytes read
    size_t total_bytes_read = 0;

    // Read data in chunks
    for (hsize_t i = 0; i < dims[0]; i += chunk_dims[0]) {
        hsize_t start[1] = {i};
        hsize_t count[1] = {
            std::min(chunk_dims[0], dims[0] - i)
        };

        //std::cerr << "Selecting hyperslab: start=(" << start[0] << ") count=(" << count[0] << ")" << std::endl;
        herr_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
        if (status < 0) {
            std::cerr << "Error selecting hyperslab" << std::endl;
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return;
        }

        // Create memory dataspace for chunk
        hid_t memspace_id = H5Screate_simple(1, count, NULL);
        if (memspace_id < 0) {
            std::cerr << "Error creating memory dataspace" << std::endl;
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return;
        }

        // Allocate local buffer
        size_t buf_size = count[0] * elem_size_bytes;
        void* local_buffer = malloc(buf_size);
        if (local_buffer == NULL) {
            std::cerr << "Error allocating memory for local buffer" << std::endl;
            H5Sclose(memspace_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return;
        }

        total_bytes_read += buf_size;
        //std::cerr << "Reading chunk: start=(" << start[0] << ") size=" << buf_size << " bytes, total bytes read=" << total_bytes_read << " bytes" << std::endl;

        // Read data
        if (elem_size_bytes == sizeof(double)) {
            status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
        } else if (elem_size_bytes == sizeof(int)) {
            status = H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
        } else {
            std::cerr << "Unsupported element size" << std::endl;
            free(local_buffer);
            H5Sclose(memspace_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return;
        }

        if (status < 0) {
            std::cerr << "Error reading data from dataset" << std::endl;
            free(local_buffer);
            H5Sclose(memspace_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return;
        }

        // Copy data to the main buffer
        for (hsize_t x = 0; x < count[0]; ++x) {
            size_t buf_index = i + x;
            size_t local_index = x;
            if (elem_size_bytes == sizeof(double)) {
                ((double*)buffer)[buf_index] = ((double*)local_buffer)[local_index];
            } else if (elem_size_bytes == sizeof(int)) {
                ((int*)buffer)[buf_index] = ((int*)local_buffer)[local_index];
            }
        }

        free(local_buffer);
        H5Sclose(memspace_id);
    }

    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    //std::cerr << "Finished reading 1D dataset: " << name << std::endl;
}
// Function to read a 2D dataset from HDF5
void read_dataset_2D(hid_t file_id, const char* name, void* buffer, const hsize_t* dims, size_t elem_size_bytes, const hsize_t* chunk_dims) {

    //std::cerr << "Reading 2D dataset: " << name << std::endl;

    // Open dataset
    hid_t dataset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "Error opening dataset" << std::endl;
        return;
    }

    // Get dataspace
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        std::cerr << "Error getting dataspace for dataset" << std::endl;
        H5Dclose(dataset_id);
        return;
    }

    // Verify dimensions match
    hsize_t actual_dims[2] = {0};
    int rank = H5Sget_simple_extent_dims(dataspace_id, actual_dims, NULL);
    if (rank != 2) {
        std::cerr << "Dimension mismatch: expected 2, got " << rank << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return;
    }

    //std::cerr << "Expected dimensions: ";
    //for (size_t i = 0; i < 2; ++i) {
    //    std::cerr << dims[i] << " ";
    //}
    //std::cerr << std::endl;

    //std::cerr << "Actual dimensions: ";
    //for (size_t i = 0; i < 2; ++i) {
    //    std::cerr << actual_dims[i] << " ";
    //}
    //std::cerr << std::endl;

    // Initialize total bytes read
    size_t total_bytes_read = 0;

    // Read data in chunks
    for (hsize_t i = 0; i < dims[0]; i += chunk_dims[0]) {
        for (hsize_t j = 0; j < dims[1]; j += chunk_dims[1]) {
            hsize_t start[2] = {i, j};
            hsize_t count[2] = {
                std::min(chunk_dims[0], dims[0] - i),
                std::min(chunk_dims[1], dims[1] - j)
            };

            //std::cerr << "Selecting hyperslab: start=(" << start[0] << "," << start[1] << ") count=(" << count[0] << "," << count[1] << ")" << std::endl;
            herr_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
            if (status < 0) {
                std::cerr << "Error selecting hyperslab" << std::endl;
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return;
            }

            // Create memory dataspace for chunk
            hid_t memspace_id = H5Screate_simple(2, count, NULL);
            if (memspace_id < 0) {
                std::cerr << "Error creating memory dataspace" << std::endl;
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return;
            }

            // Allocate local buffer
            size_t buf_size = count[0] * count[1] * elem_size_bytes;
            void* local_buffer = malloc(buf_size);
            if (local_buffer == NULL) {
                std::cerr << "Error allocating memory for local buffer" << std::endl;
                H5Sclose(memspace_id);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return;
            }

            total_bytes_read += buf_size;
            //std::cerr << "Reading chunk: start=(" << start[0] << "," << start[1] << ") size=" << buf_size << " bytes, total bytes read=" << total_bytes_read << " bytes" << std::endl;

            // Read data
            if (elem_size_bytes == sizeof(double)) {
                status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
            } else if (elem_size_bytes == sizeof(int)) {
                status = H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
            } else {
                std::cerr << "Unsupported element size" << std::endl;
                free(local_buffer);
                H5Sclose(memspace_id);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return;
            }

            if (status < 0) {
                std::cerr << "Error reading data from dataset" << std::endl;
                free(local_buffer);
                H5Sclose(memspace_id);
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return;
            }

            // Copy data to the main buffer
            for (hsize_t x = 0; x < count[0]; ++x) {
                for (hsize_t y = 0; y < count[1]; ++y) {
                    size_t buf_index = (i + x) * dims[1] + (j + y);
                    size_t local_index = x * count[1] + y;
                    if (elem_size_bytes == sizeof(double)) {
                        ((double*)buffer)[buf_index] = ((double*)local_buffer)[local_index];
                    } else if (elem_size_bytes == sizeof(int)) {
                        ((int*)buffer)[buf_index] = ((int*)local_buffer)[local_index];
                    }
                }
            }

            free(local_buffer);
            H5Sclose(memspace_id);
        }
    }

    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    //std::cerr << "Finished reading 2D dataset: " << name << std::endl;
}
void read_dataset_3D(hid_t file_id, const char* name, void* buffer, const hsize_t* dims, size_t elem_size_bytes, const hsize_t* chunk_dims) {
    // Open dataset
    hid_t dataset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dataset_id < 0) {
        //std::cerr << "Error opening dataset" << std::endl;
        return;
    }

    // Get dataspace
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        //std::cerr << "Error getting dataspace for dataset" << std::endl;
        H5Dclose(dataset_id);
        return;
    }

    // Verify dimensions match
    hsize_t actual_dims[3] = {0};
    int rank = H5Sget_simple_extent_dims(dataspace_id, actual_dims, NULL);
    if (rank != 3) {
        //std::cerr << "Dimension mismatch: expected 3, got " << rank << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return;
    }

    // Initialize total bytes read
    size_t total_bytes_read = 0;

    // Read data in chunks
    for (hsize_t i = 0; i < dims[0]; i += chunk_dims[0]) {
        for (hsize_t j = 0; j < dims[1]; j += chunk_dims[1]) {
            for (hsize_t k = 0; k < dims[2]; k += chunk_dims[2]) {
                hsize_t start[3] = {i, j, k};
                hsize_t count[3] = {
                    std::min(chunk_dims[0], dims[0] - i),
                    std::min(chunk_dims[1], dims[1] - j),
                    std::min(chunk_dims[2], dims[2] - k)
                };

                herr_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
                if (status < 0) {
                    //std::cerr << "Error selecting hyperslab" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return;
                }

                // Create memory dataspace for chunk
                hid_t memspace_id = H5Screate_simple(3, count, NULL);
                if (memspace_id < 0) {
                    //std::cerr << "Error creating memory dataspace" << std::endl;
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return;
                }

                // Allocate local buffer
                size_t buf_size = count[0] * count[1] * count[2] * elem_size_bytes;
                void* local_buffer = malloc(buf_size);
                if (local_buffer == NULL) {
                    //std::cerr << "Error allocating memory for local buffer" << std::endl;
                    H5Sclose(memspace_id);
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return;
                }

                total_bytes_read += buf_size;

                //std::cerr << "Reading chunk: start=(" << start[0] << "," << start[1] << "," << start[2] << ") size=" << buf_size << " bytes, total bytes read=" << total_bytes_read << " bytes" << std::endl;

                // Read data
                if (elem_size_bytes == sizeof(double)) {
                    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
                } else if (elem_size_bytes == sizeof(int)) {
                    status = H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, local_buffer);
                } else {
                    std::cerr << "Unsupported element size" << std::endl;
                    free(local_buffer);
                    H5Sclose(memspace_id);
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return;
                }

                if (status < 0) {
                    std::cerr << "Error reading data from dataset" << std::endl;
                    free(local_buffer);
                    H5Sclose(memspace_id);
                    H5Sclose(dataspace_id);
                    H5Dclose(dataset_id);
                    return;
                }

                //std::cout << "Copying data to main buffer" <<  std::endl;

                // Copy data to the main buffer
                for (hsize_t x = 0; x < count[0]; ++x) {
                    for (hsize_t y = 0; y < count[1]; ++y) {
                        for (hsize_t z = 0; z < count[2]; ++z) {
                            size_t buf_index = (i + x) * dims[1] * dims[2] + (j + y) * dims[2] + (k + z);
                            size_t local_index = x * count[1] * count[2] + y * count[2] + z;

                            // Ensure index is within bounds
                            if (buf_index >= dims[0] * dims[1] * dims[2] || local_index >= count[0] * count[1] * count[2]) {
                                std::cerr << "Index out of bounds: buf_index=" << buf_index << ", local_index=" << local_index << std::endl;
                                continue;
                            }

                            // Print out data type size and value being copied
                            if (elem_size_bytes == sizeof(double)) {
                                if (buf_index == 0) {
                                    std::cout << "Copying data from local_buffer to buffer:" << std::endl;
                                    std::cout << "Data type size (double): " << sizeof(double) << " bytes" << std::endl;
                                    std::cout << "Index in buffer: " << buf_index << std::endl;
                                    std::cout << "Index in local_buffer: " << local_index << std::endl;

                                    std::cout << "Buffer address: " << buffer << std::endl;
                                    std::cout << "Buffer size (in bytes): " << dims[0] * dims[1] * dims[2] * elem_size_bytes << std::endl;
                                    std::cout << "buf_index: " << buf_index << std::endl;

                                    // Ensure that buf_index is within bounds
                                    if (buf_index >= 0 && buf_index < (dims[0] * dims[1] * dims[2])) {
                                        std::cout << "Destination buffer index before copy: " << ((double*)buffer)[buf_index] << std::endl;
                                    } else {
                                        std::cerr << "Index out of bounds: " << buf_index << std::endl;
                                    }

                                    std::cout << "Value to copy: " << ((double*)local_buffer)[local_index] << std::endl;
                                    std::cout << "Destination buffer index before copy: " << ((double*)buffer)[buf_index] << std::endl;

                                    ((double*)buffer)[buf_index] = ((double*)local_buffer)[local_index];
                                    std::cout << "Destination buffer index after copy: " << ((double*)buffer)[buf_index] << std::endl;
                                }
                            ((double*)buffer)[buf_index] = ((double*)local_buffer)[local_index];
                            } else if (elem_size_bytes == sizeof(int)) {
                                ((int*)buffer)[buf_index] = ((int*)local_buffer)[local_index];
                            }
                        }
                    }
                }

                free(local_buffer);
                H5Sclose(memspace_id);
            }
        }
    }

    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}
void read_hdf5_data(const char* filename, void* iz_f_ptr, void* ite_ptr, void* jne_ptr,
                        void* epatom_ptr, void* amamu_f_ptr, void* eumesh_ptr, void* sig_ptr, int* ierr) {
 
    
    std::cout << "Address of sig in C++: " << sig_ptr << std::endl;
    
	//uintptr_t address = reinterpret_cast<uintptr_t>(&sig[0]);
	//std::cout << "Hexadecimal Address of sig: 0x" << std::hex << address << std::endl;



    //std::cerr << "Loading OP mono data from HDF5 file..." << std::endl;

    int nel = 17;
    int np_mesh = 1648;
    int nptot = 10000;
    
    size_t expected_size_3D = static_cast<size_t>(nel) * static_cast<size_t>(np_mesh) * static_cast<size_t>(nptot);
    std::cerr << "Expected size of sig vector: " << expected_size_3D << std::endl;
    size_t expected_size_2D = static_cast<size_t>(nel) * static_cast<size_t>(np_mesh);
    std::cerr << "Expected size of sig vector: " << expected_size_2D << std::endl;
    size_t expected_size_1D = static_cast<size_t>(nel) ;
    std::cerr << "Expected size of sig vector: " << expected_size_1D << std::endl;

    hid_t file_id;
    hsize_t iz_f_dims[2] = {17, 1648};
    hsize_t ite_dims[1] = {1648};
    hsize_t jne_dims[1] = {1648};
    hsize_t epatom_dims[2] = {17, 1648};
    hsize_t amamu_f_dims[2] = {17, 1648};
    hsize_t eumesh_dims[3] = {17, 1648,10000};
    hsize_t sig_dims[3] = {17, 1648, 10000};
    
    // Cast pointers to appropriate types
    int* iz_f = static_cast<int*>(iz_f_ptr);
    int* ite = static_cast<int*>(ite_ptr);
    int* jne = static_cast<int*>(jne_ptr);

    double* epatom = static_cast<double*>(epatom_ptr);
    double* amamu_f = static_cast<double*>(amamu_f_ptr);
    double* eumesh = static_cast<double*>(eumesh_ptr);
    double* sig = static_cast<double*>(sig_ptr);
   
    std::cerr << "Buffer address for iz_f after allocation: " << iz_f_ptr << std::endl;    
    std::cerr << "Buffer address for ite after allocation: " << ite_ptr << std::endl;
    std::cerr << "Buffer address for jne after allocation: " << jne_ptr << std::endl;
    std::cerr << "Buffer address for epatom after allocation: " << epatom_ptr << std::endl;
    std::cerr << "Buffer address for amamu_f after allocation: " << amamu_f_ptr << std::endl;    
    std::cerr << "Buffer address for eumesh after allocation: " << eumesh_ptr << std::endl;
    std::cerr << "Buffer address for sig after allocation: " << sig_ptr << std::endl;
    size_t total_size = 17 * 1648 * 10000 * sizeof(double);
    size_t total_size_overflow_test = nel * np_mesh * nptot* sizeof(double);
    
    //std::cerr << "Total buffer size: " << total_size << " bytes" << std::endl;
    //std::cerr << "Total (overflow test) buffer size: " << total_size_overflow_test << " bytes" << std::endl;
    
    hsize_t chunk_dims_2[2] = {1, 1};
    hsize_t chunk_dims_1D[1] = {1};
    hsize_t chunk_dims_3[3] = {1, 1, 1};



    // Open HDF5 file
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    //std::cerr << "Successfully opened file: " << filename << std::endl;

    // Read datasets
    
    
    std::cerr << "Reading dataset /iz_f" << std::endl;
    //read_dataset(file_id, "/iz_f", iz_f, iz_f_dims, 2, sizeof(int), chunk_dims_2);
    //read_dataset(file_id, "/iz_f", iz_f.data(), iz_f_dims, 2, sizeof(int), chunk_dims_2);
    read_dataset_2D(file_id,"/iz_f",iz_f,iz_f_dims,sizeof(int),chunk_dims_2) ;
    
    std::cerr << "Reading dataset /ite" << std::endl;
    //read_dataset(file_id, "/ite", ite, ite_dims, 1, sizeof(int), chunk_dims_1D);
    //read_dataset_1D(file_id,"/ite",ite.data(),ite_dims,sizeof(int),chunk_dims_1D) ;
    read_dataset_1D(file_id,"/ite",ite,ite_dims,sizeof(int),chunk_dims_1D) ;
        
    std::cerr << "Reading dataset /jne" << std::endl;
    //read_dataset(file_id, "/jne", jne, jne_dims, 1, sizeof(int), chunk_dims_1D);
    //read_dataset_1D(file_id,"/jne",jne.data(),jne_dims,sizeof(int),chunk_dims_1D) ;
    read_dataset_1D(file_id,"/jne",jne,jne_dims,sizeof(int),chunk_dims_1D) ;    
    
    std::cerr << "Reading dataset /epatom" << std::endl;
    //read_dataset(file_id, "/epatom", epatom, epatom_dims, 2, sizeof(double), chunk_dims_2);
    //read_dataset_2D(file_id,"/epatom",epatom.data(),epatom_dims,sizeof(double),chunk_dims_2) ;
    read_dataset_2D(file_id,"/epatom",epatom,epatom_dims,sizeof(double),chunk_dims_2) ;
        
    std::cerr << "Reading dataset /amamu_f" << std::endl;
    //read_dataset(file_id, "/amamu_f", amamu_f, amamu_f_dims, 2, sizeof(double), chunk_dims_2);
    //read_dataset_2D(file_id,"/amamu_f",amamu_f.data(),amamu_f_dims,sizeof(double),chunk_dims_2) ;
    read_dataset_2D(file_id,"/amamu_f",amamu_f,amamu_f_dims,sizeof(double),chunk_dims_2) ;
        
    std::cerr << "Reading dataset /eumesh" << std::endl;
    //read_dataset(file_id, "/eumesh", eumesh, eumesh_dims, 2, sizeof(double), chunk_dims_2);
    //read_dataset_3D(file_id,"/eumesh",eumesh.data(),eumesh_dims,sizeof(double),chunk_dims_3) ;
    read_dataset_3D(file_id,"/eumesh",eumesh,eumesh_dims,sizeof(double),chunk_dims_3) ;    
    
    std::cerr << "Reading dataset /sig" << std::endl;
    //read_dataset(file_id, "/sig", sig, sig_dims, 3, sizeof(double), chunk_dims_3);
    //read_dataset_3D(file_id,"/sig",sig.data(),sig_dims,sizeof(double),chunk_dims_3) ;

    read_dataset_3D(file_id,"/sig",sig_ptr,sig_dims,sizeof(double),chunk_dims_3) ;
   
    
    // Close the file
    H5Fclose(file_id);
    std::cerr << "Closed HDF5 file: " << filename << std::endl;
}
}
   

