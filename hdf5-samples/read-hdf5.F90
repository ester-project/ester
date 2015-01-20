program read_hdf5
    use hdf5
    implicit none

    character*100 :: file_name = "star.hdf5"
    integer status, error
    integer(hsize_t) :: dims(2)
    integer :: nr, nth
    integer(HID_T) :: fid, gid, aid, did
    double precision, allocatable :: T(:,:)

    ! init interface
    call h5open_f(status)

    ! open the HDF5 file
    call h5fopen_f(file_name, H5F_ACC_RDWR_F, fid, status)

    ! open the `star' group
    call h5gopen_f(fid, "star", gid, error)

    ! open the `nr' attribute
    call h5aopen_f(gid, "nr", aid, error)

    ! read the attribute
    dims(1) = 1
    dims(2) = 0
    call h5aread_f(aid, H5T_NATIVE_INTEGER, nr, dims, error)

    ! close the attribute
    call h5aclose_f(aid, error)

    ! open the `nth' attribute
    call h5aopen_f(gid, "nth", aid, error)

    ! read the attribute
    dims(1) = 1
    dims(2) = 0
    call h5aread_f(aid, H5T_NATIVE_INTEGER, nth, dims, error)

    ! close the attribute
    call h5aclose_f(aid, error)

    print *, "nr: ", nr
    print *, "nth:", nth


    ! allocate memory for the temperature field
    allocate(T(nr, nth))

    ! open the `T' dataset
    call h5dopen_f(gid, "T", did, error)

    ! read the field
    dims(1) = nr
    dims(2) = nth
    call h5dread_f(did, H5T_NATIVE_DOUBLE, T, dims, error)

    print *, "T at the center: ", T(1, 1)
    print *, "T at the equator:", T(nr, 1)

    deallocate(T)

    ! close dataset, group and file
    call h5dclose_f(did, error)
    call h5gclose_f(gid, error)
    call h5fclose_f(fid, error)
end program
