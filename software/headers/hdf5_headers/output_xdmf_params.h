#ifndef OUTPUT_XDMF_PARAMETERS_H
#define OUTPUT_XDMF_PARAMETERS_H

/*
 * Include output_xdmf.h first so that:
 *
 *   1. HAVE_HDF5 is detected consistently.
 *   2. hdf5.h is included.
 *   3. MPI/Basilisk definitions are available.
 */
#include "output_xdmf.h"

#ifdef HAVE_HDF5

/**
 * Open an existing HDF5 file for read/write access.
 *
 * In MPI simulations, all ranks must call this function because
 * H5Fopen() is collective when using the MPI-IO file driver.
 */
static hid_t open_hdf5_file_append (const char *filename)
{
  hid_t file_id = -1;
  hid_t file_access_plist = H5P_DEFAULT;

#if _MPI
  file_access_plist = H5Pcreate (H5P_FILE_ACCESS);
  if (file_access_plist < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "open_hdf5_file_append: could not create file-access "
               "property list for '%s'\n",
               filename);
    return -1;
  }

  if (H5Pset_fapl_mpio (file_access_plist,
                        MPI_COMM_WORLD,
                        MPI_INFO_NULL) < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "open_hdf5_file_append: could not enable MPI-IO for '%s'\n",
               filename);

    H5Pclose (file_access_plist);
    return -1;
  }

  /*
   * These functions are available in HDF5 1.10 and newer.
   * They make metadata operations collective.
   */
#if (H5_VERS_MAJOR > 1) || \
    (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 10)
  H5Pset_coll_metadata_write (file_access_plist, 1);
  H5Pset_all_coll_metadata_ops (file_access_plist, 1);
#endif
#endif // _MPI

  file_id = H5Fopen (filename, H5F_ACC_RDWR, file_access_plist);

#if _MPI
  H5Pclose (file_access_plist);
#endif

  if (file_id < 0 && pid() == 0)
    fprintf (stderr,
             "open_hdf5_file_append: could not open '%s'\n",
             filename);

  return file_id;
}


/**
 * Open an existing group or create it when it does not exist.
 *
 * Under MPI, all ranks must call this function using the same group name.
 */
static hid_t open_or_create_hdf5_group (hid_t file_id,
                                        const char *group_name)
{
  htri_t exists =
    H5Lexists (file_id, group_name, H5P_DEFAULT);

  if (exists < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "open_or_create_hdf5_group: could not test group '%s'\n",
               group_name);
    return -1;
  }

  hid_t group_id;

  if (exists > 0)
    group_id = H5Gopen2 (file_id,
                         group_name,
                         H5P_DEFAULT);
  else
    group_id = H5Gcreate2 (file_id,
                           group_name,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

  if (group_id < 0 && pid() == 0)
    fprintf (stderr,
             "open_or_create_hdf5_group: could not open or create '%s'\n",
             group_name);

  return group_id;
}


/**
 * Write one double-precision scalar dataset.
 *
 * Dataset creation/opening is performed by all MPI ranks. Rank zero writes
 * the actual scalar value because every MPI rank has the same simulation
 * parameter and only one copy is needed.
 *
 * Existing datasets are overwritten rather than recreated.
 */
static int write_hdf5_scalar_double (hid_t group_id,
                                     const char *dataset_name,
                                     double value)
{
  hid_t dataset_id = -1;
  hid_t dataspace_id = -1;

  htri_t exists =
    H5Lexists (group_id, dataset_name, H5P_DEFAULT);

  if (exists < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "write_hdf5_scalar_double: could not test dataset '%s'\n",
               dataset_name);
    return -1;
  }

  if (exists > 0) {
    dataset_id = H5Dopen2 (group_id,
                           dataset_name,
                           H5P_DEFAULT);
  }
  else {
    /*
     * H5S_SCALAR represents a true scalar dataset. Its HDF5 shape is (),
     * rather than (1,) or (1,1).
     */
    dataspace_id = H5Screate (H5S_SCALAR);

    if (dataspace_id < 0) {
      if (pid() == 0)
        fprintf (stderr,
                 "write_hdf5_scalar_double: could not create dataspace "
                 "for '%s'\n",
                 dataset_name);
      return -1;
    }

    dataset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             H5T_IEEE_F64LE,
                             dataspace_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);

    H5Sclose (dataspace_id);
  }

  if (dataset_id < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "write_hdf5_scalar_double: could not open or create '%s'\n",
               dataset_name);
    return -1;
  }

  int status = 0;

  /*
   * The datasets are tiny, so rank zero performs the data write.
   * Dataset creation and opening are still called collectively above.
   */
  if (pid() == 0) {
    herr_t write_status =
      H5Dwrite (dataset_id,
                H5T_NATIVE_DOUBLE,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &value);

    if (write_status < 0) {
      fprintf (stderr,
               "write_hdf5_scalar_double: could not write '%s'\n",
               dataset_name);
      status = -1;
    }
  }

#if _MPI
  /*
   * Ensure that all ranks leave this routine consistently before the
   * dataset is closed.
   */
  MPI_Bcast (&status,
             1,
             MPI_INT,
             0,
             MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);
#endif

  H5Dclose (dataset_id);

  return status;
}


/**
 * Append scalar simulation parameters to an existing XDMF/HDF5 file.
 *
 * Arguments
 * ---------
 *
 * subname:
 *     The same base name passed to output_xmf(). For example, if
 *
 *         subname = "output/snapshot_0.0100"
 *
 *     the routine opens:
 *
 *         output/snapshot_0.0100.h5
 *
 * parameter_names:
 *     Array containing dataset names.
 *
 * parameter_values:
 *     Array containing the corresponding values.
 *
 * number_of_parameters:
 *     Number of entries in both arrays.
 */
trace
int append_xmf_parameters (const char *subname,
                           const char **parameter_names,
                           const double *parameter_values,
                           int number_of_parameters)
{
  if (!subname ||
      !parameter_names ||
      !parameter_values ||
      number_of_parameters <= 0) {
    if (pid() == 0)
      fprintf (stderr,
               "append_xmf_parameters: invalid input arguments\n");
    return -1;
  }

  char filename[260];

  int required =
    snprintf (filename,
              sizeof(filename),
              "%s.h5",
              subname);

  if (required < 0 || required >= (int) sizeof(filename)) {
    if (pid() == 0)
      fprintf (stderr,
               "append_xmf_parameters: filename is too long\n");
    return -1;
  }

  /*
   * output_xmf() closes the HDF5 file before returning. The barrier makes
   * that completion explicit before the same file is collectively reopened.
   */
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  hid_t file_id =
    open_hdf5_file_append (filename);

  if (file_id < 0)
    return -1;

  hid_t group_id =
    open_or_create_hdf5_group (file_id,
                               "/Parameters");

  if (group_id < 0) {
    H5Fclose (file_id);
    return -1;
  }

  int status = 0;

  for (int i = 0; i < number_of_parameters; i++) {
    if (!parameter_names[i] ||
        parameter_names[i][0] == '\0') {
      if (pid() == 0)
        fprintf (stderr,
                 "append_xmf_parameters: parameter %d has no name\n",
                 i);

      status = -1;
      continue;
    }

    if (write_hdf5_scalar_double (group_id,
                                  parameter_names[i],
                                  parameter_values[i]) < 0)
      status = -1;
  }

  H5Gclose (group_id);

  /*
   * Ensure that written data reaches the file before it is closed.
   */
  H5Fflush (file_id, H5F_SCOPE_GLOBAL);
  H5Fclose (file_id);

#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  return status;
}

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_PARAMETERS_H
