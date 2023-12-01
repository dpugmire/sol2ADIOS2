#include <algorithm> // std::transform
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread> // sleep_for
#include <vector>

#include <hdf5.h>
#include <adios2.h>

#include <mpi.h>
MPI_Comm comm;
int rank = 0, nproc = 1;
std::vector<std::string> fnames;
size_t nZones = 0;
std::string outFile = "output.bp";

std::vector<std::string> FlowVariables = {"Density",   "Pressure",  "Q_CRITERIA", "SCHLIEREN", "VelocityX",
                                           "VelocityY", "VelocityZ", "vort_x",     "vort_y",    "vort_z"};

// std::vector<std::string> FlowVariables = {"Density", "Pressure", "CFL", "PHI", "VelocityX", "VelocityY", "VelocityZ"};

//std::vector<std::string> FlowVariables = {"P_aver", "P_ms", "Rho_aver", "Rho_ms", "U_aver", "V_aver", "W_aver",
//                                          "uu",     "uv",   "uw",       "vv",     "vw",     "ww"};


template <typename T>
void
GetADIOSVar(adios2::IO& io, const std::string& varNm, adios2::Variable<T>& var, const std::size_t& sz)
{
  auto tmp = io.InquireVariable<T>(varNm);
  if (tmp)
    var = tmp;
  else
    var = io.DefineVariable<T>(varNm, {}, {}, {sz});

  std::cout<<"GET_VAR: "<<varNm<<" sz= "<<sz<<std::endl;
}

bool ReadVariable(int rank, const std::string &name, hid_t &hfile, hid_t memtype, void *data)
{
  //std::cout << "  Open Dataset " << name << std::endl;
  hid_t dataset_id = H5Dopen2(hfile, name.c_str(), H5P_DEFAULT);
  herr_t status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(dataset_id);
  return true;
}

void ProcessArgs(int rank, int argc, char *argv[])
{
  if (argc > 2)
  {
    nZones = atoll(argv[1]);
    int n = 2;
    while (n < argc)
    {
      fnames.push_back(argv[n++]);
    }
  }
  else
  {
    std::cout << "Usage: solread-adios N <sol1.bp> [<sol2.bp> .. <solk.bp>]\n"
              << "    where N is the number of zones" << std::endl;
    MPI_Abort(comm, 1);
  }
}

std::string GetDatasetName(std::string &path1, size_t zone, std::string &path2)
{
  return path1 + " " + std::to_string(zone) + path2;
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  hid_t file_id, dataset_id; /* identifiers */
  herr_t status;

  ProcessArgs(rank, argc, argv);

  // nz: number of zones process by one process
  // process zones: nzStart .. nzStart + nz - 1
  size_t nz = nZones / nproc;
  size_t rem = nZones - (nz * nproc);
  size_t nzStart = rank * nz;
  if (rank < rem)
  {
    nzStart += rank;
    ++nz;
  }
  else
  {
    nzStart += rem;
  }
  std::cout << "Rank " << rank << " reads zones " << nzStart << ".." << nzStart + nz - 1 << std::endl;

  double tstart = MPI_Wtime();
  size_t nBytesRead = 0;

  adios2::ADIOS adios(MPI_COMM_WORLD);
  adios2::IO io = adios.DeclareIO("io");
  adios2::Engine engine = io.Open(outFile, adios2::Mode::Write);
  engine.BeginStep();

  int timeStep = 0;
  for (auto &fname : fnames)
  {
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (!rank)
    {
      std::cout << "File :" << fname << std::endl;
      std::cout << "  Zones = " << nZones << std::endl;
    }

    nBytesRead += nZones * 2 * sizeof(int32_t);

    for (size_t zone = nzStart + 1; zone < nzStart + nz + 1; ++zone)
    {
      std::string zonePath = "/hpMusic_base/hpMusic_Zone " + std::to_string(zone);
      std::cout << "  Open Group " << zonePath << std::endl;
      hid_t zoneGroupID = H5Gopen2(file_id, zonePath.c_str(), H5P_DEFAULT);
      int32_t data[3];
      std::cout << "  Read ' data' " << std::endl;
      ReadVariable(rank, " data", zoneGroupID, H5T_NATIVE_INT32, data);
      int32_t elemdata[2];
      std::cout << "  Read 'Elem/ data' " << std::endl;
      ReadVariable(rank, "Elem/ data", zoneGroupID, H5T_NATIVE_INT32, elemdata);

      const int32_t nNodes = data[0];
      const int32_t nElems = data[1];

      std::cout << "Rank " << rank << " Zone " << zone << " nElems = " << nElems << " nNodes = " << nNodes << std::endl;

      int64_t *dataEC = static_cast<int64_t *>(malloc(nElems * 8 * sizeof(int64_t)));
      nBytesRead += nElems * 8 * sizeof(int64_t);
      int64_t *dataER = static_cast<int64_t *>(malloc(2 * sizeof(int64_t)));
      nBytesRead += 2 * sizeof(int64_t);

      std::vector<double *> ptrs;
      for (auto &fv : FlowVariables)
      {
        double *ptr = static_cast<double *>(malloc(nNodes * sizeof(double)));
        ptrs.push_back(ptr);
      }
      double *gcx = static_cast<double *>(malloc(nNodes * sizeof(double)));
      double *gcy = static_cast<double *>(malloc(nNodes * sizeof(double)));
      double *gcz = static_cast<double *>(malloc(nNodes * sizeof(double)));

      nBytesRead += (ptrs.size() + 3) * nNodes * sizeof(int64_t);

      ReadVariable(rank, "Elem/ElementConnectivity/ data", zoneGroupID, H5T_NATIVE_INT64, dataEC);
      ReadVariable(rank, "Elem/ElementRange/ data", zoneGroupID, H5T_NATIVE_INT64, dataER);
      ReadVariable(rank, "GridCoordinates/CoordinateX/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcx);
      ReadVariable(rank, "GridCoordinates/CoordinateY/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcy);
      ReadVariable(rank, "GridCoordinates/CoordinateZ/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcz);
      for (int i = 0; i < FlowVariables.size(); ++i)
      {
        ReadVariable(rank, "FlowSolution/" + FlowVariables[i] + "/ data", zoneGroupID, H5T_NATIVE_DOUBLE, ptrs[i]);
      }

      //convert connectivity to 0-based indexing
      const std::size_t nConn = 8*nElems;
      for (std::size_t i = 0; i < nConn; i++)
        dataEC[i] -= 1;

      //Write out the ADIOS.
      adios2::Variable<int64_t> varConn;
      adios2::Variable<double> varCoordsX, varCoordsY, varCoordsZ;
      GetADIOSVar(io, "ElementConnectivity", varConn, static_cast<std::size_t>(8*nElems));
      GetADIOSVar(io, "GridCoordinates/CoordinateX", varCoordsX, static_cast<std::size_t>(nNodes));
      GetADIOSVar(io, "GridCoordinates/CoordinateY", varCoordsY, static_cast<std::size_t>(nNodes));
      GetADIOSVar(io, "GridCoordinates/CoordinateZ", varCoordsZ, static_cast<std::size_t>(nNodes));

      engine.Put<int64_t>(varConn, dataEC);
      engine.Put<double>(varCoordsX, gcx);
      engine.Put<double>(varCoordsY, gcy);
      engine.Put<double>(varCoordsZ, gcz);

      for (int i = 0; i < FlowVariables.size(); i++)
      {
        adios2::Variable<double> var;
        std::string varNm = "FlowSolution/" + FlowVariables[i];
        GetADIOSVar(io, varNm, var, nNodes);
        engine.Put<double>(var, ptrs[i]);
      }

      if (rank == 0)
      {
        adios2::Variable<int> varTime;
        GetADIOSVar(io, "time", varTime, 1);
        engine.Put<int>(varTime, &timeStep);
        timeStep++;
      }

      engine.EndStep();

      //free up memory.
      for (auto p : ptrs)
        free(p);
      free(dataEC);
      free(dataER);
      free(gcx);
      free(gcy);
      free(gcz);

      status = H5Gclose(zoneGroupID);
    }
    status = H5Fclose(file_id);
  }

  MPI_Barrier(comm);
  double tend = MPI_Wtime();
  size_t nbytes[nproc];
  MPI_Gather(&nBytesRead, 1, MPI_LONG_LONG, nbytes, 1, MPI_LONG_LONG, 0, comm);
  if (!rank)
  {
    std::cout << "Total time to read = " << tend - tstart << "s" << std::endl;
    size_t n = 0;
    for (int i = 0; i < nproc; ++i)
    {
      std::cout << "Rank " << i << " bytes read = " << nbytes[i] << std::endl;
      n += nbytes[i];
    }
    std::cout << "Bytes read = " << n << std::endl;
  }

  engine.Close();

  MPI_Finalize();

  return 0;
}
