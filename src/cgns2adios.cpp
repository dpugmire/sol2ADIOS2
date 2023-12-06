#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <hdf5.h>
#include <adios2.h>
#include <stdexcept>

#include <mpi.h>


struct State
{
  MPI_Comm comm;
  int rank = 0, nproc = 1;
  std::vector<std::string> fnames;
  size_t nZones = 0;
  std::string outFile = "output.bp";
};

template <typename T>
void
GetADIOSVar(adios2::IO& io, const std::string& varNm, adios2::Variable<T>& var, const std::size_t& sz)
{
  auto tmp = io.InquireVariable<T>(varNm);
  if (tmp)
  {
    var = tmp;
    var.SetSelection({{},{sz}});
  }
  else
    var = io.DefineVariable<T>(varNm, {}, {}, {sz});
}

bool ReadVariable(int rank, const std::string& zoneName, const std::string &name, hid_t &hfile, hid_t memtype, void *data)
{
  //std::cout << "  Open Dataset " << name << std::endl;
  hid_t dataset_id = H5Dopen2(hfile, name.c_str(), H5P_DEFAULT);
  if (dataset_id == H5I_INVALID_HID)
    throw std::runtime_error(std::string("H5Dopen2 fail: " + zoneName + " " + name));
  herr_t status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0)
    throw std::runtime_error(std::string("H5read fail: " + zoneName + " " + name));

  status = H5Dclose(dataset_id);
  return true;
}

void ProcessArgs(State& state, int argc, char *argv[])
{
  if (argc > 2)
  {
    state.nZones = atoll(argv[1]);
    state.outFile = argv[2];
    int n = 3;
    while (n < argc)
      state.fnames.push_back(argv[n++]);
  }
  else
  {
    std::cout<<"Usage: solread-adios N outputFile <sol1.cgns> [<sol2.cgns> .. <solk.cgns>]"<<std::endl;
    std::cout<<"  where N is the number of zones" << std::endl;
    std::cout<<"  outputFile is the ADIOS file for output"<<std::endl;
    std::cout<<"  sol1.cgns ... are the names of the input CGNS files"<<std::endl;
    MPI_Abort(state.comm, 1);
  }
}

std::vector<std::string>
ReadVariableNames(hid_t& fileID)
{
  std::vector<std::string> varNames;

  //query all of the groups in zone 1 get to the variables names.
  hid_t zoneGroupID = H5Gopen2(fileID, "/hpMusic_base/hpMusic_Zone 1/FlowSolution", H5P_DEFAULT);
  hsize_t num;
  H5Gget_num_objs(zoneGroupID, &num);
  //std::cout<<"ReadVariableNames: "<<zoneGroupID<<" num= "<<num<<std::endl;
  for (hsize_t i = 0; i < num; i++)
  {
    char buff[128];
    size_t n;
    H5Gget_objname_by_idx(zoneGroupID, i, buff, n);
    varNames.push_back(buff);
    //std::cout<<"  "<<i<<" "<<buff<<std::endl;
  }

  return varNames;
}

int main(int argc, char *argv[])
{
  State state;

  MPI_Init(&argc, &argv);
  state.comm = MPI_COMM_WORLD;
  MPI_Comm_rank(state.comm, &state.rank);
  MPI_Comm_size(state.comm, &state.nproc);

  ProcessArgs(state, argc, argv);

  // nz: number of zones process by one process
  // process zones: nzStart .. nzStart + nz - 1
  size_t nz = state.nZones / state.nproc;
  size_t rem = state.nZones - (nz * state.nproc);
  size_t nzStart = state.rank * nz;
  if (state.rank < rem)
  {
    nzStart += state.rank;
    ++nz;
  }
  else
    nzStart += rem;

  std::cout << "Rank " << state.rank << " reads zones " << nzStart << ".." << nzStart + nz - 1 << std::endl;

  double tstart = MPI_Wtime();
  size_t nBytesRead = 0;

  adios2::ADIOS adios(MPI_COMM_WORLD);
  adios2::IO io = adios.DeclareIO("io");
  adios2::Engine engine = io.Open(state.outFile, adios2::Mode::Write);

  int timeStep = 0;
  std::vector<std::string> flowVariables;
  for (auto &fname : state.fnames)
  {
    engine.BeginStep();

    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == H5I_INVALID_HID)
      throw std::runtime_error("File not found: " + fname);

    if (timeStep == 0)
      flowVariables = ReadVariableNames(file_id);

    if (!state.rank)
      std::cout<<"File :" << fname << " Zones= "<<state.nZones<<std::endl;

    nBytesRead += state.nZones * 2 * sizeof(int32_t);

    bool firstIteration = true;
    for (size_t zone = nzStart + 1; zone < nzStart + nz + 1; ++zone)
    {
      std::string zonePath = "/hpMusic_base/hpMusic_Zone " + std::to_string(zone);
      //std::cout << "  Open Group " << zonePath << std::endl;
      hid_t zoneGroupID = H5Gopen2(file_id, zonePath.c_str(), H5P_DEFAULT);
      if (zoneGroupID == H5I_INVALID_HID)
        throw std::runtime_error(std::string("Zone H5Dopen2 fail: " + zonePath));

      int32_t data[3];
      ReadVariable(state.rank, zonePath, " data", zoneGroupID, H5T_NATIVE_INT32, data);
      int32_t elemdata[2];
      ReadVariable(state.rank, zonePath, "Elem/ data", zoneGroupID, H5T_NATIVE_INT32, elemdata);

      const int32_t nNodes = data[0];
      const int32_t nElems = data[1];

      std::cout << "Rank " << state.rank << " Zone " << zone << " nElems = " << nElems << " nNodes = " << nNodes << std::endl;

      int64_t *dataEC = static_cast<int64_t *>(malloc(nElems * 8 * sizeof(int64_t)));
      nBytesRead += nElems * 8 * sizeof(int64_t);
      int64_t *dataER = static_cast<int64_t *>(malloc(2 * sizeof(int64_t)));
      nBytesRead += 2 * sizeof(int64_t);

      std::vector<double *> ptrs;
      for (auto &fv : flowVariables)
      {
        double *ptr = static_cast<double *>(malloc(nNodes * sizeof(double)));
        ptrs.push_back(ptr);
      }
      double *gcx = static_cast<double *>(malloc(nNodes * sizeof(double)));
      double *gcy = static_cast<double *>(malloc(nNodes * sizeof(double)));
      double *gcz = static_cast<double *>(malloc(nNodes * sizeof(double)));

      nBytesRead += (ptrs.size() + 3) * nNodes * sizeof(int64_t);

      ReadVariable(state.rank, zonePath, "Elem/ElementConnectivity/ data", zoneGroupID, H5T_NATIVE_INT64, dataEC);
      ReadVariable(state.rank, zonePath, "Elem/ElementRange/ data", zoneGroupID, H5T_NATIVE_INT64, dataER);
      ReadVariable(state.rank, zonePath, "GridCoordinates/CoordinateX/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcx);
      ReadVariable(state.rank, zonePath, "GridCoordinates/CoordinateY/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcy);
      ReadVariable(state.rank, zonePath, "GridCoordinates/CoordinateZ/ data", zoneGroupID, H5T_NATIVE_DOUBLE, gcz);
      for (int i = 0; i < flowVariables.size(); ++i)
      {
        ReadVariable(state.rank, zonePath, "FlowSolution/" + flowVariables[i] + "/ data", zoneGroupID, H5T_NATIVE_DOUBLE, ptrs[i]);
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

      engine.Put<int64_t>(varConn, dataEC, adios2::Mode::Sync);
      engine.Put<double>(varCoordsX, gcx, adios2::Mode::Sync);
      engine.Put<double>(varCoordsY, gcy, adios2::Mode::Sync);
      engine.Put<double>(varCoordsZ, gcz, adios2::Mode::Sync);

      for (int i = 0; i < flowVariables.size(); i++)
      {
        adios2::Variable<double> var;
        std::string varNm = "FlowSolution/" + flowVariables[i];
        GetADIOSVar(io, varNm, var, nNodes);
        engine.Put<double>(var, ptrs[i], adios2::Mode::Sync);
      }

      if (state.rank == 0 && firstIteration)
      {
        adios2::Variable<int> varTime;
        GetADIOSVar(io, "time", varTime, 1);
        engine.Put<int>(varTime, &timeStep, adios2::Mode::Sync);
        timeStep++;
      }

      //free up memory.
      for (auto p : ptrs)
        free(p);
      free(dataEC);
      free(dataER);
      free(gcx);
      free(gcy);
      free(gcz);

      auto status = H5Gclose(zoneGroupID);
      firstIteration = false;
    }
    engine.EndStep();
    auto status = H5Fclose(file_id);
  }

  MPI_Barrier(state.comm);
  double tend = MPI_Wtime();
  size_t nbytes[state.nproc];
  MPI_Gather(&nBytesRead, 1, MPI_LONG_LONG, nbytes, 1, MPI_LONG_LONG, 0, state.comm);
  if (!state.rank)
  {
    std::cout << "Total time to read = " << tend - tstart << "s" << std::endl;
    size_t n = 0;
    for (int i = 0; i < state.nproc; ++i)
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
