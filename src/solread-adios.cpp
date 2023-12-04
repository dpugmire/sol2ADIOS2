#include <algorithm> // std::transform
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread> // sleep_for
#include <vector>

#include <adios2.h>

#include <mpi.h>
MPI_Comm comm;
int rank = 0, nproc = 1;
std::vector<std::string> fnames;
size_t nZones = 0;

std::vector<std::string> FlowVariables = {"Density",   "Pressure",  "Q_CRITERIA", "SCHLIEREN", "VelocityX",
                                           "VelocityY", "VelocityZ", "vort_x",     "vort_y",    "vort_z"};

// std::vector<std::string> FlowVariables = {"Density", "Pressure", "CFL", "PHI", "VelocityX", "VelocityY", "VelocityZ"};

//std::vector<std::string> FlowVariables = {"P_aver", "P_ms", "Rho_aver", "Rho_ms", "U_aver", "V_aver", "W_aver",
//                                          "uu",     "uv",   "uw",       "vv",     "vw",     "ww"};

/* Find the block written by rank for a given variable
 * Appends the block info to the passed Varinfo vector if found
 */
template <class T>
bool ReadVariableBlock(int rank, const std::string &name, adios2::Engine &reader, adios2::IO &io, const size_t step,
                       const size_t blockIdx, T *data)
{
    std::cout << "    Rank " << rank << " zone " << blockIdx << " Variable " << name << std::endl;
    adios2::Variable<T> variable = io.InquireVariable<T>(name);
    if (!variable)
    {
        std::cout << "    ERROR: variable " << name << "not found" << std::endl;
        return false;
    }

    if (variable.ShapeID() != adios2::ShapeID::GlobalArray && variable.ShapeID() != adios2::ShapeID::LocalArray)
    {
        if (!rank)
            std::cout << "    not an array. skip" << std::endl;
        return false;
    }

    variable.SetStepSelection({step, 1});
    const std::vector<typename adios2::Variable<T>::Info> &blocks = reader.BlocksInfo(variable, step);
    auto &bi = blocks[blockIdx];
    variable.SetBlockSelection(blockIdx);
    reader.Get(variable, data, adios2::Mode::Deferred);
    return true;
}

template <class T>
bool ReadVariable(int rank, const std::string &name, adios2::Engine &reader, adios2::IO &io, const size_t step, T *data)
{
    adios2::Variable<T> variable = io.InquireVariable<T>(name);
    if (!variable)
    {
        std::cout << "    ERROR: variable " << name << "not found" << std::endl;
        return false;
    }
    std::cout << "    Rank " << rank << " Variable " << name << " count = " << variable.Count()[0] << std::endl;

    if (variable.ShapeID() != adios2::ShapeID::GlobalArray && variable.ShapeID() != adios2::ShapeID::LocalArray)
    {
        if (!rank)
            std::cout << "    not an array. skip" << std::endl;
        return false;
    }

    variable.SetStepSelection({step, 1});
    reader.Get(variable, data, adios2::Mode::Deferred);
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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    adios2::ADIOS adios(comm);

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

    int32_t nElems[nZones];
    int32_t nNodes[nZones];

    double tstart = MPI_Wtime();
    size_t nBytesRead = 0;

    for (auto &fname : fnames)
    {
        adios2::IO io = adios.DeclareIO(fname);
        adios2::Engine reader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        adios2::Variable<int32_t> vStep = io.InquireVariable<int32_t>("/hpMusic_base/step");
        size_t steps = vStep.Steps();
        if (!rank)
        {
            std::cout << "File :" << fname << std::endl;
            std::cout << "  Steps:   " << steps << std::endl;
            std::cout << "  Zones = " << nZones << std::endl;
        }

        for (size_t step = 0; step < steps; ++step)
        {
            if (!rank)
            {
                std::cout << "Step = " << step << std::endl;
            }
            ReadVariable<int32_t>(rank, "/hpMusic_base/hpMusic_Zone/n_total_node", reader, io, step, nNodes);
            ReadVariable<int32_t>(rank, "/hpMusic_base/hpMusic_Zone/n_total_elem", reader, io, step, nElems);
            reader.PerformGets();
            nBytesRead += nZones * 2 * sizeof(int32_t);

            for (size_t zone = nzStart; zone < nzStart + nz; ++zone)
            {
                int64_t *dataEC = static_cast<int64_t *>(malloc(nElems[zone] * 8 * sizeof(int64_t)));
                nBytesRead += nElems[zone] * 8 * sizeof(int64_t);
                int64_t *dataER = static_cast<int64_t *>(malloc(2 * sizeof(int64_t)));
                nBytesRead += 2 * sizeof(int64_t);

                std::vector<double *> ptrs;
                for (auto &fv : FlowVariables)
                {
                    double *ptr = static_cast<double *>(malloc(nNodes[zone] * sizeof(double)));
                    ptrs.push_back(ptr);
                }

                double *gcx = static_cast<double *>(malloc(nNodes[zone] * sizeof(double)));
                double *gcy = static_cast<double *>(malloc(nNodes[zone] * sizeof(double)));
                double *gcz = static_cast<double *>(malloc(nNodes[zone] * sizeof(double)));

                nBytesRead += (ptrs.size() + 3) * nNodes[zone] * sizeof(int64_t);

                ReadVariableBlock<int64_t>(rank, "/hpMusic_base/hpMusic_Zone/Elem/ElementConnectivity", reader, io, step, zone,
                                           dataEC);
                ReadVariableBlock<int64_t>(rank, "/hpMusic_base/hpMusic_Zone/Elem/ElementRange", reader, io, step, zone, dataER);
                ReadVariableBlock<double>(rank, "/hpMusic_base/hpMusic_Zone/GridCoordinates/CoordinateX", reader, io, step, zone,
                                          gcx);
                ReadVariableBlock<double>(rank, "/hpMusic_base/hpMusic_Zone/GridCoordinates/CoordinateY", reader, io, step, zone,
                                          gcy);
                ReadVariableBlock<double>(rank, "/hpMusic_base/hpMusic_Zone/GridCoordinates/CoordinateZ", reader, io, step, zone,
                                          gcz);

                for (int i = 0; i < FlowVariables.size(); ++i)
                {
                    ReadVariableBlock<double>(rank, "/hpMusic_base/hpMusic_Zone/FlowSolution/" + FlowVariables[i], reader, io, step,
                                              zone, ptrs[i]);
                }

                reader.PerformGets();

                for (auto p : ptrs)
                {
                    free(p);
                }
                free(dataEC);
                free(dataER);
                free(gcx);
                free(gcy);
                free(gcz);
            }
        }
        reader.Close();
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

    MPI_Finalize();

    return 0;
}
