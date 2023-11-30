#include <algorithm> // std::transform
#include <assert.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <thread> // sleep_for
#include <utility>
#include <vector>

#include <mpi.h>

#include <adios2.h>
#include <hdf5.h>

/* Command line arguments */
std::string bpname;   // input file
std::string cgnsname; // output file
size_t step = 0;      // which step to convert

const std::string flowprefix = "/hpMusic_base/hpMusic_Zone/FlowSolution";
const std::string zoneprefix = "/hpMusic_base/hpMusic_Zone";

MPI_Comm comm;
int rank = 0, nproc = 1;
size_t nZones = 0;
size_t nSteps = 1;
std::vector<int32_t> nElems; // mesh: # of elements per zone
std::vector<int32_t> nNodes; // mesh: # of grid points per zone
std::string ElementType;     // QUAD_4 or HEXA_8
std::vector<std::string> FlowVariables;

struct BoundaryInfo
{
    size_t nZones;
    std::string ElementType;
    std::vector<size_t> nConnectivity; // mesh: size of connectivity per zone
};

std::map<std::string, BoundaryInfo> Boundaries;

/* HDF5 variables */
constexpr int FAIL = -1;
hid_t file_id;
hid_t groupBase;
struct ZoneHierarchy
{
    hid_t groupZone;
    hid_t groupElem;
    hid_t groupElementConnectivity;
    hid_t groupElementRange;
    hid_t groupGridCoordinates;
    hid_t groupCoordinateX;
    hid_t groupCoordinateY;
    hid_t groupCoordinateZ;
    hid_t groupZoneType;
    hid_t groupFlowSolution;
    std::vector<hid_t> groupFlowSolutionVariables;

    hid_t dsZone_data;
    hid_t dsElem_data; // ElemType, e.g. quad or hexa
    hid_t dsElementConnectivity_data;
    hid_t dsElementRange_data;
    hid_t dsCoordinateX_data;
    hid_t dsCoordinateY_data;
    hid_t dsCoordinateZ_data;
    hid_t dsZoneType_data;
    std::vector<hid_t> dsFlowSolutionVariables_data;
};

std::vector<ZoneHierarchy> zoneHierarchy;

struct BoundaryHierarchy
{
    hid_t groupBoundary;
    hid_t groupElementConnectivity;
    hid_t groupElementRange;

    hid_t dsBoundary_data; // ElemType, e.g. quad or hexa
    hid_t dsElementConnectivity_data;
    hid_t dsElementRange_data;
};

struct BoundaryGroup
{
    std::string boundaryName;
    BoundaryInfo *boundaryInfo;
    std::vector<BoundaryHierarchy> boundaryHierarchy;
};

std::vector<BoundaryGroup> boundaryGroups;

int ElementTypeStringToInt(const std::string &et)
{
    if (et == "QUAD_4")
    {
        return 7;
    }
    else if (et == "HEXA_8")
    {
        return 17;
    }
    return 99999999;
}

int ElementTypeStringToNodes(const std::string &et)
{
    if (et == "QUAD_4")
    {
        return 4;
    }
    else if (et == "HEXA_8")
    {
        return 8;
    }
    return 1;
}

/* Collective call */
hid_t CreateDataset(const std::string name, hid_t group, hid_t type, std::vector<hsize_t> sizes)
{
    const hsize_t *dims = sizes.data();
    int ndim = sizes.size();
    hid_t dspace = H5Screate_simple(ndim, dims, NULL);
    hid_t ds = H5Dcreate2(group, name.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Sclose(dspace);
    return ds;
}

/* Single process call */
void WriteDataset(hid_t id, hid_t type, void *data)
{
    herr_t status = H5Dwrite(id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    assert(status != FAIL);
}

/* Collective call */
void CloseDataset(hid_t id)
{
    herr_t status = H5Dclose(id);
    assert(status != FAIL);
}

void CreateAttribute(hid_t id, const std::string name, hsize_t maxsize, const std::string value)
{
    hid_t acpl = H5Pcreate(H5P_ATTRIBUTE_CREATE);
    assert(acpl != H5I_INVALID_HID);
    herr_t status = H5Pset_char_encoding(acpl, H5T_CSET_ASCII);
    assert(status >= 0);
    hid_t fspace = H5Screate(H5S_SCALAR);
    assert(fspace != H5I_INVALID_HID);
    // create a fixed size string type
    hid_t attr_type = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(attr_type, maxsize);
    assert(status >= 0);

    // std::cout << "Create attribute " << name << std::endl;
    hid_t attr = H5Acreate2(id, name.c_str(), attr_type, fspace, acpl, H5P_DEFAULT);
    assert(attr != H5I_INVALID_HID);
    if (!rank)
    {
        // std::cout << "Write attribute value " << value << std::endl;
        status = H5Awrite(attr, attr_type, value.c_str());
        assert(status >= 0);
    }
    H5Sclose(fspace);
    H5Aclose(attr);
}

void CreateAttribute(hid_t id, const std::string name, hid_t memtype, hid_t filetype, hsize_t size, const void *data)
{
    herr_t status;
    hsize_t dims[] = {size};
    hid_t dspace = H5Screate_simple(1, dims, NULL);
    // std::cout << "Create attribute " << name << std::endl;
    hid_t attribute_id = H5Acreate2(id, name.c_str(), filetype, dspace, H5P_DEFAULT, H5P_DEFAULT);
    if (!rank)
    {
        status = H5Awrite(attribute_id, memtype, data);
        status = H5Aclose(attribute_id);
    }
    status = H5Sclose(dspace);
}

hid_t CreateGroup(hid_t parent, const std::string name, const std::string labelAttribute, const std::string typeAttribute,
                  const int *flags = nullptr)
{
    hid_t group = H5Gcreate2(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (flags)
    {
        CreateAttribute(group, "flags", H5T_NATIVE_INT, H5T_STD_I32LE, 1, flags);
    }
    CreateAttribute(group, "label", 33, labelAttribute.c_str());
    CreateAttribute(group, "name", 33, name.c_str());
    CreateAttribute(group, "type", 3, typeAttribute.c_str());
    return group;
}

void CreateHDF5OutputStructure()
{
    hid_t acc_tpl1;
    herr_t status;
    MPI_Info info = MPI_INFO_NULL;
    hid_t ds_format, ds_hdf5version, ds_CGNSLibraryVersion, ds_Base;
    hid_t groupCGNSLibraryVersion;

    if (!rank)
    {
        std::cout << "Create CGNS output hierarchy " << std::endl;
    }

    /* setup file access template with parallel IO access. */
    acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
    assert(acc_tpl1 != FAIL);
    std::cout << "H5Pcreate access succeed\n";
    /* set Parallel access with communicator */
    status = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    assert(status != FAIL);
    std::cout << "H5Pset_fapl_mpio succeed\n";
    /* create the file collectively */
    file_id = H5Fcreate(cgnsname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
    assert(file_id != FAIL);
    std::cout << "H5Fcreate succeed\n";
    /* Release file-access template */
    status = H5Pclose(acc_tpl1);
    assert(status != FAIL);

    int groupFlags[1] = {1};

    std::cout << "rank " << rank << ": Create group hierarchy collectively\n";
    groupBase = CreateGroup(file_id, "hpMusic_base", "CGNSBase_t", "I4", groupFlags);
    ds_Base = CreateDataset(" data", groupBase, H5T_STD_I32LE, {2});

    groupCGNSLibraryVersion = CreateGroup(file_id, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", groupFlags);
    ds_CGNSLibraryVersion = CreateDataset(" data", groupCGNSLibraryVersion, H5T_IEEE_F32LE, {1});
    ds_format = CreateDataset(" format", file_id, H5T_STD_U8LE, {15});
    ds_hdf5version = CreateDataset(" hdf5version", file_id, H5T_STD_U8LE, {33});

    zoneHierarchy.resize(nZones);
    for (int zone = 0; zone < nZones; ++zone)
    {
        std::string name = "hpMusic_Zone " + std::to_string(zone + 1);
        std::vector<size_t> nConnectivitySize = {(size_t)(ElementTypeStringToNodes(ElementType) * nElems[zone])};
        std::vector<size_t> nNodeSize = {(size_t)nNodes[zone]};

        zoneHierarchy[zone].groupZone = CreateGroup(groupBase, name, "Zone_t", "I8", groupFlags);
        zoneHierarchy[zone].dsZone_data = CreateDataset(" data", zoneHierarchy[zone].groupZone, H5T_STD_I64LE, {3, 1});
        zoneHierarchy[zone].groupElem = CreateGroup(zoneHierarchy[zone].groupZone, "Elem", "Elements_t", "I4", groupFlags);
        zoneHierarchy[zone].dsElem_data = CreateDataset(" data", zoneHierarchy[zone].groupElem, H5T_STD_I32LE, {2});
        zoneHierarchy[zone].groupElementConnectivity =
            CreateGroup(zoneHierarchy[zone].groupElem, "ElementConnectivity", "DataArray_t", "I8", groupFlags);
        zoneHierarchy[zone].dsElementConnectivity_data =
            CreateDataset(" data", zoneHierarchy[zone].groupElementConnectivity, H5T_STD_I64LE, nConnectivitySize);
        zoneHierarchy[zone].groupElementRange =
            CreateGroup(zoneHierarchy[zone].groupElem, "ElementRange", "IndexRange_t", "I8", groupFlags);
        zoneHierarchy[zone].dsElementRange_data = CreateDataset(" data", zoneHierarchy[zone].groupElementRange, H5T_STD_I64LE, {2});
        zoneHierarchy[zone].groupGridCoordinates =
            CreateGroup(zoneHierarchy[zone].groupZone, "GridCoordinates", "GridCoordinates_t", "MT", groupFlags);
        zoneHierarchy[zone].groupCoordinateX =
            CreateGroup(zoneHierarchy[zone].groupGridCoordinates, "CoordinateX", "DataArray_t", "R8", groupFlags);
        zoneHierarchy[zone].dsCoordinateX_data =
            CreateDataset(" data", zoneHierarchy[zone].groupCoordinateX, H5T_IEEE_F64LE, nNodeSize);
        zoneHierarchy[zone].groupCoordinateY =
            CreateGroup(zoneHierarchy[zone].groupGridCoordinates, "CoordinateY", "DataArray_t", "R8", groupFlags);
        zoneHierarchy[zone].dsCoordinateY_data =
            CreateDataset(" data", zoneHierarchy[zone].groupCoordinateY, H5T_IEEE_F64LE, nNodeSize);
        zoneHierarchy[zone].groupCoordinateZ =
            CreateGroup(zoneHierarchy[zone].groupGridCoordinates, "CoordinateZ", "DataArray_t", "R8", groupFlags);
        zoneHierarchy[zone].dsCoordinateZ_data =
            CreateDataset(" data", zoneHierarchy[zone].groupCoordinateZ, H5T_IEEE_F64LE, nNodeSize);
        zoneHierarchy[zone].groupZoneType = CreateGroup(zoneHierarchy[zone].groupZone, "ZoneType", "ZoneType_t", "C1", groupFlags);
        zoneHierarchy[zone].dsZoneType_data = CreateDataset(" data", zoneHierarchy[zone].groupZoneType, H5T_STD_U8LE, {12});

        zoneHierarchy[zone].groupFlowSolution =
            CreateGroup(zoneHierarchy[zone].groupZone, "FlowSolution", "FlowSolution_t", "MT", groupFlags);
        for (auto &fv : FlowVariables)
        {
            std::string name = fv.substr(flowprefix.length() + 1);
            // std::cout << "Create group for flow variable " << name << std::endl;
            hid_t gid = CreateGroup(zoneHierarchy[zone].groupFlowSolution, name.c_str(), "DataArray_t", "R8", groupFlags);
            assert(gid != FAIL);
            zoneHierarchy[zone].groupFlowSolutionVariables.push_back(gid);
            hid_t dsid = CreateDataset(" data", gid, H5T_IEEE_F64LE, nNodeSize);
            zoneHierarchy[zone].dsFlowSolutionVariables_data.push_back(dsid);
        }
    }

    boundaryGroups.reserve(Boundaries.size());
    if (!rank)
        std::cout << "Create hierarchies for " << Boundaries.size() << " boundaries\n";
    for (auto &bv : Boundaries)
    {
        const std::string &BoundaryName = bv.first;
        BoundaryInfo &bi = bv.second;
        BoundaryGroup bg;
        bg.boundaryInfo = &bi;
        bg.boundaryName = BoundaryName;
        bg.boundaryHierarchy.reserve(bi.nZones);
        if (!rank)
            std::cout << "    Boundary " << BoundaryName << " zones = " << bi.nZones << "\n";
        for (int bzone = 0; bzone < bi.nZones; ++bzone)
        {
            BoundaryHierarchy bh;
            if (!rank)
                std::cout << "        zone " << bzone << " nelems = " << bi.nConnectivity[bzone] << "\n";

            bh.groupBoundary = CreateGroup(zoneHierarchy[bzone].groupZone, BoundaryName.c_str(), "Elements_t", "I4", groupFlags);
            bh.dsBoundary_data = CreateDataset(" data", bh.groupBoundary, H5T_STD_I32LE, {2});
            bh.groupElementConnectivity = CreateGroup(bh.groupBoundary, "ElementConnectivity", "DataArray_t", "I8", groupFlags);
            bh.dsElementConnectivity_data =
                CreateDataset(" data", bh.groupElementConnectivity, H5T_STD_I64LE, {bi.nConnectivity[bzone]});
            bh.groupElementRange = CreateGroup(bh.groupBoundary, "ElementRange", "IndexRange_t", "I8", groupFlags);
            bh.dsElementRange_data = CreateDataset(" data", bh.groupElementRange, H5T_STD_I64LE, {2});
            bg.boundaryHierarchy.push_back(bh);
        }
        boundaryGroups.push_back(bg);
    }

    /* rank 0 write informational data */
    if (!rank)
    {
        int basedata[2] = {3, 3};
        WriteDataset(ds_Base, H5T_NATIVE_INT, basedata);
        float cgnsversion[1] = {3.3};
        WriteDataset(ds_CGNSLibraryVersion, H5T_NATIVE_FLOAT, cgnsversion);
        char format[] = "IEEE_LITTLE_32";
        WriteDataset(ds_format, H5T_NATIVE_CHAR, format);
        char h5ver[] = H5_VERS_INFO;
        WriteDataset(ds_hdf5version, H5T_NATIVE_CHAR, h5ver);
        for (int zone = 0; zone < nZones; ++zone)
        {
            uint64_t data[3] = {(uint64_t)nNodes[zone], (uint64_t)nElems[zone], 0};
            WriteDataset(zoneHierarchy[zone].dsZone_data, H5T_NATIVE_UINT64, data);
            int32_t elemtype[2] = {ElementTypeStringToInt(ElementType), 0};
            WriteDataset(zoneHierarchy[zone].dsElem_data, H5T_NATIVE_INT32, elemtype);
            char zonetype[] = "Unstructured";
            WriteDataset(zoneHierarchy[zone].dsZoneType_data, H5T_NATIVE_CHAR, zonetype);
        }
        for (auto &bg : boundaryGroups)
        {
            for (auto &bh : bg.boundaryHierarchy)
            {
                int32_t elemtype[2] = {ElementTypeStringToInt(bg.boundaryInfo->ElementType), 0};
                WriteDataset(bh.dsBoundary_data, H5T_NATIVE_INT32, elemtype);
            }
        }
    }
    MPI_Barrier(comm);
    CloseDataset(ds_Base);
    CloseDataset(ds_CGNSLibraryVersion);

    status = H5Gclose(groupCGNSLibraryVersion);
    assert(status != FAIL);

    CloseDataset(ds_format);
    CloseDataset(ds_hdf5version);

    /* Close some datasets */
    for (int zone = 0; zone < nZones; ++zone)
    {
        CloseDataset(zoneHierarchy[zone].dsZone_data);
        CloseDataset(zoneHierarchy[zone].dsElem_data);
        CloseDataset(zoneHierarchy[zone].dsZoneType_data);
    }
    for (auto &bg : boundaryGroups)
    {
        for (auto &bh : bg.boundaryHierarchy)
        {
            CloseDataset(bh.dsBoundary_data);
        }
    }
}

void CloseHDF5()
{
    herr_t status;
    if (!rank)
    {
        std::cout << "Close CGNS output " << std::endl;
    }
    /* Close the datasets */
    for (int zone = 0; zone < nZones; ++zone)
    {
        CloseDataset(zoneHierarchy[zone].dsElementConnectivity_data);
        CloseDataset(zoneHierarchy[zone].dsElementRange_data);
        CloseDataset(zoneHierarchy[zone].dsCoordinateX_data);
        CloseDataset(zoneHierarchy[zone].dsCoordinateY_data);
        CloseDataset(zoneHierarchy[zone].dsCoordinateZ_data);
        for (auto &dsfv : zoneHierarchy[zone].dsFlowSolutionVariables_data)
        {
            CloseDataset(dsfv);
        }
    }
    for (auto &bg : boundaryGroups)
    {
        for (auto &bh : bg.boundaryHierarchy)
        {
            CloseDataset(bh.dsElementConnectivity_data);
            CloseDataset(bh.dsElementRange_data);
        }
    }

    /* Close the groups */
    for (int zone = 0; zone < nZones; ++zone)
    {
        status = H5Gclose(zoneHierarchy[zone].groupElementConnectivity);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupElementRange);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupElem);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupGridCoordinates);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupCoordinateX);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupCoordinateY);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupCoordinateZ);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupZoneType);
        assert(status != FAIL);
        for (auto &gfv : zoneHierarchy[zone].groupFlowSolutionVariables)
        {
            status = H5Gclose(gfv);
            assert(status != FAIL);
        }
        status = H5Gclose(zoneHierarchy[zone].groupFlowSolution);
        assert(status != FAIL);
        status = H5Gclose(zoneHierarchy[zone].groupZone);
        assert(status != FAIL);
    }

    for (auto &bg : boundaryGroups)
    {
        for (auto &bh : bg.boundaryHierarchy)
        {
            status = H5Gclose(bh.groupElementConnectivity);
            assert(status != FAIL);
            status = H5Gclose(bh.groupElementRange);
            assert(status != FAIL);
            status = H5Gclose(bh.groupBoundary);
            assert(status != FAIL);
        }
    }

    status = H5Gclose(groupBase);
    assert(status != FAIL);

    /* Close the file. */
    status = H5Fclose(file_id);
    assert(status != FAIL);
}

/* Find the block written by rank for a given variable
 * Appends the block info to the passed Varinfo vector if found
 */
template <class T>
bool ReadVariableBlock(int rank, const std::string &name, adios2::Engine &reader, adios2::IO &io, const size_t step,
                       const size_t blockIdx, T *data)
{
    // std::cout << "    Rank " << rank << " zone " << blockIdx << " Variable " << name << std::endl;
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
        bpname = std::string(argv[1]);
        cgnsname = std::string(argv[2]);
        if (argc > 3)
        {
            step = atoll(argv[3]);
        }
        else
        {
            step = 0;
        }
    }
    else
    {
        std::cout << "Usage: solread-adios sol.bp sol.cgns [step]\n"
                  << "    Only one step is converted at a time" << std::endl;
        MPI_Abort(comm, 1);
    }
}

/* Determine number of steps, zones, and list of Flow variables */
void PreProcess(adios2::Engine &reader, adios2::IO &io, size_t step)
{
    auto v = io.InquireVariable<double>("/hpMusic_base/hpMusic_Zone/GridCoordinates/CoordinateX");
    nSteps = v.Steps();
    const std::vector<typename adios2::Variable<double>::Info> &blocks = reader.BlocksInfo(v, step);
    nZones = blocks.size();

    auto varlist = io.AvailableVariables(true);

    for (auto &v : varlist)
    {
        if (!v.first.compare(0, flowprefix.length(), flowprefix))
        {
            FlowVariables.push_back(v.first);
        }
        else if (!v.first.compare(0, zoneprefix.length(), zoneprefix))
        {
            std::string name = v.first.substr(zoneprefix.length() + 1);
            if (name != "n_total_elem" && name != "n_total_node" && name.substr(0, 4).compare("Elem") &&
                name.substr(0, 15).compare("GridCoordinates"))
            {
                std::string b = name.substr(0, name.find_first_of('/'));
                auto ret = Boundaries.insert({b, BoundaryInfo()});
                // at first find of the boundary, get the number of zones as well
                if (ret.second)
                {
                    BoundaryInfo &bi = ret.first->second;
                    auto v = io.InquireVariable<int64_t>(zoneprefix + "/" + b + "/ElementConnectivity");
                    const std::vector<typename adios2::Variable<int64_t>::Info> &blocks = reader.BlocksInfo(v, step);
                    bi.nZones = blocks.size();
                    for (auto &b : blocks)
                    {
                        bi.nConnectivity.push_back(b.Count[0]);
                    }
                    auto a = io.InquireAttribute<std::string>(zoneprefix + "/" + b + "/ElementType");
                    bi.ElementType = a.Data()[0];
                }
            }
        }
    }

    auto a = io.InquireAttribute<std::string>("/hpMusic_base/hpMusic_Zone/Elem/ElementType");
    ElementType = a.Data()[0];

    if (!rank)
    {
        std::cout << "BP File Info: "
                  << "\n   steps = " << nSteps << "\n   zones = " << nZones << "\n   element type = " << ElementType
                  << "\n   flow variables = {";
        for (auto &fv : FlowVariables)
        {
            std::cout << " " << fv.substr(flowprefix.length() + 1);
        }
        std::cout << "}\n   boundaries = {";
        for (auto &bv : Boundaries)
        {
            std::cout << "\n        " << bv.first << ", with zones = " << bv.second.nZones
                      << " element type = " << bv.second.ElementType;
        }
        std::cout << "\n   }" << std::endl;
    }
}

// return: start, n  for the rank
std::pair<size_t, size_t> decomposeZones(int nproc, int rank, size_t nZones)
{
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
    return {nzStart, nz};
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    adios2::ADIOS adios(comm);

    ProcessArgs(rank, argc, argv);

    adios2::IO io = adios.DeclareIO(bpname);
    adios2::Engine reader = io.Open(bpname, adios2::Mode::ReadRandomAccess);

    PreProcess(reader, io, step);

    if (step >= nSteps)
    {
        if (!rank)
        {
            std::cout << "ERROR: selected step is not available: " << step << ". Pick a step from 0 to " << nSteps - 1 << std::endl;
            reader.Close();
            MPI_Abort(comm, 1);
        }
    }

    // nz: number of zones process by one process
    // process zones: nzStart .. nzStart + nz - 1
    auto p = decomposeZones(nproc, rank, nZones);
    size_t nz = p.second;
    size_t nzStart = p.first;
    std::cout << "Rank " << rank << " reads zones " << nzStart << ".." << nzStart + nz - 1 << std::endl;

    double tstart = MPI_Wtime();
    size_t nBytesRead = 0;

    adios2::Variable<int32_t> vStep = io.InquireVariable<int32_t>("/hpMusic_base/step");

    if (!rank)
    {
        std::cout << "Step = " << step << std::endl;
    }

    nElems.resize(nZones);
    nNodes.resize(nZones);
    ReadVariable<int32_t>(rank, zoneprefix + "/n_total_node", reader, io, step, nNodes.data());
    ReadVariable<int32_t>(rank, zoneprefix + "/n_total_elem", reader, io, step, nElems.data());
    reader.PerformGets();
    nBytesRead += nZones * 2 * sizeof(int32_t);

    CreateHDF5OutputStructure();

    for (size_t zone = nzStart; zone < nzStart + nz; ++zone)
    {
        size_t nConnectivitySize = (size_t)(ElementTypeStringToNodes(ElementType) * nElems[zone]);
        int64_t *dataEC = static_cast<int64_t *>(malloc(nConnectivitySize * sizeof(int64_t)));
        nBytesRead += nConnectivitySize * sizeof(int64_t);
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

        /* Read zone variables */
        ReadVariableBlock<int64_t>(rank, zoneprefix + "/Elem/ElementConnectivity", reader, io, step, zone, dataEC);
        ReadVariableBlock<int64_t>(rank, zoneprefix + "/Elem/ElementRange", reader, io, step, zone, dataER);
        ReadVariableBlock<double>(rank, zoneprefix + "/GridCoordinates/CoordinateX", reader, io, step, zone, gcx);
        ReadVariableBlock<double>(rank, zoneprefix + "/GridCoordinates/CoordinateY", reader, io, step, zone, gcy);
        ReadVariableBlock<double>(rank, zoneprefix + "/GridCoordinates/CoordinateZ", reader, io, step, zone, gcz);

        for (int i = 0; i < FlowVariables.size(); ++i)
        {
            ReadVariableBlock<double>(rank, FlowVariables[i], reader, io, step, zone, ptrs[i]);
        }

        reader.PerformGets();

        /* CGNS connectivity is 1-based, ADIOS data is 0 based. Fix it here */
        for (size_t k = 0; k < nConnectivitySize; ++k)
        {
            ++dataEC[k];
        }
        ++dataER[0];
        ++dataER[1];

        /* Write zone variables */
        WriteDataset(zoneHierarchy[zone].dsElementConnectivity_data, H5T_NATIVE_INT64, dataEC);
        WriteDataset(zoneHierarchy[zone].dsElementRange_data, H5T_NATIVE_INT64, dataER);
        WriteDataset(zoneHierarchy[zone].dsCoordinateX_data, H5T_NATIVE_DOUBLE, gcx);
        WriteDataset(zoneHierarchy[zone].dsCoordinateY_data, H5T_NATIVE_DOUBLE, gcy);
        WriteDataset(zoneHierarchy[zone].dsCoordinateZ_data, H5T_NATIVE_DOUBLE, gcz);
        for (int i = 0; i < FlowVariables.size(); ++i)
        {
            WriteDataset(zoneHierarchy[zone].dsFlowSolutionVariables_data[i], H5T_NATIVE_DOUBLE, ptrs[i]);
        }

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

    for (auto &bg : boundaryGroups)
    {
        const std::string &BoundaryName = bg.boundaryName;
        size_t nBoundaryZones = bg.boundaryHierarchy.size();

        // BoundaryInfo *bi = bg.boundaryInfo;
        if (!rank)
        {
            std::cout << "Process Boundary " << BoundaryName << ", with zones = " << nBoundaryZones << std::endl;
        }

        auto p = decomposeZones(nproc, rank, nBoundaryZones);
        size_t nz = p.second;
        size_t nzStart = p.first;
        std::cout << "Rank " << rank << " reads boundary zones " << nzStart << ".." << nzStart + nz - 1 << std::endl;

        for (size_t zone = nzStart; zone < nzStart + nz; ++zone)
        {
            auto &bh = bg.boundaryHierarchy[zone];

            std::string bname = zoneprefix + "/" + BoundaryName;
            auto attNBoundaries = io.InquireAttribute<int32_t>(bname + "/nBoundaries");
            size_t nConnectivitySize = bg.boundaryInfo->nConnectivity[zone];
            int64_t *dataEC = static_cast<int64_t *>(malloc(nConnectivitySize * sizeof(int64_t)));
            nBytesRead += nConnectivitySize * sizeof(int64_t);
            int64_t *dataER = static_cast<int64_t *>(malloc(2 * sizeof(int64_t)));
            nBytesRead += 2 * sizeof(int64_t);

            ReadVariableBlock<int64_t>(rank, bname + "/ElementConnectivity", reader, io, step, zone, dataEC);
            ReadVariableBlock<int64_t>(rank, bname + "/ElementRange", reader, io, step, zone, dataER);

            reader.PerformGets();

            /* CGNS connectivity is 1-based, ADIOS data is 0 based. Fix it here */
            for (size_t k = 0; k < nConnectivitySize; ++k)
            {
                ++dataEC[k];
            }
            ++dataER[0];
            ++dataER[1];

            WriteDataset(bh.dsElementConnectivity_data, H5T_NATIVE_INT64, dataEC);
            WriteDataset(bh.dsElementRange_data, H5T_NATIVE_INT64, dataER);

            free(dataEC);
            free(dataER);
        }
    }

    reader.Close();
    CloseHDF5();

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
