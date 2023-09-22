#include <nanovdb/util/IO.h>
#include <nanovdb/util/OpenToNanoVDB.h>
#include <openvdb/openvdb.h>

#include <iostream>
int main(int argc, char **argv) {
    openvdb::initialize();

    openvdb::io::File file(argv[1]);
    // Open the file.  This reads the file header, but not any grids.
    file.open();
    // Loop over all grids in the file and retrieve a shared pointer
    // to the one named "LevelSetSphere".  (This can also be done
    // more simply by calling file.readGrid("LevelSetSphere").)
    std::vector<nanovdb::GridHandle<nanovdb::HostBuffer>> handles;

    for (openvdb::io::File::NameIterator nameIter = file.beginName();
         nameIter != file.endName(); ++nameIter) {
        // Read in only the grid we are interested in.
        auto baseGrid = file.readGrid(nameIter.gridName());
        printf("Convert grid %s?\n", nameIter.gridName().c_str());

        char c;

        std::cin >> c;

        if (c == 'y') {
            auto floatGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
            auto handle = nanovdb::openToNanoVDB(*floatGrid);
            handles.emplace_back(std::move(handle));
        } else if (c == 'n')
            continue;
        else {
            printf("Default skip this grid!\n");
        }
    }
    nanovdb::io::writeGrids(argv[2], handles, nanovdb::io::Codec::NONE, 1);

    file.close();

    auto infos = nanovdb::io::readGridMetaData(argv[2]);

    for (const auto &info : infos) {
        printf("%s\n", info.gridName.c_str());
    }
}