CC6_GENERATE_LIBRARY(
    TARGET cc6_common
    SOURCES
        Sources/MultiPointSource.cc
        Sources/PlanarSource.cc
        Sources/PointSource.cc

        BinnedGeometry.cc
        DetPlane.cc
        DR_GenerallStructs.cc
        G4Input.cc
        InputReader.cc
        InputReaderEI.cc
        InputReaderGeant.cc
        InputReaderNN.cc
        InputReaderSimple.cc
        InputReaderPMI.cc
        IsectionPoint.cc
        Source.cc
        Track.cc

    HEADERS
        Sources/MultiPointSource.hh
        Sources/PlanarSource.hh
        Sources/PointSource.hh

        BinnedGeometry.hh
        DetPlane.hh
        DR_GenerallStructs.hh
        G4Input.hh
        InputReader.hh
        InputReaderEI.hh
        InputReaderGeant.hh
        InputReaderNN.hh
        InputReaderSimple.hh
        InputReaderPMI.hh
        IsectionPoint.hh
        Source.hh
        Track.hh

    LIBRARIES
        cc6_tools
        ROOT::Core
        ROOT::Physics
        ROOT::Geom
        ROOT::Tree
        ROOT::Gpad
)
