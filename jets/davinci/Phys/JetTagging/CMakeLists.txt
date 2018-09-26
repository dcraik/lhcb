################################################################################
# Package: JetTagging
################################################################################
gaudi_subdir(JetTagging v1r10)

gaudi_depends_on_subdirs(Phys/LoKiJets
			 Phys/LoKiAlgo
			 Kernel/LHCbMath
			 Phys/JetAccessories
			 Det/VeloDet
			 Phys/DaVinciInterfaces
			 )

find_package(NeuroBayesExpert)
find_package(ROOT COMPONENTS TMVA Physics)

gaudi_install_headers(Kernel)

find_package(Boost)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

gaudi_add_module(JetTagging
                 src/*.cpp
                 INCLUDE_DIRS NeuroBayesExpert
                 LINK_LIBRARIES LoKiAlgo LHCbMathLib VeloDetLib DaVinciInterfacesLib NeuroBayesExpert ROOT)

gaudi_add_dictionary(JetTagging
                     dict/JetTaggingInterface.h
                     dict/JetTaggingInterface.xml
                     LINK_LIBRARIES GaudiKernel TrackEvent RecEvent PhysEvent
                     OPTIONS " -U__MINGW32__ ")

gaudi_install_python_modules()

