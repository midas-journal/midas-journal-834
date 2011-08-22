INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

PROJECT(TrustRegionGradientDescentOptimizer)
ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages ITK-3.20.0 )

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)


#FIND_PACKAGE(ITK)
#IF(ITK_FOUND)
#  INCLUDE(${ITK_USE_FILE})
#ELSE(ITK_FOUND)
#  MESSAGE(FATAL_ERROR
#          "Cannot build without ITK.  Please set ITK_DIR.")
#ENDIF(ITK_FOUND)


SET(TRUST_SRC itkGradientDescentTrustRegionOptimizer.cxx itkGradientDescentTrustRegionOptimizer.h)

ADD_LIBRARY(gdtr ${TRUST_SRC})
TARGET_LINK_LIBRARIES(gdtr ${ITK_LIBRARIES} )

ADD_EXECUTABLE(itkGradientDescentTrustRegionOptimizerTest itkGradientDescentTrustRegionOptimizerTest.cxx )
TARGET_LINK_LIBRARIES(itkGradientDescentTrustRegionOptimizerTest ${ITK_LIBRARIES} gdtr )

ADD_TEST(GradientDescentTrustRegionOptimizerTest itkGradientDescentTrustRegionOptimizerTest)

ADD_EXECUTABLE(ImageRegistration1 ImageRegistration1.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration1 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration3 ImageRegistration3.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration3 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration4 ImageRegistration4.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration4 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration5 ImageRegistration5.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration5 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration6 ImageRegistration6.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration6 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration7 ImageRegistration7.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration7 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration9 ImageRegistration9.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration9 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration12 ImageRegistration12.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration12 ${ITK_LIBRARIES} gdtr)

ADD_EXECUTABLE(ImageRegistration13 ImageRegistration13.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration13 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration18 ImageRegistration18.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration18 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(ImageRegistration20 ImageRegistration20.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration20 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(MultiResImageRegistration1 MultiResImageRegistration1.cxx )
TARGET_LINK_LIBRARIES(MultiResImageRegistration1 ${ITK_LIBRARIES} gdtr )

ADD_EXECUTABLE(MultiResImageRegistration2 MultiResImageRegistration2.cxx )
TARGET_LINK_LIBRARIES(MultiResImageRegistration2 ${ITK_LIBRARIES} gdtr )


