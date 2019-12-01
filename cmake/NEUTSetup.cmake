# Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

################################################################################
#    This file is part of NUISANCE.
#
#    NUISANCE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NUISANCE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

include(cmake/parseConfigApp.cmake)

find_program(NEUTCONFIG NAMES neut-config)

LIST(APPEND EXTRA_CXX_FLAGS -DNEED_FILL_NEUT_COMMONS)

SET(HAVENEUTCONFIG FALSE)
# We are dealing with shiny NEUT
if(NOT "${NEUTCONFIG}" STREQUAL "NEUTCONFIG-NOTFOUND")
  SET(HAVENEUTCONFIG TRUE)
  cmessage(STATUS "Found neut-config, using it to determine configuration.")
else()
  cmessage(STATUS "Failed to find neut-config, assuming older NEUT build.")
endif()

if(HAVENEUTCONFIG)
  execute_process (COMMAND neut-config
    --version OUTPUT_VARIABLE NEUT_VERSION
             OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process (COMMAND neut-config
    --incdir OUTPUT_VARIABLE NEUT_INCLUDE_DIRS
             OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process (COMMAND neut-config
    --libdir OUTPUT_VARIABLE NEUT_LINK_DIRS
             OUTPUT_STRIP_TRAILING_WHITESPACE)

  GETLIBDIRS(neut-config --cernflags CERN_LIB_DIR)
  LIST(APPEND NEUT_LINK_DIRS ${CERN_LIB_DIR})
  GETLIBS(neut-config --cernflags CERN_LIBS)

  if(USE_REWEIGHT)
    execute_process (COMMAND neut-config
      --rwlibflags OUTPUT_VARIABLE NEUT_RWLIBS
               OUTPUT_STRIP_TRAILING_WHITESPACE)
    GETLIBS(neut-config --rwlibflags NEUT_RWLIBS)
    LIST(APPEND NEUT_LIBS ${NEUT_RWLIBS})
  else()
      GETLIBS(neut-config --libflags NEUT_GENLIBS)
      GETLIBS(neut-config --iolibflags NEUT_LIBS)
      LIST(APPEND NEUT_LIBS ${NEUT_IOLIBS})
      LIST(APPEND NEUT_LIBS ${NEUT_GENLIBS})
  endif()

  LIST(APPEND NEUT_LIBS ${CERN_LIBS};gfortran)
  LIST(APPEND EXTRA_LIBS ${NEUT_LIBS})

  PrefixList(NEUT_INCLUDE_DIRS "-I" ${NEUT_INCLUDE_DIRS})
  LIST(APPEND EXTRA_CXX_FLAGS ${NEUT_INCLUDE_DIRS} -D__NEUT_ENABLED__ -DNEUT_VERSION=${NEUT_VERSION})

  LIST(APPEND EXTRA_LINK_DIRS ${NEUT_LINK_DIRS})

  cmessage(STATUS "NEUT")
  cmessage(STATUS "     Version   : ${NEUT_VERSION}")
  cmessage(STATUS "     Flags     : ${NEUT_CXX_FLAGS}")
  cmessage(STATUS "     Includes  : ${NEUT_INCLUDE_DIRS}")
  cmessage(STATUS "     Link Dirs : ${NEUT_LINK_DIRS}")
  cmessage(STATUS "     Libs      : ${NEUT_LIBS}")


else() # Everything better be set up already

  if(NEUT_ROOT STREQUAL "")
    cmessage(FATAL_ERROR "Variable NEUT_ROOT is not defined. Please export environment variable NEUT_ROOT or configure with -DNEUT_ROOT=/path/to/NEUT. This must be set to point to a prebuilt NEUT instance.")
  endif()

  if(CERN STREQUAL "")
    cmessage(FATAL_ERROR "Variable CERN is not defined. Please export environment variable CERN or configure with -DCERN=/path/to/CERNLIB. This must be set to point to a prebuilt CERNLIB instance.")
  endif()

  if(CERN_LEVEL STREQUAL "")
      cmessage(FATAL_ERROR "Variable CERN_LEVEL is not defined. Please export environment variable CERN_LEVEL or configure with -DCERN_LEVEL=XXXX (likely to be 2005).")
  endif()

  if(${NEUT_VERSION} VERSION_LESS 5.4.0)
    set(NEUT_LIB_DIR ${NEUT_ROOT}/lib/Linux_pc)
  else()
    set(NEUT_LIB_DIR ${NEUT_ROOT}/lib)
  endif()

  set(NEUT_CLASS ${NEUT_ROOT}/src/neutclass)

  LIST(APPEND EXTRA_CXX_FLAGS -D__NEUT_ENABLED__ -DNEUT_VERSION=${NEUT_VERSION})

  LIST(APPEND EXTRA_CXX_FLAGS
    -I${NEUT_ROOT}/include
    -I${NEUT_ROOT}/src/neutclass)
  LIST(APPEND EXTRA_LINK_DIRS
    ${NEUT_LIB_DIR}
    ${CERN}/${CERN_LEVEL}/lib)

  if(USE_REWEIGHT)
    LIST(APPEND EXTRA_CXX_FLAGS
      -I${NEUT_ROOT}/src/reweight)
    LIST(APPEND EXTRA_LINK_DIRS
      ${NEUT_ROOT}/src/reweight)
  endif()

  if(${NEUT_VERSION} VERSION_EQUAL 5.4.2)
    LIST(APPEND EXTRA_LIBS
      -Wl,--as-needed)
    if(USE_REWEIGHT)
      LIST(APPEND EXTRA_LIBS
        NReWeight)
    endif()
    LIST(APPEND EXTRA_LIBS
      -Wl,--start-group
      neutcore_5.4.2
      nuccorspl_5.4.2 #typo in NEUT, may hopefully disappear
      nuceff_5.4.2
      partnuck_5.4.2
      skmcsvc_5.4.2
      tauola_5.4.2
      HT2p2h_5.4.0
      N1p1h_5.4.0
      -Wl,--end-group
      jetset74
      pdflib804
      mathlib
      packlib
      pawlib)

    LIST(APPEND EXTRA_CXX_FLAGS -DNEUT_COMMON_QEAV)
  elseif(${NEUT_VERSION} VERSION_EQUAL 5.4.0)
    LIST(APPEND EXTRA_LIBS
      -Wl,--as-needed)
    if(USE_REWEIGHT)
      LIST(APPEND EXTRA_LIBS
        NReWeight)
    endif()
    LIST(APPEND EXTRA_LIBS
      -Wl,--start-group
      neutcore_5.4.0
      nuccorspl_5.4.0 #typo in NEUT, may hopefully disappear
      nuceff_5.4.0
      partnuck_5.4.0
      skmcsvc_5.4.0
      tauola_5.4.0
      HT2p2h_5.4.0
      N1p1h_5.4.0
      specfunc_5.4.0
      radcorr_5.4.0
      gfortran
      -Wl,--end-group
      jetset74
      pdflib804
      mathlib
      packlib
      pawlib)
cmessage(STATUS "set the right neut 5.4.0")
  else()
    LIST(APPEND EXTRA_LIBS
      -Wl,--as-needed)
    if(USE_REWEIGHT)
      LIST(APPEND EXTRA_LIBS
        NReWeight)
    endif()
    LIST(APPEND EXTRA_LIBS
      -Wl,--start-group
      neutcore
      nuccorrspl
      nuceff
      partnuck
      skmcsvc
      tauola
      -Wl,--end-group
      jetset74
      pdflib804
      mathlib
      packlib
      pawlib)
  endif()

  set(NEUT_ROOT_LIBS)

  LIST(APPEND NEUT_ROOT_LIBS
    ${NEUT_CLASS}/neutctrl.so
    ${NEUT_CLASS}/neutfsivert.so)

  # Check for new versions of NEUT with NUCLEON FSI
  if(EXISTS "${NEUT_CLASS}/neutnucfsistep.so")
    set(NEUT_NUCFSI 1)
    LIST(APPEND EXTRA_CXX_FLAGS -DNEUT_NUCFSI_ENABLED)

    LIST(APPEND NEUT_ROOT_LIBS
      ${NEUT_CLASS}/neutnucfsistep.so
      ${NEUT_CLASS}/neutnucfsivert.so
      )
  endif()

  if(${NEUT_VERSION} VERSION_LESS 5.4.0)
    LIST(APPEND NEUT_ROOT_LIBS
      ${NEUT_CLASS}/neutrootTreeSingleton.so)
  endif()

  LIST(APPEND NEUT_ROOT_LIBS
    ${NEUT_CLASS}/neutvtx.so
    ${NEUT_CLASS}/neutfsipart.so
    ${NEUT_CLASS}/neutpart.so
    ${NEUT_CLASS}/neutvect.so
    )

  foreach(OBJ ${NEUT_ROOT_LIBS})
    LIST(APPEND EXTRA_SHAREDOBJS ${OBJ})
  endforeach()
endif()
