# External Data Fit Makefile
#   Modified from the Makefile in T2KReWeight
#
# For help, type ./configure --help
#
# Authors:
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#  STFC - Rutherford Appleton Laboratory, UK
#
# Jim Dobson <j.dobson07 \at imperial.ac.uk>
#  Imperial College London, UK
#
# Patrick de Perio <pdeperio \at physics.utoronto.ca>
#  University of Toronto, Canada
#
# Callum Wilkinson <callum.wilkinson \at sheffield.ac.uk>
#   University of Sheffield, UK
#
# Patrick Stowell <patrick.stowell \at sheffield.ac.uk>
#   University of Sheffield, UK
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(EXT_FIT)/make/Make.include

BUILD_TARGETS = print \
	make-bin-lib-dir \
	src \
	app;

all: $(BUILD_TARGETS)

make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating externalDataFitter lib and bin directories..."
	cd ${EXT_FIT}; \
	[ -d bin ] || mkdir bin; chmod 755 bin;\
	[ -d lib ] || mkdir lib; chmod 755 lib;

print: FORCE
	@echo " "
	@echo " "
	@echo "***** Building externalDataFitter from source trees at: $(EXT_FIT)"
	@echo " "

src:	FORCE
	@echo " "
	@echo "** Building all in src directory..."
	cd ${EXT_FIT}/src;\
	make; \
	cd ${EXT_FIT}

app: FORCE
	@echo " "
	@echo "** Building all in app directory..."
	cd ${EXT_FIT}/app;\
	make; \
	cd ${EXT_FIT}

purge: FORCE
	@echo " "
	@echo "** Purging..."
	cd ${EXT_FIT}/src;\
	make purge;\
	cd ${EXT_FIT}

clean: clean-files clean-etc

clean-files: FORCE
	@echo " "
	@echo "** Cleaning..."
	cd ${EXT_FIT}/src; \
	make clean; cd ..; \
	cd ${EXT_FIT}/app; \
	make clean; cd ..;

clean-dir: FORCE
	@echo "Deleting the lib and bin directories...";\
	cd ${EXT_FIT};\
	[ ! -d ./bin ] || rm -r ./bin/*;\
	[ ! -d ./lib ] || rm -r ./lib/* 

clean-etc: FORCE
	cd ${EXT_FIT};\
	rm -f ./*log





FORCE:

