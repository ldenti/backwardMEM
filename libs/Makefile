# Makefile for building the third-party libraries.

# Define a destination directory (local) if it is not already defined (as it should be)
LOC_DIR ?=$(CURDIR)/../local

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders sdsl fasta

.PHONY: clean-prerequisites
clean-prerequisites: clean-sdsl clean-fasta

.PHONY: folders
folders:
	@echo "* Checking folders ..." ; \
	[ -d $(LOC_DIR) ] || mkdir -p $(LOC_DIR) ; \
	[ -d $(LOC_DIR)/include ] || mkdir -p $(LOC_DIR)/include ; \
	[ -d $(LOC_DIR)/lib ] || mkdir -p $(LOC_DIR)/lib ; \
	[ -d $(LOC_DIR)/sparseMEM ] || mkdir -p $(LOC_DIR)/sparseMEM ; \
	[ -f $(LOC_DIR)/sparseMEM/fasta.o ] || echo "* Fasta"

.PHONY: sdsl
sdsl: folders $(LOC_DIR)/include/divsufsort64.h

SDSL_NAME:=sdsl-lite
$(LOC_DIR)/include/divsufsort64.h:
	@echo "* SDSL-lite ..." ;	\
	git clone https://github.com/simongog/sdsl-lite.git;	\
	if [ -d $(SDSL_NAME)/ ]; then	\
		cd $(SDSL_NAME)/ &&	\
		./install.sh $(LOC_DIR)/;	\
	else \
		echo "SDSL library not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-sdsl
clean-sdsl:
	@echo "* Cleaning SDSL library ..." ; \
	rm -rf $(SDSL_NAME)/ \
	rm -rf $(LOC_DIR)/include \
	rm -rf $(LOC_DIR)/lib \

.PHONY: fasta
fasta: folders $(LOC_DIR)/sparseMEM/fasta.o

FASTA_NAME:=sparseMEM/src/
$(LOC_DIR)/sparseMEM/fasta.o:
	@echo "* fasta.cpp ..." ;	\
	git clone https://github.com/zia1138/sparseMEM.git;	\
	if [ -d $(FASTA_NAME)/ ]; then	\
		cd $(FASTA_NAME)/ &&	\
		g++ -std=c++11 -Wall -c fasta.cpp -o $(LOC_DIR)/sparseMEM/fasta.o;	\
	else \
		echo "sparseMEM not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-fasta
clean-fasta:
	@echo "* Cleaning fasta.cpp ..." ; \
	rm -rf $(LOC_DIR)/sparseMEM \
	rm -rf $(CURDIR)/sparseMEM \
	rm -rf $(LOC_DIR)
