# Makefile for building the third-party libraries.

# Define a destination directory (libs) if it is not already defined (as it should be)
3RD_DIR ?=$(CURDIR)

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders sdsl fasta

.PHONY: clean-prerequisites
clean-prerequisites: clean-pre

.PHONY: folders
folders:
	@echo "* Checking folders..." ; \
	[ -d $(3RD_DIR) ] || mkdir -p $(3RD_DIR) ; \
	[ -d $(3RD_DIR)/sdsl-lite/include ] || mkdir -p $(3RD_DIR)/sdsl-lite/include ; \
	[ -d $(3RD_DIR)/sdsl-lite/lib ] || mkdir -p $(3RD_DIR)/sdsl-lite/lib ; \
	[ -f $(3RD_DIR)/fasta.o ] || echo "* Fasta"

.PHONY: sdsl
sdsl: folders $(3RD_DIR)/sdsl-lite/include/divsufsort64.h

SDSL_NAME:=sdsl-lite
$(3RD_DIR)/sdsl-lite/include/divsufsort64.h:
	@echo "* SDSL-lite..." ;	\
	cd $(3RD_DIR)/../$(SDSL_NAME)/ &&	\
	./install.sh $(3RD_DIR)/$(SDSL_NAME)/;	\

.PHONY: fasta
fasta: folders $(3RD_DIR)/fasta.o

FASTA_NAME:=sparseMEM/src/
$(3RD_DIR)/fasta.o:
	@echo "* fasta.cpp..." ;	\
	cd $(3RD_DIR)/../$(FASTA_NAME) &&	\
	g++ -std=c++11 -Wall -c fasta.cpp -o $(3RD_DIR)/fasta.o;	\

.PHONY: clean-pre
clean-pre:
	@echo "* Cleaning prerequisites..." ; \
	rm -rf $(3RD_DIR)/fasta.o \
	rm -rf $(3RD_DIR)/$(SDSL_NAME) \
