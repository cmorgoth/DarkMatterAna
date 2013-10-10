CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

OS_NAME:=$(shell uname -s | tr A-Z a-z)
ifeq ($(OS_NAME),darwin)
STDINCDIR := -I/opt/local/include
STDLIBDIR := -L/opt/local/lib
else
STDINCDIR :=-I/media/data/cmorgoth/git/DarkMatterAna/
STDLIBDIR := 
endif

CPPFLAGS := -Wl,--no-as-needed $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g

TARGET = DM_ANA_Err

SRC = src/test.cc src/DM_StackPlots.cc src/DM_RatioPlots.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc src/DM_DY_HTBins.cc src/DM_TT_LSLH.cc src/DM_Data.cc src/DM_ZJetsNuNu.cc src/DM_WJetsHTBins.cc src/DM_2D_MR_RSQ_Dist.cc src/DM_T2CC.cc src/DM_METPlots.cc

#TARGET = bkg_pre_1mu
#SRC = bkg_pre_1mu.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc

OBJ = $(SRC:.cc=.o)

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) *~

