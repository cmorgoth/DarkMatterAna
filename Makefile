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

#SRC = src/testV2.cc src/DM_BaseV2.cc src/DM_DataV2.cc 
SRC = src/testV3.cc src/DM_BaseV2.cc src/DM_WJetsHTBinsV3.cc src/DM_ZJetsNuNuV3.cc src/DM_DY_HTBinsV3.cc src/DM_TT_LSLHV3.cc

#SRC = src/testV3.cc src/DM_StackPlots.cc src/DM_RatioPlots.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc src/DM_DY_HTBins.cc src/DM_TT_LSLH.cc src/DM_Data.cc src/DM_ZJetsNuNu.cc src/DM_WJetsHTBinsV3.cc src/DM_2D_MR_RSQ_Dist.cc src/DM_METPlots.cc src/DM_BaseV2.cc

#SRC = src/error_study.cc src/DM_StackPlots.cc src/DM_RatioPlots.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc src/DM_DY_HTBins.cc src/DM_TT_LSLH.cc src/DM_Data.cc src/DM_ZJetsNuNu.cc src/DM_WJetsHTBins.cc src/DM_2D_MR_RSQ_Dist.cc src/DM_T2CC.cc src/DM_METPlots.cc

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

