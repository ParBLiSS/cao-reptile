LIB_DIR = /home/srirampc/lib
INC_DIR = /home/srirampc/include
#MY_CXX_FLAGS= -std=c++11 -Wall -Wextra  -DNDEBUG $(CODE_COVER)
MY_CXX_FLAGS= -std=c++11 -Wall -Wextra  $(CODE_COVER)
MY_CXX_OPT_FLAGS= -O3 -ffast-math -funroll-loops
MY_CXX=/usr/bin/c++
MY_CC=/usr/bin/cc
CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)
#CCLIB=-lsdsl -ldivsufsort -ldivsufsort64
CCLIB=
SOURCES=$(wildcard *.cpp)
EXECS=$(SOURCES:.cpp=.x)

all: $(EXECS)

build-test: $(EXECS)

%.x:%.cpp
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(CCLIB)
clean:
	rm -f $(EXECS)
	rm -rf *.dSYM
