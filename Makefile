export CXX          = g++
export CXXFLAGS     = -Wno-write-strings -W -Wall -std=c++0x -O3 -lgsl -lgslcblas

#PROTBINFO_MODULE 	= common seq align msaanalysis
#PROTBINFO_SRC_DIR 	= $(addprefix protbinfo/src/, $(PROTBINFO_MODULE))

LBFGS_MODULE 	    = lbfgs
LBFGS_SRC_DIR 	    = $(LBFGS_MODULE)

SRC_DIR 			= src $(addprefix src/, seq util)
#DIRS 				= $(SRC_DIR) $(PROTBINFO_SRC_DIR) $(LBFGS_SRC_DIR)
DIRS 				= $(SRC_DIR) $(LBFGS_SRC_DIR)

PROG_MRFMAIN		= pmrf
PROG_MRFMAIN_OBJECT	= src/$(PROG_MRFMAIN).o

PROG_OBJECTS		= $(PROG_MRFMAIN_OBJECT)

SRCS   				= $(foreach dir,. $(DIRS), $(wildcard $(dir)/*.cpp))
OBJECTS				= $(filter-out $(PROG_OBJECTS), $(SRCS:.cpp=.o))

all : prog

prog : objs
	$(CXX) $(CXXFLAGS) -o $(PROG_MRFMAIN) $(PROG_MRFMAIN_OBJECT) $(OBJECTS)

objs :
	@for dir in $(DIRS); do\
		make -C $$dir || exit $?;\
	done

clean : test_clean prog_clean
	@for dir in $(DIRS); do\
		make -C $$dir clean;\
	done

prog_clean:
	rm -rf $(PROG_MRFMAIN)

###############
# Test runner #
###############

TESTFLAGS           = -lpthread -lgtest

TESTRUNNER 			= runtest
TESTRUNNER_DIR		= tests
TESTRUNNER_OBJECT 	= $(TESTRUNNER_DIR)/$(TESTRUNNER).o

#TEST_DIRS 			= $(addsuffix /tests, $(SRC_DIR)) $(TESTRUNNER_DIR) $(addsuffix /tests, $(PROTBINFO_SRC_DIR))
TEST_DIRS 			= $(addsuffix /tests, $(SRC_DIR)) $(TESTRUNNER_DIR)
TEST_SRCS			= $(foreach dir,. $(TEST_DIRS), $(wildcard $(dir)/*.cpp))
TEST_OBJECTS		= $(filter-out $(TESTRUNNER_OBJECT), $(TEST_SRCS:.cpp=.o))

test : objs test_objs
	$(CXX) $(CXXFLAGS) -o $(TESTRUNNER) $(TESTRUNNER_OBJECT) $(TEST_OBJECTS) $(OBJECTS) $(TESTFLAGS)
	./$(TESTRUNNER) --gtest_throw_on_failure

test_objs :
	@for dir in $(TEST_DIRS); do\
		make -C $$dir || exit $?;\
	done

test_clean :
	@for dir in $(TEST_DIRS); do\
		make -C $$dir clean;\
	done
	rm -rf $(TESTRUNNER)
