export CXX          = g++
export CXXFLAGS     = -Wno-write-strings -W -Wall -std=c++0x -O3 -lgsl -lgslcblas

MODULE 				= common seq align msaanalysis machlearn
SRC_DIR 			= $(addprefix src/, $(MODULE))
DIRS 				= $(SRC_DIR)

SRCS   				= $(foreach dir,. $(DIRS), $(wildcard $(dir)/*.cpp))
OBJECTS				= $(SRCS:.cpp=.o)

all : test

objs :
	@for dir in $(DIRS); do\
		make -C $$dir || exit $?;\
	done

clean : test_clean
	@for dir in $(DIRS); do\
		make -C $$dir clean;\
	done

###############
# Test runner #
###############

TESTFLAGS           = -lpthread -lgtest

TESTRUNNER 			= runtest
TESTRUNNER_DIR		= tests
TESTRUNNER_OBJECT 	= $(TESTRUNNER_DIR)/$(TESTRUNNER).o

TEST_DIRS 			= $(addsuffix /tests, $(SRC_DIR)) $(TESTRUNNER_DIR)
TEST_SRCS			= $(foreach dir,. $(TEST_DIRS), $(wildcard $(dir)/*.cpp))
TEST_OBJECTS		= $(filter-out $(TESTRUNNER_OBJECT), $(TEST_SRCS:.cpp=.o))

test : objs test_objs
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) -o $(TESTRUNNER) $(TESTRUNNER_OBJECT) $(TEST_OBJECTS) $(OBJECTS)
	./$(TESTRUNNER)

test_objs :
	@for dir in $(TEST_DIRS); do\
		make -C $$dir || exit $?;\
	done

test_clean :
	@for dir in $(TEST_DIRS); do\
		make -C $$dir clean;\
	done
	rm -rf $(TESTRUNNER)
