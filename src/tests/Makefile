OBJECTS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
INCPATH     = -I ../ -I ../../ -I ../../protbinfo/src/

.SUFFIXES : .o .cpp
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<

all : $(OBJECTS)

clean :
	rm -rf *.o 
