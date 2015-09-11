OBJECTS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
INCPATH     = -I ../ -I ../../

.SUFFIXES : .o .cpp
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<

all : $(OBJECTS)

clean :
	rm -rf *.o 
