OBJECTS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

.SUFFIXES : .o .cpp
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all : $(OBJECTS)

clean :
	rm -f *.o 
