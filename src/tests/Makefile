OBJECTS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
INCPATH     = -I ../ -I ../../thirdparty/

.SUFFIXES : .o .cpp
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<

all : $(OBJECTS)

clean :
	rm -f *.o 
