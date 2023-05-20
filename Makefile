EXENAME = main
OBJS = main.o Matrix.o

CXX = clang++
CXXFLAGS = -std=c++14 -c -g -O0 -Wall -Wextra -pedantic
LD = clang++
LDFLAGS = -std=c++14 -lpthread -lm

all : $(EXENAME)

$(EXENAME) : $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXENAME)

main.o : main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

Matrix.o : Matrix.cpp Matrix.h
	$(CXX) $(CXXFLAGS) Matrix.cpp

clean :
	-rm -f *.o $(EXENAME) test
