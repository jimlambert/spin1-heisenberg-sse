# Makefile for Black_Jack project
MAIN = HESSE_S1
CXX = g++ 
CXX_FLAGS = -std=c++14 -O0 -Wall -Wextra -g
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(addprefix obj/, $(notdir $(CPP_FILES:.cpp=.o)))
INCLUDE = -I./inc/

# This section makes the full black_jack executable

$(MAIN): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -o $@ $^

obj/%.o: src/%.cpp  
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

clean:	
	$(RM) $(OBJ_FILES) $(MAIN)


