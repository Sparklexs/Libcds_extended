
CPP=g++

#CPPFLAGS=-g3 -Wall -O0 -I../includes/ 
#CPPFLAGS=-O3 -Wall -DNDEBUG -I../includes/
INCS=-I../includes/
LIB=../lib/libcds.a

OBJECTS= testWM.o
BIN= testWM

%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(CPPFLAGS) $(INCS) -c $< -o $@

all: clean $(OBJECTS) $(BIN)
	@echo " [MSG] Done compiling tests"
	@echo " [FLG] $(CPPFLAGS)"

testWM:
	@echo " [LNK] Building testWM"
	@$(CPP) $(CPPFLAGS) -o testWM testWM.o $(LIB)

clean:
	@echo " [CLN] Cleaning object and binary files"
	@rm -f $(OBJECTS) $(BIN)
