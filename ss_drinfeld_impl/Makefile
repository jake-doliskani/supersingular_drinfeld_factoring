BUILD_DIR = build

SOURCES = BaseChange.cpp\
HasseLift.cpp\
HasseLiftExt.cpp\
MultiComposeMod.cpp\
FrobComp.cpp\
BalancedMul.cpp\
MultipointEval.cpp\
Util.cpp\
SSFactoring.cpp

OBJECTS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

GXX = g++
GXXFLAGS = -Wall -O3 -march=native -std=gnu++11
LIBS = -lntl -lgmp -lpthread

.PHONY = clean execs

all:  $(OBJECTS) execs

clean:
	rm -rf $(BUILD_DIR)
	rm -f test
	
$(BUILD_DIR)/%.o:	%.cpp
	$(GXX) -c $< -o $@ $(GXXFLAGS)

execs:  $(OBJECTS)
	$(GXX) test.cpp -o test $(OBJECTS) $(GXXFLAGS) $(LIBS)

# make the build dir if not exist
$(OBJECTS) : | $(BUILD_DIR)
$(BUILD_DIR) :
	mkdir $(BUILD_DIR)