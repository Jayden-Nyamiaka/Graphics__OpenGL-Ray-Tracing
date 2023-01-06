CC 		:= g++
FLAGS 	:= -O2 -Wall -g -std=c++14
LDLIBS 	:= -lpng -lyaml-cpp -lGLEW -lGL -lGLU -lglut -lm -lpthread
INCLUDE := -I./include -I/usr/X11R6/include -I/usr/include/GL -I/usr/include/

CCFLAGS := $(FLAGS) $(INCLUDE) -c $(LDLIBS)
LDFLAGS := $(FLAGS)

BIN_PATH := bin
OBJ_PATH := obj
SRC_PATH := src

TARGET_NAME := renderer
TARGET 		:= $(BIN_PATH)/$(TARGET_NAME)

SRC := $(foreach x, $(SRC_PATH), $(wildcard $(addprefix $(x)/*,.c*)))
OBJ := $(addprefix $(OBJ_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC)))))

CLEAN_LIST := $(TARGET) $(OBJ)

default: makedir all

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c*
	$(CC) $(CCFLAGS) -o $@ $< 

# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(BIN_PATH) $(OBJ_PATH)

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(CLEAN_LIST)