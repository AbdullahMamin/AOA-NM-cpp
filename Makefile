# == Only change these ======================
CC := g++
ERRFLAGS := -Wall -Wextra -Wpedantic -Werror
LIBS :=
OPT := -O3
CCFLAGS := -std=c++14
OUTPUT := out
# ===========================================

SOURCE_DIR := ./src

SOURCES := $(wildcard $(SOURCE_DIR)/*.cpp)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))
DEPENDS := $(patsubst %.cpp, %.d, $(SOURCES))

.PHONY: all clean install

all: $(OUTPUT) $(SHADER_OBJECTS)

clean:
	rm -rf $(OUTPUT) $(OBJECTS) $(DEPENDS) $(SHADER_OBJECTS)

$(OUTPUT): $(OBJECTS)
	$(CC) $(CCFLAGS) $(ERRFLAGS) $(OPT) $^ -o $@ $(LIBS)

install: $(OUTPUT)
	cp $(OUTPUT) /usr/bin/

-include $(DEPENDS)

%.o: %.cpp Makefile
	$(CC) $(CCFLAGS) $(ERRFLAGS) $(OPT) -MMD -MP -c $< -o $@ $(LIBS)
