##
# source directory
##
SRC_DIR := src

##
# output directory
##
OUT_DIR := bin

##
# sources
##
SRCS := $(wildcard $(SRC_DIR)/facet/*.java $(SRC_DIR)/opal/*.java $(SRC_DIR)/opal/*/*.java $(SRC_DIR)/opal/*/*/*.java $(SRC_DIR)/com/traviswheeler/libs/*java $(SRC_DIR)/com/bluemarsh/graphmaker/core/util/*java $(SRC_DIR)/gnu/getopt/*java)
MBS := $(wildcard $(SRC_DIR)/gnu/getopt/MessagesBundle* $(SRC_DIR)/opal/IO/usage.txt)

##
# classes
## 
CLS := $(SRCS:$(SRC_DIR)/%.java=$(OUT_DIR)/%.class)
MBO := $(MBS:$(SRC_DIR)/%=$(OUT_DIR)/%)

##
# compiler and compiler flags
##
JC := javac
JCFLAGS := -d $(OUT_DIR)/ -cp $(SRC_DIR)/

##
# suffixes
##
.SUFFIXES: .java

##
# targets that do not produce output files
##
.PHONY: all clean

##
# default target(s)
##
all: $(CLS) $(MBO) $(OUT_DIR)/Opal.jar

$(CLS): $(OUT_DIR)/%.class: $(SRC_DIR)/%.java
	$(JC) $(JCFLAGS) $<

$(MBO): $(OUT_DIR)/%: $(SRC_DIR)/%
	cp $< $@

$(OUT_DIR)/Opal.jar: $(CLS)
	mkdir -p $(OUT_DIR)/META-INF
	./createManifest.pl $(OUT_DIR)
	cd $(OUT_DIR) ; jar -c -f Opal.jar -m META-INF/MANIFEST.MF *

##
# clean up any output files
##
clean:
	rm -rf $(OUT_DIR)/*
