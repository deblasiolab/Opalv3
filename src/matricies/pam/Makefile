SHELL = /bin/sh
CFLAGS =
CC = cc

all: pam 

pam: pam.o
	$(CC) -o $@ $^ -lm

.PHONY: clean veryclean

clean:
	-rm -f pam.o

veryclean: clean
	-rm -f pam
