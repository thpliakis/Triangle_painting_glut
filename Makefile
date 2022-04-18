CC=gcc
CFLAGS = -lglut -lGLU -lGL

default: all

main_glut: main_glut.c
	$(CC)  -o $@ $^ $(CFLAGS)

.PHONY: clean

all: main_glut

test: main_glut
	./main_glut
	
clean:
	rm -f main_glut
