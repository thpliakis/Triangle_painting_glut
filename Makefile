CC=gcc
CFLAGS = -lglut -lGLU -lGL -lm

default: all

main_glut: main_glut.c
	$(CC)  -o $@ $^ $(CFLAGS)

random_triangles: random_triangles.c
	$(CC)  -o $@ $^ $(CFLAGS)

.PHONY: clean

all: main_glut
rt: random_triangles
	./random_triangles

test: main_glut
	./main_glut
	
clean:
	rm -f main_glut
