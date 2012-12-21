CPP = gcc -Wall 
SRCS = main.c solver3d.c

all:
	$(CPP) $(SRCS) -o fluid -lGL -lGLU -lglut

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
