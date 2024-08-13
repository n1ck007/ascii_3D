SRC=donut
CC=clang++

${SRC}: ${SRC}.cpp
	${CC} -g -lm $< -o $@
