# -*- Makefile -*-

# Compiler and flags
CC = gcc
CFLAGS = -Wall -g

# Source and object directories
SRC_DIR = .
OBJ_DIR = o_files

SRCS = $(wildcard $(SRC_DIR)/*.c)

# List of object files
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

main: $(OBJS) h_files/nmdm.h
	$(CC) $(OBJS) -o main -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c h_files/nmdm.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f o_files/*.o
