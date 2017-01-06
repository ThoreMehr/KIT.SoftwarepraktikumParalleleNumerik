CC ?= gcc

CFLAGS = -Wall -O3
LDFLAGS =

PROGRAM = do_lbm
SRCS = 		main.c \
		allocate.c \
		collide_and_stream.c \
		free_memory.c \
		init.c \
		set_u.c 

OBJS = $(SRCS:.c=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f $(PROGRAM) $(OBJS)

