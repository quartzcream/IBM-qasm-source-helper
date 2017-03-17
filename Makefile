CC = gcc
CCFLAGS = -Wall
LDFLAGS = -lm
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = analyze

debug: CCFLAGS += -g -O0
debug: all

devel: CCFLAGS += -g -O2
devel: all

release: CCFLAGS += -O2
release: all

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	$(CC) $(CCFLAGS) -c $<

%.o: %.c
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -f *.o $(TARGET)
