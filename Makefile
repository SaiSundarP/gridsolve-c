CC ?= gcc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS ?= -Iinclude
LDLIBS = -lm
SOURCES = src/common.c src/case.c src/ybus.c src/lu.c src/pf.c src/engine.c

.PHONY: all test clean sanitize

all: gridsolve

gridsolve: $(SOURCES) src/main.c include/gridsolve.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) src/main.c -o $@ $(LDLIBS)

tests/test_gridsolve: $(SOURCES) tests/test_gridsolve.c include/gridsolve.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) tests/test_gridsolve.c -o $@ $(LDLIBS)

test: tests/test_gridsolve
	./tests/test_gridsolve

sanitize:
	$(CC) $(CPPFLAGS) -std=c11 -O1 -g -Wall -Wextra -Wpedantic -Werror -fsanitize=address,undefined -fno-omit-frame-pointer $(SOURCES) tests/test_gridsolve.c -o tests/test_gridsolve_san $(LDLIBS)
	./tests/test_gridsolve_san

clean:
	rm -f gridsolve tests/test_gridsolve tests/test_gridsolve_san