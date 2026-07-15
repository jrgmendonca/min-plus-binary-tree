CC      ?= clang
CSTD    := -std=c11
WARN    := -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-function
# Tuned for Apple silicon; override on other hardware, e.g.:
#   make ARCH=-march=native
ARCH    ?= -mcpu=apple-m1
OPT     := -O3 -flto -fomit-frame-pointer -funroll-loops
CFLAGS  ?= $(CSTD) $(WARN) $(ARCH) $(OPT)
LDFLAGS ?= -flto
LDLIBS  := -lpthread

DBGFLAGS := $(CSTD) $(WARN) -O1 -g -fsanitize=address,undefined

all: minplus minplus_ref precompute_leaves

minplus: minplus.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

minplus_ref: minplus_ref.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

precompute_leaves: precompute_leaves.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lm

debug: minplus.c minplus_ref.c
	$(CC) $(DBGFLAGS) minplus.c     -o minplus_dbg     $(LDLIBS)
	$(CC) $(DBGFLAGS) minplus_ref.c -o minplus_ref_dbg

clean:
	rm -f minplus minplus_ref minplus_dbg minplus_ref_dbg precompute_leaves

.PHONY: all clean debug
