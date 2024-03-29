PROGNAME = genetics
SRCDIR = src/
OBJDIR = obj/
EXECUTS = $(PROGNAME) $(PROGNAME).debug $(PROGNAME).profile
HELPERS = lint.out check.out

SOURCES = $(wildcard $(SRCDIR)*.c)
DEPENDS	= $(wildcard $(SRCDIR)*.h)
OBJECTS = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(SOURCES:.c=.o))

CC = gcc
CFLAGS = -Wall -Wextra -Wshadow -std=c11 -pedantic -Ofast -fopenmp
LFLAGS = -lm
CDEBUGFLAGS = -D DEBUG -ggdb -g3 -O0
CPROFILEFLAGS = -pg

$(PROGNAME): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPENDS)
	@mkdir -pv $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: release
release: release.tar.gz

release.tar.gz: $(SOURCES) $(DEPENDS) Makefile README.rst
	tar -czvf $@ $^

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
	$(RM) $(EXECUTS)
	$(RM) $(HELPERS) gmon.out
	$(RM) release.tar.gz

.PHONY: bin
bin: $(EXECUTS)

.PHONY: help
help: $(HELPERS)

.PHONY: all
all: bin help

.PHONY: debug
debug: $(PROGNAME).debug

$(PROGNAME).debug: $(SOURCES) $(DEPENDS)
	$(CC) $(CFLAGS) $(CDEBUGFLAGS) -o $@ $(SOURCES) $(LFLAGS)

.PHONY: profile
profile: $(PROGNAME).profile

$(PROGNAME).profile: $(SOURCES) $(DEPENDS)
	$(CC) $(CFLAGS) $(CPROFILEFLAGS) -o $@ $(SOURCES) $(LFLAGS)

.PHONY: lint
lint: lint.out
lint.out: $(SOURCES) $(DEPENDS)
	splint -booltype Bool -boolfalse FALSE -booltrue TRUE -D__uint128_t=int -D__linux__ $^ > $@ || true

.PHONY: check
check: check.out
check.out: $(SOURCES) $(DEPENDS)
	cppcheck --std=c11 --enable=all --suppress=missingIncludeSystem $^ &> $@
