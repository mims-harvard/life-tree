#
#	Makefile for all programs in this directory
#
#	implements:
#		all (default), create, analyze, graphlet, clean
#

# Global configuration for SNAP makefiles
GLIB = glib-core
SNAP = snap-core
GLIBADV = glib-adv
SNAPADV = snap-adv
SNAPEXP = snap-exp

CGLIB = Snap-4.0/$(GLIB)
CSNAP = Snap-4.0/$(SNAP)

EXGLIB = Snap-4.0/$(GLIB)
EXSNAP = Snap-4.0/$(SNAP)
EXGLIBADV = Snap-4.0/$(GLIBADV)
EXSNAPADV = Snap-4.0/$(SNAPADV)
EXSNAPEXP = Snap-4.0/$(SNAPEXP)

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC = g++
	CXXFLAGS += -std=c++0x -Wall
	CXXFLAGS += -O3 -DNDEBUG -fopenmp
	LDFLAGS +=
	LIBS += -lrt -lm
endif

ANALYZE = analyze
GRAPHLET = graphlet

all: $(ANALYZE) $(GRAPHLET)

# COMPILE

$(ANALYZE): $(ANALYZE).cpp $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(ANALYZE) $(ANALYZE).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(GRAPHLET): $(GRAPHLET).cpp $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(GRAPHLET) $(GRAPHLET).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(EXSNAP)/Snap.o: 
	make -C $(EXSNAP)

clean:
	rm -f *.o  $(ANALYZE) $(GRAPHLET)
