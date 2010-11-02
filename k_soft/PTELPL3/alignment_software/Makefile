CXXFLAGS=-g -O2 -Wall -Wextra -W -O3 -DNDEBUG
CFLAGS=-g -O2 -Wall -Wextra -W -O3 -DNDEBUG
LDFLAGS= -lcvd -lGVars3
CXX=g++

TARGETS=composite ncc_coarse overlay_pair ncc_fine weighted_mean
ALL=$(TARGETS) hd_imview

targets:$(TARGETS)
all:$(ALL)

hd_imview: hd_imview.o
	$(CXX) -o $@ $<  $(LDFLAGS)

composite: composite.o
	$(CXX) -o $@ $<  $(LDFLAGS)

ncc_coarse: ncc_coarse.o
	$(CXX) -o $@ $<  $(LDFLAGS)

overlay_pair: overlay_pair.o
	$(CXX) -o $@ $<  $(LDFLAGS)

ncc_fine: ncc_fine.o
	$(CXX) -o $@ $<  $(LDFLAGS)

weighted_mean: weighted_mean.o
	$(CXX) -o $@ $<  $(LDFLAGS)

clean:
	rm -f *.o $(ALL)

distclean: clean
	rm -fr autom4te.cache  config.status config.log Makefile

.PHONY: all clean distclean
