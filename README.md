# CACTI-6.0
ASIC simulation of Multi-ported Memory Module. And it can offer SRAM-based dual-port basic building block to support multiple read/write ports, and it can easily get the simulation result (like frequency, power consumption, and area... etc.). 

CACTI is an analytical tool that takes a set of cache para-
meters as input and calculates access time, power, cycle 
time and area of a cache.
CACTI was originally developed by Dr. Jouppi and Dr. Wilton
in 1993 and since then it has undergone five major 
revisions.

List of features (version 1-6):
===============================
The following are the list of features supported by the tool. 

* Power, delay, area, and cycle time model for 
                  direct mapped caches
                  set-associative caches
                  fully associative caches
                  
* Support for modeling multi-ported uniform cache access (UCA)
  and multi-banked, multi-ported non-uniform cache access (NUCA).

* Leakage power calculation that also considers the operating
  temperature of the cache.
  
* Router power model.

* Interconnect model with different delay, power, and area 
  properties including low-swing wire model.

* API to perform trade-off analysis involving power, delay,
  area, and bandwidth.

* All process specific values used by the tool are obtained
  from ITRS and currently, the tool supports 90nm, 65nm, 45nm, 
  and 32nm technology nodes.

How to use the tool?
====================
Prior versions of CACTI take input parameters such as cache
size, technology node etc. as a set of command line arguments. 
In addition to this method, CACTI 6 also lets users specify
their cache model in a more detailed manner by using a config
file (cache.cfg).

-> define the cache model using cache.cfg
-> run the "cacti" binary <./cacti>

As specified earlier, CACTI also supports command line arguments.
Please use "-h" argument <./cacti -h > to get more details.

List of files in CACTI:
=======================
def.h -> carries a list of process specific variables 

technology.c -> has ITRS values for different technology nodes

area.c -> defines area model for different cache modules

time.c -> defines power and delay model for different cache modules

leakage.c -> carries definition of basic primitives for leakage
             power calculation

basic_circuit.c -> has a set of building blocks required by
            the timing and power models

router.c -> defines the router power model

io.c -> carries different i/o functions.


For complete documentation of the tool, please refer CACTI-6
technical report and the following paper,
"Optimizing NUCA Organizations and Wiring Alternatives for 
Large Caches With CACTI 6.0" that appear in MICRO 07.
