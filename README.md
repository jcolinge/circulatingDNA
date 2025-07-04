# circulatingDNA

This repository contains various programs that we developed for
our manuscript *Circulating DNA reveals nucleosome occupancy patterns
that are associated with nucleosome-DNA affinity and are affected in cancer*.
The different codes are found in folders corresponding to the main steps of 
this study.

There is a folder ``referenceFile`` with parameter files and
a large text file ``masterScript_v3`` that illustrates how these codes were 
executed (in our setting).

In a few cases, some programs were used in several parts of our work. We put 
such programs in one folder only.

The Savitzky-Golay filter code (``sgfilter.c`` and ``sgfilter.h``) in the ``WPNAcomputation`` folder
is a slightly modified version of the excellent implementation proposed by Fredrik
Jonsson. Instructions to compile the code are in the ``sgfilter.pdf`` file.

The C++ codes are compiled with very basic options, for instance

``g++ -O2 -o normalize normalize-1.cpp``

