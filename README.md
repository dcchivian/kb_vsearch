[![Build Status](https://travis-ci.org/dylan/kb_vsearch.svg?branch=master)](https://travis-ci.org/dylan/kb_vsearch)

# kb_vsearch
---

VSearch will have several methods offered on KBase.

* Basic Search: one sequence against many
* Cluster: build groups of related sequences
* Chimera: identify sequences containing chimeras
* Dereplicate: remove redunant sequences

Basic Search Mode:

Initially, only SingleEndLibrary sequence sets will be operated on.  Subsequent versions will accept the following types for

  One sequence argument:
        read_seq_set
        feature
  Many sequences argument:
        read_seq_set
        feature_set
        genome
        genome_set
