## repeats
#### scripts for repeat analysis

**gff_merger.py** : script that takes a gff file, a cutoff and merges overlapping entries and those within the cutoff.

*e.g.: ./gff_merger.py annotation1.gff 10 source*


**repeat_analysis.py** : script that compares the gff of an annotation to a reference and outputs the percent annotated for each entry, a global confusion matrix and the Matthews correlation coefficient

*e.g.: ./repeat_analysis.py reference.gff repeats.gff*

---

#### scripts for simulating sequences with TEs
**random_sequence_TEs.py** : script needs config.yml and the files specified there (repeats_list and repeats.fa). Outputs a fasta sequences wih random TEs specified in repeat_list and a GFF file.

**random_nest_TEs.py** : script needs config.yml and the files specified there (repeats_list and repeats.fa). Outputs a fasta sequences with nested TEs as specified in repeat_list and an updated GFF file. Needs to be run after random_sequences_TEs.py, uses its output.

**repeats_list** : configuration file with a row for each TE and 9 columns: id_of_TE, number_of_repeats, %_identity, std_dev, %_indels, has_TSD_?, real_length, %_fragmentation, %_nested 

**config.yml** : configuration file with prefixes for the output, size of base sequence, and location of repeats_list and fasta files.
