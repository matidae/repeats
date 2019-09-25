## repeats
#### scripts for repeat analysis

**gff_merger.py** : script that takes a gff file, a cutoff and merges overlapping entries and those within the cutoff.

*e.g.: ./gff_merger.py annotation1.gff 10 source*


**repeat_analysis.py** : script that compares the gff of an annotation to a reference and outputs the percent annotated for each entry, a global confusion matrix and the Matthews correlation coefficient

*e.g.: ./repeat_analysis.py reference.gff repeats.gff*
