# Arabidopsis thaliana TE Transcripts
TE transcript-based annotation layered over Arabidopsis TAIR10 TE annotation

TAIR10 TE annotation provides 6 columns of information for each Transposable Element (TE):
1)	TE-ID
2)	If the orientation is 5’
3)	TE Start
4)	TE End
5)	TE Family
6)	TE superfamily and Type

This information is available in the first 9 columns of the new annotation bed file. 

A typical bed file format is followed where the first 6 columns are: Chromosome, Start, Stop, TE-ID, Score, Strand. This format is chosen for ease of data analysis using BEDTools. 

In addition, the following columns are noted for each TE:

7: Subfamily of the TE
8: Family of the TE
9: Super Family of the TE

10: Length of the TE
11: Length Category of the TE as described in Panda et. al., Genome Biology 2016.
12: Copy# of the TE subfamily.
13: Copy# category of the TE subfamily as described in Panda et. al., Genome Biology 2016.
14: Distance from the centromere
15: Position category based on distance from centromere as described in Panda et. al., Genome Biology 2016.
16: Distance from the nearest gene
17: Type of RNA-directed DNA methylation (RdDM) acting on the TE when TEs are silent (wt Col) as reported in Panda et. al., Genome Biology 2016.
18: Type of RdDM acting on the TE when TEs are active (ddm1) as reported in Panda et. al., Genome Biology 2016.

19: Evidence of Expression: 
	“.” if no expression found
	“Low Expressed” if at least 1 read reported expression but not enough to annotate transcripts
	“ExpressedAndAnnotated” if TE is expressed high enough for transcript annotation.
20: Transcriptional Start Site (TSS). Comma separated values note multiple TSSs.
21: Transcriptional Stop Site or more accurately polyadenylation site.
22: Transcriptional Strand. Comma separated valued not direction of each transcript.
23: Transcript ID. This transcript ID can be cross-referenced in the GFF file for investigating exon intron structure.

References
Panda, K., Ji, L., Neumann, D. A., Daron, J., Schmitz, R. J., & Slotkin, R. K. (2016). Full-length autonomous transposable elements are preferentially targeted by expression-dependent forms of RNA-directed DNA methylation. Genome Biology, 17(1), 170. http://doi.org/10.1186/s13059-016-1032-y
