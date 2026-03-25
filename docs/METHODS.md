# HelixHarbor Methods

## HelixHarbor: transmembrane helix and topology analysis

HelixHarbor is a tool for analyzing transmembrane helices (TMHs) and associated topological regions in transmembrane proteins. The tool supports:

1. single amino-acid sequence analysis
2. comparison of a user-supplied UniProt accession list against curated background sets
3. comparison of two UniProt accession lists
4. position-specific amino-acid composition (AAC) analysis of TMHs

## TMH/topology definition and sequence extraction

For single-sequence analysis, HelixHarbor uses TMbed to assign per-residue classes such as transmembrane helix or strand, signal peptide, or other topology states [Bernhofer2022TMbed].

Contiguous runs of the same class are collapsed into segments and reported with:

- 1-based start coordinates
- 1-based end coordinates
- segment type
- segment length

For analyses based on UniProt accessions, HelixHarbor retrieves protein sequences and annotated features from UniProt records [UniProtConsortium2025].

Segment sequences are extracted using the annotated coordinate ranges, including transmembrane and topological domain features. Where available, experimentally derived topology constraints can be incorporated via TOPDB-derived annotations [Tusnady2008TOPDB].

## Physicochemical features

For each TMH segment, HelixHarbor computes per-segment averages, defined as mean per-residue values across the segment.

### Surface area

Surface area is derived from the maximum allowed solvent accessibility reference values reported by Tien et al. [Tien2013ASA], reported in Å² per residue, and averaged across residues in the segment.

### Bulkiness

Bulkiness is calculated from the Zimmerman bulkiness index [Zimmerman1968] and averaged across residues in the segment.

### Hydrophobicity

Hydrophobicity is computed using the Kyte-Doolittle hydropathy scale [KyteDoolittle1982] and averaged across residues in the segment.

### Custom scales

The tool also supports user-supplied custom residue scales as two-column amino-acid-to-value tables. These are summarized as the mean value per residue per segment.

## Comparisons, statistics, and visualization

Distributions of per-helix feature values are visualized as kernel density estimates (KDEs). KDE curves are probability densities whose integrals equal 1. Each data point corresponds to the mean per-residue score across one transmembrane helix.

Group differences can be assessed using standard two-sample tests:

- Kolmogorov-Smirnov test
- Welch's t-test
- Mann-Whitney U test

Scientific computing and plotting are performed in the SciPy ecosystem [Virtanen2020SciPy], with figures rendered using Matplotlib [Hunter2007Matplotlib] and, where applicable, Seaborn [Waskom2021Seaborn].

## References

- [Bernhofer2022TMbed] Bernhofer M, Rost B. TMbed: transmembrane proteins predicted through language model embeddings. *BMC Bioinformatics*. 2022;23:326. doi:10.1186/s12859-022-04873-x
- [UniProtConsortium2025] The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Research*. 2025;53(D1):D609-D617. doi:10.1093/nar/gkae1010
- [Tusnady2008TOPDB] Tusnády GE, Kalmár L, Simon I. TOPDB: topology data bank of transmembrane proteins. *Nucleic Acids Research*. 2008;36(suppl_1):D234-D239. doi:10.1093/nar/gkm751
- [Tien2013ASA] Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum Allowed Solvent Accessibilities of Residues in Proteins. *PLOS ONE*. 2013;8(11):e80635. doi:10.1371/journal.pone.0080635
- [Zimmerman1968] Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences in proteins by statistical methods. *Journal of Theoretical Biology*. 1968;21(2):170-201. doi:10.1016/0022-5193(68)90069-6
- [KyteDoolittle1982] Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*. 1982;157(1):105-132. doi:10.1016/0022-2836(82)90515-0
- [Virtanen2020SciPy] Virtanen P, et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nature Methods*. 2020;17:261-272. doi:10.1038/s41592-019-0686-2
- [Hunter2007Matplotlib] Hunter JD. Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*. 2007;9(3):90-95. doi:10.1109/MCSE.2007.55
- [Waskom2021Seaborn] Waskom ML. seaborn: statistical data visualization. *Journal of Open Source Software*. 2021;6(60):3021. doi:10.21105/joss.03021

## Appendix: BibTeX entries

```bibtex
@article{Bernhofer2022TMbed,
  author  = {Bernhofer, Michael and Rost, Burkhard},
  title   = {{TMbed}: transmembrane proteins predicted through language model embeddings},
  journal = {BMC Bioinformatics},
  year    = {2022},
  volume  = {23},
  pages   = {326},
  doi     = {10.1186/s12859-022-04873-x}
}

@article{UniProtConsortium2025,
  author  = {{The UniProt Consortium}},
  title   = {{UniProt}: the Universal Protein Knowledgebase in 2025},
  journal = {Nucleic Acids Research},
  year    = {2025},
  volume  = {53},
  number  = {D1},
  pages   = {D609--D617},
  doi     = {10.1093/nar/gkae1010}
}

@article{Tusnady2008TOPDB,
  author  = {Tusn{\'a}dy, G{\'a}bor E. and Kalm{\'a}r, Lajos and Simon, Istv{\'a}n},
  title   = {{TOPDB}: topology data bank of transmembrane proteins},
  journal = {Nucleic Acids Research},
  year    = {2008},
  volume  = {36},
  number  = {suppl_1},
  pages   = {D234--D239},
  doi     = {10.1093/nar/gkm751}
}

@article{Tien2013ASA,
  author  = {Tien, Matthew Z. and Meyer, Austin G. and Sydykova, Dariya K. and Spielman, Stephanie J. and Wilke, Claus O.},
  title   = {Maximum Allowed Solvent Accessibilites of Residues in Proteins},
  journal = {PLOS ONE},
  year    = {2013},
  volume  = {8},
  number  = {11},
  pages   = {e80635},
  doi     = {10.1371/journal.pone.0080635}
}

@article{Zimmerman1968,
  author  = {Zimmerman, J. M. and Eliezer, Naomi and Simha, R.},
  title   = {The characterization of amino acid sequences in proteins by statistical methods},
  journal = {Journal of Theoretical Biology},
  year    = {1968},
  volume  = {21},
  number  = {2},
  pages   = {170--201},
  doi     = {10.1016/0022-5193(68)90069-6}
}

@article{KyteDoolittle1982,
  author  = {Kyte, Jack and Doolittle, Russell F.},
  title   = {A simple method for displaying the hydropathic character of a protein},
  journal = {Journal of Molecular Biology},
  year    = {1982},
  volume  = {157},
  number  = {1},
  pages   = {105--132},
  doi     = {10.1016/0022-2836(82)90515-0}
}

@article{Virtanen2020SciPy,
  author  = {Virtanen, Pauli and {et al.}},
  title   = {{SciPy} 1.0: Fundamental Algorithms for Scientific Computing in Python},
  journal = {Nature Methods},
  year    = {2020},
  volume  = {17},
  pages   = {261--272},
  doi     = {10.1038/s41592-019-0686-2}
}

@article{Hunter2007Matplotlib,
  author  = {Hunter, John D.},
  title   = {{Matplotlib}: A 2D Graphics Environment},
  journal = {Computing in Science \& Engineering},
  year    = {2007},
  volume  = {9},
  number  = {3},
  pages   = {90--95},
  doi     = {10.1109/MCSE.2007.55}
}

@article{Waskom2021Seaborn,
  author  = {Waskom, Michael L.},
  title   = {seaborn: statistical data visualization},
  journal = {Journal of Open Source Software},
  year    = {2021},
  volume  = {6},
  number  = {60},
  pages   = {3021},
  doi     = {10.21105/joss.03021}
}
```
