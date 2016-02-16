# ARNquant

One method to study transcription is to perform a smFISH (single-molecule RNA fluorescence in situ hybridization).
smFISH differs from conventional approaches by using many short (about 20 base pairs long) oligonucleotide probes to target different regions of a RNA transcript of interest [1]. Individual RNA molecules appear as bright, diffraction-limited spots and can be detected using a gaussian filtering process [2].
2D Spot detection programs such as Localize, return a list of (x,y) coordinate and intensity for each spot. However there is still a need to reliably quantify RNA expression by distinguishing between nascent RNA, mature RNA and transcription site. 
Automatic counting of transcripts in 3D FISH images can be done with FishQuant (great tool btw!). Here I present ARNquant, a simpler tool for RNA quantification in 2D.

2D quantification is usually done by first identifying the single-RNA molecule intensity, which will later be used to find transcription site. Indeed, transcription sites are expected to be an aggregate of nascent and mature RNA. They will thus be detected as spots with high intensity.
Since transcription is a stochastic process [3], taking the mean of all spots intensity should give the single-RNA intensity ( _Si_ ) . By dividing each spot intensity by _Si_ , the number of RNA in this spot can be estimated. 
Although this sound simple and straightforward, spots intensity distribution are often positively skew, so taking the mean intensity to compute single-RNA intensity is actually a bad idea. 
It's necessary to remove the outlier from the normal distribution before computing _Si_. ARNquant offer a list of method to automatically perform this.

## Tutorial




Ji, N., & van Oudenaarden, A. (2005). Single molecule fluorescent in situ hybridization (smFISH) of C. elegans worms and embryos.
Mueller, F., Senecal, A., Tantale, K., Marie-Nelly, H., Ly, N., Collin, O., ... & Zimmer, C. (2013). FISH-quant: automatic counting of transcripts in 3D FISH images. Nature Methods, 10(4), 277-278.
Li, G. W., & Xie, X. S. (2011). Central dogma at the single-molecule level in living cells. Nature, 475(7356), 308-315.
