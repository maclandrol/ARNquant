# ARNquant

One method to study transcription is to perform a smFISH (single-molecule RNA fluorescence in situ hybridization).
smFISH differs from conventional approaches by using many short (about 20 base pairs long) oligonucleotide probes to target different regions of a RNA transcript of interest [1]. Individual RNA molecules appear as bright, diffraction-limited spots and can be detected using a gaussian filtering process [2].
2D Spot detection programs such as Localize, return a list of (x,y) coordinate and intensity for each spot. However there is still a need to reliably quantify RNA expression by distinguishing between nascent RNA, mature RNA and transcription site. Here I present ARNquant, a matlab tool for RNA quantification in 2D.

2D quantification is usually done by first identifying the single-RNA molecule intensity, which will later be used to find transcription site. Indeed, transcription sites are expected to be an aggregate of nascent and mature RNA. They will thus be detected as spots with high intensity.
Since transcription is a stochastic process [3], taking the mean of all spots intensity should give the single-RNA intensity ( _Si_ ) . By dividing each spot intensity by _Si_ , the number of RNA in this spot can be estimated. 
Although this sound simple and straightforward, spots intensity distribution are often positively skew, so taking the mean intensity to compute single-RNA intensity is actually a bad idea. It's then necessary to remove the outlier from the normal distribution before computing _Si_. ARNquant offer a list of method, including 1D clustering of spots' intensity, to automatically perform this.

## Getting started with ARNquant

![main windows](http://i.imgur.com/xAEqnh8.png)

#### 1. Nucleus segmentation mask
  Load your segmentation mask with the `Load mask` button and it will show in this panel.


#### 2. **Spots intensity distribution**
  Load the output of your spot detection and the intensity distribution will be shown on this panel. The spot detection output file should be in Localize format (a tsv file with each line representing a spot information where the first three columns are the following : `x-coordinate y-coordinate intensity`.  Each spot position will also automatically be shown on the mask as a red dot.


#### 3. **Intensity range selection** 
  Manually remove outliers from the intensity distribution by fixing a `min` and a `max` intensity. The red vertical line represent the max intensity and the green one the min intensity. Those two lines are updated when you move the sliders or enter a numerical value to provide immediate feedback. Discarded spots will also be removed from the mask image. Manual selection of intensity range is useful if you want to use the `mean algorithm` to compute single RNA intensity without introducing bias due to the skewness of the distribution.


#### 4. Global parameters
  







## See also

FishQuant [2] for the automatic counting of transcripts in 3D FISH images.



Ji, N., & van Oudenaarden, A. (2005). Single molecule fluorescent in situ hybridization (smFISH) of C. elegans worms and embryos.
Mueller, F., Senecal, A., Tantale, K., Marie-Nelly, H., Ly, N., Collin, O., ... & Zimmer, C. (2013). FISH-quant: automatic counting of transcripts in 3D FISH images. Nature Methods, 10(4), 277-278.
Li, G. W., & Xie, X. S. (2011). Central dogma at the single-molecule level in living cells. Nature, 475(7356), 308-315.
