# ARNquant

One method to study transcription is to perform a smFISH (single-molecule RNA fluorescence in situ hybridization).
smFISH differs from conventional approaches by using many short (about 20 base pairs long) oligonucleotide probes to target different regions of a RNA transcript of interest [1]. Individual RNA molecules appear as bright, diffraction-limited spots and can be detected using a gaussian filtering process [2].
2D Spot detection programs such as Localize, return a list of (x,y) coordinate and intensity for each spot. However there is still a need to reliably quantify RNA expression by distinguishing between nascent RNA, mature RNA and transcription site. Here I present ARNquant, a matlab tool for RNA quantification in 2D.

2D quantification is usually done by first identifying the single-RNA molecule intensity, which will later be used to find transcription site. Indeed, transcription sites are expected to be an aggregate of nascent and mature RNA. They will thus be detected as spots with high intensity.
Since transcription is a stochastic process [3], taking the mean of all spots intensity should give the single-RNA intensity ( _Si_ ) . By dividing each spot intensity by _Si_ , the number of RNA in this spot can be estimated. 
Although this sound simple and straightforward, spots intensity distribution are often positively skew, so taking the mean intensity to compute single-RNA intensity is actually a bad idea. It's then necessary to remove the outlier from the normal distribution before computing _Si_. ARNquant offer a list of method, including 1D clustering of spots' intensity, to automatically perform this.

## Installation

Clone or download this repository, then [add the folder in your matlab search path](http://www.mathworks.com/help/matlab/matlab_env/add-folders-to-search-path-upon-startup-on-unix-or-macintosh.html).

You can launch the gui by running `ARNquant` in matlab.

## Getting started with ARNquant

![main windows](http://i.imgur.com/xAEqnh8.png)

#### 1. Nucleus segmentation mask
  Load your segmentation mask with the `Load mask` button and it will show in this panel.


#### 2. **Spots intensity distribution**
  Load the output of your spot detection with the `Load spot` button and the intensity distribution will be shown on this panel. The spot detection output file should be in Localize format (a tsv file with each line representing a spot information where the first three columns are the following : `x-coordinate y-coordinate intensity`.  Each spot position will also automatically be shown on the mask as a red dot. You can select each dot with the data cursor tool, ![data cursor tool](http://www.mathworks.com/help/matlab/creating_plots/datatip_icon.jpg) in order to get the position and intensity of the corresponding spot.  This is particularly useful if you find a really bright spot in your original image and would like to compare its intensity to the other spots.


#### 3. **Intensity range selection** 
  Manually remove outliers from the intensity distribution by fixing a `min` and a `max` intensity. The red vertical line represent the max intensity and the green one the min intensity. Those two lines are updated when you move the sliders or enter a numerical value to provide immediate feedback. Discarded spots will also be removed from the mask image. Manual selection of intensity range is useful if you want to use the `mean algorithm` to compute single RNA intensity without introducing bias due to the skewness of the distribution.


#### 4. Global parameters
  You can set global parameters for your quantification here

  1. **Trans coeff (c)** : let you choose the coeff to use in order to determine transcription site. A transcription site will then have an intensity ```I > c*MeanInt ```, Where  `MeanInt` is the single spot intensity.
  2. **Datatype** :  the type of RNA you're using. Some RNA are only localised in the nucleus while other in both cytoplasm and nucleus. You should set this parameters, because if your RNA is not exported to the cytoplasm, any spot found in the cytoplasm should be considered as noise. By default  both nuclear and cytoplasmic spots are expected.
  3. **BinSize** : use this to determine the binsize use for the intensity distribution. You should set a value that   will allow an acccurate representation of the intensity distribution of your spot. Failure to set the correct binsize could impact your analysis. [This section in the wikipedia page of Histogram](https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width) provide some method to choose the number of bins. I recommend playing with the values if you're not sure.
  
  
#### 5. 6. 7. Single spot intensity computation
 
This section let you choose the algorithm to compute the single spot intensity.
 
##### Mean and Median 
 1. **Mean** : This just use the mean intensity of all spots while excluding the spots that are outside the interval set by the user.
 2. **Median** : This is similar to the previous method, but it use the median instead. In a truly normal distribution the mean and the median are the same. However, the median is more resistant to outliers than the mean, which is a one of the reason to prefer it over the mean. By using deviation around the median, it's actually  possible to remove outliers. The ``"Std from median"`` parameters in  **SECTION 7** let you choose the coefficient `k` that should be used to exclude outliers according to the following equation : 

				isOutlier = abs(spotINT - median(INT)) > k*std(INT);
         
    where `spotINT` is the spot intensity and `INT` is a vector with all spot intensity.
    
    
If you click on `Apply`, a new window with the nascents RNA distribution in your spot list will be shown. The thick black line represent the single spot intensity and is expected to fit well as the mean of the intensity distribution of spots with only one nascent RNA. If you're are not satisfied, you can close this window and change the parameters. Results are only saved when you click on the `OK` button.
  
    
   ![](https://i.imgur.com/bj3ERBm.png)
             

##### Classification

  Select this option, if you want to use alternative methods (unsupervised data clustering in 1D) to group your spots into several clusters then use the best cluster to compute single spot intensity. For this particular problem, grouping your data into 2 or 3 classes should be adequate. If you're are not sure, use the automatic clustering, selected by default. The following methods are available  :
  * **kmeans** is a clustering algorithm that group a given data set into k fixed clusters (a priori). It's a two-phase iterative  algorithm that minimize the sum of point-to-centroid distances, summed over all k clusters. It often falls in local minima and is not really suited for 1D clustering.
  * **ckmeans** is a dynamic programming algorithm for optimal one-dimensional kmeans clustering that minimizes the sum of squares of within cluster distances. As an alternative to the standard heuristic kmeans algorithm, this algorithm guarantees optimality and repeatability [4].
  * **Gaussian mixture model**. This approach fit a Gaussian mixture models (GMM) to the data and then assign each spot to the multivariate normal components that maximise the posterior probability, given the data.
  * **KDE ( not implemented ).** This approach use kernel density estimation and split the data at the local minima of the density. It's however difficult to specify a number of cluster here.

**I recommend using ckmeans or in extreme case, the gaussian mixture model.**
 
If you click on `Apply`, the window shown here will be a little different. The upper subplot will show the different clusters and the intensity distribution inside those clusters, while the lower subplot is similar to the previous figure (intensity distribution when mean or median is used). ARNquant try to guess the best cluster to use in order to compute the single spot intensity. However there is an option at the bottom of the figure that enable correction if the wrong cluster has been chosen. Click `ok` when you're done.

![classification window](https://i.imgur.com/g73zoxf.jpg)

## See also

FishQuant [2] for the automatic counting of transcripts in 3D FISH images.


## References
Ji, N., & van Oudenaarden, A. (2005). Single molecule fluorescent in situ hybridization (smFISH) of C. elegans worms and embryos.

Mueller, F., Senecal, A., Tantale, K., Marie-Nelly, H., Ly, N., Collin, O., ... & Zimmer, C. (2013). FISH-quant: automatic counting of transcripts in 3D FISH images. Nature Methods, 10(4), 277-278.

Li, G. W., & Xie, X. S. (2011). Central dogma at the single-molecule level in living cells. Nature, 475(7356), 308-315.

https://cran.r-project.org/web/packages/Ckmeans.1d.dp/Ckmeans.1d.dp.pdf
