
<img src="https://github.com/michbur/AmpGram/blob/master/inst/AmpGram/AMP_log.png" alt="AmpGram" style="height: 200px;"/>

Identify antimicrobial peptides
-------------------------

AmpGram identifies antimicrobial peptides using n-gram encoding and random forests. It can be also accessed as a web-based service http://www.smorfland.uni.wroc.pl/shiny/AmpGram/. 

Local instance of AmpGram
------------------------
AmpGram can be used installed from CRAN as the R package:

```R
install.packages("AmpGram")
```

You can install the latest development version of the package:

```R
source("https://install-github.me/michbur/AmpGram")
```

After installation GUI can be accessed locally:

```R
library(AmpGram)
AmpGram_gui()
```


Installing dependency: AmpGramModel
------------------------
To be able to use AmpGram properly, you should have installed 'AmpGramModel' package available via GitHub. 
AmpGramModel contains stacked random forest model and informative n-grams required for prediction of antimicrobial peptides.
Due to the large size of a model, it needs to be stored in the external repository, as CRAN do not allow upload of files
larger than 5 MB. 

You can install AmpGramModel using the AmpGram function:

```R
install_AmpGramModel()
```

Antimicrobial peptides might be also identified in the batch mode:

```R
library(AmpGram)
library(biogram)
install_AmpGramModel()
sequences <- read_fasta(system.file("AmpGram/prots.txt", package = "AmpGram"))
predict(sequences,"new_file")
```
Unix/macOS: curl
------------------------
Function that installs AmpGramModel uses devtools library to install package from GitHub. 
On Unix and macOS systems you may encounter error concerning curl library. To install devtools 
dependencies run the following command:

```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev
```
