## General information
This repo contains all the scripts and data for the following paper: 

**Toro-Serey, C., Kane, G., & McGuire, J. (2021, May 31). Apparent preferences for cognitive effort fade when multiple forms of effort and delay are interleaved in a foraging environment. [Preprint](psyarxiv.com/5ygwh)**

## Reproduce the preprint

The paper can be reproduced in pdf format by running (well, Knitting) `Manuscript.Rmd` in R (blogdown required). The markdown file is organized so that code for the analyses precedes each written section. The files `sessionInfo_iMac.txt` and `sessionInfo_m1.txt` contain information about the environments used and tested to produce this manuscript. 

### Docker environment 

#### *This method does not work on M1 Macs yet because of Rstudio incompatibility in rocker/verse*

You can also use a [Docker](https://www.docker.com/) container to run `Manuscript.Rmd` and produce the preprint pdf [(Image produced based on this tutorial.)](https://ropenscilabs.github.io/r-docker-tutorial/)

Once you install Docker:

- Clone this repo

- Run the following command, where `<yourpath>` is the path to cloned repository (this command will download a ~3GB image to your machine):

```
docker run --rm -p 8787:8787 -e PASSWORD=foraging -v <yourpath>:/home/rstudio/Cost_studies ctoroserey/cost_studies:preprintenv
```

- In your browser of choice, go to `http://localhost:8787` (user: rstudio, password: foraging), then make sure you're in the Cost_studies working directory (`setwd(‘Cost_studies’)`) and that you see all the repo files. Open `Manuscript.Rmd`. This will allow you to load the data from the repo directory, and everything produced within the Docker container will be stored locally within the repo directory. 

- To produce the pdf, click on the `Knit` button towards the top of the window. This will take ~5 mins. If you want the figures at the end of the manuscript, set figsEnd to `TRUE` (within the Setups chunk).

Of course, loading this environment will also let you run sections of the code that you might be interested in. As long as you run the `Setups` and both `Load data *` sections, you should be able to run any section by itself. *Note that the models are not run within the manuscript file (results are loaded from previous fits).* To examine the model in detail, run `modeling_btw/wth.R` withi the Models directory.





