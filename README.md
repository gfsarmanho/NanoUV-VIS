# NanoUV-VIS

> An interactive visualization tool for monitoring the evolution of optical properties of nanoparticles throughout synthesis reactions.

<img src="https://github.com/gfsarmanho/NanoUV-VIS/tree/master/www/hex_logo-R.png" width="150" height="150" align="right"/>

<p align="center">
<img src="www/Fig_plots.png" alt="Drawing" style="width:80%;"/>
</p>

# Installation

To use *NanoUV-VIS* in your PC, you will need first to install [R](https://cloud.r-project.org/) statistical software. We recommend to use [RStudio](https://www.rstudio.com/products/rstudio/download/) editor, in order to have a better experience with this [Shiny](https://shiny.rstudio.com/) application.

Second, you will need to install the following R packages: `shiny`, `shinyjs`, `knitr`, `rmarkdown`, `shinycssloaders`, `plotly`, `plot3D`, `DT`, `crosstalk` and `htmlwidgets`. To install them from CRAN, use the following R code:

```r
install.packages(c("shiny", "shinyjs", "knitr", "rmarkdown", "shinycssloaders", "plotly", "plot3D", "DT", "crosstalk", "htmlwidgets"))
```

After the local installation, you will be ready to run the application! Just type on R/RStudio console:

```r
shiny::runGitHub("usnistgov/NanoUV-VIS")
```
If there are no error messages, you are ready to enjoy *NanoUV-VIS*! 

> **NOTE:** Remember not to close R/RStudio after initiate the web application!

# Getting Started

After start the application, you should see the following screen on your web browser:

<p align="center">
<img src="www/screenshot.png" alt="Drawing" style="width:80%;"/>
</p>

Click in **About** tab to see an overview of the application, describing the way you should upload your UV-VIS data and the available interactive plots. For more detailed information and description of the tool using a real experiment, please check [Calderon-Jimenez et. al, 2017](http://nvlpubs.nist.gov/nistpubs/jres/122/jres.122.037.pdf).

In the following gif you will find some quick instructions to run *NanoUV-VIS* using an example data set.

![](www/demo.gif)

# Contributions
Any suggestions or contributions to this project (comments, questions, typo corrections, etc.) are more than welcome! You can do it directly on this GitHub repository or mail to <gabriel.sarmanho@nist.gov> and <bryan.calderonjimenez@nist.gov>.

# Citation

Please cite the following article when using *NanoUV-VIS* web application:

**Jimenez**, B. C., **Sarmanho**, G. F., Murphy, K. E., Bustos, A. R. M., & Baudrit, J. R. V. (2017). NanoUV-VIS: an interactive visualization tool for monitoring the evolution of optical properties of nanoparticles throughout synthesis reactions. **Journal of Research (NIST JRES)**, 122. DOI: https://doi.org/10.18434/M3T952


# License and Disclaimer Statement
This software is in the public domain. Click [here](DISCLAIMER.md) to access the disclaimer statement.
