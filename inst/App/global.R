#===================#
# App:  NanoUV-VIS  #
#                   #
# File: global.R    #
#===================#

#================#==============================================================
# ggplot2 colors #
plot_colors <- c("#619CFF", "#F8766D", "#00BA38")
# gg_blue  <- "#619CFF"
# gg_red   <- "#F8766D"
# gg_green <- "#00BA38"
# gg_gray  <- "#CCCCCC80"

#===============================================================================
# Transform ".rmd" about data to ".md"
about_file_rmd <- "./www/ABOUT.Rmd"
sapply(X = about_file_rmd, FUN = knitr::knit, quiet = TRUE)

#==============================================================================#
#                             HTML/CSS functions 
#==============================================================================#

#-------------------------------------------------------------------------------
# Create a little question mark link that shows a help popup on hover
helpPopup <- function(content, title = NULL) {
  a(href = "#",
    class = "popover-link",
    `data-toggle`  = "popover",
    `data-title`   = title,
    `data-content` = content,
    `data-html`    = "true",
    `data-trigger` = "hover",
    icon("question-circle")
  )
}

#-------------------------------------------------------------------------------
# Create an upper red asterisk, indicating a mandatory input
redAsterisk <- function(label) {
  tagList(label, span("*", class="mandatory_star"))
}

#-------------------------------------------------------------------------------
# spinner for "load" plot status
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    span(
      class="btn-loading-container",
      hidden(
        img(src="ajax-loader-bar.gif", class="btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

withBusyIndicatorServer <- function(buttonId, expr) {
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl    <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl     <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector=loadingEl)
  })
  
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

#-------------------------------------------------------------------------------
# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl      <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg   <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

#==============================================================================#
#                             R functions 
#==============================================================================#

#-------------------------------------------------------------------------------
# Function: "join_zz_table()" - merge single files and extract unities
#-------------------------------------------------------------------------------
# Data to test
# files_names <- arq_names
#-------------------------------------------------------------------------------
join_zz_table <- function(files_names){
  
  # Create NULL vectors to store objects from all data
  zz <- NULL
  label_x <- label_y <- NULL
  unities <- NULL
  
  # FOR: read each ".csv" file, merge spectrum data and evaluate FWHM statistics
  for(i in 1:length(files_names)){
    #i=1
    dados_0 <- read.csv(files_names[i], header=FALSE, skip=0, colClasses="character")
    #head(dados_0)
    
    #----------------#
    # Pre-processing # Clean data and save information
    #----------------#
    # Save "times"
    line1      <- unlist(strsplit(dados_0[1, 2], " "))
    label_x[i] <- as.numeric(line1[1])
    
    # Save "unities" and wavelength (label_y) (from first file only)
    if(i==1){
      unit_x  <- line1[2]
      unit_y  <- dados_0[2, 1]
      unit_z  <- dados_0[2, 2]
      label_y <- as.numeric(dados_0[-(1:2), 1])
    }
    unities <-c(unit_x, unit_y, unit_z)
    
    # Save matrix with spectrum and wavelength vectors only
    dados <- matrix(c(as.numeric(dados_0[-(1:2), 1]),
                      as.numeric(dados_0[-(1:2), 2])), ncol=2)
    colnames(dados) <- c("wvlen", "spec")
    #head(dados)
    
    # Append data
    zz <- cbind(zz, dados[, 2])
    #head(zz)
  }
  
  # Order vectors and matrices, according to the "time" (x axis)
  ind_ord_x <- order(label_x)
  label_x   <- label_x[ind_ord_x]
  zz        <- zz[, ind_ord_x]
  
  zz_table <- data.frame(label_y, zz)
  names(zz_table) <- c("wavelength", label_x)
  
  # Create a list with objects to be returned
  RES <- list(
    zz_table = zz_table,
    label_x = label_x,
    label_y = label_y,
    unities = unities
  )
  return(RES)
}

#-------------------------------------------------------------------------------
# Function: "peak_uvvis()" - evaluate FWHM statistics for each wavelength data
#-------------------------------------------------------------------------------
# Data to test
# spec=zz[, 35] ; wvlen=label_y
#-------------------------------------------------------------------------------
peak_uvvis <- function(spec, wvlen=NULL){
  # Handle NULL wavelengths
  if(is.null(wvlen)) wvlen <- 1:length(spec)
  ind_wl <- order(wvlen)
  wvlen <- wvlen[ind_wl]
  spec <- spec[ind_wl]
  names(spec) <- as.character(wvlen)
  
  # Find maximun value of spectrum (peak or mode) 
  ind_max_y <- which(spec==max(spec, na.rm=TRUE))
  max_x     <- mean(wvlen[ind_max_y])
  max_y     <- max(spec, na.rm=TRUE)
  
  #------------------------------------------------#
  # Split spectrum according to it's maximun value #
  #------------------------------------------------#
  # BEFORE peak (group 1 or "g1") #---------------------------------------------
  ind_g1 <- which(wvlen <= max_x)
  wl_g1   <- wvlen[ind_g1]
  spec_g1 <- spec[ind_g1]
  
  # Find minimun spectrum value (x-y coordinates) of "g1"
  ind_min_y_g1 <- which(spec_g1==min(spec_g1, na.rm=TRUE))
  min_x_g1     <- max(wl_g1[ind_min_y_g1])
  min_y_g1     <- min(spec_g1, na.rm=TRUE)
  
  # AFTER peak (group 2 or "g2") #----------------------------------------------
  ind_g2  <- which(wvlen > max_x)
  wl_g2   <- wvlen[ind_g2]
  spec_g2 <- spec[ind_g2]
  
  # Find minimun spectrum value (x-y coordinates) of "g2"
  ind_min_y_g2 <- which(spec_g2==min(spec_g2, na.rm=TRUE))
  min_x_g2     <- min(wl_g2[ind_min_y_g2])
  min_y_g2     <- min(spec_g2, na.rm=TRUE)
  
  #----------------------------------------------------------------------------#
  # Find x-y coordinates of "minimun absolute" point, given by the intersec.   #
  # of vertical line (wavelenght) of maximun absorbance and the line which     #
  # conects the minimuns points of of g1 and g2                                #
  #----------------------------------------------------------------------------#
  m     <- (min_y_g2 - min_y_g1)/(min_x_g2 - min_x_g1)
  min_x <- max_x
  min_y <- min_y_g1 + m*(min_x - min_x_g1)
  
  # Evaluate "half height" of spectrum (median point)
  half_h <- (max_y + min_y)/2
  
  #----------------------------------------------------------------------------#
  # Evaluate intersection point between "half height" and the spectrum values
  #----------------------------------------------------------------------------#
  # BEFORE peak (group 1 or "g1") #---------------------------------------------
  if(any(spec_g1==half_h)){
    ind_diff_spec_g1 <- which(spec_g1==half_h)
    intersec_g1      <- wl_g1[ind_diff_spec_g1][1]
  } else {
    diff_g1 <- spec_g1 - half_h
    # Linear interpolation
    ind_M <- as.numeric(names(which.min(diff_g1[which(diff_g1>0)])))
    ind_m <- as.numeric(names(which.max(diff_g1[which(diff_g1<0)])))
    M <- spec[which(wvlen==ind_M)]
    m <- spec[which(wvlen==ind_m)]
    intersec_g1 <- approx(y=c(ind_m, ind_M), x=c(m, M),
                          xout=half_h, method="linear", n=50)$y
  }
  
  # AFTER peak (group 2 or "g2") #----------------------------------------------
  if(any(spec_g2==half_h)){
    ind_diff_spec_g2 <- which(spec_g2==half_h)
    intersec_g2      <- wl_g2[ind_diff_spec_g2][1]
  } else {
    diff_g2 <- spec_g2 - half_h
    # Linear interpolation
    ind_M <- as.numeric(names(which.min(diff_g2[which(diff_g2>0)])))
    ind_m <- as.numeric(names(which.max(diff_g2[which(diff_g2<0)])))
    M <- spec[which(wvlen==ind_M)]
    m <- spec[which(wvlen==ind_m)]
    intersec_g2 <- approx(y=c(ind_m, ind_M), x=c(m, M),
                          xout=half_h, method="linear", n=50)$y
  }
  
  # Evaluate FWHM
  intersec_x <- c(intersec_g1, intersec_g2)
  FWHM       <- diff(intersec_x)
  
  # Output
  out <- list()
  out <- c(max_y, FWHM, max_x,                     # important statistics
           min_y, min_x, half_h,                   # coord. of min + halh-height
           intersec_g1, intersec_g2,               # intersec(s) of half_h
           min_y_g1, min_x_g1, min_y_g2, min_x_g2) # coord. of min (g1 and g2)
  names(out) <- c("Abs_max", "FWHM", "Wl_max",
                  "Abs_min", "Wl_min", "Half_height",
                  "FWHM_L1", "FWHM_L2",
                  "Abs_min_1", "Wl_min_1", "Abs_min_2", "Wl_min_2")
  return(out)
} 

#-------------------------------------------------------------------------------
# Function: "plot_spec()" - create plot_ly graphs of spectrums
#-------------------------------------------------------------------------------
# Data to test
# my_i=35
# ZZ=zz[, my_i]
# UV=uvvis_data[my_i, ]
# labY=label_y
# labX=label_x[my_i]
#-------------------------------------------------------------------------------
plot_spec <- function(ZZ, UV, labX, labY, unX, unY){
  # inputs
  spec       <- ZZ
  min_y      <- as.numeric(UV["Abs_min"])
  max_y      <- as.numeric(UV["Abs_max"])
  max_x      <- as.numeric(UV["Wl_max"])
  intersec_x <- as.numeric(UV[c("FWHM_L1", "FWHM_L2")])
  half_h     <- as.numeric(UV["Half_height"])
  FWHM       <- as.numeric(UV["FWHM"])
  min_y_g1   <- as.numeric(UV["Abs_min_1"])
  min_x_g1   <- as.numeric(UV["Wl_min_1"])
  min_y_g2   <- as.numeric(UV["Abs_min_2"])
  min_x_g2   <- as.numeric(UV["Wl_min_2"])
  
  # Define plot colors
  plot_colors <- c("#619CFF", "#F8766D", "#00BA38")
  
  # Order data
  dados <- data.frame(spec=as.numeric(spec),
                      wl=as.numeric(labY))[order(labY), ]
  
  # Text in legentes (hover)
  hover_text0 <- paste("<b>Wavelength:</b> ", dados$wl, " ", unY, "<br>",
                       "<b>Absorbance:</b> ", dados$spec, sep="")
  hover_text1 <- paste("<b>Wavelength:</b> ", round(intersec_x, 2), " ", unY, "<br>",
                       "<b>Absorbance:</b> ", round(c(half_h, half_h), 4), sep="")
  hover_text2 <- paste("<b>Wavelength:</b> ", max_x, " ", unY, "<br>",
                       "<b>Absorbance:</b> ", max_y, sep="")
  hover_text3 <- paste("<b>FWHM:</b> ", round(FWHM, 2), " ", unY, sep="")
  
  # Margins and fonts
  marg        <- list(b = 87, l = 68, t = 65, r = 28, pad = 0)
  font_title  <- list(size = 20, family = "Arial", color = "black")
  font_labels <- list(size = 14, family = "Arial", color = "black")
  font_ticks  <- list(size = 14, family = "Arial", color = "black")
  
  # Plot
  plot_ly(dados, x=~wl) %>%
    # Wavelength data
    add_lines(y=~spec, line=list(color=plot_colors[1], width=4),
              hoverinfo="text", text=hover_text0, name="Spetrum") %>%
    add_markers(x=round(intersec_x, 2), y=c(half_h, half_h),
                marker=list(size=10, color=plot_colors[3]),
                hoverinfo="text", text=hover_text1, showlegend=FALSE) %>%
    add_markers(x=max_x, y=max_y,  marker=list(size=10, color=plot_colors[2]),
                hoverinfo="text", text=hover_text2, showlegend=FALSE) %>%
    add_trace(x=seq(from=intersec_x[1], to=intersec_x[2], length=100),
              y=rep(half_h, 100),
              mode="lines", line=list(color=plot_colors[3], width=4),
              hoverinfo="text", text=hover_text3, showlegend=FALSE) %>%
    layout(title       = paste("<b>Time", " = ", labX, " (", unX, ")</b>", sep=""),
           titlefont   = font_title,
           showlegend  = FALSE,
           xaxis       = list(title=paste("Wavelength (", unY, ")", sep=""),
                              titlefont=font_labels, tickfont=font_ticks,
                              linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                              tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE),
           yaxis       = list(title="Abs", titlefont=font_labels, tickfont=font_ticks,
                              linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                              tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE),
           margin      = marg,
           annotations = list(
             list(xref="x", yref="y",
                  x=intersec_x[2]*1.4, y=half_h, showarrow=FALSE, align="right",
                  text=hover_text3, font=list(size=15, color=plot_colors[3]),
                  showarrow=TRUE),
             list(xref="x", yref="y", x=max_x, y=max_y, text="SPR peak",
                  align="right", showarrow = TRUE, arrowhead=100,
                  col=plot_colors[2])
           )
    )
}

#==============================================================================#
#================================= END ========================================#
#==============================================================================#