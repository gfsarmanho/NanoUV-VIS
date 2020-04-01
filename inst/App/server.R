#===================#
# App:  NanoUV-VIS  #
#                   #
# File: server.R    #
#===================#

#====================#
# Begin Shiny Server #
#====================#
shinyServer(function(input, output) {
  
  #-----------------------------------------------------------------------------
  # Manage Example data set to download (Unique csv)
  output$downloadExampleData <- downloadHandler(
    filename = function() { "example-to-download.csv" },
    content = function(file) {
      data_example <- read.csv("./www/example-to-download.csv")
      colnames(data_example)[-1] <- gsub(pattern="X", replacement="",
                                         x=colnames(data_example)[-1])
      write.csv(x=data_example, file, row.names=FALSE)
    }
  )
  
  #-----------------------------------------------------------------------------
  # Manage Example data set to download
  output$downloadExampleMult <- downloadHandler(
    filename = function() { "example-multiple-data.zip" },
    content = function(file) {
      data_example <- read.csv("./www/example-to-download.csv")
      colnames(data_example)[-1] <- gsub(pattern="X", replacement="",
                                         x=colnames(data_example)[-1])
      write.csv(x=data_example, file, row.names=FALSE)
    }
  )
  
  #-----------------------------------------------------------------------------
  # Create reactive variables - temporary clones of input$file1 and input$filesN
  rv_file1  <- reactive(x={ if(is.null(input$file1)) return(NULL) else return(1) })
  rv_filesN <- reactive(x={ if(is.null(input$filesN)) return(NULL) else return(1) })

  # Activate "create table" button whenever there is any data
  observe(x={
    toggleState(id="loadBtn",   condition={ !is.null(rv_file1())  })
    toggleState(id="createBtn", condition={ !is.null(rv_filesN()) })
  })
  
  #-----------------------------------------------------------------------------
  # Reset buttons and variables whenever press "reset" button

  # Case 1: Single data
  observeEvent(eventExpr={input$resetBtnSINGLE}, handlerExpr={
    RV$zz_table <- NULL
    RV$zz_rev   <- NULL
    shinyjs::hide(id="down_button_data2")
    shinyjs::hide(id="down_button_FWHM")
    shinyjs::hide(id="visualTab2")
    shinyjs::show(id="message_dat")
    shinyjs::show(id="message_vis")
    shinyjs::show(id="message_res")
    
    shinyjs::reset(id="file1")
    shinyjs::disable(id="loadBtn")
    shinyjs::disable(id="createBtn")
    shinyjs::disable(id="resetBtnSINGLE")
  })
  
  #Case 2: Multiple data
  observeEvent(eventExpr={input$resetBtnMULT}, handlerExpr={
    RV$zz_table <- NULL
    RV$zz_rev   <- NULL
    shinyjs::hide(id="down_button_data2")
    shinyjs::hide(id="down_button_FWHM")
    shinyjs::hide(id="visualTab2")
    shinyjs::show(id="message_dat")
    shinyjs::show(id="message_vis")
    shinyjs::show(id="message_res")
    
    shinyjs::reset(id="filesN")
    shinyjs::disable(id="createBtn")
    shinyjs::disable(id="loadBtn")
    shinyjs::disable(id="resetBtnMULT")
  })
  
  #-----------------------------------------------------------------------------
  # Create reactive objects to store values to be used in graphs
  FWHM <- reactiveValues(
    FWHM_stats=NULL,
    uvvis_data=NULL
  )
  RV <- reactiveValues(
    zz_table    = NULL,
    unities     = NULL,
    label_x     = NULL,
    label_y     = NULL,
    nX          = NULL,
    nY          = NULL,
    zz          = NULL,
    zz_rev      = NULL,
    unit_name   = NULL,
    zz_tab_name = NULL,
    arq_names   = NULL,
    dados       = NULL
  )
  
  # Reactive values for plots
  Rplot <- reactiveValues(
    p_surface = NULL,
    p_contour = NULL,
    p_FWHM     = NULL,
    p_spec_mix = NULL
  )
  #-----------------------------------------------------------------------------
  # Main routine (after press any of "create/load table" buttons
  observeEvent(eventExpr={
    input$createBtn | input$loadBtn
  },
  handlerExpr={
    if(input$createBtn==0 && input$loadBtn==0){
      return()
    } else {
      
      shinyjs::hide(id="message_dat", time=0.2)
      shinyjs::hide(id="message_vis", time=0.1)
      shinyjs::hide(id="message_res", time=0.1)
      shinyjs::show(id="visualTab2")
      
      if( !is.null(rv_file1()) ) {
        withBusyIndicatorServer(buttonId="loadBtn", expr={
          #-----------------------# 
          # Option 1: Unique file #
          #-----------------------#
          RV$arq_names      <- input$file1$datapath
          RV$dados          <- read.csv(RV$arq_names, sep=",", dec=".",
                                        header=TRUE) 
          RV$unities        <- c(input$unit_x, input$unit_y, input$unit_z) 
          names(RV$unities) <- c("Time", "Wavelength", "Absorbance")
          
          # Create label vectors
          RV$label_x <- as.numeric(gsub(pattern="X", replacement="",
                                        x=names(RV$dados)[-1]))
          RV$label_y <- RV$dados[, 1]
          
          # Create data matrix
          RV$zz           <- as.matrix(RV$dados[, -1])
          colnames(RV$zz) <- RV$label_x
          rownames(RV$zz) <- as.numeric(RV$label_y)
          
          # Create reversal zz matrix (for plots)
          RV$zz_rev <- RV$zz[order(RV$label_y), ]
          
          # Customize "zz_table"
          RV$zz_table <- RV$dados
          names(RV$zz_table)[-1] <- RV$label_x
          
          # Constants
          RV$nX <- length(RV$label_x)
          RV$nY <- length(RV$label_y)
          
          # show "reset" button
          shinyjs::show(id="resetBtnSINGLE")

        }) #end of withBusyIndicatorServer()
      }
      
      if( !is.null(rv_filesN()) ) {
        withBusyIndicatorServer(buttonId="createBtn", expr={
          
          #--------------------------# 
          # Option 2: Multiple Files #
          #--------------------------# 
          RV$arq_names <- input$filesN$datapath
          
          # create the zz table using "join_zz_table()" function
          RV$dados          <- join_zz_table(RV$arq_names)
          RV$zz_table       <- RV$dados$zz_table
          RV$unities        <- RV$dados$unities
          names(RV$unities) <- c("Time", "Wavelength", "Absorbance")
          
          # Create label vectors
          RV$label_x <- RV$dados$label_x
          RV$label_y <- RV$dados$label_y
          
          # Create data matrix
          RV$zz           <- as.matrix(RV$zz_table[, -1])
          colnames(RV$zz) <- RV$label_x
          rownames(RV$zz) <- RV$label_y
          
          # Create reversal zz matrix (for plots)
          RV$zz_rev <- RV$zz[order(RV$label_y), ]
          
          # Customize "zz_table"
          RV$unit_name   <- paste(c("time-", "wave-", "abs-"), RV$unities,
                                  sep="", collapse="_")
          RV$zz_tab_name <- paste("zz_table_(", RV$unit_name, ").csv", sep="")
          
          # Constants
          RV$nX <- length(RV$label_x)
          RV$nY <- length(RV$label_y)
          
          # show "reset" and "download table" buttons
          shinyjs::show(id="resetBtnMULT")
          shinyjs::show(id="down_button_data2")
          
        }) #end of withBusyIndicatorServer()
      }
      
      ##########################################################################
      ########################## PART 2 (FWHM stats) ###########################
      ##########################################################################
      # Apply "peak_uvvis()" function and generate table with statistics
      FWHM$FWHM_stats <- apply(RV$zz, 2, FUN=peak_uvvis, wvlen=RV$label_y)
      FWHM$uvvis_data <- data.frame(Time=RV$label_x, t(FWHM$FWHM_stats))
      
      # Create a list with spectrum plots
      for(i in 1:RV$nX) {
        local({
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          output[[plotname]] <- renderPlotly({
            p <- plot_spec(ZZ=RV$zz[, my_i], UV=FWHM$uvvis_data[my_i, ],
                           labX=RV$label_x[my_i], labY=RV$label_y,
                           unX=RV$unities[1], unY=RV$unities[2])
            
            # Ignoring explicitly provided widget ID, since Shiny doesn't use them
            p$elementId <- NULL
            p
          })
        })
      }
      
      # Show tabs and download buttons
      shinyjs::show(id="tab_zzTable")
      
      shinyjs::show(id="down_button_FWHM")
      shinyjs::show(id="tab_resultsFWHM")

    }# end of if()
  }) #end of observeEvent()
  
  ##########################################################################
  ################################# TABS ###################################
  ##########################################################################
  
  #-----------------------------------------------------------------------------
  # TAB PANEL: Table with data #
  #----------------------------#
  output$tab_zzTable <- DT::renderDataTable(expr={
    if(is.null(RV$zz_rev)){
      return()
    } else {
      tabela <- RV$zz_table
      names(tabela)[1] <- "Wavelength/Time"
      DT::datatable(tabela, options=list(searching=TRUE, searchHighlight=TRUE,
                                         pageLength=100, paging=TRUE))
    }
  })
  
  # Button to townload the created joint table 
  output$down_button_data <- downloadHandler(
    filename = function() { "tab-dataset.csv" },
    content  = function(file) {
      write.csv(x=RV$zz_table, file, row.names=FALSE)
    }
  )
  
  #-----------------------------------------------------------------------------
  # TAB PANEL: Table with results of FWMH calculations #
  #----------------------------------------------------#
  output$tab_resultsFWHM <- DT::renderDataTable(expr={
    if(is.null(RV$zz_rev)){
      return()
    } else {
      tabelaUVVIS <- as.data.frame(FWHM$uvvis_data)
      names(tabelaUVVIS) <- c("Time(min)",
                              "Abs.max", "FWHM(nm)", "SPR.peak(nm)",
                              "P3.absorbance", "P3.wavelength(nm)", "HH",
                              "FWHM.L1(nm)", "FWHM.L2(nm)",
                              "P1.absorbance", "P1.wavelength(nm)",
                              "P2.absorbance", "P2.wavelength(nm)")
      DT::datatable(data=tabelaUVVIS, rownames=FALSE,
                    options=list(searching=TRUE, searchHighlight=TRUE,
                                 autoWidth=TRUE, paging=FALSE) ) %>%
        DT::formatRound(columns=c("FWHM(nm)",
                                  "P3.absorbance", "HH",
                                  "FWHM.L1(nm)", "FWHM.L2(nm)",
                                  "P1.absorbance", "P2.absorbance"),
                        digits=4)
                        
    }
  })
  
  # Button to townload Table with data
  output$downloadBtnFWHM <- downloadHandler(
    filename = function() { "tab_results.csv" },
    content = function(file) {
      tabelaUVVIS <- as.data.frame(FWHM$uvvis_data)
      names(tabelaUVVIS) <- c("Time(min)",
                              "Abs.max", "FWHM(nm)", "SPR.peak",
                              "Abs.min", "SPR.min", "Half.height",
                              "FWHM.L1(nm)", "FWHM.L2(nm)",
                              "Abs.min.01", "SPR.min.01",
                              "Abs.min.02", "SPR.min.02")
      write.csv(x=tabelaUVVIS, file, row.names=FALSE)
    }
  )
  
  ##########################################################################
  ################################ PLOTS ###################################
  ##########################################################################
  
  #-----------------------------------------------------------------------------
  # TAB: SURFACE (3-D surface plot) #
  #---------------------------------#
  output$plot_surface <- renderPlotly({
    if(is.null(RV$zz_rev)){
      return()
    } else {
      # Plot color scheme
      plot_colors0  <- colorRampPalette(RColorBrewer::brewer.pal(5, "Spectral"))(30)
      
      # Create matrix with showed names when click (hover)
      f_pasteXY <- function(x, y){
        paste("<b>Time:<b> ", x, " ", RV$unities[1],
              "<br>Wavelength:</b> ", y, " ", RV$unities[2], sep="")
      }
      hover_surface1 <- outer(RV$label_x, rev(RV$label_y), FUN=f_pasteXY)
      hover_surface2 <- paste("<br>Absorbance: ", t(RV$zz_rev), sep="")
      hover_surface  <- matrix(paste(hover_surface1, hover_surface2, sep=""), ncol=length(RV$label_y))

      # Scene parameters for all axis
      marg <- list(t=65, pad=0) #b=87, l=68, r=28
      
      font_title  <- list(size=20, family="Arial", color="black")
      font_labels <- list(size=11, family="Arial", color="black")
      font_ticks  <- list(size=11, family="Arial", color="black")
      
      x_breaks <- hist(RV$label_x, plot=FALSE)$breaks
      x_axis   <- list(title=paste("Time (", RV$unities[1], ")", sep=""),
                       tickmode="array", 
                       tickvals=c(0, which(RV$label_x %in% x_breaks)-1), 
                       ticktext=c(RV$label_x[1], RV$label_x[which(RV$label_x %in% x_breaks)]),
                       tickwidth=0.2, ticklen=10, tickangle=0, zeroline=FALSE,
                       titlefont=font_labels, tickfont=font_ticks, showgrid=TRUE)
      
      y_breaks  <- hist(rev(RV$label_y), plot=FALSE)$breaks
      y_axis    <- list(title=paste("Wavelength (", RV$unities[2], ")", sep=""),
                        tickmode="array", ticktext=y_breaks,
                        tickvals=which(rev(RV$label_y) %in% y_breaks)-1,
                        tickwidth=0.2, ticklen=2, tickangle=0,
                        titlefont=font_labels, tickfont=font_ticks, showgrid=TRUE) 
      
      z_axis <- list(title="Abs", titlefont=font_labels, tickfont=font_ticks, showgrid=TRUE)
      scene_surface_invertXY <- list(xaxis=y_axis, yaxis=x_axis, zaxis=z_axis,
                                     camera=list(eye=list(x=-0.7, y=-2.5, z=0.1)))
      
      # Plot graph and save
      Absorbance <- t(RV$zz_rev)
      p_surface  <- plot_ly(z=~Absorbance, type="surface", colors=plot_colors0,
                           hoverinfo="text", text=hover_surface,
                           #height="1000px", width="1000px",
                           showscale=TRUE, source="source_surface") %>%
        layout(title="<b>3D plot</b>", titlefont=font_title, margin=marg,
               scene=scene_surface_invertXY, autosize=TRUE, showlegend=FALSE)
      
      # Ignoring explicitly provided widget ID, since Shiny doesn't use them
      p_surface$elementId <- NULL
      p_surface
      
      # Save plot in a reactive object
      Rplot$p_surface <- p_surface
    }
  })

  # Button to townload plot in html
  output$downloadBtnSurf <- downloadHandler(
    filename = function() { "plot-surface.html" },
    content  = function(file) {
      p_surface2 <- Rplot$p_surface
      htmlwidgets::saveWidget(widget=p_surface2, file)
    }
  )

  #-----------------------------------------------------------------------------
  # TAB: "CONTOUR" (contour plot) #
  #-------------------------------#
  output$plot_contour <- renderPlotly({
    if(is.null(RV$zz_rev)){
      return()
    } else {
      # Plot color scheme
      plot_colors0  <- colorRampPalette(RColorBrewer::brewer.pal(5, "Spectral"))(30)
      
      # Create matrix with showed names when click (hover)
      f_pasteXY <- function(x, y){
        paste("Time: ", x, " ", RV$unities[1],
              "<br>Wavelength: ", y, " ", RV$unities[2], sep="")
      }
      hover_surface1 <- outer(RV$label_x, rev(RV$label_y), FUN=f_pasteXY)
      hover_surface2 <- paste("<br>Absorbance: ", t(RV$zz_rev), sep="")
      hover_surface  <- matrix(paste(hover_surface1, hover_surface2, sep=""),
                               ncol=length(RV$label_y))
      
      # axis parameters
      marg        <- list(b=87, l=68, t=65, r=28, pad=0)
      
      font_title  <- list(size=20, family="Arial", color="black")
      font_labels <- list(size=16, family="Arial", color="black")
      font_ticks  <- list(size=16, family="Arial", color="black")
      
      x_breaks <- hist(RV$label_x, plot=FALSE)$breaks
      x_axis   <- list(title=paste("Time (", RV$unities[1], ")", sep=""),
                       tickmode="array", 
                       tickvals=c(0, which(RV$label_x %in% x_breaks)-1), 
                       ticktext=c(RV$label_x[1], RV$label_x[which(RV$label_x %in% x_breaks)]),
                       tickwidth=0.2, ticklen=10, tickangle=0, zeroline=FALSE,
                       titlefont=font_labels, tickfont=font_ticks, showgrid=TRUE) 
      
      y_breaks <- hist(rev(RV$label_y), plot=FALSE)$breaks
      y_axis   <- list(title=paste("Wavelength (", RV$unities[2], ")", sep=""),
                       tickmode="array", ticktext=y_breaks,
                       tickvals=which(rev(RV$label_y) %in% y_breaks)-1,
                       tickwidth=0.2, ticklen=2, tickangle=0,
                       titlefont=font_labels, tickfont=font_ticks, showgrid=TRUE)
      
      # Plot
      Absorbance <- t(RV$zz_rev)
      p_contour  <- plot_ly(z=~Absorbance, type="contour", colors=plot_colors0,
                           hoverinfo="text", text=hover_surface,
                           showscale=TRUE, source="source_contour") %>% 
        layout(title="<b>Contour plot</b>", titlefont=font_title, margin=marg,
               xaxis=y_axis, yaxis=x_axis, autosize=TRUE, showlegend=FALSE)
      
      # Ignoring explicitly provided widget ID, since Shiny doesn't use them
      p_contour$elementId <- NULL
      p_contour
      
      # Save plot in a reactive object
      Rplot$p_contour <- p_contour
    }
  })

  # Button to townload plot in html
  output$downloadBtnCont <- downloadHandler(
    filename = function() { "plot-contour.html" },
    content  = function(file) {
      p_contour2 <- Rplot$p_contour
      htmlwidgets::saveWidget(widget=p_contour2, file)
    }
  )
  
  #----------------------------------------------------------------------------#
  # TAB: "Optical" (Line graphs of FWHM statistics on same x-axis) #
  #----------------------------------------------------------------#
  output$plot_FWHM <- renderPlotly({
    if(is.null(RV$zz_rev)){
      return()
    } else { 
      plot_colors <- c("#619CFF", "#F8766D", "#00BA38")
      
      # Plot "Max Absorbance" data (blue circles) ------------------------------
      hover_p1 <- paste("<b>Max Absorbance:</b> ", FWHM$uvvis_data$Abs_max, "<br>", 
                        "<b>Time:</b> ", FWHM$uvvis_data$Time, " ", RV$unities[1], sep="")
      p1 <- plot_ly(y=FWHM$uvvis_data$Abs_max, x=FWHM$uvvis_data$Time,
                    mode="lines", type="scatter",
                    hoverinfo="text", text=hover_p1,
                    line=list(color=plot_colors[1], width=2), showlegend=FALSE)
      
      # Add "Max Wavelength" curve
      hover_p2 <- paste("<b>Max Wavelength:</b> ", FWHM$uvvis_data$Wl_max, " ", RV$unities[2], "<br>", 
                        "<b>Time:</b> ", FWHM$uvvis_data$Time, " ", RV$unities[1], sep="")
      p2 <- plot_ly(y=FWHM$uvvis_data$Wl_max, x=FWHM$uvvis_data$Time,
                    mode="lines", type="scatter", hoverinfo="text", text=hover_p2,
                    line=list(color=plot_colors[2], width=2), showlegend=FALSE)
      
      # Add "FWHM" curve
      hover_p3 <- paste("<b>FWHM:</b> ", round(FWHM$uvvis_data$FWHM, 2), " ", RV$unities[2], "<br>", 
                        "<b>Time:</b> ", FWHM$uvvis_data$Time, " ", RV$unities[1], sep="")
      p3 <- plot_ly(y=FWHM$uvvis_data$FWHM, x=FWHM$uvvis_data$Time,
                    mode="lines", type="scatter",
                    hoverinfo="text", text=hover_p3,
                    line=list(color=plot_colors[3], width=2), showlegend=FALSE)
      
      marg   <- list(b=87, l=68, t=65, r=28, pad=0)
      xrange <- c(min(RV$label_x), max(RV$label_x)+5)
      font_title  <- list(size=20, family="Arial", color="black")
      font_labels <- list(family="Arial", color="black")
      font_ticks  <- list(family="Arial", color="black")
      
      p_FWHM <- subplot(p1, p2, p3, nrows=3) %>%
        layout(
          title="<b>Optical Properties</b>", titlefont=font_title,
          margin=marg,
          yaxis=list(title="Absorbance"),
          yaxis2=list(title=paste("SPR peak (", RV$unities[2], ")", sep=""),
                      titlefont=font_labels, tickfont=font_ticks),
          yaxis3=list(title=paste("FWHM (", RV$unities[2], ")", sep=""),
                      titlefont=font_labels, tickfont=font_ticks),
          xaxis=list(range=xrange, titlefont=font_labels, tickfont=font_ticks),
          xaxis2=list(range=xrange, titlefont=font_labels, tickfont=font_ticks),
          xaxis3=list(title=paste("Time (", RV$unities[1], ")", sep=""), range=xrange,
                      titlefont=font_labels, tickfont=font_ticks)
        )
      
      # Ignoring explicitly provided widget ID, since Shiny doesn't use them
      p_FWHM$elementId <- NULL
      p_FWHM
      
      # Save plot in a reactive object
      Rplot$p_FWHM <- p_FWHM
    }
  })
  # Button to townload plot in html
  output$downloadBtnStats <- downloadHandler(
    filename = function() { "plot-optical.html" },
    content  = function(file) {
      p_FWHM2 <- Rplot$p_FWHM
      htmlwidgets::saveWidget(widget=p_FWHM2, file)
    }
  )
  #---------------------------------------------------------------------------#
  # TAB: SPECTRUMS (Indiv. plots of Absorbance vs Wavelengths, for each time) #
  #---------------------------------------------------------------------------#
  # Define number of columns to be displayed
  n_col <- 3
  
  # Define UI for plot spectrum
  output$plot_spectrum_1 <- renderUI({
    if(is.null(RV$zz_rev)){
      return()
    } else {
      # Define general parameters of plot
      col_width <- round(12/n_col)        # Calculate bootstrap column width
      n_row     <- ceiling(RV$nX/n_col) # calculate number of rows
      counter   <<- 0                     # Counter variable
      
      # Create matrix layout scheme of plot tags, based on the number of colums
      rows  <- lapply(1:n_row, function(row_num){
        cols  <- lapply(1:n_col, function(i) {
          counter   <<- counter + 1
          plot_name <- paste("plot", counter, sep="")
          column(width=col_width, plotlyOutput(plot_name))  # , height=280, width=250
        })
        fluidRow( do.call(tagList, cols) )
      })
      do.call(tagList, rows)
    }
  })
  
  
  #----------------------------------------------------------------------------#
  # TAB: "MixSpectrums" (All spectrum together in one interactive plot) #
  #---------------------------------------------------------------------#
  # Define UI for plot spectrum
  output$plot_spectrum_2 <- renderPlotly({
    if(is.null(RV$zz_rev)){
      return()
    }
    if(!is.null(RV$zz_rev) | input$reloadBtnSpec > 0){  
      # Plot color scheme
      plot_colors  <- c("#619CFF", "#F8766D", "#00BA38")
      
      # Create merged data-frame with all information available 
      ZZ_ALL_0 <- data.frame(Time=rep(RV$label_x, each=nrow(RV$zz_rev)),
                             wl=rep(rev(RV$label_y), ncol=RV$zz_rev),
                             absor=c(RV$zz_rev))
      ZZ_ALL   <- merge(ZZ_ALL_0, FWHM$uvvis_data, by="Time")
      head(ZZ_ALL)
      
      # Plot first graph (LEFT)
      sd   <- SharedData$new(ZZ_ALL, ~Time)
      base <- plot_ly(sd, color=I(plot_colors[1])) %>%
        group_by(Time)
      
      #--------------------------------
      # Plot first graph (on the LEFT)
      #--------------------------------
      hover_p_left <- paste("<b>Time:</b> ", FWHM$uvvis_data$Time, " ", RV$unities[1], "<br>",
                            "<b>Max Absorbance:</b> ", FWHM$uvvis_data$Abs_max, "<br>", 
                            "<b>Max Wavelength:</b> ", FWHM$uvvis_data$Wl_max, " ", RV$unities[2], "<br>",
                            "<b>FWHM:</b> ", round(FWHM$uvvis_data$FWHM, 2), " ", RV$unities[2], "", sep="")
      p_left <- base %>%
        summarise(Time2=unique(Time), Abs_max2=unique(Abs_max)) %>%
        add_markers(x=~Abs_max2, y=~Time2, showlegend=FALSE,
                    hoverinfo="text", text=hover_p_left) %>%
        layout(xaxis=list(title="Absorbance"),
               yaxis=list(title=paste("Time (", RV$unities[1], ")", sep=""),
                          linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                          tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE)
        )
      p_left
      
      #----------------------------------
      # Plot second graph (on the RIGHT)
      #----------------------------------
      hover_p_right <- paste("<b>Time:</b> ", ZZ_ALL$Time, " ", RV$unities[1], "<br>",
                             "<b>Wavelength:</b> ", ZZ_ALL$wl, " ", RV$unities[2], "<br>",
                             "<b>Absorbance:</b> ", ZZ_ALL$absor, sep="")
      p_right <- base %>%
        add_lines(x=~wl, y=~absor, alpha=0.3, type='scatter', xaxis="x", yaxis="y",
                  hoverinfo="text", text=hover_p_right,
                  showlegend=TRUE) %>%
        layout(xaxis=list(title=paste("Wavelength (", RV$unities[2], ")", sep=""),
                          linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                          tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE),
               yaxis=list(title="Absorbance", linewidth=1, showline=TRUE)
        )
      p_right
      
      #------------------
      # Plot JOINT graph
      #------------------
      marg <- list(b=87, l=68, t=65, r=28, pad=0)
      
      font_title  <- list(size=20, family="Arial", color="black")
      font_labels <- list(family="Arial", color="black")
      font_ticks  <- list(family="Arial", color="black")
      
      p_joint <- subplot(p_left, p_right, titleX=TRUE, titleY=TRUE,
                         widths=c(0.2, 0.8)) %>%
        layout(title="<b>Mix Spectrum</b>",
               titlefont=font_title,
               margin=marg,
               dragmode="select", showlegend=TRUE,
               yaxis=list(title=paste("Time (", RV$unities[1], ")", sep=""),
                          titlefont=font_labels, tickfont=font_ticks),
               xaxis=list(title="Absorbance",
                          linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                          tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE,
                          titlefont=font_labels, tickfont=font_ticks),
               yaxis2=list(title="Absorbance", side="right",
                           linewidth=1, showline=TRUE, showgrid=TRUE, showticklabels=TRUE,
                           tickangle=0, tickwidth=1, ticklen=6, zeroline=FALSE,
                           titlefont=font_labels, tickfont=font_ticks),
               xaxis2=list(title=paste("Wavelength (", RV$unities[2], ")", sep=""),
                           titlefont=font_labels, tickfont=font_ticks)
        )
      p_spec_mix <- highlight(p_joint, dynamic=FALSE, persistent=TRUE)
      
      # Ignoring explicitly provided widget ID, since Shiny doesn't use them
      p_spec_mix$elementId <- NULL
      p_spec_mix
      
      # Save plot in a reactive object
      Rplot$p_spec_mix <- p_spec_mix
    }
  })
  # Button to townload plot in html
  output$downloadBtnSpec <- downloadHandler(
    filename = function() { "plot-mixSpectrum.html" },
    content  = function(file) {
      p_spec_mix2 <- Rplot$p_spec_mix
      htmlwidgets::saveWidget(widget=p_spec_mix2, file)
    }
  )
}) #endof server{()}

#==============================================================================#
#================================= END ========================================#
#==============================================================================#