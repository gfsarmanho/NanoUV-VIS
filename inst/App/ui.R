#===================#
# App:  NanoUV-VIS  #
#                   #
# File: ui.R        #
#===================#

# Install required packages
#install.packages(c("shiny", "shinyjs", "knitr", "rmarkdown", "shinycssloaders"))
#install.packages(c("plotly", "plot3D", "DT", "crosstalk", "htmlwidgets"))

# Load packages
# library(shiny)
# library(shinyjs)
# library(knitr)
# library(rmarkdown)
# library(plotly)           # interactive plots (surface, contour, curves, etc.)
# library(plot3D)           # create fancy surface plots
# library(DT)               # create interactive tables
# library(crosstalk)        # generate interactions between widgets
# library(shinycssloaders)  # create "loading signals" when processing plots
# library(htmlwidgets)      # save plots in ".html"

if (!(require(shiny))){install.packages("shiny"); library(shiny, quietly = TRUE)} 
if (!(require(shinyjs))){install.packages("shinyjs"); library(shinyjs, quietly = TRUE)} 
if (!(require(knitr))){install.packages("knitr"); library(knitr, quietly = TRUE)} 
if (!(require(rmarkdown))){install.packages("rmarkdown"); library(rmarkdown, quietly = TRUE)} 
if (!(require(plotly))){install.packages("plotly"); library(plotly, quietly = TRUE)} 
if (!(require(plot3D))){install.packages("plot3D"); library(plot3D, quietly = TRUE)} 
if (!(require(DT))){install.packages("DT"); library(DT, quietly = TRUE)} 
if (!(require(crosstalk))){install.packages("crosstalk"); library(crosstalk, quietly = TRUE)} 
if (!(require(shinycssloaders))){install.packages("shinycssloaders"); library(shinycssloaders, quietly = TRUE)} 
if (!(require(htmlwidgets))){install.packages("htmlwidgets"); library(htmlwidgets, quietly = TRUE)} 

#================#
# Begin Shiny UI #
#================#
shinyUI(navbarPage(
  # Some options
  title       = "NanoUV-VIS",
  windowTitle = "NanoUV-VIS",
  id          = "mainNavBar",
  selected    = "dataTab",
  fluid       = TRUE,
  inverse     = TRUE,             # Changes the (default) navbar color

  # "Load" shinyjs library inside application
  shinyjs::useShinyjs(),
  
  # Load JS code 
  includeScript("www/shinyjs-funcs.js"), # button to stop computations

  # Load CSS functions
  includeCSS("www/style.css"),
  
  #==================#
  # NAV PANEL: About #
  #==================#
  tabPanel(
    title = "About",
    icon  = icon("info-circle"),
    id    = "aboutTab",
    value = "aboutTab",
    name  = "aboutTab",
    
    # Load Rmarkdown data with "about" info
    tags$iframe(src = "ABOUT.html",       # put ".html" file on "www" folder
                width="100%", height="800px", seamless="yes",
                frameborder=0, scrolling="auto"
    )
    
  ), # end of "NAV PANEL: About"
  
  #=================#
  # NAV PANEL: Data # Contains two panels on sidebarPanel and a table on mainPanel
  #=================#
  tabPanel(
    title = "Data",
    icon  = icon(name="file-text", lib="font-awesome"),
    id    = "dataTab",
    value = "dataTab",
    
    # Panel layout (left/sidebarPanel=inputs, right/mainPanel = outputs) 
    sidebarLayout(
      sidebarPanel(
        
        # Organize "Data" navbar as two tabs
        tabsetPanel(
          id   = "dataTabs",
          type = "tabs",
          
          #-------------------#
          # Tab "Single File" # upload a pre-formated data base 
          #-------------------#
          tabPanel(
            title = "Single File",
            id    = "singleFileTab",
            br(),
            downloadLink('downloadExampleData', "Example - Single data file"),
            br(),
            # Button: "load table"
            fileInput(inputId="file1", multiple=FALSE,
                      label=div(h4("Choose unique (.csv) file",
                                   helpPopup("Choose a pre-formated file."))),
                      accept=c("text/csv", "text/comma-separated-values",
                               "text/plain", ".csv")),
            fluidRow(
              column(width=5, 
                     disabled(
                       withBusyIndicatorUI(
                         button=actionButton(inputId="loadBtn",
                                             label="Load table",
                                             class="btn-primary")
                       )
                     )
              ),
              column(width=5,
                     shinyjs::hidden(
                       div(id="resetBtnSINGLE", style="display:inline-block",
                           actionButton(inputId="resetBtnSINGLE",
                                        label="Reset data"),
                           style="float:right")
                     )
              )
            ),
            br(), br(),
            h4("Define coordinates unities labels:"),
            fluidRow(
              column(width=4,
                     textInput(inputId="unit_x", value="min",
                               label=div(h5("Time", helpPopup("X-axis (abscissa)")))
                     )
              ),
              column(width=4, 
                     textInput(inputId="unit_y", value="nm",
                               label=div(h5("Wavelength", helpPopup("Y-axis (ordinate)")))
                     )
              ),
              column(width=4, 
                     textInput(inputId="unit_z", value="",
                               label=div(h5("Absorbance", helpPopup("Z-axis (applicate)")))
                     )
              )
            ),
            # NOT USED
            shinyjs::hidden(
              div(id="RunBtn",
                  h4("Press button to generate graphs and tables"),
                  actionButton(inputId="submitBtn", label="Submit calculations", class="btn-primary")
              )
            )
          ), 
          
          #----------------------#
          # Tab "Multiple Files" # upload files from equipaments 
          #----------------------#
          tabPanel(
            title = "Multiple Files",
            id    = "multipleFilesTab",
            
            br(),
            downloadLink(outputId="downloadExampleMult",
                         label = "Example - Multiple data files"),
            br(),
            
            # Button: "create table"
            fileInput(inputId="filesN", multiple=TRUE,
                      label=div(h4("Choose multiple (.csv) files",
                                   helpPopup("Those files will create the final dataset"))),
                      accept=c("text/csv", "text/comma-separated-values",
                               "text/plain", ".csv")),
            fluidRow(
              column(width=4, 
                     disabled(
                       withBusyIndicatorUI(
                         button=actionButton(inputId="createBtn",
                                             label="Create table",
                                             class="btn-primary"))
                     )
              ),
              column(width=4,
                     shinyjs::hidden(
                       div(id="resetBtnMULT", style="display:inline-block",
                           actionButton(inputId="resetBtnMULT",
                                        label="Reset data"),
                           style="float:right"
                       )
                     )
              ),
              column(width=4,
                     shinyjs::hidden(
                       div(id="down_button_data2",
                           downloadButton(outputId="down_button_data",
                                          label="Download"))
                     )
              )
            ),
            br(), br(),
            h4("Define coordinates unities labels:"),
            fluidRow(
              column(width=4,
                     textInput(inputId="unit_x", value="min",
                               label=div(h5("Time", helpPopup("X-coor.")))
                     )
              ),
              column(width=4, 
                     textInput(inputId="unit_y", value="nm",
                               label=div(h5("Wavelength", helpPopup("Y-coor.")))
                     )
              ),
              column(width=4, 
                     textInput(inputId="unit_z", value="",
                               label=div(h5("Absorbance", helpPopup("Z-coord.")))
                     )
              )
            )
          )
          
        )
      ), #end of sidebarPanel()
      mainPanel(
        div(id="message_dat", h3("Please load a valid data set."),
            style="text-align: center"),
        DT::dataTableOutput(outputId="tab_zzTable")
      ) #end of mainPanel()
    ) #end of sidebarLayout()
  ),
  
  #==========================#
  # NAV PANEL: Visualization #
  #==========================#
  tabPanel(
    title = "Visualization",
    icon  = icon(name="bar-chart", lib="font-awesome"),
    id    = "visualTab",
    value = "visualTab",
    
    div(id="message_vis", h3("Please load a valid data set."),
        style="text-align: center"),
    
    shinyjs::hidden(div(id="visualTab2",
      tabsetPanel(
        tabPanel(title="Surface",
                 div(id="down_button_surface",
                     downloadButton(outputId="downloadBtnSurf",
                                    label="Download Plot in html")),
                 withSpinner(color=plot_colors[1], color.background="white", type=3,
                             plotlyOutput(outputId="plot_surface", height="600px", width="auto"))
                 #verbatimTextOutput("hover_surface"),
                 #verbatimTextOutput("click_surface")
        ),
        tabPanel(title="Contour",
                 div(id="down_button_contour",
                     downloadButton(outputId="downloadBtnCont",
                                    label="Download Plot in html")),
                 withSpinner(color=plot_colors[1], color.background="white", type=3,
                   plotlyOutput(outputId="plot_contour", height="800px", width="auto")
                 )
                 #verbatimTextOutput("hover_contour"),
                 #verbatimTextOutput("click_contour")
        ),
        tabPanel(title="Spectrum",
                 withSpinner(color=plot_colors[1], color.background="white", type=3,
                             uiOutput(outputId="plot_spectrum_1")
                 )
        ),
        tabPanel(title="MixSpectrum", 
                 div(id="down_button_spec",
                     downloadButton(outputId="downloadBtnSpec",
                                    label="Download Plot in html")),
                 div(id="reload_button_spec",
                     actionButton(inputId="reloadBtnSpec",
                                    label="Reload Plot")),
                 withSpinner(color=plot_colors[1], color.background="white", type=3,
                             plotlyOutput(outputId="plot_spectrum_2", height="800px", width="auto")
                 )
        ),
        tabPanel(title="Optical",
                 div(id="down_button_stats",
                     downloadButton(outputId="downloadBtnStats",
                                    label="Download Plot in html")),
                 withSpinner(color=plot_colors[1], color.background="white", type=3,
                             plotlyOutput(outputId="plot_FWHM", height="800px", width="auto")
                 )
        )
      )
    ))
  ), 
  
  #====================#
  # NAV PANEL: Results #
  #====================#
  tabPanel(
    title = "Results",
    icon  = icon(name="table", lib="font-awesome"),
    id    = "resultsTab",
    value = "resultsTab",
    
    div(id="message_res", h3("Please load a valid data set."),
        style="text-align: center"),
    
    shinyjs::hidden(
      div(id="down_button_FWHM",
          downloadButton(outputId="downloadBtnFWHM",
                         label="Download Table")
      )
    ),
    br(), 
    DT::dataTableOutput(outputId="tab_resultsFWHM")
  )
)) #endof navbarPage() and shinyUI()

#==============================================================================#
#================================= END ========================================#
#==============================================================================#