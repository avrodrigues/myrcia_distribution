map_clade_richness <- function(
    tree_df,
    occ_df,
    sf_polygon,
    color_pkg = "viridis", 
    palette = "viridis", 
    line_color_sf = "#505045", 
    background_color = "#FAF8F4",
    title = NULL,
    plot_theme = NULL
){
  
  require(rcartocolor)
  require(viridis)
  require(ggplot2)
  require(ggmap)
  
  
  # CHECK ARGUMENTS --------------------------------------------------------
  
  
  match.arg(color_pkg, c("viridis", "rcartocolors"))
  
  # |- palettes ----
  if(color_pkg == "viridis"){
    
  viridis_names <- c("cividis", "inferno", "magma", "mako", "plasma", 
                     "rocket", "turbo", "viridis")
  
  match.arg(palette, viridis_names)
  color_idx <- pmatch(palette, viridis_names)
  }
  
  if(color_pkg == "rcartocolors"){
    
    carto_names <- c("TealGrn", "Teal", "SunsetDark", "Sunset", "RedOr",     
    "PurpOr", "Purp", "PinkYl", "Peach", "OrYel", "Mint", "Magenta", 
    "Emrld", "DarkMint" , "BurgYl", "Burg", "BrwnYl", "BluYl", "BluGrn") 
    
    palette = "Emrld"
    
    match.arg(palette, carto_names)
    color_idx <- pmatch(palette, carto_names)
  }
  
  if(is.null(line_color_sf)|is.na(line_color_sf)) {
    stop("Please provide a color to 'line_color'")
    }
  if(is.null(background_color)|is.na(background_color)){
    stop("Please provide a color to 'line_color'")
  }
  
  
  # |- plot theme
  if(is.null(plot_theme)){
    greys <- c(
      "#040400",
      "#1F1F1B",
      "#3B3B37",
      "#575753",
      "#73736F",
      "#8E8E8A",
      "#AAAAA6",
      "#C6C6C2",
      "#E2E2DE",
      "#FEFEFA"
    )
    
    plot_theme <- list(
      theme(
        panel.background = element_rect(fill = background_color), 
        panel.grid = element_blank(), 
        text = element_text(color = greys[1]), 
        title = element_text(color = greys[2]),
        axis.text = element_text(color = greys[2]), 
        axis.ticks = element_line(color = greys[3]), 
        panel.border = element_rect(color = greys[4], fill = NA), 
        axis.title = element_blank()
      )
    )
  }
  else{
    if(!inherits(plot_theme, "list")) stop("'plot_theme' must be of class 'list'")
  }
  
  
  # |- prepare data ----
  rich_df <- calc_richness(tree_df, occ_df)
  max_rich <- max(rich_df$richness)
  limits <- xy_limits(occ_df)
  
  # CREATE GGPLOT OBJECT ----------------------------------------------------
  
  # |- base layer ----
  gg_base <- 
    ggplot() +
    geom_raster(data = rich_df, aes(x, y, fill = richness)) 
  
  # |- fill options ----
  
  if(color_pkg == "viridis"){
  gg_fill <- 
    scale_fill_viridis_c(
      name = "Richness",
      option = palette, 
      limits = c(0, max_rich), 
      direction = 1) 
  }
  
  if(color_pkg == "rcartocolors"){
    gg_fill <- 
      scale_fill_carto_c(
        name = "Richness",
        palette  = palette, 
        limits = c(0, max_rich),
        direction = 1) 
  }
  
  # |- add sf layer ----
  
  gg_sf <- 
    gg_base + gg_fill +
    geom_sf(data = sf_polygon, fill = NA, colour  = line_color_sf) +
    coord_sf(xlim = limits$x, ylim = limits$y) 
  
  # |- theme options
  gg_out <-gg_sf + 
    theme_bw() +
    ggtitle(title) +
    plot_theme
  
  gg_out
  
}
  
