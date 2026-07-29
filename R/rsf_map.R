rsf_map <- function(rsfmap, aoi, classify,
                         class_breaks, title = "Habitat Suitability") {
  stopifnot("`rsfmap` must be a SpatRaster object"
            = class(rsfmap) == "SpatRaster",
            "`rsfmap` layer must have values between 0-1"
            = (round(terra::minmax(rsfmap)[1], 0) >= 0
               && round(terra::minmax(rsfmap)[2], 0) <= 1),
            "`rsfmap` layer must have a CRS defined"
            = terra::crs(rsfmap) != "")

  names(rsfmap) <- "rsfmap"
  exp <- rsfmap

  if (!missing(aoi)) {
    stopifnot("`aoi` must be a SpatVector object"
              = class(aoi) == "SpatVector",
              "`aoi` extent must be within `rsfmap` extent"
              = terra::relate(aoi, exp, "within"),
              "`rsfmap` and `aoi` must have same CRS"
              = terra::same.crs(exp, aoi))

    aoi_r <- terra::rasterize(aoi, exp)

    exp_aoi <- exp * aoi_r

    exp <- terra::trim(exp_aoi, padding = 10)
  }


  if (!missing(classify)) {
    stopifnot("`classify` must be one of: 'local', 'landscape', or 'custom'"
              = classify %in% c("local", "landscape", "custom"))

    if (classify == "custom") {
      stopifnot("must provide 'class_breaks' if `classify` = 'custom'"
                = !missing(class_breaks))
    }

    if (classify == "landscape") {
      class_breaks <- c(0.2, 0.4, 0.6, 0.8, 1)
    }

    if (classify == "local") {
      class_breaks <- c(0.15, 0.3, 0.45, 1)
    }

    class_breaks <- sort(class_breaks)

    # class_breaks checks
    stopifnot("`class_breaks` must be a vector of numbers"
              = class(class_breaks) == "numeric",
              "`class_breaks` must have 1 as the maximum value"
              = max(class_breaks) == 1,
              "`class_breaks` must be greater than 0"
              = class_breaks > 0)

    class_labels <- character()

    label_breaks <- c(0, class_breaks)
    for (i in seq_along(label_breaks)) {
      class_labels[i] <- paste(label_breaks[i], "-", label_breaks[i + 1])
    }

    class_labels <- c("Nil", utils::head(class_labels, -1))

    lut <- data.frame(start = c(0, 0, utils::head(class_breaks, -1)),
                      end = c(0, class_breaks),
                      factor = 0:length(class_breaks),
                      label = class_labels)

    rcmats <- as.matrix(lut[, 1:3])
    exp <- terra::classify(exp, rcmats, include.lowest = TRUE)

    levels(exp) <- lut[, 3:4]

    n_color <- length(class_breaks)

    cols <- c("grey40", tidyterra::whitebox.colors(n_color,
                                                   palette = "bl_yl_rd"))
    names(cols) <- lut$label


    plt <- tmap::tm_shape(exp) +
      tmap::tm_raster(
        col.scale = tmap::tm_scale_categorical(values = cols),
        col_alpha = 0.7,
        col.legend = tmap::tm_legend(title = "Probability")
      )
  } else {
    cols <- tidyterra::whitebox.colors(n = 100, palette = "bl_yl_rd")

    plt <- tmap::tm_shape(exp) +
      tmap::tm_raster(
        col.scale = tmap::tm_scale_continuous(values = cols, limits = c(0, 1)),
        col_alpha = 0.7,
        col.legend = tmap::tm_legend(title = "Probability",
                                     position = tmap::tm_pos_out(
                                       "right"))
      )
  }

  if (!missing(aoi)) {
    plt <- plt + tmap::tm_shape(sf::st_as_sf(aoi)) +
      tmap::tm_borders(lwd = 2)
  }


  plt <- plt + tmap::tm_basemap("Esri.WorldImagery") +
    tmap::tm_credits("Basemap: Esri World Imagery",
                     color = "white") +
    tmap::tm_title(title) +
    tmap::tm_compass(
      type = "arrow",
      position = c("LEFT", "BOTTOM"),
      color.light = "black",
      color.dark = "white",
      text.color = "white"
    ) +
    tmap::tm_scalebar(position = tmap::tm_pos_out("center", "bottom"),
                      text.size = 0.9) +
    tmap::tm_layout(inner.margins = 0.05) +
    tmap::tm_crs("auto")

  return(plt)
}