
# Packages

library(terra)   
library(sf)      
library(dplyr)  
library(tidyr)
library(stringr)
library(spatialEco)


############################################################
#  Function to process a single DEM and RGB image
############################################################

process_rgb_dem <- function(rgb_path, dem_path, grid_path, target_res, out_dir) {
  
  message("Processing: ", basename(dem_path))
  message("Processing: ", basename(rgb_path))
  
  #Load RGB
  rgb_data <- rast(rgb_path)
  plotRGB(rgb_data, r = 1, g = 2, b = 3)
  
  # Load DEM
  dem <- rast(dem_path)
  plot(dem)
  
  #Load the grid
  grid_sf   <- st_read(grid_files[5], quiet = TRUE)
  grid_vect <- vect(grid_sf)
  plot(grid_vect)
  text(grid_vect, label = grid_vect$id)
  
  # Project grid if needed
  # if (print(are_equal) == FALSE) {
  #   #print('yes')
  #   grid_vect <- project(grid_vect, crs(dem))
  # }
  
  # CRS check
  # rgb_crs <- sf::st_crs(rgb_data)
  # dem_crs <- sf::st_crs(dem)
  # vector_crs <- sf::st_crs(grid_vect)
  # are_equal <- rgb_crs == vector_crs
  # are_equal2 <- dem_crs == vector_crs
  # 
  # #print
  # print("RGB and grids are in the same CRS: ", are_equal)
  # print("DEM and grids are in the same CRS: ", are_equal2)
  # 
  
  # Create target template raster for RGB
  # rgb_template <- rast(
  #   ext(rgb_data),
  #   resolution = target_res,
  #   crs = crs(rgb_data)
  # )
  # 
  # # Resample RGB
  # rgb_rs <- resample(rgb_data, rgb_template, method = "bilinear")
  # plot(rgb_rs)
  
  
  ########################################################################################################## Orthomosaic Analysis
  
  blue  <- rgb_data[[1]]
  green <- rgb_data[[2]]
  red   <- rgb_data[[3]]
  #nir   <- rgb_path[[4]]
  
  # ---- Indices ----
  rb <- red - blue
  rg <- red - green
  gr = (green / red)
  rgb_intensity = (red + green + blue) / 3
  nrg = (red - green) / (red + green)
  nrb = (red - blue) / (red + blue)
  exg  = (green * 2) - (red - blue)
  vari = (green - red) / (green + red - blue)
  ci = (red - blue) / red
  rr = red / (red + green + blue)
  gr = green / (red + green + blue)
  br = blue / (red + green + blue)
  
  #ndvi <- (nir - red) / (nir + red)
  #gndvi <- (nir - green) / (nir + green)
  #ndwi <- (green - nir) / (green + nir)
  
  message("Finished calculating RGB Indices ")
  
  # Extraction
  blue_stats <- terra::extract(blue, grid_vect, fun = c("mean", "sd"),na.rm = TRUE)
  green_stats <- terra::extract(green, grid_vect, fun = c("mean", "sd"),na.rm = TRUE)
  red_stats <- terra::extract(red, grid_vect, fun = c("mean", "sd"),na.rm = TRUE)
  rgb_stats  <- terra::extract(rgb_intensity, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  nrg_stats  <- terra::extract(nrg, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  nrb_stats <- terra::extract(nrb, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  exg_stats <- terra::extract(exg, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  vari_stats <- terra::extract(vari, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  ci_stats <- terra::extract(ci, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  rr_stats <- terra::extract(rr, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  rr_stats <- terra::extract(rr, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  gr_stats <- terra::extract(gr, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  br_stats <- terra::extract(br, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  
  
  ######################################################################################################### DEM Analysis 
  
  # Create target template raster for DEM
  dem_template <- rast(
    ext(dem),
    resolution = target_res,
    crs = crs(dem)
  )
  
  # Resample DEM
  dem_rs <- resample(dem, dem_template, method = "bilinear")
  plot(dem_rs)
  
  
  ##########################################################################################################
  
  # Terrain metrics
  # Circular roughness: SD of elevation within 10-cell radius
  # custom_w <- focalMat(x = dem_rs, 1.5, type = "circle" )
  # rough <- focal(dem_rs, w = custom_w, fun = sd, na.rm = TRUE)
  # plot(rough)
  ##########################################################################################################
  
  slope <- terrain( dem_rs, neighbors=8, v = "slope",unit = "degrees")   # Slope (degrees)
  aspect <- terrain(dem_rs, neighbors=8, v = "aspect", unit = "degrees") # Aspect
  triRMSD <- terrain(dem_rs, neighbors=8, v = "TRIrmsd")   # TRIrmsd
  flowdir <- terrain(dem_rs, neighbors=8, v = "flowdir")
  tpi <- terrain(dem_rs, neighbors=8, v = "TPI")
  #vrm_raster <- vrm(dem_rs, s = 7)
  ##########################################################################################################
  
  ecoplan_curv <- curvature(x = dem_rs, type = "planform")
  ecoprof_curv <- curvature(x = dem_rs, type = "profile")
  ecototal_curv <- curvature(x = dem_rs, type = "total")
  
  ##########################################################################################################
  # Zonal statistics
  elev_stats  <- terra::extract(dem_rs, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  tpi_stats <- terra::extract(tpi,  grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  slope_stats <- terra::extract(slope, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  aspect_stats <- terra::extract(aspect, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  triRMSD_stats <- terra::extract(triRMSD, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  flowdir_stats <- terra::extract(flowdir, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  #total_curv_stats <- terra::extract(profile_curv, grid_vect, fun = "mean", na.rm = TRUE)
  ecoplan_curv <- terra::extract(ecoplan_curv, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  ecoprof_curv <- terra::extract(ecoprof_curv, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  ecototal_curv <- terra::extract(ecototal_curv, grid_vect, fun = c("mean", "sd"), na.rm = TRUE)
  
  
  message("Processing: Extracting values ", basename(dem_path))
  
  # Build output table
  out_df <- data.frame(
    grid_id     = elev_stats$ID,
    mean_elev   = elev_stats$mean,
    mean_flowdir  = flowdir_stats$mean,
    mean_slope  = slope_stats$mean,
    mean_aspect  = aspect_stats$mean,
    mean_tpi  = tpi_stats$mean,
    mean_triRMSD  = triRMSD_stats$mean,
    mean_ecoplan  = ecoplan_curv$mean,
    mean_eco_prof = ecoprof_curv$mean,
    mean_eco_total = ecoprof_curv$mean,
    mean_blue = blue_stats$mean,
    mean_green = green_stats$mean,
    mean_reds = red_stats$mean,
    mean_rgb = rgb_stats$mean,
    mean_nrg = nrg_stats$mean,
    mean_nrb = nrb_stats$mean,
    mean_exg = exg_stats$mean,
    mean_vari = vari_stats$mean,
    mean_ci = ci_stats$mean,
    mean_rr = rr_stats$mean,
    mean_rr = rr_stats$mean,
    mean_gr = gr_stats$mean,
    mean_br = br_stats$mean,
    # sd
    sd_blue = blue_stats$sd,
    sd_greens = green_stats$sd,
    sd_red = green_stats$sd,
    sd_rgb  = rgb_stats$sd,
    sd_nrg = nrg_stats$sd,
    sd_nrb = nrb_stats$sd,
    sd_exg = exg_stats$sd,
    sd_vari = vari_stats$sd,
    sd_ci = ci_stats$sd,
    sd_rr = rr_stats$sd,
    sd_rr = rr_stats$sd,
    sd_gr = gr_stats$sd,
    sd_br = br_stats$sd,
    sd_elev     = elev_stats$sd,
    sd_flowdir  = flowdir_stats$sd,
    sd_slope  = slope_stats$sd,
    sd_triRMSD  = triRMSD_stats$sd,
    sd_ecoplan  = ecoplan_curv$sd,
    sd_eco_prof = ecoprof_curv$sd,
    sd_eco_total = ecoprof_curv$sd,
    sd_aspect  = aspect_stats$sd,
    sd_tpi  = tpi_stats$sd)
  
  # Add DEM identifier
  out_df$dem_name <- tools::file_path_sans_ext(basename(dem_path))
  
  # Save DEM output
  out_csv <- file.path(
    out_dir,
    paste0(out_df$dem_name[1] , " _grid_metrics.csv")
  )
  
  message("Processing: Saving file ", basename(dem_path))
  
  write.csv(out_df, out_csv)
  
  return(out_df) 
  
  }


############################################################


# Dataset links
rgb_folder <- "C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/chapter_2_data_by_structure/all_drone_data/drone_data_2024/10by10_drone_data_2024"

dem_folder <- "C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/chapter_2_data_by_structure/all_drone_data/drone_data_2024/10by10_drone_data_2024"

# grid_folder <- "C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/chapter_new_analysis/roughness/new_roughness_analysis_jan15/"

grid_folder <- "C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/ipad_data_2024/Amabalis/SEG_results/"

#  Output folder
output_dir <- 'C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/chapter_2_data_analysis/chapter_2_obj2/data_extract'

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
# Extracting DEM, RGB, and GRIDS

############################################################
#  List RGB
############################################################

rgb_files <- list.files(path = rgb_folder, pattern = "\\_RGB.tif$",full.names = TRUE)
rgb_files

# Testing: Check the data
# rgb_raster <- rast(rgb_files[1])
# plotRGB(rgb_raster, r = 1, g = 2, b = 3)


############################################################
#  List DEMs
############################################################

dem_files <- list.files(path = dem_folder, 
                        pattern = "\\DEM.tif$", 
                        full.names = TRUE)
dem_files
# Testing: check the dem
# plot(rast(dem_files[1]))

############################################################
#  List Grids
############################################################

grid_files <- list.files(path = grid_folder, 
                         pattern = "\\_5m_Grids.shp$", 
                         full.names = TRUE)
grid_files

#check the data
grid_sf   <- st_read(grid_files[1], quiet = TRUE)
grid_vect <- vect(grid_sf)
plot(grid_vect)
############################################################
#  Apply the function
############################################################

process_rgb_dem(rgb_path = rgb_files[3], 
            dem_path = dem_files[3], 
            grid_path = grid_files[3], 
            target_res = 0.15, 
            out_dir = output_dir)

process_rgb_dem(rgb_path = rgb_files[5], 
                dem_path = dem_files[5], 
                grid_path = grid_files[5], 
                target_res = 0.15, 
                out_dir = output_dir)









