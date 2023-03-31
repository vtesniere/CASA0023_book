library(tidyverse)
library(terra)
library(fs)
# install.packages("Rtools")
library(devtools)
library(glcm)
library(raster)
library(here)
library(dplyr)
library(raster)
# all images were taken between Jan 24 and Feb 22 2022 (7 frames total)
test <- list.files("Week_3/Landsat/")
test[1:5]
list
func <- function(x){
  dir_info(here("Week_3", "Landsat", x))%>%
    filter(str_detect(path, "[B123456790].TIF"))%>%
    select(path)%>%
    pull()%>%
    as.character()%>%
    rast()
}
landsat8_list <- lapply(test[1:5], func)
unlist(landsat8_list)
for (i in length(test[1:5])) {
  landsat8_[i] <- func(test[i])
  i+1
}
getwd()
landsat8_1 <- func( test[1] )
landsat8_2 <- func( test[2] )
landsat8_3 <- func( test[3] )
landsat8_4 <- func( test[4] )
landsat8_5 <- func( test[5] )
landsat9_6 <- func( test[6] )
landsat9_7 <- func( test[7] )

# making a mosaic with all the images

m1 <- mosaic(landsat8_1,landsat8_2, landsat8_3, landsat8_4, landsat8_5,
             landsat9_6, landsat9_7, fun = "mean")
m1_samedate <- mosaic(landsat8_1,landsat8_2, landsat8_3, landsat8_5,
                      landsat9_6, landsat9_7, fun = "mean")
plot(m1, axes=0, sub = NULL)
plot(m1_samedate, axes=0, sub = NULL)
title(main = NULL, sub = NULL)
typeof(m1)

# Using the spectral traits to use the NDVI index
m1_NDVI <- (m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5 - m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4 ) /
  (m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5 + m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4)
# so there are only 8 bands that we chose initially and which i guess were merged during the creation of the mosaic

m1_NDVI %>%
  plot(., axes = 0)
# reclassyfying where NDVI is above 0.1, removing Low NDVI values
veg_no_low <- m1_NDVI %>%
  classify(., cbind(-Inf, 0.1, NA))

veg_medium <- m1_NDVI %>%
  classify(., cbind(-Inf, 0.2, NA))

veg_med_high <- m1_NDVI %>%
  classify(., cbind(-Inf, 0.4, NA))

veg_low_high <- m1_NDVI %>%
  classify(., cbind(-Inf, 0.6, NA))

veg_high <- m1_NDVI %>%
  classify(., cbind(-Inf, 0.9, NA))
par(mfrow = c(2,2))
m1_NDVI %>%
  plot(., axes = 0)
veg_no_low %>%
  plot(., axes = 0)
veg_medium %>%
  plot(., axes = 0)
veg_med_high %>%
  plot(., axes = 0)

max(veg_med_high)

# Using the spectral traits to use the NDMI (moisture) index
# uses bands B5 and B6 on Landsat 8-9
m1_NDMI <- (m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5 - m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B6 ) / 
  (m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5 + m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B6)

m1_NDMI %>%
  plot(., axes = 0)
moisture_positive %>%
  plot(., axes = 0)
moisture_low %>%
  plot(., axes = 0)
moisture_medium %>%
  plot(., axes = 0)

moisture_positive <- m1_NDMI %>%
  classify(., cbind(-Inf, 0, NA))
moisture_low <- m1_NDMI %>%
  classify(., cbind(-Inf, 0.1, NA))
moisture_medium <- m1_NDMI %>%
  classify(., cbind(-Inf, 0.2, NA))
plot(moisture_low, axes=0)
plot(moisture_medium, axes=0)
par(mfrow = c(1,1))
par(mfrow = c(1,2))
plot(moisture_medium, axes=0)
plot(veg_medium, axes=0)
par(mfrow = c(2,1))
plot(m1_filter_SWIR, axes = 0)
plot(m1_filter_SWIR_same, axes = 0, main = "Comparison Between filtering with 6 (below) and 7 (above) frames, \nwhich allows to check for effects of temporal differences")
plot(m1_filter_SWIR, axes = 0)
plot(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B, axes = 0)
# now on to filtering (on a single)
# this is for band 4 (aka red)
m1_filter <- focal(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4, w=matrix(3,3))
m1_filter_green <- focal(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B3, w=matrix(3,3))
m1_filter_blue <- focal(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B2, w=matrix(3,3))
m1_filter_SWIR <- focal(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B6, w=matrix(3,3))
plot(m1_filter_green, axes = 0)
plot(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B3)
# we see that one of the images has outliers (image set 4 which was a couple weeks later than the other 6)
# here we run it again to see if there are different results without said outlier
m1_filter_same <- focal(m1_samedate$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4, w=matrix(3,3))
m1_filter_green_same <- focal(m1_samedate$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B3, w=matrix(3,3))
m1_filter_blue_same <- focal(m1_samedate$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B2, w=matrix(3,3))
m1_filter_SWIR_same <- focal(m1_samedate$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B6, w=matrix(3,3))
plot(m1_filter_green_same, axes = 0)
plot(m1_samedate$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B3)

# it seems mac does not need Rtools as already uses UNIX
library(glcm)
library(raster)

red_raster <- raster(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4)
green_raster <- raster(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B3)
blue_raster <- raster(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B2)
NIR_raster <- raster(m1$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5)

red_texture <- glcm(red_raster,
                    window = c(7,7),
                    # window = c(9,9),
                    # shift = list(c(0,1),c(1,1), c(1,0), c(1,-1)),
                    statistics = c("homogeneity"))
green_texture <- glcm(green_raster,
                      window = c(7,7),
                      # window = c(9,9),
                      # shift = list(c(0,1),c(1,1), c(1,0), c(1,-1)),
                      statistics = c("homogeneity"))

NIR_texture <- glcm(NIR_raster,
                    window = c(7,7),
                    # window = c(9,9),
                    # shift = list(c(0,1),c(1,1), c(1,0), c(1,-1)),
                    statistics = c("mean",
                                   "variance",
                                   "homogeneity", #,
                                   "contrast",
                                   "entropy", 
                                   "dissimilarity",
                                   "second_moment", 
                                   "correlation"))

library(pkgbuild)

red_texture$glcm_homogeneity %>% plot(., axes=0)
NIR_texture%>% plot(., axes=0)
green_texture$glcm_homogeneity %>% plot(., axes = 0)


# changing M1 to raster but not through terra
m1_raster <- stack(m1)
fuse <- stack(m1_raster, red_texture)
Fuse_3_bands <- stack(fuse$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B4,
                      fuse$LC08_L2SP_090088_20220121_20220128_02_T1_SR_B5,
                      fuse$glcm_homogeneity)
# scaling of the pixel values for standardisation
scale_fuse <- scale(Fuse_3_bands)

.toRaster <- function(x) {
  if (inherits(x, "SpatRaster")) {
    p <- crs(x)
    s <- stack(x)
    crs(s) <- p
    return(s)
  } else {
    return(x)
  }
}

paraRasterFun <- function(raster, rasterFun, args = list(), wrArgs = list()){
  if (isTRUE( getOption('rasterCluster'))) {
    do.call("clusterR", args = c(list(x = raster, fun = rasterFun, args=args), wrArgs))
  } else {
    do.call("rasterFun", args=c(raster, args, wrArgs))
  }
}
# trying to locally create rasterPCA
rasterPCA <- function(img, nSamples = NULL, nComp = nlayers(img), spca = FALSE,  maskCheck = TRUE, ...){      
  img <- toRaster(img)
  if(nlayers(img) <= 1) stop("Need at least two layers to calculate PCA.")   
  ellip <- list(...)
  
  ## Deprecate norm, as it has the same effect as spca
  if("norm" %in% names(ellip)) {
    warning("Argument 'norm' has been deprecated. Use argument 'spca' instead.\nFormer 'norm=TRUE' corresponds to 'spca=TRUE'.", call. = FALSE)
    ellip[["norm"]] <- NULL
  }
  
  if(nComp > nlayers(img)) nComp <- nlayers(img)
  
  if(!is.null(nSamples)){    
    trainData <- sampleRandom(img, size = nSamples, na.rm = TRUE)
    if(nrow(trainData) < nlayers(img)) stop("nSamples too small or img contains a layer with NAs only")
    model <- princomp(trainData, scores = FALSE, cor = spca)
  } else {
    if(maskCheck) {
      totalMask <- !sum(calc(img, is.na))
      if(cellStats(totalMask, sum) == 0) stop("img contains either a layer with NAs only or no single pixel with valid values across all layers")
      img <- mask(img, totalMask , maskvalue = 0) ## NA areas must be masked from all layers, otherwise the covariance matrix is not non-negative definite   
    }
    covMat <- layerStats(img, stat = "cov", na.rm = TRUE)
    model  <- princomp(covmat = covMat[[1]], cor=spca)
    model$center <- covMat$mean
    model$n.obs  <- cellStats(!any(is.na(img)), sum)
    if(spca) {    
      ## Calculate scale as population sd like in in princomp
      S <- diag(covMat$covariance)
      model$scale <- sqrt(S * (model$n.obs-1)/model$n.obs)
    }
  }
  ## Predict
  out   <- paraRasterFun(img, rasterFun=raster::predict, args = list(model = model, na.rm = TRUE, index = 1:nComp), wrArgs = ellip)  
  names(out) <- paste0("PC", 1:nComp)
  structure(list(call = match.call(), model = model, map = out), class = c("rasterPCA", "RStoolbox"))  
  
}

pca <- rasterPCA(fuse, 
                 nSamples =100,
                 spca = TRUE)
summary(pca$model)
par(mfrow = c(2,2))
plot(pca$map$PC1, axes = 0)
plot(pca$map$PC2, axes = 0)
plot(pca$map$PC3, axes = 0)
plot(pca$map$PC4, axes = 0)
par(mfrow = c(1,1))
plot(pca$map$PC1, axes = 0)

install.packages("RStoolbox")
install.packages("devtools")
library(devtools)
install_github("bleutner/RStoolbox")
library(RStoolbox)


# below this was just some training/playing around with the data

landsat9_list <- lapply(test[6:7], func)
typeof(landsat8_1)
typeof(practical_way)
m1 <- mosaic(practical_way, practical_way2, fun = "mean")

practical_way <- dir_info(here("Week_3", "Landsat","LC08_L2SP_090088_20220121_20220128_02_T1"))%>%
  filter(str_detect(path, "[B123456790].TIF"))%>%
  select(path)%>%
  pull()%>%
  as.character()%>%
  rast()

practical_way2 <- dir_info(here("Week_3", "Landsat","LC08_L2SP_090089_20220121_20220128_02_T1"))%>%
  filter(str_detect(path, "[B123456790].TIF"))%>%
  select(path)%>%
  pull()%>%
  as.character()%>%
  rast()

# List your raster files excluding band 8 using the patter argument
listlandsat_8<-dir_info(here::here("Week_3", "Landsat", "Lsat8"))%>%
  dplyr::filter(str_detect(path, "[B123456790].TIF")) %>%
  dplyr::select(path)%>%
  pull()%>%
  as.character()%>%
  # Load our raster layers into a stack
  terra::rast()

# List your raster files excluding band 8 using the patter argument
listlandsat_9<-dir_info(here::here("prac_3", "Landsat", "Lsat9"))%>%
  dplyr::filter(str_detect(path, "[1B23456790].TIF")) %>%
  dplyr::select(path)%>%
  pull()%>%
  as.character()%>%
  # Load our raster layers into a stack
  terra::rast()