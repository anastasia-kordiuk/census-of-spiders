brd <- st_read("D:/ural_border.gpkg", "border")

# landcover ---------------------------------------------------------------
nm <- dir(pattern = "*.tif")
i <- 1

for(i in 1:12){
cat(i, "\n")
rst <- raster(nm[i])
rst_cropped <- crop(rst, extent(brd))
rst_final <- aggregate(rst_cropped, fact=2)
raster::writeRaster(
    rst_final, 
    paste0("export/landcover_", str_extract(nm[i], "class_[:digit:]+"), ".tif")
)
rm(rst_cropped, rst_final)
}
    
# elevation ---------------------------------------------------------------


library(geodata)

geodata::geodata_path("D:/geodata.gis.tmp")

srtm <- geodata::elevation_30s("RU")
kz <- geodata::elevation_30s("KZ")
ru <- srtm


elev <- raster::merge(kz, ru)
plot(elev)
rst_cropped <- crop(elev, extent(brd))
rst_final <- aggregate(rst_cropped, fact=2)
raster::writeRaster(
    rst_final, 
    paste0("tmp.tif")
)



ru <- geodata::worldclim_country("RU", 1)


# WorldClim ---------------------------------------------------------------
nm <- dir(pattern = "*.tif")
i <- 1

for(i in 2:19){
    cat(i, "\n")
    rst <- raster(nm[i])
    rst_cropped <- crop(rst, extent(brd))
    rst_final <- aggregate(rst_cropped, fact=3)
    raster::writeRaster(
        rst_final, 
        paste0("export/", str_remove_all(nm[i], "30s_"))
    )
    rm(rst_cropped, rst_final)
}

