Variations to Consider
======================

Effect of projections on metrics
Effect of simplificaton on metrics
Effect of simple vs non-simple vs holes on metrics


http://www.electiondataarchive.org/datacenter-gred.html

ogr2ogr -f 'ESRI Shapefile' -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs' tl_2016_us_cd115.shp /z/conv.shp

./dgfinder.exe /z/conv.shp conv /z/out

ogr2ogr -f "ESRI Shapefile" -s_srs '+proj=lcc +lat_1=41.71666666666667 +lat_2=42.68333333333333 +lat_0=41 +lon_0=-71.5 +x_0=200000 +y_0=750000 +datum=NAD83 +units=m +no_defs ' -t_srs "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"  /z/conv.shp CONGRESSMA_POLY.shp
