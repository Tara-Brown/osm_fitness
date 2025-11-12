import geopandas as gpd
from shapely.validation import make_valid

gpkg_path = "unified_geopackages/state_06_unified.gpkg"
gdf = gpd.read_file(gpkg_path, layer="unified_park_area")

# Ensure CRS is equal-area for area calculation
gdf = gdf.to_crs("EPSG:5070")

# Fix invalid geometries
gdf["geometry"] = gdf.geometry.apply(lambda g: make_valid(g) if g else g)

# Compute areas (square meters)
gdf["area_m2"] = gdf.geometry.area

print(gdf[["city_name", "area_m2_parks"]].head())
print("Total area:", gdf["area_m2_parks"].sum())
print(gdf.columns)
