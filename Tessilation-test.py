import math
import geopandas as gpd
import pandas as pd
import osmnx as ox
from shapely.geometry import Polygon
from pyproj import CRS
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, vectors
import rpy2.robjects as robj
import os
import requests

# --- CONFIG ---
PLACE_BASE_URL = "https://www2.census.gov/geo/tiger/TIGER2025/PLACE/tl_2025_{statefp}_place.zip"
CACHE_DIR = "cache/place"
os.makedirs(CACHE_DIR, exist_ok=True)

STATE_FIPS = {
    "01": "Alabama", "02": "Alaska", "04": "Arizona", "05": "Arkansas",
    "06": "California", "08": "Colorado", "09": "Connecticut", "10": "Delaware",
    "12": "Florida", "13": "Georgia", "16": "Idaho", "17": "Illinois",
    "18": "Indiana", "19": "Iowa", "20": "Kansas", "21": "Kentucky",
    "22": "Louisiana", "23": "Maine", "24": "Maryland", "25": "Massachusetts",
    "26": "Michigan", "27": "Minnesota", "28": "Mississippi", "29": "Missouri",
    "30": "Montana", "31": "Nebraska", "32": "Nevada", "33": "New Hampshire",
    "34": "New Jersey", "35": "New Mexico", "36": "New York", "37": "North Carolina",
    "38": "North Dakota", "39": "Ohio", "40": "Oklahoma", "41": "Oregon",
    "42": "Pennsylvania", "44": "Rhode Island", "45": "South Carolina",
    "46": "South Dakota", "47": "Tennessee", "48": "Texas", "49": "Utah",
    "50": "Vermont", "51": "Virginia", "53": "Washington", "54": "West Virginia",
    "55": "Wisconsin", "56": "Wyoming"
}

# -------------------------
# Load city boundaries
# -------------------------
def load_state_places(statefp):
    zip_path = os.path.join(CACHE_DIR, f"tl_2025_{statefp}_place.zip")
    if not os.path.exists(zip_path):
        url = PLACE_BASE_URL.format(statefp=statefp)
        print(f"Downloading {STATE_FIPS[statefp]} PLACE file...")
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with open(zip_path, "wb") as f:
            f.write(r.content)
    else:
        print(f"Using cached file for {STATE_FIPS[statefp]}")
    gdf = gpd.read_file(f"zip://{zip_path}")
    print(f"State shapefile loaded: {gdf.shape[0]} rows")
    return gdf

def get_city_boundary(city_name, state):
    city_name_lower = city_name.lower()
    state_lookup = {v.lower(): k for k, v in STATE_FIPS.items()}
    if state.isdigit():
        statefp = state.zfill(2)
    else:
        statefp = state_lookup.get(state.lower())
        if not statefp:
            raise ValueError(f"Unknown state: {state}")
    gdf = load_state_places(statefp)
    city_rows = gdf[gdf["NAME"].str.lower() == city_name_lower].copy()
    if city_rows.empty:
        raise ValueError(f"City '{city_name}' not found in {STATE_FIPS[statefp]}")
    print(f"City '{city_name}' found: {len(city_rows)} rows")
    return city_rows

# -------------------------
# UTM CRS helper
# -------------------------
def _get_utm_crs_from_lonlat(lon, lat):
    zone = int((lon + 180) / 6) + 1
    south = lat < 0
    crs_str = f"+proj=utm +zone={zone}" + (" +south" if south else "") + " +datum=WGS84 +units=m +no_defs"
    crs = CRS.from_proj4(crs_str)
    print(f"Computed UTM CRS: {crs}")
    return crs

# -------------------------
# Hex grid helpers
# -------------------------
def make_hexagon(cx, cy, s):
    angles = [0, 60, 120, 180, 240, 300]
    return Polygon([(cx + s*math.cos(math.radians(a)), cy + s*math.sin(math.radians(a))) for a in angles])

def make_hex_grid(polygon_gdf, hex_area_m2=1_000_000, buffer_factor=0.05):
    s = math.sqrt((2*hex_area_m2)/(3*math.sqrt(3)))
    hex_height = math.sqrt(3)*s
    hex_width = 2*s
    x_step = 1.5*s
    y_step = hex_height
    minx, miny, maxx, maxy = polygon_gdf.total_bounds
    minx -= (maxx-minx)*buffer_factor + hex_width
    maxx += (maxx-minx)*buffer_factor + hex_width
    miny -= (maxy-miny)*buffer_factor + hex_height
    maxy += (maxy-miny)*buffer_factor + hex_height

    hexes = []
    col = 0
    x = minx
    while x <= maxx + x_step:
        y_offset = 0.5*y_step if col % 2 == 1 else 0.0
        y = miny - y_offset
        while y <= maxy + y_step:
            hexes.append(make_hexagon(x, y, s))
            y += y_step
        x += x_step
        col += 1

    hex_gdf = gpd.GeoDataFrame({"geometry": hexes}, crs=polygon_gdf.crs)
    union_poly = polygon_gdf.unary_union
    hex_gdf = hex_gdf[hex_gdf.intersects(union_poly)].reset_index(drop=True)
    hex_gdf["hex_id"] = hex_gdf.index
    print(f"Hex grid created: {len(hex_gdf)} hexes after clipping")
    return hex_gdf

# -------------------------
# Master workflow with debug
# -------------------------
def hex_tessellate_city_with_grts(city_name, state, hex_area_km2=1.0, sample_frac=0.2, random_seed=42):
    print(f"--- Processing {city_name}, {state} ---")

    # Load city
    city = get_city_boundary(city_name, state).to_crs(epsg=4326)
    centroid = city.geometry.unary_union.centroid
    print(f"City centroid: {centroid.x}, {centroid.y}")
    utm_crs = _get_utm_crs_from_lonlat(centroid.x, centroid.y)
    city_utm = city.to_crs(utm_crs)

    # Hex tessellation
    hex_area_m2 = hex_area_km2*1_000_000
    hexes_utm = make_hex_grid(city_utm, hex_area_m2=hex_area_m2)
    hexes_clipped_utm = hexes_utm[hexes_utm.intersects(city_utm.unary_union)].reset_index(drop=True)

    # Prepare data for R
    hexes_r = hexes_clipped_utm[["hex_id"]].copy()
    hexes_r["geometry_wkt"] = hexes_clipped_utm.geometry.apply(lambda g: g.wkt)
    print("Sample hexes DF for R:", hexes_r.head())

    n_base_val = max(1, int(len(hexes_clipped_utm)*sample_frac))
    print(f"n_base={n_base_val}, seed={random_seed}, utm_epsg={utm_crs.to_epsg()}")

    with pandas2ri.converter.context():
        ro.globalenv["hexes_df"] = pandas2ri.py2rpy(hexes_r)
        ro.globalenv["n_base"] = vectors.IntVector([n_base_val])
        ro.globalenv["seed"] = vectors.IntVector([random_seed])
        ro.globalenv["utm_epsg"] = vectors.IntVector([utm_crs.to_epsg()])
        print("utm_crs EPSG:", utm_crs.to_epsg())


    # --- R code with tryCatch ---
    r_code = """
    library(sf)
library(spsurvey)

set.seed(seed[1])

tryCatch({
    cat("R: Converting geometry WKT to sfc...\n")
    hexes_df$geometry <- st_as_sfc(hexes_df$geometry_wkt)
    
    hexes_sf <- st_sf(hexes_df, crs = utm_epsg[1])
    cat("R: Checking hexes_sf...\n")
    print(st_geometry_type(hexes_sf))
    print(st_crs(hexes_sf))
    
    # Run GRTS sampling
    samples <- grts(sframe = hexes_sf, n_base = n_base[1], DesignID = "CityGRTS")
    
    # Convert to plain data.frame (without sfc geometry class)
    samples_base <- as.data.frame(samples$sites_base)
    samples_df <- data.frame(lapply(samples_base, as.character), stringsAsFactors = FALSE)
    
    cat("R: grts sampling done\n")
}, error=function(e){
    cat("R error:", e$message, "\n")
    samples_df <- data.frame()
})

    """

    ro.r(r_code)

    # Retrieve sampled hexes

    with pandas2ri.converter.context():
        try:
            sampled_r = ro.globalenv["samples_df"]
            if sampled_r is None or len(sampled_r) == 0:
                sampled_df = pd.DataFrame()
            else:
                # If it's still an rpy2 object
                if isinstance(sampled_r, robj.vectors.DataFrame):
                    sampled_df = pd.DataFrame({col: list(sampled_r.rx2(col)) for col in sampled_r.names})
                else:
                    # Already a pandas DataFrame
                    sampled_df = pd.DataFrame(sampled_r)
            print("Sampled hexes:", len(sampled_df))
        except Exception as e:
            print("Error retrieving samples from R:", e)
            sampled_df = pd.DataFrame()


    # Merge sampled hexes
    # Ensure 'hex_id' is int for merging
    if not sampled_df.empty:
        sampled_df["hex_id"] = sampled_df["hex_id"].astype(int)
        sampled_hexes_utm = hexes_clipped_utm.merge(
            sampled_df,
            on="hex_id"
        )

        # Fix geometry after merge
        if "geometry_x" in sampled_hexes_utm.columns:
            sampled_hexes_utm = sampled_hexes_utm.set_geometry("geometry_x")
            sampled_hexes_utm = sampled_hexes_utm.rename_geometry("geometry")
        else:
            sampled_hexes_utm = sampled_hexes_utm.set_geometry("geometry")
        sampled_hexes_utm.crs = hexes_clipped_utm.crs
    else:
        sampled_hexes_utm = gpd.GeoDataFrame(columns=hexes_clipped_utm.columns)
        sampled_hexes_utm = sampled_hexes_utm.set_geometry("geometry")
        sampled_hexes_utm.crs = hexes_clipped_utm.crs


    # Determine which geometry column exists after merge
    geom_col = "geometry"
    if "geometry_x" in sampled_hexes_utm.columns:
        geom_col = "geometry_x"
    elif "geometry_y" in sampled_hexes_utm.columns:
        geom_col = "geometry_y"

    # Fetch OSM buildings using itertuples
    # Use itertuples to safely access hex_id and geometry
    buildings_list = []
    for row in sampled_hexes_utm.itertuples():
        hex_id = row.hex_id
        # Use geometry from hexes_clipped_utm by index
        hex_geom = hexes_clipped_utm.loc[hexes_clipped_utm.hex_id == hex_id, "geometry"].iloc[0]

        # Convert to WGS84
        hex_geom_wgs = gpd.GeoSeries([hex_geom], crs=utm_crs).to_crs(epsg=4326).iloc[0]

        print(f"Fetching OSM buildings for hex_id={hex_id}")
        try:
            result = ox.features_from_polygon(hex_geom_wgs, {"building": True})
            if not result.empty:
                result = gpd.GeoDataFrame(result).reset_index(drop=True)
                result["hex_id"] = hex_id
                if result.crs is None:
                    result.set_crs(epsg=4326, inplace=True)
                else:
                    result = result.to_crs(epsg=4326)
                buildings_list.append(result)
        except Exception as e:
            print(f"OSM fetch error for hex_id={hex_id}: {e}")
            continue

    # Combine all OSM buildings
    buildings = gpd.GeoDataFrame(
        pd.concat(buildings_list, ignore_index=True) if buildings_list else [],
        columns=["geometry","hex_id"],
        crs="EPSG:4326"
    )





    return {
        "city": city,
        "city_utm": city_utm,
        "hexes_utm": hexes_utm,
        "hexes_clipped_utm": hexes_clipped_utm,
        "sampled_hexes_utm": sampled_hexes_utm,
        "buildings": buildings
    }

# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    results = hex_tessellate_city_with_grts("Missoula", "Montana",
                                            hex_area_km2=1.0,
                                            sample_frac=0.2,
                                            random_seed=42)
    print("Hexes created:", len(results["hexes_clipped_utm"]))
    print("Sampled hexes:", len(results["sampled_hexes_utm"]))
    print("OSM buildings retrieved:", len(results["buildings"]))

    # --- Visualization ---
    import matplotlib.pyplot as plt
    import folium
    import os
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Static plot
    fig, ax = plt.subplots(figsize=(10, 10))
    results["hexes_clipped_utm"].plot(ax=ax, facecolor="none", edgecolor="lightgray", alpha=0.5)
    results["sampled_hexes_utm"].plot(ax=ax, facecolor="none", edgecolor="red", linewidth=2, alpha=0.8)
    if not results["buildings"].empty:
        results["buildings"].to_crs(results["hexes_clipped_utm"].crs).plot(ax=ax, color="blue", alpha=0.6, edgecolor="black", linewidth=0.5)
    plt.title("Missoula Hex Grid: Sampled Hexes & OSM Buildings")
    plt.axis("equal")
    plt.show()

    # Interactive map
    city_center = [results["city"].geometry.centroid.y.mean(),
                   results["city"].geometry.centroid.x.mean()]
    m = folium.Map(location=city_center, zoom_start=13, tiles="cartodbpositron")
    for geom in results["hexes_clipped_utm"].to_crs(epsg=4326).geometry:
        folium.GeoJson(geom, style_function=lambda x: {"color": "lightgray", "weight": 1, "fillOpacity": 0}).add_to(m)
    for geom in results["sampled_hexes_utm"].to_crs(epsg=4326).geometry:
        folium.GeoJson(geom, style_function=lambda x: {"color": "red", "weight": 2, "fillOpacity": 0.1}).add_to(m)
    for geom in results["buildings"].to_crs(epsg=4326).geometry:
        folium.GeoJson(geom, style_function=lambda x: {"color": "blue", "weight": 0.5, "fillOpacity": 0.5}).add_to(m)

    m.save(os.path.join(output_dir, "missoula_map.html"))
    print("Interactive map saved: missoula_map.html")
