"""
Understanding GIS: Assessment 1
@author [10906064]
Calculate the length of the World's Shortest Border, as fast as possible
"""

from time import time
from pyproj import Geod
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from geopandas import read_file, GeoSeries
from matplotlib.pyplot import subplots, savefig
from matplotlib_scalebar.scalebar import ScaleBar
from rtree import index

# set start time
start_time = time()	# NO CODE ABOVE HERE

'''ALL CODE MUST BE INSIDE HERE'''

# Constants (Best Practice)
CRS_WGS84 = 'EPSG:4326'
ELLIPSOID = 'WGS84'
EA_PROJ = '+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs'

def find_shortest_border(world):
      
    """
    Find the shortest border between neighboring countries using an R-tree for spatial indexing.
    
    Args:
        world (a GeoDataFrame): a dataframe containing geoseries' for countries geometry, ISO_A3 and name.
    
    Returns:
        a tuple: Shortest border length, country_a, country_b, and the border geometry as a shapley object.
    """
        
    # Initialise variables for tracking the shortest border, bordering countries and associated geometry
    shortest_border_length = float('inf')
    border_geom = None
    border_country_a = None
    border_country_b = None
    
    # Create an R-tree index for spatial queries
    rtree_idx = index.Index()
    
    # Insert each country's box into the R-tree with the index of the country in the world GeoDataFrame
    for idx, country_a in world.iterrows():
        rtree_idx.insert(idx, country_a.geometry.bounds)
    
    # Loop through each country in GeoDataFrame world
    for idx_a, country_a in world.iterrows():
        
        # Find neighboring countries based on bounding box intersections using R-tree spatial index
        possible_neighbors = list(rtree_idx.intersection(country_a.geometry.bounds))
        
        # Loop through each possible neigbouring country 
        for idx_b in possible_neighbors:
            
            # Avoid redundant and self comparisons
            if idx_a >= idx_b: 
                continue  
            
            #find GeoSeries from GeoDataFrame
            country_b = world.iloc[idx_b]
            
            # Check for and return a shared border between potentially neighbouring countries
            border = check_border(country_a, country_b)
            
            # Handle empty border returns, for robustness
            if border: 
                
                border_length = calculate_length(border, shortest_border_length)
                    
                # Update the knowledge of shortest border if border is valid, and smaller than current shortest
                if border_length and shortest_border_length > border_length:
                    
                    shortest_border_length = border_length
                    border_country_a = country_a
                    border_country_b = country_b
                    border_geom = border
                    
    return shortest_border_length, border_country_a, border_country_b, border_geom

def check_border(country_a, country_b) :
    
    """
    Check if two countries share a border, and return the geometry of the shared border.
    
    Args:
        country_a (GeoSeries: First country geometry.
        country_b (GeoSeries: Second country geometry.
    
    Returns:
        shapely geometry: LineString or MultiLineString of the shared border. 
        Avoids it being calculated twice.
    """
     
    geom_a = country_a.geometry
    geom_b = country_b.geometry
    
    # Calculate intersection (potential shared border)
    border = geom_a.intersection(geom_b)
    
    # Return the border if it is a LineString or MultiLineString
    # Also validates that border is not empty
    if border and border.geom_type in ['LineString', 'MultiLineString']:
        
        # Returns the border geometry, to avoid it being calculated twice 
        return border
    
    # Handle cases where one country is fully enclosed in another (which returns a polygon rather than LineString/MultiLineString)
    elif border:
        
        # Ensures working with LineString or MultiLineString
        return border.boundary
    
    return None

def calculate_length(border, max_length) :
    
    """
    Calculate the length of a given border (LineString or MultiLineString).
    
    Args:
        border (shapely geometry): Border geometry.
    
    Returns:
        float: Total length of the border.
    """
    
    cumulative_length = 0.0
    
    # Check if the geometry is a MultiLineString or LineString
    if border.geom_type == 'MultiLineString':
        
        # Iterate over each line segment in the MultiLineString
        for segment in border.geoms:
            distance = g.inv(segment.coords[0][0], segment.coords[0][1], 
                             segment.coords[1][0], segment.coords[1][1])[2]
            if (cumulative_length + distance) > max_length :
                return None
        
            else :
                cumulative_length += distance

    elif border.geom_type == 'LineString':
        
    # Handle simple LineString
        segment = border
        distance = g.inv(segment.coords[0][0], segment.coords[0][1], 
                     segment.coords[1][0], segment.coords[1][1])[2]
        cumulative_length += distance
        
        if (cumulative_length + distance) > max_length :
            return None
    
        else :
            cumulative_length += distance

    else : 
        
        # Handle cases where the geometry is a Point, or any other type (no length to calculate)
        return None

    return cumulative_length

def draw_map(border_length, country_a, country_b, border, graticule) :
    """
    Draw a map of the shared border between two countries and their surrounding areas,
    with a small inset context map in the bottom left corner.

    Args:
        border_length (float): Length of the shared border.
        country_a (GeoDataFrame): First country.
        country_b (GeoDataFrame): Second country.
        border (shapely geometry): Geometry of the shared border.
        graticule (GeoDataFrame): Graticule lines for the map.
    """

    # Create a GeoSeries for the border in WGS84
    border_series = GeoSeries([border], crs=CRS_WGS84)

    # Estimate a suitable UTM CRS for accurate plotting
    utm_crs = border_series.estimate_utm_crs()

    # Reproject all layers to UTM for plotting
    border_series = border_series.to_crs(utm_crs)
    country_a_series = GeoSeries([country_a.geometry], crs=CRS_WGS84).to_crs(utm_crs)
    country_b_series = GeoSeries([country_b.geometry], crs=CRS_WGS84).to_crs(utm_crs)
    graticule = graticule.to_crs(utm_crs)

    # Create main figure and axis
    my_fig, my_ax = subplots(1, 1, figsize=(16, 10))

    # Remove axis for main map
    my_ax.axis('off')

    # Set title with country names and border length
    my_ax.set_title(f"{country_a.loc['NAME']} and {country_b.loc['NAME']} share a border of {round(border_length, 2)}m")

    # Calculate bounds of the zoomed area, with buffer
    minx, miny, maxx, maxy = border_series.total_bounds
    buffer = 50  # 50 meter buffer around border
    my_ax.set_xlim([minx - buffer, maxx + buffer])
    my_ax.set_ylim([miny - buffer, maxy + buffer])

    # Plot countries and border with consistent styling
    country_a_series.plot(ax=my_ax, color='#ccebc5', edgecolor='#4daf4a', linewidth=0.5)
    country_b_series.plot(ax=my_ax, color='#fed9a6', edgecolor='#ff7f00', linewidth=0.5)
    border_series.plot(ax=my_ax, color='#984ea3', linewidth=2)
    graticule.plot(ax=my_ax, color='grey', linewidth=1)

    # Add legend
    my_ax.legend(handles=[
        Patch(facecolor='#ccebc5', edgecolor='#4daf4a', label=country_a.loc['NAME']),
        Patch(facecolor='#fed9a6', edgecolor='#ff7f00', label=country_b.loc['NAME']),
        Line2D([0], [0], color='#984ea3', lw=2, label='Border')
    ], loc='lower right')

    # Add north arrow
    x, y, arrow_length = 0.98, 0.99, 0.1
    my_ax.annotate('N', xy=(x, y), xytext=(x, y - arrow_length),
                   arrowprops=dict(facecolor='black', width=5, headwidth=15),
                   ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

    # Add scale bar
    my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.25))
    
    # Save the result
    try:
        savefig('out/2.png', bbox_inches='tight')

    # Handle issues saving map
    except Exception as e:
        print(f"Error saving the map: {e}")

if __name__ == "__main__":
    
    # Load country data geodatabase into a geodataframe 
    try:
       # Load country data geodatabase into a geodataframe
       world = read_file("./data/natural-earth/ne_10m_admin_0_countries.shp")
       graticule = read_file("./data/natural-earth/ne_110m_graticules_15.shp")
       
       # Catch normal errors, like an inccorect file path, or another error. Return status 1
    except FileNotFoundError as e :
        print(f"Error: {e}. Please check the file paths.")
        exit(1)
        
    except Exception as e:
        print(f"An unexpected error occurred while loading files: {e}")
        exit(1)

    # Set Geod for use in distance calculations
    g = Geod(ellps=ELLIPSOID)

    # Reproject all series in world to WGS84
    world = world.to_crs(CRS_WGS84)
    
    # Validation for robustness
    if graticule.empty:
        print("The 'graticule' GeoDataFrame is empty. Please provide a valid graticule shapefile.")
        exit(1)
    
    # Validation for robustness 
    if world.empty:
        raise ValueError("The GeoDataFrame 'world' is empty.")

    # Run algorithm and plot the map
    shortest_border, country_a, country_b, border = find_shortest_border(world)
    draw_map(shortest_border, country_a, country_b, border, graticule)

'''NO CODE BELOW HERE'''

# report runtime
print(f"completed in: {time() - start_time} seconds")	# NO CODE BELOW HERE

