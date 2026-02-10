"""
Understanding GIS: Assessment 2
@author 10906064

An Implementation Weighted Redistribution Algorithm (Huck et al. 2015)
"""

''' 
Loosely styled on PEP8 :
    
    "Know when to be inconsistent -- sometimes the style guide just doesn't apply.
    When in doubt, use your best judgment." - PEP8

'''

# Library Imports
from math import pi, sqrt
from time import time
from sys import exit
import numpy as np
from numpy import column_stack, zeros
from numpy.random import uniform, randint, default_rng
from geopandas import read_file
from rasterio import open as rio_open
from rasterio.features import geometry_mask
from rasterio.plot import show as rio_show
from shapely.geometry import Point
from skimage.draw import disk
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# set start time
start_time = time()  # NO CODE ABOVE HERE

'''ALL CODE MUST BE INSIDE HERE'''

# =============================================================================
#                                   Constants
# =============================================================================

# Spatial Ambiguity Parameter
s = 0.3

# Number of potential seed points to generate for each administration area (W)
sample_number = 25

# set British National Grid projection for validation
crs = 'EPSG:27700'

'''Seed and rng object for reproducible results : Used in Figure 3 generation'''
seed = randint(0, 5000)
rng = default_rng(seed)

# =============================================================================
#                               Helper Functions
# =============================================================================
    
def get_weight(seed, weighting_surface, weight_data):

    """
    Retrieve the value of the weighting surface at the location of a given seed

    Args:
        seed (tuple): Geographic coordinates (x, y).
        weighting_surface (rasterio DatasetReader): Weighting raster.
        weight_data (numpy.ndarray): Raster data as numpy array.

    Returns:
        float: Weight of the seed. None if out of bounds.
    """

    try:
        
        # Convert geographic coordinates to row and column image space
        row, col = weighting_surface.index(seed[0], seed[1])

        # Obtain value from seed location
        value = weight_data[row, col]

        return value

    except IndexError:

        # Handle unexpected conversion out of bounds errors 
        print("Potential seed is out of bounds")

        # Move on to next seed rather than exiting the program
        return None

def generate_seed(geom, weighting_surface, weight_data):

    """
    Generate random seeds within an admin area. Find the weight for each corresponding 
    point in the weighting layer. Return the seed with the greatest corresponding weight

    Args:
        geom (shapely.geometry): Geometry defining bounds.
        sample_number (int): Number of seeds to generate.

    Returns:
        seed (tuple) : Geographic coordinates (x, y)
    """

    # Create empty list to store seed points within the administrative area
    viable_seeds = []

    # Get bounding box co-ordinates
    minx, miny, maxx, maxy = geom.bounds
    
    # Generate 'sample_number' of seeds within the administrative area
    while len(viable_seeds) < sample_number:
        
        # Generate array of pseudorandom x, y co-ordinates within the admin geoms
        # bounding box, the length of the remaining number of seeds required
        x_coords = rng.uniform(low=minx, high=maxx, size=(sample_number - len(viable_seeds)))
        y_coords = rng.uniform(low=miny, high=maxy, size=(sample_number - len(viable_seeds)))

        # List comprehend the zip of these arrays into a list of co-ordinates
        potential_seeds = [(x, y) for x, y in zip(x_coords, y_coords)]
        
        # Filter the seeds outside the admin area
        # Add the list of potential seeds remaining, to the list of viable seeds 
        viable_seeds.extend([seed for seed in potential_seeds if geom.contains(Point(seed))])

    # Calculate the weight of each viable seed. 
    # Add each weight and corresponding seed co-ordinate as a tuple, to a list 
    weights = [(seed, get_weight(seed, weighting_surface, weight_data)) 
               for seed in viable_seeds]

    # Find the maximum weight in the list, and return the corresponding seed
    # More efficient than iterative conditional comparison of the list
    return max(weights, key=lambda item: item[1])[0]

def calculate_distribution(seed, output_surface, r, weighting_surface):

    """
   Distribute values around a seed point using a Euclidean decay function.
   Applies these values to an output surface.

   Args:
       seed (tuple): Geographic coordinates of the seed (x, y).
       output_surface (numpy.ndarray): Raster array to update.
       radius (float): Radius in geographic distance.
       weighting_surface (rasterio.DatasetReader): Raster for coordinate
       transformations.

   Returns:
       numpy.ndarray: Updated output surface.
   """

    # Calculate the seed's geographic co-ordinate in image space for continuity
    seed = weighting_surface.index(seed[0], seed[1])
    
    # Euclidean distance function, declared for reusability 
    euclidean_distance = lambda row, col: sqrt((row - seed[0])**2+ (col - seed[1])**2)

    # Iterate through each cell in a circle of radius 'R' around centerpoint 'seed'
    # Column stack of disk function returns an iterable 2D-array of only relevant cells
    for row, col in column_stack(disk(seed, r)):

        # Out of bounds error handling for when a cell is out of bounds of the 
        # output surface
        try:
            output_surface[row, col] += 1 - ((euclidean_distance(row, col) / r))

        # Handle error and continue to next cell
        except IndexError:
            continue

    # Return output surface with updated distribuution radius around seed
    return output_surface

def filter_points(geom, spatial_index, point_data):
    
    """
    Filters points from a GeoDataFrame that are located within the admin geometry.

    Args:
        geom (shapely.geometry): The geometry of an admin area
        spatial_index (geopandas.sindex): A spatial index createn in core algorithm
        point_data (GeoDataFrame): A GeoDataFrame containing point geometries 

    Returns:
        GeoDataFrame: A subset of `point_data` located within admin area 'geom'
    """
    # Use spatial index to filter points within the bounding box
    potential_points = list(spatial_index.intersection(geom.bounds))

    # Retrieve point data from the GeoDataFrame
    potential_points = point_data.iloc[potential_points]

    # Find points that are  within the geometry of administration area
    points = potential_points[potential_points.geometry.within(geom)]

    return points
    
# =============================================================================
#                               Core Function
# =============================================================================

def weighted_redistribution(admin_areas, point_data, weighting_surface,
                            weight_data):
    """
    Redistribute point data within administrative areas based 
    on a weighting surface (e.g., population density). For each administrative 
    area, the function generates a distribution radius, identifies tweet locatin 
    points within the area, generates a new location point, then applies a euclidian 
    decay of distribution of weight radiating up to a certain distance of ambiguity
    away.

    Args:
        admin_areas (GeoDataFrame): Geospatial DataFrame containing 
            administrative boundaries.
        point_data (GeoDataFrame): Geospatial DataFrame containing point data 
            (e.g., spatially ambiguous points to redistribute).
        weighting_surface (rasterio.DatasetReader): Raster dataset representing 
            the weighting surface (e.g., population density).
        weight_data (numpy.ndarray): Numpy array of raster data representing 
            the weighting surface.

    Returns:
        numpy.ndarray: A 2D array (same shape as `weight_data`) containing the 
            redistributed output values.
    """
    
    # Create a blank output surface
    output_surface = zeros(weight_data.shape)

    # Create a spatial index of point data
    spatial_index = point_data.sindex

    # Iterate through each administrative area in level 3
    for idx, admin in admin_areas.iterrows():

        # Avoid overhead from repeatedly accessing the geometry from the GeoDataFrame
        geom = admin.geometry
 
        # Validate geometry
        if not geom.is_valid:
            print(f"Skipping invalid geometry at index {idx}")
            continue

        # Obtain points in administrative area
        points = filter_points(geom, spatial_index, point_data)

        # Skip admin areas with no points, avoid redundant computation
        if points.empty:
            continue

        # Calculate the distribution radius for administration area
        distribution_radius = sqrt(geom.area * s / pi)

        # Convert length from geographic to image-space by dividing by pixel resolution
        r = round(distribution_radius / abs(weighting_surface.transform.a))

        # Iterate through each point within the administrative area
        for _, point in points.iterrows():

            # Generate seed and add distribution to output surface
            seed = generate_seed(geom, weighting_surface, weight_data)
            output_surface = calculate_distribution(seed, output_surface, r, weighting_surface)

    return output_surface

# =============================================================================
#                               Visualise Results
# =============================================================================

def visualise_results(admin_areas, point_data, raster, output_surface, ax=None, title=None, show=True):
    """
    Visualise the redistributed tweet data and weighted surface.

    Args:
        admin_areas (GeoDataFrame): Administrative boundaries.
        point_data (GeoDataFrame): Original tweet point data.
        raster (rasterio.DatasetReader): Raster used for redistribution.
        output_surface (np.ndarray): The 2D redistribution output array.
        ax (matplotlib.axes.Axes, optional): Axes object for subplot. If None, creates new figure.
        title (str, optional): Title for the plot. Defaults to a static string if not provided.
        show (bool): Whether to display the figure. Set False when plotting subplots.
    """

    if raster.crs != admin_areas.crs:
        raise ValueError("CRS mismatch: Ensure raster and admin areas have the same CRS.")

    # Create new figure and axis if not provided
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 16))

    # Remove axes
    ax.axis('off')

    # Add title if given, otherwise use default
    plot_title = title or "Re-distrubted Spatially Ambiguous Tweet Data (Huck et al., 2015)"
    ax.set_title(plot_title, fontsize=20, fontweight="bold", ha="center")

    # Plot point data
    point_data.plot(
        ax=ax,
        color="red",
        markersize=10,
        alpha=0.8,
        zorder=3
    )

    # Combine all geometries into one
    admin_union = admin_areas.geometry.unary_union

    # Create a mask for pixels outside the GM polygon
    mask = geometry_mask(
        [admin_union],
        transform=raster.transform,
        invert=True,
        out_shape=output_surface.shape
    )

    # Apply the mask: Set outside pixels to NaN
    clipped_raster = np.where(mask, output_surface, np.nan)

    # Use the cividis colormap suitable for colourblind people
    cividis_cmap = plt.cm.cividis
    cividis_colors = cividis_cmap(np.linspace(0, 1, 256))
    cividis_colors[0, -1] = 0  # Set alpha of the first color to transparent
    transparent_cividis = ListedColormap(cividis_colors)

    # Plot the administrative boundaries as the base layer
    admin_areas.plot(
        ax=ax,
        edgecolor='black',
        facecolor='lightblue',
        alpha=1,
        zorder=1
    )

    # Plot the weighted redistribution output surface
    rio_show(
        clipped_raster,
        ax=ax,
        transform=raster.transform,
        cmap=transparent_cividis,
        alpha=0.9,
        zorder=2
    )

    # Only add north arrow, scalebar, and legend if this is a standalone figure
    if fig is not None:

        # Add north arrow
        x, y, arrow_length = 0.98, 0.95, 0.05
        ax.annotate(
            "N",
            xy=(x, y),
            xytext=(x, y - arrow_length),
            arrowprops=dict(facecolor="black", width=5, headwidth=15),
            ha="center",
            va="center",
            fontsize=20,
            xycoords=ax.transAxes
        )

        # Add scalebar
        ax.add_artist(ScaleBar(
            dx=1,
            units="m",
            location="lower left",
            length_fraction=0.35
        ))

        # Add colorbar and label
        cbar = fig.colorbar(
            ScalarMappable(cmap=transparent_cividis, norm=plt.Normalize(0, 1)),
            ax=ax,
            shrink=0.5,
            aspect=20,
            pad=0.02
        )  
        cbar.set_label("Likelihood of Tweet Location")

        # Create gradient patch for legend
        gradient_patch = Patch(
            facecolor=plt.cm.cividis(0.7),
            edgecolor='none',
            label='Redistribution Likelihood'
        )

        # Create handles for legend
        handles = [
            Line2D([0], [0], marker='o', color='w', label='Tweet Locations',
                   markerfacecolor='red', markersize=10),
            gradient_patch
        ]

        ax.legend(
            handles=handles,
            loc='lower right',
            title="Legend",
            bbox_to_anchor=(1.1, 0.0)
        )

    # Show figure if not inside a subplot
    if fig is not None and show:
        plt.show()

def visualise_multiple_distributions(admin_areas, point_data, raster, s_values, w_values):
    """
    Generate a 12-panel subplot visualisation of tweet redistributions
    across combinations of spatial ambiguity (s) and seed influence (w).

    Args:
        admin_areas (GeoDataFrame): Administrative boundaries.
        point_data (GeoDataFrame): Original tweet point data.
        raster (rasterio.DatasetReader): Population weighting raster.
        s_values (list): List of spatial ambiguity values.
        w_values (list): List of weighting values.
    """
    fig, axes = plt.subplots(len(s_values), len(w_values), figsize=(18, 16))

    for i, s_val in enumerate(s_values):
        for j, w_val in enumerate(w_values):
            global s, sample_number
            s = s_val
            sample_number = w_val
            weight_data = raster.read(1)
            output_surface = weighted_redistribution(admin_areas, point_data, raster, weight_data)
            subplot_title = f"s = {s_val}, w = {w_val}"
            visualise_results(admin_areas, point_data, raster, output_surface,
                              ax=axes[i, j], title=subplot_title, show=False)
            axes[i, j].text(0.02, 0.95, chr(65 + i * len(w_values) + j) + '.',
                            transform=axes[i, j].transAxes,
                            fontsize=12, verticalalignment='top', weight='bold')

    # Axis labels
    for i, s_val in enumerate(s_values):
        axes[i, 0].text(-0.08, 0.5, f"{s_val}", transform=axes[i, 0].transAxes,
                       fontsize=14, va='center', ha='right', weight='bold')

    for j, w_val in enumerate(w_values):
        axes[-1, j].text(0.5, -0.15, f"{w_val}", transform=axes[-1, j].transAxes,
                        fontsize=14, ha='center', weight='bold')

    fig.text(0.06, 0.5, 'Spatial Ambiguity (s):', va='center', rotation='vertical', fontsize=16, weight='bold')
    fig.text(0.5, 0.06, 'Influence of the weighting layer:', ha='center', fontsize=16, weight='bold')

    plt.tight_layout(rect=[0.08, 0.08, 0.95, 0.98])
    plt.show()

if __name__ == "__main__":

    # Load and validate data to work with
    try:

        # Variables to be passed as arguments, not accessed as global variables
        admin_areas = read_file("./data/gm-districts.shp")
        
        if admin_areas.empty:
            raise ValueError("Admin areas dataset is empty. Please check the input file.")

        point_data = read_file("./data/level3-tweets-subset.shp")

        if point_data.empty:
            raise ValueError("Point data dataset is empty. Please check the input file.")

        weighting_surface = "./data/100m_pop_2019.tif"

    # Catch nornal errors : File path issues, or unexpected errors
    except FileNotFoundError as e:
        print(f"Error: {e}. Please check the file paths.")
        exit(1)

    except Exception as e:
        print(f"An unexpected error occurred while loading files: {e}")
        exit(1)
   
    # Reproject data to match raster CRS
    with rio_open(weighting_surface) as raster:
        
        # Confirm that weighting surface is in correct CRS (BNG in this case)
        if raster.crs != crs :
            raster.to_crs(crs)
            
        # Ensure all working data is in the same CRS
        admin_areas = admin_areas.to_crs(raster.crs)
        point_data = point_data.to_crs(raster.crs)

        # Obtain the first band of population density values
        weight_data = raster.read(1)

        # Run weighted redistribution algorithm
        output_surface = weighted_redistribution(admin_areas, point_data,raster, weight_data)

        # Visualise the 12 panel comparative figure
        s_values = [0.05, 0.1, 0.2, 0.5]
        w_values = [10, 25, 50]
        visualise_multiple_distributions(admin_areas, point_data, raster, s_values, w_values)
        
        # Visualise the single figure
        #visualise_results(admin_areas, point_data, raster, output_surface)
        
# report runtime
print(f"completed in: {time() - start_time} seconds")  # NO CODE BELOW HERE