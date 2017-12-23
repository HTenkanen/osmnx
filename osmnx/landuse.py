################################################################################
# Module: landuse.py
# Description: Download and plot landuse features from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

from shapely.geometry import Point, LineString, Polygon, box
import geopandas as gpd
from .core import overpass_request, bbox_from_point, gdf_from_place
from .utils import log, geocode, bbox_to_poly

# List of available amenities
available_landuse = []

def parse_landuse_query(west, south, east, north, landuse_classes=None, timeout=180, maxsize=''):
    """
    Parse the Overpass QL query based on the list of amenities.

    Parameters
    ----------
    west : float
        Westernmost coordinate of the bounding box of the search area.
    south : float
        Southernmost coordinate from bounding box of the search area.
    east : float
        Easternmost coordinate from bounding box of the search area.
    north : float
        Northernmost coordinate from bounding box of the search area.
    landuse_classes : list
        List of landuse classes that will be used for finding the landuse features from the selected area.
    timeout : int
        Timeout for the API request. 
    """
    if landuse_classes:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["landuse"~"{landuse_classes}"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(way["landuse"~"{landuse_classes}"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["landuse"~"{landuse_classes}"]'
                          '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')

        # Parse amenties
        query_str = query_template.format(landuse_classes="|".join(landuse_classes), north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)
    else:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["landuse"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(way["landuse"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["landuse"]'
                          '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')

        # Parse amenties
        query_str = query_template.format(north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)

    return query_str


def validate_landuse(landuse_classes=None):
    """
    Validate the list of amenities. Warn users about unrecognized amenities. 
    """
    # Count valid amenities
    valid_cnt = 0
    for amenity in landuse_classes:
        # Check known amenities
        if amenity in available_landuse:
            valid_cnt += 1
        else:
            Warning("Amenity {} was not recognized.".format(amenity))
    if valid_cnt == 0:
        Warning(
            "There were not any recognized amenities. You might not get any results. Check available amenities by running: ox.pois.available")


def osm_landuse_download(polygon=None, landuse_classes=None, north=None, south=None, east=None, west=None,
                     timeout=180, max_query_area_size=50 * 1000 * 50 * 1000):
    """
    Get landuse features from OpenStreetMap based on selected landuse types.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Polygon that will be used to limit the POI search. 
    landuse_classes : list
        List of landuse classes that will be used for finding the landuse features from the selected area.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Landuse features and the tags associated to them as Geopandas GeoDataFrame.
    """

    # Validate amenities
    if landuse_classes:
        validate_landuse(landuse_classes)

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = parse_landuse_query(landuse_classes=landuse_classes, west=west, south=south, east=east, north=north)

    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        # polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = parse_landuse_query(landuse_classes=landuse_classes, west=west, south=south, east=east, north=north)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the POIs
    responses = overpass_request(data={'data': query}, timeout=timeout)

    return responses


def create_landuse_gdf(polygon=None, landuse_classes=None, north=None, south=None, east=None, west=None, retain_invalid=False):
    """
    Parse GeoDataFrames from POI json that was returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    landuse_classes: list
        List of landuse classes that will be used for finding the landuse features from the selected area.
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    Geopandas GeoDataFrame with POIs and the associated attributes. 
    """
    responses = osm_landuse_download(polygon=polygon, landuse_classes=landuse_classes, north=north, south=south, east=east, west=west)

    # Parse Landuse features
    landuse_nodes = {}
    landuse_ways = {}

    for result in responses['elements']:
        if result['type'] == 'node' and 'tags' in result:
            try:
                point = Point(result['lon'], result['lat'])

                node = {
                    'osmid': result['id'],
                    'geometry': point
                }
                if 'tags' in result:
                    for tag in result['tags']:
                        node[tag] = result['tags'][tag]
                landuse_nodes[result['id']] = node

            except Exception:
                log('Point has invalid geometry: {}'.format(result['id']))

        elif result['type'] == 'relation':
            # TODO: Add functionalities to parse 'relation' tags.
            pass
        elif result['type'] == 'way':
            vertices = {}
            # Get the vertices
            for result in responses['elements']:
                if 'type' in result and result['type'] == 'node':
                    vertices[result['id']] = {'lat': result['lat'],
                                              'lon': result['lon']}
            # Parse the Polygons
            for result in responses['elements']:
                if result['type'] == 'way' and 'tags' in result:
                    nodes = result['nodes']
                    try:
                        polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])
                    except Exception:
                        log('Polygon has invalid geometry: {}'.format(nodes))
                    landuse_feat = {'nodes': nodes,
                                'geometry': polygon}

                    if 'tags' in result:
                        for tag in result['tags']:
                            landuse_feat[tag] = result['tags'][tag]

                    landuse_ways[result['id']] = landuse_feat

    gdf = gpd.GeoDataFrame(landuse_ways).T
    gdf.crs = {'init': 'epsg:4326'}
    return gdf


def landuse_from_point(point, distance=None, landuse_classes=None, retain_invalid=False):
    """
    Get landuse features within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    landuse_classes : list
        List of landuse classes that will be used for finding the landuse features from the selected area.
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame 
    """
    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_landuse_gdf(landuse_classes=landuse_classes, north=north, south=south, east=east, west=west)


def landuse_from_address(address, distance, landuse_classes=None, retain_invalid=False):
    """
    Get OSM Points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    landuse_classes : list
        List of landuse classes that will be used for finding the landuse features from the selected area.
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return landuse_from_point(point=point, landuse_classes=landuse_classes, distance=distance)


def landuse_from_polygon(polygon, landuse_classes=None, retain_invalid=False):
    """
    Get landuse features within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the POIs are search from. 
    landuse_classes : list
        List of landuse classes that will be used for finding the landuse features from the selected area.
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    return create_landuse_gdf(polygon=polygon, landuse_classes=landuse_classes)


def landuse_from_place(place, landuse_classes=None, retain_invalid=False):
    """
    Get landuse features within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    landuse_classes : list
        List of landuse classes that will be used for finding the POIs from the selected area.
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_landuse_gdf(polygon=polygon, landuse_classes=landuse_classes, retain_invalid=retain_invalid)
