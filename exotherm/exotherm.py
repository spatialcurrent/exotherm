# -*- coding: utf-8 -*-
#########################################################################
#
# Copyright (C) 2018 Spatial Current, Inc.
#
#########################################################################
import math
import re

import numpy as np

from pyextract.extract import extract

from enumerations import resolutions, webmercator_bbox, D2R, R2D


def to_wkt(value, value_type=None):

    if value is None:
        raise Exception("value missing")

    if not isinstance(value, dict):
        raise Exception("value is not a dict")

    if value_type is None:
        value_type = value.get("type")

    if value_type is None:
        raise Exception("value_type is missing and not given by value.")

    wkt = None

    value_type_lc = value_type.lower()

    if value_type_lc == "lonlat":

        coords = value.get("value") or value.get("coords")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        wkt = lonlat_to_wkt(coords)

    elif value_type_lc == "feature":

        wkt = feature_to_wkt(value)

    elif value_type_lc == "bbox":

        coords = value.get("value") or value.get("coords")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        wkt = bbox_to_wkt(*coords)

    else:
        raise Exception("Unknown value type "+value_type_lc+".")

    return wkt


def bbox_to_wkt(x0, y0, x1, y1):
    if None not in [x0, y0, x1, y1]:
        return 'POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s))' % (x0, y0, x0, y1, x1, y1, x1, y0, x0, y0)
        # wkt = 'SRID=%s;POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))' % (srid, x0, y0, x0, y1, x1, y1, x1, y0, x0, y0)
    else:
        return None

def lonlat_to_wkt(lonlat):
    return "POINT ("+str(lonlat[0])+" "+str(lonlat[1])+")"


def polygon_to_wkt(coords):
    return "POLYGON ("+(", ".join([("("+(", ".join([str(p[0])+" "+str(p[1]) for p in ring]))+")") for ring in coords]))+")"


def feature_to_wkt(feature):

    if feature is None:
        raise Exception("missing feature")

    geom = feature.get("geometry")

    if geom is None:
        raise Exception("missing geometry")

    return geom_to_wkt(geom)

def geom_to_wkt(geom):

    wkt = None

    geom_type_lc = (geom.get("type") or "").lower()
    if geom_type_lc == "point":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        wkt = lonlat_to_wkt(coords)

    elif geom_type_lc == "polygon":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        wkt = polygon_to_wkt(coords)

    elif geom_type_lc == "geometrycollection":

        geometries = geom.get("geometries")

        if geometries is None or len(geometries) == 0:
            raise Exception("missing geometries")

        wkt = "GEOMETRYCOLLECTION ("+(",".join([geom_to_wkt(g) for g in geometries]))+")"

    else:
        raise Exception("unknown geometry type "+geom_type_lc+".")

    return wkt


def wkt_to_geom(wkt):

    if wkt is None or len(wkt) == 0:
        raise Exception("missing wkt")

    m = re.match("(?P<type>\\w+)\\s*\\((?P<value>.+)\\)", wkt)
    geometry_type = m.group("type")
    value = m.group("value")

    if geometry_type is None or len(geometry_type) == 0:
        raise Exception("missing geometry_type")

    if value is None or len(value) == 0:
        raise Exception("missing value")

    geom_type_lc = (m.group("type") or "").lower()

    geom = None
    if geom_type_lc == "point":

        coords = [float(x.strip()) for x in value.split(" ")]

        if len(coords) != 2:
            raise Exception("invalid number of coordinates")

        geom = {
            "type": "Point",
            "coordinates": coords
        }

    elif geom_type_lc == "polygon":

        coordinates = [[[float(z.strip()) for z in y.strip().split(" ")] for y in x.group("ring").split(",")] for x in re.finditer(",?(\\s*)\\((?P<ring>[^)]+)\\)(\\s*)", value)]

        if len(coords) == 0:
            raise Exception("invalid number of coordinates")

        geom = {
            "type": "Polygon",
            "coordinates": coords
        }

    else:
        raise Exception("invalid geometry type")

    return geom

def point_to_bbox(xy, radius):
    return [
        xy[0] - radius,
        xy[1] - radius,
        xy[0] + radius,
        xy[1] + radius
    ]

def llbbox_to_mercator(llbbox):
    minlonlat = forward_mercator([llbbox[0], llbbox[1]])
    maxlonlat = forward_mercator([llbbox[2], llbbox[3]])
    return [minlonlat[0], minlonlat[1], maxlonlat[0], maxlonlat[1]]


def mercator_to_llbbox(bbox):
    minlonlat = inverse_mercator([bbox[0], bbox[1]])
    maxlonlat = inverse_mercator([bbox[2], bbox[3]])
    return [minlonlat[0], minlonlat[1], maxlonlat[0], maxlonlat[1]]


def inverse_mercator(xy):
    lon = (xy[0] / 20037508.34) * 180
    lat = (xy[1] / 20037508.34) * 180
    lat = 180 / math.pi * (2 * math.atan(math.exp(lat * math.pi / 180)) - math.pi / 2)
    return [lon, lat]


def forward_mercator(lonlat):
    x = lonlat[0] * 20037508.34 / 180
    try:
        n = math.tan((90 + lonlat[1]) * math.pi / 360)
    except ValueError:
        n = 0
    if n <= 0:
        y = float("-inf")
    else:
        y = math.log(n) / math.pi * 20037508.34
    return (x, y)


def bbox_intersects(a,b):
    return ( a[0] < b[2] and a[2] > b[0] ) and ( a[1] < b[3] and a[3] > b[1] )


def getMaxX(res, size, bbox):
    maxX = int(
        round(
            (bbox[2] - bbox[0]) /
            (res * size)
        )
    ) - 1
    return maxX


def getMaxY(res, size, bbox):
    maxY = int(
        round(
            (bbox[3] - bbox[1]) /
            (res * size)
        )
    ) - 1
    return maxY


def flip_y(y, z, size=256, bbox=[-20037508.34,-20037508.34,20037508.34,20037508.34]):
    res = resolutions[int(z)]
    maxY = int(
        round(
            (bbox[3] - bbox[1]) /
            (res * size)
        )
    ) - 1
    return maxY - y


def tms_to_bbox(x,y,z):
    e = tile_to_lon(x+1,z)
    w = tile_to_lon(x,z)
    s = tile_to_lat(y+1,z)
    n = tile_to_lat(y,z)
    return [w, s, e, n]


def tile_to_lon(x, z):
    return (x/math.pow(2,z)*360-180);


def tile_to_lat(y, z):
    n = math.pi - 2 * math.pi * y / math.pow(2,z);
    return ( R2D * math.atan(0.5*(math.exp(n)-math.exp(-n))));


def lon_to_tile(lon, z):
    return (180 + lon) *(math.pow(2,z)/360)


def lat_to_tile(lat, z):
    lat_rad = math.radians(lat)
    return (1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * math.pow(2,z)


def tms_to_quadkey(x,y,z):
    quadKey = []
    for i in range(z,0,-1):
        digit = 0
        mask = 1 << ( i - 1)
        if ((x & mask) != 0):
            digit += 1
        if ((y & mask) != 0):
            digit += 1
            digit += 1
        quadKey.append(str(digit));

    return ''.join(quadKey);


def quadkey_to_tms(u):
    iz = len(u)
    ix = quadkey_to_x(u)
    iy = quadkey_to_y(u)
    return iz, ix, iy


def quadkey_to_x(u):
    x = 0
    for i in range(0,len(u)):
        x = x * 2
        if ( int(u[i]) == 1 ) or ( int(u[i]) == 3 ):
            x += 1
    return x


def quadkey_to_y(u):
    y = 0
    for i in range(0,len(u)):
        y = y * 2
        if ( int(u[i]) == 2 ) or ( int(u[i]) == 3 ):
            y += 1
    return y



def to_bbox(value, buffer_distance=0, value_type=None):

    if value is None:
        raise Exception("value missing")

    if not isinstance(value, dict):
        raise Exception("value is not a dict")

    if value_type is None:
        value_type = value.get("type")

    if value_type is None:
        raise Exception("value_type is missing and not given by value.")

    bbox = None

    value_type_lc = value_type.lower()

    if value_type_lc == "lonlat":

        coords = value.get("value") or value.get("coords")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        bbox = mercator_to_llbbox(point_to_bbox(forward_mercator(coords), buffer_distance))

    elif value_type_lc == "feature":

        bbox = feature_to_bbox(value, buffer_distance=buffer_distance)

    else:
        raise Exception("Unknown value type "+value_type_lc+".")

    return bbox

def feature_to_bbox(feature, buffer_distance=None):

    if feature is None:
        raise Exception("missing feature")

    geom = feature.get("geometry")

    if geom is None:
        raise Exception("missing geometry")

    return geom_to_bbox(geom, buffer_distance=buffer_distance)

def geom_to_bbox(geom, buffer_distance=None):

    geom_type_lc = (geom.get("type") or "").lower()
    if geom_type_lc == "point":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        bbox = mercator_to_llbbox(point_to_bbox(forward_mercator(coords), buffer_distance))

    elif geom_type_lc == "polygon":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        points = [forward_mercator(p) for p in coords[0]]
        minx = points[0][0] - buffer_distance
        miny = points[0][1] - buffer_distance
        maxx = points[0][0] + buffer_distance
        maxy = points[0][1] + buffer_distance
        for p in points[1:]:
            if p[0] - buffer_distance < minx:
                minx = p[0] - buffer_distance
            if p[0] + buffer_distance > maxx:
                maxx = p[0] + buffer_distance
            if p[1] - buffer_distance < miny:
                miny = p[1] - buffer_distance
            if p[1] + buffer_distance > maxy:
                maxy = p[1] + buffer_distance
        bbox = mercator_to_llbbox([minx, miny, maxx, maxy])

    elif geom_type_lc == "geometrycollection":

        geometries = geom.get("geometries")

        if geometries is None or len(geometries) == 0:
            raise Exception("missing geometries")

        bounding_boxes = [geom_to_bbox(g, buffer_distance=buffer_distance) for g in geometries]
        minx = bounding_boxes[0][0]
        miny = bounding_boxes[0][1]
        maxx = bounding_boxes[0][2]
        maxy = bounding_boxes[0][3]

        for bb in bounding_boxes[1:]:
            if bb[0] < minx:
                minx = bb[0]
            if bb[1] < miny:
                miny = bb[1]
            if bb[2] > maxx:
                maxx = bb[2]
            if bb[3] > maxy:
                maxy = bb[3]

        bbox = [minx, miny, maxx, maxy]
    else:
        raise Exception("unknown geometry type")

    return bbox



def zxy_to_lonlat(zxy):
    z, x, y = [int(n) for n in zxy.split("/")]
    return (tile_to_lon(x, z), tile_to_lat(y, z))


def to_centroid(value, value_type=None):

    if value is None:
        raise Exception("value missing")

    if not isinstance(value, dict):
        raise Exception("value is not a dict")

    if value_type is None:
        value_type = value.get("type")

    if value_type is None:
        raise Exception("value_type is missing and not given by value.")

    centroid = None

    value_type_lc = value_type.lower()

    if value_type_lc == "lonlat":

        coords = value.get("value") or value.get("coords")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        centroid = coords

    elif value_type_lc == "feature":

        centroid = feature_to_centroid(value)

    elif value_type_lc == "bbox":

        coords = value.get("value") or value.get("coords")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 4:
            raise Exception("invalid number of coordinates")

        centroid = [
           ((coords[0]+coords[2])/2.0),
           ((coords[1]+coords[3])/2.0)
        ]

    else:
        raise Exception("Unknown value type "+value_type_lc+".")

    return centroid


def feature_to_centroid(feature):

    if feature is None:
        raise Exception("missing feature")

    geom = feature.get("geometry")

    if geom is None:
        raise Exception("missing geometry")

    return geom_to_centroid(geom)

def geom_to_centroid(geom):

    centroid = None

    geom_type_lc = (geom.get("type") or "").lower()
    if geom_type_lc == "point":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        if len(coords) != 2:
            raise Exception("Invalid number of coordinates")

        centroid = coords

    elif geom_type_lc == "polygon":

        coords = geom.get("coordinates")

        if coords is None or len(coords) == 0:
            raise Exception("missing coordinates")

        centroid = centroid_of_rings(coords)

    elif geom_type_lc == "geometrycollection":

        geometries = geom.get("geometries")

        if geometries is None or len(geometries) == 0:
            raise Exception("missing geometries")

        centroid = centroid_of_points([geom_to_centroid(g) for g in geometries])

    else:
        raise Exception("unknown geometry type "+geom_type_lc+".")

    return centroid

def centroid_of_rings(rings):
    return centroid_of_points([centroid_of_points(r) for r in rings])

def centroid_of_points(points):

    lon = np.mean(np.array([p[0] for p in points][0:2], dtype=np.float64))
    lat = np.mean(np.array([p[1] for p in points][0:2], dtype=np.float64))

    return [lon, lat]
