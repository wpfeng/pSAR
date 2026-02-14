#!/usr/bin/env python
#
#
import numpy as np
from pyproj import CRS, Transformer

class CustomUTMConverter:
    """
    This was firstly given by Doubao based on my request to convert lon/lat to local coordinates in km.
    Then I slightly modify the codes to allow a reference site given. Below are two new functions to 
    return easting and northing, and back to lon/lat...
    def ll2local(self, lats,lons)
    def local2ll(self, easting, northing)
    
    Provided by Wanpeng Feng, @SYSU, Guangzhou, 2025/06/21
    #

    Custom UTM projection converter that allows specifying the central meridian directly
    instead of relying on standard UTM zones. This provides more flexibility for
    regional projections with minimal distortion.
    """

    def __init__(self, central_longitude=0, latitude_of_origin=0, false_easting=500000,
                 false_northing=0, scale_factor=0.9996, datum='WGS84',origin=None):
        """
        Initialize the custom UTM projection converter.
        
        origin, [lat,lon]

        Args:
            central_longitude: Central meridian for the projection (degrees)
            latitude_of_origin: Origin latitude (degrees), default 0
            false_easting: Easting value assigned to the central meridian (meters), default 500000
            false_northing: Northing value assigned to the origin latitude (meters),
                            default 0 for Northern Hemisphere, 10000000 for Southern Hemisphere
            scale_factor: Scale factor applied to the central meridian, default 0.9996 (UTM standard)
            datum: Geodetic datum, default 'WGS84'
        """
        self.central_longitude  = central_longitude
        self.latitude_of_origin = latitude_of_origin
        self.false_easting      = false_easting
        self.false_northing     = false_northing
        self.scale_factor       = scale_factor
        self.datum              = datum
        self.origin             = origin
        self.utm_zone           = None
        self.utm_letter         = None
        #
        if origin is not None:
            self.latitude_of_origin = origin[0]
            self.central_longitude  = origin[1]
        #
        # Create the custom UTM coordinate reference system
        self.crs = self._create_custom_crs()

        # Create transformers for coordinate conversions
        self.geodetic_crs = CRS("EPSG:4326")  # WGS84 geographic coordinate system
        self.to_utm_transformer = Transformer.from_crs(self.geodetic_crs, self.crs)
        self.from_utm_transformer = Transformer.from_crs(self.crs, self.geodetic_crs)
        #
        if self.origin is not None:
            e0,n0 = self.to_utm(self.origin[0],self.origin[1])
        else:
            e0,n0 = 0,0
        #
        self.origin_utm = [e0,n0]
    #
    def get_utm_zone(self):
        #
        if self.origin is not None:
            #
            lon = self.origin[1]
            lat = self.origin[0]
            #
            lon = (lon + 180) % 360 - 180
            #
            # 基本带号计算（每6°一个带）
            zone_number = int((lon + 180) // 6) + 1

            # 处理特殊区域：挪威和Svalbard
            # 挪威北部（北纬56°~64°，东经3°~12°）使用32V
            if (56 <= lat < 64) and (3 <= lon < 12):
               zone_number = 32

            # Svalbard群岛特殊带号
            if lat >= 72:
               if 0 <= lon < 9:
                  zone_number = 31
               elif 9 <= lon < 21:
                  zone_number = 33
               elif 21 <= lon < 33:
                  zone_number = 35
               elif 33 <= lon < 42:
                  zone_number = 37

            # 确定半球字母（北半球N，南半球S，特殊区域可能有其他字母）
            if lat >= 0:
               # 北半球：除特殊区域外一般为N
               # 这里简化处理，完整字母表可参考UTM规范
               zone_letter = 'N'
            else:
               # 南半球：除特殊区域外一般为S
               zone_letter = 'S'
        else:
            #
            zone_number = None
            zone_letter = None
        return zone_number, zone_letter
    def _create_custom_crs(self):
        """Create a custom Transverse Mercator CRS based on the specified parameters"""
        # Build PROJ string for the custom projection
        proj_str = (
            f"+proj=tmerc +lat_0={self.latitude_of_origin} "
            f"+lon_0={self.central_longitude} +k={self.scale_factor} "
            f"+x_0={self.false_easting} +y_0={self.false_northing} "
            f"+datum={self.datum} +units=m +no_defs"
        )

        # Adjust false northing for Southern Hemisphere if not already set
        if self.latitude_of_origin < 0 and self.false_northing == 0:
            proj_str = proj_str.replace("+y_0=0", "+y_0=10000000")

        return CRS.from_proj4(proj_str)

    def to_utm(self, lat, lon):
        """
        Convert geographic coordinates (lat, lon) to custom UTM projection coordinates.

        Args:
            lat: Latitude in degrees
            lon: Longitude in degrees

        Returns:
            tuple: (easting, northing) coordinates in meters
        """
        #
        easting, northing = self.to_utm_transformer.transform(lat, lon)
        return easting, northing

    def from_utm(self, easting, northing):
        """
        Convert custom UTM projection coordinates back to geographic coordinates.

        Args:
            easting: Easting coordinate in meters
            northing: Northing coordinate in meters

        Returns:
            tuple: (latitude, longitude) in degrees
        """
        lat, lon = self.from_utm_transformer.transform(easting, northing)
        return lat, lon

    def get_proj_info(self):
        """Get the projection information as a WKT string"""
        return self.crs.to_wkt(pretty=True)
    #
    def ll2local(self, lats,lons):
        #
        east_m  = np.copy(lons)*0
        north_m = np.copy(lats)*0
        #
        #
        for i in range(lats.shape[0]):
            ce,cn = self.to_utm(lats[i],lons[i])
            #
            east_m[i] = ce - self.origin_utm[0]
            north_m[i] = cn - self.origin_utm[1]
            #
        #
        return east_m, north_m
    #
    #
    def local2ll(self, east_m, north_m):
        #
        lats  = np.copy(east_m)*0
        lons  = np.copy(east_m)*0
        #
        #
        for i in range(lats.shape[0]):
            clat,clon = self.from_utm(east_m[i] + self.origin_utm[0], north_m[i] + self.origin_utm[1])
            #
            lats[i] = clat
            lons[i] = clon
        #
        return lats,lons 


# Example usage
if __name__ == "__main__":

    # Example 1: Create a custom UTM projection for New York area
    # (central meridian approximately 75°W)
    origin = [33.65,96.25]  # Note: negative for west longitude
    converter = CustomUTMConverter(origin=origin)

    # Convert coordinates of Times Square
    lons = np.linspace(96.2,96.3,num=6)
    lats = np.linspace(33.6,33.7,num=6)
    #
    easting, northing = converter.ll2local(lats,lons)
    olats, olons      = converter.local2ll(easting, northing)
    #
    print(easting)
    print(northing)
    #
    print("Lons:")
    print(lons)
    print(olons)
    print("Lats:")
    #
    print(lats)
    print(olats)
