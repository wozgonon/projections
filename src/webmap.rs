//    projections - provide a set of map projections
//
//    Copyright (C) 2018  Wozgonon
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published
//    by the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


use std::f64::consts::PI;
use std::fmt::Display;
use std::fmt;

pub const GREENWICH_MERIDIAN: f64 = 0.0;
pub const INTERNATIONAL_DATE_LINE: f64 = 180.0;
pub const LATITUDE_OF_EQUATOR: f64 = 0.0;
pub const NORTH_MOST_LATITUDE: f64 = 85.051129; // FIXME
pub const SOUTH_MOST_LATITUDE: f64 = -NORTH_MOST_LATITUDE; // FIXME
pub const ONE_EIGHTY_OVER_PI : f64 = 57.29577951308;

///
///  A coordinate, in [Radians](https://en.wikipedia.org/wiki/Radian), on the surface of the earth (or another such body),
///
#[derive(Clone, Copy)]
#[derive(Debug)]
pub struct LonLat {
    lon : f64,
    lat : f64
}

impl LonLat {

    ///
    /// Creates a new reference point.
    ///
    /// The Longitude is -PI..PI radians.
    /// the Latitudue is -PI/2..PI/2 radians.
    ///
    pub fn new (lon : f64, lat : f64) -> LonLat {
        return LonLat { lon : lon % 180.0, lat : lat % 90.0 };
    }
    ///
    ///  Converts this coordinate to one in degrees.
    ///  ```
    ///  use libwebmap::webmap::LonLat;
    ///  use libwebmap::webmap::LonLatD;
    ///  assert_eq!(LonLat::new(PI,-PI).to_degrees (), LonLatD(180.0, -180.0))
    ///  ```
    ///
    pub fn to_degrees (&self) -> LonLatD {
        return LonLatD { lon : LonLat::radians_to_degrees (self.lon), lat : LonLat::radians_to_degrees (self.lat) };
    }

    ///
    /// Convert an angle in [Radians](https://en.wikipedia.org/wiki/Radian) to one in degrees.
    /// ```
    /// use libwebmap::webmap::radians_to_Degrees;
    /// assert_eq(radians_to_Degrees (0.0), 0.0);
    /// assert_eq(radians_to_Degrees (PI), 180.0);
    /// assert_eq(radians_to_Degrees (PI/2.0), 90.0);
    /// assert_eq(radians_to_Degrees (-PI), -180.0);
    /// ```
    ///
    #[inline]
    fn radians_to_degrees(radians : f64) -> f64 {
       return radians / 180.0 * PI;
    }

    ///
    /// Google Maps cannot show the poles because Mercator projects the poles at infinity.
    /// It cuts off at 85.051129 north and south (the latitudes at which the full map becomes a square).
    ///  ```
    ///  assert!((max_latitude() - north_most_latitude) < .0001);
    /// ```
    ///
    pub fn max_latitude () {
       2.0 * PI.exp().atan() - PI/2.0;
    }

    ///
    ///  Makes a [regular convex polygon](https://en.wikipedia.org/wiki/Regular_polygon) of 'n' sides.
    ///  * n=3     => triangle
    ///  * n=4     => square
    ///  * n=5     => pentagram
    ///  * n>~20   => approximates a cirle to the naked eye when plotted on a screen.
    ///  * n-> inf => circle
    ///
    /// ## Examples
    ///
    /// ```
    /// use libwebmap::webmap;
    /// use libwebmap::webmap::LonLat;
    /// assert_eq!(LonLat::make_regular_convex_polygon(0, 10.).len(), 0);
    /// assert_eq!(LonLat::make_regular_convex_polygon(1, 10.).len(), 1);
    /// assert_eq!(LonLat::make_regular_convex_polygon(3, 10.).len(), 3);
    /// assert_eq!(LonLat::make_regular_convex_polygon(11, 10.).len(), 11);
    /// ```
    ///
   pub fn make_regular_convex_polygon (sides : u32, radius : f64) -> Vec<LonLat> {
       let mut result = Vec::new ();
       let angle = 2.0 * PI / sides as f64;
       for _ in 0.. sides {
           let x = radius * angle.cos ();
           let y = radius * angle.cos ();
           result.push (LonLat::new (x,y));
       }
       result
   }
}

///
///  A coordinate, in degrees, on the surface of the earth (or another such body).
///
#[derive(Debug)]
pub struct LonLatD {
    lon : f64,  // Degrees
    lat : f64   // Degrees
}

///
///  A coordinate, in [Degrees](https://en.wikipedia.org/wiki/Degree_(angle)), on the surface of the earth (or another planet,
///
impl LonLatD {
    ///
    /// Creates a new reference point.
    /// The Longitude is modified to be in the range -180..180 degrees and
    /// the Latitudue is modified to be in the range -90..90 degrees.
    ///
    pub fn new (lon : f64, lat : f64) -> LonLatD {
        return LonLatD { lon : lon % 180.0, lat : lat % 90.0 };
    }
    ///
    ///  Converts this coordinate to one in radians.
    ///  ```
    ///  assert_eq!(LonLatD::new(90.0,-90.0).to_radians (), LonLat(PI/2.0, -PI/2.0))
    ///  ```
    ///
    pub fn to_radians (&self) -> LonLat {
        return LonLat { lon : LonLatD::degrees_to_radians (self.lon), lat : LonLatD::degrees_to_radians (self.lat) };
    }

    ///
    /// Convert an angle in degrees to one in [Radians](https://en.wikipedia.org/wiki/Radian).
    /// ```
    /// use libwebmap::webmap::degrees_to_Radians;
    /// assert_eq(degrees_to_Radians (0.0), 0.0);
    /// assert_eq(degrees_to_Radians (180.0), PI);
    /// assert_eq(degrees_to_Radians (90.0), PI/2.0);
    /// assert_eq(degrees_to_Radians (-180), -PI);
    /// ```
    ///
    #[inline]
    fn degrees_to_radians(degrees : f64) -> f64 {
       return degrees * 180.0 / PI;
    }
}

///
///  A point of a 2D surface.
///
#[derive(Clone, Copy)]
#[derive(Debug)]
#[derive(PartialEq)]
pub struct PointXY {
    x : f64,
    y : f64
}

impl Display for PointXY {
     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
     }
}

impl PointXY {
    pub fn new(x: f64, y: f64) -> PointXY {
        return PointXY { x: x, y: y };
    }
    ///
    /// Find the difference between two points.
    ///
    /// ## Examples
    ///
    /// ```
    /// use libwebmap::webmap::PointXY;
    /// let (dx, dy) = PointXY::new(50., 25.).minus(PointXY::new(25., 50.));
    /// assert_eq!(dx, 25.);
    /// assert_eq!(dy, -25.);
    /// ```
    ///
    pub fn minus (&self, end : PointXY) -> (f64,f64) {
       let dy = self.y - end.y;
       let dx = self.x - end.x;
       (dx, dy)
    }
}

///
///  A list of points of a 2D surface.
///
///struct PointsXY {
///    points : [PointXY]
///}

///
///  Returns some notion of 'distance' between a point and a line,
///  see [Metric or distance function](https://en.wikipedia.org/wiki/Metric_(mathematics)).
///
///  Note that the distance between two points on the Earth's surface is not the same as the
///  familiar straight line or euclidean distance.
///  Nor, strictly speaking, is it the [Spherical Distance](http://mathworld.wolfram.com/SphericalDistance.html)
///  since the Earth is not actually a sphere.
///
///  This 'distance' could be defined in a variety of ways: it could for simplicity be perpidicular euclidean distance
///  between the point and line or the distance between the point and the mid-point of the line.   It could defined
///  as being related to the area of the triangle formed by the point and the two end-points of the line.
///

pub trait Metric {
    fn distance (&self, point : &PointXY, start : &PointXY, end : &PointXY) -> f64;
}

///
///  A particular simple implementation of the Metric [trait](https://doc.rust-lang.org/book/traits.html) that
///  defines the distance between a point and a line as the perpidicular euclidean distance.
///
pub struct PerpendicularDistance {}

impl Metric for PerpendicularDistance {

   ///
   /// Distance between a point and the [closest point a line segement](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line).
   /// (Note that this formula will only be accurate for small distances on a sphere).
   ///
   ///  # Examples
   ///
   /// ```
   /// use libwebmap::webmap::PointXY;
   /// use libwebmap::webmap::PerpendicularDistance;
   /// use libwebmap::webmap::Metric;
   /// let metric = PerpendicularDistance {};
   /// assert_eq!(metric.distance(&PointXY::new(1.,0.),&PointXY::new(0.,0.),&PointXY::new(0.,1.)), 1.0);
   /// assert_eq!(metric.distance(&PointXY::new(0.,0.),&PointXY::new(0.,0.),&PointXY::new(0.,1.)), 0.0);
   /// assert_eq!(metric.distance(&PointXY::new(2.,2.),&PointXY::new(0.,0.),&PointXY::new(4.,0.)), 2.0);  // -2.0
   /// ```

    fn distance (&self, point : &PointXY, start : &PointXY, end : &PointXY) -> f64 {
        let dy = end.y - start.y;
        let dx = end.x - start.x;
        let dxy = (dy * point.x - dx * point.y - end.x*start.y - end.y*start.x) / ((dy*dy + dx*dx).sqrt());
	dxy.abs ()
    }
}

///
///  A boundary or a river may be represented by a sequence of coordinates,
///  which can be seen as a Polygon or [PolyLine](https://en.wikipedia.org/wiki/Polygonal_chain).
///  For some situations,
///  such as the coastline of large countries or continents, these sequences
///  can run into hundreds of thousands or even millions of points which may
///  have to be copied and translated even when such detail is not required.
///
///  For this reason it is very useful to simplify such a sequence by removing
///  unnecessarilly detailed points.
///

pub trait Simplifier {
    ///
    ///  Return a subset of the input points.
    ///  Typically a simplifier should not try to smooth or modify the set of
    ///  inputs, it should return a subset of the inputs so that maps appear to
    ///  align properly on-screen.
    ///
    fn simplify (&self, points : &[PointXY], result :  &mut Vec<PointXY>);
}

///
///  An implmentation of the [Ramer Douglas Peucker](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm)
///  algorithm for simplifying a Polygon or Polyline.
///

pub struct RamerDouglasPeuker<'a> {
    epsilon : f64,
    metric : Box<Metric + 'a>
}

impl<'a> RamerDouglasPeuker<'a> {
    pub fn new (epsilon : f64, metric : Box<Metric + 'a>) -> RamerDouglasPeuker<'a> {
    	RamerDouglasPeuker { epsilon : epsilon, metric : metric }
    }
}

impl<'a> Simplifier for RamerDouglasPeuker<'a> {

    ///
    ///  Returns a subset of the input coordinates.   This is based on divide and conquer recursion,
    ///  splitting the input array in half, then half again and so on,
    ///  (somewhat reministant of the quick sort algorithm), with an algorithmic time complexity of O(nlog(n)).
    ///
    fn simplify (&self, points : &[PointXY], result :  &mut Vec<PointXY>) {
        let mut max_distance = 0.0;
        let mut index = 0;
        let end = points.len() - 1;
        for i in 1..end {
            let distance = self.metric.distance(&points[i], &points[1], &points[end]);
            if distance > max_distance {
                index = i;
                max_distance = distance;
            }
        }
        if max_distance as f64 > self.epsilon {
            self.simplify(&points[0..index], result);
            self.simplify(&points[index..end], result);
        } else {
            result.push(points[1]);
            result.push(points[end]);
        }
    }
}


///
/// A [map projection](https://en.wikipedia.org/wiki/Map_projection) provides a systematic transformation
/// from coordinates on the Earth (or any other such body, such as sphere or elipsoid) to the plane (a flat surface).
/// Esssentially it converts latitude, longitude coordinates to XY coordinates.
///
/// There are many such map projections, [see here for a list](https://en.wikipedia.org/wiki/List_of_map_projections).
///

pub trait MapProjection {
    fn to_point_xy(&self, lonlat : LonLat) -> PointXY;
    fn to_lonlat (&self, point : PointXY) -> LonLat;
    fn centre (&self) -> LonLat { LonLat::new (PI, 0.0) }
//       let centre = LonLat::new(PI, 0);    // FIXME PI is centre.x or lambda_0  - MapProjection::centre().lon
}

///
///  The [Mercator projection](https://en.wikipedia.org/wiki/Mercator_projection) is perhaps the most familiar
///  projection, it models the earth is an elipsoid (rather than a sphere).
///
///  It is a [conformal map](https://en.wikipedia.org/wiki/Conformal_map)
///  which means that angles remain unchanged which makes it useful for navigation though shapes and areas in particular are grossly distorted
///  particularly towards the poles.
///
///  It is a cylindrical projection, so:
///  * horizontal lines on the map represent lines of constant latitude, or [parallels](https://en.wikipedia.org/wiki/Circle_of_latitude).
///  * vertical lines on the map represent lines of constant longitude, or [meridians]) (https://en.wikipedia.org/wiki/Meridian_(geography)).
///
pub struct Mercator {}

impl MapProjection for Mercator {

    ///
    ///   The [Mercator projection](http://mathworld.wolfram.com/MercatorProjection.html)
    ///
    fn to_point_xy(&self, lonlat : LonLat) -> PointXY {
        //let e = 0.006.sqrt();
        //let y = ((lonlat.lat/2.0 + PI/4.0).tan()*(1-e*lonlat.lat.sin())/(1+e*lonlat.lat.sin())).ln();
        let x = lonlat.lon-self.centre().lon;
	let y = (lonlat.lat.tan () + 1.0/lonlat.lat.cos()).ln();
        PointXY::new(x, y)
    }
    ///
    ///  The [inverse mercator projection](https://www.johndcook.com/blog/2009/09/21/gudermannian/).
    ///  Formula from: [mercator projection](http://mathworld.wolfram.com/MercatorProjection.html)
    ///
    ///  The formula for the latitude is actually the [Gudermannian](http://mathworld.wolfram.com/Gudermannian.html).
    ///
    fn to_lonlat (&self, point : PointXY) -> LonLat {
        let lon = point.x+self.centre().lon;
	let lat = point.y.sinh().atan();
	return LonLat::new (lon, lat);
    }
}

///
///  The WebMercator projection has become standard for online maps, the formulas are based on those for the mercator
///  but are simpler because they model the earth as a Sphere rather than as an Elipsoid.
///  On screen they look very similar and because the formulas are simpler they will be quicker to render,
/// however they should not be used for practical navigation because of the [discrepencies](https://thatsmaths.files.wordpress.com/2015/05/wm-vs-merc-detail.jpg).
///

pub struct WebMercator {}

impl MapProjection for WebMercator {

    ///
    ///  Performs the basic web mercator transformation, not taking zoom or tiling into account.
    ///
    fn to_point_xy(&self, lonlat : LonLat) -> PointXY {
        let x = lonlat.lon + PI;  // FIXME PI is centre.x or lambda_0
        let y = (lonlat.lat/2.0 + PI/4.0).tan().ln();
        PointXY::new(x, y)
    }
    ///
    ///  Performs the basic inverse web mercator transformation, not taking zoom or tiling into account.
    ///
    fn to_lonlat (&self, point : PointXY) -> LonLat {
        let lon = point.x - PI;  // FIXME PI is centre.x or lambda_0
        let lat = 2.0*(point.y.exp().atan() - PI/4.0);
        return LonLat::new (lon, lat);
    }
}

///
///  The [Sinusoidal_projection](https://en.wikipedia.org/wiki/Sinusoidal_projection).
///
///  * As on the earth, the length of each parallel (line of latitude) is proportional to the cosine of the latitude.
///  * Each meridian (line of longitude) is half of a sine wave.
///
pub struct Sinusoidal{}

impl MapProjection for Sinusoidal {
    fn to_point_xy(&self, lonlat : LonLat) -> PointXY {
        PointXY::new (lonlat.lon - lonlat.lon-self.centre().lon, lonlat.lat)
    }
    fn to_lonlat (&self, point : PointXY) -> LonLat {
        LonLat::new (point.y + self.centre().lon, point.y)
    }
}

/// ///
/// /// https://en.wikipedia.org/wiki/Rectangular_polyconic_projection
/// ///
/// struct RectangularPolyconicProjection  {}
/// 
/// impl MapProjection for RectangularPolyconicProjection {
///     fn to_point_xy(&self, lonlat : LonLat) -> PointXY {
///         let A = (0.5*(lonlat.lon - MapProjection::centre().lon) * MapProjection::centre.lat.sin()).tan()/MapProjection::centre.lat.sin();
///         let E = 2.0 * (A * lonlat.lat.sin ()).atan ();
///         let x = lonlat.lat.tan().recip()*E.sin ();
///         let y = lonlat.lat - MapProjection::centre.lat () + (1.0 + E.cos ()) / lonlat.lat;
///         PointXY::new (x, y)
///     }
///     fn to_lonlat (&self, point : PointXY) -> LonLat {
///         LonLat::new (point.y + MapProjection::centre().lon, point.y)
///     }
/// }

///
/// Before applying zoom, the "world coordinates" are adjusted such that the upper left corner
///  is (0, 0) and the lower right corner is (256, 256)  
///
///
/// http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/

pub struct WebMap<'a> {
    projection : Box<MapProjection + 'a>,
    zoom_scaling : f64
}

impl<'a> WebMap<'a> {

    ///
    ///  Creates a new WebMap projection.
    ///
    pub fn new (zoom : u32, projection : Box<MapProjection + 'a>) -> WebMap<'a> {
        let zoom_scaling =  2.0_f64.powi(zoom as i32) / 2.0 / PI;
        WebMap { projection : projection, zoom_scaling : zoom_scaling }
    }

    pub  fn simplify_lonlats (&self, lonlats : &Vec<LonLat>, simplifier : &Simplifier) -> Box<Vec<LonLat>> {
        let xys = self.lonlats_to_points_xy (lonlats);
        let mut simplified_xys = Vec::new (); //PointXY; len + recResults2.len()];
        simplifier.simplify(&xys, &mut simplified_xys);
        return self.points_xy_to_lonlats (&simplified_xys);
    }

    ///
    /// Project sequence of points to Web Mercator
    ///
    fn lonlats_to_points_xy (&self, points : &Vec<LonLat>) -> Box<Vec<PointXY>> {
        let len = points.len();
        let mut result = Vec::new();
        for ii in 0..len {
            let point = &points[ii];
            let xy = self.to_point_xy (*point);
            result.push(xy);
        }
        return Box::new(result);
    }

    ///
    ///  Convert a sequence of points projected web mercator back to
    ///  a sequence of Longitude, Latutudes.
    ///
    fn points_xy_to_lonlats (&self, points : &Vec<PointXY>) -> Box<Vec<LonLat>> {
        let len = points.len ();
        let mut result = Vec::new();
        for xx in 0..len {
            let point = &points[xx];
            let lonlat = self.projection.to_lonlat (*point);
            result.push (lonlat);
        }
        return Box::new (result);
    }

    ///
    ///  Project the point to a Web Mercator (Pyramid) Coordinate.
    ///
    ///  # Examples
    ///
    /// TODO
    ///  use std::f64::INFINITY;
    ///  use std::f64::NAN;
    ///  use libwebmap::webmap::PointXY;
    ///  use libwebmap::webmap::LonLat;
    ///  use libwebmap::webmap::GREENWICH_MERIDIAN;
    ///  use libwebmap::webmap::INTERNATIONAL_DATE_LINE;
    ///  use libwebmap::webmap::NORTH_MOST_LATITUDE;
    ///  use libwebmap::webmap::SOUTH_MOST_LATITUDE;
    ///  use libwebmap::webmap::LATITUDE_OF_EQUATOR;
    ///  let zoom = 0;
    ///  let projection = WebMercator {};
    ///  let webmap = WebMap::new (zoom, Box::new(projection));
    ///  assert_eq!(webmap.to_point_xy (LonLat::new (GREENWICH_MERIDIAN, equator))), PointXY::new(0.5,0.));
    ///  assert_eq!(webmap.to_point_xy (LonLat::new (90., LATITUDE_OF_EQUATOR).to_point_xy ()), PointXY::new(INFINITY,0.));
    ///  assert_eq!(webmap.to_point_xy (LonLat::new (NORTH_MOST_LATITUDE, INTERNATIONAL_DATE_LINE).to_point_xy ()), PointXY::new(0.,0.));
    ///  assert_eq!(webmap.to_point_xy (LonLat::new (SOUTH_MOST_LATITUDE, INTERNATIONAL_DATE_LINE).to_point_xy ()), PointXY::new(1.,1.));
    /// 
    pub fn to_point_xy (&self, point : LonLat) -> PointXY {
        let xy = self.projection.to_point_xy (point);
        PointXY::new (xy.x * self.zoom_scaling, xy.y * self.zoom_scaling)
    }

    pub fn to_lonlat (&self, point : PointXY) -> LonLat {
        let xy = PointXY::new (point.x /self.zoom_scaling, point.y / self.zoom_scaling);
        self.projection.to_lonlat (xy)
   }
}


#[cfg(test)]
mod tests {

    use std::f64::INFINITY;
    use std::f64::consts::PI;
    use webmap::PointXY;
    use webmap::LonLat;
    use webmap::LonLatD;
    use webmap::WebMap;
    use webmap::WebMercator;
    use webmap::GREENWICH_MERIDIAN;
    use webmap::INTERNATIONAL_DATE_LINE;
    use webmap::NORTH_MOST_LATITUDE;
    use webmap::SOUTH_MOST_LATITUDE;
    use webmap::LATITUDE_OF_EQUATOR;

    macro_rules! assert_float_eq {
       ($a:expr, $b:expr) => {
           let (dx, dy) = $a.minus($b);
	   eprintln!("a={} b={} dx={} dy={}", $a, $b, dx, dy);
           assert_eq!(dx.abs() < 0.000001, true);
           assert_eq!(dy.abs() < 0.000001, true);
       }
    }
    #[test]
    fn should_project_to_web_mercator1 () {
        let zoom = 0;
    	let projection = WebMercator {};
	let webmap = WebMap::new (zoom, Box::new(projection));

        let atlantic : LonLat = LonLatD::new (GREENWICH_MERIDIAN, LATITUDE_OF_EQUATOR).to_radians ();
        let pacific : LonLat = LonLatD::new (INTERNATIONAL_DATE_LINE, LATITUDE_OF_EQUATOR).to_radians ();
	
	assert_float_eq!(webmap.to_point_xy (atlantic), PointXY::new(0.5,0.));
	assert_float_eq!(webmap.to_point_xy (pacific), PointXY::new(0.5,0.));  
    }
    #[test]
    fn should_project_to_web_mercator2 () {
        let zoom = 0;
    	let projection = WebMercator {};
	let webmap = WebMap::new (zoom, Box::new(projection));

	let east_equator : LonLat  = LonLatD::new (-90.0, LATITUDE_OF_EQUATOR).to_radians ();
	let west_equator : LonLat  = LonLatD::new (90.0, LATITUDE_OF_EQUATOR).to_radians ();
	
        assert_float_eq!(webmap.to_point_xy (east_equator), PointXY::new(0.25,0.));  // 0.25,0.5
        assert_float_eq!(webmap.to_point_xy (west_equator), PointXY::new(0.75,0.));  // 0.75,0.5
    }
    #[test]
    fn should_project_to_web_mercator3_extreme_latitudes () {
        let zoom = 0;
    	let projection = WebMercator {};
	let webmap = WebMap::new (zoom, Box::new(projection));

	let north_pacific : LonLat  = LonLatD::new (INTERNATIONAL_DATE_LINE, NORTH_MOST_LATITUDE).to_radians ();
	let south_pacific : LonLat  = LonLatD::new (INTERNATIONAL_DATE_LINE, SOUTH_MOST_LATITUDE).to_radians ();

	use std::f64::NAN;
        assert_float_eq!(webmap.to_point_xy (south_pacific), PointXY::new(1.,1.));  // 1.,1.
        assert_float_eq!(webmap.to_point_xy (north_pacific), PointXY::new(0.5, NAN));  // 1.,0.
    }
    #[test]
    fn should_project_to_web_mercator4_poles () {
        let zoom = 0;
    	let projection = WebMercator {};
	let webmap = WebMap::new (zoom, Box::new(projection));

	let north_pole : LonLat  = LonLatD::new (GREENWICH_MERIDIAN, PI/2.).to_radians ();
	let south_pole : LonLat  = LonLatD::new (GREENWICH_MERIDIAN, PI/2.).to_radians ();

	assert_float_eq!(webmap.to_point_xy (south_pole), PointXY::new(0.5, -INFINITY));
	assert_float_eq!(webmap.to_point_xy (north_pole), PointXY::new(0.5, INFINITY));
    }
}
