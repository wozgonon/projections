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

extern crate serde_json;
extern crate serde;
extern crate libwebmap;

//use libwebmap::webmap::WebMap;
use libwebmap::webmap::LonLatD;
use libwebmap::webmap::WebMercator;
use libwebmap::webmap::round7;
use std::env::args;

///
/// For documentation see:
/// * https://docs.serde.rs/serde/index.html
/// * https://crates.io/crates/serde_json
/// * https://crates.io/crates/serde
///
use serde_json::{Value, Error};


use std::fs::File;
use std::f64::INFINITY;

fn main() {

    let mut count = 0;
    for arg in  args () {
        if count > 0 {
            let path = arg;  // Name of a file containing Geojson
            let mut reader = GeoJsonReader::new ();

            match reader.read_from_file(&path) {
                Ok(_) => eprintln!("Read file successfully: {}", path),
                Err(message) => eprintln!("Error while reading geojson from file {} : {}", path, message),
            }
        }
        count = count + 1;
    }
}


struct GeoJsonReader {
    features: u32,
    max_x : f64,
    max_y : f64,
    min_x : f64,
    min_y : f64,
    projection :  WebMercator
}

impl GeoJsonReader {
    pub fn new() -> GeoJsonReader {
        GeoJsonReader { features: 0,
            max_x : -INFINITY,
            max_y : -INFINITY,
            min_x : INFINITY,
            min_y : INFINITY,
            projection: WebMercator {} }
    }
    pub fn read_from_file(&mut self, path: &String) -> Result<(), Error> {
        eprintln!("Reading file: {}", path);
        let mut file = File::open(path).unwrap();
        let mut data = String::new();
        use std::io::Read;
        file.read_to_string(&mut data).unwrap();
        self.read_from_string(&data)
    }
    pub fn read_from_string(&mut self, data: &String) -> Result<(), Error> {
        let v: Value = serde_json::from_str(&data)?;
        let typ = GeoJsonReader::get_type(&v);
        eprintln!("Type '{}' {}", typ, typ == "FeatureCollection");

        match &typ as &str {
            "FeatureCollection" => {
                let features = &v["features"];
                let feature_array = features.as_array().unwrap();
                self.features = 0;
                for feature in feature_array {
                    self.read_feature_collection(feature);
                }
                self.print ("topleft", self.min_x, self.min_y);
                self.print ("bottomright", self.max_x, self.max_y);
            },
            _ => eprintln!("Expected a FeatureCollection, not: {}", typ)
        }
        Ok(())
    }
    ///
    /// Pulls out the value of the 'type' attribute and removes
    /// unnecessary quotes from either end of the vslue.
    ///
    fn get_type(v: &Value) -> String {
        let t = v["type"].to_string();
        t.trim_matches('"').to_string()
    }

    fn read_feature_collection(&mut self, feature: &Value) {
        self.features = self.features + 1;
        let geometry = &feature["geometry"];
        let typ = GeoJsonReader::get_type(&geometry);
        eprintln!("Feature: {}  geometry type={}", self.features, typ);
        match typ.as_ref() {
            "MultiPolygon" => {
                self.read_multi_polygon(geometry)
            },
            "MultiLineString" => {
                // For type "MultiLineString", the "coordinates" member must be an array of LineString coordinate arrays.
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            "Point" => {
                //  For type "Point", the "coordinates" member must be a single position.
            },
            "MultiPoint" => {
                //For type "MultiPoint", the "coordinates" member must be an array of positions.
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            "GeometryCollection" => {
                // A GeoJSON object with type "GeometryCollection" is a geometry object which represents a collection of geometry objects.
                // A geometry collection must have a member with the name "geometries".
                // The value corresponding to "geometries" is an array.
                // Each element in this array is a GeoJSON geometry object.
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            "Polygon" => {
                self.read_polygon(geometry)
            },
            "LineString" => {
                // For type "LineString", the "coordinates" member must be an array of two or more positions.
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            _ => eprintln!("Unexpected type: {}", typ)
        }
    }
    fn read_multi_polygon(&mut self, geometry: &Value) {
        let mut pp = 0;
        let polygons = geometry["coordinates"].as_array().unwrap();
        for polygon in polygons {
            pp = pp + 1;
            // eprintln!("Feature: {} Polygon: {}", self.ff, pp);
            let ring_array = polygon.as_array().unwrap();
            self.read_rings(ring_array);
        }
    }
    fn read_polygon(&mut self, geometry: &Value) {
        let ring_array = geometry["coordinates"].as_array().unwrap();
        self.read_rings(ring_array);
    }
    fn read_rings(&mut self, ring_array: &Vec<Value>) {
        let mut rr = 0;
        for ring in ring_array {
            rr = rr + 1;
            // eprintln!("Feature: {} Polygon: {} Ring: {}", self.ff, pp, rr);
            let mut cc = 0;
            for coordinate in ring.as_array().unwrap() {
                let lon = coordinate[0].as_f64().unwrap();
                let lat = coordinate[1].as_f64().unwrap();
                //if lon.abs() > 360.0 || lat.abs() > 90.0 {
                    // From the http://geojson.org/geojson-spec.html
                    // * x, y, z order (easting, northing, altitude for coordinates in a projected coordinate reference system,
                //} else {
                    // * or longitude, latitude, altitude for coordinates in a geographic coordinate reference system).
                //}
                let command = match cc == 0 {
                    true => "moveto",
                    false => "drawto"
                };
                self.project_to_xy (command, lon, lat);
                cc = cc + 1;
            }
        }
    }
    fn project_to_xy (&mut self, command : &str, lon : f64, lat : f64) {
        let lonlat = LonLatD::new(lon, lat).to_radians();
        use libwebmap::webmap::MapProjection;
        let point = self.projection.to_point_xy(lonlat);
        if point.x.is_infinite() || point.y.is_infinite() {
            eprintln!("Ignoring point projected to infinity: {}", point);
        } else {
            self.max_x = point.x.max(self.max_x);
            self.min_x = point.x.min(self.min_x);
            self.max_y = point.y.max(self.max_y);
            self.min_y = point.y.min(self.min_y);
            self.print(command, point.x, point.y);
        }
    }
    fn print(&self, command : &str, x : f64, y : f64) {
        let tab = "	";
        println!("{}{}{}{}{}", command, tab, round7(x), tab, round7(y));
        // eprintln!("Coordinate {} ->  lonlat={}   ->    point={}", coordinate, lonlat.round7(), point.round7());
    }
}

