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
    projection :  WebMercator
}

impl GeoJsonReader {
    pub fn new() -> GeoJsonReader {
        GeoJsonReader { features: 0, projection: WebMercator {} }
    }
    pub fn read_from_file(&mut self, path: &String) -> Result<(), Error> {
        println!("Reading file: {}", path);
        let mut file = File::open(path).unwrap();
        let mut data = String::new();
        use std::io::Read;
        file.read_to_string(&mut data).unwrap();
        self.read_from_string(&data)
    }
    fn get_type  (v :  &Value) -> String {
        let t = v["type"].to_string();
        t.trim_matches('"').to_string()
    }
    pub fn read_from_string(&mut self, data: &String) -> Result<(), Error> {

        let v : Value = serde_json::from_str(&data)?;
        let typ = GeoJsonReader::get_type (&v);
        eprintln!("Type '{}' {}", typ, typ == "FeatureCollection");

        match &typ as &str {
            "FeatureCollection" => {
                let features = &v["features"];
                let feature_array = features.as_array().unwrap();
                self.features = 0;
                for feature in feature_array {
                    self.read_feature_collection(feature);
                }
            },
            _ => eprintln!("Expected a FeatureCollection, not: {}", typ)
        }
        Ok(())
    }

    fn read_feature_collection(&mut self, feature : &Value) {
        self.features = self.features + 1;
        let geometry = &feature["geometry"];
        let typ = GeoJsonReader::get_type (&geometry);
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
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            "LineString" => {
                // For type "LineString", the "coordinates" member must be an array of two or more positions.
                eprint!("Handling of GeoJson type '{}' not yet implemented", typ);
            },
            _ => eprintln!("Unexpected type: {}", typ)
        }
    }
    fn read_multi_polygon(&self, geometry: &Value) {
        let mut pp = 0;
        let polygons = geometry["coordinates"].as_array().unwrap();
        for polygon in polygons {
            pp = pp + 1;
            // eprintln!("Feature: {} Polygon: {}", self.ff, pp);
            let ring_array = polygon.as_array().unwrap();
            let mut rr = 0;
            for ring in ring_array {
                rr = rr + 1;
                // eprintln!("Feature: {} Polygon: {} Ring: {}", self.ff, pp, rr);
                let mut cc = 0;
                for coordinate in ring.as_array().unwrap() {
                    // From the http://geojson.org/geojson-spec.html
                    // * x, y, z order (easting, northing, altitude for coordinates in a projected coordinate reference system,
                    // * or longitude, latitude, altitude for coordinates in a geographic coordinate reference system).
                    //
                    let lon = coordinate[0].as_f64().unwrap();
                    let lat = coordinate[1].as_f64().unwrap();
                    let lonlat = LonLatD::new(lon, lat).to_radians();
                    use libwebmap::webmap::MapProjection;
                    let point = self.projection.to_point_xy(lonlat);
                    let command = match cc == 1 { true => "moveto", false => "drawto" };
                    println!("{}    {}  {}", command, round7(point.x), round7(point.y));
                    // eprintln!("Coordinate {} ->  lonlat={}   ->    point={}", coordinate, lonlat.round7(), point.round7());
                    cc = cc + 1;
                }
            }
        }
    }
}
