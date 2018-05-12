extern crate serde_json;
extern crate serde;
extern crate libwebmap;

//use libwebmap::webmap::WebMap;
use libwebmap::webmap::LonLatD;
use libwebmap::webmap::WebMercator;
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
            let path = arg;  /// Name of a file containing Geojson
            println!("Reading file: {}", path);
            let mut file = File::open(path).unwrap();
            let mut data = String::new();
            use std::io::Read;
            file.read_to_string(&mut data).unwrap();

            match read_geojson(&data) {
                Ok(v) => println!("ok {:?}", v),
                Err(message) => println!("err {}", message),
            }
        }
        count = count + 1;
    }
}

fn read_geojson(data : &String) -> Result<(), Error> {


    // Parse the string of data into serde_json::Value.
    let v: Value = serde_json::from_str(&data)?;

    let feature_collection = "FeatureCollection";
    let typ = v["type"].to_string();
    println!("Type {}", typ);

    let projection = WebMercator {};

    match typ {
        feature_collection => {
            let features = &v["features"];
            //use serde::de::Deserializer;
            let feature_array = features.as_array().unwrap();
            //for xx = 0..feature_array.length()
            let mut ff = 0;
            for feature in feature_array{
                ff = ff + 1;
                let geometry = &feature["geometry"];
                let typ = geometry["type"].to_string();
                println!("Feature: {}  geometry type={}", ff, typ);
                let multipolygon = "MultiPolygon";
                match typ {
                    multipolygon => {
                        let mut pp = 0;
                        let polygons = geometry["coordinates"].as_array().unwrap();
                        for polygon in polygons {
                            pp = pp + 1;
                            println!("Feature: {} Polygon: {}", ff, pp);
                            let ring_array= polygon.as_array().unwrap();
                            let mut rr = 0;
                            for ring in ring_array {
                                rr = rr + 1;
                                println!("Feature: {} Polygon: {} Ring: {}", ff, pp, rr);
                                for coordinate in ring.as_array().unwrap() {
                                    let lon = coordinate[0].as_f64().unwrap();
                                    let lat = coordinate[1].as_f64().unwrap();
                                    let lonlat = LonLatD::new(lon, lat).to_radians ();
                                    use libwebmap::webmap::MapProjection;
                                    let point = projection.to_point_xy (lonlat);
                                    println!("Coordinate {} ->  lonlat={}   ->    point={}", coordinate, lonlat.round7(), point.round7());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}