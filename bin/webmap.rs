extern crate libwebmap;

use libwebmap::webmap::LonLat;
use libwebmap::webmap::WebMap;
use libwebmap::webmap::WebMercator;
use libwebmap::webmap::RamerDouglasPeuker;
use libwebmap::webmap::PerpendicularDistance;

fn main() {

    let projection = WebMercator {};
    let zoom = 18;
    let webmap = WebMap::new (zoom, Box::new(projection));
    let epsilon = 0.5;
    let metric = PerpendicularDistance {};
    let simplifier = RamerDouglasPeuker::new (epsilon, Box::new (metric));
    for nn in 0..20 {
        let sides = 1 << nn;
	let radius = 100.;
        let polygon = LonLat::make_regular_convex_polygon (sides, radius);
    	let simplified_polygon = webmap.simplify_lonlats (&polygon, &simplifier);
    	println!("Polygon points: {} simplified to: {}", polygon.len (), simplified_polygon.len ());
    }
}


//  GeoJson classes

//struct Polygon {
//    outer : LonLats,
//    inner : LonLats []
//}

//struct MultiLine {
    //   points : LonLats []
//}

//struct MultiPolygon {
    //polygons : Polygon [];
//}
