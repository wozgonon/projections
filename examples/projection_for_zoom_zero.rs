extern crate libwebmap;

use libwebmap::webmap::WebMap;
use libwebmap::webmap::LonLatD;
use libwebmap::webmap::WebMercator;

///
/// For zoom level one displays a table of the range of longitude and latitudes to x,y
///
fn main() {

    let projection = WebMercator {};
    let zoom = 0;
    let webmap = WebMap::new (zoom, Box::new(projection));
    let steps = 20;
    let lon_increment = 360.0 / steps as f64;
    let lat_increment = 180.0 / steps as f64;

    for step in 0..steps+1 {
        let lon = -180. + lon_increment * step as f64;
        let lat = 90. - lat_increment * step as f64;
        let lonlat = LonLatD::new(lon, lat).to_radians();
        let point = webmap.to_point_xy(lonlat);
        let (x, y) = point.to_pixel_coordinates ();
        println!("{:2}.   lon:{:4} -> x:{:.4} ({:3})  lat:{:4} -> y:{:.4} ({:3})", step, lon, point.x, x, lat, point.y, y);
    }
}
