extern crate libwebmap;

use libwebmap::webmap::LonLatD;
use libwebmap::webmap::WebMap;
use libwebmap::webmap::WebMercator;

///
/// For zoom level one displays a table of the range of longitude and latitudes to x,y
///
fn main() {

    let steps = 20;
    let levels = 10;
    let lon_increment = 360.0 / steps as f64;
    let lat_increment = 180.0 / steps as f64;

    print!("Step  Lon   Lat ");
    for zoom in 0..levels {
        print!("    {:4}x   {:4}y", zoom, zoom);
    }
    println!();
    for step in 0..steps+1 {
        let lon = -180. + lon_increment * step as f64;
        let lat = 90. - lat_increment * step as f64;
        print!("{:4}  {:4}  {:4}", step, lon, lat);
        let lonlat = LonLatD::new(lon, lat).to_radians();
        for zoom in 0..levels {
            let projection = WebMercator {};
            let webmap = WebMap::new(zoom, Box::new(projection));
            let point = webmap.to_point_xy(lonlat);
            let (x, y) = point.to_pixel_coordinates();
            print!("  {:5}    {:5}", x, y);
        }
        println!();
    }
}
