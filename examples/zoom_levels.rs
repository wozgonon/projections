extern crate libwebmap;

use libwebmap::webmap::WebMap;
use libwebmap::webmap::WebMercator;

///
/// Displays information about each zoom level.
/// (note: the output looks likes a pyramind, hence the term pyramidal coordinates)
///
fn main() {

    let levels = 23;

    println!("{}    {:^43}   {:^}", "Zoom",  "Tiles",   "Pixels");
    for zoom in 0..levels+1 {
        let projection = WebMercator {};
        let webmap = WebMap::new(zoom, Box::new(projection));
        let tiles = webmap.number_of_tiles();
        let pixels = webmap.number_of_pixels();
        println!("{:02}  {:>9} x {:<9}={:16}  {:12} x {:<12}", zoom, tiles.0, tiles.1, tiles.0 as f64 * tiles.1 as f64, pixels.0, pixels.1);
    }
    println!();
}

