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
