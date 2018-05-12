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
