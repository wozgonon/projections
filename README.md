# Projections

Provides a number of map projection.


## Build

The build uses the standard [Rust](https://rust-lang.org)  [Cargo](https://doc.rust-lang.org/cargo/) package manager.

```
$ cargo build               # To build the library
$ cargo build --examples    # To just build the examples
$ cargo test                # To build and run all tests
$ cargo test --doc          # To just run the doctests
$ cargo install
```


## Web Mercator


### A Geojson example

To parse a [GeoJson](http://geojson.org/) file and project the coordinates to a web mercator projection:

```
$ cargo build --examples && ./target/debug/examples/geojson_to_web_mercator <file name>[.geo].json
```

### Javascript


How to render a picture of [Scotland](doc/scotland1.png):

```
$ wget https://github.com/martinjc/UK-GeoJSON/blob/master/json/electoral/sco/wpc.json
$ cargo build --examples && ./target/debug/examples/geojson_to_web_mercator data/wpc.json | ./examples/draw_to_js.py > examples/_scotland.js
```

Then place the following URL in your web browser:
* file:///<...your path...>/projections/examples/scotland_canvas.html


## Mercator

