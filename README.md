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

## Mercator

