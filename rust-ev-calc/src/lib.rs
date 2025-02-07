use wasm_bindgen::prelude::*;
use serde::Deserialize;
use serde_json::{Map, Value as JsonValue};
use log::Level;
use wasm_bindgen_futures::future_to_promise;
use wasm_bindgen_futures::js_sys;

use geo::{Coord, CoordsIter, MultiPoint};
use geo::algorithm::concave_hull::ConcaveHull;
use geojson::{Feature, FeatureCollection, GeoJson, Geometry, Value};
use geo::MapCoords;
use delaunator::{triangulate, Point};
use std::collections::HashSet;
use serde_json::json;

pub mod graph;
pub mod osm_fetcher;
pub mod types;

use self::graph::reachable_points;
use self::types::{Location, BoundingCoordinates};

#[derive(Deserialize)]
pub struct JsLocation {
    pub lat: f64,
    pub lng: f64,
}

#[derive(Deserialize)]
pub struct EVParams {
    pub consumption_wh_per_km: f64,
    pub battery_size_kwh: f64,
    pub current_charge_kwh: f64,
    pub regen_efficiency: f64,
    pub weight: f64,
    pub start_location: JsLocation,
    pub view_area: BoundingCoordinates,
}

#[wasm_bindgen]
pub fn rust_init() {
    console_log::init_with_level(Level::Error).expect("error initializing logger");
    log::info!("Logger initialized from library");
}

impl Location {
    fn to_tuple(&self) -> [f64; 2] {
        [self.longitude, self.latitude]
    }
    fn to_vec(&self) -> Vec<f64> {
        vec![self.longitude, self.latitude]
    }
    fn to_point(&self) -> Value {
        Value::Point(self.to_vec())
    }
}
/*
fn expand(tree: Tree) -> Vec<[f64; 2]> {
    let mut coordinates = Vec::new();
    coordinates.push(tree.location.to_tuple());
    for child in tree.children {
        coordinates.push(tree.location.to_tuple());
        coordinates.extend(expand(child));
    }
    coordinates
}
*/

// Helper: compute the circumradius of a triangle given three points.
fn circumradius(a: &Point, b: &Point, c: &Point) -> f64 {
    let ab = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
    let bc = ((c.x - b.x).powi(2) + (c.y - b.y).powi(2)).sqrt();
    let ca = ((a.x - c.x).powi(2) + (a.y - c.y).powi(2)).sqrt();
    let s = (ab + bc + ca) / 2.0;
    let area = (s * (s - ab) * (s - bc) * (s - ca)).sqrt();
    if area == 0.0 {
        return f64::INFINITY;
    }
    ab * bc * ca / (4.0 * area)
}

// Computes the signed area of a ring of coordinates.
// The ring is assumed to be closed (first point equals last).
fn signed_area(ring: &[Vec<f64>]) -> f64 {
    let mut area = 0.0;
    // Iterate over each segment except the last duplicate point.
    for i in 0..ring.len()-1 {
        let (x1, y1) = (ring[i][0], ring[i][1]);
        let (x2, y2) = (ring[i+1][0], ring[i+1][1]);
        area += (x1 * y2) - (x2 * y1);
    }
    area / 2.0
}

// Ensure the ring is counterclockwise.
// If not, reverse the order (except keep the first point at the start and end).
fn enforce_counterclockwise(mut ring: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    if signed_area(&ring) < 0.0 {
        // Remove the duplicate last point, reverse, then close again.
        ring.pop();
        ring.reverse();
        // Close the ring by adding the first coordinate.
        ring.push(ring[0].clone());
    }
    ring
}

pub async fn compute_range_async(params_json: &str) -> Result<String, String> {
    let params: EVParams = serde_json::from_str(params_json)
        .map_err(|e| e.to_string())?;

    let remaining_energy_kwh = params.current_charge_kwh;
    let consumption_wh_per_km = params.consumption_wh_per_km;
    let origin = Location{latitude: params.start_location.lat, longitude: params.start_location.lng};
    
    let (radius_m, locations) = reachable_points(&origin, params.view_area, remaining_energy_kwh * 1000.0, consumption_wh_per_km / 1000.0)
        .await
        .map_err(|e| e.to_string())?;

    let lines: Vec<Vec<Vec<f64>>> = locations
        .iter()
        .map(|pair| pair.iter().map(|loc| loc.to_vec()).collect())
        .collect();

    let earth_radius_m = 6371000.0;
    let lat1 = origin.latitude.to_radians();
    let lon1 = origin.longitude.to_radians();
    let d = radius_m / earth_radius_m;
    let steps = 36;
    let mut circle_coords = Vec::with_capacity(steps + 1);
    // Change direction to counterclockwise (right-hand rule) by iterating in reverse
    for i in (0..=steps).rev() {
        let bearing = 2.0 * std::f64::consts::PI * (i as f64) / (steps as f64);
        let lat2 = (lat1.sin() * d.cos() + lat1.cos() * d.sin() * bearing.cos()).asin();
        let lon2 = lon1
            + (bearing.sin() * d.sin() * lat1.cos())
                .atan2(d.cos() - lat1.sin() * lat2.sin());
        circle_coords.push(vec![lon2.to_degrees(), lat2.to_degrees()]);
    }
    circle_coords.push(circle_coords[0].clone());


    // Calculate bounding box [min_x, min_y, max_x, max_y]
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    // Update bounds from lines
    for line in &lines {
        for coord in line {
            min_x = min_x.min(coord[0]);
            min_y = min_y.min(coord[1]);
            max_x = max_x.max(coord[0]);
            max_y = max_y.max(coord[1]);
        }
    }

    // Also update bounds from circle coordinates
    for coord in &circle_coords {
        min_x = min_x.min(coord[0]);
        min_y = min_y.min(coord[1]);
        max_x = max_x.max(coord[0]);
        max_y = max_y.max(coord[1]);
    }

    let line_feature = Feature {
        geometry: Some(Geometry::new(Value::MultiLineString(lines))),
        ..Default::default()
    };
    let circle_feature = Feature {
        geometry: Some(Geometry::new(Value::Polygon(vec![circle_coords]))),
        ..Default::default()
    };

    let geojson = GeoJson::FeatureCollection(FeatureCollection {
        bbox: Some(vec![min_x, min_y, max_x, max_y]),
        features: vec![line_feature, circle_feature],
        foreign_members: None,
    });

    let geojson_str = serde_json::to_string(&geojson)
        .map_err(|e| e.to_string())?;
    log::info!("GeoJson created: {}", geojson_str);

    Ok(geojson_str)
}

#[wasm_bindgen]
pub fn compute_range(params: String) -> js_sys::Promise {
    future_to_promise(async move {
        match compute_range_async(&params).await {
            Ok(geojson_str) => Ok(JsValue::from_str(&geojson_str)),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    })
}

