use serde::Deserialize;
use std::collections::HashMap;
// ** Coordinates **

#[derive(Clone)]
pub struct Location {
    pub latitude: f64,
    pub longitude: f64,
}

#[derive(Deserialize)]
pub struct BoundingCoordinates {
    pub north_latitude: f64,
    pub south_latitude: f64,
    pub east_longitude: f64,
    pub west_longitude: f64,
}

// ** OSM data types **

// Represents a single segment going out from a point to a nearby one.
// (this is what we want to use for the graph)
pub type OsmNodeId = u64;

pub struct SegmentInfo {
    pub end: OsmNodeId,
    pub maxspeed: Option<u64>,
}

pub struct NodeInfo {
    pub location: Location,
    pub connections : Vec<SegmentInfo>,
}

pub type OsmGraph = HashMap<u64, NodeInfo>;

pub struct Tree {
    pub location: Location,
    pub children: Vec<Tree>,
}

