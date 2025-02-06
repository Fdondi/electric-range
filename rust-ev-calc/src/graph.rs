use priority_queue::PriorityQueue;

use std::cmp::{Ord, Ordering, PartialOrd};
use std::fmt;
use std::ops::Add;

use crate::osm_fetcher::fetch_zurich_graph;
use std::error::Error;
use crate::types::{Location, BoundingCoordinates, OsmGraph, OsmNodeId};

struct ReachedNode {
    elevation: f64,
    charge_wh: f64,
}

#[derive(Clone, Copy)]
struct OrderedFloat(f64);

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).expect("NaN encountered")
    }
}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.0.is_nan() || other.0.is_nan() {
            None
        } else {
            Some(self.0.partial_cmp(&other.0)?)
        }
    }
}

impl PartialEq for OrderedFloat {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 && !self.0.is_nan()
    }
}

impl Eq for OrderedFloat {}

impl fmt::Debug for OrderedFloat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "OrderedFloat({})", self.0)
    }
}

// Implement arithmetic operations
impl Add for OrderedFloat {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.0.is_nan() || other.0.is_nan() {
            panic!("Cannot perform addition with NaN values");
        }
        OrderedFloat(self.0 + other.0)
    }
}

impl From<f64> for OrderedFloat {
    fn from(value: f64) -> Self {
        if value.is_nan() {
            panic!("Cannot create OrderedFloat with NaN value");
        }
        OrderedFloat(value)
    }
}

impl Add<f64> for OrderedFloat {
    type Output = Self;

    fn add(self, other: f64) -> Self::Output {
        self.add(OrderedFloat::from(other))
    }
}

#[derive(Eq, PartialEq, Hash)]
struct Edge {
    start_id: OsmNodeId,
    end_id: OsmNodeId,
}

fn distance_meters(a: &Location, b: &Location) -> f64 {
    let ap = haversine_rs::point::Point{ latitude: a.latitude, longitude: a.longitude };
    let bp = haversine_rs::point::Point{ latitude: b.latitude, longitude: b.longitude };
    return haversine_rs::distance(ap, bp, haversine_rs::units::Unit::Meters);
}

fn closest_point(point: &Location, graph: &OsmGraph) -> OsmNodeId {
    let mut closest_node_id = None;
    let mut closest_distance = std::f64::INFINITY;
    for (node_id, node_info) in graph.iter() {
        let distance = distance_meters(&point, &node_info.location);
        if distance < closest_distance {
            closest_distance = distance;
            closest_node_id = Some(node_id);
        }
    }
    *closest_node_id.unwrap()
}

pub async fn reachable_points(origin: &Location, view_area: BoundingCoordinates, initial_charge_wh: f64, consumption_wh_per_m: f64) -> Result<(f64, Vec<[Location; 2]>), Box<dyn Error>> {
    log::info!("Called reachable_points");
    
    // Example fetch function for elevation (stub)
    fn get_elevation(_origin: &Location) -> f64 {
        // In a real implementation, call an API here
        100.0
    }
    
    let graph = fetch_zurich_graph(view_area).await?;

    log::info!("Fetch succesful");

    let root_node_id = closest_point(&origin, &graph);
    log::info!("Closest found");

    /* I can't just expand a node, because I might get to the same node earlier from another one.
    Suppose I can get from root to A in 10, and from root to B in 5, but from A to C in 2 and from B to C in 8.
    Then expanding first B would get us to C in 13, but expanding A would get us to C in 12.
    So we need to add the expansions to a priority queue, and expand the one with the lowest time to reach.
    */

    let mut res = Vec::new();

    let mut reached_nodes = std::collections::HashMap::new();
    let origin_node_id = 0;  // Invalid ID 
    reached_nodes.insert(origin_node_id, ReachedNode{elevation: get_elevation(&origin), charge_wh: initial_charge_wh});
    // Priority queue for next points to visit
    let mut candidate_ids = PriorityQueue::new();
    candidate_ids.push(Edge{start_id: origin_node_id, end_id: root_node_id}, OrderedFloat(0.0));

    let mut closest_end_distance_m = std::f64::INFINITY;

    while let Some((Edge { start_id: source_id, end_id: candidate_id }, reach_time)) = candidate_ids.pop() {
        log::info!("Edge from {} to {} at time {}", source_id, candidate_id, reach_time.0);
        if reached_nodes.contains_key(&candidate_id) {
            log::info!("Already reached node {}", candidate_id);
            continue;
        }
        
        let target_id = candidate_id;
        let source_node_location = if source_id != origin_node_id { &graph.get(&source_id).unwrap().location } else {origin};
        let source_node_info = reached_nodes.get(&source_id).unwrap();
        let target_node = graph.get(&target_id).unwrap();
        let distance_m = distance_meters(&source_node_location, &target_node.location);
        let elevation = get_elevation(&target_node.location);
        let consumption_wh = distance_m * consumption_wh_per_m; // TODO consider elevation
        if source_node_info.charge_wh < consumption_wh {
            let reached_distance_m = distance_meters(&source_node_location, &origin);
            if  reached_distance_m < closest_end_distance_m {
                closest_end_distance_m = reached_distance_m;
                log::info!("New closest end distance: {}", closest_end_distance_m);
            }
            log::info!("Charge {} not enough for consumption {}", source_node_info.charge_wh, consumption_wh);
            continue;
        }
        res.push([source_node_location.clone(), target_node.location.clone()]);
        let remaining_charge_wh = source_node_info.charge_wh - consumption_wh;
        reached_nodes.insert(candidate_id, ReachedNode{elevation, charge_wh: remaining_charge_wh});
        log::info!("Saved node {} with charge {}; now {}", target_id, remaining_charge_wh, res.len());
        for connection in &target_node.connections {
            let connected_node_id = connection.end;
            let new_distance_m = distance_meters(&target_node.location, &graph.get(&connected_node_id).unwrap().location);
            let speed_km_h = connection.maxspeed.unwrap_or(30) as f64;
            let speed_m_h = if speed_km_h < 30.0 {
                log::info!("Speed {} too low", speed_km_h);
                30000.0 // 30 km/h in m/h
            } else {
                speed_km_h * 1000.0
            };
            let time_to_reach = new_distance_m / speed_m_h;
            candidate_ids.push(Edge{start_id: target_id, end_id: connected_node_id}, reach_time + time_to_reach);
            log::info!("Added candidate {} with time {}", connected_node_id, time_to_reach);
        }
    }

    log::info!("Completed reachable_points: {}; all points before {}", res.len(), closest_end_distance_m);

    // Filter all points that are closer to the origin than the closest end point
    let filtered_res: Vec<_> = res.into_iter().filter(|[start, end]| f64::max(distance_meters(&origin, start), distance_meters(&origin, end)) > closest_end_distance_m).collect();

    log::info!("Remaining reachable_points outside the radius of {}: {}",closest_end_distance_m, filtered_res.len());
    Ok((closest_end_distance_m, filtered_res))
}

#[test]
fn test_reachable_points() {
    let origin = Point{lat: 0, lng: 0};
    let initial_charge = 100.0;
    let movement_step_m = 200.0;
    let points = reachable_points(origin, initial_charge, movement_step_m);
    assert_eq!(points.len(), 1);
    assert_eq!(points[0].lat, 0);
    assert_eq!(points[0].lng, 0);
}