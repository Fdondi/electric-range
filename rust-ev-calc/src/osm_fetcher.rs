use reqwest::Client;
use serde_json::Value;
use std::error::Error;
use std::collections::HashMap;
use crate::types::{Location, BoundingCoordinates, OsmGraph, NodeInfo, SegmentInfo};

#[derive(PartialEq)]
enum Direction {
    ForwardOnly,
    BackwardOnly,
    BothWays,
}

// Represents a whole road (this is what the API gives us, only used here)
struct WayInfo {
    pub node_ids: Vec<u64>,
    pub maxspeed: Option<u64>,
    pub direction: Direction,
}

fn parse_ways(json: Value) -> Result<OsmGraph, Box<dyn Error>> {
    // Prase both ways and nodes independently in case ways refer to nodes before they are specified. 
    let mut nodes = HashMap::new();
    let mut ways = HashMap::new();
    if let Value::Object(map) = &json {
        let elements = &map["elements"];    
        for element in elements.as_array().unwrap() {
            match element["type"].as_str().unwrap() {
                "node" => {
                    let id = element["id"].as_u64().unwrap();
                    let latitude = element["lat"].as_f64().unwrap();
                    let longitude = element["lon"].as_f64().unwrap();
                    nodes.insert(id, NodeInfo{location: Location{latitude, longitude}, connections: Vec::new()});
                }
                "way" => {
                    let id = element["id"].as_u64().unwrap();
                    let connected_nodes = element["nodes"].as_array().unwrap();
                    let node_ids = connected_nodes.iter().map(|n| n.as_u64().unwrap()).collect::<Vec<u64>>();
                    let direction: Direction = match element.get("tags").and_then(|tags| tags.get("oneway")) {
                        Some(Value::String(ref s)) => if s == "yes" { Direction::ForwardOnly } else if s == "-1" { Direction::BackwardOnly } else { Direction::BothWays },
                        _ => Direction::BothWays,
                    };
                    let maxspeed = element["tags"]["maxspeed"].as_u64();
                    let way_info = WayInfo {
                        node_ids: node_ids,
                        maxspeed: maxspeed,
                        direction: direction,
                    };
                    ways.insert(id, way_info);
                }
                _ => {
                    // Return error
                    println!("Unknown element type: {}", element["type"]);
                    return Err("Unknown element type.".into());
                }
            }
        }
        
        // Now that we have all the nodes and ways, we can incorporate the ways into the nodes
        // Ways go both ways... usually?
        for way in ways.values() {
            for i in 0..way.node_ids.len() - 1 {
                let start = &way.node_ids[i];
                let end = &way.node_ids[i + 1];
                if !(way.direction == Direction::BackwardOnly) {
                    let segment_forward = SegmentInfo {
                        end: end.clone(),
                        maxspeed: way.maxspeed,
                    };
                    nodes.get_mut(&start).unwrap().connections.push(segment_forward);
                }
                if !(way.direction == Direction::ForwardOnly) {
                    let segment_backward = SegmentInfo {
                        end: start.clone(),
                        maxspeed: way.maxspeed,
                    };
                    nodes.get_mut(&end).unwrap().connections.push(segment_backward);
                }
            }
        }

        return Ok(nodes)
    } else {
        println!("Response is not a JSON object.");
        return Err("Response is not a JSON object.".into());
    }
}

/// Fetches OSM road data for Zurich using Overpass API
async fn fetch_zurich_roads(view_area: BoundingCoordinates) -> Result<Value, Box<dyn Error>> {
    let overpass_url = "http://overpass-api.de/api/interpreter";
    // Download the roads passable by car
    // TODO consider downloading separately through and local roads
    // TODO cache instead of re-downloading
    let query = format!(r#"
    [out:json][timeout:60];
    (
       way["highway"~"^(motorway|trunk|primary|secondary|tertiary|unclassified|residential|service)$"]({}, {}, {}, {});
    );
    out body;
    >;
    out skel qt;
    "#, 
    view_area.south_latitude,
    view_area.west_longitude,
    view_area.north_latitude,
    view_area.east_longitude);

    let client = Client::new();
    let response = client
        .post(overpass_url)
        .body(query)
        .send().await?
        .json::<Value>()
        .await?;
        
    Ok(response)
}

pub async fn fetch_zurich_graph(view_area: BoundingCoordinates) -> Result<OsmGraph, Box<dyn Error>> {
    let json = fetch_zurich_roads(view_area).await?;
    parse_ways(json)
}

// test that the fetcher works

#[test]
fn test_fetch_zurich_roads() {
    let response = fetch_zurich_roads().unwrap();
    print!("Response: {:?}", response);
    assert!(response.is_object());
}