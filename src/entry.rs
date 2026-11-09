//! Fitted helical entries and checked routes through existing pocket guides.
use crate::{
    Point,
    motion::{Curve, MoveKind, Plane, Vec3},
    navigation::Guide,
    pocket::{PocketPlan, Stroke},
    stock::CancelSignal,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    f64::consts::{PI, TAU},
};

// Pocket entry fitting uses one copy of each admitted graph edge, plus any
// leaf endpoint candidates, rather than a doubled traversal. Actual pending
// routes are not length-capped.
/// Tours longer than this are compacted (duplicate guide edges removed) to
/// keep entry search efficient. A performance threshold, not a limit.
pub(crate) const COMPACT_TOURS_ABOVE: usize = 8192;
fn angle_default() -> f64 {
    PI / 60.
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PocketEntry {
    Helix(HelicalEntry),
    StockEdgeOrHelix(HelicalEntry),
    StockEdge(HelicalEntry),
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum EntryPolicy {
    Helix,
    StockEdgeOrHelix,
    StockEdge,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HelicalEntry {
    #[serde(default = "angle_default")]
    pub angle_rad: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radius_mm: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feed_mm_min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clockwise: Option<bool>,
}
impl Default for HelicalEntry {
    fn default() -> Self {
        Self {
            angle_rad: angle_default(),
            radius_mm: None,
            feed_mm_min: None,
            clockwise: None,
        }
    }
}
impl PocketEntry {
    pub(crate) fn helix(&self) -> &HelicalEntry {
        match self {
            Self::Helix(h) | Self::StockEdgeOrHelix(h) | Self::StockEdge(h) => h,
        }
    }
    pub(crate) fn policy(&self) -> EntryPolicy {
        match self {
            Self::Helix(_) => EntryPolicy::Helix,
            Self::StockEdgeOrHelix(_) => EntryPolicy::StockEdgeOrHelix,
            Self::StockEdge(_) => EntryPolicy::StockEdge,
        }
    }
}
impl HelicalEntry {
    pub(crate) fn validate(&self, r: f64, tol: f64, default_feed: f64) -> Result<(), String> {
        if !self.angle_rad.is_finite()
            || !(0. ..PI / 2.).contains(&self.angle_rad)
            || self.angle_rad == 0.
            || self
                .radius_mm
                .is_some_and(|h| !h.is_finite() || h <= 64. * tol || h > r)
            || self.feed_mm_min.unwrap_or(default_feed) <= 0.
            || !self.feed_mm_min.unwrap_or(default_feed).is_finite()
        {
            return Err("invalid helical entry angle, radius or feed".into());
        }
        Ok(())
    }
    pub(crate) fn radii(&self, r: f64, preferred_cap: f64, minimum: f64) -> Vec<f64> {
        if let Some(h) = self.radius_mm {
            return vec![h];
        }
        let preferred = (r / 2.).min(preferred_cap).max(minimum);
        let floor = preferred.min(0.1).max(minimum);
        let mut radii = vec![preferred];
        while radii.len() < 7 && *radii.last().unwrap() / 2. > floor {
            radii.push(*radii.last().unwrap() / 2.);
        }
        if *radii.last().unwrap() > floor {
            radii.push(floor);
        }
        radii
    }
    pub(crate) fn turns(&self, h: f64, depth: f64) -> Result<usize, String> {
        let advance = TAU * h * self.angle_rad.tan() * (1. - 128. * f64::EPSILON);
        let count = (depth / advance).ceil().max(1.);
        if !count.is_finite() || count > 1000. {
            return Err(
                "helical entry exceeds 1000-turn bound; increase radius/angle or reduce depth"
                    .into(),
            );
        }
        Ok(count as usize)
    }
    pub(crate) fn curves(
        &self,
        center: Point,
        h: f64,
        top: f64,
        bottom: f64,
        clockwise: bool,
        opposite: bool,
    ) -> Result<Vec<Curve>, String> {
        let turns = self.turns(h, top - bottom)?;
        let start = [center[0] + if opposite { -h } else { h }, center[1], top];
        let end = [start[0], start[1], bottom];
        let sweep = if self.clockwise.unwrap_or(clockwise) {
            -TAU
        } else {
            TAU
        };
        // One motion per full turn: every controller can represent a
        // <=one-turn helix as a single G2/G3 block.
        let mut curves = Vec::with_capacity(turns + 1);
        let mut z = top;
        for turn in 1..=turns {
            let next = if turn == turns {
                bottom
            } else {
                top - (top - bottom) * turn as f64 / turns as f64
            };
            curves.push(Curve::Helix {
                start: [start[0], start[1], z],
                end: [start[0], start[1], next],
                center: [center[0], center[1], z],
                plane: Plane::Xy,
                sweep_rad: sweep,
            });
            z = next;
        }
        curves.push(Curve::Arc {
            start: end,
            end,
            center: [center[0], center[1], bottom],
            plane: Plane::Xy,
            sweep_rad: sweep,
        });
        Ok(curves)
    }
}
fn dist(a: Point, b: Point) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}
pub(crate) fn point_key(p: Point) -> [u64; 2] {
    p.map(|x| if x == 0. { 0 } else { x.to_bits() })
}
#[derive(Clone)]
pub(crate) struct EntryLocation {
    pub center: Point,
    pub radius: f64,
    pub attachment: usize,
    pub max_turns: usize,
    pub stock_edge: bool,
}
pub(crate) struct EntryFit {
    pub locations: Vec<EntryLocation>,
    pub work: usize,
}
/// Compact only oversized entry-search walks. Connectivity is carried by the
/// admitted edges, not by their traversal order; PocketEntries already builds
/// an undirected graph. Keep every edge and every endpoint as an entry candidate.
/// Adaptive traversal consumers retain their original continuous tours.
pub(crate) fn compact_entry_tours(
    tours: &mut [Vec<Guide>],
    work: &mut usize,
    cancel: &impl CancelSignal,
) -> Result<(), String> {
    for tour in tours {
        if tour.len() <= COMPACT_TOURS_ABOVE {
            continue;
        }
        let mut edges = std::collections::HashSet::new();
        let mut starts = std::collections::HashSet::new();
        let mut compact = Vec::new();
        for &g in tour.iter() {
            if cancel.is_cancelled() {
                return Err("entry graph compaction cancelled".into());
            }
            *work += 1;
            if edges.insert(g.edge) {
                starts.insert(point_key(g.stroke.a));
                compact.push(g);
            }
        }
        let unique = compact.len();
        for i in 0..unique {
            if cancel.is_cancelled() {
                return Err("entry graph compaction cancelled".into());
            }
            *work += 1;
            let g = compact[i];
            if starts.insert(point_key(g.stroke.b)) {
                compact.push(Guide {
                    stroke: g.stroke.reversed(),
                    ..g
                });
            }
        }
        *tour = compact;
    }
    Ok(())
}
/// Every candidate is only a proposal. Original-region disk and complete-line
/// tests admit it; all later navigation follows already admitted guide strokes.
pub(crate) struct EntryFitOptions {
    pub cutter_radius: f64,
    pub minimum_orbit: f64,
    pub preferred_cap: f64,
    pub max_depth: f64,
    pub starting_work: usize,
    pub policy: EntryPolicy,
}
/// The start of `fit_entries`'s refusal when a component has no entry.
pub(crate) const NO_ENTRY: &str = "no source-admitted entry/route";

pub(crate) fn fit_entries(
    pocket: &PocketPlan,
    tours: &[Vec<Guide>],
    h: &HelicalEntry,
    options: EntryFitOptions,
    cancel: &impl CancelSignal,
) -> Result<EntryFit, String> {
    let EntryFitOptions {
        cutter_radius: cutter,
        minimum_orbit: minimum,
        preferred_cap,
        max_depth,
        starting_work,
        policy,
    } = options;
    let mut chosen = vec![None; tours.len()];
    let mut work = starting_work;
    // Each connected tour consists of already source-admitted positive-radius
    // guides/links. Its nodes cannot cross original filled-component boundaries.
    let labels: Vec<_> = tours
        .iter()
        .map(|tour| {
            let first = tour
                .first()
                .ok_or("empty helical entry navigation component")?;
            pocket.source_component(first.stroke.a, &mut work, usize::MAX, cancel)
        })
        .collect::<Result<_, String>>()?;
    let mut seed_labels = HashMap::new();
    let stages: &[bool] = match policy {
        EntryPolicy::Helix => &[false],
        EntryPolicy::StockEdge => &[true],
        EntryPolicy::StockEdgeOrHelix => &[true, false],
    };
    for &stock_edge in stages {
        for radius in h.radii(cutter, preferred_cap, minimum) {
            if chosen.iter().all(Option::is_some) {
                break;
            }
            if radius <= 64. * pocket.operation_tolerance() {
                return Err("automatic helical radius is below numerical resolution".into());
            }
            let max_turns = if stock_edge {
                0
            } else {
                h.turns(radius, max_depth)?
            };
            let candidates = pocket.entry_centers(cutter + radius.max(minimum), cancel)?;
            for (component, tour) in tours.iter().enumerate() {
                if chosen[component].is_some() {
                    continue;
                }
                let Some(first) = tour.first() else {
                    return Err("empty helical entry navigation component".into());
                };
                let mut candidates = candidates.clone();
                candidates
                    .sort_by(|a, b| dist(*a, first.stroke.a).total_cmp(&dist(*b, first.stroke.a)));
                'candidate: for center in candidates {
                    if stock_edge
                        && !pocket.outside_stock_disk(center, cutter + radius.max(minimum))
                    {
                        continue;
                    }
                    if cancel.is_cancelled() {
                        return Err("helical entry fit cancelled".into());
                    }
                    if let Some(label) = labels[component] {
                        let seed = if let Some(label) = seed_labels.get(&point_key(center)) {
                            *label
                        } else {
                            let label =
                                pocket.source_component(center, &mut work, usize::MAX, cancel)?;
                            seed_labels.insert(point_key(center), label);
                            label
                        };
                        if seed.is_some_and(|seed| seed != label) {
                            continue;
                        }
                    }
                    work += pocket.link_work_per_query();
                    if !pocket.link_geometry_fits_radius(
                        center,
                        center,
                        cutter + radius.max(minimum),
                    ) {
                        continue;
                    }
                    let mut nodes: Vec<_> = tour
                        .iter()
                        .enumerate()
                        .map(|(i, g)| (dist(center, g.stroke.a), i))
                        .collect();
                    nodes.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                    for (_, attachment) in nodes {
                        if cancel.is_cancelled() {
                            return Err("helical entry fit cancelled".into());
                        }
                        work += pocket.link_work_per_query();
                        let mut center = center;
                        let node = tour[attachment].stroke.a;
                        if dist(center, node) < 1e-9 && center != node {
                            work += pocket.link_work_per_query();
                            if !pocket.link_geometry_fits_radius(
                                node,
                                node,
                                cutter + radius.max(minimum),
                            ) {
                                continue;
                            }
                            center = node;
                        }
                        if stock_edge
                            && !pocket.outside_stock_disk(center, cutter + radius.max(minimum))
                        {
                            continue;
                        }
                        if pocket.link_geometry_fits_radius(center, node, cutter + minimum) {
                            chosen[component] = Some(EntryLocation {
                                center,
                                radius,
                                attachment,
                                max_turns,
                                stock_edge,
                            });
                            break 'candidate;
                        }
                    }
                }
            }
        }
    }
    let locations=chosen.into_iter().enumerate().map(|(i,c)|c.ok_or_else(||format!("{NO_ENTRY} for component {i}; choose a smaller radius or another entry strategy"))).collect::<Result<_,_>>()?;
    Ok(EntryFit { locations, work })
}
type RouteEdge = (usize, Stroke, f64, Option<usize>);
struct Walk {
    root: Point,
    nodes: Vec<Point>,
    adjacent: Vec<Vec<RouteEdge>>,
}
#[derive(PartialEq)]
struct RouteState {
    cost: f64,
    node: usize,
}
impl Eq for RouteState {}
impl PartialOrd for RouteState {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RouteState {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.node.cmp(&self.node))
    }
}
pub(crate) struct PocketEntries {
    pub locations: Vec<EntryLocation>,
    pub work: usize,
    walks: Vec<Walk>,
    points: HashMap<[u64; 2], (usize, usize)>,
}
impl PocketEntries {
    pub(crate) fn new(tours: Vec<Vec<Guide>>, fit: EntryFit, tol: f64) -> Result<Self, String> {
        let mut walks = vec![];
        let mut points = HashMap::new();
        for (component, mut tour) in tours.into_iter().enumerate() {
            tour.rotate_left(fit.locations[component].attachment);
            let root = tour[0].stroke.a;
            let mut nodes = vec![];
            let mut adjacent: Vec<Vec<RouteEdge>> = vec![];
            let mut used = std::collections::HashSet::new();
            for g in tour {
                if !used.insert(g.edge) {
                    continue;
                }
                let mut ends = [0; 2];
                for (slot, point) in [g.stroke.a, g.stroke.b].into_iter().enumerate() {
                    let key = point_key(point);
                    ends[slot] = points
                        .entry(key)
                        .or_insert_with(|| {
                            let id = nodes.len();
                            nodes.push(point);
                            adjacent.push(vec![]);
                            (component, id)
                        })
                        .1;
                }
                let length = g.stroke.curve(0.).prepare(tol)?.length_mm();
                adjacent[ends[0]].push((ends[1], g.stroke, length, None));
                adjacent[ends[1]].push((ends[0], g.stroke.reversed(), length, None));
            }
            walks.push(Walk {
                root,
                nodes,
                adjacent,
            });
        }
        Ok(Self {
            locations: fit.locations,
            work: fit.work,
            walks,
            points,
        })
    }
    /// Install known analytic subdivisions AFTER fitting the compact guide graph.
    /// The caller creates these pieces from each original curve, in order. This
    /// changes graph resolution only: no new connecting chord is introduced.
    pub(crate) fn subdivide(
        &mut self,
        originals: &[Stroke],
        pieces: &[Stroke],
        ends: &[usize],
        tolerance: f64,
        cancel: &impl CancelSignal,
    ) -> Result<(), String> {
        let same = |a: Stroke, b: Stroke| a.a == b.a && a.b == b.b && a.arc == b.arc;
        let mut start = 0;
        for (group_id, &original) in originals.iter().enumerate() {
            if cancel.is_cancelled() {
                return Err("rest entry subdivision cancelled".into());
            }
            let end = *ends.get(start).ok_or("missing rest guide group")?;
            let group = pieces.get(start..end).ok_or("invalid rest guide group")?;
            let &(component, a) = self
                .points
                .get(&point_key(original.a))
                .ok_or("rest guide start absent")?;
            let &(other, b) = self
                .points
                .get(&point_key(original.b))
                .ok_or("rest guide end absent")?;
            if component != other
                || group.first().is_none_or(|s| s.a != original.a)
                || group.last().is_none_or(|s| s.b != original.b)
            {
                return Err("rest guide endpoints changed".into());
            }
            let walk = &mut self.walks[component];
            self.work += walk.adjacent[a].len() + walk.adjacent[b].len();
            let forward = walk.adjacent[a]
                .iter()
                .position(|(n, s, _, _)| *n == b && same(*s, original))
                .ok_or("rest original forward guide absent")?;
            walk.adjacent[a].remove(forward);
            let reverse = walk.adjacent[b]
                .iter()
                .position(|(n, s, _, _)| *n == a && same(*s, original.reversed()))
                .ok_or("rest original reverse guide absent")?;
            walk.adjacent[b].remove(reverse);
            let mut from = a;
            let mut previous = original.a;
            for &piece in group {
                if cancel.is_cancelled() {
                    return Err("rest entry subdivision cancelled".into());
                }
                self.work += 1;
                if piece.a != previous {
                    return Err("discontinuous rest guide subdivision".into());
                }
                let (c, to) = *self.points.entry(point_key(piece.b)).or_insert_with(|| {
                    let node = walk.nodes.len();
                    walk.nodes.push(piece.b);
                    walk.adjacent.push(vec![]);
                    (component, node)
                });
                if c != component {
                    return Err("rest subdivision crosses components".into());
                }
                let length = piece.curve(0.).prepare(tolerance)?.length_mm();
                walk.adjacent[from].push((to, piece, length, Some(group_id)));
                walk.adjacent[to].push((from, piece.reversed(), length, Some(group_id)));
                previous = piece.b;
                from = to;
            }
            start = end;
        }
        if start != pieces.len() {
            return Err("unattached rest guide pieces".into());
        }
        Ok(())
    }
    pub(crate) fn component(&self, p: Point) -> Result<usize, String> {
        self.points
            .get(&point_key(p))
            .map(|x| x.0)
            .ok_or("entry navigation point absent".into())
    }
    pub(crate) fn root(&self, component: usize) -> Point {
        self.walks[component].root
    }
    pub(crate) fn distances(
        &self,
        from: Point,
        work: &mut usize,
        cancel: &impl CancelSignal,
    ) -> Result<RouteTree, String> {
        self.shortest_tree(from, None, work, cancel)
    }
    fn shortest_tree(
        &self,
        from: Point,
        to: Option<Point>,
        work: &mut usize,
        cancel: &impl CancelSignal,
    ) -> Result<RouteTree, String> {
        let &(component, start) = self
            .points
            .get(&point_key(from))
            .ok_or("entry route start absent")?;
        let stop = if let Some(to) = to {
            let &(c, n) = self
                .points
                .get(&point_key(to))
                .ok_or("entry route end absent")?;
            if c != component {
                return Err("entry route crosses disconnected components".into());
            }
            Some(n)
        } else {
            None
        };
        let w = &self.walks[component];
        let mut distances = vec![f64::INFINITY; w.nodes.len()];
        let mut previous = vec![None; w.nodes.len()];
        let mut queue = std::collections::BinaryHeap::new();
        distances[start] = 0.;
        queue.push(RouteState {
            cost: 0.,
            node: start,
        });
        while let Some(RouteState { cost, node }) = queue.pop() {
            if cancel.is_cancelled() {
                return Err("helical entry route cancelled".into());
            }
            if cost > distances[node] {
                continue;
            }
            if stop == Some(node) {
                break;
            }
            for &(next, stroke, length, group) in &w.adjacent[node] {
                *work = work.checked_add(1).ok_or("entry route work overflow")?;
                let cost = cost + length;
                if cost < distances[next] {
                    distances[next] = cost;
                    previous[next] = Some((node, stroke, group));
                    queue.push(RouteState { cost, node: next });
                }
            }
        }
        Ok(RouteTree {
            component,
            start,
            distances,
            previous,
        })
    }
    pub(crate) fn route(
        &self,
        from: Point,
        to: Point,
        z: f64,
        cancel: &impl CancelSignal,
    ) -> Result<Vec<Curve>, String> {
        let tree = self.shortest_tree(from, Some(to), &mut 0, cancel)?;
        let &(_, mut b) = self
            .points
            .get(&point_key(to))
            .ok_or("entry route end absent")?;
        let mut result = vec![];
        while tree.start != b {
            if cancel.is_cancelled() {
                return Err("helical entry route cancelled".into());
            }
            let (node, stroke, group) =
                tree.previous[b].ok_or("no admitted entry navigation route")?;
            result.push((stroke, group));
            b = node;
        }
        result.reverse();
        // Only recombine subdivisions of the SAME admitted analytic guide.
        // No proximity-based simplification or new corner-cutting chord.
        let mut compact: Vec<(Stroke, Option<usize>)> = vec![];
        for (stroke, group) in result {
            if let Some((previous, id)) = compact.last_mut() {
                let same_curve = match (previous.arc, stroke.arc) {
                    (None, None) => true,
                    (Some((a, u)), Some((b, v))) => a == b && u.signum() == v.signum(),
                    _ => false,
                };
                if group.is_some() && *id == group && previous.b == stroke.a && same_curve {
                    previous.b = stroke.b;
                    if let (Some((_, u)), Some((_, v))) = (&mut previous.arc, stroke.arc) {
                        *u += v;
                    }
                    continue;
                }
            }
            compact.push((stroke, group));
        }
        Ok(compact.into_iter().map(|(s, _)| s.curve(z)).collect())
    }
}
pub(crate) struct RouteTree {
    component: usize,
    start: usize,
    distances: Vec<f64>,
    previous: Vec<Option<(usize, Stroke, Option<usize>)>>,
}
impl RouteTree {
    pub(crate) fn distance(&self, entries: &PocketEntries, to: Point) -> Result<f64, String> {
        let &(component, node) = entries
            .points
            .get(&point_key(to))
            .ok_or("entry distance target absent")?;
        if component != self.component {
            return Err("entry distance crosses disconnected components".into());
        }
        let d = self.distances[node];
        if !d.is_finite() {
            return Err("no admitted entry distance".into());
        }
        Ok(d)
    }
}

pub(crate) type PendingMove = (Curve, MoveKind, Option<f64>);
pub(crate) fn entry_line(
    out: &mut Vec<PendingMove>,
    from: Vec3,
    to: Vec3,
    kind: MoveKind,
    feed: Option<f64>,
) {
    if from != to {
        out.push((
            Curve::Line {
                start: from,
                end: to,
            },
            kind,
            feed,
        ));
    }
}

#[cfg(test)]
mod compact_tests {
    use super::*;
    #[test]
    fn compact_entry_graph_keeps_every_admitted_edge_and_candidate_with_bounded_work() {
        let mut forward = vec![];
        // Large closed loop, plus a leaf branch whose end has no outgoing edge
        // after removing reverse traversals. Neither its edge nor seed may vanish.
        for i in 0..5000 {
            let point = |j: usize| {
                let a = TAU * (j % 5000) as f64 / 5000.;
                [a.cos(), a.sin()]
            };
            forward.push(Guide {
                stroke: Stroke::line(point(i), point(i + 1)),
                radius: 0.,
                coverage: 0.,
                edge: i,
            });
        }
        forward.push(Guide {
            stroke: Stroke::line([1., 0.], [2., 0.]),
            radius: 0.,
            coverage: 0.,
            edge: 5000,
        });
        let mut doubled = forward.clone();
        doubled.extend(forward.iter().rev().map(|&g| Guide {
            stroke: g.stroke.reversed(),
            ..g
        }));
        let mut tours = vec![doubled.clone()];
        let mut work = 0;
        compact_entry_tours(&mut tours, &mut work, &|| false).unwrap();
        assert_eq!(tours[0].len(), 5002);
        assert_eq!(work, doubled.len() + forward.len());
        let original: std::collections::HashSet<_> =
            doubled.iter().map(|g| point_key(g.stroke.a)).collect();
        let compact: std::collections::HashSet<_> =
            tours[0].iter().map(|g| point_key(g.stroke.a)).collect();
        assert_eq!(original, compact);
        for g in forward {
            assert!(
                tours[0].iter().any(|v| v.edge == g.edge
                    && v.stroke.a == g.stroke.a
                    && v.stroke.b == g.stroke.b)
            );
        }
        let attachment = tours[0]
            .iter()
            .position(|g| g.stroke.a == [2., 0.])
            .unwrap();
        let fit = EntryFit {
            locations: vec![EntryLocation {
                center: [2., 0.],
                radius: 1.,
                attachment,
                max_turns: 1,
                stock_edge: false,
            }],
            work,
        };
        let entries = PocketEntries::new(tours, fit, 1e-5).unwrap();
        let route = entries
            .route([2., 0.], [-1., PI.sin()], 3., &|| false)
            .unwrap();
        assert_eq!(route.first().unwrap().start(), [2., 0., 3.]);
        assert_eq!(route.last().unwrap().end(), [-1., PI.sin(), 3.]);
        let length: f64 = route
            .iter()
            .map(|c| c.prepare(1e-5).unwrap().length_mm())
            .sum();
        assert!((length - (1. + PI)).abs() < 1e-5, "{length}");
        let original_points: Vec<_> = doubled
            .iter()
            .map(|g| (g.edge, g.stroke.a, g.stroke.b))
            .collect();
        // Cancellation leaves the tour untouched (no work cap exists).
        {
            let mut tours = vec![doubled.clone()];
            let mut work = 0;
            assert!(compact_entry_tours(&mut tours, &mut work, &|| true).is_err());
            assert_eq!(
                original_points,
                tours[0]
                    .iter()
                    .map(|g| (g.edge, g.stroke.a, g.stroke.b))
                    .collect::<Vec<_>>()
            );
        }
    }
}
