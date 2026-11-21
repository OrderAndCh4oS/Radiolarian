//! Trace ink strokes across the visible silica. Both the nib and the strip of
//! paper beside it are sized in specimen units, independent of the export.

use super::{
    LIGHT, Point,
    guide::{Guide, GuideSample},
    ink::Plate,
    lineprint::{Ink, Nib, Segment},
};
use crate::growth::Random;

const STEP: f32 = 0.003;
type Axis = [f32; 2];

#[derive(Clone, Copy)]
enum Layer {
    Along,
    Undertone,
    Across,
}

#[derive(Clone, Copy)]
struct Flow {
    sweep: Axis,
    front: f32,
    turn: f32,
}

impl Layer {
    fn nib_radius(self) -> f32 {
        match self {
            Self::Along => 0.00118,
            Self::Undertone => 0.00081,
            Self::Across => 0.00092,
        }
    }

    fn return_turn(self, ink: Ink, shaft: f32, variation: f32) -> f32 {
        let angle = match (self, ink) {
            (Self::Across, _) => -1.13,
            (_, Ink::Silver) => 1.02,
            (_, Ink::Graphite) => -0.83,
            (_, Ink::Charcoal) => 0.67,
        };
        (angle + 0.20 * (variation - 0.5)) * (1.0 - 0.28 * shaft)
    }

    fn direction(self, s: GuideSample, flow: Flow) -> Axis {
        // On the broad face, marks follow the turn of the entire chamber. As
        // a rib turns edge-on or rises into a spine, its own axis takes over.
        // Double angles let either end of a strut guide a mark at a weld.
        let edge = 1.0 - s.normal.z;
        let along =
            s.coherence * (0.34 + edge * edge + 0.64 * flow.front * (edge + 0.42) + 1.5 * s.shaft);
        let oblique = 0.86 * (1.0 - s.shaft) * (1.0 - 0.32 * flow.front * edge);
        let xx = along * (s.tangent_x * s.tangent_x - s.tangent_y * s.tangent_y)
            + flow.sweep[0] * oblique;
        let xy = along * 2.0 * s.tangent_x * s.tangent_y + flow.sweep[1] * oblique;
        let (ty, tx) = (0.5 * xy.atan2(xx)).sin_cos();
        // Oblique strokes cross the first family in the deeper half-tones.
        // Foreshorten the crossing on a steep wall. On a shaft the principal
        // strokes follow its long axis all the way into the tip.
        let turn = match self {
            Self::Along => 0.0_f32,
            Self::Undertone => -0.04 * (1.0 - s.shaft),
            Self::Across => 0.94 - flow.front * (0.10 + 0.23 * s.shaft),
        };
        let (sin, cos) = (turn + flow.turn).sin_cos();
        let x = cos * tx - sin * ty;
        let y = sin * tx + cos * ty;
        let tilt = (x * s.normal.x + y * s.normal.y) * (0.64 - 0.48 * s.shaft);
        let u = x - tilt * s.normal.x;
        let v = y - tilt * s.normal.y;
        let length = u.hypot(v);
        [u / length, v / length]
    }
}

fn following(axis: Axis, previous: Axis) -> Option<Axis> {
    let agreement = axis[0] * previous[0] + axis[1] * previous[1];
    (agreement.abs() > 0.56).then(|| [axis[0] * agreement.signum(), axis[1] * agreement.signum()])
}

fn chamber_sweep(normal: Point) -> Axis {
    let hand = Point {
        x: 0.68,
        y: -0.70,
        z: -0.22,
    };
    let tilt = 0.78 * hand.dot(normal);
    let x = hand.x - tilt * normal.x;
    let y = hand.y - tilt * normal.y;
    let length2 = x * x + y * y;
    [(x * x - y * y) / length2, 2.0 * x * y / length2]
}

// A stroke reserves a little surrounding paper. The wider clearance in the
// lights and the narrower clearance in the shadows set the local ink density.
struct Spacing {
    left: f32,
    top: f32,
    width: usize,
    height: usize,
    pitch: f32,
    used: Vec<bool>,
}

impl Spacing {
    fn new([left, top, right, bottom]: [f32; 4]) -> Self {
        let pitch = 0.0016_f32.max((right - left).max(bottom - top) / 2400.0);
        let width = ((right - left) / pitch).ceil() as usize + 2;
        let height = ((bottom - top) / pitch).ceil() as usize + 2;
        Self {
            left,
            top,
            width,
            height,
            pitch,
            used: vec![false; width * height],
        }
    }

    fn occupied(&self, p: Point) -> bool {
        let x = ((p.x - self.left) / self.pitch) as usize;
        let y = ((p.y - self.top) / self.pitch) as usize;
        x >= self.width || y >= self.height || self.used[y * self.width + x]
    }

    fn reserve(&mut self, p: Point, radius: f32) {
        let cx = (p.x - self.left) / self.pitch - 0.5;
        let cy = (p.y - self.top) / self.pitch - 0.5;
        let r = radius / self.pitch;
        let right = ((cx + r) as usize).min(self.width - 1);
        let bottom = ((cy + r) as usize).min(self.height - 1);
        for y in (cy - r).ceil() as usize..=bottom {
            let dy = y as f32 - cy;
            for x in (cx - r).ceil() as usize..=right {
                let dx = x as f32 - cx;
                if dx * dx + dy * dy <= r * r {
                    self.used[y * self.width + x] = true;
                }
            }
        }
    }
}

struct Pen {
    ink: Ink,
    length: f32,
    gap: f32,
    opacity: f32,
}

// The second direction enters only as a rib rolls into a sufficiently dark
// tone. Guide normals then keep the packet on the shaded side of the crest.
fn crossing_pressure(wash: f32, surface: GuideSample) -> f32 {
    let value = ((0.75 - wash) / 0.36).clamp(0.0, 1.0);
    let turn = (0.80 - 0.62 * surface.normal.dot(LIGHT)).clamp(0.20, 1.0);
    value * value * (3.0 - 2.0 * value) * turn
}

impl Layer {
    fn pen(self, wash: f32, front: f32, surface: GuideSample, rng: &mut Random) -> Option<Pen> {
        let shaft = surface.shaft;
        let shade = ((0.93 - wash) / 0.46).clamp(0.0, 1.0);
        let length = (0.055 + 0.065 * rng.unit()) * (1.0 + 1.6 * shaft);
        let pressure = 0.88 + 0.12 * rng.unit();
        // Use a lighter hand on the distant wall, leaving the deepest
        // crevices room to carry a full charcoal stroke.
        let t = ((wash - 0.12) / 0.17).clamp(0.0, 1.0);
        let distant = 1.0 - 0.30 * t * t * (3.0 - 2.0 * t);
        // The transition from soft charcoal in the depth to graphite on the
        // face is shared between neighbouring pulls of the two pencils.
        let near = ((front - 0.12) / 0.66).clamp(0.0, 1.0);
        let near = near * near * (3.0 - 2.0 * near);
        match self {
            Self::Along => Some(Pen {
                ink: if rng.unit() < 0.74 * (1.0 - near) {
                    Ink::Charcoal
                } else {
                    Ink::Graphite
                },
                length: length * (1.0 + 0.38 * front),
                gap: 0.0082 - 0.0043 * shade,
                opacity: pressure * (distant + (0.83 - distant) * front),
            }),
            Self::Undertone => (rng.unit() < 0.32 + 0.68 * near).then_some(Pen {
                ink: Ink::Silver,
                length: length * 0.93,
                gap: 0.0055 - 0.0017 * shade,
                opacity: pressure * (0.46 + 0.28 * shade) * (0.76 + 0.24 * near),
            }),
            // The rear wall retains its subdued, fine charcoal crossing.
            Self::Across if front < 0.40 => {
                let load = ((0.72 + 0.14 * front - wash) / 0.33).clamp(0.0, 1.0);
                (rng.unit() < load * (0.75 + 0.22 * front - 0.65 * shaft)).then_some(Pen {
                    ink: Ink::Charcoal,
                    length: length * (0.56 - 0.11 * front),
                    gap: 0.0120 - 0.0035 * load - 0.0030 * front,
                    opacity: pressure * distant,
                })
            }
            Self::Across => {
                let load = crossing_pressure(wash, surface);
                (load > 0.07).then_some(Pen {
                    ink: Ink::Charcoal,
                    length: (0.027 + 0.044 * load) * (1.0 - 0.25 * shaft),
                    gap: 0.0058 - 0.0023 * load,
                    opacity: pressure,
                })
            }
        }
    }
}

struct Drawing<'a> {
    plate: &'a Plate,
    guide: &'a Guide,
    layer: Layer,
    spacing: Spacing,
}

impl Drawing<'_> {
    fn extend(
        &self,
        mut p: Point,
        mut dir: Axis,
        flow: Flow,
        count: usize,
        points: &mut Vec<Point>,
    ) {
        for _ in 0..count {
            let Some(mid) = self
                .guide
                .sample(p.x + STEP * 0.5 * dir[0], p.y + STEP * 0.5 * dir[1])
            else {
                break;
            };
            let Some(tangent) = following(self.layer.direction(mid, flow), dir) else {
                break;
            };
            let x = p.x + STEP * tangent[0];
            let y = p.y + STEP * tangent[1];
            let Some(next) = self.guide.sample(x, y) else {
                break;
            };
            let Some(out) = following(self.layer.direction(next, flow), tangent) else {
                break;
            };
            let tolerance = 0.006 + 0.22 * next.radius.min(mid.radius);
            if (mid.depth - p.z).abs() > tolerance || (next.depth - mid.depth).abs() > tolerance {
                break;
            }
            p = Point {
                x,
                y,
                z: next.depth,
            };
            if self.spacing.occupied(p) {
                break;
            }
            points.push(p);
            dir = out;
        }
    }

    fn stroke(
        &mut self,
        start: Point,
        sample: GuideSample,
        pen: &Pen,
        flow: Flow,
        output: &mut Vec<Segment>,
    ) {
        let dir = self.layer.direction(sample, flow);
        let count = (0.5 * pen.length / STEP).round() as usize;
        let mut points = Vec::with_capacity(2 * count + 1);
        self.extend(start, [-dir[0], -dir[1]], flow, count, &mut points);
        points.reverse();
        points.push(start);
        self.extend(start, dir, flow, count, &mut points);
        if points.len() < 5 {
            return;
        }

        let last = points.len() - 1;
        let mut ink = Random(u64::from(start.x.to_bits()) << 32 | u64::from(start.y.to_bits()));
        // A pencil settles into the stroke, carrying a different load on
        // each pull. Equal-wash neighbours can leave quite different greys.
        let grip = 0.34 + 0.65 * ink.unit();
        let attack = 0.43 + 0.55 * ink.unit();
        let release = 0.30 + 0.58 * ink.unit();
        let strongest = 0.28 + 0.46 * ink.unit();
        let dryness = 0.21 + 0.23 * ink.unit();
        let pressure: Vec<_> = (0..=last / 2 + 1)
            .map(|_| 1.20 - dryness * ink.unit().powi(2) * 1.7)
            .collect();
        let nib = |i: usize| {
            let p = points[i];
            let opacity = if matches!(self.layer, Layer::Across) && flow.front >= 0.40 {
                // Lift the nib at a light-facing lip even if the opposite end
                // of this stroke begins well inside the shadow.
                self.guide.sample(p.x, p.y).map_or(0.0, |s| {
                    let (wash, _, _) = self.plate.tone(p, s.normal, s.shaft);
                    pen.opacity * crossing_pressure(wash, s)
                })
            } else {
                pen.opacity
            };
            let t = i as f32 / last as f32;
            // Set the lead down, bear into the page, then lift. The darkest
            // part of each sweep travels with its own hand movement.
            let (u, end_load) = if t < strongest {
                (t / strongest, attack)
            } else {
                ((1.0 - t) / (1.0 - strongest), release)
            };
            let load = grip * (end_load + (1.0 - end_load) * u * u * (3.0 - 2.0 * u));
            let tip = (t / 0.14).min((1.0 - t) / 0.26).min(1.0);
            let contact = if i.is_multiple_of(2) {
                pressure[i / 2]
            } else {
                0.5 * (pressure[i / 2] + pressure[i / 2 + 1])
            };
            Nib {
                p,
                opacity: (opacity * load * contact * tip * (2.0 - tip)).clamp(0.0, 1.0),
            }
        };
        let depth_tolerance = (0.006 + 0.32 * sample.radius).min(0.028);
        let stroke_id = output.len();
        let mut a = nib(0);
        for i in 1..=last {
            let b = nib(i);
            output.push(Segment {
                a,
                b,
                depth_tolerance,
                ink: pen.ink,
                stroke_id,
                radius: self.layer.nib_radius(),
            });
            a = b;
        }
        for p in points {
            self.spacing.reserve(p, 0.78 * pen.gap);
        }
    }

    fn return_strokes(&mut self, first: usize, end: usize, seed: u64, output: &mut Vec<Segment>) {
        // Place a short return stroke on an actual first impression. Separate
        // spacing for the second pass lets it cross the first and keeps the
        // return strokes apart from one another. Long pulls can carry a small
        // group of crossings, short ones a single accent.
        let mut start = first;
        while start < end {
            let id = output[start].stroke_id;
            let count = output[start..end].partition_point(|s| s.stroke_id == id);
            let offset_in_pass = (start - first) as u64;
            let mut rng = Random(seed ^ offset_in_pass.wrapping_mul(0x9e37_79b9_7f4a_7c15));
            let interval = 4 + (3.0 * rng.unit()) as usize;
            let offset = ((count.min(interval) - 1) as f32 * (0.30 + 0.4 * rng.unit())) as usize;
            for i in (start + offset..start + count).step_by(interval) {
                let segment = output[i];
                let p = segment.b.p;
                if self.spacing.occupied(p) {
                    continue;
                }
                let Some(s) = self.guide.sample(p.x, p.y) else {
                    continue;
                };
                let (wash, normal, front) = self.plate.tone(p, s.normal, s.shaft);
                let darkness = ((0.92 - wash) / 0.54).clamp(0.0, 1.0);
                let acceptance = match segment.ink {
                    Ink::Silver => 0.40 + 0.30 * darkness,
                    Ink::Graphite => 0.28 + 0.38 * darkness,
                    Ink::Charcoal => 0.12 + 0.40 * darkness,
                };
                if rng.unit() >= acceptance {
                    continue;
                }
                let flow = Flow {
                    sweep: chamber_sweep(normal),
                    front,
                    turn: self.layer.return_turn(segment.ink, s.shaft, rng.unit()),
                };
                let pen = Pen {
                    ink: segment.ink,
                    length: (0.030 + 0.055 * darkness) * (0.70 + 0.55 * rng.unit()),
                    gap: 0.0086 - 0.0030 * darkness,
                    opacity: (0.59 + 0.35 * darkness)
                        * (0.82 + 0.18 * front)
                        * (0.88 + 0.12 * rng.unit()),
                };
                self.stroke(p, s, &pen, flow, output);
            }
            start += count;
        }
    }

    fn shadow_packet(
        &mut self,
        p: Point,
        s: GuideSample,
        pen: &Pen,
        flow: Flow,
        output: &mut Vec<Segment>,
    ) {
        let [tx, ty] = self.layer.direction(s, flow);
        // Short parallel pulls share a start along a shadow band. Neighbouring
        // pulls are displaced along the band, staggered slightly at the ends.
        for lane in [0_i32, -1, 1, -2, 2, -3, 3] {
            let offset = lane as f32 * pen.gap;
            let stagger = offset * 0.30;
            let x = p.x - ty * offset + tx * stagger;
            let y = p.y + tx * offset + ty * stagger;
            let Some(side) = self.guide.sample(x, y) else {
                continue;
            };
            let q = Point {
                x,
                y,
                z: side.depth,
            };
            if (q.z - p.z).abs() > 0.007 + 1.5 * offset.abs() || self.spacing.occupied(q) {
                continue;
            }
            let (wash, _, front) = self.plate.tone(q, side.normal, side.shaft);
            if front >= 0.40
                && crossing_pressure(wash, side) > 0.05
                && following(self.layer.direction(side, flow), [tx, ty]).is_some()
            {
                self.stroke(q, side, pen, flow, output);
            }
        }
    }

    fn fill(&mut self, seed: u64, gaps: Option<&Spacing>, output: &mut Vec<Segment>) {
        let [left, top, right, bottom] = self.guide.bounds();
        let patch = match self.layer {
            Layer::Along => 0.0080,
            Layer::Undertone => 0.0063,
            Layer::Across => 0.014,
        };
        let columns = (left / (2.0 * patch)).floor() as i32..(right / (2.0 * patch)).ceil() as i32;
        let rows = (top / (2.0 * patch)).floor() as i32..(bottom / (2.0 * patch)).ceil() as i32;
        // The first, widely spaced seeds have room for a full sweep. The other
        // three interleaved subgrids fill the strips of paper between them.
        for (dr, dc) in [(0, 0), (1, 1), (1, 0), (0, 1)] {
            for r in rows.clone() {
                let row = 2 * r + dr;
                for c in columns.clone() {
                    let col = 2 * c + dc;
                    let cell = ((row as u32 as u64) << 32) | col as u32 as u64;
                    let mut rng = Random(seed ^ cell);
                    let x = patch * (col as f32 + rng.unit());
                    let y = patch * (row as f32 + rng.unit());
                    let Some(s) = self.guide.sample(x, y) else {
                        continue;
                    };
                    let p = Point { x, y, z: s.depth };
                    if self.spacing.occupied(p) || gaps.is_some_and(|primary| primary.occupied(p)) {
                        continue;
                    }
                    let (wash, normal, front) = self.plate.tone(p, s.normal, s.shaft);
                    if let Some(pen) = self.layer.pen(wash, front, s, &mut rng) {
                        let flow = Flow {
                            sweep: chamber_sweep(normal),
                            front,
                            turn: 0.13 * (rng.unit() - 0.5),
                        };
                        if matches!(self.layer, Layer::Across) && front >= 0.40 {
                            self.shadow_packet(p, s, &pen, flow, output);
                        } else {
                            self.stroke(p, s, &pen, flow, output);
                        }
                    }
                }
            }
        }
    }
}

// Lay short, overlapping pulls in the surface's own direction. Staggered
// seeds give the rib a fairly even density as its direction turns. Looking
// either side of a seed also reaches the edges of a narrow or occluded wall.
fn pencil_rows(guide: &Guide, plate: &Plate, ink: Ink, seed: u64, output: &mut Vec<Segment>) {
    let [left, top, right, bottom] = guide.bounds();
    let (pitch, radius) = match ink {
        Ink::Silver => (0.00365_f32, 0.00110),
        Ink::Graphite => (0.00420, 0.00097),
        Ink::Charcoal => unreachable!(),
    };
    let first_col = (left / pitch).floor() as i32 - 1;
    let last_col = (right / pitch).ceil() as i32;
    for row in (top / pitch).floor() as i32 - 1..=(bottom / pitch).ceil() as i32 {
        for col in first_col..=last_col {
            let cell = ((row as u32 as u64) << 32) | col as u32 as u64;
            let mut hand = Random(seed ^ cell);
            let x = (col as f32 + 0.5 * (row & 1) as f32 + 0.82 * hand.unit()) * pitch;
            let y = (row as f32 + 0.82 * hand.unit()) * pitch;
            if !guide.near_wall(Point { x, y, z: 0.0 }) {
                continue;
            }
            // The same square may contain the near wall, a spine, and the far
            // wall behind it. Give each depth an independent pencil pull.
            let mut depths = [f32::INFINITY; 6];
            let mut count = 0;
            for (dx, dy) in [
                (0.0, 0.0),
                (-0.009, 0.0),
                (0.009, 0.0),
                (0.0, -0.009),
                (0.0, 0.009),
                (-0.014, -0.014),
                (0.014, 0.014),
                (-0.014, 0.014),
                (0.014, -0.014),
            ] {
                let Some(z) = guide.depth(x + dx, y + dy) else {
                    continue;
                };
                if depths[..count]
                    .iter()
                    .any(|&depth| (depth - z).abs() < 0.046)
                {
                    continue;
                }
                let Some(s) = guide.sample(x + dx, y + dy) else {
                    continue;
                };
                depths[count] = z;
                count += 1;
                let (_, normal, front) = plate.tone(
                    Point {
                        x: x + dx,
                        y: y + dy,
                        z,
                    },
                    s.normal,
                    s.shaft,
                );
                let skew = if ink == Ink::Silver { -0.14 } else { 0.88 };
                let turn = skew * (1.0 - 0.94 * s.shaft) + 0.16 * (hand.unit() - 0.5);
                let flow = Flow {
                    sweep: chamber_sweep(normal),
                    front,
                    turn,
                };
                let axis = Layer::Along.direction(s, flow);
                let slope = |s: GuideSample, d: Axis| {
                    (-(s.normal.x * d[0] + s.normal.y * d[1]) / s.normal.z.max(0.22))
                        .clamp(-1.4, 1.4)
                };
                let center = Point {
                    x,
                    y,
                    z: s.depth
                        + ((dx * s.normal.x + dy * s.normal.y) / s.normal.z.max(0.22))
                            .clamp(-0.020, 0.020),
                };
                let step = (0.032 + 0.022 * hand.unit()) / 4.0;
                let mut positions = [center; 5];
                for sign in [-1_i32, 1] {
                    let mut p = center;
                    let mut direction = axis;
                    let mut dz = slope(s, direction);
                    for n in 1..=2 {
                        let advance = step * sign as f32;
                        p.x += direction[0] * advance;
                        p.y += direction[1] * advance;
                        p.z += dz * advance;
                        if let Some(next) = guide.sample(p.x, p.y)
                            && (next.depth - p.z).abs() < 0.033
                        {
                            p.z = next.depth;
                            if let Some(tangent) =
                                following(Layer::Along.direction(next, flow), direction)
                            {
                                direction = tangent;
                            }
                            dz = slope(next, direction);
                        }
                        positions[(2 + sign * n) as usize] = p;
                    }
                }
                let load = 0.16 + 0.22 * hand.unit();
                let pressure = [
                    0.13,
                    0.75 + 0.2 * hand.unit(),
                    1.0,
                    0.67 + 0.25 * hand.unit(),
                    0.10,
                ];
                let id = output.len();
                for i in 1..positions.len() {
                    let nib = |j: usize| Nib {
                        p: positions[j],
                        opacity: load * pressure[j],
                    };
                    output.push(Segment {
                        a: nib(i - 1),
                        b: nib(i),
                        ink,
                        radius,
                        stroke_id: id,
                        depth_tolerance: 0.034 + 0.28 * s.radius.min(0.05),
                    });
                }
                if count == depths.len() {
                    break;
                }
            }
        }
    }
}

pub(super) fn trace(guide: &Guide, plate: &Plate, seed: u64) -> Vec<Segment> {
    let mut output = Vec::new();
    pencil_rows(guide, plate, Ink::Silver, seed ^ 0x22bb_d35f, &mut output);
    pencil_rows(guide, plate, Ink::Graphite, seed ^ 0xc47a_7a53, &mut output);
    let mut primary = Spacing::new(guide.bounds());
    for (layer, salt) in [
        (Layer::Along, 0),
        (Layer::Undertone, 0x51b0_d799_c37f_818e),
        (Layer::Across, 0x5ee1_d8b6_d75a_2d1f),
    ] {
        let mut drawing = Drawing {
            plate,
            guide,
            layer,
            spacing: Spacing::new(guide.bounds()),
        };
        let gaps = matches!(layer, Layer::Undertone).then_some(&primary);
        let first = output.len();
        drawing.fill(seed ^ salt, gaps, &mut output);
        let end = output.len();
        if matches!(layer, Layer::Along) {
            // Reserve the actual nib cores for choosing undertone seeds. The
            // original wider spacing remains open to this second impression.
            for segment in &output[first..end] {
                primary.reserve(segment.a.p, segment.radius + 0.0008);
            }
        }
        drawing.spacing = Spacing::new(guide.bounds());
        drawing.return_strokes(first, end, seed ^ salt ^ 0xe772_ac9b, &mut output);
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::growth::{Material, Stroke, Vec3};
    use crate::illustration::{Viewport, envelope::Envelope, eye};

    fn rib() -> GuideSample {
        GuideSample {
            normal: Point {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            coherence: 1.0,
            tangent_x: 1.0,
            tangent_y: 0.0,
            depth: 0.0,
            radius: 0.02,
            shaft: 0.0,
        }
    }

    #[test]
    fn the_underdrawing_follows_two_shafts_at_distinct_depths() {
        // Two shafts cross in projection. The near one runs across the page
        // and the rear one runs down it.
        let model = |x, y, z| {
            let p = Point { x, y, z };
            Vec3 {
                x: p.dot(eye(Vec3 {
                    x: 1.0,
                    ..Vec3::default()
                })),
                y: p.dot(eye(Vec3 {
                    y: 1.0,
                    ..Vec3::default()
                })),
                z: p.dot(eye(Vec3 {
                    z: 1.0,
                    ..Vec3::default()
                })),
            }
        };
        let horizontal = Stroke {
            a: model(-0.11, 0.0, 0.25),
            b: model(0.11, 0.0, 0.25),
            ra: 0.012,
            rb: 0.012,
            tone: 0.5,
            material: Material::Spine,
            shaft: [1.0; 2],
        };
        let vertical = Stroke {
            a: model(0.0, -0.11, 0.0),
            b: model(0.0, 0.11, 0.0),
            ..horizontal
        };
        let frame = [horizontal, vertical];
        let guide = Guide::new(&frame, &frame);
        let plate = Plate {
            envelope: Envelope::from_frame(&frame, crate::BodyPlan::Bell),
            view: Viewport {
                x: 0.0,
                y: 0.0,
                scale: 1.0,
            },
        };
        for ink in [Ink::Silver, Ink::Graphite] {
            let mut strokes = Vec::new();
            pencil_rows(&guide, &plate, ink, 45, &mut strokes);
            let (mut near, mut far) = (0, 0);
            for segment in strokes {
                let p = segment.a.p;
                if p.x.abs().max(p.y.abs()) > 0.065 {
                    continue;
                }
                let d = segment.b.p.minus(segment.a.p);
                let length = d.x.hypot(d.y);
                if p.z > 0.2 {
                    assert!(d.x.abs() > 0.90 * length);
                    near += 1;
                } else {
                    assert!(p.z < 0.05 && d.y.abs() > 0.90 * length);
                    far += 1;
                }
            }
            assert!(near > 100 && far > 100);
        }
    }

    #[test]
    fn crossing_strokes_ease_off_toward_a_lit_crest() {
        let lit = GuideSample {
            normal: LIGHT,
            ..rib()
        };
        let turned = GuideSample {
            normal: Point {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            ..rib()
        };
        let dark = crossing_pressure(0.40, turned);
        let middle = crossing_pressure(0.63, turned);
        assert!(dark > middle && middle > 0.0);
        assert!(crossing_pressure(0.40, lit) < dark);
        assert_eq!(crossing_pressure(0.80, turned), 0.0);
    }

    #[test]
    fn each_pencil_has_a_return_across_its_own_first_direction() {
        for (layer, ink) in [
            (Layer::Along, Ink::Graphite),
            (Layer::Undertone, Ink::Silver),
            (Layer::Along, Ink::Charcoal),
            (Layer::Across, Ink::Charcoal),
        ] {
            let flow = Flow {
                sweep: [1.0, 0.0],
                front: 0.8,
                turn: 0.0,
            };
            let a = layer.direction(rib(), flow);
            let b = layer.direction(
                rib(),
                Flow {
                    turn: layer.return_turn(ink, 0.0, 0.5),
                    ..flow
                },
            );
            assert!((a[0] * b[1] - a[1] * b[0]).abs() > 0.60);
        }
    }

    #[test]
    fn the_rear_wall_receives_both_pencils_and_a_lighter_underdrawing() {
        let mut graphite = 0;
        let mut silver = 0;
        for seed in 0..128 {
            let pen = Layer::Along
                .pen(0.32, 0.05, rib(), &mut Random(seed))
                .unwrap();
            graphite += usize::from(pen.ink == Ink::Graphite);
            silver += usize::from(
                Layer::Undertone
                    .pen(0.32, 0.05, rib(), &mut Random(seed))
                    .is_some(),
            );
        }
        // Charcoal leads on the far wall, with other pencils still visible
        // at exactly the same wash value and depth.
        assert!((1..64).contains(&graphite));
        assert!((1..64).contains(&silver));
    }

    #[test]
    fn primary_strokes_reach_the_lights_and_pack_closer_in_shade() {
        for seed in 0..64 {
            let lit = Layer::Along
                .pen(0.80, 1.0, rib(), &mut Random(seed))
                .unwrap();
            let shaded = Layer::Along
                .pen(0.48, 1.0, rib(), &mut Random(seed))
                .unwrap();
            assert!(lit.gap > shaded.gap);
            assert!(
                Layer::Across
                    .pen(0.90, 1.0, rib(), &mut Random(seed))
                    .is_none()
            );
        }
        assert!((0..64).any(|seed| {
            Layer::Across
                .pen(0.48, 1.0, rib(), &mut Random(seed))
                .is_some()
        }));
    }
}
