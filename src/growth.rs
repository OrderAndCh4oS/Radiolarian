//! An accreting lattice on a three-chambered shell of revolution.
//! Distances are in specimen units. The axial profile and spine placement are
//! artistic choices. The pore spacing and deposition schedule model local growth.

use crate::{BodyPlan, Settings};
use std::collections::BTreeMap;
use std::f32::consts::{PI, TAU};

mod sphere;

pub(crate) fn specimen(settings: &Settings, age: f32) -> Vec<Stroke> {
    let generate = match settings.body_plan {
        BodyPlan::Bell => grow,
        BodyPlan::Sphere => sphere::grow,
    };
    generate(
        settings.seed,
        settings.flare,
        settings.spacing,
        settings.thickness,
        age,
        settings.fine_spines,
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    fn plus_scaled(self, v: Self, t: f32) -> Self {
        Self {
            x: self.x + t * v.x,
            y: self.y + t * v.y,
            z: self.z + t * v.z,
        }
    }

    fn length(self) -> f32 {
        self.x.hypot(self.z).hypot(self.y)
    }

    fn lerp(self, other: Self, t: f32) -> Self {
        Self {
            x: self.x + (other.x - self.x) * t,
            y: self.y + (other.y - self.y) * t,
            z: self.z + (other.z - self.z) * t,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Material {
    Lattice,
    Spine,
}

/// A length of silica skeleton with endpoint radii and reflectance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stroke {
    pub a: Vec3,
    pub b: Vec3,
    pub ra: f32,
    pub rb: f32,
    pub tone: f32,
    pub material: Material,
    /// Fraction of exposed shaft at each end. A mineral foot shares the pencil tone of the shell.
    pub shaft: [f32; 2],
}

pub(crate) struct Random(pub(crate) u64);
impl Random {
    pub(crate) fn unit(&mut self) -> f32 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((z ^ (z >> 31)) >> 40) as f32 / 16_777_216.0
    }
    fn signed(&mut self) -> f32 {
        self.unit() * 2.0 - 1.0
    }
}

#[derive(Clone, Copy)]
struct Surface {
    // The azimuth can run past one full turn while joining the seam.
    angle: f32,
    depth: f32,
}

fn angle_difference(a: f32, b: f32) -> f32 {
    (b - a + PI).rem_euclid(TAU) - PI
}

impl Surface {
    fn toward(self, other: Self, t: f32) -> Self {
        Self {
            angle: self.angle + angle_difference(self.angle, other.angle) * t,
            depth: self.depth + (other.depth - self.depth) * t,
        }
    }
}

struct Shell {
    flare: f32,
    spacing: f32,
    thickness: f32,
    phase: f32,
    lobe: f32,
    shoulder: f32,
    tilt: f32,
}

impl Shell {
    fn radius(&self, t: f32) -> f32 {
        // Monotone, smooth chamber profiles. The small constrictions correspond
        // to the collar below the cephalis and the base of the thorax.
        let knots = [
            (0.0, 0.053),
            (0.055, 0.135),
            (0.135, 0.184),
            (0.200, 0.156),
            (0.29, 0.303),
            (0.425, self.shoulder),
            (0.525, self.shoulder + 0.077),
            (0.605, self.shoulder + 0.047),
            (0.690, self.shoulder + 0.14),
            (0.835, 0.818),
            (0.92, 0.890),
            (1.0, 0.871),
        ];
        let i = knots
            .partition_point(|&(x, _)| x < t)
            .saturating_sub(1)
            .min(knots.len() - 2);
        let (a, r) = knots[i];
        let (b, s) = knots[i + 1];
        let slope = |j: usize| {
            let k = j.saturating_sub(1);
            let l = (j + 1).min(knots.len() - 1);
            if k != j && j != l && (knots[j].1 - knots[k].1) * (knots[l].1 - knots[j].1) <= 0.0 {
                0.0
            } else {
                (knots[l].1 - knots[k].1) / (knots[l].0 - knots[k].0)
            }
        };
        let u = ((t - a) / (b - a)).clamp(0.0, 1.0);
        let v = 1.0 - u;
        let radius = v * v * ((1.0 + 2.0 * u) * r + u * (b - a) * slope(i))
            + u * u * ((3.0 - 2.0 * u) * s - v * (b - a) * slope(i + 1));
        radius * (1.0 + (self.flare - 1.0) * t)
    }

    fn pitch(&self, t: f32) -> f32 {
        (0.058 + 0.096 * t + 0.031 * t * t) * self.spacing
    }

    fn position(&self, site: Surface) -> Vec3 {
        let t = site.depth;
        let a = site.angle;
        let r = self.radius(t)
            * (1.0
                + 0.029 * (3.0 * a + self.phase + 1.3 * t).cos()
                + 0.015 * (7.0 * a + self.lobe - 1.8 * t).cos());
        Vec3 {
            x: r * a.cos() + self.tilt * t * t,
            y: 0.94 - 1.60 * t
                + (0.005 + 0.022 * t * t) * (3.0 * a + self.phase + 2.3 * t).sin()
                + 0.018 * t * t * (7.0 * a + self.lobe).cos(),
            z: r * a.sin() + 0.034 * (4.0 * t + self.phase).sin() * t,
        }
    }

    fn row_step(&self, t: f32) -> f32 {
        // Measure spacing along the sloped shell as well as around each ring.
        let lo = (t - 0.003).max(0.0);
        let hi = (t + 0.003).min(1.0);
        let dr = (self.radius(hi) - self.radius(lo)) / (hi - lo);
        self.pitch(t) * 0.88 / dr.hypot(1.60)
    }

    fn sites(&self, rng: &mut Random) -> (Vec<Surface>, Vec<Vec<usize>>) {
        let mut depths = vec![0.0];
        loop {
            let t = *depths.last().unwrap();
            let step = self.row_step(t);
            if t + step * 1.45 >= 1.0 {
                depths.push(1.0);
                break;
            }
            depths.push(t + step);
        }
        let mut sites = Vec::new();
        let mut rows = Vec::new();
        for (row, &depth) in depths.iter().enumerate() {
            let count = (TAU * self.radius(depth) / self.pitch(depth))
                .round()
                .max(5.0) as usize;
            let shift = (row % 2) as f32 * 0.5 + 0.11 * rng.signed();
            let mut indices = Vec::with_capacity(count);
            for column in 0..count {
                let angle = (column as f32 + shift + 0.135 * rng.signed()) * TAU / count as f32;
                let jitter = if depth > 0.0 && depth < 1.0 {
                    self.row_step(depth) * 0.125
                } else {
                    0.0
                };
                indices.push(sites.len());
                sites.push(Surface {
                    angle,
                    depth: depth + jitter * rng.signed(),
                });
            }
            rows.push(indices);
        }
        (sites, rows)
    }
}

// Each triangle of pore centres contributes one three-way junction. Walking
// both rings by their angular bisectors inserts additional cells as the shell
// flares, and closes the seam without privileging a particular meridian.
fn join_rows(upper: &[usize], lower: &[usize], sites: &[Surface]) -> Vec<[usize; 3]> {
    let mut faces = Vec::with_capacity(upper.len() + lower.len());
    let midpoint = |row: &[usize], i: usize| {
        if i == row.len() {
            return f32::INFINITY;
        }
        let a = sites[row[i]].angle;
        let b = sites[row[(i + 1) % row.len()]].angle;
        a + 0.5 * (b - a).rem_euclid(TAU)
    };
    let (mut i, mut j) = (0, 0);
    while i < upper.len() || j < lower.len() {
        let a = upper[i % upper.len()];
        let b = lower[j % lower.len()];
        if midpoint(upper, i) < midpoint(lower, j) {
            i += 1;
            faces.push([a, upper[i % upper.len()], b]);
        } else {
            j += 1;
            faces.push([a, lower[j % lower.len()], b]);
        }
    }
    faces
}

fn centroid(a: Surface, b: Surface, c: Surface) -> Surface {
    Surface {
        angle: a.angle
            + (angle_difference(a.angle, b.angle) + angle_difference(a.angle, c.angle)) / 3.0,
        depth: (a.depth + b.depth + c.depth) / 3.0,
    }
}

fn lattice(shell: &Shell, rng: &mut Random) -> (Vec<Surface>, Vec<[usize; 2]>, Vec<usize>) {
    let (sites, rows) = shell.sites(rng);
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut boundary = BTreeMap::new();
    for pair in rows.windows(2) {
        for [a, b, c] in join_rows(&pair[0], &pair[1], &sites) {
            let index = nodes.len();
            nodes.push(centroid(sites[a], sites[b], sites[c]));
            for (u, v) in [(a, b), (b, c), (c, a)] {
                let key = (u.min(v), u.max(v));
                if let Some(prev) = boundary.remove(&key) {
                    edges.push([prev, index]);
                } else {
                    boundary.insert(key, index);
                }
            }
        }
    }
    let mut ends = [Vec::new(), Vec::new()];
    for ((a, b), inner) in boundary {
        let site = sites[a].toward(sites[b], 0.5);
        let rim = nodes.len();
        nodes.push(site);
        ends[usize::from(site.depth > 0.5)].push(rim);
        edges.push([inner, rim]);
    }
    for ring in &mut ends {
        ring.sort_by(|&a, &b| {
            nodes[a]
                .angle
                .rem_euclid(TAU)
                .total_cmp(&nodes[b].angle.rem_euclid(TAU))
        });
        for i in 0..ring.len() {
            edges.push([ring[i], ring[(i + 1) % ring.len()]]);
        }
    }
    (nodes, edges, ends[1].clone())
}

/// The surface at the end of a rib. At a budding node all incident ribs
/// share a lifted endpoint and the radius of the emerging shaft.
#[derive(Clone, Copy, Default)]
struct Bud {
    lift: Vec3,
    radius: f32,
}

#[derive(Clone, Copy)]
struct Node {
    site: Surface,
    bud: Bud,
}

#[derive(Clone, Copy)]
struct Station {
    pos: Vec3,
    radius: f32,
    shaft: f32,
}

// Exposed shaft relief where the supporting ribs first meet.
const ROOT_RELIEF: f32 = 0.12;

fn ease(t: f32) -> f32 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// Slow the tangential motion at a lifted junction: the incoming bars turn
// together into the growing axis. With two plain ends this is the shell curve.
fn rib(shell: &Shell, a: Node, b: Node, t: f32) -> Station {
    let v = 1.0 - t;
    let start = f32::from(a.bud.radius > 0.0);
    let end = f32::from(b.bud.radius > 0.0);
    let along = t - start * t * v * v + end * v * t * t;
    let site = a.site.toward(b.site, along);
    let mut pos = shell.position(site);
    if a.site.depth == 1.0 && b.site.depth == 1.0 {
        pos.y += 0.23 * shell.pitch(1.0) * (PI * along).sin();
    }
    let wa = start * v * v * v;
    let wb = end * t * t * t;
    pos = pos.plus_scaled(a.bud.lift, wa).plus_scaled(b.bud.lift, wb);
    let depth = 0.5 * (a.site.depth + b.site.depth);
    let bar = shell.pitch(depth) / shell.spacing * shell.thickness * (0.116 - 0.017 * depth);
    let radius = bar * (1.0 + 0.45 * (2.0 * t - 1.0).powi(2));
    Station {
        pos,
        radius: radius + wa * (a.bud.radius - radius) + wb * (b.bud.radius - radius),
        shaft: ROOT_RELIEF * (wa + wb),
    }
}

struct Growth {
    age: f32,
    material: Material,
    strokes: Vec<Stroke>,
}

impl Growth {
    fn segment(&mut self, a: Station, b: Station, tone: f32, start: f32, end: f32) {
        if self.age < start {
            return;
        }
        let progress = if end <= start {
            1.0
        } else {
            ((self.age - start) / (end - start)).min(1.0)
        };
        if progress <= 0.0 {
            return;
        }
        let maturity = |born: f32| 0.61 + 0.39 * ((self.age - born) / 0.095).clamp(0.0, 1.0);
        self.strokes.push(Stroke {
            a: a.pos,
            b: a.pos.lerp(b.pos, progress),
            ra: a.radius * maturity(start),
            rb: (a.radius + (b.radius - a.radius) * progress)
                * maturity(start + (end - start) * progress),
            tone,
            material: self.material,
            shaft: [a.shaft, a.shaft + progress * (b.shaft - a.shaft)],
        });
    }

    fn strut(&mut self, shell: &Shell, mut a: Node, mut b: Node, tone: f32) {
        if b.site.depth < a.site.depth {
            std::mem::swap(&mut a, &mut b);
        }
        let steps = if a.bud.radius + b.bud.radius > 0.0 {
            10
        } else {
            4
        };
        let time = |t: f32| (a.site.depth + (b.site.depth - a.site.depth) * t) * 0.76;
        let mut prev = rib(shell, a, b, 0.0);
        for k in 1..=steps {
            let t = k as f32 / steps as f32;
            let next = rib(shell, a, b, t);
            self.segment(
                prev,
                next,
                tone,
                time((k - 1) as f32 / steps as f32),
                time(t),
            );
            prev = next;
        }
    }

    // A quadratic fillet follows both pore walls tangentially. Sampling both
    // ribs at equal distance also carries the fillet up a developing bud.
    fn fillet(&mut self, shell: &Shell, center: Node, a: Node, b: Node) {
        let reach = (0.62 * shell.thickness / shell.spacing).clamp(0.22, 0.67);
        let radius = 0.100 * shell.pitch(center.site.depth) / shell.spacing * shell.thickness;
        let point = |u: f32| {
            let left = (1.0 - u).powi(2);
            let right = u * u;
            let distance = reach * (left + right);
            let mix = right / (left + right);
            let start = rib(shell, center, a, distance);
            let end = rib(shell, center, b, distance);
            let width = start.radius + (end.radius - start.radius) * mix;
            Station {
                pos: start.pos.lerp(end.pos, mix),
                radius: radius.min(0.91 * width),
                shaft: start.shaft + mix * (end.shaft - start.shaft),
            }
        };
        let born = 0.76
            * (center.site.depth
                + reach * (a.site.depth.max(b.site.depth) - center.site.depth).max(0.0));
        let mut prev = point(0.0);
        for k in 1..=6 {
            let u = k as f32 / 6.0;
            let next = point(u);
            self.segment(
                prev,
                next,
                0.96,
                born + 0.016 * (k - 1) as f32 / 6.0,
                born + 0.016 * u,
            );
            prev = next;
        }
    }
}

struct Branch {
    span: Vec3,
    bend: Vec3,
    radius: f32,
    tone: f32,
    times: [f32; 2],
}

impl Branch {
    fn offset(&self, t: f32) -> Vec3 {
        Vec3::default()
            .plus_scaled(self.span, t)
            .plus_scaled(self.bend, t * t)
    }

    fn radius(&self, t: f32) -> f32 {
        self.radius * ((1.0 - t).powf(1.35) + 0.018)
    }

    fn takeoff(&self, pitch: f32) -> f32 {
        (0.69 * pitch / self.span.length()).clamp(0.12, 0.50)
    }

    fn extend(&self, growth: &mut Growth, root: Vec3, start: f32) {
        self.draw_shaft(growth, root, start, |t| self.radius(t));
    }

    fn extend_fine(&self, growth: &mut Growth, root: Vec3) {
        // The first part carries the width of the welded ribs. Beyond that
        // foot the walls converge into one much finer, uninterrupted needle.
        self.draw_shaft(growth, root, 0.0, |t| {
            let v = 1.0 - t;
            self.radius * v * (0.21 + 0.79 * v.powi(6))
        });
    }

    fn draw_shaft(&self, growth: &mut Growth, root: Vec3, start: f32, radius: impl Fn(f32) -> f32) {
        let sample = |u: f32| {
            let t = start + u * (1.0 - start);
            Station {
                pos: root.plus_scaled(self.offset(t), 1.0),
                radius: radius(t),
                shaft: ROOT_RELIEF + (1.0 - ROOT_RELIEF) * ease(u / 0.48),
            }
        };
        let mut prev = sample(0.0);
        for i in 1..=10 {
            let u = i as f32 / 10.0;
            let next = sample(u);
            growth.segment(
                prev,
                next,
                self.tone,
                self.times[0] + self.times[1] * (i - 1) as f32 / 10.0,
                self.times[0] + self.times[1] * u,
            );
            prev = next;
        }
    }
}

// Each whorl visits its sites in strides around the circumference. Increasing
// density reveals more of the same specimen, keeping all earlier roots fixed.
fn fine_bell_spines(
    shell: &Shell,
    nodes: &mut [Node],
    neighbors: &[Vec<usize>],
    seed: u64,
    density: f32,
) -> Vec<(usize, Branch)> {
    let mut rng = Random(seed ^ 0x5116_25bf_3bc4_2113);
    let whorls = [(0.28, 0.25), (0.39, 0.33), (0.66, 0.44), (0.79, 0.39)];
    let per_whorl = 14;
    let count = (density * (whorls.len() * per_whorl) as f32).round() as usize;
    let mut branches = Vec::new();
    for i in 0..count {
        let row = i % whorls.len();
        let (depth, reach) = whorls[row];
        let k = (5 * (i / whorls.len())) % per_whorl;
        let angle = shell.phase
            + TAU * (k as f32 + 0.39 * row as f32 + 0.08 * rng.signed()) / per_whorl as f32;
        let depth = depth + 0.013 * rng.signed();
        let shell_radius = shell.radius(depth);
        let distance = |site: Surface| {
            (shell_radius * angle_difference(site.angle, angle)).hypot((site.depth - depth) * 1.60)
        };
        // Leave a complete pore wall between adjacent feet. On coarse
        // lattices this naturally reduces the number of available junctions.
        let root = nodes
            .iter()
            .enumerate()
            .filter(|(i, n)| {
                n.bud.radius == 0.0 && neighbors[*i].iter().all(|&j| nodes[j].bud.radius == 0.0)
            })
            .map(|(i, n)| (i, distance(n.site)))
            .min_by(|a, b| a.1.total_cmp(&b.1));
        let Some((index, _)) = root.filter(|&(_, d)| d < 0.85 * shell.pitch(depth)) else {
            continue;
        };
        let site = nodes[index].site;
        let (s, c) = site.angle.sin_cos();
        // Lean with the chamber profile, measuring its slope along the axis.
        let slope = (shell.radius(site.depth + 0.005) - shell.radius(site.depth - 0.005)) / 0.016;
        let rise = (0.65 * slope + 0.08).clamp(-0.26, 0.85);
        let direction = Vec3 {
            x: c,
            y: rise,
            z: s,
        };
        let direction = Vec3::default().plus_scaled(direction, direction.length().recip());
        let length = reach * (0.91 + 0.25 * rng.unit());
        let pitch = shell.pitch(site.depth);
        let radius = 0.088 * pitch / shell.spacing * shell.thickness;
        nodes[index].bud = Bud {
            lift: Vec3::default().plus_scaled(direction, (0.52 * pitch).min(0.24 * length)),
            radius,
        };
        let sweep = length * (0.05 + 0.06 * rng.signed());
        // Give the three incoming walls and their corner fillets time to join.
        let last_wall = neighbors[index]
            .iter()
            .fold(site.depth, |d, &n| d.max(nodes[n].site.depth));
        branches.push((
            index,
            Branch {
                span: Vec3::default().plus_scaled(direction, length),
                bend: Vec3 {
                    x: -s * sweep,
                    y: 0.035 * length,
                    z: c * sweep,
                },
                radius,
                tone: 0.94,
                times: [0.76 * last_wall + 0.02, 0.17],
            },
        ));
    }
    branches
}

/// Generate a single silica skeleton, with a growth front travelling down its axis.
pub fn grow(
    seed: u64,
    flare: f32,
    spacing: f32,
    thickness: f32,
    age: f32,
    fine_spines: f32,
) -> Vec<Stroke> {
    let mut rng = Random(seed);
    let shell = Shell {
        flare,
        spacing,
        thickness,
        phase: rng.unit() * TAU,
        lobe: rng.unit() * TAU,
        shoulder: 0.444 + rng.unit() * 0.092,
        tilt: rng.signed() * 0.061,
    };
    let (sites, edges, basal) = lattice(&shell, &mut rng);
    let mut nodes: Vec<_> = sites
        .iter()
        .map(|&site| Node {
            site,
            bud: Bud::default(),
        })
        .collect();
    let mut spines = Random(seed ^ 0x30282eb9a529d3a8);
    let mut occupied: Vec<_> = sites
        .iter()
        .map(|p| p.depth == 0.0 || p.depth == 1.0)
        .collect();
    let mut nearest = |angle, depth| {
        let r = shell.radius(depth);
        let distance =
            |p: &Surface| (r * angle_difference(p.angle, angle)).hypot((p.depth - depth) * 1.6);
        let index = sites
            .iter()
            .enumerate()
            .filter(|(i, _)| !occupied[*i])
            .min_by(|(_, a), (_, b)| distance(a).total_cmp(&distance(b)))
            .unwrap()
            .0;
        occupied[index] = true;
        index
    };

    // The small top ring draws inward and upward before its ribs merge into the horn.
    let apex = Vec3 {
        x: 0.0,
        y: 1.011,
        z: 0.0,
    };
    let horn = Branch {
        span: Vec3 {
            x: 0.064,
            y: 0.40 + spines.unit() * 0.16,
            z: -0.028,
        },
        bend: Vec3 {
            x: -0.042 + 0.050 * spines.signed(),
            y: 0.04,
            z: 0.08 * spines.signed(),
        },
        radius: thickness * 0.037,
        tone: 0.985,
        times: [0.034, 0.30],
    };
    for node in nodes.iter_mut().filter(|n| n.site.depth == 0.0) {
        let root = shell.position(node.site);
        let lift = apex.plus_scaled(root, -1.0);
        node.bud = Bud {
            lift: Vec3::default().plus_scaled(lift, 0.42),
            radius: 0.014 * thickness,
        };
    }
    let mut branches = Vec::new();
    // Three shell zones bear differently inclined spines. Attach to the actual
    // junctions so partial growth stages never leave an ornament suspended.
    for (count, depth, reach, rise, width) in [
        (6, 0.15, 0.14, 0.25, 0.014),
        (17, 0.44, 0.35, 0.24, 0.013),
        (11, 0.58, 0.28, 0.07, 0.009),
    ] {
        for k in 0..count {
            let a = shell.phase + (k as f32 + 0.17 * spines.signed()) * TAU / count as f32;
            let index = nearest(a, depth + 0.055 * spines.signed());
            let site = sites[index];
            let len = 0.65 + 0.85 * spines.unit();
            let twist = 0.1 * spines.signed();
            let (s, c) = site.angle.sin_cos();
            branches.push((
                index,
                Branch {
                    span: Vec3 {
                        x: reach * len * c,
                        y: rise * len,
                        z: reach * len * s,
                    },
                    bend: Vec3 {
                        x: -s * twist,
                        y: 0.05 * spines.signed(),
                        z: c * twist,
                    },
                    radius: width * thickness,
                    tone: 0.89 + 0.1 * spines.unit(),
                    times: [site.depth * 0.76, 0.20],
                },
            ));
        }
    }
    for (k, &idx) in basal.iter().enumerate() {
        let site = sites[idx];
        let (s, c) = site.angle.sin_cos();
        let big = k % 5 == 2 || k % 9 == 4;
        let len = if big {
            0.38 + 0.20 * spines.unit()
        } else {
            0.085 + 0.175 * spines.unit()
        };
        let twist = if big { 0.14 } else { 0.065 } * spines.signed();
        let lean = if big { 0.16 } else { 0.085 } * (0.65 + spines.unit());
        let span = Vec3 {
            x: c * lean,
            y: -len,
            z: s * lean,
        };
        let inward = if big { -0.23 } else { -0.093 } * (0.7 + spines.unit() * 0.6);
        let bend = Vec3 {
            x: c * inward - s * twist,
            y: -0.025,
            z: s * inward + c * twist,
        };
        let width = if big { 0.043 } else { 0.024 } * thickness;
        let times = [0.76, 0.16 + 0.08 * spines.unit()];
        branches.push((
            idx,
            Branch {
                span,
                bend,
                radius: width,
                times,
                tone: if big { 0.98 } else { 0.97 },
            },
        ));
    }

    for (index, branch) in &branches {
        let node = &mut nodes[*index];
        let t = branch.takeoff(shell.pitch(node.site.depth));
        node.bud = Bud {
            lift: branch.offset(t),
            radius: branch.radius(t),
        };
    }

    let mut neighbors = vec![Vec::with_capacity(3); nodes.len()];
    for &[a, b] in &edges {
        neighbors[a].push(b);
        neighbors[b].push(a);
    }
    let fine = fine_bell_spines(&shell, &mut nodes, &neighbors, seed, fine_spines);
    let mut growth = Growth {
        age,
        material: Material::Lattice,
        strokes: Vec::with_capacity(edges.len() * 10),
    };
    for [a, b] in edges {
        growth.strut(&shell, nodes[a], nodes[b], 0.94 + 0.06 * rng.unit());
    }
    for (node, around) in nodes.iter().zip(&neighbors) {
        for (i, &a) in around.iter().enumerate() {
            for &b in &around[i + 1..] {
                // Around an open end, the two circumferential bars meet across
                // the open mouth. Only soften corners bounding a pore.
                if (node.site.depth == 0.0 || node.site.depth == 1.0)
                    && nodes[a].site.depth == node.site.depth
                    && nodes[b].site.depth == node.site.depth
                {
                    continue;
                }
                growth.fillet(&shell, *node, nodes[a], nodes[b]);
            }
        }
    }

    growth.material = Material::Spine;
    for node in nodes.iter().filter(|n| n.site.depth == 0.0) {
        let foot = shell.position(node.site).plus_scaled(node.bud.lift, 1.0);
        let station = |u: f32| {
            let curve = foot.lerp(apex, 1.0 - (1.0 - u).powi(2));
            Station {
                pos: Vec3 {
                    y: foot.y + (apex.y - foot.y) * u,
                    ..curve
                },
                radius: node.bud.radius + (horn.radius(0.0) - node.bud.radius) * ease(u),
                shaft: ROOT_RELIEF,
            }
        };
        let mut prev = station(0.0);
        for i in 1..=5 {
            let u = i as f32 / 5.0;
            let next = station(u);
            growth.segment(prev, next, 0.985, 0.034 * (i - 1) as f32 / 5.0, 0.034 * u);
            prev = next;
        }
    }
    horn.extend(&mut growth, apex, 0.0);
    for (index, branch) in branches {
        let site = nodes[index].site;
        branch.extend(
            &mut growth,
            shell.position(site),
            branch.takeoff(shell.pitch(site.depth)),
        );
    }
    for (index, branch) in fine {
        let node = nodes[index];
        let foot = shell.position(node.site).plus_scaled(node.bud.lift, 1.0);
        branch.extend_fine(&mut growth, foot);
    }

    growth.strokes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn increasing_density_keeps_existing_needles_and_they_grow_from_completed_walls() {
        let mature = |density| grow(31, 1.0, 1.0, 1.0, 1.0, density);
        let bare = mature(0.0);
        let primary_count = bare
            .iter()
            .filter(|s| s.material == Material::Spine)
            .count();
        let needles = |strokes: Vec<Stroke>| {
            strokes
                .into_iter()
                .filter(|s| s.material == Material::Spine)
                .skip(primary_count)
                .collect::<Vec<_>>()
        };
        let half = needles(mature(0.5));
        let full = needles(mature(1.0));
        assert!(!half.is_empty() && full.len() > half.len());
        assert!(full.starts_with(&half));
        for age in [0.25, 0.38, 0.55, 0.68] {
            let developing = grow(31, 1.0, 1.0, 1.0, age, 1.0);
            for foot in full.iter().step_by(10) {
                if developing
                    .iter()
                    .any(|s| s.material == Material::Spine && s.a == foot.a)
                {
                    let supports = developing
                        .iter()
                        .filter(|s| {
                            s.material == Material::Lattice
                                && (s.a.plus_scaled(foot.a, -1.0).length() < 1e-5
                                    || s.b.plus_scaled(foot.a, -1.0).length() < 1e-5)
                        })
                        .count();
                    assert_eq!(supports, 3, "incomplete junction at age {age}");
                }
            }
        }
    }

    #[test]
    fn three_supporting_ribs_turn_into_each_shaft() {
        let unit = |v: Vec3| Vec3::default().plus_scaled(v, 1.0 / v.length());
        for (seed, flare, spacing, thickness) in [
            (31, 1.0, 1.0, 1.0),
            (71, 1.6, 0.6, 0.6),
            (12, 0.6, 1.6, 1.6),
        ] {
            let mature = grow(seed, flare, spacing, thickness, 1.0, 1.0);
            let feet: Vec<_> = mature
                .iter()
                .filter(|s| {
                    s.material == Material::Spine
                        && s.shaft[0] == ROOT_RELIEF
                        && s.shaft[1] > s.shaft[0]
                })
                .collect();
            assert!(feet.len() > 34);
            for foot in feet {
                let axis = unit(foot.b.plus_scaled(foot.a, -1.0));
                let supporting = mature
                    .iter()
                    .filter(|s| {
                        if s.shaft[0] > ROOT_RELIEF || s.shaft[1] > ROOT_RELIEF {
                            return false;
                        }
                        let approach = if s.a.plus_scaled(foot.a, -1.0).length() < 1e-5
                            && (s.ra - foot.ra).abs() < 1e-5
                        {
                            s.a.plus_scaled(s.b, -1.0)
                        } else if s.b.plus_scaled(foot.a, -1.0).length() < 1e-5
                            && (s.rb - foot.ra).abs() < 1e-5
                        {
                            s.b.plus_scaled(s.a, -1.0)
                        } else {
                            return false;
                        };
                        let tangent = unit(approach);
                        tangent.x * axis.x + tangent.y * axis.y + tangent.z * axis.z > 0.75
                    })
                    .count();
                assert!(
                    supporting >= 3,
                    "seed {seed} has {supporting} supporting ribs at {:?}",
                    foot.a
                );
            }
        }
    }

    #[test]
    fn widening_surface_inserts_more_pore_sites() {
        let shell = Shell {
            flare: 1.0,
            spacing: 1.0,
            thickness: 1.0,
            phase: 0.0,
            lobe: 0.0,
            shoulder: 0.49,
            tilt: 0.0,
        };
        let (_, rows) = shell.sites(&mut Random(7));
        assert!(rows.last().unwrap().len() > 4 * rows[0].len());
        let coarse = Shell {
            spacing: 1.5,
            ..shell
        };
        assert!(coarse.sites(&mut Random(7)).1.last().unwrap().len() < rows.last().unwrap().len());
    }

    #[test]
    fn early_struts_remain_in_place_and_thicken() {
        let early = grow(12, 1.0, 1.0, 1.0, 0.3, 1.0);
        let mature = grow(12, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(early.len() < mature.len() / 2);
        assert!(early.iter().all(|s| s.b.y > 0.25));
        assert!(mature.iter().any(|s| s.b.y < -0.8));
        let fixed = early
            .iter()
            .filter(|s| {
                mature
                    .iter()
                    .any(|m| m.a == s.a && m.b == s.b && m.ra >= s.ra)
            })
            .count();
        assert!(fixed > early.len() * 8 / 10);
    }
}
