//! Two concentric pore lattices joined by six radial beams. The beams continue
//! beyond the outer shell in the manner of six-spined spumellarians.

use super::{Branch, Growth, Material, ROOT_RELIEF, Random, Station, Stroke, Vec3};
use std::collections::{BTreeMap, VecDeque};

fn scale(v: Vec3, factor: f32) -> Vec3 {
    Vec3::default().plus_scaled(v, factor)
}

fn unit(v: Vec3) -> Vec3 {
    scale(v, v.length().recip())
}

fn dot(a: Vec3, b: Vec3) -> f32 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

struct Cage {
    // Each face of the geodesic triangulation gives one three-way weld.
    directions: Vec<Vec3>,
    neighbors: Vec<Vec<usize>>,
}

impl Cage {
    fn new(frequency: usize, rng: &mut Random) -> Self {
        let phi = 0.5 * (1.0 + 5.0_f32.sqrt());
        let vertices = [
            [-1.0, phi, 0.0],
            [1.0, phi, 0.0],
            [-1.0, -phi, 0.0],
            [1.0, -phi, 0.0],
            [0.0, -1.0, phi],
            [0.0, 1.0, phi],
            [0.0, -1.0, -phi],
            [0.0, 1.0, -phi],
            [phi, 0.0, -1.0],
            [phi, 0.0, 1.0],
            [-phi, 0.0, -1.0],
            [-phi, 0.0, 1.0],
        ]
        .map(|[x, y, z]| unit(Vec3 { x, y, z }));
        let faces = [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        let mut sites = Vec::new();
        let mut site_ids = BTreeMap::new();
        let mut triangles = Vec::with_capacity(20 * frequency * frequency);
        for face in faces {
            let mut previous = Vec::new();
            for row in 0..=frequency {
                let mut current = Vec::with_capacity(row + 1);
                for col in 0..=row {
                    let weights = [frequency - row, row - col, col];
                    // Barycentric addresses agree exactly on shared edges.
                    let mut key = std::array::from_fn::<_, 3, _>(|i| {
                        if weights[i] == 0 {
                            (0, 0)
                        } else {
                            (face[i] + 1, weights[i])
                        }
                    });
                    key.sort_unstable();
                    let site = *site_ids.entry(key).or_insert_with(|| {
                        let mut p = Vec3::default();
                        for i in 0..3 {
                            p = p.plus_scaled(vertices[face[i]], weights[i] as f32);
                        }
                        let jitter = Vec3 {
                            x: rng.signed(),
                            y: rng.signed(),
                            z: rng.signed(),
                        };
                        sites.push(unit(p.plus_scaled(jitter, 0.15)));
                        sites.len() - 1
                    });
                    current.push(site);
                }
                for c in 0..row {
                    triangles.push([previous[c], current[c], current[c + 1]]);
                    if c + 1 < row {
                        triangles.push([previous[c], current[c + 1], previous[c + 1]]);
                    }
                }
                previous = current;
            }
        }
        let directions = triangles
            .iter()
            .map(|&[a, b, c]| {
                unit(
                    sites[a]
                        .plus_scaled(sites[b], 1.0)
                        .plus_scaled(sites[c], 1.0),
                )
            })
            .collect();
        let mut neighbors = vec![Vec::with_capacity(3); triangles.len()];
        let mut boundary: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for (i, [a, b, c]) in triangles.into_iter().enumerate() {
            for (u, v) in [(a, b), (b, c), (c, a)] {
                let key = (u.min(v), u.max(v));
                if let Some(other) = boundary.remove(&key) {
                    neighbors[i].push(other);
                    neighbors[other].push(i);
                } else {
                    boundary.insert(key, i);
                }
            }
        }
        debug_assert!(boundary.is_empty());
        Self {
            directions,
            neighbors,
        }
    }

    fn nearest(&self, axis: Vec3) -> usize {
        (0..self.directions.len())
            .max_by(|&a, &b| {
                dot(self.directions[a], axis).total_cmp(&dot(self.directions[b], axis))
            })
            .unwrap()
    }

    fn distances(&self, roots: &[usize]) -> Vec<usize> {
        let mut distance = vec![usize::MAX; self.neighbors.len()];
        let mut queue = VecDeque::from(roots.to_vec());
        for &root in roots {
            distance[root] = 0;
        }
        while let Some(i) = queue.pop_front() {
            for &j in &self.neighbors[i] {
                if distance[j] == usize::MAX {
                    distance[j] = distance[i] + 1;
                    queue.push_back(j);
                }
            }
        }
        distance
    }
}

struct Wall {
    radius: f32,
    rib: f32,
    tone: f32,
    // The three ribs feeding each spine bend outwards to form its foot.
    buds: Vec<Bud>,
}

#[derive(Clone, Copy, Default)]
struct Bud {
    lift: f32,
    // The foot of an accessory spine narrows as the three ribs merge.
    neck_reduction: f32,
}

impl Wall {
    fn station(&self, cage: &Cage, a: usize, b: usize, t: f32) -> Station {
        let lift_a = self.buds[a].lift;
        let lift_b = self.buds[b].lift;
        let v = 1.0 - t;
        let along = t - f32::from(lift_a > 0.0) * t * v * v + f32::from(lift_b > 0.0) * v * t * t;
        let direction = unit(cage.directions[a].lerp(cage.directions[b], along));
        let lift = lift_a * v.powi(3) + lift_b * t.powi(3);
        Station {
            pos: scale(direction, self.radius + lift),
            radius: self.rib * (1.0 + 0.30 * (2.0 * t - 1.0).powi(2)) + 0.33 * lift
                - self.buds[a].neck_reduction * v.powi(3)
                - self.buds[b].neck_reduction * t.powi(3),
            shaft: ROOT_RELIEF
                * (f32::from(lift_a > 0.0) * v.powi(3) + f32::from(lift_b > 0.0) * t.powi(3)),
        }
    }

    fn deposit(
        &self,
        growth: &mut Growth,
        cage: &Cage,
        roots: &[usize],
        start: f32,
        duration: f32,
    ) -> Vec<f32> {
        let distance = cage.distances(roots);
        let step = duration / (*distance.iter().max().unwrap() as f32 + 1.0);
        let joints: Vec<_> = distance
            .iter()
            .map(|&d| start + step * (d as f32 + 1.0))
            .collect();
        for (i, around) in cage.neighbors.iter().enumerate() {
            for &j in around.iter().filter(|&&j| j > i) {
                let (a, b) = if distance[i] <= distance[j] {
                    (i, j)
                } else {
                    (j, i)
                };
                let born = start + step * distance[a] as f32;
                let mut prev = self.station(cage, a, b, 0.0);
                let count = if self.buds[a].lift + self.buds[b].lift > 0.0 {
                    8
                } else {
                    3
                };
                for k in 1..=count {
                    let t = k as f32 / count as f32;
                    let next = self.station(cage, a, b, t);
                    growth.segment(
                        prev,
                        next,
                        self.tone,
                        born + step * (k - 1) as f32 / count as f32,
                        born + step * t,
                    );
                    prev = next;
                }
            }
            // Round the three pore corners by growing short bridges between
            // the bars. Fillets are added after all three bars reach the weld.
            let born = joints[i];
            for (index, &a) in around.iter().enumerate() {
                let b = around[(index + 1) % around.len()];
                let mut prev = self.station(cage, i, a, 0.36);
                for k in 1..=3 {
                    let t = k as f32 / 3.0;
                    let left = (1.0 - t).powi(2);
                    let right = t * t;
                    let reach = 0.36 * (left + right);
                    let u = self.station(cage, i, a, reach);
                    let v = self.station(cage, i, b, reach);
                    let next = Station {
                        pos: u.pos.lerp(v.pos, right / (left + right)),
                        radius: self.rib * 0.77,
                        shaft: u.shaft,
                    };
                    growth.segment(
                        Station {
                            radius: self.rib * 0.77,
                            ..prev
                        },
                        next,
                        self.tone,
                        born + 0.004 * (k - 1) as f32,
                        born + 0.004 * k as f32,
                    );
                    prev = next;
                }
            }
        }
        joints
    }
}

// Opposite pairs occupy the eight octants first, then the twelve gaps halfway
// between adjacent principal axes. Progressive filling preserves both poles.
const ACCESSORY_AXES: [[f32; 3]; 10] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, -1.0],
    [1.0, 0.0, 1.0],
    [1.0, -1.0, 0.0],
    [0.0, 1.0, 1.0],
    [1.0, 0.0, -1.0],
];

pub(super) fn grow(
    seed: u64,
    flare: f32,
    spacing: f32,
    thickness: f32,
    age: f32,
    fine_spines: f32,
) -> Vec<Stroke> {
    let mut rng = Random(seed ^ 0x94d0_c7e5_4520_70f7);
    let frequency = ((4.6 + 0.65 * rng.unit()) / spacing).round() as usize;
    // Changing pore spacing retessellates the walls while the specimen's
    // proportions, orientation and six polar axes stay fixed.
    let mut inner = Cage::new(frequency.div_ceil(2).max(2), &mut Random(seed ^ 0x5268));
    let mut outer = Cage::new(frequency, &mut Random(seed ^ 0x8aec));
    let outer_radius = 0.70 + 0.06 * rng.unit();
    let inner_radius = outer_radius * (0.37 + 0.09 * rng.unit());
    let rib = 0.0145 * thickness;
    let (sp, cp) = (-0.25 + 0.14 * rng.signed()).sin_cos();
    let (sy, cy) = (0.90 + 0.14 * rng.signed()).sin_cos();
    let (sr, cr) = (0.05 + 0.18 * rng.signed()).sin_cos();
    let rotate = |v: Vec3| {
        let (y, z) = (cp * v.y - sp * v.z, sp * v.y + cp * v.z);
        let (x, z) = (cy * v.x + sy * z, -sy * v.x + cy * z);
        Vec3 {
            x: cr * x - sr * y,
            y: sr * x + cr * y,
            z,
        }
    };
    let axes = [
        Vec3 {
            x: 1.0,
            ..Default::default()
        },
        Vec3 {
            y: 1.0,
            ..Default::default()
        },
        Vec3 {
            z: 1.0,
            ..Default::default()
        },
    ];
    // Avoid centring the six spines on the twelve pentagonal pores.
    let bias = |v: Vec3| Vec3 {
        x: 0.80 * v.x - 0.48 * v.y + 0.36 * v.z,
        y: 0.60 * v.x + 0.64 * v.y - 0.48 * v.z,
        z: 0.60 * v.y + 0.80 * v.z,
    };
    let rays: Vec<_> = axes
        .into_iter()
        .flat_map(|v| {
            let v = bias(v);
            [v, scale(v, -1.0)]
        })
        .collect();
    let root_inner: Vec<_> = rays.iter().map(|&r| inner.nearest(r)).collect();
    let root_outer: Vec<_> = rays.iter().map(|&r| outer.nearest(r)).collect();
    for ((&a, &b), &ray) in root_inner.iter().zip(&root_outer).zip(&rays) {
        inner.directions[a] = ray;
        outer.directions[b] = ray;
    }
    for d in inner.directions.iter_mut().chain(&mut outer.directions) {
        *d = rotate(*d);
    }
    let inner_wall = Wall {
        radius: inner_radius,
        rib: 0.68 * rib,
        tone: 0.79,
        buds: vec![Bud::default(); inner.directions.len()],
    };
    let mut outer_wall = Wall {
        radius: outer_radius,
        rib,
        tone: 0.98,
        buds: vec![Bud::default(); outer.directions.len()],
    };
    for &node in &root_outer {
        outer_wall.buds[node].lift = 0.062 * thickness.sqrt();
    }
    let pairs = (ACCESSORY_AXES.len() as f32 * fine_spines).round() as usize;
    let mut accessory = Vec::with_capacity(2 * pairs);
    let mut accessory_rng = Random(seed ^ 0xcc16);
    for &[x, y, z] in &ACCESSORY_AXES[..pairs] {
        let direction = rotate(unit(
            scale(rays[0], x)
                .plus_scaled(rays[2], y)
                .plus_scaled(rays[4], z),
        ));
        let in_octant = x * y * z != 0.0;
        let span = if in_octant { 0.30 } else { 0.22 };
        for sign in [1.0, -1.0] {
            let i = outer.nearest(scale(direction, sign));
            let lift = 0.030 * thickness.sqrt();
            debug_assert_eq!(outer_wall.buds[i].lift, 0.0);
            outer_wall.buds[i] = Bud {
                lift,
                neck_reduction: 0.74 * rib + 0.33 * lift,
            };
            let length = outer_radius * span * (0.92 + 0.16 * accessory_rng.unit()) * flare;
            accessory.push((i, length));
        }
    }
    let mut growth = Growth {
        age,
        material: Material::Lattice,
        strokes: Vec::new(),
    };
    inner_wall.deposit(&mut growth, &inner, &[root_inner[0]], 0.0, 0.235);
    for ((&a, &b), &ray) in root_inner.iter().zip(&root_outer).zip(&rays) {
        let axis = rotate(ray);
        let start = inner_wall.station(&inner, a, inner.neighbors[a][0], 0.0);
        let end = outer_wall.station(&outer, b, outer.neighbors[b][0], 0.0);
        growth.segment(
            start,
            Station {
                pos: scale(axis, outer_radius),
                radius: 0.84 * rib,
                shaft: 0.0,
            },
            0.86,
            0.25,
            0.46,
        );
        growth.segment(
            Station {
                pos: scale(axis, outer_radius),
                radius: 0.84 * rib,
                shaft: 0.0,
            },
            end,
            0.91,
            0.46,
            0.5,
        );
    }
    let joints = outer_wall.deposit(&mut growth, &outer, &root_outer, 0.5, 0.34);
    growth.material = Material::Spine;
    for &i in &root_outer {
        let base = outer_wall.station(&outer, i, outer.neighbors[i][0], 0.0);
        let axis = outer.directions[i];
        let length = outer_radius * (0.72 + 0.21 * rng.unit()) * flare;
        let branch = Branch {
            span: scale(axis, length),
            bend: Vec3::default(),
            radius: base.radius / 1.018,
            tone: 0.96,
            times: [0.56, 0.39],
        };
        branch.extend(&mut growth, base.pos, 0.0);
    }
    for (i, length) in accessory {
        let base = outer_wall.station(&outer, i, outer.neighbors[i][0], 0.0);
        let branch = Branch {
            span: scale(outer.directions[i], length),
            bend: Vec3::default(),
            radius: base.radius,
            tone: outer_wall.tone,
            times: [joints[i] + 0.016, 0.14],
        };
        branch.extend_fine(&mut growth, base.pos);
    }
    growth.strokes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neighboring_faces_form_one_closed_trivalent_cage() {
        for n in [2, 5, 9] {
            let cage = Cage::new(n, &mut Random(31));
            assert_eq!(cage.directions.len(), 20 * n * n);
            assert!(cage.distances(&[0]).iter().all(|&d| d != usize::MAX));
            for (i, neighbors) in cage.neighbors.iter().enumerate() {
                assert_eq!(neighbors.len(), 3);
                assert!(
                    neighbors
                        .iter()
                        .all(|&j| i != j && cage.neighbors[j].contains(&i))
                );
            }
        }
    }

    #[test]
    fn beams_grow_out_of_inner_shell_and_reach_six_spine_roots() {
        let grow_to = |age| grow(31, 1.0, 1.0, 1.0, age, 0.0);
        assert!(grow_to(0.0).is_empty());
        let core = grow_to(0.245);
        assert!(!core.is_empty() && core.iter().all(|s| s.b.length() < 0.36));
        let beams = grow_to(0.5);
        let poles: Vec<_> = beams.iter().filter(|s| s.b.length() > 0.76).collect();
        assert_eq!(poles.len(), 6);
        let shell = grow_to(1.0);
        let spines: Vec<_> = shell
            .iter()
            .filter(|s| s.material == Material::Spine)
            .collect();
        assert_eq!(spines.len(), 60);
        for spine in spines.as_chunks::<10>().0 {
            assert!(
                poles
                    .iter()
                    .any(|p| p.b.plus_scaled(spine[0].a, -1.0).length() < 1e-5)
            );
            assert!(spine[9].b.length() > 1.2);
        }
    }

    fn spine_strokes(skeleton: &[Stroke]) -> Vec<Stroke> {
        skeleton
            .iter()
            .copied()
            .filter(|s| s.material == Material::Spine)
            .collect()
    }

    fn touching(a: Vec3, b: Vec3) -> bool {
        a.plus_scaled(b, -1.0).length() < 2e-5
    }

    fn feeds(skeleton: &[Stroke], axis: Vec3) -> usize {
        skeleton
            .iter()
            .filter(|s| {
                s.material == Material::Lattice && (touching(s.a, axis) || touching(s.b, axis))
            })
            .count()
    }

    #[test]
    fn denser_accessory_patterns_preserve_the_roots_and_shaft_of_each_previous_pair() {
        for (seed, spacing, thickness) in [(2, 0.6, 1.6), (31, 1.0, 1.0), (791, 1.6, 0.6)] {
            let grow_with = |density| grow(seed, 1.0, spacing, thickness, 1.0, density);
            let bare = spine_strokes(&grow_with(0.0));
            let partial = spine_strokes(&grow_with(0.4));
            let skeleton = grow_with(1.0);
            let all = spine_strokes(&skeleton);
            assert!(partial.starts_with(&bare) && all.starts_with(&partial));
            let mut roots = Vec::new();
            let mut previous = None;
            for s in &all[bare.len()..] {
                if previous.is_none_or(|p| !touching(p, s.a)) {
                    assert_eq!(feeds(&skeleton, s.a), 3);
                    assert!(roots.iter().all(|&p| !touching(p, s.a)));
                    roots.push(s.a);
                }
                assert!(s.rb <= s.ra);
                previous = Some(s.b);
            }
            assert_eq!(roots.len(), 20);
            assert_eq!(
                all[bare.len()..].iter().filter(|s| s.rb == 0.0).count(),
                roots.len()
            );
            assert_eq!(
                partial[bare.len()..].iter().filter(|s| s.rb == 0.0).count(),
                8
            );
        }
    }

    #[test]
    fn the_first_fine_spine_grows_after_its_three_ribs_have_reached_the_joint() {
        let mut first = None;
        for tick in 50..=90 {
            let age = tick as f32 * 0.01;
            let bare = spine_strokes(&grow(31, 1.0, 1.0, 1.0, age, 0.0));
            let shell = grow(31, 1.0, 1.0, 1.0, age, 1.0);
            let all = spine_strokes(&shell);
            if let Some(first_fine) = all.get(bare.len()) {
                first = Some((age, first_fine.a));
                break;
            }
        }
        let (age, root) = first.expect("some cortical ribs should have sprouted fine spines");
        assert!(age > 0.5 && age < 0.9);
        let preceding = grow(31, 1.0, 1.0, 1.0, age - 0.01, 1.0);
        assert_eq!(feeds(&preceding, root), 3);
    }
}
