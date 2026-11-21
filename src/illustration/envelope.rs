//! The mature lattice guides the direction and weight of the pencil across
//! adjacent pores, including when only part of the skeleton has grown.

use super::{Material, Point, Stroke, Vec3, eye};
use crate::BodyPlan;

const SECTIONS: usize = 49;

pub(super) struct ShellSample {
    pub(super) normal: Point,
    pub(super) front: f32,
    pub(super) collar_slope: f32,
}

struct Chambers {
    // Each section stores centre X, mean radius, centre Z in specimen space.
    rings: [Point; SECTIONS],
    slopes: [Point; SECTIONS],
    bottom: f32,
    per_y: f32,
    axes: [Point; 3],
    horizontal_view: f32,
}

// Fit a circle to a section, allowing for uneven angular spacing of the
// lattice intersections. Subtracting the means leaves a 2x2 least-squares
// system for the centre. A sparse or collinear section uses the model axis.
fn section(points: &[Point]) -> Point {
    let count = points.len().max(1) as f32;
    let mut mean = Point::default();
    for p in points {
        mean.x += p.x / count;
        mean.y += p.y / count; // squared distance to the model axis
        mean.z += p.z / count;
    }
    let (mut xx, mut xz, mut zz, mut xr, mut zr) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for p in points {
        let d = p.minus(mean);
        xx += d.x * d.x;
        xz += d.x * d.z;
        zz += d.z * d.z;
        xr += d.x * d.y;
        zr += d.z * d.y;
    }
    let det = xx * zz - xz * xz;
    let (cx, cz) = if det > (xx + zz).powi(2) * 0.02 {
        (
            (zz * xr - xz * zr) / (2.0 * det),
            (xx * zr - xz * xr) / (2.0 * det),
        )
    } else {
        (0.0, 0.0)
    };
    let radius = (mean.y + cx * (cx - 2.0 * mean.x) + cz * (cz - 2.0 * mean.z))
        .max(0.0)
        .sqrt();
    Point {
        x: cx,
        y: radius,
        z: cz,
    }
}

fn mix(a: Point, b: Point, t: f32) -> Point {
    Point {
        x: a.x + t * (b.x - a.x),
        y: a.y + t * (b.y - a.y),
        z: a.z + t * (b.z - a.z),
    }
}

fn world_axes() -> [Point; 3] {
    [
        eye(Vec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }),
        eye(Vec3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }),
        eye(Vec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }),
    ]
}

impl Chambers {
    pub(super) fn from_frame(frame: &[Stroke]) -> Self {
        let lattice: Vec<_> = frame
            .iter()
            .filter(|s| s.material == Material::Lattice)
            .collect();
        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        for s in &lattice {
            lo = lo.min(s.a.y).min(s.b.y);
            hi = hi.max(s.a.y).max(s.b.y);
        }
        if lattice.is_empty() {
            (lo, hi) = (0.0, 1.0);
        }
        // A horizontal cut at the very tip of a scalloped rim can hit just
        // one lobe. Keep both end cuts inside the shell for a stable fit.
        let span = (hi - lo).max(0.001);
        let bottom = lo + 0.026 * span;
        let per_y = (SECTIONS - 1) as f32 / (0.948 * span);
        let mut points = Vec::new();
        let rings: [_; SECTIONS] = std::array::from_fn(|i| {
            let y = bottom + i as f32 / per_y;
            points.clear();
            for s in &lattice {
                let dy = s.b.y - s.a.y;
                if dy.abs() < 1e-7 || y < s.a.y.min(s.b.y) || y > s.a.y.max(s.b.y) {
                    continue;
                }
                let t = (y - s.a.y) / dy;
                let x = s.a.x + t * (s.b.x - s.a.x);
                let z = s.a.z + t * (s.b.z - s.a.z);
                points.push(Point {
                    x,
                    y: x * x + z * z,
                    z,
                });
            }
            section(&points)
        });
        // A centred secant spans several cuts. This averages out individual
        // pores while retaining the slope reversal at a chamber neck.
        let slopes = std::array::from_fn(|i| {
            let a = i.saturating_sub(2);
            let b = (i + 2).min(SECTIONS - 1);
            mix(
                Point::default(),
                rings[b].minus(rings[a]),
                per_y / (b - a) as f32,
            )
        });
        let axes = world_axes();
        let horizontal_view = 1.0 / axes[0].z.hypot(axes[2].z);
        Self {
            rings,
            slopes,
            bottom,
            per_y,
            axes,
            horizontal_view,
        }
    }

    pub(super) fn sample(&self, p: Point) -> ShellSample {
        let [x, y, z] = self.axes.map(|axis| p.dot(axis));
        let h = ((y - self.bottom) * self.per_y).clamp(0.0, (SECTIONS - 1) as f32);
        let i = (h as usize).min(SECTIONS - 2);
        let t = h - i as f32;
        let ring = mix(self.rings[i], self.rings[i + 1], t);
        let slope = mix(self.slopes[i], self.slopes[i + 1], t);
        let (dx, dz) = (x - ring.x, z - ring.z);
        let r = (dx * dx + dz * dz).sqrt();
        let (rx, rz) = if r > 1e-6 {
            (dx / r, dz / r)
        } else {
            (0.0, 1.0)
        };
        // The axial component includes the tilt of the local centreline.
        let ny = -(slope.y + slope.x * rx + slope.z * rz);
        let weight = 1.0 / (1.0 + ny * ny).sqrt();
        let normal = eye(Vec3 {
            x: rx * weight,
            y: ny * weight,
            z: rz * weight,
        });
        let front = (rx * self.axes[0].z + rz * self.axes[2].z) * self.horizontal_view;
        ShellSample {
            normal,
            front,
            // The open basal lip has no adjoining chamber below it.
            collar_slope: slope.y * ((h - 1.0) / 3.0).clamp(0.0, 1.0),
        }
    }
}

enum Form {
    Bell(Box<Chambers>),
    Sphere { radius: f32 },
}

pub(super) struct Envelope(Form);

impl Envelope {
    pub(super) fn from_frame(frame: &[Stroke], plan: BodyPlan) -> Self {
        Self(match plan {
            BodyPlan::Bell => Form::Bell(Box::new(Chambers::from_frame(frame))),
            BodyPlan::Sphere => {
                let radius = frame
                    .iter()
                    .filter(|s| s.material == Material::Lattice)
                    .flat_map(|s| [s.a, s.b])
                    .map(|p| p.x * p.x + p.y * p.y + p.z * p.z)
                    .fold(0.001_f32, f32::max)
                    .sqrt();
                Form::Sphere { radius }
            }
        })
    }

    pub(super) fn sample(&self, p: Point) -> ShellSample {
        match &self.0 {
            Form::Bell(chambers) => chambers.sample(p),
            Form::Sphere { radius } => {
                let distance = p.dot(p).sqrt().max(0.0001);
                let normal = mix(Point::default(), p, 1.0 / distance);
                // Let the medullary cage recede behind the lit cortical lattice.
                let outer = ((distance / radius - 0.50) / 0.42).clamp(0.0, 1.0);
                let front = -0.15 + outer * (normal.z + 0.15);
                ShellSample {
                    normal,
                    front,
                    collar_slope: 0.0,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inner_cage_recedes_behind_the_near_hemisphere() {
        let frame = crate::growth::specimen(
            &crate::Settings {
                body_plan: BodyPlan::Sphere,
                ..Default::default()
            },
            1.0,
        );
        let shell = Envelope::from_frame(&frame, BodyPlan::Sphere);
        let front_at = |z| {
            shell
                .sample(Point {
                    z,
                    ..Default::default()
                })
                .front
        };
        assert!(front_at(0.74) > front_at(0.30) + 0.7);
        assert!(front_at(0.30) > front_at(-0.74) + 0.4);
    }

    #[test]
    fn recovers_an_offset_cones_normal_and_ignores_appendages() {
        let mut cage = Vec::new();
        for i in 0..12 {
            let (z, x) = (i as f32 * std::f32::consts::TAU / 12.0).sin_cos();
            let ring = |y: f32| Vec3 {
                x: 0.2 + x * (1.0 - 0.5 * y),
                y,
                z: -0.1 + z * (1.0 - 0.5 * y),
            };
            cage.push(Stroke {
                a: ring(0.0),
                b: ring(1.0),
                ra: 0.01,
                rb: 0.01,
                tone: 1.0,
                material: Material::Lattice,
                shaft: [0.0; 2],
            });
        }
        cage.push(Stroke {
            material: Material::Spine,
            b: Vec3 {
                x: 10.0,
                y: 20.0,
                z: 10.0,
            },
            ..cage[0]
        });
        let envelope = Envelope::from_frame(&cage, BodyPlan::Bell);
        for sign in [1.0, -1.0] {
            let point = eye(Vec3 {
                x: 0.2,
                y: 0.5,
                z: -0.1 + sign * 0.75,
            });
            let sample = envelope.sample(point);
            let expected = eye(Vec3 {
                x: 0.0,
                y: 0.5 / 1.25_f32.sqrt(),
                z: sign / 1.25_f32.sqrt(),
            });
            assert!(sample.normal.dot(expected) > 0.999);
            assert!(sign * sample.front > 0.9);
        }
        assert!(
            Envelope::from_frame(&[], BodyPlan::Bell)
                .sample(Point::default())
                .front
                .is_finite()
        );
    }
}
