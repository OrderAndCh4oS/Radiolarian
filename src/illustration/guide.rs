//! A camera-space line field on the visible silica. Its sampling grid stays
//! registered to the specimen through changes of output size and growth.

use super::{Hit, Point, Stroke, Tube, Viewport, eye};

const DENSITY: f32 = 290.0;
const MAX_SIDE: f32 = 1536.0;
const PAD: f32 = 2.0;
const TILE: usize = 16;

#[derive(Clone, Copy)]
pub(super) struct GuideSample {
    pub(super) tangent_x: f32,
    pub(super) tangent_y: f32,
    pub(super) depth: f32,
    pub(super) radius: f32,
    pub(super) shaft: f32,
    pub(super) normal: Point,
    pub(super) coherence: f32,
}

// Sixteen bytes per guide sample. The double angle represents a line without
// a preferred traversal direction. Opposing endpoints therefore agree at a weld.
#[derive(Clone, Copy, Default)]
struct Cell {
    depth: f32,
    radius: f32,
    axis: [i8; 2],
    normal: [i8; 3],
    shaft: u8,
    weight: u16,
}

fn signed(v: f32) -> i8 {
    (127.0 * v).round().clamp(-127.0, 127.0) as i8
}

impl Cell {
    fn deposit(&mut self, hit: Hit, tube: &Tube, axis: [f32; 2]) {
        let radius = tube.ra + hit.along * (tube.rb - tube.ra);
        let shaft = tube.shaft[0] + hit.along * (tube.shaft[1] - tube.shaft[0]);
        let dz = hit.depth - self.depth;
        let weld = (0.43 * self.radius.min(radius)).max(0.12);
        if self.weight == 0 || dz >= 5.0 * weld {
            *self = Self {
                depth: hit.depth,
                radius,
                axis: axis.map(signed),
                normal: [hit.normal.x, hit.normal.y, hit.normal.z].map(signed),
                shaft: (255.0 * shaft).round() as u8,
                weight: 256,
            };
        } else if dz > -5.0 * weld {
            let old = (self.weight as f32 / 256.0) * (-dz.max(0.0) / weld).exp2();
            let new = (dz.min(0.0) / weld).exp2();
            let t = new / (old + new);
            self.depth = self.depth.max(hit.depth);
            self.radius += t * (radius - self.radius);
            let blend = |a: i8, b: f32| (a as f32 + t * (127.0 * b - a as f32)).round() as i8;
            self.axis = [blend(self.axis[0], axis[0]), blend(self.axis[1], axis[1])];
            self.normal = [
                blend(self.normal[0], hit.normal.x),
                blend(self.normal[1], hit.normal.y),
                blend(self.normal[2], hit.normal.z),
            ];
            self.shaft =
                (self.shaft as f32 + t * (255.0 * shaft - self.shaft as f32)).round() as u8;
            self.weight = (256.0 * (old + new)).round() as u16;
        }
    }
}

pub(super) struct Guide {
    view: Viewport,
    width: usize,
    height: usize,
    cells: Vec<Cell>,
    visible_tiles: Vec<bool>,
}

impl Guide {
    pub(super) fn new(strokes: &[Stroke], frame: &[Stroke]) -> Self {
        let mut left = f32::INFINITY;
        let mut top = f32::INFINITY;
        let mut right = f32::NEG_INFINITY;
        let mut bottom = f32::NEG_INFINITY;
        for s in if frame.is_empty() { strokes } else { frame } {
            for (p, r) in [(eye(s.a), s.ra), (eye(s.b), s.rb)] {
                let r = r.max(0.0);
                left = left.min(p.x - r);
                top = top.min(p.y - r);
                right = right.max(p.x + r);
                bottom = bottom.max(p.y + r);
            }
        }

        let extent = (right - left).max(bottom - top);
        if !extent.is_finite() {
            return Self {
                width: 0,
                height: 0,
                cells: Vec::new(),
                visible_tiles: Vec::new(),
                view: Viewport {
                    x: 0.0,
                    y: 0.0,
                    scale: DENSITY,
                },
            };
        }
        let scale = DENSITY.min((MAX_SIDE - 2.0 * PAD - 2.0) / extent.max(0.001));
        let x = PAD - (left * scale).floor();
        let y = PAD - (top * scale).floor();
        let width = ((right * scale).ceil() + x + PAD) as usize;
        let height = ((bottom * scale).ceil() + y + PAD) as usize;
        let view = Viewport { scale, x, y };
        let mut cells = vec![Cell::default(); width * height];
        for s in strokes {
            let tube = Tube::new(s, view, eye);
            let axis = [
                (tube.ux * tube.ux - tube.uy * tube.uy) * tube.side_slope,
                (2.0 * tube.ux * tube.uy) * tube.side_slope,
            ];
            tube.raster(width, 0, height, |x, y, hit| {
                cells[y * width + x].deposit(hit, &tube, axis);
            });
        }
        // Only a fraction of the bounding rectangle contains silica. Keep
        // a coarse silhouette to skip pencil paths well outside the shell.
        let columns = width.div_ceil(TILE);
        let mut visible_tiles = vec![false; columns * height.div_ceil(TILE)];
        for (y, row) in cells.chunks_exact(width).enumerate() {
            for (x, tile) in row.chunks(TILE).enumerate() {
                if tile.iter().any(|c| c.weight > 0) {
                    visible_tiles[y / TILE * columns + x] = true;
                }
            }
        }
        Self {
            view,
            width,
            height,
            cells,
            visible_tiles,
        }
    }

    // Include adjacent tiles. The printer clips each short section of the
    // stroke against its full-resolution silhouette, including narrow lips.
    pub(super) fn near_wall(&self, p: Point) -> bool {
        let columns = self.width.div_ceil(TILE);
        let rows = self.height.div_ceil(TILE);
        let col = ((self.view.x + p.x * self.view.scale) / TILE as f32).floor() as i32;
        let row = ((self.view.y + p.y * self.view.scale) / TILE as f32).floor() as i32;
        for y in (row - 1).max(0)..(row + 2).min(rows as i32) {
            for x in (col - 1).max(0)..(col + 2).min(columns as i32) {
                if self.visible_tiles[y as usize * columns + x as usize] {
                    return true;
                }
            }
        }
        false
    }

    pub(super) fn bounds(&self) -> [f32; 4] {
        [
            -self.view.x / self.view.scale,
            -self.view.y / self.view.scale,
            (self.width as f32 - self.view.x) / self.view.scale,
            (self.height as f32 - self.view.y) / self.view.scale,
        ]
    }

    /// Fast rejection before interpolating a new pencil starting surface.
    pub(super) fn depth(&self, x: f32, y: f32) -> Option<f32> {
        let x = self.view.x + x * self.view.scale;
        let y = self.view.y + y * self.view.scale;
        if !(0.0..self.width as f32).contains(&x) || !(0.0..self.height as f32).contains(&y) {
            return None;
        }
        let cell = self.cells[y as usize * self.width + x as usize];
        (cell.weight > 0).then_some(cell.depth / self.view.scale)
    }

    pub(super) fn sample(&self, x: f32, y: f32) -> Option<GuideSample> {
        let x = self.view.x + x * self.view.scale;
        let y = self.view.y + y * self.view.scale;
        if !(0.0..self.width as f32).contains(&x) || !(0.0..self.height as f32).contains(&y) {
            return None;
        }
        let center = self.cells[y as usize * self.width + x as usize];
        if center.weight == 0 {
            return None;
        }

        let u = x - 0.5;
        let v = y - 0.5;
        let i = u.floor() as isize;
        let j = v.floor() as isize;
        let fx = u - i as f32;
        let fy = v - j as f32;
        let mut axis = [0.0; 2];
        let mut normal = Point::default();
        let mut depth = 0.0;
        let mut radius = 0.0;
        let mut shaft = 0.0;
        let mut sum = 0.0;
        for (row, wy) in [(j, 1.0 - fy), (j + 1, fy)] {
            for (col, wx) in [(i, 1.0 - fx), (i + 1, fx)] {
                if col < 0 || row < 0 || col as usize >= self.width || row as usize >= self.height {
                    continue;
                }
                let cell = self.cells[row as usize * self.width + col as usize];
                // Avoid mixing the far wall into the near wall across a pore.
                let dz = (cell.depth - center.depth).abs();
                if cell.weight == 0 || dz > 1.0 + 0.65 * cell.radius.min(center.radius) {
                    continue;
                }
                let w = wx * wy;
                axis[0] += w * cell.axis[0] as f32;
                axis[1] += w * cell.axis[1] as f32;
                normal.x += w * cell.normal[0] as f32;
                normal.y += w * cell.normal[1] as f32;
                normal.z += w * cell.normal[2] as f32;
                depth += w * cell.depth;
                radius += w * cell.radius;
                shaft += w * cell.shaft as f32;
                sum += w;
            }
        }
        let inv_length = 1.0 / normal.dot(normal).sqrt().max(0.00001);
        normal.x *= inv_length;
        normal.y *= inv_length;
        normal.z *= inv_length;
        let theta = 0.5 * axis[1].atan2(axis[0]);
        let (tangent_y, tangent_x) = theta.sin_cos();
        Some(GuideSample {
            tangent_x,
            tangent_y,
            depth: depth / (sum * self.view.scale),
            radius: radius / (sum * self.view.scale),
            shaft: shaft / (255.0 * sum),
            normal,
            coherence: (axis[0].hypot(axis[1]) / (127.0 * sum)).min(1.0),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::growth::{Material, Vec3};

    #[test]
    fn a_reversed_rib_keeps_its_direction_in_front_of_a_crossing() {
        // Offsetting along the camera ray leaves the two midpoints in register.
        let z = |x, y, z| eye(Vec3 { x, y, z }).z;
        let near = Vec3 {
            x: z(1.0, 0.0, 0.0) * 0.2,
            y: z(0.0, 1.0, 0.0) * 0.2,
            z: z(0.0, 0.0, 1.0) * 0.2,
        };
        let rib = Stroke {
            a: Vec3 {
                x: near.x - 0.2,
                ..near
            },
            b: Vec3 {
                x: near.x + 0.2,
                ..near
            },
            ra: 0.025,
            rb: 0.025,
            tone: 0.8,
            material: Material::Lattice,
            shaft: [0.0, 1.0],
        };
        let reversed = Stroke {
            a: rib.b,
            b: rib.a,
            shaft: [1.0, 0.0],
            ..rib
        };
        let rear = Stroke {
            a: Vec3 {
                y: -0.2,
                ..Vec3::default()
            },
            b: Vec3 {
                y: 0.2,
                ..Vec3::default()
            },
            shaft: [0.0, 0.0],
            ..rib
        };
        for strokes in [[rear, rib, reversed], [reversed, rib, rear]] {
            let guide = Guide::new(&strokes, &strokes);
            let c = guide.sample(0.0, 0.0).unwrap();
            assert!(c.tangent_x.abs() > 0.96 && c.coherence > 0.85);
            assert!(c.depth > 0.2 && (c.shaft - 0.5).abs() < 0.08);
            let [l, t, ..] = guide.bounds();
            assert!(guide.sample(l, t).is_none());
        }
        assert!(Guide::new(&[], &[]).sample(0.0, 0.0).is_none());
    }
}
