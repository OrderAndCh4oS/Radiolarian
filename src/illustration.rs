//! Orthographic rendering of the silica skeleton in reflected light.
//! The camera and the pen's scale stay fixed throughout a growth sequence.

mod contour;
mod envelope;
mod guide;
mod hatching;
mod ink;
mod lineprint;

use crate::Drawing;
use crate::growth::{Material, Stroke, Vec3};

#[derive(Clone, Copy, Default)]
struct Point {
    x: f32,
    y: f32,
    z: f32,
}

impl Point {
    fn minus(self, p: Self) -> Self {
        Self {
            x: self.x - p.x,
            y: self.y - p.y,
            z: self.z - p.z,
        }
    }

    fn dot(self, p: Self) -> f32 {
        self.x * p.x + self.y * p.y + self.z * p.z
    }
}

// Camera-space coordinates are right, down, and towards the observer.
// Pitch below the equator exposes the underside of the apertural rim.
fn eye(p: Vec3) -> Point {
    let side = 0.951057 * p.x - 0.309017 * p.z;
    let near = 0.309017 * p.x + 0.951057 * p.z;
    let down = -0.970296 * p.y - 0.241922 * near;
    Point {
        x: 0.999391 * side - 0.034899 * down,
        y: 0.034899 * side + 0.999391 * down,
        z: -0.241922 * p.y + 0.970296 * near,
    }
}

const LIGHT: Point = Point {
    x: -0.526782,
    y: -0.622561,
    z: 0.578713,
};
#[derive(Clone, Copy)]
struct Viewport {
    scale: f32,
    x: f32,
    y: f32,
}

impl Viewport {
    fn fit(
        frame: &[Stroke],
        width: usize,
        height: usize,
        transform: impl Fn(Vec3) -> Point,
        fill: f32,
    ) -> Self {
        let mut left = f32::INFINITY;
        let mut top = f32::INFINITY;
        let mut right = f32::NEG_INFINITY;
        let mut bottom = f32::NEG_INFINITY;
        for s in frame {
            for (end, r) in [(s.a, s.ra), (s.b, s.rb)] {
                let p = transform(end);
                left = left.min(p.x - r);
                right = right.max(p.x + r);
                top = top.min(p.y - r);
                bottom = bottom.max(p.y + r);
            }
        }
        let scale = fill
            * (width as f32 / (right - left).max(0.01))
                .min(height as f32 / (bottom - top).max(0.01));
        Self {
            scale,
            x: (width as f32 - scale * (left + right)) * 0.5,
            y: (height as f32 - scale * (top + bottom)) * 0.5,
        }
    }

    fn screen(self, p: Point) -> Point {
        Point {
            x: self.x + self.scale * p.x,
            y: self.y + self.scale * p.y,
            z: self.scale * p.z,
        }
    }
}

/// A tapered silica strut. Each end has a round cap so adjoining struts fuse.
struct Tube {
    a: Point,
    b: Point,
    delta: Point,
    ra: f32,
    rb: f32,
    ux: f32,
    uy: f32,
    depth_slope: f32,
    side_slope: f32,
    projected_length: f32,
    tone: f32,
    shaft: [f32; 2],
    spine: bool,
    // Radius of one lengthwise graphite stroke, in raster samples.
    pencil_radius: f32,
    branch: usize,
    filter: f32,
}

#[derive(Clone, Copy)]
struct Hit {
    depth: f32,
    along: f32,
    coverage: f32,
    normal: Point,
    axial_contact: f32,
}

struct TipFragment {
    index: usize,
    branch: usize,
    depth: f32,
    coverage: f32,
    graphite: u8,
}

impl Tube {
    fn new(stroke: &Stroke, view: Viewport, transform: impl Fn(Vec3) -> Point) -> Self {
        let a = view.screen(transform(stroke.a));
        let b = view.screen(transform(stroke.b));
        let delta = b.minus(a);
        let projected_length = delta.x.hypot(delta.y);
        let inv_length = 1.0 / projected_length.hypot(delta.z).max(0.00001);
        let inv_projected = 1.0 / projected_length.max(0.00001);
        Self {
            a,
            b,
            delta,
            projected_length,
            ra: stroke.ra * view.scale,
            rb: stroke.rb * view.scale,
            ux: delta.x * inv_projected,
            uy: delta.y * inv_projected,
            depth_slope: delta.z * inv_length,
            side_slope: projected_length * inv_length,
            tone: stroke.tone,
            shaft: stroke.shaft,
            spine: stroke.material == Material::Spine,
            pencil_radius: view.scale * 0.00165,
            branch: 0,
            filter: 0.0,
        }
    }

    // Integrate the coverage of a mineral tip narrower than one sample.
    fn coverage(&self, radius: f32, distance: f32) -> f32 {
        if self.filter == 0.0 {
            return 1.0;
        }
        ((radius + self.filter - distance).min(2.0 * radius) / (2.0 * self.filter)).clamp(0.0, 1.0)
    }

    fn side_normal(&self, side: f32, chord: f32, radius: f32) -> Point {
        Point {
            x: (-self.uy * side - self.ux * self.depth_slope * chord) / radius,
            y: (self.ux * side - self.uy * self.depth_slope * chord) / radius,
            z: self.side_slope * chord / radius,
        }
    }

    // A single fixed-width longitudinal mark takes over as the point
    // gets too narrow to accommodate the separate cross-hatching strokes.
    fn axial_contact(&self, radius: f32, side: f32) -> f32 {
        if !self.spine {
            return 0.0;
        }
        let narrow = (1.55 - 0.45 * radius / self.pencil_radius).clamp(0.0, 1.0);
        // When even the pencil is subpixel, its source colour belongs to
        // the covered shaft. The filter separately accounts for its area.
        let distance = if self.pencil_radius < self.filter {
            side.abs().min(radius)
        } else {
            side.abs()
        };
        let edge = (1.4 * self.pencil_radius - distance) / self.pencil_radius;
        narrow * narrow * (3.0 - 2.0 * narrow) * edge.clamp(0.0, 1.0)
    }

    fn hit(&self, x: f32, y: f32) -> Option<Hit> {
        let dx = x - self.a.x;
        let dy = y - self.a.y;
        let along = dx * self.ux + dy * self.uy;
        let side = dy * self.ux - dx * self.uy;
        let t = (along / self.projected_length.max(0.00001)).clamp(0.0, 1.0);
        let physical_radius = self.ra + t * (self.rb - self.ra);
        let radius = physical_radius + self.filter;
        let chord2 = radius * radius - side * side;
        if chord2 < 0.0 || radius <= 0.0 {
            return None;
        }
        let chord = chord2.sqrt();
        // The closest visible point lies nearer the eye than the projected axis.
        let axis = along + chord * self.depth_slope;
        if self.projected_length > 0.00001 && (0.0..=self.projected_length).contains(&axis) {
            let along = axis / self.projected_length;
            return Some(Hit {
                along,
                coverage: self.coverage(radius - self.filter, side.abs()),
                depth: self.a.z + along * self.delta.z + chord * self.side_slope,
                normal: self.side_normal(side, chord, radius),
                axial_contact: self.axial_contact(physical_radius, side),
            });
        }
        let (end, radius, along) =
            if axis < 0.0 || (self.projected_length <= 0.00001 && self.a.z > self.b.z) {
                (self.a, self.ra, 0.0)
            } else {
                (self.b, self.rb, 1.0)
            };
        let dx = x - end.x;
        let dy = y - end.y;
        let radius = radius + self.filter;
        let dz2 = radius * radius - dx * dx - dy * dy;
        if dz2 < 0.0 || radius <= 0.0 {
            return None;
        }
        let dz = dz2.sqrt();
        Some(Hit {
            along,
            depth: end.z + dz,
            coverage: if self.filter > 0.0 {
                self.coverage(radius - self.filter, dx.hypot(dy))
            } else {
                1.0
            },
            // Consecutive stations share these spherical caps geometrically.
            // On a shaft their pencil guides follow the sidewall through the
            // joint instead of picking up a separate end-on highlight.
            normal: if self.spine && self.projected_length > 0.00001 {
                self.side_normal(
                    side,
                    (radius * radius - side * side).max(0.0).sqrt(),
                    radius,
                )
            } else {
                Point {
                    x: dx / radius,
                    y: dy / radius,
                    z: dz / radius,
                }
            },
            axial_contact: self.axial_contact(radius - self.filter, side),
        })
    }

    fn raster(
        &self,
        width: usize,
        first: usize,
        last: usize,
        mut visit: impl FnMut(usize, usize, Hit),
    ) {
        let r = self.ra.max(self.rb) + self.filter;
        let y0 = ((self.a.y.min(self.b.y) - r).floor() as usize).max(first);
        let y1 = ((self.a.y.max(self.b.y) + r).ceil() as usize).min(last);
        for y in y0..y1 {
            let py = y as f32 + 0.5;
            // Clip the axis to this scanline's strip before visiting columns.
            // This keeps the finest long diagonal spines inexpensive to render.
            let (lo, hi) = if self.delta.y.abs() > 0.0001 {
                let u = ((py - r - self.a.y) / self.delta.y).clamp(0.0, 1.0);
                let v = ((py + r - self.a.y) / self.delta.y).clamp(0.0, 1.0);
                let x0 = self.a.x + u * self.delta.x;
                let x1 = self.a.x + v * self.delta.x;
                (x0.min(x1), x0.max(x1))
            } else {
                (self.a.x.min(self.b.x), self.a.x.max(self.b.x))
            };
            let x0 = (lo - r).floor() as usize;
            let x1 = ((hi + r).ceil() as usize).min(width);
            for x in x0..x1 {
                if let Some(hit) = self.hit(x as f32 + 0.5, py) {
                    visit(x, y, hit);
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
struct Surface {
    depth: f32,
    normal: Point,
    radius: f32,
    tone: f32,
    weight: f32,
    shaft: f32,
    axial_contact: f32,
}

impl Surface {
    const EMPTY: Self = Self {
        depth: f32::NEG_INFINITY,
        normal: Point {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        radius: 0.0,
        tone: 0.0,
        weight: 0.0,
        shaft: 0.0,
        axial_contact: 0.0,
    };

    fn new(tube: &Tube, hit: Hit) -> Self {
        Self {
            depth: hit.depth,
            normal: hit.normal,
            radius: tube.ra + hit.along * (tube.rb - tube.ra),
            tone: tube.tone,
            weight: 1.0,
            shaft: tube.shaft[0] + hit.along * (tube.shaft[1] - tube.shaft[0]),
            axial_contact: hit.axial_contact,
        }
    }

    fn axial_graphite(&self, load: u8) -> u8 {
        ((40.0 + 1.2 * f32::from(load)).min(255.0) * self.axial_contact).round() as u8
    }

    // Smooth the normals at silica joins. Both sides of a joint contribute
    // near a contact; distant surfaces in the same camera ray stay occluded.
    fn deposit(&mut self, other: Self) {
        let depth_delta = other.depth - self.depth;
        let blend = (self.radius.min(other.radius) * 0.43).max(0.12);
        if depth_delta >= 5.0 * blend {
            *self = other;
        } else if depth_delta > -5.0 * blend {
            let prior_weight = self.weight * (-depth_delta.max(0.0) / blend).exp2();
            let incoming_weight = (depth_delta.min(0.0) / blend).exp2();
            self.weight = prior_weight + incoming_weight;
            let f = incoming_weight / self.weight;
            self.normal.x += f * (other.normal.x - self.normal.x);
            self.normal.y += f * (other.normal.y - self.normal.y);
            self.normal.z += f * (other.normal.z - self.normal.z);
            self.tone += f * (other.tone - self.tone);
            self.shaft += f * (other.shaft - self.shaft);
            self.axial_contact += f * (other.axial_contact - self.axial_contact);
            self.radius += f * (other.radius - self.radius);
            self.depth = self.depth.max(other.depth);
        }
    }
}

// Supersample rows available to the normal and contour filters on either
// side of a render strip. This is also the maximum normal filter radius.
const FILTER_MARGIN: usize = 12;

// Round the shading across adjacent silica planes. Each spoke walks the
// visible surface one sample at a time, stopping at a hole or a depth jump.
// The symmetric tangent-plane test rejects a second wall glimpsed through a
// pore and accepts the gradual change of depth on an inclined rib.
fn local_surface(samples: &[Surface], index: usize, width: usize, x: usize) -> Surface {
    let mut center = samples[index];
    let reach = (0.55 * center.radius).min(FILTER_MARGIN as f32);
    let limit = reach.ceil() as usize;
    if reach < 1.2
        || index / width < limit
        || index / width + limit >= samples.len() / width
        || x < limit
        || x + limit >= width
    {
        return center;
    }

    let mut normal = center.normal;
    let mut shaft = center.shaft;
    let mut weight = 1.0;
    for (dx, dy) in [
        (-1, 0),
        (1, 0),
        (0, -1),
        (0, 1),
        (-1, -1),
        (1, -1),
        (-1, 1),
        (1, 1),
    ] {
        let stride = if dx != 0 && dy != 0 {
            std::f32::consts::SQRT_2
        } else {
            1.0
        };
        let offset = dx + dy * width as isize;
        let mut previous = center.depth;
        for step in 1..=(reach / stride).floor() as isize {
            let other = samples[(index as isize + step * offset) as usize];
            let radius = center.radius.min(other.radius);
            if other.weight == 0.0 || (other.depth - previous).abs() > 0.35 * radius + 2.0 * stride
            {
                break;
            }
            let delta = Point {
                x: (step * dx) as f32,
                y: (step * dy) as f32,
                z: other.depth - center.depth,
            };
            let departure = center
                .normal
                .dot(delta)
                .abs()
                .max(other.normal.dot(delta).abs());
            let continuity = 1.0 - departure / (0.34 * radius + 0.85);
            if continuity <= 0.0 {
                break;
            }
            let distance = step as f32 * stride / (reach + 0.5);
            let spatial = 1.0 - distance * distance;
            let w = continuity * spatial * spatial;
            normal.x += w * other.normal.x;
            normal.y += w * other.normal.y;
            normal.z += w * other.normal.z;
            shaft += w * other.shaft;
            weight += w;
            previous = other.depth;
        }
    }
    center.normal = Point {
        x: normal.x / weight,
        y: normal.y / weight,
        z: normal.z / weight,
    };
    center.shaft = shaft / weight;
    center
}

/// Render one specimen and its supersampled ink layers for the live shader.
/// `frame` is the mature specimen and keeps the camera still as `strokes` grows.
pub fn draw(
    strokes: &[Stroke],
    frame: &[Stroke],
    width: u32,
    height: u32,
    seed: u64,
    body_plan: crate::BodyPlan,
) -> Drawing {
    let (width, height) = (width as usize, height as usize);
    let mut drawing = Drawing {
        pixels: vec![0; width * height],
        samples: vec![[0; 4]; width * height * 4],
        tips: Vec::new(),
    };
    if width == 0 || height == 0 || strokes.is_empty() {
        return drawing;
    }
    let frame = if frame.is_empty() { strokes } else { frame };
    let (sw, sh) = (2 * width, 2 * height);
    let view = Viewport::fit(frame, sw, sh, eye, 0.856);
    let mut previous_end = None;
    let mut branch = 0;
    let mut tubes: Vec<_> = strokes
        .iter()
        .map(|s| {
            let mut tube = Tube::new(s, view, eye);
            if previous_end != Some(s.a) {
                branch += 1;
            }
            previous_end = Some(s.b);
            tube.branch = branch;
            if tube.spine && tube.ra.max(tube.rb) < 1.5_f32.min(0.0014 * view.scale) {
                tube.filter = 0.7;
            }
            tube
        })
        .collect();
    tubes.sort_by(|a, b| a.filter.total_cmp(&b.filter));
    let plate = ink::Plate {
        view,
        envelope: envelope::Envelope::from_frame(frame, body_plan),
    };
    let rows = 32.min(sh);
    let lines = {
        let guide = guide::Guide::new(strokes, frame);
        let segments = hatching::trace(&guide, &plate, seed);
        lineprint::PrintLines::new(segments, view, sw, sh, rows)
    };
    let pen_width = (view.scale * 0.0035).clamp(0.7, 5.4);
    let opaque = tubes.partition_point(|t| t.filter == 0.0);
    let mut samples = vec![Surface::EMPTY; (rows + 2 * FILTER_MARGIN) * sw];
    let mut pigment = vec![0u8; rows * sw];
    let mut tips = Vec::new();
    for top in (0..sh).step_by(rows) {
        let bottom = (top + rows).min(sh);
        let padded_top = top.saturating_sub(FILTER_MARGIN);
        let padded_bottom = (bottom + FILTER_MARGIN).min(sh);
        let patch = &mut samples[..(padded_bottom - padded_top) * sw];
        patch.fill(Surface::EMPTY);
        for tube in &tubes[..opaque] {
            tube.raster(sw, padded_top, padded_bottom, |x, y, hit| {
                patch[(y - padded_top) * sw + x].deposit(Surface::new(tube, hit));
            });
        }
        for y in top..bottom {
            for x in 0..sw {
                let index = (y - padded_top) * sw + x;
                pigment[(y - top) * sw + x] = if patch[index].weight > 0.0 {
                    let rim = contour::edge(patch, index, sw, pen_width);
                    let surface = local_surface(patch, index, sw, x);
                    plate.pencil_load(surface, x, y, rim)
                } else {
                    0
                };
            }
        }
        lines.draw(
            top,
            bottom,
            padded_top,
            patch,
            &mut pigment,
            &mut drawing.samples[top * sw..bottom * sw],
        );
        for tube in &tubes[opaque..] {
            tube.raster(sw, top, bottom, |x, y, hit| {
                let wall = patch[(y - padded_top) * sw + x];
                if hit.depth > wall.depth && hit.coverage > 0.0 {
                    let index = (y - top) * sw + x;
                    let surface = Surface::new(tube, hit);
                    let load = plate.pencil_load(surface, x, y, contour::Rim::default());
                    tips.push(TipFragment {
                        index,
                        branch: tube.branch,
                        depth: hit.depth,
                        coverage: hit.coverage,
                        graphite: surface.axial_graphite(load),
                    });
                }
            });
        }
        // Two successive pieces of a taper share a cap. Their silhouette
        // covers this sample once. Distinct crossing spines still layer.
        tips.sort_unstable_by(|a, b| {
            (a.index, a.branch)
                .cmp(&(b.index, b.branch))
                .then_with(|| b.coverage.total_cmp(&a.coverage))
                .then_with(|| b.depth.total_cmp(&a.depth))
        });
        tips.dedup_by_key(|tip| (tip.index, tip.branch));
        tips.sort_unstable_by(|a, b| {
            a.index
                .cmp(&b.index)
                .then_with(|| a.depth.total_cmp(&b.depth))
        });
        for tip in tips.drain(..) {
            let inks = [0, tip.graphite, 0];
            let source = f32::from(lineprint::print_inks(255, inks));
            let old = f32::from(pigment[tip.index]);
            pigment[tip.index] = (old + tip.coverage * (source - old)).round() as u8;
            // The shader reads a u32 sample address, four ink bytes and f32 coverage.
            let address = (top * sw + tip.index) as u32;
            drawing.tips.extend_from_slice(&address.to_le_bytes());
            drawing.tips.extend_from_slice(&[255, 0, tip.graphite, 0]);
            drawing.tips.extend_from_slice(&tip.coverage.to_le_bytes());
        }
        for y in (top..bottom).step_by(2) {
            let first = (y - top) * sw;
            let dest = (y / 2) * width;
            for x in 0..width {
                let i = first + x * 2;
                let sum = u16::from(pigment[i])
                    + u16::from(pigment[i + 1])
                    + u16::from(pigment[i + sw])
                    + u16::from(pigment[i + sw + 1]);
                drawing.pixels[dest + x] = ((sum + 2) / 4) as u8;
            }
        }
    }
    drawing
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shader_ink_layers_reproduce_the_reference_drawing_with_filtered_tips() {
        let settings = crate::Settings {
            width: 121,
            height: 153,
            ..Default::default()
        };
        let drawing = crate::render(&settings).unwrap();
        let sw = 2 * settings.width as usize;
        assert_eq!(drawing.samples.len(), sw * 2 * settings.height as usize);
        let mut samples: Vec<u8> = drawing
            .samples
            .iter()
            .map(|s| {
                assert!(s[0] == 0 || s[0] == 255);
                if s[0] == 0 {
                    assert_eq!(s[1..], [0; 3]);
                }
                lineprint::print_inks(s[0], [s[1], s[2], s[3]])
            })
            .collect();
        let (tips, remainder) = drawing.tips.as_chunks::<12>();
        assert!(!tips.is_empty() && remainder.is_empty());
        let mut last_index = 0;
        let mut graphite_tip = false;
        for event in tips {
            let index = u32::from_le_bytes(event[..4].try_into().unwrap()) as usize;
            assert!(index >= last_index);
            last_index = index;
            let source: [u8; 4] = event[4..8].try_into().unwrap();
            assert_eq!(source[0], 255);
            graphite_tip |= source[2] > 0;
            let grey = f32::from(lineprint::print_inks(
                source[0],
                [source[1], source[2], source[3]],
            ));
            let coverage = f32::from_le_bytes(event[8..].try_into().unwrap());
            let old = f32::from(samples[index]);
            samples[index] = (old + coverage * (grey - old)).round() as u8;
        }
        assert!(graphite_tip);
        let mut resolved = Vec::new();
        for rows in samples.chunks_exact(2 * sw) {
            for x in (0..sw).step_by(2) {
                let sum: u16 = [x, x + 1, x + sw, x + sw + 1]
                    .into_iter()
                    .map(|i| u16::from(rows[i]))
                    .sum();
                resolved.push(((sum + 2) / 4) as u8);
            }
        }
        assert_eq!(resolved, drawing.pixels);
    }

    #[test]
    fn shaft_relief_varies_along_shaft_and_remains_constant_on_each_cap() {
        let stroke = Stroke {
            a: Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            b: Vec3 {
                x: 10.0,
                y: 0.0,
                z: 0.0,
            },
            ra: 1.0,
            rb: 1.0,
            tone: 1.0,
            material: Material::Spine,
            shaft: [0.2, 0.8],
        };
        let mut tube = Tube::new(
            &stroke,
            Viewport {
                x: 0.0,
                y: 0.0,
                scale: 1.0,
            },
            |p| Point {
                x: p.x,
                y: p.y,
                z: p.z,
            },
        );
        for filter in [0.0, 0.7] {
            tube.filter = filter;
            for (x, expected) in [(-0.5, 0.2), (5.0, 0.5), (10.5, 0.8)] {
                let surface = Surface::new(&tube, tube.hit(x, 0.0).unwrap());
                assert!((surface.shaft - expected).abs() < 1e-6);
            }
            // End caps overlap at every construction station. At the same
            // offset across a shaft the pencil sees the same sidewall normal.
            let interior = tube.hit(5.0, 0.25).unwrap().normal;
            for x in [-0.5, 10.5] {
                let cap = tube.hit(x, 0.25).unwrap().normal;
                assert!(interior.dot(cap) > 0.9999);
            }
        }
    }

    #[test]
    fn touching_wall_and_shaft_share_relief_and_a_distant_wall_stays_occluded() {
        let wall = Surface {
            depth: 1.0,
            normal: Point {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            weight: 1.0,
            radius: 4.0,
            ..Surface::EMPTY
        };
        let shaft = Surface {
            shaft: 1.0,
            normal: Point {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            ..wall
        };
        let mut joined = wall;
        joined.deposit(shaft);
        assert_eq!(
            (joined.shaft, joined.normal.x, joined.normal.z),
            (0.5, 0.5, 0.5)
        );

        let near = Surface {
            depth: 50.0,
            ..shaft
        };
        joined.deposit(near);
        joined.deposit(wall);
        assert_eq!(
            (joined.shaft, joined.normal.x, joined.normal.z),
            (1.0, 1.0, 0.0)
        );
    }

    #[test]
    fn smooths_adjoining_planes_and_stops_at_pores_and_occlusions() {
        let wall = Surface {
            depth: 12.0,
            radius: 12.0,
            weight: 1.0,
            normal: Point {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            ..Surface::EMPTY
        };
        let relief = Surface {
            normal: Point {
                x: 0.6,
                y: 0.0,
                z: 0.8,
            },
            shaft: 1.0,
            ..wall
        };
        let width = 42;
        let mut patch = vec![wall; width * 60];
        for row in patch.chunks_mut(width) {
            row[20..].fill(relief);
        }
        let index = 32 * width + 19;
        let rounded = local_surface(&patch, index, width, 19);
        assert!(rounded.normal.x > 0.05 && rounded.shaft > 0.05);
        assert_eq!(rounded.depth, wall.depth);
        let cropped = &patch[(32 - FILTER_MARGIN) * width..];
        assert_eq!(
            local_surface(cropped, FILTER_MARGIN * width + 19, width, 19)
                .normal
                .x,
            rounded.normal.x
        );

        for barrier in [
            Surface::EMPTY,
            Surface {
                depth: -100.0,
                ..wall
            },
        ] {
            for row in patch.chunks_mut(width) {
                row[20] = barrier;
            }
            let edge = local_surface(&patch, index, width, 19);
            assert_eq!((edge.normal.x, edge.shaft), (0.0, 0.0));
        }
    }
}
