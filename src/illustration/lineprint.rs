//! Deposit traced strokes on the visible wall of the shell.

use super::{Point, Surface, Viewport};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Ink {
    Silver,
    Graphite,
    Charcoal,
}

#[derive(Clone, Copy)]
pub(super) struct Nib {
    pub p: Point,
    /// Fraction of the ink transferred to paper at this end of the stroke.
    pub opacity: f32,
}

#[derive(Clone, Copy)]
pub(super) struct Segment {
    /// Segments of one pencil stroke are emitted consecutively with the same ID.
    pub stroke_id: usize,
    /// Fixed pencil point radius in specimen units.
    pub radius: f32,
    pub a: Nib,
    pub b: Nib,
    pub depth_tolerance: f32,
    pub ink: Ink,
}

struct Ribbon {
    a: Point,
    delta: Point,
    inv_length_sq: f32,
    radius: f32,
    opacity: f32,
    opacity_change: f32,
    depth_tolerance: f32,
    stroke: (usize, Ink),
    left: usize,
    right: usize,
    top: usize,
    bottom: usize,
}

impl Ribbon {
    fn raster(
        &self,
        width: usize,
        rows: std::ops::Range<usize>,
        padded_top: usize,
        patch: &[Surface],
        mut deposit: impl FnMut(usize, u8),
    ) {
        for y in rows.start.max(self.top)..rows.end.min(self.bottom) {
            let dy = y as f32 + 0.5 - self.a.y;
            let row = (y - rows.start) * width;
            let wall = (y - padded_top) * width;
            for x in self.left..self.right {
                let surface = &patch[wall + x];
                if surface.weight == 0.0 {
                    continue;
                }
                let dx = x as f32 + 0.5 - self.a.x;
                let t =
                    ((dx * self.delta.x + dy * self.delta.y) * self.inv_length_sq).clamp(0.0, 1.0);
                let z = self.a.z + t * self.delta.z;
                if (z - surface.depth).abs() > self.depth_tolerance {
                    continue;
                }
                let side_x = dx - t * self.delta.x;
                let side_y = dy - t * self.delta.y;
                let distance = (side_x * side_x + side_y * side_y).sqrt();
                // The lead keeps the same footprint as pressure rises and falls.
                let contact = |rad: f32| (rad + 0.5 - distance).min(2.0 * rad).clamp(0.0, 1.0);
                let amount =
                    0.8 * contact(0.88 * self.radius) + 0.2 * contact(SHOULDER * self.radius);
                let opacity = self.opacity + t * self.opacity_change;
                deposit(row + x, (amount * opacity * 255.0).round() as u8);
            }
        }
    }
}

pub(super) struct PrintLines {
    ribbons: Vec<Ribbon>,
    bands: Vec<Vec<usize>>,
    width: usize,
    rows: usize,
}

const SHOULDER: f32 = 1.5;

// Keep the print and the live shader on the same side of half-grey steps.
const ONE: u32 = 4096;
const fn multiply(a: u32, b: u32) -> u32 {
    (a * b + ONE / 2) / ONE
}

// The reference drawing uses the same three pass strengths as the shader, at
// 100% overall pencil strength. A second pass takes up some remaining paper.
const PRINT_STRENGTHS: [u32; 3] = [176, 130, 108];
const fn pencil_transmission(percent: u32) -> [u32; 256] {
    let weight = (percent * ONE + 50) / 100;
    let first = if weight < ONE { weight } else { ONE };
    let mut table = [0; 256];
    let mut ink = 0;
    while ink < 256 {
        let coverage = (ink as u32 * ONE + 127) / 255;
        let second = multiply(coverage, ONE - coverage);
        let deposit = (coverage * first + second * (weight - first) + ONE / 2) / ONE;
        table[ink] = ONE - deposit;
        ink += 1;
    }
    table
}
const TRANSMISSION: [[u32; 256]; 3] = [
    pencil_transmission(PRINT_STRENGTHS[0]),
    pencil_transmission(PRINT_STRENGTHS[1]),
    pencil_transmission(PRINT_STRENGTHS[2]),
];

pub(super) fn print_inks(paper: u8, inks: [u8; 3]) -> u8 {
    let mut remaining = ONE;
    for (pass, ink) in TRANSMISSION.iter().zip(inks) {
        remaining = multiply(remaining, pass[usize::from(ink)]);
    }
    ((u32::from(paper) * remaining + ONE / 2) / ONE) as u8
}

impl PrintLines {
    pub(super) fn new(
        segments: Vec<Segment>,
        view: Viewport,
        width: usize,
        height: usize,
        rows: usize,
    ) -> Self {
        let mut bands = vec![Vec::new(); height.div_ceil(rows)];
        let mut ribbons = Vec::with_capacity(segments.len());
        for segment in segments {
            let mut a = view.screen(segment.a.p);
            let b = view.screen(segment.b.p);
            let radius = (segment.radius * view.scale).max(0.0);
            let mut opacity = segment.a.opacity.clamp(0.0, 1.0);
            let end_opacity = segment.b.opacity.clamp(0.0, 1.0);
            if opacity == 0.0 && end_opacity == 0.0 {
                continue;
            }
            let delta = b.minus(a);
            let length_sq = delta.x * delta.x + delta.y * delta.y;
            let inv_length_sq = if length_sq > 1.0e-8 {
                length_sq.recip()
            } else {
                // Draw coincident endpoints with the stronger nib impression.
                opacity = opacity.max(end_opacity);
                a.z = a.z.max(b.z);
                0.0
            };
            let extent = SHOULDER * radius + 0.5;
            let left = ((a.x.min(b.x) - extent).floor() as usize).min(width);
            let right = ((a.x.max(b.x) + extent).ceil() as usize).min(width);
            let top = ((a.y.min(b.y) - extent).floor() as usize).min(height);
            let bottom = ((a.y.max(b.y) + extent).ceil() as usize).min(height);
            if left >= right || top >= bottom || extent <= 0.5 {
                continue;
            }
            for band in &mut bands[top / rows..=(bottom - 1) / rows] {
                band.push(ribbons.len());
            }
            ribbons.push(Ribbon {
                a,
                delta,
                inv_length_sq,
                radius,
                opacity,
                opacity_change: end_opacity - opacity,
                depth_tolerance: (segment.depth_tolerance * view.scale).max(0.0),
                stroke: (segment.stroke_id, segment.ink),
                left,
                right,
                top,
                bottom,
            });
        }
        Self {
            ribbons,
            bands,
            width,
            rows,
        }
    }

    pub(super) fn draw(
        &self,
        top: usize,
        bottom: usize,
        padded_top: usize,
        patch: &[Surface],
        pigment: &mut [u8],
        proof: &mut [[u8; 4]],
    ) {
        if top >= bottom {
            return;
        }
        for band_id in top / self.rows..=(bottom - 1) / self.rows {
            let band = &self.bands[band_id];
            let first = (band_id * self.rows).max(top);
            let last = ((band_id + 1) * self.rows).min(bottom);
            let area = (first - top) * self.width..(last - top) * self.width;
            let pixels = &mut pigment[area.clone()];
            let wall = &patch[(first - padded_top) * self.width..(last - padded_top) * self.width];
            // Until the final resolve, pixels stores the local graphite load.
            // The only way that value can reach the print is through a touched nib.
            let mut coverage = vec![[0u8; 3]; (last - first) * self.width];
            let mut nib = vec![0u8; coverage.len()];
            let mut touched = Vec::new();
            for stroke in band.chunk_by(|&i, &j| self.ribbons[i].stroke == self.ribbons[j].stroke) {
                let ink = self.ribbons[stroke[0]].stroke.1 as usize;
                // Segments share the footprint of one continuous nib, including
                // caps and joins. Collect the strongest impression at each sample.
                for &id in stroke {
                    self.ribbons[id].raster(
                        self.width,
                        first..last,
                        padded_top,
                        patch,
                        |i, impression| {
                            if impression > nib[i] {
                                if nib[i] == 0 {
                                    touched.push(i);
                                }
                                nib[i] = impression;
                            }
                        },
                    );
                }
                // Each new pencil pass covers a fraction of the remaining paper.
                // Reset just the samples touched by this stroke for the next one.
                for i in touched.drain(..) {
                    let contact = u16::from(std::mem::take(&mut nib[i]));
                    let deposit = (contact * u16::from(pixels[i]) + 127) / 255;
                    let previous = &mut coverage[i][ink];
                    let paper = u16::from(255 - *previous);
                    *previous += ((paper * deposit + 127) / 255) as u8;
                }
            }
            let samples = &mut proof[area];
            for (i, mut ink) in coverage.into_iter().enumerate() {
                let axial = u16::from(wall[i].axial_graphite(pixels[i]));
                let graphite = &mut ink[Ink::Graphite as usize];
                *graphite += ((u16::from(255 - *graphite) * axial + 127) / 255) as u8;
                let paper = if wall[i].weight > 0.0 { 255 } else { 0 };
                pixels[i] = print_inks(paper, ink);
                samples[i] = [paper, ink[0], ink[1], ink[2]];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_deposited_graphite_can_change_the_paper() {
        for channel in 0..3 {
            let mut previous = 255;
            for amount in 0..=255 {
                let mut ink = [0; 3];
                ink[channel] = amount;
                let value = print_inks(255, ink);
                assert!(value <= previous && previous - value <= 2);
                assert!(value <= 255 - amount);
                assert_eq!(print_inks(0, ink), 0);
                previous = value;
            }
            assert_eq!(previous, 0);
        }
        assert_eq!(print_inks(255, [0; 3]), 255);
        assert_eq!(
            [
                print_inks(255, [128, 0, 0]),
                print_inks(255, [0, 128, 0]),
                print_inks(255, [0, 0, 128])
            ],
            [79, 108, 122]
        );
    }

    #[test]
    fn a_pencil_carries_the_local_value_over_bare_white_paper() {
        let nib = |x| Nib {
            p: Point { x, y: 0.5, z: 0.0 },
            opacity: 1.0,
        };
        let stroke = Segment {
            stroke_id: 0,
            radius: 0.8,
            a: nib(0.5),
            b: nib(255.5),
            ink: Ink::Graphite,
            depth_tolerance: 0.0,
        };
        let view = Viewport {
            x: 0.0,
            y: 0.0,
            scale: 1.0,
        };
        let print = PrintLines::new(vec![stroke], view, 256, 3, 3);
        let wall = [Surface {
            depth: 0.0,
            weight: 1.0,
            ..Surface::EMPTY
        }; 768];
        let mut pixels: [u8; 768] = std::array::from_fn(|i| (i % 256) as u8);
        print.draw(0, 3, 0, &wall, &mut pixels, &mut [[0; 4]; 768]);
        for x in 0..256 {
            assert_eq!(pixels[x], print_inks(255, [0, x as u8, 0]));
            assert_eq!(pixels[512 + x], 255);
        }
    }

    #[test]
    fn a_coincident_nib_leaves_a_soft_round_cap_and_preserves_paper() {
        let a = Nib {
            p: Point {
                x: 1.5,
                y: 0.5,
                z: 0.0,
            },
            opacity: 0.0,
        };
        let b = Nib { opacity: 1.0, ..a };
        let cap = Segment {
            stroke_id: 0,
            radius: 0.8,
            a,
            b,
            depth_tolerance: 0.0,
            ink: Ink::Graphite,
        };
        let view = Viewport {
            x: 0.0,
            y: 0.0,
            scale: 1.0,
        };
        let inks = [Ink::Charcoal, Ink::Graphite, Ink::Silver];
        let print = PrintLines::new(
            inks.map(|ink| Segment { ink, ..cap }).to_vec(),
            view,
            3,
            1,
            1,
        );
        let patch = [Surface {
            depth: 0.0,
            weight: 1.0,
            ..Surface::EMPTY
        }; 3];
        let mut pigment = [140, 140, 0];
        print.draw(0, 1, 0, &patch, &mut pigment, &mut [[0; 4]; 3]);
        assert!(pigment[1] < pigment[0] && pigment[0] < 255);
        assert_eq!(pigment[1], print_inks(255, [140; 3]));
        assert_eq!(pigment[2], 255);
        let mut highlight = [0; 3];
        print.draw(0, 1, 0, &patch, &mut highlight, &mut [[0; 4]; 3]);
        assert_eq!(highlight, [255; 3]);
    }

    #[test]
    fn pressure_changes_graphite_density_across_the_same_nib_footprint() {
        let nib = |x, opacity| Nib {
            p: Point { x, y: 2.5, z: 0.0 },
            opacity,
        };
        let stroke = Segment {
            stroke_id: 0,
            radius: 1.0,
            a: nib(0.5, 0.0),
            b: nib(4.5, 1.0),
            ink: Ink::Charcoal,
            depth_tolerance: 0.0,
        };
        let view = Viewport {
            scale: 1.0,
            x: 0.0,
            y: 0.0,
        };
        let lines = PrintLines::new(vec![stroke, stroke], view, 5, 5, 1);
        let load = 255;
        let mut pigment = [load; 25];
        let patch = [Surface {
            depth: 0.0,
            weight: 1.0,
            ..Surface::EMPTY
        }; 25];
        let mut samples = [[0; 4]; 25];
        lines.draw(0, 5, 0, &patch, &mut pigment, &mut samples);
        assert_eq!(
            pigment[10..15],
            [0, 64, 128, 191, 255].map(|alpha| print_inks(255, [0, 0, alpha]))
        );
        // The softer end and the harder end have identical cross-sections
        // after normalizing for the quantity of graphite deposited.
        let full_shoulder = samples[9][3];
        assert!(full_shoulder > 0 && full_shoulder < 255);
        for x in 0..5 {
            let core = f32::from(samples[10 + x][3]);
            let shoulder = samples[5 + x][3];
            assert_eq!(samples[x][3], 0);
            assert!((f32::from(shoulder) - core * f32::from(full_shoulder) / 255.0).abs() <= 1.0);
            assert_eq!(samples[15 + x][3], shoulder);
            assert_eq!(samples[20 + x][3], 0);
        }
    }

    #[test]
    fn crossing_pencil_strokes_deepen_each_separation_without_darkening_its_joins() {
        let nib = |x, y, opacity| Nib {
            p: Point { x, y, z: 0.0 },
            opacity,
        };
        let view = Viewport {
            x: 0.0,
            y: 0.0,
            scale: 1.0,
        };
        let wall = [Surface {
            depth: 0.0,
            weight: 1.0,
            ..Surface::EMPTY
        }; 49];
        for ink in [Ink::Silver, Ink::Graphite, Ink::Charcoal] {
            let pull = Segment {
                stroke_id: 0,
                radius: 0.8,
                a: nib(0.5, 3.5, 0.5),
                b: nib(3.5, 3.5, 0.5),
                ink,
                depth_tolerance: 0.0,
            };
            let joined = Segment {
                a: pull.b,
                b: nib(6.5, 3.5, 0.5),
                ..pull
            };
            let crossed = Segment {
                stroke_id: 2,
                a: nib(3.5, 0.5, 0.25),
                b: nib(3.5, 6.5, 0.25),
                ..pull
            };
            let lines = PrintLines::new(vec![pull, joined, crossed], view, 7, 7, 3);
            let mut samples = [[0; 4]; 49];
            let mut shade = [255; 49];
            lines.draw(0, 7, 0, &wall, &mut shade, &mut samples);
            // Two joined halves of a 50% stroke meet exactly under a 25% stroke:
            // the intersection covers 62.5% of the paper (160/255 after rounding).
            for (i, alpha) in [(0, 0), (21, 128), (3, 64), (24, 160)] {
                let mut expected = [0; 3];
                expected[ink as usize] = alpha;
                assert_eq!(&samples[i][1..], &expected);
                assert_eq!(shade[i], print_inks(255, expected));
            }
            assert!(shade[24] < shade[21] && shade[21] < shade[3]);
        }
    }

    #[test]
    fn joins_and_render_strip_boundaries_leave_no_extra_ink() {
        let nib = |x, y| Nib {
            p: Point { x, y, z: 1.0 },
            opacity: 1.0,
        };
        let stroke = Segment {
            stroke_id: 0,
            radius: 0.8,
            a: nib(0.5, 0.5),
            b: nib(5.5, 5.5),
            depth_tolerance: 0.1,
            ink: Ink::Graphite,
        };
        let view = Viewport {
            x: 0.0,
            y: 0.0,
            scale: 1.0,
        };
        let printed = PrintLines::new(vec![stroke], view, 6, 6, 3);
        let repeated = PrintLines::new(vec![stroke, stroke], view, 6, 6, 3);
        let mut patch = vec![
            Surface {
                depth: 1.0,
                weight: 1.0,
                ..Surface::EMPTY
            };
            36
        ];
        patch[7] = Surface::EMPTY;
        patch[14].depth = -3.0;
        let mut whole = [160; 36];
        let mut whole_samples = [[0; 4]; 36];
        printed.draw(0, 6, 0, &patch, &mut whole, &mut whole_samples);
        let mut strips = [160; 36];
        let mut strip_samples = [[0; 4]; 36];
        repeated.draw(0, 3, 0, &patch, &mut strips[..18], &mut strip_samples[..18]);
        repeated.draw(
            3,
            6,
            2,
            &patch[12..],
            &mut strips[18..],
            &mut strip_samples[18..],
        );
        assert_eq!(whole, strips);
        assert_eq!(whole_samples, strip_samples);
        assert_eq!((whole[7], whole[14]), (0, 255));
        assert!(whole[0] < 160 && whole[35] < 160);
        assert_eq!(whole[5], 255);
    }
}
