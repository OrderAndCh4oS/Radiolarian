use super::Surface;
use std::f32::consts::SQRT_2;

#[derive(Clone, Copy, Default)]
pub(super) struct Rim {
    pub ink: f32,
}

/// The near wall owns each depth break. A fine line describes its shaded
/// lip and trails off as the perimeter turns towards the light.
pub(super) fn edge(samples: &[Surface], i: usize, width: usize, pen_width: f32) -> Rim {
    let center = samples[i];
    if center.weight == 0.0 || center.radius < 1.5 || pen_width <= 0.0 {
        return Rim::default();
    }

    let x = (i % width) as isize;
    let y = (i / width) as isize;
    let height = (samples.len() / width) as isize;
    // The rasterizer pads its strips with adjacent rows. Cap the search inside
    // that margin and leave narrow spines enough room for an unmarked face.
    let pen = (0.78 + 0.10 * center.shaft) * pen_width;
    let reach = pen.min(0.32 * center.radius).min(6.0);
    let mut outward = (0.0, 0.0);
    let mut ink_width = 0.0_f32;
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
        let stride = if dx != 0 && dy != 0 { SQRT_2 } else { 1.0 };
        let steps = (reach / stride).ceil() as isize;
        for step in 1..=steps {
            let nx = x + step * dx;
            let ny = y + step * dy;
            if nx < 0 || nx >= width as isize || ny < 0 || ny >= height {
                break;
            }
            let other = samples[ny as usize * width + nx as usize];
            if other.weight > 0.0 {
                // An inclined continuous surface can recede substantially over
                // a few samples. Only mark a step beyond its expected slope.
                let slope = (center.normal.x * dx as f32 + center.normal.y * dy as f32)
                    / center.normal.z.max(0.12);
                let tolerance = 0.58 * center.radius + step as f32 * (stride + slope.max(0.0));
                let drop = center.depth - other.depth;
                if drop < -tolerance {
                    break;
                }
                if drop <= tolerance {
                    continue;
                }
            }
            let offset = ((step - 1) as f32) * stride;
            let coverage = ((reach - offset) / stride).clamp(0.0, 1.0);
            outward.0 += coverage * dx as f32 / stride;
            outward.1 += coverage * dy as f32 / stride;
            ink_width = ink_width.max(coverage);
            break;
        }
    }
    let length = outward.0.hypot(outward.1);
    if length < 0.001 {
        return Rim::default();
    }
    let shadow_side = (0.65 * outward.0 + 0.76 * outward.1) / length;
    let pressure = ((center.radius - 1.5) / 3.0).clamp(0.0, 1.0);
    Rim {
        ink: ink_width * pressure * shadow_side.clamp(0.0, 1.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::illustration::Point;

    #[test]
    fn near_wall_owns_the_mark_on_its_shaded_lip() {
        let near = Surface {
            depth: 30.0,
            radius: 8.0,
            weight: 1.0,
            normal: Point {
                z: 1.0,
                ..Point::default()
            },
            ..Surface::EMPTY
        };
        let far = Surface { depth: 0.0, ..near };
        let mut samples = [far; 15 * 15];
        for row in samples.as_chunks_mut::<15>().0 {
            row[4..11].fill(near);
        }
        let at = |x: usize| edge(&samples, 7 * 15 + x, 15, 2.0);
        let left = at(4);
        let right = at(10);
        assert_eq!(left.ink, 0.0);
        assert!(right.ink > 0.5);
        for x in [3, 11] {
            let mark = at(x);
            assert_eq!(mark.ink, 0.0);
        }

        for row in samples.as_chunks_mut::<15>().0 {
            row[7..11].fill(Surface { shaft: 1.0, ..near });
        }
        let root = edge(&samples, 7 * 15 + 7, 15, 2.0);
        assert_eq!(root.ink, 0.0);
    }
}
