//! Use chamber and rib relief to guide the pressure of a pencil in contact with the shell.

use super::{LIGHT, Point, Surface, Viewport, contour::Rim, envelope::Envelope};

fn smooth(lo: f32, hi: f32, x: f32) -> f32 {
    let t = ((x - lo) / (hi - lo)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// Increase the graphite on the turning face, easing the transition at both ends.
fn form_value(tone: f32, front: f32) -> f32 {
    let halftone = smooth(0.12, 0.52, tone) * (1.0 - smooth(0.66, 0.86, tone));
    tone - 0.18 * front * halftone
}

pub(super) struct Plate {
    pub view: Viewport,
    pub envelope: Envelope,
}

impl Plate {
    /// Return the target pencil tone, chamber normal and near-wall weight.
    pub fn tone(&self, p: Point, normal: Point, exposure: f32) -> (f32, Point, f32) {
        let shell = self.envelope.sample(p);
        let length = normal.dot(normal).max(0.01).sqrt();
        let detail = normal.dot(LIGHT) / length;
        let front = smooth(-0.27, 0.29, shell.front);
        let field = shell.normal.dot(Point {
            x: -0.79,
            y: -0.49,
            z: 0.37,
        });
        let chamber_light = (0.5 + 0.5 * field).clamp(0.0, 1.0);
        let rib_light = (0.5 + 0.5 * detail).clamp(0.0, 1.0);
        // Each rib modulates the chamber value used in spacing and loading strokes.
        let back = 0.29 + 0.115 * chamber_light + 0.075 * rib_light;
        // Let the side turning away from the light fall into a broad half-tone.
        // This is shared by neighbouring bars across an entire chamber.
        let face = 0.45 + 0.40 * chamber_light + 0.075 * rib_light;
        let value = back + (face - back) * front;
        let alignment = normal.dot(shell.normal) / length;
        let turned = 1.0 - smooth(-0.08, 0.61, alignment);
        let bevel = 1.0 - smooth(-0.66, 0.13, detail);
        let recess = (0.025 + 0.095 * front) * turned * bevel;
        let collar = 0.042 * smooth(0.035, 0.33, shell.collar_slope);
        let wall = value - recess - collar;

        // The neck of each outgrowth shares the surrounding wall tone. Use
        // more of the local relief as the exposed shaft rises from the ribs.
        let shaft = 0.17 + 0.54 * value + 0.25 * rib_light;
        let surface = wall + exposure * (shaft - wall);
        (form_value(surface, front), shell.normal, front)
    }

    pub fn pencil_load(&self, s: Surface, x: usize, y: usize, rim: Rim) -> u8 {
        let view = self.view;
        let p = Point {
            x: (x as f32 + 0.5 - view.x) / view.scale,
            y: (y as f32 + 0.5 - view.y) / view.scale,
            z: s.depth / view.scale,
        };
        let (value, _, front) = self.tone(p, s.normal, s.shaft);
        let tone = ((value - 0.12 * (1.0 - s.tone)) * (1.0 - 0.27 * rim.ink)).clamp(0.0, 1.0);
        pencil_load(tone, front)
    }
}

// Fraction of graphite a full-pressure stroke may transfer at this position.
// A smooth shoulder eases the pencil off the lit crests.
fn pencil_load(tone: f32, front: f32) -> u8 {
    let lift = front * smooth(0.755, 0.86, tone);
    let transfer = 0.83 - 0.13 * smooth(0.34, 0.65, tone);
    (255.0 * (1.0 - tone) * transfer * (1.0 - lift)).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graphite_load_eases_to_zero_at_the_lit_crest() {
        let mut previous = 255;
        for i in 0..=255 {
            let load = pencil_load(i as f32 / 255.0, 1.0);
            assert!(load <= previous);
            previous = load;
        }
        assert_eq!(pencil_load(0.90, 1.0), 0);
        assert!(pencil_load(0.70, 0.3) < pencil_load(0.35, 0.3));
    }
}
