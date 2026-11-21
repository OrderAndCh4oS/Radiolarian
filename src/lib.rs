//! Growth and monochrome rendering for a single radiolarian specimen.

mod growth;
mod illustration;
mod wasm;

/// Body plan indices used by the browser worker.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u32)]
enum BodyPlan {
    #[default]
    Bell = 0,
    Sphere = 1,
}

impl TryFrom<u32> for BodyPlan {
    type Error = &'static str;

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Bell),
            1 => Ok(Self::Sphere),
            _ => Err("unknown body plan"),
        }
    }
}

/// Flare adjusts bell width or sphere spine length for the selected plan.
#[derive(Clone, Copy, Debug)]
struct Settings {
    body_plan: BodyPlan,
    width: u32,
    height: u32,
    seed: u64,
    flare: f32,
    spacing: f32,
    thickness: f32,
    fine_spines: f32,
    age: f32,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            body_plan: BodyPlan::Bell,
            width: 2400,
            height: 3000,
            seed: 31,
            flare: 1.0,
            spacing: 1.0,
            thickness: 1.0,
            fine_spines: 1.0,
            age: 1.0,
        }
    }
}

impl Settings {
    fn validate(&self) -> Result<(), String> {
        if !(1..=16_384).contains(&self.width) || !(1..=16_384).contains(&self.height) {
            return Err("width and height must each be between 1 and 16,384 pixels".into());
        }
        if u64::from(self.width) * u64::from(self.height) > 16_777_216 {
            return Err("choose dimensions of at most 16,777,216 pixels in total".into());
        }
        for (name, value) in [
            ("flare", self.flare),
            ("spacing", self.spacing),
            ("thickness", self.thickness),
        ] {
            if !(0.6..=1.6).contains(&value) {
                return Err(format!("{name} must be between 0.6 and 1.6"));
            }
        }
        if !(0.0..=1.0).contains(&self.fine_spines) {
            return Err("fine-spines must be between 0 and 1".into());
        }
        if !(0.0..=1.0).contains(&self.age) {
            return Err("age must be between 0 and 1".into());
        }
        Ok(())
    }
}

struct Drawing {
    /// Reference drawing, one grayscale byte per pixel.
    pixels: Vec<u8>,
    /// Four supersamples per pixel, each storing paper, silver, graphite and charcoal.
    samples: Vec<[u8; 4]>,
    /// Spine tip records consumed by the shader, in sample and draw order.
    tips: Vec<u8>,
}

/// Growth stages share the mature specimen's camera framing.
fn render(settings: &Settings) -> Result<Drawing, String> {
    settings.validate()?;
    let at_age = |age| growth::specimen(settings, age);
    let mature = at_age(1.0);
    let developing = (settings.age < 1.0).then(|| at_age(settings.age));
    let specimen = developing.as_deref().unwrap_or(&mature);
    Ok(illustration::draw(
        specimen,
        &mature,
        settings.width,
        settings.height,
        settings.seed,
        settings.body_plan,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_a_lit_specimen_on_a_black_field() {
        let settings = Settings {
            width: 180,
            height: 220,
            ..Settings::default()
        };
        let image = render(&settings).unwrap().pixels;
        assert_eq!(image.len(), 180 * 220);
        assert!(image[..180].iter().all(|&v| v == 0));
        assert!(image.iter().any(|&v| v > 128));
        assert_eq!(render(&settings).unwrap().pixels, image);
    }

    #[test]
    fn rejects_nonfinite_parameters_and_oversized_images_before_allocating() {
        for settings in [
            Settings {
                age: f32::NAN,
                ..Settings::default()
            },
            Settings {
                spacing: f32::INFINITY,
                ..Settings::default()
            },
            Settings {
                width: 0,
                ..Settings::default()
            },
            Settings {
                width: 8192,
                height: 8192,
                ..Settings::default()
            },
        ] {
            assert!(render(&settings).is_err());
        }
    }
}
