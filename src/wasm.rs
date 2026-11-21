//! A worker owns one instance. Its output slices live until the next render or release.

use crate::{BodyPlan, Drawing, Settings, render};
use std::cell::RefCell;

thread_local! {
    static OUTPUT: RefCell<Option<Result<Drawing, String>>> = const { RefCell::new(None) };
}

#[unsafe(no_mangle)]
pub extern "C" fn render_specimen(
    width: u32,
    height: u32,
    seed: u64,
    flare: f32,
    spacing: f32,
    thickness: f32,
    age: f32,
    body_plan: u32,
    fine_spines: f32,
) -> bool {
    release_specimen();
    // JavaScript converts a missing float argument from older clients to NaN.
    let fine_spines = if fine_spines.is_nan() {
        Settings::default().fine_spines
    } else {
        fine_spines
    };
    let print = BodyPlan::try_from(body_plan)
        .map_err(str::to_owned)
        .and_then(|body_plan| {
            render(&Settings {
                width,
                height,
                seed,
                flare,
                spacing,
                thickness,
                age,
                body_plan,
                fine_spines,
            })
        });
    let success = print.is_ok();
    OUTPUT.set(Some(print));
    success
}

fn read_buffer<T>(kind: u32, read: impl FnOnce(&[u8]) -> T) -> T {
    OUTPUT.with_borrow(|output| {
        read(match (output.as_ref(), kind) {
            (Some(Ok(print)), 0) => &print.pixels,
            (Some(Ok(print)), 1) => print.samples.as_flattened(),
            (Some(Ok(print)), 2) => &print.tips,
            (Some(Err(message)), 3) => message.as_bytes(),
            _ => &[],
        })
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn specimen_buffer(kind: u32) -> *const u8 {
    read_buffer(kind, <[u8]>::as_ptr)
}

#[unsafe(no_mangle)]
pub extern "C" fn specimen_len(kind: u32) -> usize {
    read_buffer(kind, <[u8]>::len)
}

#[unsafe(no_mangle)]
pub extern "C" fn release_specimen() {
    OUTPUT.set(None);
}
