# Radiolarian

Radiolarian is a browser-based drawing tool for exploring the silica skeletons of radiolarians. A Rust WebAssembly worker grows a chambered bell or concentric sphere and traces the pencil strokes along its surface. A WebGPU shader combines the hatching layers so you can adjust their strength on the finished drawing.

## Run locally

From the project directory, build the renderer with Rust and serve the viewer with Python 3.

```sh
rustup target add wasm32-unknown-unknown
cargo build --release --locked --target wasm32-unknown-unknown
python3 -m http.server 8000 --bind 127.0.0.1
```

Open the [local viewer](http://127.0.0.1:8000/). Specimen 31 appears immediately while the worker prepares the interactive drawing. A browser with WebGPU enables live shading. Browsers with WebAssembly and web workers can also display the reference drawing and generate new specimens.

## Draw a specimen

- Choose a body plan, adjust its proportions or growth, and select **Apply form**. The camera stays fixed through the growth stages of a specimen.
- Enter a seed to reproduce a drawing with the same settings. **Randomise seed** starts a new preview using the current form. A progress indicator shows which specimen is rendering.
- Adjust **Fine hatching**, **Rib hatching**, **Deep crosshatching**, or **Overall pencil strength**. Select **Reference drawing** to compare the result with the original pencil settings.
- The preview is 960 × 1200 pixels. **Render at 2400 × 3000** produces the full-size drawing. **Pore detail** lets you inspect the strokes by dragging or using the arrow keys.

## Test the renderer

```sh
cargo test --lib --release --locked
```
