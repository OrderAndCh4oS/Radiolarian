#!/bin/sh
set -eu

export PATH="${CARGO_HOME:-$HOME/.cargo}/bin:$PATH"
if ! command -v rustup >/dev/null 2>&1; then
    mkdir -p target
    curl --proto '=https' --tlsv1.2 -fsSL https://sh.rustup.rs -o target/rustup-init.sh
    sh target/rustup-init.sh -y --profile minimal --no-modify-path
fi

wasm_target=wasm32-unknown-unknown
rustup target add "$wasm_target"
cargo build --release --locked --target "$wasm_target"

site_dir=target/site
wasm_dir="target/$wasm_target/release"
mkdir -p "$site_dir/$wasm_dir"
cp index.html "$site_dir/"
cp "$wasm_dir/haeckel.wasm" "$site_dir/$wasm_dir/"
