# Install needletail dependency

```bash
# Install maturin >=0.14 & <0.15
sudo dnf install patchelf
pip install maturin==0.14
maturin build --features=python --release --strip
pip install ./target/wheels/needletail-0.5.1-cp312-cp312-manylinux_2_28_x86_64.whl
```
