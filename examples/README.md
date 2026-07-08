# Examples

## Demo GIF (`docs/demo.gif`)

The engine itself is headless — it computes physics, it does not draw. The
demo GIF is produced by a tiny two-step offline pipeline that lives here,
so it stays out of the core library and adds no dependencies to it:

1. **`dump_scene.cpp`** links `melon3d`, simulates the showcase scene
   (a wrecking ball on a chain of ball joints demolishing a box pyramid,
   then a rain of mixed shapes settling over the debris) and writes every
   frame's body transforms to a binary file `scene.bin`.
2. **`render.py`** reads `scene.bin`, renders each frame with a small
   perspective rasterizer (shaded box faces, painter's-algorithm sorting,
   soft blob shadows, a floor grid, 3x supersampling) and writes an
   optimized `melon3d.gif` (shared palette + inter-frame transparency
   diffing, so the static background costs almost nothing).

### Build and run

```sh
# from the repo root, after building the library into build/
g++ -O2 -std=c++17 -Iinclude examples/dump_scene.cpp -Lbuild -lmelon3d -o dump_scene
./dump_scene              # writes scene.bin

python -m pip install pillow numpy
python examples/render.py preview   # a few PNG stills to check the look
python examples/render.py           # writes melon3d.gif
```

The renderer needs only Pillow and NumPy (Python 3.8+). Nothing here is
part of the shipped library; it is purely a visualization aid.
