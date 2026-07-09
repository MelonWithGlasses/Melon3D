import struct, math, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

SS = 3                      # supersample factor
W, H = 760, 560
WS, HS = W * SS, H * SS
FOV = math.radians(46)
FOCAL = (HS / 2) / math.tan(FOV / 2)

EYE = np.array([13.6, 8.2, 17.6])
TARGET = np.array([0.0, 3.9, 0.0])
UP = np.array([0.0, 1.0, 0.0])
LIGHT = np.array([-0.35, 1.0, 0.55]); LIGHT = LIGHT / np.linalg.norm(LIGHT)
AMBIENT = 0.40

# ---- camera basis (look along -z) ----
fwd = TARGET - EYE; fwd = fwd / np.linalg.norm(fwd)
right = np.cross(fwd, UP); right = right / np.linalg.norm(right)
up = np.cross(right, fwd)
VIEW = np.stack([right, up, -fwd])   # rows map world->view

def to_view(pw):
    return VIEW @ (pw - EYE)

def project(pv):
    z = -pv[2]
    if z <= 0.05:
        return None, z
    sx = WS / 2 + FOCAL * pv[0] / z
    sy = HS / 2 - FOCAL * pv[1] / z
    return (sx, sy), z

def qrot(q, v):
    x, y, z, w = q
    u = np.array([x, y, z])
    t = 2.0 * np.cross(u, v)
    return v + w * t + np.cross(u, t)

PALETTE = {
    1: (86, 196, 138),    # pyramid - melon green
    2: (242, 90, 78),     # wrecking ball - coral
    3: (84, 182, 255),    # rain sphere - sky
    4: (255, 188, 84),    # rain capsule - amber
    5: (172, 126, 255),   # rain box - purple
    6: (150, 160, 176),   # chain link - steel
    7: (98, 106, 122),    # anchor - dark steel
}

# screen-space light direction for the fake sphere highlight
_lv = VIEW @ LIGHT
_llen = math.hypot(_lv[0], _lv[1])
LX, LY = (_lv[0] / _llen, -_lv[1] / _llen) if _llen > 1e-6 else (-0.5, -0.5)

def hue_shift(c, i):
    # small per-index variation so identical bodies read as distinct
    f = 0.88 + 0.24 * ((i * 0.61803398875) % 1.0)
    return tuple(min(255, int(v * f)) for v in c)

def shade(base, n):
    d = AMBIENT + (1 - AMBIENT) * max(0.0, float(np.dot(n, LIGHT)))
    return tuple(min(255, int(v * d)) for v in base)

BOX_CORNERS = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                        [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], dtype=float)
BOX_FACES = [([1,2,6,5],[1,0,0]),([0,4,7,3],[-1,0,0]),
             ([3,7,6,2],[0,1,0]),([0,1,5,4],[0,-1,0]),
             ([4,5,6,7],[0,0,1]),([0,3,2,1],[0,0,-1])]

def read_scene(path):
    with open(path, 'rb') as f:
        numFrames, numBodies = struct.unpack('ii', f.read(8))
        bodies = []
        for _ in range(numBodies):
            st, r, hh, ex, ey, ez, col, isst = struct.unpack('ifffffii', f.read(32))
            bodies.append(dict(st=st, r=r, hh=hh, e=np.array([ex,ey,ez]), col=col, static=isst))
        frames = np.frombuffer(f.read(numFrames*numBodies*7*4), dtype=np.float32)
        frames = frames.reshape(numFrames, numBodies, 7)
    return bodies, frames

# background with a vertical gradient and a soft radial vignette,
# rendered once and reused for every frame
def make_background():
    yy = (np.arange(HS) / HS)[:, None]
    xx = (np.arange(WS) / WS)[None, :]
    top = np.array([24, 28, 38]); bot = np.array([46, 52, 66])
    col = top[None, None, :] * (1 - yy)[:, :, None] + bot[None, None, :] * yy[:, :, None]
    # vignette: darken the corners a touch
    r2 = ((xx - 0.5) * 1.15) ** 2 + ((yy - 0.45) * 1.15) ** 2
    vig = 1.0 - 0.35 * np.clip(r2 * 1.6, 0, 1)
    col = col * vig[:, :, None]
    return Image.fromarray(col.astype(np.uint8))

BACKGROUND = None

def draw_floor(draw):
    # filled floor quad + grid lines at y=0, fading away from the origin
    ext = 19.0
    quad = []
    ok = True
    for (x, z) in [(-ext,-ext),(ext,-ext),(ext,ext),(-ext,ext)]:
        p, _ = project(to_view(np.array([x, 0.0, z])))
        if p is None: ok = False; break
        quad.append(p)
    if ok:
        draw.polygon(quad, fill=(37, 43, 55))
    step = 2.0
    n = int(ext / step)
    for i in range(-n, n+1):
        c = i * step
        fade = 1.0 - 0.55 * abs(i) / n
        col = tuple(int(v * fade + 37 * (1 - fade)) for v in (60, 68, 83))
        for a, b in [((c,-ext),(c,ext)), ((-ext,c),(ext,c))]:
            pa, _ = project(to_view(np.array([a[0],0.0,a[1]])))
            pb, _ = project(to_view(np.array([b[0],0.0,b[1]])))
            if pa and pb:
                draw.line([pa, pb], fill=col, width=SS)

def draw_shadows(img, bodies, frame):
    # soft blob shadows: project each body onto y=0, alpha falls with height
    ov = Image.new('L', (WS, HS), 0)
    d = ImageDraw.Draw(ov)
    for bi, b in enumerate(bodies):
        if b['static']:
            continue
        p = frame[bi, 0:3].astype(float)
        if p[1] < -0.5 or p[1] > 14.0:
            continue
        if b['st'] == 0:
            rad = b['r']
        elif b['st'] == 1:
            rad = b['r'] + 0.5 * b['hh']
        else:
            rad = float(np.mean(b['e'])) * 1.15
        h = max(0.0, p[1] - rad)
        alpha = 95.0 / (1.0 + 0.55 * h)
        spread = 1.0 + 0.16 * h
        # shadow ellipse in world space, projected as a polygon
        cx, cz = p[0], p[2]
        pts = []
        okp = True
        for k in range(14):
            a = 2 * math.pi * k / 14
            wp = np.array([cx + rad * spread * math.cos(a), 0.02, cz + rad * spread * math.sin(a)])
            sp, _ = project(to_view(wp))
            if sp is None:
                okp = False
                break
            pts.append(sp)
        if okp:
            d.polygon(pts, fill=int(alpha))
    ov = ov.filter(ImageFilter.GaussianBlur(3 * SS))
    black = Image.new('RGB', (WS, HS), (10, 12, 16))
    img.paste(black, (0, 0), ov)

def render_frame(bodies, frame, idx):
    global BACKGROUND
    if BACKGROUND is None:
        BACKGROUND = make_background()
    img = BACKGROUND.copy()
    draw = ImageDraw.Draw(img)
    draw_floor(draw)
    draw_shadows(img, bodies, frame)
    draw = ImageDraw.Draw(img)

    prims = []  # (depth, kind, payload)
    for bi, b in enumerate(bodies):
        if b['col'] == 0:
            continue  # ground drawn as floor grid
        p = frame[bi, 0:3].astype(float)
        q = frame[bi, 3:7].astype(float)
        base = hue_shift(PALETTE[b['col']], bi)
        if b['st'] == 0:  # sphere
            pv = to_view(p); sc, z = project(pv)
            if sc is None: continue
            rad = FOCAL * b['r'] / max(0.05, -pv[2])
            prims.append((-pv[2], 'sphere', (sc, rad, base)))
        elif b['st'] == 1:  # capsule (local Y axis)
            axis = qrot(q, np.array([0.0, b['hh'], 0.0]))
            c1 = p - axis; c2 = p + axis
            pv = to_view(p)
            s1, _ = project(to_view(c1)); s2, _ = project(to_view(c2))
            if s1 is None or s2 is None: continue
            rad = FOCAL * b['r'] / max(0.05, -pv[2])
            prims.append((-pv[2], 'capsule', (s1, s2, rad, base)))
        else:  # box
            corners = np.array([p + qrot(q, c * b['e']) for c in BOX_CORNERS])
            cv = np.array([to_view(c) for c in corners])
            for idxs, ln in BOX_FACES:
                nrm = qrot(q, np.array(ln, dtype=float))
                fc = corners[idxs].mean(axis=0)
                if np.dot(nrm, EYE - fc) <= 0:
                    continue
                poly = []
                zsum = 0.0; okf = True
                for k in idxs:
                    sc, z = project(cv[k])
                    if sc is None: okf = False; break
                    poly.append(sc); zsum += -cv[k][2]
                if not okf: continue
                prims.append((zsum/4, 'poly', (poly, shade(base, nrm))))

    prims.sort(key=lambda t: -t[0])  # farthest first
    for _, kind, pl in prims:
        if kind == 'poly':
            poly, col = pl
            edge = tuple(int(v * 0.45) for v in col)
            draw.polygon(poly, fill=col, outline=edge)
        elif kind == 'sphere':
            sc, rad, base = pl
            _draw_ball(draw, sc, rad, base)
        elif kind == 'capsule':
            s1, s2, rad, base = pl
            # connecting body first, then the end caps on top
            dx, dy = s2[0]-s1[0], s2[1]-s1[1]
            L = math.hypot(dx, dy)
            if L > 1e-3:
                nx, ny = -dy/L*rad, dx/L*rad
                bodycol = shade(base, np.array([0, 0.35, 0.94]))
                quad = [(s1[0]+nx,s1[1]+ny),(s2[0]+nx,s2[1]+ny),(s2[0]-nx,s2[1]-ny),(s1[0]-nx,s1[1]-ny)]
                draw.polygon(quad, fill=bodycol)
                # simple two-tone: a lighter core strip along the axis
                lt = tuple(min(255, int(v * 1.25)) for v in bodycol)
                q2 = [(s1[0]+nx*0.45,s1[1]+ny*0.45),(s2[0]+nx*0.45,s2[1]+ny*0.45),
                      (s2[0]-nx*0.1,s2[1]-ny*0.1),(s1[0]-nx*0.1,s1[1]-ny*0.1)]
                draw.polygon(q2, fill=lt)
            _draw_ball(draw, s1, rad, base)
            _draw_ball(draw, s2, rad, base)

    img = img.resize((W, H), Image.LANCZOS)
    # small caption, drawn after the downscale so it stays crisp
    d2 = ImageDraw.Draw(img)
    try:
        font = ImageFont.load_default()
    except Exception:
        font = None
    d2.text((12, H - 22), "Melon3D - XPBD rigid body physics", fill=(150, 158, 172), font=font)
    return img

def _draw_ball(draw, sc, rad, base):
    # radial gradient faked with concentric circles drifting toward the light
    cx, cy = sc
    # dark rim first: adjacent round bodies then read as occlusion instead of
    # merging into one blob (a resting sphere next to a capsule looked like
    # interpenetration without it)
    rim = tuple(int(v * 0.35) for v in base)
    draw.ellipse([cx - rad - SS, cy - rad - SS, cx + rad + SS, cy + rad + SS], fill=rim)
    steps = 14
    for i in range(steps):
        t = i / (steps - 1)
        r = rad * (1 - 0.60 * t)
        ox = LX * rad * 0.32 * t
        oy = LY * rad * 0.32 * t
        d = AMBIENT + (1 - AMBIENT) * (0.20 + 0.80 * t)
        col = tuple(min(255, int(v * d)) for v in base)
        draw.ellipse([cx+ox-r, cy+oy-r, cx+ox+r, cy+oy+r], fill=col)
    # specular dot
    r = rad * 0.13
    ox, oy = LX * rad * 0.40, LY * rad * 0.40
    col = tuple(min(255, int(v * 0.6 + 120)) for v in base)
    draw.ellipse([cx+ox-r, cy+oy-r, cx+ox+r, cy+oy+r], fill=col)

if __name__ == '__main__':
    bodies, frames = read_scene('scene.bin')
    if len(sys.argv) > 1 and sys.argv[1] == 'preview':
        for fr in [0, 60, 100, 160, 260, 420]:
            fr = min(fr, len(frames) - 1)
            im = render_frame(bodies, frames[fr], fr)
            im.save(f'preview_{fr}.png')
            print('wrote', f'preview_{fr}.png')
    else:
        every = 3
        imgs = []
        for fr in range(0, len(frames), every):
            imgs.append(render_frame(bodies, frames[fr], fr))
            if fr % 60 == 0: print('frame', fr)

        # One shared 255-color palette (index 255 reserved for transparency).
        # Auto palettes (mediancut/octree) starve the small colorful bodies
        # because 85% of the pixels are dark background, so the palette is
        # built explicitly: dark scene colors sampled from real frames plus
        # a deterministic swatch ramp for every object hue at every shading
        # level. No dithering, so static pixels quantize identically every
        # frame and the inter-frame diff turns them transparent.
        n = len(imgs)
        picks = [imgs[0], imgs[n // 3], imgs[2 * n // 3], imgs[-1]]
        montage = Image.new('RGB', (W * 2, H * 2))
        for k, im in enumerate(picks):
            montage.paste(im, ((k % 2) * W, (k // 2) * H))
        arr = np.array(montage).reshape(-1, 3)
        darks = arr[arr.max(axis=1) < 95]
        dark_img = Image.fromarray(darks.reshape(-1, 1, 3))
        dark_q = dark_img.quantize(colors=84, method=Image.MEDIANCUT)
        colors = dark_q.getpalette()[: 84 * 3]

        swatch = []
        for hue in PALETTE.values():
            for f in (0.88, 0.96, 1.04, 1.12):
                base = np.minimum(255, np.array(hue, dtype=float) * f)
                for d in np.linspace(0.36, 1.0, 5):
                    swatch.append(base * d)
                swatch.append(np.minimum(255, base * 0.6 + 120))  # specular
        swatch.append(np.array([150.0, 158, 172]))                # caption text
        for c in swatch[: 255 - len(colors) // 3]:
            colors.extend(int(v) for v in c)
        colors = colors[: 255 * 3]
        colors.extend([0] * (768 - len(colors)))
        pal = Image.new('P', (1, 1))
        pal.putpalette(colors)
        qs = [im.quantize(palette=pal, dither=Image.Dither.NONE) for im in imgs]
        out = [qs[0]]
        prev = np.array(qs[0])
        for q in qs[1:]:
            cur = np.array(q)
            diff = cur.copy()
            diff[cur == prev] = 255
            im = Image.fromarray(diff, 'P')
            im.putpalette(q.getpalette())
            out.append(im)
            prev = cur
        out[0].save('melon3d.gif', save_all=True, append_images=out[1:],
                    duration=int(1000 * every / 60), loop=0,
                    transparency=255, disposal=1, optimize=False)
        print('wrote melon3d.gif', len(out), 'frames')
