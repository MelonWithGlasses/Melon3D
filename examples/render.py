import struct, math, sys
import numpy as np
from PIL import Image, ImageDraw

SS = 2                      # supersample factor
W, H = 760, 560
WS, HS = W * SS, H * SS
FOV = math.radians(46)
FOCAL = (HS / 2) / math.tan(FOV / 2)

EYE = np.array([11.5, 7.0, 15.5])
TARGET = np.array([0.0, 2.2, 0.0])
UP = np.array([0.0, 1.0, 0.0])
LIGHT = np.array([-0.35, 1.0, 0.55]); LIGHT = LIGHT / np.linalg.norm(LIGHT)
AMBIENT = 0.42

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
    1: (78, 190, 130),    # pyramid - melon green
    2: (240, 84, 74),     # projectile - coral
    3: (77, 179, 255),    # rain sphere - sky
    4: (255, 184, 77),    # rain capsule - amber
    5: (168, 120, 255),   # rain box - purple
}

def hue_shift(c, i):
    # small per-index variation so identical bodies read as distinct
    f = 0.85 + 0.30 * ((i * 0.61803398875) % 1.0)
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

def draw_floor(draw):
    # filled floor quad + grid lines at y=0
    ext = 17.0
    quad = []
    ok = True
    for (x, z) in [(-ext,-ext),(ext,-ext),(ext,ext),(-ext,ext)]:
        p, _ = project(to_view(np.array([x, 0.0, z])))
        if p is None: ok = False; break
        quad.append(p)
    if ok:
        draw.polygon(quad, fill=(38, 44, 56))
    step = 2.0
    n = int(ext / step)
    for i in range(-n, n+1):
        c = i * step
        for a, b in [((c,-ext),(c,ext)), ((-ext,c),(ext,c))]:
            pa, _ = project(to_view(np.array([a[0],0.0,a[1]])))
            pb, _ = project(to_view(np.array([b[0],0.0,b[1]])))
            if pa and pb:
                draw.line([pa, pb], fill=(58, 66, 80), width=SS)

def render_frame(bodies, frame, idx):
    # background vertical gradient
    img = Image.new('RGB', (WS, HS))
    top = np.array([26, 30, 40]); bot = np.array([44, 50, 64])
    col = (top[None,:] * (1 - (np.arange(HS)/HS))[:,None] + bot[None,:] * (np.arange(HS)/HS)[:,None]).astype(np.uint8)
    img = Image.fromarray(np.repeat(col[:,None,:], WS, axis=1))
    draw = ImageDraw.Draw(img)
    draw_floor(draw)

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
            draw.polygon(poly, fill=col, outline=(20, 24, 30))
        elif kind == 'sphere':
            sc, rad, base = pl
            _draw_ball(draw, sc, rad, base)
        elif kind == 'capsule':
            s1, s2, rad, base = pl
            _draw_ball(draw, s1, rad, base)
            _draw_ball(draw, s2, rad, base)
            # connecting body: quad between the two circle edges
            dx, dy = s2[0]-s1[0], s2[1]-s1[1]
            L = math.hypot(dx, dy)
            if L > 1e-3:
                nx, ny = -dy/L*rad, dx/L*rad
                quad = [(s1[0]+nx,s1[1]+ny),(s2[0]+nx,s2[1]+ny),(s2[0]-nx,s2[1]-ny),(s1[0]-nx,s1[1]-ny)]
                draw.polygon(quad, fill=shade(base, np.array([0,0.4,0.9])))
    img = img.resize((W, H), Image.LANCZOS)
    return img

def _draw_ball(draw, sc, rad, base):
    # a few concentric offset circles to fake a lit sphere
    cx, cy = sc
    steps = 5
    for i in range(steps):
        t = i / (steps - 1)
        r = rad * (1 - 0.55 * t)
        ox = -rad * 0.32 * t; oy = -rad * 0.32 * t
        d = AMBIENT + (1 - AMBIENT) * (0.25 + 0.75 * t)
        col = tuple(min(255, int(v * d)) for v in base)
        draw.ellipse([cx+ox-r, cy+oy-r, cx+ox+r, cy+oy+r], fill=col)

if __name__ == '__main__':
    bodies, frames = read_scene('scene.bin')
    if len(sys.argv) > 1 and sys.argv[1] == 'preview':
        for fr in [30, 110, 200, 320]:
            im = render_frame(bodies, frames[fr], fr)
            im.save(f'preview_{fr}.png')
            print('wrote', f'preview_{fr}.png')
    else:
        every = 3
        imgs = []
        for fr in range(0, len(frames), every):
            imgs.append(render_frame(bodies, frames[fr], fr))
            if fr % 60 == 0: print('frame', fr)
        imgs[0].save('melon3d.gif', save_all=True, append_images=imgs[1:],
                     duration=int(1000 * every / 60), loop=0, optimize=True)
        print('wrote melon3d.gif', len(imgs), 'frames')
