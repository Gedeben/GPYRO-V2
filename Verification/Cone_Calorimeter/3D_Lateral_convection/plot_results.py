import glob
import re
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.patches import Rectangle, Polygon
import sys


# Paramètres de style pour l’article
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 22
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['legend.fontsize'] = 22
mpl.rcParams['xtick.labelsize'] = 22
mpl.rcParams['ytick.labelsize'] = 22

# -------- Extraction temps / images depuis fichier .ssf --------
def extract_time_image_pairs(ssf_file):
    with open(ssf_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    time_image_pairs = []
    current_time = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if line == "SETTIMEVAL":
            i += 1
            if i < len(lines):
                try:
                    current_time = int(lines[i])
                except ValueError:
                    current_time = None
        elif line.startswith("RENDERONCE"):
            i += 1
            if i < len(lines):
                filename = f"{lines[i]}.png"
                if current_time is not None:
                    time_image_pairs.append((filename, current_time))
                    current_time = None
        i += 1
    return time_image_pairs

# -------- Lecture images --------
pairs = extract_time_image_pairs("reference_cc_3D_conv_01.ssf")  
images = []
image_times = []

# Load images
for fname, t in pairs:
    try:
        img = Image.open(fname)
        images.append(img)
        image_times.append(t)
    except FileNotFoundError:
        print(f"[!] Image {fname} not found, skipping.")

if not images:
    print("[!] No valid image found. Exiting.")
    sys.exit(1)



# -------- Création image composite --------
width, height = images[0].size
images = [img.resize((width, height)) for img in images]
total_height = height * len(images)
composite = Image.new('RGB', (width, total_height))
for idx, img in enumerate(images):
    composite.paste(img, (0, idx * height))
composite.save("composite_raw.png")

# -------- Création figure matplotlib --------

# -------- Création de la figure --------
fig, ax = plt.subplots(figsize=(15, len(images)*3))
img = plt.imread("composite_raw.png")
ax.imshow(img)


ax.axis('off')

# -------- Ajouter le texte des temps --------
for i, t in enumerate(image_times):
    ypos = (len(image_times)-i - 0.5) / len(images)
    ax.text(1.01, ypos, f"t = {t} s", transform=ax.transAxes, va='center')

# -------- Ajouter une flèche verticale --------




# --- Dégradé vertical ---
x_start = 1.17
width = 0.04
y_start = 0.07
y_end = 0.97
n_steps = 100

# Inverser le dégradé : gris clair (haut) → gris foncé (bas)
grays = np.linspace(0.6, 0.2, n_steps)[::-1]

for i in range(n_steps):
    y0 = y_start + (y_end - y_start) * (i / n_steps)
    y1 = y_start + (y_end - y_start) * ((i + 1) / n_steps)
    height = y1 - y0
    color = str(grays[i])

    rect = Rectangle(
        (x_start, y0), width, height,
        transform=ax.transAxes,
        facecolor=color,
        edgecolor='none',
        clip_on=False,
        zorder=5
    )
    ax.add_patch(rect)

# --- Triangle pour la tête de flèche ---
triangle = Polygon(
    [
        (x_start + width / 2, y_start - 0.06),        # sommet
        (x_start - 0.015, y_start+0.01),                   # coin gauche
        (x_start + width + 0.015, y_start+0.01)            # coin droit
    ],
    closed=True,
    facecolor=str(grays[1]),
    edgecolor='none',
    transform=ax.transAxes,
    clip_on=False,
    zorder=23
)
ax.add_patch(triangle)




# -------- Ajouter le texte "Time" dans la flèche --------
ax.text(
    x_start+width/3, 0.5, "Time",
    rotation=270,
    va='center', ha='center',
    color='white',
    transform=ax.transAxes,
    fontsize=30,
    zorder=11
)


# -------- Couleurs de la colorbar --------
colorbar_colors = [
    (0.000000, 0.000000, 0.000000),
    (0.169118, 0.155025, 0.140931),
    (0.338235, 0.310049, 0.281863),
    (0.507353, 0.465074, 0.422794),
    (0.676471, 0.620098, 0.563725),
    (0.845588, 0.775123, 0.704657),
    (1.000000, 0.913386, 0.000000),
    (1.000000, 0.732283, 0.000000),
    (1.000000, 0.551181, 0.000000),
    (1.000000, 0.370079, 0.000000),
    (1.000000, 0.188976, 0.000000),
    (1.000000, 0.000000, 0.000000)
]

n_top = 250
n_bot = 250
colors_top = colorbar_colors[:6]    # noir à jaune
colors_bot = colorbar_colors[6:]    # jaune à rouge foncé

# -------- Lecture du fichier pour récupérer min et max --------
with open("reference_cc_3D_conv_01_001.sf.bnd", "r") as f:
    line = f.readline().strip()
    parts = list(map(float, line.split()))
    val_min = parts[1]
    val_max = parts[2]


def interpolate_gradient(colors, n):
    gradient = np.linspace(0, 1, n)
    interpolated = np.zeros((n, 1, 3))
    base = np.linspace(0, 1, len(colors))
    for i in range(3):  # R, G, B
        channel = [c[i] for c in colors]
        interpolated[:, 0, i] = np.interp(gradient, base, channel)
    return interpolated

grad_top = interpolate_gradient(colors_top, n_top)
grad_bot = interpolate_gradient(colors_bot, n_bot)

# Inverser le sens pour que rouge soit en haut, noir en bas
full_gradient = np.vstack([grad_bot[::-1], grad_top[::-1]])

# -------- Ajouter la colorbar à gauche --------
fig.subplots_adjust(left=0.22)
cbar_ax = fig.add_axes([0.15, 0.36, 0.04, 0.29])
cbar_ax.imshow(full_gradient, aspect='auto')

# Supprimer les ticks sur l'axe x
cbar_ax.set_xticks([])

# Supprimer le cadre (spines)
for spine in cbar_ax.spines.values():
    spine.set_visible(False)

# Afficher ticks sur l'axe y uniquement
cbar_ax.yaxis.tick_left()

nb_ticks = 9
ticks_pos = np.linspace(0, full_gradient.shape[0]-1, nb_ticks)
ticks_labels = np.linspace(val_min, val_max, nb_ticks)
ticks_labels = [f"{v:.2f}" for v in ticks_labels[::-1]]

cbar_ax.set_yticks(ticks_pos)
cbar_ax.set_yticklabels(ticks_labels)

# Ajouter la légende
cbar_ax.set_ylabel(r"Reaction Rate (Kg.m$^{⁻3}$.s$^{⁻1}$)", rotation=90, labelpad=15)


# -------- Sauvegarde --------


fig.savefig("smokeview_3D_conv.png", dpi=300, bbox_inches='tight')

sys.exit(0)
