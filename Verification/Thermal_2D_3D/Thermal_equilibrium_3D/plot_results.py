import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob
import re
from PIL import Image
from matplotlib.patches import Rectangle, Polygon
import sys


import os


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))



#%% color and style parameters for the plots 
param_sets = [
    { 'color':'blue', 'linewidth': 3, 'linestyle': '-', "label" :'Point 1'},
    { 'color':'red' ,'linewidth': 3, 'linestyle': '-', "label" :"Point 2"},
    { 'color':'green','linewidth': 3, 'linestyle': '-', "label" : "Point 3"},
    ]

#%%


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 26  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)



#%%

Ntime=-1
# Load data
try:
    data= pd.read_csv("thermal_equilibrium_summary_01_0001.csv")
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)



#%%





DIM_TOTAL=0.20 #20 cm*20 cm *20 cm
DIM_HOT= 0.10 # 10 cm *10 cm *10 cm
THOT= 673.15 # in K
TCOLD= 293.15 # in K

RHO_HOT= 1000
CP_HOT = 1000

RHO_COLD= 700
CP_COLD = 600


VTOT=DIM_TOTAL**3
VHOT=DIM_HOT**3
VCOLD=VTOT-VHOT

H_HOT= CP_HOT*RHO_HOT*VHOT*THOT
H_COLD= CP_COLD*RHO_COLD*VCOLD*TCOLD

H_TOTAL= H_HOT+H_COLD

T_final = data['001_TEMPERATURE( 0.0500_ 0.0500_ 0.0500)'].values[-1]


#H_TOTAL= (CP_HOT*RHO_HOT*VHOT+CP_COLD*RHO_COLD*VCOLD)*T_EQ

T_EQ= H_TOTAL/(CP_HOT*RHO_HOT*VHOT+CP_COLD*RHO_COLD*VCOLD) -273.15




plt.figure()
plt.axhline(y=T_EQ, color='k', lw=3,linestyle='--', label=f'Analytical = {T_EQ:.2f} °C' )

plt.plot(data['t'].values,data['001_TEMPERATURE( 0.0500_ 0.0500_ 0.0500)'].values, **param_sets[0])
plt.plot(data['t'].values,data['002_TEMPERATURE( 0.0500_ 0.1500_ 0.1500)'].values, **param_sets[1])
plt.plot(data['t'].values,data['003_TEMPERATURE( 0.1300_ 0.0500_ 0.0500)'].values, **param_sets[2])


plt.xlabel('Time (s)')
plt.ylabel("Temperature [°C]")
plt.legend(frameon=False)
plt.tick_params()
plt.tight_layout()
plt.savefig("point_temperature_3D_Teq", dpi=300, bbox_inches='tight')  # Enregistrement de la figure



#%%



# Function to compute absolute error between simulation and analytical solution
def compute_absolute_error_uniform_grid(sim_depths, sim_temps, ana_depths, ana_temps, num_points=100):
    # Define a uniform grid over the common depth range
    min_depth = max(min(sim_depths), min(ana_depths))
    max_depth = min(max(sim_depths), max(ana_depths))
    uniform_depths = np.linspace(min_depth, max_depth, num_points)

    # Interpolate both datasets on the uniform grid
    sim_interp = np.interp(uniform_depths, sim_depths, sim_temps)
    ana_interp = np.interp(uniform_depths, ana_depths, ana_temps)

    # Compute sum of absolute differences
    error = np.mean(np.abs(sim_interp - ana_interp))
    return error

# Compute error
error_value = abs(T_final- T_EQ)






#%%



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
pairs = extract_time_image_pairs("thermal_equilibrium_01.ssf")  
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

if images :  
        
    # -------- Création image composite --------
    width, height = images[0].size
    images = [img.resize((width, height)) for img in images]
    total_width  = width * len(images)
    composite = Image.new('RGB', (total_width, height))
    for idx, img in enumerate(images):
        composite.paste(img, ( idx * width,0))
    composite.save("composite_raw.png")
    
    
    # -------- Création figure matplotlib --------
    
    # -------- Création de la figure --------
    fig, ax = plt.subplots(figsize=(len(images)*5, 5))
    img = plt.imread("composite_raw.png")
    ax.imshow(img)
    
    
    ax.axis('off')
    
    
    # -------- Ajouter le texte des temps --------
    for i, t in enumerate(image_times):
        xpos = (i + 0.5) / len(images)
        ax.text(
            xpos, 1.05, f"t = {t} s",
            transform=ax.transAxes,
            va='center',  # vertical centering
            ha='center'   # horizontal centering
        )
    
    
    # -------- Couleurs de la colorbar --------
    colorbar_colors = [
     (0.186347, 0.198469, 0.970456),
     (0.186347, 0.198469, 0.970456),
     (0.282280, 0.397470, 0.907141),
     (0.274480, 0.569688, 0.837588),
     (0.114252, 0.727861, 0.765207),
     (0.222277, 0.838418, 0.572903),
     (0.194393, 0.938675, 0.329380),
     (0.400266, 0.944289, 0.000000),
     (0.706114, 0.790856, 0.000000),
     (0.858515, 0.606036, 0.000000),
     (0.962672, 0.335775, 0.000000),
     (0.962672, 0.335775, 0.000000)
    ]
    
    
    
    n_color_bar=10
    
    
    # -------- Lecture du fichier pour récupérer min et max --------
    with open("thermal_equilibrium_01_001.sf.bnd", "r") as f:
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
    
    gradp = interpolate_gradient(colorbar_colors, n_color_bar)
    
    grad=gradp[::-1]
    # Inverser le sens pour que rouge soit en haut, noir en bas
    
    
    
    # -------- Ajouter la colorbar à gauche --------
    fig.subplots_adjust(left=0.22)
    cbar_ax = fig.add_axes([0.18, 0.17, 0.03, 0.65])
    cbar_ax.imshow(grad, aspect='auto')
    
    
    nb_ticks = 5
    ticks_pos = np.linspace(0, grad.shape[0]-1, nb_ticks)
    ticks_labels = np.linspace(val_min, val_max, nb_ticks)
    ticks_labels = [f"{v:.0f}" for v in ticks_labels[::-1]]
    
    cbar_ax.set_yticks(ticks_pos)
    cbar_ax.set_yticklabels(ticks_labels)
    
    # Supprimer les ticks sur l'axe x
    cbar_ax.set_xticks([])
    
    # Supprimer le cadre (spines)
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    
    # Afficher ticks sur l'axe y uniquement
    cbar_ax.yaxis.tick_left()
    
    
    # -------- Ajouter axes Y et X en bas à gauche --------
    arrow_len = 0.1   # Longueur relative des flèches (en fraction de la largeur/hauteur)
    offset_x = 0.05   # Position horizontale de l'origine
    offset_y = 0.15   # Position verticale de l'origine
    
    # Flèche Y (horizontale)
    ax.annotate("", xy=(offset_x + arrow_len, offset_y), xytext=(offset_x, offset_y),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=3, color='red'))
    ax.text(offset_x + arrow_len + 0.015, offset_y, "y", transform=ax.transAxes,
            va='center', ha='left', fontsize=30, color='red')
    
    # Flèche X (verticale)
    ax.annotate("", xy=(offset_x, offset_y + 3.5*arrow_len), xytext=(offset_x, offset_y),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=3, color='red'))
    ax.text(offset_x, offset_y + 3.5*arrow_len + 0.015, "z", transform=ax.transAxes,
            va='bottom', ha='center', fontsize=30, color='red')
    
    ax.axis('off')
    
    
    
    # -------- Sauvegarde --------
    
    
    fig.savefig("slice_temperature_3D_Teq.png", dpi=300, bbox_inches='tight')


#%%

# Define a threshold for validation
threshold = 0.3  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error = {error_value:.6f} °C")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
