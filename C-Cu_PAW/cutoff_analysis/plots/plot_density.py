import numpy as np
import matplotlib.pyplot as plt
from ase.io.cube import read_cube_data

# === Carica dati dal file .cube ===
cube_file = "density_to_be_plotted_80_Ry.cube"  # <-- Cambia nome se necessario
data, atoms = read_cube_data(cube_file)

# === Info ===
print("Forma della densità elettronica (Nx, Ny, Nz):", data.shape)
print("Atomi nel sistema:", atoms)

# === Seleziona un piano centrale (taglio lungo Y) ===
y_index = data.shape[1] // 2
slice = data[:, y_index, :]

# === Plot ===
plt.figure(figsize=(6, 5))
plt.imshow(slice.T, origin='lower', cmap='grey')
plt.colorbar(label='Densità elettronica (a.u.)')
plt.title(f'Densità elettronica - slice')
plt.xlabel('x grid index')
plt.ylabel('y grid index')
plt.tight_layout()
plt.show()
