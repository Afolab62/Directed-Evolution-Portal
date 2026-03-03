import numpy as np
import matplotlib.pyplot as plt

# only rows where generation and activity_score exist utilised

df = df.dropna(subset=['Directed_Evolution_Generation', 'activity_score'])

# generations sorted in order

gens = sorted(df['Directed_Evolution_Generation'].unique())

# activity scores per generation extraction
data = []
for g in gens:
    v = df[df['Directed_Evolution_Generation'] == g]['activity_score'].values.astype(float)
    data.append(v)

# figure
fig, ax = plt.subplots(figsize=(10, 7))

# violin plots

violins = ax.violinplot(
    data,
    widths = 0.9,
    points = 200,
)

# whiskers and mean/median activity score per generation

for i, v in enumerate(data, start = 1):

    vmin = np.min(v)
    vmax = np.max(v)
    mean = np.mean(v)
    median = np.median(v)

    ax.vlines(i, vmin, vmax, linewidth = 2)
    ax.hlines(vmin, i - 0.12, i + 0.12, linewidth = 2)
    ax.hlines(vmax, i - 0.12, i + 0.12, linewidth = 2)
    ax.hlines(mean, i - 0.18, i + 0.18, linewidth = 3)
    ax.hlines(median, i - 0.18, i + 0.18, linewidth = 3)

# plot labels
ax.set_xticks(range(1, len(gens) + 1))
ax.set_xticklabels([f'{g}' for g in gens])

ax.set_xlabel('Generation')
ax.set_ylabel('Activity Score')
ax.set_title('Activity Score Distribution by Generation')
ax.grid(True)

plt.tight_layout()
plt.show()