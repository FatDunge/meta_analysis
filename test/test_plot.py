#%%
from meta_analysis import main

def test_plot():
    result_model = main.csv_meta_analysis('../test_data/test.csv',)
    result_model.plot_forest()

test_plot()

# %%
import matplotlib.pyplot as plt
dpi = 100
width = 14
height = 5
grid_width = 1 / width
grid_height = 1 / height
fig, ax = plt.subplots(figsize=(width, height),
                        dpi=dpi)
# set default to invisible
ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.hlines(1, 0, 1)
ax.hlines(0, 0, 1)
ax.hlines(1/14, 0, 1)

plt.show()

# %%
