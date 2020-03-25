#%%
from meta_analysis import main

def test_plot():
    result_model = main.csv_meta_analysis('../test_data/test.csv',)
    result_model.plot_forest()

test_plot()