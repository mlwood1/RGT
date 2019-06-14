'''docstring'''
import matplotlib.pyplot as plt

class CountsPlotter():

    @static
    def plot_counts_table(counts_table, export_directory):
        sorted_table = sorted(counts_table.items(),reverse=True)
        x, y = zip(*sorted_table)
        graph = plt.plot(x,y)
        plt.savefig(export_directory)
        plt.clf()