'''docstring'''
import matplotlib.pyplot as plt

class CountsPlotter():

    @staticmethod
    def plot_counts_table(counts_table, export_directory):
        sorted_table = sorted(counts_table.items(),reverse=True)
        x, y = zip(*sorted_table)
        
        graph = plt.bar(x,y, align='center', alpha=0.5)
        plt.xticks(list(range(1,max(x))), list(range(1,max(x))),   fontsize=4.5)

        #plt.show()
        plt.savefig(export_directory, dpi=600)
        plt.clf()
