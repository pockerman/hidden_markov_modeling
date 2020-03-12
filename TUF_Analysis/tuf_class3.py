import pymysql, sqlite3, os, tk, matplotlib.pyplot as plt, numpy as np, sys, itertools
from matplotlib import style
from scipy.stats import norm
style.use('fivethirtyeight')

print(sys.version)

class TUF_graph:

    print("Class initiating...")

    def __init__(self, sqlite_db, *args): # args = sqlite3 tables to input
        #connecting to sqlite3 database from python
        con = sqlite3.connect(sqlite_db)
        cur = con.cursor()
        self.tables = []
        self.table_names = [arg for arg in args]
        # extracting data from each table for plotting
        for arg in args: 
            cur.execute(f"""select position from {arg} order by position asc limit 1""")
            start = cur.fetchall()
            cur.execute(f"""select position from {arg} order by position desc limit 1""")
            end = cur.fetchall()
            cur.execute(f"""select position, read_depth, poisson_transformation, 
            negative_binomial from {arg} where position between {start[0][0]} and {end[0][0]}""")
            self.tables.append(cur.fetchall())

    print("Class initiated, data loaded into python")

# if plotting m585 and m605, check which one has higher max_var.

            
    def rd_plot(self, start, end, display, output): # start and end = positions from the sqlite tables
        table_rds = []
        x = 0 # x is an index for the table_names list.
        print(start, end)
        for table, table_name in zip(self.tables, self.table_names):
            rd = []
            rd.append(self.table_names[x])
            x += 1
            for row in table:
                if row[0] >= start and row[0] <= end:
                    #print(row, table_name)
                    if row[1] != 0:
                        rd.append(row[1])
            table_rds.append(rd)

        cat_chroms = [read for chrom in table_rds for read in chrom[1:]] # concatenating the reads for m585 and m605 datasets.
        upper_limit = np.percentile(cat_chroms, 99.99) # returns the 9.99369th percentile of both datasets combined, e.g. the value corresponding
                                                         # to four standard deviations from the mean to use for plotting range.
        lower_limit = np.percentile(cat_chroms, 0.001)

        print("rd function ", upper_limit, lower_limit)

        clr = 0
        #graph settings
        bar_colours = ["r", "b", "g"]
        line_colours = ["crimson", "deepskyblue"] 
        for rd_tab in table_rds:
            #plotting the rd
            plt.hist(rd_tab[1:], bins=50 ,range=(lower_limit, upper_limit), density=True,
            alpha=0.6, color=bar_colours[clr], label=f"{rd_tab[0]}\n mu = {np.round(np.mean(rd_tab[1:]), 2)}, std = {np.round(np.std(rd_tab[1:]), 2)}") # density=True - shows the PDF 
            #plotting the pdf
            a = np.linspace(lower_limit, upper_limit, (lower_limit)) # set hist range/linspace end as a variable == highest rd value for specified region
            pdf = norm.pdf(a, np.mean(rd_tab[1:]), np.std(rd_tab[1:]))
            plt.plot(a, pdf, line_colours[clr], linewidth=2)
            clr += 1
        plt.ylabel("Frequency")
        plt.xlabel("Coverage")
        plt.legend()
        plt.title(f"read_depth\n Chr1:{start}_{end}")
        #plt.title("Mean read depth per 100bp window")
        path = f"TUF_plots/Chr1_{start}:{end}"  #setting path to write graphs to later
        if output == True:
            if os.path.exists(path):
               plt.savefig(f"TUF_plots/Chr1_{start}:{end}/rd_hist.png", bbox_inches='tight')
            else:
                os.mkdir(path)
                plt.savefig(f"TUF_plots/Chr1_{start}:{end}/rd_hist.png", bbox_inches='tight')
        if display == True:
            plt.show()
        plt.close()

        return

    def pt_plot(self, start, end, display, output):

        table_pts = []
        x = 0 # x is an index for the table_names list.
        print(start, end)
        for table, table_name in zip(self.tables, self.table_names):
            pt = []
            pt.append(self.table_names[x])
            x += 1
            for row in table:
                if row[0] >= start and row[0] <= end:
                    #print(row, table_name)
                    if row[2] > 0.1:
                        pt.append(row[2])
            table_pts.append(pt)

        cat_chroms = [read for chrom in table_pts for read in chrom[1:]] # concatenating the reads for m585 and m605 datasets.
        upper_limit = np.percentile(cat_chroms, 99.9999) # returns the 99.99th percentile of both datasets combined, e.g. the value corresponding
                                                       # to five standard deviations from the mean to use for plotting range.
        lower_limit = np.percentile(cat_chroms, 0.00001)

        print("pt function ", upper_limit, lower_limit)

        clr = 0
        #graph settings
        bar_colours = ["r", "b", "g"]
        line_colours = ["crimson", "deepskyblue"] 
        for pt_tab in table_pts:
            #plotting the rd
            plt.hist(pt_tab[1:], bins=50 ,range=(lower_limit, upper_limit), density=True, alpha=0.6, 
                     color=bar_colours[clr], label=f"{pt_tab[0]}\n mu = {np.round(np.mean(pt_tab[1:]), 2)}, std = {np.round(np.std(pt_tab[1:]), 2)}") # density=True - shows the PDF 
            #plotting the pdf
            a = np.linspace(lower_limit, upper_limit, (lower_limit)) # set hist range/linspace end as a variable == highest rd value for specified region
            pdf = norm.pdf(a, np.mean(pt_tab[1:]), np.std(pt_tab[1:]))
            plt.plot(a, pdf, line_colours[clr], linewidth=2)
            clr += 1
        plt.ylabel("Frequency")
        plt.xlabel("coverage")
        plt.legend()
        plt.title(f"poisson_transformation\n Chr1:{start}_{end}")
        path = f"TUF_plots/Chr1_{start}:{end}"  #setting path to write graphs to later
        if output == True:
            if os.path.exists(path):
               plt.savefig(f"TUF_plots/Chr1_{start}:{end}/pt_hist.png", bbox_inches='tight') 
            else:
                os.mkdir(path)
                plt.savefig(f"TUF_plots/Chr1_{start}:{end}/pt_hist.png", bbox_inches='tight')
        if display == True:
            plt.show()
        plt.close()
        return

    def nb_plot(self, start, end, display, output):

        table_nbs = []
        x = 0 # x is an index for the table_names list.
        print(start, end)
        for table, table_name in zip(self.tables, self.table_names):
            nb = []
            nb.append(self.table_names[x])
            x += 1
            for row in table:
                if row[0] >= start and row[0] <= end:
                    #print(row, table_name)
                    if row[3] > 0.1:
                        nb.append(row[3])
            table_nbs.append(nb)

        cat_chroms = [read for chrom in table_nbs for read in chrom[1:]] # concatenating the reads for m585 and m605 datasets.
        upper_limit = np.percentile(cat_chroms, 99.9999) # returns the 99.99th percentile of both datasets combined, e.g. the value corresponding
                                                       # to five standard deviations from the mean to use for plotting range.
        lower_limit = np.percentile(cat_chroms, 0.00001)

        print("nb function ", upper_limit, lower_limit)

        clr = 0
        #graph settings
        bar_colours = ["r", "b", "g"]
        line_colours = ["crimson", "deepskyblue"] 
        for nb_tab in table_nbs:
            #plotting the rd
            plt.hist(nb_tab[1:], bins=50 ,range=(lower_limit, upper_limit), density=True, alpha=0.6, color=bar_colours[clr], label=f"{nb_tab[0]}\n mu = {np.round(np.mean(nb_tab[1:]), 2)}, std = {np.round(np.std(nb_tab[1:]), 2)}") # density=True - shows the PDF 
            #plotting the pdf
            a = np.linspace(lower_limit, upper_limit, (lower_limit)) # set hist range/linspace end as a variable == highest rd value for specified region
            pdf = norm.pdf(a, np.mean(nb_tab[1:]), np.std(nb_tab[1:]))
            plt.plot(a, pdf, line_colours[clr], linewidth=2)
            clr += 1
        plt.ylabel("Frequency")
        plt.xlabel("coverage")
        plt.legend()
        plt.title(f"negative_binomial\n Chr1:{start}_{end}")
        path = f"TUF_plots/Chr1_{start}:{end}"  #setting path to write graphs to later
        if output == True:
            if os.path.exists(path):
               plt.savefig(f"TUF_plots/Chr1_{start}:{end}/nb_hist.png", bbox_inches='tight')
            else:
                os.mkdir(path)
                plt.savefig(f"TUF_plots/Chr1_{start}:{end}/nb_hist.png", bbox_inches='tight')
        if display == True:
            plt.show()
        plt.close()
        return