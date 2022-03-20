import time, os
from modules.dbmodel import *
import datetime
from datetime import timedelta
#math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import squarify    # pip install squarify (algorithm for treemap)
import calmap
#special library for imports
import importlib
dbmodel = importlib.import_module('modules.dbmodel')


def facts_to_str(user_data):
    facts = list()

    for key, value in user_data.items():
        facts.append('{} - {}'.format(key, value))

    return "\n".join(facts).join(['\n', '\n'])

#doughnut plot
def draw_douhnut_plot(user_id):
    plt.clf()
    #DATA
    time_netto_min = dbmodel.netto_time_today(user_id)
    time_brutto_min = 24*60
    wasted_time = time_brutto_min - time_netto_min
    categories_list = dbmodel.list_categories(user_id)
    ts = time.localtime()
    ts = time.strftime("%Y-%m-%d", ts)
    # print(ts)
    # print (categories_list, "\n", time_netto_min)

    sizes = []
    sizes.append(wasted_time)
    labels = []
    labels.append("time wasted, {}h.{}min.".format(int(wasted_time//60), int(wasted_time%60)))
    sizes.append(time_netto_min)
    labels.append("time netto, {}min.".format(int(time_netto_min)))

    labels = tuple(labels)

    explode = (0, 0.1)
    fig, ax = plt.subplots(figsize=(15, 10), subplot_kw=dict(aspect="equal"))

    wedges, texts = ax.pie(sizes, explode=explode, wedgeprops=dict(width=0.5), shadow=True, startangle=-40)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, fontsize='30', **kw)

    ax.set_title("потраченное время, \n{}".format(ts), fontsize=25)

    # generation of name
    save_name = "doughnut_{}_{}".format(user_id, ts)

    #saving picture
    plt.savefig(os.path.join('./src', save_name + '.png'), bbox_inches='tight', transparent=True, pad_inches=0.1)
    plt.close()
    return True

def draw_barplot(user_id):

    #DATA
    time_netto_min = dbmodel.netto_time_today(user_id)
    time_brutto_min = 24*60
    wasted_time = time_brutto_min - time_netto_min
    categories_list = dbmodel.list_categories(user_id)
    ts = time.localtime()
    ts = time.strftime("%Y-%m-%d", ts)
    # print(ts)
    # print (categories_list, "\n", time_netto_min)

    if time_netto_min == 0:
        return None
    else:
        sizes = []
        sizes.append(wasted_time)
        labels = []
        labels.append("time brutto, {}".format(time_brutto_min))
        for category in categories_list:
            sizes.append(dbmodel.total_time_per_category_today(user_id, category)["total_time"])
            labels.append(dbmodel.total_time_per_category_today(user_id, category)["category"])

        labels = tuple(labels)

        r = [0]
        # raw_data = {'greenBars': [20], 'orangeBars': [5],'blueBars': [2]}

        raw_data = {}
        raw_data[ts] = []
        for category in categories_list:
            raw_data[ts].append(dbmodel.total_time_per_category_today(user_id, category)["total_time"])

        def survey(results, category_names):
            """
            Parameters
            ----------
            results : dict
                A mapping from question labels to a list of answers per category.
                It is assumed all lists contain the same number of entries and that
                it matches the length of *category_names*.
            category_names : list of str
                The category labels.
            """
            labels = list(results.keys())
            data = np.array(list(results.values()))
            data_cum = data.cumsum(axis=1)
            category_colors = plt.get_cmap('twilight_shifted')(
                np.linspace(0.15, 0.85, data.shape[1]))

            fig, ax = plt.subplots(figsize=(15, 3))
            ax.invert_yaxis()
            ax.xaxis.set_visible(False)
            ax.set_xlim(0, np.sum(data, axis=1).max())

            for i, (colname, color) in enumerate(zip(category_names, category_colors)):
                widths = data[:, i]
                starts = data_cum[:, i] - widths
                ax.barh(labels, widths, left=starts, height=0.5,
                        label=colname, color=color)
                xcenters = starts + widths / 2

                r, g, b, _ = color
                text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
                for y, (x, c) in enumerate(zip(xcenters, widths)):
                    if c != 0:
                        ax.text(x, y, str(int(c)), ha='center', va='center',
                                color=text_color, fontsize='15')
            ax.legend(ncol=1, bbox_to_anchor=(1, 1), loc="upper left", fontsize='15')
            txt = ax.yaxis.get_ticklabels()[-1]
            txt.set_fontstretch('condensed')
            txt.set_fontsize(20)
            return fig, ax

        survey(raw_data, categories_list)
        # generation of name
        save_name = "barplot_{}_{}".format(user_id, ts)
        #saving picture
        plt.savefig(os.path.join('./src', save_name + '.png'), bbox_inches='tight', transparent=True, pad_inches=0.1)
        plt.close()
        return True

def draw_barplot_by_date(user_id, date):
    #DATA
    time_netto_min = dbmodel.netto_time_by_date(user_id, date)
    print(time_netto_min)
    time_brutto_min = 24*60
    wasted_time = time_brutto_min - time_netto_min
    categories_list = dbmodel.list_categories(user_id)
    #date format "%Y-%m-%d"

    if time_netto_min == 0:
        return None
    else:
        raw_data = {}
        raw_data[date] = []
        for category in categories_list:
            raw_data[date].append(dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"])
        print(raw_data)
        def survey(results, category_names):
            """
            Parameters
            ----------
            results : dict
                A mapping from question labels to a list of answers per category.
                It is assumed all lists contain the same number of entries and that
                it matches the length of *category_names*.
            category_names : list of str
                The category labels.
            """
            labels = list(results.keys())
            data = np.array(list(results.values()))
            data_cum = data.cumsum(axis=1)
            category_colors = plt.get_cmap('twilight_shifted')(
                np.linspace(0.15, 0.85, data.shape[1]))

            fig, ax = plt.subplots(figsize=(15, 3))
            ax.invert_yaxis()
            ax.xaxis.set_visible(False)
            ax.set_xlim(0, np.sum(data, axis=1).max())

            for i, (colname, color) in enumerate(zip(category_names, category_colors)):
                widths = data[:, i]
                starts = data_cum[:, i] - widths
                ax.barh(labels, widths, left=starts, height=0.5,
                        label=colname, color=color)
                xcenters = starts + widths / 2

                r, g, b, _ = color
                text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
                for y, (x, c) in enumerate(zip(xcenters, widths)):
                    if c != 0:
                        ax.text(x, y, str(int(c)), ha='center', va='center',
                                color=text_color, fontsize='15')
            ax.legend(ncol=1, bbox_to_anchor=(1, 1), loc="upper left", fontsize='15')
            txt = ax.yaxis.get_ticklabels()[-1]
            txt.set_fontstretch('condensed')
            txt.set_fontsize(20)
            return fig, ax

        survey(raw_data, categories_list)
        # generation of name
        save_name = "barplot_{}_{}".format(user_id, date)
        #saving picture
        plt.savefig(os.path.join('./src', save_name + '.png'), bbox_inches='tight', transparent=True, pad_inches=0.1)
        plt.close()
        return True

def draw_barplot_by_date_list(user_id, date_list, figuresize):
    raw_data = {}
    for date in date_list:
        categories_list = dbmodel.list_categories(user_id)
        #date format "%Y-%m-%d"

        raw_data[date] = []
        for category in categories_list:
            raw_data[date].append(dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"])
    print(raw_data)
    def survey(results, category_names):
        """
        Parameters
        ----------
        results : dict
            A mapping from question labels to a list of answers per category.
            It is assumed all lists contain the same number of entries and that
            it matches the length of *category_names*.
        category_names : list of str
            The category labels.
        """
        labels = list(results.keys())
        data = np.array(list(results.values()))
        data_cum = data.cumsum(axis=1)
        category_colors = plt.get_cmap('twilight_shifted')(
            np.linspace(0.15, 0.85, data.shape[1]))
        if figuresize == "small":
            fig, ax = plt.subplots(figsize=(15, 7))
        else:
            fig, ax = plt.subplots(figsize=(15, 20))
        ax.invert_yaxis()
        ax.xaxis.set_visible(True)
        ax.set_xlim(0, np.sum(data, axis=1).max())

        for i, (colname, color) in enumerate(zip(category_names, category_colors)):
            widths = data[:, i]
            starts = data_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.85,
                    label=colname, color=color)
            xcenters = starts + widths / 2

            r, g, b, _ = color
            text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
            for y, (x, c) in enumerate(zip(xcenters, widths)):
                if c != 0:
                    ax.text(x, y, str(int(c)), ha='center', va='center',
                            color=text_color, fontsize='15')
        # ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
        #           loc='lower left', fontsize='15')
        ax.legend(ncol=1, bbox_to_anchor=(1, 1), loc="upper left", fontsize='15')
        return fig, ax

    survey(raw_data, categories_list)
    # generation of name
    save_name = "barplot_{}_week".format(user_id)
    #saving picture
    plt.savefig(os.path.join('./src', save_name + '.png'), bbox_inches='tight', transparent=True, pad_inches=0.1)
    plt.close()
    return True

def treemap_chart_by_date(user_id, date):

    #DATA
    time_netto_min = dbmodel.netto_time_by_date(user_id, date)
    time_brutto_min = 24*60
    wasted_time = time_brutto_min - time_netto_min
    categories_list = dbmodel.list_categories(user_id)
    #date format "%Y-%m-%d"

    if time_netto_min == 0:
        return None
    else:
        raw_data = {}
        raw_data[date] = []
        for category in categories_list:
            if dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"] != 0:
                raw_data[date].append(dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"])
    labels = []
    for category in categories_list:
        if dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"] != 0:
            labels.append("{}\n{:.0f} min.".format(category, dbmodel.total_time_per_category_by_date(user_id, category, date)["total_time"]))

    data = np.array(categories_list)

    #choosing colours
    category_colors = [plt.cm.Spectral(i/float(data.shape[0])) for i in range(data.shape[0])]
    # category_colors = plt.get_cmap('rainbow_r')(np.linspace(0.15, 0.85, data.shape[0]))

    squarify.plot(sizes=raw_data[date], label=labels, alpha=.7, color=category_colors) #, text_kwargs={'fontsize':15}
    plt.axis('off')

    # generation of name
    save_name = "treemap_{}_{}".format(user_id, date)
    #saving picture
    plt.savefig(os.path.join('./src', save_name + '.png'), bbox_inches='tight', transparent=True, pad_inches=0.1)
    plt.close()
    return True

def get_date_format(number):
    today = datetime.date.today()
    date = today - timedelta(days=int(number))
    date = date.strftime('%Y-%m-%d')
    return date

def draw_year_heatmap(user_id, category):

    date = str(datetime.date.today().replace(day=2))

    past_thirty_days = pd.date_range(end=datetime.date.today(), periods=30, freq='D')
    current_month = datetime.date.today().month
    current_month_days = pd.date_range(start=datetime.date.today().replace(day=1), periods=30, freq='D')
    current_year_days = pd.date_range(start=datetime.date.today().replace(day=1, month=1), end=datetime.date.today().replace(day=31, month=12), freq='D')

    date_list = list(current_year_days)
    for i in range(len(date_list)):
        date_list[i] = date_list[i].to_pydatetime().strftime("%Y-%m-%d")

    total_time_per_category_per_day_list =[]
    for date in date_list:
        total_time = dbmodel.total_time_per_category_by_date(user_id, category, date)['total_time']
        total_time_per_category_per_day_list.append(total_time)

    time_series = np.array(total_time_per_category_per_day_list)
    data = pd.Series(time_series, index=current_year_days)

    calmap.calendarplot(data, daylabels='MTWTFSS', linewidth=5, yearlabel_kws={'color':'black', 'fontsize':20}, dayticks=True, fig_kws=dict(figsize=(25, 4)))

    save_name = "heatmap_{}".format(user_id)
    calmap.plt.savefig(os.path.join('./src', save_name + '.jpg'), bbox_inches='tight', transparent=True, pad_inches=0.1)
    plt.close()
    return True

def build_menu(buttons,
               n_cols,
               header_buttons=None,
               footer_buttons=None):
    menu = [buttons[i:i + n_cols] for i in range(0, len(buttons), n_cols)]
    if header_buttons:
        menu.insert(0, [header_buttons])
    if footer_buttons:
        menu.append([footer_buttons])
    return menu
