#!/usr/bin/env python
"""Plotting COVID-19 for my own curiosity
"""
import csv
import datetime as dt
import urllib.request
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

URL_BASE = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/'
URL_BASE += 'csse_covid_19_data/csse_covid_19_time_series/'
URL_CONFIRMED = 'time_series_covid19_confirmed_US.csv'
URL_DEATHS = 'time_series_covid19_deaths_US.csv'
SEARCH = ('Province_State', 'Michigan')
# SEARCH = ('Province_State', 'Massachusetts')


def _parse_url_csv(url):
    """Parse the url csv file into a dict"""
    return_data = {}
    with urllib.request.urlopen(url) as datafile:
        lines = [line.decode('ascii') for line in list(datafile)]
        headers = lines[0].strip().split(',')
        reader = csv.DictReader(lines)
        for header in headers:
            return_data[header] = []
        for row in reader:
            for key, value in row.items():
                return_data[key].append(value)
    return return_data


def get_covid_data(url_confirmed=URL_BASE+URL_CONFIRMED,
                   url_deaths=URL_BASE+URL_DEATHS):
    """Returns dictionary of COVID-19 data from John Hopkins University.

    Args:
        url_confirmed (str): url for csv with confirmed cases.
        url_deaths (str): url for csv with death cases.

    Returns:
        (dict) of 'confirmed' and 'removed'
    """
    data_confirmed = _parse_url_csv(url_confirmed)
    data_deaths = _parse_url_csv(url_deaths)
    return {'confirmed': data_confirmed, 'deaths': data_deaths}


def basic_reproduction_ratio(time_series):
    """Returns the basic reproduction ratio as well as gamma and beta"""
    d_t = np.diff(time_series['mtimes'])
    population = float(time_series['population'])
    d_s = np.diff(time_series['susceptible'])
    ds_dt = d_s/d_t
    d_i = np.diff(time_series['infected'])
    di_dt = d_i/d_t
    dr_dt = -ds_dt-di_dt
    infected = time_series['infected'][1:].astype(float)
    susceptible = time_series['susceptible'][1:].astype(float)
    gamma = dr_dt/infected
    beta = -population*d_s/d_t/infected/susceptible
    return beta/gamma, beta, gamma


def gather_data(data, search):
    """Returns data as SIR values

    Args:
        data (dict): Dictionary of data from #get_covid_data()
        search (tuple): Tuple of strings to search ('header', 'value')

    Returns:
        (list(datetime.datetime), np.array, np.array, np.array) of
        times, susceptible, infected, removed
    """
    search_header, search_value = search[0], search[1]

    population = 0
    local_pop = (num for header_value, num in zip(data['deaths'][search_header],
                                                  data['deaths']['Population'])
                 if header_value == search_value)
    for num in local_pop:
        population += int(num)

    # build time series
    time_series = {
        'times': [],
        'mtimes': [],
        'infected': [],
        'susceptible': [],
        'deaths': [],  # figure this out from removed
        'removed': [],  # figure this out from deaths
        'population': population
        }

    # Go through days
    dates = (date for date in data['confirmed']
             if date.endswith('20') or date.endswith('19'))
    for date in dates:

        day = dt.datetime.strptime(date, '%m/%d/%y')
        time_series['times'].append(day)
        time_series['mtimes'].append(mdates.date2num(day))

        # count
        val_deaths = (num for header_value, num in zip(
            data['deaths'][search_header], data['deaths'][date])
                      if header_value == search_value)
        num_of_deaths = 0
        for num in val_deaths:
            num_of_deaths += int(num)
        time_series['deaths'].append(num_of_deaths)
        time_series['removed'].append(num_of_deaths)

        val_infected = (num for header_value, num in zip(
            data['confirmed'][search_header], data['confirmed'][date])
                        if header_value == search_value)
        num_of_infected = 0
        for num in val_infected:
            num_of_infected += int(num)
        time_series['susceptible'].append(population-num_of_infected-num_of_deaths)
        time_series['infected'].append(num_of_infected)

    time_series['susceptible'] = np.array(time_series['susceptible'],
                                          dtype=int)
    time_series['infected'] = np.array(time_series['infected'],
                                       dtype=int)
    time_series['removed'] = np.array(time_series['removed'],
                                      dtype=int)
    return time_series


def make_plots(data, search=SEARCH):
    """Makes plots of COVID-19 data.

    Args:
        data (dict): Dictionary of data from #get_covid_data()
        search (tuple): Tuple of strings to search ('header', 'value')
    """

    figure, axis = plt.subplots(2, 1, sharex=True)
    plot_options = {}
    time_series = gather_data(data, search)
    mtimes = [mdates.date2num(time) for time in time_series['times']]
    population = time_series['population']
    axis[0].plot_date(mtimes, time_series['susceptible'],
                      color='grey',
                      label=f'susceptible ({str(population/1e6)[:4]})M',
                      **plot_options)
    axis[0].plot_date(mtimes, time_series['infected'],
                      color='blue', label='infected', **plot_options)
    axis[0].plot_date(mtimes, time_series['removed'],
                      color='orange', label='removed', **plot_options)
    axis[0].fill_between(mtimes,
                         time_series['infected'],
                         time_series['susceptible'],
                         facecolor='grey', alpha=0.2)
    axis[0].fill_between(mtimes,
                         time_series['infected'],
                         time_series['removed'],
                         facecolor='blue', alpha=0.5)
    axis[0].fill_between(mtimes,
                         0,
                         time_series['removed'],
                         facecolor='orange', alpha=0.5)
    axis[0].set_ylim((0, 1.25*max(time_series['infected'])))
    axis[0].set_xlim((mtimes[0], mtimes[-1]))
    axis[0].set_title(SEARCH[0] + ': ' + SEARCH[1])
    axis[0].legend()
    ratio, beta, gamma = basic_reproduction_ratio(time_series)

    axis[1].plot_date(mtimes[1:], ratio, color='green', label=r'$R_0$')

    axis[1].set_title(r'Basic reproduction ratio $R_0$')
    axis[1].legend()
    figure.show()
    # plt.show(figure)
    # plt.tight_layout()


if __name__ == '__main__':
    DATA = get_covid_data()
    make_plots(DATA)
