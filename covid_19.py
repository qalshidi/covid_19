#!/usr/bin/env python
"""Plotting COVID-19 for my own curiosity
"""
import csv
import datetime as dt
import urllib
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
        (dict) of 'confirmed' and 'deaths'
    """
    data_confirmed = _parse_url_csv(url_confirmed)
    data_deaths = _parse_url_csv(url_deaths)
    return {'confirmed': data_confirmed, 'deaths': data_deaths}


def basic_reproduction_ratio(times, susceptible, infected, removed):
    """Returns the basic reproduction ratio as well as gamma and beta"""
    time = np.array([time.timestamp() for time in times], dtype=float)
    d_t = np.diff(time)
    population = susceptible[0]
    d_s = np.diff(susceptible)
    d_r = np.diff(removed)
    gamma = d_r/d_t/infected[1:]
    beta = -population*d_s/d_t/infected[1:]/susceptible[1:]
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
    for header, local_pop in zip(data['deaths'][search_header],
                                 data['deaths']['Population']):
        if header == search_value:
            population += int(local_pop)

    # build time series
    times = []
    infected = []
    susceptible = []
    deaths = []
    for date in data['confirmed']:
        if date.endswith('20'):
            times.append(dt.datetime.strptime(date, '%m/%d/%y'))

            # count
            num_of_deaths = 0
            for header, local_pop in zip(data['deaths'][search_header],
                                         data['deaths'][date]):
                if header == search_value:
                    num_of_deaths += int(local_pop)
            deaths.append(num_of_deaths)
            num_of_infected = 0
            for header, local_pop in zip(data['confirmed'][search_header],
                                         data['confirmed'][date]):
                if header == search_value:
                    num_of_infected += int(local_pop)
            susceptible.append(population-num_of_infected-num_of_deaths)
            infected.append(num_of_infected)

    return (times,
            np.array(susceptible, dtype=int),
            np.array(infected, dtype=int),
            np.array(deaths, dtype=int))


def make_plots(data, search=SEARCH):
    """Makes plots of COVID-19 data.

    Args:
        data (dict): Dictionary of data from #get_covid_data()
        search (tuple): Tuple of strings to search ('header', 'value')
    """

    _, axis = plt.subplots(2, 1, sharex=True)
    plot_options = {}
    times, susceptible, infected, deaths = gather_data(data, search)
    mtimes = [mdates.date1num(time) for time in times]
    axis[0].plot_date(mtimes, susceptible,
                      color='grey',
                      label=f'susceptible ({str(susceptible[0]/1e6)[:4]})M',
                      **plot_options)
    axis[0].plot_date(mtimes, infected,
                      color='blue', label='infected', **plot_options)
    axis[0].plot_date(mtimes, deaths,
                      color='orange', label='deaths', **plot_options)
    axis[0].fill_between(mtimes, infected, susceptible,
                         facecolor='grey', alpha=0.2)
    axis[0].fill_between(mtimes, infected, deaths, facecolor='blue', alpha=0.5)
    axis[0].fill_between(mtimes, 0, deaths, facecolor='orange', alpha=0.5)
    axis[0].set_ylim((0, 1.25*max(infected)))
    axis[0].set_xlim((mtimes[0], mtimes[-1]))
    axis[0].set_title(SEARCH[0] + ': ' + SEARCH[1])
    axis[0].legend()
    ratio, beta, gamma = basic_reproduction_ratio(times,
                                                  susceptible,
                                                  infected,
                                                  deaths)

    # coeff, const = np.polyfit(mtimes[1:],
    #                           np.log(ratio),
    #                           1,
    #                           w=np.sqrt(ratio))
    # axis[1].plot_date(mtimes[1:], np.exp(const)*np.exp(coeff*mdates[1:]))

    coeff, y_interc = np.polyfit(np.log(mtimes[1:]), ratio, 1)
    axis[1].plot_date(mtimes[1:], y_interc+coeff*np.log(mtimes[1:]),
                      color='black', label=r'$R_0$')

    axis[1].plot_date(mtimes[1:], ratio, color='green', label=r'$R_0$')

    axis[1].set_title(r'Basic reproduction ratio $R_0$')
    axis[1].legend()
    plt.show()
    plt.tight_layout()


if __name__ == '__main__':
    DATA = get_covid_data()
    make_plots(DATA)
