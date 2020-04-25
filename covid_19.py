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
# SEARCH = ('code3', '840')  # NYC
# SEARCH = ('Province_State', 'New York')


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
    """Returns the basic reproduction ratio as well as gamma and beta
    """
    population = float(time_series['population'])
    infected = time_series['infected'].astype(float)
    susceptible = time_series['susceptible'].astype(float)
    ds_dt, _, dr_dt = derivatives(time_series)
    gamma = dr_dt/infected
    beta = -population*ds_dt/infected/susceptible
    return beta/gamma, beta, gamma


def gather_data(data, search):
    """Returns data as SIR values

    Args:
        data (dict): Dictionary of data from #get_covid_data()
        search (tuple): Tuple of strings to search ('header', 'value')

    Returns:
        (dict) of 'times', 'susceptible', 'infected', 'removed'
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
    time_series['mtimes'] = np.array(time_series['mtimes'],
                                     dtype=int)
    return time_series


def derivatives(time_series):
    """Returns the derivatives with respect to time
    """
    ds_dt = np.gradient(time_series['susceptible'], time_series['mtimes'])
    di_dt = np.gradient(time_series['infected'], time_series['mtimes'])
    dr_dt = -ds_dt-di_dt
    return ds_dt, di_dt, dr_dt


def make_sir_plot(axis, time_series, **kwargs):
    """Makes the SIR plot
    """
    population = time_series['population']
    mtimes = time_series['mtimes']
    s_color = kwargs.get('s_color', 'grey')
    i_color = kwargs.get('i_color', 'blue')
    r_color = kwargs.get('r_color', 'orange')
    axis.plot_date(mtimes, time_series['susceptible'],
                   color=s_color,
                   label=f'susceptible ({_pretty_num_str(population)})',
                   **kwargs)
    axis.plot_date(mtimes, time_series['infected'],
                   color=i_color, label='infected', **kwargs)
    axis.plot_date(mtimes, time_series['removed'],
                   color=r_color, label='removed', **kwargs)
    axis.fill_between(mtimes,
                      time_series['infected'],
                      time_series['susceptible'],
                      facecolor=s_color, alpha=0.2)
    axis.fill_between(mtimes,
                      time_series['infected'],
                      time_series['removed'],
                      facecolor=i_color, alpha=0.5)
    axis.fill_between(mtimes,
                      0,
                      time_series['removed'],
                      facecolor=r_color, alpha=0.5)
    axis.set_ylim((0, 1.25*max(time_series['infected'])))
    axis.set_xlim((mtimes[0], mtimes[-1]))
    axis.set_title(SEARCH[0] + ': ' + SEARCH[1])
    axis.legend()

def make_derivative_plot(axis, time_series, **kwargs):
    """Makes plots of dS/dt dI/dt dR/dt
    """
    ds_dt, di_dt, dr_dt = derivatives(time_series)
    mtimes = time_series['mtimes']
    s_color = kwargs.get('s_color', 'grey')
    i_color = kwargs.get('i_color', 'blue')
    r_color = kwargs.get('r_color', 'orange')
    axis.plot_date(mtimes, ds_dt,
                   color=s_color,
                   label=r'$\frac{dS}{dt}$',
                   **kwargs)
    axis.plot_date(mtimes, di_dt,
                   color=i_color,
                   label=r'$\frac{dI}{dt}$',
                   **kwargs)
    axis.plot_date(mtimes, dr_dt,
                   color=r_color,
                   label=r'$\frac{dR}{dt}$',
                   **kwargs)
    axis.legend()

def make_r0_plot(axis, time_series, **kwargs):
    """Makes basic reproduction ratio plot
    """
    mtimes = time_series['mtimes']
    ratio, *_ = basic_reproduction_ratio(time_series)

    poly_y = ratio[np.argwhere(np.isfinite(ratio))].flatten()
    poly_x = mtimes[np.argwhere(np.isfinite(ratio))].flatten()
    coeff, y_interc = np.polyfit(np.log(poly_x), poly_y, 1)

    r0_color = kwargs.get('r0_color', 'green')
    r0_fit_color = kwargs.get('r0_fit_color', 'black')
    axis.plot_date(mtimes, ratio,
                   color=r0_color,
                   label=r'$R_0$',
                   **kwargs)
    axis.plot_date(poly_x, y_interc+coeff*np.log(poly_x),
                   label=f'fit ({str(int(coeff))} ' + r'$\log(t)$)',
                   color=r0_fit_color,
                   marker='',
                   linestyle='-',
                   **kwargs)

    axis.set_title(r'Basic reproduction ratio $R_0$')
    axis.set_xlim((poly_x[0], poly_x[-1]))
    axis.legend()

def make_gamma_beta_plot(axis, time_series, **kwargs):
    """Make plots for gamma and beta
    """
    _, gamma, beta = basic_reproduction_ratio(time_series)
    mtimes = time_series['mtimes']
    axis.plot_date(mtimes, gamma, label=r'$\gamma$', **kwargs)
    axis.plot_date(mtimes, beta, label=r'$\beta$', **kwargs)
    axis.legend()

def make_plots(data, search=SEARCH):
    """Makes plots of COVID-19 data.

    Args:
        data (dict): Dictionary of data from #get_covid_data()
        search (tuple): Tuple of strings to search ('header', 'value')
    """

    _, axis = plt.subplots(2, 2, sharex='col')
    time_series = gather_data(data, search)

    make_sir_plot(axis[0][0], time_series)
    make_derivative_plot(axis[1][0], time_series)
    make_r0_plot(axis[0][1], time_series)
    make_gamma_beta_plot(axis[1][1], time_series)


    # figure.show()
    plt.show()
    plt.tight_layout()


# HIDDEN FUNCTIONS
def _nearest(pivot, items):
    """Returns nearest point
    """
    return min(items, key=lambda x: abs(x - pivot))


def _parse_url_csv(url):
    """Parse the url csv file into a dict
    """
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


def _pretty_num_str(number, length=1):
    """Returns a string of the number in a prettier way.

    Args:
        number (int, float): Number you want to prettify into string.
        length (int, default=1): Number of decimal places. 0 means round to
                                 integer.

    Examples:
        ```python
        _pretty_num_str(1000, length=0)
        # 1k
        ```
    """
    num_str = str(number)
    suffixes = {
        2: (1, '', ''),
        4: (1000, 'k', 'thousand'),
        7: (1e6, 'M', 'million'),
        10: (1e9, 'B', 'billion'),
        }
    nearest_place = _nearest(len(num_str), suffixes)
    suffix = suffixes[nearest_place]
    return str(round(number/suffix[0], length)) + suffix[1]


# MAIN
if __name__ == '__main__':
    DATA = get_covid_data()
    make_plots(DATA)
