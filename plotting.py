import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from datetime import timedelta


def plot_ic(p, name, t0_date):

    # plot the ic occupation and add axes labels
    dates = [t0_date + timedelta(days=i) for i in range(len(p.ic_t))]
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%b"))
    if len(dates) < 300:
        ax.xaxis.set_major_locator(mdates.WeekdayLocator())
    else:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_tick_params(rotation=90)
    ax.set_xlabel("Date")
    ax.set_ylabel("Number of IC beds")
    (ic_oc,) = ax.plot(dates, p.ic_t, color="b", label="required IC capacity")

    # add hzline for ic capacity
    ic_max_kwargs = {"color": "black", "linestyle": "--"}
    ax.axhline(2000, **ic_max_kwargs)
    ax.plot([], [], label="IC capacity", **ic_max_kwargs)

    # plot vertical lines for lockdowns and openups if present
    if p.lockdowns:
        lockdown_kwargs = {"color": "r", "linestyle": "--"}
        for x in p.lockdowns:
            ax.axvline(dates[x], **lockdown_kwargs)
        ax.plot([], [], **lockdown_kwargs, label="start lockdown")

    if p.openups:
        openups_kwargs = {"color": "g", "linestyle": "--"}
        for x in p.openups:
            ax.axvline(dates[x], **openups_kwargs)
        ax.plot([], [], **openups_kwargs, label="end lockdown")

    # add legend
    ax.legend(loc="upper right")
    ax.set_xlim([dates[0], dates[-1]])

    # save fig
    fig.set_figwidth(12)
    fig.set_figheight(10)
    fig.savefig(name)
