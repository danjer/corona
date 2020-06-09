from datetime import date
from scipy.stats import norm
from plotting import plot_ic


def day_to_p(dist, day):
    return dist.cdf(day) - dist.cdf(day - 1)


class Population:
    def __init__(
        self,
        size=17000000,
        infected=100,
        si_distribution=norm(7.5, 3.5),
        ic_admission_distribution=norm(10, 3.5),
        max_treatment_duration=30,
        r0=2.2,
        ic_p=0.005,
    ):

        # describe the current state
        self.size = size
        self.infected = infected
        self.ic = 0

        # probability distributions used to update the states
        self.si_d = si_distribution
        self.ic_d = ic_admission_distribution
        self.max_treatment_duration = max_treatment_duration

        # estimates for update parameters
        self.ic_p = ic_p
        self.r0 = r0

        # the number of (new) infected and ic patients after updating the state
        self.nw_infected = infected
        self.nw_ic = 0

        # to keep track of the history of the population
        self.infected_t = []
        self.ic_t = []
        self.susceptible_t = []
        self.nw_infected_t = []
        self.nw_ic_t = []

        # store the the lockdown and openup events
        self.lockdowns = []
        self.openups = []

        # write day zero to history
        self.record_state()

    @property
    def susceptible(self):
        return self.size - self.infected

    @property
    def r_effective(self):
        susceptible = self.susceptible_t[-1]
        return self.r0 * (susceptible / self.size)

    def update_infected(self):
        nw_infected = 0
        days_back = 15
        for day in range(1, days_back):
            try:

                # the number of new infected on day(s) back
                n_nwinf = self.nw_infected_t[-day]

                # the current reproduction
                re = self.r_effective

                # the probability that a given transmission occurs after day(s)
                si_p = day_to_p(self.si_d, day)

                # add the contribution to nw infected
                nw_infected += re * si_p * n_nwinf
            except IndexError:
                break

        self.nw_infected = nw_infected
        self.infected += nw_infected

    def update_ic(self):

        # remove dismissed patients
        dismissed_ic = 0
        for day in range(self.max_treatment_duration):
            try:
                index = self.max_treatment_duration - day
                n = self.nw_ic_t[-index]
                dismissed_ic += n * (1 / self.max_treatment_duration)

            except IndexError:
                pass

        # add new ic patients
        days_back = 20
        nw_ic = 0

        for day in range(1, days_back):
            try:

                # the number of new infected on days(s) back
                n_nwinf = self.nw_infected_t[-day]

                # the probability that an infected needs ic
                ic_p = self.ic_p

                # the probability that, given that ic treatment is required, it
                # is required after day(s)
                ic_p2 = day_to_p(self.ic_d, day)

                # add the new ic patients that got infected day days back
                nw_ic += n_nwinf * ic_p * ic_p2
            except IndexError:
                break

        # update state
        balance_ic = nw_ic - dismissed_ic
        self.nw_ic = nw_ic
        self.ic = self.ic_t[-1] + balance_ic
        
    def record_state(self):
        self.infected_t.append(self.infected)
        self.ic_t.append(self.ic)
        self.susceptible_t.append(self.susceptible)
        self.nw_infected_t.append(self.nw_infected)
        self.nw_ic_t.append(self.nw_ic)

    def run(self, days=1):
        for _ in range(days):
            self.update_infected()
            self.update_ic()
            self.record_state()


def main():

    EXIT_STRATEGY = "open_up"
    NUMBER_OF_DAYS = 365

    # initialize population
    p = Population()

    # synchronize based on the number of new ic patients on day 20
    while True:
        p.run()

        if p.nw_ic > 80:  # on 20 march the number of new ic patients was 82

            t0 = -20  # 1th the of march
            t1 = -5  # 15th of march, start of lockdown

            # remove history, keep 1th to 15th of march
            p.infected_t = p.infected_t[t0:t1]
            p.ic_t = p.ic_t[t0:t1]
            p.susceptible_t = p.susceptible_t[t0:t1]
            p.nw_ic_t = p.nw_ic_t[t0:t1]
            p.nw_infected_t = p.nw_infected_t[t0:t1]

            t0_date = date(2020, 3, 1)
            break

    # 15 march inteligent lockdown started..
    p.r0 = 0.9
    p.lockdowns.append(len(p.ic_t))
    p.run(60)

    if EXIT_STRATEGY == None:
        plot_ic(p, "initial_outbreak.png", t0_date)

    if EXIT_STRATEGY == "open_up":
        p.r0 = 2.2
        p.openups.append(len(p.ic_t))
        p.run(NUMBER_OF_DAYS)
        plot_ic(p, "open_up.png", t0_date)

    if EXIT_STRATEGY == "interrupted_braking":
        day = 0
        while day < NUMBER_OF_DAYS:
            if p.nw_ic > 80 and p.r0 != 0.9:
                p.r0 = 0.9
                p.lockdowns.append(len(p.ic_t))
            elif p.ic < 300 and p.r0 != 2.2:
                p.r0 = 2.2
                p.openups.append(len(p.ic_t))
            p.run()
            day += 1
        plot_ic(p, "interrupted_braking.png", t0_date)

    if EXIT_STRATEGY == "titration":
        day = 0
        while day < NUMBER_OF_DAYS:
            if p.ic < 300 and p.r_effective < 1:
                p.r0 += 0.4
            p.run()
            day += 1
        plot_ic(p, "titration.png", t0_date)

    if EXIT_STRATEGY == "extensive_testing":
        p.r0 = 2.2
        day = 0
        while day < NUMBER_OF_DAYS:
            p.nw_infected_t[-3] *= 0.5
            p.run()
            day += 1
        plot_ic(p, "extensive_testing.png", t0_date)


if __name__ == "__main__":
    main()
