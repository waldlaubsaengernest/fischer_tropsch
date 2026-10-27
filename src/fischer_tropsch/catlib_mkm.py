import numpy as np
import math
from copy import deepcopy

def reaction_names():
     return ["CO_oxidation","NO_reduction","FT"]

def mkm_class(reaction_name):
    if reaction_name =="ft": from fischer_tropsch.ft import FT; return FT
    if reaction_name =="CO_oxidation": from coo import COO; return COO
    if reaction_name =="NO_reduction": from nor import NOR; return NOR

class Reaction:
    def __init__(self,a,b,ab,reverse,c=None,d=None,nu=None):
        self.a  = a
        self.b  = b
        self.ab = ab
        self.nu = nu  # dictionary of stoichiometric coefficients
        self.reverse = reverse

        if self.ab is not None:
            reactants = [self.species_string(ab)]
            products  = [self.species_string(a), self.species_string(b)]
        else:
            self.c = c
            self.d = d
            reactants = [self.species_string(c), self.species_string(d)]
            products  = [self.species_string(a), self.species_string(b)]

        if self.reverse:
            reactants, products = products, reactants

        self.name = "+".join(reactants) + "->" + "+".join(products)

        self.er = None
        self.ea = None
        self.k  = None
        self.k_cap = None

    def species_string(self,species):
        if species is None:
            return None
        coeff = ""
        if self.nu is not None and self.nu.get(species, 1) > 1:
            coeff = str(self.nu[species])
        return coeff + species

    def get_er_ea(self,reverse):
        if self.er is None or self.ea is None: return 0,0
        if reverse:
            er_new = -1*self.er
            ea_new = max([0,self.ea - self.er,-self.er])
        else:
            er_new = self.er
            ea_new = max([0,self.ea,self.er])
        return er_new,ea_new


class MKM:
    def __init__(self,eads,temperature=None,gas_pressures=None,steady_state=True,sites=None,nguess=200,max_cost=1e-8,initial_guess=None,optcyc=2000,adsorption_steady_state=False,cov_dep_eads=False,cov_ads="CO",min_rates=None):
        self.eads_database = eads
        self.eads          = eads

        self.cov_dep_eads  = cov_dep_eads
        self.cov_ads       = None
        self.eads_fit      = None

        if self.cov_dep_eads:
            self.cov_ads = cov_ads + "*"
            self.eads_fit = {}

            # Collect all adsorbates contained in the database.
            adsorbates = set()
            for coverage in self.eads_database:
                adsorbates.update(self.eads_database[coverage])

            # Build one coverage fit for each adsorbate.
            for adsorbate in adsorbates:
                x = []
                y = []

                for coverage in sorted(self.eads_database):
                    if adsorbate in self.eads_database[coverage]:
                        x.append(coverage)
                        y.append(self.eads_database[coverage][adsorbate])

                if len(x) >= 2:
                    self.eads_fit[adsorbate] = piecewise_polyfit_curve(x, y)

            if initial_guess is not None:
                theta = initial_guess
            else:
                theta = {self.cov_ads: 0.0}

            self.update_coverage_eads(theta)

        self.temperature   = temperature
        self.get_gas_pressures(gas_pressures)
        self.steady_state  = steady_state
        self.sites         = sites
        self.nguess        = nguess
        self.max_cost      = max_cost
        self.initial_guess = initial_guess
        self.optcyc        = optcyc
        self.min_rates     = {} if min_rates is None else min_rates
        self.adsorption_steady_state = adsorption_steady_state
        self.constants     = {"R" : 8.314*10**(-3)/96.485, "h":4.1356677 *10**(-15)}
        self.initialize_reactions()
        if self.temperature is not None and self.gas_pressures is not None:
            self.run_model()

    def update_coverage_eads(self, theta):
        coverage = theta[self.cov_ads]

        self.eads = {
            adsorbate: real_eads(eads_fit, coverage)
            for adsorbate, eads_fit in self.eads_fit.items()
        }

    def initialize_reactions(self):
        self.reactions = []
        self.read_reaction_energies()
    
    def read_reaction_energies(self):
        self.reactions = []
    
    def get_gas_pressures(self,gas_pressures):
        self.gas_pressures = gas_pressures
    
    def get_main_barriers(self):
        return {}
    
    def get_catmap_adsorption_energies(self):
        return {}
    
    def prepare_mkm(self):
        self.adsorbates            = []
        self.in_steady_state       = []
        self.include_initial_guess = []

    def get_final_rates_and_selectivities(self):
        return
    
    def sabotier_analysis(self):
        from catlib_addons import Rate_FT
        rate1, p_C2H4       = Rate_FT([reac.k[0] for reac in self.reactions],[reac.k[1] for reac in self.reactions],[reac.k_cap for reac in self.reactions],self.gas_pressures["CO"],self.gas_pressures["H2"],self.gas_pressures["H2O"],self.gas_pressures["CH4"],RDS=3)
        self.reaction_rates = [rate1]
        self.theta          = None
        self.selectivity    = None
    
    def run_model(self):
        self.get_reaction_k()
        steady_state = True
        if not steady_state:
            self.sabotier_analysis()
        else:
            self.prepare_mkm()
            if self.adsorption_steady_state: self.in_steady_state = deepcopy(self.adsorbates)
            theta, theta_star, sol = self.solve_mkm()
            self.reaction_rates = compute_reaction_rates(self.reactions, theta, self.gas_pressures, self.sites)

            self.get_final_rates_and_selectivities()
            self.theta = theta
            for ads in self.ads_ignore:
                self.theta[ads] = theta_fixed(ads,self.reactions,self.gas_pressures)
            self.get_print_strings()

    def get_reaction_k(self):
        r = self.constants["R"]
        h = self.constants["h"]
        t = self.temperature
        for reaction in self.reactions:
            er_reverse, ea_reverse = reaction.get_er_ea(reverse=True)
            er_normal, ea_normal   = reaction.get_er_ea(reverse=False)
            ea = [ea_normal,ea_reverse]
            k = [r*t/h * np.exp(-ea[i]/(r*t)) for i in [0,1]]
            reaction.k = k
            eps = 10**(-14)
            if k[0] == 0 and k[1]==0: reaction.k_cap=0
            else: reaction.k_cap = reaction.k[0]/max(eps,reaction.k[1])

    def solve_mkm(
        self,
        enforce_nonzero=5,
    ):
        adsorbates = self.adsorbates
        steady_state = self.in_steady_state
        ads_ignore = self.ads_ignore
        gas_pressures = self.gas_pressures
        sites = self.sites
        include_initial_guess = self.include_initial_guess

        # Build once to determine the order of the unknown coverages.
        _, unknowns = build_residual(
            self.reactions,
            adsorbates,
            steady_state,
            ads_ignore,
            gas_pressures,
            sites,
            self.adsorption_steady_state,
        )

        if self.cov_dep_eads:
            def residual(theta_vec):
                # theta changes during every least_squares iteration.
                theta = dict(zip(unknowns, theta_vec))

                # Update all coverage-dependent adsorption energies.
                self.update_coverage_eads(theta)

                # Rebuild reactions and rate constants for the current theta.
                self.initialize_reactions()
                self.get_reaction_k()

                # initialize_reactions() replaces self.reactions, so always use
                # the current self.reactions rather than a stale local list.
                current_residual, current_unknowns = build_residual(
                    self.reactions,
                    adsorbates,
                    steady_state,
                    ads_ignore,
                    gas_pressures,
                    sites,
                    self.adsorption_steady_state,
                )

                if current_unknowns != unknowns:
                    raise RuntimeError(
                        "Coverage-dependent reaction update changed the MKM unknowns."
                    )

                return current_residual(theta_vec)
        else:
            residual, unknowns = build_residual(
                self.reactions,
                adsorbates,
                steady_state,
                ads_ignore,
                gas_pressures,
                sites,
                self.adsorption_steady_state,
            )

        N = len(unknowns)

        lb = np.zeros(N)
        ub = np.ones(N)

        best = None
        best_cost = np.inf

        def solution_species_rate(sol, species):
            # min_rates is defined through the surface formation rate of the
            # corresponding adsorbed species, matching surface_species_rate().
            #
            # Example:
            #     min_rates={"CH4": 1e-8}
            # checks the net surface formation rate of "CH4*".
            #
            # Adsorption/desorption reactions are deliberately ignored.
            if self.cov_dep_eads:
                residual(sol.x)

            surface_species = species if "*" in species else species + "*"
            theta_sol = dict(zip(unknowns, sol.x))

            reaction_rates = compute_reaction_rates(
                self.reactions,
                theta_sol,
                gas_pressures,
                sites,
            )

            total = 0.0

            for r in self.reactions:
                if is_adsorption_reaction(r):
                    continue

                net = reaction_rates[r.name]["net"]

                if r.ab is not None:
                    reactants = expand_species([r.a, r.b], r.nu)
                    products = expand_species([r.ab], r.nu)
                else:
                    reactants = expand_species([r.a, r.b], r.nu)
                    products = expand_species([r.c, r.d], r.nu)

                if not r.reverse:
                    reactants, products = products, reactants

                nu = (
                    products.count(surface_species)
                    - reactants.count(surface_species)
                )
                total += nu * net

            return total

        def acceptable(sol):
            if sol is None or not sol.success:
                return False

            for species_group, min_rate in self.min_rates.items():
                if isinstance(species_group, str):
                    species_group = (species_group,)

                # Tuple/list means OR: at least one species in the group must
                # have a sufficiently large positive net production rate.
                if not any(
                    solution_species_rate(sol, species) >= min_rate
                    for species in species_group
                ):
                    return False

            return True

        if include_initial_guess is not None:
            theta_d = {s: 0 for s in unknowns}
            initial_guesses = []

            for i in range(self.nguess):
                for s in include_initial_guess:
                    theta_d[s] = np.random.rand()
                initial_guesses.append(
                    np.array([theta_d.get(s, 0) for s in unknowns])
                )

            if self.initial_guess is not None:
                theta_ig = np.array(
                    [self.initial_guess.get(s, 0) for s in unknowns]
                )

                sol = least_squares(
                    residual,
                    theta_ig,
                    # bounds=(lb, ub),
                    xtol=1e-12,
                    ftol=1e-12,
                    gtol=1e-12,
                    max_nfev=self.optcyc,
                    verbose=0,
                )

                if acceptable(sol) and sol.cost < self.max_cost:
                    best_cost = sol.cost
                    theta = dict(zip(unknowns, sol.x))

                    # Make sure self.reactions and reaction.k correspond exactly
                    # to the converged coverage before returning.
                    if self.cov_dep_eads:
                        residual(sol.x)

                    theta_star = get_theta_star(theta, sites)
                    return theta, theta_star, sol

            cycle = 0

            for theta0 in initial_guesses:
                sol = least_squares(
                    residual,
                    theta0,
                    # bounds=(lb, ub),
                    xtol=1e-12,
                    ftol=1e-12,
                    gtol=1e-12,
                    max_nfev=self.optcyc,
                    verbose=0,
                )

                if acceptable(sol) and sol.cost < best_cost:
                    best = sol
                    best_cost = sol.cost
                    print(best_cost, flush=True)

                    if best_cost < self.max_cost:
                        break

                cycle += 1

            print(best_cost, flush=True)

            if best is None:
                raise RuntimeError(
                    f"No acceptable MKM solution found after {self.nguess} guesses."
                )

            sol = best

        else:
            sol = least_squares(
                residual,
                theta0,
                # bounds=(lb, ub),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                max_nfev=2000,
                verbose=2,
            )

        theta = dict(zip(unknowns, sol.x))

        # Make sure the final reaction energies and rate constants correspond
        # exactly to the returned coverage.
        if self.cov_dep_eads:
            residual(sol.x)

        theta_star = get_theta_star(theta, sites)

        return theta, theta_star, sol

    def surface_species_rate(self,species):
        total = 0.0
        for r in self.reactions:
            # Ignore adsorption/desorption reactions.
            if is_adsorption_reaction(r):
                continue

            net = self.reaction_rates[r.name]["net"]

            if r.ab is not None:
                reactants = expand_species([r.a, r.b], r.nu)
                products  = expand_species([r.ab], r.nu)
            else:
                reactants = expand_species([r.a, r.b], r.nu)
                products  = expand_species([r.c, r.d], r.nu)

            if not r.reverse:
                reactants, products = products, reactants

            nu = products.count(species) - reactants.count(species)
            total += nu * net

        return total

    def get_print_strings(self):
        self.final_rates = dict()
        for rate in self.reaction_rates:
            if not "+" in rate: self.final_rates[rate] = self.reaction_rates[rate]
        sorted_theta = dict(sorted(self.theta.items(), key=lambda item: item[1],reverse=True))
        self.print_head = "   ".join(["N"] + ["theta_"+pt for pt in sorted_theta]+ ["rate("+r+")" for r in self.final_rates]) + " selectivity\n"
        self.print_line =" ".join([str(format_zahl(x)) for x in ([self.theta[tn] for tn in sorted_theta]+
            [self.final_rates[r] for r in self.final_rates] +
            [self.selectivity*100])])+"\n"

def to_dict(reactions):
    reac = dict()
    if reactions is None: return None
    for reaction in reactions:
        reac[reaction.name] = reaction
    return reac

def to_list(reactions):
    reac = []
    if reactions is None: return None
    for name in reactions:
        reac.append(reactions[name])
    return reac


from scipy.optimize import least_squares

def get_theta_star(theta,sites):
    for a in theta:
        if theta[a]>1:
            theta[a]=1
        if theta[a]<0:
            theta[a]=0

    if sites is not None:
            theta_star = dict()
            for s in sites:
                theta_star[s] = max(0,1.0-sum(theta[a] for a in sites[s] if a in theta))
    else:
        theta_star = max(0, 1.0 - sum(theta.values()))
    return theta_star

def build_residual(
    reactions,
    adsorbates,
    steady_state,
    ads_ignore,
    gas_pressures,
    sites,
    adsorption_steady_state=False
):
    unknowns = [
        s for s in adsorbates
        if s not in ads_ignore and is_surface(s)
    ]

    def residual(theta_vec):
        unknowns = [s for s in adsorbates if s not in ads_ignore and is_surface(s)]

        theta = dict(zip(unknowns, theta_vec))
        #if theta["CO*"] > 0.95:
        #    print(theta["CO*"])

        theta_star = get_theta_star(theta,sites)

        res = []
        eq_ss =[]

        for r in reactions:
            if not is_adsorption_reaction(r) or adsorption_steady_state:
                continue
            
            # check if at least one reactant/product is an unknown
            involved = [r.a, r.b, r.ab]
            if r.ab is None:
                involved = [r.a, r.b, r.c,r.d]
            unk2 = [s for s in adsorbates if s not in ads_ignore and s not in steady_state and is_surface(s)]

            unknown_involved = any(s in unk2 for s in involved if is_surface(s))
            
            if not unknown_involved:
                continue  # skip, alles steady_state oder ads_ignore
            #(r.a,r.b,r.ab)
            rf, rr = rate_forward_backward(r, theta, gas_pressures,reactions,sites)
            res.append(rf - rr)
            #print(r.name,rf-rr)
        #print("f")

        # B) Steady State nur für reine Oberflächenreaktionen
        for ads in steady_state:
            eq = 0.0
            for r in reactions:
                if is_adsorption_reaction(r) and not adsorption_steady_state:
                    continue
                nu = stoich_coeff(ads, r)
                if nu == 0:
                    continue
                rf, rr = rate_forward_backward(r, theta, gas_pressures,reactions,sites)
                eq += nu * (rf - rr)
                #print(r.name,rf-rr)
            res.append(eq)
            eq_ss.append(eq)

        #print(len(res),len(eq_ss),len(unknowns),unknowns)
        return np.array(res)

    return residual, unknowns

def is_gas(s):
    return "*" not in s

def is_surface(s):
    return "*" in s and s != "*"

def n_sites(s):
    if s == "*":
        return 1
    if "*" in s:
        return s.count("*")
    return 0

def is_adsorption_reaction(r):
    if r.ab is not None:
        return any(is_gas(s) for s in (r.a, r.b, r.ab))
    else:
        return any(is_gas(s) for s in (r.a, r.b, r.c,r.d))

def expand_species(species,nu):
    expanded = []
    for s in species:
        coeff = 1 if nu is None else nu.get(s, 1)
        expanded.extend([s] * coeff)
    return expanded

def stoich_coeff(ads, r):
    if not r.reverse:
        products = [r.a, r.b]
        reactants  = [r.ab]
        if r.ab is None:
            reactants  = [r.c,r.d]
            products   = [r.a,r.b]
    else:
        products   = [r.ab]
        reactants  = [r.a, r.b]
        if r.ab is None:
            reactants  = [r.a,r.b]
            products   = [r.c,r.d]

    reactants = expand_species(reactants,r.nu)
    products = expand_species(products,r.nu)

    return products.count(ads) - reactants.count(ads)

def get_site_of_star(reactants,products,sites):
    count = dict()
    for r in reactants:
        for site in sites:
            if site not in count:
                count[site]=0
            if r in sites[site]:
                count[site]-=1
    for r in products:
        for site in sites:
            if r in sites[site]:
                count[site]+=1
    sites_not_zero = [site for site in count if count[site] != 0]
    #if len(sites_not_zero)!=1:
    return {site:count[site] for site in sites_not_zero}
    #return sites_not_zero[0]

def rate_forward_backward(r, theta, gas_pressures,reactions,sites):
    
    theta_star = get_theta_star(theta,sites)

    def val(s):
        if s == "*":
            if sites is None:
                return theta_star
            else:
                return 1
        if is_gas(s):
            return gas_pressures[s]
        if s not in theta:
            return theta_fixed(s,reactions,gas_pressures)
        return theta[s]

    if not r.reverse:
        products = [r.a, r.b]
        reactants  = [r.ab]
        if r.ab is None: reactants  = [r.c,r.d]
        kf, kr = r.k
    else:
        products = [r.ab]
        reactants  = [r.a, r.b]
        if r.ab is None: products  = [r.c,r.d]
        kf, kr = r.k

    reactants = expand_species(reactants,r.nu)
    products  = expand_species(products,r.nu)

    rf = kf
    for s in reactants:
        rf *= val(s)

    rr = kr
    for s in products:
        rr *= val(s)

    site_diff = (
        sum(n_sites(p) for p in products)
        - sum(n_sites(r_) for r_ in reactants)
    )

    #if r.ab=="C2H2*" and r.a=="C2H2" and r.b=="*":
    #    print(products,theta_star,rf,rr,rf-rr,site_diff)

    if sites is None:
        if site_diff < 0:
            rr *= theta_star ** abs(site_diff)
        elif site_diff > 0:
            rf *= theta_star ** site_diff

    else:
        site_count = get_site_of_star(reactants,products,sites)
        if site_count!={}:
            theta_star_ed = math.prod(theta_star[s]**site_count[s] for s in site_count if site_count[s]>0)
            theta_star_prod   = math.prod(theta_star[s]**abs(site_count[s]) for s in site_count if site_count[s]<0)
            #if round(theta_star_ed,1)!=1 or round(theta_star_prod,1)!=1:
                #print("check")
            rr *= theta_star_prod
            rf *= theta_star_ed

    return rf, rr

def get_zeros(res):
    zeros = 0
    for r in res:
        if r==0:
            zeros+=1
    return zeros

def reaction_string(r):
    return r.name

def compute_reaction_rates(reactions, theta, gas_pressures,sites):
    reaction_rates = {}

    for r in reactions:
        rf, rr = rate_forward_backward(r, theta, gas_pressures,reactions,sites)

        key = r.name

        reaction_rates[key] = {
            "forward": rf,
            "backward": rr,
            "net": rf - rr
        }
    return reaction_rates

def theta_fixed(species,reactions,gas_pressures):
    reactions = to_dict(reactions)
    if species=="H*":
        kp_05 = np.sqrt(reactions["H2->H*+H*"].k_cap * gas_pressures["H2"])
        return kp_05/(1+kp_05)
    if species=="CO*":
        k = reactions["CO+*->CO*"].k_cap
        p = gas_pressures["CO"]
        return k*p/(1+k*p)
    if species=="H2O*":
        k = 1/reactions["H2O*->H2O+*"].k_cap
        p = gas_pressures["H2O"]
        return k*p/(1+k*p)


def format_zahl(n):
    if n == 0:
        return "0.00"
    if abs(n) < 0.01:
        return f"{n:.1e}"
    else:
        return f"{n:.2f}"
    

import numpy as np


def piecewise_polyfit_curve(x, y, n_points=300):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one-dimensional.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) < 2:
        raise ValueError("At least two finite points are required.")

    sort_index = np.argsort(x)
    x = x[sort_index]
    y = y[sort_index]

    unique_x = np.unique(x)

    if len(unique_x) < 2:
        raise ValueError("At least two different x values are required.")

    if len(unique_x) != len(x):
        y = np.array([np.mean(y[x == xv]) for xv in unique_x], dtype=float)
        x = unique_x

    total_width = float(x[-1] - x[0])
    if total_width <= 0.0:
        raise ValueError("x values must span a finite interval.")

    x_parts = []
    y_parts = []
    coefficients = []

    for i in range(len(x) - 1):
        x0 = float(x[i])
        x1 = float(x[i + 1])
        y0 = float(y[i])
        y1 = float(y[i + 1])

        coefficients_i = np.polyfit([x0, x1], [y0, y1], 1)
        coefficients.append(coefficients_i)

        interval_fraction = (x1 - x0) / total_width
        n_interval = max(2, int(round(n_points * interval_fraction)))

        x_interval = np.linspace(x0, x1, n_interval)
        y_interval = np.polyval(coefficients_i, x_interval)

        # Avoid duplicate points at interval boundaries.
        if i > 0:
            x_interval = x_interval[1:]
            y_interval = y_interval[1:]

        x_parts.append(x_interval)
        y_parts.append(y_interval)

    return (
        np.concatenate(x_parts),
        np.concatenate(y_parts),
        coefficients,
        1,
    )

def real_eads(eads_fit, theta):
    """
    Return the adsorption energy at a given coverage theta.
    """
    x_fit, y_fit, _, _ = eads_fit

    return float(np.interp(theta, x_fit, y_fit))

def get_eads_from_database(
    gads_database,
    co_coverage,
    o_coverage=0.0,
):
    """
    Build eads[adsorbate] for the requested coverage.

    Supported database formats:

        old:
            database[co_coverage][adsorbate]

        new:
            database[(co_coverage, o_coverage)][adsorbate]
    """

    def split_coverage(coverage):
        if isinstance(coverage, tuple):
            if len(coverage) != 2:
                raise ValueError(
                    f"Invalid coverage key: {coverage!r}"
                )
            return float(coverage[0]), float(coverage[1])

        # Old database format: float means CO coverage with O coverage = 0.
        return float(coverage), 0.0

    adsorbates = set()

    for coverage, values in gads_database.items():
        _, current_o_cov = split_coverage(coverage)

        if abs(current_o_cov - o_coverage) < 1e-12:
            adsorbates.update(values)

    eads = {}

    for adsorbate in adsorbates:
        x = []
        y = []

        for coverage in sorted(
            gads_database,
            key=lambda cov: split_coverage(cov),
        ):
            current_co_cov, current_o_cov = split_coverage(coverage)

            if abs(current_o_cov - o_coverage) >= 1e-12:
                continue

            if adsorbate not in gads_database[coverage]:
                continue

            x.append(current_co_cov)
            y.append(gads_database[coverage][adsorbate])

        if len(x) >= 2:
            fit = piecewise_polyfit_curve(x, y)
            eads[adsorbate] = real_eads(fit, co_coverage)

        elif len(x) == 1:
            eads[adsorbate] = float(y[0])

    return eads
