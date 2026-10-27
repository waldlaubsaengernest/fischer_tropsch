import numpy as np
import math
import os
import signal
import socket
from copy import deepcopy
from fischer_tropsch.catlib import format_tabelle

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
    def __init__(self,eads,temperature=None,gas_pressures=None,steady_state=True,sites=None,nguess=200,max_cost=1e-8,initial_guess=None,optcyc=2000,adsorption_steady_state=False,cov_dep_eads=False,cov_ads="CO",parallel=False,min_rates=None):
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
        self.parallel      = parallel
        self.min_rates = {}
        if min_rates is not None:
            items = min_rates.items() if hasattr(min_rates, "items") else min_rates
            for species_group, min_rate in items:
                if isinstance(species_group, list):
                    species_group = tuple(species_group)
                self.min_rates[species_group] = min_rate
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
        enforce_nonzero=None,
    ):
        adsorbates = self.adsorbates
        steady_state = self.in_steady_state
        ads_ignore = self.ads_ignore
        gas_pressures = self.gas_pressures
        sites = self.sites
        include_initial_guess = self.include_initial_guess

        if self.parallel:
            try:
                from mpi4py import MPI
            except ImportError as exc:
                raise ImportError(
                    "parallel=True requires mpi4py. Install/load mpi4py or use parallel=False."
                ) from exc
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
        else:
            comm = None
            rank = 0
            size = 1

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
                theta = dict(zip(unknowns, theta_vec))
                self.update_coverage_eads(theta)
                self.initialize_reactions()
                self.get_reaction_k()

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

            surface_exists = any(surface_species in (r.a, r.b, r.ab, getattr(r, "c", None), getattr(r, "d", None)) for r in self.reactions)

            total = 0.0

            for r in self.reactions:
                if is_adsorption_reaction(r) and surface_exists:
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

                if surface_exists:
                    nu = ( products.count(surface_species) - reactants.count(surface_species))
                else:
                    nu = ( products.count(species) - reactants.count(species))
                total += nu * net

            return total

        def acceptable(sol):
            if sol is None or not sol.success:
                return False

            # Optional legacy coverage filter.  Disabled by default: a physically
            # valid steady state may legitimately contain many zero coverages.
            if enforce_nonzero is not None:
                zeros = get_zeros(sol.x)
                if zeros >= len(sol.x) - enforce_nonzero:
                    return False

            # Physically useful steady state: requested product/species formation
            # rates must be at least their user-specified minimum values.
            #
            # A tuple/list of species is interpreted as an OR condition:
            # at least one species in the group must reach min_rate.
            for species_group, min_rate in self.min_rates.items():
                if isinstance(species_group, str):
                    species_group = (species_group,)

                if not any(
                    solution_species_rate(sol, species) >= min_rate
                    for species in species_group
                ):
                    return False

            return True

        def make_random_guess():
            theta_d = {s: 0.0 for s in unknowns}
            for s in include_initial_guess:
                theta_d[s] = np.random.rand()
            return np.array([theta_d.get(s, 0.0) for s in unknowns])

        # ------------------------------------------------------------------
        # SERIAL PATH: preserve the old behaviour as closely as possible.
        # ------------------------------------------------------------------
        if not self.parallel or size == 1:
            best = None
            best_cost = np.inf

            if self.initial_guess is not None:
                theta_ig = np.array(
                    [self.initial_guess.get(s, 0) for s in unknowns]
                )
                sol = least_squares(
                    residual,
                    theta_ig,
                    xtol=1e-12,
                    ftol=1e-12,
                    gtol=1e-12,
                    max_nfev=self.optcyc,
                    verbose=0,
                )
                if acceptable(sol) and sol.cost < self.max_cost:
                    best = sol
                    best_cost = sol.cost

            if best is None and include_initial_guess is not None:
                for _ in range(self.nguess):
                    theta0 = make_random_guess()
                    sol = least_squares(
                        residual,
                        theta0,
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
                sol = best
            elif best is not None:
                sol = best
            else:
                theta0 = np.zeros(len(unknowns))
                sol = least_squares(
                    residual,
                    theta0,
                    xtol=1e-12,
                    ftol=1e-12,
                    gtol=1e-12,
                    max_nfev=self.optcyc,
                    verbose=0,
                )

            if sol is None:
                raise RuntimeError(
                    f"MKM solver found no acceptable solution in {self.nguess} guesses."
                )

        # ------------------------------------------------------------------
        # PARALLEL PATH
        #
        # Rank 0 is only the monitor. Workers 1..size-1 solve independent
        # guesses continuously.
        #
        # HARD STOP INSIDE MKM:
        # On the one-node Slurm setup, rank 0 knows every worker PID. As soon
        # as one worker finds an acceptable solution below max_cost, rank 0
        # sends SIGUSR1 to every other worker. Their signal handler raises
        # _MKMStop out of the currently running least_squares() call. Thus the
        # expensive solver is interrupted immediately, but the MPI processes
        # themselves stay alive and the communicator remains valid.
        #
        # The winning OptimizeResult is sent directly to rank 0. Rank 0 waits
        # until every worker has left the search, then broadcasts exactly that
        # winning solution back to all ranks. No gather of all trial results is
        # needed after convergence.
        # ------------------------------------------------------------------
        else:
            TAG_BEST       = 301
            TAG_FOUND      = 302
            TAG_DONE       = 303
            TAG_RELEASE    = 304
            TAG_INITIAL_GO = 305

            nworkers = size - 1
            if nworkers < 1:
                raise RuntimeError("parallel=True requires at least 2 MPI ranks.")

            # All ranks run on one node in the intended Slurm setup. Collect
            # host/PID information once so rank 0 can interrupt workers without
            # polling MPI from inside scipy.
            process_info = comm.gather(
                (socket.gethostname(), os.getpid()),
                root=0,
            )

            if rank == 0:
                worker_hosts = {r: process_info[r][0] for r in range(1, size)}
                worker_pids  = {r: process_info[r][1] for r in range(1, size)}

                hosts = set(worker_hosts.values())
                if len(hosts) != 1:
                    raise RuntimeError(
                        "The signal-stop MPI implementation currently requires "
                        "all MKM ranks to run on one node."
                    )

                print(
                    f"MKM MPI: {size} ranks, {nworkers} solver workers",
                    flush=True,
                )

                # If an initial guess from a previous T/p point is available,
                # all workers evaluate one guess around it in parallel first.
                # No worker is allowed to trigger STOP during this first wave.
                initial_candidates = []

                if self.initial_guess is not None:
                    initial_results = comm.gather(None, root=0)

                    initial_candidates = [
                        candidate
                        for candidate in initial_results[1:]
                        if candidate is not None
                    ]

                    if initial_candidates:
                        initial_best = min(
                            initial_candidates,
                            key=lambda candidate: candidate.cost,
                        )
                        print(initial_best.cost, flush=True)
                    else:
                        initial_best = None

                    initial_converged = (
                        initial_best is not None
                        and initial_best.cost < self.max_cost
                    )

                    for worker in range(1, size):
                        comm.send(initial_converged, dest=worker, tag=TAG_INITIAL_GO)

                    if initial_converged:
                        # The full initial wave has finished, so this is the
                        # best acceptable continuation solution from all workers.
                        # No random search is needed.
                        sol = initial_best

                        for worker in range(1, size):
                            comm.send(None, dest=worker, tag=TAG_RELEASE)

                        print(
                            f"MKM final best cost: {sol.cost}",
                            flush=True,
                        )

                        sol = comm.bcast(sol, root=0)

                        if not acceptable(sol):
                            raise RuntimeError(
                                "Internal MKM error: selected MPI solution is not acceptable."
                            )

                        theta = dict(zip(unknowns, sol.x))

                        if self.cov_dep_eads:
                            residual(sol.x)

                        theta_star = get_theta_star(theta, sites)
                        return theta, theta_star, sol

                finished_workers = set()
                best_from_workers = {}
                global_best_cost = initial_best.cost if self.initial_guess is not None and initial_best is not None else np.inf
                winning_sol = None
                winning_rank = None
                stopping = False

                while len(finished_workers) < nworkers:
                    status = MPI.Status()
                    payload = comm.recv(
                        source=MPI.ANY_SOURCE,
                        tag=MPI.ANY_TAG,
                        status=status,
                    )
                    worker = status.Get_source()
                    tag = status.Get_tag()

                    if tag == TAG_BEST:
                        cost = float(payload)

                        if cost < global_best_cost:
                            global_best_cost = cost
                            print(global_best_cost, flush=True)

                    elif tag == TAG_FOUND and not stopping:
                        winning_rank = worker
                        winning_sol = payload
                        stopping = True

                        if winning_sol.cost < global_best_cost:
                            global_best_cost = float(winning_sol.cost)
                            print(global_best_cost, flush=True)

                        print(
                            f"MKM converged on rank {worker}: "
                            f"cost={winning_sol.cost} -> stopping all ranks",
                            flush=True,
                        )

                        # Interrupt every non-winning worker immediately inside
                        # its current least_squares() call.
                        for other, pid in worker_pids.items():
                            if other == winning_rank:
                                continue
                            try:
                                os.kill(pid, signal.SIGUSR1)
                            except ProcessLookupError:
                                pass

                    elif tag == TAG_DONE:
                        # Payload is the best acceptable local solution (or None).
                        finished_workers.add(worker)
                        best_from_workers[worker] = payload

                # Select the globally best acceptable solution known by ANY
                # worker, not simply the first rank that triggered STOP.
                #
                # winning_sol is the first solution below max_cost and is used
                # only to stop the search quickly. Every worker then reports
                # its local_best through TAG_DONE. A different worker may
                # already have found an even better acceptable solution before
                # receiving the stop signal.
                candidates = [
                    candidate
                    for candidate in best_from_workers.values()
                    if candidate is not None
                ]

                if winning_sol is not None:
                    candidates.append(winning_sol)

                if not candidates:
                    raise RuntimeError(
                        f"MKM solver found no acceptable solution in "
                        f"{self.nguess} guesses."
                    )

                sol = min(
                    candidates,
                    key=lambda candidate: candidate.cost,
                )

                # Nobody is still inside MKM search now. Release workers into
                # the common final broadcast.
                for worker in range(1, size):
                    comm.send(None, dest=worker, tag=TAG_RELEASE)

                print(
                    f"MKM final best cost: {sol.cost}",
                    flush=True,
                )

            else:
                class _MKMStop(Exception):
                    pass

                local_best = None
                local_best_cost = np.inf
                local_reported_cost = np.inf
                stop_requested = False
                inside_solver = False

                # First parallel continuation wave around self.initial_guess.
                # Rank 1 uses the exact previous theta; all other workers use
                # small multiplicative perturbations so that the neighbourhood
                # is explored in parallel instead of serializing onto rank 1.
                if self.initial_guess is not None:
                    theta_ig = np.array(
                        [self.initial_guess.get(s, 0) for s in unknowns],
                        dtype=float,
                    )

                    if rank != 1:
                        theta_ig = theta_ig * np.exp(
                            np.random.normal(0.0, 0.20, len(theta_ig))
                        )

                    initial_trial = least_squares(
                        residual,
                        theta_ig,
                        xtol=1e-12,
                        ftol=1e-12,
                        gtol=1e-12,
                        max_nfev=self.optcyc,
                        verbose=0,
                    )

                    if acceptable(initial_trial):
                        local_best = initial_trial
                        local_best_cost = initial_trial.cost
                        local_reported_cost = initial_trial.cost
                        initial_candidate = initial_trial
                    else:
                        initial_candidate = None

                    # Full first wave: every rank waits here until all workers
                    # have completed their continuation trial.
                    comm.gather(initial_candidate, root=0)

                    initial_converged = comm.recv(source=0, tag=TAG_INITIAL_GO)

                    if initial_converged:
                        # Rank 0 has already selected the globally best
                        # acceptable result from the entire first wave.
                        comm.recv(source=0, tag=TAG_RELEASE)
                        sol = None
                        sol = comm.bcast(sol, root=0)

                        if not acceptable(sol):
                            raise RuntimeError(
                                "Internal MKM error: selected MPI solution is not acceptable."
                            )

                        theta = dict(zip(unknowns, sol.x))

                        if self.cov_dep_eads:
                            residual(sol.x)

                        theta_star = get_theta_star(theta, sites)
                        return theta, theta_star, sol

                def _signal_stop(signum, frame):
                    nonlocal stop_requested, inside_solver
                    stop_requested = True
                    if inside_solver:
                        raise _MKMStop()

                old_signal_handler = signal.signal(signal.SIGUSR1, _signal_stop)

                def run_trial(theta0):
                    nonlocal inside_solver
                    inside_solver = True
                    try:
                        return least_squares(
                            residual,
                            theta0,
                            xtol=1e-12,
                            ftol=1e-12,
                            gtol=1e-12,
                            max_nfev=self.optcyc,
                            verbose=0,
                        )
                    except _MKMStop:
                        return None
                    finally:
                        inside_solver = False

                def process_trial(trial):
                    nonlocal local_best, local_best_cost
                    nonlocal local_reported_cost, stop_requested

                    if trial is None:
                        return False

                    trial_acceptable = acceptable(trial)

                    # Ignore physically invalid solutions completely.
                    # Only acceptable MKM solutions are allowed to update the
                    # reported/global best cost.
                    if trial_acceptable and trial.cost < local_best_cost:
                        local_best = trial
                        local_best_cost = trial.cost

                        if trial.cost < local_reported_cost:
                            local_reported_cost = trial.cost
                            comm.send(
                                float(trial.cost),
                                dest=0,
                                tag=TAG_BEST,
                            )

                        if trial.cost < self.max_cost:
                            # Send the complete converged solution to rank 0.
                            # Rank 0 will interrupt every other worker.
                            comm.send(trial, dest=0, tag=TAG_FOUND)
                            stop_requested = True
                            return True

                    return False

                try:
                    # Static distribution over workers. There is no scheduler
                    # and no MPI polling inside least_squares().
                    if not stop_requested and include_initial_guess is not None:
                        worker_index = rank - 1

                        for guess_index in range(
                            worker_index,
                            self.nguess,
                            nworkers,
                        ):
                            if stop_requested:
                                break

                            theta0 = make_random_guess()
                            #print( f"rank {rank}: starting guess {guess_index}",flush=True,)
                            trial = run_trial(theta0)
                            #print(f"rank {rank}: finished guess {guess_index}, "f"cost={None if trial is None else trial.cost}",flush=True,)

                            # SIGUSR1 may have interrupted this solve.
                            if stop_requested and trial is None:
                                break

                            if process_trial(trial):
                                break

                finally:
                    signal.signal(signal.SIGUSR1, old_signal_handler)

                # Worker is completely out of least_squares at this point.
                comm.send(local_best, dest=0, tag=TAG_DONE)

                # Rank 0 only releases workers after every worker has stopped.
                comm.recv(source=0, tag=TAG_RELEASE)
                sol = None

            # Every rank receives the exact globally best solution found before
            # termination. No rank leaves solve_mkm() while any worker is still
            # solving.
            sol = comm.bcast(sol, root=0)

            if not acceptable(sol):
                raise RuntimeError(
                    "Internal MKM error: selected MPI solution is not acceptable."
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
