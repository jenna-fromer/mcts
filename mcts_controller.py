import itertools
import networkx as nx
import numpy as np
import os
import time
from api.expand_one_api import ExpandOneAPI
from api.historian_api import HistorianAPI
from api.pathway_ranker_api import PathwayRankerAPI
from api.pricer_api import PricerAPI
from api.reaction_classification_api import ReactionClassificationAPI
from api.scscorer_api import SCScorerAPI
from options import ExpandOneOptions, BuildTreeOptions, EnumeratePathsOptions
from rdkit import Chem
from typing import List, Optional, Set, Tuple
from utils import get_graph_from_tree, is_terminal, nx_graph_to_paths, nx_paths_to_json

GATEWAY_URL = os.environ.get("GATEWAY_URL", "http://0.0.0.0:9100")
expand_one = ExpandOneAPI(
    default_url=f"{GATEWAY_URL}/api/tree_search/expand_one/call_sync_without_token",
)
historian = HistorianAPI(
    default_url=f"{GATEWAY_URL}/api/historian/lookup_smiles"
)
pathway_ranker = PathwayRankerAPI(
    url=f"{GATEWAY_URL}/api/pathway_ranker/call_sync"
)
pricer = PricerAPI(
    default_url=f"{GATEWAY_URL}/api/pricer/lookup_smiles"
)
reaction_classifier = ReactionClassificationAPI(
    url=f"{GATEWAY_URL}/api/get_top_class_batch/call_sync"
)
scscorer = SCScorerAPI(
    default_url=f"{GATEWAY_URL}/api/scscore/call_sync"
)


class MCTS:
    def __init__(self):
        self.expand_one_options = None
        self.build_tree_options = None
        self.enumerate_paths_options = None

        self.tree = nx.DiGraph()        # directed graph
        self.target = None              # the target compound
        self.target_uuid = None         # unique identifier for the target in paths
        self.paths = None               # pathway results as nx graphs
        self.chemicals = []
        self.reactions = []
        self.iterations = 0
        self.time_to_solve = 0

        self.expand_one = expand_one
        self.historian = historian
        self.pathway_ranker = pathway_ranker
        self.pricer = pricer
        self.scscorer = scscorer
        self.reaction_classifier = reaction_classifier

    @property
    def num_unique_chemicals(self):
        """Number of unique chemicals explored."""
        return len(self.chemicals)

    @property
    def num_unique_reactions(self):
        """Number of unique reactions explored."""
        return len(self.reactions)

    @property
    def num_total_reactions(self):
        """Total number of reactions explored."""
        return sum(self.tree.out_degree(chem) for chem in self.chemicals)

    @property
    def done(self):
        """Determine if we're done expanding the tree."""
        return (
            self.tree.nodes[self.target]["done"]
            or (
                self.build_tree_options.max_iterations is not None
                and self.iterations >= self.build_tree_options.max_iterations
            )
            or (
                self.build_tree_options.max_chemicals is not None
                and self.num_unique_chemicals >= self.build_tree_options.max_chemicals
            )
            or (
                self.build_tree_options.max_reactions is not None
                and self.num_unique_reactions >= self.build_tree_options.max_reactions
            )
            or (
                self.build_tree_options.max_templates is not None
                and self.num_total_reactions >= self.build_tree_options.max_templates
            )
        )

    def get_buyable_paths(
        self,
        target: str,
        expand_one_options: ExpandOneOptions = ExpandOneOptions(),
        build_tree_options: BuildTreeOptions = BuildTreeOptions(),
        enumerate_paths_options: EnumeratePathsOptions = EnumeratePathsOptions(),
    ) -> Tuple[List[dict], dict, dict]:
        """
        Build retrosynthesis tree and return paths to buyable precursors.
        *Based on v2 tree builder*
        Args:
            target (str): SMILES of target chemical
            expand_one_options (ExpandOneOptions object): options for one-step retro
            build_tree_options (BuildTreeOptions object): options for build_tree
            enumerate_paths_options (EnumeratePathsOptions object):
                options for enumerate_paths

        Returns:
            trees (list of dict): List of synthetic routes as networkx json
            stats (dict): Various statistics about the expansion
            graph (dict): Full explored graph as networkx node link json
        """
        self.expand_one_options = expand_one_options
        self.build_tree_options = build_tree_options
        self.enumerate_paths_options = enumerate_paths_options

        start = time.time()
        self.build_tree(target=target)
        build_time = time.time() - start

        start = time.time()
        paths = self.enumerate_paths()
        path_time = time.time() - start

        graph = nx.node_link_data(get_graph_from_tree(self.tree))
        stats = {
            "total_iterations": self.iterations,
            "total_chemicals": self.num_unique_chemicals,
            "total_reactions": self.num_unique_reactions,
            "total_templates": self.num_total_reactions,
            "total_paths": len(paths),
            "first_path_time": self.time_to_solve,
            "build_time": build_time,
            "path_time": path_time,
        }

        return paths, stats, graph

    def _initialize(self, target: str) -> None:
        """
        Initialize the tree by with the target chemical.
        """
        self.target = Chem.MolToSmiles(
            Chem.MolFromSmiles(target),
            isomericSmiles=True
        )           # Canonicalize SMILES
        self.create_chemical_node(smiles=self.target)
        self.tree.nodes[self.target]["terminal"] = False
        self.tree.nodes[self.target]["solved"] = False
        self.tree.nodes[self.target]["done"] = False

    def create_chemical_node(self, smiles: str) -> None:
        """
        Create a new chemical node from the provided SMILES and populate node
        properties with chemical data.

        Includes purchase price and *no* template info
        """
        purchase_price, properties = self.pricer(
            smiles=smiles,
            source=self.build_tree_options.buyable_source,
            canonicalize=False
        )

        template_sets = [option.retro_model_name for option
                         in self.expand_one_options.retro_backend_options]
        hist = self.historian(
            smiles=smiles,
            template_sets=template_sets,
            canonicalize=False
        )

        terminal = is_terminal(
            smiles=smiles,
            build_tree_options=self.build_tree_options,
            scscorer=self.scscorer,
            ppg=purchase_price,
            hist=hist,
            properties=properties
        )
        est_value = 1.0 if terminal else 0.0

        self.chemicals.append(smiles)
        # *terminal* is like a static property, e.g., buyable
        # *expanded* indicates whether _expand() has been called on this node
        # *solved* indicates whether this node falls on *a* buyable path
        # *done* is similar to "proven", indicating whether all subtrees have been expanded
        # Of course, a node can only be "done" if it has been "expanded",
        # unless it's "terminal", in which case it'd been "done" at creation
        self.tree.add_node(
            smiles,
            as_reactant=hist["as_reactant"],
            as_product=hist["as_product"],
            est_value=est_value,    # total value of node
            min_depth=None,         # minimum depth at which this chemical appears
            properties=properties,  # properties from buyables database if any
            purchase_price=purchase_price,
            solved=terminal,        # whether a path to terminal leaves has been found from this node
            terminal=terminal,      # whether this chemical meets terminal criteria
            done=terminal,          # simplified update logic from is_chemical_done()
            expanded=False,
            type="chemical",
            visit_count=1
        )

    def create_reaction_node(
        self,
        smiles: str,
        template_tuple: Optional[Tuple[str, str]],
        rxn_score_from_model: float,
        plausibility: float
    ):
        """Create a new reaction node from the provided smiles and data."""
        self.reactions.append(smiles)
        self.tree.add_node(
            smiles,
            est_value=0.0,      # score for how feasible a route is, based on whether precursors are terminal
            plausibility=plausibility,
            solved=False,       # whether a path to terminal leaves has been found from this node
            rxn_score_from_model=rxn_score_from_model,
            templates=[template_tuple] if template_tuple is not None else [],
            type="reaction",
            visit_count=1
        )

    def is_reaction_done(self, smiles: str):
        """
        Determine if the specified reaction node should be expanded further.

        Reaction nodes are done when all of its children chemicals are done.
        """
        return all(c["done"] for c in self.tree.successors(smiles))

    def build_tree(
        self,
        target: str
    ) -> None:
        """
        Build retrosynthesis tree by iterative expansion of precursor nodes.
        """
        print("Initializing tree...")
        self._initialize(target)

        print("Starting tree expansion...")
        start_time = time.time()
        elapsed_time = time.time() - start_time

        while elapsed_time < self.build_tree_options.expansion_time and not self.done:
            chem_path, rxn_path = self._select()
            self._expand(chem_path)
            self._update(chem_path, rxn_path)

            elapsed_time = time.time() - start_time

            self.iterations += 1
            if self.iterations % 100 == 0:
                print(f"Iteration {self.iterations} ({elapsed_time: .2f}s): "
                      f"|C| = {len(self.chemicals)} "
                      f"|R| = {len(self.reactions)}")

            if not self.time_to_solve and self.tree.nodes[self.target]["solved"]:
                self.time_to_solve = elapsed_time
                print(f"Found first pathway after {elapsed_time:.2f} seconds.")
                if self.build_tree_options.return_first:
                    print("Stopping expansion to return first pathway.")

        print("Tree expansion complete.")
        self.print_stats()

    def _select(self) -> Tuple[List[str], List[str]]:
        """
        Select next unexpanded leaf node to be expanded.

        This starts at the root node (target chemical), and at each level,
        use UCB to score each of the "reaction" options which can be taken.
        It will take the optimal option, which will now be an already explored
        reaction. It will descend to the next level and repeat the process until
        reaching an unexpanded node.
        """
        chem_path = [self.target]
        rxn_path = []
        invalid_options = set()

        while True:
            leaf = chem_path[-1]
            if not self.tree.nodes[leaf]["expanded"]:
                # Termination criteria for selection; found an unexpanded node
                # "done" nodes (including "terminal" ones) would have been excluded
                # when choosing the min precursor
                break

            options = self.ucb(
                node=leaf,
                chem_path=chem_path,
                invalid_options=invalid_options,
                exploration_weight=self.build_tree_options.exploration_weight
            )
            if not options:
                # There are no valid options from this chemical node, need to backtrack
                invalid_options.add(leaf)
                del chem_path[-1]
                del rxn_path[-1]
                continue

            score, reaction = options[0]
            # With ASKCOSv2 refactor, a reaction would always have been *explored*
            precursor = min(
                (
                    c for c in self.tree.successors(reaction)
                    if not c["done"] and c not in invalid_options
                ),
                key=lambda _node: self.tree.nodes[_node]["visit_count"],
                default=None
            )
            if precursor is None:
                # There are no valid options from this reaction node, need to backtrack
                invalid_options.add(reaction)
                continue
            else:
                chem_path.append(precursor)
                rxn_path.append(reaction)

        return chem_path, rxn_path

    def ucb(
        self,
        node: str,
        chem_path: List[str],
        invalid_options: Set[str],
        exploration_weight: float
    ) -> List[Tuple[float, str]]:
        """
        Calculate UCB score for all exploration options from the specified node. Only
        considers *explored reactions*, as reactions would have always been explored.

        Returns a list of (score, option) tuples sorted by score.
        """
        options = []
        product_visits = self.tree.nodes[node]["visit_count"]

        # Get scores for explored templates (reaction node exists)
        for rxn in self.tree.successors(node):
            rxn_data = self.tree.nodes[rxn]

            if (
                rxn in invalid_options
                # Simplified is_reaction_done
                or all(c["done"] for c in self.tree.successors(rxn))
                or len(set(self.tree.successors(rxn)) & set(chem_path)) > 0     # FIXME: why this condition
            ):
                continue

            est_value = rxn_data["est_value"]
            node_visits = rxn_data["visit_count"]
            # normalized_score is a generalized version of template score,
            # which would have been handled/computed by the retro controller
            rxn_probability = rxn_data["normalized_score"]

            # Q represents how good a move is
            q_sa = rxn_probability * est_value / node_visits
            # U represents how many times this move has been explored
            u_sa = np.sqrt(np.log(product_visits) / node_visits)

            score = q_sa + exploration_weight * u_sa

            # The options here are to follow a reaction down one level
            options.append((score, rxn))

        # Sort options from highest to lowest score
        options.sort(key=lambda x: x[0], reverse=True)

        return options

    def _expand(self, chem_path: List[str]) -> None:
        """
        Expand the tree by running one-step retro prediction to a chemical node
        """
        leaf = chem_path[-1]
        retro_results = self.expand_one(
            smiles=leaf,
            expand_one_options=self.expand_one_options
        )

        for result in retro_results:
            precursor_smiles = result["outcome"]
            reaction_smiles = precursor_smiles + ">>" + leaf
            reactant_list = precursor_smiles.split(".")

            try:
                template = result["template"]
                template_tuple = (template["index"], template["template_set"])
            except KeyError:
                # try-except just for backward compatibility
                template_tuple = None

            if reaction_smiles in self.reactions:
                # This Reaction node already exists
                rxn_data = self.tree.nodes[reaction_smiles]
                if template_tuple is not None:
                    rxn_data["templates"].append(template_tuple)

                # retro controller now computes normalized_model_score
                rxn_data["rxn_score_from_model"] = max(
                    rxn_data["rxn_score_from_model"], result["normalized_model_score"]
                )
            else:
                # This is new, so create a Reaction node
                self.create_reaction_node(
                    smiles=reaction_smiles,
                    template_tuple=template_tuple,
                    rxn_score_from_model=result["normalized_model_score"],
                    plausibility=result["plausibility"]
                )

            # Add edges to connect target -> reaction -> precursors
            self.tree.add_edge(leaf, reaction_smiles)
            for reactant in reactant_list:
                if reactant not in self.chemicals:
                    # This is new, so create a Chemical node
                    self.create_chemical_node(smiles=reactant)
                self.tree.add_edge(reaction_smiles, reactant)
            # This _update_value only updates reactions *below* leaf
            self._update_value(reaction_smiles)

        self.tree.nodes[leaf]["expanded"] = True

    def _update(self, chem_path: List[str], rxn_path: List[str]) -> None:
        """
        Update status and reward for nodes in this path.

        Reaction nodes are guaranteed to only have a single parent. Thus, the
        status of its parent chemical will always be updated appropriately in
        ``_update`` and will not change until the next time the chemical is
        in the selected path. Thus, the done state of the chemical can be saved.

        However, chemical nodes can have multiple parents (i.e. can be reached
        via multiple reactions), so a given update cycle may only pass through
        one of multiple parent reactions. Thus, the done state of a reaction
        must be determined dynamically and cannot be saved.
        """
        assert (
            chem_path[0] == self.target
        ), "Chemical path should start at the root node."

        # Iterate over the full path in reverse
        # On each iteration, rxn will be the parent reaction of chem
        # For the root (target) node, rxn will be None
        for i, chem, rxn in itertools.zip_longest(
            range(len(chem_path) - 1, -1, -1), reversed(chem_path), reversed(rxn_path)
        ):
            chem_data = self.tree.nodes[chem]
            chem_data["visit_count"] += 1
            chem_data["min_depth"] = (
                min(chem_data["min_depth"], i)
                if chem_data["min_depth"] is not None
                else i
            )
            # simplified update logic from is_chemical_done()
            if not chem_data["done"]:
                chem_data["done"] = (
                    chem_data["min_depth"] < self.build_tree_options.max_depth
                    and
                    chem_data["expanded"]
                    and
                    all(self.is_reaction_done(r) for r in self.tree.successors(chem))
                )

            if rxn is not None:
                rxn_data = self.tree.nodes[rxn]
                rxn_data["visit_count"] += 1
                # This _update_value only updates reactions *above* the expanded leaf
                self._update_value(rxn)

    def _update_value(self, smiles: str):
        """
        Update the value of the specified reaction node and its parent.
        """
        rxn_data = self.tree.nodes[smiles]

        if rxn_data["type"] == "reaction":
            # Calculate value as the sum of the values of all precursors
            est_value = sum(
                self.tree.nodes[c]["est_value"] for c in self.tree.successors(smiles)
            )

            # Update estimated value of reaction
            rxn_data["est_value"] += est_value

            # Update estimated value of parent chemical
            chem_data = self.tree.nodes[next(self.tree.predecessors(smiles))]
            chem_data["est_value"] += est_value

            # Check if this node is solved
            solved = rxn_data["solved"] or all(
                self.tree.nodes[c]["solved"] for c in self.tree.successors(smiles)
            )
            chem_data["solved"] = rxn_data["solved"] = solved

    def enumerate_paths(self) -> List:
        """
        Return list of paths to buyables starting from the target node.
        """
        if self.build_tree_options.return_first:
            self.enumerate_paths_options.score_trees = False
            self.enumerate_paths_options.cluster_trees = False

        self.paths, self.target_uuid = nx_graph_to_paths(
            self.tree,
            self.target,
            max_depth=self.build_tree_options.max_depth,
            max_trees=self.build_tree_options.max_trees,
            sorting_metric=self.enumerate_paths_options.sorting_metric,
            validate_paths=self.enumerate_paths_options.validate_paths,
            pathway_ranker=self.pathway_ranker,
            cluster_method=self.enumerate_paths_options.cluster_method,
            min_samples=self.enumerate_paths_options.min_samples,
            min_cluster_size=self.enumerate_paths_options.min_cluster_size
        )

        print(f"Found {len(self.paths)} paths to buyable chemicals.")

        path_format = self.enumerate_paths_options.path_format
        json_format = self.enumerate_paths_options.json_format

        if path_format == "graph":
            paths = self.paths
        elif path_format == "json":
            paths = nx_paths_to_json(
                paths=self.paths,
                root_uuid=self.target_uuid,
                json_format=json_format
            )
        else:
            raise ValueError(f"Unrecognized format type {path_format}")

        return paths

    def print_stats(self) -> None:
        """
        Print tree statistics.
        """
        info = "\n"
        info += f"Number of iterations: {self.iterations}\n"
        num_nodes = self.tree.number_of_nodes()
        info += f"Number of nodes: {num_nodes:d}\n"
        info += f"    Chemical nodes: {len(self.chemicals):d}\n"
        info += f"    Reaction nodes: {len(self.reactions):d}\n"
        info += f"Number of edges: {self.tree.number_of_edges():d}\n"
        if num_nodes > 0:
            info += f"Average in degree: " \
                    f"{sum(d for _, d in self.tree.in_degree()) / num_nodes:.4f}\n"
            info += f"Average out degree: " \
                    f"{sum(d for _, d in self.tree.out_degree()) / num_nodes:.4f}"

        print(info)
