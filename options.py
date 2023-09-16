import pydantic
from pydantic import BaseModel
from typing import Any, Dict, List, Literal, Optional, Union
from utils import canonicalize


class RetroBackendOption(BaseModel):
    retro_backend: str = "template_relevance"
    retro_model_name: str = "reaxys"
    max_num_templates: int = 100
    max_cum_prob: float = 0.995
    attribute_filter: List[Dict[str, Any]] = []


class ClusterSetting(BaseModel):
    feature: str = "original"
    cluster_method: str = "rxn_class"
    fp_type: str = "morgan"
    fp_length: int = 512
    fp_radius: int = 1
    classification_threshold: float = 0.2


class ExpandOneOptions(BaseModel):
    # v1 tree builder fields for backward compatibility
    template_count: int = 100
    max_cum_template_prob: float = 0.995
    forbidden_molecules: List[str] = []
    known_bad_reactions: List[str] = []

    # v2 tree builder fields
    template_max_count: int | None
    template_max_cum_prob: float | None
    banned_chemicals: List[str] | None
    banned_reactions: List[str] | None

    retro_backend_options: List[RetroBackendOption] = [RetroBackendOption()]
    use_fast_filter: bool = True
    filter_threshold: float = 0.75
    retro_rerank_backend: str | None = None
    cluster_precursors: bool = False
    cluster_setting: ClusterSetting = None
    extract_template: bool = False
    return_reacting_atoms: bool = True
    selectivity_check: bool = False

    @pydantic.validator("template_max_count", pre=True, always=True)
    def template_max_count_alias(cls, v: int, *, values: Dict[str, Any]) -> int:
        return v or values["template_count"]

    @pydantic.validator("template_max_cum_prob", pre=True, always=True)
    def template_max_cum_prob_alias(cls, v: int, *, values: Dict[str, Any]) -> float:
        return v or values["max_cum_template_prob"]

    @pydantic.validator("banned_chemicals", pre=True, always=True)
    def banned_chemicals_alias(cls, v: List[str], *, values: Dict[str, Any]) -> List[str]:
        banned_chemicals = v or values["forbidden_molecules"]
        banned_chemicals = list(set(canonicalize(smi) for smi in banned_chemicals))

        return banned_chemicals

    @pydantic.validator("banned_reactions", pre=True, always=True)
    def banned_reactions_alias(cls, v: List[str], *, values: Dict[str, Any]) -> List[str]:
        banned_reactions = v or values["known_bad_reactions"]
        banned_reactions = list(set(canonicalize(smi) for smi in banned_reactions))

        return banned_reactions


class BuildTreeOptions(BaseModel):
    expansion_time: int = 30
    max_iterations: Optional[int] = None
    max_chemicals: Optional[int] = None
    max_reactions: Optional[int] = None
    max_templates: Optional[int] = None
    max_branching: int = 25
    max_depth: int = 5
    exploration_weight: float = 1.0
    return_first: bool = False
    max_trees: Optional[int] = None
    max_ppg: Optional[float] = None
    max_scscore: Optional[float] = None
    max_elements: Optional[Dict[str, int]] = None
    min_history: Optional[Dict[str, int]] = None
    property_criteria: Optional[List[Dict[str, Any]]] = None
    termination_logic: Optional[Dict[str, List[str]]] = {"and": ["buyable"]}
    buyables_source: Optional[Union[str, List[str]]] = None
    custom_buyables: Optional[List[str]] = None


class EnumeratePathsOptions(BaseModel):
    path_format: Literal["json", "graph"] = "json"
    json_format: Literal["treedata", "nodelink"] = "treedata"
    sorting_metric: Literal[
        "plausibility",
        "number_of_starting_materials",
        "number_of_reactions",
        "score"
    ] = "plausibility"
    validate_paths: bool = True
    score_trees: bool = False
    cluster_trees: bool = False
    cluster_method: Literal["hdbscan", "kmeans"] = "hdbscan"
    min_samples: int = 5
    min_cluster_size: int = 5
    paths_only: bool = False
