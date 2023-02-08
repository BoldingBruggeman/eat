from typing import Any, List, MutableMapping
import datetime

import numpy as np

from .. import shared


class Log(shared.Plugin):
    def __init__(self, *variable_names, transform_obs: bool=True):
        self.variable_names = frozenset(variable_names)
        self.variable_metadata: List[Any] = []
        self.transform_obs = transform_obs

    def initialize(self, variables: MutableMapping[str, Any], *args, **kwargs):
        for name in self.variable_names:
            self.variable_metadata.append(variables[name])

    def before_analysis(
        self,
        time: datetime.datetime,
        state: np.ndarray,
        iobs: np.ndarray,
        obs: np.ndarray,
        *args,
        **kwargs
    ):
        for metadata in self.variable_metadata:
            affected_obs = (iobs >= metadata["start"]) & (iobs < metadata["stop"])
            if self.transform_obs and affected_obs.any():
                obs[affected_obs] = np.log10(obs[affected_obs])
            metadata["data"][...] = np.log10(metadata["data"])

    def after_analysis(self, *args, **kwargs):
        for metadata in self.variable_metadata:
            metadata["data"][...] = 10.0 ** metadata["data"]
